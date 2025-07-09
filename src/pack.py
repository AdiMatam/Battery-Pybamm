import json
import traceback
import pybamm
import numpy as np
from src.cell import Cell
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS, T, THEORETICAL_CAPACITY
import pandas as pd
import os
import pickle
import time

from src.variator import Variator
import concurrent.futures
import warnings
from enum import IntEnum, unique

@unique
class Protocol(IntEnum):
    CC_Discharge = 0
    CC_Charge = 1
    CV_Charge = 2
    Rest = 3

PROTOCOL_NAMES = {
    Protocol.CC_Discharge: "CC_Discharge",
    Protocol.CC_Charge   : "CC_Charge",
    Protocol.CV_Charge   : "CV_Charge",
    Protocol.Rest        : "Rest"
}

class Pack:
    def __init__(self, experiment: str, parallel, series, discrete_pts: int,
        model:pybamm.BaseModel, geo:dict, parameters:dict, aging: bool
    ):

        self.experiment = experiment
        if os.path.exists(f"data/{self.experiment}"):
            a = input("Experiment already exists. Data will be overwritten! 'Y' to proceed anyway: ")
            if (a != 'Y'):
                raise ValueError("Experiment already exists!")
        else:
            os.makedirs(f"data/{self.experiment}")

        self.data_path = f"data/{self.experiment}/data.csv"

        self.parallel = parallel
        self.series = series
        self.temperature = T

        self.model = model
        self.geo = geo
        self.parameters = parameters
        self.i_total = pybamm.Variable("Pack Current")

        self.iapps = [
            pybamm.Variable(f"String {i+1} Iapp") for i in range(parallel)
        ]

        self.aging = aging
        self.charging = pybamm.Parameter("Pack Charging?")
        self.p_aging = pybamm.Parameter("Aging On?")

        BIND_VALUES(parameters, 
            {
                self.charging: "[input]",
                self.p_aging: 1 if aging else 0
            }
        )

        self.shape = (series, parallel)

        cells = np.empty(self.shape, dtype=Cell)
        for i in range(series):
            for j in range(parallel):
                cells[i, j] = Cell(f"Cell {i + 1},{j + 1}", self.iapps[j], self.charging, self.p_aging, model, geo, parameters)


        self.cells = cells
        self.flat_cells = self.cells.flatten()

        self.discrete_pts = discrete_pts

        self.particles = [] 
        for cell in self.flat_cells:
            self.particles.append(cell.pos)
            self.particles.append(cell.neg)

        self.voltage = 0
        for i in range(self.series):
            self.voltage += self.cells[i, 0].vvolt

        self.outputs = ["Pack Voltage", "Pack Current"]
        for cell in self.flat_cells:
            for var in (cell.pos.c, cell.neg.c, cell.voltage, cell.neg.sei_L, cell.capacity):
                self.outputs.append(var.name)

        self.inps = {}
        for cell in self.flat_cells:
            binder = {param: param.value for param in (cell.pos.c0, cell.neg.c0, cell.pos.phi0, cell.neg.phi0, cell.neg.sei0)}
            BIND_VALUES(self.inps, binder)

        self.model.variables.update({
            "Pack Voltage": self.voltage,
            "Pack Current": self.i_total
        })


        self.reset()

    def simulate(self, protocol: Protocol, final_time: float, c_rate=None, iapp=None, until=None):
        assert(isinstance(protocol, Protocol))

        # TODO: determined by the input arguments to this function
        if c_rate and iapp:
            warnings.warn("Either c_rate and iapp should be given (but not both). C_rate argument will take precedence by default", UserWarning)
        elif (c_rate or iapp) and protocol == Protocol.CV_Charge:
            warnings.warn("For CV charge, c_rate and iapp will be ignored. Current voltage will be held until final time or until=<current_cut_ratio>, whichever is reached first", UserWarning)
        elif (not c_rate and not iapp) and (protocol == Protocol.CC_Charge or protocol == Protocol.CC_Discharge):
            raise ValueError("Please provide either a c_rate or iapp to perform a CC_charge/discharge protocol!")

        if until is None:
            warnings.warn("No stop condition given. Simulation will progress until the provided 'final_time'", UserWarning)
            until = 0.0

        iappt = -1.0
        if protocol == Protocol.CC_Charge or protocol == Protocol.CC_Discharge:
            print("[NOTICE]: 'until' parameter interpreted as voltage cutoff")
            if c_rate:
                iappt = THEORETICAL_CAPACITY * c_rate * self.parallel
            else:
                iappt = iapp

            self.latest_current = iappt

            if protocol == Protocol.CC_Discharge:
                self.model.algebraic.update({
                    self.i_total: (-iappt - self.i_total)
                })
                self.__setupDAE()
                BIND_VALUES(self.inps, { self.charging: 0 })

                self.model.initial_conditions.update({
                    self.i_total: -iappt
                })
                self.model.initial_conditions.update({
                    **{ self.iapps[i]: -iappt / self.parallel for i in range(self.parallel) },
                })

                self.model.events = [
                    pybamm.Event("Min Voltage Cutoff", (self.voltage - until)),
                ]

            else:
                self.model.algebraic.update({
                    self.i_total: (iappt - self.i_total)
                })
                self.__setupDAE()
                BIND_VALUES(self.inps, { self.charging: 1 })

                self.model.initial_conditions.update({
                    self.i_total: +iappt
                })
                self.model.initial_conditions.update({
                    **{ self.iapps[i]: +iappt / self.parallel for i in range(self.parallel) },
                })

                self.model.events = [
                    pybamm.Event("Max Voltage Cutoff", (until - self.voltage)),
                ]

        elif protocol == Protocol.CV_Charge:
            print("[NOTICE]: 'until' parameter interpreted as min current")
            self.model.algebraic.update({
                self.i_total: (self.latest_voltage - self.voltage)
            })
            self.__setupDAE()
            BIND_VALUES(self.inps, { self.charging: 1 })

            self.model.initial_conditions.update({
                self.i_total: self.latest_current
            })
            self.model.initial_conditions.update({
                **{ self.iapps[i]: self.latest_current / self.parallel for i in range(self.parallel) },
            })

            self.model.events = [
                pybamm.Event("Min Current Cutoff", pybamm.AbsoluteValue(self.i_total) - until),
            ]

        else:
            raise ValueError("Given protocol has not been implemented")


        self.param_ob = pybamm.ParameterValues(self.parameters)
        self.param_ob.process_model(self.model)
        self.param_ob.process_geometry(self.geo)

        mesh = pybamm.Mesh(self.geo, 
            { p.domain: pybamm.Uniform1DSubMesh for p in self.particles },
            { p.r: self.discrete_pts for p in self.particles }
        )

        self.disc = pybamm.Discretisation(mesh, 
            { p.domain: pybamm.FiniteVolume() for p in self.particles }
        )

        temp = self.disc.process_model(self.model, inplace=False)

        time_steps = np.linspace(0, final_time, 100)

        solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-5, root_tol=1e-10, dt_max=1e-10, root_method='lm', extra_options_setup={"max_num_steps": 100000})
        solution = solver.solve(temp, time_steps, inputs=self.inps)
        
        print(f"[TERMINATED BY]: {solution.termination}")

        # 1) Start a dict with the time‐vector
        data = {"Time": solution.t, "Clock": solution.t + self.prev_time}

        # 2) Loop over every cell and pull out the five variables you asked for
        for var in self.outputs:
            # use the variable’s .name as the column header
            entries = solution[var].entries
            if (len(entries.shape) == 2):
                data[var] = entries[-1]
            else:
                data[var] = entries

        # inject your cycle & protocol columns and set MultiIndex
        results_df = pd.DataFrame(data)

        results_df["#"] = self.cycle_number
        results_df["protocol"] = PROTOCOL_NAMES[protocol]
        results_df.set_index(["#", "protocol"], inplace=True)

        # self.protocol_chain.append(PROTOCOL_NAMES[protocol])
        # self.iapp_chain.append(iappt)
        # if ('final time' in solution.termination):
            # self.cutoff_chain.append(-1.0)
        # else:
            # self.cutoff_chain.append(until)
        # self.cycle_number_chain.append(self.cycle_number)

        # write the column header ONLY on first attempt
        write_header = not os.path.exists(self.data_path) or os.stat(self.data_path).st_size == 0
        results_df.to_csv(
            self.data_path,
            mode="a",
            header=write_header,
            index=True,
            index_label=["#", "protocol"],
        )

        self.__update_pack_state(solution)
        self.prev_time += solution.t[-1]

    def reset(self):
        if os.path.exists(self.data_path):
            open(self.data_path, 'w').close()

        self.prev_time = 0
        self.cycle_number = 1

    def next_cycle(self):
        self.cycle_number += 1

    def __update_pack_state(self, solution: pybamm.Solution):
        for cell in self.flat_cells:
            BIND_VALUES(self.inps, 
                {
                    cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                    cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                    cell.pos.phi0: solution[cell.pos.phi.name].entries[-1],
                    cell.neg.phi0: solution[cell.neg.phi.name].entries[-1],
                    cell.neg.sei0: solution[cell.neg.sei_L.name].entries[-1],
                }
            )

        self.latest_voltage = 0
        for i in range(self.series):
            c = self.cells[i, 0]
            self.latest_voltage += solution[c.voltage.name].entries[-1]
        

    def __setupDAE(self):

        # cutoffs[1] (max V-cut is effectively the vlock)
        # 'boolean algebra' to switch state from CC <-> CV
        self.model.algebraic.update({
            self.iapps[0]: self.i_total - sum(self.iapps),
        })

        for i in range(1, self.parallel):
            vbalance = 0
            ## V{str{n}} - V{str{n-1}} = 0 from n=[1, num-strings]
            for j in range(self.series):
                vbalance += self.cells[j, i].vvolt
                vbalance -= self.cells[j, i-1].vvolt

            #expr = cells[0, i].vvolt - cells[0, i-1].vvolt
            self.model.algebraic[self.iapps[i]] = vbalance #expr
    

    def export_profile(self):
        data = {
            'Experiment': self.experiment,
            'Parallel': self.parallel,
            'Series': self.series,
            'Temperature': self.temperature,
            'Aging': self.aging,

            'Protocols': self.protocol_chain,
            'Cutoffs': self.cutoff_chain,
            'I_apps': self.iapp_chain,
            'Cycles': self.cycle_number_chain,
        }

        data.update(Variator.JSON())

        file_path = f"data/{self.experiment}/profile.json"
        with open(file_path, 'w') as json_file:
            json.dump(data, json_file, indent=4)
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

class Pack:
    STATEMAP = {
        0: 'CC-discharge',
        1: 'CC-charge',
        2: 'CV-charge'
    }


    def __init__(self, experiment: str, parallel, series,
        model:pybamm.BaseModel, geo:dict, parameters:dict
    ):

        self.experiment = experiment
        if os.path.exists(f"data/{self.experiment}"):
            a = input("Experiment already exists. Data will be overwritten! 'Y' to proceed anyway: ")
            if (a != 'Y'):
                raise ValueError("Experiment already exists!")
        else:
            os.makedirs(f"data/{self.experiment}")

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

        self.cv_mode = pybamm.Parameter("CV Mode")
        self.cc_mode = pybamm.Negate(self.cv_mode - 1)
        self.charging = pybamm.Parameter("Pack Charging?")
        self.ilock = pybamm.Parameter("Current Lock")

        BIND_VALUES(parameters, 
            {
                self.ilock: "[input]",
                self.charging: "[input]",
                self.cv_mode: "[input]"
            }
        )

        self.shape = (series, parallel)

        cells = np.empty(self.shape, dtype=Cell)
        for i in range(series):
            for j in range(parallel):
                cells[i, j] = Cell(f"Cell {i + 1},{j + 1}", self.iapps[j], self.charging, model, geo, parameters)

        self.cells = cells


    def set_charge_protocol(self, cycles, crate_or_current, use_c_rate=True):
        self.cycles = cycles
        if use_c_rate:
            self.c_rate = crate_or_current
            self.iappt = THEORETICAL_CAPACITY * self.c_rate * self.parallel
        else:
            self.iappt = crate_or_current
            self.c_rate = self.iappt / (THEORETICAL_CAPACITY * self.parallel)

    def set_cutoffs(self, voltage_window: tuple, current_cut, capacity_cut):
        self.voltage_window = voltage_window
        self.current_cut = current_cut
        self.capacity_cut = capacity_cut

    
    # ------------

    def export_profile(self):
        data = {
            'Experiment': self.experiment,
            'Parallel': self.parallel,
            'Series': self.series,
            'Temperature': self.temperature,
            'Voltage Window': self.voltage_window,
            "C-rate": self.c_rate,
            'I-app': self.iappt,
            'I-app Cut Factor': self.current_cut,
            'Cycles': self.cycles,
        }

        data.update(Variator.JSON())

        file_path = f"data/{self.experiment}/profile.json"
        with open(file_path, 'w') as json_file:
            json.dump(data, json_file, indent=4)

    def build(self, discrete_pts):
        
        self.__setupDAE()
        self.__IC_and_StopC()

        self.param_ob = pybamm.ParameterValues(self.parameters)
        self.param_ob.process_model(self.model)
        self.param_ob.process_geometry(self.geo)

        self.model.variables.update({
            "Pack Voltage": self.voltage,
            "Pack Current": self.i_total
        })
        SET_MODEL_VARS(self.model, self.iapps)


        particles = [] 
        for cell in self.flat_cells:
            particles.append(cell.pos)
            particles.append(cell.neg)

        mesh = pybamm.Mesh(self.geo, 
            { p.domain: pybamm.Uniform1DSubMesh for p in particles },
            { p.r: discrete_pts for p in particles }
        )

        disc = pybamm.Discretisation(mesh, 
            { p.domain: pybamm.FiniteVolume() for p in particles }
        )
        disc.process_model(self.model)

    def cycler(self, hours, time_pts):
        solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-5, root_tol=1e-10, dt_max=1e-10, root_method='lm', extra_options_setup={"max_num_steps": 100000}, return_solution_if_failed_early=True)
        time_steps = np.linspace(0, 3600 * hours, time_pts)
        
        inps = {}
        self.capacity_dict = {}
        csv_cols= self.__setup_initialization(inps)

        cycle_storage = {col: [] for col in csv_cols}
        pd.DataFrame(columns=csv_cols)                         .to_csv(f"data/{self.experiment}/data.csv", index=False)
        pd.DataFrame(columns=list(self.capacity_dict.keys()))  .to_csv(f"data/{self.experiment}/capacities.csv", index=False)

        prev_time = 0
        state = 0
        i = 0

        def cycle_dump(storage: dict, i: int, state: int):
            subdf = pd.DataFrame(storage)
            subdf = pd.concat({(i+1, Pack.STATEMAP[state]): subdf})
            print(subdf.iloc[:, :2])

            subdf = subdf.rename_axis(index=["Cycle", "Protocol", "Stamps"])            
            subdf.to_csv(f"data/{self.experiment}/data.csv", mode='a', header=False, index=True)

        def cap_dump(i: int):
            with open(f"data/{self.experiment}/capacities.csv", mode='a') as f:
                f.write(str(i+1))
                for vals in self.capacity_dict.values():
                    f.write(f",{str(vals[1])}")
                f.write('\n')

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            try:
                while i < self.cycles:
                    solution = solver.solve(self.model, time_steps, inputs=inps)

                    print(f"{i}) {solution.termination}")
                    print(f"Completed cycle {i+1}, {Pack.STATEMAP[state]}")                

                    cycle_storage['Time'] = solution.t
                    cycle_storage['Global Time'] = solution.t + prev_time
                    prev_time += solution.t[-1]

                    ## KEYS ARE SOLVED VARIABLES
                    for col in csv_cols:
                        data = solution[col].entries
                        if len(data.shape) == 2:
                            data = data[-1]
                        cycle_storage[col].extend(data)

                    ## 1) set initial conditions for the next cycle (with 'last' data from this cycle)
                    ## 2) Store discharge capacity in sep capacity_dict
                    terminate = self.__update_pack_state(inps, solution, state, i)
                    
                    executor.submit(cycle_dump, cycle_storage, i, state)
                    if (state == 0):
                        executor.submit(cap_dump, i)

                    state = self.__next_protocol(inps, state)

                    if (terminate):
                        break

                    cycle_storage = {col: [] for col in csv_cols}

                    if (state == 0):
                        i += 1
                    
            except Exception as e:
                print(traceback.format_exc())
                print (f"FAILED AT CYCLE # {i+1}. Dumping collected data so far")
                self.cycles = i

            finally:
                with open(f"data/{self.experiment}/model.pkl", 'wb') as f:
                    pickle.dump(self, f)


    def __setup_initialization(self, inps: dict):
        outputs = ["Pack Voltage", "Pack Current"]
        SET_OUTPUTS(outputs, self.iapps)

        BIND_VALUES(inps, 
            {
                self.ilock: -self.iappt,
                self.cv_mode: 0,
                self.charging: 0,
            }
        )
        
        for cell in self.flat_cells:
            SET_OUTPUTS(outputs, [cell.pos.c, cell.neg.c, cell.sei, cell.voltage, cell.capacity])
            BIND_VALUES(inps, 
                {
                    cell.pos.c0: cell.pos.c0.value,
                    cell.neg.c0: cell.neg.c0.value,
                    cell.pos.phi0: cell.pos.phi0.value,
                    cell.neg.phi0: cell.neg.phi0.value,
                    cell.neg.sei0: 5.e-9,
                }
            )
            self.capacity_dict[cell.name] = [0, 0]

        return outputs


    def __update_pack_state(self, inps: dict, solution: pybamm.Solution, state: int, cycle_num: int):
        a = 0
        for cell in self.flat_cells:
            BIND_VALUES(inps, 
                {
                    cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                    cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                    cell.pos.phi0: solution[cell.pos.phi.name].entries[-1],
                    cell.neg.phi0: solution[cell.neg.phi.name].entries[-1],
                    cell.neg.sei0: solution[cell.neg.sei_L.name].entries[-1],
                }
            )
            if (state == 0):
                self.capacity_dict[cell.name][1] = solution[cell.capacity.name].entries[-1] 
                if (cycle_num == 1):
                    ## store reference capacity at index 0
                    self.capacity_dict[cell.name][0] = solution[cell.capacity.name].entries[-1] 

                elif (cycle_num > 1):
                    ## degradation check
                    cur = self.capacity_dict[cell.name][1]
                    orig = self.capacity_dict[cell.name][0]

                    if (cur <= orig*self.capacity_cut):
                        print(f"Capacity cut triggered: {cell.name}")
                        a = 1

        return a
        
    def __next_protocol(self, inps: dict, state: int):
        # CC charge up next
        if (state == 0):
            BIND_VALUES(inps, 
                {
                    self.ilock: +self.iappt,
                    self.charging: 1,
                    self.cv_mode: 0 
                }
            )

        # CV charge up next
        elif (state == 1):
            BIND_VALUES(inps, 
                {
                    self.charging: 1,
                    self.cv_mode: 1 
                }
            )

        # Discharge next
        else:
            BIND_VALUES(inps, 
                {
                    self.ilock: -self.iappt,
                    self.charging: 0,
                    self.cv_mode: 0 
                }

            )
    
        nstate = (state + 1) % 3

        return nstate



    def __setupDAE(self):
        self.voltage = 0
        for i in range(self.series):
            self.voltage += self.cells[i, 0].vvolt

        # cutoffs[1] (max V-cut is effectively the vlock)
        self.model.algebraic.update({
            self.i_total: (self.ilock - self.i_total)*self.cc_mode + 
            (self.voltage_window[1] - self.voltage)*self.cv_mode
        })

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
    

    def __IC_and_StopC(self):
        self.model.initial_conditions.update({
            self.i_total: self.ilock
        })

        self.model.initial_conditions.update({
            **{ self.iapps[i]: self.ilock / self.parallel for i in range(self.parallel) },
        })

        self.flat_cells = self.cells.flatten()

        min_current = self.iappt * self.current_cut

        self.model.events += [
            pybamm.Event("Min Voltage Cutoff", self.voltage - self.voltage_window[0]),
            pybamm.Event("Max Voltage Cutoff", (self.voltage_window[1] - self.voltage)*self.cc_mode + 1*self.cv_mode),
            pybamm.Event("Min Current Cutoff", pybamm.AbsoluteValue(self.i_total) - min_current),
        ]

        # for cell in self.flat_cells:
            # self.model.events.extend([
                # pybamm.Event(f"{cell.name} Min Anode Concentration Cutoff", cell.neg.surf_c - 10),
                # pybamm.Event(f"{cell.name} Max Cathode Concentration Cutoff", cell.pos.cmax - cell.pos.surf_c),

                # pybamm.Event(f"{cell.name} Max Anode Concentration Cutoff", cell.neg.cmax - cell.neg.surf_c),
            # ])



if __name__ == '__main__':

    NUM_PARALLEL = 2
    NUM_SERIES = 2
    BASE_CURRENT = 13.6319183090575
    I_INPUT = BASE_CURRENT * NUM_PARALLEL

    VOLTAGE_LOW_CUT = 2.0
    VOLTAGE_HIGH_CUT =4.1

    #--------------------
    model = pybamm.BaseModel()
    geo = {}
    parameters = {}

    pack = Pack(NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), I_INPUT / 10, model, geo, parameters)

    print(pack[1, 1])
import json
import traceback
import pybamm
import numpy as np
from cell import Cell
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS, T, THEORETICAL_CAPACITY
import pandas as pd
import os
import pickle
import time

from variator import Variator
import concurrent.futures

class Pack:
    STATEMAP = {
        0: 'CC-discharge',
        1: 'CC-charge',
        2: 'CV-charge'
    }

    def __init__(self, experiment: str, parallel, series):

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

        # self.model = model
        # self.geo = geo
        # self.parameters = parameters
        self.i_total = pybamm.Variable("Pack Current")

        self.iapps = [
            pybamm.Variable(f"String {i+1} Iapp") for i in range(parallel)
        ]

        # self.cv_mode = pybamm.Parameter("CV Mode")
        # self.cc_mode = pybamm.Negate(self.cv_mode - 1)
        self.charging = pybamm.Parameter("Pack Charging?")
        self.discharging = pybamm.Negate(self.charging - 1)
        self.ilock = pybamm.Parameter("Current Lock")


        self.shape = (self.series, self.parallel)


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

    def export_profile(self, i):
        data = {
            'Experiment': self.experiment,
            'Parallel': self.parallel,
            'Series': self.series,
            'Temperature': self.temperature,
            'Voltage Window': self.voltage_window,
            "C-rate": self.c_rate,
            'I-app': self.iappt,
            'I-app Cut Factor': self.current_cut,
            'Capacity Cut Factor': self.capacity_cut,
            'Cycles': f"{i}/{self.cycles}",
        }

        data.update(Variator.JSON())

        file_path = f"data/{self.experiment}/profile.json"
        with open(file_path, 'w') as json_file:
            json.dump(data, json_file, indent=4)

    def build(self, discrete_pts):
        
        # self.__setupDAE(cc=True)
        # self.__IC_and_StopC()

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

        self.solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-5, root_tol=1e-8, dt_max=1e-12, max_step_decrease_count=15,
                    root_method='casadi')
        # self.solver =pybamm.IDAKLUSolver(root_tol=1e-10)

    def cycler(self, hours, time_pts):
        time_steps = np.linspace(0, 3600 * hours, time_pts)
        
        inps = {}
        outputs = self.__setup_initialization_and_outputs(inps)

        ## insert at front
        cycle_columns = ['Time', 'Global Time'] + outputs
        cycle_data = {col: [] for col in cycle_columns}

        self.__create_dataframe_files(cycle_columns, ["Pack Capacity"] + [c.name for c in self.cells[0]])

        prev_time = 0
        prev_cycle_time = 0
        state = 0
        i = 0

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            try:
                futures = []
                while i < self.cycles:
                    if i == 1:
                        time_steps = time_steps[1:]

                    solution = self.solver.solve(self.model, time_steps, inputs=inps)

                    print(f"Completed cycle {i+1}, {Pack.STATEMAP[state]} -- HIT {solution.termination}")                

                    cycle_data['Time'] = solution.t + prev_cycle_time
                    cycle_data['Global Time'] = solution.t + prev_time
                    prev_cycle_time += solution.t[-1]
                    prev_time += solution.t[-1]

                    ## KEYS ARE SOLVED VARIABLES
                    for var in outputs:
                        data = solution[var].entries
                        if len(data.shape) == 2:
                            data = data[-1]
                        cycle_data[var].extend(data)
                    
                    if len(futures) != 0: 
                        concurrent.futures.wait(futures)
                        futures.clear()

                    ## 1) set initial conditions for the next cycle (with 'last' data from this cycle)
                    ## 2) Store discharge capacity in sep capacity_dict
                    capcut = self.__update_pack_state(inps, solution, i, state)
                    
                    futures.append(executor.submit(self.__cycle_dump, cycle_data, i, state))
                    if (state == 0):
                        futures.append(executor.submit(self.__cap_dump, i))

                    cycle_data = {col: [] for col in cycle_columns}

                    state = self.__next_protocol(inps, state)
                    if (state == 0):
                        prev_cycle_time = 0
                        i += 1
                    
            except Exception as e:
                print(traceback.format_exc())
                print (f"FAILED AT CYCLE # {i+1}. Dumping collected data so far")

            finally:
                self.cycles = i
                self.export_profile(i)
                with open(f"data/{self.experiment}/model.pkl", 'wb') as f:
                    pickle.dump(self, f)


    def __cycle_dump(self, data: dict, i: int, state: int):
        subdf = pd.DataFrame(data)
        subdf = pd.concat({(i+1, Pack.STATEMAP[state]): subdf})
        subdf.to_csv(f"data/{self.experiment}/data.csv", mode='a', header=False, index=True)

    def __cap_dump(self, i: int):
        pass
        # with open(f"data/{self.experiment}/capacities.csv", mode='a') as f:
        #     f.write(str(i+1))
        #     # f.write(f",{self.capacity_value}")

        #     ## Only need to look at the first row of cells (first cell in each parallel branch)
        #     ## Each cell in a branch will have the same 'real' capacity (same current integrated over time)
        #     for cell in self.cells[0]:
        #         f.write(f",{cell.capacity_value}")

        #     f.write('\n')

    def __create_dataframe_files(self, cycle_columns, cell_names):
        pd.DataFrame(
            columns=cycle_columns, 
            index=pd.MultiIndex.from_product([[], [], []], names=["Cycle", "Protocol", "Stamps"])
        ).to_csv(f"data/{self.experiment}/data.csv", index=True)

        # pd.DataFrame(
        #     columns=cell_names,
        #     index=pd.MultiIndex.from_product([[]], names=["Cycle"])
        # ).to_csv(f"data/{self.experiment}/capacities.csv", index=True)
    

    def __reconstruct(self):
        self.model = pybamm.BaseModel()
        self.geo = {}
        self.parameters = {}

        BIND_VALUES(self.parameters, 
            {
                self.ilock: "[input]",
                self.charging: "[input]",
                # self.cv_mode: "[input]"
            }
        )

        cells = np.empty(self.shape, dtype=Cell)
        for i in range(self.series):
            for j in range(self.parallel):
                cells[i, j] = Cell(f"Cell {i + 1},{j + 1}", self.iapps[j], self.charging, self.model, self.geo, self.parameters)

        self.cells = cells
        self.flat_cells = self.cells.flatten()

    def __setup_initialization_and_outputs(self, inps: dict):

        self.__reconstruct()

        outputs = ["Pack Voltage", "Pack Current"]
        SET_OUTPUTS(outputs, self.iapps)

        BIND_VALUES(inps, 
            {
                self.ilock: -self.iappt,
                # self.cv_mode: 0,
                self.charging: 0,
            }
        )
        
        for cell in self.flat_cells:
            # SET_OUTPUTS(outputs, [cell.pos.c, cell.neg.c, cell.sei, cell.voltage, cell.neg.i_sei, cell.neg.i_int])
            # SET_OUTPUTS(outputs, [cell.pos.c, cell.neg.c, cell.voltage, cell.capacity])
            SET_OUTPUTS(outputs, [cell.pos.c, cell.neg.c, cell.sei, cell.voltage, cell.capacity])
            BIND_VALUES(inps, 
                {
                    cell.pos.c0: cell.pos.c0.value,
                    cell.neg.c0: cell.neg.c0.value,
                    cell.neg.sei0: cell.neg.sei0.value, # 0 #5.e-9,
                }
            )

        self.__setupDAE(cc=True)
        self.__IC_and_StopC()
        self.build(100)

        return outputs


    def __update_pack_state(self, inps: dict, solution: pybamm.Solution, i: int, state: int):
        for cell in self.flat_cells:
            BIND_VALUES(inps, 
                {
                    cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                    cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                    cell.neg.sei0: solution[cell.sei.name].entries[-1],
                }
            )

        return False
        
    def __next_protocol(self, inps: dict, state: int):
        # CC charge up next
        if (state == 0):
            BIND_VALUES(inps, 
                {
                    self.ilock: +self.iappt,
                    self.charging: 1,
                    # self.cv_mode: 0 
                }
            )
            # REDISCRETIZATION
            self.__reconstruct()
            self.__setupDAE(cc=True)
            self.__IC_and_StopC()
            self.build(100)


        # CV charge up next
        elif (state == 1):
            BIND_VALUES(inps, 
                {
                    self.ilock: self.iappt,
                    self.charging: 1,
                    # self.cv_mode: 0
                }
            )
            self.__reconstruct()
            self.__setupDAE(cc=False)
            self.__IC_and_StopC()
            self.build(100)

        # Discharge next
        else:
            BIND_VALUES(inps, 
                {
                    self.ilock: -self.iappt,
                    self.charging: 0,
                    # self.cv_mode: 0 
                }
            )
            # REDISCRETIZATION
            self.__reconstruct()
            self.__setupDAE(cc=True)
            self.__IC_and_StopC()
            self.build(100)

        nstate = (state + 1) % 3

        return nstate



    def __setupDAE(self, cc=True):
        self.voltage = 0
        for i in range(self.series):
            self.voltage += self.cells[i, 0].vvolt

        # cutoffs[1] (max V-cut is effectively the vlock)
        # 'boolean algebra' to switch state from CC <-> CV
        if cc:
            self.model.algebraic.update({
                self.i_total: self.ilock - self.i_total
            })
        else:
            self.model.algebraic.update({
                self.i_total: (self.voltage_window[1] - self.voltage)
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

        if cc:
            self.model.events = [
                pybamm.Event("Min Voltage Cutoff", (self.voltage - self.voltage_window[0])*self.discharging + 1*self.charging),
                pybamm.Event("Max Voltage Cutoff", (self.voltage_window[1] - self.voltage)*self.charging + 1*self.discharging),
            ]
        else:
            min_current = self.iappt * self.current_cut
            self.model.events = [
                pybamm.Event("Min Current Cutoff", (pybamm.AbsoluteValue(self.i_total) - min_current)*self.charging + 1*self.discharging),
            ]  
    
    def __IC_and_StopC(self):
        self.model.initial_conditions.update({
            self.i_total: self.ilock
        })

        self.model.initial_conditions.update({
            **{ self.iapps[i]: self.ilock / self.parallel for i in range(self.parallel) },
        })

        # self.model.events += [
        #     pybamm.Event("Min Voltage Cutoff", (self.voltage - self.voltage_window[0])*self.discharging + 1*self.charging),
        #     pybamm.Event("Max Voltage Cutoff", (self.voltage_window[1] - self.voltage)*self.charging + 1*self.discharging),
        #     # pybamm.Event("Min Current Cutoff", (pybamm.AbsoluteValue(self.i_total) - min_current)*self.charging + 1*self.discharging),
        # ]


if __name__ == '__main__':
    pass

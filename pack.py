import pybamm
import numpy as np
from cell import Cell
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS
import params as p
import pandas as pd
from typing import List

class Pack:
    def __init__(self,i_total: pybamm.Parameter, parallel: int, series: int, cutoffs: tuple, 
        model:pybamm.BaseModel, geo:dict, parameters:dict
    ):

        self.parallel = parallel
        self.series = series

        self.model = model
        self.geo = geo
        self.parameters = parameters
        self.i_total = i_total

        self.iapps = [
            pybamm.Variable(f"String {i+1} Iapp") for i in range(parallel)
        ]


        size = (series, parallel)
        cells = np.empty(size, dtype=Cell)
        for i in range(series):
            for j in range(parallel):
                cells[i, j] = Cell(f"Cell {i + 1},{j + 1}", self.iapps[j], model, geo, parameters)
                # c = cells[i,j]
                # model.algebraic[c.volt] = c.pos.phi - c.neg.phi - c.volt
        #c[0, 1] - c[0, 0] + c[1, 1] - c[1, 0]

        self.cells = cells

        # TODO: temporary (for parallel-only pack)
        self.voltage = cells[0, 0].vvolt

        model.algebraic.update({
            self.iapps[0]: i_total - sum(self.iapps[1:]) - self.iapps[0],
        })

        for i in range(1, parallel):
            #vbalance = 0
            #for j in range(series):
                #vbalance += cells[j, i].volt
                #vbalance -= cells[j, i-1].volt

            expr = cells[0, i].vvolt - cells[0, i-1].vvolt
            model.algebraic[self.iapps[i]] = expr
    

        model.initial_conditions.update({
            **{ self.iapps[i]: -i_total / parallel for i in range(parallel) },
        })

        self.flat_cells = self.cells.flatten()

        model.events += [
            pybamm.Event("Min Voltage Cutoff", self.voltage - cutoffs[0]*series),
            pybamm.Event("Max Voltage Cutoff", cutoffs[1]*series - self.voltage),
        ]

        for cell in self.flat_cells:
            model.events.extend([
                pybamm.Event(f"{cell.name} Min Anode Concentration Cutoff", cell.neg.surf_c - 10),
                pybamm.Event(f"{cell.name} Max Cathode Concentration Cutoff", cell.pos.cmax - cell.pos.surf_c),

                pybamm.Event(f"{cell.name} Max Anode Concentration Cutoff", cell.neg.cmax - cell.neg.surf_c),
            ])

        self.param_ob = pybamm.ParameterValues(parameters)
        self.param_ob.process_model(model)
        self.param_ob.process_geometry(geo)


        SET_MODEL_VARS(model,
            self.iapps
        )

        
        # for i in range(parallel):
            # model.variables[self.iapps[i].name] = self.iapps[i]
            # model.variables[self.volts[i].name] = 0

            # for j in range(series):
                # c = self.cells[j, i]
                # model.variables[self.volts[i].name] += c.volt
                # model.variables[c.volt.name] = c.volt

                # if i == 0:
                    # # capture symbolic sum of one of the strings
                    # self.voltsum += c.volt


        # self.capacity = 0 
        # for i in range(self.parallel):
            # cap = min(map(lambda c: c.capacity, self.cells[:, i]))
            # self.capacity += cap
    
    def build(self, discrete_pts):
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

    def cycler(self, i_input, cycles, hours, time_pts, output_path=""):
        solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, dt_max=1e-10, extra_options_setup={"max_num_steps": 100000})

        time_steps = np.linspace(0, 3600 * hours, time_pts)
        
        sign = -1
        inps = {}
        outputs = SET_OUTPUTS(self.iapps)

        BIND_VALUES(inps, 
            {
                self.i_total: sign * i_input,
            }
        )
        for cell in self.flat_cells:
            outputs.extend(
                SET_OUTPUTS([cell.pos.c, cell.neg.c, cell.neg.sei_L, cell.voltage])
            )
            BIND_VALUES(inps, 
                {
                    cell.pos.c0: cell.GET[cell.pos.c0.name],
                    cell.neg.c0: cell.GET[cell.neg.c0.name],
                    cell.neg.sei0: 5.e-9,
                    cell.neg.iflag: 0 if (sign == -1) else 1
                }
            )

        ### EVERYTHING BELOW THIS IS JUST RUNNING / CAPTURING SIMULATION DATA.
        ### NO PARAMETER-RELEVANT CODE BELOW

        caps = []
        subdfs = []

        solution = None
        prev_time = 0

        for i in range(cycles):

            solution = solver.solve(self.model, time_steps, inputs=inps)

            subdf = pd.DataFrame(columns=['Global Time', 'Time'] + outputs)
            subdf['Time'] = solution.t
            subdf['Global Time'] = solution.t + prev_time
            prev_time += solution.t[-1]

            if (i % 2 == 0):
                caps.append( i_input * solution.t[-1] / 3600 )

            ## KEYS ARE VARIABLES
            for key in outputs:
                data = solution[key].entries
                if len(data.shape) == 2:
                    # TODO: see if there's a fix for this
                    data = data[-1] # last node (all nodes 'equal' due to broadcasted surface concentration)

                subdf[key] = data

            subdf = pd.concat({f'C{i+1}': subdf})
            subdfs.append(subdf)
            print(f"Finished Cycle #{i}")

            sign *= -1

            BIND_VALUES(inps, 
                {
                    self.i_total: sign * i_input,
                }
            )
            for cell in self.flat_cells:
                BIND_VALUES(inps, 
                    {
                        cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                        cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                        cell.neg.sei0: solution[cell.neg.sei_L.name].entries[-1],
                        cell.neg.iflag: 0 if (sign == -1) else 1
                    }
                )

        df = pd.concat(subdfs)
        
        print(df)
        print(caps)

        df.to_csv(output_path)
        return df, caps


"""

        solver = pybamm.CasadiSolver()
        seconds = 3600 * runtime_hours
        time_steps = np.linspace(0, seconds, time_pts)

        inps = {
            self.i_total.name: -iapp,
            ## intial conditions for particle concentrations are set here!!
            ** {cell.pos.c0.name: cell.pos_csn_ival for cell in self.flat_cells},
            ** {cell.neg.c0.name: cell.neg_csn_ival for cell in self.flat_cells},
        }

        caps = []
        subdfs = []

        solution = None
        prev_time = 0
        for _ in range(cycles):
            solution = solver.solve(self.model, time_steps, inputs=inps)

            subdf = pd.DataFrame(columns=['Time'] + variables)
            subdf['Time'] = solution.t + prev_time
            prev_time += solution.t[-1]

            caps.append( iapp * solution.t[-1] / 3600 )

            ## KEYS ARE VARIABLES
            for key in variables:
                data = solution[key].entries
                if len(data.shape) == 2:
                    data = data[-1] # last node (all nodes 'equal' due to broadcast)

                subdf[key] = data

            subdfs.append(subdf)

            inps[self.i_total.name] *= -1
            for cell in self.flat_cells:
                inps.update({
                    cell.pos.c0.name: solution[cell.pos.csn.name].entries[-1][-1],
                    cell.neg.c0.name: solution[cell.neg.csn.name].entries[-1][-1]
                })

        df = pd.concat(subdfs, ignore_index=True)
        if output_path:
            df.to_csv(output_path, index=False)
        return df, caps

"""


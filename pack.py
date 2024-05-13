import pybamm
import numpy as np
from cell import Cell
import params as p
import pandas as pd
from typing import List

class Pack:
    def __init__(self, parallel: int, series: int, model:pybamm.BaseModel, geo:dict, parameters:dict, 
                    i_total: pybamm.Parameter
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
        self.volts = [
            pybamm.Variable(f"String {i+1} Voltage") for i in range(parallel)
        ]
        self.voltsum = 0

        size = (series, parallel)
        cells = np.empty(size, dtype=Cell)
        for i in range(series):
            for j in range(parallel):
                cells[i, j] = Cell(f"Cell {i + 1},{j + 1}", model, geo, parameters, iapp=self.iapps[j])
                c = cells[i,j]
                model.algebraic[c.volt] = c.pos.phi - c.neg.phi - c.volt
        #c[0, 1] - c[0, 0] + c[1, 1] - c[1, 0]

        self.cells = cells

        model.algebraic.update({
            self.iapps[0]: i_total - sum(self.iapps[1:]) - self.iapps[0],
        })
        for i in range(1, parallel):
            vbalance = 0
            for j in range(series):
                vbalance += cells[j, i].volt
                vbalance -= cells[j, i-1].volt

            model.algebraic[self.iapps[i]] = vbalance
        
        for i in range(parallel):
            model.variables[self.iapps[i].name] = self.iapps[i]
            model.variables[self.volts[i].name] = 0

            for j in range(series):
                c = self.cells[j, i]
                model.variables[self.volts[i].name] += c.volt
                model.variables[c.volt.name] = c.volt

                if i == 0:
                    # capture symbolic sum of one of the strings
                    self.voltsum += c.volt

        self.flat_cells = self.cells.flatten()

        self.param_ob = pybamm.ParameterValues(parameters)
        self.param_ob.process_model(model)
        self.param_ob.process_geometry(geo)

        self.capacity = 0 #np.inf
        for i in range(self.parallel):
            cap = min(map(lambda c: c.capacity, self.cells[:, i]))
            self.capacity += cap
    
    def get_capacity(self):
        return self.capacity

    def get_cells(self) -> List[Cell]:
        return self.cells

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

    def cycler(self, iapp, cycles, runtime_hours, time_pts, variables, output_path=""):
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


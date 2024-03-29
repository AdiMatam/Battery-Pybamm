import pybamm
import numpy as np
from cell import Cell
import params as p
import pandas as pd
from typing import List

class Pack:
    def __init__(self, num_cells:int, model:pybamm.BaseModel, geo:dict, parameters:dict, 
                    i_param:pybamm.Parameter, voltage_cutoff: tuple
        ):
        self.num_cells = num_cells
        self.model = model
        self.geo = geo
        self.parameters = parameters
        self.i_total = i_param

        size = (2, 2)
        cells = np.empty(size, dtype=Cell)
        for i in range(2):
            for j in range(2):
                cells[i, j] = Cell(f"Cell {i + 1}{j + 1}", model, geo, parameters, voltage_cutoff)

        self.cells = cells
        

        str_voltages = [0 for _ in range(size[1])] # np.empty((size[1], 0), dtype=pybamm.Variable)
        for j in range(size[1]):
            str_voltages[j] = pybamm.Variable(f"String Voltage {j}")

        for j in range(size[1]-1):        
            model.algebraic[str_voltages[j]] = str_voltages[j+1] - str_voltages[j]

        last_col_v = 0
        for i in range(size[0]):
            last_col_v += cells[i, -1].pos.phival - cells[i, -1].neg.phival

        model.algebraic[str_voltages[-1]] = (last_col_v) - str_voltages[-1]


        for j in range(size[1]):
            v = 0
            for i in range(size[0]-1):
                model.algebraic[cells[i, j].iapp] = cells[i+1, j].iapp - cells[i, j].iapp
                v += (cells[i, j].pos.phival - cells[i, j].neg.phival)

            model.algebraic[cells[-1, j].iapp] = v


        i_branches = 0
        for j in range(size[1]):
            i_branches += cells[0, j].iapp

        model.algebraic[cells[-1, -1].iapp] = i_param - i_branches


        model.variables.update({
            **{ str_voltages[i].name: str_voltages[i] for i in range(len(str_voltages)) }
        })

        model.initial_conditions.update({
            **{ str_voltages[i]: p.POS_OCP(cells[0,i].pos_csn_ival / cells[0,i].pos_csn_maxval) 
                                - p.NEG_OCP(cells[-1,i].neg_csn_ival / cells[-1,i].pos_csn_maxval) 
            
                for i in range(size[1])
            }, 
        })


        self.flat_cells = self.cells.flatten()
        self.capacity = self.num_cells * min(cell.capacity for cell in self.flat_cells) # this is in Ah/m^2

        self.param_ob = pybamm.ParameterValues(parameters)
        self.param_ob.process_model(model)
        self.param_ob.process_geometry(geo)


    def get_num_cells(self):
        return self.num_cells

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
        cycles = cycles
        solver = pybamm.CasadiSolver()
        time_steps = np.linspace(0, 3600 * runtime_hours, time_pts)

        inps = {
            self.i_total.name: -iapp,
            ## intial conditions for particle concentrations are set here!!
            ** {cell.pos.c_0.name: cell.pos_csn_ival for cell in self.flat_cells},
            ** {cell.neg.c_0.name: cell.neg_csn_ival for cell in self.flat_cells},
        }

        subdfs = []

        solution = None
        for _ in range(cycles):
            solution = solver.solve(self.model, time_steps, inputs=inps)

            subdf = pd.DataFrame(columns=['Time'] + variables)
            subdf['Time'] = list(solution.t)

            ## KEYS ARE VARIABLES
            for key in variables:
                data = solution[key].entries
                if "Concentration" in key:
                    data = data[-1] # last node (all nodes 'equal' due to broadcast)

                subdf[key] = data

            subdfs.append(subdf)

            inps[self.i_total.name] *= -1
            for cell in self.flat_cells:
                inps.update({
                    cell.pos.c_0_name: solution[cell.pos.surf_csn_name].entries[-1][-1],
                    cell.neg.c_0_name: solution[cell.neg.surf_csn_name].entries[-1][-1]
                })

        df = pd.concat(subdfs, ignore_index=True)
        if output_path:
            df.to_csv(output_path, index=False)
        return df
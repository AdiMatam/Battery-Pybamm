import pybamm
import numpy as np
from cell import Cell
import params as p
import pandas as pd
from typing import List

class Pack:
    def __init__(self, num_cells:int, model:pybamm.BaseModel, geo:dict, parameters:dict, i_total:pybamm.Parameter):
        self.num_cells = num_cells
        self.model = model
        self.geo = geo
        self.parameters = parameters
        self.i_total = i_total

        cells = [Cell(f"Cell {i + 1}", model, geo, parameters) for i in range(num_cells)]
        self.cells = cells

        # Vcell1 - Vcell2 = 0
        # (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
        for i in range(len(cells) - 1):
            model.algebraic[cells[i].pos.phi] = cells[i + 1].pos.phi - cells[i].pos.phi
            model.algebraic[cells[i].neg.phi] = cells[i + 1].neg.phi - cells[i].neg.phi

            model.algebraic[cells[i].iapp] = cells[i].pos.bv_term - cells[i].neg.bv_term - (cells[i].pos.phi - cells[i].neg.phi)

        model.algebraic[cells[-1].pos.phi] = cells[-1].pos.bv_term - cells[-1].pos.phi
        model.algebraic[cells[-1].neg.phi] = cells[-1].neg.bv_term - cells[-1].neg.phi
        model.algebraic[cells[-1].iapp] = i_total - sum(cell.iapp for cell in cells)


        self.capacity = self.num_cells * min(cell.capacity for cell in cells) # this is in Ah/m^2

        self.param_ob = pybamm.ParameterValues(parameters)
        self.param_ob.process_model(model)
        self.param_ob.process_geometry(geo)


    def get_num_cells(self):
        return self.num_cells

    def get_cells(self) -> List[Cell]:
        return self.cells

    def build(self, discrete_pts):
        particles = [] 
        for cell in self.cells:
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
        total_time_steps = np.linspace(0, 3600 * runtime_hours * cycles, time_pts * cycles)

        inps = {
            self.i_total.name: -iapp,
            ## intial conditions for particle concentrations are set here!!
            ** {cell.pos.c_0.name: cell.pos_csn_ival for cell in self.cells},
            ** {cell.neg.c_0.name: cell.neg_csn_ival for cell in self.cells},
        }

        subdfs = []

        solution = None
        for i in range(cycles):
            solution = solver.solve(self.model, time_steps, inputs=inps)

            subdf = pd.DataFrame(columns=variables)

            ## KEYS ARE VARIABLES
            for key in variables:
                data = solution[key].entries
                if "Concentration" in key:
                    data = data[-1] # last node (all nodes 'equal' due to broadcast)

                subdf[key] = data

            subdfs.append(subdf)

            inps[self.i_total.name] *= -1
            for cell in self.cells:
                inps.update({
                    cell.pos.c_0_name: solution[cell.pos.surf_csn_name].entries[-1][-1],
                    cell.neg.c_0_name: solution[cell.neg.surf_csn_name].entries[-1][-1]
                })

        df = pd.concat(subdfs, ignore_index=True)
        df.insert(0, 'Time', total_time_steps)
        if output_path:
            df.to_csv(output_path, index=False)
        return df
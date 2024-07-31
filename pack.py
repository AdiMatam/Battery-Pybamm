import pybamm
import numpy as np
from cell import Cell
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS
import pandas as pd

class Pack:
    def __init__(self, parallel: int, series: int, cutoffs: tuple, min_current: float, 
        model:pybamm.BaseModel, geo:dict, parameters:dict
    ):

        self.parallel = parallel
        self.series = series
        self.cutoffs = cutoffs
        self.min_current = min_current

        self.model = model
        self.geo = geo
        self.parameters = parameters
        self.i_total = pybamm.Variable("Solved Total Current")

        self.iapps = [
            pybamm.Variable(f"String {i+1} Iapp") for i in range(parallel)
        ]

        self.cv_mode = pybamm.Parameter("CV Mode")
        self.cc_mode = pybamm.Negate(pybamm.Subtraction(self.cv_mode, 1))
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
                c = cells[i,j]
                model.algebraic[c.voltage] = c.pos.phi - c.neg.phi - c.voltage
        #c[0, 1] - c[0, 0] + c[1, 1] - c[1, 0]

        self.cells = cells

        self.voltage = 0
        for i in range(series):
            self.voltage += cells[i, 0].vvolt

        # cutoffs[1] (max V-cut is effectively the vlock)
        model.algebraic.update({
            self.i_total: (self.ilock - self.i_total)*self.cc_mode + (cutoffs[1]*series - self.voltage)*self.cv_mode
        })

        model.algebraic.update({
            self.iapps[0]: self.i_total - sum(self.iapps),
        })

        for i in range(1, parallel):
            vbalance = 0
            for j in range(series):
                vbalance += cells[j, i].vvolt
                vbalance -= cells[j, i-1].vvolt

            #expr = cells[0, i].vvolt - cells[0, i-1].vvolt
            model.algebraic[self.iapps[i]] = vbalance #expr
    

        model.initial_conditions.update({
            self.i_total: self.ilock
        })

        model.initial_conditions.update({
            **{ self.iapps[i]: self.ilock / parallel for i in range(parallel) },
        })

        self.flat_cells = self.cells.flatten()

        model.events += [
            pybamm.Event("Min Voltage Cutoff", self.voltage - cutoffs[0]*series),
            pybamm.Event("Max Voltage Cutoff", (cutoffs[1]*series - self.voltage)*self.cc_mode + 1*self.cv_mode),
            pybamm.Event("Min Current Cutoff", pybamm.AbsoluteValue(self.i_total) - min_current*parallel),
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

        self.volt_var = pybamm.Variable("Pack Voltage")
        model.variables.update({
            self.volt_var.name: self.voltage
        })
        SET_MODEL_VARS(model,
            self.iapps
        )
    
    def __getitem__(self, index):
        if (index[0] < 0 or index[0] >= self.series or index[1] < 0 or index[1] >= self.series):
            raise IndexError(f"{index} invalid for {self.shape} Pack")

        return self.cells[index]

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

    def cycler(self, iappt, cycles, hours, time_pts, output_path=""):
        self.cycles = cycles
        self.iappt = iappt

        solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-5, dt_max=1e-10, extra_options_setup={"max_num_steps": 100000})

        time_steps = np.linspace(0, 3600 * hours, time_pts)
        
        inps = {}
        outputs = [self.volt_var.name]
        SET_OUTPUTS(outputs, self.iapps)

        BIND_VALUES(inps, 
            {
                self.ilock: -iappt,
                self.cv_mode: 0,
                self.charging: 0,
            }
        )
        for cell in self.flat_cells:
            SET_OUTPUTS(outputs, [cell.pos.c, cell.neg.c, cell.neg.sei_L, cell.voltage, cell.capacity])

            BIND_VALUES(inps, 
                {
                    cell.pos.c0: cell.pos.c0.value,
                    cell.neg.c0: cell.neg.c0.value,
                    cell.neg.sei0: 5.e-9,
                }
            )

        ### EVERYTHING BELOW THIS IS JUST RUNNING / CAPTURING SIMULATION DATA.
        ### NO PARAMETER-RELEVANT CODE BELOW

        subdfs = []

        solution = None
        prev_time = 0

        state = 0
        i = 0
        while i < cycles:
            try:
                solution = solver.solve(self.model, time_steps, inputs=inps)
            except Exception as e:
                print(e)
                print (f"FAILED AT CYCLE # {i}. Dumping collected data so far")
                break

            subdf = pd.DataFrame(columns=['Global Time', 'Time'] + outputs)
            subdf['Time'] = solution.t
            subdf['Global Time'] = solution.t + prev_time
            prev_time += solution.t[-1]

            ## KEYS ARE VARIABLES
            for key in outputs:
                data = solution[key].entries
                if len(data.shape) == 2:
                    # TODO: see if there's a fix for this
                    data = data[-1] # last node (all nodes 'equal' due to broadcasted surface concentration)

                subdf[key] = data

            for cell in self.flat_cells:
                BIND_VALUES(inps, 
                    {
                        cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                        cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                        cell.neg.sei0: solution[cell.neg.sei_L.name].entries[-1],
                    }
                )

            just_finished = 0
            
            # CC charge up next
            if (state == 0):
                just_finished = "CC-discharge"
                BIND_VALUES(inps, 
                    {
                        self.ilock: +iappt,
                        self.charging: 1,
                        self.cv_mode: 0 
                    }
                )

            # CV charge up next
            elif (state == 1):
                just_finished = "CC-charge"
                BIND_VALUES(inps, 
                    {
                        self.charging: 1,
                        self.cv_mode: 1 
                    }
                )

            # Discharge next
            else:
                just_finished = "CV-charge"
                BIND_VALUES(inps, 
                    {
                        self.ilock: -iappt,
                        self.charging: 0,
                        self.cv_mode: 0 
                    }
                )

            print(f"Finished Cycle #{i} -- {just_finished}")

            subdf = pd.concat({(i+1, just_finished): subdf})
            subdfs.append(subdf)
            state = (state + 1) % 3
            if (state == 0):
                i += 1
            

        df = pd.concat(subdfs)
        df.to_csv(output_path)
        return df


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
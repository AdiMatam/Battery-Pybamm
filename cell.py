import pybamm
from particle_anode import Anode
from particle_cathode import Cathode
import consts as c
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS
import params as p

class Cell:
    CELLS = list()
    def __init__(self, name: str,iapp: pybamm.Variable, model: pybamm.BaseModel, geo:dict, parameters:dict):
        if name in self.CELLS:
            raise ValueError("Must have unique cell names/IDs")

        self.CELLS.append(name)
        self.name = name
        self.model = model

        self.iapp = iapp
        
        ## JUST NAMESAKE (not used)
        self.voltage = pybamm.Variable(name + " Voltage")

        self.pos = Cathode(name + " Cathode", iapp)
        self.neg = Anode(name + " Anode", iapp)

        self.vvolt = self.pos.phi - self.neg.phi

        self.pos.process_model(model)
        self.neg.process_model(model)

        self.GET = {}
        self.__attach_parameters(parameters)
        
        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)

        self.capacity = self.compute_capacity(self.GET)

        model.variables.update({
            self.voltage.name: self.vvolt
        })

    def compute_capacity(self, G: dict):
        # (csn_max - csn_min) * L * (1-eps_n) * (mols->Ah)
        # cathode
        pos_cap = (G[self.pos.cmax.name] - G[self.pos.c0.name])
        pos_cap *= G[self.pos.L.name] * (1-G[self.pos.eps_n.name]) * (c.F / 3600)

        # anode
        neg_cap = (G[self.neg.c0.name])
        neg_cap *= G[self.neg.L.name] * (1-G[self.neg.eps_n.name]) * (c.F / 3600)

        # TODO: anode is coming limited... wrong?

        return min(pos_cap, neg_cap)

    def __attach_parameters(self, param_dict: dict):

        BIND_VALUES(param_dict, {
            self.pos.c0:               "[input]",
            self.pos.L:                p.POS_ELEC_THICKNESS.rand_sample(),
            self.pos.eps_n:            p.POS_ELEC_POROSITY.rand_sample(),
            self.pos.cmax:             p.POS_CSN_MAX.rand_sample(),

            self.pos.ocp:              p.POS_OCP,
            self.pos.D:                p.POS_DIFFUSION.rand_sample(),
            self.pos.R:                p.PARTICLE_RADIUS.rand_sample(),
        })

        BIND_VALUES(param_dict, {
            self.neg.c0:               "[input]",
            self.neg.L:                p.NEG_ELEC_THICKNESS.rand_sample(),
            self.neg.eps_n:            p.NEG_ELEC_POROSITY.rand_sample(),
            self.neg.cmax:             p.NEG_CSN_MAX.rand_sample(),

            self.neg.ocp:              p.NEG_OCP,
            self.neg.D:                p.NEG_DIFFUSION.rand_sample(),
            self.neg.R:                p.PARTICLE_RADIUS.rand_sample(),
            self.neg.sei0:             "[input]",
            self.neg.iflag:            "[input]",
        })

        self.GET = param_dict.copy()
        BIND_VALUES(self.GET, {
            self.pos.c0: p.POS_CSN_INITIAL.rand_sample(),
            self.neg.c0: p.NEG_CSN_INITIAL.rand_sample(),
        })

if __name__ == '__main__':
    import pandas as pd
    import numpy as np

    HOURS = 2 
    I_INPUT = 13.6319183090575 #2.4
    DISCRETE_PTS = 30
    TIME_PTS = 100

    geo = {}
    parameters = {}
    model = pybamm.BaseModel()
    iapp = pybamm.Parameter("Input Current")

    parameters[iapp.name] = "[input]"

    cell = Cell("Cell1", iapp, model, geo, parameters)
    print(cell.capacity / 2)

    model.events += [
        pybamm.Event("Min Anode Concentration Cutoff", cell.neg.surf_c - 10),
        pybamm.Event("Max Cathode Concentration Cutoff", cell.pos.cmax - cell.pos.surf_c),

        pybamm.Event("Max Anode Concentration Cutoff", cell.neg.cmax - cell.neg.surf_c),

        pybamm.Event("Min Voltage Cutoff", (cell.pos.phi - cell.neg.phi) - 2.0),
        pybamm.Event("Max Voltage Cutoff", 4.1 - (cell.pos.phi - cell.neg.phi)),
    ]

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    particles = [cell.pos, cell.neg]
    mesh = pybamm.Mesh(geo, 
        { d.domain: pybamm.Uniform1DSubMesh for d in particles },
        { d.r: DISCRETE_PTS for d in particles }
    )

    disc = pybamm.Discretisation(mesh, 
        { d.domain: pybamm.FiniteVolume() for d in particles }
    )

    disc.process_model(model)

    CYCLES = 10
    solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, dt_max=1e-10, extra_options_setup={"max_num_steps": 100000})

    time_steps = np.linspace(0, 3600 * HOURS, TIME_PTS)
    
    sign = -1
    inps = {}
    BIND_VALUES(inps, 
        {
            iapp: sign * I_INPUT,
            cell.pos.c0: p.POS_CSN_INITIAL.get_value(),

            cell.neg.c0: p.NEG_CSN_INITIAL.get_value(),
            cell.neg.sei0: 5.e-9,
            cell.neg.iflag: 0 if (sign == -1) else 1
        }
    )

    ### EVERYTHING BELOW THIS IS JUST RUNNING / CAPTURING SIMULATION DATA.
    ### NO PARAMETER-RELEVANT CODE BELOW

    outputs = SET_OUTPUTS(
        [cell.pos.c, cell.pos.phi, cell.neg.c, cell.neg.phi, cell.neg.sei_L]
    )
    caps = []
    subdfs = []

    solution = None
    prev_time = 0

    for i in range(CYCLES):

        solution = solver.solve(model, time_steps, inputs=inps)

        subdf = pd.DataFrame(columns=['Time'] + outputs)
        subdf['Time'] = solution.t #+ prev_time
        #prev_time += solution.t[-1]

        caps.append( I_INPUT * solution.t[-1] / 3600 )

        ## KEYS ARE VARIABLES
        for key in outputs:
            data = solution[key].entries
            if len(data.shape) == 2:
                data = data[-1] # last node (all nodes 'equal' due to broadcasted surface concentration)

            subdf[key] = data

        subdf = pd.concat({f'C{i+1}': subdf})
        subdfs.append(subdf)
        print(f"Finished Cycle #{i}")

        sign *= -1
        BIND_VALUES(inps, 
            {
                iapp: sign * I_INPUT,
                cell.pos.c0: solution[cell.pos.c.name].entries[-1][-1],
                cell.neg.c0: solution[cell.neg.c.name].entries[-1][-1],
                cell.neg.sei0: solution[cell.neg.sei_L.name].entries[-1],
                cell.neg.iflag: 0 if (sign == -1) else 1
            }
        )

    df = pd.concat(subdfs)
    
    print(df)
    print(caps)

    df.to_csv(f"CELL.csv")

    # from zvalidate import plotter
    # plotter(df)
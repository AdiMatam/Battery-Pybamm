import pybamm
import consts as c
from params import POS_OCP
import params as p
from single_particle import SingleParticle
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS

class Cathode(SingleParticle):
    OCP_INIT = 4.234963004675769

    def __init__(self, name: str, 
            iapp: pybamm.Variable):

        super().__init__(name, +1, iapp)
    
    def process_model(self, model: pybamm.BaseModel):
        # dc/dt = d^2c/dr^2
        flux = self.D * -pybamm.grad(self.c)
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.c: dcdt,
        })

        KINT = 1.04e-11
        x = c.F / (2 * c.R_GAS * c.T) * (self.phi - self.ocp)
        j0 = c.F * KINT * self.surf_c**0.5 * (self.cmax - self.surf_c)**0.5 

        model.algebraic.update({
            self.phi: j0 * 2 * pybamm.sinh(x) - self.j,
        })

        model.initial_conditions.update({
            self.c: self.c0,
            self.phi: self.OCP_INIT
        }) 

        model.boundary_conditions.update({
            self.c: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * self.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        # model.variables.update{}
        model.variables.update({
            self.c.name: self.c, # pybamm.PrimaryBroadcast(self.surf_c, self.domain),
            self.phi.name: self.phi
        })

    def attach_parameters(self, parameters: dict):
        BIND_VALUES(parameters, {
            self.c0:               "[input]",
            self.L:                p.POS_ELEC_THICKNESS.sample(),
            self.eps_n:            p.POS_ELEC_POROSITY.sample(),
            self.cmax:             p.POS_CSN_MAX.sample(),

            self.ocp:              p.POS_OCP,
            self.D:                p.POS_DIFFUSION.sample(),
            self.R:                p.PARTICLE_RADIUS.sample(),
        })

        self.c0.set_value(p.POS_CSN_INITIAL.sample()) 

if __name__ == '__main__':
    import params as p
    import numpy as np
    import pandas as pd

    C_RATE = 1.0
    I_INPUT = 27.263836618115 * C_RATE
    HOURS = (1./C_RATE)
    DISCRETE_PTS = 100
    TIME_PTS = 100

    model = pybamm.BaseModel()
    iapp = pybamm.Parameter("Input Current") 
    geo = {}
    parameters = {}

    cath = Cathode("Cathode", iapp)
    cath.process_model(model)
    cath.process_geometry(geo)
    cath.attach_parameters(parameters)

    BIND_VALUES(parameters, {
        iapp: "[input]",
    })

    # model.events += [
    #     pybamm.Event("Min Concentration", cath.surf_c - 500),
    #     pybamm.Event("Max Concentration", cath.cmax.value + 100 - cath.surf_c)
    # ]

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    particles = [cath]
    mesh = pybamm.Mesh(geo, 
        { d.domain: pybamm.Uniform1DSubMesh for d in particles },
        { d.r: DISCRETE_PTS for d in particles }
    )

    disc = pybamm.Discretisation(mesh, 
        { d.domain: pybamm.FiniteVolume() for d in particles }
    )

    disc.process_model(model)

    cycles = 1
    solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, extra_options_setup={"max_num_steps": 100000})

    time_steps = np.linspace(0, 3600 * HOURS, TIME_PTS)
    total_time_steps = np.linspace(0, 3600 * HOURS * cycles, TIME_PTS * cycles)

    inps = {
        iapp.name: -1 * I_INPUT,
        cath.c0.name: p.POS_CSN_INITIAL.sample(),
    }

    solution = solver.solve(model, time_steps, inputs=inps)
    # solution.plot([cath.c.name])

    subdf = pd.DataFrame(columns=['Time', 'Cathode Concentration', 'Cathode Potential'])
    subdf['Time'] = solution.t
    subdf['Cathode Concentration'] = solution[cath.c.name].entries[-1]
    subdf['Cathode Potential'] = solution[cath.phi.name].entries

    t = []
    conc = []
    pot = []

    with open("sundata.txt") as f:
        for line in f:
            data = line.split("|")[:-1]
            t.append(float(data[0].strip()))
            conc.append(float(data[1].strip()))
            pot.append(float(data[3].strip()))

    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(2)
    
    ax[0].plot(solution.t, subdf['Cathode Concentration'], label="Pybamm")
    ax[0].plot(t, conc, label="Sundials")
    ax[1].plot(solution.t, subdf['Cathode Potential'], label="Pybamm")
    ax[1].plot(t, pot, label="Pybamm")
    ax[0].legend()
    ax[1].legend()
    
    plt.show()
    
    # sign = -1

    # ### EVERYTHING BELOW THIS IS JUST RUNNING / CAPTURING SIMULATION DATA.
    # ### NO PARAMETER-RELEVANT CODE BELOW

    # outputs = SET_OUTPUTS([cc.c, cc.phi])
    # caps = []
    # subdfs = []

    # solution = None
    # prev_time = 0

    # for _ in range(cycles):

        # solution = solver.solve(model, time_steps, inputs=inps)

        # subdf = pd.DataFrame(columns=['Time'] + outputs)
        # subdf['Time'] = solution.t + prev_time
        # prev_time += solution.t[-1]

        # caps.append( I_INPUT * solution.t[-1] / 3600 )

        # ## KEYS ARE VARIABLES
        # for key in outputs:
            # data = solution[key].entries
            # if len(data.shape) == 2:
                # data = data[-1] # last node (all nodes 'equal' due to broadcasted surface concentration)

            # subdf[key] = data

        # subdfs.append(subdf)

        # sign *= -1
        # BIND_VALUES(inps, 
            # {
                # iapp: sign * I_INPUT,
                # cc.c0: solution[cc.c.name].entries[-1][-1],
            # }
        # )

    # df = pd.concat(subdfs, ignore_index=True)
    
    # print(df)
    # print(caps)

    # df.to_csv(f"CATHODE_{cycles}.csv", index=False)
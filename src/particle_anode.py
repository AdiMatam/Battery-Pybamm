import pybamm
import consts as cc
from consts import SET_MODEL_VARS, SET_OUTPUTS, BIND_VALUES
from params import NEG_OCP
from src.single_particle import SingleParticle
import params as p

#pybamm.set_logging_level("DEBUG")

class Anode(SingleParticle): 
    OCP_INIT = 0.08352811644995728

    def __init__(self, name: str, iapp: pybamm.Variable):

        super().__init__(name, -1, iapp)

        self.i_sei = pybamm.Variable(name + " Side Current")
        self.i_int = pybamm.Variable(name + " Intercalation Current")
        self.sei_L = pybamm.Variable(name + " SEI Length")
        self.sei0 = pybamm.Parameter(name + " Initial SEI Length")

    def process_model(self, model: pybamm.BaseModel, charging):
        flux = self.D * -pybamm.grad(self.c)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

        ## see params.py
        # self.eps_n <- NEG_ELEC_POROSITY
        # self.L     <- NEG_ELEC_THICKNESS

        KSEI = 5.0e-6
        M_SEI = 0.162
        RHO_SEI = 1690
        KINT = 2.07e-11

        ## -- SEI START -- 
        dLdt = (-self.i_sei / (2*cc.F)) * (M_SEI / RHO_SEI)

        # solve the ODEs -- diffusion equation (del * del(c))
        model.rhs.update({
            self.c: dcdt,
            self.sei_L: dLdt
        })

        ## anode (SEI case)
        x = cc.F / (2 * cc.R_GAS * cc.T) * (self.phi - self.ocp - (self.sei_L/KSEI)*self.j)

        ## SEE PAPER
        kfs = 1.36e-12 #* 10
        cec_init = 0.05 * 4541
        is_rhs = charging * -cc.F*kfs*cec_init * pybamm.exp( (-0.5*cc.F)/(cc.R_GAS*cc.T) * (self.phi - (self.sei_L/KSEI)*self.j) ) 

        j0 = cc.F * KINT * self.surf_c**0.5 * (self.cmax - self.surf_c)**0.5 

        # algebraic equations. Equation AFTER the colon is relevant, 
        # ( self.XX BEFORE the colon can be ignored as it's just a syntactical requirement )

        model.algebraic.update({
            self.phi: j0 * 2*pybamm.sinh(x) - self.i_int,
            self.i_sei: is_rhs - self.i_sei,
            self.i_int: -self.i_int - self.i_sei + self.j
        })

        model.initial_conditions.update({
            self.c: self.c0,
            self.phi: self.OCP_INIT,
            self.i_sei: 1e-8,
            self.i_int: 1e-2,
            self.sei_L: self.sei0,
        }) 

        # TODO: Sign check on surface boundary condition
        model.boundary_conditions.update({
            self.c: {
                "left":  (0, "Neumann"),
                "right": (-self.i_int / (cc.F * self.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        # model.variables.update{}
        model.variables.update({
            self.c.name: pybamm.PrimaryBroadcast(self.surf_c, self.domain),
        })
        SET_MODEL_VARS(model,
            [
                self.phi, 
                self.i_int, 
                self.i_sei,
                self.sei_L
            ]
        )

    def attach_parameters(self, parameters: dict):
        BIND_VALUES(parameters, {
            self.c0:               "[input]",
            self.L:                p.NEG_ELEC_THICKNESS.sample(),
            self.eps_n:            p.NEG_ELEC_POROSITY.sample(),
            self.cmax:             p.NEG_CSN_MAX.sample(),

            self.ocp:              p.NEG_OCP2,
            self.D:                p.NEG_DIFFUSION.sample(),
            self.R:                p.PARTICLE_RADIUS.sample(),
            self.sei0:             "[input]",
        })

        self.c0.set_value(p.NEG_CSN_INITIAL.sample()) 

if __name__ == '__main__':
    import params as p
    import numpy as np
    import pandas as pd

    HOURS = 10 
    I_INPUT = 2.4
    DISCRETE_PTS = 30
    TIME_PTS = 250

    geo = {}
    model = pybamm.BaseModel()

    iapp = pybamm.Parameter("Input Current") 
    parameters = {}

    a = Anode("Anode", iapp)
    a.process_model(model)
    a.process_geometry(geo)

    BIND_VALUES(parameters, {
        iapp:               "[input]",
        a.c0:               "[input]",
        a.L:                p.NEG_ELEC_THICKNESS.sample(),
        a.eps_n:            p.NEG_ELEC_POROSITY.sample(),
        a.cmax:             p.NEG_CSN_MAX.sample(),

        a.ocp:              p.NEG_OCP,
        a.D:                p.NEG_DIFFUSION.sample(),
        a.R:                p.PARTICLE_RADIUS.sample(),
        a.sei0:             "[input]",
        a.charging:            "[input]",
    })

    model.events += [
        pybamm.Event("Min Concentration", a.surf_c - 500),
        pybamm.Event("Max Concentration", p.NEG_CSN_INITIAL.get_value() + 100 - a.surf_c)
    ]

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    particles = [a]
    mesh = pybamm.Mesh(geo, 
        { d.domain: pybamm.Uniform1DSubMesh for d in particles },
        { d.r: DISCRETE_PTS for d in particles }
    )

    disc = pybamm.Discretisation(mesh, 
        { d.domain: pybamm.FiniteVolume() for d in particles }
    )


    disc.process_model(model)

    cycles = 3
    solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, extra_options_setup={"max_num_steps": 100000})

    time_steps = np.linspace(0, 3600 * HOURS, TIME_PTS)
    total_time_steps = np.linspace(0, 3600 * HOURS * cycles, TIME_PTS * cycles)
    
    sign = -1
    inps = {
        iapp.name: sign * I_INPUT,
        a.c0.name: p.NEG_CSN_INITIAL.get_value(),
        a.sei0.name: 5.e-9,
        a.charging.name: 0 if (sign == -1) else 1
    }

    ### EVERYTHING BELOW THIS IS JUST RUNNING / CAPTURING SIMULATION DATA.
    ### NO PARAMETER-RELEVANT CODE BELOW

    outputs = SET_OUTPUTS([a.c, a.phi, a.sei_L])
    caps = []
    subdfs = []

    solution = None
    prev_time = 0

    for _ in range(cycles):

        solution = solver.solve(model, time_steps, inputs=inps)

        subdf = pd.DataFrame(columns=['Time'] + outputs)
        subdf['Time'] = solution.t + prev_time
        prev_time += solution.t[-1]

        caps.append( I_INPUT * solution.t[-1] / 3600 )

        ## KEYS ARE VARIABLES
        for key in outputs:
            data = solution[key].entries
            if len(data.shape) == 2:
                data = data[-1] # last node (all nodes 'equal' due to broadcasted surface concentration)

            subdf[key] = data

        subdfs.append(subdf)

        sign *= -1
        BIND_VALUES(inps, 
            {
                iapp: sign * I_INPUT,
                a.c0: solution[a.c.name].entries[-1][-1],
                a.sei0: solution[a.sei_L.name].entries[-1],
                a.charging: 0 if (sign == -1) else 1
            }
        )

    df = pd.concat(subdfs, ignore_index=True)
    
    print(df)
    print(caps)

    df.to_csv(f"ANODE_{cycles}.csv", index=False)
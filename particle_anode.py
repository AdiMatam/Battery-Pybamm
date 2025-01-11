import pybamm
import consts as cc
from consts import SET_MODEL_VARS, SET_OUTPUTS, BIND_VALUES
from params import NEG_OCP
from single_particle import SingleParticle
from wrapped_parameter import WrappedParameter
import params as p

#pybamm.set_logging_level("DEBUG")

class Anode(SingleParticle): 
    OCP_INIT = 0.09280796340471076

    def __init__(self, name: str, iapp: pybamm.Variable):

        super().__init__(name, -1, iapp)

        self.i_sei = pybamm.Variable(name + " Side Current")
        self.i_int = pybamm.Variable(name + " Intercalation Current")
        self.sei_L = pybamm.Variable(name + " SEI Length")
        self.sei0 = WrappedParameter(name + " Initial SEI Length")

    def process_model(self, model: pybamm.BaseModel, charging):
    # def process_model(self, model: pybamm.BaseModel):
        flux = self.D * -pybamm.grad(self.c)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

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
        # x = cc.F / (2 * cc.R_GAS * cc.T) * (self.phi - self.ocp)

        # j = 2*j0*sinh(F/2RT * (V - U - PSEI))
        # j/(2*j0) = sinh(F/2RT * (V - U - PSEI))
        # Asinh(X) = F/2RT * (V-U-PSEI)
        # 2RT*Asinh(X)/F + U + PSEI = V


        ## SEE PAPER
        kfs = p.AGING * 1.36e-12 #* 10
        cec_init = 0.05 * 4541
        is_rhs = charging * -cc.F*kfs*cec_init * pybamm.exp( (-0.5*cc.F)/(cc.R_GAS*cc.T) * (self.phi - (self.sei_L/KSEI)*self.j) ) 

        j0 = cc.F * KINT * self.surf_c**0.5 * (self.cmax - self.surf_c)**0.5 

        # algebraic equations. Equation AFTER the colon is relevant, 
        # ( self.XX BEFORE the colon can be ignored as it's just a syntactical requirement )

        model.algebraic.update({
            self.phi: j0 * 2*pybamm.sinh(x) - self.i_int,
            # self.phi: j0 * 2 * pybamm.sinh(x) - self.j,
            self.i_sei: is_rhs - self.i_sei,
            self.i_int: -self.i_int - self.i_sei + self.j
        })

        model.initial_conditions.update({
            self.c: self.c0,
            self.phi: self.OCP_INIT,
            self.i_sei: 0,
            self.i_int: 1e-2,
            self.sei_L: self.sei0,
        }) 

        # TODO: Sign check on surface boundary condition
        model.boundary_conditions.update({
            self.c: {
                "left":  (0, "Neumann"),
                "right": (-self.i_int / (cc.F * self.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
                # "right": (-self.j / (cc.F * self.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        # model.variables.update{}
        model.variables.update({
            self.c.name: self.c # pybamm.PrimaryBroadcast(self.surf_c, self.domain),
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
        self.sei0.set_value(p.SEI_INITIAL.sample()) 


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

    ano = Anode("Anode", iapp)
    ano.process_model(model)
    ano.process_geometry(geo)
    ano.attach_parameters(parameters)

    BIND_VALUES(parameters, {
        iapp: "[input]",
    })

    # model.events += [
    #     pybamm.Event("Min Concentration", ano.surf_c - 500),
    #     pybamm.Event("Max Concentration", ano.cmax.value + 100 - ano.surf_c)
    # ]

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    particles = [ano]
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
        ano.c0.name: p.NEG_CSN_INITIAL.sample(),
    }

    solution = solver.solve(model, time_steps, inputs=inps)
    # solution.plot([ano.c.name])

    subdf = pd.DataFrame(columns=['Time', 'Anode Concentration', 'Anode Potential'])
    subdf['Time'] = solution.t
    subdf['Anode Concentration'] = solution[ano.c.name].entries[-1]
    subdf['Anode Potential'] = solution[ano.phi.name].entries

    t = []
    conc = []
    pot = []

    with open("sundata.txt") as f:
        for line in f:
            data = line.split("|")[:-1]
            t.append(float(data[0].strip()))
            conc.append(float(data[2].strip()))
            pot.append(float(data[4].strip()))

    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(2)
    
    ax[0].plot(solution.t, subdf['Anode Concentration'], label="Pybamm")
    ax[0].plot(t, conc, label="Sundials")
    ax[1].plot(solution.t, subdf['Anode Potential'], label="Pybamm")
    ax[1].plot(t, pot, label="Pybamm")
    ax[0].legend()
    ax[1].legend()
    plt.show()


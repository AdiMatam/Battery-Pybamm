import pybamm
import consts as c

pybamm.set_logging_level("DEBUG")

POS_OCP_INIT = 4.08138601219583

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Variable, r_val: float):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.c0 = pybamm.Parameter(name + " pInitial Concentration")
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.cmax = pybamm.Parameter(name + " pMax Concentration")
        self.D = pybamm.Parameter(name + " pDiffusion Coefficient")

        self.iflag = pybamm.Parameter(name + " pCharge?")

        self.iapp = iapp

        self.i_s = pybamm.Variable(name + " vSide Current")
        self.i_int = pybamm.Variable(name + " vIntercalation Current")
        self.seiL = pybamm.Variable(name + " vSEI Thickness")
        self.sei0 = pybamm.Parameter(name + " pInitial SEI Length")

        self.phi = pybamm.Variable(name + " vPhi")
        self.csn = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.surf_csn = pybamm.surf(self.csn)

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")
        self.r_val = r_val

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.name + " fParticle Voltage", 
            {
                "Stoichiometry": sto
            }
        )

    
    def process_model(self, model: pybamm.BaseModel, clear=False):
        flux = self.D * -pybamm.grad(self.csn)
        # dc/dt = d^2c/dr^2
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.csn: dcdt,
        })

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.a_term = (3 * (1 - self.eps_n)) / self.r_val
        self.j = (self.charge * self.iapp) / (self.L * self.a_term)

        self.ocp = self.u_func(self.surf_csn / self.cmax)

        cathode_kint = 1.04e-11
        j0 = c.F * cathode_kint * self.surf_csn**0.5 * (self.cmax - self.surf_csn)**0.5 
        x = c.F / (2 * c.R_GAS * c.T) * (self.phi - self.ocp)

        model.algebraic.update({
            self.phi: j0 * 2 * pybamm.sinh(x) - self.j,
        })

        model.initial_conditions.update({
            self.csn: pybamm.x_average(self.c0),
            self.phi: POS_OCP_INIT
        }) 

        model.boundary_conditions.update({
            self.csn: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * self.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        model.variables.update({
            self.csn.name: self.csn,
            self.phi.name: self.phi,
            "JTerm": self.j

        })

    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

        ## TODO: need another 'domain' for SEI
        geo.update({
            self.domain: {self.r: {"min": 0, "max": self.r_val}}
        })
    
    def process_parameters(self, all_params: dict, particle_params: dict, clear=False):
        if clear:
            all_params.clear()

        all_params.update(
            {key.name : value for key, value in particle_params.items()}
        )









if __name__ == '__main__':
    import params as p
    import numpy as np


    HOURS = 11 # 20 / 2
    I_TOTAL = +1.2027659291666666 * 2
    # HOURS = 20
    # I_TOTAL = +1.2027659291666666
    DISCRETE_PTS = 30
    TIME_PTS = 250

    pos_cap = p.POS_CSN_MAX.get_value() - p.POS_CSN_INITIAL.get_value()
    pos_cap *= p.POS_ELEC_THICKNESS.get_value() * (1-p.POS_ELEC_POROSITY.get_value()) * (c.F / 3600)

    neg_cap = p.NEG_CSN_INITIAL.get_value()
    neg_cap *= p.NEG_ELEC_THICKNESS.get_value() * (1-p.NEG_ELEC_POROSITY.get_value()) * (c.F / 3600)
    
    print(pos_cap, neg_cap)
    print(neg_cap / 20)

    geo = {}
    model = pybamm.BaseModel()

    i_total = pybamm.Parameter("Input Current") 
    parameters = {}

    a = SingleParticle("Cathode", +1, i_total, p.PARTICLE_RADIUS.get_value())
    a.process_model(model)
    a.process_geometry(geo)

    a.process_parameters(parameters, {
        i_total:            "[input]",
        a.c0:               "[input]",
        a.L:                p.POS_ELEC_THICKNESS.get_value(),
        a.eps_n:            p.POS_ELEC_POROSITY.get_value(),
        a.cmax:             p.POS_CSN_MAX.get_value(),

        a.ocp:              p.pos_ocp,
        a.D:                p.POS_DIFFUSION.get_value(),
    })

    particles = [a]
    mesh = pybamm.Mesh(geo, 
        { d.domain: pybamm.Uniform1DSubMesh for d in particles },
        { d.r: DISCRETE_PTS for d in particles }
    )

    disc = pybamm.Discretisation(mesh, 
        { d.domain: pybamm.FiniteVolume() for d in particles }
    )

    param_ob = pybamm.ParameterValues(parameters)
    param_ob.process_model(model)
    param_ob.process_geometry(geo)

    disc.process_model(model)

    cycles = 2
    solver = pybamm.CasadiSolver(mode='safe', atol=1e-6, rtol=1e-5, extra_options_setup={"max_num_steps": 100000})

    time_steps = np.linspace(0, 3600 * HOURS, TIME_PTS)
    total_time_steps = np.linspace(0, 3600 * HOURS * cycles, TIME_PTS * cycles)
    
    sign = -1
    last_conc = p.POS_CSN_INITIAL.get_value()

    conc_array = np.empty((TIME_PTS*cycles,))
    volt_array = np.empty((TIME_PTS*cycles,))
    J_array = np.empty((TIME_PTS*cycles,))

    for i in range(cycles):
        inps = {
            i_total.name: sign * I_TOTAL,
            a.c0.name: last_conc,
        }
        solution = solver.solve(model, time_steps, inputs=inps)

        arr = solution[a.csn.name].entries[-1]
        conc_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = arr
        volt_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution[a.phi.name].entries
        J_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = solution['JTerm'].entries

        last_conc = arr[-1]
        sign *= -1
        print(f"Passed {i}")

    print(len(conc_array))
    print(len(total_time_steps))

    import pandas as pd
    df = pd.DataFrame({
        'Time': total_time_steps,
        'Concentration': conc_array,
        "Voltage": volt_array,
        "JTerm": J_array,
    })

    df.to_csv(f"CATHODE_SEI_{cycles}.csv", index=False)


        ## Butler volmer calculation
        # self.ocp = self.u_func(self.surf_csn / self.cmax)
        # self.overp = (2*c.RTF*pybamm.arcsinh(self.j / (2 * self.j0)))
        # self.phi = self.ocp + self.overp
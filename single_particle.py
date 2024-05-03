import pybamm
import consts as c

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Variable, r_val: float):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.c_0_name = name + " pInitial Concentration"
        self.c_0 = pybamm.Parameter(self.c_0_name)
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.c_max = pybamm.Parameter(name + " pMax Concentration")

        self.iapp = iapp

        self.csn_name = name + " vConcentration"
        self.csn = pybamm.Variable(self.csn_name, domain=self.domain)
        self.surf_csn_name = name + " vSurface Concentration"
        self.surf_csn = pybamm.surf(self.csn)

        self.j0_name = name + " fExchange Current Density"

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")
        self.r_val = r_val

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.name + " fParticle Voltage", 
            {
                "Stoichiometry": sto
            }
        )

    def j0_func(self, c_e, c_s_surf, c_s_max):
        return pybamm.FunctionParameter(
            self.j0_name,
            {
                "Electrolyte Concentration": c_e,
                "Electrode Surface Concentration": c_s_surf,
                "Electrode Max Concentration": c_s_max
            }
        )

    def process_model(self, model: pybamm.BaseModel, electrolyte_conc: float, clear=False):
        if clear:
            model = pybamm.BaseModel()

        flux = c.D * -pybamm.grad(self.csn)
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.csn: dcdt,
        })

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.a_term = (3 * (1 - self.eps_n)) / self.r_val
        self.j = (self.charge * self.iapp) / (self.L * self.a_term)
        self.j0 = self.j0_func(pybamm.Scalar(electrolyte_conc), self.surf_csn, self.c_max)

        self.ocp = self.u_func(self.surf_csn / self.c_max)
        self.bv = self.ocp + (2*c.RTF*pybamm.arcsinh(self.j / (2 * self.j0)))

        model.initial_conditions.update({
            self.csn: pybamm.x_average(self.c_0),
        }) 

        model.boundary_conditions.update({
            self.csn: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * c.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })

        model.variables.update({
            self.surf_csn_name: pybamm.PrimaryBroadcast(
                self.surf_csn, self.domain
            ),
        })
    
    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

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

    I_TOTAL = +1.2027659291666666
    DISCRETE_PTS = 30
    TIME_PTS = 250

    # iapp = pybamm.Variable("Iapp")
    geo = {}
    model = pybamm.BaseModel()

    i_total = pybamm.Parameter("Input Current / Area") 
    parameters = {}

    a = SingleParticle("Positive", +1, i_total, p.PARTICLE_RADIUS.get_value())
    a.process_model(model, p.ELECTROLYTE_CONC.get_value())
    a.process_geometry(geo)

    a.process_parameters(parameters, {
        i_total:            "[input]",
        a.c_0:               "[input]",
        a.L:                p.POS_ELEC_THICKNESS.get_value(),
        a.eps_n:            p.POS_ELEC_POROSITY.get_value(),
        a.c_max:          p.POS_CSN_MAX.get_value(),

        a.j0:               p.POS_J0,
        a.ocp:              p.POS_OCP
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

    cycles = 1
    solver = pybamm.CasadiSolver()
    time_steps = np.linspace(0, 3600 * 20, TIME_PTS)
    total_time_steps = np.linspace(0, 3600 * 20 * cycles, TIME_PTS * cycles)
    
    sign = -1
    last_conc = p.POS_CSN_INITIAL.get_value()

    conc_array = np.empty((250*cycles,))
    volt_array = np.empty((250*cycles,))

    for i in range(cycles):
        inps = {
            i_total.name: sign * I_TOTAL,
            a.c_0.name: last_conc
        }
        solution = solver.solve(model, time_steps, inputs=inps)

        arr = solution[a.surf_csn_name].entries[-1]
        conc_array[i*TIME_PTS: (i*TIME_PTS)+TIME_PTS] = arr

        last_conc = arr[-1]
        sign *= -1

    print(len(conc_array))
    print(len(total_time_steps))

    from  matplotlib import pyplot as plt

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (s)')

    color = 'tab:red'
    ax1.set_ylabel('Concentration (mol / m2)', color=color)
    ax1.plot(total_time_steps, conc_array, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    plt.tight_layout()
    plt.show()


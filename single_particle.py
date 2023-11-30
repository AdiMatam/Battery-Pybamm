import pybamm
import consts as c

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Variable):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.conc_0 = pybamm.Parameter(name + " pInitial Concentration")
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.conc_max = pybamm.Parameter(name + " pMax Concentration")

        self.iapp = iapp
        self.conc = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.surf_conc = pybamm.surf(self.conc)
        self.surf_conc_name = name + " vSurface Concentration"

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")
        self.particle_voltage_name = name + " fParticle Voltage"
        self.j0_name = name + " fExchange Current Density"

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.particle_voltage_name,
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

    def process_model(self, model: pybamm.BaseModel, clear=False):
        if clear:
            model = pybamm.BaseModel()

        flux = c.D * -pybamm.grad(self.conc)
        dcdt = -pybamm.div(flux)

        # solve the diffusion equation (del * del(c))
        model.rhs.update({
            self.conc: dcdt,
        })

        model.initial_conditions.update({
            self.conc: pybamm.x_average(self.conc_0),
        }) 
        
        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        a_term = (3 * (1 - self.eps_n)) / c.R
        self.j = (self.charge * self.iapp) / (self.L * a_term)
        self.j0 = self.j0_func(pybamm.Scalar(c.ELECTROLYTE_CONC), self.surf_conc, self.conc_max)
        
        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * c.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })
        model.variables.update({
            self.surf_conc_name: pybamm.PrimaryBroadcast(
                self.surf_conc, self.domain
            )
        })
    
    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

        geo.update({
            self.domain: {self.r: {"min": 0, "max": c.R}}
        })
    
    def process_parameters(self, all_params: dict, particle_params: dict, clear=False):
        if clear:
            all_params.clear()

        all_params.update(
            {key.name : value for key, value in particle_params.items()}
        )
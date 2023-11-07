import pybamm
import consts as c

from marquis import lico2_ocp_Dualfoil1998 as Up
from marquis import graphite_mcmb2528_ocp_Dualfoil1998 as Un

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Parameter):
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

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")
        self.voltage_name = self.name + " fVoltage"

    def j0_func(self, c_e, c_s_surf, c_s_max):
        return pybamm.FunctionParameter(
            self.name + " fExchange Current Density",
            {
                "Electrolyte Concentration": c_e,
                "Electrode Surface Concentration": c_s_surf,
                "Electrode Max Concentration": c_s_max
            }
        )

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.voltage_name,
            {
                "Stoichiometry": sto
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
        
        a_term = (3 * (1 - self.eps_n)) / c.R

        # self.charge = +1 for positive electrode
        #               -1 for negative electrode
        self.j = (self.charge * self.iapp) / (self.L * a_term)

        RTF = c.R_GAS * c.T / c.F
        
        self.j0 = self.j0_func(pybamm.Scalar(c.ELECTROLYTE_CONC), self.surf_conc, self.conc_max)
        overpotential = 2 * RTF * pybamm.arcsinh(self.j / (2 * self.j0))

        self.ocp = self.u_func(self.surf_conc / self.conc_max) 
        self.voltage = overpotential + self.ocp
        
        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * c.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })
        model.variables.update({
            self.conc.name: pybamm.PrimaryBroadcast(
                self.surf_conc, self.domain #self.conc, 
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
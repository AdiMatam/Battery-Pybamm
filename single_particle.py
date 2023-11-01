import pybamm
import consts as c
import numpy as np

class SingleParticle:
    def __init__(self, name: str, charge: int, current: pybamm.Parameter):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.conc_0 = pybamm.Parameter(name + " pInitial Concentration")
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.conc_max = pybamm.Parameter(name + " pMax Concentration")
        self.current = current

        self.conc_name = name + " vConcentration"
        self.surf_conc_name = name + " vSurface Concentration"

        self.conc = pybamm.Variable(self.conc_name, domain=self.domain)
        self.surf_conc = pybamm.surf(self.conc)
        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")

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
        self.j = (self.charge * self.current) / (self.L * a_term)
        
        ## TEMPORARY, TAKEN FROM PYBAMM CODE
        if self.charge == 1:
            m_ref = 6 * 10 ** (-7)  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
            self.j0 = m_ref * 1000**0.5 * self.surf_conc**0.5 * (self.conc_max - self.surf_conc) ** 0.5
        else:
            m_ref = 2 * 10 ** (-5)  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
            self.j0 = m_ref * 1000**0.5 * self.surf_conc**0.5 * (self.conc_max - self.surf_conc) ** 0.5

        self.j_name = self.name + " vCurrent Density"
        self.j0_name = self.name + " vExchange Current Density"

        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (-self.j / (c.F * c.D), "Neumann") # outer boundary condition (dc/dr behavior @r=R)
            },
        })
        model.variables.update({
            self.conc_name: self.conc, 
            self.surf_conc_name: self.surf_conc,
            self.j_name: self.j,
            self.j0_name: self.j0,
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
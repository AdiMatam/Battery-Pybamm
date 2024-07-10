import pybamm
import consts as cc
from abc import ABC, abstractmethod

class SingleParticle:
    def __init__(self, name: str, charge: int, iapp: pybamm.Variable):
        self.name = name
        self.domain = name + " dDomain"
        self.charge = charge

        self.c0 = pybamm.Parameter(name + " pInitial Concentration")
        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")
        self.cmax = pybamm.Parameter(name + " pMax Concentration")
        self.D = pybamm.Parameter(name + " pDiffusion Coefficient")
        self.R = pybamm.Parameter(name + " pParticle Radius")


        self.iapp = iapp

        self.phi = pybamm.Variable(name + " vPhi")
        self.c = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.surf_c = pybamm.surf(self.c)

        self.ocp = self.u_func(self.surf_c / self.cmax)

        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")


        a_term = (3 * (1 - self.eps_n)) / self.R
        self.j = (self.charge * self.iapp) / (self.L * a_term)

    def u_func(self, sto):
        return pybamm.FunctionParameter(
            self.name + " fParticle Voltage", 
            {
                "Stoichiometry": sto
            }
        )

    @abstractmethod    
    def process_model(self, model: pybamm.BaseModel):
        pass

    def process_geometry(self, geo: dict):
        geo.update({
            self.domain: {self.r: {"min": 0, "max": self.R}}
        })
    

if __name__ == '__main__':
    pass
import pybamm
from single_particle import SingleParticle
import consts as c

class Cell:
    CELLS = list()
    def __init__(self, name: str, model: pybamm.BaseModel, geo:dict):
        # if self in self.CELL_NAMES:
            # raise ValueError("Must have unique cell names/IDs")

        self.CELLS.append(self)
        self.name = name
        self.model = model

        self.iapp = pybamm.Variable(name + " Iapp")

        self.pos = SingleParticle(name + " Pos Particle", +1, self.iapp)
        self.neg = SingleParticle(name + " Neg Particle", -1, self.iapp)

        self.pos.process_model(model)
        self.neg.process_model(model)

        ### FIGURE THIS OUT
        self.voltage = pybamm.Variable(name + " Voltage")

        model.variables.update({
            self.pos.phi.name: self.pos.phi,
            self.neg.phi.name: self.neg.phi,
            self.iapp.name: self.iapp
        })

        # ## PARAMETERIZE
        # model.events += [
            # pybamm.Event("Voltage Min Cutoff", self.voltage - 3.0)
        # ]

        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)


    def get_meshed_objects(self):
        return { self.pos.domain: pybamm.Uniform1DSubMesh, self.neg.domain: pybamm.Uniform1DSubMesh } 

    def get_radial_mesh(self, pts: int):
        return { self.pos.r: pts, self.neg.r: pts }

    def get_discretized(self):
        return { self.pos.domain: pybamm.FiniteVolume(), self.neg.domain: pybamm.FiniteVolume() }

    def set_parameters(self, param_dict: dict):
        self.pos.process_parameters(param_dict, {
            self.pos.conc_0:    c.POS_CSN_INITIAL,
            self.pos.L:         c.POS_ELEC_THICKNESS,
            self.pos.eps_n:     c.POS_ELEC_POROSITY,
            self.pos.conc_max:  c.POS_CSN_MAX,
            self.pos.j0:        c.POS_EXCHANGE_CURRENT_DENSITY,
            self.pos.ocp:       c.POS_OPEN_CIRCUIT_POTENTIAL
        })

        self.neg.process_parameters(param_dict, {
            self.neg.conc_0:    c.NEG_CSN_INITIAL,
            self.neg.L:         c.NEG_ELEC_THICKNESS,
            self.neg.eps_n:     c.NEG_ELEC_POROSITY,
            self.neg.conc_max:  c.NEG_CSN_MAX,
            self.neg.j0:        c.NEG_EXCHANGE_CURRENT_DENSITY,
            self.neg.ocp:       c.NEG_OPEN_CIRCUIT_POTENTIAL
            
        })

    
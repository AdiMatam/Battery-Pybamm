import pybamm
from single_particle import SingleParticle
import consts as c
import params as p

class Cell:
    CELLS = list()
    def __init__(self, name: str, model: pybamm.BaseModel, geo:dict, parameters:dict):
        # if self in self.CELL_NAMES:
            # raise ValueError("Must have unique cell names/IDs")

        self.CELLS.append(self)
        self.name = name
        self.model = model

        self.iapp = pybamm.Variable(name + " Iapp")

        self.pos = SingleParticle(name + " Pos Particle", +1, self.iapp)
        self.neg = SingleParticle(name + " Neg Particle", -1, self.iapp)

        self.electrolyte_conc    = p.ELECTROLYTE_CONC.rand_sample()
        self.__create_parameter_samples()

        self.pos.process_model(model, self.electrolyte_conc)
        self.neg.process_model(model, self.electrolyte_conc)

        self.__attach_parameters(parameters)

        model.variables.update({
            self.pos.phi.name: self.pos.phi,
            self.neg.phi.name: self.neg.phi,
            self.iapp.name: self.iapp
        })

        ## PARAMETERIZE
        model.events += [
            pybamm.Event("Min Concentration Cutoff", self.neg.surf_conc - self.neg_csn_min)
        ]

        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)

        ## pos_phi = p_OCP() @ t=0
        ## neg_phi = n_OCP() @ t=0
        self.pos_phi_init = p.POS_OPEN_CIRCUIT_POTENTIAL(self.pos_csn_initial / self.pos_csn_max)
        self.neg_phi_init = p.NEG_OPEN_CIRCUIT_POTENTIAL(self.neg_csn_initial / self.neg_csn_max)

        self.capacity = self.__compute_capacity()

    def __compute_capacity(self):
        # (csn_max - csn_min) * L * (1-eps_n) * (mols->Ah)
        pos_cap = (self.pos_csn_max - self.pos_csn_initial) * self.pos_elec_thickness * (1-self.pos_elec_porosity) * (c.F / 3600)
        neg_cap = (self.neg_csn_initial - self.neg_csn_min) * self.neg_elec_thickness * (1-self.neg_elec_porosity) * (c.F / 3600)

        return min(pos_cap, neg_cap)

    def __create_parameter_samples(self):    

        self.pos_csn_max         = p.POS_CSN_MAX.rand_sample()
        self.pos_csn_initial     = p.POS_CSN_INITIAL.rand_sample()
        self.pos_elec_thickness  = p.POS_ELEC_THICKNESS.rand_sample()
        self.pos_elec_porosity   = p.POS_ELEC_POROSITY.rand_sample()

        self.neg_csn_min         = p.NEG_CSN_MIN.rand_sample()
        self.neg_csn_initial     = p.NEG_CSN_INITIAL.rand_sample()
        self.neg_csn_max         = p.NEG_CSN_MAX.rand_sample()
        self.neg_elec_thickness  = p.NEG_ELEC_THICKNESS.rand_sample()
        self.neg_elec_porosity   = p.NEG_ELEC_POROSITY.rand_sample()


    def __attach_parameters(self, param_dict: dict):
        self.pos.process_parameters(param_dict, {
            self.pos.conc_0:    self.pos_csn_initial,
            self.pos.L:         self.pos_elec_thickness,
            self.pos.eps_n:     self.pos_elec_porosity,
            self.pos.conc_max:  self.pos_csn_max,

            self.pos.j0:        p.POS_EXCHANGE_CURRENT_DENSITY,
            self.pos.ocp:       p.POS_OPEN_CIRCUIT_POTENTIAL
        })

        self.neg.process_parameters(param_dict, {
            self.neg.conc_0:    self.neg_csn_initial,
            self.neg.L:         self.neg_elec_thickness,
            self.neg.eps_n:     self.neg_elec_porosity,
            self.neg.conc_max:  self.neg_csn_max,

            self.neg.j0:        p.NEG_EXCHANGE_CURRENT_DENSITY,
            self.neg.ocp:       p.NEG_OPEN_CIRCUIT_POTENTIAL
            
        })

    

    # def get_meshed_objects(self):
        # return { self.pos.domain: pybamm.Uniform1DSubMesh, self.neg.domain: pybamm.Uniform1DSubMesh } 

    # def get_radial_mesh(self, pts: int):
        # return { self.pos.r: pts, self.neg.r: pts }

    # def get_discretized(self):
        # return { self.pos.domain: pybamm.FiniteVolume(), self.neg.domain: pybamm.FiniteVolume() }
import pybamm
from single_particle import SingleParticle
import consts as c
import params as p

class Cell:
    CELLS = list()
    def __init__(self, name: str, model: pybamm.BaseModel, geo:dict, parameters:dict, iapp: pybamm.Variable):
        if name in self.CELLS:
            raise ValueError("Must have unique cell names/IDs")

        self.CELLS.append(name)
        self.name = name
        self.model = model

        self.iapp = iapp
        self.volt = pybamm.Variable(name + " Voltage")

        self.particle_radius = p.PARTICLE_RADIUS.rand_sample()
        self.pos = SingleParticle(name + " Pos Particle", +1, iapp, self.particle_radius)
        self.neg = SingleParticle(name + " Neg Particle", -1, iapp, self.particle_radius)

        ## TEMPORARILY ELECTROLYLTE HANDLED DONE THIS WAY
        self.electrolyte_conc    = p.ELECTROLYTE_CONC.rand_sample()

        self.pos.process_model(model, self.electrolyte_conc)
        self.neg.process_model(model, self.electrolyte_conc)

        self.__create_parameter_samples()
        self.__attach_parameters(parameters)
        
        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)

        self.capacity = self.__compute_capacity()

    def __compute_capacity(self):
        # (csn_max - csn_min) * L * (1-eps_n) * (mols->Ah)
        pos_cap = (self.pos_csn_maxval - self.pos_csn_ival) * self.pos_L * (1-self.pos_eps_n) * (c.F / 3600)
        neg_cap = (self.neg_csn_ival - self.neg_csn_min) * self.neg_L * (1-self.neg_eps_n) * (c.F / 3600)

        return min(pos_cap, neg_cap)

    
    ### RETHINK PARAMETER GENERATION AND ATTACHMENT
    ### THIS KIND OF SUCKS (too many names, identifiers...)
    def __create_parameter_samples(self):    

        self.pos_csn_maxval         = p.POS_CSN_MAX.rand_sample()
        self.pos_csn_ival     = p.POS_CSN_INITIAL.rand_sample()
        self.pos_L  = p.POS_ELEC_THICKNESS.rand_sample()
        self.pos_eps_n   = p.POS_ELEC_POROSITY.rand_sample()

        self.neg_csn_min         = p.NEG_CSN_MIN.rand_sample()
        self.neg_csn_ival     = p.NEG_CSN_INITIAL.rand_sample()
        self.neg_csn_maxval         = p.NEG_CSN_MAX.rand_sample()
        self.neg_L  = p.NEG_ELEC_THICKNESS.rand_sample()
        self.neg_eps_n   = p.NEG_ELEC_POROSITY.rand_sample()


    def __attach_parameters(self, param_dict: dict):
        self.pos.process_parameters(param_dict, {
            self.pos.c_0:       "[input]",
            self.pos.L:         self.pos_L,
            self.pos.eps_n:     self.pos_eps_n,
            self.pos.c_max:  self.pos_csn_maxval,

            self.pos.j0:        p.POS_J0,
            self.pos.ocp:       p.POS_OCP
        })

        self.neg.process_parameters(param_dict, {
            self.neg.c_0:       "[input]",
            self.neg.L:         self.neg_L,
            self.neg.eps_n:     self.neg_eps_n,
            self.neg.c_max:  self.neg_csn_maxval,

            self.neg.j0:        p.NEG_J0,
            self.neg.ocp:       p.NEG_OCP
            
        })

        




        # self.voltage_name = f"{self.name} Voltage"
        # self.voltage = pybamm.Variable(self.voltage_name)   # self.pos.phi - self.neg.phi
        # model.variables.update({
            # # self.voltage_name: self.voltage,
            # self.iapp.name: self.iapp
        # })

        ## PARAMETERIZE
        # model.events += [
            # pybamm.Event("Min Concentration Cutoff", self.neg.surf_csn - self.neg_csn_min),
            # pybamm.Event("Max Concentration Cutoff", self.pos_csn_maxval - self.pos.surf_csn),
            # pybamm.Event("Min Voltage Cutoff", (self.voltage) - voltage_cutoff[0]),
            # pybamm.Event("Max Voltage Cutoff", voltage_cutoff[1] - (self.voltage)),
        # ]


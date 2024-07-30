import pybamm
from ocp import NEG_OCP, POS_OCP
from particle_anode import Anode
from particle_cathode import Cathode
from consts import BIND_VALUES, SET_MODEL_VARS, SET_OUTPUTS
import params as p

class Cell:
    CELLS = list()
    def __init__(self, name: str,iapp: pybamm.Variable, charging: pybamm.Parameter,
            model: pybamm.BaseModel, geo:dict, parameters:dict
    ):

        if name in self.CELLS:
            raise ValueError("Must have unique cell names/IDs")

        self.CELLS.append(name)
        self.name = name
        self.model = model

        self.iapp = iapp
        
        self.voltage = pybamm.Variable(name + " Voltage")

        self.pos = Cathode(name + " Cathode", iapp)
        self.neg = Anode(name + " Anode", iapp)

        self.sei = self.neg.sei_L

        self.vvolt = self.pos.phi - self.neg.phi

        self.capacity = pybamm.Variable(name + " Capacity by Area")

        discharging = pybamm.Negate(pybamm.Subtraction(charging, 1))
        model.rhs.update({
            self.capacity: discharging * pybamm.AbsoluteValue(iapp / 3600)
        })

        model.initial_conditions.update({
            self.capacity: 0
        })

        self.pos.process_model(model)
        self.neg.process_model(model, charging)

        self.__attach_parameters(parameters)
        
        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)

        model.initial_conditions.update({
            self.voltage: POS_OCP(self.pos.c0.value / self.pos.cmax.value) 
                - NEG_OCP(self.neg.c0.value / self.neg.cmax.value)
        })

        model.variables.update({
            self.voltage.name: self.vvolt,
            self.capacity.name: self.capacity
        })

    def __attach_parameters(self, param_dict: dict):

        rad = p.PARTICLE_RADIUS.rand_sample()
        BIND_VALUES(param_dict, {
            self.pos.c0:               "[input]",
            self.pos.L:                p.POS_ELEC_THICKNESS.rand_sample(),
            self.pos.eps_n:            p.POS_ELEC_POROSITY.rand_sample(),
            self.pos.cmax:             p.POS_CSN_MAX.rand_sample(),

            self.pos.ocp:              p.POS_OCP,
            self.pos.D:                p.POS_DIFFUSION.rand_sample(),
            self.pos.R:                rad,
        })

        BIND_VALUES(param_dict, {
            self.neg.c0:               "[input]",
            self.neg.L:                p.NEG_ELEC_THICKNESS.rand_sample(),
            self.neg.eps_n:            p.NEG_ELEC_POROSITY.rand_sample(),
            self.neg.cmax:             p.NEG_CSN_MAX.rand_sample(),

            self.neg.ocp:              p.NEG_OCP,
            self.neg.D:                p.NEG_DIFFUSION.rand_sample(),
            self.neg.R:                rad,
            self.neg.sei0:             "[input]",
        })

        self.pos.c0.set_value(p.POS_CSN_INITIAL.rand_sample()) 
        self.neg.c0.set_value(p.NEG_CSN_INITIAL.rand_sample()) 

        # self.GET = param_dict.copy()

        # BIND_VALUES(self.GET, {
            # self.pos.c0: p.POS_CSN_INITIAL.rand_sample(),
            # self.neg.c0: p.NEG_CSN_INITIAL.rand_sample(),
        # })

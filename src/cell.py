import pybamm
from src.particle_anode import Anode
from src.particle_cathode import Cathode
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
        
        self.pos = Cathode(name + " Cathode", iapp)
        self.neg = Anode(name + " Anode", iapp)

        ## cell-level 'reference' to sei length 
        self.sei = self.neg.sei_L

        self.voltage = pybamm.Variable(name + " Voltage")
        self.vvolt = self.pos.phi - self.neg.phi

        self.capacity = pybamm.Variable(name + " Capacity by Area")
        discharging = pybamm.Negate(pybamm.Subtraction(charging, 1))
        model.rhs.update({
            self.capacity: discharging * pybamm.AbsoluteValue(iapp / 3600)
        })

        model.initial_conditions.update({
            self.capacity: 0
        })

        self.capacity_value = 0

        self.pos.process_model(model)
        self.neg.process_model(model, charging)

        self.pos.attach_parameters(parameters)
        self.neg.attach_parameters(parameters)
        # self.__attach_parameters(parameters)
        
        self.pos.process_geometry(geo)
        self.neg.process_geometry(geo)

        model.variables.update({
            self.voltage.name: self.vvolt,
            self.capacity.name: self.capacity
        })

        self.__stopC()

    def __stopC(self):
        self.model.events.extend([
            pybamm.Event(f"{self.name} Min Anode Concentration Cutoff", self.neg.surf_c - 10),
            pybamm.Event(f"{self.name} Max Cathode Concentration Cutoff", self.pos.cmax - self.pos.surf_c),

            pybamm.Event(f"{self.name} Max Anode Concentration Cutoff", self.neg.cmax - self.neg.surf_c),
        ])

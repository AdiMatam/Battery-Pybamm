import pybamm 

class WrappedParameter(pybamm.Parameter):
    def __init__(self, name: str):
        super().__init__(name)
        self.value = None
    
    def set_value(self, value):
        self.value = value
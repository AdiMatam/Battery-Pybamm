
import random

class VariatedParameter:
    OVERRIDE_VARIATON = False

    def __init__(self, value: float, low_high: tuple):
        self.value = value
        self.low_high = low_high

    def __float__(self):
        return float(self.get_value())

    @classmethod
    def from_percent(cls, value: float, percent: float=1):
        percent *= int(not cls.OVERRIDE_VARIATON)
        offset = value * (percent / 100)
        return cls(value, (value - offset, value + offset))

    @classmethod
    def from_gaussian(cls, value: float, stddev: float):
        raise NotImplementedError("TBD")

    def rand_sample(self):
        return random.uniform(self.low_high[0], self.low_high[1]) 

    def get_value(self):
        return self.value


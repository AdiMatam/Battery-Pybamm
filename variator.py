import random

class Variator:
    def __init__(self, mean_value: float, func):
        self.mean_value = mean_value
        self.func = func

    @classmethod
    def from_percent(cls, mean: float, percent: float):
        offset = mean * (percent / 100)
        func = lambda: random.uniform(mean - offset, mean + offset)
        return cls(mean, func)

    @classmethod
    def from_gaussian_percent(cls, mean: float, percent: float):
        stddev = mean * (percent / 100)
        func = lambda: random.gauss(mean, stddev)
        return cls(mean, func)

    @classmethod
    def from_gaussian_stddev(cls, mean: float, stddev: float):
        func = lambda: random.gauss(mean, stddev)
        return cls(mean, func)

    def sample(self):
        return self.func()

    def get_mean_value(self):
        return self.mean_value


import random

class Variator:
    ALL = []    


    def __init__(self, name: str, mean_value: float, func, string: str):
        self.name = name
        self.mean_value = mean_value
        self.func = func
        self.string = string
        Variator.ALL.append(self)

    @classmethod
    def from_percent(cls, name: str, mean: float, percent: float):
        offset = mean * (percent / 100)
        func = lambda: random.uniform(mean - offset, mean + offset)
        return cls(name, mean, func, f"Uniform: {percent:.3f}%")

    @classmethod
    def from_gaussian_percent(cls, name: str, mean: float, percent: float):
        stddev = mean * (percent / 100)
        func = lambda: random.gauss(mean, stddev)
        return cls(name, mean, func, f"Gaussian: {percent:.3f}% stddev")

    @classmethod
    def from_gaussian_stddev(cls, name: str, mean: float, stddev: float):
        func = lambda: random.gauss(mean, stddev)
        return cls(name, mean, func, f"Gaussian: {stddev:.3f} stddev")

    def __str__(self):
        return self.string

    def sample(self):
        return self.func()

    def get_mean_value(self):
        return self.mean_value

    @classmethod
    def JSON(cls):
        master = {}
        for variator in cls.ALL:
            key = f"{variator.name}"
            val = f"({variator.mean_value}, {variator.string})"
            master[key] = val

        d = {}
        d['Parameter Variations'] = master
        return d;

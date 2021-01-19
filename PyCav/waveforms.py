import numpy as np

class waveforms:

    def __init__(self,config={}):

        if "form" in config:
            self.form = config["form"]
        else:
            self.form = "constant"

        if "amplitude" in config:
            self.amplitude = config["amplitude"]
        else:
            raise Exception("No amplitude")

        if "period" in config:
            self.period = config["period"]
        else:
            if self.form == "sine" or self.form == "square":
                raise Exception("Need period")
            else:
                self.period = 1.

        if "cycles" in config:
            self.cycles = config["cycles"]
        else:
            self.cycles = 1.

        if "ambient" in config:
            self.ambient = config["ambient"]
        else:
            self.ambient = 1.

        if self.form == "constant":
            self.p = self.p_constant
        elif self.form == "sine":
            self.p = self.p_sine
        elif self.form == "square":
            self.p = self.p_square
        else:
            NotImplementedError

    def p_constant(self,t):
        return self.amplitude

    def p_sine(self,t):
        if t <= self.period * self.cycles:
            return self.amplitude * np.sin(2.*np.pi*t/self.period)
        else:
            return self.ambient

    def p_square(self,t):
        if t <= self.period:
            return self.amplitude
        else:
            return self.ambient

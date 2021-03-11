import numpy as np
from . import Label

class Binary(Label):
    """
    Class that handles the manipulation of binary molecular labels.
    """
    def __init__(self,name,
                      values,
                      context='whole',
                      sensitivity=0.1):
        self.name        = name
        self.values      = values
        self.context     = context
        self.sensitivity = sensitivity
        self.av          = np.mean(values)
        self.entropy     = self.entropy()

    def summary(self,verbose=False):
        if self.entropy < self.sensitivity:
            if self.context == "whole":
                return f"Molecules are {self.name}" if self.av > 0.5 else f"Molecules are not {self.name}"
            elif self.context == "part":
                return f"Molecules contain {self.name}" if self.av > 0.5 else f"Molecules do not contain {self.name}"
            else:
                return f"{self.context} {self.name}" if self.av > 0.5 else f"Not {self.context} {self.name}"

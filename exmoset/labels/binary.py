import numpy as np
from . import Label

class Binary(Label):
    """
    Class that handles the manipulation of binary molecular labels.
    """
    def __init__(self,property,
                      values,
                      verb="are",
                      noun='Molecules',
                      sensitivity=0.1):
        self.property    = property
        self.values      = values
        self.noun        = noun
        self.verb        = verb
        self.sensitivity = sensitivity
        self.av          = np.round(np.mean(values))
        self.entropy     = self.entropy()

    def summary(self,verbose=False):
        if self.entropy < self.sensitivity:
            description = f"{self.noun} {self.verb} {self.property}"
            return description if self.av > 0.5 else self.neg(description)

import numpy as np
from exmoset.abstract import Label

class Continuous(Label):
    """
    Class that handles the manipulation of continuous molecular labels.
    """
    def __init__(self,name,
                      values,
                      context='atom',
                      sensitivity=0.1):
        self.name    = name
        self.values  = values
        self.context = context
        self.av      = np.mean(values)
        self.sensitivity = sensitivity
        self.entropy     = self.entropy()

    def summary(self,verbose=False):
        if self.entropy < self.sensitivity:
            if self.context == "atom":
                return f"Atoms share a similar {self.name} centred around {self.av}"
            elif self.context == "bonds":
                return f"Bonds share a similar {self.name} centred around {self.av}"
            elif self.context == "molecule":
                return f"Molecules share a similar {self.name} centred around {self.av}"
            else:
                return f"{self.context} {self.name} at {self.av}"

from . import Label
import numpy as np

class Multiclass(Label):
    """
    Class that handles the manipulation of multiclass molecular labels.
    """
    def __init__(self,name,
                      values,
                      context='atom',
                      sensitivity=0.1):
        self.name    = name
        self.values  = values
        self.context = context
        self.sensitivity = sensitivity
        self.av      = np.int(np.mean(values))
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

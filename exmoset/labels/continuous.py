import numpy as np
from . import Label

class Continuous(Label):
    """
    Class that handles the manipulation of continuous molecular labels.
    """
    def __init__(self,name,
                      values,
                      context='atom',
                      sensitivity=0.1):
        self.name        = name
        self.values      = values
        self.context     = context
        self.sensitivity = sensitivity
        self.av          = np.mean(values)
        self.entropy     = self.entropy()

    def summary(self):
        """
        A continuous label summary.

        Grammatical Structure
        <noun> share a similar <property> centred around <average value>
        """
        if self.entropy < self.sensitivity:
            return f"{self.context} share a similar {self.name} centred around {self.av:.4f}"

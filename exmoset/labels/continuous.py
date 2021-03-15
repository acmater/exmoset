import numpy as np
from . import Label
from entropy_estimators import continuous

class Continuous(Label):
    """
    Class that handles the manipulation of continuous molecular labels.
    """
    def __init__(self,name,
                      values,
                      verb="share",
                      noun='Atoms',
                      sensitivity=0.1):
        self.name        = name
        self.values      = values
        self.noun        = noun
        self.sensitivity = sensitivity
        self.av          = np.mean(values)
        self.entropy     = self.entropy()

    def entropy(self,k=10,norm="euclidean",min_dist=0.001):
        """
        Estimates the entropy of a continuous distribution.
        """
        ent = continuous.get_h(self.values,k=k,norm=norm,min_dist=min_dist)
        return ent

    def summary(self,sensitivity=None):
        """
        A continuous label summary.

        Grammatical Structure
        <noun> share a similar <property> centred around <average value>
        """
        if not sensitivity:
            sensitivity = self.sensitivity
        if self.entropy < sensitivity:
            return f"{self.noun} share a similar {self.name} centred around {self.av:.4f}"

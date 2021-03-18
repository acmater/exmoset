from . import Label
import numpy as np
import matplotlib.pyplot as plt

class Multiclass(Label):
    """
    Class that handles the manipulation of multiclass molecular labels.
    """
    def __init__(self,property,
                      values,
                      verb="contain",
                      noun="Molecules",
                      sensitivity=0.1):
        self.property    = property
        self.noun        = noun
        self.verb        = verb
        self.values      = values
        self.sensitivity = sensitivity
        self.av          = np.int(np.mean(values))
        self.entropy     = self.entropy()

    def plot(self):
        plt.hist(self.values,width=1,alpha=0.5,align="mid",bins=len(np.unique(self.values)))
        plt.title(self.property)
        plt.xlim(min(self.values)+0.5,max(self.values)+0.5)
        plt.tight_layout()
        return plt.gcf()

    def summary(self,sensitivity=None,unimportant_label=False):
        """
        A multiclass label summary.

        Grammatical Structure
        <noun> <verb> <average> <property>
        """
        if not sensitivity:
            sensitivity = self.sensitivity
        if self.entropy < sensitivity:
            description = f"{self.noun} {self.verb} {self.av} {self.property}"
            return description
        elif unimportant_label:
            return f"{self.property} not meaningful"

import numpy as np
from . import Label
import matplotlib.pyplot as plt

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

    def plot(self):
        """
        Uses a customized matplotlib histogram to plot the frequency of the two class labels.
        """
        plt.hist(self.values,width=1,alpha=0.5,align="mid",bins=2)
        ax = plt.gca()
        ax.set_xticks([0.5,1.5])
        plt.xlim([0,2])
        ax.set_xticklabels(["No (0)","Yes (1)"])
        plt.title(self.property)
        plt.tight_layout()
        return plt.gcf()

    def summary(self,sensitivity=None,unimportant_label=False):
        """
        A binary label summary.

        Grammatical Structure
        <noun> <verb> <property>

        The description is generated in accordance with the different language components.
        If the value is zero, then the statement is negated using a method in the label superclass.
        """
        if not sensitivity:
            sensitivity = self.sensitivity
        if self.entropy < sensitivity:
            description = f"{self.noun} {self.verb} {self.property}"
            return description if self.av > 0.5 else self.neg(description)
        elif unimportant_label:
            return f"{self.property} not meaningful"

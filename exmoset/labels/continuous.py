import numpy as np
from . import Label
from entropy_estimators import continuous
from scipy import stats
import matplotlib.pyplot as plt

class Continuous(Label):
    """
    Class that handles the manipulation of continuous molecular labels.
    """
    def __init__(self,property,
                      values,
                      verb="share",
                      noun='Atoms',
                      sensitivity=0.1):
        self.property    = property
        self.values      = values
        self.noun        = noun
        self.sensitivity = sensitivity
        self.av          = np.mean(values)
        self.entropy     = self.entropy()

    def plot(self,bw_method=None,weights=None):
        """
        Returns a matplotlib continuous plot for this particular property.
        Density estimation is performed using the guassian_kde from scipy with default parameters.

        Parameters
        ----------
        bw_method : str, scalar, or callable, optional
            The method used to calculate the estimator bandwidth. See scipy.stats.gaussian_kde for more information.

        weights : array_like, option
            Weights of datapoints. See scipy.stats.gaussian_kde for more information.
        """
        kernel = stats.gaussian_kde(self.values,bw_method=bw_method,weights=weights)
        x = np.linspace(min(self.values),max(self.values),num=500)
        y = kernel(x)
        fig, ax = plt.subplots()
        ax.plot(x,y)
        ax.fill_between(x,y,0,alpha=0.1)

        plt.title(self.property)
        plt.tight_layout()
        return plt.gcf()

    def entropy(self,k=10,norm="euclidean",min_dist=0.001):
        """
        Estimates the entropy of a continuous distribution.
        """
        ent = continuous.get_h(self.values,k=k,norm=norm,min_dist=min_dist)
        return ent

    def summary(self,sensitivity=None,unimportant_label=False):
        """
        A continuous label summary.

        Grammatical Structure
        <noun> share a similar <property> centred around <average value>
        """
        if not sensitivity:
            sensitivity = self.sensitivity
        if self.entropy < sensitivity:
            return f"{self.noun} share a similar {self.property} centred around {self.av:.4f}"
        elif unimportant_label:
            return f"{self.property} not meaningful"

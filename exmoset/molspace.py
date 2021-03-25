import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm

from .molset import MolSet
from .molecule import Molecule
from .labels import Binary, Multiclass, Continuous

from math import log, e
from collections import namedtuple

from entropy_estimators import continuous

import matplotlib.pyplot as plt
plt.style.use("exmoset/utils/matplotlibrc")
from scipy.stats import gaussian_kde
from sklearn.preprocessing import normalize
from scipy.spatial import cKDTree
from scipy.special import gamma, digamma

class MolSpace():
    """
    A class that handles a set of chemical data and the associated MolSets. This method
    will delegate to the methods of those classes when determining things such as calculating
    mutual information for a given label.

    The most intuitive way to think of MolSpace as is providing the context in which each
    MolSet sits. If the MolSet is singular (for example an entire dataset), then MolSpace
    will contain only a single MolSet object.

    Attributes
    ----------
    fingerprints : list
            A list of property fingerprints that will be calculated for each system

    molecules : list, default=None
        A list of molecules that will be parsed into the Molecule class

    mol_converters : dict, default = {}
        A list of molecule converters that will be passed as kwargs to the Molecule object

        Example: {"rdkit" : Chem.MolFromSmiles}
        This will then be provided as keyword arguments to Molecule and given the particular mol as an argument.
        Molecule(mol, rdkit=Chem.MolFromSmiles(mol))

    significance : float, default = 0.1
        The default signifiance threshold used when calculating whether or not a particular label is significant.

    file : str, default=None
        An optional file (.csv) that will be imported by pandas and can be accessed by the fingerprints.

        # TODO make this file import flexible so multiple type formats can be specified.

    index_col : str, default=None
        The column name for which column in the file contains the base molecule representation that will then be assigned to self.Molecules

    clusters : {str : np.array(dtype=np.int)}, default={}
        A dictionary which stores cluster information. The str is used to define the name of the clustering approach
        and the numpy array is the indices which will define each MolSet cluster given a particular clustering approach.

    label_types : {str : <Label Class>}, default = {"binary" : Binary, "multiclass" : Multiclass, "continuous" : Continuous}
        A dictionary of possible label types for indexing.
    """
    def __init__(self,fingerprints,
                      molecules=None,
                      file=None,
                      mol_converters={},
                      significance=0.1,
                      index_col=None,
                      clusters={},
                      label_types = {"binary"     : Binary,
                                     "multiclass" : Multiclass,
                                     "continuous" : Continuous}):
        assert file is not None or molecules is not None, "Either a file or a list of molecules must be provided."

        if file is not None :
            assert index_col is not None, "An index column must be provided to determine what column of the file corresponds to the molecules."
            print(f"Importing {file}")
            self.df        = pd.read_csv(file,index_col=index_col)
            if molecules is None:
                molecules = self.df.index.to_numpy()

        print("Converting Molecules")
        self.Molecules = []
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            self.Molecules.append(Molecule(mol, **formats))
        self.Molecules = np.array(self.Molecules)

        print("Calculating Properties")
        labels = {fp.property : np.zeros(len(self.Molecules)) for fp in fingerprints}
        for i, molecule in enumerate(tqdm.tqdm(self.Molecules)):
            for fp in fingerprints:
                if fp.file is not None:
                    labels[fp.property][i] = fp.calculator(molecule[fp.mol_format],file=self.df)
                    # I want to make sure that this isn't passing an entire copy of the dataframe, as that would suck
                else:
                    labels[fp.property][i] = fp.calculator(molecule[fp.mol_format])

        print("Generating sets of molecules")
        self.data         = pd.DataFrame(labels)
        self.indices      = np.arange(len(self.data))
        self.fingerprints = {fp.property : fp for fp in fingerprints}
        if clusters:
            self.clusters     = {key : self.gen_clusters(value) for key, value in clusters.items()}

        else:
            self.clusters = {"Full" : MolSet(fingerprints=self.fingerprints.values(),
                                             indices=np.arange(len(self.Molecules)),
                                             context=self)}

    @staticmethod
    def c_H(values,k=10,norm="euclidean",min_dist=0.001):
        """
        Computes entropy of a continuous label distribution.

        Parameters
        ----------
        values : np.ndarray
            The array of label values.
        """
        ent = continuous.get_h(values,
                               k=k,
                               norm=norm,
                               min_dist=min_dist)
        return ent

    @staticmethod
    def d_H(values,base=None):
        """ Computes entropy of a discrete label distribution.

        Parameters
        ----------
        values : np.ndarray
            The array of label values.

        base: float
            Floating number used as base for entropy calculation.
        """
        n_labels = len(values)

        if n_labels <= 1:
            return 0

        value, counts = np.unique(values, return_counts=True)
        probs = counts / n_labels
        n_classes = np.count_nonzero(probs)
        if n_classes <= 1:
            return 0
        ent = 0.
        # Compute entropy
        base = e if base is None else base
        for i in probs:
            ent -= i * log(i, base)

        return ent

    def entropy(self,prop,set):
        """
        Determines the entropy for a particular property and set.

        Offloads the calculation to specific methods depending on the label type.

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        set : np.array(np.int)
            A numpy array of indexes that is used to determine which values should be used to calculate the entropy.
        """
        if self.fingerprints[prop].label_type == "continuous":
            values = self.data[prop].loc[set].to_numpy().reshape(-1,1)
            return self.c_H(values)
        else:
            return self.d_H(self.data[prop].loc[set])

    def mi_dd(self,prop,sets):
        """
        Mutual information for a discrete - discrete mixture.

        Adapted from - TODO (Add citation)

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.
        """
        max_val = int(max([max(self.data[prop].loc[set_]) for set_ in sets]))
        if max_val < 1:
            max_val = 1
        contingency = np.array([np.bincount(self.data[prop].loc[set_],minlength=max_val+1) for set_ in sets])

        total = 0
        N = np.sum(contingency)
        for x in range(contingency.shape[0]):
            for y in range(contingency.shape[1]):
                if contingency[x,y] == 0:
                    total += 0
                else:
                    total += (contingency[x,y]/N)*np.log2((N*contingency[x,y])/
                    (np.sum(contingency[x,:])*np.sum(contingency[:,y])))
        return total

    def mi_dc(self,prop,sets,k=3):
        """
        Calculates the mutual information for a discrete number of sets and a continuous
        label, hence mutual information - discrete continuous.

        Method adapted from - TODO (Add citation).

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        k : int, default=3
            The number of k nearest neighbours to consider for each point within the set.
        """
        assert prop in self.fingerprints.keys(), "Not a valid prop, must be a fingerprint name."
        N = len(self.data)
        full = cKDTree(self.data[prop].to_numpy().reshape(-1,1))
        Nxs = []
        ms = []
        for set_ in sets:
            label_data = self.data[prop].loc[set_].to_numpy().reshape(-1,1)
            N_xi = len(label_data)
            labeltree = cKDTree(label_data)
            distances, _ = labeltree.query(label_data,k=k)
            m_i = full.query_ball_point(label_data,distances[:,-1],return_length=True)
            Nxs.append(digamma(N_xi))
            ms.append(np.mean(digamma(m_i)))
        return digamma(N) + digamma(k) - np.mean(ms) - np.mean(Nxs)

    def plot_kdes(self,prop,sets):
        """
        Generates the density plots for a particular property given different integers sets.

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        TODO - Added legend to KDE plot.
        """
        assert prop in self.fingerprints.keys(), "Not a valid prop, must be a fingerprint name."
        data = [self.data[prop].loc[set_].to_numpy() for set_ in sets] #set_ to avoid shadowing set() method.
        kernels = [gaussian_kde(datum) for datum in data]
        positions = np.linspace(min([min(x) for x in data]),max([max(x) for x in data]),1000)
        for kernel in kernels:
            plt.plot(positions,kernel(positions))
            plt.fill_between(positions,kernel(positions),alpha=0.3)
        return plt.gcf()

    def plot_entropy(self,set,props="all"):
        """
        Helper function to plot the entropy (and its associated sensitivity) for each property of interest within the code.

        Parameters
        ----------
        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        prop : str, default="all"
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.
            If "all" is provided it will use every fingerprint in the code.

        Returns
        -------
        plt.fig
            The matplotlib figure that can then be plotted using plt.show() on the subsequent line.
        """
        if props == "all":
            props = self.fingerprints.keys()

        Summary = namedtuple("Summary",["val","entropy","sensitivity"])
        label_dict_sorted    = {key : Summary(np.mean(self.data[key].loc[set]),self.entropy(key,set),self.fingerprints[key].sensitivity) for key in self.fingerprints}
        label_dict_sorted    = {key : val for key, val in sorted(label_dict_sorted.items(), key=lambda item: item[1][1])}

        plt.bar(range(len(label_dict_sorted)), [x.entropy for x in label_dict_sorted.values()], align="center",alpha=0.5,color="r")
        labels = [self.fingerprints[key].summary(*label_dict_sorted[key],unimportant_label=True) for key in label_dict_sorted]
        colors = ["#404040" if "not meaningful" in label else "#FFFFFF" for label in labels ]
        plt.xticks(range(len(label_dict_sorted)), labels, rotation=45, ha='right')
        for label, color in zip(plt.gca().get_xticklabels(),colors):
            label.set_color(color)
        plt.ylabel("Entropy",fontsize=30)
        plt.plot(range(len(label_dict_sorted)), [x[2] for x in label_dict_sorted.values()], dashes=[6,2],color='w')
        plt.tight_layout()
        return plt.gcf()

    def calc_vector(self,set):
        """
        Represents the meaningful properties as a vector. Uses numpy's masking feature to only include meaningful
        values within it.

        Returns
        -------
        np.MaskedArray
            This array masks all non important values and uses the averages in other positions. Is used in
            generate outliers to determine the distance for each label.

        np.MaskedArray(Struct)
            A masked struct that has each value accessible by the fingerprints assocaited name. Is used for various dunder methods
            such as __and__.
        """
        data = self.data.loc[set]
        vector = np.zeros((len(self.fingerprints),))
        mask   = np.zeros((len(self.fingerprints),))
        for i,prop in enumerate(self.fingerprints.keys()):
            vector[i] = np.mean(data[prop])
            if self.entropy(prop,set) < self.fingerprints[prop].sensitivity:
                mask[i] = 0
            else:
                mask[i] = 1

        struct = np.array(tuple(vector), dtype=[(key,"<f4") for key in self.fingerprints.keys()])

        return ma.array(vector, mask=mask), ma.array(struct, mask=tuple(mask))

    def gen_clusters(self,indices):
        """
        Uses an array of indexes to generate the clusters by either breaking them into individual
        classes, or by creating the set and its complement.

        Parameters
        ----------
        indices : np.array[np.int]
            The numpy array of integer indexes that will be used to assign membership to different sets.
        """
        if len(indices) == len(self.indices):
            return {val : np.where(indices==val)[0] for val in np.unique(indices)}
        elif len(indices) < len(self.indices):
            return {"Set" : indices, "Complement" : np.setdiff1d(self.indices,indices)}
        else:
            raise ValueError("Invalid index array.")

        return clusters

    def query(self,label,condition):
        """
        This method is intended to allow to user to probe the dataset in question.
        """
        pass

    def __getitem__(self,idx):
        return self.clusters[idx]

import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm
import multiprocessing as mp

from .molecule import Molecule
from .fingerprints import gen_fp

from math import log, e
from collections import namedtuple

from entropy_estimators import continuous

import matplotlib.pyplot as plt
plt.style.use("exmoset/utils/light")

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

    pandas_func : function , default=pd.read_csv
        The pandas function that will be used to read the data file in.

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
                      pandas_func=pd.read_csv,
                      mol_converters={},
                      significance=0.1,
                      index_col=None,
                      clusters={},
                      label_types = ["binary","multiclass","continuous"]):

        assert file is not None or molecules is not None, "Either a file or a list of molecules must be provided."

        if file is not None :
            assert index_col is not None, "An index column must be provided to determine what column of the file corresponds to the molecules."
            print(f"Importing {file}")
            self.data     = pandas_func(file,index_col=index_col)
            if molecules is None:
                molecules = self.data.index.to_numpy()
            indexing_fingerprints = {prop : gen_fp(prop) for prop in list(self.data)}
        else:
            self.data = pd.DataFrame()
            indexing_fingerprints = {}

        print("Converting Molecules")
        mols = []
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            mols.append(Molecule(mol, **formats))
        # I may want to replace this section with a pandas dataframe or something.
        self.Molecules = np.array(mols)
        # Generates the necessary iters for the fingerprint methods ahead of time and caches them.
        self.mol_iters = {mol_type : [mol[mol_type] for mol in self.Molecules] for mol_type in self.Molecules[0].types}

        print("Calculating Properties")
        labels = {}
        for fp in tqdm.tqdm(fingerprints):
            print(fp.property)
            self.data[fp.property] = self.map_fingerprint(fp.calculator,self.mol_iters[fp.mol_format])

        print("Generating sets of molecules")
        self.indices      = np.arange(len(self.data))
        self.data         = self.data.set_index(self.indices)
        self.fingerprints = {**indexing_fingerprints, **{fp.property : fp for fp in fingerprints}}
        if clusters:
            self.clusters     = {key : self.gen_clusters(value) for key, value in clusters.items()}
        else:
            self.clusters = {"Full" : np.arange(len(self.data))}

    @staticmethod
    def map_fingerprint(func,mol_iter):
        return list(map(func,mol_iter))

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

    def plot(self,set,prop,title=None):
        """
        Uses a customized matplotlib histogram to plot the frequency of the two class labels.

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.
        """
        assert prop in self.fingerprints.keys(), "Not a valid property."

        values = self.data[prop].loc[set]
        label_type = self.fingerprints[prop].label_type

        if label_type == "binary":
            plt.hist(values,width=1,alpha=0.5,align="mid",bins=2)
            ax = plt.gca()
            ax.set_xticks([0.5,1.5])
            plt.xlim([0,2])
            ax.set_xticklabels(["No (0)","Yes (1)"])

        elif label_type == "multiclass":
            plt.hist(values,alpha=0.5,align="mid")
            plt.xlim(min(values)-0.5,max(values)+0.5)

        elif label_type == "continuous":
            kernel = gaussian_kde(values)
            x = np.linspace(min(values),max(values),num=500)
            y = kernel(x)
            fig, ax = plt.subplots()
            ax.plot(x,y)
            ax.fill_between(x,y,0,alpha=0.1)

        if title is None:
            title = prop
        plt.title(title)
        plt.tight_layout()
        return plt.gcf()


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
        colors = ["#404040" if "not meaningful" in label else "#3BB2E2" for label in labels ]
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

    def get_outliers(self,set):
        """
        Function that compares the meaningul vector description of a class and identifiers species that are not correctly
        described by it and returns the indices of these species.

        Returns
        -------
        np.BooleanArray
            An array that can be used to index self.Molecules to return the outlier species.
        """
        return np.where(np.sum((self.data.loc[set].to_numpy() - self.calc_vector(set)[0]),axis=1) > 0.5) # Need to update distance formulation


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

    def add_cluster(self,cluster):
        """
        Helper function to update the clusters.

        Parameters
        ----------
        cluster : (str, np.array)
            Takes in a tuple with the string key denoting the clusters' name and the nd.array its associated indices.
        """
        self.clusters[cluster[0]] = self.gen_clusters(cluster[1])

    def add_fingerprint(self,fp):
        """
        Helper function to add a new fingerprint to the fingerprints object and map its values into self.data.

        Parameters
        ----------
        fp : Fingerprint
            The fingerprint that will be used, must be provided using the fingerprint API.
        """
        self.fingerprints[fp.property] = fp
        self.data[fp.property] = list(self.map_fingerprint(fp.calculator, self.mol_iters[fp.mol_format]))

    def __repr__(self):
        return self.data.__repr__()

    def __getitem__(self,idx):
        return self.clusters[idx]

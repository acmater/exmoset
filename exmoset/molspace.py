import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm
import multiprocessing as mp

from .molecule import Molecule
from .fingerprints import gen_fp, template_fp

from math import log, e
from collections import namedtuple

from entropy_estimators import continuous

import matplotlib.pyplot as plt
plt.style.use("exmoset/utils/light")

from scipy.stats import gaussian_kde
from sklearn.preprocessing import normalize
from scipy.spatial import cKDTree
from scipy.special import gamma, digamma, kl_div

# Colour Configuration
meaningful = "#3BB2E2"
not_meaningful = "#808080"

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

    file : str, default=None
        An optional file (.csv) that will be imported by pandas and can be accessed by the fingerprints.

    pandas_func : function , default=pd.read_csv
        The pandas function that will be used to read the data file in.

    mol_converters : dict, default = {}
        A list of molecule converters that will be passed as kwargs to the Molecule object

        Example: {"rdkit" : Chem.MolFromSmiles}
        This will then be provided as keyword arguments to Molecule and given the particular mol as an argument.
        Molecule(mol, rdkit=Chem.MolFromSmiles(mol))

    sensitivity : float, default = 0.1
        The default sensitivity threshold used when calculating whether or not a particular label is meaningful.

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
                      sensitivity=0.1,
                      index_col=None,
                      clusters={},
                      template_fp=template_fp,
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
        self.mol_iters = {mol_type : np.array([mol[mol_type] for mol in self.Molecules]) for mol_type in self.Molecules[0].types}

        print("Calculating Properties")
        labels = {}
        for fp in tqdm.tqdm(fingerprints):
            self.data[fp.property] = self.map_fingerprint(fp.calculator,self.mol_iters[fp.mol_format])

        print("Generating Clusters")
        self.indices      = np.arange(len(self.data))
        self.sensitivity  = sensitivity
        self.data         = self.data.set_index(self.indices)
        self.fingerprints = {**indexing_fingerprints, **{fp.property : fp for fp in fingerprints}}
        if clusters:
            self.clusters     = {key : self.gen_clusters(value) for key, value in clusters.items()}
        else:
            self.clusters = {"Full" : {True : np.arange(len(self.data))}}

    @staticmethod
    def map_fingerprint(func,mol_iter):
        """ Helper function that maps a fingerprint across its associated mol_iter """
        return list(map(func,tqdm.tqdm(mol_iter)))

    @staticmethod
    def c_H(values,k=10,norm="euclidean",min_dist=0.001):
        """
        Computes entropy of a continuous label distribution.

        There are additional precautions taken to assure that values that are negative or

        Parameters
        ----------
        values : np.ndarray
            The array of label values.
        """
        if np.mean(values) < 1:
            values *= 5
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

    def entropy(self,set_idxs,prop):
        """
        Determines the entropy for a particular property and set.

        Offloads the calculation to specific methods depending on the label type.

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        set_idxs : np.array(np.int)
            A numpy array of indexes that is used to determine which values should be used to calculate the entropy.
        """
        if self.fingerprints[prop].label_type == "continuous":
            values = self.data[prop].loc[set_idxs].to_numpy().reshape(-1,1)
            return self.c_H(values)
        else:
            return self.d_H(self.data[prop].loc[set_idxs])

    def mi(self,set_, prop,set_val=0,k=3,*args,**kwargs):
        """
        Helper function that calculates mutual information using the available methods
        depending on information available in the fingerprint.
        """
        assert prop in self.fingerprints.keys(), "Not a valid prop, must be a fingerprint name."

        idxs = self[set_][set_val]
        complement = np.setdiff1d(self.indices,idxs)
        set_complement = [idxs,complement]

        if self.fingerprints[prop].label_type == "continuous":
            return self.mi_dc(set_complement,prop,k=k,*args,**kwargs)
        else:
            return self.mi_dd(set_complement,prop,*args,**kwargs)

    def mi_dd(self,sets,prop,labels=None,return_contingency=False):
        """
        Mutual information for a discrete - discrete mixture.

        Adapted from - TODO (Add citation)

        Parameters
        ----------
        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        labels : np.array, default=None
            A custom label array to use when calculating the mutual information.

        return_contingency : bool, default=False
            Whether or not to return the contingency matrix.
        """
        if labels is None:
            labels = self.data[prop]
        max_val = int(np.max([np.max(labels[set_]) for set_ in sets]))
        if max_val < 1:
            max_val = 1
        contingency = np.array([np.bincount(labels[set_],minlength=max_val+1) for set_ in sets])

        total = 0
        N = np.sum(contingency)
        for x in range(contingency.shape[0]):
            for y in range(contingency.shape[1]):
                if contingency[x,y] == 0:
                    total += 0
                else:
                    total += (contingency[x,y]/N)*np.log2((N*contingency[x,y])/
                    (np.sum(contingency[x,:])*np.sum(contingency[:,y])))

        if return_contingency:
            return total, contingency
        else:
            return total

    def mi_dc(self,sets,prop,labels=None,k=3):
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
            The number of k nearest neighbours to consider for each point within the set.\

        # Check that there are no errors in this plot.
        """
        if labels is not None:
            return self.mi_dd(sets,prop,labels)

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

    def plot(self,set_,prop,set_val=0,title=None):
        """
        Uses a customized matplotlib histogram to plot the frequency of the two class labels.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

        set_val : int, default = 0 (False)
            The value of the set to be isolated, either an integer to multiclass sets (i.e multiple clusters),
            or a numeric bool (0 for False, 1 for True)

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.
        """
        assert prop in self.fingerprints.keys(), "Not a valid property."

        idxs = self[set_][set_val]
        values = self.data[prop].loc[idxs]
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
            title = f"{set_} ({set_val})"
        plt.xlabel(prop)
        plt.title(title)
        plt.tight_layout()
        return plt.gcf()

    def plot_kdes(self,set_,prop,title=False,ax=None,*args,bins=None):
        """
        Generates the density plots for a particular property given different integers sets.

        # TODO Need to fix the summary generation so that new clusters can be used.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.
        """
        assert prop in self.fingerprints.keys(), "Not a valid prop, must be a fingerprint name."

        if ax is None:
            ax = plt.gca()

        sets = {set_val : self.data[prop].loc[idxs] for set_val,idxs in self[set_].items()} #set_ to avoid shadowing set() method.
        if bins is None:
            kernels = []
            positions = np.linspace(min([min(x) for x in sets.values()]),max([max(x) for x in sets.values()]),1000)
            for set_val in sets.keys():
                kernel = gaussian_kde(sets[set_val])
                ax.plot(positions,kernel(positions),label=self.fingerprints[set_].summary(val=set_val,entropy=0))
                ax.fill_between(positions,kernel(positions),alpha=0.3)
                kernels.append(kernel)
        else:
            for datum in sets.values():
                ax.hist(datum,bins=bins)

        if title:
            ax.title = f"{prop} Values for {set_} Sets"
        ax.legend()
        ax.set_xlabel(prop)
        return ax

    def kl_div(self,set_,prop,symmetrised=True):
        """
        Calculates the kullbach Leibler divergence between two a set and its surrounding context for a particular property.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.
        """
        assert prop in self.fingerprints.keys(), "Not a valid prop, must be a fingerprint name."

        x = self.data[prop].loc[self[set_][True]]
        y = self.data[prop].loc[self[set_][False]]

        kernel1 = gaussian_kde(x)
        kernel2 = gaussian_kde(y)

        positions = np.linspace(min(min(x),min(y)),10,1000)

        if symmetrised:
            return 0.5 * np.sum(kl_div(kernel1(positions),kernel2(positions))) + 0.5 * np.sum(kl_div(kernel2(positions),kernel1(positions)))
        else:
            return kl_div(kernel1(positions),kernel2(positions))

    def plot_mi(self,set_, props="all",ax=None,set_val=0):
        """
        Helper function to plot the mutual information (and its associated sensitivity) for each property of interest within the code.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

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

        Summary = namedtuple("Summary",["mi","sensitivity"])

        label_dict_sorted    = {key : Summary(self.mi(set_,key,set_val=set_val),self.fingerprints[key].sensitivity) for key in props}
        label_dict_sorted    = {key : val for key, val in sorted(label_dict_sorted.items(), key=lambda item: item[1][0])}

        mis = [x.mi for x in label_dict_sorted.values()]
        mis = [x if x != np.inf else 1 for x in mis]

        if ax is None:
            ax = plt.gca()

        ax.bar(range(len(label_dict_sorted)), mis, align="center",alpha=0.5,color="r")
        labels = [key for key in label_dict_sorted]
        colors = [not_meaningful if label_dict_sorted[label].sensitivity > label_dict_sorted[label].mi else meaningful for label in labels ]
        ax.set_xticks(range(len(label_dict_sorted)))
        ax.set_xticklabels(labels,rotation=45, ha='right',color=colors)
        for label, color in zip(ax.get_xticklabels(),colors):
            label.set_color(color)
        ax.set_ylabel("Mutual Information",fontsize=30)
        #plt.plot(np.linspace(-0.5,len(label_dict_sorted)+0.5,num=100), [x[1] for x in label_dict_sorted.values()], dashes=[6,2],color='k',label="Sensitivity")
        #plt.axhline(self.sensitivity,dashes=[6,2],color="k",label="Sensitivity")
        y = [x.sensitivity for x in label_dict_sorted.values()]
        xmin = [x-0.4 for x in range(len(label_dict_sorted))]
        xmax = [x+0.4 for x in range(len(label_dict_sorted))]
        ax.hlines(y,xmin,xmax,label="Sensitivity",color=meaningful)
        #ax.legend()
        #ax.title(f"Mutual Information Analysis for {set_} Set")
        return ax

    def plot_entropy(self,set_ ,props="all",set_val=True,ax=None):
        """
        Helper function to plot the entropy (and its associated sensitivity) for each property of interest within the code.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

        prop : str, default="all"
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.
            If "all" is provided it will use every fingerprint in the code.

        set_val : str, default=None
            The specific value for the set that you want to consider, defaults to True, i.e consider the set and not its complement.
            Note that for sets with multiple integers this is equivalent to defaulting to 1.

        Returns
        -------
        plt.fig
            The matplotlib figure that can then be plotted using plt.show() on the subsequent line.
        """
        if props == "all":
            props = self.fingerprints.keys()

        if ax is None:
            ax = plt.gca()

        idxs = self[set_][set_val]

        Summary = namedtuple("Summary",["val","entropy","sensitivity"])
        label_dict_sorted    = {key : Summary(np.mean(self.data[key].loc[idxs]),self.entropy(idxs,key),self.fingerprints[key].sensitivity) for key in props}
        label_dict_sorted    = {key : val for key, val in sorted(label_dict_sorted.items(), key=lambda item: item[1][1])}

        entropies = [x.entropy for x in label_dict_sorted.values()]
        entropies = [x if x != np.inf else 1 for x in entropies]

        ax.bar(range(len(label_dict_sorted)), entropies, align="center",alpha=0.5,color="r")
        labels = [self.fingerprints[key].summary(*label_dict_sorted[key],unimportant_label=True) for key in label_dict_sorted]
        colors = [not_meaningful if label_dict_sorted[key].sensitivity < label_dict_sorted[key].entropy else meaningful for key in label_dict_sorted]
        ax.set_xticks(range(len(label_dict_sorted)))
        ax.set_xticklabels(labels,rotation=45, ha='right',color=colors)
        for label, color in zip(ax.get_xticklabels(),colors):
            label.set_color(color)
        ax.set_ylabel("Entropy",fontsize=30)
        #plt.plot(np.linspace(-0.5,len(label_dict_sorted)+0.5,num=100), [x[2] for x in label_dict_sorted.values()], dashes=[6,2],color='k',label="Sensitivity")
        y = [x[2] for x in label_dict_sorted.values()]
        xmin = [x-0.4 for x in range(len(label_dict_sorted))]
        xmax = [x+0.4 for x in range(len(label_dict_sorted))]
        ax.hlines(y,xmin,xmax,label="Sensitivity",color=meaningful)
        #plt.axhline(self.sensitivity,dashes=[6,2],color="k",label="Sensitivity")
        ax.legend()
        #plt.title(f"Entropy Analysis for {set_} with value {set_val}")
        return ax

    def calc_vector(self,set_,set_val=True):
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
        data = self.data.loc[self[set_][set_val]]
        vector = np.zeros((len(self.fingerprints),))
        mask   = np.zeros((len(self.fingerprints),))
        for i,prop in enumerate(self.fingerprints.keys()):
            vector[i] = np.mean(data[prop])
            if self.entropy(self[set_][set_val],prop) < self.fingerprints[prop].sensitivity:
                mask[i] = 0
            else:
                mask[i] = 1

        struct = np.array(tuple(vector), dtype=[(key,"<f4") for key in self.fingerprints.keys()])

        return ma.array(vector, mask=mask), ma.array(struct, mask=tuple(mask))

    def get_outliers(self,set_,set_val=True):
        """
        Function that compares the meaningul vector description of a class and identifiers species that are not correctly
        described by it and returns the indices of these species.

        Returns
        -------
        np.BooleanArray
            An array that can be used to index self.Molecules to return the outlier species.
        """
        return np.where(np.sum((self.data.loc[self[set_][set_val]].to_numpy() - self.calc_vector(set_)[0]),axis=1) > 0.5) # Need to update distance formulation

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
            return {True : indices, False : np.setdiff1d(self.indices,indices)}
        else:
            raise ValueError("Invalid index array.")

        return clusters

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

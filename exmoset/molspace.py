import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm
import multiprocessing as mp

from exmoset.molecule import Molecule
from exmoset.fingerprints import gen_fp, template_fp

from math import log, e
import operator
from collections import namedtuple
from functools import reduce
from itertools import combinations

from entropy_estimators import continuous

import matplotlib.pyplot as plt
#plt.style.use("exmoset/utils/light")

from scipy.stats import gaussian_kde
from sklearn.preprocessing import normalize
from scipy.spatial import cKDTree
from scipy.special import gamma, digamma, kl_div
from scipy.optimize import dual_annealing

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

    template_fp : dict, default=template_fp (imported from template_fingerprints.py)
        The default settings for the fingerprint.

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
        self.Molecules = np.array(mols,dtype=np.object)
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
        """
        Helper function that maps a fingerprint across its associated mol_iter

        Parameters
        ----------
        mol_iter : iterable
            An iterable of molecules that the function will be mapped across.
        """
        return list(map(func,mol_iter))

    @staticmethod
    def c_H(values,k=10,norm="euclidean",min_dist=0.001):
        """
        Computes entropy of a continuous label distribution.

        This method utilises the continuous entropy estimator from Paul Broderson (https://github.com/paulbrodersen/entropy_estimators)
        which is adapted from Kraskov, A.; Stögbauer, H.; Grassberger, P. Estimating mutual information. Phys. Rev. E 2004, 69, 66138.

        Parameters
        ----------
        values : np.ndarray
            The array of label values.

        k : int, default=10
            The number of k nearest neighbours that are used to calculate the entropy.

        norm : str, default="euclidean"
            The norm to use when calculating the entropy.

        min_dist : 0.001,
            The minimum distance for the entropy calculator.
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

        Parameters
        ----------
        set_ : str
            The string identifier of the set of interest.

        prop : str
            The property of interest for the mi calculation.

        set_val : int, default=0
            The integer identifier of the subset that will define the indices of interest.
        """
        assert prop in self.fingerprints.keys() or prop == "Labels_Provided", "Not a valid prop, must be a fingerprint name."

        idxs = self[set_][set_val]
        complement = np.setdiff1d(self.indices,idxs)
        set_complement = [idxs,complement]
        if prop == "Labels_Provided":
            return self.mi_dd(set_complement,prop,*args,**kwargs)
        if self.fingerprints[prop].label_type == "continuous":
            return self.mi_dc(set_complement,prop,k=k,*args,**kwargs)
        else:
            return self.mi_dd(set_complement,prop,*args,**kwargs)

    def mi_dd(self,sets,prop,labels=None,return_contingency=False):
        """
        Mutual information for a discrete - discrete mixture.

        Adapted from https://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html

        Parameters
        ----------
        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

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

        Method adapted from - Ross, B. C. Mutual Information between Discrete and Continuous Data Sets. PLOS ONE 2014, 9, 1–5.

        Parameters
        ----------
        sets : [np.array(np.int)]
            An iterable (typically a list) of numpy index arrays.

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

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
            label_data   = self.data[prop].loc[set_].to_numpy().reshape(-1,1)
            N_xi         = len(label_data)
            labeltree    = cKDTree(label_data)
            distances, _ = labeltree.query(label_data,k=k)
            m_i          = full.query_ball_point(label_data,distances[:,-1],return_length=True)
            Nxs.append(digamma(N_xi))
            ms.append(digamma(m_i))
        ms = np.concatenate(ms)
        return digamma(N) + digamma(k) - np.mean(ms) - np.mean(Nxs)

    def plot(self,set_,prop,set_val=0,title=None):
        """
        Uses a customized matplotlib histogram to plot the frequency of the two class labels.

        Parameters
        ----------
        set_ : str
            A string to identify the cluster of interest.

        prop : str
            A property that will be analysed - must be a fingerprint name and thus a column in the dataframe.

        set_val : int, default = 0 (False)
            The value of the set to be isolated, either an integer to multiclass sets (i.e multiple clusters),
            or a numeric bool (0 for False, 1 for True)
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

    def plot_kdes(self,set_,prop,title=False,ax=None,*args,alpha=(1,0.3),bins=None,**kwargs):
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
                ax.plot(positions,kernel(positions),label=self.fingerprints[prop].summary(val=set_val,entropy=0),alpha=alpha[0])
                ax.fill_between(positions,kernel(positions),alpha=alpha[1])
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

        symmetrised : bool, default=True
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

    def plot_mi(self,set_, props="all",ax=None,set_val=0,k=3):
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

        label_dict_sorted    = {key : Summary(self.mi(set_,key,set_val=set_val,k=k),self.fingerprints[key].sensitivity) for key in props}
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

    def gen_labels(self,set_,mi_cutoff=0.005):
        """
        Method that generates labels for every subcluster within a particular clustering of the data

        Parameters
        ----------
        set_ : str
            The string identifier for a particular clustering of the space.

        mi_cutoff, float, default=0.005
            The mutual information cutoff that is used to remove unimportant labels.
        """
        summaries = {}
        for set_val,idxs in tqdm.tqdm(self[set_].items()):
            fingerprints = []
            binary_labels = []
            mis = np.zeros((len(self.fingerprints,)))
            ents = np.zeros_like(mis)
            for i,fp in enumerate(self.fingerprints.values()):
                mis[i] = self.mi(set_,fp.property,set_val=set_val)
                ents[i] = self.entropy(idxs,fp.property)
                labels = self.data[fp.property].loc[idxs]
                bounds = np.array([np.min(labels),np.max(labels)])
                mi = self.mi(set_,fp.property,set_val=set_val)
                ent = self.entropy(idxs,fp.property)

                if mi > mi_cutoff:
                    if ent < (fp.sensitivity):
                        if fp.label_type == "continuous":
                            u = np.mean(labels)
                            var = np.var(labels)
                            binary_labels.append(np.logical_and((u - 0.1*var) < self.data[fp.property], self.data[fp.property] < (u + 0.1*var)))
                            fingerprints.append(fp.summary(u,entropy=ent))
                        else:
                            binary_labels.append((self.data[fp.property] == round(np.mean(labels))).to_numpy())
                            fingerprints.append(fp.summary(round(np.mean(labels)),ent,fp.sensitivity+0.1))
                    else:
                        if fp.label_type == "binary":
                            binary_labels.append((self.data[fp.property] == round(np.mean(labels))).to_numpy())
                            fingerprints.append(fp.summary(np.mean(labels),ent,sensitivity=1))
                            continue

                        elif fp.label_type == "continuous":
                            split = float(dual_annealing(self.cost_generator(set_,fp.property,set_val=set_val),bounds=bounds.reshape(1,2),maxiter=250)["x"])
                        else:
                            iterable = np.array(range(*bounds))
                            split = iterable[np.argmin(np.apply_along_axis(self.cost_generator(set_,fp.property,set_val=set_val), 0, iterable))]

                        if split > np.mean(self.data.loc[idxs][fp.property]):
                            fingerprints.append(fp.to_binary(split,"<"))
                            binary_labels.append((self.data[fp.property] < split).to_numpy())
                        else:
                            fingerprints.append(fp.to_binary(split,">"))
                            binary_labels.append((self.data[fp.property] >= split).to_numpy())

            # To do - alter this optimizer so that it doesn't brute force it.
            if len(binary_labels) > 1:
                num_options = len(binary_labels)
                binary_labels = np.array(binary_labels)

                idxs = []
                for i in range(1,num_options+1):
                    idxs.append(list(combinations(list(range(num_options)),i))) # Generate all possible permutations
                idxs = [item for sublist in idxs for item in sublist] # Flatten the list of lists
                mut_infs = []
                for idx in idxs: # Iterate over possible label combinations
                    current_label = reduce(np.logical_and, binary_labels[tuple([idx])])
                    mut_infs.append(self.mi(set_,"Labels_Provided",set_val=set_val,labels=current_label))

                summary = [fingerprints[x] for x in idxs[np.argmax(mut_infs)]] # Isolate the fingerprint summaries for each index is the argmax collection
                summaries[set_val] = "\n\t".join(summary) + (f"\n\tMutual Information {np.max(mut_infs)}")

            elif len(binary_labels) == 1:
                summaries[set_val] = "".join(fingerprints) + (f"\n\tMutual Information {np.max(mis)}") # Don't try and iterate over label comnbinations if there is only one.

            else:
                summaries[set_val] = ("No Meaningful Label")

        return summaries

    def label_summary(self,set_,*args):
        """
        Helper function that prints out the summary labels generated by gen_labels.
        """
        summaries = self.gen_labels(set_, *args)
        for set_val,summary in enumerate(summaries.items()):
            print(f"Cluster {set_val}: \n\t{summary}")
        return None

    def cost_generator(self,set_,prop,set_val,func="mi",op=operator.ge):
        """
        Helper function that returns a function that can be used to optimise the information split.
        The cost functions return a lambda function that calculates the mutual information of a new label vector
        produced by comparing the label vector (self.data[prop]) to a cutoff point (x) that is going to be optimized.

        The floating value x determines the binary split through the multiclass or continuous label, transforming it
        into a dichotomous formulation. These lambda functions are then used to optimize the position of this split
        either by brute force (multiclass) or simulated annealing (continuous).

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

        func : str, ["mi", "ent", "entropy"], default="mi"
            The function that the arguments will be passed to when calculating the cost (usually mutual information).

        op : function _operator, default=operator.ge
            The operator that will be used to compare the decision boundary (float value) and the numpy label array.
        """
        if func == "mi":
            return np.vectorize(lambda x: -self.mi(set_, prop, set_val, labels=np.array(op(self.data[prop],float(x)))))
        else:
            return np.vectorize(lambda x: -self.ent(set_, prop, set_val, labels=np.array(op(self.data[prop],float(x)))))

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

    def sample_mols(self,set_,set_val=0,num_mols=16,molsPerRow=4,subImgSize=(200,200),**kwargs):
        """
        Function that randomly generates a sample of molecules from a particular cluster.
        """
        total_mols = len(self[set_][set_val])
        mols = self.mol_iters["rd"][self[set_][set_val]][np.random.randint(0,total_mols,size=(num_mols,))]
        return Chem.Draw.MolsToGridImage(mols,molsPerRow=molsPerRow,subImgSize=subImgSize,returnPNG=False,**kwargs) # ReturnPNG=False to allow image to be saved

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

    def __len__(self):
        return len(self.Molecules)

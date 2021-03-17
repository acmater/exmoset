import numpy as np
from rdkit import Chem
import pandas as pd
import tqdm

from .molset import MolSet
from .molecule import Molecule
from .labels import Binary, Multiclass, Continuous


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
                      mol_converters={},
                      significance=0.1,
                      file=None,
                      index_col=None,
                      clusters={},
                      label_types = {"binary"     : Binary,
                                     "multiclass" : Multiclass,
                                     "continuous" : Continuous}):
        assert  file is not None or molecules is not None, "Either a file or a list of molecules must be provided."

        if file is not None :
            print(f"Importing {file}")
            self.df        = pd.read_csv(file,index_col=index_col)
            if molecules is None:
                assert index_col is not None, "An index column must be provided to determine what column of the file corresponds to the molecules."
                molecules = self.df.index.to_numpy()

        self.Molecules = []
        print("Converting Molecules")
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            self.Molecules.append(Molecule(mol, **formats))
        self.Molecules = np.array(self.Molecules)

        self.prop_values = np.zeros((len(fingerprints),len(self.Molecules)))

        print("Calculating Properties")
        self.labels = {fp.property : np.zeros(len(self.Molecules)) for fp in fingerprints}
        for i, molecule in enumerate(tqdm.tqdm(self.Molecules)):
            for fp in fingerprints:
                if fp.file is not None:
                    self.labels[fp.property][i] = fp.calculator(molecule[fp.mol_format],file=self.df)
                    # I want to make sure that this isn't passing an entire copy of the dataframe, as that would suck
                else:
                    self.labels[fp.property][i] = fp.calculator(molecule[fp.mol_format])

        self.fingerprints = fingerprints
        self.clusters     = {key : self.gen_clusters(val) for key, val in clusters.items()}

    def gen_clusters(self,indices,fingerprints=None):
        if fingerprints is None:
            fingerprints = self.fingerprints
        clusters = []
        for val in np.unique(indices):
            clusters.append(MolSet(fingerprints=fingerprints,
                                   indices=np.where(indices==val)[0],
                                   context=self))
        return clusters

    def query(self):
        """
        This method is intended to allow to user to probe the dataset in question.
        """
        pass

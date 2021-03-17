import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm

from ..molecule import Molecule
from ..labels import Binary, Multiclass, Continuous

class MolSet():
    """
    A molecular set to be analysed.

    Attributes
    ----------
    molecules : list
        A list of molecules that will be parsed into the Molecule class

    properties : list
        A list of property fingerprints that will be calculated for each system

    indices : np.array(dtype=np.int), default=None
        A numpy array of indexes which define the molset within the context of a larger molecular space.

    significance : float, default = 0.1
        The default signifiance threshold used when calculating whether or not a particular label is significant.

    file : str, default=None
        An optional file (dataframe) that will be imported by pandas and can be accessed by the fingerprints.

    label_types : {str : <Label Class>}, default = {"binary" : Binary, "multiclass" : Multiclass, "continuous" : Continuous}
        A dictionary of possible label types for indexing.
    """
    def __init__(self,molecules=None,
                      fingerprints=None,
                      mol_converters={},
                      significance=0.1,
                      file=None,
                      indices=None,
                      label_types = {"binary"     : Binary,
                                     "multiclass" : Multiclass,
                                     "continuous" : Continuous}):

        if indices is not None:
            self.indices = indices
        else:

            if mol_converters:
                self.Molecules = []
                print("Converting Molecules")
                for mol in tqdm.tqdm(molecules):
                    formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
                    self.Molecules.append(Molecule(mol, **formats))
                # Step below is to allow array indexing.
                self.Molecules = np.array(self.Molecules)
            else:
                self.Molecules = molecules



            if file is not None:
                print(f"Importing {file}")
                df     = pd.read_csv(file,index_col="SMILES")
                sub_df = df.loc[[mol.smiles for mol in self.Molecules]]

            self.prop_values = np.zeros((len(fingerprints),len(self.Molecules)))

            print("Calculating Properties")
            for i,molecule in enumerate(tqdm.tqdm(self.Molecules)):
                for j,fp in enumerate(fingerprints):
                    if fp.file is not None:
                        self.prop_values[j,i] = fp.calculator(molecule[fp.mol_format],file=sub_df)
                    else:
                        self.prop_values[j,i] = fp.calculator(molecule[fp.mol_format])

            self.label_dict = {}
            for i,fp in enumerate(fingerprints):
                self.label_dict[fp.property] = label_types[fp.label_type](fp.property,self.prop_values[i],fp.verb,fp.context)

            self.significance       = significance
            self.vector,self.struct = self.calc_vector()

    def calc_vector(self):
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
        vector = np.zeros((len(self.label_dict),))
        mask   = np.zeros((len(self.label_dict),))
        test = []
        for i,prop in enumerate(self.label_dict.values()):
            vector[i] = prop.av
            if prop.entropy < prop.sensitivity:
                mask[i] = 0
            else:
                mask[i] = 1

        struct = np.array(tuple(vector), dtype=[(key,"<f4") for key in self.label_dict.keys()])

        return ma.array(vector, mask=mask), ma.array(struct, mask=tuple(mask))

    def get_outliers(self):
        # This is clunky af
        """
        Function that compares the meaningul vector description of a class and identifiers species that are not correctly
        described by it and returns the indices of these species.

        Returns
        -------
        np.BooleanArray
            An array that can be used to index self.Molecules to return the outlier species.
        """
        full_mask = np.zeros((len(self.vector.mask),len(self)))
        for i in range(len(self)):
            full_mask[:,i] = self.vector.mask
        masked_vals = ma.array(self.prop_values,mask=full_mask)
        return np.where(np.sum((self.vector.reshape(-1,1) - masked_vals),axis=0) != 0) # Need to update distance formulation

    def plot_entropy(self):
        """
        Helper function to plot the entropy (and its associated sensitivity) for each fingerprint within the code.

        Returns
        -------
        plt.fig
            The matplotlib figure that can then be plotted using plt.show() on the subsequent line.
        """
        import matplotlib.pyplot as plt
        from matplotlib import cm
        label_dict_sorted    = {key : (self.label_dict[key].entropy,self.label_dict[key].sensitivity) for key in self.label_dict}
        label_dict_sorted    = {key : val for key, val in sorted(label_dict_sorted.items(), key=lambda item: item[1])}
        plt.style.use("exmoset/utils/matplotlibrc")
        plt.bar(range(len(label_dict_sorted)), [x[0] for x in label_dict_sorted.values()], align="center",alpha=0.5,color="r")
        labels = [self.label_dict[key].summary(unimportant_label=True) for key in label_dict_sorted]
        colors = ["#404040" if "not meaningful" in label else "#FFFFFF" for label in labels ]
        plt.xticks(range(len(label_dict_sorted)), labels, rotation=45, ha='right')
        for label, color in zip(plt.gca().get_xticklabels(),colors):
            label.set_color(color)
        plt.ylabel("Entropy",fontsize=30)
        plt.plot(range(len(label_dict_sorted)), [x[1] for x in label_dict_sorted.values()], dashes=[6,2],color='w')
        plt.tight_layout()
        return plt.gcf()

    def __len__(self):
        """
        Helper function to determine length. Length is entirely determined by how many molecules are in the set.
        """
        return len(self.Molecules)

    def __str__(self):
        header  = ["Subset Description"]
        labels = [label.summary() for label in self.label_dict.values()]
        results = ["\t{:<60}".format(label) for label in labels if label is not None]

        if len(results) == 0:
            return "No Discernible pattern"
        else:
            final = "\n".join(header + results)
            return final

    def __and__(self, other):
        """
        Overrides the boolean `and` operation to identify meaningful labels that are shared between two MolSets.
        """
        assert isinstance(other, MolSet)
        keys = []
        for i,key in enumerate(self.label_dict.keys()):
            if self.struct[key] == other.struct[key]:
                keys.append(key)
        return keys

    def __getitem__(self,idx):
        return self.Molecules[idx]

    def __bool__(self):
        # Need to update this bool method
        if self.indices is not None:
            return True
        else:
            return False

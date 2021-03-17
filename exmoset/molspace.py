class MolSpace():
    """
    A class that handles a set of chemical data and the associated MolSets. This method
    will delegate to the methods of those classes when determining things such as calculating
    mutual information for a given label.

    Attributes
    ----------
    molecules : list
        A list of molecules that will be parsed into the Molecule class

    properties : list
        A list of property fingerprints that will be calculated for each system

    mol_converters : dict, default = {}
        A list of molecule converters that will be passed as kwargs to the Molecule object

        Example: {"rdkit" : Chem.MolFromSmiles}
        This will then be provided as keyword arguments to Molecule and given the particular mol as an argument.
        Molecule(mol, rdkit=Chem.MolFromSmiles(mol))

    significance : float, default = 0.1
        The default signifiance threshold used when calculating whether or not a particular label is significant.

    file : str, default=None
        An optional file (dataframe) that will be imported by pandas and can be accessed by the fingerprints.

    clusters : np.array(dtype=np.int), default=None
        A numpy array of indexes which will define each MolSet cluster given a particular clustering approach.

    label_types : {str : <Label Class>}, default = {"binary" : Binary, "multiclass" : Multiclass, "continuous" : Continuous}
        A dictionary of possible label types for indexing.
    """
    def __init__(self,molecules,
                      fingerprints,
                      mol_converters={},
                      significance=0.1,
                      file=None,
                      clusters=None
                      label_types = {"binary"     : Binary,
                                     "multiclass" : Multiclass,
                                     "continuous" : Continuous}):
        self.Molecules = []
        print("Converting Molecules")
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            self.Molecules.append(Molecule(mol, **formats))
        self.Molecules = np.array(self.Molecules)


    def query(self):
        """
        This method is intended to allow to user to probe the dataset in question.
        """
        pass

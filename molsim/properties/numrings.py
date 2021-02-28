import numpy as np
from rdkit import Chem
from property import Property


class NumRings(Property):
    """
    How many rings the molecule contains
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values = self.calc_property(molecules)

    @staticmethod
    def calc_property(molecules):
        return [x.GetRingInfo().NumRings() for x in molecules]


if __name__ == "__main__":
    a = "CCNC"
    b = "c1ccccc1"
    mols = list(map(Chem.MolFromSmiles,[a,b]))
    print(NumRings(mols).values)

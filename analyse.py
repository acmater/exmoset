from exmoset.data import *
from exmoset.fingerprints import *
from exmoset import MolSet

from rdkit import Chem

if __name__ == "__main__":
    fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints

    analysis = MolSet(molecules7,
                    fingerprints = fingerprints,
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    significance=0.1,
                    file="exmoset/data/QM9_Data.csv")

    analysis2 = MolSet(molecules4,
                fingerprints = fingerprints,
                mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                significance=0.1,
                file="exmoset/data/QM9_Data.csv")

    import matplotlib.pyplot as plt
    fig = analysis.plot_entropy()
    #plt.show()
    print(analysis[analysis.get_outliers()])
    print(analysis & analysis2)
    print(analysis)

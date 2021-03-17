from exmoset.data import *
from exmoset.fingerprints import *
from exmoset.molset import MolSet
from exmoset import MolSpace


from rdkit import Chem

if __name__ == "__main__":
    fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints

    analysis = MolSet(molecules3,
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
    plt.savefig("Molecules6_Entropy.png",transparent=True,dpi=300)
    #plt.show()
    print(analysis[analysis.get_outliers()])
    print(analysis & analysis2)
    print(analysis)

    #### MolSpace testing
    space = MolSpace(molecules3,
                    fingerprints = fingerprints,
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    significance=0.1,
                    file="exmoset/data/QM9_Data.csv")

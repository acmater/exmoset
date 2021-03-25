from exmoset.data import *
from exmoset.fingerprints import *
from exmoset.molset import MolSet
from exmoset import MolSpace

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

if __name__ == "__main__":
    fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints

    """analysis = MolSet(molecules3,
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
    print(analysis)"""

    #### MolSpace testing
    space = MolSpace(fingerprints = fingerprints,
                    molecules=molecules3,
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    significance=0.1,
                    file="exmoset/data/QM9_Data.csv",
                    index_col="SMILES",
                    clusters={"Full" : np.concatenate([np.zeros((200,)).reshape(-1,1),np.ones((382,)).reshape(-1,1)])})

    print(space.mi_dd("Aromatic",space.clusters["Full"].values()))
    print(space.mi_dc("Dipole Moment",space.clusters["Full"].values()))
    #fig = space.plot_kdes("Dipole Moment",[space.clusters["Full"][0],space.clusters["Full"][1],np.arange(200,300)])
    #plt.show()

    print(space.entropy("Aromatic",np.array([0,1,2,3,4])))
    space.plot_entropy(np.array([0,1,2,3,4]))
    plt.savefig("Test.png",transparent=True)

    print(space.calc_vector(np.array([0,1,2,3,4,5])))

    del fingerprints[3] # Remove the Dipole Moment fingerprint
    space2 = MolSpace(fingerprints=fingerprints,
                      file="exmoset/data/FreeSolv.csv",
                      mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                      index_col="SMILES",
                      clusters={})
    import matplotlib.pyplot as plt
    space2.clusters["Full"]["Atoms"].entropy

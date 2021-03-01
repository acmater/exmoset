import numpy as np
from rdkit import Chem

from properties import *





molecules = ['COC(CN)C#N', 'CC(C)(OC=N)C#C', 'CC(C)(OC=N)C#N', 'CN(CC#C)CC#N', 'CN(CC#N)CC#N', 'CC(C#C)C(C)C#C', 'CC(C#C)C(C)C#N', 'CC(C#C)C(N)C#N', 'CC(C#C)C(O)C#C', 'CC(C#C)C(O)C#N', 'CC(C#N)C(C)C#N', 'CC(C#N)C(N)C#C', 'CC(C#N)C(N)C#N', 'CC(C#N)C(O)C#C', 'CC(C#N)C(O)C#N', 'NC(C#C)C(O)C#N', 'NC(C#N)C(O)C#C', 'NC(C#N)C(O)C#N', 'OC(C#C)C(O)C#C', 'OC(C#C)C(O)C#N', 'OC(C#N)C(O)C#N', 'CC(NCC#C)C#N', 'CC(NCC#N)C#C', 'CC(NCC#N)C#N', 'CNC(C#N)C#CC', 'CNC(CC#C)C#N', 'CNC(CC#N)C#C', 'CNC(CC#N)C#N', 'CCNC(C#C)C#N', 'CCNC(C#N)C#N', 'CC(C)C(C#C)C#C', 'CC(C)C(C#C)C#N', 'CC(C)C(C#N)C#N', 'CC(O)C(C#C)C#C', 'CC(O)C(C#C)C#N', 'C#CCC#CCC1CN1', 'N#CCC#CCC1CN1', 'C#CCCC#CC1CN1', 'N#CCCC#CC1CN1', 'CC#CC(C#C)C(C)C', 'CC#CC(C#C)C(C)O', 'CC#CC(C#N)C(C)C', 'CC#CC(C#N)C(C)O', 'CC#CC(C)(C)OC=N', 'CC#CC(C)C(C)C#C', 'CC#CC(C)C(C)C#N', 'CC#CC(C)C(N)C#N', 'CC#CC(C)C(O)C#C', 'CC#CC(C)C(O)C#N', 'CC#CC(N)C(C)C#N', 'CC#CC(N)C(N)C#N', 'CC#CC(N)C(O)C#N', 'CC#CC(O)C(C)C#C', 'CC#CC(O)C(C)C#N', 'CC#CC(O)C(N)C#N', 'CC#CC(O)C(O)C#C', 'CC#CC(O)C(O)C#N', 'CC#CC(C)NCC#N', 'CC#CCN(C)CC#N', 'CC#CCNC(C)C#N', 'CC(=O)OC(C)(C)C#C', 'CC(=O)OC(C)(C)C#N', 'CC(C)(OC(N)=O)C#C', 'CC(C)(OC(N)=O)C#N', 'CC(C)(O)C(C#C)C#C', 'CC(C)(CC#C)OC=N', 'CC(C)(CC#N)OC=N', 'CC(C)(CC#N)OC=O', 'CC(C)(NCC#C)C#N', 'CC(C)(NCC#N)C#C', 'CC(C)(NCC#N)C#N', 'CC(C)(COC=N)C#N', 'CC(C)C(CC#C)C#C', 'CC(C)C(CC#C)C#N', 'CC(C)C(CC#N)C#C', 'CC(C)C(CC#N)C#N', 'CC(N)C(CC#N)C#N', 'CC(O)C(CC#C)C#C', 'CC(O)C(CC#C)C#N', 'CC(O)C(CC#N)C#C', 'CC(O)C(CC#N)C#N', 'CN(C)C(CC#C)C#N', 'CN(C)C(C#C)CC#N', 'CN(C)C(CC#N)C#N', 'CC(C)CC(C#C)C#C', 'CC(C)CC(C#C)C#N', 'CC(C)CC(C#N)C#N', 'CC(C)OC(C#C)C#N', 'CC(C)OC(C#N)C#N', 'CC(O)CC(C#C)C#C', 'CC(O)CC(C#C)C#N', 'CN(C)CC(C#C)C#N', 'CC(C#C)C(O)CC#C', 'CC(C#C)C(O)CC#N', 'CC(C#C)N(C)CC#N', 'CC(C#N)C(N)CC#N', 'CC(C#N)C(O)CC#C', 'CC(C#N)C(O)CC#N', 'CC(C#N)N(C)CC#C', 'CC(C#N)N(C)CC#N', 'CC(CC#C)C(C)C#C', 'CC(CC#C)C(C)C#N', 'CC(CC#C)C(N)C#N', 'CC(CC#C)C(O)C#C', 'CC(CC#C)C(O)C#N', 'CC(CC#N)C(C)C#C', 'CC(CC#N)C(C)C#N', 'CC(CC#N)C(N)C#N', 'CC(CC#N)C(O)C#C', 'CC(CC#N)C(O)C#N', 'NC(C#N)C(O)CC#C', 'NC(C#N)C(O)CC#N', 'NC(CC#C)C(N)C#N', 'NC(CC#C)C(O)C#N', 'NC(CC#N)C(N)C#N', 'NC(CC#N)C(O)C#C', 'NC(CC#N)C(O)C#N', 'OC(CC#C)C(O)C#C', 'OC(CC#C)C(O)C#N', 'OC(CC#N)C(O)C#C', 'OC(CC#N)C(O)C#N', 'CC(CC(C)C#C)C#C', 'CC(CC(C)C#N)C#C', 'CC(CC(C)C#N)C#N', 'CC(CC(N)C#N)C#C', 'CC(CC(N)C#N)C#N', 'CC(CC(O)C#C)C#C', 'CC(CC(O)C#C)C#N', 'CC(CC(O)C#N)C#C', 'CC(CC(O)C#N)C#N', 'CC(OC(C)C#C)C#C', 'CC(OC(C)C#N)C#C', 'CC(OC(C)C#N)C#N', 'NC(CC(N)C#N)C#N', 'NC(CC(O)C#C)C#N', 'NC(CC(O)C#N)C#N', 'OC(CC(O)C#C)C#C', 'OC(CC(O)C#N)C#C', 'OC(CC(O)C#N)C#N', 'CC(CC#C)NCC#N', 'CC(CC#N)NCC#C', 'CC(CC#N)NCC#N', 'CN(CCC#C)CC#N', 'CN(CCC#N)CC#C', 'CN(CCC#N)CC#N', 'CC(CNCC#C)C#N', 'CC(CNCC#N)C#C', 'CC(CNCC#N)C#N', 'CC(NCCC#C)C#N', 'CC(NCCC#N)C#C', 'CC(NCCC#N)C#N', 'NC(CNCC#C)C#N', 'NC(CNCC#N)C#N', 'OC(CNCC#C)C#N', 'OC(CNCC#N)C#C', 'OC(CNCC#N)C#N', 'CCC#CC(NC)C#N', 'CNC(C#N)C#CCO', 'CC(C#C)C(CO)C#C', 'CC(C#C)C(CO)C#N', 'CC(C#N)C(CN)C#N', 'CC(C#N)C(CO)C#C', 'CC(C#N)C(CO)C#N', 'CCC(C#C)C(C)C#C', 'CCC(C#C)C(C)C#N', 'CCC(C#C)C(N)C#N', 'CCC(C#C)C(O)C#C', 'CCC(C#C)C(O)C#N', 'CCC(C#N)C(C)C#C', 'CCC(C#N)C(C)C#N', 'CCC(C#N)C(N)C#C', 'CCC(C#N)C(N)C#N', 'CCC(C#N)C(O)C#C', 'CCC(C#N)C(O)C#N', 'COC(C#C)C(C)C#C', 'COC(C#C)C(C)C#N', 'COC(C#C)C(N)C#N', 'COC(C#C)C(O)C#C', 'COC(C#C)C(O)C#N', 'COC(C#N)C(C)C#C', 'COC(C#N)C(C)C#N', 'COC(C#N)C(N)C#C', 'COC(C#N)C(N)C#N', 'COC(C#N)C(O)C#C', 'COC(C#N)C(O)C#N', 'NC(C#C)C(CO)C#N', 'NC(C#N)C(CO)C#C', 'NC(C#N)C(CO)C#N', 'NCC(C#N)C(N)C#N', 'NCC(C#N)C(O)C#N', 'OCC(C#C)C(O)C#C', 'OCC(C#C)C(O)C#N', 'OCC(C#N)C(O)C#C', 'OCC(C#N)C(O)C#N', 'CCC(C)(OC=N)C#C', 'CCC(C)(OC=N)C#N', 'CC(CO)C(C#C)C#C', 'CC(CO)C(C#C)C#N', 'CCC(C)C(C#C)C#C', 'CCC(C)C(C#C)C#N', 'CCC(C)C(C#N)C#N', 'CCC(O)C(C#C)C#C', 'CCC(O)C(C#C)C#N', 'COC(C)C(C#C)C#C', 'COC(C)C(C#C)C#N', 'NC(CO)C(C#C)C#N', 'OCC(O)C(C#C)C#C', 'OCC(O)C(C#C)C#N', 'CNC(CC#N)C#CC', 'CCN(CC#C)CC#N', 'CCN(CC#N)CC#N', 'CNC(CC#N)CC#N', 'CNC(CC#CC)C#N', 'CCC(NCC#C)C#N', 'CCC(NCC#N)C#C', 'CCC(NCC#N)C#N', 'CNC(CCC#C)C#N', 'CNC(CCC#N)C#N', 'NCC(NCC#C)C#N', 'NCC(NCC#N)C#N', 'OCC(NCC#C)C#N', 'OCC(NCC#N)C#C', 'OCC(NCC#N)C#N', 'CCNC(C#N)C#CC', 'CCNC(CC#C)C#N', 'CCNC(CC#N)C#C', 'CCNC(CC#N)C#N', 'CNCC(CC#N)C#N', 'CCCNC(C#C)C#N', 'CCCNC(C#N)C#N', 'OCCNC(C#C)C#N', 'OCCNC(C#N)C#N']

properties = [Aromatic,NumRings,NumAtoms]



class Similarity_Analysis():
    def __init__(self,molecules,properties):
        self.molecules = molecules
        self.properties = [prop(molecules) for prop in properties]
        self.labels = [prop.values for prop in self.properties]
        self.mean = np.mean(np.vstack(self.labels),axis=1)
        self.variance = np.var(np.vstack(self.labels),axis=1)
        self.entropy = [x.entropy() for x in self.properties]



if __name__ == "__main__":
    print(Similarity_Analysis(molecules,properties).mean)
    print(Similarity_Analysis(molecules,properties).labels)
    print(Similarity_Analysis(molecules,properties).entropy)

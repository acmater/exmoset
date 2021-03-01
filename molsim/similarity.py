import numpy as np
from rdkit import Chem
import pprint
from collections import namedtuple

from molecule import *
from atom import *


molecules = ['COC(CN)C#N', 'CC(C)(OC=N)C#C', 'CC(C)(OC=N)C#N', 'CN(CC#C)CC#N', 'CN(CC#N)CC#N', 'CC(C#C)C(C)C#C', 'CC(C#C)C(C)C#N', 'CC(C#C)C(N)C#N', 'CC(C#C)C(O)C#C', 'CC(C#C)C(O)C#N', 'CC(C#N)C(C)C#N', 'CC(C#N)C(N)C#C', 'CC(C#N)C(N)C#N', 'CC(C#N)C(O)C#C', 'CC(C#N)C(O)C#N', 'NC(C#C)C(O)C#N', 'NC(C#N)C(O)C#C', 'NC(C#N)C(O)C#N', 'OC(C#C)C(O)C#C', 'OC(C#C)C(O)C#N', 'OC(C#N)C(O)C#N', 'CC(NCC#C)C#N', 'CC(NCC#N)C#C', 'CC(NCC#N)C#N', 'CNC(C#N)C#CC', 'CNC(CC#C)C#N', 'CNC(CC#N)C#C', 'CNC(CC#N)C#N', 'CCNC(C#C)C#N', 'CCNC(C#N)C#N', 'CC(C)C(C#C)C#C', 'CC(C)C(C#C)C#N', 'CC(C)C(C#N)C#N', 'CC(O)C(C#C)C#C', 'CC(O)C(C#C)C#N', 'C#CCC#CCC1CN1', 'N#CCC#CCC1CN1', 'C#CCCC#CC1CN1', 'N#CCCC#CC1CN1', 'CC#CC(C#C)C(C)C', 'CC#CC(C#C)C(C)O', 'CC#CC(C#N)C(C)C', 'CC#CC(C#N)C(C)O', 'CC#CC(C)(C)OC=N', 'CC#CC(C)C(C)C#C', 'CC#CC(C)C(C)C#N', 'CC#CC(C)C(N)C#N', 'CC#CC(C)C(O)C#C', 'CC#CC(C)C(O)C#N', 'CC#CC(N)C(C)C#N', 'CC#CC(N)C(N)C#N', 'CC#CC(N)C(O)C#N', 'CC#CC(O)C(C)C#C', 'CC#CC(O)C(C)C#N', 'CC#CC(O)C(N)C#N', 'CC#CC(O)C(O)C#C', 'CC#CC(O)C(O)C#N', 'CC#CC(C)NCC#N', 'CC#CCN(C)CC#N', 'CC#CCNC(C)C#N', 'CC(=O)OC(C)(C)C#C', 'CC(=O)OC(C)(C)C#N', 'CC(C)(OC(N)=O)C#C', 'CC(C)(OC(N)=O)C#N', 'CC(C)(O)C(C#C)C#C', 'CC(C)(CC#C)OC=N', 'CC(C)(CC#N)OC=N', 'CC(C)(CC#N)OC=O', 'CC(C)(NCC#C)C#N', 'CC(C)(NCC#N)C#C', 'CC(C)(NCC#N)C#N', 'CC(C)(COC=N)C#N', 'CC(C)C(CC#C)C#C', 'CC(C)C(CC#C)C#N', 'CC(C)C(CC#N)C#C', 'CC(C)C(CC#N)C#N', 'CC(N)C(CC#N)C#N', 'CC(O)C(CC#C)C#C', 'CC(O)C(CC#C)C#N', 'CC(O)C(CC#N)C#C', 'CC(O)C(CC#N)C#N', 'CN(C)C(CC#C)C#N', 'CN(C)C(C#C)CC#N', 'CN(C)C(CC#N)C#N', 'CC(C)CC(C#C)C#C', 'CC(C)CC(C#C)C#N', 'CC(C)CC(C#N)C#N', 'CC(C)OC(C#C)C#N', 'CC(C)OC(C#N)C#N', 'CC(O)CC(C#C)C#C', 'CC(O)CC(C#C)C#N', 'CN(C)CC(C#C)C#N', 'CC(C#C)C(O)CC#C', 'CC(C#C)C(O)CC#N', 'CC(C#C)N(C)CC#N', 'CC(C#N)C(N)CC#N', 'CC(C#N)C(O)CC#C', 'CC(C#N)C(O)CC#N', 'CC(C#N)N(C)CC#C', 'CC(C#N)N(C)CC#N', 'CC(CC#C)C(C)C#C', 'CC(CC#C)C(C)C#N', 'CC(CC#C)C(N)C#N', 'CC(CC#C)C(O)C#C', 'CC(CC#C)C(O)C#N', 'CC(CC#N)C(C)C#C', 'CC(CC#N)C(C)C#N', 'CC(CC#N)C(N)C#N', 'CC(CC#N)C(O)C#C', 'CC(CC#N)C(O)C#N', 'NC(C#N)C(O)CC#C', 'NC(C#N)C(O)CC#N', 'NC(CC#C)C(N)C#N', 'NC(CC#C)C(O)C#N', 'NC(CC#N)C(N)C#N', 'NC(CC#N)C(O)C#C', 'NC(CC#N)C(O)C#N', 'OC(CC#C)C(O)C#C', 'OC(CC#C)C(O)C#N', 'OC(CC#N)C(O)C#C', 'OC(CC#N)C(O)C#N', 'CC(CC(C)C#C)C#C', 'CC(CC(C)C#N)C#C', 'CC(CC(C)C#N)C#N', 'CC(CC(N)C#N)C#C', 'CC(CC(N)C#N)C#N', 'CC(CC(O)C#C)C#C', 'CC(CC(O)C#C)C#N', 'CC(CC(O)C#N)C#C', 'CC(CC(O)C#N)C#N', 'CC(OC(C)C#C)C#C', 'CC(OC(C)C#N)C#C', 'CC(OC(C)C#N)C#N', 'NC(CC(N)C#N)C#N', 'NC(CC(O)C#C)C#N', 'NC(CC(O)C#N)C#N', 'OC(CC(O)C#C)C#C', 'OC(CC(O)C#N)C#C', 'OC(CC(O)C#N)C#N', 'CC(CC#C)NCC#N', 'CC(CC#N)NCC#C', 'CC(CC#N)NCC#N', 'CN(CCC#C)CC#N', 'CN(CCC#N)CC#C', 'CN(CCC#N)CC#N', 'CC(CNCC#C)C#N', 'CC(CNCC#N)C#C', 'CC(CNCC#N)C#N', 'CC(NCCC#C)C#N', 'CC(NCCC#N)C#C', 'CC(NCCC#N)C#N', 'NC(CNCC#C)C#N', 'NC(CNCC#N)C#N', 'OC(CNCC#C)C#N', 'OC(CNCC#N)C#C', 'OC(CNCC#N)C#N', 'CCC#CC(NC)C#N', 'CNC(C#N)C#CCO', 'CC(C#C)C(CO)C#C', 'CC(C#C)C(CO)C#N', 'CC(C#N)C(CN)C#N', 'CC(C#N)C(CO)C#C', 'CC(C#N)C(CO)C#N', 'CCC(C#C)C(C)C#C', 'CCC(C#C)C(C)C#N', 'CCC(C#C)C(N)C#N', 'CCC(C#C)C(O)C#C', 'CCC(C#C)C(O)C#N', 'CCC(C#N)C(C)C#C', 'CCC(C#N)C(C)C#N', 'CCC(C#N)C(N)C#C', 'CCC(C#N)C(N)C#N', 'CCC(C#N)C(O)C#C', 'CCC(C#N)C(O)C#N', 'COC(C#C)C(C)C#C', 'COC(C#C)C(C)C#N', 'COC(C#C)C(N)C#N', 'COC(C#C)C(O)C#C', 'COC(C#C)C(O)C#N', 'COC(C#N)C(C)C#C', 'COC(C#N)C(C)C#N', 'COC(C#N)C(N)C#C', 'COC(C#N)C(N)C#N', 'COC(C#N)C(O)C#C', 'COC(C#N)C(O)C#N', 'NC(C#C)C(CO)C#N', 'NC(C#N)C(CO)C#C', 'NC(C#N)C(CO)C#N', 'NCC(C#N)C(N)C#N', 'NCC(C#N)C(O)C#N', 'OCC(C#C)C(O)C#C', 'OCC(C#C)C(O)C#N', 'OCC(C#N)C(O)C#C', 'OCC(C#N)C(O)C#N', 'CCC(C)(OC=N)C#C', 'CCC(C)(OC=N)C#N', 'CC(CO)C(C#C)C#C', 'CC(CO)C(C#C)C#N', 'CCC(C)C(C#C)C#C', 'CCC(C)C(C#C)C#N', 'CCC(C)C(C#N)C#N', 'CCC(O)C(C#C)C#C', 'CCC(O)C(C#C)C#N', 'COC(C)C(C#C)C#C', 'COC(C)C(C#C)C#N', 'NC(CO)C(C#C)C#N', 'OCC(O)C(C#C)C#C', 'OCC(O)C(C#C)C#N', 'CNC(CC#N)C#CC', 'CCN(CC#C)CC#N', 'CCN(CC#N)CC#N', 'CNC(CC#N)CC#N', 'CNC(CC#CC)C#N', 'CCC(NCC#C)C#N', 'CCC(NCC#N)C#C', 'CCC(NCC#N)C#N', 'CNC(CCC#C)C#N', 'CNC(CCC#N)C#N', 'NCC(NCC#C)C#N', 'NCC(NCC#N)C#N', 'OCC(NCC#C)C#N', 'OCC(NCC#N)C#C', 'OCC(NCC#N)C#N', 'CCNC(C#N)C#CC', 'CCNC(CC#C)C#N', 'CCNC(CC#N)C#C', 'CCNC(CC#N)C#N', 'CNCC(CC#N)C#N', 'CCCNC(C#C)C#N', 'CCCNC(C#N)C#N', 'OCCNC(C#C)C#N', 'OCCNC(C#N)C#N']

properties = [Aromatic,NumRings,NumAtoms,ContainsTriple,ContainsNitrogen,ContainsCarbon,ContainsFluorine,ContainsOxygen]

molecules2 = ['CC#C', 'CC#N', 'CC#CC', 'CCC#C', 'CCC#N', 'NCC#N', 'OCC#C', 'OCC#N', 'CC#CCO', 'CCC#CC', 'CCCC#C', 'CCCC#N', 'COCC#C', 'COCC#N', 'OCCC#C', 'OCCC#N', 'CCC#CCC', 'CCC#CCO', 'OCC#CCO', 'CC#CCCO', 'CCCC#CC', 'COCC#CC', 'CCCCC#C', 'CCCCC#N', 'CCOCC#C', 'CCOCC#N', 'COCCC#C', 'COCCC#N', 'OCCCC#C', 'OCCCC#N', 'CCC#CCCO', 'CCC#CCOC', 'CCCC#CCC', 'CCCC#CCO', 'COCC#CCO', 'OCCC#CCO', 'CC#CCCCO', 'CCCCC#CC', 'CCOCC#CC', 'COCCC#CC', 'CCCCCC#C', 'CCCCCC#N', 'CCCOCC#C', 'CCCOCC#N', 'CCOCCC#C', 'CCOCCC#N', 'COCCCC#C', 'COCCCC#N', 'OCCCCC#C', 'OCCCCC#N', 'OCCOCC#C', 'OCCOCC#N', 'CCCC#CCCC', 'CCCC#CCCO', 'CCCC#CCOC', 'COCC#CCCO', 'COCC#CCOC', 'OCCC#CCCO', 'CCC#CCCCO', 'CCC#CCCOC', 'CCCCC#CCC', 'CCCCC#CCO', 'CCOCC#CCC', 'CCOCC#CCO', 'COCCC#CCO', 'OCCCC#CCO', 'CC#CCCCCO', 'CC#CCOCCO', 'CCCCCC#CC', 'CCCOCC#CC', 'CCOCCC#CC', 'COCCCC#CC', 'CCCCCCC#C', 'CCCCCCC#N', 'CCCCOCC#C', 'CCCCOCC#N', 'CCCOCCC#C', 'CCCOCCC#N', 'CCOCCCC#C', 'CCOCCCC#N', 'COCCCCC#C', 'COCCCCC#N', 'COCCOCC#C', 'COCCOCC#N', 'OCCCCCC#C', 'OCCCCCC#N', 'OCCCOCC#C', 'OCCCOCC#N', 'OCCOCCC#C', 'OCCOCCC#N', 'CCCC#CCCCO', 'CCCC#CCCOC', 'CCCC#CCOCC', 'CCCCC#CCCC', 'CCCCC#CCCO', 'CCCCC#CCOC', 'CCOCC#CCCO', 'CCOCC#CCOC', 'COCC#CCCCO', 'COCCC#CCCO', 'COCCC#CCOC', 'OCCCC#CCCO', 'CCC#CCCCCO', 'CCC#CCCCOC', 'CCC#CCOCCO', 'CCCCCC#CCC', 'CCCCCC#CCO', 'CCCOCC#CCC', 'CCCOCC#CCO', 'CCOCCC#CCC', 'CCOCCC#CCO', 'COCCCC#CCO', 'OCCCCC#CCO', 'OCCOCC#CCO', 'CC#CCCCCCO', 'CC#CCCOCCO', 'CC#CCOCCCO', 'CCCCCCC#CC', 'CCCCOCC#CC', 'CCCOCCC#CC', 'CCOCCCC#CC', 'COCCCCC#CC', 'COCCOCC#CC', 'CCCCCCCC#C', 'CCCCCCCC#N', 'CCCCCOCC#C', 'CCCCCOCC#N', 'CCCCOCCC#C', 'CCCCOCCC#N', 'CCCOCCCC#C', 'CCCOCCCC#N', 'CCOCCCCC#C', 'CCOCCCCC#N', 'CCOCCOCC#C', 'CCOCCOCC#N', 'COCCCCCC#C', 'COCCCCCC#N', 'COCCCOCC#C', 'COCCCOCC#N', 'COCCOCCC#C', 'COCCOCCC#N', 'OCCCCCCC#C', 'OCCCCCCC#N', 'OCCCCOCC#C', 'OCCCCOCC#N', 'OCCCOCCC#C', 'OCCCOCCC#N', 'OCCOCCCC#C', 'OCCOCCCC#N']

molecules3 = ['O=C1C2C3CN1C23', 'O=C1C2CC3C2N13', 'ON=C1C2C3CC1C23', 'ON=C1C2C3OC1C23', 'O=C1NC23CC(C2)C13', 'O=C1OC23CN(C2)C13', 'CC12CC3C1N2C3=O', 'CC12CN3C1C2C3=O', 'OC12CN3C1C2C3=O', 'CC12C3CC1C(=O)N23', 'CC12C3CN1C(=O)C23', 'CC12CC3C=CC1C23', 'OC12CN3C=NC1C23', 'CC12CC3C1C=CC23', 'CC12C=CC3C1CC23', 'CC12C3CN(C13)C2=O', 'OC12C3CN(C13)C2=N', 'CC12C=CC3CC1C23', 'CC12CC3C1N3C2=O', 'OC12CC3C1N3C2=N', 'OC12CC3C1N3C2=O', 'CC12CC3C(C=C1)C23', 'NC12C3CN=C1OC23', 'OC12C3CN=C1NC23', 'CC12C3CC1C=CC23', 'CC12C3OC1C=CC23', 'O=C1NC2C3CC1C23', 'O=C1NC2C3OC1C23', 'O=C1OC2C3CC1C23', 'O=C1OC2C3NC1C23', 'O=C1OC2C3OC1C23', 'CC1C2C3C=CC1C23', 'O=CC1C2C3CC2C13', 'O=CN1C2C3NC2C13', 'O=CN1C2C3OC2C13', 'O=C1NC2C3CC2C13', 'O=C1NC2C3OC2C13', 'O=C1OC2C3CC2C13', 'O=C1OC2C3NC2C13', 'O=C1OC2C3OC2C13', 'O=C1NC2CC3C2C13', 'O=C1OC2CC3C2C13', 'O=C1N2CC2C11CC1', 'CC1C2C3C1C(=O)N23', 'CC1C2C3C2C(=O)N13', 'OC1C2C3C1C(=O)N23', 'O=CC1C2C3CC1C23', 'CC1C2C3C=CC2C13', 'O=C1NC2C3CC3C12', 'O=C1NC2C3NC3C12', 'O=C1NC2C3OC3C12', 'O=C1OC2C3CC3C12', 'O=C1OC2C3NC3C12', 'O=C1OC2C3OC3C12', 'O=CC1C2C3CC3C12', 'O=CN1C2C3NC3C12', 'O=CN1C2C3OC3C12', 'C1N=C2C3C4C(C13)N24', 'CC12CC3C1C2C3=NO', 'CC12OC3C1C2C3=NO', 'ON=C1C2C3C1CC23O', 'CC12C3CC1C(=NO)C23', 'CC12C3OC1C(=NO)C23', 'ON=C1C2C3CC1C23O', 'ON=C1C2C3NC1C23O', 'CC12C3CC(C13)C2=NO', 'CC12C3OC(C13)C2=NO', 'ON=C1C2CC3C2C13O', 'CC12CC3C(C13)C2=NO', 'CC12OC3C(C13)C2=NO', 'ON=C1C2C3CC1(O)C23', 'CC1C2C3C1C2C3=NO', 'CC1C2C3N1C2C3=NO', 'CN1C2C3C1C2C3=NO', 'ON=C1C2C3C(O)C2C13', 'ON=C1C2CC2C11CC1', 'ON=C1C2NC2C11CC1', 'ON=C1C2OC2C11CC1', 'CC1C2C3C(N13)C2=NO', 'CC1C2C3C2C(=NO)C13', 'CC1C2C3N2C1C3=NO', 'CN1C2C3C2C(=NO)C13', 'ON=C1C2C3C2C1C3O', 'C#CC12C3CN=C1OC23', 'N#CC12C3CN=C1OC23', 'C#CC12OC3=NCC1C23', 'N#CC12OC3=NCC1C23', 'C#CC12CN=C3OC1C23', 'N#CC12CN=C3OC1C23', 'N#CC1C2C3C=CC1C23', 'N#CC1C2C3C=CC2C13', 'C1N=C2C3NC3C3C1N23', 'C1N=C2C3OC3C3C1N23', 'C1C2C1C1C3C=CC1C23', 'O1C2C1C1C3C=CC1C23', 'C1C2C1C1C=CC3C2C13', 'N1C2C1C1N=CN3C2C13', 'O1C2C1C1C=CC3C2C13', 'O1C2C1C1N=CN3C2C13', 'N1C2C3C1C1N=C3OC21', 'O1C2C3C1C1N=C3OC21', 'C1C2C3C4C=CC(C13)C24', 'N1C2C3C1C1N=CN3C21', 'O1C2C3C4C=CC(C13)C24', 'O1C2C3C1C1N=CN3C21', 'C1N=CN2C3C4NC3C124', 'C1C2C3C2C2C=CC3C12', 'O1C2C3C2C2C=CC3C12', 'C1C2C3C4C3C1C=CC24', 'O1C2C3C4C3C1C=CC24', 'C1C2C3C4C=CC2C1C34', 'N1C2C3C4C1C2N=CN34', 'O1C2C3C4C1C2N=CN34', 'O1C2C3C4C=CC2C1C34', 'C1C2C3N=C4OCC13C24', 'C1OC2=NC3C4NC13C24', 'C1OC2=NC3C4OC13C24', 'C1C2C3C=CC4C1C2C34', 'O1C2C3C=CC4C1C2C34', 'C1C2C3C=CC4C1C4C23', 'O1C2C3C=CC4C1C4C23', 'C1C2C3CC4=CCC13C24', 'C1C2C3OC4=NCC13C24', 'C1N=C2NC3C4OC13C24', 'C1N=C2OC3C4NC13C24', 'C1N=C2OC3C4OC13C24', 'C1N=C2C3OC4C1N2C34', 'C1C2N=C3NC4C2N1C34', 'C1C2N=C3OC4C2N1C34', 'C1C2OC3=NC1C1C2N31', 'O1C=NC2C3C4C1C2N34', 'C1OC11C2C3CN=C1N23', 'C1CC11N=C2CC3C1N23', 'C1CC11N=C2NC3C2C13', 'C1CC11N=C2OC3C2C13', 'C1CC2C3C4C2C(=N1)N34', 'N=C1C2CC3(C#N)C2N13', 'O=C1C2C3C(C2C#C)N13', 'O=C1C2C3C(C2C#N)N13', 'O=C1C2C3C2N1C3C#C', 'O=C1C2C3C2N1C3C#N', 'N=C1N2C3CC1(C#N)C23', 'N=C1C2CC3N1C23C#C', 'N=C1C2CC3N1C23C#N', 'N=C1N2CC3C2C13C#N', 'O=C1NC2(CC2)C11CO1', 'O=C1OC2(CC2)C11CN1', 'O=C1OC2(CC2)C11CO1', 'O=C1NC23CC12COC3', 'O=C1NC23CC12OCC3', 'O=C1OC23CC12COC3', 'O=C1OC23CC12NCC3', 'O=C1OC23CC12OCC3', 'O=C1NC23CCC12OC3', 'O=C1OC23CCC12NC3', 'O=C1OC23CCC12OC3', 'O=C1NC2CC12C1CN1', 'O=C1OC2CC12N1CC1', 'O=C1NCC23CC12CO3', 'O=C1NCC23CC12OC3', 'O=C1NCC23CCC12O3', 'O=C1OCC23CC12CO3', 'O=C1OCC23CC12NC3', 'O=C1OCC23CC12OC3', 'O=C1OCC23CCC12N3', 'O=C1OCC23CCC12O3', 'O=CN1C2CC1C21CO1', 'O=CN1C2CC3(CO3)C12', 'O=CN1CC11C2CC1O2', 'O=CN1CC11CC2OC12', 'O=CN1CC23CC(O2)C13', 'O=COC12CC1N1CC21', 'O=COC12CC3C1N3C2', 'O=COC12CN3C1CC23', 'O=COC12CN3CC1C23', 'NC(=O)C12C3CC1OC23', 'NC(=O)C12CC1C1NC21', 'NC(=O)C12CC1C1OC21', 'NC(=O)C12NC1C1CC21', 'NC(=O)C12NC1C1NC21', 'NC(=O)C12OC1C1CC21', 'NC(=O)C12CC3C(O1)C23', 'NC(=O)C12CC3C1OC23', 'NC(=O)C12NC3C1NC23', 'NC(=O)C12OC3C1CC23', 'NC(=O)C12CC3OC1C23', 'NC(=O)C12OC3CC1C23', 'CC(=O)C1C2C3CC1C23', 'NC(=O)C1C2C3C1CN23', 'NC(=O)C1C2C3CC1N23', 'NC(=O)C1C2C3CN1C23', 'NC(=O)C1C2CC3C2N13', '[O-]C(=O)C1C2C3[NH2+]C1C23', 'CC(=O)C1C2C3CC2C13', 'CC(=O)N1C2C3NC2C13', 'NC(=O)C1C2C3CC2N13', 'NC(=O)C1C2C3CN2C13', '[O-]C(=O)C1C2C3[NH2+]C2C13', 'CC(=O)C1C2C3CC3C12', 'CC(=O)N1C2C3NC3C12', 'NC(=O)C1C2C3CC3N12', 'NC(=O)C1C2C3CN3C12', 'NC(=O)C1C2C3NC3C12', 'CC1(C)C2C(=O)CC12C', 'CC1(C)C2C3C1C(=O)N23', 'CC1(C)C2C3C2C(=O)N13', 'CC1(O)C2C3C1C(=O)N23', 'CC1(C)C2C3C=CC1C23', 'CC1(C)C2C3C=CC1N23', 'CC1(O)C2C3C1N=CN23', 'CC1(O)C2C3C=CC1C23', 'CC1(C)C2C3C=CC2C13', 'CC1(C)C2C3C=CC2N13', 'CC1(O)C2C3C=CC2C13', 'CC1(C)C2CC(=O)C12C', 'CC1(O)C2NC(=O)C12N', 'CC12C3C4C1CN=C4N23', 'OC12C3C1C1=NCC2N31', 'CC12C3C4CN=C(C14)N23', 'OC12C3C1N1C3CN=C21', 'CC12C3C4C1CN=C2N34', 'OC12C3C4C1CN=C2N34', 'NC12C3C4N=C1OC2C34', 'OC12C3N=C4CC1N4C23', 'CC12C3C(CN13)NC2=O', 'CC12C3C(CN13)OC2=O', 'NC12C3CC(OC1=O)C23', 'OC12C3CC(NC1=O)C23', 'CC12C3CC1(C)C(=O)N23', 'CC12C3CC1(O)C(=N)N23', 'CC12C3CC1(O)C(=O)N23', 'CC12C3CC1(C)C=CC23', 'CC12C3CC1(C)N=CN23', 'CC12C3CC1(O)C=CC23', 'CC12C3OC1(C)C=CC23', 'CC12CC3C(C=C1)C23O', 'CC12C3CC1(C3)NC2=O', 'CC12N3CC1(C3)OC2=N', 'CC12N3CC1(C3)OC2=O', 'NC12C3CC1(C3)OC2=O', 'OC12C3CC1(C3)NC2=N', 'OC12C3CC1(C3)NC2=O', 'OC12C3CC1(C3)OC2=O', 'CC12C3CN1C(=N)C23O', 'CC12C3CN1C(=O)C23C', 'CC12C3CN1C(=O)C23O', 'CC12C3CN1C(=O)C2N3', 'CC12C3CC1C(=O)CC23', 'CC12C3NC(=O)C1OC23', 'CC12C3NC1C(=O)OC23', 'NC12C3NC(=O)C1OC23', 'NC12C3NC1C(=O)NC23', 'NC12C3NC1C(=O)OC23', 'NC12C3OC1C(=O)OC23', 'OC12C3CC1C(=O)NC23', 'OC12C3NC(=N)C1OC23', 'OC12C3NC(=O)C1OC23', 'OC12C3NC1C(=N)OC23', 'OC12C3NC1C(=O)NC23', 'OC12C3NC1C(=O)OC23', 'OC12C3OC1C(=O)OC23', 'OC12C3CC1N(C=O)C23', 'CC12C3NC1C3N2C=O', 'CC12C3OC1C3N2C=O', 'CC12C3NC1C3NC2=O', 'CC12C3NC1C3OC2=O', 'CC12C3OC1C3NC2=O', 'NC12C3CC1C3OC2=O', 'NC12C3NC1C3NC2=O', 'NC12C3NC1C3OC2=O', 'NC12C3OC1C3NC2=O', 'NC12C3OC1C3OC2=O', 'OC12C3CC1C3NC2=O', 'OC12C3NC1C3NC2=O', 'OC12C3NC1C3OC2=N', 'OC12C3NC1C3OC2=O', 'OC12C3OC1C3NC2=N', 'OC12C3OC1C3NC2=O', 'OC12C3OC1C3OC2=O', 'CC12C3CC1CC(=O)C23', 'CC12C3N1CC2NC3=O', 'CC12C3N1CC2OC3=O', 'OC12C3CC1NC(=O)C23', 'CC12C3NC3C1N2C=O', 'CC12C3OC3C1N2C=O', 'CC12C(NC1=O)C1CN21', 'CC12C(OC1=O)C1CN21', 'CC12C3NC3C1NC2=O', 'CC12C3NC3C1OC2=O', 'CC12C3OC3C1NC2=O', 'NC12C3CC3C1OC2=O', 'NC12C3NC3C1NC2=O', 'NC12C3NC3C1OC2=O', 'NC12C3OC3C1NC2=O', 'NC12C3OC3C1OC2=O', 'OC12C3CC3C1NC2=O', 'OC12C3NC3C1NC2=N', 'OC12C3NC3C1NC2=O', 'OC12C3NC3C1OC2=N', 'OC12C3NC3C1OC2=O', 'OC12C3OC3C1NC2=N', 'OC12C3OC3C1NC2=O', 'OC12C3OC3C1OC2=O', 'CC12C3CCN(C13)C2=O', 'OC12C3CCN(C13)C2=O', 'CC12C3CCN1C(=O)C23', 'CC12C=CC3(C)C1CC23', 'CC12C=CC3(C)C1OC23', 'CC12C=CC3(C)CC1C23', 'CC12C=CC3(C)OC1C23', 'CC12C3CC1C=CC23O', 'CC12C=CC3CC1C23C', 'CC12C=CC3CC1C23O', 'CC12C=CC3OC1C23C', 'CC12CC(=O)C1C1CC21', 'CC12NC(=O)C1C1NC21', 'CC12NC(=O)C1C1OC21', 'CC12NC(=O)C1N1CC21', 'CC12OC(=O)C1C1NC21', 'CC12CC(=O)C3C1CC23', 'CC12NC(=O)C3C1NC23', 'CC12NC(=O)C3C1OC23', 'CC12OC(=O)C3C1NC23', 'CC12CC(=O)C3CC1C23', 'CC12NC(=O)C3OC1C23', 'CC12OC(=O)C3NC1C23', 'CC12CN1C(=O)C1NC21', 'CC12NC1C(=O)N1CC21', 'CC12CN1C(=N)C21CO1', 'CC12CN1C(=O)C21CC1', 'CC12CN1C(=O)C21CN1', 'CC12CN1C(=O)C21CO1', 'CC12NC1C1C2N1C=O', 'CC12OC1C1C2N1C=O', 'OC12CC1C1C2N1C=O', 'CC12CC1C1C2CC1=O', 'CC12CN1C1C2NC1=O', 'CC12NC1C1C2NC1=O', 'CC12NC1C1C2OC1=O', 'CC12OC1C1C2NC1=O', 'OC12CC1C1C2NC1=O', 'CC12CC1C1CC(=O)C21', 'CC12NC1C1NC(=O)C21', 'CC12NC1C1OC(=O)C21', 'CC12OC1C1NC(=O)C21', 'OC12CC1C1NC(=O)C21', 'CC12NC1C1CN1C2=O', 'CC12CC3(C)C1N2C3=O', 'CC12CC3(O)C1N2C3=O', 'CC12CC3(C)C1C=CC23', 'CC12CC3(O)C1C=CC23', 'CC12OC3(C)C1C=CC23', 'CC12CC3(C)C=CC1C23', 'CC12CC3(O)C(C=C1)C23', 'CC12CC3(O)C=CC1C23', 'CC12OC3(C)C=CC1C23', 'CC12CC3(CC(=O)C13)C2', 'CC12CC3=CCC1(C)C23', 'CC12CN=C3NC1(C)C23', 'CC12NC3=NCC1(O)C23', 'CC12OC3=NCC1(O)C23', 'CC12CC3=CCC1C23C', 'CC12NC3=NCC1C23O', 'CC12OC3=NCC1C23N', 'CC12OC3=NCC1C23O', 'CC12C3C(=O)N1CC23O', 'CC12CC3C(=O)N1C23C', 'CC12CN3C(=O)C1C23C', 'CC12CC3C(C13)C(=O)C2', 'CC12CN3C(C13)C(=O)N2', 'CC12NC3C(OC1=O)C23', 'CC12OC3C(NC1=O)C23', 'OC12CC3C(NC1=O)C23', 'CC12NC3C1N(C=O)C23', 'CC12OC3C1N(C=O)C23', 'OC12CC3C1N(C=O)C23', 'CC12C3N(CC13O)C2=O', 'CC12CN3C1C2(C)C3=O', 'CC12CN3C1C2(O)C3=N', 'CC12CN3C1C2(O)C3=O', 'OC12CN3C1C2(O)C3=N', 'OC12CC3C1C2N3C=O', 'CC12CC3C1C2CC3=O', 'CC12NC3C1C2OC3=O', 'CC12OC3C1C2NC3=O', 'OC12CC3C1C2NC3=O', 'CC12C=CC3C1CC23O', 'CC12CC3C1C=CC23C', 'CC12CC3C1C=CC23O', 'CC12OC3C1C=CC23C', 'CC12CC3C1CC(=O)C23', 'CC12NC3C1NC(=O)C23', 'CC12NC3C1OC(=O)C23', 'CC12OC3C1NC(=O)C23', 'OC12CC3C1NC(=O)C23', 'CC12C=CC3CC1(O)C23', 'CC12CC3C=CC1(C)C23', 'CC12CC3C=CC1(O)C23', 'CC12N=CN3CC1(O)C23', 'CC12OC3C=CC1(C)C23', 'CC12C3C=CC1CC23O', 'CC12C3N=CN1CC23O', 'CC12CC3C=CC1C23C', 'CC12CC3C=CC1C23O', 'CC12CC3N=CN1C23C', 'CC12OC3C=CC1C23C', 'CC12CC3CC(=O)C1C23', 'OC12CC3NC(=O)C1C23', 'CC12NC3CN(C13)C2=O', 'CC12CN=C3C1C1C2N31', 'OC12CN=C3C1C1C2N31', 'CC12C3NC1=NCC23O', 'CC12C3OC1=NCC23O', 'CC12CC=C3CC1C23C', 'CC12CN=C3NC1C23O', 'CC12CN=C3OC1C23N', 'CC12CN=C3OC1C23O', 'NC12C3OC1=NCC23O', 'NC12CN=C3OC1C23O', 'OC12CN=C3OC1C23O', 'CC12N=C(N)C3OC1C23', 'CC1=CC2C3C1CC23C', 'CC1=CC2C3CC1C23C', 'NC1=NC2C3OC1C23O', 'CC1=CC2C3CC2(C)C13', 'NC1=NC2C3NC2C13O', 'NC1=NC2C3OC2C13O', 'CC1=CC2CC3(C)C2C13', 'CC1=CC2CC3C1C23C', 'CN1C(=O)C2C3CC12C3', 'CC12NC(=O)C1(O)C2O', 'CC12OC(=O)C1(N)C2N', 'CC1C2(C)CC(=O)C12C', 'OC1C23CCC12N=CN3', 'OC1C23CCC12N=CO3', 'CC1C2N1C(=O)C21CN1', 'CC1C2N1C(=O)C21CO1', 'CC1C2C1C1C2CC1=O', 'CC1C2C3NC(=O)C3N12', 'CC1C2C3OC(=O)C3N12', 'CN1C2C1C1C2NC1=O', 'CN1C2C1C1C2OC1=O', 'OC1C2C1C1C2NC1=O', 'OC1C2C1C1C2OC1=O', 'CC12C3C(O)C1C(=O)N23', 'CC1C2C3C(=O)N1C23C', 'CC1C2N3C(=O)C1C23C', 'CC12C3C(C1O)C(=O)N23', 'CC1C2C3N(C2=O)C13C', 'CC1N2C3C(C2=O)C13C', 'CC1N2C3C(C2=O)C13O', 'CC1C2C3C(C=O)C2C13', 'CN1C2C3C1C2N3C=O', 'OC1C2C3C(C=O)C2C13', 'CN1C2C3C1=NCC23O', 'CC1C2C3N1C(=O)C23C', 'CC1C2C3N1C(=O)C23O', 'CC12C=CC3C1C2C3O', 'CC1C2C3C1C=CC23C', 'CC1C2C3C1C=CC23O', 'CC12C3C(C1O)N3C2=O', 'CC1C2C3N2C(=O)C13C', 'CC1C2C3N2C(=O)C13O', 'NC12C3C(C1O)N3C2=O', 'NC1C2C3N2C(=N)C13O', 'OC1C2C3N2C(=N)C13O', 'OC1C2C3N2C(=O)C13O', 'CC1C2C3C2C1C3C=O', 'OC1C2C3C2C1C3C=O', 'OC1C2C3C4N=C1N2C34', 'CC1=CC2C3C2C1C3O', 'CC1C2C3C=C(C)C1C23', 'CC1C2C3C=C(C)C2C13', 'CC12C=CC3C(C13)C2O', 'CC12N=CN3C(C13)C2O', 'CC1C2C3C=CC1(C)C23', 'CC1C2C3C=CC1(O)C23', 'CC1C2C3N2C=NC13C', 'CC12C3C(N)C1N=CN23', 'CC12C3C(O)C1N=CN23', 'CC12C3C=CC1C(O)C23', 'CC1C2C3C=CC1C23C', 'CC1C2C3C=CC1C23O', 'CC1C2N3C=NC1C23C', 'CN1C2C3C=CC1C23C', 'CC12C=CC3C1C(O)C23', 'CC1C2C3C=CC2(C)C13', 'CC1C2C3C=CC2(O)C13', 'CC12C(O)C3C1C=CC23', 'CC1C2C3C=CC2C13C', 'CC1C2C3C=CC2C13O', 'CC1N2C3C=CC2C13C', 'CC1C2C3CC(=O)C1C23', 'CN1C2C3NC(=O)C1C23', 'CN1C2C3OC(=O)C1C23', 'OC1C2C3NC(=O)C1C23', 'OC1C2C3OC(=O)C1C23', 'CC1C2C3CC(=O)C2C13', 'CN1C2C3NC(=O)C2C13', 'CN1C2C3OC(=O)C2C13', 'OC1C2C3NC(=O)C2C13', 'OC1C2C3OC(=O)C2C13', 'CN1C2C3OC(C23)C1=O', 'CN1C2C3CC2(C3)C1=O', 'CN1C2C3CC2(N3)C1=O', 'CN1C2C3CC2(O3)C1=O', 'CN1C2C3CC2C3C1=O', 'CN1C2C3NC2C3C1=O', 'CN1C2C3OC2C3C1=O', 'CC1C2C3NC3C(=O)N12', 'CN1C2C3NC3C2C1=O', 'CN1C2C3OC3C2C1=O', 'CN1C2C3CN=C1C23O', 'CC12C3C(N=CN13)C2O', 'CC12C3C=CC(C13)C2O', 'CC1C2C=CC3C2C13C', 'CC1C2C=CC3C2C13O', 'CC1C2C=CC3N2C13C', 'CC1N2C=NC3C2C13N', 'CC1N2C=NC3C2C13O', 'CC1C2CC11CC(=O)C21', 'CC1N2CC11NC(=O)C21', 'OC1C2CC11NC(=O)C21', 'CC1C2NC3C2N1C3=O', 'CC1N=C2C3C4C(C13)N24', 'CC1N=C2NC3(C)C1C23', 'CC1N=C2NC3C1C23O', 'CC1N=C2NC3C2C13C', 'CC1N=C2NC3C2C13O', 'CC1CC11C2CN2C1=O', 'CC1NC11C2CN2C1=O', 'CC1OC11C2CN2C1=O', 'OC1CC11C2CN2C1=O', 'OC1CC23CC12N=CO3', 'OC1CC23CC12OC=N3', 'CN=C1N2C3CC1(O)C23', 'CN=C1N2CC3C2C13O', 'CCC12C3CN(C13)C2=O', 'COC12C3CN(C13)C2=N', 'COC12C3CN(C13)C2=O', 'OCC12C3CN(C13)C2=O', 'COC12C3CN=C1NC23', 'OCC12C3CN=C1NC23', 'CCC12C3CC1C(=O)N23', 'CCC12C3CN1C(=O)C23', 'OCC12C3CC1C(=O)N23', 'OCC12C3CN1C(=O)C23', 'CCC12C3CC1C=CC23', 'CCC12C3OC1C=CC23', 'OCC12C3CC1C=CC23', 'OCC12C3CC1N=CN23', 'CCC12C=CC3C1CC23', 'CCC12C=CC3CC1C23', 'CCC12C=CC3OC1C23', 'OCC12NC3=NCC1C23', 'CCC12CC3C1N3C2=O', 'COC12CC3C1N3C2=O', 'OCC12CC3C1N3C2=O', 'CCC12CC3C(C=C1)C23', 'CCC12CC3C1N2C3=O', 'CCC12CN3C1C2C3=O', 'COC12CN3C1C2C3=O', 'OCC12CC3C1N2C3=O', 'OCC12CN3C1C2C3=O', 'CCC12CC3C1C=CC23', 'CCC12OC3C1C=CC23', 'OCC12CC3C1C=CC23', 'CCC12CC3C=CC1C23', 'CCC12OC3C=CC1C23', 'OCC12CC3C=CC1C23', 'COC12CN=C3NC1C23', 'CCC1C2C3C1C(=O)N23', 'CCC1C2C3C2C(=O)N13', 'COC1C2C3C1C(=O)N23', 'OCC1C2C3C1C(=O)N23', 'OCC1C2C3C2C(=O)N13', 'OCC1C2C3C=CC1C23', 'OCC1C2C3C=CC2C13', 'N1C2C3OC4C3C1=NC24', 'CC12C3C4CC1N=C4N23', 'CC12C3C4CC1C(=N4)N23', 'CC12C3C1N1CC3N=C21']

molecules4 =[1,2,3]

class Similarity_Analysis():
    def __init__(self,molecules,properties,significance=0.1):
        self.molecules = molecules
        self.properties = [prop(molecules) for prop in properties]
        self.entropy = [x.entropy() for x in self.properties]
        self.significance = significance

    def __str__(self):
        header  = ["Subset Description"]
        results = ["\t{:<60}".format(prop.summative_label(significance=self.significance)) \
                  for prop in self.properties if prop.summative_label(significance=self.significance) is not None]
        return "\n".join(header + results)


if __name__ == "__main__":
    analy = Similarity_Analysis(molecules,properties)
    print(analy.entropy)

    print(analy)

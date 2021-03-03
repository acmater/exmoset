import unittest
import numpy as np


from molecule import *
from atom import *
from bond import *
from substructure import *
from data import *

test_mols = molecules4

class TestSubstructure(unittest.TestCase):
    def test_substructures(self):
        substructures = Substructure(test_mols,substructures=["[NH2]"])
        assert sum(substructures.values["[NH2]"]) == 1, "Incorrectly identifying substructures"

if __name__ == "__main__":
    unittest.main()

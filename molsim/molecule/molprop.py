import numpy as np
from rdkit import Chem
from property import Property

class MolProp(Property):
    """
    Given a list of property labels and a dataframe, will extract meaningful labels.
    """
    def __init__(self,molecules,properties,df):
        super().__init__(molecules)
        assert isinstance(molecules[0],str), "Molecules for this method must be provided as immutable, string representations"

        self.ent_type    = "Continuous"
        self.properties  = properties
        self.values      = self.calc_property(molecules,df)
        print(self.values)

    def calc_property(self,molecules,df):
        props = {}
        for prop in self.properties:
            props[prop] = df[prop].to_numpy()
        return props

    def summative_label(self,significance=0.1,verbose=True):
        summary = []
        for prop in self.properties:
            if verbose:
                print(f"Working on property: {prop}")
                if self.entropy(self.values[prop],self.ent_type) < significance:
                    print(f"Inside the inner loop. Entropy is {self.entropy(self.values[prop],self.ent_type)}")
                    print(f"Due to signifiance, calculating average value of {prop}")
                    print(f"Average value of {prop} is {np.mean(self.values[prop])}")
                    print()

                    summary.append(f"These molecules share similar values of {prop} centred around {np.mean(self.values[prop])}")
            else:
                if self.entropy(self.values[prop],self.ent_type) < significance:
                    summary.append(f"These molecules share similar values of {prop} centred around {np.mean(self.values[prop])}")
        if summary:
            return "\n".join(summary)
import numpy as np
from molsim.property import Property
from molsim.labels import Continuous

class MolProp(Property):
    """
    Given a list of property labels and a dataframe, will extract meaningful labels.
    """
    def __init__(self,molecules,properties,df):
        super().__init__(molecules)
        assert isinstance(molecules[0],str), "Molecules for this method must be provided as immutable, string representations"

        self.ent_type    = "Continuous"
        self.values      = self.calc_property(molecules,properties,df)

    def calc_property(self,molecules,properties,df):
        props = {}
        for prop in properties:
            props[prop] = Continuous(prop,
                                     df[prop].to_numpy(),
                                     context="molecule")
        return props

    def summative_label(self,significance=0.1,verbose=True):
        summary = []
        for prop, label in self.values.items():
            if self.entropy(label) < significance:
                print(f"Working on property: {prop}")
                if verbose:
                    print(f"Inside the inner loop. Entropy is {self.entropy(label)}")
                    print(f"Due to signifiance, calculating average value of {prop}")
                    print(f"Average value of {prop} is {(self[prop].av)}")
                    print()

                summary.append(label.summary())

        if summary:
            return "\n".join(summary)

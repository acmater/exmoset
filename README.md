# EXplainable MOlecular SETs

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties.

This method is guided by the underlying principle of finding labels with minimum entropy
across the provided set of molecules. Intuitively these labels represent a "pure" description
of the space of interest, with the label adopting almost exclusively a single value.

## API
The MolSet Class handles the analysis of a given molecular set in accordance with the list of fingerprints provided. The molecules can be passed to MolSet in any format, with additional conversions specified by the mol_converter argument.

```python
analysis = MolSet(molecules,
                  fingerprints = fingerprints,
                  mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                  significance=0.1,
                  file="data/QM9_Data.csv")
```

### Fingerprints
Fingerprints are a standardized way for Molset to calculate the properties for each molecule it is analysing. Its arguments determine how the property will be calculated, custom sensititivies, and through a variety of arguments how exactly the descriptor sentence will be constructed.

```python
def contains_C(mol):
      return 1 if C in mol else 0

contains_carbon = Fingerprint(name="Contains C",
                  context="molecule",
                  label_type="binary",
                  calculator=contains_C,
                  mol_format="smiles"
                  file=None)
```

### Molecular Converters
The mol_converters argument provides the means to transform each molecule into alternate representations. The argument is a dictionary with the following structure {Identifier : Function_that_will_convert} that is expanded in the following way:

```python
formats = {key : mol_converters[key](mol) for key in mol_converters.keys()} # Assigns each identifier to its assocaited representation by
self.Molecules.append(Molecule(mol, **formats)) # Unpacks the new formats as kwargs into the Molecule object
```

## Label Types
### Binary
Binary labels indicate the presence of absence of a particular element, bond type, or molecular feature (such as aromaticity). Simplest to calculate and best behaved with respect to the entropy estimators. Uses a discrete entropy estimator.

### Multiclass
Discrete labels where the value can be any integer. Examples include number of rings, number of atoms, or number of each type of bond. Uses a discrete entropy estimator.

### Continuous
Continuous labels where the value can be any real number. Examples include electronic spatial extent, dipole moment, and free energy. Uses the continuous entropy estimator

The properties were built with extensability in mind. The API for a property is
defined by the Property abstract class with its associated methods.

The majority of the work is handled by the property metaclass, which contains the key entropy functionality
for both continuous and discrete cases. This is called by all of its subclasses to determine the importance
of each of the properties that they are assessing.

Property instances are further decomposed into three separate subclasses - atom, bond, and molecule
properties. All of these utilize the API provided by property, and the information provided
by similarity analysis is sorted in accordance with that.

## Entropy Estimation
Key to the functioning of this approach is the estimation of entropy for a given random variable. Both discrete (for properties such as aromaticity) and continuous (for properties like electronic spatial extent) have to be considered, and different approaches are used to estimate the entropy.

The discrete entropy is estimated using the plug in estimator that assumes a uniform probability over the classes represented in the vector. A simple implenetation was adapted from -

The continuous entropy estimator uses the Kozachenko and Leonenko (1987) estimator in which k-nearest neighbour distances are used to approximate the entropy of the underlying distribution. The implementation is provided by the entropy esitmators package of Paul Broderson (https://github.com/paulbrodersen/entropy_estimators).

### TODO
So the final major code refactorization that I want to consider now is a further decomposition of the problem.
Currently molprop, substructure, atom, and bond are the four main classes, but their functionality is restricted. Having these as their own superclasses that inherit from property and then spawning particular instances like atom_types would allow me to extend it further to other atomic properties such as charge, radical, hybridisation, etc.

0. Update docstrings throughout codebase.
1. Write new tests for all of the new functionalities
2. Another nice feature would be the capacity to select two distributions and identify what properties are different between them.
   1. There's actually more scope here than looking at differences. The labels can be viewed as logical descriptions of the group, with logical operations defined on them.
   2. To do this I need to break down the summative labels into properties that get assigned to the similarity analysis class (with an associated print statement) and then when operations like and or addition are applied, it performs labelwise comparison.
   3. This would be nice as it would readily enable comparison between groups.
   4. I think this is quite a bit simpler than above, just make a set of the labels that two particular groups share.
   5. Of importance is the possibility that I may need to consider mutual information labelling. I need to understand this better.
3. Add plotting functionality.
4. Add feature to represent label of cluster as a vector, can then strongest outliers for a given set.
5. Add feature to identify outliers from a group
   1. This one involves

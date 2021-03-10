# EXplainable MOlecular SETs

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties.

This method is guided by the underlying principle of finding labels with minimum entropy
across the provided set of molecules. Intuitively these labels represent a "pure" description
of the space of interest, with the label adopting almost exclusively a single value.

## API
The MolSet Class handles the analysis of the molecular set with a provided set of properties. These properties are provided as fingerprints which provide information about how to calculate the property, how it  will be described (defined by context and name).

```python
a = Fingerprint(name="Contains C",
                  context="molecule",
                  label_type="binary",
                  calculator="add")
```

## Context
There are four different contexts that will be considered by MolSet (atom, bond, molecule, and substructure). Each context dictates what the label describes and determines internally when it will be calculated.

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

1. Another nice feature would be the capacity to select two distributions and identify what properties are different between them.
   1. There's actually more scope here than looking at differences. The labels can be viewed as logical descriptions of the group, with logical operations defined on them.
   2. To do this I need to break down the summative labels into properties that get assigned to the similarity analysis class (with an associated print statement) and then when operations like and or addition are applied, it performs labelwise comparison.
   3. This would be nice as it would readily enable comparison between groups.
   4. I think this is quite a bit simpler than above, just make a set of the labels that two particular groups share.
   5. Of importance is the possibility that I may need to consider mutual information labelling. I need to understand this better.
2. Add plotting functionality.
3. Speed it up. The code was written to get a working prototype as quickly as possible, but the current codebase is slow to execute and their should be some ways to dramatically accelerate it.
4. Add code snippet examples for how to run
5. Add feature to represent label of cluster as a vector, can then strongest outliers for a given set.
6. Add feature to identify outliers from a group
   1. This one involves


#### Notes
There are some portions of the code that look like they could be offloaded to the abstract class. In particular
```python
def summative_label
```
appears to duplicated. The reason for its repeated use is that each one is slightly different in the print statements that it provides. There may be a way to automate this by providing template strings, but currently that has not been examined.

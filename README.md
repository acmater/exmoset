# EXplainable MOlecular SETs

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties.

This method is guided by the underlying principle of finding labels with minimum entropy
across the provided set of molecules. Intuitively these labels represent a "pure" description
of the space of interest, with the label adopting almost exclusively a single value.

## API
The `MolSet` Class handles the analysis of a given molecular set in accordance with the list of fingerprints provided. The molecules can be passed to MolSet in any format, with additional conversions specified by the mol_converter argument.

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

## Entropy Estimation
Key to the functioning of this approach is the estimation of entropy for a given random variable. Both discrete (for properties such as aromaticity) and continuous (for properties like electronic spatial extent) have to be considered, and different approaches are used to estimate the entropy.

The discrete entropy is estimated using the plug in estimator that assumes a uniform probability over the classes represented in the vector. A simple implenetation was adapted from -

The continuous entropy estimator uses the Kozachenko and Leonenko (1987) estimator in which k-nearest neighbour distances are used to approximate the entropy of the underlying distribution. The implementation is provided by the entropy esitmators package of Paul Broderson (https://github.com/paulbrodersen/entropy_estimators).

## Vector Description
Internal operations utilize an information vector to describe the MolSet. This information vector is produced by the `calc_vector` method and returns a masked array showing the average meaningful (meaning it is rounded to an integer for binary and multiclass labels) label with values that exceed the entropy theshold masked. The comparison operations between sets are expedited through the use of these vector descriptions.

## Grammatical Structure
All nouns are pluran

## TODO
The big task is to get mutual information working.

I also need to decide on a clustering problem to sink this codebase's teeth into.

Critically, MolSet has to be able to function in an entirely self-contained manner and as a sub-component of a molspace. It makes the most sense to precompute all label values in the molspace and then pass the information explicitly to the MolSets (either by passing them indexes or by copying the portions of the dataframe that are relevant).

The following are the steps:
0. I need to work out exactly how to phrase the problem of mutual information and deal with what seems to be a disagreement between the Stanford NLP guide and the wikipage on mutual information regarding the meaning of the two terms.
   1. I think that the resolution of the disagreement is in my own misunderstanding of independence. I assumed independence meant (in addition to its typical probability meaning) that the probability distributions were entirely separate in their outcome space, however there is no reason that the must be. Using the KLDiv definition and the fact that two distributions have a KLDiv of 0 only if they are identical mean that more positive values for mutual information indicate a greater degree of separation between the two distributions.
1. Add MolSpace class that handles and entire dataset and offloads work to MolSet class.
   1. This means that MolSet will need its index label so that it can specify which cluster it is referring to.
   2. I might want to change MolSet back to Subspace in terms of class name, although Molset does make some sense.
   3. I need to ensure that MolSet functions in an entirely self-contained manner.
      1. Actually, maybe it shouldn't. You could achieve equivalent behaviour by calling MolSpace on whatever single set you're interested in and then just specify either no clusters, or a cluster where it is all one big cluster.
2. Estimating the mutual information for discrete labels is straightforward, basically adapts the Stanford NLP approach.
   1. The first goal should be to implement code which works for binary labels using the Stanford NLP approach (https://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html).
3. Estimating the mutual information in the continuous case in a open problem, but there is what looks like a good solution in the work of Krakov (https://arxiv.org/pdf/cond-mat/0305641.pdf).
4. There is a github implementation of the Kraskov estimator (https://github.com/mutualinfo/mutual_info/tree/main).
5. Will also need a way to visualise the contingency matrix for discrete cases.
6. I need to figure out how to deal with the calculation of mutual information when a particular class has 0 elements, using this as an example (https://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html).
7. Need to resolve how to deal with the negative entropy values that are possible when using k-means clustering approaches.

# Molsim

Package to automate the identification of molecular similarity given an arbitrary set
of molecular properties.

### TODO

Figure out how to properly set abstract method in properties class to require that all
subclasses have the .values property without creating and AttributeError.

I would like to extend this code to continuous descriptors such as electronic spatial extent or something similar.
I think however that doing that would require implementing https://en.wikipedia.org/wiki/Limiting_density_of_discrete_points.
Which would be very tricky.

Additionally, these continuous decisions need to be taken against some standard I believe. What exactly
this standard should be is difficult.

The reason for this is that as they are not discrete, they will never have multiple instances of the same class, so the notion
of "purity" that is intuitive for discrete distributions is lost in the continuous case.

There are some alternatives https://thomasberrett.github.io/EntropyEstimation.pdf and https://github.com/paulbrodersen/entropy_estimators.

But the best formulation may be using KL-Div between the original distribution (in the entire dataset), and the distribution in the subset of interest.

1. Add the ability to print out molecules that could not be converted in Chem.mol objects

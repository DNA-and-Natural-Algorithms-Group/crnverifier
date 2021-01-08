# CRN verification module
Are two chemical reaction networks (CRNs) the same? This package provides code
to verify the correctness of an implementation CRN with respect to a formal CRN
using the stochastic
trajectory equivalence notions **CRN bisimulation** [[Johnson et al. (2019)]], 
**pathway decomposition** [[Shin et al.  (2019)]] and preliminary implementations of 
**compositional & integrated hybrid** [[Shin et al.  (2019)]].

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/DNA-and-Natural-Algorithms-Group/crnverifier)](https://github.com/DNA-and-Natural-Algorithms-Group/crnverifier/tags)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/DNA-and-Natural-Algorithms-Group/crnverifier?include_prereleases)](https://github.com/DNA-and-Natural-Algorithms-Group/crnverifier/releases)
[![PyPI version](https://badge.fury.io/py/crnverifier.svg)](https://badge.fury.io/py/crnverifier)
[![PyPI - License](https://img.shields.io/pypi/l/crnverifier)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/DNA-and-Natural-Algorithms-Group/crnverifier.svg?branch=master)](https://travis-ci.com/github/DNA-and-Natural-Algorithms-Group/crnverifier)
[![Codecov branch](https://img.shields.io/codecov/c/github/DNA-and-Natural-Algorithms-Group/crnverifier/master)](https://codecov.io/gh/DNA-and-Natural-Algorithms-Group/crnverifier)


### Installation
```
  $ python setup.py install
```

### Library examples

The following verification functions are currently available:
```py
from crnverifier import (pathway_decomposition_eq,
                         crn_bisimulation_test,
                         modular_crn_bisimulation_test,
                         integrated_hybrid_test, 
                         compositional_hybrid_test)
from crnverifier.utils import parse_crn

# A formal CRN and a corresponding implementation CRN
fcrn = "A -> B"
icrn = "a1 -> b1; x -> a1; x -> b1; y -> b1; y -> a1; x -> a0; a0 -> a1"

# A quick interface to the internal list representation of CRNs, 
# the first CRN contains only formal species (fs), the species
# of the second CRN are not important right now.
fcrn, fs = parse_crn(fcrn, is_file = False)
icrn, _ = parse_crn(icrn, is_file = False)

# Verify whether the two CRNs are pathway decomposition equivalent 
# given the formal species fs:
v = pathway_decomposition_eq([fcrn, icrn], fs)
print(v)

# For the other notions, we may need a (partial) interpretation:
inter = {'a0': 'A',
         'a1': 'A',
         'b1': 'B'}

# Test if there exists a correct CRN bisimulation that interprets 
# the reactions of icrn as reactions of fcrn.
# The (partial) interpretation is optional here.
v, i = crn_bisimulation_test(fcrn, icrn, fs, interpretation = inter)
print(v, i)

# Verify whether icrn is a correct implementation of fcrn using the 
# two supported hybrid notions.
# The (partial) interpretation is required here.
v, i = integrated_hybrid_test(fcrn, icrn, fs, inter)
print(v, i)
v, i = compositional_hybrid_test(fcrn, icrn, fs, inter)
print(v, i)

```

### Commandline examples
For the format *.crn files, see [crnsimulator].

Verify whether two CRNs f.crn and i.crn are pathway decomposition equivalent:
```
  $ crnverifier pathway-decomposition --crns f.crn i.crn --formal-species A B C
```
Compute the formal basis of a single CRN:
```
  $ crnverifier formal-basis --crn i.crn --formal-species A B C
```

Verify whether a correct CRN bisimulation exists to interpret *i.crn* as *f.crn*:
```
  $ crnverifier crn-bisimulation --formal-crn f.crn --implementation-crn i.crn
```
For options, e.g. to provide a partial interpretation, or to choose a more
suitable algorithm for the permissive condition, use 
```
  $ crnverifier --help
```

Verify whether *i.crn* is a correct implementation of *f.crn* using the two supported hybrid notions.
```
  $ crnverifier integrated-hybrid --formal-crn f.crn --implementation-crn i.crn --interpretation itof.crn
  $ crnverifier compositional-hybrid --formal-crn f.crn --implementation-crn i.crn --interpretation itof.crn
```

## Version
0.3

## License
MIT

## Cite
The equivalence notions:
 - Seung Woo Shin, Chris Thachuk, and Erik Winfree (2019) 
    "Verifying chemical reaction network implementations: A pathway decomposition approach"
    [[Shin et al. (2019)]].
 - Robert F. Johnson, Qing Dong, and Erik Winfree (2019)
    "Verifying chemical reaction network implementations: A bisimulation approach"
    [[Johnson et al. (2019)]].

The implementation (a part of the [Nuskell] compiler project):
 - Stefan Badelt, Seung Woo Shin, Robert F. Johnson, Qing Dong, Chris Thachuk, and Erik Winfree (2017)
    "A General-Purpose CRN-to-DSD Compiler with Formal Verification, Optimization, and Simulation Capabilities"
    [[Badelt et al. (2017)]].


[//]: References
[Shin et al. (2019)]: <https://doi.org/10.1016/j.tcs.2017.10.011>
[Johnson et al. (2019)]: <https://doi.org/10.1016/j.tcs.2018.01.002>
[Badelt et al. (2017)]: <https://doi.org/10.1007/978-3-319-66799-7_15>
[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
[Nuskell]: <https://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell>
[crnsimulator]: <https://www.github.com/bad-ants-fleet/crnsimulator>

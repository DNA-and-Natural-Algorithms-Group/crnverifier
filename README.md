# CRN verification tests
Are two chemical reaction networks (CRNs) the same? This package provides code
to verify the correctness of an implementation CRN with respect to a formal CRN
using the stochastic
trajectory equivalence notions **CRN bisimulation** [Johnson et al. (2019)], 
**pathway decomposition** [Shin et al.  (2019)] and preliminary implementations of 
**compositional & integrated hybrid** [Shin et al.  (2019)].

### Installation
```
  $ python setup.py install
```

### Examples

Verify whether *i.crn* is a CRN bisimulation of *f.crn*:

```
  $ crnverifier crn-bisimulation --formal-crn f.crn --implementation-crn i.crn
```
For options, e.g. to provide a partial interpretation, or to choose a more
suitable algorithm for the permissive condition, use 
```
  $ crnverifier --help
```

Verify whether two CRNs f.crn and i.crn are pathway decomposition equivalent:
```
  $ crnverifier pathway-decomposition --crns f.crn i.crn --formal-species A B C
```

Compute the formal basis of a single CRN:
```
  $ crnverifier formal-basis --crn i.crn --formal-species A B C
```

Verify whether *i.crn* is a correct implementation of *f.crn* using hybrid notions (given the partial interpretation *itof.crn*):
```
  $ crnverifier integrated-hybrid --formal-crn f.crn --implementation-crn i.crn --interpretation itof.crn
  $ crnverifier compositional-hybrid --formal-crn f.crn --implementation-crn i.crn --interpretation itof.crn
```

## Version
0.1

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

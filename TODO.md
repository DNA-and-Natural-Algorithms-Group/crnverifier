# Issues:

- The CRN bisimulation loopsearch algorithm is *not* the same as described in the paper.

- The *new* CRN bisimulation reactionsearch algorithm can produce duplicate
  results whenever multiple distinct sequences of interpreting reactions leads
  to the same "ordering" of implementation species. 
  (See unittest for (fcrn = " -> A"; icrn = " -> y; y <=> z; z -> a".)

- The modular CRN bisimulation algorithm cannot handle cases where shared implementation
  species are not provided as partial interpretation. This raises a NotImplementedError.

# Feature requests:

- Provide some interesting examples in the corresponding directory.

- Provide a space-efficient mode for pathway decomposition.

- A complete "compositional hybrid" and "integrated hybrid" algorithm. (The current algorithms cover only edge cases.)

- Provide an algorithm for 'soft constraints' for interpretations. E.g. by generating a set of possible partial interpretations?


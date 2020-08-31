#
#  crnverifier/__init__.py
#  CRN verification
#
__version__ = "v0.2"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .verifier import verify, modular_bisimulation
from .basis_finder import get_formal_basis
from .crn_pathway_equivalence import test as crn_pathway_equivalence_test
from .crn_bisimulation import test as crn_bisimulation_test_

# crn 
# from nuskell.verifier import (get_formal_basis,
#                               crn_pathway_equivalence_test,
#                               crn_bisimulation_test,
#                               crn_bisimulation_test,
#                               crnutils)
# from nuskell.verifier import crn_bisimulation_test
# correct, interpretation = crn_bisimulation_test(fcrn, icrn, interpretation = None, algorithm = 'whole-graph')

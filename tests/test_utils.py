#!/usr/bin/env python
#
#  tests/test_utils.py
#  Original source from the Nuskell compiler project
#
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

import unittest

from crnverifier.utils import assign_crn_species, parse_crn

class Test_utils(unittest.TestCase):
    def test_assign_species(self):
        crn = "A + B -> C + D"
        (crn, fs) = parse_crn(crn, is_file = False)

        crn = "A + B -> B + w1; w1 -> w1 + w2"
        (crn, _) = parse_crn(crn, is_file = False)
        I, W, Wn = assign_crn_species(crn, fs)

        assert I == set()
        assert W == set(['w1', 'w2'])
        assert Wn == set(['w1'])

        crn = """ A + B -> i1 + C + w1; i1 -> D + w2; w1 + w2 -> w3 """
        (crn, _) = parse_crn(crn, is_file = False)
        I, W, Wn = assign_crn_species(crn, fs)

        assert I == set(['i1'])
        assert W == set(['w1', 'w2', 'w3'])
        assert Wn == set(['w1', 'w2'])

    def test_parse_crn(self): 
        crn = """
            A <=> i
            i + B1 <=> j1
            i + B2 <=> j2
            j1 -> C
            j2 -> C
            """
        icrn, isp = parse_crn(crn)
        ecrn = [[['A'], ['i']], 
                [['i'], ['A']], 
                [['i', 'B1'], ['j1']], 
                [['j1'], ['i', 'B1']], 
                [['i', 'B2'], ['j2']],
                [['j2'], ['i', 'B2']],
                [['j1'], ['C']],
                [['j2'], ['C']]]
        esp = set(['B2', 'j1', 'C', 'A', 'B1', 'j2', 'i'])
        assert icrn == ecrn
        assert isp == esp



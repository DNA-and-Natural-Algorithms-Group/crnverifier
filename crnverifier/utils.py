#
#  crnverifier/utils.py
#  Original source from the Nuskell compiler project
#
import logging
log = logging.getLogger(__name__)

import re
import math
from itertools import chain, product
from .crn_parser import parse_crn_file, parse_crn_string

def parse_crn(string, is_file = False, modular = False):
    crn, species = parse_crn_file(string, modular = modular) if is_file else parse_crn_string(string, modular = modular)
    if modular:
        crn = [split_reversible_reactions(m) for m in crn]
    else:
        crn = split_reversible_reactions(crn)
    return crn, set(species.keys())

def split_reversible_reactions(crn):
    """
    Replace every reversible reaction with the two corresponding irreversible
    reactions.
    """
    new = []
    for [r, p, k] in crn:
        assert len(k) == 1 or len(k) == 2
        new.append([r, p])
        if len(k) == 2:
            new.append([p, r])
    return new

def natural_sort(l):
    """
    Sorts a collection in the order humans would expect. Implementation from
    http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    """
    def convert(text): 
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): 
        return [convert(c) for c in re.split('([0-9]+)', str(key))]

    return sorted(l, key=alphanum_key)

def pretty_rxn(rxn):
    return '{} -> {}'.format(' + '.join(natural_sort(rxn[0])), 
                             ' + '.join(natural_sort(rxn[1])))

def pretty_crn(crn):
    for rxn in natural_sort(crn):
        yield pretty_rxn(rxn)

def clean_crn(crn, duplicates = True, trivial = True, inter = None):
    """Takes a crn and removes trivial / duplicate reactions. """
    new = []
    seen = set()
    for [R, P] in crn:
        lR = sorted(interpret(R, inter)) if inter else sorted(R)
        lP = sorted(interpret(P, inter)) if inter else sorted(P)
        tR = tuple(lR)
        tP = tuple(lP)
        if trivial and tR == tP:
            continue
        if duplicates and (tR, tP) in seen:
            continue
        new.append([lR, lP])
        seen.add((tR, tP))
    return new

def crnsize(crn):
    """ The size of a CRN as defined in JDW 2019.
    """
    sp = set().union(*[set().union(*rxn[:2]) for rxn in crn])
    tot = 0
    for r, p in crn:
        for s in sp:
            tot += math.ceil(math.log2(len([x for x in r if x == s]) + 1)) + \
                   math.ceil(math.log2(len([x for x in p if x == s]) + 1))
    return tot + len(sp)

def intersize(inter, fsp, isp):
    """ The size of a CRN interpretation as defined in JDW 2019.
    """
    # TODO: untested
    return sum([math.ceil(math.log2(len([x for x in inter[fsp] if x == s]) + 1)) for f, s in product(fsp, isp)])

def interpret(l, inter):
    """ Replace species with their interpretation. """
    return list(chain(*[inter.get(x, [x]) for x in l]))

def remove_species(crn, const):
    for rxn in crn:
        rxn[0] = [x for x in rxn[0] if x not in const]
        rxn[1] = [x for x in rxn[1] if x not in const]
    return crn

def assign_crn_species(crn, signals):
    """ Returns types of species in a given CRN.

    On the types of species in an implementation CRN:
        - Signal species are implementation species that (are supposed to)
          correspond to formal species in a formal CRN.
        - Fuel species are implementation species that are required for some
          reactions and are assumed to be present, always.
        - Waste species are chemically inert byproducts that never react.
        - Reactive waste species are byproducts of a reaction that can react,
          but only with fuels or other reactive waste species to produce
          exclusively (reactive) waste species. A typical example for reactive
          wastes are so-called garbage collection mechanisms to turn an
          undesired waste into a desired waste.
        - Intermediate species are all other species that do not fall in any
          above category.

    Returns:
        set(): intermediates
        set(): wastes
        set(): reactive wastes
    """
    species = set().union(*[set().union(*rxn[:2]) for rxn in crn])
    # A signal species cannot be considered waste.
    assert isinstance(signals, set)
    nonwastes = signals.copy()
    wastes = set()

    while True:
        # Add x to non-waste if in any reaction x is an reactant while there
        # are non-wastes taking part in the reaction. Reset the outer loop if
        # a new non-waste species is found.
        flag = False
        for x in species:
            if x in nonwastes:
                continue
            for rxn in crn:
                if x in rxn[0] and len(nonwastes & set(rxn[0] + rxn[1])):
                    nonwastes.add(x)
                    flag = True
                    break
        if not flag:
            break

     # Inert and reactive waste species.
    wastes = species - nonwastes

    # An intermediate species 
    intermediates = species - signals - wastes

    # Let's assert here, to find the problems.
    reactive_waste = set()
    for w in list(wastes):
        if any(w in rxn[0] for rxn in crn):
            reactive_waste.add(w)

    return intermediates, wastes, reactive_waste

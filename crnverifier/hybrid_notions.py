#
#  crnverifier/hybrid_notions.py
#  Original source from the Nuskell compiler project
#
#  Authors:
#   - Seung Woo Shin (seungwoo.theory@gmail.com)
#   - Stefan Badelt (bad-ants-fleet@posteo.eu)
#
import logging
log = logging.getLogger(__name__)

from itertools import chain
from collections import Counter

from .utils import assign_crn_species, pretty_crn, natural_sort
from .pathway_decomposition import get_formal_basis, NoFormalBasisError, clean_crn
from .crn_bisimulation import crn_bisimulation_test

def integrated_hybrid_test(*kargs):
    return integrated_hybrid_dev1(*kargs)

def compositional_hybrid_test(*kargs):
    return compositional_hybrid_dev1(*kargs)

def integrated_hybrid_dev1(fcrn, icrn, fs, inter, modular):
    # SWS implementation
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    assert all(x in fs2 for x in fs) 

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(reactive_waste):
        log.warning(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.warning(f'{len(wastes)} waste species are treated as formal: ({", ".join(wastes)})')
        for x in wastes: 
            inter[x] = []

    log.debug('Formal CRN with formal species: {}\n  {}'.format(
        ", ".join(natural_sort(fs)),
        "\n  ".join(pretty_crn(fcrn))))
    log.debug('Implementation CRN with formal species: {}\n  {}'.format(
        ", ".join(natural_sort(inter.keys())),
        "\n  ".join(pretty_crn(icrn))))
    log.debug('Implementation CRN after interpretateion:\n  {}'.format(
        "\n  ".join(pretty_crn(clean_crn(icrn, inter = inter)))))

    try:
        log.debug(f'Formal species to find basis: {set(inter.keys())}')
        fbasis_raw, fbasis_int = get_formal_basis(icrn, set(inter.keys()), modular = modular, interpretation = inter)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False

    log.debug('Raw formal basis:\n  {}'.format("\n  ".join(pretty_crn(fbasis_raw))))
    log.debug('Interpreted formal basis:\n  {}'.format("\n  ".join(pretty_crn(fbasis_int))))
    return sorted(fcrn) == sorted(clean_crn(fbasis_int))

def compositional_hybrid_dev1(fcrn, icrn, fs, inter):
    # SWS implementation
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    assert all(x in fs2 for x in fs) 

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(reactive_waste):
        log.warning(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.warning(f'{len(wastes)} waste species are treated as formal: ({", ".join(wastes)})')
        for x in wastes: 
            inter[x] = []

    log.debug('Formal CRN with formal species: {}\n  {}'.format(
        ", ".join(natural_sort(fs)),
        "\n  ".join(pretty_crn(fcrn))))
    log.debug('Implementation CRN with formal species: {}\n  {}'.format(
        ", ".join(natural_sort(inter.keys())),
        "\n  ".join(pretty_crn(icrn))))

    try:
        log.debug(f'Formal species to find basis: {set(inter.keys())}')
        fbasis_raw, _ = get_formal_basis(icrn, set(inter.keys()), modular = modular)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False

    log.debug('Raw formal basis:\n  {}'.format("\n  ".join(pretty_crn(fbasis_raw))))
    return sorted(fcrn) == sorted(clean_crn(fbasis_raw, inter = inter))


def hybrid(c1, c2, inter, compositional = False, integrated = False):
    """ Test two CRNs for pathway equivalence.

    Args:
        c1 (list, list): Tuple of formal CRN and formal species.
        c2 (list, list): Tuple of implementation CRN and signal species.
        inter (dict): An interpretation of fs2 in terms of fs1.
        compositional (bool, optional): Use compositional hybrid notion.
            Defaults to False.
        integrated (bool, optional): Use integrated hybrid notion.
            Defaults to False.

    Returns:
        True: if formal basis of crn1 and formal basis of crn2 are equivalent.
        False: otherwise.
    """
    DEVELMODE = True
    # Use the *interpretation upfront* mode of integrated hybrid and the new
    # compositional hybrid implementation.
    if DEVELMODE:
        (integrated, integrated2) = (False, integrated) if integrated else (False, False)
        (compositional, compositional2) = (False, compositional) if compositional else (False, False)

    # Preprocess arguments, just to make sure ...
    crn1, fs1 = clean_crn([rxn[0:2] for rxn in c1[0]]), set(c1[1])
    if integrated2: 
        log.warning('Using integrated "interpretation upfront" mode.')
        crn2, fs2 = clean_crn([rxn[0:2] for rxn in c2[0]], inter = inter), set(c1[1])
    else:
        crn2, fs2 = clean_crn([rxn[0:2] for rxn in c2[0]]), set(c2[1])
        # Note: if inter provides interpretations for non-signal species, that's ok here.
        # They will be used as formal species when finding the formal basis.
        assert all(x in inter for x in fs2)

    assert isinstance(inter, dict)
    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(crn2, fs2)
    if len(reactive_waste):
        log.warning(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.warning(f'{len(wastes)} waste species are treated as formal: ({", ".join(wastes)})')
        for x in wastes: 
            inter[x] = []

    log.debug("Formal CRN:")
    [log.debug('    {}'.format(r)) for r in genCRN(crn1, rates = False)]
    log.debug("")
    log.debug(f"Implementation CRN with {len(set().union(*[set().union(*rxn) for rxn in crn2]))} species and {len(crn2)} reactions.")
    [log.debug(f'    {r}') for r in genCRN(crn2, rates = False)]
    log.debug("")
    if integrated:
        log.debug(f"Implementation CRN after partial interpretation:")
        [log.debug(f'    {r}') for r in genCRN(crn2, interpretation = inter, rates = False)]
        log.debug("")

    try:
        log.debug(f'Formal species to find basis: {fs2 | wastes}')
        fbasis_raw, fbasis_int = get_formal_basis(crn2, fs2 | wastes, modular = True,
                                            interpretation = inter if integrated else None)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False

    log.debug(f"Raw formal basis:")
    [log.debug(f'    {r}') for r in genCRN(fbasis_raw, rates = False)]
    log.debug(f"Interpreted formal basis:")
    [log.debug(f'    {r}') for r in genCRN(fbasis_int, rates = False)]

    if compositional:
        return sorted(crn1) == sorted(clean_crn(fbasis_raw, inter = inter))
    elif compositional2:
        # Currently, compositional and integrated use the same raw basis,
        # because the raw basis after interpretation is the interpreted basis.
        fcrn = [[Counter(part) for part in rxn] for rxn in crn1]
        icrn = [[Counter(part) for part in rxn] for rxn in fbasis_raw]

        interC = dict()
        for k, v in inter.items():
            interC[k] = Counter(v)

        v, i = get_crn_bisimulation(fcrn, icrn, fs1, 
                                    interpretation = interC, 
                                    permissive = 'whole-graph')
        return v
    elif integrated:
        # If you use the standard integrated implementation, then you need
        # to use the interpreted formal base.
        return sorted(crn1) == sorted(clean_crn(fbasis_int))

    # Pure CRN pathway decomposition or *new* integrated implementation.
    return sorted(crn1) == sorted(fbasis_raw)



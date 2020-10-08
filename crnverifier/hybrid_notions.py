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

from .utils import assign_crn_species, pretty_crn, natural_sort
from .pathway_decomposition import get_formal_basis, NoFormalBasisError, clean_crn
from .crn_bisimulation import crn_bisimulation_test

class HybridNotionError(Exception):
    pass

def integrated_hybrid_test(*kargs):
    # Wrapper function to use the recommended standard implementation
    return integrated_hybrid_dev1(*kargs)

def compositional_hybrid_test(*kargs):
    # Wrapper function to use the recommended standard implementation
    return compositional_hybrid_dev1(*kargs)

def integrated_hybrid_dev1(fcrn, icrn, fs, inter, modular = True):
    """ Integrated hybrid test - using icrn during path enumeration.
    """
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    if not all(x in fs2 for x in fs):
        raise HybridNotionError('Some formal species cannot be found in the interpretation.')

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(reactive_waste):
        log.debug(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.info(f'{len(wastes)} waste species are treated as formal: ' + \
                 f'({", ".join(wastes)})')
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
        fbasis_raw, fbasis_int = get_formal_basis(icrn, 
                                                  set(inter.keys()),
                                                  modular = modular,
                                                  interpretation = inter)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False, inter

    log.debug('Raw formal basis:\n  {}'.format(
        "\n  ".join(pretty_crn(fbasis_raw))))
    log.debug('Interpreted formal basis:\n  {}'.format(
        "\n  ".join(pretty_crn(fbasis_int))))
    return sorted(fcrn) == sorted(clean_crn(fbasis_int)), inter

def integrated_hybrid_dev2(fcrn, icrn, fs, inter, modular = True):
    """ Integrated hybrid test - interpreted icrn during path enumeration.
    """
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn, inter = inter)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    if not all(x in fs2 for x in fs):
        raise HybridNotionError('Some formal species cannot be found in the interpretation.')

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(wastes):
        log.warning(f'{len(wastes)} waste species are treated as formal: ' + \
                    f'({", ".join(wastes)})')
        for x in wastes: 
            inter[x] = []
        # Now update icrn so that wastes are removed!
        icrn = clean_crn(icrn, inter = inter)

    log.debug('Formal CRN with formal species: {}\n  {}'.format(
        ", ".join(natural_sort(fs)),
        "\n  ".join(pretty_crn(fcrn))))
    log.debug('Implementation CRN after interpretateion:\n  {}'.format(
        "\n  ".join(pretty_crn(icrn))))

    try:
        log.debug(f'Formal species to find basis: {set(inter.keys())}')
        fbasis_raw, fbasis_int = get_formal_basis(icrn,
                                                  set(inter.keys()),
                                                  modular = modular)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False, inter

    log.debug('Interpreted formal basis:\n  {}'.format(
        "\n  ".join(pretty_crn(fbasis_raw))))
    return sorted(fcrn) == sorted(fbasis_raw), inter

def compositional_hybrid_dev1(fcrn, icrn, fs, inter, modular = True):
    """ Limited compositional hybrid test - using strong bisimulation.
    """
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    if not all(x in fs2 for x in fs):
        raise HybridNotionError('Some formal species cannot be found in the interpretation.')

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(reactive_waste):
        log.debug(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.info(f'{len(wastes)} waste species are treated as formal: ' + \
                 f'({", ".join(wastes)})')
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
        return False, inter

    log.debug('Raw formal basis:\n  {}'.format("\n  ".join(pretty_crn(fbasis_raw))))
    return sorted(fcrn) == sorted(clean_crn(fbasis_raw, inter = inter)), inter

def compositional_hybrid_dev2(fcrn, icrn, fs, inter, modular = True):
    """ Limited compositional hybrid test - using crn_bisimulation_test.
    """
    fcrn = clean_crn(fcrn)
    icrn = clean_crn(icrn)

    # Note: if inter provides interpretations for non-signal species, that's ok here.
    # They will be used as formal species when finding the formal basis.
    fs2 = list(chain(*inter.values()))
    if not all(x in fs2 for x in fs):
        raise HybridNotionError('Some formal species cannot be found in the interpretation.')

    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(icrn, set(inter.keys()))
    if len(reactive_waste):
        log.debug(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.info(f'{len(wastes)} waste species are treated as formal: ' + \
                 f'({", ".join(wastes)})')
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
        return False, inter

    log.debug('Raw formal basis:\n  {}'.format("\n  ".join(pretty_crn(fbasis_raw))))

    v, i = crn_bisimulation_test(fcrn, icrn, fs, 
                                 interpretation = inter, 
                                 permissive = 'whole-graph')
    return v, inter


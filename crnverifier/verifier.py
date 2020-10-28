#!/usr/bin/env python
#
#  crnverifier/verifier.py
#  Original source from the Nuskell compiler project
#
import logging
log = logging.getLogger(__name__)

import sys
import signal
import argparse

from crnverifier import __version__
from crnverifier.utils import (parse_crn, 
                               pretty_crn,
                               clean_crn,
                               crnsize,
                               remove_species, 
                               natural_sort)
from crnverifier.crn_bisimulation import (crn_bisimulation_test, 
                                          modular_crn_bisimulation_test)
from crnverifier.pathway_decomposition import (NoFormalBasisError, 
                                               get_formal_basis, 
                                               pathway_decomposition_eq)
from crnverifier.hybrid_notions import (integrated_hybrid_test,
                                        compositional_hybrid_test)

class TimeoutError(Exception):
    pass

def handler(signum, frame):
    raise TimeoutError('Time over')

def limit_runtime(timeout, func, *kargs, **kwargs):
    """ Helper function to exit a function and return None after a specified time.
    """
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    try:
        if timeout:
            log.info(f'Function will abort after {timeout} seconds.')
        ret = func(*kargs, **kwargs)
    except TimeoutError:
        ret = None
    finally:
        tleft = signal.alarm(0)
        if timeout:
            log.info(f'Function needed {int(timeout-tleft)} seconds.')
    return ret

def get_bisimulation_inputs_interactive(constants):
    raise NotImplementedError('no interactive mode support.')

def get_pathway_decomposition_inputs_interactive(constants):
    raise NotImplementedError('no interactive mode support.')

def get_bisimulation_inputs(fCRN, iCRN, interpretation, constants):
    """ Helper function to parse commandline input files.
    """
    fcrn, fs = parse_crn(fCRN, is_file = True)
    fcrn = remove_species(fcrn, constants)
    log.info('Input formal CRN of size {} with formal species: {}\n  {}'.format(
        crnsize(fcrn),
        ", ".join(natural_sort(fs)),
        "\n  ".join(pretty_crn(fcrn))))

    icrn, _ = parse_crn(iCRN, is_file = True)
    icrn = remove_species(icrn, constants)
    log.info('Input implementation CRN:\n  {}'.format("\n  ".join(pretty_crn(icrn))))

    inter = None
    if interpretation:
        inte, _ = parse_crn(interpretation, is_file = True)
        assert all(len(rxn[0]) == 1 for rxn in inte)
        inter = {rxn[0][0]: rxn[1] for rxn in inte}
        log.info('Input interpretation:\n  {}'.format(
                 '\n  '.join(f"{k} -> {' + '.join(v)}" for k, v in natural_sort(inter.items()))))
    return fcrn, icrn, fs, inter

def get_modular_bisimulation_inputs(fCRN, iCRN, interpretation, constants):
    """ Helper function to parse commandline input files.
    """
    fcrns, fs = parse_crn(fCRN, is_file = True, modular = True)
    fcrns = [remove_species(m, constants) for m in fcrns]
    for e, m in enumerate(fcrns, 1):
        log.info('Input formal CRN module: {}\n  {}'.format(e,
            "\n  ".join(pretty_crn(m))))

    icrns, _ = parse_crn(iCRN, is_file = True, modular = True)
    icrn = [remove_species(m, constants) for m in icrns]
    for e, m in enumerate(icrns, 1):
        log.info('Input implementation CRN module: {}\n  {}'.format(e,
            "\n  ".join(pretty_crn(m))))

    inter = None
    if interpretation:
        inte, _ = parse_crn(interpretation, is_file = True)
        assert all(len(rxn[0]) == 1 for rxn in inte)
        inter = {rxn[0][0]: rxn[1] for rxn in inte}
        log.info('Input interpretation:\n  {}'.format(
                 '\n  '.join(f"{k} -> {' + '.join(v)}" for k, v in natural_sort(inter.items()))))
    return fcrns, icrns, fs, inter

def print_bisimulation_outputs(v, i, method, verify_timeout):
    if v is True:
        print(f"Verification result for {method} = {v}.",
              f"The implementation CRN is a correct implementation of the formal CRN.")
    elif v is False:
        print(f"Verification result for {method} = {v}.",
              f"The implementation CRN is not a correct implementation of the formal CRN.")
    elif v is None:
        print(f"No verification result for {method}.",
              f"Verification did not terminate within {verify_timeout} seconds.")

def get_pathway_decomposition_inputs(crn_files, formals, constants):
    """ Helper function to parse commandline input files.
    """
    log.info(f'Formal species: {{{", ".join(natural_sort(formals))}}}')
    log.info(f'Constant species: {{{", ".join(natural_sort(constants))}}}')

    crns = []
    for e, cf in enumerate(crn_files, 1):
        crn, _ = parse_crn(cf, is_file = True)
        crn = remove_species(crn, constants)
        log.info('Input CRN {}:\n  {}'.format(e, "\n  ".join(pretty_crn(crn))))
        crns.append(crn)
    return crns, set(formals)

def add_commandline_arguments(parser):
    correctness = parser.add_argument_group(
            'Arguments for crn-bisimulation & hybrid notions')
    equivalence = parser.add_argument_group(
            'Arguments for pathway decomposition equivalence')

    parser.add_argument("method", 
            action = 'store', 
            choices = ('crn-bisimulation', 'modular-crn-bisimulation', 
                       'pathway-decomposition', 'formal-basis', 
                       'integrated-hybrid', 
                       'compositional-hybrid'),
            help = """Specify verification method.  The method
            *modular-crn-bisimulation* requires a more constrained text input
            format for CRNs: each module must be written on one line. For
            example: "A + B -> C; C -> D + E" written on one line is interpeted
            as one module composed of two reactions. If the formal CRN has n
            modules, and the implementation CRN n+1 modules, then the last
            implementation CRN module is assumed to implement crosstalk, i.e.
            the corresponding formal CRN module is empty. """)

    parser.add_argument("--version", action = 'version', 
                        version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = "Print verbose output. -vv increases verbosity level.")
    parser.add_argument("--verify-timeout", type = int, default = 30, metavar = '<int>',
            help = "Specify time in seconds to wait for verification to complete.")
    parser.add_argument("-c", "--constant-species", nargs = '+', default = [], 
            action = 'store', metavar = '<str>', 
            help = """Provide a list of constant species that can be 
            removed before testing equivalence, 
            e.g. fuel species and waste species.""")
    parser.add_argument("--non-modular", action = 'store_true',
            help = """Do not split CRNs into smaller modules when using the
            pathway decomosition method.""")
    parser.add_argument("--profile", action = 'store_true',
            help = """Get some code profiling information (requires statprof-smarkets).""")
 
    correctness.add_argument("-f", "--formal-crn", 
            action = 'store', metavar = '</path/to/file>',
            help = "Read a formal CRN from a file.")
    correctness.add_argument("-i", "--implementation-crn", 
            action = 'store', metavar = '</path/to/file>',
            help = "Read an implementation CRN from a file.")
    correctness.add_argument("-m", "--interpretation", 
            action = 'store', metavar = '</path/to/file>',
            help = """Read an interpretation (a mapping between species in
            typical CRN format) from a file.""")
    correctness.add_argument("--permissive-check",
            action = 'store', metavar = '<str>', default = 'graphsearch',
            choices = ('loopsearch', 'reactionsearch', 'graphsearch'),
            help = "Choose an algorithm to check the permissive condition.")

    equivalence.add_argument("-b", "--crn-files", nargs = '+', default = [],
            action='store', metavar='</path/to/file>',
            help = "Read one or more files that contain a CRN.")
    equivalence.add_argument("-s", "--formal-species", nargs = '+', 
            default = [], action = 'store', metavar = '<str>', 
            help = "List formal species names.")

def main():
    """ Test if an implementation CRN is a correct implementation of a formal CRN.

    An example from Robert's CRN bisimulation code:
    fCRN = "a -> b"
    icrn = "a1 -> b1; x -> a1; x -> b1; y -> b1; y -> a1; x -> a0; a0 -> a1"


    """
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    add_commandline_arguments(parser)
    args = parser.parse_args()

    logger = logging.getLogger('crnverifier')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s %(filename)s: %(message)s')
    ch.setFormatter(formatter)
    if args.verbose == 0:
        ch.setLevel(logging.WARNING)
    elif args.verbose == 1:
        ch.setLevel(logging.INFO)
    elif args.verbose == 2:
        ch.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
        ch.setLevel(logging.NOTSET)
    logger.addHandler(ch)

    if args.profile:
        try:
            import statprof
            statprof.start()
        except ImportError as err:
            log.warning('Cannot import statprof module. (pip install statprof-smarkets.)')
            args.profile = False

    if args.method in ('crn-bisimulation', 'integrated-hybrid', 'compositional-hybrid'):
        # TODO: Interactive mode asks for input if fcrn isn't given, etc.
        assert args.formal_crn is not None
        assert args.implementation_crn is not None
        fcrn, icrn, fs, inter = get_bisimulation_inputs(args.formal_crn,
                                                    args.implementation_crn,
                                                    args.interpretation,
                                                    args.constant_species)
        if args.method == 'crn-bisimulation':
            v = limit_runtime(args.verify_timeout,
                              crn_bisimulation_test,
                              fcrn, icrn, fs, 
                              interpretation = inter, 
                              permissive = args.permissive_check)
        elif args.method == 'compositional-hybrid':
            v = limit_runtime(args.verify_timeout,
                              compositional_hybrid_test,
                              fcrn, icrn, fs, inter, 
                              not args.non_modular)
 
        else:
            assert args.method == 'integrated-hybrid'
            v = limit_runtime(args.verify_timeout,
                              integrated_hybrid_test,
                              fcrn, icrn, fs, inter, 
                              not args.non_modular)

        # Unpack the verification results if it wasn't a timeout.
        (v, i) = (None, None) if v is None else v
        if v:
            log.info('Returned interpretation:\n  ' + \
                     '\n  '.join(f"{k} -> {' + '.join(v)}" for k, v in natural_sort(i.items())))
            log.info('Interpreted CRN:\n  {}'.format( 
                     '\n  '.join(pretty_crn(clean_crn(icrn, inter = i)))))
        print_bisimulation_outputs(v, i, args.method, args.verify_timeout)

    elif args.method in ('modular-crn-bisimulation'):
        # TODO: Interactive mode asks for input if fcrn isn't given, etc.
        assert args.formal_crn is not None
        assert args.implementation_crn is not None
        fcrns, icrns, fs, inter = get_modular_bisimulation_inputs(args.formal_crn,
                                                                  args.implementation_crn,
                                                                  args.interpretation,
                                                                  args.constant_species)
        v = limit_runtime(args.verify_timeout,
                          modular_crn_bisimulation_test,
                          fcrns, icrns, fs, 
                          interpretation = inter, 
                          permissive = args.permissive_check)

        # Unpack the verification results if it wasn't a timeout.
        (v, i) = (None, None) if v is None else v
        if v:
            log.info("Returned interpretation:\n  " + \
                     '\n  '.join(f"{k} -> {' + '.join(v)}" for k, v in natural_sort(i.items())))
            for e, icrn in enumerate(icrns, 1):
                log.info('Interpreted CRN module {}:\n  {}'.format(e, 
                         '\n  '.join(pretty_crn(clean_crn(icrn, inter = i)))))
        print_bisimulation_outputs(v, i, args.method, args.verify_timeout)

    elif args.method in ('formal-basis', 'pathway-decomposition'):
        # TODO: Interactive mode asks for input if fcrn isn't given, etc.
        if args.method == 'pathway-decomposition':
            assert len(args.crn_files) > 1
        assert len(args.formal_species)
        crns, fs = get_pathway_decomposition_inputs(args.crn_files,
                                                    args.formal_species,
                                                    args.constant_species)
        if args.method == 'pathway-decomposition':
            v = limit_runtime(args.verify_timeout,
                              pathway_decomposition_eq,
                              crns, fs, not args.non_modular)
            if v is True:
                print(f"Verification result for {args.method} = {v}. The CRNs are equivalent.")
            elif v is False:
                print(f"Verification result for {args.method} = {v}. The CRNs are not equivalent.")
            elif v is None:
                print(f"No verification result for {args.method}.",
                      f"Verification did not terminate within {args.verify_timeout} seconds.")
        elif args.method == 'formal-basis':
            for e, crn in enumerate(crns, 1):
                try:
                    v = limit_runtime(args.verify_timeout,
                                      get_formal_basis,
                                      crn, fs, 
                                      not args.non_modular)

                    if v is None: 
                        print('Timeout, no formal basis found.')
                    else:
                        fbasis, _ = v
                        print('Formal basis {}:\n  {}'.format(e, "\n  ".join(pretty_crn(fbasis))))
                except NoFormalBasisError as err:
                    print("Could not find formal basis {}: {}".format(e, err))
    else:
        raise NotImplementedError(args.method)

    if args.profile:
        statprof.stop()
        statprof.display()

if __name__ == "__main__":
    main()

#!/usr/bin/env python

import random
import argparse

def rand_species(n, prefix = 'X'):
    return prefix + str(random.choice(range(1, n + 1)))

def rand_state(s, n):
    n = random.choice(range(1, n + 1))
    sp = []
    for i in range(n):
        sp.append(rand_species(s))
    return ' + '.join(sp)

def arrow(mode):
    if mode == 'irreversible':
        return '->'
    elif mode == 'reversible':
        return '<=>'
    return '->' if random.choice([1, 2]) == 2 else '<=>'

def generate_crn(num_species,
                 num_reactions,
                 max_reactants = 2,
                 max_products = 2,
                 arrow_choice = 'random'):
    for i in range(0, num_reactions):
        yield "{} {} {}".format(
            rand_state(num_species, max_reactants),
            arrow(arrow_choice),
            rand_state(num_species, max_products))

def main():
    """ Generate a random CRN """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--species", type=int, default=5,
                        help="Total number of species.")
    parser.add_argument("-r", "--reactions", type=int, default=10,
                        help="Total number of reactions.")
    parser.add_argument("-e", "--max-reactants", type=int, default=2,
                        help="Maximum number of reactants per reaction.")
    parser.add_argument("-p", "--max-products", type=int, default=2,
                        help="Maximum number of products per reaction.")
    parser.add_argument("-a", "--arrows", action = "store", default = 'random',
                        choices = ('irreversible', 'reversible', 'random'),
                        help="Specify types of reaction arrows.")
    args = parser.parse_args()

    for rxn in generate_crn(args.species, 
                            args.reactions,
                            args.max_reactants, 
                            args.max_products, 
                            args.arrows):
        print(rxn)

if __name__ == '__main__':
    main()

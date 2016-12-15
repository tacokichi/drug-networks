"""
Project 4 - tanimoto.py
November 6, 2016
@author: Hataka
"""

import sys as sys
import getopt as opt
import chemoUtils as utils
import time

start_time = time.time()


# function to initialize options and arguments
# outputs:
#   n: number of iterations
#   r: seed for random number generator
#   drug_file: name of file containing drug data
#   target_file: name of file containing target data
#   protein_a: protein A
#   protein_b: protein B
def initialize():
    # default options
    n = 100
    r = 214

    opts, args = opt.getopt(sys.argv[1:], "n:r:")
    
    # store non-default options
    for option in opts:
        if option[0] == "-n":
            n = option[1]
        if option[0] == "-r":
            r = option[1]

    # store arguments
    drug_file = args[0]
    target_file = args[1]
    protein_a = args[2]
    protein_b = args[3]

    return int(n), int(r), drug_file, target_file, protein_a, protein_b

    
# main method  
def main():
    # initialize options and arguments
    n, r, drug_file, target_file, protein_a, protein_b = initialize()

    # read drug and target files
    drugs = utils.read_file(drug_file)
    targets = utils.reformate_targets(drugs, utils.read_file(target_file))
    
    # determine drugs that bind to proteins A and B
    drugs_a = utils.get_binding_drugs(targets, protein_a)
    drugs_b = utils.get_binding_drugs(targets, protein_b)
    
    # construct pairs of drugs
    pairs = utils.get_pairs(drugs_a, drugs_b)
    
    # calculate the non-random tanimoto summary   
    tanimotos_sum = utils.calculate_tanimoto_summary(drugs, pairs, 0.5)
    
    # calculate the p-value from pseudo-random tanimoto summaries
    print utils.calculate_bootstrap_p(drugs, len(drugs_a), len(drugs_b), tanimotos_sum, r, n, 0.5)
    print "run time:", time.time()-start_time

if  __name__ == '__main__':
    main()

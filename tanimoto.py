"""
Project 4 - tanimoto.py
November 6, 2016
@author: Hataka
"""

import sys as sys
import chemoUtils as utils
import time

start_time = time.time()

    
# function to write to output file
# inputs:
#   output_file: name of output file
#   drugs: 2D numpy array containing drug IDs, generic names, and fingerprints
#   targets: 2D numpy array containing drug IDs, uniprot accessions, and target proteins
#   pairs: 2D numpy array of pairs of drugs
#   tanimotos: 1D numpy array of tanimoto coefficients for pairs of drugs
def write_file(output_file, drugs, targets, pairs, tanimotos):
    with open(output_file, "w") as f:
        for i in range(len(pairs)):
            f.write(drugs[pairs[i,0]][0])
            f.write(",")
            f.write(drugs[pairs[i,1]][0])
            f.write(",")
            f.write("%.6f" %tanimotos[i])
            f.write(",")    
            f.write(str(utils.contains_same_target(targets, pairs[i,0], pairs[i,1])))
            f.write("\n")


# main method  
def main():
    # read input arguments
    drug_file = sys.argv[1]
    target_file = sys.argv[2]
    output_file = sys.argv[3]

    # read drug and target files
    drugs = utils.read_file(drug_file)
    targets = utils.reformate_targets(drugs, utils.read_file(target_file))
    
    # find all unique pairs of drugs
    pairs = utils.get_pairs_to_self(drugs)
    
    # calculate all tanimoto coefficients for pairs of drugs
    tanimotos = utils.calculate_tanimotos(drugs[:,2], pairs)
    
    # write to output file
    write_file(output_file, drugs, targets, pairs, tanimotos)
    print "run time:", time.time()-start_time

if  __name__ == '__main__':
    main()

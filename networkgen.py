"""
Project 4 - networkgen.py
November 6, 2016
@author: Hataka
"""

import sys as sys
import numpy as np
import chemoUtils as utils
import time

start_time = time.time()


# function to calculate bootstrap p-values for every protein pair
# inputs:
#   drugs: 2D array of all drugs
#   targets: 2D array of all targets
#   nodes: 2D array of protein nodes
#   protein_pairs: 2D array of pairs of proteins
#   n: number of iterations
#   r: seed for random number generator
#   cutoff: threshold tanimoto summary
# outputs:
#   returns: 1D array of p-values
def calculate_p_values(drugs, targets, nodes, protein_pairs, n, r, cutoff):
    p_values = []
    for pair in protein_pairs:
        # determine drugs that bind to proteins A and B
        drugs_a = utils.get_binding_drugs(targets, nodes[pair[0]][0])
        drugs_b = utils.get_binding_drugs(targets, nodes[pair[1]][0])

        # construct pairs of drugs
        drug_pairs = utils.get_pairs(drugs_a, drugs_b)
        
        # calculate the non-random tanimoto summary
        tanimotos_sum = utils.calculate_tanimoto_summary(drugs, drug_pairs, cutoff)
        
        # calculate the p-value from pseudo-random tanimoto summaries
        p_values.append(utils.calculate_bootstrap_p(drugs, len(drugs_a), len(drugs_b), tanimotos_sum, r, n, cutoff))
    return np.array(p_values)


# function to sort protein pairs within and between rows
# inputs:
#   nodes: 2D array of protein nodes
#   protein_pairs: 2D array of pairs of proteins
# outputs:
#   returns: pairs of proteins sorted by their names
def sort_to_sif_pairs(nodes, protein_pairs):
    sif_pairs = []
    for pair in protein_pairs:
        sif_pairs.append([nodes[pair[0]][0], nodes[pair[1]][0]])
    sif_pairs = np.sort(np.array(sif_pairs), axis = 1)
    return np.sort(sif_pairs, axis = 0)
    
    
# function to create and initialize output files
def create_files():
    with open("network.sif", "w") as f:
        f.write("")
    with open("name.nodeAttr", "w") as f:
        f.write("name\n")
    with open("indication.nodeAttr", "w") as f:
        f.write("indication\n")


# function to write to output files
# inputs:
#   nodes: 2D array of protein nodes
#   protein_pairs: 2D array of pairs of proteins
def write_files(nodes, protein_pairs):
    sif_pairs = sort_to_sif_pairs(nodes, protein_pairs)
    for pair in sif_pairs:
        with open("network.sif", "a") as f:
            f.write(pair[0])
            f.write(" edge ")
            f.write(pair[1])
            f.write("\n")
            
    unique_proteins = np.unique(protein_pairs)
    for protein in unique_proteins:
        with open("name.nodeAttr", "a") as f:
            f.write(nodes[protein][0])
            f.write(" = ")
            f.write(nodes[protein][1])
            f.write("\n")
        with open("indication.nodeAttr", "a") as f:
            f.write(nodes[protein][0])
            f.write(" = ")
            f.write(nodes[protein][2])
            f.write("\n")
        

# main method
def main():
    # read input arguments
    drug_file = sys.argv[1]
    target_file = sys.argv[2]
    node_file = sys.argv[3]

    # read drug and target files
    drugs = utils.read_file(drug_file)
    targets = utils.reformate_targets(drugs, utils.read_file(target_file))
    nodes = utils.read_file(node_file)
    
    # construct pairs of proteins
    protein_pairs = utils.get_pairs_to_self(nodes)
    
    # calculate p-values
    p_values = calculate_p_values(drugs, targets, nodes, protein_pairs, 100, 214, 0.5)
    p_indices = np.argwhere(np.array(p_values) <= 0.05).flatten()
    
    # create and write output files
    create_files()
    write_files(nodes, protein_pairs[p_indices])
    print "run time:", time.time()-start_time

if  __name__ == '__main__':
    main()

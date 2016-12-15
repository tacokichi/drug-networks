"""
Project 4 - chemoUtils.py
November 6, 2016
@author: Hataka
"""

import numpy as np


# function to calculate tanimoto coefficient between two fingerprints
# inputs:
#   fingerprint1: 1D array of drug 1's fingerprint
#   fingerprint2: 1D array of drug 2's fingerprint
# outputs:
#   returns: tanimoto coefficient
def calculate_tanimoto(fingerprint1, fingerprint2):
    set1 = set(fingerprint1.split())
    set2 = set(fingerprint2.split())
    return float(len(set1.intersection(set2))) / len(set1.union(set2))


# function to calculate tanimoto coefficients for pairs of drugs
# inputs:
#   fingerprints: 2D array of fingerprints for all drugs
#   pairs: 2D array of pairs of drugs
# outputs:
#   returns: 1D array of tanimoto coefficients
def calculate_tanimotos(fingerprints, pairs):
    tanimotos = []
    for pair in pairs:
        tanimotos.append(calculate_tanimoto(fingerprints[pair[0]], fingerprints[pair[1]]))
    return np.array(tanimotos)
    

# function to calculate a tanimoto summary for pairs of drugs
# inputs:
#   drugs: 2D array of all drugs
#   pairs: 2D array of pairs of drugs
#   cutoff: threshold tanimoto coefficient
# outputs:
#   returns: tanimoto summary
def calculate_tanimoto_summary(drugs, pairs, cutoff):
    tanimotos = calculate_tanimotos(drugs[:,2], pairs)  
    return np.sum(tanimotos[np.argwhere(tanimotos > cutoff)])


# function to determine whether two drugs bind to the same target
# inputs:
#   targets: 2D array of all targets
#   drug1: integer of drug 1
#   drug2: integer of drug 2
# outputs:
#   returns: 1 (True) or 0 (False)
def contains_same_target(targets, drug1, drug2):
    if len(set(targets[drug1][:,2]).intersection(set(targets[drug2][:,2]))) != 0:
        return 1
    else:
        return 0


# function to calculate a bootstrap p-value
# inputs:
#   drugs: 2D array of all drugs
#   na: number of drugs that bind protein A
#   nb: number of drugs that bind protein B
#   r: seed for the random number generator
#   n: number of iterations to calculate the p-value
#   cutoff: threshold tanimoto coefficient
# outputs:
#   returns: bootstrap p-value
def calculate_bootstrap_p(drugs, na, nb, real_sum, r, n, cutoff):
    np.random.seed(seed = r)
    ndrugs = len(drugs)
    
    p_count = 0
    for i in range(n):
        drugs_a = np.random.randint(0, ndrugs, (na))
        drugs_b = np.random.randint(0, ndrugs, (nb))
        pairs = get_pairs(drugs_a, drugs_b)
        tanimotos_sum = calculate_tanimoto_summary(drugs, pairs, cutoff)
        if tanimotos_sum > real_sum:
            p_count += 1

    return float(p_count) / n
    

# function to return drugs that bind to a target protein
# inputs:
#   targets: 2D array of all targets
#   target: target protein
# outputs:
#   returns: 1D array of drugs
def get_binding_drugs(targets, target):
    drugs = []
    for i in range(len(targets)):
        if target in targets[i][:,1]:
            drugs.append(i)
    return np.array(drugs)


# function to return non-redundant pairs of items from one list
# inputs:
#   lst: 1D array of items
# outputs:
#   returns: 2D array of paired items
def get_pairs_to_self(lst):
    pairs = []
    for item1 in range(len(lst)):
        for item2 in range(len(lst)):
            if item1 != item2 and item1 < item2:
                pairs.append([item1, item2])
    return np.array(pairs)


# function to return non-redundant pairs of items from two lists
# inputs:
#   lst1: 1D array of items
#   lst2: 1D array of other items
# outputs:
#   returns: 2D array of paired items
def get_pairs(lst1, lst2):
    return np.array(np.meshgrid(lst1, lst2)).T.reshape(-1,2)


# function to read CSV file
# inputs:
#   input_file: name of input file
# outputs:
#   returns: 2D array of data
def read_file(input_file):
    data = []
    line_num = 1
    with open(input_file) as f:
        for line in f:
            if line_num != 1:
                line = line.strip()
                if line == "":
                    break
                data.append(line.split(","))
            line_num += 1
    return np.array(data)


# function to reformate and reshape target data
# inputs:
#   drugs: 2D array of all drugs
#   targets: original 2D array of all targets (one drug can span multiple rows)
# outputs:
#   returns: new 2D array of all targets (one drug can span one row only)
def reformate_targets(drugs, targets):
    temp_targets = []
    for drug in drugs:
        temp_targets.append(targets[np.argwhere(targets[:,0] == drug[0]).flatten()])
    return np.array(temp_targets)

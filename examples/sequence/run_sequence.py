#!/usr/bin/python
'''
Example: sequence module
Input: FASTA sequence with natural amino acids
Output: File with properties for a group of sequences and a library of peptides
'''

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa"]
__license__ = "MIT"
__version__ = "2.0"
__email__ = "rodrigo.ochoa@udea.edu.co, rodrigo.ochoa@boehringer-ingelheim.com"

########################################################################################
# Modules
########################################################################################

# PepFun modules
from pepfun.sequence import Sequence
from pepfun.sequence import Library

# List of sequences
sequences = ["KMGRLFR", "CAAAAAC", "PEPTIDES", "SDVAFRGNLLD", "ITFEDLLDYYG"]

report = open("properties_sequences.txt","w")
report.write("Sequence\tNet_charge\tHydrophobicity\tIsoelectric_point\tMolecular_weigth\n")

for sequence in sequences:
    pep = Sequence(sequence)
    print("The main peptide sequence is: {}".format(pep.sequence))

    # Sequence tables
    pep.compute_peptide_charges()
    print("Net charge at pH 7: {}".format(pep.netCharge))

    pep.calculate_properties_from_sequence()
    print("Average Hydrophobicity: {}".format(pep.avg_hydro))
    print("Isoelectric Point: {}".format(pep.isoelectric_point))

    # RDKit descriptors
    pep.calculate_properties_from_mol()
    print("Molecular weigth: {}".format(pep.mol_weight))
    print()
    report.write(f"{pep.sequence}\t{pep.netCharge}\t{pep.avg_hydro}\t{pep.isoelectric_point}\t{pep.mol_weight}\n")
report.close()

# Using a similarity matrix
pep = Sequence(sequences[0])
pep.align_position_matrix("KAGRSFR")
print("Alignment score by position with peptide KAGRSFR is: {}".format(pep.score_matrix))

# Counting matches
pep.align_position_local("KAGRSFR")
print("The number of dismatches with peptide KAGRSFR are: {}".format(pep.dismatch))

# Using the SMILES and fingerprint similarity
pep.similar_smiles("KAGRSFR")
print("The SMILES similarity is: {}".format(pep.smiles_similarity))
print()


# Create a library with all natural AAs and a size of 4 and frequency of 3 AAs per position
lib = Library(4)
peptides=lib.combinatorial_library(3)
library=open("library_sequences.txt", "w")
for p in peptides:
    library.write(p+"\n")
library.close()

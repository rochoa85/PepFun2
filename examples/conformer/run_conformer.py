#!/usr/bin/python
'''
Example: conformer module
Input: FASTA sequence with natural amino acids
Output: Conformers of the sequences depending on the method
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
from pepfun.conformer import Conformer

# List of sequences
sequences = ["PEPTIDES", "SDVAFRGNLLD"]

for sequence in sequences:
    # The auxiliar path folder can be changed depending on your system
    pep = Conformer(sequence, auxiliar_path="../../auxiliar")
    # Generate a PDB file with the conformer
    pep.generate_conformer()

# List of sequences
sequences = ["KMGRLFR", "ITFEDLLDYYG"]

for sequence in sequences:
    # The auxiliar path folder can be changed depending on your system
    pep = Conformer(sequence, auxiliar_path="../../auxiliar")
    # Predict the secondary structure
    ss=pep.run_psipred()
    # Generate the peptide templates, as well as the secondary structure restrictions
    pepT, pepM, rangeHelix, rangeBeta, rangeCycle = pep.prepare_modeller(ss)
    # Run Modeller
    pep.modelling_modeller(pepT, pepM, rangeHelix, rangeBeta, rangeCycle)
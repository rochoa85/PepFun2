#!/usr/bin/python
'''
Example: modifications module
Input: PDB files of the peptides
Output: Modified PDB files with filling amino acids, capping groups and NNAA mutations
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
from pepfun.modifications import Filling
from pepfun.modifications import Capping
from pepfun.modifications import Mutation

# Provide the template sequence and the sequence with the new amino acids
template="AKAFIAWLVRG"
model="LKAKAFIAWLVRGKL"
# The function accepts as input the peptide sequences and the code of the PDB file
code = "peptide_alone"
Filling.fill_peptide_alone(code,template,model)

# Template and sequence to be model
template="PTSYAGDDS"
model="AAPTSYALGDDSAA"
# Code of the complex and function to fill the remaining amino acids
code="peptide_complex"
Filling.fill_peptide_complex(code,template,model)

# Peptide template that should be the same in the PDB file
peptide="AKAFIAWLVRG"
# Calling the function with the pdb file and the mode to cap
pdb_file="peptide_alone.pdb"
Capping.capping_peptide(pdb_file,peptide,mode='N-term')

# Path to the input structure
pdb_file = "example_structure.pdb"
# Chain of the peptide in the structure
chainID = "C"
# Position that will be mutated
pepPosition = 2
# PDB of the NNAA that will be included
pdbNNAA = "Aib.pdb"
# Chain of the NNAA in the PDB file
chainNNAA = "A"
# 3-letter code of the NNAA
nnaa_name = "AIB"

# Call the function with the input variables
mut = Mutation(pdb_file, chainID, pepPosition, pdbNNAA, chainNNAA, nnaa_name)
# Get the backbone dihedrals from the reference structure
mut.get_bb_dihedral()
# Get dictionary with bond, angle and dihedral information
mut.get_NNAA_basic_info()
# Generate the PDB file with the mutation
mut.assign_mutation()
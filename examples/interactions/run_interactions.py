#!/usr/bin/python
'''
Example: interactions module
Input: PDB files of the protein-peptide complex
Output: Plot of the HBs, list of specific HBs, and summary report with observables per residue
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
from pepfun.interactions import Interactions

# PDB file of the complex
pdb_file = "example_structure.pdb"
# Chain of the peptide
chain = "C"
# Path to the auxiliar folder
aux_path = "../../auxiliar"
# Type of peptide conformation (linear or cyclic)
pep_conformation = "linear"
# Distance threshold to define contacts
contact_threshold = 4.0

pepStr = Interactions(pdb_file, chain, auxiliar_path=aux_path)
print("Peptide sequence based on the PDB file is: {}".format(pepStr.sequence))

pepStr.get_secondary_structure()
print("The predicted secondary structure is: {}".format(pepStr.total_dssp))

# Functions to detect the different types of interactions
pepStr.get_hydrogen_bonds()
pepStr.plot_hydrogen_bonds(pep_conformation)
pepStr.get_heavy_atom_contacts(contact_threshold)

print("The total number of contacts are: {}".format(pepStr.total_contacts))
print("The total number of hydrogen bonds are: {}".format(pepStr.number_hydrogen_bonds))

# Report with the summary per peptide position, including accesible surface area values
report = open("observables_residues.txt", "w")
header = "position"
for key in pepStr.positions[1]:
    header+=f" {key}"
report.write(header+"\n")

for key in pepStr.positions:
    line=str(key)
    for fields in pepStr.positions[key]:
        line+=f" {str(pepStr.positions[key][fields])}"
    report.write(line+"\n")
report.close()
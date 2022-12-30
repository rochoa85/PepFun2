"""
PepFun 2.0: improved protocols for the analysis of natural and modified peptides

From publication "PepFun 2.0: improved protocols for the analysis of natural and modified peptides"
Author: Rodrigo Ochoa
Year: 2023
"""

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

import unittest
import os
from pepfun.interactions import Interactions

##########################################################################
# Functions and classes
##########################################################################

def check_interactions(pdb):
    """
    Check secondary structure and interaction for a peptide in a complex

    :param pdb: Name of the PDB file
    :return: secondary structure, number of contacts, number of hydrogen bonds
    """

    try:
        pdb_file = "structures/interactions/{}".format(pdb)
        chain = "B"
        aux_path = "../auxiliar"
        contact_threshold = 4.0

        # Create the object
        pepStr = Interactions(pdb_file, chain, auxiliar_path=aux_path)
    except:
        pdb_file = "tests/structures/interactions/{}".format(pdb)
        chain = "B"
        aux_path = "auxiliar"
        contact_threshold = 4.0

        # Create the object
        pepStr = Interactions(pdb_file, chain, auxiliar_path=aux_path)

    # Get secondary structure
    pepStr.get_secondary_structure()
    ss=pepStr.total_dssp

    # Interactions
    pepStr.get_hydrogen_bonds()
    pepStr.get_heavy_atom_contacts(contact_threshold)
    contacts=pepStr.total_contacts
    hbs=pepStr.number_hydrogen_bonds

    hb_file="predicted_hbs_{}.txt".format(pepStr.sequence)
    os.unlink(hb_file)

    return ss,contacts,hbs

########################################################################################

class TestInteractions(unittest.TestCase):

    def test_interactions(self):
        """
        Check interactions for a set of PDB files with protein-peptide complexes
        """

        self.assertEqual(check_interactions('1xn2_complex.pdb'),('-----B-B--', 167, 23), )
        self.assertEqual(check_interactions('2axi_complex.pdb'),('-EEETTEEE-', 77, 5), )
        self.assertEqual(check_interactions('2qlf_complex.pdb'),('----', 35, 7), )
        self.assertEqual(check_interactions('2r5b_complex.pdb'),('--GGG-GGGHHHHH-', 50, 1), )
        self.assertEqual(check_interactions('2w6t_complex.pdb'),('----TTSS--', 115, 9), )

########################################################################################

if __name__ == "__main__":
    unittest.main()
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
from pepfun.conformer import Conformer

##########################################################################
# Functions and classes
##########################################################################

def generate_conformer(sequence):
    """
    Generate conformers from the sequence

    :param sequence: Peptide sequence
    :return: First line of the PDB files
    """

    pep = Conformer(sequence, auxiliar_path="../auxiliar")
    pep.generate_conformer()
    outFile = 'structure_{}.pdb'.format(sequence)

    with open(outFile, 'r') as f:
        data = f.read()
    os.unlink(outFile)

    line = data.split('\n')[4]

    return line

########################################################################################

class TestConformer(unittest.TestCase):

    def test_conformer(self):
        """
        Generate basic conformers to test RDKit functionalities
        """

        self.assertEqual(generate_conformer('GAANDENY'),'ATOM      1  N   GLY A   1      15.502  -0.915  -1.139  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('EAPPSYAEV'),'ATOM      1  N   GLU A   1      12.874  -0.568  -1.208  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('SDVAFRGNLLD'),'ATOM      1  N   SER A   1     -11.941   4.662  -1.609  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('GVLKEYGV'),'ATOM      1  N   GLY A   1     -11.575   3.737   2.619  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('MCLRMTAVM'),'ATOM      1  N   MET A   1      13.534   2.297   0.936  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('EEFELLISNS'),'ATOM      1  N   GLU A   1     -11.564   3.543  -0.795  1.00  0.00           N  ')
        self.assertEqual(generate_conformer('SQFDLSTRRLK'),'ATOM      1  N   SER A   1     -13.906   5.066   2.200  1.00  0.00           N  ')


if __name__ == "__main__":
    unittest.main()
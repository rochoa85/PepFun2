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
from pepfun.modifications import Filling
from pepfun.modifications import Capping
from pepfun.modifications import Mutation

##########################################################################
# Functions and classes
##########################################################################

def run_filling(template, model, pdb, mode):
    """
    Fill a peptide based on a PDB file

    :param template: Peptide sequence as template
    :param model: Peptide sequence with new amino acids to fill
    :param pdb: Name of the PDB file
    :return: First line of the PDB file with the filled peptide
    """

    try:
        path='structures/modifications/'
        if mode=="alone":
            Filling.fill_peptide_alone(pdb, template, model, path=path)
        elif mode=="complex":
            Filling.fill_peptide_complex(pdb, template, model, path=path)
    except:
        path = 'tests/structures/modifications/'
        if mode == "alone":
            Filling.fill_peptide_alone(pdb, template, model, path=path)
        elif mode == "complex":
            Filling.fill_peptide_complex(pdb, template, model, path=path)

    outFile = 'filled_{}.pdb'.format(pdb)

    with open(outFile, 'r') as f:
        data = f.read()
    os.unlink(outFile)

    line = data.split('\n')[11]

    return line

##########################################################################

def run_capping(sequence, pdb, mode):
    """
    Fill a peptide based on a PDB file

    :param sequence: Peptide sequence as template
    :param pdb: Name of the PDB file
    :param mode: Mode to cap both flanking or only one
    :return: First line of the PDB file with the filled peptide
    """

    try:
        path = 'structures/modifications/'
        Capping.capping_peptide(pdb, sequence, path=path, mode=mode)
    except:
        path = 'tests/structures/modifications/'
        Capping.capping_peptide(pdb, sequence, path=path, mode=mode)

    outFile = 'capped_{}'.format(pdb)

    with open(outFile, 'r') as f:
        data = f.read()
    os.unlink(outFile)

    line = data.split('\n')[0]

    return line

########################################################################################

def run_mutation(pdb_ref, pdb_nnaa, chain_ref, chain_nnaa, nnaa_name):
    """
    Mutate a NNAA in the position of a natural

    :param pdb_ref: PDB file of the peptide
    :param pdb_nnaa: PDB of the NNAA
    :param chain_ref: Chain of the peptide
    :param chain_nnaa: Chain of the NNAA
    :param nnaa_name: Name of the NNAA pdb code
    :return: Line of the mutated PDB file
    """

    try:
        pdb_file = "structures/modifications/{}".format(pdb_ref)
        pdbNNAA = "structures/modifications/{}".format(pdb_nnaa)
        pepPosition = 3

        # Create the object
        mut = Mutation(pdb_file, chain_ref, pepPosition, pdbNNAA, chain_nnaa, nnaa_name)
        mut.get_bb_dihedral()
        mut.get_NNAA_basic_info()
        mut.assign_mutation()
    except:
        pdb_file = "tests/structures/modifications/{}".format(pdb_ref)
        pdbNNAA = "tests/structures/modifications/{}".format(pdb_nnaa)
        pepPosition = 3

        # Create the object
        mut = Mutation(pdb_file, chain_ref, pepPosition, pdbNNAA, chain_nnaa, nnaa_name)
        mut.get_bb_dihedral()
        mut.get_NNAA_basic_info()
        mut.assign_mutation()

    outFile = 'mutated_{}.pdb'.format(nnaa_name)

    with open(outFile, 'r') as f:
        data = f.read()
    os.unlink(outFile)

    line = data.split('\n')[0]

    return line

########################################################################################

class TestModifications(unittest.TestCase):

    def test_filling(self):
        """
        Check filling approches
        """

        self.assertEqual(run_filling("AKAFIAWLVRG", "LKAKAFIAWLVRGKL", "peptide_1", mode="alone"),'ATOM      5  CD1 LEU A   1       8.408  23.959  -3.290  1.00 39.83           C')
        self.assertEqual(run_filling("PTSYAGDDS", "AAPTSYALGDDSAA", "peptide_complex_1", mode="complex"),'ATOM      1  N   ILE A   1      12.408 -35.108   3.957  1.00155.37           N')
        self.assertEqual(run_filling("YSNTLPVRK", "TYSNTIILPVRKS", "peptide_2", mode="alone"),'ATOM      5  CG2 THR A   1      32.630  -0.556  39.151  1.00 31.64           C')
        self.assertEqual(run_filling("NPVVHFFKNIVTPRTPPPSQ", "GNPVVHFFKNIVTPRTLLPPPSQG", "peptide_complex_2", mode="complex"),'ATOM      3  CB  GLU A   1      19.686   7.645 -14.427  1.00 34.44           C')

    ########################################################################################

    def test_capping(self):
        """
        Check capping approch for a single peptide
        """

        self.assertEqual(run_capping("YSNTLPVRK","peptide_2.pdb", mode="both"),'ATOM      1  CH3 ACE A   1      46.635  24.343  55.770  1.00  6.93           C  ')
        self.assertEqual(run_capping("AKAFIAWLVRG","peptide_1.pdb", mode="N-term"),'ATOM      1  CH3 ACE A   1      10.193  13.568  12.078  1.00 28.32           C  ')

    ########################################################################################

    def test_mutation(self):
        """
        Check mutation function
        """

        self.assertEqual(run_mutation("peptide_1.pdb", "Aib.pdb", "A", "A", "AIB"),'ATOM      1  N   ALA A   1      -1.777   4.394   5.174  1.00  0.00           N  ')
        self.assertEqual(run_mutation("peptide_2.pdb", "Nle.pdb", "A", "A", "NLE"),'ATOM      1  N   TYR A   1      27.561   4.601  48.452  1.00 65.91           N  ')
        self.assertEqual(run_mutation("peptide_complex_1.pdb", "Aib.pdb", "B", "A", "AIB"),'ATOM      1  N   ILE A  16      14.152 -32.013   4.101  1.00 11.52           N  ')
        self.assertEqual(run_mutation("peptide_complex_2.pdb", "Nle.pdb", "B", "A", "NLE"),'ATOM      1  N   GLU A   1      20.364   9.466 -16.009  1.00  0.00           N  ')

########################################################################################

if __name__ == "__main__":
    unittest.main()
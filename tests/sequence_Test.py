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
from pepfun.sequence import Sequence
from pepfun.sequence import Library

##########################################################################
# Functions and classes
##########################################################################

def get_properties(sequence):
    """
    Calculate properties from the sequence

    :param sequence: Peptide sequence
    :return: net charge, average hydrophobicity, isoelectric point and molecular weigth
    """
    # Create object
    pep = Sequence(sequence)

    # Properties from the sequence
    pep.compute_peptide_charges()
    netCharge=pep.netCharge
    pep.calculate_properties_from_sequence()
    avgHydro=pep.avg_hydro
    isoPoint=pep.isoelectric_point

    # RDKit descriptors
    pep.calculate_properties_from_mol()
    molWeight=pep.mol_weight

    return netCharge,avgHydro,isoPoint,molWeight

def get_rules(sequence):
    """
    Calculate similarity and some empirical rules

    :param sequence: Peptide sequence
    :return: similarity, failed solubility rules and failed synthesis rules
    """
    # Create object
    pep = Sequence(sequence)

    # Similarity
    pep.similar_smiles("KAGRSFR")
    sim=pep.smiles_similarity

    # Empirical rules
    pep.solubility_rules()
    sol=pep.solubility_rules_failed
    pep.synthesis_rules()
    syn=pep.synthesis_rules_failed


    return sim,sol,syn

##########################################################################

def library_functions(size,list_aa):
    """
    Test some of the Library functions

    :param size: size of the peptide
    :param list_aa: list of amino acids included in the library
    :return: first 5 elements of the list

    """
    lib = Library(size, list_aa=list_aa)
    frags = lib.generate_sequences()
    res_frags=frags[0:5]

    return res_frags

##########################################################################

class TestSequence(unittest.TestCase):

    def test_sequenceProperties(self):
        """
        Test the properties predicted for some peptides
        """
        # list_peps=['GAANDENY',
        # 'EAPPSYAEV',
        # 'SDVAFRGNLLD',
        # 'RRNLKGLNLNLH',
        # 'GVLKEYGV',
        # 'MCLRMTAVM',
        # 'EEFELLISNS',
        # 'SQFDLSTRRLK',
        # 'KLMFKTEGPDSD']
        #
        # for p in list_peps:
        #     a,b,c,d=get_properties(p)
        #     print(f"{a},{b},{c},{d}")

        self.assertEqual(get_properties('GAANDENY'), (-2.0008411239265196, -1.22, 4.0500284194946286, 852.8120000000001), )
        self.assertEqual(get_properties('EAPPSYAEV'), (-1.999412677161022, 1.1600000000000001, 4.0500284194946286, 962.0240000000001), )
        self.assertEqual(get_properties('SDVAFRGNLLD'), (-1.0055728771384773, 0.2000000000000003, 4.533500862121582, 1206.3229999999994), )
        self.assertEqual(get_properties('RRNLKGLNLNLH'), (3.0896133167191, -4.58, 11.999967765808105, 1447.7139999999997), )
        self.assertEqual(get_properties('GVLKEYGV'), (-0.001722146326308235, 2.2, 6.000685691833495, 864.0110000000002), )
        self.assertEqual(get_properties('MCLRMTAVM'), (0.9492040695067787, 2.39, 7.999833488464356, 1055.427), )
        self.assertEqual(get_properties('EEFELLISNS'), (-2.99679181488375, 1.3299999999999996, 4.0500284194946286, 1180.277), )
        self.assertEqual(get_properties('SQFDLSTRRLK'), (1.9935405781644697, -5.41, 10.834379386901855, 1350.545), )
        self.assertEqual(get_properties('KLMFKTEGPDSD'), (-1.0081847499417282, -2.28, 4.780807304382324, 1367.545), )

    ##########################################################################

    def test_sequenceRules(self):
        """
        Test similarity and some empirical rules for peptides
        """

        self.assertEqual(get_rules('GAANDENY'), (0.3595505617977528, 0, 0), )
        self.assertEqual(get_rules('EAPPSYAEV'), (0.3425925925925926, 2, 2), )
        self.assertEqual(get_rules('SDVAFRGNLLD'), (0.6153846153846154, 2, 0), )
        self.assertEqual(get_rules('RRNLKGLNLNLH'), (0.5434782608695652, 4, 2), )
        self.assertEqual(get_rules('GVLKEYGV'), (0.41935483870967744, 2, 0), )
        self.assertEqual(get_rules('MCLRMTAVM'), (0.40625, 2, 2), )
        self.assertEqual(get_rules('EEFELLISNS'), (0.42391304347826086, 3, 2), )
        self.assertEqual(get_rules('SQFDLSTRRLK'), (0.6043956043956044, 3, 0), )
        self.assertEqual(get_rules('KLMFKTEGPDSD'), (0.4260869565217391, 3, 1), )

    ##########################################################################

    def test_library(self):
        """
        Test some library functions
        """

        self.assertEqual(library_functions(3,['A','L','S']), (['AAA', 'AAL', 'AAS', 'ALA', 'ALL']), )
        self.assertEqual(library_functions(4,['M','T','C','L','L']), (['MMMM', 'MMMT', 'MMMC', 'MMML', 'MMML']), )
        self.assertEqual(library_functions(6,['K','P','R']), (['KKKKKK', 'KKKKKP', 'KKKKKR', 'KKKKPK', 'KKKKPP']), )
        self.assertEqual(library_functions(5,[]), (['AAAAA', 'AAAAC', 'AAAAD', 'AAAAE', 'AAAAF']), )

##########################################################################

if __name__ == "__main__":
    unittest.main()

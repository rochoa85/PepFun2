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
from pepfun.extra import Alignment
from pepfun.extra import pepDescriptors

##########################################################################
# Functions and classes
##########################################################################

def run_similarity(seq1, seq2):
    """
    Run similarity between two peptides

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :return: similarity value
    """

    similarity = Alignment.similarity_pair(seq1, seq2)

    return  similarity

##########################################################################

def run_score(seq1, seq2, mode="unweighted"):
    """
    Run alignment between two peptides (BILN format)

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :return: score
    """

    score, start, end = Alignment.score_modified(seq1, seq2, mode=mode)

    return score

##########################################################################

def generate_descriptors(sequence):
    """
    Generate descriptors for a peptide (BILN format)

    :param sequence: Peptide sequence
    :return: First 5 descriptors
    """

    # Create the dictionary
    desc = pepDescriptors(sequence, generate_properties=False)
    descriptors = desc.moranCorrelation()

    desc_pairs=[]
    counter=1
    for key, val in descriptors.items():
        if counter > 5:
            break
        desc_pairs.append((key,val))
        counter+=1

    return desc_pairs

########################################################################################

class TestExtra(unittest.TestCase):

    def test_similarity(self):
        """
        Run similarity analysis between 2 sequences (equal length)
        """

        self.assertEqual(run_similarity('AFTGYW','AGTGYL'),0.4079015827236455)
        self.assertEqual(run_similarity('NPVVHFFKNIVTPRTPPPSQ', 'AAAAAFFKNIVAAAAAAAAA'),0.22966770070528583)
        self.assertEqual(run_similarity('LLSHYTSY', 'LLSHYTSY'),1.0)

    ########################################################################################

    def test_score(self):
        """
        Run alignment analysis between 2 sequences of different lenght and having NNAAs too
        """

        self.assertEqual(run_score("W-W-S-E-V-N-10L-A-E-F", "K-T-E-E-I-S-E-V-N-STA-V-A-E-F"), 7.0)
        self.assertEqual(run_score("K-Aib-M-P","S-Aib-Iva-P", mode="weighted"), 8.0)
        self.assertEqual(run_score("L-L-S-H-Y-T-S-Y", "A-F-T-G-Y-W"), 2.0)

    ########################################################################################

    def test_descriptors(self):
        """
        Generate descriptors for a modified peptide in BILN format. Monomers should be part of the dictionary
        """

        self.assertEqual(generate_descriptors("K-Aib-M-P-C-A"),[('nrot-lag1', -0.5565217391304348), ('nrot-lag2', 0.7282608695652173), ('nrot-lag3', -0.6956521739130433),('nrot-lag4', 0.1739130434782609), ('nrot-lag5', -1.391304347826087)])
        self.assertEqual(generate_descriptors("S-Aib-A-L-L-R-E-Iva-P"),[('nrot-lag1', 0.5256849315068494), ('nrot-lag2', -0.1634050880626223), ('nrot-lag3', -0.7226027397260273),('nrot-lag4', -0.8027397260273974), ('nrot-lag5', -0.49143835616438347)])
        self.assertEqual(generate_descriptors("A-A-L-Nva-C-W-E-Orn"),[('nrot-lag1', 0.47558770343580475), ('nrot-lag2', -0.22362869198312232), ('nrot-lag3', 0.03291139240506332),('nrot-lag4', 0.2911392405063291), ('nrot-lag5', -0.628691983122363)])


########################################################################################

if __name__ == "__main__":
    unittest.main()
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
# Modules to import
########################################################################################

# External modules
import re
import pickle
from igraph import *
from random import randint
import itertools

# BioPython
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors


########################################################################################
# Classes and Functions
########################################################################################

def aminoacidSMILES(amino):
    """
    Obtain the SMILES representation of a particular amino acid

    :param amino: One-letter code of the amino acid

    :return smiles: 1D Chemical representation of the amino acid
    """

    # Dictionary with the SMILES per amino acid
    aminoacids = {'G': {'SMILES': 'NCC(=O)O'},
                  'A': {'SMILES': 'N[C@@]([H])(C)C(=O)O'},
                  'R': {'SMILES': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O'},
                  'N': {'SMILES': 'N[C@@]([H])(CC(=O)N)C(=O)O'},
                  'D': {'SMILES': 'N[C@@]([H])(CC(=O)O)C(=O)O'},
                  'C': {'SMILES': 'N[C@@]([H])(CS)C(=O)O'},
                  'E': {'SMILES': 'N[C@@]([H])(CCC(=O)O)C(=O)O'},
                  'Q': {'SMILES': 'N[C@@]([H])(CCC(=O)N)C(=O)O'},
                  'H': {'SMILES': 'N[C@@]([H])(CC1=CN=C-N1)C(=O)O'},
                  'I': {'SMILES': 'N[C@@]([H])(C(CC)C)C(=O)O'},
                  'L': {'SMILES': 'N[C@@]([H])(CC(C)C)C(=O)O'},
                  'K': {'SMILES': 'N[C@@]([H])(CCCCN)C(=O)O'},
                  'M': {'SMILES': 'N[C@@]([H])(CCSC)C(=O)O'},
                  'F': {'SMILES': 'N[C@@]([H])(Cc1ccccc1)C(=O)O'},
                  'P': {'SMILES': 'N1[C@@]([H])(CCC1)C(=O)O'},
                  'S': {'SMILES': 'N[C@@]([H])(CO)C(=O)O'},
                  'T': {'SMILES': 'N[C@@]([H])(C(O)C)C(=O)O'},
                  'W': {'SMILES': 'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O'},
                  'Y': {'SMILES': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O'},
                  'V': {'SMILES': 'N[C@@]([H])(C(C)C)C(=O)O'}}

    # Store the SMILES in a variable
    smiles = aminoacids[amino]['SMILES']
    return smiles

########################################################################################

class Sequence:
    """
    Class with functions to perform different type of analysis using a peptide sequence as an object
    """

    def __init__(self, sequence):
        """
        Inititalize the class calculating some basic properties

        :param sequence: Peptide sequence

        :return Based on the sequence, the class start with some counters, a neutral pH, the peptide length and the SMILES representation
        """
        self.sequence = sequence
        self.pH = 7
        self.solubility_rules_failed = 0
        self.set_sol_rules = []
        self.length_peptide = len(self.sequence)
        self.synthesis_rules_failed = 0
        self.set_syn_rules = []

        # Loop to create the SMILES
        connect_smiles = 'O'
        for res in sequence:
            connect_smiles = connect_smiles[:-1]
            smiles = aminoacidSMILES(res)
            connect_smiles = connect_smiles + smiles
        self.smiles = connect_smiles

    ############################################################################

    def align_position_matrix(self, peptide_to_match, matrixInput='matrix.pickle'):
        """
        Align position by position a peptide of reference with another one

        :param peptide_to_match: String with the sequence of a peptide that will be compared
        :param matrixInput: A dictionary of biopython with substitution scores.

        :return score: A numerical value to rank the alignment
        """

        try:
            import pkg_resources
            stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(matrixInput))
            matrix = pickle.load(stream)
        except:
            with open('data/{}'.format(matrixInput), 'rb') as handle:
                matrix = pickle.load(handle)

        self.score_matrix = 0

        for i, p in enumerate(self.sequence):
            # Generate the tuples with the pair of amino acids to compare
            pair1 = (p, peptide_to_match[i])
            pair2 = (peptide_to_match[i], p)

            # Obtain the score from the matrix and sum up the values
            if pair1 in matrix:
                value = matrix[pair1]
            else:
                value = matrix[pair2]
            self.score_matrix += value

    ############################################################################
    def align_position_local(self, peptide_to_match):
        """
        Align position by position a peptide of reference with another one

        :param peptide_to_match: String with the sequence of a peptide that will be compared

        :return match: Number of matches
        """
        self.dismatch = 0

        for i, p in enumerate(self.sequence):
            # Generate the tuples with the pair of amino acids to compare
            if p != peptide_to_match[i]:
                self.dismatch += 1

        self.match = len(self.sequence) - self.dismatch

    ############################################################################
    def similar_smiles(self, peptide_to_match):
        """
        Calculate similarity but using SMILES representations of the peptides

        :param peptide_to_match: peptide sequence that will be compared

        :return smiles_similarity: based on Morgan Fingerprints and Tanimoto coefficient
        """

        # Generate molecule from sequence
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol1.SetProp("_Name", self.sequence)

        connect_smiles = 'O'
        for res in peptide_to_match:
            connect_smiles = connect_smiles[:-1]
            smiles = aminoacidSMILES(res)
            connect_smiles = connect_smiles + smiles

        mol2 = Chem.MolFromSmiles(connect_smiles)
        mol2.SetProp("_Name", peptide_to_match)

        # Calculate the fingerprints and the similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, 2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, 2048)

        self.smiles_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

    ############################################################################
    def compute_peptide_charges(self, pH_internal=7):
        """
        Function to calculate the average net charge based on pka values

        :param pH_internal -- By default is 7

        :return net_charge: The net charge based on reported pka values
        """
        # Set the general variables and pka terms
        self.pH = pH_internal
        self.netCharge = 0.0
        pka_alpha_amino = {'G': 9.60, 'A': 9.69, 'V': 9.62, 'L': 9.60, 'I': 9.68, 'M': 9.21, 'F': 9.13, 'W': 9.39,
                           'P': 10.60, 'S': 9.15,
                           'T': 9.10, 'C': 10.78, 'Y': 9.11, 'N': 8.84, 'Q': 9.13, 'D': 9.82, 'E': 9.67, 'K': 8.95,
                           'R': 9.04, 'H': 9.17}
        pka_alpha_carboxy = {'G': 2.34, 'A': 2.34, 'V': 2.32, 'L': 2.36, 'I': 2.36, 'M': 2.28, 'F': 1.83, 'W': 2.38,
                             'P': 1.99, 'S': 2.21,
                             'T': 2.63, 'C': 1.71, 'Y': 2.2, 'N': 2.02, 'Q': 2.17, 'D': 2.09, 'E': 2.19, 'K': 2.18,
                             'R': 2.17, 'H': 1.82}
        pka_sidechain_positive = {'K': 10.79, 'R': 12.48, 'H': 6.04}
        pka_sidechain_negative = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10.07}

        # Calculate the net charge for the extreme groups (without modifications)
        amino = self.sequence[0]
        carboxy = self.sequence[-1]
        self.netCharge = self.netCharge + (math.pow(10, pka_alpha_amino[amino]) / (
                    math.pow(10, pka_alpha_amino[amino]) + math.pow(10, self.pH)))
        self.netCharge = self.netCharge - (
                    math.pow(10, self.pH) / (math.pow(10, pka_alpha_carboxy[carboxy]) + math.pow(10, self.pH)))

        # Calculate the net charge for the charged amino acid side chains
        for aa in self.sequence:
            if aa in pka_sidechain_positive:
                self.netCharge = self.netCharge + (math.pow(10, pka_sidechain_positive[aa]) / (
                            math.pow(10, pka_sidechain_positive[aa]) + math.pow(10, self.pH)))
            if aa in pka_sidechain_negative:
                self.netCharge = self.netCharge - (
                            math.pow(10, self.pH) / (math.pow(10, pka_sidechain_negative[aa]) + math.pow(10, self.pH)))

    ############################################################################
    def calculate_properties_from_mol(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities

        :return: Static physico-chemical properties: molecular weight, crippen logP, number of hydrogen bond acceptors and donors
        """

        # Generate molecule from sequence
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name", self.sequence)

        # Calculate the descriptors
        self.num_hdonors = Lipinski.NumHDonors(mol)
        self.num_hacceptors = Lipinski.NumHAcceptors(mol)
        self.mol_weight = Descriptors.MolWt(mol)
        self.mol_logp = Crippen.MolLogP(mol)

    ############################################################################
    def calculate_properties_from_sequence(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities

        :return Average Eisenberg hydrophobicity, ProtParam parameters: Isolectric point, aromaticity, instability index, amino acid percentage
        """

        # Hydrophobicity -> Eisenberg scale
        hydrophobicity = {'A': 0.620, 'R': -2.530, 'N': -0.780, 'D': -0.900, 'C': 0.290, 'Q': -0.850, 'E': -0.740,
                          'G': 0.480, 'H': -0.400, 'Y': 0.260,
                          'I': 1.380, 'L': 1.060, 'K': -1.500, 'M': 0.640, 'F': 1.190, 'P': 0.120, 'S': -0.180,
                          'T': -0.050, 'W': 0.810, 'V': 1.080}
        self.avg_hydro = sum([hydrophobicity[resi] for resi in self.sequence])

        # ProParam properties
        prot_parameters = ProteinAnalysis(self.sequence)
        self.aromaticity = prot_parameters.aromaticity()
        self.aa_percent = prot_parameters.get_amino_acids_percent()
        self.instability_index = prot_parameters.instability_index()
        self.isoelectric_point = prot_parameters.isoelectric_point()

    ############################################################################
    def solubility_rules(self):
        """
        Function to calculate some solubility rules based on recommendations of http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/

        :return solubility_rules_failed: return the number of rules failed based on the criteria
        """
        # Rule N1. Number of hydrophobic or charged residues
        hydro_residues = ['V', 'I', 'L', 'M', 'F', 'W', 'C']
        charged_residues = ['H', 'R', 'K', 'D', 'E']

        count_hydro_charged = 0
        for aa in self.sequence:
            if aa in hydro_residues or aa in charged_residues: count_hydro_charged += 1

        # This condition should change depending on the sequence length
        hydro_char_threshold = float(self.length_peptide) * 0.45
        if count_hydro_charged > hydro_char_threshold:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N2. Computed peptide charge
        charge_threshold = 1
        self.compute_peptide_charges()
        if self.netCharge > 1:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N3. Glycine or Proline content in the sequence
        count_gly_pro = 0
        for aa in self.sequence:
            if aa == "G" or aa == "P": count_gly_pro += 1
        # Check threshold
        if count_gly_pro > 1:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N4. First or last amino acid charged
        count_charge = 0
        if self.sequence[0] in charged_residues:
            count_charge += 1
        if self.sequence[-1] in charged_residues:
            count_charge += 1
        # Check threshold
        if count_charge > 0:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N5. Any amino acid represent more than 25% of the total sequence
        prot_parameters = ProteinAnalysis(self.sequence)
        aa_content = prot_parameters.get_amino_acids_percent()
        flag5 = 0
        for aa in aa_content:
            if aa_content[aa] >= 0.3:
                self.solubility_rules_failed += 1
                self.set_sol_rules.append(1)
                flag5 = 1
                break
        if flag5 == 0: self.set_sol_rules.append(0)

    ############################################################################
    def synthesis_rules(self):
        """
        Function to check some synthesis rules based on empirical recommendations

        :return synthesis_rules_failed: return the number of rules failed based on the criteria
        """
        # Presence of forbiden motifs
        forbidden_motifs = {'2-prolines': r'[P]{3,}', 'DG-DP': r'D[GP]', 'N-Q-Nterminal': r'^[NQ]', }
        for motif in forbidden_motifs:
            if re.search(forbidden_motifs[motif], self.sequence):
                self.synthesis_rules_failed += 1
                self.set_syn_rules.append(1)
            else:
                self.set_syn_rules.append(0)

        # test if there are charged residues every 5 amino acids
        charged_residues = ['H', 'R', 'K', 'D', 'E']
        counter_charged = 0
        for residue in self.sequence:
            counter_charged += 1
            if residue in charged_residues:
                counter_charged = 0
            if counter_charged >= 5:
                self.synthesis_rules_failed += 1
                self.set_syn_rules.append(1)
            else:
                self.set_syn_rules.append(0)

        # Check if there are oxidation-sensitive amino acids
        aa_oxidation = ['M', 'C', 'W']
        flag5 = 0
        for aa in self.sequence:
            if aa in aa_oxidation:
                self.synthesis_rules_failed += 1
                self.set_syn_rules.append(1)
                flag5 = 1
                break
        if flag5 == 0: self.set_syn_rules.append(0)

    ############################################################################
    def blast_online(self):
        """
        Function to run online blast configured with parameters suitable to compare peptides

        :return hits: List of hits with dictionary containing fields from the alignment result
        """
        # Create a temporal fasta file with the sequence
        fasta_file = open("{}.fasta".format(self.sequence), "w")
        fasta_file.write(">sequence\n{}".format(self.sequence))
        fasta_file.close()

        record = SeqIO.read("{}.fasta".format(self.sequence), format="fasta")
        result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), word_size=2, expect=20000.0,
                                       matrix_name="PAM30", gapcosts="9 1", format_object="Alignment")
        b_record = NCBIXML.read(result_handle)

        # Parse the results
        hits = []
        for alignment in b_record.alignments:
            for hsp in alignment.hsps:
                dict_hits = {}
                dict_hits["identities"] = hsp.identities
                dict_hits["positives"] = hsp.positives
                dict_hits["gaps"] = hsp.gaps
                dict_hits["align_length"] = hsp.align_length
                dict_hits["query_start"] = hsp.query_start
                dict_hits["e-value"] = hsp.expect
                dict_hits["query_sequence"] = hsp.query[0:75]
                dict_hits["match_id"] = alignment.title[:100]
                dict_hits["subject_sequence"] = hsp.sbjct[0:75]
                hits.append(dict_hits)

        return hits

    # End of Sequence class
    ############################################################

########################################################################################

class Library:

    def __init__(self, len_peptide, list_aa=[]):
        """
        Initialize with the len of the fragments and the list of AAs that will be used
        By default the 20 natural AAs are used

        :param len_peptide: Length of the peptides that are going to be generated
        :param list_aa: List of AAs that will be included in the library
        """
        # List of the amino acids that can be included
        self.list_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        if list_aa:
            self.list_aa=list_aa

        # Length of the fragments that will be generated
        self.len_peptide=len_peptide

    ########################################################################################
    def generate_sequences(self):
        """
        Method to generate a combinatorial library of peptide sequences using natural amino acids

        :return frags: List with the sequences generated
        """

        # List where the sequences will be stored
        frags = []
        # Generate the combinatorial library
        for i in itertools.product(self.list_aa, repeat=self.len_peptide):
            frags.append(''.join(map(str, i)))

        # Return the list of peptide sequences
        return frags

    ########################################################################################
    def generate_peptide_pattern(self, pattern):
        """
        Method to generate a peptide sequences following a pattern using only natural amino acids

        :param pattern: Pattern based on XX to generate the list of sequences

        :return list_peptides: List with the sequences generated
        """

        list_peptides = []

        # Combination of possible sequences following the given pattern
        for pep in itertools.product(self.list_aa, repeat=pattern.count('X')):
            pep = list(pep)
            mut_pep = []
            for res in pattern:
                # check if any amino acid can be included in the sequence
                if res != 'X':
                    mut_pep.append(res)
                else:
                    mut_pep.append(pep.pop(0))
            list_peptides.append("".join(mut_pep))

        # Return the list of peptide sequences
        return list_peptides

    ########################################################################################
    def combinatorial_library(self, frequency):
        """
        Method to generate a combinatorial library of peptide sequences but limiting the frequency of amino acids in certain positions

        :param frequency: maximum number of repetitions per amino acid

        :return peptides: List with the sequences generated
        """

        # List where the sequences will be stored
        peptides = []

        # Counter of the peptide length
        for i in range(0, self.len_peptide):
            for aa in self.list_aa:
                # Run the process 3 times per amino acids
                for j in range(0, 3):
                    sequence = [None] * self.len_peptide
                    for k in range(0, self.len_peptide):
                        if k == i:
                            sequence[k] = aa
                        else:
                            # Loop to check that the sequence does not have more than 3 equal amino acids
                            accept_change = 0
                            while accept_change == 0:
                                posRand = randint(0, len(self.list_aa) - 1)
                                new_aa = self.list_aa[posRand]
                                if sequence.count(new_aa) <= frequency:
                                    accept_change = 1
                            # Add the amino acid after checking the restraint
                            sequence[k] = new_aa

                    # Join the amino acids and append the new sequence
                    seq_string = ''.join(sequence)
                    peptides.append(seq_string)

        # Return the list
        return peptides

    ########################################################################################
    def frequencies_library(self, peptides):
        """
        Calculate the frequency of each amino acid in a particular peptide library

        :param peptides: List of peptides in the library

        :return count_dict: Dictionary with the numbers per each natural amino acid
        """

        # Create counter dictionary
        count_dict = {}
        for i in range(1, len(peptides[0]) + 1):
            count_dict[i] = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0,
                             "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0}

        # Read the sequences and count each amino acid
        for sequence in peptides:
            for pos, aa in enumerate(sequence):
                count_dict[pos + 1][aa] += 1

        # Return the dictionary
        return count_dict

    # End of Library class
    ############################################################

############################################################
## End of sequence.py
############################################################
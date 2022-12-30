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
import math
import pandas as pd
import re
import os
import pickle
import numpy as np
from itertools import combinations

# BioPython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

# Internal
from pepfun.sequence import Sequence

########################################################################################
# Classes and functions
########################################################################################

class Alignment:

    """Class to allow the alignment of two sequences, including those with natural or unnatural amino acids"""

    @staticmethod
    def generate_matrix(threshold=60, monomers='monomers.csv'):
        """
        Class to generate a similarity matrix for a set of non-natural monomers

        :param threshold: Value to define a similarity threshold for the matrix
        :return: A pickle file with the generated matrix. This matrix is provided by default in the package
        """

        # Read the CSV file with the monomer information
        try:
            import pkg_resources
            stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(monomers))
            df = pd.read_csv(stream)
        except:
            with open('data/{}'.format(monomers), 'rb') as handle:
                df = pd.read_csv(handle)

        # Extract the smiles from the dataframe
        totalNames=[]
        totalSmiles=[]
        for i, idx in enumerate(df.index):
            name = df.at[idx, 'Symbol']
            smiles = df.at[idx, 'SMILES']
            totalNames.append(name)
            totalSmiles.append(smiles)

        # Run pair alignments and create the matrix
        matrixP = {}
        print("Generating the matrix ...")
        for i,name1 in enumerate(totalNames):
            smiles1=totalSmiles[i]
            totalValues=[]
            for j,name2 in enumerate(totalNames):
                smiles2=totalSmiles[j]
                try:
                     mol1= Chem.MolFromSmiles(smiles1)
                     mol2 = Chem.MolFromSmiles(smiles2)
                     fp1 = FingerprintMols.FingerprintMol(mol1)
                     fp2 = FingerprintMols.FingerprintMol(mol2)
                     smiles_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                     totalValues.append(smiles_similarity)
                except:
                    totalValues.append(0.0)

            for i,ele in enumerate(totalValues):
                name2=totalNames[i]
                if (name1,name2) not in matrixP and (name2,name1) not in matrixP:
                    matrixP[(name1,name2)]=10*ele-(threshold/10)

        # Output file in a local data folder
        print("The matrix with a threshold of {}% has been generated".format(threshold))
        if not os.path.exists("data"):
            os.makedirs("data")
        pickle.dump(matrixP, open("data/simMatrix_{}.pickle".format(threshold), "wb"))

    ########################################################################################
    @staticmethod
    def __convert_list(peptide):
        """
        Function to convert a peptide separated by dashes in a list

        :param peptide: Peptide sequence in BILN format (separated by dashes)
        :return: List with the peptide monomers
        """
        peptideList = []
        monomers = peptide.split('-')
        for res in monomers:
            mon = re.sub("\(\d+,\d+\)", "", res)
            peptideList.append(str(mon))
        return peptideList

    ########################################################################################
    @staticmethod
    def score_modified(peptide1, peptide2, mode='unweighted', matrix='simMatrix_60.pickle'):
        """
        Method to compare two sequences with non-natural amino acids. The input sequences should have the monomers separated by dashes.

        :param peptide1: Peptide input 1 in BILN format (monomers separated by dashes)
        :param peptide2: Peptide input 2 in BILN format (monomers separated by dashes)
        :param mode: Use the matrix or not in the alignment. Options between unweighted (default) or weighted
        :param matrix: Name of the matrix file that will be used to run the alignment

        :return: score of the alignment
        """

        # Check if the peptides are in BILN format. This can be checked with the pyPept package.
        if "-" not in peptide1 or "-" not in peptide2:
            print("The monomers should be separated by dash symbols")
            exit(1)
        if mode not in ['unweighted','weighted']:
            print("Please select one mode between weighted or unweighted (default)")
            exit(1)

        # Convert the sequences to list
        peptideList1 = Alignment.__convert_list(peptide1)
        peptideList2 = Alignment.__convert_list(peptide2)

        # Without matrix
        if mode == "unweighted":
            print("Only matching alignment")
            for a in pairwise2.align.globalxx(peptideList1, peptideList2, gap_char=["-"]):
                score=a.score
                start=a.start
                end=a.end
                print(format_alignment(*a))
                break

        # With matrix
        elif mode == "weighted":
            # Read the CSV file with the matrix information
            try:
                import pkg_resources
                stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(matrix))
                simMatrix = pickle.load(stream)
            except:
                with open('data/{}'.format(matrix), 'rb') as handle:
                    simMatrix = pickle.load(handle)

            print("Alignment with matrix")
            for a in pairwise2.align.globaldx(peptideList1, peptideList2, simMatrix, gap_char=["-"]):
                score = a.score
                start = a.start
                end = a.end
                print(format_alignment(*a))
                break

        return score,start,end

    ############################################################################
    @staticmethod
    def similarity_pair(peptide1, peptide2, matrix='matrix.pickle'):
        """
        Function to calculate similarity between two peptide sequences

        :param peptide1: Sequence of one of the input peptides
        :param peptide2: Sequence of the second input peptide
        :param matrix: BioPython matrix chosen for running the analysis

        :return sim_val: similarity based on the alignments between the peptides and themselves - Max value is 1
        """

        # Alignment between peptide1 and peptide2
        pep = Sequence(peptide1)
        pep.align_position_matrix(peptide2, matrix)
        score1_2 = pep.score_matrix

        # Alignment between peptide1 with itself
        pep.align_position_matrix(peptide1)
        score1_1 = pep.score_matrix

        # Alignment between peptide2 with itself
        pep2 = Sequence(peptide2)
        pep2.align_position_matrix(peptide2, matrix)
        score2_2 = pep2.score_matrix

        # Calculate similarity value
        sim_val = float(score1_2) / math.sqrt(float(score1_1 * score2_2))

        # Return similarity
        return sim_val

    # End of Alignment class
    ############################################################

##########################################################################

class pepDescriptors:
    """
    Class to calculate descriptors for peptide containing non-natural amino acids
    """

    def __init__(self, sequence, monomers='monomers.csv', properties="monomers_properties.pkl", generate_properties=False):
        """
        Start class by assigning general values

        :param sequence: Sequence of the peptide in BILN format
        :param generate_properties:Generate or not a pickle file with the pre-calculated properties. By default the file is provided
        """

        # Check if the peptide is in BILN format. This can be checked with the pyPept package.
        if "-" not in sequence:
            print("The monomers should be separated by dash symbols")
            exit(1)
        self.sequence = sequence

        if generate_properties==True:
            print("Calculating the properties of the monomers ...")
            # Read the CSV file with the monomer information
            try:
                import pkg_resources
                stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(monomers))
                df_initial = pd.read_csv(stream)
            except:
                with open('data/{}'.format(monomers), 'rb') as handle:
                    df_initial = pd.read_csv(handle)

            dict_df = {'name':[], 'smiles':[], 'mw':[], 'logp':[], 'nrot':[], 'tpsa':[]}
            for i, idx in enumerate(df_initial.index):
                name = df_initial.at[idx, 'Symbol']
                smiles = df_initial.at[idx, 'SMILES']
                mol = Chem.MolFromSmiles(smiles)

                # Calculate the descriptors
                mol_weight = Descriptors.MolWt(mol)
                mol_logp = Descriptors.MolLogP(mol)
                mol_tpsa = Descriptors.TPSA(mol)
                mol_nrot = Descriptors.NumRotatableBonds(mol)

                dict_df['name'].append(name)
                dict_df['smiles'].append(smiles)
                dict_df['mw'].append(mol_weight)
                dict_df['logp'].append(mol_logp)
                dict_df['nrot'].append(mol_nrot)
                dict_df['tpsa'].append(mol_tpsa)

            self.df = pd.DataFrame(dict_df)
            if not os.path.exists("data"):
                os.makedirs("data")
            self.df.to_pickle("data/{}".format(properties))
        else:
            try:
                import pkg_resources
                stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(properties))
                self.df = pickle.load(stream)
            except:
                with open('data/{}'.format(properties), 'rb') as handle:
                    self.df = pickle.load(handle)

    ##########################################################################
    def moranCorrelation(self, nlag=5):
        """
        Calculate amino-acid based moran correlation descriptors

        :param nlag: number of neighbors accounted for the descriptors. The values should be smaller than the full peptide length
        :return descriptors: a dictionary with the numerical descriptors
        """
        fullMonomers = re.split('[.-]', self.sequence)

        monomers_raw = self.sequence.split('-')
        monomers = []
        for m in monomers_raw:
            mon = re.sub("\(\d+,\d+\)", "", m)
            monomers.append(mon)

        totalMon = len(monomers)

        if totalMon == 1:
            print('Warning: The sequence only have one amino acid and no Moran correlation descriptors will be calculated')

        if totalMon < nlag + 1:
            print('Warning: the sequence should be larger than nlag+1: ' + str(nlag + 1) + '. The nlag was refactored based on the chain length')
            nlag = totalMon -1

        AAidxName = ['nrot', 'logp', 'tpsa', 'mw']

        AAidx = []
        for j in AAidxName:
            property_values = []
            for i, res in enumerate(fullMonomers):
                mon = re.sub("\(\d+,\d+\)", "", res)
                property_values.append(self.df.loc[self.df['name'] == mon, j].item())
            AAidx.append(property_values)

        AAidx1 = np.array([float(j) for i in AAidx for j in i])
        AAidx = AAidx1.reshape((len(AAidx), len(fullMonomers)))

        propMean = np.mean(AAidx, axis=1)
        propStd = np.std(AAidx, axis=1)

        for i in range(len(AAidx)):
            for j in range(len(AAidx[i])):
                AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]

        index = {}
        for i in range(len(fullMonomers)):
            mon = re.sub("\(\d+,\d+\)", "", fullMonomers[i])
            index[mon] = i

        # Calculate the descriptors and include them in the dictionary
        descriptors = {}
        for prop in range(len(AAidxName)):
            xmean = sum([AAidx[prop][index[aa]] for aa in monomers]) / totalMon
            for n in range(1, nlag + 1):
                if totalMon > nlag:
                    num = sum([(AAidx[prop][index.get(monomers[j], 0)] - xmean) *
                                 (AAidx[prop][index.get(monomers[j + n], 0)] - xmean)
                                 for j in range(totalMon - n)]) / (totalMon - n)
                    den = sum(
                        [(AAidx[prop][index.get(monomers[j], 0)] - xmean) ** 2 for j in range(totalMon)]) / totalMon
                    rn = num / den
                else:
                    rn = 'NA'
                descriptors[AAidxName[prop] + '-lag' + str(n)] = rn

        return descriptors

    # End of pepDescriptors class
    ############################################################

########################################################################################
# Additional function
########################################################################################

def readProperties(properties="monomers_properties.pkl"):
    """
    Read the pickle file with the monomer properties

    :param properties: Name of the file
    :return df: dataframe with the monomer properties
    """
    try:
        import pkg_resources
        stream = pkg_resources.resource_stream(__name__, 'data/{}'.format(properties))
        df = pickle.load(stream)
    except:
        with open('data/{}'.format(properties), 'rb') as handle:
            df = pickle.load(handle)
    return df

########################################################################################
def peptideFromSMILES(smiles,df):
    """
    Class to convert from smiles to peptide sequence. Only for natural amino acids

    :param smiles: SMILES to convert into sequence
    :param df: dataframe with monomer properties

    :return final_seq: Sequence of amino acids
    """

    list_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    m = Chem.MolFromSmiles(smiles)

    # Pattern of the AA backbone
    peptide_bond_representation = Chem.MolFromSmarts('NCC(=O)N')
    am = np.array(Chem.GetAdjacencyMatrix(m))

    masses = []
    for peptide_bond in m.GetSubstructMatches(peptide_bond_representation):

        alpha = peptide_bond[1]
        nitrogens = set([peptide_bond[0],peptide_bond[-1]])

        aa_atom_idx = set([alpha])
        set2 = set()
        while aa_atom_idx != set2:
            set2 = aa_atom_idx.copy()
            temp_am = am[:, list(aa_atom_idx)]
            aa_atom_idx = set(np.where(temp_am==1)[0]) | aa_atom_idx
            aa_atom_idx -= nitrogens
        aa_atom_idx.add(peptide_bond[0])

        bonds = []
        for i,j in combinations(aa_atom_idx, 2):
            b = m.GetBondBetweenAtoms(int(i),int(j))
            if b: bonds.append(b.GetIdx())

        tempMol = Chem.PathToSubmol(m, bonds)
        mass = Descriptors.MolWt(tempMol)
        masses.append(mass)

    # Generate the final sequence based on the masses
    final_pep=[]
    for i,m in enumerate(masses):
        if i == 0:
            new_mass = m + 17.01
            names=[]
            for i, idx in enumerate(df.index):
                name = df.at[idx, 'name']
                mw = df.at[idx, 'mw']
                if new_mass - 0.03 < mw < new_mass + 0.03:
                    if name in list_aa:
                        names.append(name)
            final_pep.append(names[0])
        else:
            new_mass = m + 18.01
            names=[]
            for i, idx in enumerate(df.index):
                name = df.at[idx, 'name']
                mw = df.at[idx, 'mw']
                if new_mass - 0.03 < mw < new_mass + 0.03:
                    if name in list_aa:
                        names.append(name)
            final_pep.append(names[0])

    final_seq = ''.join(final_pep)
    return final_seq

    # End of additional functions
    ############################################################

############################################################
## End of extra.py
############################################################
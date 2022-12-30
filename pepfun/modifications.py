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

# System dependent
import os
import math
import warnings
import numpy as np

# Modeller
import modeller
from modeller.automodel import *

# BioPython
from Bio.Align import PairwiseAligner
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, calc_angle

########################################################################################
# Classes and functions
########################################################################################

class Filling:

    @staticmethod
    def fill_peptide_complex(code, template, model, path="", output_name=""):
        """
        Function to model new monomers of a bound peptide. It is assumed the peptide is the second smaller chain in the PDB

        Arguments:
        :param code: code of the structure used to model the peptide
        :param template: sequence of the peptide in the complex structure
        :param model: new sequence we want to model
        :param output_name: name of the final modelled structure

        :return: PDB structure of the protein-peptide complex with new amino acids
        """
        print("Fill a peptide in complex with a protein target")

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1]=="/":
                os.system("cp {}{}.pdb .".format(path,code))
            else:
                os.system("cp {}/{}.pdb .".format(path,code))

        structure = parser.get_structure('PEP', code + ".pdb")
        chains=[]
        for i,chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 1:
                if chSeq != template:
                    print("The second chain is not the peptide sequence. Please verify the input")
                    exit(1)
        if len(chains) != 2:
            print("The complex does not have exactly 2 chains (protein and peptide). Please verify the input")
            exit(1)

        aligner = PairwiseAligner()
        alignments = aligner.align(template, model)
        alignment = alignments[0]
        pepT=str(alignment).split()[0]
        pepM=str(alignment).split()[2]
        print("The template sequence is: {}, and the sequence to model is: {}".format(pepT, pepM))

        # Start Modeller environment and write the sequence file
        e = modeller.Environ()
        m = modeller.Model(e, file=code)
        aln = modeller.Alignment(e)
        aln.append_model(m, align_codes=code)
        aln.write(file=code +'.seq')

        # Obtain the protein sequence with the peptide
        infoSeq =[x.strip() for x in open(code +'.seq')]
        header =[]
        sequenceLine =''
        for info in infoSeq:
            if ">" not in info and ":" not in info:
                if info:
                    sequenceLine += info
            else:
                if info: header.append(info)

        # Split the sequence information as required by Modeller
        last_line =sequenceLine.split("/")
        sequenceTemp =last_line[0 ] +"/ " +pepT +"*"
        sequenceMod =last_line[0 ] +"/ " +pepM +"*"

        seqTempList = [sequenceTemp[i: i +75] for i in range(0, len(sequenceTemp), 75)]
        seqModList = [sequenceMod[i: i +75] for i in range(0, len(sequenceMod), 75)]

        # Create the alignment file
        alignmentFile =open("alignment.ali" ,"w")
        for h in header: alignmentFile.write( h +"\n")
        for s in seqTempList: alignmentFile.write( s +"\n")
        alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
        for s in seqModList: alignmentFile.write( s +"\n")
        alignmentFile.close()

        # Directories for input atom files
        e.io.atom_files_directory = ['.', '../atom_files']
        a = AutoModel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
        a.starting_model= 1
        a.ending_model  = 1
        a.make()

        # Export the generated model
        if output_name:
            if ".pdb" in output_name:
                os.system('mv {}_fill.B99990001.pdb {}'.format(code, output_name))
            else:
                os.system('mv {}_fill.B99990001.pdb {}.pdb'.format(code, output_name))
        else:
            os.system('mv {}_fill.B99990001.pdb filled_{}.pdb'.format(code, code))

        os.system("rm {c}.seq {c}_fill.D00000001 {c}_fill.ini {c}_fill.rsr {c}_fill.sch {c}_fill.V99990001 alignment.ali".format(c=code))
        if path:
            os.system("rm {}.pdb".format(code))


    ########################################################################################
    @staticmethod
    def fill_peptide_alone(code, template, model, path="", output_name=""):
        """
        Function to model new monomers of a single peptide.

        Arguments:
        :param code: code of the structure used to model the peptide
        :param template: sequence of the peptide
        :param model: new sequence we want to model
        :param output_name: name of the final modelled structure

        :return: PDB structure of the peptide with new amino acids
        """
        print("Fill a single peptide")

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1]=="/":
                os.system("cp {}{}.pdb .".format(path,code))
            else:
                os.system("cp {}/{}.pdb .".format(path,code))

        structure = parser.get_structure('PEP', code + ".pdb")
        chains = []
        for i, chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 0:
                if chSeq != template:
                    print("The template peptide is not alone in the structure. Please verify the input")
                    exit(1)
        if len(chains) > 1:
            print("The PDB must have only 1 chain with the peptide")
            exit(1)

        aligner = PairwiseAligner()
        alignments = aligner.align(template, model)
        alignment = alignments[0]
        pepT = str(alignment).split()[0]
        pepM = str(alignment).split()[2]
        print("The template sequence is: {}, and the sequence to model is: {}".format(pepT, pepM))

        # Start the Modeller environment
        e = modeller.Environ()
        m = modeller.Model(e, file=code)
        aln = modeller.Alignment(e)
        aln.append_model(m, align_codes=code)
        aln.write(file=code+'.seq')

        # Edit the information of the sequence to store the sequences in the Modeller format
        infoSeq=[x.strip() for x in open(code+'.seq')]
        header=[]
        sequenceLine=''
        for info in infoSeq:
            if ">" not in info and ":" not in info:
                if info:
                    sequenceLine+=info
            else:
                if info: header.append(info)

        # Store the sequences in variables according to Modeller format
        sequenceTemp=pepT+"*"
        sequenceMod=pepM+"*"

        seqTempList = [sequenceTemp]
        seqModList = [sequenceMod]

        # Create the alignment file
        alignmentFile=open("alignment.ali","w")
        for h in header: alignmentFile.write(h+"\n")
        for s in seqTempList: alignmentFile.write(s+"\n")
        alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
        for s in seqModList: alignmentFile.write(s+"\n")
        alignmentFile.close()

        # Directories for input atom files
        e.io.atom_files_directory = ['.', '../atom_files']
        a = AutoModel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
        a.starting_model= 1
        a.ending_model  = 1
        a.make()

        if output_name:
            if ".pdb" in output_name:
                os.system('mv {}_fill.B99990001.pdb {}'.format(code, output_name))
            else:
                os.system('mv {}_fill.B99990001.pdb {}.pdb'.format(code, output_name))
        else:
            os.system('mv {}_fill.B99990001.pdb filled_{}.pdb'.format(code, code))

        os.system("rm {c}.seq {c}_fill.D00000001 {c}_fill.ini {c}_fill.rsr {c}_fill.sch {c}_fill.V99990001 alignment.ali".format(c=code))
        if path:
            os.system("rm {}.pdb".format(code))

    # End of Filling class
    ############################################################

########################################################################################
class Capping:

    ########################################################################################
    @staticmethod
    def capping_peptide(pdb_file, peptide, path="", output_name = "", mode='both'):
        '''
        Function to cap an existing PDB file

        :param pdb_file: PDB of the peptide
        :param peptide: sequence of the original peptide
        :param output_name: name of the output capped structure
        :param mode: Mode to select which positions to cap. Choose between, N-term, C-term, or both

        :return: a peptide capped structure
        '''

        print("Capping a peptide with ACE and/or NME groups")

        template=peptide
        if mode in ['N-term','C-term','both']:
            if mode=='N-term':
                model='G' + peptide
            elif mode=='C-term':
                model = peptide + 'G'
            elif mode=='both':
                model = 'G' + peptide + 'G'
        else:
            print("The mode option should be N-term, C-term or both. Please check the option")
            exit(1)
        code=pdb_file.split('.')[0]

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1] == "/":
                os.system("cp {}{} .".format(path, pdb_file))
            else:
                os.system("cp {}/{} .".format(path, pdb_file))
        structure = parser.get_structure('PEP', pdb_file)
        chains = []
        for i, chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 0:
                if chSeq != peptide:
                    print("The peptide sequence does not coincide with the PDB file. Please verify input")
                    exit(1)
        try:
            Filling.fill_peptide_alone(code, template, model, output_name='modeller_peptide')
        except Exception as e:
            print("Failed to fill the peptide with cap groups, Error: {}".format(e))
            exit(1)

        parser = PDBParser()

        reference = parser.get_structure('REF', "modeller_peptide.pdb")
        model = reference[0]
        chain = model['A']

        for i, res in enumerate(chain):
            if i == 0:
                if mode == 'both' or mode == 'N-term':
                    # Delete the other atoms leaving only the atoms of the backbone
                    ids = []
                    for a in res:
                        atomId = a.id
                        if atomId not in ("CA", "O", "C"): ids.append(atomId)
                    for i in ids: res.detach_child(i)

                    res.resname = 'ACE'

                    for a in res:
                        atomId = a.id
                        if atomId == 'CA':
                            a.fullname = 'CH3'
                            a.name = 'CH3'
                            a.element = 'C'

            if i == len(chain)-1:
                if mode == 'both' or mode == 'C-term':
                    # Delete the other atoms leaving only the atoms of the backbone
                    ids = []
                    for a in res:
                        atomId = a.id
                        if atomId not in ("CA", "N"): ids.append(atomId)
                    for i in ids: res.detach_child(i)

                    res.resname = 'NME'

                    for a in res:
                        atomId = a.id
                        if atomId == 'CA':
                            a.fullname = ' C'
                            a.name = ' C'
                            a.element = 'C'

        # Saving the new structure
        io = PDBIO()
        io.set_structure(reference)
        if output_name:
            if ".pdb" in output_name:
                io.save("{}".format(output_name))
            else:
                io.save("{}.pdb".format(output_name))
        else:
            io.save("capped_{}".format(pdb_file))

        os.system("rm modeller_peptide.pdb")
        if path:
            os.system("rm {}".format(pdb_file))

    # End of Capping class
    ############################################################

########################################################################################

class Mutation:

    def __init__(self, pdb_file, chainID, pepPosition, pdbNNAA, chainNNAA, nnaa_name):
        '''
        Class to mutate a NNAA in a position of a natural AA

        :param pdb_file: PDB file of the peptide that will be mutated
        :param chainID: chain ID of the peptide
        :param pepPosition: position on the peptide that will be mutated
        :param pdbNNAA: PDB file of the NNAA used in the mutation
        :param chainNNAA: chain of the NNAA PDB file
        :param nnaa_name: abbreviated name of the NNAA
        '''

        self.pdb_file = pdb_file
        self.chainID = chainID
        self.pepPosition = pepPosition
        self.pdbNNAA = pdbNNAA
        self.chainNNAA = chainNNAA
        self.nnaa_name = nnaa_name

    ########################################################################################

    def get_bb_dihedral(self):
        '''
        Function to get the backbone dihedral from the structure that will be mutated

        :return: The backbone dihedrals of the input PDB file
        '''

        # Read the reference PDB file
        parser = PDBParser()
        reference = parser.get_structure('REF',self.pdb_file)
        model = reference[0]
        chain = model[self.chainID]

        # Inititalize variables
        rad = 180.0 / math.pi
        flagCap = False
        for i,res in enumerate(chain):
            # First residue
            if i == 0:
                try:
                    phi = 0
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()

                    N_next = chain[i+2]["N"]
                    n_next = N_next.get_vector()

                    psi = calc_dihedral(n, ca, c, n_next)

                    # Store the data for the position of interest
                    if self.pepPosition == i + 1:
                        print(c,n,ca)
                        print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(),self.pepPosition, phi * rad, psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad

                    # Update the previous coordinates
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                except:
                    flagCap=True

            # Last residue
            elif i == len(chain)-1:
                try:
                    n1 = N_prev.get_vector()
                    ca1 = CA_prev.get_vector()
                    c1 = C_prev.get_vector()

                    # Get current AA coordinates
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()

                    # Calculate dihedral
                    psi = 0
                    phi = calc_dihedral(c1, n, ca, c)

                    # Store the data for the position of interest
                    if self.pepPosition == i + 1:
                        print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(),self.pepPosition, phi * rad,psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad
                except:
                    flagCap=True
            else:
                if flagCap is False:
                    # Get previous AA coordinates
                    n1 = N_prev.get_vector()
                    ca1 = CA_prev.get_vector()
                    c1 = C_prev.get_vector()

                    # Get current AA coordinates
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()

                    N_next = chain[i + 2]["N"]
                    n_next = N_next.get_vector()

                    # Calculate dihedral
                    psi = calc_dihedral(n, ca, c, n_next)
                    phi = calc_dihedral(c1, n, ca, c)

                    # Store the data for the position of interest
                    if self.pepPosition == i+1:
                        print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(), self.pepPosition, phi * rad, psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad

                    # Update the previous coordinates
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                else:
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                    flagCap=False

        try:
            self.near_phi = int(self.ref_phi / 10) * 10
            self.near_psi = int(self.ref_psi / 10) * 10
        except:
            print("Failed to map the backbone dihedrals. Check the input variables")

    ########################################################################################

    def get_NNAA_basic_info(self):
        '''
        Function to capture all possible information from the NNAA pdb, including atoms, bonds, and structure parameters

        :return: Lists and dictionaries with all the required parameters
        '''
        # Input PDB information
        rad = 180.0 / math.pi
        parser = PDBParser()
        reference = parser.get_structure('REF', self.pdbNNAA)
        model = reference[0]
        chain = model[self.chainNNAA]

        # Get all the atoms in the NNAA and store those that are not hydrogens
        backbone_atoms = ("N", "O", "C", "OXT")
        atom_list = [x for x in model.get_atoms()]
        key_atoms = []
        atom_names = []
        self.atom_order = []
        for atom in atom_list:
            if atom.get_id()[0] != 'H':
                key_atoms.append(atom)
                atom_names.append(atom.get_id())
                if atom.get_id() not in backbone_atoms:
                    self.atom_order.append(atom.get_id())
        # Create an engine to look for bonds
        ns = NeighborSearch(key_atoms)
        _cutoff_dist = 2

        # Loop to store the bonds

        self.bonds = []
        self.assoc_bonds = {}
        for target in key_atoms:
            name = target.get_id()
            if name not in backbone_atoms:
                close_atoms = ns.search(target.coord, _cutoff_dist)
                for close_atom in close_atoms:
                    name2 = close_atom.get_id()
                    if name2 not in backbone_atoms and target != close_atom:
                        if (name, name2) not in self.bonds and (name2, name) not in self.bonds:
                            self.bonds.append((name,name2))
                            self.assoc_bonds[name2]=name

        self.dict_NNAA = {}

        for bond in self.bonds:

            self.dict_NNAA[bond]={}
            ele1 = bond[0]
            ele2 = bond[1]
            ind1 = self.atom_order.index(ele1)
            # Exception
            if ele1 != "CA" and self.assoc_bonds[ele1] == "CA":
                ind1 = 1

            for i,res in enumerate(chain):
                if ind1 == 0:
                    atom_1 = res["N"]
                    atom_2 = res["C"]
                    atom_3 = res[ele1]
                    atom_4 = res[ele2]
                elif ind1 == 1:
                    atom_1 = res["N"]
                    atom_2 = res[self.assoc_bonds[ele1]]
                    atom_3 = res[ele1]
                    atom_4 = res[ele2]
                else:
                    atom_1 = res[self.assoc_bonds[self.assoc_bonds[ele1]]]
                    atom_2 = res[self.assoc_bonds[ele1]]
                    atom_3 = res[ele1]
                    atom_4 = res[ele2]

                atom1_vec = atom_1.get_vector()
                atom2_vec = atom_2.get_vector()
                atom3_vec = atom_3.get_vector()
                atom4_vec = atom_4.get_vector()

                length = atom_4 - atom_3
                angle = calc_angle(atom2_vec, atom3_vec, atom4_vec) * rad
                diangle = calc_dihedral(atom1_vec, atom2_vec, atom3_vec, atom4_vec) * rad

                self.dict_NNAA[bond]={'length': length, 'angle': angle, 'diangle': diangle}

    ########################################################################################

    def assign_mutation(self):
        '''
        Function to assign the mutation based on the stored parameters

        :return: A mutated PDB file that can be minimized using external tools
        '''

        parser = PDBParser()
        reference = parser.get_structure('REF', self.pdb_file)
        model = reference[0]
        chain = model[self.chainID]

        for i, res in enumerate(chain):
            if self.pepPosition == i + 1:
                C = res["C"]
                N = res["N"]
                CA = res["CA"]

                # Delete the other atoms leaving only the atoms of the backbone
                ids = []
                for a in res:
                    atomId = a.id
                    if atomId not in ("N", "CA", "O", "C"): ids.append(atomId)
                for i in ids: res.detach_child(i)

                res.resname = self.nnaa_name
                atom_to_select = {'CA': CA}
                for bond in self.bonds:

                    # dict_NNAA[bond] = {}
                    ele1 = bond[0]
                    ele2 = bond[1]
                    ind1 = self.atom_order.index(ele1)
                    if ele1 != "CA" and self.assoc_bonds[ele1] == "CA":
                        ind1 = 1

                    if ind1 == 0:

                        atom_coord = calculateCoordinates(
                                N, C, atom_to_select[ele1], self.dict_NNAA[bond]['length'], self.dict_NNAA[bond]['angle'],
                                self.dict_NNAA[bond]['diangle']
                        )
                        atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                        if ele2 not in atom_to_select:
                            atom_to_select[ele2] = atom_object
                            res.add(atom_object)

                    elif ind1 == 1:

                        atom_coord = calculateCoordinates(
                                N, atom_to_select[self.assoc_bonds[ele1]], atom_to_select[ele1],
                                self.dict_NNAA[bond]['length'],
                                self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                        )
                        atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                        if ele2 not in atom_to_select:
                            atom_to_select[ele2] = atom_object
                            res.add(atom_object)

                    else:
                        atom_coord = calculateCoordinates(
                                atom_to_select[self.assoc_bonds[self.assoc_bonds[ele1]]], atom_to_select[self.assoc_bonds[ele1]],
                                atom_to_select[ele1], self.dict_NNAA[bond]['length'],
                                self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                        )
                        if ele2[0] not in ('C', 'N', 'O', 'S'):
                            if ele2 == 'BRZ':
                                ele2mod = 'BrZ'
                            else:
                                ele2mod = ele2
                            atom_object = Atom(ele2mod, atom_coord, 0.0, 1.0, " ", " {}".format(ele2mod), 0,
                                               ele2[0:len(ele2) - 1])
                        else:
                            atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                        if ele2 not in atom_to_select:
                            atom_to_select[ele2] = atom_object
                            res.add(atom_object)

        # Saving the new structure
        io = PDBIO()
        io.set_structure(reference)
        io.save("mutated_{}.pdb".format(self.nnaa_name))
        print("The mutation has been generated for the NNAA: {}".format(self.nnaa_name))

    # End of Mutation class
    ############################################################

########################################################################################
# Function derived from the PeptideBuilder module
########################################################################################

def calculateCoordinates(
    refA: Residue, refB: Residue, refC: Residue, L: float, ang: float, di: float
) -> np.ndarray:
    '''
    Function from PeptideBuilder to generate the coordinates of the NNAA based on the stored parameters
    '''
    AV = refA.get_vector()
    BV = refB.get_vector()
    CV = refC.get_vector()

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(
        math.pow((B * BZ - BY * G), 2)
        * (
            -(F * F) * (A * A + B * B + G * G)
            + (
                B * B * (BX * BX + BZ * BZ)
                + A * A * (BY * BY + BZ * BZ)
                - (2 * A * BX * BZ * G)
                + (BX * BX + BY * BY) * G * G
                - (2 * B * BY) * (A * BX + BZ * G)
            )
            * L
            * L
        )
    )
    denom = (
        (B * B) * (BX * BX + BZ * BZ)
        + (A * A) * (BY * BY + BZ * BZ)
        - (2 * A * BX * BZ * G)
        + (BX * BX + BY * BY) * (G * G)
        - (2 * B * BY) * (A * BX + BZ * G)
    )

    X = (
        (B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const
    ) / denom

    if (B == 0 or BZ == 0) and (BY == 0 or G == 0):
        const1 = math.sqrt(
            G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X))
        )
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = (
            (A * A * BY * F) * (B * BZ - BY * G)
            + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const)
            - A * (B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)
        ) / ((B * BZ - BY * G) * denom)
        Z = (
            (A * A * BZ * F) * (B * BZ - BY * G)
            + (B * F) * math.pow(B * BZ - BY * G, 2)
            + (A * BX * F * G) * (-B * BZ + BY * G)
            - B * BX * const
            + A * BY * const
        ) / ((B * BZ - BY * G) * denom)

    # Get the new Vector from the origin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()

############################################################
## End of modifications.py
############################################################
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
import subprocess
import os

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# Modeller
import modeller
from modeller import *
from modeller.automodel import *

########################################################################################
# Classes and Functions
########################################################################################

class Conformer:

    def __init__(self, sequence, auxiliar_path):
        """
        Inititalize the class calculate conformers using RDKit or Modeller

        :param sequence: Peptide sequence
        :param auxiliar_path: Path of the auxiliar folder

        :return Based on the sequence, the class will generate a PDB file of the peptide structure
        """
        self.sequence = sequence
        self.length_peptide = len(self.sequence)
        self.path = auxiliar_path

    ############################################################################
    def generate_conformer(self):
        """
        Function to generate basic conformer using RDKit

        :return Structure predicted in PDB format
        """

        # Generate molecule using the sequence and HELM notation
        helm = ".".join(list(self.sequence))
        mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
        mol.SetProp("_Name", self.sequence)
        mol = Chem.AddHs(mol)

        # Generate the conformer using UFF force field
        print("Generating the basic conformer for peptide {}".format(self.sequence))
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xf00d
        AllChem.EmbedMolecule(mol, ps)

        # Write PDB file
        molfile = open('structure.pdb', 'w')
        molfile.write(Chem.MolToPDBBlock(mol))
        molfile.close()

        # Add Remarks
        remarks = open('remarks.pdb', 'w')
        remarks.write('REMARK    Conformer generated with RDKit\n')
        remarks.write('REMARK    The conformer is created using a distance geometry approach.\n')
        remarks.write('REMARK    The method is explained in the PepFun2 paper.\n')
        remarks.close()

        # Concatenate final structure
        with open("structure_{}.pdb".format(self.sequence), "w") as outfile:
            for filename in ["remarks.pdb", "structure.pdb"]:
                with open(filename) as infile:
                    contents = infile.read()
                    outfile.write(contents)

        # Delete files
        os.remove("remarks.pdb")
        os.remove("structure.pdb")

    ########################################################################################

    def run_psipred(self):
        '''
        Run PSIPRED locally to assign the secondary structure

        :return: A string with the predicted SS. H: Helix, E: Strand,Sheet, C: Coil
        '''
        # Run PSIPRED
        fastaFile = open('{}.fasta'.format(self.sequence), 'w')
        fastaFile.write('>sequence\n')
        fastaFile.write('{}\n'.format(self.sequence))
        fastaFile.close()

        if self.path[-1] == "/":
            os.system('{}runpsipred_single {}.fasta {}'.format(self.path, self.sequence,self.path))
        else:
            os.system('{}/runpsipred_single {}.fasta {}'.format(self.path, self.sequence, self.path))
        bash = 'grep Pred {}.horiz | cut -f 2 -d " "'.format(self.sequence)
        ss_PSI = str(subprocess.check_output(['bash', '-c', bash]).strip().decode('utf-8'))
        os.system('rm {}*'.format(self.sequence))

        print(f"The predicted secondary structure is: {ss_PSI}")
        return ss_PSI

    ############################################################################
    def prepare_modeller(self, ss_ref):
        '''
        Protocol to generate all the parameters required for Modeller, including the RDKit template

        :param ss_ref: Secondary structure assigned to the peptide
        :return: Input parameters to run Modeller
        '''
        pos = int(len(self.sequence) / 2)
        pepT = ""
        pepM = ""
        for i, ch in enumerate(self.sequence):
            if i < pos:
                pepT += "-"
                pepM += ch
            elif i == pos:
                pepT += ch
                pepM += ch
            else:
                pepT += "-"
                pepM += ch

        amino = self.sequence[pos:pos + 1]
        baseMol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % amino)
        #baseMol = Chem.AddHs(baseMol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xf00d
        AllChem.EmbedMolecule(baseMol, ps)

        # Write PDB file
        molfile = open('template.pdb', 'w')
        molfile.write(Chem.MolToPDBBlock(baseMol))
        molfile.close()

        rangeHelix = []
        rangeBeta = []
        if ss_ref.count('H') > 2 or ss_ref.count('G') > 2:
            positions = {}
            count_groups = 0
            prev = '-'
            for i, ele in enumerate(ss_ref):
                if ele == 'H' or ele == 'G':
                    if ele != prev:
                        count_groups += 1
                        positions[count_groups] = [i + 1]
                    elif ele == prev:
                        positions[count_groups].append(i + 1)
                prev = ele
            rangeHelix = []
            for group in positions:
                rangeHelix.append((min(positions[group]), max(positions[group])))

        if ss_ref.count('E') > 2 or ss_ref.count('B') > 2:
            positions = {}
            count_groups = 0
            prev = '-'
            for i, ele in enumerate(ss_ref):
                if ele == 'E' or ele == 'B':
                    if ele != prev:
                        count_groups += 1
                        positions[count_groups] = [i + 1]
                    elif ele == prev:
                        positions[count_groups].append(i + 1)
                prev = ele
            rangeBeta = []
            for group in positions:
                rangeBeta.append((min(positions[group]), max(positions[group])))

        # To add a cycle constraint
        rangeCycle = ()
        positions = []
        for i, aa in enumerate(self.sequence):
            if aa == 'C':
                if i <= 1 or i >= len(self.sequence) - 2:
                    positions.append(i + 1)
        if len(positions) == 2:
            rangeCycle = (positions[0], positions[1])

        return pepT, pepM, rangeHelix, rangeBeta, rangeCycle

    ########################################################################################

    def modelling_modeller(self, pepT, pepM, rangeHelix, rangeBeta, rangeCycle):
        '''
        Function to use modeller in ab-initio mode based on a peptide template

        :param pepT: Sequence used as template
        :param pepM: Sequence to model
        :param rangeHelix: List of tuples with regions containing alpha-helix
        :param rangeBeta: List of tuples with regions containing beta sheets/strands
        :param rangeCycle: Tuple with the residues creating the cycle
        :return: A PDB file with the predicted peptide structure
        '''

        # Start the Modeller environment
        code = 'template'
        e = modeller.environ()
        m = modeller.model(e, file=code)
        aln = modeller.alignment(e)
        aln.append_model(m, align_codes=code)
        aln.write(file=code + '.seq')

        # Edit the information of the sequence to store the sequences in the Modeller format
        infoSeq = [x.strip() for x in open(code + '.seq')]
        header = []
        sequenceLine = ''
        for info in infoSeq:
            if ">" not in info and ":" not in info:
                if info:
                    sequenceLine += info
            else:
                if info: header.append(info)

        # Store the sequences in variables according to Modeller format
        sequenceTemp = pepT + "*"
        sequenceMod = pepM + "*"

        seqTempList = [sequenceTemp]
        seqModList = [sequenceMod]

        # Create the alignment file
        alignmentFile = open("alignment.ali", "w")
        for h in header: alignmentFile.write(h + "\n")
        for s in seqTempList: alignmentFile.write(s + "\n")
        alignmentFile.write("\n>P1;template_fill\nsequence:::::::::\n")
        for s in seqModList: alignmentFile.write(s + "\n")
        alignmentFile.close()

        # Directories for input atom files
        e.io.atom_files_directory = ['.', '../atom_files']

        if len(rangeHelix) >= 1 or len(rangeBeta) >= 1:
            class MyModel(automodel):
                def special_patches(self, aln):
                    if len(rangeCycle) >= 1:
                        self.patch(residue_type='DISU', residues=(self.residues['{}:A'.format(rangeCycle[0])],
                                                                  self.residues['{}:A'.format(rangeCycle[1])]))

                def special_restraints(self, aln):
                    rsr = self.restraints
                    at = self.atoms
                    if len(rangeHelix) >= 1:
                        for pairHelix in rangeHelix:
                            rsr.add(secondary_structure.alpha(
                                self.residue_range('{}:A'.format(pairHelix[0]), '{}:A'.format(pairHelix[1]))))

                    if len(rangeBeta) >= 1:
                        extremes = []
                        for j, pairBeta in enumerate(rangeBeta):
                            rsr.add(secondary_structure.strand(
                                self.residue_range('{}:A'.format(pairBeta[0]), '{}:A'.format(pairBeta[1]))))
                            if j == 0: extremes.append(pairBeta[0])
                            if j == 1: extremes.append(pairBeta[1])
                        rsr.add(secondary_structure.sheet(at['N:{}:A'.format(extremes[0])],
                                                          at['O:{}:A'.format(extremes[1])],
                                                          sheet_h_bonds=-5))
                    if len(rangeCycle) >= 1:
                        rsr.add(forms.gaussian(group=physical.xy_distance,
                                               feature=features.distance(at['SG:{}:A'.format(rangeCycle[0])],
                                                                         at['SG:{}:A'.format(rangeCycle[1])]),
                                               mean=2.0, stdev=0.1))

            a = MyModel(e, alnfile='alignment.ali', knowns='template', sequence='template_fill')
            a.starting_model = 1
            a.ending_model = 1
            a.make()
        else:
            a = automodel(e, alnfile='alignment.ali', knowns='template', sequence='template_fill')
            a.starting_model = 1
            a.ending_model = 1
            a.make()

        os.system("mv template_fill.B99990001.pdb modeller.pdb")
        os.system("rm template* alignment.ali")

        # Add Remarks
        remarks = open('remarks.pdb', 'w')
        remarks.write('REMARK    Conformer generated with Modeller\n')
        remarks.write('REMARK    The conformer uses secondary structure restraints obtained with PSIPred.\n')
        remarks.write('REMARK    The method is explained in the PepFun2 paper. Both Modeller and PSIPred can be used for academic purposes.\n')
        remarks.close()

        # Concatenate final structure
        with open("modeller_{}.pdb".format(self.sequence), "w") as outfile:
            for filename in ["remarks.pdb", "modeller.pdb"]:
                with open(filename) as infile:
                    contents = infile.read()
                    outfile.write(contents)

        # Delete files
        os.remove("remarks.pdb")
        os.remove("modeller.pdb")

    # End of Conformer class
    ############################################################

############################################################
## End of conformer.py
############################################################
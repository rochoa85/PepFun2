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
import numpy as np
import os
from igraph import *

# BioPython
from Bio.PDB import *

########################################################################################
# Additional function
########################################################################################

def renum_structure(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb_file)

    for i,chain in enumerate(structure[0]):
        for i,res in enumerate(structure[0][chain.id]):
            res.id = (' ', i+1, ' ')

    # Saving the new structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)
    os.system("head -n -1 {} > temp; mv temp {}".format(pdb_file,pdb_file))

########################################################################################
# Classes and functions
########################################################################################

class Interactions:

    """
    Class with functions to perform different type of analysis using a peptide structure alone or in complex with a protein
    """

    def __init__(self ,pdb_file ,chain, auxiliar_path):
        """
        Inititalize the class calculating some basic properties

        :param pdb_file: PDB file with the required information
        :param chain: chain containing the peptide in the file
        :param auxiliar_path: path of the auxiliar folder with the software

        :return Based on the structure, the class start with some counters, the sequence derived from the structure,
                its length and a dictionary with the aa positions
        """

        self.pdb_file =pdb_file
        renum_structure(self.pdb_file)
        self.chain =chain
        self.aminoacids_back ={"ALA" :"A" ,"ASP" :"D" ,"GLU" :"E" ,"PHE" :"F" ,"HIS" :"H" ,"ILE" :"I" ,"LYS" :"K"
                                ,"LEU" :"L" ,"MET" :"M" ,"GLY" :"G",
                              "ASN" :"N" ,"PRO" :"P" ,"GLN" :"Q" ,"ARG" :"R" ,"SER" :"S" ,"THR" :"T" ,"VAL" :"V"
                                ,"TRP" :"W" ,"TYR" :"Y" ,"CYS" :"C"}
        self.aminoacids ={"A" :"ALA" ,"D" :"ASP" ,"E" :"GLU" ,"F" :"PHE" ,"H" :"HIS" ,"I" :"ILE" ,"K" :"LYS"
                           ,"L" :"LEU" ,"M" :"MET" ,"G" :"GLY",
                         "N" :"ASN" ,"P" :"PRO" ,"Q" :"GLN" ,"R" :"ARG" ,"S" :"SER" ,"T" :"THR" ,"V" :"VAL"
                           ,"W" :"TRP" ,"Y" :"TYR" ,"C" :"CYS"}
        self.sequence =""
        self.initial_id_residue =0
        self.positions ={}
        self.path=auxiliar_path

        # Read the structure in BioPython format
        parser = PDBParser()
        self.reference = parser.get_structure('REF' ,self.pdb_file)
        for ch in self.reference[0]:
            if ch.get_id() == chain:
                for i ,residue in enumerate(ch):
                    try:
                        seq = self.aminoacids_back[residue.get_resname()]
                    except:
                        seq = "["+residue.get_resname()+"]"
                    # Save the sequence
                    self.sequence =self.sequence +str(seq)
                    if i== 0: self.initial_id_residue = residue.get_full_id()[3][1]
                    # Save the positions
                    self.positions[i + 1] = {"aa": str(seq), "full_aa": residue.get_resname()}

        # Store sequence length
        self.len_sequence = len(self.sequence)

    ############################################################################
    def get_secondary_structure(self):
        """
        Function to calculate the secondary structure and the accessible surface area using the auxiliary program mkdssp

        :return
        The dictionary positions will store the calculated data of dssp and the asa values
        total_dssp will contain the complete predicted secondary structure with the following conventions:
        B - beta bridge
        H - alpha helix
        E - beta strand
        S - bend
        T - turn
        G - 3/10 helix
        """

        model = self.reference[0]
        if self.path[-1] == "/":
            dssp = DSSP(model, self.pdb_file, dssp=self.path+"mkdssp")
        else:
            dssp = DSSP(model, self.pdb_file, dssp=self.path + "/mkdssp")
        self.total_dssp = ""

        # Loop over the keys from the dssp response to store ss and asa values
        counter = 1
        for keys in list(dssp.keys()):
            if keys[0] == self.chain:
                #position = counter
                position = int(keys[1][1])
                self.positions[position]["dssp"] = dssp[keys][2]
                self.positions[position]["asa"] = dssp[keys][3]
                self.total_dssp = self.total_dssp + dssp[keys][2]
                counter += 1

    ############################################################################

    def get_hydrogen_bonds(self, output=True):
        """
        Function to calculate hydrogen bonds with the other chains

        :return
        Dictionary containing the peptide amino acids and the protein amino acids forming hydrogen bonds
        File with the predicted hydrogen bonds
        """

        # Read the protein structure
        model = self.reference[0]
        if self.path[-1] == "/":
            dssp = DSSP(model, self.pdb_file, dssp=self.path + "mkdssp")
        else:
            dssp = DSSP(model, self.pdb_file, dssp=self.path + "/mkdssp")

        # Loop over the keys from the dssp response to store ss and asa values
        total_index = {}
        list_chains = []
        reference_position = 0
        list_chains.append(list(dssp.keys())[0][0])
        self.hbonds_peptide = {}

        # iterate over the dssp keys to store the corresponding hydrogen bonds based on the atoms positions
        for keys in list(dssp.keys()):
            if keys[0] == list_chains[-1]:
                if keys[0] == self.chain:
                    res=self.positions[keys[1][1]]['aa']
                    total_index[reference_position + keys[1][1]] = (keys[0], res, keys[1][1])
                else:
                    total_index[reference_position + keys[1][1]] = (keys[0], dssp[keys][1], keys[1][1])
                last_position = reference_position + keys[1][1]
            else:
                list_chains.append(keys[0])
                reference_position = last_position
                if keys[0] == self.chain:
                    res=self.positions[keys[1][1]]['aa']
                    total_index[reference_position + keys[1][1]] = (keys[0], res, keys[1][1])
                else:
                    total_index[reference_position + keys[1][1]] = (keys[0], dssp[keys][1], keys[1][1])
                last_position = reference_position + keys[1][1]

            if keys[0] == self.chain:

                amino_name = self.positions[keys[1][1]]['aa'] + str(keys[1][1])
                if amino_name not in self.hbonds_peptide: self.hbonds_peptide[amino_name] = []
                if dssp[keys][6] < -5:
                    interactions_pos = last_position + dssp[keys][6]
                    if total_index[interactions_pos][0]!=self.chain:
                        self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][8] < -5:
                    interactions_pos = last_position + dssp[keys][8]
                    if total_index[interactions_pos][0] != self.chain:
                        self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][10] < -5:
                    interactions_pos = last_position + dssp[keys][10]
                    if total_index[interactions_pos][0] != self.chain:
                        self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][12] < -5:
                    interactions_pos = last_position + dssp[keys][12]
                    if total_index[interactions_pos][0] != self.chain:
                        self.hbonds_peptide[amino_name].append(total_index[interactions_pos])

        temporal_hbonds={}
        for pos_ref in self.positions:
            code=self.positions[pos_ref]['aa']+str(pos_ref)
            if code not in self.hbonds_peptide:
                self.hbonds_peptide[code]=[]
            temporal_hbonds[code]=self.hbonds_peptide[code]
            self.positions[pos_ref]['n_hbonds']=len(self.hbonds_peptide[code])

        self.hbonds_peptide=temporal_hbonds
        # Iterate over the peptide residues to show the hydrogen bonds
        self.number_hydrogen_bonds = 0
        if output:
            output_hydrogen_bonds = open("predicted_hbs_{}.txt".format(self.sequence), "w")
            print("These are the hydrogen bonds detected:")
            for residue in self.hbonds_peptide:
                for partners in self.hbonds_peptide[residue]:
                    print("{} interacts with residue {}{} from chain {}".format(residue, partners[1], partners[2],
                                                                                partners[0]))
                    output_hydrogen_bonds.write(
                        "{} interacts with residue {}{} from chain {}\n".format(residue, partners[1], partners[2],
                                                                                partners[0]))
                    self.number_hydrogen_bonds += 1
            output_hydrogen_bonds.close()

    ############################################################################

    def plot_hydrogen_bonds(self, type_layout="linear"):
        """
        Function to plot the hydrogen bonds using the igraph module of Python.
        NOTE: Installation instructions: https://igraph.org/python/

        :param type_layout: Layout used to plot the interactions

        :return PNG file with the graph of the hydrogen bonds
        """

        # Generate fragments
        fragment_chains = {}
        for amino in self.hbonds_peptide:
            for receptor in self.hbonds_peptide[amino]:
                chain = receptor[0]
                if chain not in fragment_chains: fragment_chains[chain] = []
                if receptor not in fragment_chains[chain]:
                    fragment_chains[chain].append(receptor)

        for ch in fragment_chains:
            fragment_chains[ch].sort(key=lambda x: x[2])

        # Get graph order and properties of the nodes and edges
        names = [0] * len(self.hbonds_peptide)
        pep_interactions = []
        chain_nodes = []
        chain_colors = {"PEP": "yellow"}
        edge_colors = {"5": "black", "4": "orange", "8": "orange"}
        list_colors = ["cyan", "green", "magenta", "blue"]
        type_interactions = []

        # Iterate over the predicted hydrogen bonds
        for i, amino in enumerate(self.hbonds_peptide):
            names[i] = amino
            chain_nodes.append("PEP")

            if i != len(self.hbonds_peptide) - 1:
                pep_interactions.append((i, i + 1))
                type_interactions.append(5)

        # Iterate over the fragments to assign the properties and types of interactions
        reference = len(self.hbonds_peptide)
        for num_chain, chains in enumerate(fragment_chains):
            chain_colors[chains] = list_colors[num_chain]
            for count, residues in enumerate(fragment_chains[chains]):
                chain_nodes.append(chains)
                names.append(residues[0] + "-" + residues[1] + str(residues[2]))
                number_res = reference + count

                for pep_amino in self.hbonds_peptide:
                    if residues in self.hbonds_peptide[pep_amino]:
                        pep_interactions.append((names.index(pep_amino), number_res))
                        type_interactions.append(self.hbonds_peptide[pep_amino].count(residues) * 4)

            reference = len(names)

        # Create graph
        g = Graph(pep_interactions)
        g.vs["name"] = names
        g.vs["chains"] = chain_nodes
        g.es["count_int"] = type_interactions
        color_per_chain = [chain_colors[chColor] for chColor in g.vs["chains"]]
        color_per_edge = [edge_colors[str(edColor)] for edColor in type_interactions]

        # Use a default layout
        if type_layout == "linear":
            layout = []
            ref_chain = chain_nodes[0]
            positions = [0, -1, 1]
            c_pos = 0
            c_nod = 0
            for c, ele in enumerate(chain_nodes):
                if ele == ref_chain:
                    layout.append((c_nod, positions[c_pos]))
                    c_nod += 1
                else:
                    ref_chain = ele
                    c_pos += 1

                    number_ch_ele = chain_nodes.count(ele)
                    c_nod = (len(self.hbonds_peptide) / 2) - (number_ch_ele / 2)
                    layout.append((c_nod, positions[c_pos]))
                    c_nod += 1
            bbox = (1200, 400)
        if type_layout == "cyclic":
            layout = g.layout_fruchterman_reingold()
            bbox = (800, 800)

        visual_style = {}
        visual_style["vertex_size"] = 60
        visual_style["vertex_color"] = color_per_chain
        visual_style["edge_color"] = color_per_edge
        visual_style["vertex_label"] = g.vs["name"]
        visual_style["edge_width"] = type_interactions
        visual_style["layout"] = layout
        visual_style["bbox"] = bbox
        visual_style["margin"] = 80
        plot(g, "plot_hbs_{}.png".format(self.sequence), **visual_style)

    ############################################################################

    def get_heavy_atom_contacts(self, contact_threshold):
        """
        Function to count heavy atom contacts per residue

        :param contact_threshold: threshold to define when a contact is created

        :return
        Add to the amino acids dictionary the number of contacts per amino acid
        total_contacts -- total number of contacts calculated for the full peptide
        """

        # Get the other chains present in the PDB file
        chains = []
        for ch in self.reference[0]:
            chLetter = ch.get_id()
            chains.append(chLetter)

        # Store the count per amino acid in the peptide
        countDict = {}

        # Loop over the chains
        for chValues in chains:
            if chValues != self.chain:
                for residue in self.reference[0][self.chain]:
                    # Store the distances and a counter
                    counter = 0
                    distances = []
                    for residueA in self.reference[0][chValues]:
                        for atomA in residueA:
                            idAtomA = atomA.get_id()
                            if idAtomA[0] != "H" and idAtomA[0].isdigit() == False:
                                for atom in residue:
                                    idAtom = atom.get_id()
                                    if idAtom[0] != "H" and idAtom[0].isdigit() == False:
                                        # Get the distance differences
                                        diff = atom.coord - atomA.coord
                                        diffValue = np.sqrt(np.sum(diff * diff))
                                        distances.append(diffValue)
                    # Count the number of contacts based on the defined threshold
                    for d in distances:
                        if d < contact_threshold: counter += 1

                    # Get the information per residue
                    res = residue.get_resname()
                    resPos = str(residue.get_full_id()[3][1])
                    if res + resPos not in countDict:
                        countDict[res + resPos] = counter
                    else:
                        countDict[res + resPos] += counter

        # Calculate the final contacts per aa and in total
        self.total_contacts = 0
        for position in self.positions:

            residue = self.positions[position]["full_aa"]
            if residue + str(position) in countDict:
                self.positions[position]["contacts"] = countDict[residue + str(position)]
                self.total_contacts += countDict[residue + str(position)]

    # End of Interactions class
    ############################################################

############################################################
## End of interactions.py
############################################################

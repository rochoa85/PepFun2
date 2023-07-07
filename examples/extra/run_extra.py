#!/usr/bin/python
'''
Example: extra module
Input: Sequence of peptides (including NNAAs) with BILN format
Output: List of descriptors for each peptide sequence
'''

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

# PepFun modules
from pepfun.extra import Alignment
from pepfun.extra import pepDescriptors

similarity = Alignment.similarity_pair("KMGRLFR", "KAGRSFR")
print("The similarity between peptides KMGRLFR and KAGRSFR is: {}".format(similarity))

score,start,end=Alignment.score_modified("K-Aib-M-P","S-Aib-Iva-P",mode="weighted")
print("The score between peptides K-Aib-M-P and S-Aib-Iva-P is: {}".format(score))

# Report with the summary per peptide position, including accesible surface area values
report = open("descriptors.txt", "w")

# List of sequences
sequences = ["K-Aib-Iva-P-L-C-D", "S-D-Orn-A-F-R-G-Nva-L-L-D", "I-T-F-E-D-Nle-L-D-Pyr-Y-G"]

for i,sequence in enumerate(sequences):
    desc = pepDescriptors(sequence, generate_properties=False)
    descriptors = desc.moranCorrelation()
    # Create the header
    if i==0:
        header="Sequence"
        for key in descriptors:
            header+=f" {key}"
        report.write(header+"\n")
    # Write the descriptors for the sequences
    line=sequence
    for key in descriptors:
        line+=f" {descriptors[key]}"
    report.write(line+"\n")

report.close()

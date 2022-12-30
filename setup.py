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

import setuptools

########################################################################################
# Main
########################################################################################

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name='pepfun',
    version='2.0',
    url='',
    license='MIT',
    author='Rodrigo Ochoa',
    author_email='rodrigo.ochoa@boehringer-ingelheim.com',
    description='Tools for analysis with peptides',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy >= 1.22.2',
        'pandas >= 1.4.1',
        'requests >= 2.27.1',
        'openpyxl >= 3.0.9',
        'tqdm >= 4.62.3',
        'biopython >= 1.79',
        'igraph >= 0.9.10',
        'sphinx-jsonschema'
    ],
    include_package_data=True,
    packages=['pepfun'],
    package_data={'pepfun': ['data/*']}
)

try:
   import rdkit
except ModuleNotFoundError:
   raise ValueError("No rdkit installation found! Make sure you install it first!\n"
                    "Try: 'conda install -c rdkit rdkit'.")

if rdkit.__version__ < '2019.09.04':
   raise ValueError(f"rdkit version must be at least 2019.09.04, but is {rdkit.__version__}.\n"
                    "Please update first using 'conda update rdkit'")


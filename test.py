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

from unittest import TestLoader, TestResult
from pathlib import Path

##########################################################################
# Function
##########################################################################

def run_tests():
    """
    Function to run all the available unittests
    """

    test_loader = TestLoader()
    test_result = TestResult()

    test_directory = str(Path(__file__).resolve().parent / 'tests')

    test_suite = test_loader.discover(test_directory, pattern='*_Test.py')
    test_suite.run(result=test_result)

    if test_result.wasSuccessful():
        print("\n################################\nAll of the tests were successful!\n################################")
        exit(0)
    else:
        print("\nFailures:")
        for i,fail in enumerate(test_result.failures):
            print(f"{i+1}. {fail[0]}\n{fail[1]}")
        if test_result.errors:
            print("\nErrors:")
        exit(-1)

if __name__ == '__main__':
    run_tests()
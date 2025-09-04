
import sys
import os
import argparse

#_______________________________________________________________________________

def parse_arguments():

    """
    Returns:
    - args

    Structure:
    main/parse_arguments

    Uses:
    argparse

    """

    parser = argparse.ArgumentParser(description="Calculate xas, xes, or analyse RIXS spectra")


    parser.add_argument("--calculate", choices=["xas", "xes", "rixs"], required=True,
                        help=" Choose which type of spectrum to calculate ")


    parser.add_argument("--oscstr", choices=["ipa", "bse", "cumm", "final", "from_nc"],
                        help="""
                            ---> For XAS:
                                - from_nc: reads the oscillator strength from nc file,
                                - bse: calculates it from A_l and pmat,
                                - ipa: it calculates it using only pmat
                            ---> For XES
                                - ipa will be set automatically
                            ---> For RIXS:
                                - from_nc : reads the oscillator strength from rixs.h5
                                - cumm    : calculates oscstr[lo,lc], for specific lo, to do RIXS analysis
                                - final   : calculates oscstr[lo] for all lo, to obtain RIXS cross section
                            """)


    parser.add_argument("--inputs", nargs='+', action='append', required=True,
                        help=""" Specify a code followed by its required files.
                            Example: --inputs exciting eigvec.h5 pmat.h5 eigval.txt
                                     --inputs dp exceig.nc rhotw.nc eigval.txt """)


    parser.add_argument("--win_list", nargs='+', type=float,
                        help="""" Incident photon energy (in eV) for RIXS analysis
                            Example: --win_list 173.0 173.2 173.4 173.6 173.8  """)


    parser.add_argument("--lo_list", nargs='+', type=int,
                        help=""" List of lambda_o indices corresponding to the
                            to loss energies that we want to look in the RIXS analysis
                            Example: --lo_list 2 25 53 115 """)


    return parser.parse_args()

#_______________________________________________________________________________




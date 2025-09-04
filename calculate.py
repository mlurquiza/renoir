
import numpy as np
import netCDF4 as nc
import h5py
import sys

sys.path.append(".")

from CalculateEps  import CalculateEps
from CalculateRIXS import CalculateRIXS

#_______________________________________________________________________________

def run_calculate_eps (files, code, oscstr, xes, vol, broad, omegas, points, pol):

    """
    This function computes XAS or XES spectra.

    Parameters:
    - files  : list[str] --> input file paths for eigenvectors, matrix elements and eigenvalues
    - code   : str --> "dp" or "exciting"
    - oscstr : str --> "ipa" or "bse" oscillator strength
    - xes    : bolean --> True if calculating XES
    - vol    : float --> unit cell volume
    - broad  : complex --> broadening
    - omegas : np.ndarray --> array of energy values
    - points : int --> number of points
    - pol    : int --> polarization direction

    Returns:
    - epsM   :  np.ndarray [omegas]

    Structure:
    main/run_calculate_eps

    Uses:
    CalculateEps class

    """

    nc_file   = nc.Dataset(files[0])
    pmat_file = nc.Dataset(files[1])
    eigval_file = open(files[2], 'r')

    spectrum = CalculateEps(nc_file, pmat_file, eigval_file, code, xes, oscstr) # creates an instance within the class

    if xes:
        print("Running calculation of XES for", code, "...")
        print()
        epsM = spectrum.calc_emission(vol, broad, omegas, points, pol)       # calculates the emission spectra
    else:
        print("Running calculation of XAS for", code, "...")
        print()
        epsM = spectrum.calc_eps(vol, broad, omegas, points, pol)            # calculates the absorption spectra
        #epsM = spectrum.calc_eps_trans(vol, broad, omegas, points, pol)     # calculates the spectra in the transveral gauge
        eigval_file.close()

    eigval_file.close()
    nc_file.close()
    pmat_file.close()

    return epsM


#_______________________________________________________________________________

def run_calculate_RIXS (files, outfile, code, oscstr, eta_o, eta_c, lo_list, win_list, wloss_list):

    """
    This function calculates and analyzes RIXS spectra.

    Parameters:
    - files    : list[str] --> NetCDF files: rixs.h5 and data.h5
    - oscstr   : str --> "cumm", "final" or 'from_nc'
    - eta_o    : complex --> broadening for final state (e_loss axis)
    - eta_c    : complex --> core hole lifetime (w_in axis)
    - w_in     : float or tuple/list -->  incident energies (eV)
    - lo_list  : list[float] --> lambda_o indices for analysis

    Structure:
    main/run_calculate_RIXS

    Uses:
    CalculateRIXS class

    """

    nc_rixs = h5py.File(files[0])
    nc_data = h5py.File(files[1])

    rixs = np.zeros((len(win_list), len(wloss_list)), dtype=complex)

    print('Running calculation of RIXS for =', code)
    print()


    for win_idx, win in enumerate(win_list):

        analysis = CalculateRIXS(nc_rixs, nc_data, outfile, oscstr, eta_o, eta_c, win, win_idx, lo_list)

        if oscstr == 'final':
            rixs[win_idx] = analysis.calc_rixs(wloss_list)
            analysis.write_oscstr()


        elif oscstr == 'from_nc':
            rixs[win_idx] = analysis.calc_rixs(wloss_list)
            analysis.write_rixs(wloss_list)


    nc_rixs.close()
    nc_data.close()


    return rixs








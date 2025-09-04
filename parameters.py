
import numpy as np

#___________________Function_______________________________________________________


def get_param_eps():

    """
    Returns:
    vol, broad, pol, omegas, points

    Structure:
    main/get_param_eps

    """

    vol = 246.320618 #hBN=246.320618 #al2o3=573.141297 #silicon=266.6307
    broad = 1j*0.00367493
    pol = 0   # x=0, y=1, z=2 Must agree with the polarization used for fullchi

    ### Plot spectrum in range (Ha)
    emin = 6.00
    emax = 8.00
    points = 3001#2001
    step   = (emax-emin)/points
    omegas = np.arange(emin,emax,step)
    #TODO sicsor

    print('within a range of energies:', emin,'to', emax,'Ha,', 'with an energy step of:', step, 'Ha')


    return vol, broad, pol, omegas, points


def get_param_rixs():

    """
    Returns:
    eta_o, eta_c, omegas, points

    Structure:
    main/get_param_rixs
    """

    eta_o = 1j*0.2 #eV --> 1j*0.00367493
    eta_c = 1j*0.05 #eV --> 1j*0.00367493

    ### Calculate RIXS spectrum in range (eV)
    emin = 0.00   # 2.0 #4.0 en eV
    emax = 16.00  # 3.5 #4.3 en eV
    points = 8000 # 3500
    step   = (emax-emin)/points
    omegas = np.arange(emin,emax,step)
    #TODO sicsor


    return eta_o, eta_c, omegas, points



################################################################
#---------------- Class CalculateRIXS --------------------------
################################################################

import numpy as np
import netCDF4 as nc
import h5py

class CalculateRIXS:

    """
    Class to calculate or analyze RIXS spectra using output from BRIXS.

    List of input variables:
     self.nc_rixs   -----> HDF5 file with RIXS oscillator strengths and eigenvalues
     self.nc_data   -----> HDF5 file with t1 and t2 arrays
     self.osc       -----> mode of oscillator strength evaluation: 'cumm', 'final', or 'from_nc'
     self.eta_o     -----> broadening parameter for final optical excitation [eV]
     self.eta_c     -----> broadening parameter for intermediate (core) states [eV]
     self.win       -----> incident photon energy [eV]
     self.win_idx   -----> index for the incident photon energy
     self.lo_list   -----> list of optical states indices used in 'cumm' mode

    List of variables from files:
     self.evals_c   dim(nl_c)       -----> core exciton eigenvalues [Ha]
     self.evals_o   dim(nl_o)       -----> optical exciton eigenvalues [Ha]
     self.nl_c      int             -----> number of core excitons
     self.nl_o      int             -----> number of optical excitons
     self.t1        dim(nl_c)       -----> absorption oscillator strength
     self.t2        dim(nl_c, nl_o) -----> emission oscillator strength
     self.t3        dim(nl_o)       -----> final RIXS oscillator strength

    Calculated quantities
     self.t3_cumm   dim(nl_c, nl_o) -----> cumulative RIXS oscillator strength per optical state
     self.rixs      dim(wloss)      -----> computed RIXS spectrum vs energy loss

    """

    HA2EV = 27.211396641308

    def __init__(self, nc_rixs, nc_data, outfile, osc, eta_o, eta_c, win, win_idx, lo_list):

        self.nc_rixs = nc_rixs
        self.nc_data = nc_data
        self.outfile = outfile
        self.osc = osc
        self.eta_o = eta_o
        self.eta_c = eta_c
        self.win = win
        self.win_idx = win_idx

        self.evals_c = self.nc_rixs['/cevals'][:]
        self.evals_o = self.nc_rixs['/vevals'][:]

        self.nl_c = len(self.evals_c)

        if self.osc == 'cumm':
            self.lo_list = lo_list
            self.nl_o = len(self.lo_list)

        elif self.osc == 'final':
            self.nl_o = len(self.evals_o)
            self.lo_list = list(range(self.nl_o))

        elif self.osc == 'from_nc':
            self.nl_o = len(self.evals_o)
            self.lo_list = list(range(self.nl_o))

        else:
            raise ValueError(f"Unknown osc mode '{self.osc}'. Use 'cumm', 'final', or 'from_nc'.")


        # Absorption oscstr
        self.t1 = None

        # Emission oscstr
        self.t2 = None

        # RIXS oscstr
        if self.osc == 'cumm':         ### calculate cummulative oscstr as a function of lambda_core, only for few lambda_optical

            self.cumm_t3 = np.zeros((self.nl_c, self.nl_o), dtype=complex)

            for idx, lo in enumerate(self.lo_list):
                self.__get_t1__()

                self.__get_t2__(lo)

                self.__cummul__(idx)

                self.__write_t1t2t3__(idx, lo)


        elif self.osc == 'final':      ### calculate the final oscstr for every lambda_optical

            self.t3 = np.zeros(self.nl_o, dtype=complex)

            for lo in self.lo_list:
                self.__get_t1__()

                self.__get_t2__(lo)

                self.t3 [lo] = self.__vecmul__()


        elif self.osc == 'from_nc':     ### read the oscstr as a function of lambda_optical for each win from rixs.h5

            self.t3 = np.zeros(self.nl_o, dtype=complex)

            self.__get_t3__()


        # RIXS spectrum
        self.rixs = None

#_______________________________________________________________________________
#
#   Internal Methods
#_______________________________________________________________________________


    def __get_t2__(self, lo):
        self.t2 = self.nc_data['/t(2)'][:,lo,0] + 1j * self.nc_data['/t(2)'][:,lo,1]


    def __get_t1__(self):
        self.t1 = self.nc_data['/t(1)'][:, 0] + 1j * self.nc_data['/t(1)'][:, 1]


    def __get_t3__(self):
        dataset_name = f"{self.win_idx + 1:04d}"  # Obtain group name: '0001', '0002', etc
        if dataset_name not in self.nc_rixs['oscstr']:
            raise KeyError(f"Missing dataset: {dataset_name} in 'oscstr'")

        self.t3 [:] = self.nc_rixs['oscstr'][dataset_name][:,0] + 1j * self.nc_rixs['oscstr'][dataset_name][:,1]


    def __vecmul__(self):
        tmp = 0.00 + 1j*0.00
        for lc in range(self.nl_c):
            #tmp+= self.t1[lc] * self.t2[lc] / (self.win - self.evals_c[lc] * self.HA2EV + self.eta_c)
            tmp += self.t1[lc] * self.t2[lc] / (self.win - self.evals_c[lc] * self.HA2EV - self.eta_c)

        return tmp


    def __cummul__(self, idx):
        tmp = 0.00 + 1j*0.00
        for lc in range(self.nl_c):
            #tmp+= self.t1[lc] * self.t2[lc] / (self.win - self.evals_c[lc] * self.HA2EV + self.eta_c)
            tmp += self.t1[lc] * self.t2[lc] / (self.win - self.evals_c[lc] * self.HA2EV - self.eta_c)
            self.cumm_t3 [lc, idx] = tmp

        return tmp


#_______________________________________________________________________________
#
#   Calculate the RIXS cross section
#_______________________________________________________________________________


    def calc_rixs(self, wloss_list):

        self.rixs = np.zeros((len(wloss_list)), dtype=complex)

        self.t3 = -1.0 * np.abs(self.t3)**2

        for idx, wloss in enumerate(wloss_list):

            den = 1 / (wloss - self.evals_o * self.HA2EV + self.eta_o)

            self.rixs[idx] = np.dot(self.t3, den)

        return self.rixs



#_______________________________________________________________________________
#
#   Write
#_______________________________________________________________________________


    def __write_t1t2t3__(self, idx, lo):
        self.outfile.write("\n")
        self.outfile.write("{:<s}\t{:>15.9f}\n".format('### win = ', self.win))
        self.outfile.write("{:<s}\t{:>15.9f}\n".format('### lambda_o = ', self.evals_o[lo]))
        self.outfile.write("{:<s}\t{:>10d}\n".format('### n_lambda_c = ', self.nl_c))
        self.outfile.write("{:<s}\t{:<s}\t{:<s}\t{:<s}\n".format('### lambda_c','t1','t2','t3'))

        for lc in range(self.nl_c):
            self.outfile.write("{:e}\t{:e}\t{:e}\t{:e}\n".format(
                self.evals_c[lc] * self.HA2EV,
                abs(self.t1[lc]),
                abs(self.t2[lc]),
                abs(self.cumm_t3[lc,idx])
            ))



    def write_oscstr(self):
        self.outfile.write("\n")
        self.outfile.write("{:<s}\t{:>15.9f}\n".format('### win = ', self.win))
        self.outfile.write("{:<s}\t{:<s}\t{:<s}\t{:<s}\n".format('### lambda_o','Re[t3]','Im[t3]', 'Abs[t3]'))

        for idx, lo in enumerate(self.lo_list):
            self.outfile.write("{:> 4.6f} {:> 4.6f} {:> 4.6f} {:> 4.6f}\n".format(
                    self.evals_o[lo] * self.HA2EV,
                    self.t3[lo].real,
                    self.t3[lo].imag,
                    abs(self.t3[lo])
            ))



    def write_rixs(self, omegas):
        self.outfile.write("\n")
        self.outfile.write("{:<s} {:<s}\n".format('### win =', str(self.win)))
        self.outfile.write("{:<s}\t{:<s}\t{:<s}\t{:<s}\n".format('### omega', 'Re[rixs]', 'Im[rixs]', 'Abs[rixs]'))

        for idx, w in enumerate(omegas):
            self.outfile.write("{:> 4.6f} {:> 4.6f}\n".format(
                    w,
                    self.rixs[idx].imag/max(self.rixs.imag)
            ))





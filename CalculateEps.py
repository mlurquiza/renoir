################################################################
#--------------------- Class CreateEps -------------------------
################################################################

import numpy as np
import netCDF4 as nc

class CalculateEps:

    """
    Class to compute absorption (XAS) or emission (XES) spectra from BSE or IPA data.

    List of input variables:
     self.code  -----> exciting or dp
     self.osc   -----> from_scratch, ipa or bse (with or without local fields)
     self.xes   -----> whether to compute emission (XES) or absorption (XAS)

    List of variables from HDF5/netcdf:
     self.pmatvec  dim(ntrans,3) ----------> dipole matrix elements ordered by smap
     self.oscstrr  dim(nexcitons,3) -------> numerator of the spectrum
     self.evals    dim(nexcitons) ---------> BSE eigenvalues in Ha, ordered by increasing energy
     self.evalsIP  dim(nexcitons) ---------> IP transition energies in Ha, ordered by increasing energy
     self.ensort   dim(nexcitons) ---------> index list ordering evalsIP by increasing energy from vck with order smap
     self.rvec     dim(ntrans,nexcitons) --> BSE eigenvectors ordered axis 0 by smap, axis 1 by increasing energy

    List of variables from energy file:
     self.eigval  dim(number of transitions written down)
     self.evec    dim(ntrans)

    """

    
    def __init__(self, nc_file, pmat_file, eigval_file, code, xes, osc):
        self.code = code 
        self.xes = xes
        self.osc = osc
        self.nkp = None
        self.nexcitons = None
        self.evals = None
        self.evalsIP = None
        self.ensort = None
        self.koulims = None
        self.smap = None
        self.ntrans = None
        self.rvec = None 
        self.pmatvec = None
        self.evec = None
        self.oscstrr = None

        self.get_info(nc_file, pmat_file, eigval_file)

        if  self.xes == False:
            self.calc_osc(nc_file)


    def get_info(self, nc_file, pmat_file, eigval_file):

        self.__get_param__(nc_file)
        self.__get_trans__(nc_file)

        if self.osc == 'bse':
            self.__get_rvec__(nc_file)

        if self.xes == True:
            self.__get_evec__(eigval_file)

        self.__get_pmatvec__(pmat_file)


    def __get_param__(self, nc_file):

        self.evals      = nc_file['/eigvec-singlet-TDA-BAR-full/0001/evals'][:]
        self.evalsIP    = nc_file['/eigvec-singlet-TDA-BAR-full/0001/evalsIP'][:]
        self.ensort     = nc_file['/eigvec-singlet-TDA-BAR-full/0001/parameters/ensortidx'][:]
        self.nexcitons  = nc_file['/eigvec-singlet-TDA-BAR-full/0001/parameters/nexcstored'][:][0]
        self.koulims    = nc_file['/eigvec-singlet-TDA-BAR-full/0001/parameters/koulims'][:]
        self.nkp        = self.koulims.shape[0]
        self.smap       = nc_file['/eigvec-singlet-TDA-BAR-full/0001/parameters/smap'][:]



    def __get_trans__(self, nc_file):

        self.ci  = np.zeros(self.nkp, dtype=int)
        self.cf  = np.zeros(self.nkp, dtype=int)
        self.nc  = np.zeros(self.nkp, dtype=int)
        self.ui  = np.zeros(self.nkp, dtype=int)
        self.uf  = np.zeros(self.nkp, dtype=int)
        self.nu  = np.zeros(self.nkp, dtype=int)
        self.vi  = np.zeros(self.nkp, dtype=int)
        self.vf  = np.zeros(self.nkp, dtype=int)
        self.nv  = np.zeros(self.nkp, dtype=int)
        for k in range(self.nkp):
            self.ui[k] = self.koulims[k,2] - 1
            self.uf[k] = self.koulims[k,3] - 1
            self.nu[k] = self.koulims[k,3] - self.koulims[k,2] + 1
            self.ci[k] = self.koulims[k,0] - 1
            self.cf[k] = self.koulims[k,1] - 1
            self.nc[k] = self.koulims[k,1] - self.koulims[k,0] + 1
            if self.code == 'dp':   # for dp we have to remove the core states, to start the counting from the valence states
                self.ci[k]  = self.ci[k] - self.nu[k]
                self.cf[k]  = self.cf[k] - self.nu[k]
            self.nv[k] = self.ci[k] - self.vf[k]
            self.vi[k] = 0
            self.vf[k] = self.nv[k] - 1
            #print(k, self.ui[k], self.uf[k], self.nu[k])
            #print(k, self.ci[k], self.cf[k], self.nc[k])


        if self.xes == False:
            self.ntrans = self.smap.shape[0]
            if self.osc == "bse" and self.nexcitons < self.ntrans:
                print('WARNING: not enough excitons written down')

        elif self.xes == True:
            self.ntrans = self.nv[0] * self.nu[0] * self.nkp



    def __get_rvec__(self, nc_file):

        self.rvec = np.zeros((self.ntrans, self.nexcitons), dtype=complex)
        for l in range(self.nexcitons):
            name = '{:0>8d}'.format(l+1)   # make the total string size of 8, fill with zeros to the left
            self.rvec[:,l] =    nc_file['/eigvec-singlet-TDA-BAR-full/0001/rvec/'+name][:,0] + \
                             1j*nc_file['/eigvec-singlet-TDA-BAR-full/0001/rvec/'+name][:,1]



    def __get_evec__(self, eigval_file):

        line = eigval_file.readline()
        line = eigval_file.readline()
        tok  = line.strip().split()
        nlines = int(tok[0])
        eigval = np.zeros((self.nkp, nlines))
        for k in range(self.nkp):
            line = eigval_file.readline()
            line = eigval_file.readline()
            line = eigval_file.readline()
            line = eigval_file.readline()
            for it in range(nlines):
                line = eigval_file.readline()
                tok = line.split()
                eigval[k,it] = float(tok[1])

        self.evec = np.zeros ((self.ntrans))

        idx = 0
        for k in range(self.nkp):                     #loop over kpt
            for v in range(self.vi[k], self.vf[k]+1):     #loop over valence
                for u in range(self.ui[k], self.uf[k]+1):     #loop over core
                    #should be:      E_xes[idx] = E_valence[k,v] - E_core[k,u]
                    #for the moment: E_xes[idx] = E_valence[k,v] - (E_core + E_lumo) = E_valence[k,v] + E_xas[0]
                    self.evec[idx] = eigval[k,v] + 6.104581592215 #self.evalsIP[0]  # idx ordered by ascending GW energy
                    idx += 1



    def __get_pmatvec__(self, pmat_file):

        if self.code == 'exciting':
            print('Reading pmat of exciting')
        elif self.code == 'dp':
            print('Reading pmat of dp')

        n1 = pmat_file['/pmat/00000001/pmat'][:,:,:,0].shape[0]   # valence(v) + conduction(c)
        n2 = pmat_file['/pmat/00000001/pmat'][:,:,:,0].shape[1]   # core(u)
        pmat = np.zeros((self.nkp, n1, n2, 3), dtype=complex)
        for k in range(self.nkp):
            name = '{:0>8d}'.format(k+1)
            pmat[k,:,:,:] = pmat_file['/pmat/'+name+'/pmat'][:,:,:,0] + \
                           1j*pmat_file['/pmat/'+name+'/pmat'][:,:,:,1]
        print('pmat shape:',pmat.shape)


        self.pmatvec = np.zeros((self.ntrans, 3), dtype=complex)


        # For core --> valence transitions
        if self.xes == True:
            it = 0
            for k in range(self.nkp):                     #loop over kpts
                for v in range(self.vi[k], self.vf[k]+1):    #loop over valence
                    for u in range(self.ui[k], self.uf[k]+1):  #loop over core
                        self.pmatvec[it,:] = pmat[k,v,u,:]      #loop over 3 coordinates xyz
                        it += 1


        # For core --> conduction transitions
        elif self.xes == False:
            it = 0
            for k in range(self.nkp):                      # loop over kpts
                for u in range(self.ui[k], self.uf[k]+1):     # loop over core
                    for c in range(self.ci[k], self.cf[k]+1):   #loop over conduction
                        self.pmatvec[it,:] = pmat[k,c,u,:]         #loop over 3 coordinates xyz
                        it += 1


#_______________________________________________________________________________
#
#   Calculate the oscillator strength or get it from .nc file
#_______________________________________________________________________________

    def calc_osc(self, nc_file):

        self.oscstrr = np.zeros((self.ntrans, 3), dtype=complex)

        #read oscillator strengh from a nc file (only for exciting)
        if  self.osc == "from_nc":
            for i in range(3):
              group = '{:0>1d}'.format(i+1)
              self.oscstrr[:,i] = nc_file['/excitons-singlet-TDA-BAR-full/0001/'+group+'/oscstrr'][:,0] + \
                               1j*nc_file['/excitons-singlet-TDA-BAR-full/0001/'+group+'/oscstrr'][:,1]

        #calculate bse or rpa oscillator strengh
        elif self.osc == "bse":
            for i in range(3):
              for l in range(self.nexcitons): # l are already ordered by ascending E_lambda
                temp1 = 0+0*1j
                temp2 = 0+0*1j
                for it in range(self.ntrans): # oedered by transition vck (smap)
                  temp1 += np.conjugate(self.pmatvec[it,i]) * self.rvec[it,l]
                  temp2 += self.pmatvec[it,i] * np.conjugate(self.rvec[it,l])
                self.oscstrr[l,i] = temp1 * temp2



        #calculate ipa oscillator strengh, is the default option
        elif self.osc == "ipa":
            self.evals = self.evalsIP   # overwrite energies to use IP ones in the denominator of eps
            for i in range(3):
              for idx in range(self.ntrans):   # it will run over transitions with ascending GW energy
                it = self.ensort[idx] - 1              # let's find the corresponding index in the basis of transitions vck
                self.oscstrr[idx,i] = self.pmatvec[it,i]* np.conjugate(self.pmatvec[it,i])   # oscstrr[idx] is sorted by ascending GW energy



#_______________________________________________________________________________
#
#  Calculate epsilon macroscopic
#_______________________________________________________________________________

    def calc_eps(self, vol, eta, omegas, points, pol):

        prefac = 8*np.pi/(self.nkp*vol)    # 1/q has been ommited
        epsM = np.zeros((points,), dtype=complex)
        for idx, w in enumerate(omegas):
            for l in range(self.nexcitons):
                epsM[idx] += self.oscstrr[l,pol] * \
                             (1/(self.evals[l]-w-eta) + 1/(self.evals[l]+w+eta))

        epsM = 1 + prefac * epsM
        if self.code == 'dp':
            epsM = epsM/self.nu[0]

        return(epsM)



    def calc_eps_trans(self, vol, eta, omegas, points, pol):

        prefac = 8*np.pi/(self.nkp*vol)
        epsM = np.zeros((points,), dtype=complex)
        for idx, w in enumerate(omegas):
            for l in range(self.nexcitons):
                epsM[idx] += self.oscstrr[l,pol] * \
                            (1/(self.evals[l]-w-eta) + 1/(self.evals[l]+w+eta)) * (1/w)**2

        epsM = 1 + prefac * epsM
        if self.code == 'dp':
            epsM = epsM/self.nu[0]

        return(epsM)



    def calc_emission(self, vol, eta, omegas, points, pol):

        prefac = 8*np.pi/(self.nkp*vol)    # 1/q has been ommited
        epsM = np.zeros((points,), dtype=complex)
        for idx, w in enumerate(omegas):
            for l in range(self.ntrans):
               epsM[idx] += (np.conjugate(self.pmatvec[l,pol]) * self.pmatvec[l,pol]) * \
                            (1/(self.evec[l]-w-eta) + 1/(self.evec[l]+w+eta))

        epsM = 1 + prefac * epsM

        if self.code == 'dp':
            epsM = epsM/self.nu[0]

        return(epsM)


  

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from copy import copy
from scipy import interpolate
from tqdm import tqdm
from scipy.interpolate import interp1d
from parameters import *


def vector_dot(a,b):
    """
    dot product of two vector list
    """
    if np.ndim(a)>1:
        temp = a*b 
        result = np.sum(temp,axis = 1)
    elif np.ndim(a)==1:
        result = np.dot(a,b)
    return result
def generate_1_band_model(paras):
    t1, t2, t3, tz = [paras.t1, paras.t2, paras.t3, paras.tz]
    norm = paras.norm
    c = paras.c_latt
    kz = np.pi/c
    def one_band_model(kx, ky):
        term1 = -2*t1*(np.cos(norm*kx)+ np.cos(norm*ky)) 
        term2 = -4*t2*np.cos(norm*kx)*np.cos(norm*ky)
        term3 = -2*t3*(np.cos(2*norm*kx)+np.cos(2*norm*ky))
        term4 = -2*tz*np.cos(norm*kx/2)*np.cos(norm*ky/2)*np.cos(norm*kz/2)*(np.cos(norm*kx)-np.cos(norm*ky))**2
        return term1+term2+term3+term4
    
    return one_band_model

def generate_3_band_model(paras, num_intp = 1000):
    """
    calculate the energy band (with mu = 0!!) and orbital contribution by directly diagonalizing the Hamiltonian, and Interpolate afterwards (to increase the efficiency)
    """
    def meshgrid2array(*mesh, type = 2):
        """
        convert a meshgrid element to an array. We assume the grids are retangular.
        
        key arguements:
        -----------------------------------------
        - (typexNxN array) | mesh: input mesh grids  
        - (int)       | type: if 1, only input 1 array; if 2, input 2 arraies.
        return:
        - (N^2x1 array) | return reshaped array(s)
        """
        if type == 2:
            X,Y = mesh
            numx = np.shape(X)[0]
            numy = np.shape(X)[1]
            X_temp = X.reshape((1,numx*numy))[0]
            Y_temp = Y.reshape((1,numx*numy))[0]
            output = np.array(list(zip(X_temp, Y_temp)))
            return output
        elif type == 1:
            X = mesh[0]
            numx = np.shape(X)[0]
            numy = np.shape(X)[1]
            X_temp = X.reshape((1,numx*numy))[0]
            return X_temp  
        else:
            raise TypeError('type should either be 1 or 2')  
    norm = paras.norm
    mu = 0
    def hamiltonian(k):
        """
        Hamiltonian for a 3-bands model of CuO2
        """
        def sx(k):
            if np.ndim(k)>1:
                kx = k[:,0]
                return np.sin(kx*norm/2) ## CHECK AGAIN ##
            else:
                return np.sin(k[0]*norm/2)
        def sy(k):
            if np.ndim(k)>1:
                ky = k[:,1]
                return np.sin(ky*norm/2)
            else:
                return np.sin(k[1]*norm/2)
        t_pd, ti_pp, td_pp = [1600.0, -1000.0, 0.0] # in the unit of eV
        t_pp = ti_pp + td_pp
        energy_d, energy_p = [0.0, -900.0] # in the unit of eV

        H11 = energy_d-mu
        H12 = 2*t_pd*sx(k)
        H13 = -2*t_pd*sy(k)
        H22 = energy_p-mu + 4*ti_pp*sx(k)**2
        H23 = 4 * t_pp*sx(k)*sy(k)
        H33 = energy_p-mu + 4*ti_pp*sy(k)**2
        H = np.array([[H11,H12,H13],[H12,H22,H23],[H13,H23,H33]])
        return H

    def hamiltonian_mesh(k1_mesh,k2_mesh):
        k_array = meshgrid2array(k1_mesh,k2_mesh)
        H_array = np.array(list(map(hamiltonian,k_array)))
        H_mesh = H_array.reshape((np.shape(k1_mesh)[0],np.shape(k1_mesh)[1],3,3))
        return H_mesh

    def energy_hamiltonian(k1_mesh,k2_mesh,band = 2):
        H_mesh = hamiltonian_mesh(k1_mesh,k2_mesh)
        eigs, eigvs = np.linalg.eigh(H_mesh)
        return eigs[:,:,band], np.abs(eigvs[:,:,:,band])

    ###########################
    # interpolate
    tempx = np.linspace(-2.5,2.5,num_intp)
    tempy = tempx.copy()
    Tempx, Tempy = np.meshgrid(tempx,tempy)
    Temp_energy2, Temp_eigvec = energy_hamiltonian(Tempx,Tempy,band=2)
    Energy_itp = interpolate.RectBivariateSpline(tempx,tempy,Temp_energy2)
    phid_itp = interpolate.RectBivariateSpline(tempx,tempy,Temp_eigvec[:,:,0])
    phipx_itp = interpolate.RectBivariateSpline(tempx,tempy,Temp_eigvec[:,:,1])
    phipy_itp = interpolate.RectBivariateSpline(tempx,tempy,Temp_eigvec[:,:,2])
    return Energy_itp, phid_itp, phipx_itp, phipy_itp



#####################################################################
#   Class
#####################################################################
class Energy_band():
    """ The class of energy band. It taks around 30 secs to generate the class. 
    """
    def __init__(self, paras:Paras, num_itp:int):
        # underscored parameters
        self.paras = copy(paras)
        self._mu = self.paras.mu
        self._num_rows = self.paras.num_rows
        self._t1 = self.paras.t1
        self._t2 = self.paras.t2
        self._t3 = self.paras.t3
        self._tz = self.paras.tz
        self._c = self.paras.c 
        self._chi_Delta = self.paras.chi_Delta
        self._chi_Gamma = self.paras.chi_Gamma
        self._kx = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, self.paras.num_rows)
        self._ky = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, self.paras.num_rows)
        self._Kx, self._Ky = np.meshgrid(self._kx, self._ky)
        self._RawEnergyBand = generate_1_band_model(paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self._Kx,self._Ky)
        self._energy_band = self._raw_energy_band - self._mu

        self.num_itp = num_itp
        self.OldEnergyBand_itp, self.phid_itp, self.phipx_itp, self.phipy_itp = generate_3_band_model(self.paras, num_intp = self.num_itp) 
        self.old_energy_band = self.OldEnergyBand_itp(self._ky,self._kx) - self.paras.old_mu
  

    @property
    def energy_band(self): return self._energy_band
    @property
    def mu(self): return self._mu
    @property
    def t1(self): return self._t1
    @property
    def t2(self): return self._t2
    @property
    def t3(self): return self._t3
    @property
    def tz(self): return self._tz
    @property 
    def num_rows(self): return self._num_rows
    @property
    def kx(self): return self._kx
    @property
    def Kx(self): return self._Kx
    @property
    def ky(self): return self._ky
    @property
    def Ky(self): return self._Ky
    @property
    def chi_Delta(self): return self._chi_Delta
    @property
    def chi_Gamma(self): return self._chi_Gamma
    @property
    def c(self): return self._c

    @mu.setter
    def mu(self,value):
        self._mu = value
        self.paras.mu = value
        self._energy_band = self._raw_energy_band - self._mu
    @t1.setter
    def t1(self,value):
        self._t1 = value
        self.paras.t1 = value
        self._RawEnergyBand = generate_1_band_model(self.paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self.Kx,self.Ky)
        self._energy_band = self._raw_energy_band - self._mu
    @t2.setter
    def t2(self,value):
        self._t2 = value
        self.paras.t2 = value
        self._RawEnergyBand = generate_1_band_model(self.paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self.Kx,self.Ky)
        self._energy_band = self._raw_energy_band - self._mu
    @t3.setter
    def t3(self,value):
        self._t3 = value
        self.paras.t3 = value
        self._RawEnergyBand = generate_1_band_model(self.paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self.Kx,self.Ky)
        self._energy_band = self._raw_energy_band - self._mu
    @tz.setter
    def tz(self,value):
        self._tz = value
        self.paras.tz = value
        self._RawEnergyBand = generate_1_band_model(self.paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self.Kx,self.Ky)
        self._energy_band = self._raw_energy_band - self._mu
    @num_rows.setter
    def num_rows(self,value):
        self._num_rows = value
        self.paras.num_rows = value
        self._kx = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, value)
        self._ky = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, value)
        self._Kx, self._Ky = np.meshgrid(self._kx, self._ky)
        self._RawEnergyBand = generate_1_band_model(self.paras) # raw energy band means energy band with mu = 0
        self._raw_energy_band = self._RawEnergyBand(self.Kx,self.Ky)
        self._energy_band = self._raw_energy_band - self._mu
    @chi_Delta.setter
    def chi_Delta(self,value):
        self._chi_Delta = value
        self.paras.chi_Delta = value
    @chi_Gamma.setter
    def chi_Gamma(self,value):
        self._chi_Gamma = value
        self.paras.chi_Gamma = value
    @c.setter
    def c(self,value):
        self._c = value
        self.paras.c = value

class System(Energy_band):
    def __init__(self, paras:Paras, num_itp:int) -> None:
        super().__init__(paras, num_itp)
        
    def energy(self, q:float, isold = False) -> np.array:
        """ return energy band.
        """
        if isold:
            return self.OldEnergyBand_itp(self._ky,self._kx+q) - self.paras.old_mu
        else:
            return self._RawEnergyBand(self._Kx+q,self._Ky) - self.mu

    def phi(self, q:float) -> np.array:
        """ return the phi_{px} orbital character function
        """
        kx = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, self.paras.num_rows)
        ky = kx.copy()
        return self.phipx_itp(ky,kx+q)

    def g(self, q: float, mode = 'ac') -> np.array:
        kx = np.linspace(-self.paras.knorm/2, self.paras.knorm/2, self.paras.num_rows)
        ky = kx.copy()
        k1_mesh, k2_mesh = np.meshgrid(kx,ky)
        a = self.paras.norm
        if mode == 'br':
            p1_mesh = k1_mesh + q
            p2_mesh = k2_mesh
            return np.sqrt(np.sin(q * self.paras.norm/2)**2) * self.paras.EPC[mode]
        elif mode == 'b1g':
            p1_mesh = k1_mesh + q
            p2_mesh = k2_mesh
            return self.paras.EPC[mode] * (np.sin(k1_mesh*a/2)*np.sin(p1_mesh*a/2) - np.sin(k2_mesh*a/2)*np.sin(p2_mesh*a/2)*np.cos(q*a/2))
        elif mode == 'a1g':
            p1_mesh = k1_mesh + q
            p2_mesh = k2_mesh
            return self.paras.EPC[mode] * (np.sin(k1_mesh*a/2)*np.sin(p1_mesh*a/2) + np.sin(k2_mesh*a/2)*np.sin(p2_mesh*a/2)*np.cos(q*a/2))
        elif mode == 'ac':
            return self.paras.EPC[mode] # currently it's a const function
        elif mode == 'ap':
            p1_mesh = k1_mesh + q
            p2_mesh = k2_mesh
            return self.paras.EPC[mode] * (np.cos(k1_mesh*a) - np.cos(k2_mesh*a)) * (np.cos(p1_mesh*a)- np.cos(p2_mesh*a))

    def gamma_e(self, omega:float):
        return self.paras.gamma_e

    def Omega(self, q:float , mode = 'a1g') -> float:  # BIG PROBLEM IN THIS FUNCTION, modify the ac mode dispersion a bit
        """
        dispersion of phonon of a specific mode
        """
        Omega = self.paras.OMEGA
        norm = self.paras.norm

        if mode == 'a1g':
            return Omega['a1g']
        elif mode == 'b1g':
            return  Omega['b1g']
        elif mode == 'br':
            return  Omega['br']*(1-0.2*(np.sin(q*norm/2)**2))
            # return  Omega['br']
        elif mode == 'ac':
            return Omega['ac']
        elif mode == 'ap':
            return Omega['ap']
        else:
            TypeError('the phonon mode should lie in {a1g, b1g, br, ac, ap}')

    def D0(self, q:float, omega:float, mode = 'a1g'):
        """
        BARE phonon propagator of a specific mode
        """
        delta = self.paras.delta
        temp1 = omega - self.Omega(q, mode=mode) + 1j*delta
        temp2 = omega + self.Omega(q, mode=mode) + 1j*delta
        return 1/temp1 - 1/temp2
    
    def D(self, q:float, omega:float, mode = 'a1g'):
        ALPHA = self.paras.ALPHA
        D0 = self.D0(q,omega,mode)
        return D0 / (1 + ALPHA[mode]**2 * D0 * self.chi(q, omega))

    def chi(self, q:float, omega:float):
        q_cdw = self.paras.q_cdw[0]
        Delta = self.paras.chi_Delta
        Gamma = self.paras.chi_Gamma
        c = self.paras.c
        normalization_factor = 1/2 * 1/np.sqrt(Delta**2 + c**2 * (q - q_cdw)**2)
        chi = 1/(Delta**2 + c**2 * (q - q_cdw)**2 - (omega + 1j*(Gamma))**2)
        return chi * normalization_factor
    
    def nf(self, q:float, T = 0) -> np.array:
        """
        Fermi-Dirac distribution function with T in the unit of K.

        return n_f( - energy_(k+q)), which is a meshgrid
        """
        energy = self.energy(q)
        ## CHECK THE DEFINITION AGAIN ##
        kB = 8.6 * 10**(-2)
        if T == 0:
            return np.heaviside(energy,1/2)
        else:
            return 1/(1+np.exp(-energy/(T*kB)))

    def intensity(self, q:float, Omega:float, mode = 'a1g', isfast = False, isold = False):
        """ single point intensity
        """
        paras = self.paras
        gamma_e, Gamma = paras.gamma_e, paras.Gamma
        g = self.g(q, mode = mode)
        nf1 = self.nf(0, T=self.paras.T); nf2 = self.nf(q, T = self.paras.T)
        phi1 = self.phi(0); phi2 = self.phi(q)
        #phi1 = 1; phi2 = 1
        energy1 = self.energy(0, isold=isold); energy2 = self.energy(q, isold=isold)
        matrix1 = ((nf1)/(energy1 + Omega - 1j*Gamma) - (nf2)/(energy2 - 1j*Gamma) ) * g * phi1**2 * phi2**2
        matrix2 = ((nf1)/(energy1 + Omega + 1j*Gamma) - (nf2)/(energy2 + 1j*Gamma)) * g * phi1**2 * phi2**2
        D = self.D(q,Omega, mode = mode)
        temp1 = 1/(energy1 - energy2 + Omega + 1j*self.gamma_e(Omega))
        temp2 = 1/(energy1 - energy2 + Omega - 1j*self.gamma_e(Omega))
        if (not isfast):
            result = (-1/np.pi) * np.sum(np.outer(matrix1,matrix2) * np.imag(np.outer(temp1,temp2) * D)) 
        if isfast:
            ss = np.sum(matrix1 * temp1 * g)
            result = -1/np.pi * ss * ss.conjugate() * np.imag(D) 
        return result
            
    def intensity_mesh(self, q_array:np.array, Omega_array:np.array, mode = 'a1g', full = False, isfast = False, isold = False) -> np.array:
        """ generate the entire spectrum given an array of q and Omega 
        """
        q_mesh, Omega_mesh = np.meshgrid(q_array, Omega_array)
        if not full:
            I = np.ones(np.shape(q_mesh)) * 1.0
            for j in tqdm(range(len(Omega_array))):
                omega = Omega_array[j]
                for i,q in enumerate(q_array):
                    I[j,i] = self.intensity(q, omega, mode = mode, isfast = isfast, isold = isold)
            return q_mesh, Omega_mesh, I

        else:
            I1 = np.ones(np.shape(q_mesh)) * 1.0
            I2 = I1.copy()
            I3 = I1.copy()
            I4 = I1.copy()
            for j in tqdm(range(len(Omega_array))):
                omega = Omega_array[j]
                for i,q in enumerate(q_array):
                    I1[j,i] = self.intensity(q, omega, mode = 'ac', isfast = isfast, isold=isold)
                    I2[j,i] = self.intensity(q, omega, mode = 'b1g', isfast = isfast, isold=isold)
                    I3[j,i] = self.intensity(q, omega, mode = 'a1g', isfast= isfast, isold=isold)
                    I4[j,i] = self.intensity(q, omega, mode = 'br', isfast=isfast, isold=isold) 

            return q_mesh, Omega_mesh, [I1+I2+I3+I4, I1, I2, I3, I4]
        
    def Lindhard(self, q,omega):
        nq = self.nf(q, T = self.paras.T)
        n0 = self.nf(0, T = self.paras.T)
        eq = self.energy(q)
        e0 = self.energy(0)

        term = (nq - n0)/(eq - e0 - omega - 1j*self.gamma_e(omega) )
        return np.sum(term)
    
    def intensity_without_phi(self, q:float, Omega:float, mode = 'a1g', isfast = False, isold = False):
        paras = self.paras
        gamma_e, Gamma = paras.gamma_e, paras.Gamma
        g = self.g(q, mode = mode)
        nf1 = self.nf(0, T=self.paras.T); nf2 = self.nf(q, T = self.paras.T)
        phi1 = 1; phi2 = 1
        energy1 = self.energy(0, isold=isold); energy2 = self.energy(q, isold=isold)
        matrix1 = ((nf1)/(energy1 + Omega - 1j*Gamma) - (nf2)/(energy2 - 1j*Gamma) ) * g * phi1**2 * phi2**2
        matrix2 = ((nf1)/(energy1 + Omega + 1j*Gamma) - (nf2)/(energy2 + 1j*Gamma)) * g * phi1**2 * phi2**2
        D = self.D(q,Omega, mode = mode)
        temp1 = 1/(energy1 - energy2 + Omega + 1j*self.gamma_e(Omega))
        temp2 = 1/(energy1 - energy2 + Omega - 1j*self.gamma_e(Omega))
        if (not isfast):
            result = (-1/np.pi) * np.sum(np.outer(matrix1,matrix2) * np.imag(np.outer(temp1,temp2) * D)) 
        elif isfast:
            ss = np.sum(matrix1 * temp1 * g)
            result = -1/np.pi * ss * ss.conjugate() * np.imag(D) 
        return result

    def Linhard_mesh(self, q_array, omega_array):
        q_mesh, omega_mesh = np.meshgrid(q_array, omega_array)
        I = np.ones(np.shape(q_mesh)) * (1.0+0j)
        for j in tqdm(range(len(omega_array))):
            omega = omega_array[j]
            for i,q in enumerate(q_array):
                I[j,i] = self.Lindhard(q, omega)
        return q_mesh, omega_mesh, I
    
    def intensity_without_phi_mesh(self, q_array, omega_array, mode = 'a1g', isfast = False, isold = False):
        q_mesh, omega_mesh = np.meshgrid(q_array, omega_array)
        I = np.ones(np.shape(q_mesh)) * 1.0
        for j in tqdm(range(len(omega_array))):
            Omega = omega_array[j]
            for i,q in enumerate(q_array):
                I[j,i] = self.intensity_without_phi(q, Omega, mode=mode, isfast = isfast, isold=isold)
        return q_mesh, omega_mesh, I
    
    def occupation_ratio(self, isold = False):
        E = self.energy(0,isold)
        total_num = np.sum(E/E)*1.0
        occupation = np.sum(E<=0)*1.0
        return occupation/total_num
    
    def S(self, q_array:np.array, Omega_array:np.array, T = 0):
        q_mesh, Omega_mesh = np.meshgrid(q_array, Omega_array)
        chi_mesh = self.chi(q_mesh, Omega_mesh)
        if T != 0:
            S_mesh = 2/(1-np.exp(-Omega_mesh/T))*np.imag(chi_mesh)
        else:
            S_mesh = 2*np.heaviside(Omega_mesh,0.5)*np.imag(chi_mesh) 
        return S_mesh
        



####################################################################
# visualization
def image(ax,X,Y,Z):
    xr = [np.min(X[0,:]) * norm / (2*np.pi), np.max(X[0,:])* norm / (2*np.pi)]
    yr = [np.min(Y[:,0]), np.max(Y[:,0])]
    ax.imshow(Z,origin='lower', aspect=(xr[1]-xr[0])/(yr[1]-yr[0]),extent = (*xr, *yr))
    AXIS = dict(xlabel = '$q$ ($2\pi/a$)',
                ylabel = '$\Omega (meV)$')
    ax.set(**AXIS)
    return ax
    


#####################################
# create system
def create_system():
    return System(paras, 1000)


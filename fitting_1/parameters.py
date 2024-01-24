import numpy as np
from types import SimpleNamespace
from copy import copy
# assign parameters
a = np.array([3.79,0]); b = np.array([0,3.79]) # lattice in real sapce, in the unit of A
c_latt = 13.2
norm = a[0]
temp = np.array([a,b]); Brill = 2*np.pi*np.linalg.inv(temp)
ka,kb = [Brill[:,0], Brill[:,1]] # lattice in k-space; 
knorm = ka[0]
mu = -134.3; # chemical potential
old_mu = 2114 # standard chemical potential
num_rows = 100; num_grids = num_rows**2; # simulation grids
g_ac, g_b1g, g_a1g, g_br, g_ap = [1000,1000,1000,1000,1000] # EPC coefficients
EPC_coefficient = {'a1g': g_a1g, 'b1g': g_b1g, 'br': g_br, 'ac': g_ac, 'ap': g_ap}
Omega_br, Omega_a1g, Omega_b1g, Omega_ac, Omega_ap= [80.0, 52.0, 42 ,20.0, 52] # parameters in dispersion relation of phonon
OMEGA = {'a1g': Omega_a1g, 'b1g': Omega_b1g, 'ac': Omega_ac, 'br': Omega_br, 'ap': Omega_ap}
alpha_ac, alpha_b1g, alpha_a1g, alpha_br, alpha_ap = np.array([1.2,1.45,1.,1.,1.])*0 # CDW-phonon interaction coefficient
ALPHA = {'a1g': alpha_a1g, 'b1g': alpha_b1g, 'br': alpha_br, 'ac': alpha_ac, 'ap': alpha_ap}
Delta, c, delta = [20., 300*norm, 22.] # parameters in CDW susceptibility, question in definition of c
q_cdw = np.array([0.235*2*np.pi/norm,0]) # wave_vector of CDW
sigma_q, sigma_Omega = [0.023*2*np.pi/norm, 5] # resolution funtion
Gamma_ch = 750 # core-hole lifetime parameter, meV
gamma_e = 100  # electron lifetime in meV
chi_Delta = 20 # 
chi_Gamma = 20 #
mb = 0
t1, t2, t3, tz = [190, -28.5, 14.25, 5.7]

class SettableNamespace(SimpleNamespace):
    """Contains a collection of parameters.

    Used to make a System object.

    Takes keyword arguments and stores them as attributes.
    """
    def __init__(self, namespace=None, **kwargs):
        super().__init__()
        if namespace:
            self.__dict__.update(namespace.__dict__)
        self.__dict__.update(kwargs)

    def get(self, name, default=None):
        """Look up a variable.

        name: string varname
        default: value returned if `name` is not present
        """
        try:
            return self.__getattribute__(name, default)
        except AttributeError:
            return default

    def set(self, **variables):
        """Make a copy and update the given variables.

        returns: Params
        """
        new = copy(self)
        new.__dict__.update(variables)
        return new

class Paras(SettableNamespace):
    pass


##############################################
# create parameters class
##############################################
paras_init = Paras(
    a = a, b = b,  # lattice constants
    ka = ka, kb = kb,  # reciprocal lattices
    norm = norm, knorm = knorm,
    EPC = EPC_coefficient, ALPHA = ALPHA, OMEGA = OMEGA,
    delta = delta, Gamma = Gamma_ch, gamma_e = gamma_e,
    q_cdw = q_cdw,
    sigma_q = sigma_q, sigma_Omega = sigma_Omega,
    mu = mu,
    num_rows = num_rows,
    old_mu = old_mu,
    T=0,
    c = c,
    c_latt = c_latt,
    chi_Gamma= chi_Gamma, chi_Delta = chi_Delta,
    mb=mb,
    t1 = t1, t2 = t2, t3 = t3, tz = tz
)
paras = copy(paras_init)

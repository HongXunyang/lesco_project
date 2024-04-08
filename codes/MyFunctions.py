import scipy as sp
import numpy as np
from scipy.signal import convolve
from rixs_experiment import RIXS_EXP
_temp_exp = RIXS_EXP()

def fl(x,x0,A,k):
    # lorentzian curve
    temp = A*k**2/((x-x0)**2 + k**2)
    return temp
def fq(x,a,b):
    # quadratic curve
    y=x
    temp = np.where(y<=0, 0 , -(a/b**2)*(x-b)**2 + a)
    return temp

def fg(x,x0,A,k):
    # gaussian curve
    sigma = k/(np.sqrt(2)*np.sqrt(np.log(2)))
    y = A*np.exp(-(x-x0)**2/(2*sigma**2))
    return y

def fv(x,x0,A,fG,fL):
    """
    the function of Vogit profile with peak value A, centered at x0
    
    key arguements:
    -----------------------------------------
    - (array) | x: the input x 
    - (float) | x0: the center of the profile 
    - (float) | A: the peak value of the profile 
    - (float) | fG: FWHM of the gaussian function 
    - (float) | fL: FWHM of the Lorentzian function
    
    return:
    - (array) | y: return Vogit(x) 
    """
    y = sp.special.voigt_profile((x-x0),fG/(2*np.sqrt(2*np.log(2))),fL/2)
    y = y/max(y)*A
    return y




def fun_voigt(x,x0,A,fG,fL):
    """
    the function of Vogit profile with peak value A, centered at x0
    
    key arguements:
    -----------------------------------------
    - (array) | x: the input x 
    - (float) | x0: the center of the profile 
    - (float) | A: the peak value of the profile 
    - (float) | fG: FWHM of the gaussian function 
    - (float) | fL: FWHM of the Lorentzian function
    
    return:
    - (array) | y: return Vogit(x) 
    """
    y = sp.special.voigt_profile((x-x0),fG/(2*np.sqrt(2*np.log(2))),fL/2)
    y = y/max(y)*A
    return y

def fun_pvoigt(x, x0, A, res, mu):
    """pseudo-Voigt profile function
    
    Keyword arguments:
    ----------------------------------
    - (array) | x: the input variable
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile
    - (float) | res: the FWHM of the profile
    - (float) | mu: the mixing parameter of the Lorentzian and Gaussian profile. mu = 1 -> Lorentzian, mu = 0 -> Gaussian.
    
    return:
    - (array) | y: return the pseudo-Voigt profile
    """
    
    return mu * fun_lorentzian(x, x0, A, res) + (1 - mu) * fun_gaussian(x, x0, A, res)

def fun_lorentzian(x,x0,A,res):
    """
    the function of Lorentzian profile with peak value A, centered at x0

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile
    - (float) | res: the FWHM of the profile

    return:
    - (array) | y: return Lorentzian(x)
    """
    gamma = res/2
    temp = A*gamma**2/((x-x0)**2 + gamma**2)
    return temp

def fun_quadratic(x,x0,A):
    """
    the function of quadratic profile with peak value A, centered at x0

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile

    return:
    - (array) | y: return quadratic(x)
    """
    y=x
    temp = np.where(y<=0, 0 , -(A/x0**2)*(x-x0)**2 + A)
    return temp

def fun_gaussian(x,x0,A,res):
    """
    the function of Gaussian profile with peak value A, centered at x0

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile
    - (float) | res: the FWHM of the profile

    return:
    - (array) | y: return Gaussian(x)
    """
    gamma = res/2
    sigma = gamma/(np.sqrt(2*np.log(2)))
    y = A*np.exp(-(x-x0)**2/(2*sigma**2))
    return y


def fun_DHO(x, x0, A, res):
    """
    the function of DHO profile with peak value A, centered at x0

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile
    - (float) | res: the FWHM of the profile

    return:
    - (array) | y: return DHO(x)
    """
    gamma = res/2
    nominator = A*(gamma*x)**2
    denominator = (x**2 - x0**2)**2 + (x*gamma)**2
    y = nominator/denominator
    return np.where(x<0,0,y)

def fun_DHO_original(x, x0, A, res):
    """
    the function of DHO profile with peak value A, centered at x0

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | x0: the center of the profile
    - (float) | A: the peak value of the profile
    - (float) | res: the FWHM of the profile

    return:
    - (array) | y: return DHO(x)
    """
    gamma = res/2
    nominator = A*(gamma*x)**2
    denominator = (x**2 - x0**2)**2 + (x*gamma)**2
    y = nominator/denominator
    return y

# convolve functions with a given gaussian function with given FWHM
def convoluted_fun(x, res, fun, *parameters):
    """
    convolute the given function with a Gaussian function

    Keyword arguements:
    -----------------------------------------
    - (array) | x: the input x
    - (float) | res: the FWHM of the Gaussian function
    - (function) | fun: the function to be convoluted
    - (float) | parameters: the parameters of the function

    return:
    - (array) | y: return the convoluted function
    """
    # Generate a range for the Gaussian, centered at 0 with some width
    x_left = min(x)
    x_right = max(x)
    x_left_shifted = x_left - (x_right + x_left)/2
    x_right_shifted = x_right - (x_right + x_left)/2
    x_gaussian = np.linspace(x_left_shifted, x_right_shifted, len(x))
    gaussian_curve = fun_gaussian(x_gaussian, 0, 1, res)
    
    # Calculate fun_s for a range of x values
    fun_values = fun(x, *parameters)
    max_value = max(fun_values)
    # Convolute fun_s with the Gaussian
    # Assume x is a numpy array for simplicity; adjust as necessary for your case
    convoluted_values = convolve(fun_values, gaussian_curve, mode='same')
    max_convoluted_value = max(convoluted_values)
    return convoluted_values/max_convoluted_value*max_value



def fun_chi(y,q,c,A,Gamma,Delta):   
        q_cdf = 0.24
        #normalization_factor = 1/2 * 1/np.sqrt(c**2*(x-q_cdf)**2+Delta**2) # do I need this?????????????
        normalization_factor = 1
        chi = 1/(Delta**2 + c**2*(q-q_cdf)**2- (y + 1j*(Gamma))**2)
        output = np.imag(chi * normalization_factor)
        return output/np.max(output)*A



def fun_S_old(y,q,c,A,Gamma,Delta,T):
    out = 2/(1-np.exp(-y/T))*fun_chi(y,q,c,A,Gamma,Delta)
    return out/max(out) * A

def fun_S_convoluted_old(y, q, c, A, Gamma, Delta, T):
    """
    Calculate the convoluted function of S(x, y) with a Gaussian resolution function.

    Parameters:
    x (float): Momentum (r.l.u.)
    y (array-like): Energy loss (meV)
    c (float): Slope of the dispersion
    A (float): Amplitude
    Gamma (float): Damping (meV)
    Delta (float): Energy gap (meV)
    T (float): Temperature in meV

    Returns:
    array-like: The convoluted function of S(x, y) with a Gaussian resolution function, normalized to the maximum value.

    """
    y_left = min(y)
    y_right = max(y)
    y_left_shifted = y_left - (y_right + y_left) / 2
    y_right_shifted = y_right - (y_right + y_left) / 2
    y_gaussian = np.linspace(y_left_shifted, y_right_shifted, len(y))
    res = _temp_exp._resolution
    z = fun_S_old(y,q,c, A, Gamma, Delta, T)
    z_conv = convolve(z, fun_gaussian(y_gaussian, 0, 1, res), mode='same')
    return z_conv / max(z_conv) * A


def fun_S_new(y,q,c,A,Gamma,Delta,T):
    temp = fun_chi(y,q,c,A,Gamma,Delta)*(n(y,T)+1)
    # if temp[y<0] is not empty
    if temp[y<0].size !=0:
        temp[y<0] = np.exp(y[y<0]/T)*fun_chi(-y[y<0],q,c,A,Gamma,Delta)
    return temp/max(temp) * A

def n(y,T):
    return 1/(np.exp(y/T)-1)

def fun_S_convoluted_new(y, q, c, A, Gamma, Delta, T):
    """
    Calculate the convoluted function of S(x, y) with a Gaussian resolution function.

    Parameters:
    x (float): Momentum (r.l.u.)
    y (array-like): Energy loss (meV)
    c (float): Slope of the dispersion
    A (float): Amplitude
    Gamma (float): Damping (meV)
    Delta (float): Energy gap (meV)
    T (float): Temperature in meV

    Returns:
    array-like: The convoluted function of S(x, y) with a Gaussian resolution function, normalized to the maximum value.

    """
    y_left = np.min(y)
    y_right = np.max(y)
    y_left_shifted = y_left - (y_right + y_left) / 2
    y_right_shifted = y_right - (y_right + y_left) / 2
    y_gaussian = np.linspace(y_left_shifted, y_right_shifted, len(y))
    res = _temp_exp._resolution
    z = fun_S_new(y, q,c, A, Gamma, Delta, T)
    z_conv = convolve(z, fun_gaussian(y_gaussian, 0, 1, res), mode='same')
    return z_conv / max(z_conv) * A

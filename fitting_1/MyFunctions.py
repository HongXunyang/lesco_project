import scipy as sp
import numpy as np

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
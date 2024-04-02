
import numpy as np
from scipy.optimize import curve_fit



class RIXS_EXP():
    def __init__(self):
        self.__name__ = 'rixs_experiment'
        self.__material__ = "LSCO-Eu at p=0.125"
        self.__edge__ = "Oxygen K-edge"
        self.__incident_energy__ = 529
        self._temperature_list_K = [21,62,104,155]
        self._temperature_list_meV = np.array([21,62,104,155])*0.0862  
        self._temperature_list_string = ['T21K','T62K','T104K','T155K']  
        self._q_list = np.array([0.05, 0.1, 0.15, 0.18, 0.2, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.287])
        self._q_size = len(self._q_list)
        self._temperature_size = len(self._temperature_list_K)
        self._resolution = 22.22
        self._fG_preset = 16.875
        self._fL_preset = 9
        self._Raw_data = dict(  T21K = np.genfromtxt(('./simulation_data/'+'21K'+'.csv'), delimiter=','),
                                T62K = np.genfromtxt(('./simulation_data/'+'62K'+'.csv'), delimiter=','),
                                T104K = np.genfromtxt(('./simulation_data/'+'104K'+'.csv'), delimiter=','),
                                T155K = np.genfromtxt(('./simulation_data/'+'155K'+'.csv'), delimiter=','),\
                                    info="")
        self._Raw_data['T21K'] = self._Raw_data['T21K'][0:1195:,:]
        self._Raw_data['T62K'] = self._Raw_data['T62K'][0:1195:,:]
        self._Raw_data['T104K'] = self._Raw_data['T104K'][0:1195:,:]
        self._Raw_data['T155K'] = self._Raw_data['T155K'][0:1195:,:]
        self._interesting_range_index = [246, 307]
        self._interesting_range = [-80, 150]
        self._Interesting_raw_data = dict(\
            T21K = self._Raw_data['T21K'][self._interesting_range_index[0]:self._interesting_range_index[1],:], 
            T62K = self._Raw_data['T62K'][self._interesting_range_index[0]:self._interesting_range_index[1],:], 
            T104K = self._Raw_data['T104K'][self._interesting_range_index[0]:self._interesting_range_index[1],:], 
            T155K = self._Raw_data['T155K'][self._interesting_range_index[0]:self._interesting_range_index[1],:],\
                 info="")
        
        self.Aligned_interesting_data = dict(T21K = np.copy((self._Interesting_raw_data['T21K'])),\
                                                T62K = np.copy((self._Interesting_raw_data['T62K'])),\
                                                T104K = np.copy((self._Interesting_raw_data['T104K'])),\
                                                T155K = np.copy((self._Interesting_raw_data['T155K'])),\
                                                info="")

                                            
        self.Fit_results = dict(fitting_info=[], T21K = dict(optimized_parameters=[[] for i in range(self._q_size)],errors=[[] for i in range(self._q_size)]),\
                                T62K = dict(optimized_parameters=[[] for i in range(self._q_size)],errors=[[] for i in range(self._q_size)]),\
                                T104K = dict(optimized_parameters=[[] for i in range(self._q_size)],errors=[[] for i in range(self._q_size)]),\
                                T155K = dict(optimized_parameters=[[] for i in range(self._q_size)],errors=[[] for i in range(self._q_size)]),\
                                    info="")
        self.Subtracted_realigned_interesting_data = dict(T21K = np.zeros(np.shape(self._Interesting_raw_data['T21K'])),\
                                                T62K = np.zeros(np.shape(self._Interesting_raw_data['T62K'])),\
                                                T104K = np.zeros(np.shape(self._Interesting_raw_data['T104K'])),\
                                                T155K = np.zeros(np.shape(self._Interesting_raw_data['T155K'])),\
                                                    info="")
        self.Interpolated_subtracted_realigned_interesting_data = dict(T21K = np.zeros(np.shape(self._Interesting_raw_data['T21K'])),\
                                                T62K = np.zeros(np.shape(self._Interesting_raw_data['T62K'])),\
                                                T104K = np.zeros(np.shape(self._Interesting_raw_data['T104K'])),\
                                                T155K = np.zeros(np.shape(self._Interesting_raw_data['T155K'])),\
                                                    info="")
        self.Pure_CDF_data = dict(T21K = None, T62K = None, T104K = None, T155K = None, info="Subtracted spectrum showing only CDF data. all are mesh. Z_subtracted are spectrum subtracted by elastic, a1g,b1g,breathing phonons. Z_extra_subtracted are additionally subtracted by the first column (lowest q) of the spectrum (presumingly acoustic phonon)")
        self.Phonon_simulation = dict(T21K = None, T62K = None, T104K = None, T155K = None, info="Simulation data of . All are mesh. Phonons are convoluted by a gaussian (FWHM = resolution)")

    @property
    def Raw_data(self):
        return self._Raw_data
    @property
    def q_list(self):
        return self._q_list
    @property
    def temperature_list_K(self):
        return self._temperature_list_K
    @property
    def temperature_list_eV(self):
        return self._temperature_list_eV
    @property
    def temperature_list_string(self):
        return self._temperature_list_string
    @property
    def Interesting_raw_data(self):
        return self._Interesting_raw_data

    def to_dict(self):
        """
        store all entries to a dictionary
        """
        return 0


    def info(self):
        print("\n call `q_list` for the list of q in the experiment \n \
               call `temperature_list_K` for the list of temperature in Kelvin \n \
               call `temperature_list_eV` for the list of temperature in eV \n \
               call `Raw_data` for the raw data of the experiment. \n")
    
    def align_data(self, method = None):
        if method == None:
            T_list_chars = self._temperature_list_string
            for T_char in T_list_chars:
                data_T = self.Interesting_raw_data[T_char]
                for i in range(self._q_size):
                    energy = data_T[:,2*i]
                    intensity = data_T[:,2*i+1]
                    max_energy = energy[np.argmax(intensity)]
                    aligned_energy = energy - max_energy
                    self.Aligned_interesting_data[T_char][:,2*i] = aligned_energy
        elif method == 'gaussian':
            T_list_chars = self._temperature_list_string
            for T_char in T_list_chars:
                data_T = self.Interesting_raw_data[T_char]
                for i in range(self._q_size):
                    energy = data_T[:,2*i]
                    intensity = data_T[:,2*i+1]
                    resolution = self._resolution
                    # fit the data with a gaussian
                    fit_energy = energy[(energy>-resolution/2) & (energy<resolution/2)]
                    fit_intensity = intensity[(energy>-resolution/2) & (energy<resolution/2)]
                    popt,_ = curve_fit(mf.fun_gaussian, fit_energy, fit_intensity, p0=[0, 5, resolution])
                    max_energy = popt[0]
                    aligned_energy = energy - max_energy
                    self.Aligned_interesting_data[T_char][:,2*i] = aligned_energy


    def plot_data(self, T_char, q_index, ax = None):
        """
        Plot the raw data of the experiment
        
        keyword arguments:
        -----------------------------
        - (char) | T_char : 'T21K', 'T62K', 'T104K', 'T155K'
        - (int) | q_index : from 0 to 12
        - (ax) | ax : matplotlib axis to plot the data
        """
        data_T = self.Aligned_interesting_data[T_char]
        energy = data_T[:,2*q_index]
        intensity = data_T[:,2*q_index+1]
        if ax == None:
            plt.plot(energy, intensity, '.', label = T_char+' q='+str(self.q_list[q_index]), markersize = 7)
            ax.legend()
            ax.legend(loc='upper right')
        else:
            ax.plot(energy, intensity, '.', label = T_char+' q='+str(self.q_list[q_index]), markersize = 7)
            ax.legend()
            ax.legend(loc='upper right')
 

    def basic_fit(self, temperature_char, q_index, fun, parameters_guess, parameters_bound, sigma = None):
        """
        Fit the data with a given function
        
        keyword arguments:
        -----------------------------
        - (char) | temperature_char : 'T21K', 'T62K', 'T104K', 'T155K'
        - (int) | q_index : from 0 to 12
        - (function) | fun : function to fit the data
        - (list) | parameters_guess : guess of the parameters of the function
        - (list) | parameters_bound : bounds of the parameters of the function
        """
        data_T = self.Aligned_interesting_data[temperature_char]
        energy = data_T[:,2*q_index]
        intensity = data_T[:,2*q_index+1]
        popt, pcov = curve_fit(fun, energy, intensity, p0=parameters_guess, bounds=parameters_bound, sigma=sigma)
        return popt, np.sqrt(np.diag(pcov))
    
    def plot_fit_conponent(self, fun_conponent, parameters, fill = False, ax = None, label = None, color = 'lightgrey'):
        """
        Plot the fit of the data with a given function
        
        keyword arguments:
        -----------------------------
        - (function) | fun_conponent : function to fit the data
        - (array) | parameters : parameters of the function
        """
        energy = np.linspace(self._interesting_range[0], self._interesting_range[1], 300)
        intensity = fun_conponent(energy, *parameters)
        if ax == None:
            if not fill:
                plt.plot(energy, intensity, label = label, color = color)
            else:
                plt.fill_between(energy, intensity, alpha=0.8, label = label, color = color, edgecolor = 'none')
        else:
            if not fill:
                ax.plot(energy, intensity, label = label, color = color)
            else:
                ax.fill_between(energy, intensity, alpha=0.8, label = label, color = color, edgecolor = 'none')


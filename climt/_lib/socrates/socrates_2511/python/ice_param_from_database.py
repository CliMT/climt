# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Derive parametrizations for ice optical properties from database.

Derive ice optical property parametrizations from a database of
ice optical properties as a function of wavelength and ice water 
content (e.g. code/data/cloud/scatter_ice_ens013). This was used to 
develop the type 10 parametrization, with a rescaling of 0.64039
'''

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from os.path import expandvars

class IceData:
    '''
    Class to hold ice optical property data.

    :param int ndata: number of data points
    :param array wavelength: 1D array of wavelengths
    :param array qi: 1D array of ice mass mixing ratio
    :param array kext: 1D array of mass extinction coefficients 
    :param array ksca: 1D array of mass scattering coefficients
    :param array w0: 1D array of single scattering albedos (equals kscat/kext)
    :param array g: 1D array of asymmetry, the first moment of the phase function.
    '''
    def __init__(self, filename):
        with open(filename, "r") as f:
            data = f.read().strip().split()
            data = np.array([float(i) for i in data[data.index('*BEGIN_DATA')+1:data.index('*END')]])
            self.ndata = len(data)/6
            data = data.reshape(self.ndata, 6)
            self.wavelength = data[:, 0]
            self.qi = data[:,1]
            self.kext = data[:,2]
            self.ksca = data[:,3]
            self.w0 = data[:,4]
            self.g = data[:,5]

class SpecData:
    '''
    Class to hold band limits read in from spectral file.

    :param int n_band: Number of spectral bands
    :param array lower_bound: 1D array of lower bound limits of bands
    :param array upper_bound: 1D array of upper bound limits of bands
    '''
    def __init__(self, filename):
        with open(filename, "r") as f:
            data = f.read().strip().replace(",", "").split()
            if data[0] == '&R2SWSP':
                n_band_ind = data.index('N_BAND')
                self.n_band = int(data[n_band_ind+2])
                wave_length_short_ind = data.index('WAVE_LENGTH_SHORT')
                self.lower_bound =  [float(i) for i in (data[wave_length_short_ind + 2 : wave_length_short_ind + 2 + self.n_band])]
                self.upper_bound = [float(i) for i in (data[wave_length_short_ind + 4 + self.n_band : wave_length_short_ind + 4 + 2 * self.n_band])]
            elif data[0] == '*BLOCK:':
                n_band_ind = data.index('bands')
                self.n_band = int(data[n_band_ind+2])
                wave_length_short_ind = data.index('Lower')
                wavelength_data = np.array([float(i) for i in data[wave_length_short_ind + 4 : wave_length_short_ind + 4 + 3 * self.n_band]]).reshape(self.n_band,3)
                self.lower_bound = wavelength_data[:,1]
                self.upper_bound = wavelength_data[:,2]
            else:
                raise NameError('Unable to read spectral file')

class Param_from_file:
    '''
    Class to hold parameters for ice optical properties read in from spectral file.

    :param array parameters: array of parameters for ice optical properties parametrization, size depends on parametrization
    '''
    def __init__(self, filename, SpecData, index_param):
        with open(filename, "r") as f:
            data = f.read().strip().replace(",", "").split()
            param_ind1 = data.index('ICE_PARAMETER_LIST(11'+index_param+')')
            param_ind2 = data.index('ICE_PARAMETER_LIST(12'+index_param+')')
            n_parameters = param_ind2 - param_ind1 - 2
            self.parameters = np.zeros(SpecData.n_band * n_parameters).reshape(SpecData.n_band, n_parameters)
            for i in range(SpecData.n_band):
                param_ind1 = data.index('ICE_PARAMETER_LIST(1'+str(i+1)+index_param+')')
                self.parameters[i, :] = [float(j) for j in (data[param_ind1 + 2 : param_ind1 + 2 + n_parameters])]

class filename:
    '''
    Class to hold filename for spectral files and single scattering database.

    :param str sw_spec: SW spectral file name
    :param str lw_spec: LW spectral file name
    :param str ss_ice: Ice single scattering database file name
    '''
    def __init__(self):
        self.sw_spec = get_filename('SW spectral')
        self.lw_spec = get_filename('LW spectral')
        self.ss_ice = get_filename('single scattering')


class parameters:
    '''
    Class to hold the derived parameters in each spectral band.

    :param array kext: n_band x 3 array of parameters for mass extinction
    :param array ksca: n_band x 3 array of parameters for mass scattering
    :param array g: n_band x 3 array of parameters for asymmetry
    '''
    def __init__(self, n_band):
        self.kext = np.zeros(n_band * 3).reshape(n_band, 3)
        self.ksca = np.zeros(n_band * 3).reshape(n_band, 3)
        self.g = np.zeros(n_band * 3).reshape(n_band, 3)


class spectrum_data:
    '''
    Class to hold solar spectrum
    '''
    def __init__(self, filename):
        with open(filename, "r") as f:
            data = f.read().strip().split()
            data = np.array([float(i) for i in data[data.index('*BEGIN_DATA')+1:data.index('*END')]])
            self.ndata = len(data)/2
            data = data.reshape(self.ndata, 2)
            self.wavelength = data[:, 0]
            self.irradiance = data[:, 1]
        


def get_filename(filetype):
    '''
    Request and read file name from user.
    '''
    while True:
        try:
            filename = expandvars(input("Enter file containing " + filetype + " data\n"))
            print("You entered '" + filename + "'")
            f = open(filename, "r")
            break
        except IOError:
            print("Unable to open '" + filename + "'. Please try again")
        else:
            f.close()
    return filename   


def identify_region(spectrum, database):
    '''
    Find the indeces of the data with wavelengths in the given region
    '''
    ind1 = database.wavelength < np.max(spectrum.upper_bound)
    ind2 = database.wavelength >= np.min(spectrum.lower_bound)
    region = ind1*ind2
    return region

def planck_weights(wvlen, wvlen_lower, wvlen_upper, temperature):
    '''
    Calculate weights for each spectral point from Planck function
    '''
    h = 6.626068963 * 10**(-34)
    c = 2.99792548 * 10**(8) 
    k = 1.3806504 * 10**(-23)

    n_planck = 100000
    n_wvlen = len(set(wvlen))
    n_dim2 = len(wvlen) / n_wvlen
    wvlen_unique = wvlen[0:n_wvlen]
    Planckian = []
    for i in range(n_wvlen):
        if i == 0:
            wv1 = wvlen_lower
            wv2 = wvlen_unique[i] + (wvlen_unique[i+1] - wvlen_unique[i])/2.0
        elif i == n_wvlen-1:
            wv1 = wv2
            wv2 = wvlen_upper
        else:
            wv1 = wv2
            wv2 = wvlen_unique[i] + (wvlen_unique[i+1] - wvlen_unique[i])/2.0
        delta_wv = (wv2 - wv1) / n_planck
        wv = wv1 + np.arange(n_planck)*delta_wv
        B = (2*h*(c**2)/(wv**5)) / (np.exp((h*c)/(wv*k*temperature)) - 1)
        B = np.mean(B)
        Planckian.append(B * (wv2-wv1))
    Planck_sum = np.sum(Planckian)
    weights = [i/Planck_sum for i in Planckian]
    weights = weights * n_dim2
    return weights


def spectral_weights(wvlen, wvlen_lower, wvlen_upper):
    '''
    Calculate weights for each spectral point using data from an input file
    '''
    global spectrum_file
    try:
        spectrum_file
    except NameError:
        spectrum_file = get_filename('spectrum')
    spectrum = spectrum_data(spectrum_file)
    n_wvlen = len(set(wvlen))
    n_dim2 = len(wvlen) / n_wvlen
    wvlen_unique = wvlen[0:n_wvlen]
    spectral_weight = []
    for i in range(n_wvlen):
        if i == 0:
            wv1 = wvlen_lower
            wv2 = wvlen_unique[i] + (wvlen_unique[i+1] - wvlen_unique[i])/2.0
        elif i == n_wvlen-1:
            wv1 = wv2
            wv2 = wvlen_upper
        else:
            wv1 = wv2
            wv2 = wvlen_unique[i] + (wvlen_unique[i+1] - wvlen_unique[i])/2.0
        ind1 = spectrum.wavelength < wv2
        ind2 = spectrum.wavelength >= wv1
        region = ind1*ind2
        if sum(region > 0):
            spectral_weight.append(np.mean(spectrum.irradiance[region])* (wv2-wv1))
        else:
             spectral_weight.append(np.interp(wvlen_unique[i], spectrum.wavelength, spectrum.irradiance))
    spectral_sum = np.sum(spectral_weight)
    weights = [i/spectral_sum for i in spectral_weight]
    weights = weights * n_dim2
    return weights


def linear_qi(variable, A, B):
    '''
    For variable x and constant A, calculate Ax
    '''
    return A*(variable)

def power_qi(variable, A, B):
    '''
    For variable x and constants A,B, calculate A(x^B)
    '''
    return A*(variable**B)

def constant(variable, A, B):
    '''
    For variable x and constant A, return A(x^0)
    '''
    return A*(variable**0)

def fit_opt_prop(wavelength, qi, opt_prop, spectral_data, opt_prop_name, function_name, plot=False, weighting=None):
    '''
    Fit given ice optical property.

    Fit ice optical property as a function of ice water content.
    Includes option to plot the fits.
    '''
    comb_coeffs = np.zeros(2*spectral_data.n_band).reshape(spectral_data.n_band, 2)
    comb_fit = np.zeros(len(opt_prop))
    for i in range(spectral_data.n_band):
        band_ind = list(np.where((np.array(wavelength) < spectral_data.upper_bound[i]) * (np.array(wavelength) >= spectral_data.lower_bound[i]))[0])
        if (weighting == 'Spectrum'):
            weights = spectral_weights(wavelength[band_ind], spectral_data.lower_bound[i], spectral_data.upper_bound[i])
        elif (weighting == 'Planck'):
            weights = planck_weights(wavelength[band_ind], spectral_data.lower_bound[i], spectral_data.upper_bound[i], 255.0)
        else:
            weights = None
#        comb_coeffs[i,:], comb_cov = curve_fit(eval(function_name), qi[band_ind], opt_prop[band_ind])
        comb_coeffs[i,:], comb_cov = curve_fit(eval(function_name), qi[band_ind], opt_prop[band_ind], sigma=weights)
        comb_fit[band_ind] = eval(function_name)(qi[band_ind], comb_coeffs[i, 0], comb_coeffs[i, 1])
    if plot == True:
        plot_opt_prop(wavelength, qi, opt_prop, opt_prop_name, comb_fit)
    return comb_coeffs


def plot_opt_prop(wavelength, qi, opt_prop, opt_prop_name, comb_fit):
    '''
    Plot fits to ice optical property data
    '''
    fig = plt.figure()
    fig.subplots_adjust(hspace=.4)
    fig.subplots_adjust(wspace=.4)
    ax = fig.add_subplot(221)
    unique_wavelength = np.array(sorted(list(set(list(wavelength)))))
    n_unique_wavelength = len(unique_wavelength)
    unique_qi = np.array(sorted(list(set(list(qi))), reverse = True))
    n_unique_qi = len(unique_qi)
    opt_prop2d = opt_prop.reshape(n_unique_qi, n_unique_wavelength)
    plt.pcolor((unique_wavelength), unique_qi, opt_prop2d, vmax=abs(opt_prop2d).max())
    ax.set_xscale('log')
    plt.colorbar()
    plt.title(opt_prop_name)
    plt.xlabel('wavelength (m)')
    plt.yticks(np.arange(5)*0.001)
    plt.ylabel('ice water content (kg/kg)')

    ax = fig.add_subplot(222)
    plt.scatter(wavelength, opt_prop, color='blue', marker = 'o')
    plt.scatter(wavelength, comb_fit, color='red', marker = '+')
    ax.set_xscale('log')
    plt.xlabel('wavelength (m)')
    plt.ylabel(opt_prop_name)
    minw = 10**np.floor(np.log10(np.min(wavelength)))
    maxw = 10**np.ceil(np.log10(np.max(wavelength)))
    plt.xlim([minw, maxw])
    ax.locator_params(axis='y', nbins=5) 

    ax = fig.add_subplot(223)
    minqi = 0.0
    maxqi = 0.004
    plt.scatter(qi, opt_prop, color = 'blue', marker = 'o')
    plt.scatter(qi, comb_fit, color = 'red', marker = '+')
    plt.xlim([minqi, maxqi])
    plt.xlabel('ice water content (kg/kg)')
    plt.xticks(np.arange(5)*0.001)
    plt.ylabel(opt_prop_name)
    ax.locator_params(axis='y', nbins=5) 

    ax = fig.add_subplot(224)
    plt.plot(sorted(opt_prop), sorted(opt_prop), color='blue')
    plt.scatter(opt_prop, comb_fit, color='red')
    plt.xlabel('Truth')
    plt.ylabel('Parametrization')
    ax.locator_params(axis='both', nbins=5) 
    plt.savefig('ice_'+opt_prop_name[0:2]+'_'+opt_prop_name[3:]+'.png')
    print('Figure saved to "ice_'+opt_prop_name[0:2]+'_'+opt_prop_name[3:]+'.png"')
    plt.show()


def write_opt_prop_param(output_file, Dataset, SW_specfile, LW_specfile, SW_spectrum, LW_spectrum, SW_params, LW_params):
    '''
    Write derived parametrizations to file in reader friendly format.
    '''
    file = open(output_file, 'w')
    file.write('Dataset='+Dataset+'\n')
    file.write('SW spectral file='+SW_specfile+'\n')
    file.write('LW spectral file='+LW_specfile+'\n')
    file.write('\n')
    write_spectral_region(file, 'SW', SW_spectrum, SW_params)
    file.write('\n')
    write_spectral_region(file, 'LW', LW_spectrum, LW_params)
    file.close()

def write_opt_prop_prep_spec(output_file, spectrum, params):
    '''
    Write derived parametrizations to file in format required by prep_spec.
    '''
    file = open(output_file, 'w')
    file.write('*FILE_TYPE =    24' + '\n\n' + '                    Parameters fitting single scattering' + '\n')
    file.write('\n' + 'Index of parametrization scheme =   18: ' + '\n')
    file.write('Number of fitting parameters =     4 for each band' + '\n')
    file.write('Number of terms (Phase fnc.) =   1' + '\n')
    for i in range(spectrum.n_band):
        file.write('Band =     ' + str(i+1) + '     Fitting Parameters:' + '\n')
        file.write('    ' + "{:.5e}".format(params.kext[i, 0]))
        file.write('     ' + "{:.5e}".format(params.ksca[i, 0]))
        file.write('     ' + "{:.5e}".format(params.g[i, 0]))
        file.write('     ' + "{:.5e}".format(params.g[i, 1]) + '\n')
    file.close()

def write_spectral_region(file, region, spectral_data, params):
    '''
    Write parameters for given spectral region to file.
    '''
    file.write(region+':\n\n')
    file.write('           Wavelength (m)         |        Kext      |      Ksca        |               g                 |\n')
    file.write('----------------------------------------------------------------------------------------------------------|\n')
    for i in range(spectral_data.n_band):
        file.write("{:.9e}".format(spectral_data.lower_bound[i]) + ' - ' + "{:.9e}".format(spectral_data.upper_bound[i]))
        file.write(' | ')
        file.write("{:.5e}".format(params.kext[i, 0]) + ' * qi')
        file.write(' | ')
        file.write("{:.5e}".format(params.ksca[i, 0]) + ' * qi')
        file.write(' | ')
        file.write("{:.5e}".format(params.g[i, 0]) + ' * (qi ^ ' + "{:.5e}".format(params.g[i, 1]))
        file.write(' |\n')
 

def rescale(SW_params, LW_params):
    '''
    Detemine whether rescaling of optical properties is required.

    Ask user whether resclaing is required and if yes request a
    scaling factor. Use this to rescale kext and ksca parameters.
    '''
    success = False
    while success == False:
        rescale = input("Do you wish to rescale the extinction to represent a different distribution of ice crystal shapes? (Y/N)\n")
        print("You entered '" + rescale + "'")
        success = rescale in ['y', 'Y', 'yes', 'n', 'N', 'no']
    if rescale in ['y', 'Y', 'yes']:
        while True:
            try:
                factor = input("Enter the scaling factor\n")
                factor = float(factor)
                SW_params.kext = factor * SW_params.kext
                SW_params.ksca = factor * SW_params.ksca
                LW_params.kext = factor * LW_params.kext
                LW_params.ksca = factor * LW_params.ksca
                break
            except ValueError:
                print("'" + factor + "' is not a valid number.")
    return SW_params, LW_params

if __name__ == "__main__":
#   Read spectral files and ice optical properties
    file = filename()
    SW_spectra = SpecData(file.sw_spec)
    LW_spectra = SpecData(file.lw_spec)
    database = IceData(file.ss_ice)
    SW = identify_region(SW_spectra, database)
    LW = identify_region(LW_spectra, database)
    SW_params = parameters(SW_spectra.n_band)
    LW_params = parameters(LW_spectra.n_band)

#   Fit ice optical properties.
    SW_params.kext = fit_opt_prop(database.wavelength[SW], database.qi[SW], database.kext[SW], SW_spectra, 'SW kext', 'linear_qi', True, weighting='Spectrum')
    LW_params.kext = fit_opt_prop(database.wavelength[LW], database.qi[LW], database.kext[LW], LW_spectra, 'LW kext', 'linear_qi', True, weighting='Planck')
    SW_params.ksca = fit_opt_prop(database.wavelength[SW], database.qi[SW], database.ksca[SW], SW_spectra, 'SW ksca', 'linear_qi', True, weighting='Spectrum')
    LW_params.ksca = fit_opt_prop(database.wavelength[LW], database.qi[LW], database.ksca[LW], LW_spectra, 'LW ksca', 'linear_qi', True, weighting='Planck')
    SW_params.g = fit_opt_prop(database.wavelength[SW], database.qi[SW], database.g[SW], SW_spectra, 'SW g', 'power_qi', True, weighting='Spectrum')
    LW_params.g = fit_opt_prop(database.wavelength[LW], database.qi[LW], database.g[LW], LW_spectra, 'LW g', 'power_qi', True, weighting='Planck')

#   Rescale if required
    SW_params, LW_params = rescale(SW_params, LW_params)

#   Write parameters to file in reader friendly format
    write_opt_prop_param(file.ss_ice + '_parametrization.txt', file.ss_ice, file.sw_spec, file.lw_spec, SW_spectra, LW_spectra, SW_params, LW_params)

#   Write parameters to file in format to be used by prep_spec
    write_opt_prop_prep_spec(file.ss_ice + '_sw_type10', SW_spectra, SW_params)
    write_opt_prop_prep_spec(file.ss_ice + '_lw_type10', LW_spectra, LW_params)


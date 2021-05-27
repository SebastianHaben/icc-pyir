"""Functions to fit data"""

import numpy as np
import lmfit


def gaussian(x, amplitude, mean, width):
    y = amplitude * np.exp(-np.log(2) * ((x - mean) / width) ** 2)
    y = np.nan_to_num(y)
    return y


def gaussian_fitting(wavenumbers, data, band, band_height):
    model = lmfit.Model(gaussian)
    par = lmfit.Parameters()
    par.add(lmfit.Parameter('amplitude', value=band_height, vary=False))
    par.add(lmfit.Parameter('mean', value=band, vary=False))
    par.add(lmfit.Parameter('width', value=5))
    out = model.fit(data, par, x=wavenumbers)
    return out


def create_lin_baseline_between_endpoints(wavenumbers, data):
    start_point = (wavenumbers[0], data[0])
    end_point = (wavenumbers[-1], data[-1])
    slope = (end_point[1]-start_point[1])/(end_point[0]-start_point[0])
    c = start_point[1]-slope*start_point[0]
    baseline = np.array([linear(x, slope, c) for x in wavenumbers])
    corrected_data = data-baseline
    return corrected_data, baseline


def linear(x, slope, c):
    return slope*x+c

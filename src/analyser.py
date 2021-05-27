"""Analyser"""

# Copyright (c) 2020
# Authors: Sebastian Haben
# License: MIT

import datetime
import os
import time
import warnings

from scipy import signal, integrate

from fitting import gaussian, gaussian_fitting, create_lin_baseline_between_endpoints
from errors import NormalizationError, WavenumberArrayShapeMissMatch, PositionWarning, WrongFileType, \
    NoPeaksFoundError, ExperimentNoMatch
from reader import SPAReader, CSVReader, np
from requests import Request


class Analyser:
    """Class to analyse pyridine IR experiments

        Attributes:
        -----------
        pre_adsorption_request : utils.requests.Request
            Request object of the file belonging to the pre pyridine adsorption file.
        post_adsorption_request : utils.requests.Request
            Request object of the file belonging to the post pyridine adsorption file.
        date_measured : datetime.datetime
            Date when the sample was measured.
        material_id : str
            Id of the material.
        pre_adsorption_wavenumbers : numpy.ndarray
            Array of the pre adsorption spectrum wavenumbers.
        pre_adsorption_data : numpy.ndarray
            Array of the pre adsorption spectrum intensities.
        post_adsorption_wavenumbers : numpy.ndarray
            Array of the post adsorption spectrum wavenumbers.
        post_adsorption_data : numpy.ndarray
            Array of the post adsorption spectrum intensities.
        normalization_band_guess : float
            Guess for the position of the normalization band.
        normalization_bands : dict of floats
            Dictionary containing the actual positions of the pre and post adsorptions spectra. Keys: 'pre' and 'post'.
        norm_pre_adsorption_data : numpy.ndarray
            Array of the normalized pre adsorption spectrum intensities.
        norm_post_adsorption_data : numpy.ndarray
            Array of the normalized post adsorption spectrum intensities.
        self.resolution : float
            Spectral resolution in cm-1.
        difference_spectra_data : numpy.ndarray
            Array of the difference spectrum intensities.
        difference_spectra_wavenumbers : numpy.ndarray
            Array of the difference spectrum wavenumbers.
        self.py_peaks : list
            List of peak positions.
        rois : dict
            Dictionary of ROIs sorted by peak positions.
        py_peaks_fitted : dict
            Dictionary of fitting and integration results sorted by peak positions.
        bas : dict
            Dictionary including all fitting, integration and acid sites data about the BAS.
        las : dict
            Dictionary including all fitting, integration and acid sites data about the LAS.
        sample_weight : float
            Sample weight in grams.
        git_sha : str
            Git SHA-key of the project.
        short_git_sha : str
            First 8 characters of the Git SHA-key.
        sample_disk_radius : float
                Sample disk radius in cm.

        Methods
        -------
        add_pre_adsorption_spectrum(pre_adsorption_request):
            Setter function for the pre_adsorption_request attribute.
        add_post_adsorption_spectrum(post_adsorption_request):
            Setter function for the post_adsorption_request attribute.
        normalize_spectra(normalization_band_guess):
            Normalization of the pre and post pyridine adsorption spectrum.
        calculate_difference_spectra():
            Calculates differences between pre and post pyridine adsorption spectra.
        find_py_adsorption_peaks(py_adsorption_band_guess, search_window_size):
            Detects the bands originating in pyridine adsorption.
        fit_py_adsorption_peaks():
           Fits the individual ROIs data to individual Gauss functions that can be integrated.
        integrate_fitted_peaks():
            Integrates the fitted pyridine peaks
        calculate_acid_sites(sample_weight)
            Calculation of BAS and LAS acid site density in µmol/g.
        """
    def __init__(self):
        self.pre_adsorption_request = None
        self.post_adsorption_request = None
        self.date_measured = None
        self.pre_adsorption_wavenumbers = None
        self.pre_adsorption_data = None
        self.post_adsorption_wavenumbers = None
        self.post_adsorption_data = None
        self.normalization_band_guess = None
        self.normalization_bands = {}
        self.norm_pre_adsorption_data = None
        self.norm_post_adsorption_data = None
        self.resolution = None
        self.difference_spectra_data = None
        self.difference_spectra_wavenumbers = None
        self.py_peaks = None
        self.rois = None
        self.py_peaks_fitted = {}
        self.bas = None
        self.las = None
        self.sample_weight = None
        self.sample_disk_radius = None

    def add_pre_adsorption_spectrum(self, pre_adsorption_request_uri):
        """Load data from pre pyridine adsorption request.

            Parameters:
            -----------
            pre_adsorption_request_uri : requests.Request or str
                Request object or of the pre pyridine adsorption file or path to the file.

            Raises:
            -------
            TypeError:
                If neither a Request object nor a file path was provided.
            """
        if isinstance(pre_adsorption_request_uri, Request):
            pre_adsorption_request = pre_adsorption_request_uri
        elif os.path.exists(pre_adsorption_request_uri):
            pre_adsorption_request = Request(pre_adsorption_request_uri)
        else:
            raise TypeError("Please either provide a Request object or a file path")
        self.pre_adsorption_wavenumbers, self.pre_adsorption_data = self._get_reader_get_data(pre_adsorption_request)
        self.pre_adsorption_request = pre_adsorption_request
        self.date_measured = datetime.datetime.strptime(time.ctime(
            os.path.getmtime(pre_adsorption_request.abs_file_path)), "%a %b %d %H:%M:%S %Y")

    def add_post_adsorption_spectrum(self, post_adsorption_request_uri):
        """Load data from pre pyridine adsorption file.

            Parameters:
            -----------
            post_adsorption_request_uri : requests.Request
                Request object of the post pyridine adsorption file or path to the file.

            Raises:
            -------
            TypeError:
                If neither a Request object nor a file path was provided.
            Warns:
            ------
            UserWarning:
                If date of measurement is not identical between pre- and post adsorption file.
            """
        if isinstance(post_adsorption_request_uri, Request):
            post_adsorption_request = post_adsorption_request_uri
        elif os.path.exists(post_adsorption_request_uri):
            post_adsorption_request = Request(post_adsorption_request_uri)
        else:
            raise TypeError("Please either provide a Request object or a file path")
        date_measured = datetime.datetime.strptime(time.ctime(
            os.path.getmtime(post_adsorption_request.abs_file_path)), "%a %b %d %H:%M:%S %Y").date()
        if date_measured != self.date_measured.date():
            warnings.warn('The measuremnt day for files is not identical!', UserWarning)
        self.post_adsorption_wavenumbers, self.post_adsorption_data = self._get_reader_get_data(post_adsorption_request)
        self.post_adsorption_request = post_adsorption_request

    def normalize_spectra(self, normalization_band_guess):
        """Normalization of spectra.

            Parameters:
            -----------
            normalization_band_guess : int
                Wavenumber position of vibration band which is used for normalization."""
        self.normalization_band_guess = normalization_band_guess
        self.normalization_bands['pre_adsorption'], self.norm_pre_adsorption_data = \
            self._normalize_data(normalization_band_guess, self.pre_adsorption_wavenumbers, self.pre_adsorption_data)
        self.normalization_bands['post_adsorption'], self.norm_post_adsorption_data = \
            self._normalize_data(normalization_band_guess, self.post_adsorption_wavenumbers, self.post_adsorption_data)

    def calculate_difference_spectra(self, normalization=False):
        """Calculates differences between pre and post pyridine adsorption spectra.

            Parameters:
            -----------
            normalization : bool
                Indicate if normalized data should be used.

            Raises:
            -------
            WavenumberArraysShapeMissMatch:
                If the shape of the pre- and post-adsorption spectra aren't equal.

            Warnings:
            ---------
            PositionWarning:
                If wavenumbers of pre- and post-adsorption spectra aren't equal."""
        if self.pre_adsorption_wavenumbers.size == self.post_adsorption_wavenumbers.size:
            if np.array_equal(self.pre_adsorption_wavenumbers, self.post_adsorption_wavenumbers):
                if normalization:
                    self.difference_spectra_data = self.norm_post_adsorption_data - self.norm_pre_adsorption_data
                else:
                    self.difference_spectra_data = self.post_adsorption_data - self.pre_adsorption_data
                self.difference_spectra_wavenumbers = self.pre_adsorption_wavenumbers
                self.resolution = abs(np.diff(self.difference_spectra_wavenumbers).mean())
            else:
                warnings.warn(PositionWarning(), category='UserWarning')
                if normalization:
                    self.difference_spectra_data = self.norm_post_adsorption_data - self.norm_pre_adsorption_data
                else:
                    self.difference_spectra_data = self.post_adsorption_data - self.pre_adsorption_data
                self.difference_spectra_wavenumbers = np.mean([self.pre_adsorption_wavenumbers,
                                                               self.post_adsorption_wavenumbers], axis=0)
                self.resolution = abs(np.diff(self.difference_spectra_wavenumbers).mean())
        else:
            raise WavenumberArrayShapeMissMatch()

    def find_py_adsorption_peaks(self, py_adsorption_band_guess, search_window_size):
        """Detects the bands originating in pyridine adsorption.

            Parameters:
            -----------
            py_adsorption_band_guess : list of float
                List of band positions as floats.
            search_window_size : list of floats
                List of values of how big the search area around each peak is as floats. The index of the window size is
                identical to the index of the band guess in py_adsorption_band_guess
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample.

            Raises:
            -------
            NoPeaksFoundError:
                If no Peaks could be detected.
                """

        rois = {}
        rois_data_list = []
        py_peaks = []
        for band in py_adsorption_band_guess:
            band_index = py_adsorption_band_guess.index(band)
            roi_index = (np.abs(self.difference_spectra_wavenumbers-band)).argmin()
            window_index_size = round(search_window_size[band_index]/self.resolution)
            roi_index_low = int(roi_index - window_index_size/2)
            roi_index_high = int(roi_index + window_index_size/2)
            roi_data = self.difference_spectra_data[roi_index_low:roi_index_high]
            rois_data_list.append({'wavenumbers': self.difference_spectra_wavenumbers[roi_index_low:roi_index_high],
                                   'data': roi_data})
        for roi in rois_data_list:
            peaks = signal.find_peaks(roi['data'], prominence=0.1)
            try:
                peak = roi['wavenumbers'][peaks[0][0]]
            # Try to find peaks with lower prominence
            except IndexError:
                peaks = signal.find_peaks(roi['data'], prominence=0.01)
                try:
                    peak = roi['wavenumbers'][peaks[0][0]]
                except IndexError:
                    raise NoPeaksFoundError('No peaks were found!')
            py_peaks.append(peak)
            rois[peak] = roi
        self.rois = rois
        self.py_peaks = py_peaks

    def fit_py_adsorption_peaks(self):
        """Fits the individual ROIs data to individual Gauss functions that can be integrated."""
        for band, roi in self.rois.items():
            band_index = list(roi['wavenumbers']).index(band)
            corrected_data, baseline = create_lin_baseline_between_endpoints(roi['wavenumbers'], roi['data'])
            corrected_band_height = corrected_data[band_index]
            result = gaussian_fitting(roi['wavenumbers'], corrected_data, band, corrected_band_height)
            fitting_dict = {'wavenumbers': roi['wavenumbers'], 'data': roi['data'], 'baseline': baseline,
                            'corrected data': corrected_data, 'results': result}
            self.py_peaks_fitted[band] = fitting_dict

    def integrate_fitted_peaks(self):
        """Integrates the fitted pyridine peaks"""
        for results in self.py_peaks_fitted.values():
            area = integrate.quad(gaussian, results['wavenumbers'].min(), results['wavenumbers'].max(),
                                  args=tuple(results['results'].best_values.values()))
            results['area'] = area

    def calculate_acid_sites(self, sample_weight, sample_disk_radius=1.287/2):
        """Calculation of BAS and LAS acid site density in µmol/g.

            Parameters:
            ----------
            sample_weight : float
                Sample weight in mg.
            sample_disk_radius : float
                Sample disk radius in cm. Default 0.6435 cm from standard waver stamps.

            References:
            -----------------
            Calculation is based on:
                https://doi.org/10.1016/0021-9517(63)90102-7.
            """
        self.sample_weight = sample_weight
        self.sample_disk_radius = sample_disk_radius
        for band, results in self.py_peaks_fitted.items():
            if abs(band-1455) < 10:
                results['acid sites'] = 1.88*results['area'][0]*(self.sample_disk_radius**2 / sample_weight)
                self.las = results
                self.las['position'] = band
            elif abs(band-1546) < 10:
                results['acid sites'] = 1.42*results['area'][0]*(self.sample_disk_radius**2 / sample_weight)
                self.bas = results
                self.bas['position'] = band

    @staticmethod
    def _normalize_data(normalization_band_guess, wavenumber_array, data):
        """Normalizes the data around a certain band.

            Parameters:
            -----------
            normalization_band_guess : int
                Wavenumber position of vibration band which is used for normalization.
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample.

            Returns:
            ---------
            Normalized Data
        """
        peaks = signal.find_peaks(data, prominence=0.1)
        index_guess = (np.abs(wavenumber_array[peaks[0]]-normalization_band_guess)).argmin()
        calculated_band_position = wavenumber_array[peaks[0][index_guess]]
        if abs(normalization_band_guess-calculated_band_position) < 10:
            normalized_data = data/data[index_guess]
            return calculated_band_position, normalized_data
        else:
            peaks = signal.find_peaks(data, prominence=0.01)
            index_guess = (np.abs(wavenumber_array[peaks[0]] - normalization_band_guess)).argmin()
            calculated_band_position = wavenumber_array[index_guess]
            if abs(normalization_band_guess - calculated_band_position) < 10:
                normalized_data = data / data[index_guess]
                return calculated_band_position, normalized_data
            else:
                raise NormalizationError(wavenumber_array, data, peaks)

    @staticmethod
    def _get_reader(request):
        """Finds suitable reader for the requested file

            Parameters:
            ----------
            request : requests.Request
                Request object of the file.

            Returns:
            --------
            reader : reader.Reader subclass
                Suitable reader object for the file type."""

        if request.extension in ['.spa', '.SPA', '.srs']:
            reader = SPAReader(request)
        elif request.extension in ['.csv', '.CSV']:
            reader = CSVReader(request)
        else:
            raise WrongFileType(request.extension)
        return reader

    def _get_reader_get_data(self, request):
        """Gets the data from the reader.

            Parameters:
            ----------
            request : requests.Request
                Request object of the file.

            Returns:
            --------
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample."""
        reader = self._get_reader(request)
        wavenumber_array, data = reader.read_data()
        return wavenumber_array, data

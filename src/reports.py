import datetime
import os
from xml.etree.ElementTree import Element, SubElement, Comment, ElementTree

import jinja2
from matplotlib import pyplot as plt
import numpy as np
from xhtml2pdf import pisa

from utils.errors import ExperimentNoMatch, ExperimentMissMatch


class XMLReport:
    """Class to create xml report

        Parameters:
        -----------
        analyser : utils.analyser.Analyser
            Analyser object of the experiment.
        date : datetime.date
            Date when the experiment was performed.
        experiment_id : str
            ID of the experiment.
        short_git_sha : str
            First 8 characters of the Git SHA-key.

        Attributes:
        -----------
        analysed_dt : datetime.date
            Date when the report was created.
        analyser : utils.analyser.Analyser
            Analyser object of the experiment.
        experiment_id : str
            ID of the experiment.
        recorded_dt : datetime.date
            Date when the experiment was performed.
        root : xml.etree.ElementTree.Element
            Root element of the xml structure.
        tree : xml.etree.ElementTree.ElementTree
            ElementTree object of the xml structure.
        peak_section : bol
            Indication if the peak section variable is available.
        fitting_section : bol
            Indication if the fitting section variable is available.
        short_git_sha : str
            First 8 characters of the Git SHA-key.

        Methods:
        -------
        build_report(add_raw_data=True, add_normalization_data=True):
            Builds xml file containing the report.
        save_report(file_path):
            Saves the report as .xml file.
        get_experiment_id(experiment_id_pattern, automatic_search=True, experiment_id=None):
            Automatically finds experiment ID based on automatic search in file name ore uses provided experiment_id.
        """
    def __init__(self, analyser, short_git_sha='No version specified'):
        self.analysed_dt = datetime.datetime.now()
        self.analyser = analyser
        self.material_id = analyser.material_id
        self.recorded_dt = analyser.date_measured
        self.root = None
        self.tree = None
        self.peak_section = False
        self.fitting_section = False
        self.short_git_sha = short_git_sha

    def build_report(self, add_raw_data=True, add_normalization_data=True):
        """Builds xml file containing the report.

            Parameters:
            -----------
            add_raw_data : bool
                Should the raw data of the analysed spectra be included. Default: True.
            add_normalization_data : bool
                Should the normalized data of the analysed spectra be included. Default: True.
            """
        root = Element('report')
        root.set('version', self.short_git_sha)
        root.append(Comment('Created by pyir.'))
        tree = ElementTree(root)
        # General Information section
        general_info_section = SubElement(root, 'section', {'sectionName': 'General Information'})
        material_id = SubElement(general_info_section, 'materialID')
        material_id.text = self.material_id
        sample_weight = SubElement(general_info_section, 'sampleWeight', {'unit': 'mg'})
        sample_weight.text = str(self.analyser.sample_weight)
        sample_radius = SubElement(general_info_section, 'sampleDiskRadius', {'unit': 'cm'})
        sample_radius.text = str(self.analyser.sample_disk_radius)
        recorded_date = SubElement(general_info_section, 'dateRecorded', {'format': 'DD.MM.YYYY'})
        recorded_date.text = self.recorded_dt.strftime('%d.%m.%Y')
        analysed_date = SubElement(general_info_section, 'dateAnalysed', {'format': 'DD.MM.YYYY'})
        analysed_date.text = self.analysed_dt.strftime('%d.%m.%Y')
        try:
            normalization_bands = SubElement(general_info_section, 'normalizationBands',
                                             {'initialBandGuess': str(self.analyser.normalization_band_guess)})
            pre_adsorption_norm_band = SubElement(normalization_bands, 'normalizationBand', {'type': 'pre adsorption'})
            pre_adsorption_norm_band.text = str(self.analyser.normalization_bands['pre_adsorption'])
            post_adsorption_norm_band = \
                SubElement(normalization_bands, 'normalizationBand', {'type': 'post adsorption'})
            post_adsorption_norm_band.text = str(self.analyser.normalization_bands['post_adsorption'])
        except KeyError:
            pass
        bas_sites = SubElement(general_info_section, 'bas', {'unit': 'mmol/g'})
        bas_sites.text = str(self.analyser.bas['acid sites'])
        las_sites = SubElement(general_info_section, 'las', {'unit': 'mmol/g'})
        las_sites.text = str(self.analyser.las['acid sites'])
        # Raw data section
        raw_data_section = SubElement(root, 'section', {'sectionName': 'Raw Data'})
        pre_adsorption_raw_section = SubElement(raw_data_section, 'rawData', {'type': 'pre adsorption'})
        pre_adsorption_file_name = SubElement(pre_adsorption_raw_section, 'filename', {'type': 'pre adsorption'})
        pre_adsorption_file_name.text = self.analyser.pre_adsorption_request.abs_file_path
        post_adsorption_raw_section = SubElement(raw_data_section, 'rawData', {'type': 'post adsorption'})
        post_adsorption_file_name = SubElement(post_adsorption_raw_section, 'filename', {'type': 'post adsorption'})
        post_adsorption_file_name.text = self.analyser.post_adsorption_request.abs_file_path
        if add_raw_data:
            pre_adsorption_raw_spectrum_section = SubElement(pre_adsorption_raw_section, 'spectrum',
                                                             {'type': 'pre adsorption', 'dataType': 'raw data'})
            pre_adsorption_raw_spectrum_wavenumbers = SubElement(pre_adsorption_raw_spectrum_section, 'wavenumbers',
                                                                 {'type': 'pre adsorption', 'dataType': 'raw data',
                                                                  'unit': 'cm-1'})
            self._generate_data_points(pre_adsorption_raw_spectrum_wavenumbers,
                                       self.analyser.pre_adsorption_wavenumbers)
            pre_adsorption_raw_spectrum_intensity = SubElement(pre_adsorption_raw_spectrum_section, 'intensities',
                                                               {'type': 'pre adsorption', 'dataType': 'raw data',
                                                                'unit': 'a.u.'})
            self._generate_data_points(pre_adsorption_raw_spectrum_intensity, self.analyser.pre_adsorption_data)
            post_adsorption_raw_spectrum_section = SubElement(post_adsorption_raw_section, 'spectrum',
                                                              {'type': 'post adsorption', 'dataType': 'raw data'})
            post_adsorption_raw_spectrum_wavenumbers = SubElement(post_adsorption_raw_spectrum_section, 'wavenumbers',
                                                                  {'type': 'post adsorption', 'dataType': 'raw data',
                                                                   'unit': 'cm-1'})
            self._generate_data_points(post_adsorption_raw_spectrum_wavenumbers,
                                       self.analyser.post_adsorption_wavenumbers)
            post_adsorption_raw_spectrum_intensity = SubElement(post_adsorption_raw_spectrum_section, 'intensities',
                                                                {'type': 'post adsorption', 'dataType': 'raw data',
                                                                 'unit': 'a.u.'})
            self._generate_data_points(post_adsorption_raw_spectrum_intensity,
                                       self.analyser.post_adsorption_data)
        # Normalization section
        if add_normalization_data:
            normalization_section = SubElement(root, 'section', {'sectionName': 'Normalized Data'})
            pre_adsorption_norm_section = SubElement(normalization_section, 'normalizedData',
                                                     {'type': 'pre adsorption'})
            pre_adsorption_norm_spectrum_section = SubElement(pre_adsorption_norm_section, 'spectrum',
                                                              {'type': 'pre adsorption', 'dataType': 'normalized data'})
            pre_adsorption_norm_spectrum_intensity = SubElement(pre_adsorption_norm_spectrum_section, 'intensities',
                                                                {'type': 'pre adsorption',
                                                                 'dataType': 'normalized data',
                                                                 'unit': 'a.u.'})
            self._generate_data_points(pre_adsorption_norm_spectrum_intensity, self.analyser.norm_pre_adsorption_data)
            post_adsorption_norm_section = SubElement(normalization_section, 'normalizedData',
                                                      {'type': 'post adsorption'})
            post_adsorption_norm_spectrum_section = SubElement(post_adsorption_norm_section, 'spectrum',
                                                               {'type': 'post adsorption', 'dataType': 'raw data'})
            post_adsorption_norm_spectrum_intensity = SubElement(post_adsorption_norm_spectrum_section, 'wavenumbers',
                                                                 {'type': 'post adsorption',
                                                                  'dataType': 'normalized data', 'unit': 'a.u.'})
            self._generate_data_points(post_adsorption_norm_spectrum_intensity,
                                       self.analyser.norm_post_adsorption_data)
        # Difference spectrum section
        diff_data_section = SubElement(root, 'section', {'sectionName': 'Difference spectrum data'})
        diff_spectrum_section = SubElement(diff_data_section, 'spectrum', {'dataType': 'difference spectrum data'})
        diff_spectrum_wavenumbers = SubElement(diff_spectrum_section, 'wavenumbers',
                                               {'dataType': 'difference spectrum data', 'unit': 'cm-1'})
        self._generate_data_points(diff_spectrum_wavenumbers,
                                   self.analyser.difference_spectra_wavenumbers)
        diff_spectrum_intensity = SubElement(diff_spectrum_section, 'intensities',
                                             {'dataType': 'difference spectrum data', 'unit': 'a.u.'})
        self._generate_data_points(diff_spectrum_intensity, self.analyser.difference_spectra_data)
        # Peaks section
        peak_data_section = SubElement(root, 'section', {'sectionName': 'Peaks data'})
        bas_section = SubElement(peak_data_section, 'bas',
                                 {'peakPosition': str(round(self.analyser.bas['position'], 2))})
        bas_peak_data = SubElement(bas_section, 'peakData', {'dataType': 'peak data'})
        bas_peak_wavenumbers = SubElement(bas_peak_data, 'wavenumbers', {'dataType': 'peak data', 'unit': 'cm-1'})
        self._generate_data_points(bas_peak_wavenumbers, self.analyser.bas['wavenumbers'])
        bas_peak_intensities = SubElement(bas_peak_data, 'intensities', {'dataType': 'peak data', 'unit': 'a.u.'})
        self._generate_data_points(bas_peak_intensities, self.analyser.bas['data'])
        las_section = SubElement(peak_data_section, 'las',
                                 {'peakPosition': str(round(self.analyser.las['position'], 2))})
        las_peak_data = SubElement(las_section, 'peakData', {'dataType': 'peak data'})
        las_peak_wavenumbers = SubElement(las_peak_data, 'wavenumbers', {'dataType': 'peak data', 'unit': 'cm-1'})
        self._generate_data_points(las_peak_wavenumbers, self.analyser.las['wavenumbers'])
        las_peak_intensities = SubElement(las_peak_data, 'intensities', {'dataType': 'peak data', 'unit': 'a.u.'})
        self._generate_data_points(las_peak_intensities, self.analyser.las['data'])
        peak_counter = 1
        try:
            excess_bands = (set(list(self.analyser.py_peaks_fitted.keys())) - {self.analyser.las['position'],
                                                                               self.analyser.bas['position']})
            if len(excess_bands) > 0:
                peak_section = SubElement(peak_data_section, 'peakSection')
                for band, values in self.analyser.py_peaks_fitted.items():
                    if band not in [self.analyser.las['position'], self.analyser.bas['position']]:
                        self._add_peak_results(band, values, peak_section, peak_counter)
                        peak_counter += 1

        except (KeyError, TypeError):
            peak_section = SubElement(peak_data_section, 'peakSection')
            for band, values in self.analyser.py_peaks_fitted.items():
                self._add_peak_results(band, values, peak_section, peak_counter)
                peak_counter += 1
        # Fitting section
        fitting_section = SubElement(root, 'section', {'sectionName': 'Fitting data'})
        bas_fitting_section = SubElement(fitting_section, 'bas',
                                         {'peakPosition': str(round(self.analyser.bas['position'], 2))})
        self._add_fitting_results(self.analyser.bas, bas_fitting_section)
        las_fitting_section = SubElement(fitting_section, 'las',
                                         {'peakPosition': str(round(self.analyser.las['position'], 2))})
        self._add_fitting_results(self.analyser.las, las_fitting_section)
        peak_counter = 1
        try:
            excess_bands = (set(list(self.analyser.py_peaks_fitted.keys())) - {self.analyser.las['position'],
                                                                               self.analyser.bas['position']})
            if len(excess_bands) > 0:
                peak_fitting_section = SubElement(fitting_section, 'peakSection')
                for band, values in self.analyser.py_peaks_fitted.items():
                    peak_data = peak_fitting_section.makeelement('peak' + str(peak_counter),
                                                                 {'dataType': 'fitting data',
                                                                  'peakPosition': str(round(band, 2))})
                    peak_fitting_section.append(peak_data)
                    self._add_fitting_results(values, peak_fitting_section)
                    peak_counter += 1
        except (KeyError, TypeError):
            peak_fitting_section = SubElement(fitting_section, 'peakSection')
            for band, values in self.analyser.py_peaks_fitted.items():
                peak_data = peak_fitting_section.makeelement('peak' + str(peak_counter),
                                                             {'dataType': 'fitting data',
                                                              'peakPosition': str(round(band, 2))})
                peak_fitting_section.append(peak_data)
                self._add_fitting_results(values, peak_fitting_section)
                peak_counter += 1
        self.root = root
        self.tree = tree

    def save_report(self):
        """Saves the report as .xml file"""
        if self.root is None:
            self.build_report()
        self._indent(self.root)
        dir_path = os.path.split(self.analyser.pre_adsorption_request.abs_file_path)[0]
        file_path = os.path.join(dir_path, '_'.join([self.material_id, 'report', self.short_git_sha])+'.xml')
        self.tree.write(file_path, encoding='utf-8', xml_declaration=True)

    def get_experiment_id(self, experiment_id_pattern, automatic_search=True, experiment_id=None):
        """Automatically finds experiment ID based on automatic search in file name ore uses provided experiment_id.

            Parameters:
            -----------
            experiment_id_pattern : re.Pattern
                RegEx to describe the file pattern.
            automatic_search : bool
                Turn automatic search on or off. Default: True.
            experiment_id : str
                Name of experiment or material."""
        if automatic_search:
            try:
                pre_experiment_id = experiment_id_pattern.search(self.analyser.pre_adsorption_request.abs_file_path)[0]
                post_experiment_id = experiment_id_pattern.search(
                    self.analyser.post_adsorption_request.abs_file_path)[0]
            except IndexError:
                raise ExperimentNoMatch()
            if pre_experiment_id is not post_experiment_id:
                raise ExperimentMissMatch(pre_experiment_id, post_experiment_id)
            self.material_id = pre_experiment_id
        elif experiment_id is not None and not automatic_search:
            self.material_id = experiment_id

    def _add_peak_results(self, band, values, peak_section, peak_counter):
        """Helper function to create peak section of xml report.

            Parameters:
            ----------
            band : float
                Position of peak in cm-1.
            values : dict
                Dictionary containing all information about peak finding and fitting.
                Created by utils.analyser.Analyser.
            peak_data_section : xml.etree.ElementTree.SubElement
                Section of the xml report to which the data belongs to.
            peak_counter = int
                Counter variable to sequentially label the peaks."""
        peak_data = SubElement(peak_section, 'peak' + str(peak_counter), {'dataType': 'peak data',
                                                                          'peakPosition': str(round(band, 2))})
        peak_wavenumbers = SubElement(peak_data, 'wavenumbers', {'dataType': 'peak data', 'unit': 'cm-1'})
        self._generate_data_points(peak_wavenumbers, values['wavenumbers'])
        peak_intensities = SubElement(peak_data, 'intensities', {'dataType': 'peak data', 'unit': 'a.u.'})
        self._generate_data_points(peak_intensities, values['data'])

    def _add_fitting_results(self, values, section):
        """Helper function to add fitting results of known acid sites.

            Parameters:
            -----------
            values : dict
                Dictionary containing all information about peak finding and fitting.
                Created by utils.analyser.Analyser.
            section : xml.etree.ElementTree.SubElement
                Section of the xml report to which the data belongs to.
            """
        results = values['results']
        fitting_general = section.makeelement('fittingGeneral', {'method': results.method,
                                                                 'dataPoints': str(results.ndata),
                                                                 'iterations': str(results.nfev)})
        section.append(fitting_general)
        fitting_quality = section.makeelement('fittingQuality', {'chiSqr': str(results.chisqr),
                                                                 'redChiSqr': str(results.redchi)})
        section.append(fitting_quality)
        fitting_parameters = section.makeelement('fittingParameters', {'initial': str(results.init_values),
                                                                       'best': str(results.best_values)})
        section.append(fitting_parameters)
        fitting_area = SubElement(section, 'area', {'stdDev': str(values['area'][1])})
        fitting_area.text = str(values['area'][0])
        fitting_data = SubElement(section, 'fittingData', {'dataType': 'fitting data'})
        fitting_wavenumbers = SubElement(fitting_data, 'wavenumbers', {'dataType': 'fitting data', 'unit': 'cm-1'})
        self._generate_data_points(fitting_wavenumbers, values['wavenumbers'])
        fitting_intensities = SubElement(fitting_data, 'intensities', {'dataType': 'fitting data', 'unit': 'a.u.'})
        self._generate_data_points(fitting_intensities, results.best_fit)

    @staticmethod
    def _generate_data_points(parent_element, arr):
        """"Function to create transfer arrays to xml.

            Parameters:
            -----------
            parent_element : xml.etree.ElementTree.SubElement
                Parent element where the array will be nested into.
            arr : np.array
                Array containing the data."""
        for data in list(arr):
            datapoint = SubElement(parent_element, 'datapoint')
            datapoint.text = str(data)

    def _indent(self, elem, level=0):
        """Helper function to create new lines and indents for better readability.

            Parameters:
            -----------
            elements : xml.etree.ElementTree.Element
                Root of the xml file.

            References:
            -----------
            Code copied from:
            https://stackoverflow.com/questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree
            -elementtree-in-python
        """
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self._indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i


class PDFReport:
    """Class to create pdf report

        Parameters:
        -----------
        analyser : utils.analyser.Analyser
            Analyser object of the experiment.
        date : datetime.date
            Date when the experiment was performed.
        experiment_id : str
            ID of the experiment.
        short_git_sha : str
            First 8 characters of the Git SHA-key.

        Attributes:
        -----------
        analysed_dt : datetime.date
            Date when the report was created.
        analyser : utils.analyser.Analyser
            Analyser object of the experiment.
        experiment_id : str
            ID of the experiment.
        recorded_dt : datetime.date
            Date when the experiment was performed.
        root : xml.etree.ElementTree.Element
            Root element of the xml structure.
        tree : xml.etree.ElementTree.ElementTree
            ElementTree object of the xml structure.
        peak_section : bol
            Indication if the peak section variable is available.
        fitting_section : bol
            Indication if the fitting section variable is available.
        short_git_sha : str
            First 8 characters of the Git SHA-key.

        Methods:
        -------
        build_report(add_raw_data=True, add_normalization_data=True):
            Builds xml file containing the report.
        save_report(file_path):
            Saves the report as .xml file.
        get_experiment_id(experiment_id_pattern, automatic_search=True, experiment_id=None):
            Automatically finds experiment ID based on automatic search in file name ore uses provided experiment_id.
        """

    def __init__(self, analyser, short_git_sha='No version specified'):
        self.analysed_dt = datetime.datetime.now()
        self.analyser = analyser
        self.material_id = analyser.material_id
        self.recorded_dt = analyser.date_measured
        self.html = None
        self.short_git_sha = short_git_sha

    def build_report(self, include_norm=False):
        """Builds PDF report.

            Parameters:
            -----------
            include_norm : bool
                Indicate if normalized data should be included. Default = False."""
        self._plot_data(include_norm)

        raw_data_plot = {'plot': os.path.abspath('.\\images\\raw.png'), 'caption': 'Measured FT-IR spectrum pre- and '
                                                                                   'post-pyridine adsorption'}
        fit_data_plot = {'plot': os.path.abspath('.\\images\\fit.png'), 'caption': 'Fitted LAS and BAS adsorption peaks'
                                                                                   ' with fit.'}
        bas = {'chi_sq': round(self.analyser.bas['results'].chisqr, 4),
               'red_chi_sq': round(self.analyser.bas['results'].redchi, 4)}
        las = {'chi_sq': round(self.analyser.las['results'].chisqr, 4),
               'red_chi_sq': round(self.analyser.las['results'].redchi, 4)}
        if include_norm:
            norm_data_plot = {'plot': os.path.abspath('.\\images\\norm.png'),
                              'caption': 'Normalized FT-IR spectrum pre- and post-pyridine '
                                         'adsorption'}
            html = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='.\\utils')).\
                get_template('template.html').render(
                    date_analysed=self.analysed_dt.strftime('%d.%m.%Y'),
                    date_recorded=self.recorded_dt.strftime('%d.%m.%Y'),
                    short_git_sha=self.short_git_sha, material_id=self.material_id,
                    sample_weight=str(round(self.analyser.sample_weight, 4)) + ' g',
                    sample_radius=str(round(self.analyser.sample_disk_radius, 3)) + ' g',
                    bas_sites=str(round(self.analyser.bas['acid sites'], 2)) + ' mmol/g',
                    las_sites=str(round(self.analyser.las['acid sites'])) + ' mmol/g', raw_data_plot=raw_data_plot,
                    norm_data_plot=norm_data_plot, fitingt_data_plot=fit_data_plot, fitting_numb=4, bas=bas, las=las)
        else:
            html = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='.\\utils')).get_template(
                'template.html').render(
                date_analysed=self.analysed_dt.strftime('%d.%m.%Y'),
                date_recorded=self.recorded_dt.strftime('%d.%m.%Y'),
                short_git_sha=self.short_git_sha, material_id=self.material_id,
                ssample_weight=str(round(self.analyser.sample_weight, 4)) + ' g',
                    sample_radius=str(round(self.analyser.sample_disk_radius, 3)) + ' g',
                bas_sites=str(round(self.analyser.bas['acid sites'], 2))+' mmol/g',
                las_sites=str(round(self.analyser.las['acid sites']))+' mmol/g', raw_data_plot=raw_data_plot,
                fitting_data_plot=fit_data_plot, fitting_numb=3, bas=bas, las=las)
        self.html = html

    def save_report(self):
        """Saves the report as .pdf file"""
        dir_path = os.path.split(self.analyser.pre_adsorption_request.abs_file_path)[0]
        file_path = os.path.join(dir_path, '_'.join([self.material_id, 'report', self.short_git_sha]) + '.pdf')
        with open(file_path, 'w+b') as pdf:
            pisa.CreatePDF(src=self.html, dest=pdf)

    def _plot_data(self, plot_norm=False):
        """Helper function to plot data

            Parameters:
            -----------
            plot_norm : bool
                Indicate if normalized data should be plotted.
                """
        try:
            os.mkdir('.\\images')
        except FileExistsError:
            pass
        raw_fig, raw_ax = plt.subplots()
        raw_ax.plot(self.analyser.pre_adsorption_wavenumbers, self.analyser.pre_adsorption_data, label='Pre-adsorption')
        raw_ax.plot(self.analyser.post_adsorption_wavenumbers, self.analyser.pre_adsorption_data+1,
                    label='Pre-adsorption')
        raw_ax.invert_xaxis()
        raw_ax.set_xlabel = 'Wavenumber $cm^-1$'
        raw_ax.set_ylabel = 'Absorbance [a.u.]'
        raw_ax.legend()
        raw_fig.savefig('.\\images\\raw.png')
        if plot_norm:
            norm_fig, norm_ax = plt.subplots()
            norm_ax.plot(self.analyser.pre_adsorption_wavenumbers, self.analyser.norm_pre_adsorption_data,
                         label='Pre-adsorption')
            norm_ax.plot(self.analyser.post_adsorption_wavenumbers, self.analyser.norm_post_adsorption_data + 1,
                         label='Pre-adsorption')
            norm_ax.invert_xaxis()
            norm_ax.set_xlabel = 'Wavenumber $cm^-1$'
            norm_ax.set_ylabel = 'Absorbance [a.u.]'
            norm_ax.legend()
            norm_fig.savefig('.\\images\\norm.png')
        las, bas = list(self.analyser.py_peaks_fitted.keys())
        las_lv = self.analyser.py_peaks_fitted[las]['wavenumbers'][-1]
        bas_fv = self.analyser.py_peaks_fitted[bas]['wavenumbers'][0]
        idx_f = np.where(self.analyser.pre_adsorption_wavenumbers == las_lv)[0][0]
        idx_b = np.where(self.analyser.pre_adsorption_wavenumbers == bas_fv)[0][0]
        fitting_wavenumbers = self.analyser.pre_adsorption_wavenumbers[idx_b:idx_f]
        fitting_data = self.analyser.difference_spectra_data[idx_b:idx_f]
        zeros = np.zeros(len(fitting_data) - (len(self.analyser.py_peaks_fitted[bas]['wavenumbers']) + len(
            self.analyser.py_peaks_fitted[las]['wavenumbers'])))
        fitted_data_combined = \
            self.analyser.py_peaks_fitted[bas]['results'].best_fit.copy() + \
            self.analyser.py_peaks_fitted[bas]['baseline']
        fitted_data_combined = \
            np.append(fitted_data_combined, np.append(zeros + self.analyser.py_peaks_fitted[bas]['baseline'][-1],
                                                      self.analyser.py_peaks_fitted[las]['results'].best_fit +
                                                      self.analyser.py_peaks_fitted[las]['baseline']))
        residual = fitted_data_combined - fitting_data
        fitting_fig, (fitting_ax0, fitting_ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
        fitting_ax0.plot(fitting_wavenumbers, fitting_data, label='original data')
        fitting_ax0.plot(fitting_wavenumbers, fitted_data_combined, label='fit', linestyle='--')
        fitting_ax0.fill_between(self.analyser.py_peaks_fitted[las]['wavenumbers'],
                                 self.analyser.py_peaks_fitted[las]['baseline'],
                                 self.analyser.py_peaks_fitted[las]['results'].best_fit +
                                 self.analyser.py_peaks_fitted[las]['baseline'], color='limegreen', label='LAS sites')
        fitting_ax0.fill_between(self.analyser.py_peaks_fitted[bas]['wavenumbers'],
                                 self.analyser.py_peaks_fitted[bas]['baseline'],
                                 self.analyser.py_peaks_fitted[bas]['results'].best_fit +
                                 self.analyser.py_peaks_fitted[bas]['baseline'], color='forestgreen', label='BAS sites')
        fitting_ax1.scatter(fitting_wavenumbers, residual, label='residual')
        fitting_ax0.invert_xaxis()
        fitting_ax1.invert_xaxis()
        fitting_ax0.set_xlabel = 'Wavenumber $cm^-1$'
        fitting_ax0.set_ylabel = 'Absorbance [a.u.]'
        fitting_ax0.legend(frameon=False)
        fitting_ax1.legend(frameon=False)
        # TODO: add marker for position of peaks
        fitting_fig.savefig('.\\images\\fit.png')

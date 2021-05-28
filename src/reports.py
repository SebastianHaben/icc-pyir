import datetime
import os

import jinja2
from matplotlib import pyplot as plt
import numpy as np
from xhtml2pdf import pisa

from errors import ExperimentNoMatch, ExperimentMissMatch


class PDFReport:
    """Class to create pdf report

        Parameters:
        -----------
        analyser : utils.analyser.Analyser
            Analyser object of the experiment.
        material_id : str
            ID of the experiment.

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

        Methods:
        -------
        build_report(add_raw_data=True, add_normalization_data=True):
            Builds xml file containing the report.
        save_report(file_path):
            Saves the report as .xml file.
        get_experiment_id(experiment_id_pattern, automatic_search=True, experiment_id=None):
            Automatically finds experiment ID based on automatic search in file name ore uses provided experiment_id.
        """

    def __init__(self, analyser, material_id):
        self.analysed_dt = datetime.datetime.now()
        self.analyser = analyser
        self.material_id = material_id
        self.recorded_dt = analyser.date_measured
        self.html = None

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
            html = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='.')).\
                get_template('template.html').render(
                    date_analysed=self.analysed_dt.strftime('%d.%m.%Y'),
                    date_recorded=self.recorded_dt.strftime('%d.%m.%Y'),
                    sample_weight=str(round(self.analyser.sample_weight, 4)) + ' g',
                    sample_radius=str(round(self.analyser.sample_disk_radius, 3)) + ' g',
                    bas_sites=str(round(self.analyser.bas['acid sites'], 2)) + ' mmol/g',
                    las_sites=str(round(self.analyser.las['acid sites'])) + ' mmol/g', raw_data_plot=raw_data_plot,
                    norm_data_plot=norm_data_plot, fitingt_data_plot=fit_data_plot, fitting_numb=4, bas=bas, las=las)
        else:
            html = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='.')).get_template(
                'template.html').render(
                date_analysed=self.analysed_dt.strftime('%d.%m.%Y'),
                date_recorded=self.recorded_dt.strftime('%d.%m.%Y'),
                ssample_weight=str(round(self.analyser.sample_weight, 4)) + ' g',
                    sample_radius=str(round(self.analyser.sample_disk_radius, 3)) + ' g',
                bas_sites=str(round(self.analyser.bas['acid sites'], 2))+' mmol/g',
                las_sites=str(round(self.analyser.las['acid sites']))+' mmol/g', raw_data_plot=raw_data_plot,
                fitting_data_plot=fit_data_plot, fitting_numb=3, bas=bas, las=las)
        self.html = html

    def save_report(self):
        """Saves the report as .pdf file"""
        dir_path = os.path.join("")
        file_path = os.path.join(dir_path, '_'.join([self.material_id, 'report' + '.pdf')
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
        fitting_fig.savefig('.\\images\\fit.png')

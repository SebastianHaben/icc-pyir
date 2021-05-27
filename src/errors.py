"""Exceptions"""

# Copyright (c) 2020
# Authors: Sebastian Haben
# License: MIT


from matplotlib import pyplot as plt


class NotReadableError(Exception):
    """Exception raised for errors in the input.
        Parameters:
        -----------
        message : String
            Explanation of the error.

        Attributes:
        -----------
        message : String
            Explanation of the error.
        """

    def __init__(self, message):
        self.message = message


class NoPeaksFoundError(Exception):
    """Error raised when no peaks could be detected."""
    pass


class NormalizationError(Exception):
    """Exception raised for errors in the input.

        Attributes:
        -----------
            message : String
                      Explanation of the error.
        """

    def __init__(self, wavenumber_array, data, peaks):
        self.message = "A problem occurred during the normalization process. The calculated band position is to far " \
                       "from the initial guess."
        fig, ax = plt.subplots()
        ax.plot(wavenumber_array, data)
        low_wavenumber, high_wavenumber = ax.get_xlim()
        ax.set_xlim(high_wavenumber, low_wavenumber)
        ax.set_xlabel('Wavenumber [$cm^(-1)$]')
        ax.set_ylabel('Intensity [a.u.]')
        for peak in peaks[0]:
            ax.plot(wavenumber_array[peak], data[peak], color="red", marker="o")
        self.fig = fig
        plt.show()


class WrongFileType(Exception):
    """Raised when the input file type is not supported.

        Parameters:
        -----------
        file_ext

        Attributes:
        -----------
            message : String
                      Explanation of the error."""
    def __init__(self, file_ext):
        self.message = "The used file extension: %s is not supported" % file_ext


class ExperimentNoMatch(Exception):
    """Raised when the the two experiment IDs don't match.

        Attributes:
        -----------
            message : String
                      Explanation of the error."""
    def __init__(self):
        self.message = "No match could be found for the stated RegEx pattern!"


class AcidSiteTypeError(Exception):
    """Raised when the stated acid type doesn't match the allowed types.
        Parameters:
        -----------
            site_type : str
                Used site type.
        Attributes:
        -----------
            site_type : str
                Used site type.
            message : String
                      Explanation of the error."""
    def __init__(self, site_type):
        self.message = "The provided site_type: %s is not supported!" % site_type


class ExperimentMissMatch(Exception):
    """Raised when the the two experiment IDs don't match.

        Parameters:
        -----------
        pre_experiment_id : str
            Experiment ID of the pre pyridine adsorption spectra.
        post_experiment_id : str
            Experiment ID of the post pyridine adsorption spectra.
        Attributes:
        -----------
        message : String
            Explanation of the error."""
    def __init__(self, pre_experiment_id, post_experiment_id):
        self.message = "Two different IDs where found. Pre adsorption ID: %s Post adsorption ID: %s" % \
                       (pre_experiment_id, post_experiment_id)


class WavenumberArrayShapeMissMatch(Exception):
    """Raised when the shape of the wavenumber arrays are not of the same size.

        Attributes:
        -----------
        message : String
            Explanation of the error."""
    def __init__(self):
        self.message = "The pre- and post-pyridine adsorption spectra wavenumber arrays have a different size!"


class PositionWarning(Warning):
    """Warning if the wavenumber positions are not equal.

        Attributes:
        -----------
        message : String
            Explanation of the error."""
    def __init__(self):
        self.message = "The indicated positions of the pre- and post pyridine adsorption don't match!. A difference " \
                       "spectra was calculated and a wavenumber array created!"


class VersionError(ValueError):
    """Exception to be used when a version is not supported."""
    pass

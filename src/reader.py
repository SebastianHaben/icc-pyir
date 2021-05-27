"""Reader"""


# Copyright (c) 2020
# Authors: Sebastian Haben, Joëlle Siewe
# License: MIT
import requests
import struct

import numpy as np

from errors import NotReadableError


class Reader(object):
    """Superclass for reading files

    Parameters:
    -----------
    request : requests.Request
        Request of to be read file.

    Attributes:
    -----------
    _allowed_extensions : list of str
        List of allowed file extensions in the form ['.ext', ...].
    request : requests.Request
        Request of to be read file.

    Methods:
    --------
    read_data()
        Read the file to extract the spectral data. Has to be implemented for each subclass.

    """

    def __init__(self, request: requests.Request, allowed_extensions):
        self._allowed_extensions = allowed_extensions
        self._can_read(request)
        self.request = request

    def _can_read(self, request):
        """Checks if the provided files have the right extensions.

        Parameters
        ----------
        request: requests.Request
             Request of to be read file.
        """

        if request.extension not in self._allowed_extensions:
            raise NotReadableError(('The file: %s does not have the right file type'
                                    '. Allowed file types are: {}' % request.filename).format(
                self._allowed_extensions))

    def read_data(self):
        """Function needs to be implemented by each Reader subclass."""
        pass


class SPAReader(Reader):
    """Class to read .spa files.

    Parameters:
    -----------
    request : requests.Request
        request object of file.

    Attributes:
    -----------
    filename : str
        name of the file.

    Methods:
    --------
    read_data():
        Read the file to extract the spectral data.

    References:
    -----------
    Based on code by Zack Gainsforth.
    """
    SPEC_HEADER = 2

    def __init__(self, request):
        allowed_extensions = [".spa", ".SPA", ".srs"]
        super().__init__(request, allowed_extensions)
        self.filename = request.abs_file_path
        self.sheet = None
        self.saved_sections = None
        self.type = None

    def _sections(self):
        if self.saved_sections is None:
            with open(self.filename, 'rb') as f:
                name = struct.unpack("30s", f.read(30))[0].decode("ascii")
                extended = "Exte" in name
                f.seek(288)
                ft, v, _, n = struct.unpack('<hhhh', f.read(8))
                self.saved_sections = []
                self.type = ft
                for i in range(n):
                    # Go to the section start.
                    f.seek(304 + (22 if extended else 16) * i)
                    t, offset, length = struct.unpack('<hqi' if extended else '<hii',
                                                      f.read(14 if extended else 10))
                    self.saved_sections.append((i, t, offset, length))
        return self.saved_sections

    def _find_indextype(self, t):
        for i, a, _, _ in self._sections():
            if a == t:
                return i

    def _read_spec_header(self):
        info = self._find_indextype(self.SPEC_HEADER)
        _, _, offset, length = self._sections()[info]
        with open(self.filename, 'rb') as f:
            f.seek(offset)
            data_type, num_points, x_units, y_units, first_x, last_x, noise = \
                struct.unpack('<iiiifff', f.read(28))
            return num_points, first_x, last_x,

    def _read_spectra(self):

        self._sections()

        if self.type == 1:
            type_int = 3
        else:
            type_int = 3

        num_points, first_x, last_x = self._read_spec_header()

        _, _, offset, length = self._sections()[self._find_indextype(type_int)]

        with open(self.filename, 'rb') as f:
            f.seek(offset)
            data = np.fromfile(f, dtype='float32', count=length//4)

        if len(data) == num_points:
            domvals = np.linspace(first_x, last_x, num_points)
        else:
            domvals = np.arange(len(data))

        data = np.array([data])
        return domvals, data, None

    def read_data(self):
        """Reads data from file.

            Returns:
            --------
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample.
            """
        wavenumber_array, data, _ = self._read_spectra()
        data.shape = data.size
        return wavenumber_array, data

    def _select_sheet(self, sheet):
        """Select sheet to be read
        Parameters:
        -----------
        sheet : str
            sheet name
        """
        self.sheet = sheet


class CSVReader(Reader):
    """Superclass for reading .csv files

    Parameters:
    ----------
    request : requests.Request
              Request of to be read file.

    Methods:
    --------
    read_data()
        Read the file to extract the spectral data. Has to be implemented for each subclass.

    """

    def __init__(self, request):
        self.allowed_extensions = ['.csv', '.CSV']
        super().__init__(request, self.allowed_extensions)
        self.request = request

    def read_data(self):
        """Reads data from file.

            Returns:
            --------
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample.

            Raises:
            -------
            ValueError:
                If the given .csv file is not in the correct format.
            """
        try:
            wavenumbers, data = np.genfromtxt(self.request.abs_file_path, delimiter=',', unpack=True, skip_header=1)
        except ValueError as e:
            raise ValueError("The selected .csv file is not in the correct format be sure that is has no header and "
                             "that wavenumbers and intensities are separated by ','!") from e
        return wavenumbers, data


class ASCReader(Reader):
    """Superclass for reading .asc files

    Parameters:
    ----------
    request : requests.Request
              Request of to be read file.

    Methods:
    --------
    read_data()
        Read the file to extract the spectral data. Has to be implemented for each subclass.

    """

    def __init__(self, request):
        self.allowed_extensions = ['.asc', '.ASC']
        super().__init__(request, self.allowed_extensions)
        self.request = request

    def read_data(self):
        """Reads data from file.

            Returns:
            --------
            wavenumber_array : np.array
                Array of recorded wavenumbers
            data : np.array
                Data array. Usually Absorbance of sample.

            Raises:
            -------
            ValueError:
                If the given .asc file is not in the correct format.
            """
        try:
            wavenumbers, data = np.genfromtxt(self.request.abs_file_path, unpack=True, skip_header=82)
        except ValueError as e:
            raise ValueError("The selected .asc file is not in the correct format be sure that is has no header and "
                             "that wavenumbers and intensities are separated by ' '!") from e
        return wavenumbers, data

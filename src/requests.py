"""Requests"""

# Copyright (c) 2020
# Authors: Sebastian Haben
# License: MIT

import glob
import os


class Requests(object):
    """Handles request for multiple files.

    Parameters:
    -----------
    uris : path or list
        Either path to directory or list of paths or files.
    extension : str
        File extension to search for in directory.

    Attributes:
    -----------
    requests_list : list
        List of requests.Request objects.
    """

    def __init__(self, uris, extension=None):
        self.requests_list = []
        self._extension = extension
        self._parse_uris(uris)

    def _parse_uris(self, uris):
        """Parses the given uris to either find all applicable files or checks if all the given paths/files are
        applicable.

        Parameters:
        -----------
        uris : path or list
            Either path to directory or list of paths or files.

        Raises:
        -------
        FileNotFoundError
            If no files with the correct extension could be found.
        TypeError:
            If no file extension was provided.
        """
        # List provided
        if isinstance(uris, list):
            for item in uris:
                # If list of directories was given
                if os.path.isdir(item):
                    try:
                        search_path = os.path.join(item, '*' + self._extension)
                        file_list = glob.glob(search_path)
                        # Handles case where no file was found
                        if len(file_list) == 0:
                            raise FileNotFoundError('No files could be found. Make sure to give the correct file '
                                                    'extension')
                        else:
                            for file in file_list:
                                self.requests_list.append(Request(file))
                    except TypeError:
                        raise TypeError('Please give a file extension when working with directory paths.')
                # If list of files was given
                elif os.path.isfile(item):
                    if os.path.splitext(item)[1] in self._extension:
                        self.requests_list.append(Request(item))
                    else:
                        raise FileNotFoundError('No files could be found. Make sure to give the correct file '
                                                'extension')

        # Single file given
        elif hasattr(uris, 'read') and hasattr(uris, 'close'):
            self.requests_list.append(Request(uris))

        # Single string given
        elif isinstance(uris, str):
            # Single Directory path given
            if os.path.isdir(uris):
                try:
                    search_path = os.path.join(uris, '*' + self._extension)
                    file_list = glob.glob(search_path)
                    if len(file_list) == 0:
                        try:
                            search_path = os.path.join(uris, '*', '*' + self._extension)
                            file_list = glob.glob(search_path, recursive=True)
                            # Handles case where no file was found
                            if len(file_list) == 0:
                                raise FileNotFoundError('No files could be found. Make sure to give the correct file '
                                                        'extension')
                            else:
                                for file in file_list:
                                    self.requests_list.append(Request(file))
                        except TypeError:
                            raise TypeError('Please give a file extension when working with directory paths.')
                    # If directory contains files with correct extension
                    else:
                        for file in file_list:
                            self.requests_list.append(Request(file))
                except TypeError:
                    raise TypeError('Please give a file extension when working with directory paths.')


class Request(object):
    """Handles the fetching of a single file.

    Parameters:
    -----------
    uri : path or file
        Path or file object of file.

    Attributes:
    -----------
    file_name : str
    file_extension : str
    absolute_file_path : str

    Methods:
    --------
    get_file()
        Returns file object."""

    def __init__(self, uri):
        file_name, extension, abs_file_path = self._parse_uri(uri)
        self.file_name = file_name
        self.extension = extension
        self.abs_file_path = abs_file_path

    def _parse_uri(self, uri):
        """Extracts filename , extension and makes path absolute.

        Parameters:
        -----------
        uri : file or path
            Path or file object of file.

        Returns:
        --------
        file_name : str
        file_extension : str
        absolute_file_path : str

        Raises:
        -------
        FileNotFoundError:
            If provided uri is neither a file object nor a path to a file or if the specified path does not exist.
        """

        # Check if file
        if hasattr(uri, 'read') and hasattr(uri, 'close'):
            file_path = uri.name
            file_name = os.path.basename(file_path)
            _, extension = os.path.splitext(file_name)
            abs_file_path = os.path.abspath(file_path)

        # Check if path
        elif isinstance(uri, str):
            if os.path.exists(uri):
                file_path = uri
                file_name = os.path.basename(file_path)
                _, extension = os.path.splitext(file_name)
                abs_file_path = os.path.abspath(file_path)
            else:
                raise FileNotFoundError('The requested file could not be found.')

        else:
            raise FileNotFoundError('The requested file could not be found.')

        return file_name, extension, abs_file_path

    def get_file(self):
        """Get a file for the resource

        Returns:
        --------
        file : file
            Opened file object in 'r' mode."""

        return open(self.abs_file_path, 'r')

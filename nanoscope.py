# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, unicode_literals

import io

import numpy as np
import six

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'figure.figsize': (7.2,4.45)})
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)




mpl.rcParams.update({'xtick.labelsize': 14})

mpl.rcParams.update({'ytick.labelsize': 14})

mpl.rcParams.update({'font.size': 12})

mpl.rcParams.update({'figure.autolayout': True})

mpl.rcParams.update({'axes.titlesize': 14})

mpl.rcParams.update({'axes.labelsize': 14})

mpl.rcParams.update({'lines.linewidth': 2})

mpl.rcParams.update({'lines.markersize': 6})

# mpl.rcParams.update({'legend.fontsize': 13})

mpl.rcParams.update({'legend.fontsize': 10})



def read(f, encoding='cp1252', header_only=False, check_version=True):
    """
    Reads the specified file, given as either a filename or an already opened
    file object. Passed file objects must be opened in binary mode. Meant as the
    typical entry point for loading in afm data.

    :param f: Filename of the file to read or an opened file object. File
              objects must be opened in binary mode.
    :param encoding: The encoding to use when reading the file header. Defaults
                     to cp1252.
    :param header_only: Whether to read only the header of the file. Defaults to
                        False.
    :param check_version: Whether to enforce version checking for known
                          supported versions. Defaults to True.
    :returns: A NanoscopeFile object containing the image data.
    :raises OSError: If a passed file object is not opened in binary mode.
    """
    try:
        with io.open(f, 'rb') as file_obj:
            images = NanoscopeFile(file_obj, encoding, header_only, check_version)
    except TypeError:
        if 'b' not in f.mode:
            raise OSError('File must be opened in binary mode.')
        images = NanoscopeFile(f, encoding, header_only, check_version)
    return images


class NanoscopeFile(object):
    """
    Handles reading and parsing Nanoscope files.
    """
    supported_versions = ['0x05120000', '0x05120130', '0x09300201', '0x05310001']

    def __init__(self, file_object, encoding='utf-8', header_only=False, check_version=True):
        self.images = {}
        self.config = {'_Images': {}}
        self.encoding = encoding

        self._read_header(file_object, check_version)
        if not header_only:
            for image_type in six.iterkeys(self.config['_Images']):
                self._read_image_data(file_object, image_type)

    @property
    def height(self):
        """
        Return the height image if it exists, else ``None``.
        """
        return self.image('Height')

    @property
    def amplitude(self):
        """
        Return the amplitude image if it exists, else ``None``.
        """
        return self.image('Amplitude')

    @property
    def phase(self):
        """
        Return the phase image if it exists, else ``None``.
        """
        return self.image('Phase')

    def image(self, image_type):
        """
        Returns the specified image type if it exists, else ``None``.
        """
        return self.images.get(image_type, None)

    def image_types(self):
        """
        Returns a list of names for all image types.
        """
        return list(self.images.keys())

    def describe_images(self):
        """
        Returns a list of tuples (key, info) describing the image types.
        """
        return [(k, self.image(k).description) for k in self.image_types()]

    def __iter__(self):
        for v in six.itervalues(self.images):
            yield v

    def _read_header(self, file_object, check_version=True):
        """
        Read the Nanoscope file header.

        :param file_object: Opened file
        :param check_version: Whether to enforce version checking for known
                              supported versions. Defaults to True.
        :raises UnsupportedVersion: If the version is not supported and version
                                    checking is enabled.
        """
        file_object.seek(0)
        for line in file_object:
            parameter = parse_parameter(line, self.encoding)
            if not self._validate_version(parameter) and check_version:
                raise UnsupportedVersion(parameter.hard_value)
            if self._handle_parameter(parameter, file_object):
                return

    def _read_image_data(self, file_object, image_type):
        """
        Read the raw data for the specified image type if it is in the file.

        :param image_type: String indicating which image type to read.
        :returns: A NanoscopeImage instance of the specified type
        :raises MissingImageData: If the image_type indicated is not in the file
        """
        if image_type not in self.config['_Images']:
            raise MissingImageData(image_type)

        config = self.config['_Images'][image_type]
        data_offset = config['Data offset']
        data_size = config['Bytes/pixel']
        number_lines = config['Number of lines']
        samples_per_line = config['Samps/line']

        file_object.seek(data_offset)
        number_points = number_lines * samples_per_line
        raw_data = (np.frombuffer(file_object.read(data_size * number_points),
                                  dtype='<i{}'.format(data_size),
                                  count=number_points)
                   .reshape((number_lines, samples_per_line)))

        scan_size = self._get_config_fuzzy_key(config, ['Scan size', 'Scan Size'])

        self.images[image_type] = NanoscopeImage(
            image_type,
            raw_data,
            config['Bytes/pixel'],
            config['Z magnify'],
            self._get_sensitivity_value(image_type, 'Z scale'),
            self._get_sensitivity_value(image_type, 'Z offset'),
            scan_size * scan_size,
            config['Description'],
        )
        return self.images[image_type]

    def _get_sensitivity_value(self, image_type, key):
        parameter = self.config['_Images'][image_type][key]
        sensitivity = self.config[parameter.soft_scale]
        value = parameter.hard_value
        return sensitivity * value

    def _get_config_fuzzy_key(self, config, keys):
        for k in keys:
            value = config.get(k, None)
            if value is not None:
                return value
        raise KeyError

    def _validate_version(self, parameter):
        if parameter.type == 'H' or parameter.parameter != 'Version':
            return True
        return parameter.hard_value in self.supported_versions

    def _handle_parameter(self, parameter, f):
        if parameter.type == 'H':  # header
            if parameter.header == 'File list end':
                return True
            if parameter.header == 'Ciao image list':
                return self._handle_parameter(self._read_image_header(f), f)
        elif parameter.type == 'V':
            if not parameter.soft_scale and not parameter.hard_scale:
                self.config[parameter.parameter] = parameter.hard_value
            else:
                self.config[parameter.parameter] = parameter
        elif parameter.type != 'S':
            self.config[parameter.parameter] = parameter.hard_value
        return False

    def _read_image_header(self, f):
        image_config = {}
        for line in f:
            parameter = parse_parameter(line, self.encoding)
            if parameter.type == 'H':
                return parameter
            elif parameter.type == 'S':
                if parameter.parameter == 'Image Data':
                    image_config['Image Data'] = parameter.internal
                    image_config['Description'] = parameter.external
                    self.config['_Images'][parameter.internal] = image_config
            elif parameter.type == 'V':
                if not parameter.soft_scale and not parameter.hard_scale:
                    image_config[parameter.parameter] = parameter.hard_value
                else:
                    image_config[parameter.parameter] = parameter
            else:
                image_config[parameter.parameter] = parameter.hard_value


# -*- coding: utf-8 -*-



class Error(Exception):
    """Base exception for nanoscope errors."""
    pass


class UnsupportedVersion(Error):
    """Error for unsupported SPM file version."""

    def __init__(self, version):
        self.version = version

    def __str__(self):
        return 'Unsupported file version {}'.format(self.version)


class MissingImageData(Error):
    """Error for missing data for defined image in header."""

    def __init__(self, image):
        self.image = image

    def __str__(self):
        return ('Image type {} found in header '
                'but is missing data'.format(self.image))


class InvalidParameter(Error):
    """Error for incorrectly formatted Ciao parameter."""

    def __init__(self, parameter):
        self.parameter = parameter

    def __str__(self):
        return '"{}" is not a valid Ciao parameter'.format(self.parameter)


# -*- coding: utf-8 -*-
# from __future__ import absolute_import, division, unicode_literals

import numpy as np


class NanoscopeImage(object):
    """
    Holds the data associated with a Nanoscope image.
    """
    supported_colortables = {
        12: {
            'r': (lambda p: np.clip(np.round(
                  p * (10200 / 37) - (765 / 37)), 0, 255)),
            'g': (lambda p: np.clip(np.round(
                  p * (30600 / 73) - (11985 / 73)), 0, 255)),
            'b': (lambda p: np.clip(np.round(
                  p * (6800 / 9) - (4505 / 9)), 0, 255)),
        },
    }

    def __init__(self, image_type, raw_data, bytes_per_pixel, magnify,
                 scale, offset, scan_area, description):
        self.unit = scale.unit.to_string()
        self.bytes_per_pixel = bytes_per_pixel
        self.magnify = magnify
        self.raw_data = raw_data
        self.flat_data = None
        self.converted_data = None
        self.type = image_type
        self.scale = scale.value
        self.offset = offset.value
        self.height_scale = self.scale * magnify
        self.scan_area = scan_area
        self.description = description

        self._cache = {}

    @property
    def data(self):
        """
        Returns the most processed form of the data.
        """
        if self.converted_data is None:
            if self.flat_data is None:
                return self.raw_data
            return self.flat_data
        return self.converted_data

    def process(self, order=1):
        """
        Flattens and converts the raw data. Convenience function that reduces
        the manual steps needed.

        :param order: The order of the polynomial to use when flattening.
                      Defaults to 1 (linear), which should give good results
                      for most images.
        :returns: The image with flattened and converted data for chaining
                  commands.
        """
        return self.flatten(order).convert()

    def flatten(self, order=1):
        """
        Flattens the raw data, by fitting each scanline to a polynomial with
        the order specified and subtracting that fit from the raw data.

        Typically happens prior to converting from raw data.

        :param order: The order of the polynomial to use when flattening.
                      Defaults to 1 (linear).
        :returns: The image with flattened data for chaining commands.
        """
        self.flat_data = np.round([self._flatten_scanline(line, order)
                                   for line in self.raw_data])
        self._cache.clear()
        return self

    def convert(self):
        """
        Converts the raw data into data with the proper units for that image
        type (i.e. nm for Height, V for Amplitude).

        Typically happens after flattening the data.

        :returns: The image with converted data for chaining commands.
        """
        if self.flat_data is None:
            self.flat_data = self.raw_data
        value = self.scale / pow(2, 8 * self.bytes_per_pixel)
        self.converted_data = self.flat_data * value
        self._cache.clear()
        return self

    def colorize(self, colortable=12):
        """
        Colorizes the data according to the specified height scale. Currently
        uses colorscale #12 from Nanoscope as hardcoded behavior.

        :param colortable: The Nanoscope colortable to use.
                           Only 12 is supported, and is the default.
        :returns: The pixels of the image ready for use with
                  ``Pillow.Image.fromarray``.
        :raises ValueError: If the colortable is not supported.
        """
        if colortable not in self.supported_colortables:
            raise ValueError('Colortable {} is not '
                             'currently supported'.format(colortable))

        colors = self.supported_colortables[colortable]
        get_color = (lambda v:
            np.array([colors[c]((v + (self.height_scale / 2)) /
                                self.height_scale) for c in 'rgb']))
        if self.converted_data is None:
            self.converted_data = self.data

        data = []
        for row in reversed(self.converted_data):
            data.append([])
            for col in row:
                data[-1].append(get_color(col))
        return np.array(data, dtype=np.uint8)

    def reset_height_scale(self):
        """
        Resets the height scale to the original value from the file.
        """
        self.height_scale = self.magnify * self.scale

    @property
    def mean_height(self):
        """
        Returns the mean height of the data in nm. For a typical processed
        scan, this value should be about zero.

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'mean_height' not in self._cache:
            self._cache['mean_height'] = np.mean(self.data)
        return self._cache['mean_height']

    @property
    def mean_roughness(self):
        """
        Returns the mean roughness of the data in nm.

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'mean_roughness' not in self._cache:
            self._cache['mean_roughness'] = np.mean(np.abs(self.data - self.mean_height))
        return self._cache['mean_roughness']

    @property
    def rms_roughness(self):
        """
        Returns the root mean square roughness of the data in nm.

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'rms_roughness' not in self._cache:
            self._cache['rms_roughness'] = (
                np.sqrt(np.sum(np.square(
                    self.data - self.mean_height)) / self.data.size))
        return self._cache['rms_roughness']

    @property
    def total_roughness(self):
        """
        Returns the total roughness of the data in nm. This is defined as the
        difference between the highest peak and the lowest valley.
        """
        return self.max_valley + self.max_peak

    @property
    def max_valley(self):
        """
        Returns the depth of the lowest valley in nm. A valley is defined
        relative to the mean height (typically 0nm).
        """
        return abs(self.min_height - self.mean_height)

    @property
    def max_peak(self):
        """
        Returns the height of the highest valley in nm. A peak is defined
        relative to the mean height (typically 0nm).
        """
        return self.max_height - self.mean_height

    @property
    def mean_valley(self):
        """
        Returns the depth of the average valley in nm. A valley is defined
        relative to the mean height (typically 0nm).

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'mean_valley' not in self._cache:
            valley_elems = self.data[self.data < self.mean_height]
            self._cache['mean_valley'] = (
                np.sum(np.abs(
                    valley_elems - self.mean_height)) / valley_elems.size)
        return self._cache['mean_valley']

    @property
    def mean_peak(self):
        """
        Returns the height of the average peak in nm. A peak is defined
        relative to the mean height (typically 0nm).

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'mean_peak' not in self._cache:
            peak_elems = self.data[self.data > self.mean_height]
            self._cache['mean_peak'] = (
                np.sum(peak_elems - self.mean_height) / peak_elems.size)
        return self._cache['mean_peak']

    @property
    def mean_total_roughness(self):
        """
        Returns the mean total roughness in nm. This is defined as the
        difference between the mean peak and valley.
        """
        return self.mean_peak + self.mean_valley

    @property
    def min_height(self):
        """
        Returns the minimum height in the image in nm.

        The value is calculated on first access and cached for later. Running
        convert or flatten will force a recalculation on the next access.
        """
        if 'min_height' not in self._cache:
            self._cache['min_height'] = np.min(self.data)
        return self._cache['min_height']

    @property
    def max_height(self):
        """
        Returns the maximum height in the image in nm.

        The value is calculated on first access and cached for later. Running
        conver or flatten will force a recalculation on the next access.
        """
        if 'max_height' not in self._cache:
            self._cache['max_height'] = np.max(self.data)
        return self._cache['max_height']

    def n_point_roughness(self, n=5):
        """
        Returns the average roughness in nm, defined as the mean of the n
        highest peaks and n lowest valleys.

        :param n: The number of points to take from both peaks and valleys.
        :returns: The average roughness of the n highest peaks and n lowest
                  valleys, in nm.
        """
        peak_elems = np.sort(self.data[self.data > self.mean_height])[-n:]
        valley_elems = np.sort(self.data[self.data < self.mean_height])[:n]
        return np.mean(peak_elems + valley_elems)

    def peak_count(self, threshold=None):
        """
        Calculates the total number of peaks and valleys in the image. A peak
        or valley is defined as any feature that exceeds the provided threshold.

        :param threshold: The threshold to use for defining a peak or valley.
                          Defaults to the mean roughness (Ra).
        :returns: The total number of peaks and valleys in the image.
        """
        threshold = threshold or self.mean_roughness
        return self.data[np.abs(self.data) >= abs(threshold)].size

    def peak_density(self, threshold=None):
        """
        Calclulates the number of peaks or valleys per unit area in the image,
        with units of peaks per square μm. A peak or valley is defined as any
        feature that exceeds the provided threshold.

        :param threshold: The threshold to use for defining a peak or valley.
                          Defaults to the mean roughness (Ra).
        :returns: The number of peaks and valleys in the image per square μm.
        """
        return self.peak_count(threshold) / self.scan_area

    def high_spot_count(self, threshold=None):
        threshold = threshold or self.mean_roughness
        return self.data[self.data >= threshold].size

    def low_spot_count(self, threshold=None):
        threshold = threshold or self.mean_roughness
        return self.data[self.data <= threshold].size

    def _flatten_scanline(self, data, order=1):
        coefficients = np.polyfit(range(len(data)), data, order)
        correction = np.array(
            [sum([pow(i, n) * c
            for n, c in enumerate(reversed(coefficients))])
            for i in range(len(data))])
        return data - correction

    Ra = mean_roughness
    Rq = rms_roughness
    rms = rms_roughness
    Rp = max_peak
    Rv = max_valley
    Rt = total_roughness
    zrange = total_roughness
    Rpm = mean_peak
    Rvm = mean_valley
    Rz = property(lambda self: self.n_point_roughness(n=5))
    Pc = property(lambda self: self.peak_count(self.mean_roughness))
    Pd = property(lambda self: self.peak_density(self.mean_roughness))
    HSC = property(lambda self: self.high_spot_count(self.mean_roughness))
    LSC = property(lambda self: self.low_spot_count(self.mean_roughness))


# -*- coding: utf-8 -*-
# from __future__ import absolute_import, division, unicode_literals

import datetime
import re

import six
from astropy import units as u



__all__ = ['parse_parameter']


class CiaoParameter(object):
    """
    Parent class for generic values from the header.
    """
    type = 'P'

    def __init__(self, parameter, hard_value):
        self.parameter = parameter
        self.hard_value = self._parse_value(hard_value)

    def __str__(self):
        return '{0}: {1}'.format(self.parameter, self.hard_value)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.parameter == other.parameter and
                self.hard_value == other.hard_value)

    def __ne__(self, other):
        return not self.__eq__(other)

    def _parse_value(self, value):
        """
        Try to parse and return the value as the following values:

        * ``datetime.datetime``
        * ``int``
        * ``float``
        * ``str``
        * ``None``

        Trailing whitespace is stripped and an empty string is returned
        as ``None``.
        """
        try:
            return datetime.datetime.strptime(value, '%I:%M:%S %p %a %b %d %Y')
        except:
            try:
                split_value = value.strip().split(' ')[0]
                if split_value in ('', 'None'):
                    return None
            except AttributeError:
                return value
            try:
                return int(split_value)
            except ValueError:
                try:
                    return float(split_value)
                except ValueError:
                    return value


class CiaoValue(CiaoParameter):
    """
    Represents a Ciao value object.
    """
    type = 'V'

    def __init__(self, parameter, soft_scale, hard_scale, hard_value):
        self.parameter = parameter
        self.soft_scale = super(CiaoValue, self)._parse_value(soft_scale)
        self.hard_scale = self._parse_value(hard_scale)
        self.hard_value = self._parse_value(hard_value)

    def __str__(self):
        return '{0}: [{1}] ({2}) {3}'.format(self.parameter, self.soft_scale,
                                             self.hard_scale, self.hard_value)

    def __eq__(self, other):
        return (self.parameter == other.parameter and
                self.soft_scale == other.soft_scale and
                self.hard_scale == other.hard_scale and
                self.hard_value == other.hard_value)

    def _parse_value(self, value):
        if value is None:
            return None

        try:
            return datetime.datetime.strptime(value, '%I:%M:%S %p %a %b %d %Y')
        except:
            if value.strip() in ('', 'None'):
                return None
            try:
                return u.Quantity(value.strip().replace('º', 'deg'))
            except (ValueError, TypeError):
                return value


class CiaoScale(CiaoParameter):
    """
    Represents a Ciao scale object.
    """
    type = 'C'

    def __init__(self, parameter, soft_scale, hard_value):
        self.parameter = parameter
        self.soft_scale = self._parse_value(soft_scale)
        self.hard_value = self._parse_value(hard_value)

    def __str__(self):
        return '{0}: [{1}] {2}'.format(self.parameter, self.soft_scale,
                                       self.hard_value)

    def __eq__(self, other):
        return (self.parameter == other.parameter and
                self.soft_scale == other.soft_scale and
                self.hard_value == other.hard_value)


class CiaoSelect(CiaoParameter):
    """
    Represents a Ciao select object.
    """
    type = 'S'

    def __init__(self, parameter, internal, external):
        self.parameter = parameter
        self.internal = internal
        self.external = external

    def __str__(self):
        return '{0}: [{1}] "{2}"'.format(self.parameter,
                                         self.internal,
                                         self.external)

    def __eq__(self, other):
        return (self.parameter == other.parameter and
                self.internal == other.internal and
                self.external == other.external)


class CiaoSectionHeader(CiaoParameter):
    """
    Represents a Ciao section header.
    """
    type = 'H'

    def __init__(self, header):
        self.header = header

    def __str__(self):
        return self.header

    def __eq__(self, other):
        return self.header == other.header


def decode(string, encoding='utf-8'):
    """
    Decodes the binary string with the specified encoding (or passes through
    a non-binary string) and strips newlines off the end.

    :param string: The string, may be binary or non-binary but should be text
                   data.
    :param encoding: The encoding to use for a binary string. Defaults to utf-8.
    :returns: The decoded and stripped string.
    :raises TypeError: If the passed parameter is not a valid string.
    :raises UnicodeError: When decoding a binary string if there are encoding
                          errors.
    """
    if isinstance(string, six.text_type):
        return string.rstrip('\r\n')

    try:
        string = string.decode(encoding=encoding).rstrip('\r\n')
    except AttributeError:
        raise TypeError('Invalid type {} passed.'.format(type(string)))

    return string


def parse_parameter(string, encoding='utf-8'):
    """
    Factory function that parses the parameter string and creates the
    appropriate CiaoParameter object.

    :param string: The parameter string to parse.
    :returns: The CiaoParameter that corresponds with the parameter string.
    :raises ValueError: If the string is not a valid CiaoParameter.
    """
    string = decode(string, encoding)
    header_match = re.match(r'\\\*(?P<header>.+)', string)
    if header_match is not None:
        return CiaoSectionHeader(header_match.group('header'))

    regex = re.compile(r'\\(?P<ciao>@?)(?:(?P<group>[0-9]+):)?'
                       r'(?P<parameter>[^:]+): '
                       r'(?:(?P<type>[VCS]) )?(?P<value>.*)')
    m = regex.match(string)
    if m is None:
        raise InvalidParameter(string)

    parameter_type = m.group('type')
    value = m.group('value')
    if parameter_type == 'V':  # value
        value_regex = (r'(?:\[(?P<soft_scale>[^\[\]]*)\] )?'
                       r'(?:\((?P<hard_scale>[^\(\)]*)\) )?'
                       r'(?P<hard_value>[^\[\]\(\)]*)')
        vm = re.match(value_regex, value)
        return CiaoValue(m.group('parameter'), vm.group('soft_scale'),
                         vm.group('hard_scale'), vm.group('hard_value'))
    elif parameter_type == 'C':  # scale
        scale_regex = (r'(?:\[(?P<soft_scale>[^\[\]]*)\] )?'
                       r'(?P<hard_value>[^\[\]]*)')
        cm = re.match(scale_regex, value)
        return CiaoScale(m.group('parameter'), cm.group('soft_scale'),
                         cm.group('hard_value'))
    elif parameter_type == 'S':  # select
        select_regex = (r'(?:\[(?P<internal_designation>[^\[\]]*)\] )?'
                        r'(?:\"(?P<external_designation>[^\"]*)\")?')
        sm = re.match(select_regex, value)
        return CiaoSelect(m.group('parameter'),
                          sm.group('internal_designation'),
                          sm.group('external_designation'))
    else:  # simple value
        return CiaoParameter(m.group('parameter'), value)
################################################xxxxxx######################################################
################################################xxxxxx######################################################
####################################### AFM image pipelines below ##########################################
################################################xxxxxx######################################################
################################################xxxxxx######################################################

from scipy.io import FortranFile
from scipy import interpolate as sci_interpolate
class single_lattice_image:
    '''reads the simulated lattice images
 
    '''
    import os
    import shutil,glob

    def __init__(self,filepath, filename,qz_prime):
        '''here, it prepares to read the lattice image. Reading happens, and scattering image is shown when .rdframe() is called
        filename is letters after 000 format. 
        '''
        self.filepath = filepath
        self.filename = filename
        self.qz_prime = qz_prime
        
        p = read(filepath+filename)
        p.height.process()
        image = p.height
        data = image.data
        
        #below is just to test if the fourier transform has other multiplicative factors .. keep it disabled
#         data = np.concatenate((data,data), axis = 1)
#         data = np.concatenate((data,data), axis = 0)
        
#         data = np.concatenate((data,data), axis = 1)
#         data = np.concatenate((data,data), axis = 0)
        
#         p.height.process()
# #         print(p.height.zrange, p.height.rms)
#         pixels = p.height.colorize()
#         data = Image.fromarray(pixels)
#         data = np.asarray(data)[:,:,0]
        data = np.rot90(data) #so that the image is rotated 90 degree counter clock-wise
        self.image = data
        self.lattice_size = np.shape(self.image)[0]
        self.range = p.height.Rp
        

    
    def rdframe(self,n):
        
        #spline fit
        img_shape = [self.lattice_size,self.lattice_size]
        x = np.arange(0, img_shape[1], 1)
        y = np.arange(0, img_shape[0], 1)
        xx, yy = np.meshgrid(x, y)
        image = self.image
        
        ##### enable the line below (divison should be by (latticeUnit_to_nm*img_shape[1])**2) or
        fixed_img = image
        
        ###### enable spline fit below ///// in this case division should be just (img_shape[1])**2
#         img_raw = image[yy,xx] #np.sin(xx**2+yy**2)
# #         print(np.shape(image),np.shape(img_raw))
#         spline_img = sci_interpolate.interp2d(x, y, img_raw, kind='cubic')
        
        
#         xnew = np.arange(0, latticeUnit_to_nm*img_shape[1], 1)
#         ynew = np.arange(0, latticeUnit_to_nm*img_shape[0], 1)
# #                 print(centery)
#         fixed_img = spline_img(xnew, ynew)
# #         print(xnew)
       

        #Fourier transform and then do the absolute square
#         fftimage= np.fft.fftn(np.exp(-1j *self.qz_prime*fixed_img)) ####
#         gixaxs = fftimage * np.conj(fftimage)
#         gixaxs = gixaxs/(self.qz_prime**2)/(img_shape[1])**2 #### divided by qzprime squre and then by the area
# #         print((latticeUnit_to_nm*img_shape[1])**2)
#         scattering = np.real(gixaxs[:lattice_size//2,:])
#         return(scattering)
        fftimage= np.fft.fft2(np.exp(-1j *self.qz_prime*fixed_img),axes=(0, 1)) #   image_mod  np.exp(-1j *self.qz_prime*image_mod)
        gixaxs = fftimage * np.conj(fftimage)
        gixaxs = gixaxs/((self.qz_prime)**2)/(self.lattice_size)**2 #normalized by area
        
        gixaxs= np.fft.fftshift(gixaxs)#,axes=1)
#         gixaxs= np.fft.fftshift(gixaxs,axes=0)

        scattering = np.real(gixaxs[:,:]) # may do lattice_size//2 limit if you desire. will have lattice_size x lattice_size 
        half_pt = self.lattice_size//2
        scattering[half_pt,half_pt] = np.average(scattering[half_pt-2:half_pt+3:3,half_pt-2:half_pt+3:3]) # skip the center point! it will affect average
        return(scattering)
    
    def rdlattice(self,n):
        
        return(self.image)
    
    def height_range(self):
        
        return(self.range)
# #rewrite the pixel to q conversion
# def convert_pixel_to_q(x,y):
#     x,y=np.int(x),np.int(y)
#     """for Fourier transform images, we know the conversion already."""
#     q_range = np.linspace(0, (lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm)), num=lattice_size//2 +1)
    
#     q_x,q_y,q_z,q,alpha_f,qx_prime,qy_prime,qz_prime = 0,q_range[x],q_range[y],(q_range[x]**2+q_range[x]**2)**(0.5),0,0,q_range[x],q_range[y]
    
#     return(q_x,q_y,q_z,q,alpha_f,qx_prime,qy_prime,qz_prime) 

# #plot the lattaice image
# '''
# The following functions work after loading the class

# '''
# from datetime import date, datetime
# print(datetime.now())
# def plot_lattice_image(frames=[], show_subtitle = True, save = 0):
#     image = FD.rdlattice(frames[0])
#     for i in frames[1:]:
#         image += FD.rdlattice(i)
#     plt.figure(figsize=(7,7))
#     dimension = np.shape(image)[0]
#     rms = np.sqrt(np.mean(image**2))
    
#     ax = plt.gca()
#     im = ax.imshow(image, vmin = -FD.height_range() , vmax = FD.height_range(), origin = 'lower', extent = [0,latticeUnit_to_nm*dimension,0,latticeUnit_to_nm*dimension], interpolation = 'none')#cmap=cmap_albula) 
#     ax.set_ylabel(r'$L$ [nm]')
#     ax.set_xlabel(r'$L$ [nm]')
#     # create an axes on the right side of ax. The width of cax will be 5%
#     # of ax and the padding between cax and ax will be fixed at 0.05 inch.
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     plt.colorbar(im, cax=cax)
    
#     if show_subtitle == True:
#         ax.set_title('Lattice image')
    
#     if save == 1:
#         # current date and time
#         now = datetime.now()
#         # timestamp = datetime.timestamp(now)
#         now = np.str_(now)
#         now = now.replace('-','')
#         now = now.replace(':','')
#         now = now.replace(' ','')
#         index = now.index('.')
#         now = now[:index]

#         try:
#             uid
#         except NameError:
#             print("well, uid WASN'T defined! saving in primary folder")
#             fp = '/home/pmyint/' +'AFM'  + '_'+ now + '.eps'
#         else:
#             print("uid was defined. Saving in respective folder")
#             directory_exists = os.path.isdir('/home/pmyint/' + '%s_'%(uid) +'/')
#             if directory_exists == False:
#                 os.makedirs('/home/pmyint/' + '%s_'%(uid) +'/')
#             fp = '/home/pmyint/' + '%s_'%(uid) +'/' +'AFM' +  now + '.eps'

#         plt.savefig(fp,bbox_inches='tight')
#     plt.show()
#     print('RMS of lattice image is ',np.sqrt(np.mean(image**2)))
#     print('Variance of lattice image is ',np.var(image))

# #plot scattering pattern
# def plot_scattering_image(frames=[], logscale = False, pixels= False): ####
#     scattering = FD.rdframe(frames[0])
#     for i in frames[1:]:
#         scattering += FD.rdframe(i)
#     plt.figure(figsize=(7,7))
#     dimension = np.shape(scattering)[0]
# #     if logscale:
# #         scattering = np.log(scattering)
#     rms = np.sqrt(np.mean(scattering**2))
    
#     scattering = scattering[0:dimension,0:dimension] ####
    
#     ax = plt.gca()
#     if logscale:
#         if pixels:
#             im = ax.imshow(scattering, origin = 'lower', interpolation = 'none',norm=LogNorm(vmin=1, vmax=np.amax(scattering)//2))#cmap=cmap_albula) 
#         else:
#             im = ax.imshow(scattering, origin = 'lower', extent = [0,(lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),0,(lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm))], interpolation = 'none',norm=LogNorm(vmin=1, vmax=np.amax(scattering)))#cmap=cmap_albula) 
#     else:
#         if pixels:
#             im = ax.imshow(scattering, vmin = 0 , vmax = rms//100, origin = 'lower')#cmap=cmap_albula) 
#         else:
#             im = ax.imshow(scattering, vmin = 0 , vmax = rms, origin = 'lower', extent = [0,(lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),0,(lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm))])#cmap=cmap_albula) 

#     ax.set_ylabel(r'$q_y$ $[nm^{-1}]$')
#     ax.set_xlabel(r'$q_x$ $[nm^{-1}]$')
#     # create an axes on the right side of ax. The width of cax will be 5%
#     # of ax and the padding between cax and ax will be fixed at 0.05 inch.
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     plt.colorbar(im, cax=cax)
#     ax.set_title('Scattering image')
#     plt.show()
#     print('RMS of scattering image is ',rms)
#     print('Max of scattering image is ',np.amax(scattering))
    
# def slice_scat_img(x,y,frames=[0]):
#     scattering = FD.rdframe(frames[0])
#     for i in frames[1:]:
#         scattering += FD.rdframe(i)
#     plt.figure(figsize=(7,7))
#     dimension = np.shape(scattering)[0]
# #     if logscale:
# #         scattering = np.log(scattering)
#     rms = np.sqrt(np.mean(scattering**2))
    
#     scattering = scattering[0:dimension,0:dimension] ####
    
#     if isinstance(y,list):
#         slice_avg = np.average(scattering[y[0]:y[1],x[0]:x[1]], axis = 0)
#         y_avg = y[0]
#     else:
#         slice_avg = scattering[y,x[0]:x[1]]
#         y_avg = y
#     x_pixels = np.arange(x[0],x[1])
#     qvals = []
    
#     for i in range(len(x_pixels)):
#         qvals.append(convert_pixel_to_q(i,y_avg)[1])  
        
#     qvals = np.asarray(qvals)
#     plt.figure()#figsize= (7.2*1.5,4.45*1.5))
#     plt.scatter(qvals,slice_avg)
#     # plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
#     plt.ylabel('Intensity [A.U.]')
#     plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
#     # plt.legend(loc='best')
#     plt.xlim(qvals[0],qvals[-1])
#     if True:
# #         plt.gca().set_ylim(1E5, 2E8)
#         plt.yscale('log')
#     plt.show()
#     return(slice_avg)





#rewrite the pixel to q conversion
def convert_pixel_to_q(x,y):
    x,y=np.int(x),np.int(y)
    """for Fourier transform images, we know the conversion already."""
    q_range = np.linspace(-(lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm)), (lattice_size//2)*((2*np.pi)/(lattice_size*latticeUnit_to_nm)), num=lattice_size +1)
    
    q_x,q_y,q_z,q,alpha_f,qx_prime,qy_prime,qz_prime = 0,q_range[x],q_range[y],(q_range[x]**2+q_range[x]**2)**(0.5),0,0,q_range[x],q_range[y]
    
    return(q_x,q_y,q_z,q,alpha_f,qx_prime,qy_prime,qz_prime) 

#plot the lattaice image
'''
The following functions work after loading the class
'''
def plot_lattice_image(frames=[], slice = False, x = [], y = [], save = 0):
    image = FD.rdlattice(frames[0])
    for i in frames[1:]:
        image += FD.rdlattice(i)
    fig, ax = plt.subplots(figsize=(7,7))
    dimension = np.shape(image)[0]
    rms = np.sqrt(np.mean(image**2))
    ax = plt.gca()
    print(np.shape(image))
    
    
    im = ax.imshow(image, vmin = -FD.height_range() , vmax = FD.height_range(), origin = 'lower', extent = [0,latticeUnit_to_nm*dimension,0,latticeUnit_to_nm*dimension], interpolation = 'none')#cmap=cmap_albula) 
    ax.set_ylabel(r'$L$ [nm]')
    ax.set_xlabel(r'$L$ [nm]')
    
    #plot the inset
    if slice == True:
        # These are in unitless percentages of the figure size. (0,0 is bottom left)
        left, bottom, width, height = [0.34, 0.25, 0.50, 0.25]
        ax2 = fig.add_axes([left, bottom, width, height])
        dimension = np.shape(image)[0]
    #     if logscale:
    #         scattering = np.log(scattering)
        rms = np.sqrt(np.mean(image**2))

        img = image[0:dimension,0:dimension] ####

        if isinstance(y,list):
            slice_avg = np.average(img[y[0]:y[1],x[0]:x[1]], axis = 0)
            y_avg = y[0]
        else:
            slice_avg = img[y,x[0]:x[1]]
            y_avg = y
        x_pixels = np.arange(x[0],x[1])
        L_vals = x_pixels
        ax2.plot(L_vals,slice_avg)
        # plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
        ax2.set_ylabel('Height [nm]')
        ax2.set_xlabel('L [nm]')
        ax2.set_title('Cross Section')
    
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    if save == 0:
        ax.set_title('Lattice image')
#     plt.suptitle('A (c0x) = {0}'.format(metadata[0]) +
#                  r', $\nu_y$ (c1y) = {0}'.format(metadata[1]) +
#                  r', $\nu_x$ (c1x) = {0}'.format(metadata[2]) +
#                  ', ' + r'$\kappa$ (c4) = {0}'.format(metadata[3]) + 
#                  r', $\lambda_x$ (c2x) = {0}'.format(metadata[4])+
#                  r', $\lambda_y$ (c2y) = {0}'.format(metadata[5])+
#                  ', \n ' + r'$\gamma_y$ (c3x) = {0}'.format(metadata[6]) +
#                  r', $\eta$ (c3) = {0}'.format(metadata[7]) +
#                  r', $\Delta$t (delt) = {0}'.format(metadata[8])+
#                  ', niter = {0}'.format(metadata[9])+
#                  ', nsteps = {0}'.format(number_of_images)+
#                  ', seed = {0}'.format(metadata[10]), fontsize=10)
    if save == 1:
        plt.savefig(filename+'.eps',bbox_inches='tight')
    plt.show()
    print('RMS of lattice image is ',np.sqrt(np.mean(image**2)))

def slice_lattice_img(x,y,frames=[0]):
    img = FD.rdlattice(frames[0])
    for i in frames[1:]:
        img += FD.rdlattice(i)
    dimension = np.shape(img)[0]
#     if logscale:
#         scattering = np.log(scattering)
    rms = np.sqrt(np.mean(img**2))
    
    img = img[0:dimension,0:dimension] ####
    
    if isinstance(y,list):
        slice_avg = np.average(img[y[0]:y[1],x[0]:x[1]], axis = 0)
        y_avg = y[0]
    else:
        slice_avg = img[y,x[0]:x[1]]
        y_avg = y
    x_pixels = np.arange(x[0],x[1])
    L_vals = x_pixels
    
    
    plt.figure(figsize= (7.2*1.5,4.45*1.5))
    plt.plot(L_vals,slice_avg)
    # plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
    plt.ylabel('Height [nm]')
    plt.xlabel('L [nm]')
#     plt.suptitle('A (c0x) = {0}'.format(metadata[0]) +
#                  r', $\nu_y$ (c1y) = {0}'.format(metadata[1]) +
#                  r', $\nu_x$ (c1x) = {0}'.format(metadata[2]) +
#                  ', ' + r'$\kappa$ (c4) = {0}'.format(metadata[3]) + 
#                  r', $\lambda_x$ (c2x) = {0}'.format(metadata[4])+
#                  r', $\lambda_y$ (c2y) = {0}'.format(metadata[5])+
#                  ', \n ' + r'$\gamma_y$ (c3x) = {0}'.format(metadata[6]) +
#                  r', $\eta$ (c3) = {0}'.format(metadata[7]) +
#                  r', $\Delta$t (delt) = {0}'.format(metadata[8])+
#                  ', niter = {0}'.format(metadata[9])+
#                  ', nsteps = {0}'.format(number_of_images)+
#                  ', seed = {0}'.format(metadata[10]), fontsize=10)
    # plt.legend(loc='best')

    plt.show()

#plot scattering pattern
def plot_scattering_image(frames=[], logscale = False, pixels= False,slice_pt=[0,0,0,0]): ####
    scattering = FD.rdframe(frames[0])
    
    #determine extent
    y1,y2,x1,x2 = -(lattice_size//2),(lattice_size//2),-(lattice_size//2),(lattice_size//2)
    for i in frames[1:]:
        scattering += FD.rdframe(i)
        
    if slice_pt!=[0,0,0,0]:
        scattering = scattering[slice_pt[0]:slice_pt[1],slice_pt[2]:slice_pt[3]]
        y1,y2,x1,x2 = -(lattice_size//2-slice_pt[0]),(slice_pt[1]-lattice_size//2), -(lattice_size//2-slice_pt[2]),(slice_pt[3]-lattice_size//2)

    plt.figure(figsize=(7,7)) #figsize=(10*((slice_pt[1]-slice_pt[0])/lattice_size),10*((slice_pt[3]-slice_pt[2])/lattice_size)))
    print(np.shape(scattering))
    
#     dimension = np.shape(scattering)[0]


#     if logscale:
#         scattering = np.log(scattering)
    rms = np.sqrt(np.average(scattering)**2)
    avg_img = np.average(scattering)
    maxval = np.amax(scattering)
    
#     scattering = scattering[0:dimension,0:dimension] ####
    
    ax = plt.gca()
    if logscale:
        if pixels:
            im = ax.imshow(scattering, origin = 'lower', interpolation = 'none',norm=LogNorm(vmin=0.01, vmax=maxval))#cmap=cmap_albula) 
        else:
            im = ax.imshow(scattering, origin = 'lower', \
                           extent = [y1*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     y2*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     x1*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     x2*((2*np.pi)/(lattice_size*latticeUnit_to_nm))],\
                           interpolation = 'none',norm=LogNorm(vmin=0.01, vmax=maxval))#cmap=cmap_albula) 
    else:
        if pixels:
            im = ax.imshow(scattering, vmin = 0.01 , vmax = maxval, origin = 'lower')#cmap=cmap_albula) 
        else:
            im = ax.imshow(scattering, vmin = 0.01 , vmax = maxval, origin = 'lower', \
                           extent = [y1*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     y2*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     x1*((2*np.pi)/(lattice_size*latticeUnit_to_nm)),\
                                     x2*((2*np.pi)/(lattice_size*latticeUnit_to_nm))],\
                           interpolation = 'none',norm=LogNorm(vmin=0.01, vmax=maxval))#cmap=cmap_albula) 

    ax.set_ylabel(r'$q_y$ $[nm^{-1}]$')
    ax.set_xlabel(r'$q_x$ $[nm^{-1}]$')
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Scattering image')
#     plt.suptitle('A (c0x) = {0}'.format(metadata[0]) +
#                  r', $\nu_y$ (c1y) = {0}'.format(metadata[1]) +
#                  r', $\nu_x$ (c1x) = {0}'.format(metadata[2]) +
#                  ', ' + r'$\kappa$ (c4) = {0}'.format(metadata[3]) + 
#                  r', $\lambda_x$ (c2x) = {0}'.format(metadata[4])+
#                  r', $\lambda_y$ (c2y) = {0}'.format(metadata[5])+
#                  ', \n ' + r'$\gamma_y$ (c3x) = {0}'.format(metadata[6]) +
#                  r', $\eta$ (c3) = {0}'.format(metadata[7]) +
#                  r', $\Delta$t (delt) = {0}'.format(metadata[8])+
#                  ', niter = {0}'.format(metadata[9])+
#                  ', nsteps = {0}'.format(number_of_images)+
#                  ', seed = {0}'.format(metadata[10]), fontsize=10)
#     if slice_pt!=[0,0,0,0]:
#         plt.gca().axis([slice_pt[0],slice_pt[1],slice_pt[2],slice_pt[3]])
    plt.show()
    print('RMS of scattering image is ',rms)
    print('Max of scattering image is ',np.amax(scattering))
    print('Avg of scattering image is ',np.average(scattering))
    
def slice_scat_img(x,y,frames=[0], fold=False):
    scattering = FD.rdframe(frames[0])
    for i in frames[1:]:
        scattering += FD.rdframe(i)
#     dimension = np.shape(scattering)[0]
#     if logscale:
#         scattering = np.log(scattering)
#     rms = np.sqrt(np.mean(scattering)**2)
    
#     scattering = scattering[0:dimension,0:dimension] ####
    
    if isinstance(y,list):
        slice_avg = np.average(scattering[y[0]:y[1],x[0]:x[1]], axis = 0)
        y_avg = y[0]
    else:
        slice_avg = scattering[y,x[0]:x[1]]
        y_avg = y
    x_pixels = np.arange(x[0],x[1])
    qvals = []
    rms = np.sqrt(np.mean(slice_avg)**2) # [:len(slice_avg)//2-4]
    maxval = np.amax(slice_avg)
    print(rms)
    for i in range(len(x_pixels)):
        qvals.append(convert_pixel_to_q(x_pixels[i],y_avg)[1])  
    qvals = np.asarray(qvals)
    
    if fold == True:
        size = len(slice_avg)
        slice_avg1 = slice_avg[size//2:]
        qvals1 = qvals[size//2:]
        
        slice_avg2 = slice_avg[:size//2]
        qvals2 = qvals[:size//2]*(-1)
    
    plt.figure(figsize= (7.2,4.45))
    if fold == True:
        plt.scatter(qvals1,slice_avg1, label='right')
        plt.scatter(qvals2,slice_avg2, label='left')
    else:    
        plt.scatter(qvals,slice_avg)
    # plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
    plt.ylabel('Intensity [A.U.]')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    # plt.legend(loc='best')
    if fold == True:
        plt.xlim(qvals1[0],qvals1[-1])
    else:
        plt.xlim(qvals[0],qvals[-1])
    if True:
        plt.yscale('log')
        plt.gca().set_ylim(0.01, 2*maxval)
#     plt.suptitle('A (c0x) = {0}'.format(metadata[0]) +
#                  r', $\nu_y$ (c1y) = {0}'.format(metadata[1]) +
#                  r', $\nu_x$ (c1x) = {0}'.format(metadata[2]) +
#                  ', ' + r'$\kappa$ (c4) = {0}'.format(metadata[3]) + 
#                  r', $\lambda_x$ (c2x) = {0}'.format(metadata[4])+
#                  r', $\lambda_y$ (c2y) = {0}'.format(metadata[5])+
#                  ', \n ' + r'$\gamma_y$ (c3x) = {0}'.format(metadata[6]) +
#                  r', $\eta$ (c3) = {0}'.format(metadata[7]) +
#                  r', $\Delta$t (delt) = {0}'.format(metadata[8])+
#                  ', niter = {0}'.format(metadata[9])+
#                  ', nsteps = {0}'.format(number_of_images)+
#                  ', seed = {0}'.format(metadata[10]), fontsize=10)
    plt.legend(loc='best')
    plt.show()
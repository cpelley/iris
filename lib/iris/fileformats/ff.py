# (C) British Crown Copyright 2010 - 2013, Met Office
#
# This file is part of Iris.
#
# Iris is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Iris is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Iris.  If not, see <http://www.gnu.org/licenses/>.
"""
Provides UK Met Office Fields File (FF) format specific capabilities.

"""

import os

import numpy as np

from iris.exceptions import NotYetImplementedError
from iris.fileformats.manager import DataManager
import pp


FF_HEADER_DEPTH = 256  # In words (64-bit).
FF_WORD_DEPTH = 8      # In bytes.

# UM marker to signify empty lookup table entry.
_FF_LOOKUP_TABLE_TERMINATE = -99

# UM FieldsFile fixed length header names and positions.
UM_FIXED_LENGTH_HEADER = [
        ('data_set_format_version',    (1, )),
        ('sub_model',                  (2, )),
        ('vert_coord_type',            (3, )),
        ('horiz_grid_type',            (4, )),
        ('dataset_type',               (5, )),
        ('run_identifier',             (6, )),
        ('experiment_number',          (7, )),
        ('calendar',                   (8, )),
        ('grid_staggering',            (9, )),
        ('time_type',                  (10, )),
        ('projection_number',          (11, )),
        ('model_version',              (12, )),
        ('obs_file_type',              (14, )),
        ('last_fieldop_type',          (15, )),
        ('first_validity_time',        (21, 22, 23, 24, 25, 26, 27, )),
        ('last_validity_time',         (28, 29, 30, 31, 32, 33, 34, )),
        ('misc_validity_time',         (35, 36, 37, 38, 39, 40, 41, )),
        ('integer_constants',          (100, 101, )),
        ('real_constants',             (105, 106, )),
        ('level_dependent_constants',  (110, 111, 112, )),
        ('row_dependent_constants',    (115, 116, 117, )),
        ('column_dependent_constants', (120, 121, 122, )),
        ('fields_of_constants',        (125, 126, 127, )),
        ('extra_constants',            (130, 131, )),
        ('temp_historyfile',           (135, 136, )),
        ('compressed_field_index1',    (140, 141, )), 
        ('compressed_field_index2',    (142, 143, )),
        ('compressed_field_index3',    (144, 145, )),
        ('lookup_table',               (150, 151, 152, )),
        ('total_prognostic_fields',    (153, )),
        ('data',                       (160, 161, 162, )), ]

# Offset value to convert from UM_FIXED_LENGTH_HEADER positions to
# FF_HEADER offsets.
UM_TO_FF_HEADER_OFFSET = 1
# Offset the UM_FIXED_LENGTH_HEADER positions to FF_HEADER offsets.
FF_HEADER = [
    (name, tuple(position - UM_TO_FF_HEADER_OFFSET for position in positions))
        for name, positions in UM_FIXED_LENGTH_HEADER]

# UM marker to signify a null pointer address.
_FF_HEADER_POINTER_NULL = 0
# UM FieldsFile fixed length header pointer names.
_FF_HEADER_POINTERS = [
        'integer_constants',
        'real_constants',
        'level_dependent_constants',
        'row_dependent_constants',
        'column_dependent_constants',
        'fields_of_constants',
        'extra_constants',
        'temp_historyfile',
        'compressed_field_index1',
        'compressed_field_index2',
        'compressed_field_index3',
        'lookup_table',
        'data', ]

_LBUSER_DTYPE_LOOKUP = {1: '>f{word_depth}',
                        2: '>i{word_depth}',
                        3: '>i{word_depth}',
                        'default': '>f{word_depth}', }


class FFHeader(object):
    """A class to represent the FIXED_LENGTH_HEADER section of a FieldsFile."""
    
    def __init__(self, filename, word_depth=FF_WORD_DEPTH):
        """
        Create a FieldsFile header instance by reading the
        FIXED_LENGTH_HEADER section of the FieldsFile.
        
        Args:
        
        * filename (string):
            Specify the name of the FieldsFile.
            
        Returns:
            FFHeader object.
            
        """
        
        self.ff_filename = filename
        '''File name of the FieldsFile.'''
        self._word_depth = word_depth

        # Read the FF header data
        with open(filename, 'rb') as ff_file:
            # typically 64-bit words (aka. int64 or ">i8")
            header_data = np.fromfile(ff_file,
                                      dtype='>i{0}'.format(self._word_depth),
                                      count=FF_HEADER_DEPTH)
            header_data = tuple(header_data)
            # Create FF instance attributes
            for name, offsets in FF_HEADER:
                if len(offsets) == 1:
                    value = header_data[offsets[0]]
                else:
                    value = header_data[offsets[0]:offsets[-1] + 1]
                setattr(self, name, value)

    def __str__(self):
        attributes = []
        for name, offsets in FF_HEADER:
            attributes.append('    {}: {}'.format(name, getattr(self, name)))
        return 'FF Header:\n' + '\n'.join(attributes)

    def __repr__(self):
        return '{}({!r})'.format(type(self).__name__, self.ff_filename)

    def valid(self, name):
        """
        Determine whether the FieldsFile FIXED_LENGTH_HEADER pointer attribute
        has a valid FieldsFile address.
        
        Args:
        
        * name (string):
            Specify the name of the FIXED_LENGTH_HEADER attribute.
            
        Returns:
            Boolean.
        
        """
        
        if name in _FF_HEADER_POINTERS:
            value = getattr(self, name)[0] > _FF_HEADER_POINTER_NULL
        else:
            msg = '{!r} object does not have pointer attribute {!r}'
            raise AttributeError(msg.format(self.__class__.__name__, name))
        return value

    def address(self, name):
        """
        Return the byte address of the FieldsFile FIXED_LENGTH_HEADER pointer attribute.
        
        Args:
        
        * name (string):
            Specify the name of the FIXED_LENGTH_HEADER attribute.
            
        Returns:
            int.
        
        """
        
        if name in _FF_HEADER_POINTERS:
            value = getattr(self, name)[0] * self._word_depth
        else:
            msg = '{!r} object does not have pointer attribute {!r}'
            raise AttributeError(msg.format(self.__class__.__name__, name))
        return value
    
    def shape(self, name):
        """
        Return the dimension shape of the FieldsFile FIXED_LENGTH_HEADER pointer attribute.
        
        Args:
        
        * name (string):
            Specify the name of the FIXED_LENGTH_HEADER attribute.
            
        Returns:
            Dimension tuple.
        
        """
        
        if name in _FF_HEADER_POINTERS:
            value = getattr(self, name)[1:]
        else:
            msg = '{!r} object does not have pointer address {!r}'
            raise AttributeError(msg.format(self.__class_.__name__, name))
        return value


class FF2PP(object):
    """A class to extract the individual PPFields from within a FieldsFile."""

    def __init__(self, filename, read_data=False, word_depth=FF_WORD_DEPTH):
        """
        Create a FieldsFile to Post Process instance that returns a generator
        of PPFields contained within the FieldsFile.
        
        Args:
        
        * filename (string):
            Specify the name of the FieldsFile.
            
        Kwargs:
        
        * read_data (boolean):
            Specify whether to read the associated PPField data within the FieldsFile.
            Default value is False.
            
        Returns:
            PPField generator.
        
        For example::
    
            >>> for field in ff.FF2PP(filename):
            ...     print field
            
        """
        
        self._ff_header = FFHeader(filename, word_depth=word_depth)
        self._word_depth = word_depth
        self._filename = filename
        self._read_data = read_data

    def _payload(self, field):
        """Calculate the payload data depth (in bytes) and type."""
        if field.lbpack.n1 == 0:
            # Data payload is not packed.
            data_depth = (field.lblrec - field.lbext) * self._word_depth
            # Determine PP field 64-bit payload datatype.
            lookup = _LBUSER_DTYPE_LOOKUP
            dtype_template = lookup.get(field.lbuser[0], lookup['default'])
            data_type = np.dtype(dtype_template.format(word_depth=self._word_depth))
        else:
            # Data payload is packed.
            if field.lbpack.n1 == 1:
                # Data packed using WGDOS archive method.
                data_depth = ((field.lbnrec * 2) - 1) * pp.PP_WORD_DEPTH
            elif field.lbpack.n1 == 2:
                # Data packed using CRAY 32-bit method.
                data_depth = (field.lblrec - field.lbext) * pp.PP_WORD_DEPTH
            else:
                msg = 'PP fields with LBPACK of {} are not supported.'
                raise NotYetImplementedError(msg.format(field.lbpack))

            # Determine PP field payload datatype.
            lookup = pp.LBUSER_DTYPE_LOOKUP
            data_type = lookup.get(field.lbuser[0], lookup['default'])

        return data_depth, data_type
        
    def _extract_field(self):
        # FF table pointer initialisation based on FF LOOKUP table configuration. 
        table_index, table_entry_depth, table_count = self._ff_header.lookup_table
        table_offset = (table_index - 1) * self._word_depth       # in bytes
        table_entry_depth = table_entry_depth * self._word_depth  # in bytes
        # Open the FF for processing.
        ff_file = open(self._ff_header.ff_filename, 'rb')
        ff_file_seek = ff_file.seek

        # Check for an instantaneous dump.
        if self._ff_header.dataset_type == 1:
            table_count = self._ff_header.total_prognostic_fields

        # Process each FF LOOKUP table entry.
        while table_count:
            table_count -= 1
            # Move file pointer to the start of the current FF LOOKUP table entry.
            ff_file_seek(table_offset, os.SEEK_SET)
            # Read the current PP header entry from the FF LOOKUP table.
            header_integers = np.fromfile(ff_file, dtype='>i{0}'.format(self._word_depth),
                                          count=pp.NUM_LONG_HEADERS)
            header_floats = np.fromfile(ff_file, dtype='>f{0}'.format(self._word_depth),
                                        count=pp.NUM_FLOAT_HEADERS)
            # In 64-bit words.
            header_data = tuple(header_integers) + tuple(header_floats)
            # Check whether the current FF LOOKUP table entry is valid.
            if header_data[0] == _FF_LOOKUP_TABLE_TERMINATE:
                # There are no more FF LOOKUP table entries to read. 
                break
            # Calculate next FF LOOKUP table entry.
            table_offset += table_entry_depth
            # Construct a PPField object and populate using the header_data
            # read from the current FF LOOKUP table.
            # (The PPField sub-class will depend on the header release number.)
            field = pp.make_pp_field(header_data)
            # Calculate start address of the associated PP header data.
            data_offset = field.lbegin * self._word_depth
            # Determine PP field payload depth and type.
            data_depth, data_type = self._payload(field)
            # Determine PP field data shape.
            data_shape = (field.lbrow, field.lbnpt)
            # Determine whether to read the associated PP field data.
            if self._read_data:
                # Move file pointer to the start of the current PP field data.
                ff_file_seek(data_offset, os.SEEK_SET)
                # Get the PP field data.
                data = field.read_data(ff_file, data_depth, data_shape, data_type)
                field._data = data
                field._data_manager = None
            else:
                proxy = pp.PPDataProxy(self._filename, data_offset,
                                       data_depth, field.lbpack)
                field._data = np.array(proxy)
                field._data_manager = DataManager(data_shape, data_type, field.bmdi)
            yield field
        ff_file.close()
        return
        
    def __iter__(self):
        return self._extract_field()


class FF32bit2PP(FF2PP):
    """
    Specialises the 64bit FF2PP class to handle fieldfiles
    converted to 32bit with the ieee tool.

    """
    def __init__(self, filename, read_data=False, word_depth=4):
        FF2PP.__init__(self, filename, read_data=False, word_depth=word_depth)


def load_cubes(filenames, callback):
    """
    Loads cubes from a list of fields files filenames.
    
    Args:
    
    * filenames - list of fields files filenames to load
    
    Kwargs:
    
    * callback - a function which can be passed on to :func:`iris.io.run_callback`
    
    .. note::

        The resultant cubes may not be in the order that they are in the file (order 
        is not preserved when there is a field with orography references).
         
    """
    return pp._load_cubes_variable_loader(filenames, callback, FF2PP)


def load_cubes_32bit_ieee(filenames, callback):
    """
    Loads cubes from a list of 32bit ieee converted fieldsfiles filenames.

    .. seealso::

        :func:`load_cubes` for keyword details

    """
    return pp._load_cubes_variable_loader(filenames, callback, FF32bit2PP)
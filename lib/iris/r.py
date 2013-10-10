# (C) British Crown Copyright 2013, Met Office
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
Provide conversion to and from r data structures.

See also: http://rpy.sourceforge.net/rpy2.html

"""
from collections import OrderedDict

import numpy as np
import rpy2
import rpy2.robjects as robjects


def _as_rpy_coord(coord):
    index = coord.points
    return index


def as_data_frame(cube, copy=True):
    """
    Convert a 2D cube to an rpy2 DataFrame
                                                                      
    Args:

        * cube - The cube to convert to a Pandas DataFrame.

    Kwargs:

        * copy - Whether to make a copy of the data.
                 Defaults to True.

    """
    data = cube.data
    if isinstance(data, np.ma.MaskedArray):
        if not copy:
            raise ValueError("Masked arrays must always be copied.")
        data = data.astype('f').filled(np.nan)
    elif copy:
        data = data.copy()

    rows = columns = None
    if cube.coords(dimensions=[0]):
        columns = _as_rpy_coord(cube.coord(dimensions=[0]))
    if cube.coords(dimensions=[1]):
        rows = _as_rpy_coord(cube.coord(dimensions=[1]))

    if not columns:
        columns = range(data.shape[1])

    frame_data = OrderedDict(
        {columns[ind]: rpy2.robjects.vectors.FloatVector(data[:, ind]) for
         ind in xrange(data.shape[1])})

    data_frame = robjects.DataFrame(frame_data)

    if rows:
        data_frame.rownames = rows

    return data_frame
   

# (C) British Crown Copyright 2017, Met Office
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
"""Test function :func:`iris.util.sanitise_auxcoords`."""

from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa

# import iris tests first so that some things can be initialised before
# importing anything else
import iris.tests as tests

import iris
import numpy as np

from iris.util import sanitise_auxcoords


class TestAll(tests.IrisTest):
    def test(self):
        # Ensure that multidimensional coordinates are transposed where
        # necessary so that their dimension mapping is in increasing order.
        data = np.zeros((1, 2, 3, 4, 5, 6))
        cube = iris.cube.Cube(data)
        aux_coord = iris.coords.AuxCoord(np.zeros((3, 1, 5, 2)),
                                         long_name='bing')
        cube.add_aux_coord(aux_coord, (2, 0, 4, 1))

        sanitise_auxcoords(cube)

        coord = cube.coord('bing')
        self.assertEqual(coord.shape, (1, 2, 3, 5))
        self.assertEqual(cube.coord_dims(coord), (0, 1, 2, 4))


if __name__ == '__main__':
    tests.main()

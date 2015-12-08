# (C) British Crown Copyright 2014 - 2015, Met Office
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
"""Unit tests for the :class:`iris.coord_systems.CoordSystem` class."""

from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa

# Import iris.tests first so that some things can be initialised before
# importing anything else.
import iris.tests as tests

import numpy as np

import cartopy.crs as ccrs
from iris.coord_systems import CoordSystem
import iris.util as util


class DummyCoordSystem(CoordSystem):
    def __init__(self, var1, var2):
        self.var1 = var1
        self.var2 = var2

    def as_cartopy_crs(self):
        pass

    def as_cartopy_projection(self):
        pass



class Test___eq__(tests.IrisTest):
    def eq(self, s, o):
        res = not (s.__class__ != o.__class__ or
                   s.__dict__.keys() != o.__dict__.keys())
        if res:
            for key in s.__dict__:
                if np.isreal(s.__dict__[key]):
                    res = util.approx_equal(s.__dict__[key],
                                            o.__dict__[key])
                else:
                    res = s.__dict__[key] == o.__dict__[key]
                if res is False:
                    break
        return res

    def test_relaxed_equality(self):
        cs1 = DummyCoordSystem(323 + 1e-8, 323)
        cs2 = DummyCoordSystem(323.0, 323)
        self.assertTrue(cs1.__eq__(cs2))

    def test_mixed_types(self):
        # Ensure that passing an integer and float type is accepted
        cs1 = DummyCoordSystem(323, 323)
        cs2 = DummyCoordSystem(323 + 1e-8, 323)
        self.assertTrue(cs1.__eq__(cs2))

    def test_utils_approx_called(self):
        with mock.patch('iris.coord_systems.utils.approx_equal'


if __name__ == '__main__':
    tests.main()

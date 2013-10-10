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


# import iris tests first so that some things can be initialised before
# importing anything else
import iris.tests as tests

import unittest

try:
    import rpy2
except ImportError:
    rpy2 = None

skip_rpy2 = unittest.skipIf(rpy2 is None,
                            'Test(s) require "rpy2", '
                            'which is not available.')

if rpy2:
    import numpy as np

    from iris.cube import Cube
    import iris.r


@skip_rpy2
class TestAsDataFrame(tests.IrisTest):
    """
    Test converting cubes to an rpy2 dataframe.

    """
    def setUp(self):
        self.cube = iris.cube.Cube(np.arange(4).reshape(2, 2)) 

    def test_no_dim_coords(self):
        data_frame = iris.r.as_data_frame(self.cube)
        self.assertArrayEqual(data_frame, cube.data)

    def test_no_x_coord(self):
        pass

    def test_no_y_coord(self):
        pass

    def test_cube_complete(self):
        pass

    def test_cube_masked(self):
        pass

    def test_copy_true(self):
        pass
        # ????

    def test_copy_false(self):
        pass
        # ????


@skip_rpy2
class TestDataFrameAsCube(tests.IrisTest):
    """
    Test converting rpy2 dataframes into cubes.

    """
    pass

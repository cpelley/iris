# (C) British Crown Copyright 2013 - 2014, Met Office
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
"""Test function :func:`iris.util.gen_regular_target`."""

# import iris tests first so that some things can be initialised before
# importing anything else
import iris.tests as tests

import numpy as np

import iris
import iris.tests.stock as stock
from iris.util import gen_regular_target


class TestAll(tests.IrisTest):
    def test_global_unbounded_lat_lon(self):
        cube = stock.global_lat_lon((5, 5), with_bounds=False)
        ref_cube = stock.global_lat_lon((1, 1), with_bounds=False)
        shape = (4, 4)
        result = gen_regular_target(cube, ref_cube, shape)

        target = iris.cube.Cube(np.zeros((4, 4)))
        target.add_dim_coord(cube.coord(axis='y').copy(
                             [-72., -24.,  24.,  72.]), 0)
        target.add_dim_coord(cube.coord(axis='x').copy(
                             [-144.,  -48.,   48.,  144.]), 1)

        self.assertEqual(result, target)

    def test_global_unbounded_lon_lat(self):
        cube = stock.global_lat_lon((5, 5), with_bounds=False)
        lon = cube.coord(axis='x')
        lat = cube.coord(axis='y')
        cube.remove_coord(lon)
        cube.remove_coord(lat)
        cube.add_dim_coord(lon, 0)
        cube.add_dim_coord(lat, 1)

        ref_cube = stock.global_lat_lon((1, 1), with_bounds=False)
        shape = (4, 4)
        result = gen_regular_target(cube, ref_cube, shape)

        target = iris.cube.Cube(np.zeros((4, 4)))
        target.add_dim_coord(cube.coord(axis='y').copy(
                             [-72., -24.,  24.,  72.]), 1)
        target.add_dim_coord(cube.coord(axis='x').copy(
                             [-144.,  -48.,   48.,  144.]), 0)

        self.assertEqual(result, target)

    def test_global_bounded_lat_lon(self):
        cube = stock.global_lat_lon((5, 5), with_bounds=True)
        ref_cube = stock.global_lat_lon((1, 1))
        shape = (4, 4)
        result = gen_regular_target(cube, ref_cube, shape)

        target = iris.cube.Cube(np.zeros((4, 4)))
        target.add_dim_coord(cube.coord(axis='y').copy(
                             [-72., -24.,  24.,  72.]), 0)
        target.add_dim_coord(cube.coord(axis='x').copy(
                             [-144.,  -48.,   48.,  144.]), 1)

        import ipdb; ipdb.set_trace()
        self.assertEqual(result, target)

if __name__ == '__main__':
    tests.main()

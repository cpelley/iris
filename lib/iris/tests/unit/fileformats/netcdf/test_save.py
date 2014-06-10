# (C) British Crown Copyright 2014, Met Office
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
"""Unit tests for the `iris.fileformats.netcdf.save` function."""

# Import iris.tests first so that some things can be initialised before
# importing anything else.
import iris.tests as tests

import mock
import netCDF4 as nc
import numpy as np

import iris
from iris.coord_systems import GeogCS, TransverseMercator, RotatedGeogCS
from iris.coords import DimCoord
from iris.cube import Cube
from iris.fileformats.netcdf import save
import iris.tests.stock as stock


class Test_global_attributes(tests.IrisTest):
    def _simple_cube(self, dtype):
        data = np.arange(12, dtype=dtype).reshape(3, 4)
        points = np.arange(3, dtype=dtype)
        bounds = np.arange(6, dtype=dtype).reshape(3, 2)
        cube = Cube(data, 'air_pressure_anomaly')
        coord = DimCoord(points, bounds=bounds)
        cube.add_dim_coord(coord, 0)
        return cube

    def test_global_attributes_as_list(self):
        # Ensure that we can extend a global attributes list.
        cube = self._simple_cube('>f4')
        cube.attributes['Conventions'] = 'convention1 convention2, adwd'

        print 'before: ', cube.attributes
        with self.temp_filename('.nc') as nc_path:
            save(cube, nc_path, 'NETCDF4')
            ds = nc.Dataset(nc_path)
#            import ipdb; ipdb.set_trace()
            print ds.getncattr('Conventions')
            ds.close()

if __name__ == "__main__":
    tests.main()

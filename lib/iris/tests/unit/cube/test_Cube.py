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
"""Unit tests for the `iris.cube.Cube` class."""

# Import iris.tests first so that some things can be initialised before
# importing anything else.
import iris.tests as tests

import mock
import numpy as np

from iris import FUTURE
from iris.analysis import WeightedAggregator, Aggregator
from iris.cube import Cube
from iris.coords import AuxCoord, DimCoord


class Test_xml(tests.IrisTest):
    def test_checksum_ignores_masked_values(self):
        # Mask out an single element.
        data = np.ma.arange(12).reshape(3, 4)
        data[1, 2] = np.ma.masked
        cube = Cube(data)
        self.assertCML(cube)

        # If we change the underlying value before masking it, the
        # checksum should be unaffected.
        data = np.ma.arange(12).reshape(3, 4)
        data[1, 2] = 42
        data[1, 2] = np.ma.masked
        cube = Cube(data)
        self.assertCML(cube)


class Test_collapsed__warning(tests.IrisTest):
    def setUp(self):
        self.cube = Cube([[1, 2], [1, 2]])
        lat = DimCoord([1, 2], standard_name='latitude')
        lon = DimCoord([1, 2], standard_name='longitude')
        grid_lat = AuxCoord([1, 2], standard_name='grid_latitude')
        grid_lon = AuxCoord([1, 2], standard_name='grid_longitude')
        wibble = AuxCoord([1, 2], long_name='wibble')

        self.cube.add_dim_coord(lat, 0)
        self.cube.add_dim_coord(lon, 1)
        self.cube.add_aux_coord(grid_lat, 0)
        self.cube.add_aux_coord(grid_lon, 1)
        self.cube.add_aux_coord(wibble, 1)

    def _aggregator(self, uses_weighting):
        # Returns a mock aggregator with a mocked method (uses_weighting)
        # which returns the given True/False condition.
        aggregator = mock.Mock(spec=WeightedAggregator)
        aggregator.cell_method = None
        aggregator.uses_weighting = mock.Mock(return_value=uses_weighting)

        return aggregator

    def _assert_warn_collapse_without_weight(self, coords, warn):
        # Ensure that warning is raised.
        msg = "Collapsing spatial coordinate {!r} without weighting"
        for coord in coords:
            self.assertIn(mock.call(msg.format(coord)), warn.call_args_list)

    def _assert_nowarn_collapse_without_weight(self, coords, warn):
        # Ensure that warning is not rised.
        msg = "Collapsing spatial coordinate {!r} without weighting"
        for coord in coords:
            self.assertNotIn(mock.call(msg.format(coord)), warn.call_args_list)

    def test_lat_lon_noweighted_aggregator(self):
        # Collapse latitude coordinate with unweighted aggregator.
        aggregator = mock.Mock(spec=Aggregator)
        aggregator.cell_method = None
        coords = ['latitude', 'longitude']

        with mock.patch('warnings.warn') as warn:
            self.cube.collapsed(coords, aggregator, somekeyword='bla')

        self._assert_nowarn_collapse_without_weight(coords, warn)

    def test_lat_lon_weighted_aggregator(self):
        # Collapse latitude coordinate with weighted aggregator without
        # providing weights.
        aggregator = self._aggregator(False)
        coords = ['latitude', 'longitude']

        with mock.patch('warnings.warn') as warn:
            self.cube.collapsed(coords, aggregator)

        coords = filter(lambda coord: 'latitude' in coord, coords)
        self._assert_warn_collapse_without_weight(coords, warn)

    def test_lat_lon_weighted_aggregator_with_weights(self):
        # Collapse latitude coordinate with a weighted aggregators and
        # providing suitable weights.
        weights = np.array([[0.1, 0.5], [0.3, 0.2]])
        aggregator = self._aggregator(True)
        coords = ['latitude', 'longitude']

        with mock.patch('warnings.warn') as warn:
            self.cube.collapsed(coords, aggregator, weights=weights)

        self._assert_nowarn_collapse_without_weight(coords, warn)

    def test_lat_lon_weighted_aggregator_alt(self):
        # Collapse grid_latitude coordinate with weighted aggregator without
        # providing weights.  Tests coordinate matching logic.
        aggregator = self._aggregator(False)
        coords = ['grid_latitude', 'grid_longitude']

        with mock.patch('warnings.warn') as warn:
            self.cube.collapsed(coords, aggregator)

        coords = filter(lambda coord: 'latitude' in coord, coords)
        self._assert_warn_collapse_without_weight(coords, warn)

    def test_no_lat_weighted_aggregator_mixed(self):
        # Collapse grid_latitude and an unmatched coordinate (not lat/lon)
        # with weighted aggregator without providing weights.
        # Tests coordinate matching logic.
        aggregator = self._aggregator(False)
        coords = ['wibble']

        with mock.patch('warnings.warn') as warn:
            self.cube.collapsed(coords, aggregator)

        self._assert_nowarn_collapse_without_weight(coords, warn)


class Test_summary(tests.IrisTest):
    def test_cell_datetime_objects(self):
        # Check the scalar coordinate summary still works even when
        # iris.FUTURE.cell_datetime_objects is True.
        cube = Cube(0)
        cube.add_aux_coord(AuxCoord(42, units='hours since epoch'))
        with FUTURE.context(cell_datetime_objects=True):
            summary = cube.summary()
        self.assertIn('1970-01-02 18:00:00', summary)


class Test_is_compatible(tests.IrisTest):
    def setUp(self):
        self.test_cube = Cube([1.])
        self.other_cube = self.test_cube.copy()

    def test_noncommon_array_attrs_compatible(self):
        # Non-common array attributes should be ok.
        self.test_cube.attributes['array_test'] = np.array([1.0, 2, 3])
        self.assertTrue(self.test_cube.is_compatible(self.other_cube))

    def test_matching_array_attrs_compatible(self):
        # Matching array attributes should be ok.
        self.test_cube.attributes['array_test'] = np.array([1.0, 2, 3])
        self.other_cube.attributes['array_test'] = np.array([1.0, 2, 3])
        self.assertTrue(self.test_cube.is_compatible(self.other_cube))

    def test_different_array_attrs_incompatible(self):
        # Differing array attributes should make the cubes incompatible.
        self.test_cube.attributes['array_test'] = np.array([1.0, 2, 3])
        self.other_cube.attributes['array_test'] = np.array([1.0, 2, 777.7])
        self.assertFalse(self.test_cube.is_compatible(self.other_cube))


if __name__ == "__main__":
    tests.main()

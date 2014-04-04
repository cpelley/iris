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
"""Unit tests for :class:`iris.analysis._interpolator.LinearInterpolator`."""

# Import iris.tests first so that some things can be initialised before
# importing anything else.
import iris.tests as tests

from nose.tools import assert_raises_regexp
import numpy as np

import iris
import iris.exceptions
import iris.tests.stock as stock
from iris.analysis._interpolator import LinearInterpolator


class ThreeDimCube(tests.IrisTest):
    def setUp(self):
        cube = stock.simple_3d_w_multidim_coords()
        cube.add_aux_coord(iris.coords.DimCoord(range(2), 'height'), 0)
        cube.add_dim_coord(iris.coords.DimCoord(range(3), 'latitude'), 1)
        cube.add_dim_coord(iris.coords.DimCoord(range(4), 'longitude'), 2)
        self.data = np.arange(24).reshape(2, 3, 4).astype(np.float32)
        cube.data = self.data
        self.cube = cube


class Test___init__(ThreeDimCube):
    def test_properties(self):
        interpolator = LinearInterpolator(self.cube, ['latitude'])

        # Default extrapolation mode.
        self.assertEqual(interpolator.extrapolation_mode, 'linear')

        # Access to cube property of the LinearInterpolator instance.
        self.assertEqual(interpolator.cube, self.cube)

        # Access to the resulting coordinate which we are interpolating over.
        self.assertEqual(interpolator.coords,
                         [self.cube.coord('latitude')])


class Test___init____circular(tests.IrisTest):
    # Capture the behaviour of the interpolator initialisation based on the
    # circular flag.
    def setUp(self):
        # Load some data from a file so that we have some deferred data.
        self.cube = stock.simple_pp()

    # TO-DO: Do not touch the data for this case.
    def test_not_circular(self):
        self.cube.coord('longitude').circular = False
        interpolator = LinearInterpolator(self.cube, ['longitude'])
        self.assertFalse(interpolator.cube.has_lazy_data())

    def test_circular(self):
        self.cube.coord('longitude').circular = True
        interpolator = LinearInterpolator(self.cube, ['longitude'])
        self.assertFalse(interpolator.cube.has_lazy_data())


class Test___init____validation(ThreeDimCube):
    def test_interpolator_overspecified(self):
        # Over specification by means of interpolating over two coordinates
        # mapped to the same dimension.
        msg = 'Coordinates repeat a data dimension - '\
            'the interpolation would be over-specified'
        with self.assertRaisesRegexp(ValueError, msg):
            LinearInterpolator(self.cube, ['wibble', 'height'])

    def test_interpolator_overspecified_scalar(self):
        # Over specification by means of interpolating over one dimension
        # coordinate and a scalar coordinate (not mapped to a dimension).
        self.cube.add_aux_coord(iris.coords.AuxCoord(
            1, long_name='scalar'), None)

        msg = 'Coordinates repeat a data dimension - '\
            'the interpolation would be over-specified'
        with self.assertRaisesRegexp(ValueError, msg):
            LinearInterpolator(self.cube, ['wibble', 'scalar'])

    def test_interpolate_monotonic(self):
        # Decreasing monotonic.
        self.cube = self.cube[:, ::-1]
        self.cube.data = self.data
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])
        expected = self.data[:, 0:1, :]
        result = self.interpolator([[0]])
        self.assertArrayEqual(result.data, expected)

    def test_interpolate_non_monotonic(self):
        self.cube.add_aux_coord(iris.coords.AuxCoord(
            [0, 3, 2], long_name='non-monotonic'), 1)
        msg = ('Cannot interpolate over the non-monotonic coordinate '
               'non-monotonic.')
        with self.assertRaisesRegexp(ValueError, msg):
            LinearInterpolator(self.cube, ['non-monotonic'])


class Test___call____1D(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])

    def test_interpolate_bad_coord_name(self):
        with self.assertRaises(iris.exceptions.CoordinateNotFoundError):
            LinearInterpolator(self.cube, ['doesnt exist'])

    def test_interpolate_data_single(self):
        # Single sample point.
        result = self.interpolator([[1.5]])
        expected = self.data[:, 1:, :].mean(axis=1).reshape(2, 1, 4)
        self.assertArrayEqual(result.data, expected)

        foo_res = result.coord('foo').points
        bar_res = result.coord('bar').points
        expected_foo = self.cube[:, 1:, :].coord('foo').points.mean(
            axis=0).reshape(1, 4)
        expected_bar = self.cube[:, 1:, :].coord('bar').points.mean(
            axis=0).reshape(1, 4)

        self.assertArrayEqual(foo_res, expected_foo)
        self.assertArrayEqual(bar_res, expected_bar)

    def test_interpolate_data_multiple(self):
        # Multiple sample points for a single coordinate (these points are not
        # interpolated).
        result = self.interpolator([[1, 2]])
        self.assertArrayEqual(result.data, self.data[:, 1:3, :])

        foo_res = result.coord('foo').points
        bar_res = result.coord('bar').points
        expected_foo = self.cube[:, 1:, :].coord('foo').points
        expected_bar = self.cube[:, 1:, :].coord('bar').points

        self.assertArrayEqual(foo_res, expected_foo)
        self.assertArrayEqual(bar_res, expected_bar)

    def _interpolate_data_linear_extrapolation(self, result):
        expected = self.data[:, 0:1] - (self.data[:, 1:2] - self.data[:, 0:1])
        self.assertArrayEqual(result.data, expected)

    def test_interpolate_data_linear_extrapolation(self):
        # Sample point outside the coordinate range.
        result = self.interpolator([[-1]])
        self._interpolate_data_linear_extrapolation(result)

    def test_interpolate_data_default_extrapolation(self):
        # Sample point outside the coordinate range.
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='linear')
        result = interpolator([[-1]])
        self._interpolate_data_linear_extrapolation(result)

    def _extrapolation_dtype(self, dtype):
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        result = interpolator([[-1]])
        self.assertTrue(np.all(np.isnan(result.data)))

    def test_extrapolation_nan_float32(self):
        # Ensure np.nan in a float32 array results.
        self._extrapolation_dtype(np.float32)

    def test_extrapolation_nan_float64(self):
        # Ensure np.nan in a float64 array results.
        self._extrapolation_dtype(np.float64)

    def test_interpolate_data_error_on_extrapolation(self):
        msg = 'One of the requested xi is out of bounds in dimension 0'
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='error')
        with assert_raises_regexp(ValueError, msg):
            interpolator([[-1]])

    def test_interpolate_data_unsupported_extrapolation(self):
        msg = "Extrapolation mode 'unsupported' not supported"
        with assert_raises_regexp(ValueError, msg):
            LinearInterpolator(self.cube, ['latitude'],
                               extrapolation_mode='unsupported')

    def test_multi_points_array(self):
        # Providing a multidimensional sample points for a 1D interpolation.
        # i.e. points given for two coordinates where there are only one
        # specified.
        msg = 'Expected sample points for 1 coordinates, got 2.'
        with assert_raises_regexp(ValueError, msg):
            self.interpolator([[1, 2], [1]])

    def test_interpolate_data_dtype_casting(self):
        data = self.data.astype(int)
        self.cube.data = data
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])
        result = self.interpolator([[0.125]])
        self.assertEqual(result.data.dtype, np.float64)

    def test_default_collapse_scalar(self):
        interpolator = LinearInterpolator(self.cube, ['wibble'])
        result = interpolator([0])
        self.assertEqual(result.shape, (3, 4))

    def test_collapse_scalar(self):
        interpolator = LinearInterpolator(self.cube, ['wibble'])
        result = interpolator([0], collapse_scalar=True)
        self.assertEqual(result.shape, (3, 4))

    def test_no_collapse_scalar(self):
        interpolator = LinearInterpolator(self.cube, ['wibble'])
        result = interpolator([0], collapse_scalar=False)
        self.assertEqual(result.shape, (1, 3, 4))

    def test_unsorted_datadim_mapping(self):
        # Currently unsorted data dimension mapping is not supported as the
        # indexing is not yet clever enough to remap the interpolated
        # coordinates.
        self.cube.transpose((0, 2, 1))
        interpolator = LinearInterpolator(self.cube, ['latitude'])
        msg = 'Currently only increasing data_dims is supported.'
        with self.assertRaisesRegexp(NotImplementedError, msg):
            interpolator([0])


class Test___call____1D_circular(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.cube.coord('longitude')._points = np.linspace(0, 360, 4,
                                                           endpoint=False)
        self.cube.coord('longitude').circular = True
        self.cube.coord('longitude').units = 'degrees'
        self.interpolator = LinearInterpolator(self.cube, ['longitude'],
                                               extrapolation_mode='nan')

    def test_interpolate_data_fully_wrapped(self):
        expected = self.interpolator([[180, 270]])
        result = self.interpolator([[-180, -90]])
        self.assertArrayEqual(expected.data, result.data)

    def test_interpolate_data_partially_wrapped(self):
        expected = self.interpolator([[180, 90]])
        result = self.interpolator([[-180, 90]])
        self.assertArrayEqual(expected.data, result.data)

    def test_interpolate_data_fully_wrapped_twice(self):
        xs = np.linspace(-360, 360, 100)
        xs_not_wrapped = (xs + 360) % 360
        expected = self.interpolator([xs])
        result = self.interpolator([xs_not_wrapped])
        self.assertArrayEqual(expected.data, result.data)


class Test___call____1D_singlelendim(ThreeDimCube):
    def setUp(self):
        """
        thingness / (1)                     (wibble: 2; latitude: 1)
             Dimension coordinates:
                  wibble                           x            -
                  latitude                         -            x
             Auxiliary coordinates:
                  height                           x            -
                  bar                              -            x
                  foo                              -            x
             Scalar coordinates:
                  longitude: 0
        """
        ThreeDimCube.setUp(self)
        self.cube = self.cube[:, 0:1, 0]
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])

    def test_interpolate_data_linear_extrapolation(self):
        # Linear extrapolation of a single valued element.
        result = self.interpolator([[1001]])
        self.assertArrayEqual(result.data, self.cube.data)

    def test_interpolate_data_nan_extrapolation(self):
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        result = interpolator([[1001]])
        self.assertTrue(np.all(np.isnan(result.data)))

    def test_interpolate_data_nan_extrapolation_not_needed(self):
        # No extrapolation for a single length dimension.
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        result = interpolator([[0]])
        self.assertArrayEqual(result.data, self.cube.data)


class Test___call____masked(tests.IrisTest):
    def setUp(self):
        self.cube = stock.simple_4d_with_hybrid_height()
        mask = np.isnan(self.cube.data)
        mask[::3, ::3] = True
        self.cube.data = np.ma.masked_array(self.cube.data,
                                            mask=mask)

    def test_orthogonal_cube(self):
        interpolator = LinearInterpolator(self.cube, ['grid_latitude'])
        result_cube = interpolator([1])

        # Explicit mask comparison to ensure mask retention.
        # Masked value input
        self.assertTrue(self.cube.data.mask[0, 0, 0, 0])
        # Mask retention on output
        self.assertTrue(result_cube.data.mask[0, 0, 0])

        self.assertCML(result_cube, ('experimental', 'analysis',
                                     'interpolate', 'LinearInterpolator',
                                     'orthogonal_cube_with_factory.cml'))


class Test___call____2D(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.interpolator = LinearInterpolator(self.cube,
                                               ['latitude', 'longitude'])

    def test_interpolate_data(self):
        result = self.interpolator([[1, 2], [2]])
        expected = self.data[:, 1:3, 2:3]
        self.assertArrayEqual(result.data, expected)

        index = (slice(None), slice(1, 3, 1), slice(2, 3, 1))
        for coord in self.cube.coords():
            coord_res = result.coord(coord).points
            coord_expected = self.cube[index].coord(coord).points

            self.assertArrayEqual(coord_res, coord_expected)

    def test_orthogonal_points(self):
        result = self.interpolator([[1, 2], [1, 2]])
        expected = self.data[:, 1:3, 1:3]
        self.assertArrayEqual(result.data, expected)

        index = (slice(None), slice(1, 3, 1), slice(1, 3, 1))
        for coord in self.cube.coords():
            coord_res = result.coord(coord).points
            coord_expected = self.cube[index].coord(coord).points

            self.assertArrayEqual(coord_res, coord_expected)

    def test_multi_dim_coord_interpolation(self):
        msg = 'Interpolation coords must be 1-d for rectilinear interpolation.'
        with self.assertRaisesRegexp(ValueError, msg):
            interpolator = LinearInterpolator(self.cube, ['foo', 'bar'])
            interpolator([[15], [10]])


class Test___call____2D_non_contiguous(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        coords = ['height', 'longitude']
        self.interpolator = LinearInterpolator(self.cube, coords)

    def test_interpolate_data_multiple(self):
        result = self.interpolator([[1], [1, 2]])
        expected = self.data[1:2, :, 1:3]
        self.assertArrayEqual(result.data, expected)

        index = (slice(1, 2), slice(None), slice(1, 3, 1))
        for coord in self.cube.coords():
            coord_res = result.coord(coord).points
            coord_expected = self.cube[index].coord(coord).points

            self.assertArrayEqual(coord_res, coord_expected)

    def test_orthogonal_cube(self):
        result_cube = self.interpolator([np.int64([0, 1, 1]),
                                         np.int32([0, 1])])
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator', 'basic_orthogonal_cube.cml')
        self.assertCMLApproxData(result_cube, result_path)
        self.assertEqual(result_cube.coord('longitude').dtype, np.int32)
        self.assertEqual(result_cube.coord('height').dtype, np.int64)

    def test_orthogonal_cube_squash(self):
        result_cube = self.interpolator([np.int64(0),
                                         np.int32([0, 1])])
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator', 'orthogonal_cube_1d_squashed.cml')
        self.assertCMLApproxData(result_cube, result_path)
        self.assertEqual(result_cube.coord('longitude').dtype, np.int32)
        self.assertEqual(result_cube.coord('height').dtype, np.int64)

        non_collapsed_cube = self.interpolator([[np.int64(0)],
                                                np.int32([0, 1])],
                                               collapse_scalar=False)
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator',
                       'orthogonal_cube_1d_squashed_2.cml')
        self.assertCML(non_collapsed_cube[0, ...], result_path)
        self.assertCML(result_cube, result_path)
        self.assertEqual(result_cube, non_collapsed_cube[0, ...])


if __name__ == "__main__":
    tests.main()

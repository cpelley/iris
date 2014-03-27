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
"""Unit tests for :class:`iris.analysis.interpolate.LinearInterpolator`."""

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
        self.cube = cube
        self.data = np.arange(24).reshape(2, 3, 4).astype(np.float32)

        cube.data = self.data


class Test_LinearInterpolator_1D(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])

    def test_interpolate_bad_coord_name(self):
        with self.assertRaises(iris.exceptions.CoordinateNotFoundError):
            LinearInterpolator(self.cube, ['doesnt exist'])

    def test_interpolate_data_single(self):
        sample_points = [['latitude', [1.5]]]
        result = self.interpolator.points(sample_points, self.data)
        expected = self.data[:, 1:3, :].mean(axis=1).reshape(2, 1, 4)
        self.assertArrayEqual(result, expected)

    def test_interpolate_data_multiple(self):
        sample_points = [['latitude', [1, 2]]]
        result = self.interpolator.points(sample_points, self.data)
        self.assertArrayEqual(result, self.data[:, 1:3, :])

    def test_interpolate_data_linear_extrapolation(self):
        sample_points = [['latitude', [-1]]]
        result = self.interpolator.points(sample_points, self.data)
        expected = self.data[:, 0:1] - (self.data[:, 1:2] - self.data[:, 0:1])
        self.assertArrayEqual(result, expected)

    def _extrapolation_dtype(self, dtype):
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        data = self.data.astype(dtype)
        sample_points = [['latitude', [-1]]]
        result = interpolator.points(sample_points, data)
        self.assertTrue(np.all(np.isnan(result)))

    def test_extrapolation_nan_float32(self):
        self._extrapolation_dtype(np.float32)

    def test_extrapolation_nan_float64(self):
        self._extrapolation_dtype(np.float64)

    def test_interpolate_data_error_on_extrapolation(self):
        msg = 'One of the requested xi is out of bounds in dimension 0'
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='error')
        sample_points = [['latitude', [-1]]]
        with assert_raises_regexp(ValueError, msg):
            result = interpolator.points(sample_points, self.data)

    def test_bad_sample_points_array(self):
        sample_points = [['latitude', [1, 2]],
                         ['longitude', [1]]]
        with self.assertRaises(ValueError):
            self.interpolator.points(sample_points, self.data)

    def test_interpolate_data_dtype_casting(self):
        data = self.data.astype(int)
        sample_points = [['latitude', [0.125]]]
        result = self.interpolator.points(sample_points, data)
        self.assertEqual(result.dtype, np.float64)

    def test_interpolate_data_with_data_dims(self):
        sample_points = [['latitude', [1]]]
        result = self.interpolator.points(sample_points, self.data[0],
                                          data_dims=[1, 2])
        temp = self.data[0, 1].reshape(1, 1, 4)
        expected = np.vstack((temp, temp))
        self.assertArrayEqual(result, expected)

    def teet_inetrpolate_data_with_data_dims_not_interp(self):
        # No interpolation is taking place...
        result = self.interpolator.interpolate_data([1], self.data[:, 0],
                                                    data_dims=[0, 2])
        self.assertArrayEqual(result, self.data[:, 0])


class SingleLengthDimension(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.cube = self.cube[:, 0:1, 0]
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])

    def test_interpolate_data_linear_extrapolation(self):
        # Linear extrapolation of a single valued element.
        result = self.interpolator.points([['latitude', [1001]]],
                                          self.cube.data)
        self.assertArrayEqual(result, self.cube.data)

    def test_interpolate_data_nan_extrapolation(self):
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        result = interpolator.points([['latitude', [1001]]], self.cube.data)
        self.assertTrue(np.all(np.isnan(result)))

    def test_interpolate_data_nan_extrapolation_not_needed(self):
        # No extrapolation for a single length dimension.
        interpolator = LinearInterpolator(self.cube, ['latitude'],
                                          extrapolation_mode='nan')
        result = interpolator.points([['latitude', [0]]], self.cube.data)
        self.assertArrayEqual(result, self.cube.data)

    def test_interpolator_overspecified(self):
        msg = 'Coordinates repeat a data dimension - '\
            'the interpolation would be over-specified'
        with self.assertRaisesRegexp(ValueError, msg):
            LinearInterpolator(self.cube, ['latitude', 'longitude'])


class Test_LinearInterpolator_monotonic(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.cube = self.cube[:, ::-1]
        self.cube.data = self.data
        self.interpolator = LinearInterpolator(self.cube, ['latitude'])

    def test_interpolate_data(self):
        sample_points = [['latitude', [0]]]
        expected = self.data[:, 0:1, :]
        result = self.interpolator.points(sample_points, self.data)
        self.assertArrayEqual(result, expected)


class Test_LinearInterpolator_circular(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.cube.coord('longitude').points = np.linspace(0, 360, 4,
                                                          endpoint=False)
        self.cube.coord('longitude').circular = True
        self.cube.coord('longitude').units = 'degrees'
        self.interpolator = LinearInterpolator(self.cube, ['longitude'],
                                               extrapolation_mode='nan')

    def test_interpolate_data_fully_wrapped(self):
        expected = self.interpolator.points([['longitude', [180, 270]]],
                                            self.data)
        result = self.interpolator.points([['longitude', [-180, -90]]],
                                          self.data)
        self.assertArrayEqual(expected, result)

    def test_interpolate_data_partially_wrapped(self):
        expected = self.interpolator.points([['longitude', [180, 90]]],
                                            self.data)
        result = self.interpolator.points([['longitude', [-180, 90]]],
                                          self.data)
        self.assertArrayEqual(expected, result)

    def test_interpolate_data_fully_wrapped_twice(self):
        xs = np.linspace(-360, 360, 100)
        xs_not_wrapped = (xs + 360) % 360
        expected = self.interpolator.points([['longitude', xs]], self.data)
        result = self.interpolator.points([['longitude', xs_not_wrapped]],
                                          self.data)
        self.assertArrayEqual(expected, result)

# XXX Test Masked data...


class Test_LinearInterpolator_masked_and_factory(tests.IrisTest):
    def setUp(self):
        self.cube = stock.simple_4d_with_hybrid_height()
        mask = np.isnan(self.cube.data)
        mask[::3, ::3] = True
        self.cube.data = np.ma.masked_array(self.cube.data,
                                            mask=mask)

    def test_orthogonal_cube(self):
        interpolator = LinearInterpolator(self.cube, ['grid_latitude'])
        result_cube = interpolator([['grid_latitude', 1]])
        self.assertCML(result_cube, ('experimental', 'analysis',
                                     'interpolate', 'LinearInterpolator',
                                     'orthogonal_cube_with_factory.cml'))


class Test_LinearInterpolator_2D(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        self.interpolator = LinearInterpolator(self.cube,
                                               ['latitude', 'longitude'])

    def test_interpolate_data(self):
        sample_points = [['longitude', [2]],
                         ['latitude', [1, 2]]]
        result = self.interpolator.points(sample_points,
                                          self.cube.data)
        expected = self.data[:, 1:3, 2:3]
        self.assertArrayEqual(result, expected)

    def test_orthogonal_points(self):
        sample_points = [['longitude', [1, 2]],
                         ['latitude', [1, 2]]]
        result = self.interpolator.points(sample_points,
                                          self.cube.data)
        expected = self.data[:, 1:3, 1:3]
        self.assertArrayEqual(result, expected)

    def test_interpolate_data_data_dims(self):
        sample_points = [['longitude', [0]],
                         ['latitude', [0]]]
        result = self.interpolator.points(sample_points, self.data[:, 0],
                                          data_dims=[0, 2])
        expected = self.data[:, 0, 0].reshape(2, 1, 1)
        self.assertArrayEqual(result, expected)

    def test_interpolate_data_data_dims_multi_point(self):
        sample_points = [['longitude', [0]],
                         ['latitude', [0, 1, 2]]]
        result = self.interpolator.points(sample_points, self.data[:, 0],
                                          data_dims=[0, 2])
        expected = np.repeat(self.data[:, 0, 0].reshape(2, 1, 1), 3, axis=1)
        self.assertArrayEqual(result, expected)

    def test_orthogonal_points_data_dims(self):
        sample_points = [['latitude', 0],
                         ['longitude', 0]]
        result = self.interpolator.points(sample_points,
                                          self.data[:, 0],
                                          data_dims=[0, 2])
        expected = np.array([[[0]], [[12]]])
        assert np.all(result == expected)

    def test_orthogonal_points_data_dims_multiple_points(self):
        sample_points = [['latitude', [0, 1]],
                         ['longitude', [0, 1, 1]]]
        result = self.interpolator.points(sample_points,
                                          self.data[:, 0],
                                          data_dims=[0, 2])

        expected = np.concatenate([self.data[:, 0:1, 0], self.data[:, 0:1, 0],
                                   self.data[:, 0:1, 1], self.data[:, 0:1, 1],
                                   self.data[:, 0:1, 1], self.data[:, 0:1, 1]],
                                  axis=1).reshape(2, 3, 2).transpose([0, 2, 1])
        assert np.all(result == expected)

    def test_orthogonal_points_data_dims_multiple_points_transposed(self):
        sample_points = [['latitude', [0, 1, 1]],
                         ['longitude', [0, 1]]]
        result = self.interpolator.points(sample_points,
                                          self.data[:, 0],
                                          data_dims=[0, 2])
        expected = np.array([[[0,  1], [0,  1], [0,  1]],
                             [[12, 13], [12, 13], [12, 13]]])
        assert np.all(result == expected)


class Test_LinearInterpolator_2D_non_contiguous(ThreeDimCube):
    def setUp(self):
        ThreeDimCube.setUp(self)
        coords = ['height', 'longitude']
        self.interpolator = LinearInterpolator(self.cube, coords)

    def test_interpolate_data_multiple(self):
        sample_points = [['longitude', [1, 2]],
                         ['height', [1]]]
        result = self.interpolator.points(sample_points, self.data)
        expected = self.data[1:2, :, 1:3]
        self.assertArrayEqual(result, expected)

    def test_intepolate_data_wrong_data_shape(self):
        msg = 'Data being interpolated is not consistent with ' \
            'the data passed through'
        sample_points = [['longitude', [1]],
                         ['height', [0]]]
        with assert_raises_regexp(ValueError, msg):
            self.interpolator.points(sample_points, self.data[0])

    def test_interpolate_points_data_order(self):
        sample_points1 = [['longitude', [0]],
                          ['height', [0, 1]]]
        result1 = self.interpolator.points(sample_points1, self.cube.data)
        sample_points2 = [['height', [0, 1]],
                          ['longitude', [0]]]
        result2 = self.interpolator.points(sample_points2, self.cube.data)
        expected = self.cube.data[0:2, :, 0:1]
        assert np.all(result1 == expected), 'Wrong result'
        assert np.all(result1 == result2), 'Wrong result'

        coords = ['longitude', 'height']
        inverse_interpolator = LinearInterpolator(self.cube, coords)
        sample_points3 = [['longitude', [0]],
                          ['height', [0, 1]]]
        result3 = inverse_interpolator.points(sample_points3, self.cube.data)
        assert np.all(result1 == result3), 'Wrong result'
        sample_points4 = [['height', [0, 1]],
                          ['longitude', [0]]]
        result4 = inverse_interpolator.points(sample_points4, self.cube.data)
        assert np.all(result1 == result4), 'Wrong result'

    def test_orthogonal_cube(self):
        result_cube = self.interpolator([['height', np.int64([0, 1, 1])],
                                         ['longitude', np.int32([0, 1])]])
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator', 'basic_orthogonal_cube.cml')
        self.assertCMLApproxData(result_cube, result_path)
        self.assertEqual(result_cube.coord('longitude').dtype, np.int32)
        self.assertEqual(result_cube.coord('height').dtype, np.int64)

    def test_orthogonal_cube_scalar_value(self):
        cube = self.cube[0, 0, 0]
        with self.assertRaises(ValueError):
            LinearInterpolator(cube, ['latitude'])

    def test_orthogonal_cube_squash(self):
        result_cube = self.interpolator([['height', np.int64(0)],
                                         ['longitude', np.int32([0, 1])]])
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator', 'orthogonal_cube_1d_squashed.cml')
        self.assertCMLApproxData(result_cube, result_path)
        self.assertEqual(result_cube.coord('longitude').dtype, np.int32)
        self.assertEqual(result_cube.coord('height').dtype, np.int64)

        sample_points = [['height', np.int64(0)],
                         ['longitude', np.int32([0, 1])]]
        non_collapsed_cube = self.interpolator(sample_points,
                                               collapse_scalar=False)
        result_path = ('experimental', 'analysis', 'interpolate',
                       'LinearInterpolator',
                       'orthogonal_cube_1d_squashed_2.cml')
        self.assertCML(non_collapsed_cube[0, ...], result_path)
        self.assertCML(result_cube, result_path)
        self.assertEqual(result_cube, non_collapsed_cube[0, ...])


if __name__ == "__main__":
    tests.main()

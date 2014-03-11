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
"""
Classes for regridding and interpolating.

"""
import iris.analysis.interpolate as interpolate
import iris.analysis as analysis
import iris.experimental.regrid as exregrid


class _Regridder(object):
    def __init__(self, src_cube, target_cube=None, mask_tolerance=0.0, **kwargs):
        self._src_cube = src_cube
        self._target_cube = target_cube
        self._weights = weights
        self._mdtol = mask_tolerance

        self.setup()

    def setup(self):
        # Method that contructs the sparse matrix etc. and caches this
        # information in the class instance.
        return self._setup(src_cube, target_cube)

    def regrid(self, data, **kwargs):
        return self._regrid(data, **kwargs)


class _WeightedRegridder(_Regridder):
    def __init__(self, *args, method=analysis.MEAN, weights=None, **kwargs):
        self._weights = self._weights
        super(_WeightedRegridder, self).__init__(*args, **kwargs)

    @property
    def weights(self):
        return self._weights

    def _calculate_weights(self):
        return NotImplemented

    def setup(self):
        # Method that calculates the weights and contructs the sparse matrix
        # etc. and caches this information in the class instance.
        self._setup()
        if not weights:
            self._calculate_weights()


class _Interpolator(object):
    def __init__(self, src_cube, **kwargs):
        self._src_cube = src_cube
        self.setup(src_cube)

    def setup(self, src_cube):
        # Method that performs a caching of information and associates this
        # with the class instance.
        return NotImplemented

    def interpolate(self, sample_points, **kwargs):
        return _interpolate(sample_points, **kwargs)


class _NDInterpolator(Regridder):
    def __init__(interpolator):
        pass


class LinearInterpolator(_Interpolator):
    def _interpolate(self, sample_points, extrapolation='linear', **kwargs):
        return interpolate.linear(sample_points, extrapolation, **kwargs)


class LinearRegridder(_Regridder):
    def _regrid(self, sample_points, **kwargs):
        return interpolate.regrid(self, data, mode='bilinear', **kwargs)


class NearestInterpolator(_Interpolator):
    def _interpolate(self, sample_points, **kwargs):
        return interpolate.extract_nearest_neighbour(sample_points, **kwargs)

    def value(self, sample_points, **kwargs):
        return interpolate.nearest_neighbour_data_value(sample_points, **kwargs)

    def indices(self, sample_points, **kwargs):
        return interpolate.nearest_neighbour_indices(sample_points, **kwargs)


class NearestRegridder(_Regridder):
    def _regrid(self, data, **kwargs):
        return interpolate.regrid(data, mode='nearest', **kwargs)


class AreaOverlapRegridder(_WeightedRegridder):
    def _regrid(self, data, **kwargs):
        if method != 'conservative':
            return interpolate.regrid_area_weighted_rectilinear_src_and_grid(
                data, weights=self._weights, method=self._method, **kwargs)
        else:
            return interpolate.regrid_conservative(
                data, weights=self._weights, method=self._method, **kwargs)


class PointInCellRegridder(_WeightedRegridder):
    def _regrid(self, data, **kwargs):
        return interpolate.regrid_weighted_curvilinear_to_rectilinear(
            data, weights=self._weights, method=self._method, **kwargs)


def regrid(cube, grid_cube, *args, regridder=None, cached=None, **kwargs):
    if regridder and cached:
        raise TypeError('Either a regridder class or a regridder class '
            'instance should be provided, but not both.')
    if regridder:
        regridder = regridder(cube, grid_cube, *args, **kwargs)
    else:
        regridder = cached

    return regridder.regrid(cube, grid_cube.data)


def interpolate(cube, sample_points, *args, interpolator=None, cached=None,
        **kwargs):
    if interpolator and cached_interpolator:
        raise TypeError('Either an interpolator class or an interpolator '
            'class instance should be provided, but not both.')
    if interpolator:
        interpolator = interpolator(cube, grid_cube, *args, **kwargs)
    else:
        interpolator = cached

    return interpolator.interpolate(cube, sample_points)

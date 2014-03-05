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
    def __init__(self, src_cube, target_cube=None, weights=None,
            mask_tolerance=0.0, **kwargs):
        self._src_cube = src_cube
        self._target_cube = target_cube
        self._weights = weights
        self._mdtol = mask_tolerance

        self.setup(src_cube, target_cube)

    def setup(src_cube, target_cube):
        # Method that calculates weights, sparse matrix etc. and caches this
        # information in the class instance.
        return NotImplemented

    def regrid(self, data, **kwargs):
        return _regrid(data, **kwargs)


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
    def _interpolate(self, *args, extrapolation='linear'):
        return interpolate.linear(*args, extrapolation, **kwargs)


class LinearRegridder(_Regridder):
    def _regrid(self, *args, **kwargs):
        return interpolate.regrid(self, *args, mode='bilinear', **kwargs)


class NearestInterpolator(_Interpolator):
    def _interpolate(self, *args, **kwargs):
        return interpolate.extract_nearest_neighbour(*args, **kwargs)

    def value(self, *args, **kwargs):
        return interpolate.nearest_neighbour_data_value(*args, **kwargs)

    def indices(self, *args, **kwargs):
        return interpolate.nearest_neighbour_indices(*args, **kwargs)


class NearestRegridder(_Regridder):
    def _regrid(self, *args, **kwargs):
        return interpolate.regrid(*args, mode='nearest', **kwargs)


class AreaOverlapRegridder(_Regridder):
    def _regrid(self, data, method=analysis.MEAN, **kwargs):
        if method != 'conservative':
            return interpolate.regrid_area_weighted_rectilinear_src_and_grid(
                *args, weights=weights, method=method, **kwargs)
        else:
            return interpolate.regrid_conservative(*args, weights=weights,
                **kwargs)


class PointInCellRegridder(_Regridder):
    def _regrid(self, *args, weights=None, method=analysis.MEAN, **kwargs):
        if weights is None:
            weights = self._weights
        return interpolate.regrid_weighted_curvilinear_to_rectilinear(
            *args, weights=weights, method=method, **kwargs)


def regrid(cube, grid_cube, regridder=None, cached=None):
    if regridder and cached:
        raise TypeError('Either a regridder class or a regridder class '
            'instance should be provided, but not both.')
    if regridder:
        regridder = regridder(cube, grid_cube)
    else:
        regridder = cached

    return regridder.regrid(cube, grid_cube)


def interpolate(cube, sample_points, interpolator=None, cached=None):
    if interpolator and cached_interpolator:
        raise TypeError('Either an interpolator class or an interpolator '
            'class instance should be provided, but not both.')
    if interpolator:
        interpolator = interpolator(cube, grid_cube)
    else:
        interpolator = cached

    return interpolator.interpolate(cube, grid_cube)

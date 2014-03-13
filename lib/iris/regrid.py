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
Classes for regridding.

"""
from abc import ABCMeta, abstractmethod

import iris.analysis as analysis
import iris.analysis.interpolate as interpolate
import iris.experimental.regrid as exregrid


class _Regridder(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def setup():
        # Weights calculation and/sparse matrix calculation perhaps?
        pass

    @abstractmethod
    def regrid(self):
        pass


class Linear(_Regridder):
    def __init__(self, src_cube, target_cube):
        self._src_cube = src_cube
        self._target_cube = target_cube
        self.setup()

    def regrid(self, target_cube):
        # The two algorithms need to have their setup made uniform (similar
        # restrictions)
        # regrid_bilinear_rectilinear_src_and_grid:
        #     Ignores coordinates that share the same dimension as the
        #     horizontal coordinates.
        # regrid:
        #     requires that no coordinate share the same dimension as the
        # horizontal coordinates.
        if (target_cube.coord(axis='x').coord_system !=
                self.src_cube.coord(axis='x').coord_system):
            return exregrid.regrid_bilinear_rectilinear_src_and_grid(
                self._src_cube, target_cube)
        else:
            return interpolate.regrid(self._src_cube, target_cube,
                                      mode='bilinear')


class Nearest(_Regridder):
    def __init__(self, src_cube, target_cube):
        self._src_cube = src_cube
        self._target_cube = target_cube
        self.setup()

    def setup():
        # Some setup
        pass

    def regrid(self, target_cube):
        return interpolate.regrid(self._src_cube, target_cube, mode='nearest')


class AreaOverlap(_Regridder):
    def __init__(self, src_cube, target_cube, method=analysis.SUM,
                 weights=None):
        self._method = method
        self._weights = weights
        self._src_cube = src_cube
        self._target_cube = target_cube
        self.setup()

    def setup():
        # Currently weights are calculated by the underlying regrid algorithms.
        # We should give the ability to pass a precomputed weights array.
        pass

    def weights(self):
        return self._weights

    def regrid(self, target_cube, method=analysis.SUM):
        # Currently these two algorithms differ greatly in what they expect.
        # regrid_area_weighted_rectilinear_src_and_grid:
        #    Requires both grids to have the same coordinate systems (changing
        #    this would be major if calculating the weights)
        # regrid_conservative:
        #    Grids do not have to be the same as they are transformed to
        #    lat-lon before going into ESMPY.
        #    Assumes spherical earth!
        #    Calculation made to ESMPY in slices over other dimensions (ESMPY
        #    can cope with multidimensional)
        # Additional coordinates that map to the horizontal dimensions are
        # ignroed.
        if method != 'conservative':
            return interpolate.regrid_area_weighted_rectilinear_src_and_grid(
                self._src_cube, target_cube, weights=self._weights)
        else:
            return interpolate.regrid_conservative(
                self._src_cube, target_cube, weights=self._weights)


class PointInCell(_Regridder):
    def __init__(self, src_cube, target_cube, weights=None):
        self._weights = weights
        self._src_cube = src_cube
        self._target_cube = target_cube
        self.setup()

    def setup():
        # Currently weights are provided by the underlying regrid algorithm.
        # We should give the ability to calculate this weights array if not
        # provided (if this is even accurate enough for all projections).
        # Currently restricted to lat-lon grids anyway.
        pass

    def weights(self):
        return self._weights

    def regrid(self, target_cube):
        return interpolate.regrid_weighted_curvilinear_to_rectilinear(
            self._src_cube, target_cube, weights=self._weights)

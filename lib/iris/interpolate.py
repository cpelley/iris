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
Classes for and interpolating.

"""
from abc import ABCMeta, abstractmethod

import iris.analysis.interpolate as interpolate


class _Interpolator(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def setup():
        # Weights calculation and/sparse matrix calculation perhaps?
        pass

    @abstractmethod
    def interpolate(self):
        pass


class Linear(_Interpolator):
    def __init__(self, src_cube):
        self._src_cube = src_cube
        self.setup()

    def setup(self):
        pass


class Nearest(_Interpolator):
    def __init__(self, src_cube):
        self._src_cube = src_cube
        self.setup()

    def setup(self):
        pass

    def interpolate(self, sample_points):
        return interpolate.extract_nearest_neighbour(
            self._src_cube, sample_points)

    def value(self, sample_points):
        return interpolate.nearest_neighbour_data_value(sample_points)

    def indices(self, sample_points):
        return interpolate.nearest_neighbour_indices(sample_points)

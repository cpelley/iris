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


class LINEAR(object):
    def __init__(self, mdtol=None, weights=None):
        self._mdtol = mdtol
        self._weights = weights

    def regrid(self, src_cube, target_cube):
        return interpolate.regrid(src_cube, target_cube, mode='bilinear')


def regrid(cube, grid_cube, regrid_method):
    return regrid_method.regrid(cube, grid_cube)




class Regridder(object):
    def __init__(self, mdtol, weights):
        self._mdtol = mdtol
        self._weights = weights

    def regrid(self):
        return NotImplemented


class LINEARv2(Regridder):
    def regrid(self, src_cube, target_cube):
        return interpolate.regrid(src_cube, target_cube, mode='bilinear')


class NEARESTv2(Regridder):
    def regrid(self, src_cube, target_cube):
        return interpolate.regrid(src_cube, target_cube, mode='nearest')


def regridv2(cube, grid_cube, regrid_method, mdtol=None, weights=None):
    regrid_object = regrid_method(mdtol, weights)
    return regrid_object.regrid(cube, grid_cube)

# (C) British Crown Copyright 2010 - 2013, Met Office
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
Basic mathematical and statistical operations.

"""
from __future__ import division
import warnings
import math

import numpy as np

import iris.analysis
import iris.coords
import iris.cube
import iris.exceptions


def abs(cube, in_place=False):
    """
    Calculate the absolute values of the data in the Cube provided.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _math_op_common(cube, np.abs, cube.units, in_place=in_place)


def intersection_of_cubes(cube, other_cube):
    """
    Return the two Cubes of intersection given two Cubes.

    .. note:: The intersection of cubes function will ignore all single valued
        coordinates in checking the intersection.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * other_cube:
        An instance of :class:`iris.cube.Cube`.

    Returns:
        A pair of :class:`iris.cube.Cube` instances in a tuple corresponding
        to the original cubes restricted to their intersection.

    """
    # Take references of the original cubes (which will be copied when slicing later)
    new_cube_self = cube
    new_cube_other = other_cube

    # This routine has not been written to cope with multi-dimensional coordinates.
    for coord in cube.coords() + other_cube.coords():
        if coord.ndim != 1:
            raise iris.exceptions.CoordinateMultiDimError(coord)

    coord_comp = iris.analysis.coord_comparison(cube, other_cube)

    if coord_comp['ungroupable_and_dimensioned']:
        raise ValueError('Cubes do not share all coordinates in common, cannot intersect.')

    # cubes must have matching coordinates
    for coord in cube.coords():
        other_coord = other_cube.coord(coord=coord)

        # Only intersect coordinates which are different, single values coordinates may differ.
        if coord.shape[0] > 1 and coord != other_coord:
            intersected_coord = coord.intersect(other_coord)
            new_cube_self = new_cube_self.subset(intersected_coord)
            new_cube_other = new_cube_other.subset(intersected_coord)

    return new_cube_self, new_cube_other


def _assert_is_cube(cube):
    if not isinstance(cube, iris.cube.Cube):
        raise TypeError('The "cube" argument must be an instance of '
                        'iris.cube.Cube.')


def _assert_compatible(cube, other):
    """Checks to see if cube.data and another array can be broadcast to the same shape using ``numpy.broadcast_arrays``."""
    # This code previously returned broadcasted versions of the cube data and the other array.
    # As numpy.broadcast_arrays does not work with masked arrays (it returns them as ndarrays) operations
    # involving masked arrays would be broken.

    try:
        data_view, other_view = np.broadcast_arrays(cube.data, other)
    except ValueError, err:
        # re-raise
        raise ValueError("The array was not broadcastable to the cube's data shape. The error message from numpy when broadcasting:\n%s\n"
                         "The cube's shape was %s and the array's shape was %s" % (err, cube.shape, other.shape))

    if cube.shape != data_view.shape:
        raise ValueError("The array operation would increase the dimensionality of the cube. The new cubes data would "
                         "have had to become: %s" % (data_view.shape, ))


def _assert_matching_units(cube, other, operation_noun):
    """
    Check that the units of the cube and the other item are the same, or if
    the other does not have a unit, skip this test
    """
    if cube.units != getattr(other, 'units', cube.units):
        raise iris.exceptions.NotYetImplementedError(
            'Differing units (%s & %s) %s not implemented' %
            (cube.units, other.units, operation_noun))


def add(cube, other, dim=None, ignore=True, in_place=False):
    """
    Calculate the sum of two cubes, or the sum of a cube and a coordinate or scalar
    value.

    When summing two cubes, they must both have the same coordinate systems & data resolution.

    When adding a coordinate to a cube, they must both share the same number of elements
    along a shared axis.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * other:
        An instance of :class:`iris.cube.Cube` or :class:`iris.coords.Coord`,
        or a number or :class:`numpy.ndarray`.

    Kwargs:

    * dim:
        If supplying a coord with no match on the cube, you must supply the dimension to process.
    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _add_subtract_common(np.add, 'addition', 'added', cube, other,
                                dim=dim, ignore=ignore, in_place=in_place)


def subtract(cube, other, dim=None, ignore=True, in_place=False):
    """
    Calculate the difference between two cubes, or the difference between
    a cube and a coordinate or scalar value.

    When subtracting two cubes, they must both have the same coordinate systems & data resolution.

    When subtracting a coordinate to a cube, they must both share the same number of elements
    along a shared axis.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * other:
        An instance of :class:`iris.cube.Cube` or :class:`iris.coords.Coord`,
        or a number or :class:`numpy.ndarray`.

    Kwargs:

    * dim:
        If supplying a coord with no match on the cube, you must supply the dimension to process.
    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _add_subtract_common(np.subtract, 'subtraction', 'subtracted', cube,
                                other, dim=dim, ignore=ignore,
                                in_place=in_place)


def _add_subtract_common(operation_function, operation_noun,
                         operation_past_tense, cube, other, dim=None,
                         ignore=True, in_place=False):
    """
    Function which shares common code between addition and subtraction of cubes.

    operation_function   - function which does the operation (e.g. numpy.subtract)
    operation_symbol     - the textual symbol of the operation (e.g. '-')
    operation_noun       - the noun of the operation (e.g. 'subtraction')
    operation_past_tense - the past tense of the operation (e.g. 'subtracted')

    """
    _assert_is_cube(cube)
    _assert_matching_units(cube, other, operation_noun)

    if isinstance(other, iris.cube.Cube):
        # get a coordinate comparison of this cube and the cube to do the
        # operation with
        coord_comp = iris.analysis.coord_comparison(cube, other)

        if coord_comp['transposable']:
            # User does not need to transpose their cubes if numpy
            # array broadcasting will make the dimensions match
            broadcast_padding = cube.ndim - other.ndim
            coord_dims_equal = True
            for coord_group in coord_comp['transposable']:
                cube_coord, other_coord = coord_group.coords
                cube_coord_dims = cube.coord_dims(coord=cube_coord)
                other_coord_dims = other.coord_dims(coord=other_coord)
                other_coord_dims_broadcasted = tuple(
                    [dim + broadcast_padding for dim in other_coord_dims])
                if cube_coord_dims != other_coord_dims_broadcasted:
                    coord_dims_equal = False

            if not coord_dims_equal:
                raise ValueError('Cubes cannot be %s, differing axes. '
                                 'cube.transpose() may be required to '
                                 're-order the axes.' % operation_past_tense)

        # provide a deprecation warning if the ignore keyword has been set
        if ignore is not True:
            warnings.warn('The "ignore" keyword has been deprecated in '
                          'add/subtract. This functionality is now automatic. '
                          'The provided value to "ignore" has been ignored, '
                          'and has been automatically calculated.')

        bad_coord_grps = (coord_comp['ungroupable_and_dimensioned']
                          + coord_comp['resamplable'])
        if bad_coord_grps:
            raise ValueError('This operation cannot be performed as there are '
                             'differing coordinates (%s) remaining '
                             'which cannot be ignored.'
                             % ', '.join({coord_grp.name() for coord_grp
                                          in bad_coord_grps}))
    else:
        coord_comp = None

    new_cube = _binary_op_common(operation_function, operation_noun, cube,
                                 other, cube.units, dim, in_place)

    if coord_comp:
        # If a coordinate is to be ignored - remove it
        ignore = filter(None, [coord_grp[0] for coord_grp
                               in coord_comp['ignorable']])
        for coord in ignore:
            new_cube.remove_coord(coord)

    return new_cube


def multiply(cube, other, dim=None, in_place=False):
    """
    Calculate the product of a cube and another cube or coordinate.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * other:
        An instance of :class:`iris.cube.Cube` or :class:`iris.coords.Coord`,
        or a number or :class:`numpy.ndarray`.

    Kwargs:

    * dim:
        If supplying a coord with no match on the cube, you must supply the dimension to process.

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    _assert_is_cube(cube)
    other_unit = getattr(other, 'units', '1')
    new_unit = cube.units * other_unit
    return _binary_op_common(np.multiply, 'multiplication', cube, other,
                             new_unit, dim, in_place)


def divide(cube, other, dim=None, in_place=False):
    """
    Calculate the division of a cube by a cube or coordinate.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * other:
        An instance of :class:`iris.cube.Cube` or :class:`iris.coords.Coord`,
        or a number or :class:`numpy.ndarray`.

    Kwargs:

    * dim:
        If supplying a coord with no match on the cube, you must supply the dimension to process.

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    _assert_is_cube(cube)
    other_unit = getattr(other, 'units', '1')
    new_unit = cube.units / other_unit
    return _binary_op_common(np.divide, 'divison', cube, other, new_unit, dim,
                             in_place)


def exponentiate(cube, exponent, in_place=False):
    """
    Returns the result of the given cube to the power of a scalar.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.
    * exponent:
        The integer or floating point exponent.

        .. note:: When applied to the cube's unit, the exponent must result in a unit
            that can be described using only integer powers of the basic units.

            e.g. Unit('meter^-2 kilogram second^-1')

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    _assert_is_cube(cube)
    def power(data, out=None):
        return np.power(data, exponent, out)
    return _math_op_common(cube, power, cube.units ** exponent,
                           in_place=in_place)


def exp(cube, in_place=False):
    """
    Calculate the exponential (exp(x)) of the cube.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.

    .. note::

        Taking an exponential will return a cube with dimensionless units.

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _math_op_common(cube, np.exp, iris.unit.Unit('1'),
                           in_place=in_place)


def log(cube, in_place=False):
    """
    Calculate the natural logarithm (base-e logarithm) of the cube.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _math_op_common(cube, np.log, cube.units.log(math.e),
                           in_place=in_place)


def log2(cube, in_place=False):
    """
    Calculate the base-2 logarithm of the cube.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _math_op_common(cube, np.log2, cube.units.log(2),
                           in_place=in_place)


def log10(cube, in_place=False):
    """
    Calculate the base-10 logarithm of the cube.

    Args:

    * cube:
        An instance of :class:`iris.cube.Cube`.

    Kwargs:

    * in_place:
        Whether to create a new Cube, or alter the given "cube".

    Returns:
        An instance of :class:`iris.cube.Cube`.

    """
    return _math_op_common(cube, np.log10, cube.units.log(10),
                           in_place=in_place)


def _binary_op_common(operation_function, operation_noun, cube, other,
                      new_unit, dim=None, in_place=False):
    """
    Function which shares common code between binary operations.

    operation_function   - function which does the operation (e.g. numpy.divide)
    operation_noun       - the noun of the operation (e.g. 'division')
    cube                 - the cube whose data is used as the first argument
                           to `operation_function`
    other                - the cube, coord, ndarray or number whose data is
                           used as the second argument
    new_unit             - unit for the resulting quantity
    dim                  - dimension along which to apply `other` if it's a
                           coordinate that is not found in `cube`
    in_place             - whether or not to apply the operation in place to
                           `cube` and `cube.data`
    """
    _assert_is_cube(cube)

    if isinstance(other, iris.coords.Coord):
        other = _broadcast_cube_coord_data(cube, other, operation_noun, dim)
    elif isinstance(other, iris.cube.Cube):
        # TODO: add intelligent broadcasting along coordinate dimensions for
        # all binary operators, not just + and -
        other = other.data
    # don't worry about checking for other data types (such as scalers or
    # np.ndarrays) because _assert_compatible validates that they are broadcast
    # compatible with cube.data
    _assert_compatible(cube, other)

    def unary_func(x, out=None):
        ret = operation_function(x, other, out)
        if ret is NotImplemented:
            # explicitly raise the TypeError, so it gets raised even if, for
            # example, `iris.analysis.maths.multiply(cube, other)` is called
            # directly instead of `cube * other`
            raise TypeError('cannot %s %r and %r objects' %
                            (operation_function.__name__, type(x).__name__,
                             type(other).__name__))
        return ret
    return _math_op_common(cube, unary_func, new_unit, in_place)


def _broadcast_cube_coord_data(cube, other, operation_noun, dim=None):
    # What dimension are we processing?
    data_dimension = None
    if dim is not None:
        # Ensure the given dim matches the coord
        if other in cube.coords() and cube.coord_dims(other) != [dim]:
            raise ValueError("dim provided does not match dim found for coord")
        data_dimension = dim
    else:
        # Try and get a coord dim
        if other.shape != (1,):
            try:
                coord_dims = cube.coord_dims(other)
                data_dimension = coord_dims[0] if coord_dims else None
            except iris.exceptions.CoordinateNotFoundError:
                raise ValueError("Could not determine dimension for %s. "
                                 "Use %s(cube, coord, dim=dim)"
                                 % (operation_noun, operation_noun))

    if other.ndim != 1:
        raise iris.exceptions.CoordinateMultiDimError(other)

    if other.has_bounds():
        warnings.warn('%s by a bounded coordinate not well defined, ignoring '
                      'bounds.' % operation_noun)

    points = other.points

    # If the `data_dimension` is defined then shape the provided points for
    # proper array broadcasting
    if data_dimension is not None:
        points_shape = [1] * cube.ndim
        points_shape[data_dimension] = -1
        points = points.reshape(points_shape)

    return points


def _math_op_common(cube, operation_function, new_unit, in_place=False):
    _assert_is_cube(cube)
    if in_place:
        new_cube = cube
        operation_function(new_cube.data, out=new_cube.data)
    else:
        new_cube = cube.copy(data=operation_function(cube.data))
    iris.analysis.clear_phenomenon_identity(new_cube)
    new_cube.units = new_unit
    return new_cube

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
from itertools import izip, product

import numpy as np
from numpy.lib.stride_tricks import as_strided

from iris.analysis.interpolate import (_extend_circular_coord_and_data,
                                       _extend_circular_data)
from iris.coords import Coord, DimCoord, AuxCoord
import iris.cube
from iris.experimental.regrid import _RegularGridInterpolator
from iris.util import remap_cube_dimensions


_DEFAULT_DTYPE = np.float16
_MODE_LINEAR = 'linear'
_MODE_ERROR = 'error'
_MODE_NAN = 'nan'


class LinearInterpolator(object):
    """


    """
    def __init__(self, cube, coords, extrapolation_mode=_MODE_LINEAR):
        """

        Args:

        * cube:

        * coords:

        Kwargs:

        * extrapolation_mode:


        """
        self.cube = cube
        self.coords = []
        self.mode = extrapolation_mode

        self._coord_points = []
        self._coord_dims = []
        self._circulars = []
        self._interpolator = None

        self._setup(coords)

    def _account_for_circular(self, coord_points, data):
        """
        Extend the given data array, and re-centralise coordinate points
        for circular (1D) coordinates.

        """
        # Map all the requested values into the range of the source
        # data (centred over the centre of the source data to allow
        # extrapolation where required).
        if self._circulars:
            for _, data_dim, _, _, _ in self._circulars:
                data = _extend_circular_data(data, data_dim)

            for (interp_dim, _, src_min, src_max,
                    src_modulus) in self._circulars:
                offset = (src_max + src_min - src_modulus) * 0.5
                coord_points[:, interp_dim] -= offset
                coord_points[:, interp_dim] = (coord_points[:, interp_dim] %
                                               src_modulus) + offset
        return coord_points, data

    def _interpolate(self, data, coord_points):
        """
        Create and cache the underlying interpolator instance
        before invoking it to perform interpolation over the provided
        data with the coordinate point values.

        """
        dtype = self.interpolated_dtype(data.dtype)
        if data.dtype != dtype:
            # Perform dtype promotion.
            data = data.astype(dtype)

        if self._interpolator is None:
            # Cache the interpolator instance.
            self._interpolator = _RegularGridInterpolator(self._coord_points,
                                                          data,
                                                          bounds_error=False,
                                                          fill_value=None)
        self._interpolator_mode()
        self._interpolator.values = data
        return self._interpolator(coord_points)

    def _interpolator_mode(self):
        """
        Configure the underlying interpolator given the mode of
        extrapolation.

        """
        if self.mode is None:
            self.mode = _MODE_LINEAR

        if self.mode == _MODE_LINEAR:
            self._interpolator.bounds_error = False
            self._interpolator.fill_value = None
        elif self.mode == _MODE_ERROR:
            self._interpolator.bounds_error = True
        elif self.mode == _MODE_NAN:
            self._interpolator.bounds_error = False
            self._interpolator.fill_value = np.nan
        else:
            msg = 'Extrapolation mode {!r} not supported.'
            raise ValueError(msg.format(self.mode))

    def _prepare_points(self, sample_points):
        interp_points = [None] * len(self._coord_dims)
        # The shape of the data after full interpolation of each dimension.
        interpolated_shape = list(self.cube.shape)
        cube_dim_to_result_dim = {}

        interp_coords_seen = 0
        total_non_interp_coords = self.cube.ndim - len(self._coord_dims)

        for dim in range(self.cube.ndim):
            if dim not in self._coord_dims:
                cube_dim_to_result_dim[dim] = dim - interp_coords_seen
            else:
                cube_dim_to_result_dim[dim] = total_non_interp_coords + \
                    self._coord_dims.index(dim)
                interp_coords_seen += 1

        result_to_cube_dim = {v: k for k, v in cube_dim_to_result_dim.items()}
        coord_names = [coord.name() for coord in self.coords]

        for sample_coord, points in sample_points:
            if isinstance(sample_coord, Coord):
                sample_coord = sample_coord.name()
            if sample_coord not in coord_names:
                raise ValueError('Coordinate {!r} was not one of those passed '
                                 'to the constructor.'.format(sample_coord))
            coord = self.cube.coord(sample_coord)
            points = np.array(points, ndmin=1)
            dtype = self.interpolated_dtype(points.dtype)
            points = list(np.asanyarray(points, dtype=dtype))

            index = coord_names.index(sample_coord)
            interp_points[index] = points
            interpolated_shape[self._coord_dims[index]] = len(points)

        # Given an expected shape for the final array, compute the shape of
        # the array which puts the interpolated dimensions last.
        new_dimension_order = (lambda (dim, length):
                               cube_dim_to_result_dim[dim])
        _, target_shape = zip(*sorted(enumerate(interpolated_shape),
                                      key=new_dimension_order))

        # Now compute the transpose array which needs to be applied to the
        # previously computed shape to get back to the expected interpolated
        # shape.
        old_dimension_order = lambda (dim, length): result_to_cube_dim[dim]
        transpose_order, _ = zip(*sorted(enumerate(interpolated_shape),
                                         key=old_dimension_order))

        # Convert the interpolation points into a cross-product array
        # with shape (n_interp_points, n_dims)
        interp_points = np.asarray([i for i in product(*interp_points)])

        return interp_points, target_shape, transpose_order

    def _resample_coord(self, sample_points, coord, coord_dims):
        data = self.points(sample_points, coord.points, coord_dims)
        index = tuple(0 if dim not in coord_dims else slice(None)
                      for dim in range(self.cube.ndim))
        new_points = data[index]
        # Watch out for DimCoord instances that are no longer monotonic
        # after the resampling.
        try:
            new_coord = coord.copy(new_points)
        except ValueError:
            aux_coord = AuxCoord.from_coord(coord)
            new_coord = aux_coord.copy(new_points)
        return new_coord

    def _setup(self, coords):
        """
        Perform initial start-up configuration and validation based on the
        cube and the specified coordinates to be interpolated over.

        """
        cube = self.cube
        # Triggers the loading - is that really necessary...?
        data = self.cube.data
        self.coords = [self.cube.coord(coord) for coord in coords]

        coord_points_list = []

        for interp_dim, coord in enumerate(self.coords):
            coord_dims = cube.coord_dims(coord)

            if getattr(coord, 'circular', False):
                # Only DimCoords can be circular.
                coord_points, data = _extend_circular_coord_and_data(
                    coord, data, coord_dims[0])
                modulus = getattr(coord.units, 'modulus', 0)
                self._circulars.append((interp_dim, coord_dims[0],
                                        coord_points.min(),
                                        coord_points.max(), modulus))
            else:
                coord_points = coord.points
            coord_points_list.append([coord_points, coord_dims])

        self._coord_points, coord_dims_lists = zip(*coord_points_list)

        for dim_list in coord_dims_lists:
            for dim in dim_list:
                if dim not in self._coord_dims:
                    self._coord_dims.append(dim)

        self._validate()

    def _validate(self):
        """
        Perform all sanity checks to ensure that the interpolation request
        over the cube with the specified coordinates is valid and can be
        performed.

        """
        if len(set(self._coord_dims)) != len(self.coords):
            raise ValueError('Coordinates repeat a data dimension - the '
                             'interpolation would be over-specified.')

        for coord in self.coords:
            if coord.ndim != 1:
                raise ValueError('Interpolation coords must be 1-d for '
                                 'rectilinear interpolation.')

            if not isinstance(coord, DimCoord):
                # Check monotonic.
                if not iris.util.monotonic(coord.points, strict=True):
                    msg = 'Cannot interpolate over the non-' \
                        'monotonic coordinate {}.'
                    raise ValueError(msg.format(coord.name()))

        # Force all coordinates to be monotonically increasing. Generally this
        # isn't always necessary for a rectilinear interpolator, but it is a
        # common requirement.
        self.coord_decreasing = [np.all(np.diff(points[:2]) < 0)
                                 for points in self._coord_points]
        if np.any(self.coord_decreasing):
            pairs = izip(self.coord_decreasing, self._coord_points)
            self._coord_points = [points[::-1] if is_decreasing else points
                                  for is_decreasing, points in pairs]

    def interpolated_dtype(self, dtype):
        """
        Determine the base dtype required by the underlying interpolator.

        Args:

        * dtype:
            The :class:`~numpy.dtype` of the candidate data to be interpolated.

        Returns:
            The :class:`~numpy.dtype` that the underlying interpolator
            is expecting at a minimum for its input data.

        .. note::

            The interpolator will promote the :class:`~numpy.dtype` to
            a floating point representation.

        """
        # Default to float.
        return np.result_type(_DEFAULT_DTYPE, dtype)

    def points(self, sample_points, data, data_dims=None):
        """
        Interpolate the given data values at the specified list of orthogonal
        (coord, points) pairs.

        Args:

        * sample_points:
            Sequence of (coord, points) pairs. Order of coordinates needn't be
            the same order of the coordinates passed to this interpolator's
            constructor.
        * data:
            The data to interpolate - not necessarily the data from the cube
            that was used to construct this interpolator. If the data has
            fewer dimensions, then data_dims must be defined.

        Kwargs:

        * data_dims:
            The dimensions of the given data array in terms of the original
            cube passed through to this interpolator's constructor. If None,
            the data dimensions must map one-to-one onto the increasing
            dimension order of the cube.

        Returns:
            An :class:`~numpy.ndarray` or :class:`~numpy.ma.MaskedArray`
            instance of the interpolated data.

        .. note::

            The implementation of this method means that even for small
            subsets of the original cube's data, the data to be interpolated
            will be broadcast into the orginal cube's shape - thus resulting
            in more interpolation calls than are optimally needed. This has
            been done for implementation simplification, but there is no
            fundamental reason this must be the case.

        """
        points, final_shape, final_order = self._prepare_points(sample_points)
        data_dims = data_dims or range(self.cube.ndim)

        if len(data_dims) != data.ndim:
            msg = 'Data being interpolated is not consistent with ' \
                'the data passed through.'
            raise ValueError(msg)

        if sorted(data_dims) != list(data_dims):
            # To do this, a pre & post transpose will be necessary.
            msg = 'Currently only increasing data_dims is supported.'
            raise NotImplementedError(msg)

        if data_dims != range(self.cube.ndim):
            # Broadcast the data into the shape of the original cube.
            strides = list(data.strides)
            for dim in range(self.cube.ndim):
                if dim not in data_dims:
                    strides.insert(dim, 0)

            data = as_strided(data, strides=strides, shape=self.cube.shape)

        points, data = self._account_for_circular(points, data)

        # Build up a shape suitable for passing to ndindex, inside the loop we
        # will insert slice(None) on the data indices.
        iter_shape = list(length if dim not in self._coord_dims
                          else 1
                          for dim, length in enumerate(data.shape))

        result_shape = list(length for index, length in enumerate(data.shape)
                            if index not in self._coord_dims)
        # Keep track of the dimension to put the interpolated data into.
        interpolate_dimension = len(result_shape)
        result_shape.insert(interpolate_dimension, points.shape[0])

        masked = isinstance(data, np.ma.MaskedArray)
        dtype = self.interpolated_dtype(data.dtype)
        if masked:
            result_data = np.ma.empty(result_shape, dtype=dtype)
            if not isinstance(data.mask, np.ma.MaskType):
                result_data.mask = np.zeros(result_shape, dtype=np.bool)
        else:
            result_data = np.empty(result_shape, dtype=dtype)

        # Iterate through each slice of the data, updating the interpolator
        # with the new data as we go.
        for ndindex in np.ndindex(tuple(iter_shape)):
            interpolant_index = [position
                                 for dim, position in enumerate(ndindex)
                                 if dim not in self._coord_dims]
            interpolant_index.insert(interpolate_dimension, slice(None))
            index = tuple(position if dim not in self._coord_dims
                          else slice(None)
                          for dim, position in enumerate(ndindex))
            sub_data = data[index]

            order, _ = zip(*sorted(enumerate(self._coord_dims),
                                   key=lambda (i, dim): dim))
            sub_data = np.transpose(sub_data, order).copy()

            interpolated_data = self._interpolate(sub_data, points)
            result_data[interpolant_index] = interpolated_data
            if masked and not isinstance(data.mask, np.ma.MaskType):
                interpolated_data = self._interpolate(sub_data.mask, points)
                result_data.mask[interpolant_index] = interpolated_data > 0

        # Turn the interpolated data back into the order that
        # it was given to us in the first place.
        return np.transpose(result_data.reshape(final_shape), final_order)

    def __call__(self, sample_points, collapse_scalar=True):
        """
        Construct a cube from the specified orthogonal interpolation points.

        Args:

        * sample_points:
            A sequence of (coord, points) pairs over which to interpolate.
            The order of coordinates needn't be the same order of the
            coordinates passed to this interpolator's constructor.

        Kwargs:

        * collapse_scalar:
            Whether to collapse the dimension of the scalar sample points
            in the resulting cube. Default is True.

        Returns:
            A cube interpolated at the given sample points. The dimensionality
            of the cube will be the number of original cube dimensions minus
            the number of scalar coordinates, if collapse_scalar is True.

        """
        data = self.cube.data
        # Interpolate the cube payload.
        interpolated_data = self.points(sample_points, data)
        # Ensure the data payload has non-zero dimensionality.
        if interpolated_data.ndim == 0:
            interpolated_data = np.asanyarray(interpolated_data, ndmin=1)

        # Get hold of the original interpolation coordinates in terms of the
        # given cube.
        interp_coords = self.coords
        sample_point_order = [self.cube.coord(coord)
                              for coord, _ in sample_points]

        # Keep track of the dimensions for which sample points is scalar when
        # collapse_scalar is True - we will remove these scalar dimensions
        # later on.
        _new_scalar_dims = []
        if collapse_scalar:
            for coord, points in sample_points:
                coord = self.cube.coord(coord)
                if np.array(points).ndim == 0:
                    new_dim = self._coord_dims[interp_coords.index(coord)]
                    _new_scalar_dims.append(new_dim)

        cube = self.cube
        new_cube = iris.cube.Cube(interpolated_data)
        new_cube.metadata = cube.metadata

        def construct_new_coord_given_points(coord, points):
            # Handle what was previously a DimCoord which may no longer be
            # monotonic.
            try:
                return DimCoord.from_coord(coord).copy(points)
            except ValueError:
                return AuxCoord.from_coord(coord).copy(points)

        # Keep track of id(coord) -> new_coord for aux factory construction
        # later on.
        coord_mapping = {}

        dims_with_dim_coords = []

        # Copy/interpolate the coordinates.
        for coord in cube.dim_coords:
            dim, = cube.coord_dims(coord)
            if coord in interp_coords:
                new_points = sample_points[sample_point_order.index(coord)][1]
                new_coord = construct_new_coord_given_points(coord, new_points)
            elif set([dim]).intersection(set(self._coord_dims)):
                # Interpolate the coordinate payload.
                new_coord = self._resample_coord(sample_points, coord, [dim])
            else:
                new_coord = coord.copy()

            # new_coord may no longer be a dim coord, so check we don't need
            # to add it as an aux coord (thus leaving the dim anonymous).
            if isinstance(new_coord, DimCoord) and dim is not None:
                new_cube._add_unique_dim_coord(new_coord, dim)
                dims_with_dim_coords.append(dim)
            else:
                new_cube._add_unique_aux_coord(new_coord, dim)
            coord_mapping[id(coord)] = new_coord

        for coord in cube.aux_coords:
            dims = cube.coord_dims(coord)
            if coord in interp_coords:
                new_points = sample_points[sample_point_order.index(coord)][1]
                new_coord = construct_new_coord_given_points(coord, new_points)
                dims = [self._coord_dims[interp_coords.index(coord)]]
            elif set(dims).intersection(set(self._coord_dims)):
                # Interpolate the coordinate payload.
                new_coord = self._resample_coord(sample_points, coord, dims)
            else:
                new_coord = coord.copy()

            new_dims = dims

            if (isinstance(new_coord, DimCoord) and len(new_dims) > 0
                    and new_dims[0] not in dims_with_dim_coords):
                new_cube._add_unique_dim_coord(new_coord, new_dims)
                dims_with_dim_coords.append(new_dims[0])
            else:
                new_cube._add_unique_aux_coord(new_coord, new_dims)
            coord_mapping[id(coord)] = new_coord

        for factory in self.cube.aux_factories:
            new_cube.add_aux_factory(factory.updated(coord_mapping))

        if _new_scalar_dims:
            remap_cube_dimensions(new_cube, remove_axes=_new_scalar_dims)

        return new_cube

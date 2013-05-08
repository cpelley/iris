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
High-level plotting extensions to :mod:`iris.plot`.

These routines work much like their :mod:`iris.plot` counterparts, but they
automatically add a plot title, axis titles, and a colour bar when appropriate.

See also: :ref:`matplotlib <matplotlib:users-guide-index>`.

"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import iris.config
import iris.coords
import iris.plot as iplt
import cartopy.crs as ccrs


def _title(cube_or_coord, with_units):
    if cube_or_coord is None:
        title = ''
    else:
        title = cube_or_coord.name().replace('_', ' ').capitalize()
        units = cube_or_coord.units
        if with_units and not (units.is_unknown() or
                               units.is_no_unit() or
                               units == iris.unit.Unit('1')):

            # For non-time units use the shortest unit representation e.g.
            # prefer 'K' over 'kelvin', but not '0.0174532925199433 rad'
            # over 'degrees'
            if (not units.is_time() and not units.is_time_reference() and
                len(units.symbol) < len(str(units))):
                units = units.symbol
            title += ' / {}'.format(units)

    return title


def _label(cube, mode, result=None, ndims=2, coords=None):
    """Puts labels on the current plot using the given cube."""
    
    plt.title(_title(cube, with_units=False))

    if result is not None:
        draw_edges = mode == iris.coords.POINT_MODE
        bar = plt.colorbar(result, orientation='horizontal',
                           drawedges=draw_edges)
        has_known_units = not (cube.units.is_unknown() or cube.units.is_no_unit())
        if has_known_units and cube.units != iris.unit.Unit('1'):
            # Use shortest unit representation for anything other than time
            if (not cube.units.is_time() and not cube.units.is_time_reference() and
                len(cube.units.symbol) < len(str(cube.units))):
                bar.set_label(cube.units.symbol)
            else:
                bar.set_label(cube.units)
        # Remove the tick which is put on the colorbar by default.
        bar.ax.tick_params(length=0)
    
    if coords is None:
        plot_defn = iplt._get_plot_defn(cube, mode, ndims)
    else:
        plot_defn = iplt._get_plot_defn_custom_coords_picked(cube, coords, mode, ndims=ndims)
    
    if ndims == 2:
        if not iplt._can_draw_map(plot_defn.coords):
            plt.ylabel(_title(plot_defn.coords[0], with_units=True))
            plt.xlabel(_title(plot_defn.coords[1], with_units=True))
    elif ndims == 1:
        plt.xlabel(_title(plot_defn.coords[0], with_units=True))
        plt.ylabel(_title(cube, with_units=True))
    else:
        raise ValueError('Unexpected number of dimensions (%s) given to _label.' % ndims)


def _label_with_bounds(cube, result=None, ndims=2, coords=None):
    _label(cube, iris.coords.BOUND_MODE, result, ndims, coords)


def _label_with_points(cube, result=None, ndims=2, coords=None):
    _label(cube, iris.coords.POINT_MODE, result, ndims, coords)


def animate(cube, coords, plot_func=plt.contourf, **kwargs):
    """
    Animates the given cube across a given coordinate.

    Args:

    * cube (:class:`iris.cube.Cube`):
        A :class:`iris.cube.Cube`, to be animated.

    Kwargs:

    * coords: list of :class:`~iris.coords.Coord` objects or coordinate names
        Use the given coordinates as the axes for the plot. The order of the
        given coordinates indicates which axis to use for each, where the
        first element is the horizontal axis of the plot and the second
        element is the vertical axis of the plot.

    * plot_func: plotting function to animate.

    See :func:`matplotlib.animation.FuncAnimation` for details of other valid
    keyword arguments.

    .. note::

        If no plotting function is supplied, the :class:`matplotlib.pyplot`
        plot function is used.  The raw data will be plotted without any
        coordinate system in this case.

    """
    kwargs.setdefault('interval', 100)

    def update_animation_proj(i, anim_cube, ax):
        plt.gca().cla()
        im = plot_func(anim_cube[i])
        return im,

    def update_animation(i, anim_cube, ax):
        plt.gca().cla()
        im = plot_func(anim_cube[i].data)
        return im,

    fig = plt.figure()

    # Slice cube according to the specified coordinates
    # Requires to be turned into a list to repeat anim.
    anim_cube = list(cube.slices(coords))
    frames = xrange(len(anim_cube))

    #if using matplotlib plot function, do not set-up a projection
    if plot_func.__module__ in ['iris.plot']: 
        update = update_animation_proj
    elif plot_func.__module__ in ['iris.quickplot']:
        raise LookupError(
            '{} not supported for animation'.format(plot_func.__module__))
    else:
        update = update_animation

    ani = animation.FuncAnimation(fig, update,
                                  frames=frames,
                                  fargs=(anim_cube, None),
                                  **kwargs
                                  )
    plt.show()


def contour(cube, *args, **kwargs):
    """
    Draws contour lines on a labelled plot based on the given Cube.

    With the basic call signature, contour "level" values are chosen automatically::

        contour(cube)

    Supply a number to use *N* automatically chosen levels::

        contour(cube, N)

    Supply a sequence *V* to use explicitly defined levels::
        
        contour(cube, V)
    
    See :func:`iris.plot.contour` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.contour(cube, *args, **kwargs)
    _label_with_points(cube, coords=coords)
    return result


def contourf(cube, *args, **kwargs):
    """
    Draws filled contours on a labelled plot based on the given Cube.
    
    With the basic call signature, contour "level" values are chosen automatically::

        contour(cube)

    Supply a number to use *N* automatically chosen levels::

        contour(cube, N)

    Supply a sequence *V* to use explicitly defined levels::
        
        contour(cube, V)
    
    See :func:`iris.plot.contourf` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.contourf(cube, *args, **kwargs)
    _label_with_points(cube, result, coords=coords)
    return result


def outline(cube, coords=None):
    """Draws cell outlines on a labelled plot based on the given Cube."""
    result = iplt.outline(cube, coords=coords)
    _label_with_bounds(cube, coords=coords)
    return result


def pcolor(cube, *args, **kwargs):
    """
    Draws a labelled pseudocolor plot based on the given Cube.
    
    See :func:`iris.plot.pcolor` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.pcolor(cube, *args, **kwargs)
    _label_with_bounds(cube, result, coords=coords)
    return result


def pcolormesh(cube, *args, **kwargs):
    """
    Draws a labelled pseudocolour plot based on the given Cube.
    
    See :func:`iris.plot.pcolormesh` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.pcolormesh(cube, *args, **kwargs)
    _label_with_bounds(cube, result, coords=coords)
    return result


def points(cube, *args, **kwargs):
    """
    Draws sample point positions on a labelled plot based on the given Cube.
    
    See :func:`iris.plot.points` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.points(cube, *args, **kwargs)
    _label_with_points(cube, coords=coords)
    return result


def plot(cube, *args, **kwargs):
    """
    Draws a labelled line plot based on the given Cube.
    
    See :func:`iris.plot.plot` for details of valid keyword arguments.
    
    """
    coords = kwargs.get('coords')
    result = iplt.plot(cube, *args, **kwargs)
    _label_with_points(cube, ndims=1, coords=coords)
    return result

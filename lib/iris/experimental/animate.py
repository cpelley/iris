# (C) British Crown Copyright 2013, Met Office
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
Wrapper for animating iris cubes using iris or matplotlib plotting functions

"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import iris
import iris.plot as iplt


def animate(cube, plot_func, fig=None, *args, **kwargs):
    """
    Animates the given cube across a given coordinate.

    Args:

    * cube (:class:`iris.cube.Cube`):
        A :class:`~iris.cube.Cube`, to be animated.

    * plot_func: :class:`~matplotlib.pyplot` or :class:`~iris.plot`
            plotting function
        Plotting function to animate.

    Kwargs:

    * coords: list of :class:`~iris.coords.Coord` objects or coordinate names
        Slice cube according to the specified coordinates and animate along
        chosen coordinates.

    * vmin,vmax:
        Color scaling, see :class:`matplotlib.colors.Normalize` for further
        details.  Default values are determined by the min-max animation data.

    See :func:`matplotlib.animation.FuncAnimation` for details of other valid
    keyword arguments.

    Returns:
        :class:`~matplotlib.animation.FuncAnimation` object suitable for
        saving and or plotting.

    """
    kwargs.setdefault('interval', 100)
    coords = kwargs.pop('coords', None)

    if fig is None:
        fig = plt.gcf()

    def update_animation_iris(i, anim_cube, ax, vmin, vmax, coords):
        # Update frame using iris plotting wrapper
        plt.gca().cla()
        im = plot_func(anim_cube[i], vmin=vmin, vmax=vmax, coords=coords)
        return im

    def update_animation(i, anim_cube, ax, vmin, vmax, coords):
        # Update frame indexed based
        plt.gca().cla()
        im = plot_func(anim_cube[i].data, vmin=vmin, vmax=vmax)
        return im

    suitable_plot = {'plot': 1, 'points': 1, 'contour': 2, 'contourf': 2,
                     'pcolor': 2, 'pcolormesh': 2}

    # Determine dimensionality of plot
    if plot_func.__name__ in suitable_plot:
        ndims = suitable_plot[plot_func.__name__]
    else:
        msg = ('specified plotting function not suitable: {}, '
               'those available: {}')
        msg = msg.format(plot_func.__name__, suitable_plot.keys())
        raise ValueError(msg)

    # Check cube dimensionality against plot for animating
    if cube.ndim-1 != ndims:
        msg = 'Cube must be %s-dimensional. Got %s dimensions.'
        raise ValueError(msg % (ndims+1, cube.data.ndim))

    # Set default coords
    mode = iris.coords.POINT_MODE
    if coords is not None:
        plot_defn = iplt._get_plot_defn_custom_coords_picked(
            cube[0], coords, mode, ndims=ndims)
    else:
        plot_defn = iplt._get_plot_defn(cube[0], mode, ndims=ndims)
    coords = list(reversed(plot_defn.coords))

    # Slice cube according to the specified coordinates
    # Requires to be turned into a list to repeat anim.
    anim_cube = list(cube.slices(coords))
    frames = xrange(len(anim_cube))

    # Determine plot range
    vmin = kwargs.pop('vmin', min([cc.data.min() for cc in anim_cube]))
    vmax = kwargs.pop('vmax', max([cc.data.max() for cc in anim_cube]))

    # If using matplotlib plot function, do not set-up a projection.
    if plot_func.__module__ in ['iris.plot', 'iris.quickplot']:
        update = update_animation_iris
    # Capture providing other plotting classes as input (eg. matplotlib.pyplot)
    else:
        update = update_animation

    ani = animation.FuncAnimation(fig, update,
                                  frames=frames,
                                  fargs=(anim_cube, None, vmin, vmax, coords),
                                  **kwargs
                                  )
    return ani

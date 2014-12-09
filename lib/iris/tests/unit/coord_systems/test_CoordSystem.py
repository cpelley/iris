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
"""Unit tests for the :class:`iris.coord_systems.CoordSystem` class."""

# Import iris.tests first so that some things can be initialised before
# importing anything else.
import iris.tests as tests

from iris.coord_systems import CoordSystem


class _TestCoordSystem(CoordSystem):
    # Subclass of abstract CoordSystem for testing purposes which removes
    # abstract methods.
    def as_cartopy_crs(self):
        pass

    def as_cartopy_projection(self):
        pass


class Test___eq__(tests.IrisTest):
    def test_different_class(self):
        # Ensure that we fail equality between coordinate systems of different
        # classes.
        class TestClass(_TestCoordSystem):
            pass

        self.assertFalse(_TestCoordSystem() == TestClass())

    def test_different_dict_keys(self):
        c1 = _TestCoordSystem()
        c1.attribute1 = 1
        c2 = _TestCoordSystem()
        c2.attribute2 = 1
        self.assertFalse(c1 == c2)

    def test_different_dict_values(self):
        c1 = _TestCoordSystem()
        c1.attribute = 1
        c2 = _TestCoordSystem()
        c2.attribute = 2
        self.assertFalse(c1 == c2)

    def test_equal(self):
        # Ensure that we return True for equality test between dictionary
        # values within default tolerance.
        c1 = _TestCoordSystem()
        c1.attribute = 1.
        c2 = _TestCoordSystem()
        c2.attribute = 1.0000000001
        self.assertTrue(c1 == c2)


if __name__ == '__main__':
    tests.main()

from abc import ABCMeta

from iris.analysis import interpolate


class Regridder(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def _setup(self):
        pass

    @abstractmethod
    def regrid(self):
        pass
    
    
class Interpolator(Regridder):
    @abstractmethod
    def point(self):
        pass

    @abstractmethod
    def value(self):
        pass


class Nearest(Interpolator):
    def index(self):



if __name__ == '__main__':
>>> intreg_object = Linear(src_cube, tgt_cube=None, weights=None, mask_tolerance=0.5) # Cache interpolation Setup (common to both regrid and interp)
>>> dir(intreg_object)
['value', 'points', 'regrid']
>>> intreg_object.points(sample_points)
>>> intreg_object.regrid(cube)
RuntimeError: 'tgt_cube' is not undefined

Linear(tgt_cube=other_cube, weights=weights, out=intreg_object) # Update cached regrid to existing cached interp.
>>> intreg_object.regrid(cube)


>>> intreg_object = Nearest(....)
>>> dir(intreg_object)
['index', 'value', 'points', 'regrid']

>>> intreg_object = PointInCell(....)
>>> dir(intreg_object)
['regrid']
from iris.analysis import interpolate


class Interpolator(object):
    def __init__(self, src_cube, mask_tolerance=0.0):
        """
        Setup for our interpolation and returns an object which holds any
        information that can persist between interpolation calls.
        
        """
        self._src_cube = src_cube

    def data_value(self, sample_point):
        """Returns a single interpolated data value.""" 
        return NotImplemented

    def points(self, sample_points):
        """
        Returns the interpolated points.
        """
        return NotImplemented


class Regridder(object):
    def __init__(self, source_cube, target_cube, weights=None, mask_tolerance=0.0):
        """
        Setup for our regridder and returns an object which holds any
        information that can persist between regrid calls.
        
        """
        self._src_grid = source_cube
        self._tgt_grid = target_cube
        self._weights = weights

    def regrid(self):
        return NotImplemented


class _NearestInterpolator(Interpolator):
    def index(self, sample_points):
        """
        Returns the indexes of the interpolated values.
        Applicable to extraction only (nearest neighbour).
        
        """
        return interpolate.nearest_neighbour_indices(self._src_cube,
                                                     sample_points)

    def data_value(self, sample_point):
        return interpolate.nearest_neighbour_data_value(self._src_cube,
                                                        sample_points)
    
    def points(self, sample_point):
        return interpolate.extract_nearest_neighbour(self._src_cube,
                                                     sample_points)


class _NearestRegridder(Regridder):
    def regrid(self, src_cube):
        return interpolate.regrid(self._src_grid, self_tgt_grid,
                                  weights=self._weights, mode='nearest')


class Nearest(object):
    @staticmethod
    def regridder(*args, **kwargs):
        return Regridder(*args, **kwargs)

    @staticmethod
    def interpolator(*args, **kwargs):
        return _NearestInterpolator(*args, **kwargs)
    

if __name__ == '__main__':
    import iris
    cube = iris.load_cube(iris.sample_data_path('air_temp.pp'))

    interpolator = Nearest.interpolator(cube)
    sample_points = [('latitude', 0), ('longitude', 10)]
    print interpolator.data_value(sample_points)

    regridder = Nearest.regridder(src_cube, target_cube)

    print 'completed'
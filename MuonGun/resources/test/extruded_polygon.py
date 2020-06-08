#!/usr/bin/env python

"""
Reference implementation of ExtrudedPolygon
"""

import numpy

from icecube.phys_services import Surface
from icecube.dataclasses import make_pair

def convex_hull(points):
    """Computes the convex hull of a set of 2D points.
 
    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    
    Lifted from http://code.icecube.wisc.edu/svn/sandbox/ckopper/eventinjector/python/util/__init__.py
    """
 
    # convert to a list of tuples
    points = [(p[0],p[1]) for p in points]
    
    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))
     
    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points
 
    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
 
    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
 
    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    hull = lower[:-1] + upper[:-1]
    
    # convert into numpy array
    return numpy.array(hull)

def hull_to_normals(points):
    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    vecs = points[1:]-points[:-1]
    magn = numpy.sqrt(vecs[:,0]**2 + vecs[:,1]**2)
    
    normals = numpy.array([vecs[:,1]/magn, -vecs[:,0]/magn, numpy.zeros(magn.shape)]).T
    
    return normals

def hull_to_lengths(points):
    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    vecs = points[1:]-points[:-1]
    
    return numpy.sqrt(vecs[:,0]**2 + vecs[:,1]**2)

def signed_area(points):
    """Returns the signed area of a given simple (i.e. non-intersecting) polygon.
    Positive if points are sorted counter-clockwise.
    """

    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    return numpy.sum(points[:,0][:-1]*points[:,1][1:] - points[:,0][1:]*points[:,1][:-1])/2.

class ExtrudedPolygon(Surface):
    """
    A convex polygon in the x-y plane, extruded in the z direction
    """
    def __init__(self, xy_points, z_range):
        """
        :param xy_points: a list of x-y coordinate pairs. The convex hull of
                          these points will form the edge of the surface in the
                          x-y plane
        :param z_range: a pair giving the lower and upper boundary of the
                        surface in z.
        """
        super(ExtrudedPolygon, self).__init__()
        assert len(xy_points) >= 3, "Need at least 3 points to form a closed polygon"
        
        hull = convex_hull(xy_points)
        
        # hull points, in counterclockwise order
        self._x = hull
        # next neighbor in the hull
        self._nx = numpy.roll(hull, -1, axis=0)
        # vector connecting each pair of points in the hull
        self._dx = self._nx - self._x
        
        self._z_range = z_range
        self.length = z_range[1] - z_range[0]
        self._side_lengths = hull_to_lengths(hull)
        
        side_normals = hull_to_normals(hull)
        side_areas = self._side_lengths*self.length
        cap_area = [signed_area(hull)]*2
        cap_normals = numpy.array([[0., 0., 1.], [0., 0., -1.]])
        
        self._areas = numpy.concatenate((side_areas, cap_area))
        self._normals = numpy.concatenate((side_normals, cap_normals))
        assert self._areas.size == self._normals.shape[0]
    
    def expand(self, padding):
        """
        Expand the x-y footprint by moving each edge out by a distance *padding*.
        """
        # A convex polygon can be offset by moving each vertex parallel to the
        # edges by a distance that is inversely proportional to the sine of the
        # counterclockwise angle between the edges that meet at each vertex.
        # This breaks down for edges that are [anti]parallel or, but neither
        # case should occur for maximally simplified polygons.
        
        # normalized vector connecting each vertex to the next one
        d = self._dx/self._side_lengths[:,None]
        # and the one connecting the previous vertex
        prev_d = numpy.roll(d, 1, axis=0)
        # sine of the inner angle of each vertex
        det = prev_d[:,0]*d[:,1] - prev_d[:,1]*d[:,0]
        assert (det != 0.).all(), "Edges can't be [anti]parallel"
        points = self._x + (padding/det[:,None])*(prev_d - d)
        
        z_range = [self._z_range[0]-padding, self._z_range[1]+padding]
        
        return type(self)(points, z_range)
    
    @classmethod
    def from_I3Geometry(cls, i3geo, padding=0):
        from collections import defaultdict
        strings = defaultdict(list)
        for omkey, omgeo in i3geo.omgeo:
            if omgeo.omtype != omgeo.IceTop:
                strings[omkey.string].append(list(omgeo.position))
        mean_xy = [numpy.mean(positions, axis=0)[0:2] for positions in strings.values()]
        zmax = max(max(p[2] for p in positions) for positions in strings.values())
        zmin = min(min(p[2] for p in positions) for positions in strings.values())
        
        self = cls(mean_xy, [zmin, zmax])
        if padding != 0:
            return self.expand(padding)
        else:
            return self
    
    @classmethod
    def from_file(cls, fname, padding=0):
        from icecube import icetray, dataio, dataclasses
        f = dataio.I3File(fname)
        fr = f.pop_frame(icetray.I3Frame.Geometry)
        f.close()
        return cls.from_I3Geometry(fr['I3Geometry'], padding)
    
    def area(self, dir):
        """
        Return projected area in the given direction
        
        :param dir: an I3Direction
        """
        # inner product with component normals
        inner = numpy.dot(self._normals, numpy.asarray((dir.x, dir.y, dir.z)))
        # only surfaces that face the requested direction count towards the area
        mask = inner < 0
        return -(inner*self._areas*mask).sum(axis=0)
    
    def partial_area(self, dir):
        inner = numpy.dot(self._normals, numpy.asarray((dir.x, dir.y, dir.z)))
        # only surfaces that face the requested direction count towards the area
        mask = inner < 0
        return -(inner*self._areas*mask)
    
    def azimuth_averaged_area(self, cos_theta):
        """
        Return projected area at the given zenith angle, averaged over all
        azimuth angles.
        
        :param cos_theta: cosine of the zenith angle
        """
        cap = self._areas[-1]
        sides = self._side_lengths.sum()*self.length/numpy.pi
        
        return cap*abs(cos_theta) + sides*numpy.sqrt(1-cos_theta**2)
    
    @staticmethod
    def _integrate_area(a, b, cap, sides):
        return numpy.pi*(cap*(b**2-a**2) + sides*(numpy.arccos(a) - numpy.arccos(b) - numpy.sqrt(1-a**2)*a + numpy.sqrt(1-b**2)*b))
    
    def entendue(self, cosMin=-1., cosMax=1.):
        """
        Integrate A * d\Omega over the given range of zenith angles
        
        :param cosMin: cosine of the maximum zenith angle
        :param cosMax: cosine of the minimum zenith angle
        :returns: a product of area and solid angle. Divide by
                  2*pi*(cosMax-cosMin) to obtain the average projected area in
                  this zenith angle range
        """
        
        # First, integrate over all azimuthal angles, exploiting the fact that
        # the projected area of a plane, averaged over a 2\pi rotation that
        # passes through the normal, is
        # A*\int_0^\pi \Theta(\sin\alpha)\sin\alpha d\alpha / 2\pi = A/\pi
        sides = self._side_lengths.sum()*self.length/numpy.pi
        # The projected area of the cap is independent of azimuth
        cap = self._areas[-1]
        
        if (cosMin >= 0 and cosMax >= 0):
            return self._integrate_area(cosMin, cosMax, cap, sides)
        elif (cosMin < 0 and cosMax <= 0):
            return self._integrate_area(-cosMax, -cosMin, cap, sides)
        elif (cosMin < 0 and cosMax > 0):
            return self._integrate_area(0, -cosMin, cap, sides) \
                + self._integrate_area(0, cosMax, cap, sides)
        else:
            raise ValueError("Can't deal with zenith range [%.1e, %.1e]" % (cosMin, cosMax))
        return numpy.nan
    
    def _point_in_hull(self, point):
        """
        Test whether point is inside the 2D hull by ray casting
        """
        x, y = point[0:2]
        # Find segments whose y range spans the current point
        mask = ((self._x[:,1] > y)&(self._nx[:,1] <= y))|((self._x[:,1] <= y)&(self._nx[:,1] > y))
        # Count crossings to the right of the current point
        xc = self._x[:,0] + (y-self._x[:,1])*self._dx[:,0]/self._dx[:,1]
        crossings = (x < xc[mask]).sum()
        inside = (crossings % 2) == 1
        
        return inside
    
    def _distance_to_hull(self, point, vec):
        """
        Calculate the most extreme displacements from x,y along dx,dy to points
        on the 2D hull
        """
        # calculate the distance along the ray to each line segment
        x, y = (self._x - point[:2]).T
        dx, dy = self._dx.T
        dirx, diry = vec[0:2]
        
        assert dirx+diry != 0, "Direction vector may not have zero length"
        
        # proportional distance along edge to intersection point
        # NB: if diry/dirx == dy/dx, the ray is parallel to the line segment
        nonparallel = diry*dx != dirx*dy
        alpha = numpy.where(nonparallel, (dirx*y - diry*x)/(diry*dx - dirx*dy), numpy.nan) 
        # check whether the intersection is actually in the segment
        mask = (alpha >= 0)&(alpha < 1)
        
        # distance along ray to intersection point
        if dirx != 0:
            beta = ((x + alpha*dx)/dirx)[mask]
        else:
            beta = ((y + alpha*dy)/diry)[mask]
        
        if beta.size == 0:
            return (numpy.nan,)*2
        else:
            return (numpy.nanmin(beta), numpy.nanmax(beta))
    
    def _distance_to_cap(self, point, dir, cap_z):
        return (cap_z-point[2])/dir[2]
    
    def _distance_to_caps(self, point, dir):
        return sorted((self._distance_to_cap(point, dir, cap_z) for cap_z in self._z_range))
    
    def GetIntersection(self, pos, dir):
        point = numpy.array((pos.x, pos.y, pos.z))
        vec = numpy.array((dir.x, dir.y, dir.z))
        
        no_intersection = make_pair(numpy.nan, numpy.nan)
        
        # perfectly vertical track: only check intersections with caps
        if abs(dir.z) == 1.:
            if not self._point_in_hull(point):
                return no_intersection
            else:
                return make_pair(*self._distance_to_caps(point, vec))
        # perfectly horizontal track: only check intersections with sides
        elif dir.z == 0.:
            if pos.z < self._z_range.first or pos.z > self._z_range.second:
                return no_intersection
            else:
                return make_pair(*self._distance_to_hull(point, vec))
        # general case: both rho and z components nonzero
        else:
            sides = numpy.array(self._distance_to_hull(point, vec))
            caps = self._distance_to_caps(point, vec)
            intersections = numpy.concatenate((sides, caps))
            
            if (caps[0] >= sides[1] or caps[1] <= sides[0]):
                return no_intersection
            else:
                return make_pair(numpy.max((sides[0], caps[0])), numpy.min((sides[1], caps[1])))

if __name__ == "__main__":
    
    from icecube.MuonGun import ExtrudedPolygon as CExtrudedPolygon
    from icecube.dataclasses import I3Position, I3Direction
    from icecube.phys_services import I3GSLRandomService, AxialCylinder
    
    numpy.random.seed(0)
    
    def make_surfaces(npoints=10, padding=0):
        npoints = 10
        x = numpy.random.uniform(0, 1e3, size=npoints)-5e2
        y = numpy.random.uniform(0, 1e3, size=npoints)-5e2
        z = [-500, 500]
        points = numpy.vstack((x,y)).T
    
        py_surface = ExtrudedPolygon(points, z).expand(padding)
        cpoints = []
        for z_ in z:
            cpoints += [I3Position(x_,y_,z_) for x_,y_ in zip(x,y)]
        cpp_surface = CExtrudedPolygon(cpoints, padding)
        
        return py_surface, cpp_surface
    
    def random_points(radius, npoints=10000):
        r = numpy.random.uniform(0, radius**3, size=npoints)**(1./3)
        azi = numpy.random.uniform(0, 2*numpy.pi, size=2*npoints)
        zen = numpy.arccos(numpy.random.uniform(-1, 1, size=2*npoints))
        return [(I3Position(r_, z1_, a1_, I3Position.sph), I3Direction(z2_, a2_)) for r_,z1_,a1_,z2_,a2_ in zip(r, zen[:npoints], azi[:npoints], zen[npoints:], azi[npoints:])]
    
    py_surface, cpp_surface = make_surfaces()
    numpy.testing.assert_almost_equal(py_surface._x, numpy.vstack((cpp_surface.x, cpp_surface.y)).T, err_msg="basic hull is not identical")
    numpy.testing.assert_almost_equal(py_surface._z_range, cpp_surface.z)
    
    py_surface, cpp_surface = make_surfaces(padding=50)
    numpy.testing.assert_almost_equal(py_surface._x, numpy.vstack((cpp_surface.x, cpp_surface.y)).T, err_msg="padded hull is not identical")
    numpy.testing.assert_almost_equal(py_surface._z_range, cpp_surface.z)
    
    # test intersection and projection code
    for i, (pos, dir) in enumerate(random_points(7e3)):
        i1 = py_surface.GetIntersection(pos, dir)
        i2 = cpp_surface.intersection(pos, dir)
        numpy.testing.assert_almost_equal(i1.first, i2.first)
        numpy.testing.assert_almost_equal(i1.second, i2.second)
        numpy.testing.assert_almost_equal(py_surface.area(dir), cpp_surface.area(dir))
    
    # test acceptance calculation
    ct = numpy.random.uniform(-1, 1, size=2*1000).reshape(2,1000)
    ct.sort(axis=0)
    for ct_lo, ct_hi in ct.T:
        py_acceptance = py_surface.entendue(ct_lo, ct_hi)
        cpp_acceptance = cpp_surface.acceptance(ct_lo, ct_hi)
        numpy.testing.assert_almost_equal(py_acceptance, cpp_acceptance)
    
    rng = I3GSLRandomService(1)
    # test area sampling
    pos, dir = cpp_surface.sample_impact_ray(rng)
    
    sampler = AxialCylinder(1000, 1000)
    volume = py_surface.length * py_surface.area(I3Direction(0,0))
    # compare volume and projected area to sampled versions
    for i, (pos, dir) in enumerate(random_points(0, 10)):
        nsamples = 50000
        hits = 0
        depth = 0.
        for _ in range(nsamples):
            pos = sampler.sample_impact_position(dir, rng)
            intersections = cpp_surface.intersection(pos, dir)
            if numpy.isfinite(intersections.first):
                hits += 1
                depth += (intersections.second - intersections.first)
        numpy.testing.assert_allclose(float(hits)/nsamples, py_surface.area(dir)/sampler.area(dir), rtol=2e-2, err_msg="Sampled area should match projected area")
        numpy.testing.assert_allclose(sampler.area(dir)*depth/nsamples, volume, rtol=2e-2, err_msg="Sampled volume should match exact volume")

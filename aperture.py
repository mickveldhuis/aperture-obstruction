import enum
import configparser
import numpy as np

# For coordinate transformations
from transformations import vec3, vec4, transform, rot_x, rot_z
from pytransform3d import transformations as pt

# CONSTANTS
config = configparser.ConfigParser()
config.read('config.ini')

# Telescope:
L_1 = config['mount'].getfloat('length_1') # distance floor-HA axis
L_2 = config['mount'].getfloat('length_2') # distance HA axis-Dec axis
L_3 = config['mount'].getfloat('length_3') # distance Dec axis-tube center
L_4 = config['guider'].getfloat('offset') # distance primary tube center-guider center
L_5 = config['finder'].getfloat('offset') # distance guider center-finder center

GUIDER_ANGLE = np.radians(config['guider'].getfloat('angle'))
FINDER_ANGLE =  np.radians(config['finder'].getfloat('angle'))

APERTURE_RADIUS = config['telescope'].getfloat('diameter')/2
APERTURE_SEC_RADIUS = config['telescope'].getfloat('sec_diameter')/2
GUIDER_RADIUS = config['guider'].getfloat('diameter')/2
GUIDER_SEC_RADIUS = config['guider'].getfloat('sec_diameter')/2
FINDER_RADIUS = config['finder'].getfloat('diameter')/2

# Dome:
RADIUS = config['dome'].getfloat('diameter')/2   # radius
EXTENT = config['dome'].getfloat('extent')  # extent of cylindrical dome wall
SLIT_WIDTH = config['dome'].getfloat('slit_width') # = MODEL / OLD: 1.84 (= measured)  # Slit width

# Observatory:
LAT = config['observatory'].getfloat('latitude') # degrees

class Instruments(enum.Enum):
    """
    Enum for selecting what aperture
    to use in the transformation.
    """
    TELESCOPE = enum.auto()
    GUIDER = enum.auto()
    FINDER = enum.auto()

    @classmethod
    def get_default(c):
        return c.TELESCOPE

def compute_azimuth(x, y):
    """
    Return the azimuth (rad) using
    the north clockwise convension.
    """
    az = np.arctan2(x, y) # note (x,y) rather than (y,x)
    
    if az < 0:
        az += 2*np.pi
    
    return az

def find_intersection(point, direction):
    """Find ray-capsule intersection."""
    has_intersection = False
    t = None
    
#     if direction.z < 0:
#         return has_intersection, t
    
    if np.isclose(direction[0], 0) and np.isclose(direction[1], 0):
        z = EXTENT + np.sqrt(RADIUS**2 - point[0]**2 - point[1]**2)
        t = z - point[2]
        
        has_intersection = True
        return has_intersection, t
    
    # If the direction vector is not (nearly) parallel to 
    # the z-axis of the capsule
    a2 = direction[0]**2 + direction[1]**2
    a1 = point[0]*direction[0] + point[1]*direction[1]
    a0 = point[0]**2 + point[1]**2 - RADIUS**2
    
    delta = a1**2-a0*a2
    t     = (-a1+np.sqrt(delta))/a2
    
    if point[2] + t * direction[2] >= EXTENT:
        a0 = point[0]**2 + point[1]**2 + (point[2] - EXTENT)**2 - RADIUS**2
        a1 = point[0]*direction[0] + point[1]*direction[1] + (point[2] - EXTENT)*direction[2]
        
        t = -a1+np.sqrt(a1**2-a0)
    
    if t:
        has_intersection = True
    
    return has_intersection, t

def get_ray_intersection(point, direction, t):
    """
    Return the ray intersection, based on the origin 
    (point) and direction vectors.
    """
    return point + t*direction

class Aperture:
    def __init__(self, radius, rate=100):
        self.radius = radius
        self.sample_rate = rate

        self._name = None

        # Add a vectorized instance of the _is_ray_blocked function
        self._is_blocked = np.vectorize(self._is_ray_blocked, signature='(d),(),(),(),()->()')

    def _transform(self, ha, dec):
        H_01 = transform(0, 0, L_1)
        H_12 = rot_x(90-LAT) @ rot_z(-ha) @ transform(0, 0, L_2)
        H_23 = rot_x(dec) @ transform(-L_3, 0, 0)

        H = H_01 @ H_12 @ H_23

        return H
    
    def _uniform_disk(self, r_min=0):
        """
        Sample points on a disk, uniformly.
        """
        radii  = np.random.uniform(0, 1, self.sample_rate)
        angles = np.random.uniform(0, 2*np.pi, self.sample_rate)
        
        x = np.sqrt(radii) * np.cos(angles)
        y = np.sqrt(radii) * np.sin(angles)
        
        xy = self.radius*np.column_stack([x, y])
        
        if r_min > 0:
            cond = x**2 + y**2 > r_min
            xy = xy[cond]

        return xy

    def _equidistant_disk(self, r_min=0):
        """Equidistant disk sampling based on:
            http://www.holoborodko.com/pavel/2015/07/23/generating-equidistant-points-on-unit-disk/
        """
    
        if not 0 <= r_min < 1:
            raise ValueError('r_min should be between 0 and 1...')
        
        dr = 1/self.sample_rate
        
        x = np.empty(0)
        y = np.empty(0)
        
        rs = np.linspace(r_min, 1, self.sample_rate) 
        k = np.ceil(r_min*(self.sample_rate+1))
        
        if not r_min:
            x = np.concatenate([x, [0]])
            y = np.concatenate([y, [0]])
            
            rs = np.linspace(dr, 1, self.sample_rate)
            k = 1
        
        for r in rs:
            n = int(np.round(np.pi/np.arcsin(1/(2*k))))
            
            theta = np.linspace(0, 2*np.pi, n+1)
            
            x_r = r * np.cos(theta)
            y_r = r * np.sin(theta)
            
            x = np.concatenate([x, x_r])
            y = np.concatenate([y, y_r])
            
            k += 1
        
        xy = self.radius*np.column_stack([x,y])
        
        # print('Generated {} points on a circle'.format(x.size))
        
        return xy
    
    def _sample_aperture(self, ha, dec, x, z):
        """
        Compute the position of a vector in 
        the aperture's frame.
        """
        y = np.zeros(x.size)
        dummy = np.ones(x.size)
        points = np.column_stack((x, y, z, dummy))
        
        pose_matrix = self._transform(ha, dec)

        product = pt.transform(pose_matrix, points)

        return product[:, :3]
    
    def _aperture_direction(self, ha, dec):
        H_ap = self._transform(ha, dec)
        H_unit = transform(0, 1, 0)

        H_diff = H_ap @ H_unit - H_ap

        direction = H_diff @ vec4(0, 0, 0)
        
        return vec3(direction)

    def _is_ray_blocked(self, point, ha, dec, dome_az, has_print=False):
        """
        Return True when the given ray in the aperture is blocked.
        """
        is_blocked = True
        
        direction = self._aperture_direction(ha, dec)

        try:
            has_intersection, t = find_intersection(point, direction)

            if has_intersection:
                points = get_ray_intersection(point, direction, t)
                
                ray_az = np.degrees(compute_azimuth(points[0], points[1]))

                if ray_az > dome_az + 180:
                    ray_az -= dome_az + 360
                elif ray_az > dome_az or ray_az < dome_az:
                    ray_az -= dome_az

                # Compute the dome slit - dome hemisphere intersection
                elev = np.arctan2(points[2] - EXTENT, np.sqrt(points[0]**2 + points[1]**2))
                
                y_length_sq = np.power(RADIUS*np.cos(elev), 2) - np.power(SLIT_WIDTH/2, 2)

                if y_length_sq < 0:
                    if points[2] > EXTENT:
                        # If on the other hand its below EXTENT => aperture's blocked
                        is_blocked = False
                    
                    return is_blocked

                y_length = np.sqrt(y_length_sq)
                x_length = SLIT_WIDTH/2
                slit_offset_rad = np.pi/2 - np.arctan2(y_length, x_length)
                
                # Compute off-set and compare ray az w/ dome position
                slit_offset = np.degrees(slit_offset_rad)

                slit_az_min = -slit_offset
                slit_az_max = slit_offset

                if has_print:
                    print('ray az = {:.2f} & dome az = {:.2f}'.format(np.degrees(compute_azimuth(points[0], points[1])), dome_az))
                    print('{:.2f} (min) < {:.2f} (ray) < {:.2f} (max)'.format(slit_az_min, ray_az, slit_az_max))
                
                is_ray_in_slit = ray_az > slit_az_min and ray_az < slit_az_max
                
                if points[2] > EXTENT and is_ray_in_slit:
                    is_blocked = False

                if not is_ray_in_slit and (ray_az < -90 or ray_az > 90):
                    rot = rot_z(dome_az)

                    dummy = np.ones(points[0].size)
                    pp = np.column_stack((points[0], points[1], points[2], dummy))

                    product = pt.transform(rot, pp)
                    
                    if np.abs(product[:, 0]) < SLIT_WIDTH/2 and np.abs(product[:,1]) < SLIT_WIDTH/2:
                        is_blocked = False
                
        except Exception as ex:
            print('ERROR OCCURRED DURING _is_blocked CALC...!\nERROR MSG:', str(ex))
            
        return is_blocked

    def obstruction(self, ha, dec, dome_az):
        ratio = None
        
        # Sample points in a disk; resembling the aperture
        ap_xz = self._equidistant_disk(self.radius)
        # ap_xz = self._uniform_disk(self.radius)
        
        ap_x, ap_z = ap_xz.T

        # Transfor those points to the aperture frame
        ap_pos = self._sample_aperture(ha, dec, -ap_x, ap_z)

        # Compute the no. rays, emanating from those points, blocked by the dome
        blocked = self._is_blocked(ap_pos, ha, dec, dome_az, has_print=False)
    
        ratio = blocked[blocked].size/blocked.size

        return ratio
    
    def get_name(self):
        return self._name

class TelescopeAperture(Aperture):
    def __init__(self, rate=100):
        # TODO: load & set telescope info
        super().__init__(APERTURE_RADIUS, rate=rate)
        
        self._name = 'telescope'

class GuiderAperture(Aperture):
    def __init__(self, rate=100):
        # TODO: load & set telescope info
        super().__init__(GUIDER_RADIUS, rate=rate)

        self._name = 'guider'

    def _transform(self, ha, dec):
        # Transform telescope aperture to guider aperture
        H_34 = transform(L_4*np.cos(GUIDER_ANGLE), 0, L_4*np.sin(GUIDER_ANGLE))

        # Get the telescope aperture pose
        H_telescope = super()._transform(ha, dec)

        H = H_telescope @ H_34

        return H

class FinderAperture(Aperture):
    def __init__(self, rate=100):
        # TODO: load & set telescope info
        super().__init__(FINDER_RADIUS, rate=rate)

        self._name = 'finder'

    def _transform(self, ha, dec):
        # Transform telescope aperture to guider aperture & guider to finder
        H_34 = transform(L_4*np.cos(GUIDER_ANGLE), 0, L_4*np.sin(GUIDER_ANGLE))
        H_45 = transform(-L_5*np.cos(FINDER_ANGLE), 0, L_5*np.sin(FINDER_ANGLE))

        # Get the telescope aperture pose
        H_telescope = super()._transform(ha, dec)

        H = H_telescope @ H_34 @ H_45

        return H
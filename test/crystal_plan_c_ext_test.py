import unittest
import numpy as np
import numpy.testing as npt
from collections import namedtuple

from model.crystals import PointGroup
from model.crystal_calc import getq
from model.detectors import FlatDetector
from model.numpy_utils import rotation_matrix, vector_length, normalize_vector
import crystal_plan_c_ext as ext

class CExtensionTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_q_with_different_angles(self):
        wavelength = 1.0
        azimuthal, elevation = 0., 0.
        rotation = np.eye(3)

        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_array_equal(q, np.array([0., 0., 0.]))

        azimuthal, elevation = 0., np.radians(90)
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([0., 1, -1]), 1e-3)

        azimuthal, elevation = np.radians(90), 0.
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([1, 0., -1]), 1e-3)

        azimuthal, elevation = np.radians(90), np.radians(90)
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([0., 1, -1]), 1e-3)

        azimuthal, elevation = np.radians(45), np.radians(45)
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([np.pi, np.sin(elevation)*2*np.pi, -np.pi]), 1e-3)

    def test_get_q_with_different_rotations(self):
        wavelength = 1.0
        azimuthal, elevation = 0., np.radians(45.)
        rotation = np.eye(3)

        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation), (np.cos(elevation) - 1)]), 1e-3)

        rotation = rotation_matrix(0, 0, np.radians(45))
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-1.30129027, 4.44288301, -1.30129027]), 1e-3)

        rotation = rotation_matrix(0, np.radians(45), 0)
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-np.pi, np.pi, -1.8403]), 1e-3)

        rotation = rotation_matrix(np.radians(45), 0, 0)
        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-1.30129027, 4.44288301, -1.30129027]), 1e-3)

    def test_get_q_with_different_wavelength(self):
        wavelength = 1.0
        azimuthal, elevation = 0., np.radians(45.)
        rotation = np.eye(3)

        q = ext.getq(wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation), (np.cos(elevation) - 1)]), 1e-6)

        q = ext.getq(2, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation) / 2., (np.cos(elevation) - 1)]), 1e-6)

        q = ext.getq(5, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation) / 2., (np.cos(elevation) - 1)]), 1e-6)

    def test_get_q_inelastic_with_different_angles(self):
        wavelength = 1.0
        azimuthal, elevation = 0., 0.
        rotation = np.eye(3)

        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_array_equal(q, np.array([0., 0., 0.]))

        azimuthal, elevation = 0., np.radians(90)
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([0., 1, -1]), 1e-3)

        azimuthal, elevation = np.radians(90), 0.
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([1, 0., -1]), 1e-3)

        azimuthal, elevation = np.radians(90), np.radians(90)
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi*np.array([0., 1, -1]), 1e-3)

        azimuthal, elevation = np.radians(45), np.radians(45)
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([np.pi, np.sin(elevation)*2*np.pi, -np.pi]), 1e-3)

    def test_get_q_inelastic_with_different_rotations(self):
        wavelength = 1.0
        azimuthal, elevation = 0., np.radians(45.)
        rotation = np.eye(3)

        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation), (np.cos(elevation) - 1)]), 1e-3)

        rotation = rotation_matrix(0, 0, np.radians(45))
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-1.30129027, 4.44288301, -1.30129027]), 1e-3)

        rotation = rotation_matrix(0, np.radians(45), 0)
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-np.pi, np.pi, -1.8403]), 1e-3)

        rotation = rotation_matrix(np.radians(45), 0, 0)
        q, _ = ext.getq_inelastic(wavelength, wavelength, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, np.array([-1.30129027, 4.44288301, -1.30129027]), 1e-3)

    def test_get_q_inelastic_with_different_wavelength(self):
        wavelength = 1.0
        azimuthal, elevation = 0., np.radians(45.)
        rotation = np.eye(3)

        q, _ = ext.getq_inelastic(1, 1, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation), (np.cos(elevation) - 1)]), 1e-6)

        q, _ = ext.getq_inelastic(2, 1, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation), (np.cos(elevation) - 1/2.)]), 1e-6)

        q, _ = ext.getq_inelastic(1, 2, azimuthal, elevation, np.pi, rotation)
        npt.assert_almost_equal(q, 2*np.pi * np.array([0., np.sin(elevation) / 2., (np.cos(elevation)/2. - 1)]), 1e-6)

    def test_get_detector_coordinates(self):
        detector = setup_detector()
        beam = np.array([[0, 0, 1]]).T

        _, array_size = beam.shape
        h_out = np.zeros( array_size )
        v_out = np.zeros( array_size )
        wl_out = np.zeros( array_size )
        distance_out = np.zeros( array_size )
        hits_it = np.zeros( array_size, dtype=np.bool)

        n_dot_base = np.dot(detector.base_point, detector.normal)
        wl_min, wl_max = 0, 100

        error_cout = ext.get_detector_coordinates(
            detector.base_point,
            detector.horizontal,
            detector.vertical,
            detector.normal,
            h_out, v_out, wl_out, distance_out, hits_it, beam,
            int(array_size), n_dot_base, int(detector.height), int(detector.width),
            wl_min, wl_max)

        self.assertAlmostEquals(wl_out, 6.283, 3)
        self.assertEquals(hits_it, True)

    def test_make_volume_symmetry_map(self):
        cell = (5,5,5, np.pi/2., np.pi/2., np.pi/2.)
        B = get_B_matrix(cell)
        B_inv = np.linalg.inv(B)

        point_group = PointGroup("-43m", "-43m (Cubic)", 24, 'lkh', "2;0,0,1", "2;0,1,0", "3;1,1,1", "-2;1,-1,0")
        order = len(point_group.table)
        table = np.array(point_group.table) #Turn the list of 3x3 arrays into a Nx3x3 array

        n = 10
        numpix = n**3
        symm = np.zeros( (numpix, order) , dtype=int)
        qres, qlim = 0.01, 10

        symm = ext.make_volume_symmetry_map(B, B_inv, symm, qres, qlim, n, order, table)
        self.assertEqual(len(symm[symm >= 0]), 6000)

    def test_angle_fitness(self):
        """
        Parameters:
            rot_angle_list: numpy array of the rotation of initial_rotation_matrix we are making, around the axis ending_vec
            initial_rotation_matrix: initial rotation matrix
            ending_vec: goal vector. MUST BE NORMALIZED TO 1
            starting_vec: start vector

        Return:
            fitness: metric of fitness to minimize
            best_angles: sample orientation angles corresponding to the given fitness
        """
        ending_vec = np.array([1, 0, 0])
        rot_angle_list = np.array([np.radians(0), np.radians(90)])

        size = rot_angle_list.shape[0] * 3
        fitnesses = np.zeros(size)
        chi_list = np.zeros(size)
        phi_list = np.zeros(size)
        omega_list = np.zeros(size)

        initial_rotation_matrix = np.eye(3)
        params = [0, np.pi*2, 0, np.pi*2, 0, np.pi*2]
        func_name = "general"

        ext.angle_fitness(rot_angle_list, ending_vec, initial_rotation_matrix, fitnesses, chi_list, phi_list, omega_list, func_name, params)

        npt.assert_almost_equal(fitnesses[:3], np.pi*2, 1e-3)
        npt.assert_almost_equal(fitnesses[3:], np.array([np.pi*3, np.pi*3, np.pi]), 1e-3)


    def test_calculate_coverage_stats(self):
        num = 30
        qlim = 20
        n = 100
        q_step = qlim/num

        covered_points0 = np.zeros(num)
        covered_points1 = np.zeros(num)
        covered_points2 = np.zeros(num)
        covered_points3 = np.zeros(num)
        total_points = np.zeros(num)

        qspace_radius = np.array([1])
        qspace = np.zeros((n, n, n))
        qspace = qspace.flatten()
        qspace[:100] = 1  # 1% of total volume is covered
        qspace[25:75] = 2 # 0.5% of total volume is covered more than once
        qspace_size = qspace.size

        results = ext.calculate_coverage_stats(qspace, qspace_radius, q_step, qlim, total_points, qspace_size, num, covered_points0, covered_points1, covered_points2, covered_points3)

        stats, total_points, covered_points0, covered_points1, covered_points2, covered_points3 = results
        (overall_points, overall_covered_points, overall_redundant_points) = stats

        overall_coverage = 100.0 * overall_covered_points / overall_points;
        overall_redundancy = 100.0 * overall_redundant_points / overall_points;

        self.assertEquals(overall_coverage, 0.01)
        self.assertEquals(overall_redundancy, 0.005)

    def test_calculate_coverage(self):
        det = FlatDetector('test-bank')
        det.calculate_pixel_angles()

        params = get_coverage_params(det, "elastic")
        coverage = ext.calculate_coverage(*params)


        frac = np.count_nonzero(coverage) / float(coverage.size) * 100
        self.assertAlmostEqual(frac, 0.024, places=3)


    def test_calculate_coverage_inelastic(self):
        det = FlatDetector('test-bank')
        det.calculate_pixel_angles()

        params = get_coverage_params(det, "inelastic")
        coverage = ext.calculate_coverage_inelastic(*params)


        self.assertEqual(np.min(coverage), 0.)
        self.assertEqual(np.max(coverage[coverage <1e6]), .0)
        self.assertEqual(coverage[coverage <1e6].size, 74087)

#---------------------------------------------------------------------------------------------

def get_coverage_params(det, emode):
    qlim = 10
    q_resolution = .1
    wl_input = 0.5
    wl_min, wl_max = 0.4, 0.8
    rot_matrix = np.eye(3)

    q0 = getq(det.azimuthal_angle[0, 0], det.elevation_angle[0, 0], wl_min, rot_matrix, wl_input=wl_input)
    q_xmax = getq(
        det.azimuthal_angle[0, -1], det.elevation_angle[0, -1], wl_min, rot_matrix, wl_input=wl_input)
    q_ymax = getq(
        det.azimuthal_angle[-1, 0], det.elevation_angle[-1, 0], wl_min, rot_matrix, wl_input=wl_input)

    # Make sure they aren't too long
    length = max( vector_length(q_xmax),  vector_length(q_ymax),  vector_length(q0) )
    q_xmax = normalize_vector(q_xmax, length)
    q_ymax = normalize_vector(q_ymax, length)

    qx_list = np.arange(-qlim, qlim, q_resolution)
    n = qx_list.size

    if emode == "inelastic":
        coverage = np.zeros((n,n,n), dtype=np.float) + 1e6
    elif emode == "elastic":
        coverage = np.zeros((n,n,n,1), dtype=np.int)


    nx = (vector_length(q_xmax - q0) / q_resolution) * 1.5
    ny = (vector_length(q_ymax - q0) / q_resolution) * 1.5

    azimuthal_angle = det.azimuthal_angle
    elevation_angle = det.elevation_angle

    #Dimensions of the array
    s = coverage.shape
    stride = s[0]
    max_index = s[0]*s[1]*s[2] #largest index into the array +1

    xlist = [x for x in set(np.linspace(
        0, det.xpixels - 1, nx, endpoint=True).astype(int))]
    xlist.sort()
    ylist = [x for x in set(np.linspace(
        0, det.ypixels - 1, ny, endpoint=True).astype(int))]
    ylist.sort()

    xlist = np.array(xlist)
    ylist = np.array(ylist)

    set_value1 = 0
    set_value2 = 0

    params = (xlist, ylist, azimuthal_angle,
        elevation_angle, rot_matrix, set_value1,
        set_value2, n, coverage,
        stride, max_index, wl_min, wl_max,
        qlim, q_resolution)

    if emode == "inelastic":
        params = (wl_input, ) + params

    return params

def unit_cell_from_metric_tensor(G):
    a, b, c = np.sqrt(np.diag(G))
    alpha = np.arccos(G[1, 2] / b / c)
    beta = np.arccos(G[0, 2] / a / c)
    gamma = np.arccos(G[0, 1] / a / b)
    return a, b, c, alpha, beta, gamma

def calculate_metric_tensor(unit_cell):
    a, b, c, alpha, beta, gamma = unit_cell
    G = np.eye(3)
    G[0,0] = a * a
    G[1,1] = b * b
    G[2,2] = c * c
    G[0,1] = a * b * np.cos(gamma)
    G[0,2] = a * c * np.cos(beta)
    G[1,2] = b * c * np.cos(alpha)
    G[1,0] = G[0,1]
    G[2,0] = G[0,2]
    G[2,1] = G[1,2]
    return G

def reciprocate(unit_cell):
    G = calculate_metric_tensor(unit_cell)
    G_star = np.linalg.inv(G)
    return unit_cell_from_metric_tensor(G_star)

def get_B_matrix(unit_cell):
    # B matrix using a right handed coordinate system with b1 along x and y in
    # the (b1,b2) plane.
    # This is the convention in Busing and Levy.
    # | b1 b2cos(beta3)      b3cos(beta2)        |
    # | 0  b2sin(beta3) -b3sin(beta2)cos(alpha1) |
    # | 0       0                  1/a3          |
    a1, a2, a3, alpha1, alpha2, alpha3 = unit_cell
    b1, b2, b3, beta1, beta2, beta3 = reciprocate(unit_cell)
    B = np.zeros((3,3))
    B[0,0] = b1
    B[0,1] = b2 * np.cos(beta3)
    B[0,2] = b3 * np.cos(beta2)
    B[1,0] = 0.
    B[1,1] = b1 * np.sin(beta3)
    B[1,2] = -b3 * np.sin(b2) * np.cos(alpha1)
    B[2,0] = 0.
    B[2,1] = 0.
    B[2,2] = 1. / a3
    return B

def setup_detector():
    Detector = namedtuple(
        'Detector', ['base_point', 'horizontal', 'vertical', 'normal', 'height', 'width',
                     'azimuthal_angle', 'elevation_angle', 'xpixels', 'ypixels'])
    return Detector(
        base_point = np.array([0, 0, 10]),
        horizontal = np.array([1, 0, 0]),
        vertical = np.array([0, 1, 0]),
        normal = np.array([0, 0, -1]),
        height = 10,
        width = 10,
        azimuthal_angle = [0., 0.],
        elevation_angle = [0., 0.],
        xpixels = 256,
        ypixels = 256
    )

if __name__ == "__main__":
    unittest.main()

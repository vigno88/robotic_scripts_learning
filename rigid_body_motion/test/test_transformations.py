import unittest

import numpy as np
import transformations as trans

R = np.array([[-0.263, -0.14, 0.955],
              [ 0.857,  0.42, 0.298],
              [-0.443, 0.897, 0.01 ]])
p = np.array([1,2,3])
s_a = np.array([0.254,0.381,0.889,2.53,-0.083,3.194])

class TestTransformations(unittest.TestCase):

  def test_rotation_inverse(self):
    inn = np.array([[0,-1,0],[0,0,1],[1,0,0]])    
    c_out = np.array([[0,0,1],[-1,0,0],[0,1,0]])    
    
    outt = trans.rotation_inverse(inn)
    np.testing.assert_allclose(c_out, outt)
  
  def test_so3(self):
    vec = np.array([0.54,0.78,0.32])
    so3 = np.array([[0, -0.32, 0.78],
                    [0.32, 0, -0.54],
                    [-0.78, 0.54, 0]])

    outt = trans.vec_to_so3(vec)
    np.testing.assert_allclose(so3,outt)

    outt = trans.so3_to_vec(so3)
    np.testing.assert_allclose(vec,outt)
  
  def test_exp_coords_3(self):
    axis = np.array([0.54,0.78,0.32])
    angle = 2.0
    exp_coords = angle * axis

    axis_o, angle_o = trans.axis_angle_3(exp_coords)
    np.testing.assert_allclose(axis,axis_o, rtol=0.1)
    np.testing.assert_allclose(angle,angle_o, rtol=0.1)

  def test_matrix_exp_3(self):
    omega_theta = np.array([[0, -1.098, 1.537],
                            [1.098, 0, -0.659],
                            [-1.537, 0.659, 0]])
    R_o = trans.matrix_exp_3(omega_theta)
    # The correct rotation matrix was computed by hand
    np.testing.assert_allclose(R,R_o, rtol=0.1)

  def test_matrix_log_3(self):
    # There is 3 case:
    # If R = I, we get an exception
    self.assertRaises(Exception,trans.matrix_log_3,np.identity(3))

    # If tr R == -1
    # The correct omega_theta matrix was computer by hand
    c_o_t = np.array([[0, -2.22 ,2.22],
                      [2.22, 0, 0],
                      [-2.22, 0, 0]])
    R_tr = np.array([[-1,0,0],
                  [0,0,1],
                  [0,1,0]]) 
    
    o_t = trans.matrix_log_3(R_tr)
    np.testing.assert_allclose(c_o_t,o_t, rtol=0.1)
    
    # Else
    # The correct omega_theta matrix was computer by hand
    c_o_t = np.array([[0, -1.098, 1.537],
                      [1.098, 0, -0.659],
                      [-1.537, 0.659, 0]])
    o_t = trans.matrix_log_3(R)
    np.testing.assert_allclose(c_o_t,o_t, rtol=0.1)
  
  def test_Rp_t(self):
    # Go one way
    T = trans.Rp_to_T(R,p)    
    # Go the other way
    Ro, po = trans.T_to_Rp(T)

    np.testing.assert_allclose(R,Ro)
    np.testing.assert_allclose(p,po)

  def test_T_inverse(self):
    T = trans.Rp_to_T(R,p)     
    # Go one way 
    T_inv = trans.T_inverse(T)
    # Go the other way
    To = trans.T_inverse(T_inv)
    np.testing.assert_allclose(T,To, rtol=0.1)
    
  def test_twist_to_se3(self):
    twist = np.array([0.487,0.324,0.811,0.62,0.248,0.744])
    c_se3 = np.array([[0, -0.811, 0.324, 0.62],
                      [0.811, 0, -0.487, 0.248],
                      [-0.324, 0.487, 0, 0.744],
                      [0,0,0,0]])
    
    se3_o = trans.twist_to_se3(twist)
    np.testing.assert_allclose(c_se3,se3_o)
   
  def test_T_adjoint(self):
    T = trans.Rp_to_T(R,p)
    # The adjoint as been computer by hand
    adj_T = np.array([[-0.263,  -0.14,  0.955, 0, 0, 0],
                      [ 0.857,   0.42,  0.298, 0, 0, 0],
                      [-0.443,  0.897,  0.01 , 0, 0, 0],
                      [-3.457,  0.534, -0.872,-0.263,  -0.14,  0.955],
                      [-0.345, -1.319,  2.853, 0.857,   0.42,  0.298],
                      [1.383,   0.701, -1.611,-0.443,  0.897,  0.01 ]])
    out = trans.adjoint(T)
    np.testing.assert_allclose(adj_T,out,rtol=0.1)
   
  def test_screw_to_axis(self):
    s = np.array([0.254,0.381,0.889])
    q = np.array([3,4,5])
    h = 3.45

    screw_axis_o = trans.screw_to_axis(q,s,h)
    np.testing.assert_allclose(s_a, screw_axis_o, rtol=0.1)
 
  def test_screw_mat(self):
    # Go one way
    mat = trans.screw_axis_to_mat(s_a)
    # Go back
    s_a_o = trans.mat_to_screw_axis(mat)
    np.testing.assert_allclose(s_a, s_a_o, rtol=0.1)
  
  def test_matrix_exp_log_6(self):
    # Make a screw axis matrix, 2 is theta
    mat = trans.screw_axis_to_mat(s_a * 2)

    T_o = trans.matrix_exp_6(mat)
    mat_o = trans.matrix_log_6(T_o)
    np.testing.assert_allclose(mat, mat_o, rtol=0.1)



     

if __name__ == '__main__':
  unittest.main()

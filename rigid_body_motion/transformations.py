import numpy as np

# Computes the inverse of rotation matrix R
def rotation_inverse(R):
  # Inverse of rotation matrix is equal to its transpose
  return R.T

# Returns the 3x3 skew-symmetric matrix corresponding to v
def vec_to_so3(v):
  if len(v) != 3:
    raise Exception("""Cannot create a matrix of a vector with length
    different than 3""")
  return np.array([0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0])

# Returns the 3-vector corresponding to 3x3 skew-symmetric matrix M
def so3_to_vec(M):
  if M.shape != (3,3):
    raise Exception("Cannot create a vector of a non 3X3 Matrix")
  return np.array([M[2,1],M[0,2],M[1,0]]) 

# Extracts rotation axis omega_hat and the rotation amount theta from
# exponential coordinates vector for rotation
# Return the axis and the angle
def axis_angle_3(exp_v):
  ang = np.linalg.norm(exp_v)
  return exp_v/ang, ang

# Compute the rotation matrix R corresponding to matrix exponential of M
def matrix_exp_3(M):
  ang = np.linalg.norm(so3_to_vec(M))
  axis_m = M/ang
  return np.identity(3) + np.multiply(np.sin(ang),axis_m) + 
    np.multiply(1-np.cos(ang),np.square(axis_m))
    

# Compute matrix log of rotation matrix R
def matrix_log_3(R):
  axis, angle = (0,0)
  if R == np.identity(3):
    raise Exception("axis of rotation is undefined")
  if np.trace(R) == -1:
    axis = np.array([R[0,2],R[1,2],1+R[2,2]]).multiply(1/np.sqrt(2*(1+R[2,2])))
    return vec_to_so3(axis).multiply(np.pi)
  else:
    angle = np.arccos(0.5*(np.trace(R)-1))
    M = (R-R.T)/(2*np.sin(angle))
    return M.multiply(angle)

# Return homogeneous transformation matrix T corresponding to rotation matrix R
# and position vector p
def  Rp_to_T(R,p):
  return np.array([R[0,0], R[0,1], R[0,2], p[0]],
                  [R[1,0], R[1,1], R[1,2], p[1]],
                  [R[2,0], R[2,1], R[2,2], p[2]],
                  [0     , 0     , 0     , 1])
# Extracts the rotation matrix and position vector from a homogeneous
# transformation matrix T
def T_to_Rp(T):
  return T[0:3,0:3], np.array([T[3,0],T[3,1],T[3,2]])

# Computes inverse of homogeneous transformation matrix T
def T_inverse(T):
  R,p = T_to_Rp(T)
  return Rp_to_T(R.T, -R.T.mutliply(p.T))

# Return the se(3) matrix corresponding to 6-vector twist V
def twist_to_se3(V):
  omega,v = (V[0:3],V[3:6])
  o_m = vec_to_so3(omega)
  return np.array([o_m[0,0], o_m[0,1], o_m[0,2], v[0]],
                  [o_m[1,0], o_m[1,1], o_m[1,2], v[1]],
                  [o_m[2,0], o_m[2,1], o_m[2,2], v[2]],
                  [0       , 0       , 0       , 0])

# Computes the 6x6 adjoint representation [Ad_T] of the homogeneous
# transformation matrix T
def adjoint(T):
  R,p = T_to_Rp(T)
  p_m = vec_to_so3(p)
  up = np.concatenate((R,np.zeros(3,3)),axis=1)
  down = np.concatenate((p_m.multiply(R)), R,axis=1)
  return np.concatenate((up,down), axis=0)
  

# Return a normalized screw axis represention S of a screw described by a unit
# vector s in the dir of screw axis, located a point q, with pitch h
def screw_to_axis(q,s,h):
  v = np.cross(-s,q) + np.multiply(h,s)
  return np.concatenate((s,v), axis=1)

# Extracts normalized screw axis S and the distance traveled along the screw and
# the distance traveled along the screw theta from the 6-vector of exponential
# coords S*theta
def axis_angle(exp_6):
  theta = np.linalg.norm(exp_6)
  return exp_6/theta, theta

def screw_axis_to_mat(s):
  up = np.concatenate((vec_to_so3(s[0:3]),s[3:6]), axis=1)
  return np.concatenate((up,np.zeros(1,4),axis=0)

def mat_to_screw_axis(S_m):
  w = so3_to_vec(S_m[0:3,0:3])
  v = np.array([S_m[0,3],S_m[1,3],S_m[2,3]])
  return np.concatenate(w,v)

# Computes the homogeneous transformation T corresponding to the matrix
# exponential of screw axis matrix
def matrix_exp_6(S_m):
  s_a = mat_to_screw_axis(S_m)
  angle = np.linalg.norm(s_a[0:3])
  top_left = matrix_exp_3(vec_to_so3(s_a[0:3]))
  mat_omega = vect_to_so3(s_a[0:3]/angle)
  top_right = ((np.identity(3).multiply(angle)) + 
    mat_omega.multiply(1-np.cos(angle)) +
    np.square(mat_omega).multiply(angle-np.sin(angle))).
    multiply(s_a[3:6][...,np.newaxis])
  top = np.concatenate((top_left,top_right.T),axis=1)
  return np.concatenate((top, np.array([0,0,0,1])))

# Compute matrix logarithm of the homogeneous transformation matrix T
def matrix_log_6(T):



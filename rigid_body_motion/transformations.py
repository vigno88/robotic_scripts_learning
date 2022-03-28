import numpy

# Computes the inverse of rotation matrix R
def rotation_inverse(R):

# Returns the 3x3 skew-symmetric matrix corresponding to v
def vec_to_so3(v):

# Returns the 3-vector corresponding to 3x3 skew-symmetric matrix M
def so3_to_vec(M):

# Extracts rotation axis omega_hat and the rotation amount theta from
# exponential coordinates vector for rotation
def axis_angle_3(exp_v):

# Compute the rotation matrix R corresponding to matrix exponential of M
def matrix_exp_3(M):

# Compute matrix log of rotation matrix R
def matrix_log_3(R):

# Return homogeneous transformation matrix T corresponding to rotation matrix R
# and position vector p
def  Rp_to_T(R,p):

# Extracts the rotation matrix and position vector from a homogeneous
# transformation matrix T
def T_to_Rp(T):

# Computes inverse of homogeneous transformation matrix T
def T_inverse(T):

# Return the se(3) matrix corresponding to 6-vector twist V
def twist_to_se3(V):

# Computes the 6x6 adjoint representation [Ad_T] of the homogeneous
# transformation matrix T
def adjoint(T):

# Return a normalized screw axis represention S of a screw described by a unit
# vector s in the dir of screw axis, located a point q, with pitch h
def screw_to_axis(q,s,h):

# Extracts nomralized screw axis S and the distance traveled along the screw and
# the distance traveled along the screw theta from the 6-vector of exponential
# coords S*theta
def axis_angle(exp_6):

# Computes the homogeneous transformation T corresponding to the matrix
# exponential of M
def matrix_exp_6(M):

# Compute matrix logarithm of the homogeneous transformation matrix T
def matrix_log_6(T):



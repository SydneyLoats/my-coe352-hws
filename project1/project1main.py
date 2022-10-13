import numpy as np
import sys
from scipy.linalg import svd

#take in values for number of springs/masses, spring constants (c) and masses (m)

#if len(sys.argv) == 1:
#  print(f'Please enter input values')
#  print(f'Usage:')
#  print(f'python3 project1main.py [number of springs] [number of masses]')
#  print(f'Example:')
#  print(f'python3 project1main.py 4 3')
#  exit()

#spring_count = int(sys.argv[1])
#mass_count = int(sys.argv[2])
spring_count = 4
mass_count = 3
print(f'Spring count: {spring_count}')
print(f'Mass count: {mass_count}')

spring_constant = [2, 2, 2, 2]
m_vec = [3, 3, 3]

f_vec = np.multiply(m_vec, 9.8)
print(f'f vector is {f_vec}')




def create_identity_matrix(val_vec):
  ret_mat = np.array([[0 for c in range(len(val_vec))] for r in range(len(val_vec))])

  for r in range(len(val_vec)):
    for c in range(len(val_vec)):
      if r == c:
        ret_mat[r][c] = val_vec[c]
  return ret_mat

mat = create_identity_matrix(spring_constant)
print(mat)


def a_mat_fixed_fixed():
  ret_mat = np.array([[0 for c in range(mass_count)] for r in range(spring_count)])

  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        ret_mat[r][c] = 1 
      elif r == c+1:
        ret_mat[r][c] = -1
  return ret_mat
 


def a_mat_fixed_free():  
  ret_mat = np.array([[0 for c in range(mass_count)] for r in range(spring_count)])

  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        ret_mat[r][c] = 1 
      elif c == r-1:
        ret_mat[r][c] = -1
  return ret_mat



def a_mat_free_free():
  ret_mat = np.array([[0 for c in range(mass_count)] for r in range(spring_count)])

  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        ret_mat[r][c] = -1
      elif c == r+1:
        ret_mat[r][c] = 1
  return ret_mat


#solve for displacement
def solve_displacement(mat, f_vec):

#creating A matrix
  a_mat = mat
  print(f'A matrix is \n{a_mat}')

#A transpose
  at_mat = a_mat.transpose()
  print(f'A transpose matrix is \n{at_mat}')

#A transpose inverse
  ati_mat = np.linalg.pinv(at_mat)
  print(f'Inverse A transpose matrix is \n{ati_mat}')

#finding w
  w_vec = np.matmul(ati_mat, f_vec)
  print(f'w vector is \n{w_vec}')

#identity matrix
  c_mat = create_identity_matrix(spring_constant)
  print(f'C matrix is \n{c_mat}')

#C inverse
  ci_mat = np.linalg.pinv(c_mat)
  print(f'C inverse matrix is \n{ci_mat}')

#finding e
  e_mat = np.matmul(w_vec, ci_mat)
  print(f'e matrix is \n{e_mat}')

#a inverse
  ai_mat = np.linalg.pinv(a_mat)
  print(f'a inverse is \n{ai_mat}')

#finding u (displacement)
  u_mat = np.matmul(ai_mat, e_mat)
  print(f'u matrix is \n{u_mat}')

#stiffness matrix K
#  temp_mat = np.matmul(at_mat, c_mat)
#  k_mat = np.matmul(temp_mat, a_mat)
#  print(f'K matrix is \n{k_mat}')

#inverse stiffness matrix
#  ki_mat = np.linalg.pinv(k_mat)
#  print(f'K inverse is \n{ki_mat}')

#solving for displacement
#  u_mat = np.matmul(ki_mat, f_vec)
#  print(f'Displacement u is \n{u_mat}')

  return u_mat


#def solve_svd(mat):

#  a_mat = mat
#  U, s, VT = svd(a_mat)
#  print(f'Singular values are \n{s}')
#  return s

#def solve_eigenvalues(s):

#  eigen = np.array(s) 
#  for r in range(
#  print(f'Eigenvalues are \n{eigen}')

#solve_svd(a_mat_fixed_fixed())
#solve_eigenvalues(s)

solve_displacement(a_mat_fixed_fixed(), f_vec)



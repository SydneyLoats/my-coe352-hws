import numpy as np
import sys

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

#identity matrix
  c_mat = create_identity_matrix(spring_constant)
  print(f'C matrix is \n{c_mat}')

#stiffness matrix K
  temp_mat = np.matmul(at_mat, c_mat)
  k_mat = np.matmul(temp_mat, a_mat)
  print(f'K matrix is \n{k_mat}')

#inverse stiffness matrix
  ki_mat = np.linalg.pinv(k_mat)
  print(f'K inverse is \n{ki_mat}')

#solving for displacement
  u_mat = np.matmul(ki_mat, f_vec)
  print(f'Displacement u is \n{u_mat}')



#
#A transpose inverse matrix
#  ati_mat = np.linalg.pinv(at_mat)
#  print(f'A transpose inverse matrix is \n{ati_mat}')

#solving for w
#  w_vec = np.matmul(ati_mat, f_vec)
#  print(f'w vector is \n{w_vec}')

#solving for e
#  ii_mat = np.linalg.pinv(i_mat)
#  e_vec = np.matmul(ii_mat, w_vec) 
#  print(f'e (elongation) vector is \n{e_vec}')

#solving for u (displacement)
#  ai_mat = np.linalg.pinv(a_mat)
#  print(f'A inverse matrix is \n{ai_mat}')
#  u_vec = np.matmul(ai_mat, e_vec)
#  print(f'u (displacement) vector is \n{u_vec}')


solve_displacement(a_mat_fixed_fixed(), f_vec)



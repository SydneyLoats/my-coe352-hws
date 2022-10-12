import numpy as np
import sys

#take in values for number of springs/masses, spring constants (c) and masses (m)

if len(sys.argv) == 1:
  print(f'Please enter input values')
  print(f'Usage:')
  print(f'python3 project1main.py [number of springs] [number of masses]')
  print(f'Example:')
  print(f'python3 project1main.py 4 3')
  exit()

spring_count = int(sys.argv[1])
mass_count = int(sys.argv[2])
print(f'Spring count: {spring_count}')
print(f'Mass count: {mass_count}')

spring_constant = [2, 2, 2, 2]
m_vec = [3, 3, 3]

f_vec = np.multiply(m_vec, 9.8)
print(f'f vector is {f_vec}')




def create_identity_matrix(size, val):
  ret_mat = np.array([[0 for c in range(size)] for r in range(size)])

  for r in range(size):
    for c in range(size):
      if r == c:
        ret_mat[r][c] = val
  return ret_mat



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

#test_mat = a_mat_fixed_fixed()
#print(f'a matrix fixed fixed \n{test_mat}')




def solve(mat):

#creating A matrix
  a_mat = mat
  print(f'A matrix is \n{a_mat}')

#identity matrix
  i_mat = create_identity_matrix(len(spring_constant), spring_constant[0])
  print(f'I matrix is \n{i_mat}')

#A transpose
  at_mat = a_mat.transpose()
  print(f'A transpose matrix is \n{at_mat}')

#A transpose inverse matrix
  ati_mat = np.linalg.pinv(at_mat)
  print(f'A transpose inverse matrix is \n{ati_mat}')

#solving for w
  w_vec = np.matmul(ati_mat, f_vec)
  print(f'w vector is \n{w_vec}')

#solving for e
  ii_mat = np.linalg.pinv(i_mat)
  e_vec = np.matmul(ii_mat, w_vec) 
  print(f'e (elongation) vector is \n{e_vec}')

#solving for u (displacement)
  ai_mat = np.linalg.pinv(a_mat)
  print(f'A inverse matrix is \n{ai_mat}')
  u_vec = np.matmul(ai_mat, e_vec)
  print(f'u (displacement) vector is \n{u_vec}')


solve(a_mat_fixed_fixed())



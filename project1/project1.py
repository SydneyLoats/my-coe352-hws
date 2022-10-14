import numpy as np
import sys
from scipy.linalg import svd

#assign number of springs and spring constants
global spring_count
spring_count = int(input('Input number of springs (positive integer): '))
if spring_count < 1:
  print(f'\nError: the input {spring_count} is invalid, please choose a positive integer')
  exit()
print(f'Number of springs input is {spring_count}')

global spring_const_vec
spring_const_vec = np.array([1 for i in range(spring_count)])
for i in range(spring_count):
  spring_const_vec[i] = int(input(f'Input spring constant {i+1}: '))
print(f'Spring constant vector is {spring_const_vec}')

#assign number of masses and mass values
global mass_count
mass_count = int(input('Input number of masses (positive integer): '))
if mass_count < 1:
  print(f'\nError: the input {mass_count} is invalid, please choose a positive integer')
  exit()
print(f'Number of masses input is {mass_count}')

global m_vec
m_vec = np.array([1 for i in range(mass_count)])
for i in range(mass_count):
  m_vec[i] = int(input(f'Input mass value {i+1}: '))
print(f'Mass values vector is {m_vec}')

#determine whether to use fixed-fixed or fixed-free
global system_type
system_type = int(input(f'Input boundary condition(either 0 or 1 where 0=fixed-fixed and 1=fixed-free): '))
print(f'Boundary input condition is {system_type}')
if (spring_count - mass_count) == 1 and system_type != 0:
  print('\nError: the number of masses and springs do not match the boundary condition, please choose different numbers')
  exit()
elif (spring_count - mass_count) == 0 and system_type !=1:
  print('\nError: the number of masses and springs do not match the boundary condition, please choose different numbers')
  exit()
 
#create external force vector, f
global f_vec
f_vec = np.multiply(m_vec, 9.8)

#function to create an identity matrix, uses spring constant vector
def create_identity_matrix():
  ret_mat = np.array([[0 for c in range(len(spring_const_vec))] for r in range(len(spring_const_vec))])
  for r in range(len(spring_const_vec)):
    for c in range(len(spring_const_vec)):
      if r == c:
        ret_mat[r][c] = spring_const_vec[c]
  return ret_mat

#function to create the A matrix for a fixed_fixed system
def a_mat_fixed_fixed():
  ret_mat = np.array([[0 for c in range(mass_count)] for r in range(spring_count)])
  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        ret_mat[r][c] = 1 
      elif r == c+1:
        ret_mat[r][c] = -1
  return ret_mat

#function to create the A matrix for a fixed-free system
def a_mat_fixed_free():  
  ret_mat = np.array([[0 for c in range(mass_count)] for r in range(spring_count)])
  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        ret_mat[r][c] = 1 
      elif c == r-1:
        ret_mat[r][c] = -1
  return ret_mat

#function to solve for displacement, takes in the a matrix
def solve_displacement(mat):
#creating A matrix
  a_mat = mat
#  print(f'A matrix is \n{a_mat}')

  sv = solve_svd(a_mat)
  print(f'\nSingular values for A: \n{sv}')

  eigen_sv = solve_eigenvalues(sv)
  print(f'\nEigenvalues for A: \n{eigen_sv}')

  condition_value_a = find_largest_eigen(eigen_sv)/find_smallest_eigen(eigen_sv)
  print(f'\nCondition value for A: \n{condition_value_a}')

#creating A transpose matrix
  at_mat = a_mat.transpose()
#  print(f'A transpose matrix is \n{at_mat}')

  sva = solve_svd(at_mat)
  print(f'\nSingular values for A transpose: \n{sva}')

  eigen_sva = solve_eigenvalues(sva)
  print(f'\nEigenvalues for A transpose: \n{eigen_sva}')

  condition_value_a = find_largest_eigen(eigen_sva)/find_smallest_eigen(eigen_sva)
  print(f'\nCondition value for A transpose: \n{condition_value_a}')

#creating identity matrix
  c_mat = create_identity_matrix()
#  print(f'C matrix is \n{c_mat}')

  svc = solve_svd(c_mat)
  print(f'\nSingular values for C: \n{svc}')

  eigen_svc = solve_eigenvalues(svc)
  print(f'\nEigenvalues for C: \n{eigen_svc}')
  
  condition_value_c = find_largest_eigen(eigen_svc)/find_smallest_eigen(eigen_svc)
  print(f'\nCondition value for C: \n{condition_value_c}')

#creating stiffness matrix, K
  tempk_mat = np.matmul(at_mat, c_mat)
  k_mat = np.matmul(tempk_mat, a_mat)

  svk = solve_svd(k_mat)
  print(f'\nSingular values for K: \n{svk}')

  eigen_svk = solve_eigenvalues(svk)
  print(f'\nEigenvalues for K: \n{eigen_svk}')

  condition_value_k = find_largest_eigen(eigen_svk)/find_smallest_eigen(eigen_svk)
  print(f'\nCondition value for K: \n{condition_value_k}')  

#creating K inverse matrix
  ki_mat = np.linalg.pinv(k_mat)

#solving u matrix (displacement)
  u_vec = np.matmul(ki_mat, f_vec)

#creating e vector (elongation)
  e_vec = np.matmul(a_mat, u_vec)
  print(f'\nElongation vector, e, is \n{e_vec}')

#creating w vector (internal forces)
  w_vec = np.matmul(c_mat, e_vec)
  print(f'\nInternal stress vector, w, is \n{w_vec}')

#print out the displacement
  print(f'\nu vector is \n{u_vec}')
  for i in range(len(u_vec)):
    print(f'Displacement {i+1} = {u_vec[i]}')
  return u_vec

#function to solve svd
def solve_svd(mat):
  a_mat = mat
  U, s, VT = svd(a_mat)
  return s

#function to solve for eigenvalues
def solve_eigenvalues(s):
  eigen = np.array(s) 
  for r in range(len(s)):
    eigen[r] = s[r]*s[r]
  return eigen

#function to find the largest eigenvalue
def find_largest_eigen(eigen):
  largest = eigen[0]
  for r in range(len(eigen)):
    if eigen[r] > largest:
      largest = eigen[r]
  return largest

#function to find the smallest eigenvalue
def find_smallest_eigen(eigen):
  smallest = eigen[0]
  for r in range(len(eigen)):
    if eigen[r] < smallest:
      smallest = eigen[r]
  return smallest

#solve the system
if system_type == 0:
  solve_displacement(a_mat_fixed_fixed())
else:
  solve_displacement(a_mat_fixed_free())



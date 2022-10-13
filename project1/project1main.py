import numpy as np
import sys
from scipy.linalg import svd

spring_count = int(input('Number of springs: '))
print(spring_count)
spring_const_vec = np.array([1 for i in range(spring_count)])
#print(spring_const_vec)
for i in range(spring_count):
  spring_const_vec[i] = int(input('Spring constant {i}: '))
print(spring_const_vec)

#type = 0  #where 0=fixed-fixed, 1=fixed-free, 2=free-free

mass_count = int(input('Number of masses: '))
print(mass_count)
m_vec = np.array([1 for i in range(mass_count)])
#print(m_vec)
for i in range(mass_count):
  m_vec[i] = int(input('Mass constant {i}: '))
#print(m_vec)
#def check_mass_spring_count():
#if spring_count - mass_count == 1
#  type = 0
 #   return
#elif spring_count - mass_count == 0
#  type = 1
  #  return
#elif mass_count - spring_count == 1
#  type = 2
  #  return
#else
#  print('\nInvalid number of masses, please select a different number')
#  mass_count = int(input('Number of masses: ')
  #  return
   # check_mass_spring_count()
  
#check_mass_spring_count()
print(mass_count)

m_vec = [3, 3, 3]

f_vec = np.multiply(m_vec, 9.8)
#print(f'f vector is {f_vec}')




def create_identity_matrix(val_vec):
  ret_mat = np.array([[0 for c in range(len(val_vec))] for r in range(len(val_vec))])

  for r in range(len(val_vec)):
    for c in range(len(val_vec)):
      if r == c:
        ret_mat[r][c] = val_vec[c]
  return ret_mat

mat = create_identity_matrix(spring_const_vec)
#print(mat)


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
#  print(f'A matrix is \n{a_mat}')

#A transpose
  at_mat = a_mat.transpose()
#  print(f'A transpose matrix is \n{at_mat}')

#A transpose inverse
  ati_mat = np.linalg.pinv(at_mat)
#  print(f'Inverse A transpose matrix is \n{ati_mat}')

#finding w
  w_vec = np.matmul(ati_mat, f_vec)
#  print(f'w vector is \n{w_vec}')

#identity matrix
  c_mat = create_identity_matrix(spring_constant)
#  print(f'C matrix is \n{c_mat}')


#C inverse
  ci_mat = np.linalg.pinv(c_mat)
#  print(f'C inverse matrix is \n{ci_mat}')

#finding e
  e_mat = np.matmul(w_vec, ci_mat)
#  print(f'e matrix is \n{e_mat}')

#a inverse
  ai_mat = np.linalg.pinv(a_mat)
#  print(f'a inverse is \n{ai_mat}')

#finding u (displacement)
  u_mat = np.matmul(ai_mat, e_mat)




  sv = solve_svd(a_mat)
  print(f'\nSingular values for A: \n{sv}')

  eigen_sv = solve_eigenvalues(sv)
  print(f'\nEigenvalues for A: \n{eigen_sv}')

  condition_value_a = find_largest_eigen(eigen_sv)/find_smallest_eigen(eigen_sv)
  print(f'\nCondition value for A: \n{condition_value_a}')

  svc = solve_svd(c_mat)
  print(f'\nSingular values for C: \n{svc}')

  eigen_svc = solve_eigenvalues(svc)
  print(f'\nEigenvalues for C: \n{eigen_svc}')
  
  condition_value_c = find_largest_eigen(eigen_svc)/find_smallest_eigen(eigen_svc)
  print(f'\nCondition value for A: \n{condition_value_c}')

  sva = solve_svd(at_mat)
  print(f'\nSingular values for A transpose: \n{sva}')

  eigen_sva = solve_eigenvalues(sva)
  print(f'\nEigenvalues for A transpose: \n{eigen_sva}')

  condition_value_a = find_largest_eigen(eigen_sva)/find_smallest_eigen(eigen_sva)
  print(f'\nCondition value for A: \n{condition_value_a}')



  print(f'\nu matrix is \n{u_mat}')

#stiffness matrix K
#  temp_mat = np.matmul(at_mat, c_mat)
#  k_mat = np.matmul(temp_mat, a_mat)
#  print(f'K matrix is \n{k_mat}')

#inverse stiffness matrix
#  ki_mat = np.linalg.pinv(k_mat)
#  print(f'K inverse is \n{ki_mat}')

#solving for displacement

  return u_mat


def solve_svd(mat):

  a_mat = mat
  U, s, VT = svd(a_mat)
#  print(f'Singular values are \n{s}')
  return s

def solve_eigenvalues(s):

  eigen = np.array(s) 
  for r in range(len(s)):
    eigen[r] = s[r]*2
 # print(f'Eigenvalues are \n{eigen}')
  
  return eigen

def find_largest_eigen(eigen):

  largest = eigen[0]
  for r in range(len(eigen)):
    if eigen[r] > largest:
      largest = eigen[r]
  return largest
  


def find_smallest_eigen(eigen):

  smallest = eigen[0]
  for r in range(len(eigen)):
    if eigen[r] < smallest:
      smallest = eigen[r]
  return smallest



#solve_displacement(a_mat_fixed_fixed(), f_vec)



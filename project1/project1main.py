#import numpy as np

#take in values for number of springs/masses, spring constants (c) and masses (m)

spring_count = 4
mass_count = 3

spring_constant = [2, 2, 2, 2]
mass_value = [3, 3, 3]


def fixed_fixed_function():

#creating A matrix
  a_mat = ([[0 for c in range(mass_count)] for r in range(spring_count)])

  for r in range(spring_count):
    for c in range(mass_count):
      if r == c:
        a_mat[r][c] = 1 
    #    print(f'r: {r} and c: {c}')
      elif r == c+1:
        a_mat[r][c] = -1
    #    print(f'r: {r} and c: {c}')

  print(f'A matrix is \n{a_mat}')

#creating identity matrix
  i_mat = ([[0 for c in range(spring_count)] for r in range(spring_count)])

  for r in range(spring_count):
    for c in range(spring_count):
      if r == c:
        i_mat[r][c] = spring_constant[r]

  print(f'I matrix is \n{i_mat}')


fixed_fixed_function()

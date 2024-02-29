# Hjelpefunksjoner

"""

Generelle matrise-funksjoner

"""

import numpy as np

def matrixMultiplication(A,B):
    return np.matmul(A,B)

def transposeMatrix(A):
    return np.matrix.transpose(A) # A^T

def inverseMatrix(A):
    return np.linalg.inv(A) # A^-1

# Mer spesifikke funksjoner:

def reflection1(V): # Speiler om x-aksen
    P = [[-1,0,0],
         [0,1,0],
         [0,0,1]]
    return matrixMultiplication(P, V)

def reflection2(V): # Speiler om y-aksen
    P = [[1,0,0],
         [0,-1,0],
         [0,0,1]]
    return matrixMultiplication(P, V)

def reflection3(V): # Speiler om z-aksen
    P = [[1,0,0],
         [0,1,0],
         [0,0,-1]]
    return matrixMultiplication(P, V)

def rotation1(V, theta): # Roterer om x-aksen med vinkel theta
    theta = theta * np.pi / 180
    R = [[1, 0, 0],
         [0, np.cos(theta), np.sin(theta)],
         [0, -np.sin(theta), np.cos(theta)]]
    return matrixMultiplication(R, V)

def rotation2(V, theta): # Roterer om y-aksen med vinkel theta
    theta = theta * np.pi / 180
    R = [[np.cos(theta), 0, -np.sin(theta)],
         [0, 1, 0],
         [np.sin(theta), 0, np.cos(theta)]]
    return matrixMultiplication(R, V)

def rotation3(V, theta): # Roterer om z-aksen med vinkel theta
    theta = theta * np.pi / 180
    R = [[np.cos(theta), np.sin(theta), 0],
         [-np.sin(theta), np.cos(theta), 0],
         [0, 0, 1]]
    return matrixMultiplication(R, V)

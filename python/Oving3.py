### Øving 3 ###

"""
Reliability analysis

Beregne nøyaktig posisjon av punkt med GNSS-data og beregne nøyaktigheten av disse målingene.
"""

import numpy as np

import Hjelpefunksjoner as h

"""
Hver rad i matrisen har følgende kolonner for samme punkt 'SL24':
[Nord (NTM19) [0], Øst (NTM19) [1], Høyde (NTM19)[2], xx [3], yy [4], zz [5], xy [6], xz [7], yz [8]]
"""

data = [[2281478.131, 90170.993, 50.136, 0.00003205, 0.00002769, 0.00023304, -0.00000176, 0.00002948, 0.00000301],
        [2281478.133, 90170.986, 50.139, 0.00002331, 0.00001760, 0.00014727, 0.00000073, 0.00002464, -0.00000191],
        [2281478.125, 90170.976, 50.130, 0.00003218, 0.00002182, 0.00029197, 0.00000295, 0.00004603, 0.00001266],
        [2281478.127, 90170.970, 50.104, 0.00003853, 0.00001595, 0.00013528, 0.00000508, 0.00002023, 0.00000750]] # [m] og [m^2]

C_xyz_1 = np.array([[data[0][0], data[0][6], data[0][7]], [data[0][6], data[0][4], data[0][8]], [data[0][7], data[0][8], data[0][5]]]) # [m^2]
C_xyz_2 = np.array([[data[1][0], data[1][6], data[1][7]], [data[1][6], data[1][4], data[1][8]], [data[1][7], data[1][8], data[1][5]]]) # [m^2]
C_xyz_3 = np.array([[data[2][0], data[2][6], data[2][7]], [data[2][6], data[2][4], data[2][8]], [data[2][7], data[2][8], data[2][5]]]) # [m^2]
C_xyz_4 = np.array([[data[3][0], data[3][6], data[3][7]], [data[3][6], data[3][4], data[3][8]], [data[3][7], data[3][8], data[3][5]]]) # [m^2]

lat, lon = 69.495870257, 19.248665575 # [deg]

### Oppgave 1 ###

"""
print("De fire symmetriske kovarians-matrisene:\n")
print(C_xyz_1)
print()
print(C_xyz_2)
print()
print(C_xyz_3)
print()
print(C_xyz_4)
#"""

### Oppgave 2 ###



### Oppgave 3 ###

n = 12
e = 3
f = n - e

# E(l) = L0 - li

A = np.array([[-1, 0, 0],
              [0, -1, 0],
              [0, 0, -1],
              [-1, 0, 0],
              [0, -1, 0],
              [0, 0, -1],
              [-1, 0, 0],
              [0, -1, 0],
              [0, 0, -1],
              [-1, 0, 0],
              [0, -1, 0],
              [0, 0, -1]])

"""
p0 = 1

p0*s0**2 = pn*sn**2 => pn = s0**2/sn**2

Velger at n = 1 har p0 = 1.
"""

P = np.eye(n)

f = np.array([-data[0][0],
              -data[0][1],
              -data[0][2],
              -data[1][0],
              -data[1][1],
              -data[1][2],
              -data[2][0],
              -data[2][1],
              -data[2][2],
              -data[3][0],
              -data[3][1],
              -data[3][2],])

### Oppgave 4 ###

x0 = np.array([2281480, 90170, 50])

for i in range(n):
    f[i] += x0[i % 3]

x_hat  = h.inverseMatrix(h.transposeMatrix(A) @ P @ A) @ h.transposeMatrix(A) @ P @ f

print(x_hat)
svar = x0 + x_hat
for i in range(len(svar)):
    print(str(svar[i]) + " : " + str(x0[i]))

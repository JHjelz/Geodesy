### Øving 3 ###

"""
Reliability analysis

Beregne nøyaktig posisjon av punkt med GNSS-data og beregne nøyaktigheten av disse målingene.
"""

# Importerer bibliotek:

import numpy as np

import Hjelpefunksjoner as h

# Data:

"""
Hver rad i matrisen har følgende kolonner for samme punkt 'SL24':
[Nord (NTM19) [0], Øst (NTM19) [1], Høyde (NN2000)[2], xx [3], yy [4], zz [5], xy [6], xz [7], yz [8]]
"""

data = [[2281478.131, 90170.993, 50.136, 0.00003205, 0.00002769, 0.00023304, -0.00000176, 0.00002948, 0.00000301],
        [2281478.133, 90170.986, 50.139, 0.00002331, 0.00001760, 0.00014727, 0.00000073, 0.00002464, -0.00000191],
        [2281478.125, 90170.976, 50.130, 0.00003218, 0.00002182, 0.00029197, 0.00000295, 0.00004603, 0.00001266],
        [2281478.127, 90170.970, 50.104, 0.00003853, 0.00001595, 0.00013528, 0.00000508, 0.00002023, 0.00000750]] # [m] og [m^2]

lat, lon = 69.495870257, 19.248665575 # [deg]

x0 = np.array([2281480, 90170, 50])

### Oppgave 1 ###

C_xyz_1 = np.array([[data[0][3], data[0][6], data[0][7]], [data[0][6], data[0][4], data[0][8]], [data[0][7], data[0][8], data[0][5]]]) # [m^2]
C_xyz_2 = np.array([[data[1][3], data[1][6], data[1][7]], [data[1][6], data[1][4], data[1][8]], [data[1][7], data[1][8], data[1][5]]]) # [m^2]
C_xyz_3 = np.array([[data[2][3], data[2][6], data[2][7]], [data[2][6], data[2][4], data[2][8]], [data[2][7], data[2][8], data[2][5]]]) # [m^2]
C_xyz_4 = np.array([[data[3][3], data[3][6], data[3][7]], [data[3][6], data[3][4], data[3][8]], [data[3][7], data[3][8], data[3][5]]]) # [m^2]

C_xyz = np.array([C_xyz_1, C_xyz_2, C_xyz_3, C_xyz_4])

"""
print("De fire symmetriske varians-kovarians-matrisene:\n")
for i in range(len(C_xyz)):
    print(C_xyz[i])
    print()
#"""

### Oppgave 2 ###

C_neh = np.zeros((4, 3, 3))

for i in range(len(C_xyz)):
    C_neh[i] = h.reflection2(h.rotation2(h.rotation3(C_xyz[i], lon - 180), lat - 90))

"""
print("Varians-kovarians-matrisene gjort om til (N, E, H):\n")
for i in range(len(C_neh)):
    print(C_neh[i])
    print()
#"""

### Oppgave 3 ###

n = 12
e = 3

# E(li) = li + vi = x <=> vi = x - li, der x = l0

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

Velger at n = 1 har p1 = p0 = 1.
"""

s0 = data[0][3:6]

P = np.zeros((n, n))
for i in range(n):
    P[i][i] = s0[i % 3] / data[int((i - i % 3)/3)][i % 3 + 3]

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

"""
print(A)
print()
print(P)
print()
print(f)
#"""

### Oppgave 4 ###

deltaX = np.zeros(3)
iterate = True
iteration = 0

for i in range(n):
    f[i] += x0[i % 3]

while iterate:

    x_hat  = h.inverseMatrix(h.transposeMatrix(A) @ P @ A) @ h.transposeMatrix(A) @ P @ f
    iteration += 1

    deltaX += x_hat

    for i in range(n):
        f[i] += x_hat[i % 3]
    
    count = 0
    for val in x_hat:
        if np.abs(val) < 10**(-9): #0.0005:
            count += 1
    if count == len(x_hat):
        iterate = False

svar = x0 + deltaX

"""
print("Antall iterasjoner: " + str(iteration) + "\n")

print("Estimert verdi LSM\tOppgitt verdi\tEstimert verdi gjennomsnitt:\n")

for i in range(len(svar)):    
    
    val = 0
    for j in range(len(data)): # Beregner gjennomsnitt
        val += data[j][i]
    
    if i == 1:
        s = "  "
    elif i == 2:
        s = "     "
    else:
        s = ""

    print(s + format(svar[i], '.7f') + "\t\t" + s + str(x0[i]) + "\t\t" + s + format(val / len(data), '.7f'))

print()
#"""

for i in range(n):
    f[i] -= x0[i % 3]

v = A @ svar - f

sigma0 = np.sqrt(h.transposeMatrix(v) @ P @ v / (n - e))

Q = h.inverseMatrix(h.transposeMatrix(A) @ P @ A)

C = sigma0**2 * Q

"""
print("Standardavvik:")
print("Nord: ", np.sqrt(C[0][0]), "m")
print("Øst:  ", np.sqrt(C[1][1]), "m")
print("Høyde:", np.sqrt(C[2][2]), "m")
print()
#"""

### Oppgave 5 ###

df = n - e # 12 - 3 = 9
alpha = 0.05
t = 2.262

"""
for i in range(len(svar)):
    val = np.abs((svar[i] - x0[i]) / np.sqrt(C[i][i]))
    print(val - t)
    print((val - t) > 0)
#"""

### Oppgave 6 ###

"""
for i in range(len(v)):
    s1, s2 = "", ""
    if v[i] > 0:
        s2 = " "
    if i < 9:
        s1 = " "
    print("v" + str(i+1) + s1 + " = " + s2 + str(v[i]))
#"""

### Oppgave 7 ###


### Øving 3 ###

"""
Reliability analysis

Beregne nøyaktig posisjon av punkt med GNSS-data og beregne nøyaktigheten av disse målingene.
"""

# Importerer bibliotek:

import numpy as np
from scipy.stats import t as studentT

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

# Funksjoner:

def Normalligningene(A, P):
    """
    Beregner normalligningene N
    """
    return h.inverseMatrix(h.transposeMatrix(A) @ P @ A)

def LSM(A, P, f):
    """
    Utfører LSM-prosedyren på A- og P-matrisene og f-vektoren
    """
    return Normalligningene(A, P) @ h.transposeMatrix(A) @ P @ f

def Feilligningene(A, x, f):
    """
    Beregner feilleddene v fra LSM-prosedyren
    """
    return A @ x - f

def s0(v, P, n, e):
    """
    Beregner standardavviket for enhetsvekten
    """
    return np.sqrt(h.transposeMatrix(v) @ P @ v / (n - e))

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
    rho = np.pi / 180
    Rlat, Rlon = lat * rho, lon * rho
    R = np.array([[-np.sin(Rlat)*np.cos(Rlon), -np.sin(Rlat)*np.sin(Rlon), np.cos(Rlat)],
                  [-np.sin(Rlon), np.cos(Rlon), 0],
                  [np.cos(Rlat)*np.cos(Rlon), np.cos(Rlat)*np.sin(Rlon), np.sin(Rlat)]])
    C_neh[i] = R @ C_xyz[i] @ h.transposeMatrix(R)

"""
print("Varians-kovarians-matrisene gjort om til (N, E, H):\n")
for i in range(len(C_neh)):
    print(C_neh[i])
    print()
#"""

### Oppgave 3 ###

n = 12
e = 3

# E(li) = li + vi = x <=> vi = x - li, der x er vår ukjente og li er differansen mellom observasjoner og antatt verdi = l0

A = np.array([[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1],
              [1, 0, 0],
              [0, 1, 0],
              [0, 0, 1],
              [1, 0, 0],
              [0, 1, 0],
              [0, 0, 1],
              [1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

"""
Vekter er relative, så det har ikke så mye å si hvilke verdier
en bruker så lenge de samsvarer med målingene som er gjort

P = C^(-1)
"""

P = np.zeros((n, n))

for i in range(len(C_neh)):
    P[i*3:(i+1)*3, i*3:(i+1)*3] = h.inverseMatrix(C_neh[i])

f = np.array([data[0][0]-x0[0],
              data[0][1]-x0[1],
              data[0][2]-x0[2],
              data[1][0]-x0[0],
              data[1][1]-x0[1],
              data[1][2]-x0[2],
              data[2][0]-x0[0],
              data[2][1]-x0[1],
              data[2][2]-x0[2],
              data[3][0]-x0[0],
              data[3][1]-x0[1],
              data[3][2]-x0[2]])

"""
print("Designmatrisen A:")
print(A)
print("\nVektsmatrisen P:")
print(P)
print("\nObservasjonsvektoren f:")
print(f)
print()
#"""

### Oppgave 4 ###

x_hat = LSM(A, P, f)

svar = x0 + x_hat

"""
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

v = Feilligningene(A, x_hat, f)
sigma0 = s0(v, P, n, e)

Q = Normalligningene(A, P)

C = sigma0**2 * Q

"""
print("Standardavvik av enhetsvekten:")
print(sigma0)
print("\nKofaktor-matrisen:")
print(Q)
print("\nVarians-kovarians-matrisen:")
print(C)
print()
print("Standardavvik:")
print("Nord: ", np.sqrt(C[0][0]), "m\t=\t", np.sqrt(C[0][0]) * 10**3, "mm")
print("Øst:  ", np.sqrt(C[1][1]), "m\t=\t", np.sqrt(C[1][1]) * 10**3, "mm")
print("Høyde:", np.sqrt(C[2][2]), "m\t=\t", np.sqrt(C[2][2]) * 10**3, "mm")
print()
#"""

### Oppgave 5 ###

feil = set()
k = len(feil)

while True:
    A_test = np.zeros((n - k, e))
    P_test = np.zeros((n - k, n - k))
    f_test = np.zeros(n - k)

    next1 = 0

    for i in range(n):
        if i in feil:
            continue
        
        next2 = 0
        for j in range(e):
            A_test[next1][next2] = A[i][j]
            next2 += 1
        
        next3 = 0
        for j in range(n):
            if j in feil:
                continue
            P_test[next1][next3] = P[i][j]
            next3 += 1
        
        f_test[next1] = f[i]
        
        next1 += 1

    nabla = np.array([])
    sigma_nabla = np.array([])

    for i in range(n - k):
        col = np.zeros(n - k)
        col[i] = 1
        
        A_new = np.hstack((A_test, col.reshape(-1, 1)))
        x_hat_test = LSM(A_new, P_test, f_test)

        Q_test = Normalligningene(A_new, P_test)
        v_test = Feilligningene(A_new, x_hat_test, f_test)
        sigma0_test = s0(v_test, P_test, n - k, e + 1)
        C_test = sigma0_test**2 * Q_test

        nabla = np.hstack((nabla, x_hat_test[-1]))
        sigma_nabla = np.hstack((sigma_nabla, np.sqrt(C_test[-1][-1])))

    t = np.abs(nabla / sigma_nabla)

    alpha_tot = 0.05 # 5%
    alpha = 1 - (1 - alpha_tot)**(1/(n - k))
    df = n - k - e - 1
    T = studentT.ppf(1 - alpha / 2, df)

    newFeil = []

    for i in range(len(t)):
        if t[i] > T and i not in feil:
            newFeil.append(i)

    if not newFeil:
        break

    feil.update(newFeil)
    k = len(feil)

"""
print("Endelig T-grense:", T, "\n")
for i in range(len(nabla)):
    print("Måling " + str(i+1) + " - nabla:\t", format(nabla[i], '.5f'))
    print("Måling " + str(i+1) + " - sigma:\t", format(sigma_nabla[i], '.5f'))
    print("Måling " + str(i+1) + " - t:\t\t", format(t[i], '.5f'))
    print()
#"""

### Oppgave 6 ###

I = np.eye(n - k)
R = I - A @ h.inverseMatrix(h.transposeMatrix(A) @ P @ A) @ h.transposeMatrix(A) @ P

"""
for i in range(n-k):
    print(R[i][i])
#"""

### Oppgave 7 ###

LRE = []

T = studentT.ppf(1 - 0.05 / 2, 8)

for i in range(n - k):
    low = nabla[i] - sigma_nabla[i] * T
    high = nabla[i] + sigma_nabla[i] * T
    if np.abs(low) > np.abs(high):
        LRE.append(low)
    else:
        LRE.append(high)

"""
for i in range(n-k):
    print(LRE[i])
#"""

### Oppgave 8 ###

nabla_x = np.zeros((n, e))

for i in range(n - k):
    nabla_max = np.zeros(n - k)
    nabla_max[i] = LRE[i]
    x = h.inverseMatrix(h.transposeMatrix(A) @ P @ A) @ h.transposeMatrix(A) @ P @ nabla_max
    nabla_x[i] = x * 10**3

# Målfunksjoner:

# 1D er kun høyden h = z = nabla_x[i][2] for alle i

# 2D:
målfunksjon1 = lambda x, y: np.sqrt(x**2 + y**2)

# 3D:
målfunksjon2 = lambda x, y, z: np.sqrt(x**2 + y**2 + z**2)

for i in range(n - k):
    x, y, z = nabla_x[i][0], nabla_x[i][1], nabla_x[i][2]
    nabla_x[i] = np.array([z, målfunksjon1(x, y), målfunksjon2(x, y, z)])

max1, max2, max3 = 0, 0, 0

for i in range(n - k):
    if np.abs(nabla_x[i][0]) > max1:
        max1 = np.abs(nabla_x[i][0])
    if np.abs(nabla_x[i][1]) > max2:
        max2 = np.abs(nabla_x[i][1])
    if np.abs(nabla_x[i][2]) > max3:
        max3 = np.abs(nabla_x[i][2])

#"""
print(max1)
print(max2)
print(max3)
#"""
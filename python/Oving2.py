### Øving 2 ###

"""
GPS Computations

Beregne posisjon til satellitter og mottakere ved hjelp av GNSS-signaler fra SV - satellite vehicle.
"""

# Importerer biblioteker:

import numpy as np

import Hjelpefunksjoner as h
import KoordinatKonvertering as k
from Appendix_2 import App_2 as data

# Funksjoner:

rho = lambda xj, yj, zj, x0, y0, z0: np.sqrt((xj - x0)**2 + (yj - y0)**2 + (zj - z0)**2)

# Data:

t = 129600 # Epoch [sek]
GM = 3.986005 * 10**14 # [m^3/s^2]
omega_e = 7.2921151467 * 10**(-5) # [rad/s]

# Antatt, geodetisk GPS-posisjon for mottakeren:
phi_r0 = 47.1 # [deg]
lambda_r0 = 15.5 # [deg]
h_r0 = 400 # [m]

c = 299792458 # Lyshastighet [m/s]

# Program:

### Oppgave 1 & 2 ###

koordinater = np.zeros((4, 3))
koordinater2 = np.zeros((4, 3))

j = 0

for key in data:
    t_k = t - data[key][0]
    if t_k > 302400:
        t_k -= 604800
    elif t_k < -302400:
        t_k += 604800
    
    # Dataverdiene 7 - 15 = 0 når en ser bort fra korreksjonsleddene
    
    M_k = data[key][3] + (np.sqrt(GM / data[key][1]**6) + data[key][7]) * t_k
    M_k2 = data[key][3] + np.sqrt(GM / data[key][1]**6) * t_k

    E0 = M_k
    E02 = M_k2
    for i in range(3):
        E_k = E0 + (M_k - E0 + data[key][2] * np.sin(E0)) / (1 - data[key][2] * np.cos(E0))
        E_k2 = E02 + (M_k2 - E02 + data[key][2] * np.sin(E02)) / (1 - data[key][2] * np.cos(E02))
        E0 = E_k
        E02 = E_k2
    
    f_k = 2 * np.arctan(np.sqrt((1 + data[key][2]) / (1 - data[key][2])) * np.tan(E_k / 2))
    f_k2 = 2 * np.arctan(np.sqrt((1 + data[key][2]) / (1 - data[key][2])) * np.tan(E_k2 / 2))

    u_k = data[key][4] + f_k + data[key][10] * np.cos(2 * (data[key][4] + f_k)) + data[key][11] * np.sin(2 * (data[key][4] + f_k))
    r_k = data[key][1]**2 * (1 - data[key][2] * np.cos(E_k)) + data[key][12] * np.cos(2 * (data[key][4] + f_k)) + data[key][13] * np.sin(2 * (data[key][4] + f_k))
    i_k = data[key][5] + data[key][8] * t_k + data[key][14] * np.cos(2 * (data[key][4] + f_k)) + data[key][15] * np.sin(2 * (data[key][4] + f_k))
    lambda_K = data[key][6] + (data[key][9] - omega_e) * t_k - omega_e * data[key][0]
    
    # [X^k, Y^k, Z^k] = R_3(-lambda_k) * R_1(-i_k) * R_3(-u_k) * [r_k, 0, 0]

    koor = np.array([[np.cos(-lambda_K), np.sin(-lambda_K), 0],
                     [-np.sin(-lambda_K), np.cos(-lambda_K), 0],
                     [0, 0, 1]]) \
                        @ np.array([[1, 0, 0],
                                    [0, np.cos(-i_k), np.sin(-i_k)],
                                    [0, -np.sin(-i_k), np.cos(-i_k)]]) \
                                        @ np.array([[np.cos(-u_k), np.sin(-u_k), 0],
                                                    [-np.sin(-u_k), np.cos(-u_k), 0],
                                                    [0, 0, 1]]) \
                                                        @ np.array([r_k, 0, 0])

    koordinater[j] = koor
    
    u_k2 = data[key][4] + f_k2
    r_k2 = data[key][1]**2 * (1 - data[key][2] * np.cos(E_k))
    i_k2 = data[key][5]
    lambda_K2 = data[key][6] - omega_e * t_k - omega_e * data[key][0]

    koor2 = np.array([[np.cos(-lambda_K2), np.sin(-lambda_K2), 0],
                      [-np.sin(-lambda_K2), np.cos(-lambda_K2), 0],
                      [0, 0, 1]]) \
                        @ np.array([[1, 0, 0],
                                    [0, np.cos(-i_k2), np.sin(-i_k2)],
                                    [0, -np.sin(-i_k2), np.cos(-i_k2)]]) \
                                        @ np.array([[np.cos(-u_k2), np.sin(-u_k2), 0],
                                                    [-np.sin(-u_k2), np.cos(-u_k2), 0],
                                                    [0, 0, 1]]) \
                                                        @ np.array([r_k2, 0, 0])

    koordinater2[j] = koor2

    j += 1

#"""
print("Estimerte koordinater for de fire satellittene med korreksjoner:")
print(koordinater)
print()
print("Estimerte koordinater for de fire satellittene uten korreksjoner:")
print(koordinater2)
print()
print("Differanse mellom med og uten korreksjoner:")
print(koordinater - koordinater2)
print()
#"""

### Oppgave 3 ###

mottaker_0 = np.array(k.LatLonToCartesian(phi_r0, lambda_r0, h_r0, "WGS84"))

#"""
print("De antatte mottaker-koordinatene i kartesiske koordinater:")
print(mottaker_0)
print()
#"""

### Oppgave 4 ###

A = np.zeros((4, 4))
L = np.zeros(4)
dTr = 0

iterasjon = True
iterasjoner = 0

while iterasjon:

    xr0, yr0, zr0 = mottaker_0[0], mottaker_0[1], mottaker_0[2]

    j = 0

    for key in data:
        P = data[key][16]
        Xj, Yj, Zj = koordinater[j][0], koordinater[j][1], koordinater[j][2]
        
        A[j] = np.array([-(Xj - xr0) / rho(Xj, Yj, Zj, xr0, yr0, zr0),
                         -(Yj - yr0) / rho(Xj, Yj, Zj, xr0, yr0, zr0),
                         -(Zj - zr0) / rho(Xj, Yj, Zj, xr0, yr0, zr0),
                         -c])
        
        L[j] = P - rho(Xj, Yj, Zj, xr0, yr0, zr0)

        j += 1
    
    x_hat = h.inverseMatrix(A) @ L
    iterasjoner += 1
    
    mottaker_0 += x_hat[:3]
    dTr = x_hat[3]

    teller = 0
    for val in x_hat:
        if np.abs(val) < 0.0005:
            teller += 1
    if teller == 4:
        iterasjon = False

#"""
print("Antall iterasjoner brukt i LSM: " + str(iterasjoner))
print("Estimert posisjon / kartesiske koordinat for mottakeren:")
print(mottaker_0)
print()
#"""

### Oppgave 5 ###

Q = h.inverseMatrix(h.transposeMatrix(A) @ A)

diagonal3 = 0
for i in range(3):
    diagonal3 += Q[i][i]

PDOP = np.sqrt(diagonal3)

#"""
print("Kofaktor-matrisen Q:")
print(Q)
print("Positional dilution of precision PDOP:")
print(PDOP)
print()
#"""

### Oppgave 6 ###

mottaker = k.CartesianToLatLon(mottaker_0[0], mottaker_0[1], mottaker_0[2])

#"""
print("Mottakerens posisjon i bredde- og lengdegrader:")
print(mottaker)
print()
#"""

### Oppgave 7 ###

#"""
print("Estimert klokke-feil i mottakeren:")
print(dTr * 10**9)
#"""
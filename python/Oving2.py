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
        E = E0 + (M_k - E0 + data[key][2] * np.sin(E0)) / (1 - data[key][2] * np.cos(E0))
        E2 = E02 + (M_k2 - E02 + data[key][2] * np.sin(E02)) / (1 - data[key][2] * np.cos(E02))
        E0 = E
        E02 = E2
    E_k = E
    E_k2 = E2
    
    f_k = 2 * np.arctan(np.sqrt((1 + data[key][2]) / (1 - data[key][2])) * np.tan(E_k / 2))
    f_k2 = 2 * np.arctan(np.sqrt((1 + data[key][2]) / (1 - data[key][2])) * np.tan(E_k2 / 2))

    u_k = data[key][4] + f_k + data[key][10] * np.cos(2 * (data[key][4] + f_k)) + data[key][11] * np.sin(2 * (data[key][4] + f_k))
    r_k = data[key][1]**2 * (1 - data[key][2] * np.cos(E_k)) + data[key][12] * np.cos(2  * (data[key][4] + f_k)) + data[key][13] * np.sin(2 * (data[key][4] + f_k))
    i_k = data[key][5] + data[key][8] * t_k + data[key][14] * np.cos(2 * (data[key][4] + f_k)) + data[key][15] * np.sin(2 * (data[key][4] + f_k))
    lambda_K = data[key][6] + (data[key][9] - omega_e) * t_k - omega_e * data[key][0]

    koor = h.rotation3(h.rotation1(h.rotation3(np.array([r_k, 0, 0]), -u_k), -i_k), -lambda_K)

    koordinater[j] = koor
    
    u_k2 = data[key][4] + f_k2
    r_k2 = data[key][1]**2 * (1 - data[key][2] * np.cos(E_k))
    i_k2 = data[key][5]
    lambda_K2 = data[key][6] - omega_e * t_k - omega_e * data[key][0]

    koor2 = h.rotation3(h.rotation1(h.rotation3(np.array([r_k2, 0, 0]), -u_k2), -i_k2), -lambda_K2)

    koordinater2[j] = koor2

    j+= 1

"""
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

"""
print("De antatte mottaker-koordinatene i kartesiske koordinater:")
print(mottaker_0)
print()
#"""

### Oppgave 4 ###

A = np.zeros((4, 4))
L = np.zeros(4)

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

Q = h.inverseMatrix(h.transposeMatrix(A) @ A)


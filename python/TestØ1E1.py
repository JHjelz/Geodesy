### Øving 1 ###

"""
GNSS Baselines

Fra GNSS baseline observasjoner til observasjoner i UTM kartprojeksjon og høyder i NN2000
"""

# Importerer biblioteker:

import numpy  as np

import Hjelpefunksjoner as h
import KoordinatKonvertering as koorkon

# Funksjoner:

# Generelle

def sumListe(a, b):
    """
    Legger sammen verdier av to lister hvis de er like lange
    """
    if len(a) != len(b):
        return
    ans = []
    for i in range(len(a)):
        ans.append(a[i] + b[i])
    return ans

def deltaListe(a, b):
    """
    Finner differansen mellom to lister hvis de er like lange
    """
    if len(a) != len(b):
        return
    ans = []
    for i in range(len(a)):
        ans.append(b[i] - a[i])
    return ans

def gradTilrad(v):
    return v * np.pi / 180

def radTilgrad(v):
    return v * 180 / np.pi

### Oppgave 1 ###

# a)

def CTRS_til_LG(delta, phi, Lambda):
    """
    Transformerer koordinater fra CTRS til LG
    """
    return h.reflection2(h.rotation2(h.rotation3(delta, Lambda - 180), phi - 90))

def slopeDistance(bl):
    """
    Beregner distanse i 3D-plan
    """
    return np.sqrt(bl[0]**2 + bl[1]**2 + bl[2]**2)

def horizontalDistance(bl):
    """
    Beregner distanse i 2D-plan
    """
    return np.sqrt(bl[0]**2 + bl[1]**2)

# b)

def azimuth(bl):
    """
    Beregner azimuth for en baseline
    """
    if bl[0] > 0:
        if bl[1] > 0:
            return np.arctan(bl[1] / bl[0])
        else:
            return np.arctan(bl[1] / bl[0]) + 2 * np.pi
    else:
        if bl[1] > 0:
            return np.arctan(bl[1] / bl[0]) + np.pi
        else:
            return np.arctan(bl[1] / bl[0]) + np.pi

def w(B):
    """
    Finner konstanten W for en gitt breddegrad B
    """
    B = gradTilrad(B)
    return np.sqrt(1 - e**2 * (np.sin(B))**2)

def rM(B):
    """
    Beregner krumningsradius langs 'prime vertical'
    """
    return (a * (1 - e**2))/(w(B))**3

def rN(B):
    """
    Beregner krumningsradius langs meridianen
    """
    return a/w(B)

def r(B, azi):
    """
    Beregner radius i gitt breddegrad B i retning av azimuth azi
    """
    return (rN(B) * rM(B)) / (rN(B) * (np.cos(azi))**2 + rM(B) * (np.sin(azi))**2)

def zenithAngle(bl):
    """
    Beregner zenit-vinkel for en gitt baseline
    """
    angle = np.arctan(horizontalDistance(bl) / bl[2])
    if angle < 0:
        return angle + np.pi
    return angle

def horizontalDistance_2(bl):
    """
    Beregner distanse i 2D-planet
    """
    return slopeDistance(bl) * np.sin(zenithAngle(bl))

def earthAngle(bl, R, h):
    """
    Beregner vinkelen gamma ved jordsentrum mellom to punkt i gitt baseline
    """
    return np.arctan(horizontalDistance_2(bl) / (R + h + bl[2]))

def correctedHorizontalDistance(bl, B, h, east):
    """
    Beregner korrigert distanse i kartprojeksjonsplanet
    """
    R = r(B, azimuth(bl[:2]))
    angle = earthAngle(bl, R, h)
    s0 = R * angle
    for i in range(len(east)):
        east[i] = (east[i] - 500000) / 0.9996
    return s0 + s0 / (6 * R**2) * (east[0]**2 + east[0] * east[1] + east[1]**2) - 0.0004 * s0

# c)

def correctedAzimuth(bl, B, L, FraUTM, TilUTM):
    """
    Beregner korrigert azimuth i kartprojeksjonsplanet
    """
    korC = correctionC(B, L)
    korDelta = correctionDelta(bl, B, FraUTM, TilUTM)
    return azimuth(bl) * 180 / np.pi - np.abs(korC) - np.abs(korDelta[0])

def correctionC(B, L):
    """
    Beregner korreksjonsleddet C for korrigeringen av azimuth
    """
    B = gradTilrad(B)
    l = gradTilrad(L - ((L // 6) * 6 + 3))

    NU = nu(B)

    f1 = l * np.sin(B)
    f2 = (l**3 / 3) * np.sin(B) * np.cos(B)**2 * (1 + 3 * NU**2 + 2 * NU**4)
    f3 = (l**5 / 15) * np.sin(B) * (np.cos(B))**4 * (2 - (np.tan(B))**2)
    return radTilgrad(f1 + f2 + f3)

def correctionDelta(bl, B, p1, p2):
    """
    Beregner korreksjonsleddet delta for korrigeringen av azimuth
    """
    R = r(B, azimuth(bl))
    AB = radTilgrad(1 / (6 * R**2) * (2 * (p1[1] - 500000)/0.9996 + (p2[1] - 500000)/0.9996) * (p2[0] - p1[0])) # Korrigere for UTM-fellene
    BA = radTilgrad(1 / (6 * R**2) * (2 * (p2[1] - 500000)/0.9996 + (p1[1] - 500000)/0.9996) * (p1[0] - p2[0]))
    return AB, BA

def nu(B):
    """
    Hjelpeledd for korreksjonen i C
    """
    return np.sqrt(e**4 * (np.cos(B))**2 / (1 - e**4))

# d)

def correctedHeightDifference(bl, B, NStart, NSlutt):
    """
    Beregner korrigert høydedifferanse for en baseline
    """
    return bl[2] + correctedEarthCurv(bl, B) - (NSlutt - NStart)

def correctedEarthCurv(bl, B):
    """
    Beregner korreksjoner for høydeforskjell grunnet jordas krumning
    """
    return (slopeDistance(bl) * np.sin(zenithAngle(bl)))**2 / (2 * r(B, azimuth(bl)))

# e)

# ...

### Oppgave 2 ###

# a)

# ...

# b)

def partialDerivatives(bl):
    """
    Beregner partiellderiverte verdier for gitt baseline
    Returnerer transformasjonsmatrise som brukes for å finne standardavvik for utledede observasjoner
    
    [x, y, z] = lokale koordinat lik [dx, dy, dz] i baseline 'bl'.
    
    A = arctan(y / x)
    => dA/dx = -y / (x^2 * ((y/x)^2 + 1)) og dA/dy = 1 / (x * (y/x)^2 + 1)
    => dA/dx = -sin(azi) * rho / S_h (I) og dA/dy = cos(azi) * rho / S_h (II)

    S_h = sqrt(x^2 + y^2)
    => dS_h/dx = x / sqrt(x^2 + y^2) og dS_h/dy = y / sqrt(x^2 + y^2)
    =A dS_h/dx = cos(azi) (III) og dS_h/dy = -sin(azi) (IV)

    H = z = dH/dz = 1 (V)
    """

    azi = azimuth(bl)
    S_h = horizontalDistance(bl)
    rho = 200 / np.pi # gon / rad

    ans = np.zeros(3 * 3).reshape(3, 3)

    ans[0][0] = -np.sin(azi) * rho / S_h #-bl[1] / (bl[0]**2 * ((bl[1]/bl[0])**2 + 1)) # (I)
    ans[0][1] = np.cos(azi) * rho / S_h #1 / (bl[0] * ((bl[1]/bl[0])**2 + 1)) # (II)
    ans[1][0] = np.cos(azi) #bl[0] / np.sqrt(bl[0]**2 + bl[1]**2) # (III)
    ans[1][1] = -np.sin(azi) #bl[1] / np.sqrt(bl[0]**2 + bl[1]**2) # (IV)
    ans[2][2] = 1 # (V)

    return ans

# c)

# ...

### Oppgave 3 ###

# a)

# ...

# b)

# ...

# Data:

"""
Fasit:

EUREF89: sted = [N (UTM), E (UTM), H (NN2000), h (ellipsoide)]
"""

EUREF_S = [7032158.662, 579058.469, 422.853, 462.071, 39.218]
EUREF_L = [7036209.878, 571150.913, 69.634, 108.875, 39.241]
EUREF_H = [7032450.406, 569796.193, 71.006, 110.462, 39.456]

# Jeg velger baseline 3 og 4: Moholt -> TP342 og ST46 -> Moholt

# sted = [X, Y, Z]

S = koorkon.UTMToCartesian(EUREF_S[0], EUREF_S[1], EUREF_S[2], "32N")
L = koorkon.UTMToCartesian(EUREF_L[0], EUREF_L[1], EUREF_L[2], "32N")
H = koorkon.UTMToCartesian(EUREF_H[0], EUREF_H[1], EUREF_H[2], "32N")

# baseline = [dX, dY, dZ]

LH = [3541.713, -811.048, -1666.943]
SH = [1084.518, -9212.999, -87.373]
SL = [-2457.183, -8401.958, 1579.590]

# Bredde- og lengdegrad for aktuelle punkt: sted = [B, L, h] (EUREF89 med ellipsoidisk høyde)

bl_S = koorkon.CartesianToLatLon(S[0], S[1], S[2])
bl_L = koorkon.CartesianToLatLon(L[0], L[1], L[2])
bl_H = koorkon.CartesianToLatLon(H[0], H[1], H[2])

"""
# Konvertert med skTrans:

bl_moholt   = [63.408951187, 10.431110376, 170.843]
bl_st46     = [63.431668740, 10.434433578, 151.865]
bl_tp342    = [63.427301973, 10.380349670,  44.613]
bl_havstein = [63.406818270, 10.377168783, 171.020]
"""

# UTM-koordinater for aktuelle punkt: sted = [N, E] (EUREF89)

UTM_S = koorkon.LatLonToUTM(bl_S[0], bl_S[1])
UTM_L = koorkon.LatLonToUTM(bl_L[0], bl_L[1])
UTM_H = koorkon.LatLonToUTM(bl_H[0], bl_H[1])

# EUREF-data

a = 6378137 # [m]
b = 6356752.3141 # [m]
e = np.sqrt((a**2 - b**2) / a**2)

# Program:

### Oppgave 1 ###

# a)

def O1a():
    # Konverterer til lokalt system:

    HLG = CTRS_til_LG(SH, bl_S[0], bl_S[1])
    HLG2 = CTRS_til_LG(LH, bl_L[0], bl_L[1])

    print("\n### Oppgave 1 ###\n\na)\n")
    print("Lokale koordinat for Holtermannsveien:\t" + str(HLG))
    print("Lokale koordinat for Holtermannsveien:\t" + str(HLG2))

    """
    En enkel kontrollsjekk for å sikre riktig resultat er å beregne
    distansen til baselinen i CTRS og LG - disse skal være like lange!
    """

    len1, len2 = slopeDistance(SH), slopeDistance(LH)
    len3, len4 = slopeDistance(HLG), slopeDistance(HLG2)

    print("\nKontrollsjekk:\n")
    print("Differanse S -> H:\t" + str(len1) + " m = " + str(len3) + " m,\tdiff = " + str(np.abs(len1 - len3)) + " m")
    print("Differanse L -> H:\t" + str(len2) + " m = " + str(len4) + " m,\tdiff = " + str(np.abs(len2 - len4)) + " m")

# b)

def O1b():
    #Konverterer til lokalt system:

    HLG = CTRS_til_LG(SH, bl_S[0], bl_S[1])
    HLG2 = CTRS_til_LG(LH, bl_L[0], bl_L[1])

    # Beregner korrigerte distanser:

    corrDistH = correctedHorizontalDistance(HLG, bl_S[0], bl_S[2], [UTM_S[1], UTM_H[1]])
    corrDistH2 =  correctedHorizontalDistance(HLG2, bl_L[0], bl_L[2], [UTM_L[1], UTM_H[1]])

    print("\n### Oppgave 1 ###\n\nb)\n")
    print("Korrigert distanse fra Solemsvåttan til Holtermannsveien:\t" + str(corrDistH) + " m")
    print("Korrigert distanse fra Ladehammeren til Holtermannsveien:\t\t" + str(corrDistH2) + " m")

    """
    Utfører kontroll ved å sammenligne de korrigerte distansene over med vanlig horisontal
    distanse målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = horizontalDistance(deltaListe(EUREF_S, EUREF_H)), horizontalDistance(deltaListe(EUREF_L, EUREF_H))

    print("\nKontrollsjekk:\n")
    print("Differanse S -> H:\t" + str(corrDistH) + " m = " + str(control1) + " m,\tdiff = " + str(np.abs(control1 - corrDistH)) + " m")
    print("Differanse L -> H:\t" + str(corrDistH2) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDistH2)) + " m")

# c)

def O1c():
    #Konverterer til lokalt system:
    
    HLG = CTRS_til_LG(SH, bl_S[0], bl_S[1])
    HLG2 = CTRS_til_LG(LH, bl_L[0], bl_L[1])

    # Beregner azimuth:

    aziH = radTilgrad(azimuth(HLG))
    aziH2 = radTilgrad(azimuth(HLG2))

    # Beregner korrigerte azimuth:

    corraziH = correctedAzimuth(HLG, bl_S[0], bl_S[1], UTM_S[:2], UTM_H[:2])
    corraziH2 = correctedAzimuth(HLG2, bl_L[0], bl_L[1], UTM_L[:2], UTM_H[:2])
    
    print("\n### Oppgave 1 ###\n\nc)\n")
    print("Azimuth fra Solemsvåttan til Holtermannsveien uten korrigeringer:\t" + str(aziH * 200/180) + " gon")
    print("Azimuth fra Ladehammeren til Holtermannsveien uten korrigeringer:\t\t" + str(aziH2 * 200/180) + " gon\n")
    print("Korrigert Solemsvåttan fra Moholt til Holtermannsveien:\t\t\t" + str(corraziH * 200/180) + " gon")
    print("Korrigert Ladehammeren fra ST46 til Holtermannsveien:\t\t\t" + str(corraziH2 * 200/180) + " gon")

    """
    Utfører kontroll ved å sammenligne de korrigerte azimuthene over med vanlig 
    azimuth målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = radTilgrad(azimuth(deltaListe(EUREF_S, EUREF_H))), radTilgrad(azimuth(deltaListe(EUREF_L, EUREF_H)))

    print("\nKontrollsjekk:\n")
    print("Differanse S -> H:\t" + str(corraziH * 200/180) + " gon = " + str(control1 * 200/180) + " gon,\tdiff = " + str(np.abs(control1 - corraziH) * 200/180) + " gon")
    print("Differanse L -> H:\t" + str(corraziH2 * 200/180) + " gon = " + str(control2 * 200/180) + " gon,\tdiff = " + str(np.abs(control2 - corraziH2) * 200/180) + " gon")

# d)

def O1d():
    #Konverterer til lokalt system:
    
    HLG = CTRS_til_LG(SH, bl_S[0], bl_S[1])
    HLG2 = CTRS_til_LG(LH, bl_L[0], bl_L[1])

    # Beregner korrigert høydeforskjell:

    corrDeltaHeightH = correctedHeightDifference(HLG, bl_S[0], EUREF_S[3]-EUREF_S[2], EUREF_H[3]-EUREF_H[2])
    corrDeltaHeightH2 = correctedHeightDifference(HLG2, bl_L[0], EUREF_L[3]-EUREF_L[2], EUREF_H[3]-EUREF_H[2])

    print("\n### Oppgave 1 ###\n\nd)\n")
    print("Korrigert NN2000 høydeforskjell fra Solemsvåttan til Holtermannsveien:\t" + str(corrDeltaHeightH) + " m")
    print("Korrigert NN2000 høydeforskjell fra Ladehammeren til Holtermannsveien:\t" + str(corrDeltaHeightH2) + " m")

    """
    Utfører kontroll ved å sammenligne de korrigerte høydeforskjellene over med
    høydeforskjellen beregnet ved å bruke NN2000-høydene gitt i 'Tabell 1' i oppgaveteksten
    """
    
    control1, control2 = EUREF_H[2] - EUREF_S[2], EUREF_H[2] - EUREF_L[2]

    print("\nKontrollsjekk:\n")
    print("Differanse S -> H:\t" + str(corrDeltaHeightH) + " m = " + str(control1) + " m,\t\tdiff = " + str(np.abs(control1 - corrDeltaHeightH)) + " m")
    print("Differanse L -> H:\t" + str(corrDeltaHeightH2) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDeltaHeightH2)) + " m")

# e)

# ...

### Oppgave 2 ###

# Cofactor matrix for dX, dY, dZ:
# Moholt -> TP342:

K = np.array([[1, 0.1307, 0.5335],
              [0.1307, 1, 0.1441],
              [0.5335, 0.1441, 1]])
S2 = np.array([[0.5*10**(-3), 0, 0],
              [0, 0.4*10**(-3), 0],
              [0, 0, 1.0*10**(-3)]])

# a)

def O2a():
    C1_CTRS = S2 @ K @ S2
    
    print("\n### Oppgave 2 ###\n\na)\n")
    print("Varians-kovarians-matrise for baseline S -> H:")
    print(C1_CTRS*10**6)
    
# b)

"""
I a) brukte en kofaktor-matrisen og RMS - root mean squared - fra baseline-rapporten til å finne varians-kovarians-matrisen i CTRS.

Nå må en først konvertere til LG ved bruk av error propagation law: C_LG = A * C_CTRS * A^T
Dette er fordi en bruker matematiske funksjoner og formler for å komme seg fra CTRS til LG: Y = Ax + B

Når T er funnet og en har varians og kovarians i LG kan en finne de avledede målene.

Observasjoner:
[] Distanse i kartprojeksjon
[] Azimuth i kartprojeksjon
[] Høydeforskjell i NN2000

Når uttrykkene for disse er funnet kan en partiellderivere med hensyn på x, y og z, sette de inn i en 3x3-matrise
og regne ut varians-kovarians-matrisen for avledede mål i kartprojeksjonen: C_mål = F * C_LG * F^T
"""

def O2b():
    C1_CTRS = S2 @ K @ S2

    B1, L1 = bl_S[0]*np.pi/180, bl_S[1]*np.pi/180 # Gjør om til rad

    A1 = np.array([[-np.sin(B1)*np.cos(L1), -np.sin(B1)*np.sin(L1), np.cos(B1)],
                [-np.sin(L1), np.cos(L1), 0],
                [np.cos(B1)*np.cos(L1), np.cos(B1)*np.sin(L1), np.sin(B1)]])
    C1_LG = h.matrixMultiplication(h.matrixMultiplication(A1, C1_CTRS), h.transposeMatrix(A1))
    
    """
        [[dB/dx, dB/dy, dB/dz],
    F =  [dD/dx, dD/dy, dD/dz],
         [dH/dx, dH/dy, dH/dz]]
    """

    HLG = CTRS_til_LG(SH, bl_S[0], bl_S[1])
    F1 =  partialDerivatives(HLG)
    
    C1 = h.matrixMultiplication(h.matrixMultiplication(F1, C1_LG), h.transposeMatrix(F1))

    print("\n### Oppgave 2 ###\n\nb)\n")
    print("Varians-kovarians-matrise for utledede observasjoner S -> H:")
    print(C1*10**6)
    print("\nStandardavvik azimuth i projeksjonsplanet: " + str(np.sqrt(C1[0][0]) * 10**3 * 200/180) + " mgon")
    print("Standardavvik høydedifferanse i NN2000: " + str(np.sqrt(C1[2][2]) * 10**3) + " mm")
    print("Standardavvik distanse i projeksjonsplanet: " + str(np.sqrt(C1[1][1]) * 10**3) + " mm")

O1a()
O1b()
O1c()
O1d()
O2a()
O2b()

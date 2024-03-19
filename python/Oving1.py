### Øving 1 ###

"""
GNSS Baselines

Fra GNSS baseline observasjoner til observasjoner i UTM kartprojeksjon og høyder i NN2000
"""

# Importerer biblioteker:

import numpy as np

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
    f2 = (l**3 / 3) * np.sin(B) * (np.cos(B))**2 * (1 + 3 * NU**2 + 2 * NU**4)
    f3 = (l**5 / 15) * np.sin(B) * (np.cos(B))**4 * (2 - (np.tan(B))**2)
    
    return radTilgrad(f1 + f2 + f3)

def correctionDelta(bl, B, p1, p2):
    """
    Beregner korreksjonsleddet delta for korrigeringen av azimuth
    """
    R = r(B, azimuth(bl))
    AB = radTilgrad(1 / (6 * R**2) * (2 * (p1[1] - 500000)/0.9996 + (p2[1] - 500000)/0.9996) * (p2[0] - p1[0])) # Korrigerer y for UTM-fellene
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

    ans[0][0] = -np.sin(azi) * rho / S_h # (I)
    ans[0][1] = np.cos(azi) * rho / S_h # (II)
    ans[1][0] = np.cos(azi) # (III)
    ans[1][1] = -np.sin(azi) # (IV)
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

EUREF_moholt    = [7031952.892, 571469.041, 131.404, 170.839] # Data fra oppgavetekst
EUREF_st46      = [7034487.402, 571578.304, 112.554, 151.851]
EUREF_tp342     = [7033941.628, 568890.318,   5.194,  44.618]
EUREF_havstein  = [7031656.237, 568780.646, 131.495, 171.046]

# Jeg velger baseline 3 og 4: Moholt -> TP342 og ST46 -> Moholt

# sted = [X, Y, Z]

moholt  = [2815415.137, 518305.734, 5680680.336] # Data fra vedlegg
st46    = [2813149.450, 518057.327, 5681796.430]

# baseline = [dX, dY, dZ]

dtp342  = [-1396.884, -2834.288,   802.443] # Moholt -> TP342 # Data fra vedlegg
dmoholt = [ 2265.690,   248.405, -1116.088] # ST46 -> Moholt

tp342 = sumListe(moholt, dtp342)

# Bredde- og lengdegrad for aktuelle punkt: sted = [B, L, h] (EUREF89 med ellipsoidisk høyde)

bl_moholt = koorkon.CartesianToLatLon(moholt[0], moholt[1], moholt[2])
bl_st46 = koorkon.CartesianToLatLon(st46[0], st46[1], st46[2])
bl_tp342 = koorkon.CartesianToLatLon(tp342[0], tp342[1], tp342[2])

"""
# Konvertert med skTrans:

bl_moholt   = [63.408951187, 10.431110376, 170.843]
bl_st46     = [63.431668740, 10.434433578, 151.865]
bl_tp342    = [63.427301973, 10.380349670,  44.613]
bl_havstein = [63.406818270, 10.377168783, 171.020]
"""

# UTM-koordinater for aktuelle punkt: sted = [N, E] (EUREF89)

UTM_moholt = koorkon.LatLonToUTM(bl_moholt[0], bl_moholt[1])
UTM_st46 = koorkon.LatLonToUTM(bl_st46[0], bl_st46[1])
UTM_tp342 = koorkon.LatLonToUTM(bl_tp342[0], bl_tp342[1])

# EUREF-data

a = 6378137 # [m]
b = 6356752.3141 # [m]
e = np.sqrt((a**2 - b**2) / a**2)

# Program:

### Oppgave 1 ###

# a)

def O1a():
    # Konverterer til lokalt system:

    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    print("\n### Oppgave 1 ###\n\na)\n")
    print("Lokale koordinat for TP342:\t" + str(tp342LG))
    print("Lokale koordinat for Moholt:\t" + str(moholtLG))

    """
    En enkel kontrollsjekk for å sikre riktig resultat er å beregne
    distansen til baselinen i CTRS og LG - disse skal være like lange!
    """

    len1, len2 = slopeDistance(dtp342), slopeDistance(dmoholt)
    len3, len4 = slopeDistance(tp342LG), slopeDistance(moholtLG)

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(len1) + " m = " + str(len3) + " m,\tdiff = " + str(np.abs(len1 - len3) * 10**3) + " mm")
    print("Differanse ST46 -> Moholt:\t" + str(len2) + " m = " + str(len4) + " m,\tdiff = " + str(np.abs(len2 - len4) * 10**3) + " mm")

# b)

def O1b():
    #Konverterer til lokalt system:

    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    # Beregner korrigerte distanser:

    corrDistTP342 = correctedHorizontalDistance(tp342LG, bl_moholt[0], bl_moholt[2], [UTM_moholt[1], UTM_tp342[1]])
    corrDistMoholt = correctedHorizontalDistance(moholtLG, bl_st46[0], bl_st46[2], [UTM_st46[1], UTM_moholt[1]])

    print("\n### Oppgave 1 ###\n\nb)\n")
    print("Korrigert distanse fra Moholt til TP342:\t" + str(corrDistTP342) + " m")
    print("Korrigert distanse fra ST46 til Moholt:\t\t" + str(corrDistMoholt) + " m")

    """
    Utfører kontroll ved å sammenligne de korrigerte distansene over med vanlig horisontal
    distanse målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = horizontalDistance(deltaListe(EUREF_moholt, EUREF_tp342)), horizontalDistance(deltaListe(EUREF_st46, EUREF_moholt))

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(corrDistTP342) + " m = " + str(control1) + " m,\tdiff = " + str(np.abs(control1 - corrDistTP342) * 10**3) + " mm")
    print("Differanse ST46 -> Moholt:\t" + str(corrDistMoholt) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDistMoholt) * 10**3) + " mm")

# c)

def O1c():
    #Konverterer til lokalt system:
    
    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    # Beregner azimuth:

    aziTP = radTilgrad(azimuth(tp342LG))
    aziMoh = radTilgrad(azimuth(moholtLG))

    # Beregner korrigerte azimuth:

    corraziTP = correctedAzimuth(tp342LG, bl_moholt[0], bl_moholt[1], UTM_moholt[:2], UTM_tp342[:2])
    corraziMoh = correctedAzimuth(moholtLG, bl_st46[0], bl_st46[1], UTM_st46[:2], UTM_moholt[:2])

    print("\n### Oppgave 1 ###\n\nc)\n")
    print("Azimuth fra Moholt til TP342 uten korrigeringer:\t" + str(aziTP * 200/180) + " gon")
    print("Azimuth fra ST46 til Moholt uten korrigeringer:\t\t" + str(aziMoh * 200/180) + " gon\n")
    print("Korrigert azimuth fra Moholt til TP342:\t\t\t" + str(corraziTP * 200/180) + " gon")
    print("Korrigert azimuth fra ST46 til Moholt:\t\t\t" + str(corraziMoh * 200/180) + " gon")

    """
    Utfører kontroll ved å sammenligne de korrigerte azimuthene over med vanlig 
    azimuth målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = radTilgrad(azimuth(deltaListe(EUREF_moholt, EUREF_tp342))), radTilgrad(azimuth(deltaListe(EUREF_st46, EUREF_moholt)))

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(corraziTP * 200/180) + " gon = " + str(control1 * 200/180) + " gon,\tdiff = " + str(np.abs(control1 - corraziTP) * 200/180 * 10**3) + " mgon")
    print("Differanse ST46 -> Moholt:\t" + str(corraziMoh * 200/180) + " gon = " + str(control2 * 200/180) + " gon,\tdiff = " + str(np.abs(control2 - corraziMoh) * 200/180 * 10**3) + " mgon")

# d)

def O1d():
    #Konverterer til lokalt system:
    
    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    # Beregner korrigert høydeforskjell:

    corrDeltaHeightTP = correctedHeightDifference(tp342LG, bl_moholt[0], EUREF_moholt[3]-EUREF_moholt[2], EUREF_tp342[3]-EUREF_tp342[2])
    corrDeltaHeightMoh = correctedHeightDifference(moholtLG, bl_st46[0], EUREF_st46[3]-EUREF_st46[2], EUREF_moholt[3]-EUREF_moholt[2])

    print("\n### Oppgave 1 ###\n\nd)\n")
    print("Korrigert NN2000 høydeforskjell fra Moholt til TP342:\t" + str(corrDeltaHeightTP) + " m")
    print("Korrigert NN2000 høydeforskjell fra ST46 til Moholt:\t" + str(corrDeltaHeightMoh) + " m")

    """
    Utfører kontroll ved å sammenligne de korrigerte høydeforskjellene over med
    høydeforskjellen beregnet ved å bruke NN2000-høydene gitt i 'Tabell 1' i oppgaveteksten
    """
    
    control1, control2 = EUREF_tp342[2] - EUREF_moholt[2], EUREF_moholt[2] - EUREF_st46[2]

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(corrDeltaHeightTP) + " m = " + str(control1) + " m,\t\tdiff = " + str(np.abs(control1 - corrDeltaHeightTP) * 10**3) + " mm")
    print("Differanse ST46 -> Moholt:\t" + str(corrDeltaHeightMoh) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDeltaHeightMoh) * 10**3) + " mm")

# e)

# ...

### Oppgave 2 ###

# Cofactor matrix for dX, dY, dZ:
# Moholt -> TP342:
RMS1 = 0.7645
Q1 = np.array([[0.00000008, 0.00000001, 0.00000007],
               [0.00000001, 0.00000006, 0.00000003],
               [0.00000007, 0.00000003, 0.00000039]])

# ST46 -> Moholt:
RMS2 = 0.3428
Q2 = np.array([[0.00000012, 0.00000003, 0.00000010],
               [0.00000003, 0.00000008, 0.00000005],
               [0.00000010, 0.00000005, 0.00000060]])

# a)

def O2a():
    C1_CTRS = RMS1**2 * Q1
    C2_CTRS = RMS2**2 * Q2

    print("\n### Oppgave 2 ###\n\na)\n")
    print("Varians-kovarians-matrise for baseline Moholt -> TP342 [mm]:")
    print(C1_CTRS * 10**6)
    print("\nVarians-kovarians-matrise for baseline ST46 -> Moholt [mm]:")
    print(C2_CTRS * 10**6)

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
    C1_CTRS = RMS1**2 * Q1
    C2_CTRS = RMS2**2 * Q2

    rho = np.pi / 180 # Gjør om til rad

    B1, L1 = bl_moholt[0] * rho, bl_moholt[1] * rho
    B2, L2 = bl_st46[0] * rho, bl_st46[1] * rho

    A1 = np.array([[-np.sin(B1)*np.cos(L1), -np.sin(B1)*np.sin(L1), np.cos(B1)],
                [-np.sin(L1), np.cos(L1), 0],
                [np.cos(B1)*np.cos(L1), np.cos(B1)*np.sin(L1), np.sin(B1)]])

    A2 = np.array([[-np.sin(B2)*np.cos(L2), -np.sin(B2)*np.sin(L2), np.cos(B2)],
                [-np.sin(L2), np.cos(L2), 0],
                [np.cos(B2)*np.cos(L2), np.cos(B2)*np.sin(L2), np.sin(B2)]])

    C1_LG = h.matrixMultiplication(h.matrixMultiplication(A1, C1_CTRS), h.transposeMatrix(A1))
    C2_LG = h.matrixMultiplication(h.matrixMultiplication(A2, C2_CTRS), h.transposeMatrix(A2))
    
    """
        [[dB/dx, dB/dy, dB/dz],
    F =  [dD/dx, dD/dy, dD/dz],
         [dH/dx, dH/dy, dH/dz]]
    """

    #Konverterer til lokalt system:
    
    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    F1 = partialDerivatives(tp342LG)
    F2 = partialDerivatives(moholtLG)

    C1 = h.matrixMultiplication(h.matrixMultiplication(F1, C1_LG), h.transposeMatrix(F1))
    C2 = h.matrixMultiplication(h.matrixMultiplication(F2, C2_LG), h.transposeMatrix(F2))

    print("\n### Oppgave 2 ###\n\nb)\n")
    print("Varians-kovarians-matrise for utledede observasjoner Moholt -> TP342 [milli]:")
    print(C1 * 10**6)
    print("\nStandardavvik azimuth i projeksjonsplanet: " + str(np.sqrt(C1[0][0]) * 10**3 * 200/180) + " mgon")
    print("Standardavvik høydedifferanse i NN2000: " + str(np.sqrt(C1[2][2]) * 10**3) + " mm")
    print("Standardavvik distanse i projeksjonsplanet: " + str(np.sqrt(C1[1][1]) * 10**3) + " mm")
    print("\nVarians-kovarians-matrise for utledede observasjoner ST46 -> Moholt [milli]:")
    print(C2 * 10**6)
    print("\nStandardavvik azimuth i projeksjonsplanet: " + str(np.sqrt(C2[0][0]) * 10**3 * 200/180) + " mgon")
    print("Standardavvik høydedifferanse i NN2000: " + str(np.sqrt(C2[2][2]) * 10**3) + " mm")
    print("Standardavvik distanse i projeksjonsplanet: " + str(np.sqrt(C2[1][1]) * 10**3) + " mm")

# c)

# ...

### Oppgave 3 ###

gisline = dict()

"""

dict[sted] =
[[x_LG, y_LG, z_LG],        [m]             [0]
 [sdx_LG, sdy_LG, sdz_LG],  [m]             [1]
 [rxy, rxz, ryz],           [m]             [2]
 [azi, vert, avst],         [gon, gon, m]   [3]
 [sdazi, sdvert, sdavst],   [gon, gon, m]   [4]
 [rRdH, rRD, rDdH],         [gon / m]       [5]
 [Korreksjoner]             [gon / m]       [6]
 [R, D, dH]]                [gon, m, m]     [7]

"""

gisline["mohtp"] = [[2046.5498, -2534.5363, -127.0737],
                    [0.0022, 0.0019, 0.0049],
                    [0.0138, 0.5263, 0.1381],
                    [343.24407, 102.48206, 3260.12088],
                    [0.00004, 0.00492, 0.00205],
                    [0.5064, 0.1281, 0.2575],
                    [-1.42199, -0.00011, -2.4992, -1.1065, 0.01100],
                    [341.82198, 3256.515, -126.233]]

gisline["stmoh"] = [[-2532.3500, -166.0425, 18.4962],
                    [0.0013, 0.0009, 0.0029],
                    [-0.0367, 0.6394, 0.2007],
                    [204.16826, 99.53602, 2537.85513],
                    [0.00002, 0.00290, 0.00130],
                    [-0.1416, -0.0795, -0.6493],
                    [-1.42557, 0.00014, -0.1351, -0.8561, -0.13800],
                    [202.74283, 2536.864, 18.864]]

def O3a():
    ###
    K1 = np.array([[1, gisline["mohtp"][2][0], gisline["mohtp"][2][1]],
                   [gisline["mohtp"][2][0], 1, gisline["mohtp"][2][2]],
                   [gisline["mohtp"][2][1], gisline["mohtp"][2][2], 1]])
    K2 = np.array([[1, gisline["stmoh"][2][0], gisline["stmoh"][2][1]],
                   [gisline["stmoh"][2][0], 1, gisline["stmoh"][2][2]],
                   [gisline["stmoh"][2][1], gisline["stmoh"][2][2], 1]])
    S1 = np.array([[gisline["mohtp"][1][0], 0, 0],
                   [0, gisline["mohtp"][1][1], 0],
                   [0, 0, gisline["mohtp"][1][2]]])
    S2 = np.array([[gisline["mohtp"][1][0], 0, 0],
                   [0, gisline["mohtp"][1][1], 0],
                   [0, 0, gisline["mohtp"][1][2]]])
    C1LG = S1 @ K1 @ S1
    C2LG = S2 @ K2 @ S2

    F1 = partialDerivatives(gisline["mohtp"][0])
    F2 = partialDerivatives(gisline["stmoh"][0])

    C1 = F1 @ C1LG @ h.transposeMatrix(F1)
    C2 = F2 @ C2LG @ h.transposeMatrix(F2)

    G1 = np.array([[1, gisline["mohtp"][5][1], gisline["mohtp"][5][0]],
                   [gisline["mohtp"][5][1], 1, gisline["mohtp"][5][2]],
                   [gisline["mohtp"][5][0], gisline["mohtp"][5][2], 1]])
    G2 = np.array([[1, gisline["stmoh"][5][1], gisline["stmoh"][5][0]],
                   [gisline["stmoh"][5][1], 1, gisline["stmoh"][5][2]],
                   [gisline["stmoh"][5][0], gisline["stmoh"][5][2], 1]])
    ###
    print("### Oppgave 3 ###\n\na)\n")
    print("Baseline: Moholt -> TP342\n")
    print("Lokalt koordinat for TP342: [" + str(gisline["mohtp"][0][0]) + ", " + str(gisline["mohtp"][0][1]) + ", " + str(gisline["mohtp"][0][2]) + "] [m]")
    print("Ellipsoidiske data:")
    print("Distanse:\t" + format(gisline["mohtp"][3][2], '.4f') + " m")
    print("Azimuth:\t " + format(gisline["mohtp"][3][0], '.4f') + " gon")
    print("Senit-vinkel:\t  " + format(gisline["mohtp"][3][1], '.4f') + " gon")
    print("Korreksjoner ellipsoide -> kartprojeksjon:")
    print("Distanse: ellipse -> terreng: " + format(gisline["mohtp"][6][2], '.4f') + " m, og kartprojeksjon: " + format(gisline["mohtp"][6][3], '.4f') + " m")
    print("Azimuth: konvergens meridian: " + format(gisline["mohtp"][6][0], '.4f') + " gon, og kartprojeksjon: " + format(gisline["mohtp"][6][1], '.4f') + " gon")
    print("Senit-vinkel: loddavvik: " + format(gisline["mohtp"][6][4], '.4f') + " gon")
    print("Observasjoner i kartplanet:")
    print("Distanse:\t" + format(gisline["mohtp"][7][1], '.4f') + " m")
    print("Azimuth:\t" + format(gisline["mohtp"][7][0], '.4f') + " gon")
    print("Høydeforskjell:\t" + format(gisline["mohtp"][7][2], '.4f') + " m")
    print("Nøyaktighet observasjoner:")
    print("Beregnet varians-kovarians-matrise:")
    print(C1)
    print("Beregnet korrelasjonsmatrise:")
    print(K1)
    print("Korrelasjonsmatrise fra GISLINE:")
    print(G1)
    print("Sammenlignet beregnet varians med GISLINE-varians:")
    print("Distanse:\t" + format(np.sqrt(C1[0][0]) * 10**3, '.4f') + " - " + format(gisline["mohtp"][4][2] * 10**3, '.4f') + " [mm]")
    print("Azimuth:\t" + format(np.sqrt(C1[1][1]) * 10**3, '.4f') + " - " + format(gisline["mohtp"][4][0] * 10**3, '.4f') + " [mgon]")
    print("Høydeforskjell:\t" + format(np.sqrt(C1[2][2]) * 10**3, '.4f') + " - " + format(gisline["mohtp"][4][1] * 10**3, '.4f') + " [mgon]")
    print()
    print("Baseline: ST46 -> Moholt\n")
    print("Lokalt koordinat for Moholt: [" + str(gisline["stmoh"][0][0]) + ", " + str(gisline["stmoh"][0][1]) + ", " + str(gisline["stmoh"][0][2]) + "] [m]")
    print("Ellipsoidiske data:")
    print("Distanse:\t" + format(gisline["stmoh"][3][2], '.4f') + " m")
    print("Azimuth:\t " + format(gisline["stmoh"][3][0], '.4f') + " gon")
    print("Senit-vinkel:\t  " + format(gisline["stmoh"][3][1], '.4f') + " gon")
    print("Korreksjoner ellipsoide -> kartprojeksjon:")
    print("Distanse: ellipse -> terreng: " + format(gisline["stmoh"][6][2], '.4f') + " m, og kartprojeksjon: " + format(gisline["stmoh"][6][3], '.4f') + " m")
    print("Azimuth: konvergens meridian: " + format(gisline["stmoh"][6][0], '.4f') + " gon, og kartprojeksjon: " + format(gisline["stmoh"][6][1], '.4f') + " gon")
    print("Senit-vinkel: loddavvik: " + format(gisline["stmoh"][6][4], '.4f') + " gon")
    print("Observasjoner i kartplanet:")
    print("Distanse:\t" + format(gisline["stmoh"][7][1], '.4f') + " m")
    print("Azimuth:\t" + format(gisline["stmoh"][7][0], '.4f') + " gon")
    print("Høydeforskjell:\t" + format(gisline["stmoh"][7][2], '.4f') + " m")
    print("Nøyaktighet:")
    print("Beregnet varians-kovarians-matrise:")
    print(C2)
    print("Beregnet korrelasjonsmatrise:")
    print(K2)
    print("Korrelasjonsmatrise fra GISLINE:")
    print(G2)
    print("Sammenlignet beregnet varians med GISLINE-varians:")
    print("Distanse:\t" + format(np.sqrt(C2[0][0]) * 10**3, '.4f') + " - " + format(gisline["stmoh"][4][2] * 10**3, '.4f') + " [mm]")
    print("Azimuth:\t" + format(np.sqrt(C2[1][1]) * 10**3, '.4f') + " - " + format(gisline["stmoh"][4][0] * 10**3, '.4f') + " [mgon]")
    print("Høydeforskjell:\t" + format(np.sqrt(C2[2][2]) * 10**3, '.4f') + " - " + format(gisline["stmoh"][4][1] * 10**3, '.4f') + " [mgon]")

# Kjøring av program:

"""
O1a()
O1b()
O1c()
O1d()
O2a()
O2b()
"""
O3a()

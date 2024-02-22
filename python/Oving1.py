# Øving 1

"""
GNSS Baselines

Fra GNSS baseline observasjoner til observasjoner i UTM kartprojeksjon og høyder i NN2000
"""

# Importerer biblioteker:

import Hjelpefunksjoner as h
import KoordinatKonvertering as koorkon
import numpy  as np

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

def negativListe(l1):
    """
    Snur fortegn på alle verdiene i listen
    """
    ans = []
    for i in l1:
        ans.append(-i)
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
    return (azimuth(bl) * 200 / np.pi - np.abs(korC) - np.abs(korDelta[0])) * 180 / 200

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
    AB = radTilgrad(1 / (6 * R**2) * (2 * p1[1] + p2[1]) * (p2[0] - p1[0]))
    BA = radTilgrad(1 / (6 * R**2) * (2 * p2[1] + p1[1]) * (p1[0] - p2[0]))
    return AB, BA

def nu(B):
    """
    Hjelpeledd for korreksjonen i C
    """
    return np.sqrt(e**4 * (np.cos(B))**2 / (1 - e**4))

# d)

def correctedHeightDifference(bl, B, NStart, NSlutt):
    return bl[2] + correctedEarthCurv(bl, B) - (NSlutt - NStart)

def correctedEarthCurv(bl, B):
    return (slopeDistance(bl) * np.sin(zenithAngle(bl)))**2 / (2 * r(B, azimuth(bl)))

# e)

# ...

# Data:

"""
Fasit:

EUREF89: sted = [N (UTM), E (UTM), H (NN2000), h (ellipsoide)]
"""

EUREF_moholt = [7031952.892, 571469.041, 131.404, 170.839] # Data fra oppgavetekst
EUREF_st46 = [7034487.402, 571578.304, 112.554, 151.851]
EUREF_tp342 = [7033941.628, 568890.318, 5.194, 44.618]
EUREF_havstein = [7031656.237, 568780.646, 131.495, 171.046]

# Jeg velger baseline 3 og 4: Moholt -> TP342 og ST46 -> Moholt

# sted = [X, Y, Z]

moholt = [2815415.137, 518305.734, 5680680.336] # Data fra vedlegg
st46 = [2813149.450, 518057.327, 5681796.430]

# baseline = [dX, dY, dZ]

dtp342 = [-1396.884, -2834.288, 802.443] # Moholt -> TP342 # Data fra vedlegg
dmoholt = [2265.690, 248.405, -1116.088] # ST46 -> Moholt

tp342 = sumListe(moholt, dtp342)

# Bredde- og lengdegrad for aktuelle punkt: sted = [B, L, h] (EUREF89 med ellipsoidisk høyde)

bl_moholt = koorkon.CartesianToLatLon(moholt[0], moholt[1], moholt[2])
bl_st46 = koorkon.CartesianToLatLon(st46[0], st46[1], st46[2])
bl_tp342 = koorkon.CartesianToLatLon(tp342[0], tp342[1], tp342[2])

"""
# Konvertert med skTrans:

bl_moholt = [63.408951187, 10.431110376, 170.843]
bl_st46 = [63.431668740, 10.434433578, 151.865]
bl_tp342 = [63.427301973, 10.380349670, 44.613]
bl_havstein = [63.406818270, 10.377168783, 171.020]
"""

# UTM-koordinater for aktuelle punkt: sted = [N, E] (EUREF89)

UTM_moholt = koorkon.LatLonToUTM(bl_moholt[0], bl_moholt[1])
UTM_st46 = koorkon.LatLonToUTM(bl_st46[0], bl_st46[1])
UTM_tp342 = koorkon.LatLonToUTM(bl_tp342[0], bl_tp342[1])

# EUREF-data

a = 6378137 # [m]
b = 6356752.3141 # [m]
#f = (a - b) / a # []
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
    print("Differanse Moholt -> TP342:\t" + str(len1) + " m = " + str(len3) + " m,\tdiff = " + str(np.abs(len1 - len3)) + " m")
    print("Differanse ST46 -> Moholt:\t" + str(len2) + " m = " + str(len4) + " m,\tdiff = " + str(np.abs(len2 - len4)) + " m")

# b)

def O1b():
    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    # Beregner korrigerte distanser:

    corrDistTP342 = correctedHorizontalDistance(tp342LG, bl_moholt[0], bl_moholt[2], [UTM_moholt[1], UTM_tp342[1]])
    corrDistMoholt =  correctedHorizontalDistance(moholtLG, bl_st46[0], bl_st46[2], [UTM_st46[1], UTM_moholt[1]])

    print("\n### Oppgave 1 ###\n\nb)\n")
    print("Korrigert distanse fra Moholt til TP342:\t" + str(corrDistTP342) + " m")
    print("Korrigert distanse fra ST46 til Moholt:\t\t" + str(corrDistMoholt) + " m")

    """
    Utfører kontroll ved å sammenligne de korrigerte distansene over med vanlig horisontal
    distanse målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = horizontalDistance(deltaListe(EUREF_moholt, EUREF_tp342)), horizontalDistance(deltaListe(EUREF_st46, EUREF_moholt))

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(corrDistTP342) + " m = " + str(control1) + " m,\tdiff = " + str(np.abs(control1 - corrDistTP342)) + " m")
    print("Differanse ST46 -> Moholt:\t" + str(corrDistMoholt) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDistMoholt)) + " m")

# c)

def O1c():
    tp342LG = CTRS_til_LG(dtp342, bl_moholt[0], bl_moholt[1])
    moholtLG = CTRS_til_LG(dmoholt, bl_st46[0], bl_st46[1])

    # Beregner azimuth:

    aziTP = radTilgrad(azimuth(tp342LG))
    aziMoh = radTilgrad(azimuth(moholtLG))

    # Beregner korrigerte azimuth:

    corraziTP = correctedAzimuth(tp342LG, bl_moholt[0], bl_moholt[1], UTM_moholt[:2], UTM_tp342[:2])
    corraziMoh = correctedAzimuth(moholtLG, bl_st46[0], bl_st46[1], UTM_st46[:2], UTM_moholt[:2])

    print("\n### Oppgave 1 ###\n\nc)\n")
    print("Azimuth fra Moholt til TP342 uten korrigeringer:\t" + str(aziTP) + " grader")
    print("Azimuth fra ST46 til Moholt uten korrigeringer:\t\t" + str(aziMoh) + " grader\n")
    print("Korrigert azimuth fra Moholt til TP342:\t\t\t" + str(corraziTP) + " grader")
    print("Korrigert azimuth fra ST46 til Moholt:\t\t\t" + str(corraziMoh) + " grader")

    """
    Utfører kontroll ved å sammenligne de korrigerte azimuthene over med vanlig 
    azimuth målt ved å bruke UTM-kordinatene gitt i 'Tabell 1' i oppgaveteksten
    """

    control1, control2 = radTilgrad(azimuth(deltaListe(EUREF_moholt, EUREF_tp342))), radTilgrad(azimuth(deltaListe(EUREF_st46, EUREF_moholt)))

    print("\nKontrollsjekk:\n")
    print("Differanse Moholt -> TP342:\t" + str(corraziTP) + " grader = " + str(control1) + " grader,\tdiff = " + str(np.abs(control1 - corraziTP)) + " grader")
    print("Differanse ST46 -> Moholt:\t" + str(corraziMoh) + " grader = " + str(control2) + " grader,\tdiff = " + str(np.abs(control2 - corraziMoh)) + " grader")

# d)

def O1d():
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
    print("Differanse Moholt -> TP342:\t" + str(corrDeltaHeightTP) + " m = " + str(control1) + " m,\t\tdiff = " + str(np.abs(control1 - corrDeltaHeightTP)) + " m")
    print("Differanse ST46 -> Moholt:\t" + str(corrDeltaHeightMoh) + " m = " + str(control2) + " m,\tdiff = " + str(np.abs(control2 - corrDeltaHeightMoh)) + " m")

# e)

# ...

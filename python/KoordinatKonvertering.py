# Importerer biblioteker:

import numpy as np

# Definerer konstanter:

# EUREF89:
a = 6378137 # [m]
b = 6356752.3141 # [m]
#f = (a - b) / a # []
e = np.sqrt((a**2 - b**2) / a**2) # []

# Definerer funksjoner:

def gradTilrad(v):
    """
    Gjør om vinkel v fra grader til radianer
    """
    return v * np.pi / 180

def radTilgrad(v):
    """
    Gjør om vinkel v fra radianer til grader
    """
    return v * 180 / np.pi

# Hjelpefunksjoner for UTM-konvertering:
def A(n):
    return a * (1 - n + (5 / 4) * (n**2 - n**3) + (81 / 64) * (n**4 - n**5))

def B(n):
    return (3 / 2) * a * (n - n**2 + (7 / 8) * (n**3 - n**4) + (55 / 64) * n**5)

def C(n):
    return (15 / 16) * a * (n**2 - n**3 + (3 / 4) * (n**4 - n**5))

def D(n):
    return (35 / 48) * a * (n**3 - n**4 + (11 / 16) * n**5)

def E(n):
    return (315 / 512) * a * (n**4 - n**5)

def CartesianToLatLon(x, y, z):
    """
    Funksjon som tar inn kartesiske [X,Y,Z] og returnerer bredde-
    og lengdegrader i tillegg til høyde koblet til EUREF89
    """

    # Finner lengdegraden:
    L = np.arctan(y / x)

    p = np.sqrt(x**2 + y**2) # Hjelpekonstant (horisontal distanse)

    # Finner breddegrad og høyde via iterasjon

    oldB = 0
    newB = np.arctan(z / ((1 - e**2) * p))

    while np.abs(newB - oldB) > 10**(-12):
        n = a / np.sqrt(1 - e**2 * (np.sin(newB))**2)
        h = p / (np.cos(newB)) - n # Høyde
        oldB = newB
        newB = np.arctan(z / ((1 - e**2 * n / (n + h)) * p)) # Breddegrad

    rho = 180 / np.pi # Gjør om til grader fra rad

    return [newB * rho, L * rho, h]

def  latLonToUTM(br, le):
    """
    Funksjon som tar inn bredde- og lengdegrad og
    returnerer UTM-koordinater i riktig sone
    """

    # Gjør om til rad:
    br, le = gradTilrad(br), gradTilrad(le)

    # Konverterings-konstanter:
    l0 = gradTilrad((radTilgrad(le) // 6) * 6 + 3)
    k0 = 0.9996
    e2 = e**2 / (1 - e**2)
    n = (a - b) / (a + b)
    #rho = a * (1 - e**2) / (1 - e**2 * (np.sin(br))**2)**(3 / 2)
    nu = a / np.sqrt(1 - e**2 * (np.sin(br))**2)
    p = le - l0

    S = A(n) * br - B(n) * np.sin(2 * br) + C(n) * np.sin(4 * br) - D(n) * np.sin(6 * br) + E(n) * np.sin(8 * br)
    
    k1 = S * k0
    k2 = k0 * nu * np.sin(2 * br) / 4
    k3 = (k0 * nu * np.sin(br) * (np.cos(br))**3 / 24) * (5 - (np.tan(br))**2 + 9 * e2 * (np.cos(br))**2 + 4 * e2**2 * (np.cos(br))**4)
    k4 = k0 * nu * np.cos(br)
    k5 = (k0 * nu * (np.cos(br))**3 / 6) * (1  - (np.tan(br))**2 + e2 * (np.cos(br))**2)

    # De faktiske UTM-koordinatene:
    x = k1 + k2 * p**2 + k3 * p**4
    y = k4 * p + k5 * p**3 + 500000
    
    return [x, y]

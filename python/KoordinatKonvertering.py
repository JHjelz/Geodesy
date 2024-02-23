# Importerer biblioteker:

import numpy as np

# Definerer konstanter:

def konstanter(system="EUREF89"):
    if (system == "EUREF89"):
        # EUREF89:
        a = 6378137 # [m]
        b = 6356752.3141 # [m]
        e = np.sqrt((a**2 - b**2) / a**2) # []
    elif (system ==  "WGS84"):
        # WGS84:
        a = 6378137 # [m]
        b = 6356752.3142 # [m]
        e = np.sqrt((a**2 - b**2) / a**2) # []
    elif (system == "ED50"):
        # ED50 / ED87:
        a = 6378388 # [m]
        b = 6356911.9461 # [m]
        e = np.sqrt((a**2 - b**2) / a**2) # []
    elif (system == "NGO1948"):
        # NGO1948:
        a = 6377492.0176 # [m]
        b = 6356173.5083 # [m]
        e = np.sqrt((a**2 - b**2) / a**2) # []

    #f = (a - b) / a # []

    return a, b, e

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
def A(a, n):
    return a * (1 - n + (5 / 4) * (n**2 - n**3) + (81 / 64) * (n**4 - n**5))

def B(a, n):
    return (3 / 2) * a * (n - n**2 + (7 / 8) * (n**3 - n**4) + (55 / 64) * n**5)

def C(a, n):
    return (15 / 16) * a * (n**2 - n**3 + (3 / 4) * (n**4 - n**5))

def D(a, n):
    return (35 / 48) * a * (n**3 - n**4 + (11 / 16) * n**5)

def E(a, n):
    return (315 / 512) * a * (n**4 - n**5)

def CartesianToLatLon(x, y, z, system="EUREF89"):
    """
    Funksjon som tar inn kartesiske [X,Y,Z] og returnerer bredde-
    og lengdegrader i tillegg til høyde koblet til riktig system
    """

    a, b, e = konstanter(system)

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

def LatLonToUTM(br, le, system="EUREF89"):
    """
    Funksjon som tar inn bredde- og lengdegrad og
    returnerer UTM-koordinater i riktig system
    """

    a, b, e = konstanter(system)

    # Gjør om til rad:
    br, le = gradTilrad(br), gradTilrad(le)
    
    # Konverterings-konstanter:
    l0 = gradTilrad((radTilgrad(le) // 6) * 6 + 3)
    k0 = 0.9996
    e2 = e**2 / (1 - e**2)
    n = (a - b) / (a + b)
    nu = a / np.sqrt(1 - e**2 * (np.sin(br))**2)
    p = le - l0

    S = A(a, n) * br - B(a, n) * np.sin(2 * br) + C(a, n) * np.sin(4 * br) - D(a, n) * np.sin(6 * br) + E(a, n) * np.sin(8 * br)
    
    k1 = S * k0
    k2 = k0 * nu * np.sin(2 * br) / 4
    k3 = (k0 * nu * np.sin(br) * (np.cos(br))**3 / 24) * (5 - (np.tan(br))**2 + 9 * e2 * (np.cos(br))**2 + 4 * e2**2 * (np.cos(br))**4)
    k4 = k0 * nu * np.cos(br)
    k5 = (k0 * nu * (np.cos(br))**3 / 6) * (1  - (np.tan(br))**2 + e2 * (np.cos(br))**2)

    # De faktiske UTM-koordinatene:
    x = k1 + k2 * p**2 + k3 * p**4
    y = k4 * p + k5 * p**3 + 500000
    
    return [x, y]

def CartesianToUTM(x, y, z, system="EUREF89"):
    A1 = CartesianToLatLon(x, y, z, system)
    A2 = LatLonToUTM(A1[0], A1[1], system)
    return [A2[0], A2[1], A1[2]]

# Så må en motsatt vei:

def LatLonToCartesian(br, le, h, system="EUREF89"):
    """
    Funksjon som tar inn bredde- og lengdegrad, i tillegg til høyde,
    og returnerer kartesiske koordinater i gitt system
    """

    a, b, e = konstanter(system)

    # Gjør om til rad:
    br, le = gradTilrad(br), gradTilrad(le)

    # Konverteringskonstant:

    n = a**2 / np.sqrt(a**2 * (np.cos(br))**2 + b**2 * (np.sin(br))**2)

    # Konvertering:

    x = (n + h) * np.cos(br) * np.cos(le)
    y = (n + h) * np.cos(br) * np.sin(le)
    z = ((b**2 / a**2) * n + h) * np.sin(br)

    return [x, y, z]

def UTMToLatLon(x, y, sone, system="EUREF89"):
    """
    Funksjon som tar inn nord og øst i UTM-koordinater,
    og returnerer lengde- og breddegrader i gitt system
    """
    
    a, b, e = konstanter(system)
    y = y - 500000
    
    # Konverteringskonstanter:
    k0 = 0.9996
    m = x / k0
    mu  = m / (a * (1 - e**2 / 4 - 3 * e**4 / 64 - 5 * e**6 / 256))
    e1 = (1 - np.sqrt(1 - e**2)) / (1 + np.sqrt(1 - e**2))

    # 'Fotavtrykket' til breddegraden

    J1 = 3 * e1 / 2 - 27 * e1**3 / 32
    J2 = 21 * e1**2 / 16 - 55 * e1**4 / 32
    J3 = 151 * e1**3 / 96
    J4 = 1097  * e1**4 / 512

    fp = mu + J1 * np.sin(2 * mu) + J2 * np.sin(4 * mu) + J3 * np.sin(6 * mu) +  J4 * np.sin(8 * mu)

    e2 = e**2 / (1 - e**2)
    C1 = e2 * (np.cos(fp))**2
    T1 = (np.tan(fp))**2
    R1 = a * (1 - e**2) / (1 - e**2 * (np.sin(fp))**2)**(3/2)
    N1 = a / np.sqrt(1 - e**2 * (np.sin(fp))**2)
    D = y / (N1 * k0)

    Q1 = N1 * np.tan(fp) / R1
    Q2 = D**2 / 2
    Q3 = (5 + 3 * T1 + 10 * C1 - 4 * C1**2 - 9 * e2) * D**4 / 24
    Q4 = (61 + 90 * T1 + 298 * C1 + 45 * T1**2 - 3 * C1**2 - 252 * e2) * D**6 / 720
    Q5 = D
    Q6 = (1 + 2 * T1 + C1) * D**3 / 6
    Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1**2 + 8 * e2 + 24 * T1**2) * D**5 / 120

    if (sone == "32N"):
        l0 = 9
    elif (sone == "33N"):
        l0 = 15

    br = radTilgrad(fp - Q1 * (Q2 - Q3 + Q4))
    le = l0 + radTilgrad((Q5 - Q6 + Q7) / np.cos(fp))

    return [br, le]

def UTMToCartesian(x, y, h, sone, system="EUREF89"):
    A1 = UTMToLatLon(x, y, sone, system)
    return LatLonToCartesian(A1[0], A1[1], h, system)

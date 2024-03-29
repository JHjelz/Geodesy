"""

Program som kjører de faktiske oppgavene og konverteringsmulighetene

"""

import KoordinatKonvertering as koorkon
import Oving1 as O1

# Konstanter som holder oversikt over programmets status:

side = "hoved"
valgHoved = ['1', '2', 'avslutt']
valgKoor = ['a', 'b', 'c', 'd', 'e', 'f', 'tilbake', 'avslutt']
valgO1 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'tilbake', 'avslutt']
jaNei = ['J', 'N']
systemvalg = ['WGS84', 'ED50', 'ED87', 'NGO1948']

# Hjelpefunksjoner

def avslutteEllerTilbake(bruker):
    if (bruker == 'avslutt' or bruker == 'tilbake'):
        return "hoved"
    else:
        return side

def avslutteEllerTilbake_2(bruker):
    if (bruker == 'avslutt'):
        return "hoved"
    elif (bruker == 'tilbake'):
        return "tilbake"
    else:
        return side

def gyldigTall(n):
    try:
        n = float(n)
    except:
        return False
    return True

def giInnTall(streng):
    value = input(streng + ": ")
    side = avslutteEllerTilbake_2(value)
    
    if (side == 'hoved'):
        return value, "avslutt"
    elif (side == 'tilbake'):
        return value, "tilbake"
    
    while (not gyldigTall(value)):
        value = input("Ikke gyldig tall! " + streng + ": ")
        side = avslutteEllerTilbake_2(value)
        if (side == 'hoved'):
            return value, "avslutt"
        elif (side == 'tilbake'):
            return value, "tilbake"
    
    return value, "fortsett"

def giInnJaNei():
    valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
    side = avslutteEllerTilbake_2(valg)
    
    if (side == 'hoved'):
        return valg, "avslutt"
    elif (side == 'tilbake'):
        return valg, "tilbake"
    
    while (valg not in jaNei):
        valg = input("J/N: ")
        side = avslutteEllerTilbake_2(valg)
        if (side == 'hoved'):
            return valg, "avslutt"
        elif (side == 'tilbake'):
            return valg, "tilbake"
    
    return valg, "fortsatt"

def giInnSystem():
    print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
    system = input("System: ")
    side = avslutteEllerTilbake_2(system)
    
    if (side == 'hoved'):
        return system, "avslutt"
    elif (side == 'tilbake'):
        return system, "tilbake"
    
    while (system not in systemvalg):
        system = input("Ugyldig system! System: ")
        side = avslutteEllerTilbake_2(system)
        if (side == 'hoved'):
            return system, "avslutt"
        elif (side == 'tilbake'):
            return system, "tilbake"
    
    return system, "fortsett"

# Program:

print("""
### ### ### ### ### ### ### ### ### ### ### ###

Hei og velkommen til koordinatkonvertereren!

### ### ### ### ### ### ### ### ### ### ### ###

Du er nå på hovedsiden der du har følgende valg:

[1] Se svarene for 'Øving 1'
[2] Konvertere egne koordinater til andre typer av koordinater

Du kan når som helst skrive 'avslutt' for å avslutte programmet, eller 'tilbake' for å komme tilbake til forrige meny.

Skriv et tall i følgende intervall, 1 - 2, for å velge funksjonalitet basert på valgmulighetene over.
""")

bruker = input("Jeg vil utføre handling: ")
while (bruker not in valgHoved):
    bruker = input("Ugyldig input! Utfører handling: ")

print("\n### ### ### ### ### ### ### ### ### ### ### ###")

while (bruker != 'avslutt'):
    if (bruker == '1'):
        side = "O1"
    
    while (side == 'O1'):
        print("""\nDu kan nå se svarene for Øving 1.

I denne øvinga har jeg valgt å studere følgende baselines:
[] Moholt -> TP342
[] ST46 -> Moholt

Du har følgende valg:

[a] Se lokale koordinat for TP342 og Moholt
[b] Se korrigert distanse til kartplanet
[c] Se korrigert azimuth til kartplanet
[d] Se korrigert høydeforskjell i NN2000
[e] Se varians-kovarians matrisen for [dX, dY, dZ] til baselinene
[f] Se varians-kovarians matrisen for utledede observasjoner
[g] Se sammenligning med data fra GISLINE Trans - Landmåling

Skriv inn ønsket handling som en bokstav i input-feltet.
        """)
        bruker = input("Nå vil jeg se: ")
        
        while (bruker not in valgO1):
            print("\nUgyldig input! Prøv igjen")
            bruker = input("Jeg vil utføre handling: ")

        side = avslutteEllerTilbake(bruker)

        if (bruker == 'a'):
            O1.O1a()
            input()
        elif (bruker == 'b'):
            O1.O1b()
            input()
        elif (bruker ==  'c'):
            O1.O1c()
            input()
        elif (bruker == 'd'):
            O1.O1d()
            input()
        elif (bruker == 'e'):
            O1.O2a()
            input()
        elif (bruker == 'f'):
            O1.O2b()
            input()
        elif (bruker == 'g'):
            O1.O3a()
            input()
    
    if (bruker == '2'):
        side = "Koordinatkonverterer"
    
    while (side == 'Koordinatkonverterer'):
        print("""\nHer kan du konvertere mellom ulike typer koordinater i gitte system.

Typene du kan velge mellom er: kartesiske koordinater, bredde- og lengdegrad, og UTM-koordinater.
Systemene du kan velge mellom er: EUREF89, WGS84, ED50/ED87 og NGO1948 (EUREF89 er forhåndsinnstilt)

Du har følgende valg:

[a] Kartesiske koordinater til bredde- og lengdegrad
[b] Kartesiske koordinater til UTM-koordinater
[c] Bredde- og lengdegrad til UTM-koordinater
[d] Bredde- og lengdegrad til kartesiske koordinater
[e] UTM-koordinater til bredde- og lengdegrad
[f] UTM-koordinater til kartesiske koordinater

Skriv inn ønsket handling som en bokstav i input-feltet.
        """)
        
        bruker = input("Jeg vil utføre handling: ")

        while (bruker not in valgKoor):
            print("\nUgyldig input! Prøv igjen")
            bruker = input("Jeg vil utføre handling: ")

        side = avslutteEllerTilbake(bruker)
        
        if (bruker == 'a'):
            print("""
Du skal nå konvertere kartesiske koordinater (x,y,z) til lengde- og breddegrad.
Skriv inn verdiene du vil konvertere:
""")
            x, handling = giInnTall('X-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            y, handling = giInnTall('Y-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            z, handling = giInnTall('Z-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""

            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break
            
            if (system == ""):
                resultat = koorkon.CartesianToLatLon(float(x), float(y), float(z))
            else:
                resultat = koorkon.CartesianToLatLon(float(x), float(y), float(z), system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + x + ", " + y + ", " + z + "]")
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + ", " + str(resultat[2]) + "]")
            print("### ####\n")
            input()

        elif (bruker == 'b'):
            print("""
Du skal nå konvertere kartesiske koordinater (x,y,z) til UTM-koordinater.
Skriv inn verdiene du vil konvertere:
""")
            x, handling = giInnTall('X-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            y, handling = giInnTall('Y-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            z, handling = giInnTall('Z-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""
            
            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break

            if (system == ""):
                resultat = koorkon.CartesianToUTM(float(x), float(y), float(z))
            else:
                resultat = koorkon.CartesianToUTM(float(x), float(y), float(z), system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + x + ", " + y + ", " + z + "]")
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + ", " + str(resultat[2]) + "]")
            print("### ####\n")
            input()

        elif (bruker == 'c'):
            print("""
Du skal nå konvertere bredde- og lengdegrader (B, L, h) til UTM-koordinater.
Skriv inn verdiene du vil konvertere:
""")
            B, handling = giInnTall('Breddegrad')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            L, handling = giInnTall('Lengdegrad')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""
            
            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break

            if (system == ""):
                resultat = koorkon.LatLonToUTM(float(B), float(L))
            else:
                resultat = koorkon.LatLonToUTM(float(B), float(L), system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + B + ", " + L + "]")
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + "]")
            print("### ####\n")
            input()

        elif (bruker == 'd'):
            print("""
Du skal nå konvertere bredde- og lengdegrader (B, L, h) til kartesiske koordinater.
Skriv inn verdiene du vil konvertere:
""")
            B, handling = giInnTall('Breddegrad')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            L, handling = giInnTall('Lengdegrad')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            h, handling = giInnTall('Høyde')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""
            
            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break

            if (system == ""):
                resultat = koorkon.LatLonToCartesian(float(B), float(L), float(h))
            else:
                resultat = koorkon.LatLonToCartesian(float(B), float(L), float(h), system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + B + ", " + L + ", " + h + "]")
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + ", " + str(resultat[2]) + "]")
            print("### ####\n")
            input()

        elif (bruker == 'e'):
            print("""
Du skal nå konvertere UTM-koordinater (N, E) til bredde- og lengdegrad.
Skriv inn verdiene du vil konvertere:
""")
            N, handling = giInnTall('Nord-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            E, handling = giInnTall('Øst-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            sone, handling = giInnTall('Sone')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""
            
            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break

            if (system == ""):
                resultat = koorkon.UTMToLatLon(float(N), float(E), sone+"N")
            else:
                resultat = koorkon.UTMToLatLon(float(N), float(E), sone+"N", system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + str(N) + ", " + str(E) + "] i sone " + sone)
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + "]")
            print("### ####\n")
            input()

        elif (bruker == 'f'):
            print("""
Du skal nå konvertere UTM-koordinater (N, E) til kartesiske koordinater.
Skriv inn verdiene du vil konvertere:
""")
            N, handling = giInnTall('Nord-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            E, handling = giInnTall('Øst-verdi')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            h, handling = giInnTall('Høyde')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break

            sone, handling = giInnTall('Sone')
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            valg, handling = giInnJaNei()
            if (handling == 'avslutt'):
                side, bruker = "hoved", "avslutt"
                break
            elif (handling == 'tilbake'):
                side = "Koordinatkonverterer"
                break
            
            system = ""
            
            if (valg == "J"):
                system, handling = giInnSystem()
                if (handling == 'avslutt'):
                    side, bruker = "hoved", "avslutt"
                    break
                elif (handling == 'tilbake'):
                    side = "Koordinatkonverterer"
                    break

            if (system == ""):
                resultat = koorkon.UTMToCartesian(float(N), float(E), float(h), sone + "N")
            else:
                resultat = koorkon.UTMToCartesian(float(N), float(E), float(h), sone + "N", system)
            
            print("\n### ####")
            print("Utgangspunkt: [" + str(N) + ", " + str(E) + "] i sone " + sone)
            print("Resultat: [" + str(resultat[0]) + ", " + str(resultat[1]) + "]")
            print("### ####\n")
            input()
    
    if (bruker != 'avslutt' and bruker != 'tilbake'):
        if (side == 'Koordinatkonverterer'):
            bruker = "2"

    if (bruker != "avslutt" and bruker != '2'):
        print("""
### ### ### ### ### ### ### ### ### ### ### ###

Du er nå tilbake i hovedmenyen og har følgende valg:
        
[1] Se svarene for 'Øving 1'
[2] Konvertere egne koordinater til andre system

Skriv et tall i følgende intervall, 1 - 2, for å velge funksjonalitet basert på valgmulighetene over.
""")
        
        bruker = input("Jeg vil utføre handling: ")
        while (bruker not in valgHoved):
            bruker = input("Ugyldig input! Utfører handling: ")

        print("\n### ### ### ### ### ### ### ### ### ### ### ###\n")

print("\n### ### ### ### ### ### ### ### ### ### ### ###\n\nTakk for denne gang! 8)")

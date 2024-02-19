"""

Program som kjører de faktiske oppgavene og konverteringsmulighetene

"""

import KoordinatKonvertering as koorkon

# Konstanter som holder oversikt over programmets status:

side = "hoved"
valgHoved = ['1', '2', 'avslutt']
valgKoor = ['1', '2', '3', '4', '5', '6', 'tilbake', 'avslutt']
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

# Program:

print("""
### ### ### ### ### ### ### ### ### ### ### ###

Hei og velkommen til koordinatkonvertereren!

### ### ### ### ### ### ### ### ### ### ### ###

Du er nå på hovedsiden der du har følgende valg:

[1] Se svarene for 'Øving 1'
[2] Konvertere egne koordinater til andre system

Du kan når som helst skrive 'avslutt' for å avslutte programmet, eller 'tilbake' for å komme tilbake til forrige meny.

Skriv et tall i følgende intervall, 1 - 2, for å velge funksjonalitet basert på valgmulighetene over.
""")

bruker = input("Jeg vil utføre handling: ")
while (bruker not in valgHoved):
    bruker = input("Ugyldig input! Utfører handling: ")

print("\n### ### ### ### ### ### ### ### ### ### ### ###")

while (bruker != 'avslutt'):
    if (bruker == '1'):
        side = "Oving 1"
    
    while (side == 'Oving 1'):
        print("\nDu valgte '1'.\n")
        bruker = input("Nå vil jeg se: ")
        
        ###################################
        # TODO: Legg inn funksjonalitet her
        ###################################

        side = avslutteEllerTilbake(bruker)
    
    if (bruker == '2'):
        side = "Koordinatkonverterer"
    
    ###################################
    # TODO: Les over og se om alt er riktig stavet, ønskede ord og syntaks 8)
    ###################################

    while (side == 'Koordinatkonverterer'):
        print("""\nDu valgte '2'.

Her kan du konvertere mellom ulike typer koordinater i gitte system.

Typene du kan velge mellom er: kartesiske koordinater, bredde- og lengdegrad, og UTM-koordinater.
Systemene du kan velge mellom er: EUREF89, WGS84, ED50/ED87 og NGO1948 (EUREF89 er forhåndsinnstilt)

Du har følgende valg:

[1] Kartesiske koordinater til bredde- og lengdegrad
[2] Kartesiske koordinater til UTM-koordinater
[3] Bredde- og lengdegrad til UTM-koordinater
[4] Bredde- og lengdegrad til kartesiske koordinater
[5] UTM-koordinater til bredde- og lengdegrad
[6] UTM-koordinater til kartesiske koordinater

Skriv inn ønsket handling som et tall i input-feltet.
        """)
        
        bruker = input("Jeg vil utføre handling: ")

        while (bruker not in valgKoor):
            print("\nUgyldig input! Prøv igjen")
            bruker = input("Jeg vil utføre handling: ")

        side = avslutteEllerTilbake(bruker)

        ###################################
        # TODO: Prøv å se om det går an å minske mengden kode ved å samle mye av det som går igjen i funksjoner
        ###################################

        if (bruker == '1'):
            print("""
Du skal nå konvertere kartesiske koordinater (x,y,z) til lengde- og breddegrad.
Skriv inn verdiene du vil konvertere:
""")
            x = input("x-verdi: ")
            side = avslutteEllerTilbake_2(x)
            if (side == 'hoved'):
                bruker = x
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(x)):
                x = input("Ikke gyldig tall! X-verdi: ")
                side = avslutteEllerTilbake_2(x)
                if (side == 'hoved'):
                    bruker = x
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            y = input("y-verdi: ")
            side = avslutteEllerTilbake_2(y)
            if (side == 'hoved'):
                bruker = y
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(y)):
                y = input("Ikke gyldig tall! Y-verdi: ")
                side = avslutteEllerTilbake_2(y)
                if (side == 'hoved'):
                    bruker = y
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            z = input("z-verdi: ")
            side = avslutteEllerTilbake_2(z)
            if (side == 'hoved'):
                bruker = z
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(z)):
                z = input("Ikke gyldig tall! Z-verdi: ")
                side = avslutteEllerTilbake_2(z)
                if (side == 'hoved'):
                    bruker = z
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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

        elif (bruker == '2'):
            print("""
Du skal nå konvertere kartesiske koordinater (x,y,z) til UTM-koordinater.
Skriv inn verdiene du vil konvertere:
""")
            x = input("x-verdi: ")
            side = avslutteEllerTilbake_2(x)
            if (side == 'hoved'):
                bruker = x
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(x)):
                x = input("Ikke gyldig tall! X-verdi: ")
                side = avslutteEllerTilbake_2(x)
                if (side == 'hoved'):
                    bruker = x
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            y = input("y-verdi: ")
            side = avslutteEllerTilbake_2(y)
            if (side == 'hoved'):
                bruker = y
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(y)):
                y = input("Ikke gyldig tall! Y-verdi: ")
                side = avslutteEllerTilbake_2(y)
                if (side == 'hoved'):
                    bruker = y
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            z = input("z-verdi: ")
            side = avslutteEllerTilbake_2(z)
            if (side == 'hoved'):
                bruker = z
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(z)):
                z = input("Ikke gyldig tall! Z-verdi: ")
                side = avslutteEllerTilbake_2(z)
                if (side == 'hoved'):
                    bruker = z
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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

        elif (bruker == '3'):
            print("""
Du skal nå konvertere bredde- og lengdegrader (B, L, h) til UTM-koordinater.
Skriv inn verdiene du vil konvertere:
""")
            B = input("Breddegrad: ")
            side = avslutteEllerTilbake_2(B)
            if (side == 'hoved'):
                bruker = B
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(B)):
                B = input("Ikke gyldig tall! Breddegrad: ")
                side = avslutteEllerTilbake_2(B)
                if (side == 'hoved'):
                    bruker = B
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            L = input("Lengdegrad: ")
            side = avslutteEllerTilbake_2(L)
            if (side == 'hoved'):
                bruker = L
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(L)):
                L = input("Ikke gyldig tall! Lengdegrad: ")
                side = avslutteEllerTilbake_2(L)
                if (side == 'hoved'):
                    bruker = L
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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

        elif (bruker == '4'):
            print("""
Du skal nå konvertere bredde- og lengdegrader (B, L, h) til kartesiske koordinater.
Skriv inn verdiene du vil konvertere:
""")
            B = input("Breddegrad: ")
            side = avslutteEllerTilbake_2(B)
            if (side == 'hoved'):
                bruker = B
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(B)):
                B = input("Ikke gyldig tall! Breddegrad: ")
                side = avslutteEllerTilbake_2(B)
                if (side == 'hoved'):
                    bruker = B
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            L = input("Lengdegrad: ")
            side = avslutteEllerTilbake_2(L)
            if (side == 'hoved'):
                bruker = L
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(L)):
                L = input("Ikke gyldig tall! Lengdegrad: ")
                side = avslutteEllerTilbake_2(L)
                if (side == 'hoved'):
                    bruker = L
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            h = input("Høyde: ")
            side = avslutteEllerTilbake_2(h)
            if (side == 'hoved'):
                bruker = h
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(h)):
                h = input("Ikke gyldig tall! Høyde: ")
                side = avslutteEllerTilbake_2(h)
                if (side == 'hoved'):
                    bruker = h
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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

        elif (bruker == '5'):
            print("""
Du skal nå konvertere UTM-koordinater (N, E) til bredde- og lengdegrad.
Skriv inn verdiene du vil konvertere:
""")
            N = input("Nord-verdi: ")
            side = avslutteEllerTilbake_2(N)
            if (side == 'hoved'):
                bruker = N
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(N)):
                N = input("Ikke gyldig tall! Nord-verdi: ")
                side = avslutteEllerTilbake_2(N)
                if (side == 'hoved'):
                    bruker = N
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            E = input("Øst-verdi: ")
            side = avslutteEllerTilbake_2(E)
            if (side == 'hoved'):
                bruker = E
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(E)):
                E = input("Ikke gyldig tall! Øst-verdi: ")
                side = avslutteEllerTilbake_2(E)
                if (side == 'hoved'):
                    bruker = E
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            sone = input("Hvilken sone er du i: ")
            side = avslutteEllerTilbake_2(sone)
            if (side == 'hoved'):
                bruker = sone
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(sone)):
                sone = input("Ikke gyldig tall! Sone: ")
                side = avslutteEllerTilbake_2(sone)
                if (side == 'hoved'):
                    bruker = sone
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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

        elif (bruker == '6'):
            print("""
Du skal nå konvertere UTM-koordinater (N, E) til kartesiske koordinater.
Skriv inn verdiene du vil konvertere:
""")
            N = input("Nord-verdi: ")
            side = avslutteEllerTilbake_2(N)
            if (side == 'hoved'):
                bruker = N
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(N)):
                N = input("Ikke gyldig tall! Nord-verdi: ")
                side = avslutteEllerTilbake_2(N)
                if (side == 'hoved'):
                    bruker = N
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            E = input("Øst-verdi: ")
            side = avslutteEllerTilbake_2(E)
            if (side == 'hoved'):
                bruker = E
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(E)):
                E = input("Ikke gyldig tall! Øst-verdi: ")
                side = avslutteEllerTilbake_2(E)
                if (side == 'hoved'):
                    bruker = E
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            h = input("Høyde: ")
            side = avslutteEllerTilbake_2(h)
            if (side == 'hoved'):
                bruker = h
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(h)):
                h = input("Ikke gyldig tall! Høyde: ")
                side = avslutteEllerTilbake_2(h)
                if (side == 'hoved'):
                    bruker = h
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break

            sone = input("Hvilken sone er du i: ")
            side = avslutteEllerTilbake_2(sone)
            if (side == 'hoved'):
                bruker = sone
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(sone)):
                sone = input("Ikke gyldig tall! Sone: ")
                side = avslutteEllerTilbake_2(sone)
                if (side == 'hoved'):
                    bruker = sone
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            valg = input("Ønsker du et annet system enn EUREF89? (J/N) ")
            side = avslutteEllerTilbake_2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake_2(valg)
                if (side == 'hoved'):
                    bruker = valg
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer_exit"
                    break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
                side = "Koordinatkonverterer"
                break
            
            if (valg == "J"):
                print("Velg mellom følgende system: 'WGS84', 'ED50'/'ED87' eller 'NGO1948'.")
                system = input("System: ")
                side = avslutteEllerTilbake_2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake_2(system)
                    if (side == 'hoved'):
                        bruker = system
                        break
                    elif (side == "tilbake"):
                        side = "Koordinatkonverterer_exit"
                        break
            if (side == 'hoved'):
                break
            elif (side == "Koordinatkonverterer_exit"):
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
    
    if (side == 'Koordinatkonverterer'):
        bruker = "2"

    elif (bruker != "avslutt"):
        print("""
### ### ### ### ### ### ### ### ### ### ### ###

Du er nå tilbake i hovedmenyen og har følgende valg:
        
[1] Se svarene for 'Øving 1'
[2] Konvertere egne koordinater til andre system

Skriv et tall i følgende intervall, 1 - 2, for å velge funksjonalitet basert på valgmulighetene over.
        """)
        bruker = input("Jeg vil utføre handling: ")
        print("\n### ### ### ### ### ### ### ### ### ### ### ###\n")

print("\n### ### ### ### ### ### ### ### ### ### ### ###\n\nTakk for denne gang! 8)")
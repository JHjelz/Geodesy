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

def avslutteEllerTilbake2(bruker):
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
        
        side = avslutteEllerTilbake(bruker)
    
    if (bruker == '2'):
        side = "Koordinatkonverterer"
    while (side == 'Koordinatkonverterer'):
        print("""\nDu valgte '2'.

Her kan du konvertere mellom ulike typer koordinater i gitte system.

Typene du kan velge mellom er: kartesiske koordinater, lengde- og breddegrader, og UTM-koordinater.
Systemene du kan velge mellom er: EUREF89, WGS84, ED50/ED87 og NGO1948 (EUREF89 er forhåndsinnstilt)

Du har følgende valg:

[1] Kartesiske koordinater til lengde- og breddegrad
[2] Kartesiske koordinater til UTM-koordinater
[3] Lengde- og breddegrad til UTM-koordinater
[4] Lengde- og breddegrad til kartesiske koordinater
[5] UTM-koordinater til lengde- og breddegrad
[6] UTM-koordinater til kartesiske koordinater

Skriv inn ønsket handling som et tall i input-feltet.
        """)
        
        bruker = input("Jeg vil utføre handling: ")

        while (bruker not in valgKoor):
            print("\nUgyldig input! Prøv igjen")
            bruker = input("Jeg vil utføre handling: ")

        side = avslutteEllerTilbake(bruker)

        if (bruker == '1'):
            print("""
Du skal nå konvertere kartesiske koordinater (x,y,z) til lengde- og breddegrad.
Skriv inn verdiene du vil konvertere:
""")
            x = input("x-verdi: ")
            side = avslutteEllerTilbake2(x)
            if (side == 'hoved'):
                bruker = x
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(x)):
                x = input("Ikke gyldig tall! X-verdi: ")
                side = avslutteEllerTilbake2(x)
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
            side = avslutteEllerTilbake2(y)
            if (side == 'hoved'):
                bruker = y
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(y)):
                y = input("Ikke gyldig tall! Y-verdi: ")
                side = avslutteEllerTilbake2(y)
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
            side = avslutteEllerTilbake2(z)
            if (side == 'hoved'):
                bruker = z
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            while (not gyldigTall(z)):
                z = input("Ikke gyldig tall! Z-verdi: ")
                side = avslutteEllerTilbake2(z)
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
            side = avslutteEllerTilbake2(valg)
            if (side == 'hoved'):
                bruker = valg
                break
            elif (side == "tilbake"):
                side = "Koordinatkonverterer"
                break
            system = ""
            while (valg not in jaNei):
                valg = input("J/N: ")
                side = avslutteEllerTilbake2(valg)
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
                side = avslutteEllerTilbake2(system)
                if (side == 'hoved'):
                    bruker = system
                    break
                elif (side == "tilbake"):
                    side = "Koordinatkonverterer"
                    break
                while (system not in systemvalg):
                    system = input("Ugyldig system! System: ")
                    side = avslutteEllerTilbake2(system)
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
            continue
        elif (bruker == '3'):
            continue
        elif (bruker == '4'):
            continue
        elif (bruker == '5'):
            continue
        elif (bruker == '6'):
            continue
    
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
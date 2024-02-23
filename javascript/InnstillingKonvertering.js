// Status for sida:

var system = "EUREF89";
var type = "kartesisk";

// Funksjoner for å endre status:

function endreSystem(id) {
    document.getElementById(system).style.backgroundColor = "green";
    document.getElementById(id).style.backgroundColor = "orangered";
    system = id;
}

function endreType(id) {
    document.getElementById(type).style.backgroundColor = "green";
    document.getElementById(id).style.backgroundColor = "orangered";
    type = id;
    
    if (type == "kartesisk") {
        document.getElementById("in1").placeholder = "X-verdi";
        document.getElementById("in2").placeholder = "Y-verdi";
        document.getElementById("in3").placeholder = "Z-verdi";
        document.getElementById("in4").style.display = "none";
        fjernInput()
    } else if (type == "brelen") {
        document.getElementById("in1").placeholder = "Breddegrad";
        document.getElementById("in2").placeholder = "Lengdegrad";
        document.getElementById("in3").placeholder = "Hoyde";
        document.getElementById("in4").style.display = "none";
        fjernInput();
    } else if (type == "UTM") {
        document.getElementById("in1").placeholder = "Nord-verdi";
        document.getElementById("in2").placeholder = "Ost-verdi";
        document.getElementById("in3").placeholder = "Hoyde";
        document.getElementById("in4").style.display = "inline-block";
        document.getElementById("in4").value = "";
        fjernInput();
    }
}

function fjernInput() {
    document.getElementById("in1").value = "";
    document.getElementById("in2").value = "";
    document.getElementById("in3").value = "";
    document.getElementById("resultat").innerHTML = "";
    fjernPunkt();
}

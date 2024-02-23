function reset() {
    endreSystem('EUREF89');
    endreType('kartesisk');
    document.getElementById("in1").value = "";
    document.getElementById("in2").value = "";
    document.getElementById("in3").value = "";
    document.getElementById("in4").value = "";
    document.getElementById("resultat").innerHTML = "";
    fjernPunkt();
    map.setView([63.418529, 10.40284], 13)
}

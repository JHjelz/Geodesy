function transform_1()  {
    x = parseFloat(document.getElementById("X"));
    y = parseFloat(document.getElementById("Y"));
    z = parseFloat(document.getElementById("Z"));

    var resultat = cartesianToLatLong(x, y, z);

    document.getElementById("Resultat").innerHTML("[" + resultat[0] + ", " + resultat[1] + ", " + resultat[2] + "]");
}
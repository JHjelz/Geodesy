/*

Funksjoner som konverterer koordinater mellom de ulike typene

*/

function konverter() {
    var val = hentVerdier(system);
    if (type == "kartesisk") {
        Cartesian(val);
    } else if (type == "brelen") {
        LatLong(val);
    } else if (type == "UTM") {
        UTM(val);
    }
}

function Cartesian(val) {
    var x = parseFloat(document.getElementById("in1").value);
    var y = parseFloat(document.getElementById("in2").value);
    var z = parseFloat(document.getElementById("in3").value);

    var latlon = cartesianToLatLong(x, y, z, val[0], val[2]);
    var UTM = cartesianToUTM(x, y, z, val[0], val[1], val[2]);

    plasserPunkt(latlon[0], latlon[1]);

    document.getElementById("resultat").innerHTML = `<b>Bredde- og lengdegrader:</b><br>[${latlon[0]}, ${latlon[1]}, ${latlon[2]}]<br><b>UTM-koordinater:</b><br>[${UTM[0]}, ${UTM[1]}, ${UTM[2]}]`;
}

function LatLong(val) {
    var B = parseFloat(document.getElementById("in1").value);
    var L = parseFloat(document.getElementById("in2").value);
    var h = parseFloat(document.getElementById("in3").value);

    var kartesisk = latlonToCartesian(B, L, h, val[0], val[1]);
    var UTM = latLonToUTM(B, L, val[0], val[1], val[2]);

    plasserPunkt(B, L);

    document.getElementById("resultat").innerHTML = `<b>Kartesiske koordinater:</b><br>[${kartesisk[0]}, ${kartesisk[1]}, ${kartesisk[2]}]<br><b>UTM-koordinater:</b><br>[${UTM[0]}, ${UTM[1]}, ${h}]`;
}

function UTM(val) {
    var N = parseFloat(document.getElementById("in1").value);
    var E = parseFloat(document.getElementById("in2").value);
    var h = parseFloat(document.getElementById("in3").value);
    var sone = parseFloat(document.getElementById("in4").value);

    var kartesisk = UTMToCartesian(N, E, h, sone, val[0], val[1], val[2]);
    var latlon = UTMToLatLon(N, E, sone, val[0], val[2]);

    plasserPunkt(latlon[0], latlon[1]);

    document.getElementById("resultat").innerHTML = `<b>Kartesiske koordinater:</b><br>[${kartesisk[0]}, ${kartesisk[1]}, ${kartesisk[2]}]<br><b>UTM-koordinater:</b><br>[${latlon[0]}, ${latlon[1]}, ${h}]`;
}

// Hjelpefunksjoner:

function gradTilrad(v) {
    return v * Math.PI / 180;
}

function radTilgrad(v) {
    return v * 180 / Math.PI;
}

function A(a, n) {
    return a * (1 - n + (5 / 4) * (n**2 - n**3) + (81 / 64) * (n**4 - n**5));
}

function B(a, n) {
    return (3 / 2) * a * (n - n**2 + (7 / 8) * (n**3 - n**4) + (55 / 64) * n**5);
}

function C(a, n) {
    return (15 / 16) * a * (n**2 - n**3 + (3 / 4) * (n**4 - n**5));
}

function D(a, n) {
    return (35 / 48) * a * (n**3 - n**4 + (11 / 16) * n**5);
}

function E(a, n) {
    return (315 / 512) * a * (n**4 - n**5);
}

// Hovedfunksjoner:

function cartesianToLatLong(x, y, z, a, e) {
    var L = Math.atan(y / x);
    var p = Math.sqrt(x**2 + y**2);

    var oldB = 0;
    var newB = Math.atan(z / ((1 - e**2) * p))

    while (Math.abs(newB - oldB) > 10**(-12)) {
        var n = a / Math.sqrt(1 - e**2 * (Math.sin(newB))**2);
        var h = p / (Math.cos(newB)) - n;
        oldB = newB;
        newB = Math.atan(z / ((1 - e**2 * n / (n + h)) * p));
    }

    return [radTilgrad(newB), radTilgrad(L), h];
}

function latLonToUTM(br, le, a, b, e) {
    var l0 = gradTilrad((Math.floor(le / 6) * 6) + 3);
    
    br = gradTilrad(br);
    le = gradTilrad(le);
    
    var k0 = 0.9996;
    var e2 = e**2 / (1 - e**2);
    var n = (a - b) / (a + b);
    var nu = a / Math.sqrt(1 - e**2 * (Math.sin(br))**2);
    var p = le - l0;

    var S = A(a, n) * br - B(a, n) * Math.sin(2 * br) + C(a, n) * Math.sin(4 * br) - D(a, n) * Math.sin(6 * br) + E(a, n) * Math.sin(8 * br);

    var k1 = S * k0;
    var k2 = k0 * nu * Math.sin(2 * br) / 4;
    var k3 = (k0 * nu * Math.sin(br) * (Math.cos(br))**3 / 24) * (5 - (Math.tan(br))**2 + 9 * e2 * (Math.cos(br))**2 + 4 * e2**2 * (Math.cos(br))**4);
    var k4 = k0 * nu * Math.cos(br);
    var k5 = (k0 * nu * (Math.cos(br))**3 / 6) * (1  - (Math.tan(br))**2 + e2 * (Math.cos(br))**2);

    var x = k1 + k2 * p**2 + k3 * p**4;
    var y = k4 * p + k5 * p**3 + 500000;

    return [x, y];
}

function cartesianToUTM(x, y, z, a, b, e) {
    var A1 = cartesianToLatLong(x, y, z, a, e);
    var A2 = latLonToUTM(A1[0], A1[1], a, b, e);
    return [A2[0], A2[1], A1[2]];
}

function latlonToCartesian(br, le, h, a, b) {
    br = gradTilrad(br);
    le = gradTilrad(le);

    var n = a**2 / Math.sqrt(a**2 * (Math.cos(br))**2 + b**2 * (Math.sin(br))**2);

    var x = (n + h) * Math.cos(br) * Math.cos(le);
    var y = (n + h) * Math.cos(br) * Math.sin(le);
    var z = ((b**2 /  a**2) * n + h) * Math.sin(br);

    return [x, y, z];
}

function UTMToLatLon(x, y, sone, a, e) {
    y = y - 500000;

    var k0 = 0.9996;
    var m = x / k0;
    var mu  = m / (a * (1 - e**2 / 4 - 3 * e**4 / 64 - 5 * e**6 / 256));
    var e1 = (1 - Math.sqrt(1 - e**2)) / (1 + Math.sqrt(1 - e**2));

    var J1 = 3 * e1 / 2 - 27 * e1**3 / 32;
    var J2 = 21 * e1**2 / 16 - 55 * e1**4 / 32;
    var J3 = 151 * e1**3 / 96;
    var J4 = 1097  * e1**4 / 512;

    var fp = mu + J1 * Math.sin(2 * mu) + J2 * Math.sin(4 * mu) + J3 * Math.sin(6 * mu) +  J4 * Math.sin(8 * mu);

    var e2 = e**2 / (1 - e**2);
    var C1 = e2 * (Math.cos(fp))**2;
    var T1 = (Math.tan(fp))**2;
    var R1 = a * (1 - e**2) / (1 - e**2 * (Math.sin(fp))**2)**(3/2);
    var N1 = a / Math.sqrt(1 - e**2 * (Math.sin(fp))**2);
    var D = y / (N1 * k0);

    var Q1 = N1 * Math.tan(fp) / R1;
    var Q2 = D**2 / 2;
    var Q3 = (5 + 3 * T1 + 10 * C1 - 4 * C1**2 - 9 * e2) * D**4 / 24;
    var Q4 = (61 + 90 * T1 + 298 * C1 + 45 * T1**2 - 3 * C1**2 - 252 * e2) * D**6 / 720;
    var Q5 = D;
    var Q6 = (1 + 2 * T1 + C1) * D**3 / 6;
    var Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1**2 + 8 * e2 + 24 * T1**2) * D**5 / 120;

    if (sone == 32) {
        l0 = 9;
    } else if (sone == 33) {
        l0 = 15;
    }

    var br = radTilgrad(fp - Q1 * (Q2 - Q3 + Q4));
    var le = l0 + radTilgrad((Q5 - Q6 + Q7) / Math.cos(fp));

    return [br, le];
}

function UTMToCartesian(x, y, h, sone, a, b, e) {
    var A1 = UTMToLatLon(x, y, sone, a, e);
    return latlonToCartesian(A1[0], A1[1], h, a, b);
}

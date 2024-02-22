function hentVerdier(system) {
    if (system == "EUREF89") {
        // EUREF89:
        var a = 6378137; // [m]
        var b = 6356752.3141; // [m]
    } else if (system ==  "WGS84") {
        // WGS84:
        var a = 6378137; // [m]
        var b = 6356752.3142; // [m]
    } else if (system == "ED50") {
        // ED50 / ED87:
        var a = 6378388; // [m]
        var b = 6356911.9461; // [m]
    } else if (system == "NGO1948") {
        // NGO1948:
        var a = 6377492.0176; // [m]
        var b = 6356173.5083; // [m]
    }

    var e = Math.sqrt((a**2 - b**2) / a**2); // []

    return [a, b, e]
}
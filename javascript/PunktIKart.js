var marker = L.marker();

function plasserPunkt(B, L) {
    marker.setLatLng([B, L]).addTo(map);
    map.setView([B, L], 8);
}

function fjernPunkt() {
    map.removeLayer(marker);
}
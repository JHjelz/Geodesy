// Initialiserer kartet:
const map = L.map('map', {
    maxZoom: 18, // Justerer zoom-nivået
    minZoom: 6,
    zoomControl: false, // Fjerner default zoom-knapper
}).setView([63.418529, 10.40284], 13)

// Bestemme URL til basiskartet:
var osm_map = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright"<OpenStreetMap</a> contributors'
});

osm_map.addTo(map); // Legger basiskartet til i kartet slik at det er det som vises først på nettsiden

// Layer control er lagd i 'leafletLayerControl.js'

// Legger til zoom-knapper på egnet sted
L.control.zoom({
    position:'topright'
}).addTo(map);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Ordnance Survey Grid Reference functions                           (c) Chris Veness 2005-2017  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-gridref.html                                            */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-osgridref.html                              */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

'use strict';
if (typeof module != 'undefined' && module.exports) var LatLon = require('./latlon-ellipsoidal.js'); // ≡ import LatLon from 'latlon-ellipsoidal.js'


/**
 * Convert OS grid references to/from OSGB latitude/longitude points.
 *
 * Formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is inferior
 * to Krüger as used by e.g. Karney 2011.
 *
 * www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf.
 *
 * @module   osgridref
 * @requires latlon-ellipsoidal
 */
/*
 * Converted 2015 to work with WGS84 by default, OSGB36 as option;
 * www.ordnancesurvey.co.uk/blog/2014/12/confirmation-on-changes-to-latitude-and-longitude
 */


/**
 * Creates an OsGridRef object.
 *
 * @constructor
 * @param {number} easting - Easting in metres from OS false origin.
 * @param {number} northing - Northing in metres from OS false origin.
 *
 * @example
 *   var grid = new OsGridRef(651409, 313177);
 */
function OsGridRef(easting, northing) {
    // allow instantiation without 'new'
    if (!(this instanceof OsGridRef)) return new OsGridRef(easting, northing);

    this.easting = Number(easting);
    this.northing = Number(northing);
}


/**
 * Converts latitude/longitude to Ordnance Survey grid reference easting/northing coordinate.
 *
 * Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
 * inferior to Krüger as used by e.g. Karney 2011.
 *
 * @param   {LatLon}    point - latitude/longitude.
 * @returns {OsGridRef} OS Grid Reference easting/northing.
 *
 * @example
 *   var p = new LatLon(52.65798, 1.71605);
 *   var grid = OsGridRef.latLonToOsGrid(p); // grid.toString(): TG 51409 13177
 *   // for conversion of (historical) OSGB36 latitude/longitude point:
 *   var p = new LatLon(52.65757, 1.71791, LatLon.datum.OSGB36);
 */
OsGridRef.latLonToOsGrid = function (point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');

    // if necessary convert to OSGB36 first
    if (point.datum != LatLon.datum.OSGB36) point = point.convertDatum(LatLon.datum.OSGB36);

    var φ = point.lat.toRadians();
    var λ = point.lon.toRadians();

    console.log(`φ is ${φ}`);
    console.log(`λ is ${λ}`);

    var a = 6377563.396, b = 6356256.909;              // Airy 1830 major & minor semi-axes
    var F0 = 0.9996012717;                             // NatGrid scale factor on central meridian
    var φ0 = (49).toRadians(), λ0 = (-2).toRadians();  // NatGrid true origin is 49°N,2°W
    var N0 = -100000, E0 = 400000;                     // northing & easting of true origin, metres
    var e2 = 1 - (b * b) / (a * a);                          // eccentricity squared
    var n = (a - b) / (a + b), n2 = n * n, n3 = n * n * n;         // n, n², n³

    var cosφ = Math.cos(φ), sinφ = Math.sin(φ);
    var ν = a * F0 / Math.sqrt(1 - e2 * sinφ * sinφ);            // nu = transverse radius of curvature
    var ρ = a * F0 * (1 - e2) / Math.pow(1 - e2 * sinφ * sinφ, 1.5); // rho = meridional radius of curvature
    var η2 = ν / ρ - 1;                                    // eta = ?

    var Ma = (1 + n + (5 / 4) * n2 + (5 / 4) * n3) * (φ - φ0);
    var Mb = (3 * n + 3 * n * n + (21 / 8) * n3) * Math.sin(φ - φ0) * Math.cos(φ + φ0);
    var Mc = ((15 / 8) * n2 + (15 / 8) * n3) * Math.sin(2 * (φ - φ0)) * Math.cos(2 * (φ + φ0));
    var Md = (35 / 24) * n3 * Math.sin(3 * (φ - φ0)) * Math.cos(3 * (φ + φ0));
    var M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

    console.log(`φ is ${φ}`);
    console.log(`φ0 is ${φ0}`);

    console.log(`φ - φ0 is ${φ - φ0}`);
    console.log(`φ + φ0 is ${φ + φ0}`);
    console.log(`Ma is ${Ma}`);
    console.log(`Ma is ${Ma}`);
    console.log(`Mb is ${Mb}`);
    console.log(`Mc is ${Mc}`);
    console.log(`Md is ${Md}`);

    var cos3φ = cosφ * cosφ * cosφ;
    var cos5φ = cos3φ * cosφ * cosφ;
    var tan2φ = Math.tan(φ) * Math.tan(φ);
    var tan4φ = tan2φ * tan2φ;

    var I = M + N0;
    console.log(`M is ${M}`);
    console.log(`N0 is ${N0}`);
    console.log(`I is ${I}`);
    var II = (ν / 2) * sinφ * cosφ;
    console.log(`II is ${II}`);
    var III = (ν / 24) * sinφ * cos3φ * (5 - tan2φ + 9 * η2);
    console.log(`III is ${III}`);
    var IIIA = (ν / 720) * sinφ * cos5φ * (61 - 58 * tan2φ + tan4φ);
    console.log(`IIIA is ${IIIA}`);
    var IV = ν * cosφ;
    console.log(`IV is ${IV}`);
    var V = (ν / 6) * cos3φ * (ν / ρ - tan2φ);
    console.log(`V is ${V}`);
    var VI = (ν / 120) * cos5φ * (5 - 18 * tan2φ + tan4φ + 14 * η2 - 58 * tan2φ * η2);
    console.log(`VI is ${VI}`);

    var Δλ = λ - λ0;
    var Δλ2 = Δλ * Δλ, Δλ3 = Δλ2 * Δλ, Δλ4 = Δλ3 * Δλ, Δλ5 = Δλ4 * Δλ, Δλ6 = Δλ5 * Δλ;

    var N = I + II * Δλ2 + III * Δλ4 + IIIA * Δλ6;
    var E = E0 + IV * Δλ + V * Δλ3 + VI * Δλ5;

    // N = Number(N.toFixed(3)); // round to mm precision
    // E = Number(E.toFixed(3));

    return new OsGridRef(E, N); // gets truncated to SW corner of 1m grid square
};


/**
 * Converts Ordnance Survey grid reference easting/northing coordinate to latitude/longitude
 * (SW corner of grid square).
 *
 * Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
 * inferior to Krüger as used by e.g. Karney 2011.
 *
 * @param   {OsGridRef}    gridref - Grid ref E/N to be converted to lat/long (SW corner of grid square).
 * @param   {LatLon.datum} [datum=WGS84] - Datum to convert grid reference into.
 * @returns {LatLon}       Latitude/longitude of supplied grid reference.
 *
 * @example
 *   var gridref = new OsGridRef(651409.903, 313177.270);
 *   var pWgs84 = OsGridRef.osGridToLatLon(gridref);                     // 52°39′28.723″N, 001°42′57.787″E
 *   // to obtain (historical) OSGB36 latitude/longitude point:
 *   var pOsgb = OsGridRef.osGridToLatLon(gridref, LatLon.datum.OSGB36); // 52°39′27.253″N, 001°43′04.518″E
 */
OsGridRef.osGridToLatLon = function (gridref, datum) {
    if (!(gridref instanceof OsGridRef)) throw new TypeError('gridref is not OsGridRef object');
    if (datum === undefined) datum = LatLon.datum.WGS84;

    var E = gridref.easting;
    var N = gridref.northing;

    var a = 6377563.396, b = 6356256.909;              // Airy 1830 major & minor semi-axes
    var F0 = 0.9996012717;                             // NatGrid scale factor on central meridian
    var φ0 = (49).toRadians(), λ0 = (-2).toRadians();  // NatGrid true origin is 49°N,2°W
    var N0 = -100000, E0 = 400000;                     // northing & easting of true origin, metres
    var e2 = 1 - (b * b) / (a * a);                          // eccentricity squared
    var n = (a - b) / (a + b), n2 = n * n, n3 = n * n * n;         // n, n², n³

    var φ = φ0, M = 0;
    do {
        φ = (N - N0 - M) / (a * F0) + φ;

        var Ma = (1 + n + (5 / 4) * n2 + (5 / 4) * n3) * (φ - φ0);
        var Mb = (3 * n + 3 * n * n + (21 / 8) * n3) * Math.sin(φ - φ0) * Math.cos(φ + φ0);
        var Mc = ((15 / 8) * n2 + (15 / 8) * n3) * Math.sin(2 * (φ - φ0)) * Math.cos(2 * (φ + φ0));
        var Md = (35 / 24) * n3 * Math.sin(3 * (φ - φ0)) * Math.cos(3 * (φ + φ0));
        M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

    } while (N - N0 - M >= 0.00001);  // ie until < 0.01mm

    var cosφ = Math.cos(φ), sinφ = Math.sin(φ);
    var ν = a * F0 / Math.sqrt(1 - e2 * sinφ * sinφ);            // nu = transverse radius of curvature
    var ρ = a * F0 * (1 - e2) / Math.pow(1 - e2 * sinφ * sinφ, 1.5); // rho = meridional radius of curvature
    var η2 = ν / ρ - 1;                                    // eta = ?

    var tanφ = Math.tan(φ);
    var tan2φ = tanφ * tanφ, tan4φ = tan2φ * tan2φ, tan6φ = tan4φ * tan2φ;
    var secφ = 1 / cosφ;
    var ν3 = ν * ν * ν, ν5 = ν3 * ν * ν, ν7 = ν5 * ν * ν;
    var VII = tanφ / (2 * ρ * ν);
    var VIII = tanφ / (24 * ρ * ν3) * (5 + 3 * tan2φ + η2 - 9 * tan2φ * η2);
    var IX = tanφ / (720 * ρ * ν5) * (61 + 90 * tan2φ + 45 * tan4φ);
    var X = secφ / ν;
    var XI = secφ / (6 * ν3) * (ν / ρ + 2 * tan2φ);
    var XII = secφ / (120 * ν5) * (5 + 28 * tan2φ + 24 * tan4φ);
    var XIIA = secφ / (5040 * ν7) * (61 + 662 * tan2φ + 1320 * tan4φ + 720 * tan6φ);

    var dE = (E - E0), dE2 = dE * dE, dE3 = dE2 * dE, dE4 = dE2 * dE2, dE5 = dE3 * dE2, dE6 = dE4 * dE2, dE7 = dE5 * dE2;
    φ = φ - VII * dE2 + VIII * dE4 - IX * dE6;
    var λ = λ0 + X * dE - XI * dE3 + XII * dE5 - XIIA * dE7;

    var point = new LatLon(φ.toDegrees(), λ.toDegrees(), LatLon.datum.OSGB36);
    if (datum != LatLon.datum.OSGB36) point = point.convertDatum(datum);

    return point;
};


/**
 * Parses grid reference to OsGridRef object.
 *
 * Accepts standard grid references (eg 'SU 387 148'), with or without whitespace separators, from
 * two-digit references up to 10-digit references (1m × 1m square), or fully numeric comma-separated
 * references in metres (eg '438700,114800').
 *
 * @param   {string}    gridref - Standard format OS grid reference.
 * @returns {OsGridRef} Numeric version of grid reference in metres from false origin (SW corner of
 *   supplied grid square).
 * @throws Error on Invalid grid reference.
 *
 * @example
 *   var grid = OsGridRef.parse('TG 51409 13177'); // grid: { easting: 651409, northing: 313177 }
 */
OsGridRef.parse = function (gridref) {
    gridref = String(gridref).trim();

    // check for fully numeric comma-separated gridref format
    var match = gridref.match(/^(\d+),\s*(\d+)$/);
    if (match) return new OsGridRef(match[1], match[2]);

    // validate format
    match = gridref.match(/^[A-Z]{2}\s*[0-9]+\s*[0-9]+$/i);
    if (!match) throw new Error('Invalid grid reference');

    // get numeric values of letter references, mapping A->0, B->1, C->2, etc:
    var l1 = gridref.toUpperCase().charCodeAt(0) - 'A'.charCodeAt(0);
    var l2 = gridref.toUpperCase().charCodeAt(1) - 'A'.charCodeAt(0);
    // shuffle down letters after 'I' since 'I' is not used in grid:
    if (l1 > 7) l1--;
    if (l2 > 7) l2--;

    // convert grid letters into 100km-square indexes from false origin (grid square SV):
    var e100km = ((l1 - 2) % 5) * 5 + (l2 % 5);
    var n100km = (19 - Math.floor(l1 / 5) * 5) - Math.floor(l2 / 5);

    // skip grid letters to get numeric (easting/northing) part of ref
    var en = gridref.slice(2).trim().split(/\s+/);
    // if e/n not whitespace separated, split half way
    if (en.length == 1) en = [en[0].slice(0, en[0].length / 2), en[0].slice(en[0].length / 2)];

    // validation
    if (e100km < 0 || e100km > 6 || n100km < 0 || n100km > 12) throw new Error('Invalid grid reference');
    if (en.length != 2) throw new Error('Invalid grid reference');
    if (en[0].length != en[1].length) throw new Error('Invalid grid reference');

    // standardise to 10-digit refs (metres)
    en[0] = (en[0] + '00000').slice(0, 5);
    en[1] = (en[1] + '00000').slice(0, 5);

    var e = e100km + en[0];
    var n = n100km + en[1];

    return new OsGridRef(e, n);
};


/**
 * Converts ‘this’ numeric grid reference to standard OS grid reference.
 *
 * @param   {number} [digits=10] - Precision of returned grid reference (10 digits = metres);
 *   digits=0 will return grid reference in numeric format.
 * @returns {string} This grid reference in standard format.
 *
 * @example
 *   var ref = new OsGridRef(651409, 313177).toString(); // TG 51409 13177
 */
OsGridRef.prototype.toString = function (digits) {
    digits = (digits === undefined) ? 10 : Number(digits);
    if (isNaN(digits) || digits % 2 != 0 || digits > 16) throw new RangeError('Invalid precision ‘' + digits + '’');

    var e = this.easting;
    var n = this.northing;
    if (isNaN(e) || isNaN(n)) throw new Error('Invalid grid reference');

    // use digits = 0 to return numeric format (in metres, allowing for decimals & for northing > 1e6)
    if (digits == 0) {
        var eInt = Math.floor(e), eDec = e - eInt;
        var nInt = Math.floor(n), nDec = n - nInt;
        var ePad = ('000000' + eInt).slice(-6) + (eDec > 0 ? eDec.toFixed(3).slice(1) : '');
        var nPad = (nInt < 1e6 ? ('000000' + nInt).slice(-6) : nInt) + (nDec > 0 ? nDec.toFixed(3).slice(1) : '');
        return ePad + ',' + nPad;
    }

    // get the 100km-grid indices
    var e100k = Math.floor(e / 100000), n100k = Math.floor(n / 100000);

    if (e100k < 0 || e100k > 6 || n100k < 0 || n100k > 12) return '';

    // translate those into numeric equivalents of the grid letters
    var l1 = (19 - n100k) - (19 - n100k) % 5 + Math.floor((e100k + 10) / 5);
    var l2 = (19 - n100k) * 5 % 25 + e100k % 5;

    // compensate for skipped 'I' and build grid letter-pairs
    if (l1 > 7) l1++;
    if (l2 > 7) l2++;
    var letterPair = String.fromCharCode(l1 + 'A'.charCodeAt(0), l2 + 'A'.charCodeAt(0));

    // strip 100km-grid indices from easting & northing, and reduce precision
    e = Math.floor((e % 100000) / Math.pow(10, 5 - digits / 2));
    n = Math.floor((n % 100000) / Math.pow(10, 5 - digits / 2));

    // pad eastings & northings with leading zeros (just in case, allow up to 16-digit (mm) refs)
    e = ('00000000' + e).slice(-digits / 2);
    n = ('00000000' + n).slice(-digits / 2);

    return letterPair + ' ' + e + ' ' + n;
};


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
if (typeof module != 'undefined' && module.exports) module.exports = OsGridRef; // ≡ export default OsGridRef
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Geodesy tools for an ellipsoidal earth model                       (c) Chris Veness 2005-2016  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-convert-coords.html                                     */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-latlon-ellipsoidal.html                     */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

'use strict';
if (typeof module != 'undefined' && module.exports) var Vector3d = require('./vector3d.js'); // ≡ import Vector3d from 'vector3d.js'
if (typeof module != 'undefined' && module.exports) var Dms = require('./dms.js');           // ≡ import Dms from 'dms.js'


/**
 * Library of geodesy functions for operations on an ellipsoidal earth model.
 *
 * Includes ellipsoid parameters and datums for different coordinate systems, and methods for
 * converting between them and to cartesian coordinates.
 *
 * q.v. Ordnance Survey ‘A guide to coordinate systems in Great Britain’ Section 6
 * www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf.
 *
 * @module   latlon-ellipsoidal
 * @requires dms
 */


/**
 * Creates lat/lon (polar) point with latitude & longitude values, on a specified datum.
 *
 * @constructor
 * @param {number}       lat - Geodetic latitude in degrees.
 * @param {number}       lon - Longitude in degrees.
 * @param {LatLon.datum} [datum=WGS84] - Datum this point is defined within.
 *
 * @example
 *     var p1 = new LatLon(51.4778, -0.0016, LatLon.datum.WGS84);
 */
function LatLon(lat, lon, datum) {
    // allow instantiation without 'new'
    if (!(this instanceof LatLon)) return new LatLon(lat, lon, datum);

    if (datum === undefined) datum = LatLon.datum.WGS84;

    this.lat = Number(lat);
    this.lon = Number(lon);
    this.datum = datum;
}


/**
 * Ellipsoid parameters; major axis (a), minor axis (b), and flattening (f) for each ellipsoid.
 */
LatLon.ellipsoid = {
    WGS84: { a: 6378137, b: 6356752.314245, f: 1 / 298.257223563 },
    Airy1830: { a: 6377563.396, b: 6356256.909, f: 1 / 299.3249646 },
    AiryModified: { a: 6377340.189, b: 6356034.448, f: 1 / 299.3249646 },
    Bessel1841: { a: 6377397.155, b: 6356078.962818, f: 1 / 299.1528128 },
    Clarke1866: { a: 6378206.4, b: 6356583.8, f: 1 / 294.978698214 },
    Clarke1880IGN: { a: 6378249.2, b: 6356515.0, f: 1 / 293.466021294 },
    GRS80: { a: 6378137, b: 6356752.314140, f: 1 / 298.257222101 },
    Intl1924: { a: 6378388, b: 6356911.946, f: 1 / 297 }, // aka Hayford
    WGS72: { a: 6378135, b: 6356750.5, f: 1 / 298.26 },
};

/**
 * Datums; with associated ellipsoid, and Helmert transform parameters to convert from WGS 84 into
 * given datum.
 *
 * Note that precision of various datums will vary, and WGS-84 (original) is not defined to be
 * accurate to better than ±1 metre. No transformation should be assumed to be accurate to better
 * than a meter; for many datums somewhat less.
 */
LatLon.datum = {
    // transforms: t in metres, s in ppm, r in arcseconds                    tx       ty        tz       s        rx       ry       rz
    ED50: { ellipsoid: LatLon.ellipsoid.Intl1924, transform: [89.5, 93.8, 123.1, -1.2, 0.0, 0.0, 0.156] },
    Irl1975: { ellipsoid: LatLon.ellipsoid.AiryModified, transform: [-482.530, 130.596, -564.557, -8.150, -1.042, -0.214, -0.631] },
    NAD27: { ellipsoid: LatLon.ellipsoid.Clarke1866, transform: [8, -160, -176, 0, 0, 0, 0] },
    NAD83: { ellipsoid: LatLon.ellipsoid.GRS80, transform: [1.004, -1.910, -0.515, -0.0015, 0.0267, 0.00034, 0.011] },
    NTF: { ellipsoid: LatLon.ellipsoid.Clarke1880IGN, transform: [168, 60, -320, 0, 0, 0, 0] },
    OSGB36: { ellipsoid: LatLon.ellipsoid.Airy1830, transform: [-446.448, 125.157, -542.060, 20.4894, -0.1502, -0.2470, -0.8421] },
    Potsdam: { ellipsoid: LatLon.ellipsoid.Bessel1841, transform: [-582, -105, -414, -8.3, 1.04, 0.35, -3.08] },
    TokyoJapan: { ellipsoid: LatLon.ellipsoid.Bessel1841, transform: [148, -507, -685, 0, 0, 0, 0] },
    WGS72: { ellipsoid: LatLon.ellipsoid.WGS72, transform: [0, 0, -4.5, -0.22, 0, 0, 0.554] },
    WGS84: { ellipsoid: LatLon.ellipsoid.WGS84, transform: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] },
};
/* sources:
 * - ED50:          www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4
 * - Irl1975:       www.osi.ie/wp-content/uploads/2015/05/transformations_booklet.pdf
 *   ... note: many sources have opposite sign to rotations - to be checked!
 * - NAD27:         en.wikipedia.org/wiki/Helmert_transformation
 * - NAD83: (2009); www.uvm.edu/giv/resources/WGS84_NAD83.pdf
 *   ... note: functionally ≡ WGS84 - if you *really* need to convert WGS84<->NAD83, you need more knowledge than this!
 * - NTF:           Nouvelle Triangulation Francaise geodesie.ign.fr/contenu/fichiers/Changement_systeme_geodesique.pdf
 * - OSGB36:        www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
 * - Potsdam:       kartoweb.itc.nl/geometrics/Coordinate%20transformations/coordtrans.html
 * - TokyoJapan:    www.geocachingtoolbox.com?page=datumEllipsoidDetails
 * - WGS72:         www.icao.int/safety/pbn/documentation/eurocontrol/eurocontrol wgs 84 implementation manual.pdf
 *
 * more transform parameters are available from earth-info.nga.mil/GandG/coordsys/datums/NATO_DT.pdf,
 * www.fieldenmaps.info/cconv/web/cconv_params.js
 */


/**
 * Converts ‘this’ lat/lon coordinate to new coordinate system.
 *
 * @param   {LatLon.datum} toDatum - Datum this coordinate is to be converted to.
 * @returns {LatLon} This point converted to new datum.
 *
 * @example
 *     var pWGS84 = new LatLon(51.4778, -0.0016, LatLon.datum.WGS84);
 *     var pOSGB = pWGS84.convertDatum(LatLon.datum.OSGB36); // 51.4773°N, 000.0000°E
 */
LatLon.prototype.convertDatum = function (toDatum) {
    var oldLatLon = this;
    var transform = null;

    if (oldLatLon.datum == LatLon.datum.WGS84) {
        // converting from WGS 84
        transform = toDatum.transform;
    }
    if (toDatum == LatLon.datum.WGS84) {
        // converting to WGS 84; use inverse transform (don't overwrite original!)
        transform = [];
        for (var p = 0; p < 7; p++) transform[p] = -oldLatLon.datum.transform[p];
    }
    if (transform == null) {
        // neither this.datum nor toDatum are WGS84: convert this to WGS84 first
        oldLatLon = this.convertDatum(LatLon.datum.WGS84);
        transform = toDatum.transform;
    }

    var oldCartesian = oldLatLon.toCartesian();                // convert polar to cartesian...
    console.log(`oldCartesian is ${oldCartesian.x}, ${oldCartesian.y}, ${oldCartesian.z}`);
    var newCartesian = oldCartesian.applyTransform(transform); // ...apply transform...
    console.log(`newCartesian is ${newCartesian.x}, ${newCartesian.y}, ${newCartesian.z}`);
    var newLatLon = newCartesian.toLatLonE(toDatum);           // ...and convert cartesian to polar

    console.log(`newlatlon is ${newLatLon.lat / 180 * Math.PI}, ${newLatLon.lon / 180 * Math.PI}`);

    return newLatLon;
};


/**
 * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric) cartesian
 * (x/y/z) coordinates.
 *
 * @returns {Vector3d} Vector pointing to lat/lon point, with x, y, z in metres from earth centre.
 */
LatLon.prototype.toCartesian = function () {
    var φ = this.lat.toRadians(), λ = this.lon.toRadians();
    console.log(`φ is ${φ}`);
    console.log(`λ is ${λ}`);
    var h = 0; // height above ellipsoid - not currently used
    var a = this.datum.ellipsoid.a, f = this.datum.ellipsoid.f;

    var sinφ = Math.sin(φ), cosφ = Math.cos(φ);
    var sinλ = Math.sin(λ), cosλ = Math.cos(λ);

    var eSq = 2 * f - f * f;                      // 1st eccentricity squared ≡ (a²-b²)/a²
    var ν = a / Math.sqrt(1 - eSq * sinφ * sinφ); // radius of curvature in prime vertical

    var x = (ν + h) * cosφ * cosλ;
    var y = (ν + h) * cosφ * sinλ;
    var z = (ν * (1 - eSq) + h) * sinφ;

    var point = new Vector3d(x, y, z);

    return point;
};


/**
 * Converts ‘this’ (geocentric) cartesian (x/y/z) point to (ellipsoidal geodetic) latitude/longitude
 * coordinates on specified datum.
 *
 * Uses Bowring’s (1985) formulation for μm precision in concise form.
 *
 * @param {LatLon.datum.transform} datum - Datum to use when converting point.
 */
Vector3d.prototype.toLatLonE = function (datum) {
    var x = this.x, y = this.y, z = this.z;
    var a = datum.ellipsoid.a, b = datum.ellipsoid.b, f = datum.ellipsoid.f;

    var e2 = 2 * f - f * f;   // 1st eccentricity squared ≡ (a²-b²)/a²
    var ε2 = e2 / (1 - e2); // 2nd eccentricity squared ≡ (a²-b²)/b²
    var p = Math.sqrt(x * x + y * y); // distance from minor axis
    var R = Math.sqrt(p * p + z * z); // polar radius

    // parametric latitude (Bowring eqn 17, replacing tanβ = z·a / p·b)
    var tanβ = (b * z) / (a * p) * (1 + ε2 * b / R);
    var sinβ = tanβ / Math.sqrt(1 + tanβ * tanβ);
    var cosβ = sinβ / tanβ;

    // geodetic latitude (Bowring eqn 18: tanφ = z+ε²bsin³β / p−e²cos³β)
    var φ = isNaN(cosβ) ? 0 : Math.atan2(z + ε2 * b * sinβ * sinβ * sinβ, p - e2 * a * cosβ * cosβ * cosβ);

    // longitude
    var λ = Math.atan2(y, x);

    // height above ellipsoid (Bowring eqn 7) [not currently used]
    var sinφ = Math.sin(φ), cosφ = Math.cos(φ);
    var ν = a / Math.sqrt(1 - e2 * sinφ * sinφ); // length of the normal terminated by the minor axis
    var h = p * cosφ + z * sinφ - (a * a / ν);

    var point = new LatLon(φ.toDegrees(), λ.toDegrees(), datum);

    return point;
};

/**
 * Applies Helmert transform to ‘this’ point using transform parameters t.
 *
 * @private
 * @param   {number[]} t - Transform to apply to this point.
 * @returns {Vector3} Transformed point.
 */
Vector3d.prototype.applyTransform = function (t) {

    // this point
    var x1 = this.x, y1 = this.y, z1 = this.z;

    // transform parameters
    var tx = t[0];                    // x-shift
    var ty = t[1];                    // y-shift
    var tz = t[2];                    // z-shift
    var s1 = t[3] / 1e6 + 1;            // scale: normalise parts-per-million to (s+1)
    var rx = (t[4] / 3600).toRadians(); // x-rotation: normalise arcseconds to radians
    var ry = (t[5] / 3600).toRadians(); // y-rotation: normalise arcseconds to radians
    var rz = (t[6] / 3600).toRadians(); // z-rotation: normalise arcseconds to radians

    // apply transform
    var x2 = tx + x1 * s1 - y1 * rz + z1 * ry;
    var y2 = ty + x1 * rz + y1 * s1 - z1 * rx;
    var z2 = tz - x1 * ry + y1 * rx + z1 * s1;

    return new Vector3d(x2, y2, z2);
};


/**
 * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
 * degrees+minutes+seconds.
 *
 * @param   {string} [format=dms] - Format point as 'd', 'dm', 'dms'.
 * @param   {number} [dp=0|2|4] - Number of decimal places to use - default 0 for dms, 2 for dm, 4 for d.
 * @returns {string} Comma-separated latitude/longitude.
 */
LatLon.prototype.toString = function (format, dp) {
    return Dms.toLat(this.lat, format, dp) + ', ' + Dms.toLon(this.lon, format, dp);
};


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

/** Extend Number object with method to convert numeric degrees to radians */
if (Number.prototype.toRadians === undefined) {
    Number.prototype.toRadians = function () { return this * Math.PI / 180; };
}

/** Extend Number object with method to convert radians to numeric (signed) degrees */
if (Number.prototype.toDegrees === undefined) {
    Number.prototype.toDegrees = function () { return this * 180 / Math.PI; };
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
if (typeof module != 'undefined' && module.exports) module.exports = LatLon, module.exports.Vector3d = Vector3d; // ≡ export { LatLon as default, Vector3d }

var poole = new LatLon(50.715812, -2.011931);
var osgrid = OsGridRef.latLonToOsGrid(poole);
console.log(`osgrid ref is ${osgrid.easting}, ${osgrid.northing}`)
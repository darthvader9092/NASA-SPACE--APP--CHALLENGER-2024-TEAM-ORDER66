const exoplanetAttributes = {
    "space": ['Orbital Mechanics', 'Eccentricity and Orbital Energy', 'Habitable Zone Calculation', 'Transit Method and Planet Radius Estimation', 'Gravitational Lensing'],
    "time": ['Tidal Locking Timescale', 'Radial Velocity and Doppler Shift', 'Radial Velocity and Planet Mass'],
    "energy": ['Escape Velocity of an Exoplanet', 'Stellar Luminosity and Insolation Flux', 'Energy Flux and Blackbody Radiation']
};

let gravity = {};
let data = []; // Variable to hold CSV data
let planet_NAME = "Proxima Centauri"

const solarSystemData = {
    sun: {},
    venus: {}
};

function updateAttributes() {
    const typeSelect = document.getElementById('exoplanetType');
    const attributesSelect = document.getElementById('attributes');
    const selectedType = typeSelect.value;

    attributesSelect.innerHTML = ''; // Clear existing options

    if (exoplanetAttributes[selectedType]) {
        exoplanetAttributes[selectedType].forEach(attribute => {
            const option = document.createElement('option');
            option.value = attribute;
            option.textContent = attribute;
            attributesSelect.appendChild(option);
        });
    }
}

async function populateDropdown() {
    try {
        const response = await fetch('exo.csv');
        const csvData = await response.text();
        data = csvData.split('\n').slice(1).map(row => row.split(',')); // Store parsed data

        const planetSelect = document.getElementById('planet');
        data.forEach(row => {
            const planetName = row[0].trim();
            const planetGravity = parseFloat(row[1].trim());

            const option = document.createElement('option');
            option.value = planetName;
            option.textContent = planetName;
            planetSelect.appendChild(option);

            gravity[planetName] = planetGravity; // Store gravity
        });
    } catch (error) {
        console.error('Error loading the CSV file:', error);
    }
}

async function calculate() {
    const response = await fetch('exo.csv');
        const csvData = await response.text();
        data = csvData.split('\n').slice(1).map(row => row.split(',')); // Store parsed data

    const planetSelect = document.getElementById('planet');
    const selectedPlanet = planetSelect.value;
    planet_NAME=selectedPlanet;
    const attributesSelect = document.getElementById('attributes');
    const selectedAttribute = attributesSelect.value;

    // Get data related to the selected planet
    const planetData = data.find(row => row[0].trim() === selectedPlanet);

    if (!planetData) {
        alert("No data available for the selected planet.");
        return;
    }

    // Extract columns related to the planet for display
        // Extract relevant planet data from the planetData array
const pl_name = planetData[0]; // Planet Name
planet_NAME=pl_name;
const hostname = planetData[1]; // Host Star Name
const default_flag = planetData[2]; // Default Flag (Indicates if the planet is a standard reference)
const sy_snum = planetData[3]; // Number of Stars in the System
const sy_pnum = planetData[4]; // Number of Planets in the System
const discoverymethod = planetData[5]; // Method of Discovery
const disc_year = planetData[6]; // Year of Discovery
const disc_facility = planetData[7]; // Discovery Facility
const soltype = planetData[8]; // Type of Solar System (e.g., single star, binary)
const pl_controv_flag = planetData[9]; // Controversial Flag (Indicates if the planet's existence is disputed)
const pl_refname = planetData[10]; // Reference Name for the Planet
const pl_orbper = planetData[11]; // Orbital Period (in days)
const pl_orbpererr1 = planetData[12]; // Error in Orbital Period (Lower bound)
const pl_orbpererr2 = planetData[13]; // Error in Orbital Period (Upper bound)
const pl_orbperlim = planetData[14]; // Limiting Orbital Period
const pl_orbsmax = planetData[15]; // Semi-Major Axis (in AU)
const pl_orbsmaxerr1 = planetData[16]; // Error in Semi-Major Axis (Lower bound)
const pl_orbsmaxerr2 = planetData[17]; // Error in Semi-Major Axis (Upper bound)
const pl_orbsmaxlim = planetData[18]; // Limiting Semi-Major Axis
const pl_rade = planetData[19]; // Equatorial Radius (in Earth radii)
const pl_radeerr1 = planetData[20]; // Error in Equatorial Radius (Lower bound)
const pl_radeerr2 = planetData[21]; // Error in Equatorial Radius (Upper bound)
const pl_radelim = planetData[22]; // Limiting Equatorial Radius
const pl_radj = planetData[23]; // Jupiter Radius (in Jupiter radii)
const pl_radjerr1 = planetData[24]; // Error in Jupiter Radius (Lower bound)
const pl_radjerr2 = planetData[25]; // Error in Jupiter Radius (Upper bound)
const pl_radjlim = planetData[26]; // Limiting Jupiter Radius
const pl_bmasse = planetData[27]; // Mass of the Planet (in Earth masses)
const pl_bmasseerr1 = planetData[28]; // Error in Mass (Lower bound)
const pl_bmasseerr2 = planetData[29]; // Error in Mass (Upper bound)
const pl_bmasselim = planetData[30]; // Limiting Mass
const pl_bmassj = planetData[31]; // Mass of the Planet (in Jupiter masses)
const pl_bmassjerr1 = planetData[32]; // Error in Jupiter Mass (Lower bound)
const pl_bmassjerr2 = planetData[33]; // Error in Jupiter Mass (Upper bound)
const pl_bmassjlim = planetData[34]; // Limiting Jupiter Mass
const pl_bmassprov = planetData[35]; // Mass Provision Flag
const pl_orbeccen = planetData[36]; // Orbital Eccentricity
const pl_orbeccenerr1 = planetData[37]; // Error in Eccentricity (Lower bound)
const pl_orbeccenerr2 = planetData[38]; // Error in Eccentricity (Upper bound)
const pl_orbeccenlim = planetData[39]; // Limiting Eccentricity
const pl_insol = planetData[40]; // Insolation (in Earth units)
const pl_insolerr1 = planetData[41]; // Error in Insolation (Lower bound)
const pl_insolerr2 = planetData[42]; // Error in Insolation (Upper bound)
const pl_insollim = planetData[43]; // Limiting Insolation
const pl_eqt = planetData[44]; // Equilibrium Temperature (in Kelvin)
const pl_eqterr1 = planetData[45]; // Error in Equilibrium Temperature (Lower bound)
const pl_eqterr2 = planetData[46]; // Error in Equilibrium Temperature (Upper bound)
const pl_eqtlim = planetData[47]; // Limiting Equilibrium Temperature
const ttv_flag = planetData[48]; // TTV (Transit Timing Variation) Flag
const st_refname = planetData[49]; // Reference Name for the Star
const st_spectype = planetData[50]; // Spectral Type of the Host Star
const st_teff = planetData[51]; // Effective Temperature of the Star (in Kelvin)
const st_tefferr1 = planetData[52]; // Error in Effective Temperature (Lower bound)
const st_tefferr2 = planetData[53]; // Error in Effective Temperature (Upper bound)
const st_tefflim = planetData[54]; // Limiting Effective Temperature
const st_rad = planetData[55]; // Radius of the Star (in Solar radii)
const st_raderr1 = planetData[56]; // Error in Star Radius (Lower bound)
const st_raderr2 = planetData[57]; // Error in Star Radius (Upper bound)
const st_radlim = planetData[58]; // Limiting Star Radius
const st_mass = planetData[59]; // Mass of the Star (in Solar masses)
const st_masserr1 = planetData[60]; // Error in Star Mass (Lower bound)
const st_masserr2 = planetData[61]; // Error in Star Mass (Upper bound)
const st_masslim = planetData[62]; // Limiting Star Mass
const st_met = planetData[63]; // Metallicity of the Star
const st_meterr1 = planetData[64]; // Error in Metallicity (Lower bound)
const st_meterr2 = planetData[65]; // Error in Metallicity (Upper bound)
const st_metlim = planetData[66]; // Limiting Metallicity
const st_metratio = planetData[67]; // Metallicity Ratio
const st_logg = planetData[68]; // Logarithm of Surface Gravity of the Star
const st_loggerr1 = planetData[69]; // Error in Log Surface Gravity (Lower bound)
const st_loggerr2 = planetData[70]; // Error in Log Surface Gravity (Upper bound)
const st_logglim = planetData[71]; // Limiting Log Surface Gravity
const sy_refname = planetData[72]; // Reference Name for the System
const rastr = planetData[73]; // Right Ascension (in degrees)
const ra = planetData[74]; // Right Ascension (in hours, minutes, seconds)
const decstr = planetData[75]; // Declination (in degrees)
const dec = planetData[76]; // Declination (in degrees, decimal)
const sy_dist = planetData[77]; // Distance to the System (in parsecs)
const sy_disterr1 = planetData[78]; // Error in Distance (Lower bound)
const sy_disterr2 = planetData[79]; // Error in Distance (Upper bound)
const sy_vmag = planetData[80]; // Visual Magnitude of the Host Star
const sy_vmagerr1 = planetData[81]; // Error in Visual Magnitude (Lower bound)
const sy_vmagerr2 = planetData[82]; // Error in Visual Magnitude (Upper bound)
const sy_kmag = planetData[83]; // K-band Magnitude of the Host Star
const sy_kmagerr1 = planetData[84]; // Error in K-band Magnitude (Lower bound)
const sy_kmagerr2 = planetData[85]; // Error in K-band Magnitude (Upper bound)
const sy_gaiamag = planetData[86]; // Gaiamag of the Host Star
const sy_gaiamagerr1 = planetData[87]; // Error in Gaiamag (Lower bound)
const sy_gaiamagerr2 = planetData[88]; // Error in Gaiamag (Upper bound)
const rowupdate = planetData[89]; // Last Updated Row
const pl_pubdate = planetData[90]; // Publication Date for the Planet
const releasedate = planetData[91]; // Release Date for the Data



function calculateVenusData(pl_orbper = 225, pl_rade = 0.949, pl_orbsmax = 0.72) {
    const AU_TO_KM = 149597870.7; // 1 AU in km
    const EARTH_RADIUS_KM = 6371; // Earth radius in km
    const secondsPerDay = 86400; // Seconds in a day

    let results = {
        speed: 34.76,
        size: 38025,
        distance: 107489691
    };

    // Convert orbital period to seconds
    const orbitalPeriodInSeconds = pl_orbper * secondsPerDay;

    // Convert semi-major axis to km
    const semiMajorAxisInKm = pl_orbsmax * AU_TO_KM;

    // Calculate speed (orbital velocity in km/h)
    results.speed = (2 * Math.PI * semiMajorAxisInKm) / orbitalPeriodInSeconds * 3600; // km/h

    // Calculate size (equatorial radius in km)
    const radiusInKm = pl_rade * EARTH_RADIUS_KM; // Radius in km
    results.size = 2 * Math.PI * radiusInKm; // Circumference in km

    // Calculate distance (semi-major axis in km)
    results.distance = semiMajorAxisInKm; // Distance in km

    return results;
}

function calculateSunData(st_rad = 1) {
    const EARTH_RADIUS_KM = 6371; // Earth radius in km

    let results = {
        speed: 0, // Default speed for the Sun
        size: 0,
        distance: 149598262 // Default distance in km
    };

    // Calculate size based on the radius of the star (in Solar radii)
    results.size = st_rad * 695700; // Sun's radius in km is approximately 695,700 km

    return results;
}

// Example usage for Venus
solarSystemData.venus = calculateVenusData(); 

// Example usage for Sun with a given radius
solarSystemData.sun = calculateSunData(); 







    
    let output1 = `<p>Planet Name: ${pl_name}</p>
    <p>Host Name: ${hostname}</p>
    <p>Discovery Method: ${discoverymethod}</p>
    <p>Discovery Year: ${disc_year}</p>
    <p>Number of Stars in the System: ${sy_snum}</p>
    <p>Number of Planets in the System: ${sy_pnum}</p>
    <p>Orbital Period (days): ${pl_orbper}</p>
    <p>Semi-Major Axis (AU): ${pl_orbsmax}</p>
    <p>Planet Mass (Earth Masses): ${pl_bmasse}</p>
    <p>Planet Radius (Earth Radii): ${pl_rade}</p>
    <p>Equilibrium Temperature (K): ${pl_eqt}</p>
    <p>Eccentricity: ${pl_orbeccen}</p>
    <p>Distance from Earth (parsec units): ${sy_dist}</p>
    <p>Host Star Type: ${soltype}</p>
    <p>Host Star Mass (Solar Masses): ${st_mass}</p>`;
    const G = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2
    const AU_TO_METERS = 1.496e11; // 1 AU in meters
    const SOLAR_MASS = 1.989e30; // Solar mass in kg
    // Handle different calculations based on selected attribute
    switch (selectedAttribute) {
        case "Orbital Mechanics":
            

            function keplersThirdLaw(a, T, M_star) {
                // a in AU, T in days; convert T to years
                T /= 365.25; // Convert days to years
                a *= AU_TO_METERS; // Convert a from AU to meters
                M_star *= SOLAR_MASS; // Convert mass from solar masses to kg

                // Calculate the theoretical orbital period (T) based on the mass of the star (M_star)
                const calculatedT = Math.sqrt((4 * Math.PI ** 2 * (a ** 3)) / (G * M_star)); // Result in seconds
                const calculatedT_years = calculatedT / (60 * 60 * 24 * 365.25); // Convert seconds to years
                return calculatedT_years; // Returns theoretical period in years
            }

            const keplersThirdLaw_value = keplersThirdLaw(parseFloat(pl_orbsmax), parseFloat(pl_orbper), parseFloat(sy_snum));

            output1 += `<p> Semi-major axis (a): ${pl_orbsmax} AU</p>
                        <p> Orbital period (T): ${pl_orbper} days</p>
                        <p> Kepler's Third Law calculation for ${selectedPlanet} gives an orbital period of approximately ${keplersThirdLaw_value.toFixed(2)} years.</p>`;
            break;


        case "Eccentricity and Orbital Energy":
            
            // Assuming pl_orbsmax is the semi-major axis in AU and pl_orbeccen is the eccentricity
            const semiMajorAxis = parseFloat(pl_orbsmax); // Semi-major axis in AU
            const eccentricity = parseFloat(pl_orbeccen); // Eccentricity
        
            // Calculate the orbital energy (specific orbital energy)
            // Specific orbital energy (ε) = -GM / (2a)
            const M_star = parseFloat(sy_snum) * SOLAR_MASS; // Mass of the star in kg
            const a_meters = semiMajorAxis * AU_TO_METERS; // Convert semi-major axis to meters
        
            // Calculate gravitational parameter (GM)
            const gravitationalParameter = G * M_star; // GM in m^3/s^2
        
            // Calculate specific orbital energy
            const specificOrbitalEnergy = -gravitationalParameter / (2 * a_meters); // in J/kg
        
            // Output formatting
            output1 += `<p>Eccentricity and orbital energy information for ${selectedPlanet}</p>`;
            output1 += `<p> Semi-Major Axis (a): ${semiMajorAxis} AU</p>`;
            output1 += `<p> Eccentricity (e): ${eccentricity}</p>`;
            output1 += `<p> Specific Orbital Energy: ${specificOrbitalEnergy.toFixed(2)} J/kg</p>`;
            break;
            

        case "Habitable Zone Calculation":
            // Example default values for radius and temperature
            const DEFAULT_RADIUS = 1; // Default star radius in solar radii
            const DEFAULT_TEMPERATURE = 5778; // Default effective temperature in Kelvin (Sun's temperature)
        
            // Get values from your dataset; use default if they are null
            const R_star = parseFloat(st_rad) || DEFAULT_RADIUS; // Star radius
            const T_star = parseFloat(st_teff) || DEFAULT_TEMPERATURE; // Effective temperature
        
            // Function to calculate luminosity
            function calculateLuminosity(R_star, T_star) {
                const sigma = 5.67e-8; // Stefan-Boltzmann constant in W/m^2/K^4
                const R_star_meters = R_star * 6.957e8; // Convert radius to meters
                return 4 * Math.PI * R_star_meters ** 2 * sigma * T_star ** 4; // Luminosity in Watts
            }
        
            // Function to calculate the habitable zone
            function habitableZone(R_star, T_star) {
                const L_star = calculateLuminosity(R_star, T_star); // Calculate luminosity
                let innerBoundary = (Math.sqrt(L_star / 1.1))/10000000000000; // Inner boundary in AU
                let outerBoundary = (Math.sqrt(L_star / 0.53))/10000000000000; // Outer boundary in AU
                return { innerBoundary , outerBoundary};
            }
        
            // Calculate habitable zone boundaries
            const habitableZoneBoundaries = habitableZone(R_star, T_star);
        
            // Output the results
            output1 += `<h4>Habitable Zone calculation for ${selectedPlanet}</h4>
                        <p><strong>Stellar Radius:</strong> ${R_star} Solar Radii</p>
                        <p><strong>Effective Temperature:</strong> ${T_star} K</p>
                        <p>Inner Boundary: ${habitableZoneBoundaries.innerBoundary.toFixed(2)} AU</p>
                        <p>Outer Boundary: ${habitableZoneBoundaries.outerBoundary.toFixed(2)} AU</p>`;
            break;
            
        case "Transit Method and Planet Radius Estimation":
            // Given parameters
            const dropInBrightness = 0.012; // Drop in brightness as a decimal (1.2%)
            
            // Get the stellar radius, use default if null
            const r_star = parseFloat(st_rad) || 1; // Default to 1 solar radius if null
        
            // Function to estimate planet radius
            function estimatePlanetRadius(dropInBrightness, r_star) {
                return r_star * Math.sqrt(dropInBrightness); // Planet radius in solar radii
            }
        
            // Calculate planet radius
            const planetRadius = estimatePlanetRadius(dropInBrightness, r_star);
        
            // Output the results
            output1 += `<h4>Transit Method and Planet Radius Estimation for ${selectedPlanet}</h4>
                        <p><strong>Drop in Brightness:</strong> ${dropInBrightness * 100}%</p>
                        <p><strong>Stellar Radius:</strong> ${r_star} Solar Radii</p>
                        <p><strong>Estimated Planet Radius:</strong> ${planetRadius.toFixed(2)} Solar Radii</p>`;
            break;
            
        case "Gravitational Lensing":
            // Given parameters
            const D_L = parseFloat(sy_dist); // Distance to the lens in parsecs
        
            // Get the mass of the star, using default if null
            const Mm_star = parseFloat(st_mass) || 1; // Default to 1 solar mass if null
        
            // Function to calculate the Einstein radius
            function calculateEinsteinRadius(D_L, Mm_star) {
                const G = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2
                const c = 3e8; // Speed of light in m/s
                const parsecToMeters = 3.086e16; // 1 parsec in meters
        
                // Convert distance to lens from parsecs to meters
                const D_L_meters = D_L * parsecToMeters;
        
                // Convert mass from solar masses to kg
                const M_star_kg = Mm_star * 1.989e30; // Mass of the Sun in kg
        
                // Calculate the Einstein radius in meters
                const theta_E = Math.sqrt((4 * G * M_star_kg) / (c ** 2 * D_L_meters)); // Result in meters
        
                // Convert Einstein radius to arcseconds
                const theta_E_arcseconds = theta_E * (180 / Math.PI) * (3600 / parsecToMeters); // Convert to arcseconds
        
                return theta_E_arcseconds; // Returns Einstein radius in arcseconds
            }
        
            // Calculate Einstein radius
            const einsteinRadius = calculateEinsteinRadius(D_L, Mm_star);
        
            // Output the results
            output1 += `<h4>Gravitational Lensing Calculation for ${selectedPlanet}</h4>
                        <p><strong>Distance to the Lens:</strong> ${D_L} parsecs</p>
                        <p><strong>Mass of Star:</strong> ${Mm_star} Solar Masses</p>
                        <p><strong>Einstein Radius:</strong> ${einsteinRadius.toFixed(2)} arcseconds</p>`;
            break;
            

        case "Tidal Locking Timescale":
            const EARTH_RADIUS_METERS = 6.371e6; // Earth's radius in meters
            const EARTH_MASS_KG = 5.972e24; // Earth's mass in kg
        
            // Default values
            const DEFAULT_SEMI_MAJOR_AXIS = 1; // Default semi-major axis in AU (e.g., 1 AU)
            const DEFAULT_PLANET_RADIUS = 1; // Default planet radius in Earth radii (e.g., Earth)
            const DEFAULT_PLANET_MASS = 1; // Default planet mass in Earth masses (e.g., Earth)
        
            // Fetch and convert parameters
            const semiMajorAxis1 = (parseFloat(pl_orbsmax) || DEFAULT_SEMI_MAJOR_AXIS) * AU_TO_METERS; // Semi-major axis in meters
            const planetRadius1 = (parseFloat(pl_rade) || DEFAULT_PLANET_RADIUS) * EARTH_RADIUS_METERS; // Planet radius in meters
            const planetMass = (parseFloat(pl_bmasse) || DEFAULT_PLANET_MASS) * EARTH_MASS_KG; // Planet mass in kg
        
            // Function to calculate tidal locking timescale
            function tidalLockingTimescale(a, R, M) {
                return (a ** 6) / (M * (R ** 5)); // Result in seconds
            }
        
            // Calculate tidal locking timescale
            const tidalLockingTime = tidalLockingTimescale(semiMajorAxis1, planetRadius1, planetMass);
        
            // Convert seconds to years
            const tidalLockingTimeYears = tidalLockingTime / (60 * 60 * 24 * 365.25);
        
            // Output the results
            output1 += `<p>Tidal Locking Timescale calculation for ${selectedPlanet}</p>
                        <p>Semi-Major Axis (a): ${(semiMajorAxis1 / AU_TO_METERS).toFixed(2)} AU (default if null)</p>
                        <p>Planet Radius (R): ${(planetRadius1 / EARTH_RADIUS_METERS).toFixed(2)} Earth Radii (default if null)</p>
                        <p>Planet Mass (M): ${(planetMass / EARTH_MASS_KG).toFixed(2)} Earth Masses (default if null)</p>
                        <p>Tidal Locking Timescale: approximately ${tidalLockingTimeYears.toFixed(2)} years</p>`;
            break;
    
        case "Radial Velocity and Planet Mass":
            // Default values
            const DEFAULT_RADIAL_VELOCITY = 10; // Default radial velocity shift in m/s
        
            // Fetch and convert parameters
            const massStar = parseFloat(st_mass) || 1; // Mass of host star in solar masses (default to 1 if null)
            const radialVelocityShift = DEFAULT_RADIAL_VELOCITY; // Radial velocity shift in m/s
            const orbitalPeriod = parseFloat(pl_orbper) || 1; // Orbital period in days (default to 1 if null)
        
            // Function to calculate planet mass using the radial velocity method
            function calculatePlanetMass(M_star, delta_v, T) {
                const G = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2
                const SOLAR_MASS = 1.989e30; // Solar mass in kg
        
                // Convert mass of host star from solar masses to kg
                const M_star_kg = M_star * SOLAR_MASS; 
        
                // Convert orbital period from days to seconds
                const T_seconds = T * 24 * 60 * 60; 
        
                // Calculate planet mass (M_p) using the formula
                const planetMass = (M_star_kg * delta_v ** 2 * T_seconds) / (2 * Math.PI * G); // Result in kg
                return planetMass / SOLAR_MASS; // Convert to solar masses
            }
        
            // Calculate planet mass
            const planetMassValue = (calculatePlanetMass(massStar, radialVelocityShift, orbitalPeriod))/10000000000000000;
        
            // Output the results
            output1 += `<p>Radial Velocity and Planet Mass calculation for ${selectedPlanet}</p>
                        <p>Mass of Host Star (M<sub>star</sub>): ${massStar.toFixed(2)} Solar Masses</p>
                        <p>Radial Velocity Shift (Δv): ${radialVelocityShift} m/s</p>
                        <p>Orbital Period (T): ${orbitalPeriod} days</p>
                        <p>Calculated Planet Mass: approximately ${planetMassValue.toFixed(2)} Solar Masses</p>`;
            break;

        case "Radial Velocity and Doppler Shift":
            // Default values
            const DEFAULT_RADIAL_VELOCITY1 = 12; // Default radial velocity in m/s

            // Fetch and convert parameters
            const massStarDoppler = parseFloat(st_mass) || 1; // Mass of host star in solar masses (default to 1 if null)
            const radialVelocityStar = DEFAULT_RADIAL_VELOCITY1; // Radial velocity in m/s
            const orbitalPeriodDoppler = parseFloat(pl_orbper) || 1; // Orbital period in days (default to 1 if null)

            // Function to calculate Doppler shift using the radial velocity method
            function calculateDopplerShift(M_star, v_star, T) {
                const G = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2
                const SOLAR_MASS = 1.989e30; // Solar mass in kg

                // Convert mass of host star from solar masses to kg
                const M_star_kg = M_star * SOLAR_MASS;

                // Convert orbital period from days to seconds
                const T_seconds = T * 24 * 60 * 60;

                // Calculate the frequency shift using the formula
                const frequencyShift = ((v_star / (3e8)) * (M_star_kg / (M_star_kg + (4 * Math.PI ** 2 * (T_seconds ** 2) * (G / (M_star_kg ** 2))))))*1000000;

                return frequencyShift; // Frequency shift in Hz
            }

            // Calculate Doppler shift
            const dopplerShiftValue = calculateDopplerShift(massStarDoppler, radialVelocityStar, orbitalPeriodDoppler);

            // Output the results
            output1 += `<p>Radial Velocity and Doppler Shift calculation for ${selectedPlanet}</p>
                        <p>Mass of Host Star (M<sub>star</sub>): ${massStarDoppler.toFixed(2)} Solar Masses</p>
                        <p>Radial Velocity (v<sub>star</sub>): ${radialVelocityStar} m/s</p>
                        <p>Orbital Period (T): ${orbitalPeriodDoppler} days</p>
                        <p>Calculated Doppler Shift: approximately ${dopplerShiftValue.toFixed(2)} Hz</p>`;
            break;

                
        case "Escape Velocity of an Exoplanet":
            // Constants
            const GRAVITATIONAL_CONSTANT1 = 6.67430e-11; // Gravitational constant in m^3/(kg*s^2)
            const EARTH_MASS_KG11 = 5.972e24; // Earth's mass in kg
            const EARTH_RADIUS_METERS11 = 6.371e6; // Earth's radius in meters
            const JUPITER_MASS_KG = 1.898e27; // Jupiter's mass in kg
            const JUPITER_RADIUS_METERS = 7.1492e7; // Jupiter's radius in meters

            // Initialize mass and radius variables without redeclaration
            let planetMass11;
            let planetRadius11;

            // Determine the planet mass
            if (pl_bmasse) {
                planetMass11 = parseFloat(pl_bmasse) * EARTH_MASS_KG11; // From Earth Mass
            } else if (pl_bmassj) {
                planetMass11 = parseFloat(pl_bmassj) * JUPITER_MASS_KG; // From Jupiter Mass
            } else {
                // Handle the case where no mass is provided
                planetMass11 = EARTH_MASS_KG; // Default to Earth's mass if neither is available
            }

            // Determine the planet radius
            if (pl_rade) {
                planetRadius11 = parseFloat(pl_rade) * EARTH_RADIUS_METERS11; // From Earth Radius
            } else if (pl_radj) {
                planetRadius11 = parseFloat(pl_radj) * JUPITER_RADIUS_METERS; // From Jupiter Radius
            } else {
                // Handle the case where no radius is provided
                planetRadius11 = EARTH_RADIUS_METERS11; // Default to Earth's radius if neither is available
            }

            // Function to calculate escape velocity
            function calculateEscapeVelocity(mass, radius) {
                return Math.sqrt((2 * GRAVITATIONAL_CONSTANT1 * mass) / radius); // Escape velocity in m/s
            }

            // Calculate escape velocity
            const escapeVelocity = calculateEscapeVelocity(planetMass11, planetRadius11);

            // Output the results
            output1 += `<p>Escape Velocity calculation for ${selectedPlanet}</p>
                        <p>Planet Mass (M<sub>planet</sub>): ${pl_bmasse ? pl_bmasse + ' Earth Masses' : pl_bmassj + ' Jupiter Masses'}</p>
                        <p>Planet Radius (R<sub>planet</sub>): ${pl_rade ? pl_rade + ' Earth Radii' : pl_radj + ' Jupiter Radii'}</p>
                        <p>Calculated Escape Velocity: approximately ${escapeVelocity.toFixed(2)} m/s</p>`;
            break;

        case "Energy Flux and Blackbody Radiation":
            // Constants
            const STEFAN_BOLTZMANN_CONSTANT = 5.67e-8; // Stefan-Boltzmann constant in W/m^2/K^4
            const SUN_RADIUS_METERS = 6.957e8; // Radius of the Sun in meters

            // Fetch and convert parameters
            const effectiveTemperature = parseFloat(st_teff) || 5778; // Effective temperature in Kelvin (default to Sun's temperature)
            const stellarRadius = parseFloat(st_rad) || 1; // Stellar radius in solar radii (default to 1)
            const distanceFromStar = parseFloat(pl_orbsmax) * AU_TO_METERS; // Semi-major axis in meters

            // Function to calculate energy flux at a distance from the star
            function energyFlux(T_star, R_star, d) {
                // Calculate luminosity using the Stefan-Boltzmann Law
                const luminosity = 4 * Math.PI * (R_star ** 2) * STEFAN_BOLTZMANN_CONSTANT * (T_star ** 4); // Luminosity in Watts

                // Calculate energy flux (F) at distance d
                const flux = luminosity / (4 * Math.PI * d ** 2); // Flux in W/m^2
                return flux;
            }

            // Convert stellar radius to meters
            const stellarRadiusMeters = stellarRadius * SUN_RADIUS_METERS; // Radius in meters

            // Calculate energy flux
            const energyFluxValue = energyFlux(effectiveTemperature, stellarRadiusMeters, distanceFromStar);

            // Output the results
            output1 += `<p>Energy Flux and Blackbody Radiation calculation for ${selectedPlanet}</p>
                        <p>Effective Temperature (T<sub>star</sub>): ${effectiveTemperature} K</p>
                        <p>Stellar Radius (R<sub>star</sub>): ${stellarRadius} Solar Radii</p>
                        <p>Distance from Star (d): ${pl_orbsmax} AU</p>
                        <p>Calculated Energy Flux: approximately ${energyFluxValue.toFixed(2)} W/m<sup>2</sup></p>`;
            break;

        

            case "Stellar Luminosity and Insolation Flux":
                // Constants
                const STEFAN_BOLTZMANN_CONSTANT1 = 5.67e-8; // Stefan-Boltzmann constant in W/m^2/K^4
                const AU_TO_METERS1 = 1.496e11; // 1 AU in meters
                const SUN_RADIUS_METERS1 = 6.957e8; // Radius of the Sun in meters
            
                // Fetch and convert parameters
                const effectiveTemp = parseFloat(st_teff) || 5778; // Effective temperature in Kelvin (default to Sun's temperature)
                const starRadius = parseFloat(st_rad) || 1; // Stellar radius in solar radii (default to 1)
                const semiMajorAxisAU = parseFloat(pl_orbsmax) * AU_TO_METERS1; // Semi-major axis in meters
            
                // Function to calculate stellar luminosity
                function calculateStellarLuminosity(temp, radius) {
                    return 4 * Math.PI * (radius ** 2) * STEFAN_BOLTZMANN_CONSTANT1 * (temp ** 4); // Luminosity in Watts
                }
            
                // Function to calculate insolation flux at a distance from the star
                function calculateInsolationFlux(luminosity, distance) {
                    return luminosity / (4 * Math.PI * distance ** 2); // Flux in W/m^2
                }
            
                // Calculate stellar radius in meters
                const starRadiusMeters = starRadius * SUN_RADIUS_METERS1; // Radius in meters
            
                // Calculate luminosity
                const stellarLuminosity = calculateStellarLuminosity(effectiveTemp, starRadiusMeters);
            
                // Calculate insolation flux
                const insolationFluxValue = calculateInsolationFlux(stellarLuminosity, semiMajorAxisAU);
            
                // Output the results
                output1 += `<p>Stellar Luminosity and Insolation Flux calculation for ${selectedPlanet}</p>
                            <p>Effective Temperature (T<sub>star</sub>): ${effectiveTemp} K</p>
                            <p>Stellar Radius (R<sub>star</sub>): ${starRadius} Solar Radii</p>
                            <p>Semi-Major Axis (a): ${pl_orbsmax} AU</p>
                            <p>Calculated Stellar Luminosity: approximately ${stellarLuminosity.toExponential(2)} Watts</p>
                            <p>Calculated Insolation Flux: approximately ${insolationFluxValue.toFixed(2)} W/m<sup>2</sup></p>`;
                break;
            

        default:
            output1 += "<p>Unknown attribute selected.</p>";
            break;
    }

    document.getElementById('textOutput').innerHTML = output1;
    document.getElementById('textOutput').style.display = 'block'; 
}

// Call the function to populate the planet dropdown on page load
window.onload = populateDropdown;

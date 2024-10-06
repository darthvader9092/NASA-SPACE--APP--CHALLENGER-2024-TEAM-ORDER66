// Function to handle the pendulum swing based on user input
function updatePendulum() {
    const force = parseFloat(document.getElementById("force").value);
    const weight = parseFloat(document.getElementById("weight").value);
    const planet = document.getElementById("planet").value;

    if (isNaN(force) || isNaN(weight) || weight <= 0) {
        // Invalid input handling (if needed)
        document.getElementById("decelerationInfo").innerText = `Invalid input`;
        return;
    }

    // Gravity values for different planets
    const g = {
        "Earth": 9.81,
        "Gliese 667 Cc": 9.27,
        "Kepler-186f": 9.67,
        "Kepler-22b": 9.47,
        "Proxima Centauri b": 9.81,
        "Kepler-1649c": 9.81,
        "Kepler-62f": 9.81,
        "TRAPPIST-1e": 9.81,
        "Wolf 1061c": 9.81,
        "Tau Ceti e": 9.81
    }[planet];

    // Calculate acceleration and deceleration
    const acceleration = force / weight; // a = F / m
    const deceleration = g - acceleration; // Net downward acceleration

    // Calculate time to rest
    const length = 1; // Length of the pendulum in meters
    const timeToRest = (Math.sqrt(2 * length) / Math.sqrt(deceleration)).toFixed(2); // Time to stop

    // Update deceleration and time to rest in real time
    document.getElementById("decelerationInfo").innerText =
        `Deceleration: ${deceleration.toFixed(2)} m/s²\n` +
        `Time to Rest: ${timeToRest} seconds`;

    const rod = document.getElementById("rod");

    // Set the rod's swing angle based on the applied force
    const maxSwingAngle = Math.min(30, (force / weight) * 45); // Maximum angle based on force
    const duration = 1000; // Duration for one full swing in milliseconds

    // Start the pendulum swinging
    let direction = 1;
    let angle = 0;

    function swingPendulum() {
        // Swing the pendulum left and right
        angle = maxSwingAngle * direction;
        rod.style.transform = `rotate(${angle}deg)`;

        // Change the direction after each swing
        direction *= -1;

        // Use requestAnimationFrame to schedule the next swing
        setTimeout(() => {
            requestAnimationFrame(swingPendulum);
        }, duration);
    }

    swingPendulum();

    // Stop the swinging after the time to rest
    setTimeout(() => {
        rod.style.transition = 'transform 1s ease-out';
        rod.style.transform = 'rotate(0deg)'; // Reset to original position
    }, timeToRest * 1000);
}

// Add event listeners to inputs for real-time updates
document.getElementById("force").addEventListener("input", updatePendulum);
document.getElementById("weight").addEventListener("input", updatePendulum);
document.getElementById("planet").addEventListener("change", updatePendulum);

// Initial pendulum setup
updatePendulum();

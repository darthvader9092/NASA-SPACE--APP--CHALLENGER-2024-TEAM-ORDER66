const gravity = {
    "Gliese 667 Cc": 1.13,
    "Kepler-186f": 0.61,
    "Kepler-22b": 2.4,
    "Proxima Centauri b": 1.17,
    "Kepler-1649c": 1.06,
    "Kepler-62f": 1.24,
    "TRAPPIST-1e": 0.93,
    "Wolf 1061c": 1.7,
    "Tau Ceti e": 1.55
    };
    
    let chart = null; // Variable to store chart instance
    
    function calculate() {
    // Get user input
    const height = document.getElementById("height").value;
    const weight = document.getElementById("weight").value;
    const planet = document.getElementById("planet").value;
    
    if (height === "" || weight === "") {
        alert("Please enter your height and weight.");
        return;
    }
    
    // Calculate new weight on selected exoplanet
    const newWeight = (weight * gravity[planet]).toFixed(2);
    
    // Update the text output and make it visible
    document.getElementById("textOutput").innerHTML = `
        Audios Amigos: Your meat on  <strong>${planet}</strong>
        will approximately weigh :  <strong>${newWeight} kg</strong>.
        Your height remains <strong>${height} cm</strong> (gravity does not affect height).
    `;
    
    // Show the text output container
    document.getElementById("textOutput").style.display = "block";
    
    // If a chart already exists, destroy it before creating a new one
    if (chart) {
        chart.destroy();
    }
    
    // Plot the result
    const ctx = document.getElementById("resultChart").getContext('2d');
    document.getElementById("resultChart").style.display = "block";
    
    chart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: ["Earth", planet],
            datasets: [{
                label: 'Weight (kg)',
                data: [weight, newWeight],
                backgroundColor: ['rgba(54, 162, 235, 0.6)', 'rgba(255, 99, 132, 0.6)'],
                borderColor: ['rgba(54, 162, 235, 1)', 'rgba(255, 99, 132, 1)'],
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            scales: {
                y: {
                    beginAtZero: true
                }
            }
        }
    });
    }
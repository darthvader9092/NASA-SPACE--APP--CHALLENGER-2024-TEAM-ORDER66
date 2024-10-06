// Create stars in random positions on the screen
const stars = 400;

for (let i = 0; i < stars; i++) {
    let star = document.createElement("div");
    star.className = 'stars';
    var xy = randomPosition();
    star.style.top = xy[0] + 'px';
    star.style.left = xy[1] + 'px';
    document.body.append(star);
}

// Function to generate random positions for the stars
function randomPosition() {
    var y = window.innerWidth; // Get the width of the window
    var x = window.innerHeight; // Get the height of the window
    var randomX = Math.floor(Math.random() * x); // Generate random x position
    var randomY = Math.floor(Math.random() * y); // Generate random y position
    return [randomX, randomY]; // Return the random position
}

// Wait for the window to load
window.onload = function() {
    var audio = document.getElementById('myAudio'); // Get the audio element
    var playButton = document.getElementById('playButton'); // Get the play button
    
    // Add click event to play button
    playButton.onclick = function() {
        if (audio.paused) {
            audio.play().catch(function(error) {
                console.log('Playback failed:', error); // Log playback error
                // You can add fallback behavior here if desired
            });
        } else {
            audio.pause(); // Pause the audio if it's currently playing
            audio.currentTime = 0; // Reset audio to the start
        }
    };
};


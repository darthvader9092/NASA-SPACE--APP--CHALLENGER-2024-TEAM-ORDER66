var scale = 1,
    panning = false,
    pointX = 0,
    pointY = 0,
    start = { x: 0, y: 0 };

const zoom = document.querySelector('.container');
const mercury = document.querySelector('.mercury');

function setTransform() {
    zoom.style.transform = "translate(" + pointX + "px, " + pointY + "px) scale(" + scale + ")";
}

zoom.onmousedown = function (e) {
    e.preventDefault();

    start = { x: e.clientX - pointX, y: e.clientY - pointY };
    panning = true;

    zoom.style.transition = "0s";
}

zoom.onmouseup = function (e) {
    panning = false;
}

zoom.onmousemove = function (e) {
    e.preventDefault();

    if (!panning) {
        return;
    }

    pointX = (e.clientX - start.x);
    pointY = (e.clientY - start.y);

    setTransform();
}

zoom.onwheel = function (e) {
    e.preventDefault();

    var ys = (e.clientY - pointY) / scale,
        delta = (e.wheelDelta ? e.wheelDelta : -e.deltaY);

    (delta > 0) ? (scale *= 1.05) : (scale /= 1.05);
    pointY = e.clientY - ys * scale;

    setTransform();

    zoom.style.transition = "0.3s";
}
body {
    margin: 0;
    overflow: hidden; /* Prevent scrollbars */
    background: #222;
    font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    user-select: none;
    -webkit-user-select: none;
}
body.mars {/*You can ignore the mars classes*/
    background: #821e1e;
}
body.earth,
body.slime {
    background: #3697f2;
}
body.moon { /*Moon is the end*/
    background: #1e0223;
}
body.ocean {
    background: #87ceeb;
}
.ground {
  position: fixed;   
  bottom: 0;         
  left: 0;           
  width: 100%;       
  height: 40vh;      
  z-index: 1000;     
  background-repeat: repeat-x, repeat;
  background-position: top left, top left;
  background-size:100px 100px, 100px 100px;
  image-rendering: pixelated; 
}

.ground.earth{
  background-image: url('grass.png'), url('dirt.png');
}
.ground.mars{
    background-image: url('ngrass.png'), url('netherrack.png');
}
.ground.moon{/*Moon is the end*/
    background-image: url('endstone.png'), url("endstone.png");
}
.ground.ocean{
    background-image: url('water.gif'), url("water.gif");
    opacity: 0.9;
}
.ground.slime{
    background-image: url('slime.png'), url("slime.png");
}

/* World & layers */

#world {
    position: relative;
    width: 100vw;
    height: 100vh;
    overflow: hidden;
}

#projectile-layer {
    position: absolute;
    inset: 0;
    pointer-events: none;
    z-index: 1500; /* Above ground (1000), below UI (2000) */
}

.projectile {
    position: absolute;
    transform: translate(-50%, 50%);
    image-rendering: pixelated;
    pointer-events: none;
    z-index: 1600;
}

.placed-block {
    position: absolute;
    transform: translate(0%, 0%);
    image-rendering: pixelated;
    pointer-events: none;
    z-index: 1100; /* Above ground, below projectiles */
}

.grid-highlight {
    position: absolute;
    border: 0.3vh dashed rgba(255, 255, 255, 0.7);
    background: rgba(255, 255, 255, 0.1);
    pointer-events: none;
    z-index: 1050;
    display: none;
}

.particle {
    position: absolute;
    background: #e0f7fa;
    border-radius: 50%;
    pointer-events: none;
    z-index: 1600;
    transform: translate(-50%, 50%);
}

/* Aim arrow */

.aim-line {
    position: absolute;
    transform-origin: left center;
    height: 0.5vh;
    background: rgba(255, 255, 255, 0.8);
    border-radius: 0.25vh;
    pointer-events: none;
    z-index: 1700; /* Above everything except UI */
}

.aim-line::after {
    content: "";
    position: absolute;
    right: 0;
    top: 50%;
    transform: translateY(-50%);
    border-left: 1vh solid rgba(255, 255, 255, 0.9);
    border-top: 0.6vh solid transparent;
    border-bottom: 0.6vh solid transparent;
}
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
    image-rendering: pixelated;
}

.ui-button img {
    width: 4vh;
    height: 4vh;
    object-fit: cover;
}

.ui-button.selected {
    outline: 0.3vh solid #ffe36e;
    box-shadow: 0 0 1vh rgba(255, 227, 110, 0.8);
}

#pause-button {
    color: white;
    font-size: 2vh;
    width: 4vh;
    height: 4vh;
}

#hotbar {
    position: absolute;
    left: 50%;
    transform: translateX(-50%);
    bottom: 1vh;
    display: flex;
    gap: 0.6vh;
    padding: 0.6vh 0.8vh;
    background: rgba(0, 0, 0, 0.45);
    border-radius: 1vh;
    backdrop-filter: blur(4px);
    pointer-events: auto;
}

.item-button img {
    width: 5vh;
    height: 5vh;
}

/* Aim arrow */

.aim-line {
    position: absolute;
    transform-origin: left center;
    height: 0.5vh;
    background: rgba(255, 255, 255, 0.8);
    border-radius: 0.25vh;
    pointer-events: none;
}

.aim-line::after {
    content: "";
    position: absolute;
    right: 0;
    top: 50%;
    transform: translateY(-50%);
    border-left: 1vh solid rgba(255, 255, 255, 0.9);
    border-top: 0.6vh solid transparent;
    border-bottom: 0.6vh solid transparent;
}


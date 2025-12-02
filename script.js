

const SCALE_M_PER_VH = 1 / 8; 
const VH_PER_M = 8;
const TARGET_FPS = 60;
const DT = 1 / TARGET_FPS; 
const BLOCK_SIZE = 1.0; 


const BiomeData = {
  earth: {
    name: "Earth",
    displayName: "Earth",
    gravity: 9.81,
    airDensity: 2.5,
    groundBounce: 1.0,
    isWater: false,
    waterLevelMeters: 0
  },
  moon: {
    name: "Moon",
    displayName: "Moon",
    gravity: 1.62,
    airDensity: 0.0,
    groundBounce: 1.0,
    isWater: false,
    waterLevelMeters: 0
  },
  ocean: {
    name: "Water",
    displayName: "Water",
    gravity: 9.81,
    airDensity: 2.5,
    groundBounce: 0.3,
    isWater: true,
    waterLevelMeters: 0
  },
  slime: {
    name: "Trampoline",
    displayName: "Trampoline",
    gravity: 9.81,
    airDensity: 2.5,
    groundBounce: 1.5,
    isWater: false,
    waterLevelMeters: 0
  }
};

const ItemData = {
  iron: {
    name: "Iron Cannonball",
    displayName: "Iron Cannonball",
    mass: 2354, 
    radius: 0.75,
    dragCoefficient: 0.5,
    restitution: 0.1,
    texture: "cannonball.png"
  },
  beach: {
    name: "Beach Ball",
    displayName: "Beach Ball",
    mass: 5.0,
    radius: 0.4,
    dragCoefficient: 0.5,
    restitution: 0.9,
    texture: "beachball.png"
  },
  snow: {
    name: "Normal Ball",
    displayName: "Normal Ball",
    mass: 0.5,
    radius: 0.2,
    dragCoefficient: 0.5,
    restitution: 0.8,
    texture: "snowball.png"
  },
  slime: {
    name: "Bouncy Ball",
    displayName: "Bouncy Ball",
    mass: 1.0,
    radius: 0.25,
    dragCoefficient: 0.6,
    restitution: 1.2,
    texture: "slimeball.png"
  }
};


const BlockData = {
  cobblestone: {
    name: "Cobblestone",
    texture: "cobblestone.png",
    restitution: 0.3
  },
  slime: {
    name: "Slime Block",
    texture: "slime.png",
    restitution: 1.5
  }
};

function physicsToScreen(x, y, worldHeightVh = 100, groundHeightVh = 40) {
  const pxVh = x * VH_PER_M;
  const pyVh = y * VH_PER_M;
  const bottom = groundHeightVh + pyVh;
  const left = pxVh;
  return { left, bottom };
}

let currentBiomeKey = "earth";
let currentItemKey = "iron";
let paused = false;
let timeScale = 1.0;
let launchPower = 1.2;
let integrationMethod = "semi-implicit";

let projectiles = [];
let particles = [];
let blocks = []; 
let aimState = {
  isAiming: false,
  startX: 0,
  startY: 0,
  currentX: 0,
  currentY: 0,
  lineEl: null
};
let dragState = {
  isDragging: false,
  projectile: null,
  offsetX: 0,
  offsetY: 0
};
let ctrlPressed = false;
let shiftPressed = false;
let gridHighlightEl = null;

const worldEl = document.getElementById("world");
const projectileLayerEl = document.getElementById("projectile-layer");
const groundEl = document.getElementById("ground");
const pauseButtonEl = document.getElementById("pause-button");
const velocityDisplayEl = document.getElementById("velocity-display");


const speedSlider = document.getElementById("speed-slider");
const speedVal = document.getElementById("speed-val");
const powerSlider = document.getElementById("power-slider");
const powerVal = document.getElementById("power-val");
const integratorSelect = document.getElementById("integrator-select");


function getMousePositionMeters(evt) {
  const rect = worldEl.getBoundingClientRect();
  const xVh = ((evt.clientX - rect.left) / rect.height) * 100;
  const yVhFromBottom = ((rect.bottom - evt.clientY) / rect.height) * 100;
  const xMeters = xVh / VH_PER_M;
  
  const yMeters = (yVhFromBottom - 40) / VH_PER_M;
  
  return { x: xMeters, y: yMeters };
}


function snapToGrid(x, y) {
  return {
    x: Math.floor(x / BLOCK_SIZE) * BLOCK_SIZE,
    y: Math.floor(y / BLOCK_SIZE) * BLOCK_SIZE
  };
}

function getBlockAt(gridX, gridY) {
  return blocks.find(b => b.gridX === gridX && b.gridY === gridY);
}

function placeBlock(gridX, gridY, blockType) {
  if (getBlockAt(gridX, gridY)) return; 
  
  const blockData = BlockData[blockType];
  const el = document.createElement("img");
  el.src = blockData.texture;
  el.className = "placed-block";
  el.draggable = false;
  const sizeVh = BLOCK_SIZE * VH_PER_M;
  el.style.width = `${sizeVh}vh`;
  el.style.height = `${sizeVh}vh`;
  
  const screen = physicsToScreen(gridX, gridY);
  el.style.left = screen.left + "vh";
  el.style.bottom = screen.bottom + "vh";
  
  projectileLayerEl.appendChild(el);
  
  blocks.push({
    gridX,
    gridY,
    type: blockType,
    restitution: blockData.restitution,
    element: el
  });
}

function removeBlockAt(gridX, gridY) {
  const idx = blocks.findIndex(b => b.gridX === gridX && b.gridY === gridY);
  if (idx !== -1) {
    blocks[idx].element.remove();
    blocks.splice(idx, 1);
  }
}

function findNearestProjectile(x, y, maxDist = 2) {
  let nearest = null;
  let minDistSq = maxDist * maxDist;
  
  for (const p of projectiles) {
    if (!p.alive) continue;
    const dx = p.x - x;
    const dy = p.y - y;
    const distSq = dx * dx + dy * dy;
    if (distSq < minDistSq) {
      minDistSq = distSq;
      nearest = p;
    }
  }
  return nearest;
}

function updateGridHighlight(evt) {
  if (!ctrlPressed) {
    if (gridHighlightEl) {
      gridHighlightEl.style.display = "none";
    }
    return;
  }
  
  if (!gridHighlightEl) {
    gridHighlightEl = document.createElement("div");
    gridHighlightEl.className = "grid-highlight";
    projectileLayerEl.appendChild(gridHighlightEl);
  }
  
  const pos = getMousePositionMeters(evt);
  const grid = snapToGrid(pos.x, pos.y);
  const screen = physicsToScreen(grid.x, grid.y);
  const sizeVh = BLOCK_SIZE * VH_PER_M;
  
  gridHighlightEl.style.left = screen.left + "vh";
  gridHighlightEl.style.bottom = screen.bottom + "vh";
  gridHighlightEl.style.width = sizeVh + "vh";
  gridHighlightEl.style.height = sizeVh + "vh";
  gridHighlightEl.style.display = "block";
}

function startAim(evt) {
  
  if (evt.button !== 0) return;

  
  if (evt.target.closest && evt.target.closest("#biome-controls, #hotbar, #controls-panel, #controls-help")) return;
  
  const pos = getMousePositionMeters(evt);
  
  
  if (ctrlPressed) {
    const grid = snapToGrid(pos.x, pos.y);
    if (shiftPressed) {
      placeBlock(grid.x, grid.y, "slime");
    } else {
      placeBlock(grid.x, grid.y, "cobblestone");
    }
    return;
  }
  
  
  if (shiftPressed) {
    const nearest = findNearestProjectile(pos.x, pos.y, 3);
    if (nearest) {
      dragState.isDragging = true;
      dragState.projectile = nearest;
      dragState.offsetX = nearest.x - pos.x;
      dragState.offsetY = nearest.y - pos.y;
      
      nearest.vx = 0;
      nearest.vy = 0;
      return;
    }
  }
  
  
  aimState.isAiming = true;
  aimState.startX = pos.x;
  aimState.startY = pos.y;
  aimState.currentX = pos.x;
  aimState.currentY = pos.y;

  if (!aimState.lineEl) {
    aimState.lineEl = document.createElement("div");
    aimState.lineEl.className = "aim-line";
    projectileLayerEl.appendChild(aimState.lineEl);
  }
  velocityDisplayEl.style.display = "block";
  updateAimLine();
}

function updateAim(evt) {
  updateGridHighlight(evt);
  
  
  if (dragState.isDragging && dragState.projectile) {
    const pos = getMousePositionMeters(evt);
    dragState.projectile.x = pos.x + dragState.offsetX;
    dragState.projectile.y = pos.y + dragState.offsetY;
    dragState.projectile.vx = 0;
    dragState.projectile.vy = 0;
    
    dragState.projectile.prevX = dragState.projectile.x;
    dragState.projectile.prevY = dragState.projectile.y;
    renderProjectiles();
    return;
  }
  
  if (!aimState.isAiming) return;
  const pos = getMousePositionMeters(evt);
  aimState.currentX = pos.x;
  aimState.currentY = pos.y;
  updateAimLine();
}

function endAim(evt) {
  
  if (dragState.isDragging) {
    dragState.isDragging = false;
    dragState.projectile = null;
    return;
  }
  
  if (!aimState.isAiming) return;
  const pos = getMousePositionMeters(evt);
  aimState.currentX = pos.x;
  aimState.currentY = pos.y;

  
  
  const dxMeters = aimState.currentX - aimState.startX;
  const dyMeters = aimState.currentY - aimState.startY;
  const SPEED_MULTIPLIER = launchPower; 

  const item = ItemData[currentItemKey];
  const biome = BiomeData[currentBiomeKey];
  const vx = dxMeters * SPEED_MULTIPLIER;
  const vy = dyMeters * SPEED_MULTIPLIER;

  
  
  const floorLevel = biome.isWater ? -10 : 0;
  
  
  let spawnY = aimState.startX === aimState.currentX && aimState.startY === aimState.currentY ? aimState.startY : aimState.startY;
  spawnY = Math.max(floorLevel + item.radius, spawnY);
  
  spawnProjectile(aimState.startX, spawnY, vx, vy, item, biome);

  aimState.isAiming = false;
  if (aimState.lineEl) {
    aimState.lineEl.remove();
    aimState.lineEl = null;
  }
  velocityDisplayEl.style.display = "none";
}

function updateAimLine() {
  if (!aimState.lineEl) return;
  const start = physicsToScreen(aimState.startX, aimState.startY);
  const end = physicsToScreen(aimState.currentX, aimState.currentY);
  const dx = end.left - start.left;
  const dy = end.bottom - start.bottom;
  const length = Math.sqrt(dx * dx + dy * dy);
  const angleDeg = Math.atan2(dy, dx) * (180 / Math.PI);

  
  const dxMeters = aimState.currentX - aimState.startX;
  const dyMeters = aimState.currentY - aimState.startY;
  const SPEED_MULTIPLIER = launchPower;
  const vx = dxMeters * SPEED_MULTIPLIER;
  const vy = dyMeters * SPEED_MULTIPLIER;
  const totalV = Math.hypot(vx, vy);

  velocityDisplayEl.innerText = `Velocity: ${totalV.toFixed(2)} m/s\nVx: ${vx.toFixed(2)} m/s\nVy: ${vy.toFixed(2)} m/s`;

  aimState.lineEl.style.left = start.left + "vh";
  aimState.lineEl.style.bottom = start.bottom + "vh";
  aimState.lineEl.style.width = length + "vh";
  aimState.lineEl.style.transform = `rotate(${-angleDeg}deg)`;
}

let nextProjectileId = 1;

function spawnProjectile(x, y, vx, vy, item, biome) {
  const el = document.createElement("img");
  el.src = item.texture;
  el.className = "projectile";
  el.draggable = false;
  const diameterVh = item.radius * 2 * VH_PER_M;
  el.style.width = `${diameterVh}vh`;
  el.style.height = `${diameterVh}vh`;
  projectileLayerEl.appendChild(el);

  const p = {
    id: nextProjectileId++,
    itemKey: currentItemKey,
    biomeKey: currentBiomeKey,
    mass: item.mass,
    radius: item.radius,
    restitution: item.restitution,
    dragCoefficient: item.dragCoefficient,
    x,
    y,
    prevX: x - vx * DT, 
    prevY: y - vy * DT, 
    vx,
    vy,
    angle: 0,
    element: el,
    alive: true,
    wasInWater: false
  };
  projectiles.push(p);
}

function spawnSplash(x, y, mass, speed, radius = 0.5) {
  const energy = 0.5 * mass * speed * speed;
  const count = Math.min(150, Math.max(15, Math.floor(10 * Math.sqrt(energy))));
  
  for (let i = 0; i < count; i++) {
    const el = document.createElement("div");
    el.className = "particle";
    const size = Math.random() * 0.5 + 0.2; // vh
    el.style.width = size + "vh";
    el.style.height = size + "vh";
    projectileLayerEl.appendChild(el);

    particles.push({
      x: x + (Math.random() - 0.5) * (radius * 3), 
      y: y,
      vx: (Math.random() - 0.5) * 8 * (speed / 5), 
      vy: Math.random() * 8 + 3, 
      life: 1.0 + Math.random() * 0.5, 
      element: el
    });
  }
}

function updateParticles(dt) {
  for (let i = particles.length - 1; i >= 0; i--) {
    const p = particles[i];
    p.life -= dt;
    if (p.life <= 0) {
      p.element.remove();
      particles.splice(i, 1);
      continue;
    }
    
    
    p.vy -= 9.81 * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    
    
    const screen = physicsToScreen(p.x, p.y);
    p.element.style.left = screen.left + "vh";
    p.element.style.bottom = screen.bottom + "vh";
    p.element.style.opacity = p.life;
  }
}


function updatePhysics(dt) {
  const biome = BiomeData[currentBiomeKey];
  const g = biome.gravity;
  const airDensity = biome.airDensity;
  const WATER_DENSITY = 1000; 

  let maxSpeed = 0;
  for (const p of projectiles) {
    if (!p.alive) continue;
    const spd = Math.hypot(p.vx, p.vy);
    if (spd > maxSpeed) maxSpeed = spd;
  }
  const SUBSTEPS = Math.max(8, Math.min(128, Math.ceil(maxSpeed / 8)));
  const subDt = dt / SUBSTEPS;

  for (let step = 0; step < SUBSTEPS; step++) {
    for (const p of projectiles) {
      if (!p.alive) continue;

      let Fx = 0;
      let Fy = 0;
      Fy -= p.mass * g;

      let density = airDensity;
      let submergedRatio = 0;
      let inWater = false;

      if (biome.isWater) {
        const waterLevel = biome.waterLevelMeters;
        const bottomY = p.y - p.radius;
        const topY = p.y + p.radius;

        if (bottomY < waterLevel) {
          inWater = true;
          const totalVolume = (4 / 3) * Math.PI * Math.pow(p.radius, 3);
          if (topY <= waterLevel) {
              submergedRatio = 1.0;
          } else {
              const h = waterLevel - bottomY;
              const submergedVolume = (Math.PI * h * h / 3) * (3 * p.radius - h);
              submergedRatio = submergedVolume / totalVolume;
          }
          
          const buoyancy = WATER_DENSITY * (totalVolume * submergedRatio) * g;
          Fy += buoyancy;
          
          density = airDensity * (1 - submergedRatio) + WATER_DENSITY * submergedRatio;
        }
      }

      if (inWater && !p.wasInWater) {
          const speed = Math.hypot(p.vx, p.vy);
          if (speed > 0.5) {
              spawnSplash(p.x, biome.waterLevelMeters, p.mass, speed, p.radius);
          }
      }
      p.wasInWater = inWater;

      const speed = Math.hypot(p.vx, p.vy);
      
      if (speed > 0.0001) {
        const area = Math.PI * p.radius * p.radius;
        const dragCoeff = p.dragCoefficient;
        
        const k = (0.5 * density * dragCoeff * area) / p.mass;
        
        const decayFactor = 1 / (1 + k * speed * subDt);
        
        p.vx *= decayFactor;
        p.vy *= decayFactor;
      }

      const ax = Fx / p.mass;
      const ay = Fy / p.mass;

      if (integrationMethod === "verlet") {
        const tempX = p.x;
        const tempY = p.y;
        p.vx += ax * subDt;
        p.vy += ay * subDt;
        p.x += p.vx * subDt;
        p.y += p.vy * subDt;
        p.prevX = tempX;
        p.prevY = tempY;
      } else {
        p.vx += ax * subDt;
        p.vy += ay * subDt;
        p.x += p.vx * subDt;
        p.y += p.vy * subDt;
      }

      const groundLevel = biome.isWater ? -10 : 0;
      const groundY = groundLevel + p.radius;
      
      if (p.y < groundY) {
        p.y = groundY;
        p.angle += (p.vx * subDt / p.radius) * (180 / Math.PI);
        
        if (p.vy < 0) {
          let restitution = p.restitution;
          if (currentBiomeKey === "slime") restitution *= biome.groundBounce;
          p.vy = -p.vy * restitution;
          if (Math.abs(p.vy) < 0.2) p.vy = 0;
          if (integrationMethod === "verlet") p.prevY = p.y - p.vy * subDt;
        }
        if (Math.abs(p.vx) < 0.05) p.vx = 0;
        if (integrationMethod === "verlet") p.prevX = p.x - p.vx * subDt;
      }
      
      for (const block of blocks) {
        const blockLeft = block.gridX;
        const blockRight = block.gridX + BLOCK_SIZE;
        const blockBottom = block.gridY;
        const blockTop = block.gridY + BLOCK_SIZE;
        
        const closestX = Math.max(blockLeft, Math.min(p.x, blockRight));
        const closestY = Math.max(blockBottom, Math.min(p.y, blockTop));
        
        const dx = p.x - closestX;
        const dy = p.y - closestY;
        const distSq = dx * dx + dy * dy;
        
        if (distSq < p.radius * p.radius) {
          const dist = Math.sqrt(distSq);
          if (dist === 0) {
            p.y = blockTop + p.radius;
            continue;
          }
          
          const nx = dx / dist;
          const ny = dy / dist;
          
          const overlap = p.radius - dist;
          p.x += nx * overlap;
          p.y += ny * overlap;
          
          const velAlongNormal = p.vx * nx + p.vy * ny;
          if (velAlongNormal < 0) {
            const e = Math.max(p.restitution, block.restitution);
            p.vx -= (1 + e) * velAlongNormal * nx;
            p.vy -= (1 + e) * velAlongNormal * ny;
            
            const tx = -ny;
            const ty = nx;
            const velAlongTangent = p.vx * tx + p.vy * ty;
            p.vx -= 0.02 * velAlongTangent * tx;
            p.vy -= 0.02 * velAlongTangent * ty;
          }
          
          p.angle += (p.vx * subDt / p.radius) * (180 / Math.PI);
          
          if (integrationMethod === "verlet") {
            p.prevX = p.x - p.vx * subDt;
            p.prevY = p.y - p.vy * subDt;
          }
        }
      }
    }
  }

  projectiles = projectiles.filter(p => {
    if (!p.alive) return false;
    if (p.x < -10 || p.x > 200 || p.y < -15 || p.y > 200) {
      p.element.remove();
      return false;
    }
    return true;
  });
}

function resolveCollisions() {
  for (let i = 0; i < projectiles.length; i++) {
    for (let j = i + 1; j < projectiles.length; j++) {
      const p1 = projectiles[i];
      const p2 = projectiles[j];
      if (!p1.alive || !p2.alive) continue;

      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      const distSq = dx * dx + dy * dy;
      const minDist = p1.radius + p2.radius;

      if (distSq < minDist * minDist) {
        const dist = Math.sqrt(distSq);
        if (dist === 0) continue;

        const overlap = minDist - dist;
        const nx = dx / dist;
        const ny = dy / dist;
        
        const totalMass = p1.mass + p2.mass;
        const m1Ratio = p2.mass / totalMass;
        const m2Ratio = p1.mass / totalMass;

        p1.x -= nx * overlap * m1Ratio;
        p1.y -= ny * overlap * m1Ratio;
        p2.x += nx * overlap * m2Ratio;
        p2.y += ny * overlap * m2Ratio;

        const dvx = p2.vx - p1.vx;
        const dvy = p2.vy - p1.vy;
        const velAlongNormal = dvx * nx + dvy * ny;

        if (velAlongNormal > 0) continue;

        const e = Math.max(p1.restitution, p2.restitution);

        let jVal = -(1 + e) * velAlongNormal;
        jVal /= (1 / p1.mass + 1 / p2.mass);

        const impulseX = jVal * nx;
        const impulseY = jVal * ny;

        p1.vx -= impulseX / p1.mass;
        p1.vy -= impulseY / p1.mass;
        p2.vx += impulseX / p2.mass;
        p2.vy += impulseY / p2.mass;
      }
    }
  }
}

function renderProjectiles() {
  for (const p of projectiles) {
    const screen = physicsToScreen(p.x, p.y);
    p.element.style.left = screen.left + "vh";
    p.element.style.bottom = screen.bottom + "vh";
    p.element.style.transform = `translate(-50%, 50%) rotate(${p.angle}deg)`;
  }
}

let lastTimestamp = null;

function loop(timestamp) {
  if (lastTimestamp == null) lastTimestamp = timestamp;
  const elapsedMs = timestamp - lastTimestamp;
  lastTimestamp = timestamp;
  const dt = DT * timeScale; 

  if (!paused) {
    updatePhysics(dt);
    resolveCollisions();
    updateParticles(dt);
  }
  renderProjectiles(); 

  requestAnimationFrame(loop);
}


function setBiome(key) {
  if (!BiomeData[key]) return;
  currentBiomeKey = key;
  document.body.className = key;
  groundEl.className = `ground ${key}`;
}


function setItem(key) {
  if (!ItemData[key]) return;
  currentItemKey = key;
  document.querySelectorAll(".item-button").forEach(btn => {
    btn.classList.toggle("selected", btn.dataset.item === key);
  });
}

function removeNearestProjectile(x, y) {
  let nearest = null;
  let minDistSq = Infinity;
  
  for (const p of projectiles) {
    const dx = p.x - x;
    const dy = p.y - y;
    const distSq = dx * dx + dy * dy;
    if (distSq < minDistSq) {
      minDistSq = distSq;
      nearest = p;
    }
  }
  
  
  if (nearest && minDistSq < 4) {
    nearest.alive = false;
    nearest.element.remove();
    projectiles = projectiles.filter(p => p.alive);
  }
}

function setupUI() {
  
  const imagesToPreload = [
    "grass.png", "dirt.png", "ngrass.png", "netherrack.png", 
    "endstone.png", "water.gif", "slime.png", "cobblestone.png",
    "cannonball.png", "beachball.png", "snowball.png", "slimeball.png"
  ];
  imagesToPreload.forEach(src => {
    const img = new Image();
    img.src = src;
  });

  document.querySelectorAll(".biome-button").forEach(btn => {
    btn.addEventListener("click", () => {
      const key = btn.dataset.biome;
      setBiome(key);
    });
  });

  document.querySelectorAll(".item-button").forEach(btn => {
    btn.addEventListener("click", () => {
      setItem(btn.dataset.item);
    });
  });

  pauseButtonEl.addEventListener("click", () => {
    paused = !paused;
    pauseButtonEl.textContent = paused ? "▶" : "❚❚";
  });
  
  
  speedSlider.addEventListener("input", (e) => {
    timeScale = parseFloat(e.target.value);
    speedVal.textContent = timeScale.toFixed(1);
  });
  
  powerSlider.addEventListener("input", (e) => {
    launchPower = parseFloat(e.target.value);
    powerVal.textContent = launchPower.toFixed(1);
  });
  
  integratorSelect.addEventListener("change", (e) => {
    integrationMethod = e.target.value;
  });

  worldEl.addEventListener("mousedown", startAim);
  worldEl.addEventListener("mousemove", updateAim);
  window.addEventListener("mouseup", endAim);
  
  
  window.addEventListener("contextmenu", (e) => {
    e.preventDefault();
    const pos = getMousePositionMeters(e);
    
    
    if (ctrlPressed) {
      const grid = snapToGrid(pos.x, pos.y);
      removeBlockAt(grid.x, grid.y);
      return;
    }
    
    
    removeNearestProjectile(pos.x, pos.y);
  });
  
  
  window.addEventListener("keydown", (e) => {
    if (e.key === "Control") {
      ctrlPressed = true;
      updateGridHighlight(lastMouseEvent);
    }
    if (e.key === "Shift") {
      shiftPressed = true;
    }
    
    if (e.shiftKey && e.key.toLowerCase() === 'c') {
      projectiles.forEach(p => p.element.remove());
      projectiles = [];
      particles.forEach(p => p.element.remove());
      particles = [];
      blocks.forEach(b => b.element.remove());
      blocks = [];
    }
    if (e.code === 'Space' || e.key.toLowerCase() === 'k') {
      e.preventDefault();
      paused = !paused;
      pauseButtonEl.textContent = paused ? "▶" : "❚❚";
    }
  });
  
  window.addEventListener("keyup", (e) => {
    if (e.key === "Control") {
      ctrlPressed = false;
      updateGridHighlight(lastMouseEvent);
    }
    if (e.key === "Shift") {
      shiftPressed = false;
    }
  });
  
  
  worldEl.addEventListener("mousemove", (e) => {
    lastMouseEvent = e;
  });
}

let lastMouseEvent = { clientX: 0, clientY: 0 };

window.addEventListener("load", () => {
  setupUI();
  setBiome("earth");
  setItem("iron");
  requestAnimationFrame(loop);
});

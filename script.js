"use strict";

const SCALE_M_PER_VH = 1 / 8; 
const VH_PER_M = 8;
const GROUND_HEIGHT_VH = 40;
const OCEAN_FLOOR_Y_M = -GROUND_HEIGHT_VH / VH_PER_M;
const TARGET_FPS = 60;
const DT = 1 / TARGET_FPS; 
const BLOCK_SIZE = 1.0; 
let WORLD_WIDTH_M = 100 / VH_PER_M * (window.innerWidth / window.innerHeight); 
var PXM=-1;//for fluidsim

const paramsString = window.location.search
const urlParams = new URLSearchParams(paramsString)
const watertype = urlParams.get("watertype") === "true" ? true : false

if (!watertype){
  window.FLUID_GRID_SCALE = 3;//canvas is 1/n screen size
window.FLUID_PARTICLE_RADIUS = 0.06;
window.WATER_PARTICLE_DENSITY = 1000;
window.FLUID_REST_DENSITY = 8;
window.FLUID_RENDER_PIXEL_SIZE = 4;
window.FLUID_RENDER_USE_BLOB_KERNEL = false;
}else{
   window.FLUID_GRID_SCALE = 1;//canvas is 1/n screen size
window.FLUID_PARTICLE_RADIUS = 0.06;
window.WATER_PARTICLE_DENSITY = 1000;
window.FLUID_REST_DENSITY = 8;
window.FLUID_RENDER_PIXEL_SIZE = 4;
window.FLUID_RENDER_USE_BLOB_KERNEL = true;
}
const FLUID_GRID_SCALE = window.FLUID_GRID_SCALE;
const FLUID_PARTICLE_RADIUS = window.FLUID_PARTICLE_RADIUS;
const WATER_PARTICLE_DENSITY = window.WATER_PARTICLE_DENSITY;
const FLUID_REST_DENSITY = window.FLUID_REST_DENSITY;
const FLUID_RENDER_PIXEL_SIZE = window.FLUID_RENDER_PIXEL_SIZE;
const FLUID_RENDER_USE_BLOB_KERNEL = window.FLUID_RENDER_USE_BLOB_KERNEL;
const _LE_TEST = new Uint8Array(new Uint32Array([0x0a0b0c0d]).buffer)[0] === 0x0d;
function packRGBA32(r, g, b, a) {
  return _LE_TEST
    ? ((a & 255) << 24) | ((b & 255) << 16) | ((g & 255) << 8) | (r & 255)
    : ((r & 255) << 24) | ((g & 255) << 16) | ((b & 255) << 8) | (a & 255);
}

function buildSquareKernelOffsets(size) {
  const s = Math.max(1, size | 0);
  const start = -Math.floor(s / 2);
  const xs = new Int8Array(s * s);
  const ys = new Int8Array(s * s);
  let k = 0;
  for (let y = 0; y < s; y++) {
    for (let x = 0; x < s; x++) {
      xs[k] = start + x;
      ys[k] = start + y;
      k++;
    }
  }
  return { xs, ys };
}

function buildBlobKernelOffsets() {
  const pts = [
    [-2, -1], [-2, 0], [-2, 1],
    [-1, -2], [-1, -1], [-1, 0], [-1, 1], [-1, 2],
    [0, -2], [0, -1], [0, 0], [0, 1], [0, 2],
    [1, -2], [1, -1], [1, 0], [1, 1], [1, 2],
    [2, -1], [2, 0], [2, 1]
  ];
  const xs = new Int8Array(pts.length);
  const ys = new Int8Array(pts.length);
  for (let i = 0; i < pts.length; i++) {
    xs[i] = pts[i][0];
    ys[i] = pts[i][1];
  }
  return { xs, ys };
}

const { xs: FLUID_RENDER_OFFSETS_X, ys: FLUID_RENDER_OFFSETS_Y } =
  FLUID_RENDER_USE_BLOB_KERNEL
    ? buildBlobKernelOffsets()
    : buildSquareKernelOffsets(FLUID_RENDER_PIXEL_SIZE);
window.getPixelColor=(x,y)=>{return [58,92,166,220]}
async function defGetPixelColor(){const watert = await fetch("/projectilemania/water.json")
if (!watert.ok){
  console.log("warning:/water.json not fetched")
  window.getPixelColor=(x,y)=>{return [58,92,166,220]}
} else{
  const waterdata = await watert.json()
  let wlen=waterdata.length
  console.log(waterdata,wlen)
  window.getPixelColor=(x,y)=>{
    let xmod=y%wlen //90deg
    let ymod=x%wlen
    let upr=waterdata[xmod][ymod]
    return [upr[0]*0.85,upr[1]*0.85,upr[2]*0.95,220]
  }
  try { if (fluidSim) fluidSim.updateWaterTexture(); } catch {}
}}
defGetPixelColor()
class FluidSimulation {
  constructor(canvas, ctx) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.particles = []; 
    this.gravity = 15;
    this.cellSize = 0.25; 
    this.gridWidth = 0;
    this.gridHeight = 0;
    this.u = null; 
    this.v = null; 
    this.uWeight = null;
    this.vWeight = null;
    this.uOld = null;
    this.vOld = null;
    this.cellType = null;//0-air,1=fluid,2=solid
    this.cellCount = null;
    this.pressure = null;
    this.imageData = null;
    this.imageData32 = null;
    this.waterTex32 = null;
    this.mToCanvas = 1;
    this.canvasY0 = 0;
    this._cellTypePrev = null;
    this._relaxHead = null;
    this._relaxNext = null;
    this._relaxNx = 0;
    this._relaxNy = 0;
    this.resize();
  }
  resize() {
    this.canvas.width = Math.floor(window.innerWidth / FLUID_GRID_SCALE);
    this.canvas.height = Math.floor(window.innerHeight / FLUID_GRID_SCALE);
    WORLD_WIDTH_M = 100 / VH_PER_M * (window.innerWidth / window.innerHeight);
    const worldHeightM = 100 / VH_PER_M;
    this.worldMinX = 0;
    this.worldMinY = -12;
    this.gridWidth = Math.ceil(WORLD_WIDTH_M / this.cellSize) + 4;
    this.gridHeight = Math.ceil((worldHeightM + 12) / this.cellSize) + 4; 
    const gw = this.gridWidth;
    const gh = this.gridHeight;
    this.u = new Float32Array((gw + 1) * gh);
    this.v = new Float32Array(gw * (gh + 1));
    this.uWeight = new Float32Array((gw + 1) * gh);
    this.vWeight = new Float32Array(gw * (gh + 1));
    this.uOld = new Float32Array((gw + 1) * gh);
    this.vOld = new Float32Array(gw * (gh + 1));
    this.cellType = new Uint8Array(gw * gh);
    this.cellCount = new Float32Array(gw * gh);
    this.pressure = new Float32Array(gw * gh);
    this.imageData = this.ctx.createImageData(this.canvas.width, this.canvas.height);
    this.imageData32 = new Uint32Array(this.imageData.data.buffer);
    this.waterTex32 = new Uint32Array(this.canvas.width * this.canvas.height);
    this._cellTypePrev = new Uint8Array(gw * gh);
    const vhPx = window.innerHeight / 100;
    this.mToCanvas = (VH_PER_M * vhPx) / FLUID_GRID_SCALE;
    this.canvasY0 = this.canvas.height - ((GROUND_HEIGHT_VH * vhPx) / FLUID_GRID_SCALE);
    this.updateWaterTexture();
  }
  worldToGrid(x, y) {
    return {
      i: Math.floor((x - this.worldMinX) / this.cellSize),
      j: Math.floor((y - this.worldMinY) / this.cellSize)
    };
  }
  gridToWorld(i, j) {
    return {
      x: i * this.cellSize + this.worldMinX,
      y: j * this.cellSize + this.worldMinY
    };
  }
  worldToCanvas(x, y) {
    const px = Math.floor(x * this.mToCanvas);
    const py = Math.floor(this.canvasY0 - (y * this.mToCanvas));
    return { px, py };
  }
  updateWaterTexture() {
    if (!this.waterTex32) return;
    const w = this.canvas.width;
    const h = this.canvas.height;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        const col = window.getPixelColor(x, y);
        this.waterTex32[y * w + x] = packRGBA32(col[0], col[1], col[2], col[3]);
      }
    }
  }
  addParticle(x, y, vx = 0, vy = 0) {
    this.particles.push({ x, y, vx, vy });
  }
  addParticlesInRadius(cx, cy, radius, count = 20) {
    for (let i = 0; i < count; i++) {
      const angle = Math.random() * Math.PI * 2;
      const r = Math.sqrt(Math.random()) * radius;
      const x = cx + Math.cos(angle) * r;
      const y = cy + Math.sin(angle) * r;
      const vx = (Math.random() - 0.5) * 0.5;
      const vy = (Math.random() - 0.5) * 0.5;
      this.addParticle(x, y, vx, vy);
    }
  }
  addParticlesInCell(gridX, gridY) {
    const particlesPerSide = 3;
    const spacing = BLOCK_SIZE / particlesPerSide;
    for (let i = 0; i < particlesPerSide; i++) {
      for (let j = 0; j < particlesPerSide; j++) {
        const x = gridX + (i + 0.5) * spacing;
        const y = gridY + (j + 0.5) * spacing;
        this.addParticle(x, y, 0, 0);
      }
    }
  }
  removeParticlesInRadius(cx, cy, radius) {
    const r2 = radius * radius;
    this.particles = this.particles.filter(p => {
      const dx = p.x - cx;
      const dy = p.y - cy;
      return dx * dx + dy * dy > r2;
    });
  }
  removeParticlesInCell(gridX, gridY) {
    this.particles = this.particles.filter(p => {
      return !(p.x >= gridX && p.x < gridX + BLOCK_SIZE && 
               p.y >= gridY && p.y < gridY + BLOCK_SIZE);
    });
  }
  spawnOceanWater() {
    const startY = -9.5;
    const endY = 0;
    const startX = -1;
    const endX = 20;
    const spacing = 0.2;
    for (let y = startY; y < endY; y += spacing) {
      for (let x = startX; x < endX; x += spacing) {
        const jx = x + (Math.random() - 0.5) * spacing * 0.5;
        const jy = y + (Math.random() - 0.5) * spacing * 0.5;
        this.addParticle(jx, jy, 0, 0);
      }
    }
  }
  removeParticlesBelowGround(groundY = 0) {
    this.particles = this.particles.filter(p => p.y >= groundY - 0.5);
  }
  lerp(a, b, t) { return a + (b - a) * t; }
  clamp(v, min, max) { return Math.max(min, Math.min(max, v)); }
  sampleU(x, y) {
    const i = (x - this.worldMinX) / this.cellSize;
    const j = (y - this.worldMinY) / this.cellSize - 0.5;
    const i0 = Math.floor(i), j0 = Math.floor(j);
    const s = i - i0, t = j - j0;
    const gw = this.gridWidth + 1, gh = this.gridHeight;
    const getU = (ii, jj) => (ii < 0 || ii >= gw || jj < 0 || jj >= gh) ? 0 : this.u[jj * gw + ii];
    return this.lerp(this.lerp(getU(i0, j0), getU(i0+1, j0), s), this.lerp(getU(i0, j0+1), getU(i0+1, j0+1), s), t);
  }
  sampleV(x, y) {
    const i = (x - this.worldMinX) / this.cellSize - 0.5;
    const j = (y - this.worldMinY) / this.cellSize;
    const i0 = Math.floor(i), j0 = Math.floor(j);
    const s = i - i0, t = j - j0;
    const gw = this.gridWidth, gh = this.gridHeight + 1;
    const getV = (ii, jj) => (ii < 0 || ii >= gw || jj < 0 || jj >= gh) ? 0 : this.v[jj * gw + ii];
    return this.lerp(this.lerp(getV(i0, j0), getV(i0+1, j0), s), this.lerp(getV(i0, j0+1), getV(i0+1, j0+1), s), t);
  }
  sampleUOld(x, y) {
    const i = (x - this.worldMinX) / this.cellSize;
    const j = (y - this.worldMinY) / this.cellSize - 0.5;
    const i0 = Math.floor(i), j0 = Math.floor(j);
    const s = i - i0, t = j - j0;
    const gw = this.gridWidth + 1, gh = this.gridHeight;
    const getU = (ii, jj) => (ii < 0 || ii >= gw || jj < 0 || jj >= gh) ? 0 : this.uOld[jj * gw + ii];
    return this.lerp(this.lerp(getU(i0, j0), getU(i0+1, j0), s), this.lerp(getU(i0, j0+1), getU(i0+1, j0+1), s), t);
  }
  sampleVOld(x, y) {
    const i = (x - this.worldMinX) / this.cellSize - 0.5;
    const j = (y - this.worldMinY) / this.cellSize;
    const i0 = Math.floor(i), j0 = Math.floor(j);
    const s = i - i0, t = j - j0;
    const gw = this.gridWidth, gh = this.gridHeight + 1;
    const getV = (ii, jj) => (ii < 0 || ii >= gw || jj < 0 || jj >= gh) ? 0 : this.vOld[jj * gw + ii];
    return this.lerp(this.lerp(getV(i0, j0), getV(i0+1, j0), s), this.lerp(getV(i0, j0+1), getV(i0+1, j0+1), s), t);
  }
  particlesToGrid() {
    const gw = this.gridWidth, gh = this.gridHeight;
    this.u.fill(0); this.v.fill(0);
    this.uWeight.fill(0); this.vWeight.fill(0);
    this.cellType.fill(0);
    this.cellCount.fill(0);
    for (const p of this.particles) {
      const ci = Math.floor((p.x - this.worldMinX) / this.cellSize);
      const cj = Math.floor((p.y - this.worldMinY) / this.cellSize);
      if (ci >= 0 && ci < gw && cj >= 0 && cj < gh) {
        const cidx = cj * gw + ci;
        this.cellType[cidx] = 1;
        this.cellCount[cidx] += 1;
      }
      {
        const i = (p.x - this.worldMinX) / this.cellSize;
        const j = (p.y - this.worldMinY) / this.cellSize - 0.5;
        const i0 = Math.floor(i), j0 = Math.floor(j);
        for (let di = 0; di <= 1; di++) {
          for (let dj = 0; dj <= 1; dj++) {
            const ii = i0 + di, jj = j0 + dj;
            if (ii < 0 || ii > gw || jj < 0 || jj >= gh) continue;
            const w = (1 - Math.abs(i - ii)) * (1 - Math.abs(j - jj));
            const idx = jj * (gw + 1) + ii;
            this.u[idx] += p.vx * w;
            this.uWeight[idx] += w;
          }
        }
      }
      {
        const i = (p.x - this.worldMinX) / this.cellSize - 0.5;
        const j = (p.y - this.worldMinY) / this.cellSize;
        const i0 = Math.floor(i), j0 = Math.floor(j);
        for (let di = 0; di <= 1; di++) {
          for (let dj = 0; dj <= 1; dj++) {
            const ii = i0 + di, jj = j0 + dj;
            if (ii < 0 || ii >= gw || jj < 0 || jj > gh) continue;
            const w = (1 - Math.abs(i - ii)) * (1 - Math.abs(j - jj));
            const idx = jj * gw + ii;
            this.v[idx] += p.vy * w;
            this.vWeight[idx] += w;
          }
        }
      }
    }
    for (let i = 0; i < this.u.length; i++) if (this.uWeight[i] > 0) this.u[i] /= this.uWeight[i];
    for (let i = 0; i < this.v.length; i++) if (this.vWeight[i] > 0) this.v[i] /= this.vWeight[i];
    this.uOld.set(this.u);
    this.vOld.set(this.v);
  }
  expandFluidCells(rings = 1) {
    const gw = this.gridWidth;
    const gh = this.gridHeight;
    for (let r = 0; r < rings; r++) {
      this._cellTypePrev.set(this.cellType);
      const prev = this._cellTypePrev;
      for (let j = 1; j < gh - 1; j++) {
        for (let i = 1; i < gw - 1; i++) {
          const idx = j * gw + i;
          if (prev[idx] !== 1) continue;
          const n1 = idx - 1;
          const n2 = idx + 1;
          const n3 = idx - gw;
          const n4 = idx + gw;
          if (this.cellType[n1] === 0) this.cellType[n1] = 1;
          if (this.cellType[n2] === 0) this.cellType[n2] = 1;
          if (this.cellType[n3] === 0) this.cellType[n3] = 1;
          if (this.cellType[n4] === 0) this.cellType[n4] = 1;
        }
      }
    }
  }
  markSolids(blocks, groundY) {
    const gw = this.gridWidth, gh = this.gridHeight;
    for (let j = 0; j < gh; j++) {
      for (let i = -10; i < gw; i++) {
        if (i === -10 || i === gw - 1) {
          this.cellType[j * gw + i] = 2;
          continue;
        }
        const cellY = this.gridToWorld(i, j).y;
        if (cellY < groundY) this.cellType[j * gw + i] = 2;
      }
    }
    for (const block of blocks) {
      const minGC = this.worldToGrid(block.gridX, block.gridY);
      const maxGC = this.worldToGrid(block.gridX + BLOCK_SIZE, block.gridY + BLOCK_SIZE);
      for (let j = Math.max(0, minGC.j); j <= Math.min(gh - 1, maxGC.j); j++) {
        for (let i = Math.max(0, minGC.i); i <= Math.min(gw - 1, maxGC.i); i++) {
          this.cellType[j * gw + i] = 2;
        }
      }
    }
  }
  applyGravity(dt, g) {
    const gw = this.gridWidth, gh = this.gridHeight + 1;
    for (let j = 0; j < gh; j++) {
      for (let i = 0; i < gw; i++) {
        this.v[j * gw + i] -= g * dt;
      }
    }
  }
  solvePressure(dt, iterations = 20) {
    const gw = this.gridWidth, gh = this.gridHeight, h = this.cellSize;
    const invDt = dt > 0 ? 1 / dt : 0;
    this.pressure.fill(0);
    const densityStiffness = 0.15;
    for (let iter = 0; iter < iterations; iter++) {
      for (let j = 1; j < gh - 1; j++) {
        for (let i = 1; i < gw - 1; i++) {
          const idx = j * gw + i;
          if (this.cellType[idx] !== 1) continue;
          const uRight = this.u[j * (gw + 1) + i + 1];
          const uLeft = this.u[j * (gw + 1) + i];
          const vTop = this.v[(j + 1) * gw + i];
          const vBottom = this.v[j * gw + i];
          const div = (uRight - uLeft + vTop - vBottom) / h;
          let b = div * invDt;
          const count = this.cellCount ? this.cellCount[idx] : 0;
          if (count > 0) {
            const densityErr = (count - FLUID_REST_DENSITY) / FLUID_REST_DENSITY;
            b += densityStiffness * densityErr * invDt;
          }
          let s = 0;
          if (this.cellType[idx - 1] !== 2) s++;
          if (this.cellType[idx + 1] !== 2) s++;
          if (this.cellType[idx - gw] !== 2) s++;
          if (this.cellType[idx + gw] !== 2) s++;
          if (s === 0) continue;
          const pL = this.cellType[idx - 1] !== 2 ? this.pressure[idx - 1] : 0;
          const pR = this.cellType[idx + 1] !== 2 ? this.pressure[idx + 1] : 0;
          const pB = this.cellType[idx - gw] !== 2 ? this.pressure[idx - gw] : 0;
          const pT = this.cellType[idx + gw] !== 2 ? this.pressure[idx + gw] : 0;
          this.pressure[idx] = (pL + pR + pB + pT - b * h * h) / s;
        }
      }
    }
    for (let j = 1; j < gh - 1; j++) {
      for (let i = 1; i < gw; i++) {
        const idxL = j * gw + i - 1;
        const idxR = j * gw + i;
        if ((this.cellType[idxL] === 1 || this.cellType[idxR] === 1) &&
            this.cellType[idxL] !== 2 && this.cellType[idxR] !== 2) {
          this.u[j * (gw + 1) + i] -= dt * (this.pressure[idxR] - this.pressure[idxL]) / h;
        }
      }
    }
    for (let j = 1; j < gh; j++) {
      for (let i = 1; i < gw - 1; i++) {
        const idxB = (j - 1) * gw + i;
        const idxT = j * gw + i;
        if ((this.cellType[idxB] === 1 || this.cellType[idxT] === 1) &&
            this.cellType[idxB] !== 2 && this.cellType[idxT] !== 2) {
          this.v[j * gw + i] -= dt * (this.pressure[idxT] - this.pressure[idxB]) / h;
        }
      }
    }
  }
  enforceBoundaries(blocks, groundY) {
    const gw = this.gridWidth, gh = this.gridHeight;
    for (let j = 0; j < gh; j++) {
      for (let i = 0; i < gw; i++) {
        if (this.cellType[j * gw + i] === 2) {
          this.u[j * (gw + 1) + i] = 0;
          this.u[j * (gw + 1) + i + 1] = 0;
          this.v[j * gw + i] = 0;
          if (j + 1 <= gh) this.v[(j + 1) * gw + i] = 0;
        }
      }
    }
    for (let j = 0; j < gh; j++) {
      this.u[j * (gw + 1)] = 0;
      this.u[j * (gw + 1) + gw] = 0;
    }
    for (let i = 0; i < gw; i++) {
      this.v[i] = 0;
      this.v[gh * gw + i] = 0;
    }
  }
  gridToParticles(flipRatio = 0.95) {
    for (const p of this.particles) {
      const newU = this.sampleU(p.x, p.y);
      const newV = this.sampleV(p.x, p.y);
      const oldU = this.sampleUOld(p.x, p.y);
      const oldV = this.sampleVOld(p.x, p.y);
      const flipVx = p.vx + (newU - oldU);
      const flipVy = p.vy + (newV - oldV);
      const picVx = newU;
      const picVy = newV;
      p.vx = flipRatio * flipVx + (1 - flipRatio) * picVx;
      p.vy = flipRatio * flipVy + (1 - flipRatio) * picVy;
    }
  }
  advectParticles(dt, blocks, groundY) {
    for (const p of this.particles) {
      const vx1 = p.vx, vy1 = p.vy;
      const midX = p.x + vx1 * dt * 0.5;
      const midY = p.y + vy1 * dt * 0.5;
      const vx2 = this.sampleU(midX, midY);
      const vy2 = this.sampleV(midX, midY);
      p.x += vx2 * dt;
      p.y += vy2 * dt;
      const leftWall = -10;
      const rightWall = WORLD_WIDTH_M;
      if (p.x < leftWall + FLUID_PARTICLE_RADIUS) {
        p.x = leftWall + FLUID_PARTICLE_RADIUS;
        if (p.vx < 0) p.vx = -p.vx * 0.3;
      } else if (p.x > rightWall - FLUID_PARTICLE_RADIUS) {
        p.x = rightWall - FLUID_PARTICLE_RADIUS;
        if (p.vx > 0) p.vx = -p.vx * 0.3;
      }
      if (p.y < groundY + FLUID_PARTICLE_RADIUS) {
        p.y = groundY + FLUID_PARTICLE_RADIUS;
        if (p.vy < 0) p.vy = -p.vy * 0.3;
        else p.vy = 0;
      }
      for (const block of blocks) {
        const bx = block.gridX, by = block.gridY;
        if (p.x > bx && p.x < bx + BLOCK_SIZE && p.y > by && p.y < by + BLOCK_SIZE) {
          const distLeft = p.x - bx, distRight = bx + BLOCK_SIZE - p.x;
          const distBottom = p.y - by, distTop = by + BLOCK_SIZE - p.y;
          const minDist = Math.min(distLeft, distRight, distBottom, distTop);
          if (minDist === distLeft) { p.x = bx - FLUID_PARTICLE_RADIUS; p.vx = Math.min(0, p.vx * -0.3); }
          else if (minDist === distRight) { p.x = bx + BLOCK_SIZE + FLUID_PARTICLE_RADIUS; p.vx = Math.max(0, p.vx * -0.3); }
          else if (minDist === distBottom) { p.y = by - FLUID_PARTICLE_RADIUS; p.vy = Math.min(0, p.vy * -0.3); }
          else { p.y = by + BLOCK_SIZE + FLUID_PARTICLE_RADIUS; p.vy = Math.max(0, p.vy * -0.3); }
        }
      }
      p.y = this.clamp(p.y, groundY + FLUID_PARTICLE_RADIUS, 14);
    }
  }
  relaxParticles(dt, groundY, iterations = 2) {
    if (this.particles.length === 0) return;
    const h = Math.max(FLUID_PARTICLE_RADIUS * 3.5, 0.12);
    const invH = 1 / h;
    const target = FLUID_PARTICLE_RADIUS * 2.0;
    const leftWall = 0;
    const rightWall = WORLD_WIDTH_M;
    for (let iter = 0; iter < iterations; iter++) {
      const nx = Math.max(8, Math.ceil((WORLD_WIDTH_M - this.worldMinX) * invH) + 3);
      const ny = Math.max(8, Math.ceil((16 - this.worldMinY) * invH) + 3);
      const cells = nx * ny;
      if (!this._relaxHead || this._relaxHead.length < cells || this._relaxNx !== nx || this._relaxNy !== ny) {
        this._relaxHead = new Int32Array(cells);
        this._relaxNx = nx;
        this._relaxNy = ny;
      }
      const head = this._relaxHead;
      head.fill(-1);
      if (!this._relaxNext || this._relaxNext.length < this.particles.length) {
        this._relaxNext = new Int32Array(this.particles.length);
      }
      const next = this._relaxNext;
      for (let pi = 0; pi < this.particles.length; pi++) {
        const p = this.particles[pi];
        let ix = Math.floor((p.x - this.worldMinX) * invH) + 1;
        let iy = Math.floor((p.y - this.worldMinY) * invH) + 1;
        if (ix < 1) ix = 1; else if (ix > nx - 2) ix = nx - 2;
        if (iy < 1) iy = 1; else if (iy > ny - 2) iy = ny - 2;
        const cell = iy * nx + ix;
        next[pi] = head[cell];
        head[cell] = pi;
      }
      for (let pi = 0; pi < this.particles.length; pi++) {
        const p = this.particles[pi];
        const oldX = p.x;
        const oldY = p.y;
        let ix = Math.floor((p.x - this.worldMinX) * invH) + 1;
        let iy = Math.floor((p.y - this.worldMinY) * invH) + 1;
        if (ix < 1) ix = 1; else if (ix > nx - 2) ix = nx - 2;
        if (iy < 1) iy = 1; else if (iy > ny - 2) iy = ny - 2;
        let dxSum = 0;
        let dySum = 0;
        let n = 0;
        for (let oy = -1; oy <= 1; oy++) {
          for (let ox = -1; ox <= 1; ox++) {
            const cell = (iy + oy) * nx + (ix + ox);
            for (let qi = head[cell]; qi !== -1; qi = next[qi]) {
              if (qi === pi) continue;
              const q = this.particles[qi];
              const dx = p.x - q.x;
              const dy = p.y - q.y;
              const d2 = dx * dx + dy * dy;
              if (d2 <= 0) continue;
              const d = Math.sqrt(d2);
              if (d >= target) continue;
              const push = (target - d) / target;
              dxSum += (dx / d) * push;
              dySum += (dy / d) * push;
              n++;
            }
          }
        }
        if (n > 0) {
          const strength=0.30;
          p.x += (dxSum / n) * strength;
          p.y += (dySum / n) * strength;
        }
        if (p.x < leftWall + FLUID_PARTICLE_RADIUS) p.x = leftWall + FLUID_PARTICLE_RADIUS;
        if (p.x > rightWall - FLUID_PARTICLE_RADIUS) p.x = rightWall - FLUID_PARTICLE_RADIUS;
        if (p.y < groundY + FLUID_PARTICLE_RADIUS) p.y = groundY + FLUID_PARTICLE_RADIUS;
        if (p.y > 14) p.y = 14;
        if (dt > 0) {
          p.vx += (p.x - oldX) / dt;
          p.vy += (p.y - oldY) / dt;
        }
      }
    }
  }
  interactWithProjectiles(projectiles, dt) {
    const INTERACTION_RADIUS = 0.3;
    for (const proj of projectiles) {
      if (!proj.alive) continue;
      let nearbyCount = 0;
      let avgVx = 0, avgVy = 0;
      let totalWeight = 0;
      for (const p of this.particles) {
        const dx = p.x - proj.x, dy = p.y - proj.y;
        const distSq = dx * dx + dy * dy;
        const maxDist = proj.radius + INTERACTION_RADIUS;
        if (distSq < maxDist * maxDist) {
          const dist = Math.sqrt(distSq);
          const weight = 1 - dist / maxDist;
          nearbyCount++;
          avgVx += p.vx * weight;
          avgVy += p.vy * weight;
          totalWeight += weight;
          if (dist < proj.radius + 0.05 && dist > 0.01) {
            const overlap = proj.radius + 0.05 - dist;
            const nx = dx / dist, ny = dy / dist;
            p.x += nx * overlap * 0.5;
            p.y += ny * overlap * 0.5;
            const projSpeed = Math.sqrt(proj.vx * proj.vx + proj.vy * proj.vy);
            p.vx += nx * projSpeed * 0.3 + proj.vx * 0.2;
            p.vy += ny * projSpeed * 0.3 + proj.vy * 0.2;
          }
        }
      }
      if (nearbyCount > 0 && totalWeight > 0) {
        avgVx /= totalWeight;
        avgVy /= totalWeight;
        const expectedParticles = Math.PI * Math.pow(proj.radius + INTERACTION_RADIUS, 2) * FLUID_REST_DENSITY * 4;
        const submergedRatio = Math.min(1, nearbyCount / expectedParticles);
        if (submergedRatio > 0.05) {
          const volume = (4/3) * Math.PI * Math.pow(proj.radius, 3);
          const buoyancyForce = WATER_PARTICLE_DENSITY * volume * submergedRatio * this.gravity;
          proj.vy += buoyancyForce / proj.mass * dt;
          const relVx = proj.vx - avgVx, relVy = proj.vy - avgVy;
          const relSpeed = Math.sqrt(relVx * relVx + relVy * relVy);
          if (relSpeed > 0.001) {
            const area = Math.PI * proj.radius * proj.radius;
            const dragMagnitude = 0.5 * WATER_PARTICLE_DENSITY * proj.dragCoefficient * area * relSpeed * relSpeed * submergedRatio;
            const dragAccel = Math.min(dragMagnitude / proj.mass, relSpeed / dt * 0.5); 
            proj.vx -= (relVx / relSpeed) * dragAccel * dt;
            proj.vy -= (relVy / relSpeed) * dragAccel * dt;
          }
          proj.inFluid = true;
          proj.fluidSubmergedRatio = submergedRatio;
        } else {
          proj.inFluid = false;
          proj.fluidSubmergedRatio = 0;
        }
      } else {
        proj.inFluid = false;
        proj.fluidSubmergedRatio = 0;
      }
    }
  }
  step(dt) {
    if (this.particles.length === 0) return;
    const blocksArr = typeof blocks !== 'undefined' ? blocks : [];
    const projectilesArr = typeof projectiles !== 'undefined' ? projectiles : [];
    const groundY = typeof currentBiomeKey !== 'undefined' && currentBiomeKey === 'ocean' ? OCEAN_FLOOR_Y_M : 0;
    const maxDt = 0.02;
    const steps = Math.ceil(dt / maxDt);
    const subDt = dt / steps;
    for (let s = 0; s < steps; s++) {
      this.particlesToGrid();
      this.markSolids(blocksArr, groundY);
      this.expandFluidCells(1);
      this.applyGravity(subDt, this.gravity);
      this.enforceBoundaries(blocksArr, groundY);
      this.solvePressure(subDt, 40); 
      this.enforceBoundaries(blocksArr, groundY);
      this.gridToParticles(0.5); 
      this.advectParticles(subDt, blocksArr, groundY);
      this.relaxParticles(subDt, groundY, 2);
      this.interactWithProjectiles(projectilesArr, subDt);
    }
  }
  render() {
    const w = this.canvas.width;
    const h = this.canvas.height;
    const out32 = this.imageData32;
    out32.fill(0);
    const tex32 = this.waterTex32;
    const mToCanvas = this.mToCanvas;
    const canvasY0 = this.canvasY0;
    for (const p of this.particles) {
      const px = Math.floor(p.x * mToCanvas);
      const py = Math.floor(canvasY0 - (p.y * mToCanvas));
      for (let k = 0; k < FLUID_RENDER_OFFSETS_X.length; k++) {
        const x = px + FLUID_RENDER_OFFSETS_X[k];
        const y = py + FLUID_RENDER_OFFSETS_Y[k];
        if (x < 0 || x >= w || y < 0 || y >= h) continue;
        const pi = y * w + x;
        out32[pi] = tex32[pi];
      }
    }
    this.ctx.putImageData(this.imageData, 0, 0);
  }
}
let fluidSim = null;
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
    mass: 2754, 
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
let aKeyPressed = false;
let zKeyPressed = false;
let mouseDown = false;
let gridHighlightEl = null;
let waterRadiusHighlightEl = null;
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
  const showWaterRadius = (aKeyPressed || zKeyPressed) && !shiftPressed;
  const showWaterGrid = (aKeyPressed || zKeyPressed) && shiftPressed;
  
  if (showWaterRadius) {
    if (!waterRadiusHighlightEl) {
      waterRadiusHighlightEl = document.createElement("div");
      waterRadiusHighlightEl.className = "water-radius-highlight";
      projectileLayerEl.appendChild(waterRadiusHighlightEl);
    }
    const pos = getMousePositionMeters(evt);
    const radius = launchPower * 0.5; 
    const screen = physicsToScreen(pos.x, pos.y);
    const radiusVh = radius * VH_PER_M;
    waterRadiusHighlightEl.style.left = screen.left + "vh";
    waterRadiusHighlightEl.style.bottom = screen.bottom + "vh";
    waterRadiusHighlightEl.style.width = (radiusVh * 2) + "vh";
    waterRadiusHighlightEl.style.height = (radiusVh * 2) + "vh";
    waterRadiusHighlightEl.style.display = "block";
    waterRadiusHighlightEl.style.borderColor = aKeyPressed ? "rgba(58, 92, 166, 0.8)" : "rgba(255, 100, 100, 0.8)";
  } else {
    if (waterRadiusHighlightEl) {
      waterRadiusHighlightEl.style.display = "none";
    }
  }
  const showBlockGrid = ctrlPressed || showWaterGrid;
  if (!showBlockGrid) {
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
  gridHighlightEl.style.width = (sizeVh-0.5) + "vh";
  gridHighlightEl.style.height = (sizeVh-0.4) + "vh";
  gridHighlightEl.style.display = "block";
  
  
  if (showWaterGrid) {
    gridHighlightEl.style.borderColor = aKeyPressed ? "rgba(58, 92, 166, 0.8)" : "rgba(255, 100, 100, 0.8)";
  } else {
    gridHighlightEl.style.borderColor = "rgba(255, 255, 255, 0.7)";
  }
}

function handleWaterPlacement(evt) {
  if (!mouseDown) return;
  if (!fluidSim) return;
  
  const pos = getMousePositionMeters(evt);
  
  if (aKeyPressed) {
    
    if (shiftPressed) {
     
      const grid = snapToGrid(pos.x, pos.y);
      fluidSim.addParticlesInCell(grid.x, grid.y);
    } else {
      
      const radius = launchPower * 0.5;
      fluidSim.addParticlesInRadius(pos.x, pos.y, radius, Math.floor(5 + launchPower * 3));
    }
  } else if (zKeyPressed) {
    if (shiftPressed) {
      const grid = snapToGrid(pos.x, pos.y);
      fluidSim.removeParticlesInCell(grid.x, grid.y);
    } else {
      const radius = launchPower * 0.5;
      fluidSim.removeParticlesInRadius(pos.x, pos.y, radius);
    }
  }
}
function startAim(evt) {
  if (evt.button !== 0) return;
  mouseDown = true;
  lastMouseEvent = evt;
  if (evt.target.closest && evt.target.closest("#biome-controls, #hotbar, #controls-panel, #controls-help")) return;
  
  const pos = getMousePositionMeters(evt);
  if (aKeyPressed || zKeyPressed) {
    handleWaterPlacement(evt);
    return;
  }
  if (ctrlPressed) {
    const grid = snapToGrid(pos.x, pos.y);
    if (shiftPressed) {
      placeBlock(grid.x, grid.y, "slime");
    } else {
      placeBlock(grid.x, grid.y, "cobblestone");
    }
    return;
  }
  if (shiftPressed && !aKeyPressed && !zKeyPressed) {
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
  if (mouseDown && (aKeyPressed || zKeyPressed)) {
    handleWaterPlacement(evt);
    return;
  }
  
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
  mouseDown = false;
  
  if (dragState.isDragging) {
    dragState.isDragging = false;
    dragState.projectile = null;
    return;
  }
  if (aKeyPressed || zKeyPressed) {
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
    const size = Math.random() * 0.5 + 0.2; 
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
    const minX = -2;
    const maxX = WORLD_WIDTH_M + 2;
    const minY = -25; 
    const maxY = 80;
    if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY) {
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
let waterPaintAccumulator = 0;
const WATER_PAINT_INTERVAL_S = 0.05;

function loop(timestamp) {
  if (lastTimestamp == null) lastTimestamp = timestamp;
  const elapsedMs = timestamp - lastTimestamp;
  lastTimestamp = timestamp;
  const dt = DT * timeScale; 
  if (mouseDown && (aKeyPressed || zKeyPressed) && lastMouseEvent) {
    try { updateGridHighlight(lastMouseEvent); } catch {}
    waterPaintAccumulator += dt;
    while (waterPaintAccumulator >= WATER_PAINT_INTERVAL_S) {
      waterPaintAccumulator -= WATER_PAINT_INTERVAL_S;
      handleWaterPlacement(lastMouseEvent);
    }
  } else {
    waterPaintAccumulator = 0;
  }

  if (!paused) {
    updatePhysics(dt);
    resolveCollisions();
    updateParticles(dt);
    if (fluidSim) {
      fluidSim.step(dt);
    }
  }
  renderProjectiles(); 
  if (fluidSim) {
    fluidSim.render();
  }

  requestAnimationFrame(loop);
}


function setBiome(key) {
  if (!BiomeData[key]) return;
  currentBiomeKey = key;
  document.body.className = key;
  groundEl.className = `ground ${key}`;
  if (key === "ocean" && fluidSim) {
    fluidSim.particles = [];
    const spacing = 0.12; 
    const floorY = OCEAN_FLOOR_Y_M;
    if (typeof WORLD_WIDTH_M === 'undefined' || WORLD_WIDTH_M <= 0) {
        WORLD_WIDTH_M = 100 / VH_PER_M * (window.innerWidth / window.innerHeight);
    }
    
    for (let x = spacing; x < WORLD_WIDTH_M - spacing; x += spacing) {
      for (let y = floorY + spacing; y < 0; y += spacing) {
        const px = x + (Math.random() - 0.5) * spacing * 0.3;
        const py = y + (Math.random() - 0.5) * spacing * 0.3;
        fluidSim.addParticle(px, py);
      }
    }
  } else if (fluidSim) {
    fluidSim.particles = fluidSim.particles.filter(p => p.y >= 0);
  }
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
    if (e.key.toLowerCase() === "a") {
      aKeyPressed = true;
      updateGridHighlight(lastMouseEvent);
    }
    if (e.key.toLowerCase() === "z") {
      zKeyPressed = true;
      updateGridHighlight(lastMouseEvent);
    }
    
    if (e.shiftKey && e.key.toLowerCase() === 'c') {
      projectiles.forEach(p => p.element.remove());
      projectiles = [];
      particles.forEach(p => p.element.remove());
      particles = [];
      blocks.forEach(b => b.element.remove());
      blocks = [];
      if (fluidSim) {
        fluidSim.particles = [];
      }
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
    if (e.key.toLowerCase() === "a") {
      aKeyPressed = false;
      updateGridHighlight(lastMouseEvent);
    }
    if (e.key.toLowerCase() === "z") {
      zKeyPressed = false;
      updateGridHighlight(lastMouseEvent);
    }
  });
  
  
  worldEl.addEventListener("mousemove", (e) => {
    lastMouseEvent = e;
  });
}

let lastMouseEvent = { clientX: 0, clientY: 0 };

window.addEventListener("load", () => {
  setupUI();
  PXM = window.innerHeight * (100/SCALE_M_PER_VH);
  const fluidCanvas = document.getElementById("fluid-canvas");
  const fluidCtx = fluidCanvas.getContext("2d");
  fluidSim = new FluidSimulation(fluidCanvas, fluidCtx);
  setBiome("earth");
  setItem("iron");
  requestAnimationFrame(loop);
});

window.addEventListener('mousemove', (e) => {
  try{
      ctrlPressed = e.ctrlKey; //apparently the mouse event knows this accurately
      shiftPressed = e.shiftKey;}catch{/*these vars are not defined yet*/}
    
})

window.addEventListener('blur', () => {try{
  //the user ctrl+tabbed out
    ctrlPressed = false;
    shiftPressed = false
}catch{}});

/* fluid sim resize handler */
window.addEventListener("resize", () => {
    PXM = window.innerHeight * (100/SCALE_M_PER_VH);
    if (fluidSim) {
        fluidSim.resize();
    }
});

window.switchWaterStyle = ()=>{
  let newWaterType= !watertype?"true":"false"
  //change search param "watertype" to newWaterType
  let URLu = new URL(window.location);
  URLu.searchParams.set("watertype", newWaterType);
  window.location.replace(URLu);
}
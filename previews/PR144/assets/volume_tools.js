/*
 * Volume estimates of Packmol.jl, for the browser (the "Volume guesser" and "Free volume"
 * pages of the documentation) and for Node.js (the tests of docs/test_volume_tools.js).
 * The same estimates as `guess_volume` and `free_volume` (src/tools/volume.jl).
 */
(function (root) {
  "use strict";

  const AVOGADRO = 6.02214076e23;
  // Å³ per (g/mol)/(g/mL): 1e24 / N_A
  const VOLUME_FACTOR = 1e24 / AVOGADRO;

  // Standard atomic masses (g/mol)
  const MASSES = {
    H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011, N: 14.007, O: 15.999,
    F: 18.998, Ne: 20.180, Na: 22.990, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974,
    S: 32.06, Cl: 35.45, Ar: 39.948, K: 39.098, Ca: 40.078, Ti: 47.867, Cr: 51.996,
    Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38, Ga: 69.723,
    Ge: 72.630, As: 74.922, Se: 78.971, Br: 79.904, Kr: 83.798, Rb: 85.468, Sr: 87.62,
    Ag: 107.87, Cd: 112.41, Sn: 118.71, I: 126.90, Xe: 131.29, Cs: 132.91, Ba: 137.33,
    Pt: 195.08, Au: 196.97, Hg: 200.59, Pb: 207.2,
  };

  // ---------------------------------------------------------------- volume guesser

  // The volume (Å³) of molecules of molar masses `masses` (g/mol), in the numbers
  // `numbers`, at the `density` (g/mL) of the whole system; the side of a cube and the
  // radius of a sphere of that volume.
  function guessVolume(masses, numbers, density) {
    if (!(density > 0)) throw new Error("The density must be a number greater than 0.");
    let mass = 0;
    for (let i = 0; i < masses.length; i++) {
      if (!(masses[i] > 0)) throw new Error("The molar mass of molecule type " + (i + 1) + " must be a number greater than 0.");
      if (!(numbers[i] >= 0) || !Number.isInteger(numbers[i])) throw new Error("The number of molecules of type " + (i + 1) + " must be an integer.");
      mass += masses[i] * numbers[i];
    }
    const volume = VOLUME_FACTOR * mass / density;
    return { mass, volume, side: Math.cbrt(volume), radius: Math.cbrt(3 * volume / (4 * Math.PI)) };
  }

  // The molar mass (g/mol) of a chemical formula (e.g. "C2H6O", "CH3(CH2)4CH3")
  function formulaMass(formula) {
    const s = formula.replace(/\s+/g, "");
    if (s === "") throw new Error("Empty formula.");
    let i = 0;
    function group() {
      let m = 0;
      while (i < s.length && s[i] !== ")") {
        let mi;
        if (s[i] === "(") {
          i++;
          mi = group();
          if (s[i] !== ")") throw new Error("Unbalanced parentheses in " + formula + ".");
          i++;
        } else {
          const e = /^[A-Z][a-z]?/.exec(s.slice(i));
          if (!e || !(e[0] in MASSES)) throw new Error("Unknown element at '" + s.slice(i) + "' in " + formula + ".");
          mi = MASSES[e[0]];
          i += e[0].length;
        }
        const n = /^\d+/.exec(s.slice(i));
        if (n) { mi *= parseInt(n[0]); i += n[0].length; }
        m += mi;
      }
      return m;
    }
    const m = group();
    if (i < s.length) throw new Error("Unbalanced parentheses in " + formula + ".");
    return m;
  }

  // ---------------------------------------------------------------- linear algebra

  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const norm2 = (a) => dot(a, a);
  // Matrix (rows) times vector
  const mulv = (M, v) => [dot(M[0], v), dot(M[1], v), dot(M[2], v)];
  function det3(M) {
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
      M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
      M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  }
  function inv3(M) {
    const d = det3(M);
    const c = (r1, c1, r2, c2) => M[r1][c1] * M[r2][c2] - M[r1][c2] * M[r2][c1];
    return [
      [c(1, 1, 2, 2) / d, -c(0, 1, 2, 2) / d, c(0, 1, 1, 2) / d],
      [-c(1, 0, 2, 2) / d, c(0, 0, 2, 2) / d, -c(0, 0, 1, 2) / d],
      [c(1, 0, 2, 1) / d, -c(0, 0, 2, 1) / d, c(0, 0, 1, 1) / d],
    ];
  }
  // The unit cell (columns: the cell vectors, as rows of the returned matrix of the
  // Cartesian components) from a b c α β γ (degrees), as Packmol's `unitcell_matrix`
  function unitcellMatrix(a, b, c, alpha, beta, gamma) {
    const r = Math.PI / 180;
    const [al, be, ga] = [alpha * r, beta * r, gamma * r];
    const bx = b * Math.cos(ga), by = b * Math.sin(ga);
    const cx = c * Math.cos(be);
    const cy = c * (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
    const cz = Math.sqrt(c * c - cx * cx - cy * cy);
    // H: rows are x, y, z; columns are the vectors a, b, c
    return [[a, bx, cx], [0, by, cy], [0, 0, cz]];
  }

  // Euler rotation matrix of Packmol (`eulermat`): Rx(β) Ry(γ) Rz(θ)
  function eulerMatrix(beta, gamma, theta) {
    const [c1, s1, c2, s2, c3, s3] = [Math.cos(beta), Math.sin(beta), Math.cos(gamma), Math.sin(gamma), Math.cos(theta), Math.sin(theta)];
    return [
      [c2 * c3, -c2 * s3, s2],
      [c1 * s3 + c3 * s1 * s2, c1 * c3 - s1 * s2 * s3, -c2 * s1],
      [s1 * s3 - c1 * c3 * s2, c1 * s2 * s3 + c3 * s1, c1 * c2],
    ];
  }

  // ---------------------------------------------------------------- structure files

  // The atoms of a PDB (first model) or XYZ file: positions and elements
  function readStructure(name, text) {
    const atoms = [];
    if (/\.xyz$/i.test(name)) {
      const lines = text.split(/\r?\n/);
      const n = parseInt(lines[0]);
      for (let i = 2; i < 2 + n && i < lines.length; i++) {
        const t = lines[i].trim().split(/\s+/);
        if (t.length < 4) continue;
        atoms.push({ x: [parseFloat(t[1]), parseFloat(t[2]), parseFloat(t[3])], element: elementOf(t[0]) });
      }
      return atoms;
    }
    for (const line of text.split(/\r?\n/)) {
      if (line.startsWith("ENDMDL")) break;
      if (!(line.startsWith("ATOM") || line.startsWith("HETATM"))) continue;
      const x = [parseFloat(line.slice(30, 38)), parseFloat(line.slice(38, 46)), parseFloat(line.slice(46, 54))];
      const el = line.length >= 78 ? line.slice(76, 78).trim() : "";
      atoms.push({ x, element: el !== "" ? elementOf(el) : elementOfName(line.slice(12, 16)) });
    }
    return atoms;
  }
  // The element of the element column of a PDB file (or of an XYZ file)
  function elementOf(s) {
    const letters = s.replace(/[^A-Za-z]/g, "");
    if (letters === "") return "";
    const two = letters[0].toUpperCase() + (letters[1] || "").toLowerCase();
    return letters.length === 2 && two in MASSES ? two : letters[0].toUpperCase();
  }
  // The element of an atom name, without an element column: the first letter, but for
  // the names of common ions (CA is the alpha carbon of proteins, not calcium)
  const IONS = { CL: "Cl", CLA: "Cl", NA: "Na", SOD: "Na", MG: "Mg", ZN: "Zn", BR: "Br", LI: "Li", POT: "K", CAL: "Ca" };
  function elementOfName(name) {
    const n = name.replace(/[^A-Za-z]/g, "").toUpperCase();
    if (n in IONS) return IONS[n];
    return n === "" ? "" : n[0];
  }
  function molarMass(atoms) {
    let m = 0;
    for (const a of atoms) {
      if (!(a.element in MASSES)) return NaN;
      m += MASSES[a.element];
    }
    return m;
  }

  // ---------------------------------------------------------------- input files

  const SHAPES = { box: 6, cube: 4, sphere: 4, cylinder: 8, ellipsoid: 7, plane: 4 };

  // Parses a Packmol input file: the tolerance, the periodic cell, and, for each
  // structure block, its file, number, fixed position, and constraints (those of the
  // whole molecule, and whether there are constraints of some atoms)
  function parseInput(text) {
    const lines = text.split(/\r?\n/).map((l) => l.replace(/#.*$/, "").trim());
    const input = { tolerance: 2.0, cell: null, structures: [] };
    let st = null;
    let inAtoms = false;
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      if (line === "") continue;
      const t = line.split(/\s+/);
      const key = t[0].toLowerCase();
      if (st === null) {
        if (key === "tolerance") input.tolerance = parseFloat(t[1]);
        else if (key === "pbc") input.cell = parsePbc(t.slice(1));
        else if (key === "unitcell") {
          const p = t.slice(1, 7).map(parseFloat);
          input.cell = { H: unitcellMatrix(...p), center: [0, 0, 0], wignerSeitz: false };
        } else if (key === "structure") {
          st = { file: t.slice(1).join(" "), number: 1, fixed: null, center: "none", constraints: [], someAtoms: false, unsupported: [] };
        }
        continue;
      }
      if (key === "end" && t[1] && t[1].toLowerCase() === "structure") {
        input.structures.push(st);
        st = null;
        inAtoms = false;
        continue;
      }
      if (inAtoms) {
        if (key === "end" && t[1] && t[1].toLowerCase() === "atoms") inAtoms = false;
        else if (isConstraint(key)) st.someAtoms = true;
        continue;
      }
      if (key === "constrain_rotation") {
        // A block, when followed by the lines `atoms ...` and `to ...`
        if (i + 2 < lines.length && /^atoms\b/i.test(lines[i + 1]) && /^to\b/i.test(lines[i + 2])) {
          while (i < lines.length && !/^end\s+constrain_rotation/i.test(lines[i])) i++;
        }
        continue;
      }
      if (key === "atoms") { inAtoms = true; continue; }
      if (key === "number") st.number = parseInt(t[1]);
      else if (key === "fixed") st.fixed = t.slice(1, 7).map(parseFloat);
      else if (key === "center") st.center = "geometric";
      else if (key === "centerofmass") st.center = "mass";
      else if (isConstraint(key)) {
        const shape = (t[1] || "").toLowerCase();
        if (!(shape in SHAPES)) { st.unsupported.push(line); continue; }
        const p = t.slice(2, 2 + SHAPES[shape]).map(parseFloat);
        if (p.length < SHAPES[shape] || p.some(isNaN)) throw new Error("Could not read the constraint '" + line + "'.");
        st.constraints.push(makeConstraint(key === "over" ? "above" : key, shape, p));
      }
    }
    if (st !== null) throw new Error("A structure block (" + st.file + ") is not closed with 'end structure'.");
    return input;
  }
  function isConstraint(key) {
    return key === "inside" || key === "outside" || key === "above" || key === "over" || key === "below";
  }
  function parsePbc(t) {
    const kind = (t[0] || "").toLowerCase();
    if (kind === "dodecahedral" || kind === "octahedral") {
      const [cx, cy, cz, d] = t.slice(1, 5).map(parseFloat);
      const H = kind === "dodecahedral" ? unitcellMatrix(d, d, d, 60, 60, 90) :
        unitcellMatrix(d, d, d, Math.acos(1 / 3) * 180 / Math.PI, Math.acos(-1 / 3) * 180 / Math.PI, Math.acos(1 / 3) * 180 / Math.PI);
      return { H, center: [cx, cy, cz], wignerSeitz: true };
    }
    const v = t.map(parseFloat);
    if (v.length === 3) return { H: [[v[0], 0, 0], [0, v[1], 0], [0, 0, v[2]]], center: scale(v, 0.5), wignerSeitz: false };
    if (v.length === 6) {
      const lo = v.slice(0, 3), hi = v.slice(3, 6);
      const s = sub(hi, lo);
      return { H: [[s[0], 0, 0], [0, s[1], 0], [0, 0, s[2]]], center: scale(add(lo, hi), 0.5), wignerSeitz: false };
    }
    throw new Error("Could not read the pbc line: pbc " + t.join(" "));
  }

  // A constraint: `test(x)` is true when the point x satisfies it; `box`, for the
  // inside constraints, the bounding box of their region
  function makeConstraint(kind, shape, p) {
    let inside, box = null;
    if (shape === "box") {
      const lo = p.slice(0, 3), hi = p.slice(3, 6);
      inside = (x) => x[0] >= lo[0] && x[0] <= hi[0] && x[1] >= lo[1] && x[1] <= hi[1] && x[2] >= lo[2] && x[2] <= hi[2];
      box = [lo, hi];
    } else if (shape === "cube") {
      const lo = p.slice(0, 3), hi = add(lo, [p[3], p[3], p[3]]);
      inside = (x) => x[0] >= lo[0] && x[0] <= hi[0] && x[1] >= lo[1] && x[1] <= hi[1] && x[2] >= lo[2] && x[2] <= hi[2];
      box = [lo, hi];
    } else if (shape === "sphere") {
      const c = p.slice(0, 3), r = p[3];
      inside = (x) => norm2(sub(x, c)) <= r * r;
      box = [sub(c, [r, r, r]), add(c, [r, r, r])];
    } else if (shape === "ellipsoid") {
      const c = p.slice(0, 3), [a, b, cc, s] = p.slice(3, 7);
      inside = (x) => {
        const d = sub(x, c);
        return (d[0] / a) ** 2 + (d[1] / b) ** 2 + (d[2] / cc) ** 2 - s * s <= 0;
      };
      const h = [Math.abs(s * a), Math.abs(s * b), Math.abs(s * cc)];
      box = [sub(c, h), add(c, h)];
    } else if (shape === "cylinder") {
      const c = p.slice(0, 3), axis = p.slice(3, 6), r = p[6], len = p[7];
      const u = scale(axis, 1 / Math.sqrt(norm2(axis)));
      inside = (x) => {
        const d = sub(x, c);
        const w = dot(d, u);
        return w >= 0 && w <= len && norm2(sub(d, scale(u, w))) <= r * r;
      };
      const tip = add(c, scale(u, len));
      box = [sub([Math.min(c[0], tip[0]), Math.min(c[1], tip[1]), Math.min(c[2], tip[2])], [r, r, r]),
        add([Math.max(c[0], tip[0]), Math.max(c[1], tip[1]), Math.max(c[2], tip[2])], [r, r, r])];
    } else { // plane
      const n = p.slice(0, 3), d = p[3];
      if (kind === "above") return { test: (x) => dot(n, x) - d >= 0, box: null };
      if (kind === "below") return { test: (x) => dot(n, x) - d <= 0, box: null };
      throw new Error("A plane takes 'above', 'over' or 'below', not '" + kind + "'.");
    }
    if (kind === "inside") return { test: inside, box };
    if (kind === "outside") return { test: (x) => !inside(x), box: null };
    throw new Error("A " + shape + " takes 'inside' or 'outside', not '" + kind + "'.");
  }

  // ---------------------------------------------------------------- free volume

  // A pseudo-random generator (mulberry32), for reproducible estimates
  function rng(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // The positions of the atoms of a fixed molecule, as Packmol places them
  function fixedPositions(st, atoms) {
    const [x, y, z, a, b, c] = st.fixed;
    const R = eulerMatrix(a, b, c);
    let origin = [0, 0, 0];
    if (st.center !== "none") {
      let w = 0;
      for (const at of atoms) {
        const m = st.center === "mass" ? (MASSES[at.element] || 0) : 1;
        origin = add(origin, scale(at.x, m));
        w += m;
      }
      origin = scale(origin, 1 / w);
    }
    return atoms.map((at) => add(mulv(R, sub(at.x, origin)), [x, y, z]));
  }

  // A spatial hash of points, for the test of proximity to the fixed atoms
  function spatialHash(points, h) {
    const cells = new Map();
    const key = (i, j, k) => (i * 73856093) ^ (j * 19349663) ^ (k * 83492791);
    for (const p of points) {
      const k = key(Math.floor(p[0] / h), Math.floor(p[1] / h), Math.floor(p[2] / h));
      let list = cells.get(k);
      if (list === undefined) { list = []; cells.set(k, list); }
      list.push(p);
    }
    return function near(x, tol) {
      const i0 = Math.floor(x[0] / h), j0 = Math.floor(x[1] / h), k0 = Math.floor(x[2] / h);
      const t2 = tol * tol;
      for (let i = i0 - 1; i <= i0 + 1; i++) for (let j = j0 - 1; j <= j0 + 1; j++) for (let k = k0 - 1; k <= k0 + 1; k++) {
        const list = cells.get(key(i, j, k));
        if (list === undefined) continue;
        for (const p of list) if (norm2(sub(x, p)) < t2) return true;
      }
      return false;
    };
  }

  // The image of x nearest the center of the cell (for the Wigner-Seitz cells), or the
  // image in the parallelepiped centered at the center (the other cells)
  function wrapper(cell) {
    const Hinv = inv3(cell.H);
    const col = (j) => [cell.H[0][j], cell.H[1][j], cell.H[2][j]];
    const vectors = [];
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) for (let k = -1; k <= 1; k++) {
      if (i === 0 && j === 0 && k === 0) continue;
      vectors.push(add(add(scale(col(0), i), scale(col(1), j)), scale(col(2), k)));
    }
    return function (x) {
      const d = sub(x, cell.center);
      const f = mulv(Hinv, d).map((v) => v - Math.round(v));
      let y = mulv(cell.H, f);
      if (cell.wignerSeitz) {
        // the nearest image: the lattice vector that most shortens y, repeated
        for (let iter = 0; iter < 4; iter++) {
          let best = null, bestn = norm2(y);
          for (const v of vectors) {
            const n = norm2(add(y, v));
            if (n < bestn - 1e-12) { bestn = n; best = v; }
          }
          if (best === null) break;
          y = add(y, best);
        }
      }
      return add(y, cell.center);
    };
  }

  // Estimates the volume of the region of each free structure type of an input:
  // `files` maps the names of the structure files to their texts (needed for the fixed
  // molecules, and for densities in g/mL), and `progress(fraction)` is called during the
  // sampling. Returns, for each free structure type, { file, number, volume, error,
  // density (g/mL, NaN without the file), numberDensity (molecules/nm³), note }, and
  // the warnings of the estimate.
  function freeVolume(input, files, npoints, seed) {
    const random = rng(seed === undefined ? 1 : seed);
    const warnings = [];
    const tol = input.tolerance;
    // The fixed atoms, with their periodic images near the cell
    let fixed = [];
    for (const st of input.structures) {
      if (st.fixed === null) continue;
      if (!(st.file in files)) {
        warnings.push("The structure file of the fixed molecule " + st.file + " was not given: it was ignored.");
        continue;
      }
      fixed = fixed.concat(fixedPositions(st, readStructure(st.file, files[st.file])));
    }
    const cell = input.cell;
    let wrap = null, cellVolume = 0;
    if (cell !== null) {
      wrap = wrapper(cell);
      cellVolume = Math.abs(det3(cell.H));
      if (fixed.length > 0) {
        // Images within reach of the points of the cell: within its circumradius plus
        // the tolerance of its center
        let rmax = 0;
        for (const i of [-0.5, 0.5]) for (const j of [-0.5, 0.5]) for (const k of [-0.5, 0.5]) {
          rmax = Math.max(rmax, Math.sqrt(norm2(mulv(cell.H, [i, j, k]))));
        }
        const images = [];
        for (const p of fixed) {
          for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) for (let k = -1; k <= 1; k++) {
            const q = add(wrap(p), mulv(cell.H, [i, j, k]));
            if (Math.sqrt(norm2(sub(q, cell.center))) <= rmax + tol) images.push(q);
          }
        }
        fixed = images;
      }
    }
    const near = fixed.length > 0 ? spatialHash(fixed, tol) : null;
    const results = [];
    for (const st of input.structures) {
      if (st.fixed !== null) continue;
      const r = { file: st.file, number: st.number, volume: NaN, error: NaN, density: NaN, numberDensity: NaN, note: "" };
      results.push(r);
      if (st.unsupported.length > 0) { r.note = "constraint not supported here: " + st.unsupported[0]; continue; }
      if (st.someAtoms && st.constraints.length === 0) { r.note = "region set by constraints of some atoms only"; continue; }
      // The region to sample: the periodic cell, or the bounding box of the inside constraints
      let lo = [-Infinity, -Infinity, -Infinity], hi = [Infinity, Infinity, Infinity];
      for (const c of st.constraints) {
        if (c.box === null) continue;
        lo = lo.map((v, d) => Math.max(v, c.box[0][d]));
        hi = hi.map((v, d) => Math.min(v, c.box[1][d]));
      }
      let sample, boxVolume;
      if (cell !== null) {
        boxVolume = cellVolume;
        sample = () => wrap(add(cell.center, mulv(cell.H, [random() - 0.5, random() - 0.5, random() - 0.5])));
      } else {
        if (!lo.every(isFinite) || !hi.every(isFinite)) { r.note = "region not bounded (no periodic cell or inside constraint)"; continue; }
        if (hi.some((v, d) => v <= lo[d])) { r.volume = 0; r.error = 0; r.note = "the constraints do not intersect"; continue; }
        const ext = sub(hi, lo);
        boxVolume = ext[0] * ext[1] * ext[2];
        sample = () => [lo[0] + random() * ext[0], lo[1] + random() * ext[1], lo[2] + random() * ext[2]];
      }
      let count = 0;
      for (let n = 0; n < npoints; n++) {
        const x = sample();
        let ok = true;
        for (const c of st.constraints) if (!c.test(x)) { ok = false; break; }
        if (ok && near !== null && near(x, tol)) ok = false;
        if (ok) count++;
      }
      const p = count / npoints;
      r.volume = p * boxVolume;
      r.error = boxVolume * Math.sqrt(p * (1 - p) / npoints);
      r.numberDensity = st.number / r.volume * 1000;
      if (st.file in files) {
        const m = molarMass(readStructure(st.file, files[st.file]));
        r.density = st.number * m / AVOGADRO / (r.volume * 1e-24);
      }
      if (st.someAtoms) r.note = "upper bound: constraints of some atoms ignored";
    }
    return { results, warnings };
  }

  const api = { guessVolume, formulaMass, parseInput, freeVolume, readStructure, molarMass, MASSES };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.PackmolVolume = api;
})(typeof window !== "undefined" ? window : this);

// Tests of docs/src/assets/volume_tools.js against the Julia functions:
//   node docs/test_volume_tools.js <results.json>
// <results.json>: the Julia estimates (written by docs/test_volume_tools.jl), with the
// inputs and their directories.
const fs = require("fs");
const path = require("path");
const V = require("./src/assets/volume_tools.js");
const reference = JSON.parse(fs.readFileSync(process.argv[2], "utf8"), (k, v) => (v === null ? NaN : v));
let failed = 0;
function check(name, ok, detail) {
  console.log((ok ? "  ok   " : "  FAIL ") + name + (detail ? ": " + detail : ""));
  if (!ok) failed++;
}
// guess_volume
const g = V.guessVolume([18.015], [1000], 1.0);
check("guessVolume of 1000 waters", Math.abs(g.volume - 18015 / 6.02214076e23 * 1e24) < 1e-6);
check("formulaMass(H2O)", Math.abs(V.formulaMass("H2O") - 18.015) < 1e-3);
check("formulaMass(CH3(CH2)4CH3)", Math.abs(V.formulaMass("CH3(CH2)4CH3") - 86.178) < 1e-2);
for (const ref of reference) {
  const text = fs.readFileSync(ref.input, "utf8");
  const input = V.parseInput(text);
  const files = {};
  for (const st of input.structures) {
    const f = path.join(path.dirname(ref.input), st.file);
    if (fs.existsSync(f)) files[st.file] = fs.readFileSync(f, "utf8");
  }
  const { results } = V.freeVolume(input, files, ref.npoints, 1);
  check(path.basename(ref.input) + ": structure types", results.length === ref.results.length);
  results.forEach((r, i) => {
    const j = ref.results[i];
    if (isNaN(j.volume)) {
      check(path.basename(ref.input) + " " + r.file + " (no estimate)", isNaN(r.volume), r.note);
    } else {
      const err = Math.sqrt(r.error ** 2 + j.error ** 2);
      const ok = Math.abs(r.volume - j.volume) <= 4 * err + 1e-6 * j.volume;
      check(path.basename(ref.input) + " " + r.file, ok,
        "volume " + r.volume.toFixed(0) + " (Julia " + j.volume.toFixed(0) + " ± " + err.toFixed(0) +
        "), density " + r.density.toFixed(3) + " (Julia " + j.density.toFixed(3) + ")");
      check(path.basename(ref.input) + " " + r.file + " density", Math.abs(r.density - j.density) <= 4 * err / j.volume * j.density + 2e-3);
    }
  });
}
console.log(failed === 0 ? "all tests passed" : failed + " tests failed");
process.exit(failed === 0 ? 0 : 1);

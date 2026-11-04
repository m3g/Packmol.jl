# Free volume

The free volume calculator reads a Packmol input file and estimates, by Monte Carlo, the
volume of the region available to the molecules of each structure type: the fraction of
random points of the periodic cell (or of the box of its `inside` constraints) that
satisfy the constraints of the structure type and are farther than the tolerance from the
atoms of the fixed molecules. With it comes the density at which the molecules of the
structure type fill their region, which shows when an input asks for too many molecules
(the packing is then slow, or fails) or too few. The densities of structure types that
share a region add up. Nothing leaves the page: the files are read by the browser.

Paste the input file, or open it, and select its structure files: the files of the fixed
molecules are needed to take their volume out of the regions, and those of the free
molecules give the densities in g/mL (without them, the densities are given in molecules
per nm³).

```@raw html
<div class="pm-app" id="pm-free">
  <div class="pm-row">
    <label>Input file</label>
    <input id="pm-input-file" type="file" accept=".inp,.txt,text/plain">
  </div>
  <textarea id="pm-input" rows="12" spellcheck="false" placeholder="tolerance 2.0
pbc 40. 40. 40.

structure water.pdb
  number 1000
  inside sphere 20. 20. 20. 15.
end structure"></textarea>
  <div class="pm-row">
    <label>Structure files</label>
    <input id="pm-structures" type="file" multiple accept=".pdb,.xyz">
    <span id="pm-loaded" class="pm-note"></span>
  </div>
  <div class="pm-row">
    <label for="pm-npoints">Points</label>
    <select id="pm-npoints">
      <option value="100000">100 thousand</option>
      <option value="300000" selected>300 thousand</option>
      <option value="1000000">1 million</option>
    </select>
    <button class="pm-button" id="pm-run" type="button">Estimate the volumes</button>
  </div>
  <div class="pm-result" id="pm-free-result"></div>
</div>
<script>
(function () {
  const V = window.PackmolVolume;
  const files = {};
  const textarea = document.getElementById("pm-input");
  const out = document.getElementById("pm-free-result");
  const read = (f) => new Promise((ok, fail) => { const r = new FileReader(); r.onload = () => ok(r.result); r.onerror = fail; r.readAsText(f); });
  document.getElementById("pm-input-file").addEventListener("change", async (e) => {
    if (e.target.files.length > 0) textarea.value = await read(e.target.files[0]);
  });
  document.getElementById("pm-structures").addEventListener("change", async (e) => {
    for (const f of e.target.files) files[f.name] = await read(f);
    document.getElementById("pm-loaded").textContent = Object.keys(files).length + " file(s): " + Object.keys(files).join(", ");
  });
  function fmt(x, digits) {
    if (!isFinite(x)) return "&ndash;";
    return Math.abs(x) >= 1e4 ? Math.round(x).toLocaleString("en-US") : x.toFixed(digits);
  }
  document.getElementById("pm-run").addEventListener("click", () => {
    out.innerHTML = '<div class="pm-note">Estimating&hellip;</div>';
    setTimeout(() => {
      try {
        const input = V.parseInput(textarea.value);
        if (input.structures.length === 0) throw new Error("The input has no structure blocks.");
        // the structure files are matched by their names, without the directories
        const byName = {};
        for (const st of input.structures) {
          const base = st.file.split(/[\\/]/).pop();
          if (base in files) byName[st.file] = files[base];
        }
        const npoints = parseInt(document.getElementById("pm-npoints").value);
        const { results, warnings } = V.freeVolume(input, byName, npoints, 1);
        let html = '<table class="pm-table"><thead><tr><th>Structure</th><th>Number</th><th>Volume (&#xC5;<sup>3</sup>)</th>' +
          '<th>&plusmn; error</th><th>Density (g/mL)</th><th>Molecules per nm<sup>3</sup></th><th></th></tr></thead><tbody>';
        for (const r of results) {
          html += "<tr><td>" + r.file + "</td><td>" + r.number + "</td><td>" + fmt(r.volume, 1) + "</td><td>" + fmt(r.error, 1) +
            "</td><td>" + fmt(r.density, 3) + "</td><td>" + fmt(r.numberDensity, 2) + '</td><td class="pm-note">' + r.note + "</td></tr>";
        }
        html += "</tbody></table>";
        for (const w of warnings) html += '<div class="pm-error">' + w + "</div>";
        html += '<div class="pm-note">' + npoints.toLocaleString("en-US") + " points per structure type. " +
          "Water fills its region at about 1 g/mL (33.4 molecules per nm<sup>3</sup>).</div>";
        out.innerHTML = html;
      } catch (err) {
        out.innerHTML = '<div class="pm-error">' + err.message + "</div>";
      }
    }, 20);
  });
})();
</script>
```

The region is that of the constraints that apply to the whole molecule. The constraints
of some atoms only (in `atoms` blocks) are ignored: the volume is then an upper bound, and
a structure type with constraints of some atoms only (as the lipids of a membrane, whose
heads and tails are above and below planes) has no estimate. Neither has one a structure
type whose region is not bounded, by the periodic cell or by an `inside` box, cube,
sphere, cylinder or ellipsoid. The constraints of Packmol.jl without an equivalent in
Fortran Packmol (Gaussian, wave and exponential surfaces) are not supported here.

The same estimate is given by [`free_volume`](@ref) in Julia, and by
`packmol --free-volume input.inp`, which read the structure files of the input (see
[Volume estimates](@ref)).

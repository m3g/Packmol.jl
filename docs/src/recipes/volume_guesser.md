# Volume guesser

Given the number of molecules of each type, their molar masses, and the density of the
system, the volume guesser estimates the volume of the system: each molecule occupies a
volume proportional to its molar mass, ``V = \sum_i N_i M_i / (\rho N_A)``. It gives the
side of a cubic box, or the radius of a sphere, of that volume. Unless the system is very
sensitive to the density (as a gas), this is a good starting point for the box of an input
file; a simulation at constant pressure then adjusts it. Nothing typed here leaves the
page.

```@raw html
<div class="pm-app" id="pm-guesser">
  <div class="pm-row">
    <label for="pm-density">Density of the system</label>
    <input id="pm-density" type="number" value="1.0" min="0" step="0.01"> g/mL
  </div>
  <table class="pm-table">
    <thead><tr><th>Molecule</th><th>Number of molecules</th><th>Molar mass (g/mol)</th><th>or formula</th><th></th></tr></thead>
    <tbody id="pm-types"></tbody>
  </table>
  <button class="pm-button" id="pm-add" type="button">+ Add a molecule type</button>
  <div class="pm-result" id="pm-result"></div>
</div>
<script>
(function () {
  const V = window.PackmolVolume;
  const tbody = document.getElementById("pm-types");
  const result = document.getElementById("pm-result");
  function addRow(name, number, mass) {
    const tr = document.createElement("tr");
    tr.innerHTML =
      '<td><input class="pm-name" type="text" placeholder="water" value="' + name + '"></td>' +
      '<td><input class="pm-number" type="number" min="0" step="1" value="' + number + '"></td>' +
      '<td><input class="pm-mass" type="number" min="0" step="0.001" value="' + mass + '"></td>' +
      '<td><input class="pm-formula" type="text" placeholder="H2O"></td>' +
      '<td><button class="pm-remove" type="button" title="Remove">&times;</button></td>';
    tbody.appendChild(tr);
    tr.querySelector(".pm-formula").addEventListener("input", (e) => {
      try {
        if (e.target.value.trim() !== "") tr.querySelector(".pm-mass").value = V.formulaMass(e.target.value).toFixed(3);
        e.target.classList.remove("pm-invalid");
      } catch (err) {
        e.target.classList.add("pm-invalid");
      }
      update();
    });
    tr.querySelector(".pm-remove").addEventListener("click", () => { tr.remove(); update(); });
    tr.querySelectorAll("input").forEach((i) => i.addEventListener("input", update));
  }
  function update() {
    const rows = Array.from(tbody.querySelectorAll("tr"));
    const density = parseFloat(document.getElementById("pm-density").value);
    const masses = rows.map((r) => parseFloat(r.querySelector(".pm-mass").value));
    const numbers = rows.map((r) => Number(r.querySelector(".pm-number").value));
    if (rows.length === 0) { result.innerHTML = ""; return; }
    try {
      const g = V.guessVolume(masses, numbers, density);
      const names = rows.map((r, i) => r.querySelector(".pm-name").value.trim() || ("molecule " + (i + 1)));
      const julia = "guess_volume([" + rows.map((r, i) => masses[i] + " => " + numbers[i]).join(", ") + "]; density = " + density + ")";
      result.innerHTML =
        '<div><span>The system must have a volume of</span> <b>' + g.volume.toFixed(2) + ' &#xC5;<sup>3</sup></b></div>' +
        '<div><span>a cubic box of side</span> <b>' + g.side.toFixed(2) + ' &#xC5;</b></div>' +
        '<div><span>or a sphere of radius</span> <b>' + g.radius.toFixed(2) + ' &#xC5;</b></div>' +
        '<div class="pm-note">Total mass: ' + g.mass.toFixed(2) + ' g/mol (' + names.map((n, i) => numbers[i] + " " + n).join(", ") + '). In Julia: <code>' + julia + '</code></div>';
    } catch (err) {
      result.innerHTML = '<div class="pm-error">' + err.message + '</div>';
    }
  }
  document.getElementById("pm-add").addEventListener("click", () => { addRow("", 100, ""); update(); });
  document.getElementById("pm-density").addEventListener("input", update);
  addRow("water", 1000, "18.015");
  update();
})();
</script>
```

The same estimate is given by [`guess_volume`](@ref) in Julia, which also reads the
molecules and their numbers from an input file or from structure files, and by
`packmol --guess-volume input.inp` (see [Volume estimates](@ref)). The volume of the
region available to the molecules of each structure type of an input, given its
constraints and fixed molecules, is estimated by the [free volume](@ref "Free volume")
calculator.

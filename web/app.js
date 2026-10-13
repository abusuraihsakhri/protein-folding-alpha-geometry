"use strict";

const PYODIDE_ROOT = "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/";
const EXAMPLE = [
  {name:"ALA",index:1,phi:-57,psi:-47},
  {name:"ALA",index:2,phi:-62,psi:-42},
  {name:"ALA",index:3,phi:-54,psi:-50},
  {name:"VAL",index:4,phi:-135,psi:135},
  {name:"VAL",index:5,phi:-132,psi:139},
  {name:"VAL",index:6,phi:-138,psi:131}
];
const el = id => document.getElementById(id);
let runtimePromise = null;
let latestJSON = null;

function setStatus(message, error = false) {
  el("status").textContent = message;
  el("status").classList.toggle("error", error);
}
function cutoff(id) {
  const n = Number(el(id).value);
  if (!el(id).value.trim() || !Number.isFinite(n) || n <= 0 || n > 100) {
    throw new Error("Cutoffs must be finite numbers greater than 0 and at most 100 Å.");
  }
  return n;
}
async function initializePython() {
  if (!runtimePromise) {
    runtimePromise = (async () => {
      if (typeof loadPyodide !== "function") {
        throw new Error("Python runtime unavailable. Check network access to the Pyodide CDN.");
      }
      setStatus("Loading Python runtime and geometry engine…");
      const py = await loadPyodide({indexURL: PYODIDE_ROOT});
      py.FS.mkdirTree("/home/pyodide/evofold_geometry");
      for (const name of ["__init__.py", "engine.py", "cli.py"]) {
        const url = "./evofold_geometry/" + name;
        const response = await fetch(url);
        if (!response.ok) throw new Error("Unable to load geometry file: " + url);
        py.FS.writeFile("/home/pyodide/evofold_geometry/" + name, await response.text());
      }
      py.runPython([
        "import sys, json",
        "from pathlib import Path",
        "from dataclasses import asdict",
        "sys.path.insert(0, '/home/pyodide')",
        "from evofold_geometry.cli import _parse_residues_from_json",
        "from evofold_geometry.engine import analyze_structure",
        "def run_geometry(raw_json, hbond_cutoff, contact_cutoff):",
        "    Path('/tmp/residues_input.json').write_text(raw_json, encoding='utf-8')",
        "    residues = _parse_residues_from_json('/tmp/residues_input.json')",
        "    if len(residues) > 1000:",
        "        raise ValueError('Browser analysis is limited to 1,000 residues')",
        "    result = analyze_structure(residues, float(hbond_cutoff), float(contact_cutoff))",
        "    return json.dumps(asdict(result), allow_nan=False)"
      ].join("\n"));
      return py;
    })().catch(error => { runtimePromise = null; throw error; });
  }
  return runtimePromise;
}
function drawPlot(points) {
  const canvas = el("rama-plot"), ctx = canvas.getContext("2d");
  const w = canvas.width, h = canvas.height;
  const l = 55, t = 20, innerW = w - 75, innerH = h - 66;
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = "#101a2a"; ctx.fillRect(0, 0, w, h);
  ctx.strokeStyle = "#425572"; ctx.lineWidth = 1;
  ctx.font = "13px system-ui"; ctx.fillStyle = "#baccdf";
  for (const v of [-180,-90,0,90,180]) {
    const x = l + (v+180)*innerW/360, y = t + (180-v)*innerH/360;
    ctx.beginPath(); ctx.moveTo(x,t); ctx.lineTo(x,t+innerH); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(l,y); ctx.lineTo(l+innerW,y); ctx.stroke();
    ctx.fillText(String(v), x-16, t+innerH+20);
    ctx.fillText(String(v), 5, y+4);
  }
  ctx.fillText("phi (degrees)", l+innerW/2-45, h-9);
  ctx.fillText("psi", 8, 12);
  for (const p of points) {
    const phi = ((p.phi+180)%360+360)%360-180;
    const psi = ((p.psi+180)%360+360)%360-180;
    const x = l + (phi+180)*innerW/360;
    const y = t + (180-psi)*innerH/360;
    ctx.beginPath(); ctx.arc(x,y,4,0,Math.PI*2);
    ctx.fillStyle = p.region==="favored" ? "#49cfb4" : p.region==="allowed" ? "#ffbe69" : "#fa7a83";
    ctx.fill();
  }
}
function showResult(result) {
  const rama = result.summary.ramachandran;
  el("num-residues").textContent = String(result.num_residues);
  el("favored").textContent = rama.total ? rama.favored_pct.toFixed(1) : "N/A";
  el("hbonds").textContent = String(result.hbonds.length);
  el("contacts").textContent = String(result.contact_map ? result.contact_map.contacts.length : 0);
  latestJSON = JSON.stringify(result,null,2);
  el("output").textContent = latestJSON;
  el("download").disabled = false;
  drawPlot(result.ramachandran);
}
async function analyze() {
  el("analyze").disabled = true;
  el("download").disabled = true; latestJSON = null;
  try {
    const input = el("residues").value;
    if (new Blob([input]).size > 1024*1024) throw new Error("JSON input exceeds 1 MB.");
    const hb = cutoff("hbond"), co = cutoff("contact");
    const py = await initializePython();
    setStatus("Calculating geometry…");
    py.globals.set("raw_json", input);
    py.globals.set("hbond_cutoff", hb);
    py.globals.set("contact_cutoff", co);
    showResult(JSON.parse(py.runPython("run_geometry(raw_json, hbond_cutoff, contact_cutoff)")));
    setStatus("Analysis complete. Results computed locally in this browser.");
  } catch (error) {
    el("output").textContent = "Analysis failed: " + String(error.message || error);
    setStatus(String(error.message || error), true);
  } finally { el("analyze").disabled = false; }
}
el("sample").addEventListener("click", () => {
  el("residues").value = JSON.stringify(EXAMPLE,null,2); setStatus("Example loaded.");
});
el("upload").addEventListener("change", async event => {
  const file = event.target.files && event.target.files[0];
  if (!file) return;
  if (file.size > 1024*1024) return setStatus("Selected JSON file exceeds 1 MB.",true);
  try { el("residues").value = await file.text(); setStatus("File loaded; run analysis."); }
  catch (error) { setStatus("File read failed: " + String(error.message || error), true); }
});
el("analyze").addEventListener("click",analyze);
el("download").addEventListener("click", () => {
  if (latestJSON === null) return;
  const url = URL.createObjectURL(new Blob([latestJSON],{type:"application/json"}));
  const link = document.createElement("a");
  link.href = url; link.download = "protein-geometry-analysis.json";
  document.body.appendChild(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url),1000);
});
el("residues").value = JSON.stringify(EXAMPLE,null,2);
drawPlot([]);

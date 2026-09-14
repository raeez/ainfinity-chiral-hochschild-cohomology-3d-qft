"""Record the exact source closure, compiler inputs, and local proof artifact."""

from pathlib import Path
import difflib
import hashlib
import json
import re
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[4]
REPORT = Path(__file__).resolve().parent
SOURCE = ROOT / "research-candidates/vol2_boundary026/collision005"
BUILD = REPORT / "build-pdftex"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size}


def write_json(name, value):
    (REPORT / name).write_text(json.dumps(value, indent=2) + "\n")


files = sorted(p for p in SOURCE.rglob("*") if p.is_file())
rows = [record(p) for p in files]
aggregate = hashlib.sha256("".join(f"{x['sha256']}  {Path(x['path']).relative_to(SOURCE)}\n" for x in rows).encode()).hexdigest()
with tarfile.open(REPORT / "source-freeze.tar.gz", "w:gz", dereference=True) as archive:
    for path in files:
        archive.add(path, arcname=str(path.relative_to(SOURCE)))

compiler_paths = set()
for line in (BUILD / "reader.fls").read_text().splitlines():
    if line.startswith("INPUT "):
        path = Path(line[6:])
        if not path.is_absolute():
            path = SOURCE / path
        path = path.resolve()
        if path.is_file() and BUILD not in path.parents:
            compiler_paths.add(path)
compiler_rows = [record(p) for p in sorted(compiler_paths)]
write_json("compiler-inputs.json", compiler_rows)
with tarfile.open(REPORT / "compiler-inputs.tar.gz", "w:gz") as archive:
    for index, path in enumerate(sorted(compiler_paths)):
        archive.add(path, arcname=f"{index:04d}-{path.name}")

new_source = SOURCE / "chapters/collision-chains.tex"
patch = "".join(difflib.unified_diff([], new_source.read_text().splitlines(True), fromfile="/dev/null", tofile="b/chapters/collision-chains.tex"))
(REPORT / "native-delta.patch").write_text(patch)

log = (BUILD / "reader.log").read_text(errors="replace")
issues = [line for line in log.splitlines() if any(word in line for word in ["Overfull", "Underfull", "undefined", "multiply defined", "Missing character", "Too many math"])]
aux = (BUILD / "reader.aux").read_text()
labels = [line for line in aux.splitlines() if line.startswith(r"\newlabel{eq:cc-") and "@cref" not in line]
bad_equations = [line for line in labels if not re.search(r"\{(?:equation|AMS)\.", line)]
info = subprocess.check_output(["pdfinfo", str(BUILD / "reader.pdf")], text=True)
pages = int(re.search(r"^Pages:\s*(\d+)", info, re.MULTILINE).group(1))
text = (REPORT / "reader.txt").read_text()
firewall_pattern = re.compile(r"\b(?:agent|worktree|task|candidate|reviewer|audit|acceptance|commit|workflow)\b", re.I)
new_hits = [{"line": n, "text": line} for n, line in enumerate(new_source.read_text().splitlines(), 1) if firewall_pattern.search(line)]
path_pattern = re.compile(r"(?:/Users/[^/ ]+/|~/)(?:ecosystem|centcom|kernel|moxie|momentum)(?:/|\b)")
path_hits = path_pattern.findall(text)
write_json("verification-final.json", {
    "build_diagnostics": issues, "new_equation_labels": labels,
    "incorrect_equation_destinations": bad_equations,
    "source_firewall_hits": new_hits, "reader_closed_path_hits": path_hits,
    "reader_gate": "FAIL: shared-template equation numbers and destinations",
    "visual_inspection": "New proof pages 54-59 and bibliography page 60 inspected as PNGs. No clipping, overlap, or manuscript-firewall violation identified. Equation tags are absent.",
    "calculation": json.loads((REPORT / "calculation.json").read_text()),
})
manifest = {
    "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
    "worktree": str(ROOT), "source_entrypoint": str(SOURCE / "reader.tex"),
    "files": rows, "source_aggregate_sha256": aggregate,
    "new_mathematical_body": record(new_source),
    "aggregate_diff_sha256": sha(REPORT / "native-delta.patch"),
    "compiler_inputs_manifest": record(REPORT / "compiler-inputs.json"),
    "pdf": {**record(BUILD / "reader.pdf"), "pages": pages},
    "build_command": str(REPORT / "build-command-final.json"),
    "exact_check_command": "/opt/homebrew/bin/python3 reports/research/vol2_boundary026/collision005/check_residue.py",
    "mathematical_status": "Constructed proof; fresh independent review required",
    "reader_status": "Working render; shared-template equation-reference gate fails",
    "scope": "Elementary three-current residue complex with all elementary jets; maximal unchanged test subcomplex for all three polynomial insertions; obstruction to a regular-product ternary primitive",
    "residual": "Composite collision-stratum differential and regular-coefficient comparison remain unconstructed. No nonzero higher map into bare B is claimed.",
}
write_json("manifest.json", manifest)
print(json.dumps({"manifest_sha256": sha(REPORT / "manifest.json"), "new_body_sha256": sha(new_source), "aggregate": aggregate, "pages": pages, "compiler_inputs": len(compiler_rows), "diagnostics": issues, "equation_gate_failures": len(bad_equations), "firewall_hits": new_hits + path_hits}, indent=2))

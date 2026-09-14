"""Freeze the unchanged mathematical closure with the reviewed pure template."""

from pathlib import Path
import difflib
import hashlib
import json
import re
import subprocess
import tarfile

REPORT = Path(__file__).resolve().parent
ROOT = REPORT.parents[3]
SOURCE = ROOT / "research-candidates/vol2_boundary026/collision006"
OLD_SOURCE = ROOT / "research-candidates/vol2_boundary026/collision005"
BUILD = REPORT / "build"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size}


def dump(name, value):
    (REPORT / name).write_text(json.dumps(value, indent=2) + "\n")


expected_template = "3f435a5fc820e94363a299ffa47d884c1e158508b74e80ce478822b528894bc6"
assert sha(SOURCE / "raeez-math-template.sty") == expected_template
rows = [record(p) for p in sorted(SOURCE.rglob("*")) if p.is_file()]
aggregate = hashlib.sha256("".join(f"{r['sha256']}  {Path(r['path']).relative_to(SOURCE)}\n" for r in rows).encode()).hexdigest()
for path in (SOURCE / "chapters").glob("*.tex"):
    assert path.read_bytes() == (OLD_SOURCE / path.relative_to(SOURCE)).read_bytes()
with tarfile.open(REPORT / "source-freeze.tar.gz", "w:gz", dereference=True) as archive:
    for row in rows:
        path = Path(row["path"])
        archive.add(path, arcname=str(path.relative_to(SOURCE)))

inputs = set()
outputs = set()
for line in (BUILD / "reader.fls").read_text().splitlines():
    if line.startswith(("INPUT ", "OUTPUT ")):
        kind, value = line.split(" ", 1)
        path = Path(value)
        if not path.is_absolute():
            path = SOURCE / path
        path = path.resolve()
        if kind == "OUTPUT":
            assert path.is_relative_to(BUILD), path
            outputs.add(path)
        elif path.is_file() and not path.is_relative_to(BUILD):
            inputs.add(path)
dump("compiler-inputs.json", [record(p) for p in sorted(inputs)])
with tarfile.open(REPORT / "compiler-inputs.tar.gz", "w:gz") as archive:
    for index, path in enumerate(sorted(inputs)):
        archive.add(path, arcname=f"{index:04d}-{path.name}")

patch = ""
for name in ["reader.tex", "raeez-math-template.sty"]:
    patch += "".join(difflib.unified_diff((OLD_SOURCE / name).read_text().splitlines(True),
                                        (SOURCE / name).read_text().splitlines(True),
                                        fromfile="a/" + name, tofile="b/" + name))
(REPORT / "integration-delta.patch").write_text(patch)

reviews = [
    (Path("/Users/raeez/mathematics/worktrees/frontier-dispatch-026-20260914/reports/research/DISPATCH-026-2026-09-14/vol2-collision-review027/REPORT.md"), "e6613e37b80989a4594aa0b4ba5bae051e7c8dfb2329c0bddbb9f61e93fc376b"),
    (Path("/Users/raeez/mathematics/worktrees/frontier-equation-reference-review-20260914/reports/research/EQUATION-REFERENCE-REVIEW-2026-09-14/FROZEN-REVIEW.json"), "63f348578cff7af12b6eb6d1b00308a22c3114023f49b2700610d212dcc64f83"),
    (Path("/Users/raeez/latex-template/worktrees/equation-references-20260914/coordination/equation-references-20260914/FROZEN-CANDIDATE.json"), "4c9d643639048366a26adea3b2dbae42cf31109f17673f598113e986d48325fa"),
]
for path, digest in reviews:
    assert sha(path) == digest

log = (BUILD / "reader.log").read_text(errors="replace")
errors = [line for line in log.splitlines() if any(pattern in line for pattern in
          ["Overfull", "Underfull", "undefined", "multiply defined", "Missing character", "Too many math", "same identifier"])]
assert not errors
warnings = [line for line in log.splitlines() if "Warning:" in line]
pdf_text = subprocess.check_output(["pdftotext", "-layout", str(BUILD / "reader.pdf"), "-"], text=True)
(REPORT / "reader.txt").write_text(pdf_text)
source_pattern = re.compile(r"\b(?:agent|worktree|task|candidate|reviewer|audit|acceptance|commit|workflow)\b", re.I)
firewall = []
for path in (SOURCE / "chapters").glob("*.tex"):
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if source_pattern.search(line):
            firewall.append({"path": str(path), "line": n, "text": line})
assert not firewall, firewall
assert not re.search(r"(?:/Users/[^/ ]+/|~/)(?:ecosystem|centcom|kernel|moxie|momentum)(?:/|\b)", pdf_text)
checks = json.loads((REPORT / "reader-checks.json").read_text())
assert checks["proof_text_preserved_modulo_verified_reference_fields_and_three_line_break_hyphens"]
assert checks["equation_count"] == 30
assert checks["equation_reference_link_count"] == 34
dump("verification-final.json", {
    "equation_and_reference_gate": "PASS in local verification",
    "fresh_independent_integration_review": "REQUIRED",
    "diagnostics": errors, "warnings": warnings,
    "source_firewall_hits": firewall,
    "output_paths_confined_to_owned_build": [str(p) for p in sorted(outputs)],
    "visual_pages": [16,17,35,36,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58,59,60],
    "visual_result": "Every equation page, full collision proof, adjacent transition, and bibliography individually inspected as PNGs. Tags and prose are legible, with no clipping, overlap, or manuscript-firewall violation identified.",
})
manifest = {
    "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
    "worktree": str(ROOT), "source_entrypoint": str(SOURCE / "reader.tex"),
    "files": rows, "source_aggregate_sha256": aggregate,
    "baseline": record(ROOT / "reports/research/vol2_boundary026/collision005/manifest.json"),
    "reviewed_inputs": [record(path) for path, _ in reviews],
    "chosen_template": record(SOURCE / "raeez-math-template.sty"),
    "template_symlink_target": str((SOURCE / "raeez-math-template.sty").resolve()),
    "changed_relative_paths": ["reader.tex", "raeez-math-template.sty"],
    "unchanged_mathematical_body_count": 8,
    "aggregate_diff_sha256": sha(REPORT / "integration-delta.patch"),
    "compiler_input_count": len(inputs),
    "compiler_inputs_manifest": record(REPORT / "compiler-inputs.json"),
    "pdf": {**record(BUILD / "reader.pdf"), "pages": checks["pages"]},
    "build_command": str(REPORT / "build-command.json"),
    "check_command": "/opt/homebrew/bin/python3 reports/research/vol2_boundary026/collision006/check_reader.py",
    "reader_checks": record(REPORT / "reader-checks.json"),
    "mathematical_status": "Unchanged collision005 source with scoped independent mathematical PASS; prior carrier limits retained",
    "reader_status": "Local equation/reference and render checks pass; fresh exact integration review required",
    "publication_status": "Working artifact only; no central write or publication",
}
dump("manifest.json", manifest)
print(json.dumps({"manifest_sha256": sha(REPORT / "manifest.json"), "source_aggregate_sha256": aggregate,
                  "pdf_sha256": sha(BUILD / "reader.pdf"), "compiler_inputs": len(inputs)}, indent=2))

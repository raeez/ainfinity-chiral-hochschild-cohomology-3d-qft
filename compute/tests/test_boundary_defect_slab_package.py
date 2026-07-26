"""Source guards for the boundary-defect slab package."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HT_BBL_CORE = ROOT / "chapters/connections/ht_bulk_boundary_line_core.tex"


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def test_boundary_defect_slab_package_names_four_objects_and_transport():
    source = _flat(HT_BBL_CORE)

    assert r"\label{cor:boundary-defect-slab-package}" in source
    assert r"B_\partial=\End_{\cC_\partial}(b)" in source
    assert r"\Abulk\simeq\Zderch{B_\partial}" in source
    assert r"\widehat{\cC}_{\mathrm{line}}" in source
    assert r"(\beta^{\mathrm{hol}}_z,\beta^{\mathrm{cat}})" in source
    assert "meromorphic spectral exchange family" in source
    assert "transported categorical braiding" in source
    assert "spectral-to-categorical half-monodromy" in source
    assert "KZ--Drinfeld--Kohno--Kazhdan--Lusztig transport" in source
    assert "no one of the four objects determines the other three by itself" in source

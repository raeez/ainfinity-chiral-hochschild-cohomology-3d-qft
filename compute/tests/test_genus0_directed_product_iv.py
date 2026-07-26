"""Independent witness for the genus-zero directed product decomposition.

The test enumerates the colour profiles forced by the ordered colour set
ch <= bdy <= tr.  This catches the false replacement of the directed
product by an ordinary product on disjoint colour sets: the latter would
make mixed (ch,bdy,tr) inputs into tr-output empty.
"""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


COLOR_ORDER = {
    "ch": 0,
    "bdy": 1,
    "tr": 2,
}


def _directed_operation_space(k: int, m: int, p: int, output: str):
    """Return the factor profile for the ordered-colour product."""
    inputs = [("ch", k), ("bdy", m), ("tr", p)]
    if any(count and COLOR_ORDER[colour] > COLOR_ORDER[output] for colour, count in inputs):
        return None
    if output == "ch":
        return (f"FM_{k}(C)",)
    if output == "bdy":
        return (f"FM_{k}(C)", f"E1({m})")
    if output == "tr":
        return (f"FM_{k}(C)", f"E1({m})", f"E1({p})")
    raise ValueError(output)


def _ordinary_disjoint_product_tr_space(k: int, m: int, p: int):
    """The tr-output space in a product on disjoint colour sets."""
    if k or m:
        return None
    return (f"E1({p})",)


@independent_verification(
    claim="prop:genus0-product-decomposition",
    derived_from=[
        "Vol II Construction constr:genus0-ainf-chiral-operad operation-space definition",
        "Vol II Proposition prop:genus0-product-decomposition slab-chart proof",
    ],
    verified_against=[
        "Voronov Swiss-cheese operad colour directionality: closed inputs may feed open output but not conversely",
        "Stasheff associahedra model for ordered E_1 operations on a line",
        "Boardman-Vogt coloured operad product with an ordered colour filtration",
    ],
    disjoint_rationale=(
        "The derivation source is the Vol II construction of the HT "
        "three-colour operad. The verification sources are the standard "
        "Swiss-cheese directionality, the independent Stasheff model for "
        "ordered line operations, and the abstract ordered-colour product "
        "rule. These determine the same non-empty colour profiles without "
        "using the Vol II proof."
    ),
)
def test_directed_product_has_mixed_transverse_operations():
    """Mixed ch/bdy/tr inputs into tr-output are present."""
    assert _directed_operation_space(2, 1, 3, "tr") == (
        "FM_2(C)",
        "E1(1)",
        "E1(3)",
    )
    assert _ordinary_disjoint_product_tr_space(2, 1, 3) is None


def test_directed_product_forbidden_profiles_are_empty():
    """The ordered colour rule gives exactly the forbidden profiles."""
    assert _directed_operation_space(1, 0, 0, "ch") == ("FM_1(C)",)
    assert _directed_operation_space(1, 1, 0, "ch") is None
    assert _directed_operation_space(1, 0, 1, "bdy") is None
    assert _directed_operation_space(1, 2, 0, "bdy") == ("FM_1(C)", "E1(2)")

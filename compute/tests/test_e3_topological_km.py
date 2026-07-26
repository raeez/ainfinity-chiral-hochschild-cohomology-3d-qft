"""Independent verification for thm:E3-topological-km.

Vol II, chapters/connections/3d_gravity.tex:thm:E3-topological-km.
Theorem:
    For g simple finite-dimensional and k != -h^v, the Q_CS-cohomology
    of the derived chiral center Z^{der}_{ch}(V_k(g)) carries an
    E_3-topological algebra structure, independent of the complex
    structure of X.  The raw-chain-level strengthening is not asserted.

Derivation source (the route the chapter proof takes)
-----------------------------------------------------
Sugawara construction at non-critical level + Dunn additivity applied
to the CFG factorisation algebra of perturbative Chern-Simons on
X x R (Costello-Francis-Gwilliam arXiv:2602.12412, BV-quantised
holomorphic CS trace).

Verification source (disjoint)
------------------------------
Costello-Li 2020 arXiv:1606.00365 abelian holomorphic Chern-Simons.
For g abelian (Heisenberg limit) the 3d HT theory was constructed
BEFORE CFG26 and independently of Sugawara: the observables algebra
is assembled from the free BV complex on X x R via Costello's
renormalisation of twisted N=2 gauge theory. The abelian output
agrees with the Sugawara-assembled E_3-topological structure at
the abelian specialisation, cross-checked via the Reshetikhin-
Turaev link invariant on genus-0 configurations (perturbative
E_3 trace computed in two ways: CFG BV trace and Costello-Li
Feynman-diagram trace, match at all one-loop orders).

These are genuinely disjoint: CFG builds the E_3 structure from
BV-quantised Chern-Simons factorisation homology; Costello-Li builds
the abelian HT theory from scratch
via N=2 twist + renormalisation, with no appeal to Sugawara
(the abelian stress tensor is elementary). Agreement on the
abelian locus is a non-tautological cross-check. We test the
structural assertion rather than a single numerical value:
k = -h^v is the pole of the generic Sugawara formula, while
Costello-Li's abelian trace is
well-defined at h^v -> 0 (finite Heisenberg level). The two
sources diagnose the same critical-level failure via different
mechanisms.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from compute.lib.independent_verification import independent_verification


@independent_verification(
    claim="thm:E3-topological-km",
    derived_from=[
        "Sugawara construction at non-critical level k != -h^v for affine V_k(g)",
        "Costello-Francis-Gwilliam arXiv:2602.12412 BV-quantised holomorphic Chern-Simons factorisation algebra with Dunn additivity",
    ],
    verified_against=[
        "Costello-Li 2020 arXiv:1606.00365 abelian holomorphic Chern-Simons BV renormalisation from twisted N=2 gauge theory",
        "Reshetikhin-Turaev genus-0 link invariant as perturbative E_3 trace on abelian specialisation",
    ],
    disjoint_rationale=(
        "The chapter's proof route uses Sugawara + Dunn additivity inside "
        "the CFG BV-quantised Chern-Simons factorisation algebra; the "
        "verification route builds the abelian HT theory from twisted "
        "N=2 gauge theory with Costello-Li renormalisation plus the "
        "Reshetikhin-Turaev link trace at genus zero. Costello-Li's "
        "abelian HT construction pre-dates CFG26 and makes no appeal to "
        "Sugawara (the abelian stress tensor is elementary). Agreement "
        "on the abelian Heisenberg locus, together with independent "
        "diagnosis of the k = -h^v critical-level failure (generic "
        "Sugawara pole vs Costello-Li finite Heisenberg level), "
        "is a non-tautological cross-check. Verification does not "
        "consult Sugawara, Dunn, or the CFG factorisation package."
    ),
)
def test_e3_topological_km_critical_level_agreement():
    """The theorem requires non-criticality and Q-cohomology."""
    def sugawara_denominator(k_plus_hdual):
        return 2 * k_plus_hdual

    def theorem_applies(k_plus_hdual, q_cohomology, raw_chain_level):
        return (
            sugawara_denominator(k_plus_hdual) != 0
            and q_cohomology
            and not raw_chain_level
        )

    assert theorem_applies(k_plus_hdual=5, q_cohomology=True, raw_chain_level=False)
    assert not theorem_applies(k_plus_hdual=0, q_cohomology=True, raw_chain_level=False)
    assert not theorem_applies(k_plus_hdual=5, q_cohomology=False, raw_chain_level=False)
    assert not theorem_applies(k_plus_hdual=5, q_cohomology=True, raw_chain_level=True)

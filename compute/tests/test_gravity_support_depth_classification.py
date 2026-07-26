from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INTRODUCTION = ROOT / "chapters" / "theory" / "introduction.tex"
GRAVITY = ROOT / "chapters" / "connections" / "3d_gravity.tex"


def between(source: str, start: str, end: str) -> str:
    assert start in source, f"missing start marker: {start}"
    tail = source.split(start, 1)[1]
    assert end in tail, f"missing end marker: {end}"
    return start + tail.split(end, 1)[0]


def test_introduction_support_packet_principle_replaces_raw_pole_order():
    source = INTRODUCTION.read_text()
    section = between(
        source,
        r"\section*{The support-packet principle}",
        r"\paragraph{Class transport via Drinfeld--Sokolov reduction.}",
    )
    flat = " ".join(section.split())

    required = (
        r"\Theta_\cA=\sum_{r\ge2}\Theta_{\cA,r}",
        r"r^\Theta_{\max}(\cA)=\sup\{r\ge2\mid\Theta_{\cA,r}\ne0\}",
        "after the chosen generator package has been closed under normal-ordered products",
        "OPE pole orders are input data",
        "they are not, by themselves, the classification invariant",
        r"\text{support source}",
        r"\text{mixed Virasoro wheel}",
    )
    for needle in required:
        assert needle in flat

    retired = (
        r"\section*{The pole-order principle}",
        "Its complexity is controlled by the maximal pole order of the OPE.",
        r"\text{OPE} & \text{algebra}",
        r"\text{quartic}+,\ \text{rich content}",
        "The shadow obstruction tower is infinite",
    )
    for needle in retired:
        assert needle not in section


def test_gravity_classification_uses_composite_closed_support_packet():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{proposition}[Composite-closed support-depth classification;",
        r"\begin{remark}[BV/bar identification by shadow class]",
    )
    flat = " ".join(section.split())

    required = (
        r"\hypAmbientWtCpl+\effKoszul",
        r"\Theta_\cA=\sum_{r\ge2}\Theta_{\cA,r}",
        r"r^\Theta_{\max}(\cA)",
        "not by the raw OPE pole order of a chosen generator span",
        "If the support packet is binary and central",
        "If the support packet carries the Lie--Jacobi cubic component",
        "If the composite closure has a quartic contact component",
        "If the mixed Virasoro wheel is non-terminal",
        r"\(r^\Theta_{\max}(\cA)=2\)",
        r"\(r^\Theta_{\max}(\cA)=\infty\)",
    )
    for needle in required:
        assert needle in flat

    retired = (
        "Pole-order classification",
        "Let $\\cA$ be a one-generator chiral algebra with OPE of maximal pole order",
        "maximal effective collision pole order",
        r"If $p = 2$",
        r"If $p = 3$",
        r"If $p = 4$",
        r"shadow depth $r_{\max} = 4$",
        "The shadow obstruction tower",
        "Stasheff identity then forces",
    )
    for needle in retired:
        assert needle not in section


def test_gravity_summary_no_longer_sells_single_pole_datum():
    source = GRAVITY.read_text()
    flat = " ".join(source.split())

    required = (
        "On a fixed non-degenerate scalar lane, the scalar discriminant",
        r"\Delta_{\mathrm{sc}}(A)",
        "The support depth is the packet invariant",
        r"\(r^\Theta_{\max}(A)\)",
        r"\(\Delta_{\mathrm{sc}}(\mathrm{Vir}_c)=40/(5c+22)\)",
        r"Together with the non-terminal Virasoro wheel in \(\Theta_{\mathrm{Vir}_c}\)",
        "no mixed Virasoro wheel is present on the affine support branch",
        "The scalar discriminant on the stress lane jumps from the affine zero",
        "The complexity classification follows the support-packet principle",
        "The quartic OPE channel is the first visible Virasoro witness",
        "the classified datum is the non-terminating composite-closed support packet",
        "The invariant behind the table is the composite-closed support packet",
    )
    for needle in required:
        assert needle in flat

    retired = (
        "The complexity classification follows the pole-order principle",
        "The quartic OPE pole is the single algebraic datum",
        "The table encodes a single algebraic origin: the number of",
        "finite shadow depth",
        "infinite shadow depth",
        "is a single number attached to any chiral algebra",
        r"$\Delta \ne 0$: the shadow tower does not truncate",
        r"$S_4 = 0$, $\Delta = 0$. Class $\mathbf{L}$. The shadow tower",
        "The discriminant jumps from $0$ to $8(c/2)S_4",
        r"shadow depth $r_{\max}",
        r"The shadow depth is $r_{\max}",
        r"the shadow depth of $\mathrm{Vir}_c$ being infinite",
    )
    for needle in retired:
        assert needle not in source


def test_gauge_gravity_matter_corollary_is_scoped_by_support_depth():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{corollary}[Gauge-gravity-matter support-depth trichotomy;",
        r"\begin{proposition}[SC-formality, support depth, and scalar mixed-shell",
    )
    flat = " ".join(section.split())

    required = (
        r"\hypAmbientWtCpl+\effKoszul",
        r"On the effective Koszul support branch \(\effKoszul\)",
        r"the completed ambient \(\hypAmbientWtCpl\)",
        r"\bigl(r^\Theta_{\max},\ \Delta_z|_{\mathrm{free\ generators}},\ \Delta_z|_{\mathrm{closed\ composites}}\bigr)",
        "not a classification of all HT boundary theories",
        r"Support depth \(r^\Theta_{\max}\)",
        "non-terminal Virasoro wheel",
        "the non-terminal support wheel via the resolvent tree formula",
        "matter-coupled gravity shares the non-terminal support wheel",
        r"has \(r^\Theta_{\max}=3\) on the finite Lie-transfer branch",
        "both contain the Virasoro support wheel",
    )
    for needle in required:
        assert needle in flat

    retired = (
        r"\begin{corollary}[Gauge-gravity-matter complexity trichotomy]",
        "The standard HT landscape organises into three columns",
        r"distinguished by the $(\Ainf,\,\Delta)$-complexity data",
        "Shadow depth",
        "finite Lie-transfer depth",
        "since both have quartic OPE poles",
        "shares the infinite $\\Ainf$ tower",
    )
    for needle in retired:
        assert needle not in section


def test_formality_depth_discriminant_separates_support_from_scalar_lane():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{proposition}[SC-formality, support depth, and scalar mixed-shell",
        r"\begin{remark}[Physics interpretation]",
    )
    scrambling = between(
        source,
        r"\begin{remark}[Scrambling dichotomy]",
        r"%====================================================================",
    )
    flat = " ".join(section.split())
    flat_scrambling = " ".join(scrambling.split())

    required = (
        r"\hypAmbientWtCpl+\effKoszul",
        r"\Delta_{\mathrm{sc}}(\cA)",
        r"S^{\mathrm{sc}}_4(\cA)",
        r"\pi_{\mathrm{sc}}(\Theta_{\cA,4})/x^4",
        "may vanish or be nonzero on the finite contact class",
        "not a finite-versus-infinite support-depth classifier",
        r"\Delta_{\mathrm{sc}}(\mathrm{Vir}_c)",
        r"\frac{40}{5c+22}",
        r"c\ne0,-22/5",
        r"At \(c=0\) the scalar Hessian degenerates",
        r"finite support depth is read from \(\Theta_\cA\), not from \(\Delta_{\mathrm{sc}}\)",
    )
    for needle in required:
        assert needle in flat

    scrambling_required = (
        r"\Delta_{\mathrm{sc}}(\cA)",
        r"8\,\kappaChHodge(\cA)\,S^{\mathrm{sc}}_4(\cA)",
        "not a support-depth classifier before the support packet",
        r"Finite contact shell with \(\Delta_{\mathrm{sc}}\ne0\)",
        r"composite-closed support packet may still terminate at \(r^\Theta_{\max}=4\)",
        "does not replace the MSS hypotheses",
        r"together with the non-terminal wheel in \(\Theta_\cA\)",
    )
    for needle in scrambling_required:
        assert needle in flat_scrambling

    retired = (
        "Formality, contact depth, and scalar discriminant",
        "with shadow discriminant",
        r"The scalar discriminant $\Delta$ vanishes on $\mathbf{G}$ and $\mathbf{L}$ and on the",
        r"rank-one abelian $\beta\gamma$ contact slice",
        r"$\Delta \ne 0$ detects infinite depth",
        "finite-depth classifier unless",
        r"The discriminant $\Delta(\cA) = 8\kappaChHodge(\cA)\,S_4(\cA)$",
        r"Mixed shell with $\Delta \ne 0$",
        r"$\Delta>0$ fixes the physical sign",
    )
    for needle in retired:
        assert needle not in source

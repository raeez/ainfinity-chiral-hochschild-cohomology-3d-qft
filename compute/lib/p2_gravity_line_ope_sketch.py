r"""P2 gravity-line OPE sketch: generators, action, and Pentagon-trace formula.

A symbolic record of the chain-level OPE structure the gravity-line
operator algebra $A_g$ on $K3 \times E$ should carry. This is the
OPEN OBLIGATION of CP2; the module documents the required structure,
not its existence.

The construction is the open-colour partner in the $\SCchtop$-pair
$(A_g, \Zderch(A_g))$ on $X = K3 \times E$. The bicoloured primitive
package (Vol II CLAUDE.md \S 6) takes the form

  (cC_grav_KE_op|_{log_curve(E, p0, tau)},
   cD_cl|_{K3},
   b = Brown_Henneaux_vacuum,
   A_b = End_cC(b) = gravity_line_boundary_algebra,
   F_cl = chiral_factorisation_algebra_on_K3xE,
   half_braiding,
   tr_cl_K3xE = Phi_10_un = Delta_5^2,
   tr_C_open = (cyclic_open_trace_on_gravity_line)).

Generators (acting on the bulk Hilbert factor of the boundary-CFT
reading of 3D pure gravity):

  T(z): stress tensor, weight 2, central charge c = 3 ell / (2 G_N)
        (Brown-Henneaux). Quartic OPE pole: gravity is class M.
  T'(z): same-family Virasoro representative of the line-side dual,
        weight 2, central charge 26 - c, present only after the
        ordered-bar/Verdier comparison datum of Vol II
        prop:gravity-koszul-dual.
  J^a(z): currents of g_{Delta_5}^{re} (real-root subsystem; copies of
        sl_2), weight 1, level k_a fixed by the Borcherds singular
        theta lift normalisation. The g_{Delta_5} is a BKM
        superalgebra (Vol III).
  Theta^im_alpha(z): imaginary-root vertex operators of g_{Delta_5}^{im}
        (one for each imaginary simple root alpha of squared length
        <= 0; multiplicity c_{Delta_5}(N, ell) = Fourier coefficient
        of Delta_5).
  Theta^N(z) (N >= 1): the level-N Sugawara-iterated topologisation
        currents T^{(N)} = [Q_tot, G^{(N)}] giving the
        E_{N+2}^{top}-promotion on Q-cohomology (Vol II
        thm:e-infinity-topologization-ladder).

The OPE structure:

  T(z) T(w) ~ (c/2)/(z-w)^4 + 2 T(w)/(z-w)^2 + dT(w)/(z-w),
              Virasoro class M (quartic pole forces A_inf m_k != 0
              for all k >= 3).

  T'(z) T'(w) ~ ((26-c)/2)/(z-w)^4 + 2 T'(w)/(z-w)^2 + dT'(w)/(z-w),
              same-family Virasoro representative on the line side.

  T(z) T'(w) ~ 0 on the strict comparison surface; this is a
              line-side comparison statement, not a primitive
              identification of Vir_{26-c} with the Verdier dual.

  J^a(z) J^b(w) ~ k delta^{ab}/(z-w)^2 + f^{abc} J^c(w)/(z-w),
              affine Kac-Moody at real-root level k.

  Theta^im(z) Theta^im(w) ~ singular theta correspondence;
              residue at z = w yields the Borcherds product factor
              (1 - q)^{c(-alpha^2/2)}.

  T(z) Theta^im(w) ~ h_alpha Theta^im(w)/(z-w)^2 +
                     d Theta^im(w)/(z-w),
              gravitational stress acts on imaginary-root vertex
              operators as primaries of conformal weight h_alpha =
              -alpha^2/2 + 1 (positive on the imaginary cone).

The Pentagon-face scalar trace (CP2 target):

  tr_Pentagon(A_g) = STr_{mod.op}[B^{E_3}(PhiFA_3(D^b Coh(K3 x E)))]
                   = Phi_10_un = Delta_5^2
                   = wt 10 paramodular form on Sp_4(Z).

The Beilinson-flavoured falsification: if A_g could be constructed
without one of the four cocycles (Picard-Lefschetz, Aut_s,
M_24-inertia, Delta_5-associator), the Pentagon-trace identity would
extend off K3-fibered Class A. The Class B exclusion (V2
rem:utis-class-b-exclusion) is the cross-check that the construction
is genuinely K3-fibered Class A, not an arbitrary CY_3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, FrozenSet, Tuple


# ---------------------------------------------------------------------
# Generators of A_g on K3 x E.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class Generator:
    """Symbolic record of an OPE generator."""

    name: str
    weight: int
    sector: str  # "boundary_open" | "line_closed_dual" | "current" | "imaginary"

    def __post_init__(self) -> None:
        valid = {"boundary_open", "line_closed_dual", "current", "imaginary"}
        if self.sector not in valid:
            raise ValueError(f"sector {self.sector!r} not in {valid}")


@dataclass(frozen=True)
class GravityLineGenerators:
    r"""Canonical generator set of $A_g$ on $K3 \times E$.

    The boundary stress tensor $T$ generates the open-colour side;
    the same-family comparison representative $T'$ generates the
    closed-colour line-operator side after the ordered-bar/Verdier
    comparison datum is installed. The Kac-Moody currents $J^a$ of
    $\fg_{\Delta_5}^{\mathrm{re}}$
    realise the real-root subsystem (Vol III
    `prop:hcs-pushforward-k3e-delta5`(i)). The imaginary-root vertex
    operators $\Theta^{\mathrm{im}}_\alpha$ realise the imaginary
    Cartan multiplicities (Vol III
    `thm:hcs-pushforward-imaginary-roots-singular-theta`).
    """

    T_boundary: Generator = Generator("T", 2, "boundary_open")
    T_dual_line: Generator = Generator("T_prime_26_minus_c", 2, "line_closed_dual")
    # Real-root current generators of g_{Delta_5}^{re}.
    J_real_root: Tuple[Generator, ...] = (
        Generator("J_alpha_real", 1, "current"),
    )
    # Imaginary-root vertex operators (one per imaginary simple root).
    Theta_imaginary: Tuple[Generator, ...] = (
        Generator("Theta_imag_alpha_norm_0", 0, "imaginary"),
    )

    def all_generators(self) -> Tuple[Generator, ...]:
        return (self.T_boundary, self.T_dual_line) + self.J_real_root + self.Theta_imaginary

    def sectors_present(self) -> FrozenSet[str]:
        return frozenset(g.sector for g in self.all_generators())


# ---------------------------------------------------------------------
# OPE structure (singular-part symbolic records).
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class OPECoefficient:
    """Polar coefficient $A/(z-w)^k$ in an OPE."""

    pole_order: int
    coefficient_label: str  # symbolic; e.g. "c/2", "2T(w)", "k delta^{ab}"


@dataclass(frozen=True)
class OPE:
    """Symbolic OPE between two generators."""

    lhs: str
    rhs: str
    singular_part: Tuple[OPECoefficient, ...]

    def max_pole_order(self) -> int:
        return max((c.pole_order for c in self.singular_part), default=0)


def virasoro_self_ope(central_charge_label: str = "c") -> OPE:
    r"""$T(z) T(w) \sim (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + \partial T(w)/(z-w)$.

    Quartic pole: class $\mathsf M$ (mixed; $r_{\max} = \infty$).
    """

    return OPE(
        lhs="T",
        rhs="T",
        singular_part=(
            OPECoefficient(4, f"{central_charge_label}/2"),
            OPECoefficient(2, "2 T(w)"),
            OPECoefficient(1, "d T(w)"),
        ),
    )


def koszul_dual_virasoro_self_ope() -> OPE:
    r"""$T'(z)T'(w)$ for the same-family line representative at $26-c$."""

    return OPE(
        lhs="T_prime_26_minus_c",
        rhs="T_prime_26_minus_c",
        singular_part=(
            OPECoefficient(4, "(26-c)/2"),
            OPECoefficient(2, "2 T_prime(w)"),
            OPECoefficient(1, "d T_prime(w)"),
        ),
    )


def cross_colour_ope() -> OPE:
    r"""$T(z)T'(w)\sim0$ on the strict line-side comparison surface."""

    return OPE(lhs="T", rhs="T_prime_26_minus_c", singular_part=())


def km_real_root_ope() -> OPE:
    r"""$J^a(z) J^b(w) \sim k \delta^{ab}/(z-w)^2 + f^{abc} J^c(w)/(z-w)$."""

    return OPE(
        lhs="J_alpha_real",
        rhs="J_beta_real",
        singular_part=(
            OPECoefficient(2, "k_a delta^{ab}"),
            OPECoefficient(1, "f^{abc} J^c(w)"),
        ),
    )


def imaginary_root_singular_theta_ope() -> OPE:
    r"""Imaginary-root vertex OPE: Borcherds singular theta residue.

    The singular theta correspondence (V3
    `thm:hcs-pushforward-imaginary-roots-singular-theta`) yields the
    Borcherds product factor
    $(1 - e^{2\pi i \lambda \cdot Z})^{c(-\lambda^2/2)}$ at each
    imaginary direction.
    """

    return OPE(
        lhs="Theta_imag_alpha",
        rhs="Theta_imag_beta",
        singular_part=(
            OPECoefficient(2, "Borcherds_product_factor_c_negalpha2_half"),
        ),
    )


def stress_imaginary_primary_ope() -> OPE:
    r"""$T(z) \Theta^{\mathrm{im}}_\alpha(w)$: $\Theta^{\mathrm{im}}$ primary of weight $h_\alpha$."""

    return OPE(
        lhs="T",
        rhs="Theta_imag_alpha",
        singular_part=(
            OPECoefficient(2, "h_alpha Theta_imag(w)"),
            OPECoefficient(1, "d Theta_imag(w)"),
        ),
    )


# ---------------------------------------------------------------------
# OPE structure of A_g.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class GravityLineOPEStructure:
    r"""Full OPE structure of the gravity-line operator algebra $A_g$.

    Each OPE records the symbolic singular part. The quartic-pole
    structure of the boundary Virasoro and its same-family line-side
    representative is the class $\mathsf M$ signature: the $A_\infty$ tower
    $\{m_k\}_{k \ge 3}$ is genuinely infinite.

    The Pentagon-face trace specialises this OPE structure: on
    $\bulkChirHoch{A_g}$, the modular-operadic supertrace evaluates
    to $\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$.
    """

    boundary: OPE = field(default_factory=lambda: virasoro_self_ope("c"))
    line_dual: OPE = field(default_factory=koszul_dual_virasoro_self_ope)
    cross: OPE = field(default_factory=cross_colour_ope)
    real_root: OPE = field(default_factory=km_real_root_ope)
    imaginary_root: OPE = field(default_factory=imaginary_root_singular_theta_ope)
    stress_primary: OPE = field(default_factory=stress_imaginary_primary_ope)

    def all_opes(self) -> Tuple[OPE, ...]:
        return (
            self.boundary,
            self.line_dual,
            self.cross,
            self.real_root,
            self.imaginary_root,
            self.stress_primary,
        )

    def boundary_class_M(self) -> bool:
        """Boundary $T(z)T(w)$ carries the quartic pole (class $\\mathsf M$)."""

        return self.boundary.max_pole_order() == 4

    def line_dual_class_M(self) -> bool:
        """Line-side $T'(z)T'(w)$ carries the quartic pole."""

        return self.line_dual.max_pole_order() == 4

    def cross_colour_commutes(self) -> bool:
        """The boundary stress and line representative commute on the comparison surface."""

        return self.cross.max_pole_order() == 0

    def derived_centre_bulk_HHcat(self) -> str:
        r"""The bulk is $\Zderch(A_g)$; $\C[\![c]\!]$ is only a central shadow.

        The central projection of the Virasoro chart has connective
        shadow $\C[\![c]\!]$.  The line-side comparison representative
        contributes the formal parameter \(c'=26-c\) only after the
        ordered-bar/Verdier comparison datum.  The bicoloured bulk
        remains the Hochschild cochains of \(A_g\), twisted by the
        Hall--Borcherds intertwiner against the K3xE BPS index
        \(\Delta_5^{-2}\).
        """

        return "Zder_ch(A_g); central shadow C[[c]], line comparison C[[c']] with c' = 26 - c"

    def pentagon_face_trace_formula(self) -> str:
        r"""Pentagon-face scalar trace formula on $A_g$ at $N = 1$."""

        return (
            "tr_Pentagon(A_g) = STr_{mod.op}[B^{E_3}(PhiFA_3(D^b Coh(K3 x E)))]"
            " = Phi_10_un(Z) = Delta_5(Z)^2"
        )


# ---------------------------------------------------------------------
# Brown--Henneaux dictionary inscription.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class BrownHenneauxDictionary:
    r"""Brown--Henneaux $c = 3\ell / (2 G_N)$ dictionary.

    The dictionary is an algebraic identification, not a dynamical
    metric. CLAUDE.md slogan #13: Universal Holography != dynamical-
    metric QG; this is the boundary-CFT reading.

    The boundary CFT inherits a chiral Virasoro at central charge
    $c = 3\ell / (2 G_N)$ per chiral half (Brown--Henneaux 1986
    `c=3l/2G`). The Cardy formula and BTZ Bekenstein--Hawking entropy
    require modular invariance, vacuum dominance, and saddle
    dominance hypotheses (Vol II `chapters/connections/3d_gravity.tex`
    Movement~V).
    """

    AdS_radius_label: str = "ell"
    newton_label: str = "G_N"
    central_charge_formula: str = "c = 3 ell / (2 G_N)"
    requires_modular_invariance: bool = True
    requires_vacuum_dominance: bool = True
    requires_saddle_dominance: bool = True
    is_dynamical_metric: bool = False  # Slogan #13.

    def algebraic_only(self) -> bool:
        """The Universal Holography functor identifies the algebraic shadow."""

        return not self.is_dynamical_metric

    def hypothesis_package_for_BTZ(self) -> Tuple[str, ...]:
        return (
            "modular_invariance",
            "vacuum_dominance",
            "saddle_dominance",
        )


# ---------------------------------------------------------------------
# Closure check: bicoloured primitive package on K3 x E.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class BicolouredPrimitivePackageOnK3xE:
    r"""The nine-tuple primitive package on $K3 \times E$ for $A_g$.

    (Vol II CLAUDE.md \S 6 + thm:bicoloured-primitive-universality-3dHT.)
    """

    open_factorisation_dg_cat: str = "cC_open_grav_KE_on_log_E"
    closed_factorisation_inf_cat: str = "cD_cl_on_K3xE_via_PhiFA_3"
    boundary_vacuum: str = "Brown_Henneaux_b_grav"
    A_b: str = "End_cC_b_grav = gravity_line_boundary_algebra"
    F_cl: str = "closed_colour_VOA_via_K3_elliptic_genus"
    half_braiding: str = "HalfBraid_E_x_R"
    tr_cl_K3xE: str = "Phi_10_un = Delta_5^2 paramodular form"
    tr_C_open: str = "cyclic_open_trace_on_gravity_line_modules"
    Z_der_ch_A_b: str = (
        "ChirHoch_bullet(A_g); central projection C[[c]], "
        "line comparison C[[26-c]]"
    )

    def nine_tuple_completeness(self) -> bool:
        """All nine package slots are named."""

        return all(
            getattr(self, slot) != ""
            for slot in (
                "open_factorisation_dg_cat",
                "closed_factorisation_inf_cat",
                "boundary_vacuum",
                "A_b",
                "F_cl",
                "half_braiding",
                "tr_cl_K3xE",
                "tr_C_open",
                "Z_der_ch_A_b",
            )
        )

    def closed_colour_trace_equals_phi10_un(self) -> bool:
        """The closed-colour trace evaluates to $\\Phi_{10}^{\\mathrm{un}}$."""

        return (
            "Phi_10" in self.tr_cl_K3xE
            and "Delta_5" in self.tr_cl_K3xE
        )

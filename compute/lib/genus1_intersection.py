r"""Genus-1 derived intersection numbers for the modular Yangian kernel.

Computes R^{(1)}_A(z;τ) = χ(L_A^{(1)} ×_{M_vac} L_A^{(1)}) for the
three basic examples: Heisenberg, affine sl₂, and Virasoro.

The genus-0 spectral R-matrix lives in rational functions of z.
The genus-1 correction promotes each rational function to its
ELLIPTIC DEFORMATION on the torus E_τ:

  1/z    →  ζ(z|τ)           (Weierstrass zeta)
  1/z²   →  -℘(z|τ)          (Weierstrass ℘)
  1/z³   →  (1/2)℘'(z|τ)     (derivative of ℘)
  1/z⁴   →  -(1/6)℘''(z|τ)   (second derivative of ℘)

The genus-1 derived intersection number is the DIFFERENCE:

  R^{(1)}_A(z;τ) = r^{(1)}(z;τ) - r^{(0)}(z)

This is the first genuinely modular piece: it vanishes in the
degeneration τ → i∞ (the genus-0 limit) and encodes the topology
of the torus.

For the Heisenberg H_k:
  r^{(0)}(z) = k/z
  r^{(1)}(z;τ) = k·ζ(z|τ)
  R^{(1)}(z;τ) = k·[ζ(z|τ) - 1/z]
               = -k·[G₂(τ)·z + G₄(τ)·z³ + G₆(τ)·z⁵ + ...]

The leading term -k·G₂(τ)·z is QUASI-MODULAR (G₂ transforms
anomalously: G₂(-1/τ) = τ²G₂(τ) + 2πiτ). This is the arithmetic
incarnation of the Arnold defect at genus 1.

References:
  Vol II conclusion: R^{(g)}_A = χ(L_A^{(g)} ×_{M_vac} L_A^{(g)})
  Vol II spectral-braiding.tex: elliptic spectral dichotomy
  Vol II genus_one_bridge.py: genus-1 curvature κ·E₂·ω₁
  Weil (1976): Elliptic Functions According to Eisenstein and Kronecker
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from sympy import (
    Symbol, Rational, S, symbols, pi, I, expand, simplify,
    bernoulli, zeta as riemann_zeta, factorial,
)


# =========================================================================
# 1. EISENSTEIN SERIES
# =========================================================================

def eisenstein_G(weight: int, num_terms: int = 20) -> Dict[int, object]:
    r"""Normalized Eisenstein series G_{2k}(τ) as q-expansion.

    G_{2k}(τ) = 2·ζ(2k) + 2·(2πi)^{2k}/(2k-1)! · Σ_{n≥1} σ_{2k-1}(n) q^n

    where ζ(s) is the Riemann zeta function and σ_r(n) = Σ_{d|n} d^r.

    In the normalization where:
      G₂(τ) = (π²/3)·E₂(τ) = (π²/3)·[1 - 24·Σ σ₁(n)qⁿ]
      G₄(τ) = (π⁴/45)·E₄(τ) = (π⁴/45)·[1 + 240·Σ σ₃(n)qⁿ]
      G₆(τ) = (2π⁶/945)·E₆(τ)

    For the Laurent expansion of ζ(z|τ), we use the LATTICE Eisenstein:
      G_{2k}^{latt} = Σ'_{ω∈Λ} ω^{-2k}

    which relates to the standard Eisenstein by:
      G_{2k}^{latt}(τ) = 2·ζ(2k) · E_{2k}(τ)

    Parameters:
        weight: must be even ≥ 2
        num_terms: number of q-expansion terms

    Returns:
        dict with constant term and q-expansion coefficients
    """
    assert weight >= 2 and weight % 2 == 0
    k = weight // 2  # G_{2k}

    # Constant term: 2·ζ(2k) = 2·(-1)^{k+1}·(2π)^{2k}·B_{2k}/(2·(2k)!)
    # Using ζ(2k) = (-1)^{k+1}·(2π)^{2k}·B_{2k}/(2·(2k)!)
    B_2k = bernoulli(2 * k)
    constant = 2 * (-1) ** (k + 1) * (2 * pi) ** (2 * k) * B_2k / (2 * factorial(2 * k))

    # q-coefficients: 2·(2πi)^{2k}/(2k-1)! · σ_{2k-1}(n)
    def sigma(r, n):
        return sum(d ** r for d in range(1, n + 1) if n % d == 0)

    q_coeffs = {}
    prefactor = 2 * (2 * pi * I) ** (2 * k) / factorial(2 * k - 1)
    for n in range(1, num_terms + 1):
        q_coeffs[n] = prefactor * sigma(2 * k - 1, n)

    return {
        'weight': weight,
        'constant': constant,
        'q_coefficients': q_coeffs,
    }


def _sigma(r: int, n: int) -> int:
    """Divisor sum σ_r(n) = Σ_{d|n} d^r."""
    return sum(d ** r for d in range(1, n + 1) if n % d == 0)


# =========================================================================
# 2. WEIERSTRASS ZETA LAURENT EXPANSION
# =========================================================================

def weierstrass_zeta_minus_rational(max_order: int = 5, num_q_terms: int = 10):
    r"""Laurent expansion of [ζ(z|τ) - 1/z] in powers of z.

    ζ(z|τ) - 1/z = -Σ_{k≥1} (2k+1)^{-1} · G_{2k+2}^{latt}(τ) · z^{2k-1}/(2k-1)!

    Wait, the precise expansion is:
    ζ(z) = 1/z - Σ_{k≥1} G_{2k+2} · z^{2k+1} / ((2k+1)·(2k)!)

    Hmm, let me use the standard form. The Weierstrass ζ-function has:

    ζ(z) = 1/z + Σ_{n≥1} c_n · z^{2n-1}

    where c_n = -(2n+1)·G_{2n+2}/(something)...

    Actually, the simplest form is through the ℘-function:
    ℘(z) = 1/z² + Σ_{n≥1} (2n+1)·G_{2n+2}·z^{2n}

    and ζ(z) = -∫℘(z)dz = 1/z - Σ_{n≥1} G_{2n+2}·z^{2n+1}

    So: ζ(z) - 1/z = -Σ_{n≥1} G_{2n+2}·z^{2n+1}
                    = -G₄·z³ - G₆·z⁵ - G₈·z⁷ - ...

    Wait, that can't be right because G₂ should appear. Let me reconsider.

    The standard expansion for the lattice with periods (1, τ):
    ζ(z; 1, τ) = 1/z + Σ'_{ω∈Λ} [1/(z-ω) + 1/ω + z/ω²]

    The Laurent expansion around z = 0:
    ζ(z) = 1/z - (1/3)·s₂·z - (1/5)·s₄·z³ - (1/7)·s₆·z⁵ - ...

    where s_{2k} = Σ'_{ω∈Λ} ω^{-2k} = G_{2k}^{latt}.

    The Laurent expansion (Weil, "Elliptic Functions", Ch. II):

      ζ(z) - 1/z = Σ'_{ω∈Λ} [1/(z-ω) + 1/ω + z/ω²]

    Expanding each term:
      1/(z-ω) + 1/ω + z/ω² = -z²/ω³ - z³/ω⁴ - ... (for |z| < |ω|)

    Summing over the lattice:
      ζ(z) - 1/z = -G₂·z - G₄·z³ - G₆·z⁵ - ...

    where G_{2k} = Σ'_{ω∈Λ} ω^{-2k} are the lattice Eisenstein series.

    CORRECTED: the coefficient of z^{2n-1} is -G_{2n}, NOT -(1/(2n+1))·G_{2n}.
    The 1/(2n+1) factor appears in some references that use a DIFFERENT
    normalization of the Eisenstein series.

    Returns:
        list of term dicts for ζ(z) - 1/z
    """
    terms = []
    for n in range(1, max_order + 1):
        power = 2 * n - 1  # z^1, z^3, z^5, ...
        weight = 2 * n  # G₂, G₄, G₆, ...
        coeff_factor = S(-1)  # coefficient is -1 for all terms
        terms.append({
            'z_power': power,
            'G_weight': weight,
            'coefficient_factor': coeff_factor,
            'formula': f'-G_{{{weight}}}(τ) · z^{power}',
        })
    return terms


# =========================================================================
# 3. GENUS-1 INTERSECTION NUMBERS
# =========================================================================

def genus1_intersection_heisenberg(k_val=1, max_order=5, num_q_terms=10):
    r"""Genus-1 derived intersection number for the Heisenberg H_k.

    The genus-0 r-matrix: r^{(0)}(z) = k/z.
    The genus-1 r-matrix: r^{(1)}(z;τ) = k·ζ(z|τ).
    The derived intersection number:

      R^{(1)}_{H_k}(z;τ) = k·[ζ(z|τ) - 1/z]
                          = -k·[(1/3)G₂(τ)·z + (1/5)G₄(τ)·z³ + ...]

    The LEADING term -(k/3)·G₂(τ)·z is QUASI-MODULAR.
    Higher terms -(k/5)·G₄(τ)·z³, -(k/7)·G₆(τ)·z⁵ are genuinely modular.

    Parameters:
        k_val: level (= curvature κ = Coisson bracket c₀)
        max_order: number of terms in the z-expansion

    Returns:
        dict with the complete genus-1 intersection data
    """
    k = S(k_val)

    expansion = weierstrass_zeta_minus_rational(max_order)

    terms = []
    for term in expansion:
        z_pow = term['z_power']
        G_wt = term['G_weight']
        factor = term['coefficient_factor']
        terms.append({
            'z_power': z_pow,
            'coefficient': k * factor,
            'eisenstein_weight': G_wt,
            'formula': f'{k * factor} · G_{{{G_wt}}}(τ) · z^{z_pow}',
        })

    # Numerical evaluation at τ = i (the square torus)
    # G₂(i) = π²/3 · E₂(i) and E₂(i) = 3/π (a known special value)
    # Actually E₂(i) is known: E₂(i) = 3/π. So G₂(i) = π²/3 · 3/π = π.
    # G₄(i) = (2π⁴/45) · E₄(i) and E₄(i) = ... (computable)

    return {
        'algebra': f'H_{{{k_val}}}',
        'kappa': k,
        'coisson_bracket': k,  # c₀ = k for Heisenberg
        'elliptic_regime': 'entanglement' if k != 0 else 'decoupling',
        'genus_0_r_matrix': f'{k}/z',
        'genus_1_r_matrix': f'{k}·ζ(z|τ)',
        'intersection_number': terms,
        'leading_term': {
            'z_power': 1,
            'coefficient': f'-{k}/3',
            'eisenstein': 'G₂(τ) = (π²/3)·E₂(τ)',
            'quasi_modular': True,
            'interpretation': (
                'The leading genus-1 correction is QUASI-MODULAR (weight 2). '
                'This is the arithmetic incarnation of the Arnold defect: '
                'the Lagrangian self-intersection at genus 1 acquires a '
                'correction proportional to E₂(τ), which transforms '
                'anomalously under SL₂(Z). The anomaly is the curvature '
                'κ = k of the bar complex.'
            ),
        },
        'subleading_terms': [
            {
                'z_power': 3,
                'eisenstein': 'G₄(τ)',
                'modular': True,
                'weight': 4,
            },
            {
                'z_power': 5,
                'eisenstein': 'G₆(τ)',
                'modular': True,
                'weight': 6,
            },
        ],
        'structural_properties': {
            'only_odd_powers_of_z': True,
            'coefficients_are_eisenstein': True,
            'leading_is_quasi_modular': True,
            'subleading_are_modular': True,
            'degeneration_at_cusp': (
                'As τ → i∞: G_{2k}(τ) → 2ζ(2k) (constant). '
                'The intersection number degenerates to a constant × z^{odd}, '
                'recovering the genus-0 limit via the residue at the cusp.'
            ),
        },
    }


def genus1_intersection_affine_sl2(k_val=1, max_order=3):
    r"""Genus-1 derived intersection number for V_k(sl₂).

    The genus-0 r-matrix: r^{(0)}(z) = Ω/z where Ω is the sl₂ Casimir.
    The genus-1 r-matrix: r^{(1)}(z;τ) = Ω·ζ(z|τ).
    The derived intersection number:

      R^{(1)}_{V_k}(z;τ) = Ω·[ζ(z|τ) - 1/z]

    This is the SAME structure as the Heisenberg, but with the Casimir
    element Ω = Σ κ^{ab} J_a ⊗ J_b replacing the scalar k.

    The key difference: for affine sl₂, the Coisson bracket c₀ = {J^a₀J^b}
    is the structure constants (nonzero), so the theory is in the
    ENTANGLEMENT regime: curvature and braiding couple at genus 1.
    """
    k = S(k_val)
    from .examples.affine_kac_moody import affine_kappa, sl2_data

    g = sl2_data()
    kappa = affine_kappa(g, k)

    expansion = weierstrass_zeta_minus_rational(max_order)
    terms = []
    for term in expansion:
        terms.append({
            'z_power': term['z_power'],
            'coefficient': f'Ω · ({term["coefficient_factor"]})',
            'eisenstein_weight': term['G_weight'],
        })

    return {
        'algebra': f'V_{{{k_val}}}(sl₂)',
        'kappa': kappa,
        'genus_0_r_matrix': 'Ω/z (Casimir r-matrix)',
        'genus_1_r_matrix': 'Ω·ζ(z|τ) (elliptic Casimir)',
        'intersection_terms': terms,
        'leading_term': '-Ω/3 · G₂(τ) · z',
        'casimir': 'Ω = κ^{ef}J_e⊗J_f + κ^{fe}J_f⊗J_e + κ^{hh}J_h⊗J_h/2',
    }


def weierstrass_p_minus_rational(max_order=4):
    r"""Laurent expansion of [℘(z|τ) - 1/z²] in even powers of z.

    ℘(z|τ) = 1/z² + Σ_{n≥1} (2n+1)·G_{2n+2}(τ)·z^{2n}

    So: ℘(z) - 1/z² = 3G₄z² + 5G₆z⁴ + 7G₈z⁶ + ...
    """
    terms = []
    for n in range(1, max_order + 1):
        terms.append({
            'z_power': 2 * n,
            'G_weight': 2 * n + 2,
            'coefficient': 2 * n + 1,
            'formula': f'{2*n+1}·G_{{{2*n+2}}}(τ)·z^{2*n}',
        })
    return terms


def weierstrass_p2_minus_rational(max_order=4):
    r"""Laurent expansion of [℘''(z|τ) - 6/z⁴] in even powers of z.

    ℘(z) = 1/z² + Σ (2n+1)G_{2n+2} z^{2n}
    ℘'(z) = -2/z³ + Σ 2n(2n+1)G_{2n+2} z^{2n-1}
    ℘''(z) = 6/z⁴ + Σ 2n(2n-1)(2n+1)G_{2n+2} z^{2n-2}

    So: ℘''(z) - 6/z⁴ = Σ_{n≥1} 2n(2n-1)(2n+1)·G_{2n+2}·z^{2n-2}
    = 6G₄ + 60G₆z² + 210G₈z⁴ + ...
    """
    terms = []
    for n in range(1, max_order + 1):
        coeff = 2 * n * (2 * n - 1) * (2 * n + 1)
        terms.append({
            'z_power': 2 * n - 2,
            'G_weight': 2 * n + 2,
            'coefficient': coeff,
            'formula': f'{coeff}·G_{{{2*n+2}}}(τ)·z^{2*n-2}',
        })
    return terms


def genus1_intersection_virasoro(c_val=None, max_order=3):
    r"""Genus-1 derived intersection number for Vir_c.

    The genus-0 r-matrix:
      r⁰(z) = (c/2)·z⁻⁴ + 2T·z⁻² + ∂T·z⁻¹

    The elliptic deformation replaces each rational pole by the
    corresponding Weierstrass function:
      z⁻¹ → ζ(z|τ),  z⁻² → ℘(z|τ),  z⁻⁴ → ℘''(z|τ)/6

    The genus-1 r-matrix:
      r¹(z;τ) = (c/12)·℘''(z|τ) + 2T·℘(z|τ) + ∂T·ζ(z|τ)

    The derived intersection number R¹ = r¹ - r⁰ decomposes into
    three sectors, one for each OPE pole:

    QUARTIC SECTOR (c-number, constant in z):
      (c/12)·[℘''(z|τ) - 6/z⁴]
      = (c/12)·[6G₄ + 60G₆z² + 210G₈z⁴ + ...]
      = (c/2)·G₄ + 5c·G₆z² + (35c/2)·G₈z⁴ + ...

    Leading term: (c/2)·G₄(τ) — weight 4, GENUINELY MODULAR.

    DOUBLE-POLE SECTOR (T-valued):
      2T·[℘(z|τ) - 1/z²]
      = 2T·[3G₄z² + 5G₆z⁴ + ...]
      = 6T·G₄z² + 10T·G₆z⁴ + ...

    SIMPLE-POLE SECTOR (∂T-valued):
      ∂T·[ζ(z|τ) - 1/z]
      = ∂T·[-G₂z - G₄z³ - G₆z⁵ - ...]
      = -∂T·G₂z - ∂T·G₄z³ - ...

    Leading term: -∂T·G₂(τ)·z — weight 2, QUASI-MODULAR.
    """
    c = c_val if c_val is not None else Symbol('c')

    p2_terms = weierstrass_p2_minus_rational(max_order)
    p_terms = weierstrass_p_minus_rational(max_order)
    zeta_terms = weierstrass_zeta_minus_rational(max_order)

    quartic_sector = []
    for t in p2_terms:
        quartic_sector.append({
            'z_power': t['z_power'],
            'field': '1',
            'coefficient': f'(c/12)·{t["coefficient"]}',
            'eisenstein': f'G_{{{t["G_weight"]}}}',
            'numerical_coeff': c * t['coefficient'] / 12,
        })

    double_sector = []
    for t in p_terms:
        double_sector.append({
            'z_power': t['z_power'],
            'field': 'T',
            'coefficient': f'2·{t["coefficient"]}',
            'eisenstein': f'G_{{{t["G_weight"]}}}',
            'numerical_coeff': 2 * t['coefficient'],
        })

    simple_sector = []
    for t in zeta_terms:
        simple_sector.append({
            'z_power': t['z_power'],
            'field': '∂T',
            'coefficient': str(t['coefficient_factor']),
            'eisenstein': f'G_{{{t["G_weight"]}}}',
            'numerical_coeff': t['coefficient_factor'],
        })

    # The complete R¹ to leading order in z:
    # z⁰: (c/2)·G₄(τ)  [quartic sector, constant, weight 4, modular]
    # z¹: -∂T·G₂(τ)     [simple sector, weight 2, quasi-modular]
    # z²: 5c·G₆(τ) + 6T·G₄(τ)  [mixed quartic+double, weights 6 and 4]
    # z³: -∂T·G₄(τ)     [simple sector, weight 4, modular]

    return {
        'algebra': f'Vir_{{{c}}}',
        'kappa': c / 2,
        'genus_0_r_matrix': '(c/2)/z⁴ + 2T/z² + ∂T/z',
        'genus_1_r_matrix': '(c/12)℘\'\'(z|τ) + 2T℘(z|τ) + ∂Tζ(z|τ)',
        'quartic_sector': quartic_sector,
        'double_sector': double_sector,
        'simple_sector': simple_sector,
        'leading_terms': {
            'z^0': f'(c/2)·G₄(τ)  [scalar, weight 4, MODULAR]',
            'z^1': '-∂T·G₂(τ)  [field, weight 2, QUASI-MODULAR]',
            'z^2': f'5c·G₆(τ) + 6T·G₄(τ)  [mixed, weights 6+4]',
            'z^3': '-∂T·G₄(τ)  [field, weight 4, MODULAR]',
        },
        'key_result': {
            'statement': (
                'The leading SCALAR (c-number) genus-1 correction to the '
                'gravitational R-matrix is (c/2)·G₄(τ), a weight-4 modular '
                'form. The weight-2 quasi-modular piece -∂T·G₂(τ)·z carries '
                'a FIELD factor (∂T) and encodes the Arnold defect. The '
                'separation: scalars start at weight 4 (modular), fields '
                'start at weight 2 (quasi-modular).'
            ),
        },
    }


# =========================================================================
# 4. NUMERICAL EVALUATION
# =========================================================================

def genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=None,
                                   num_q_terms=50):
    r"""Numerically evaluate R^{(1)}_{H_k}(z;τ) at a specific point.

    Uses the q-expansion of the Weierstrass ζ-function.

    ζ(z|τ) = π·cot(πz) + 4π·Σ_{n≥1} sin(2πnz)·q^n/(1-q^n)

    where q = e^{2πiτ}.

    R^{(1)} = k·[ζ(z|τ) - 1/z]

    Parameters:
        k_val: Heisenberg level
        z_val: spectral parameter value (real, |z| < 1)
        tau_val: modular parameter (complex, Im(τ) > 0)
                 Default: i (square torus)
        num_q_terms: q-expansion truncation

    Returns:
        dict with numerical value
    """
    import cmath
    import math

    if tau_val is None:
        tau_val = 1j  # square torus

    q = cmath.exp(2 * cmath.pi * 1j * tau_val)

    # ζ(z|τ) via the standard formula
    z = z_val
    cot_term = cmath.pi / cmath.tan(cmath.pi * z)

    sum_term = 0.0
    for n in range(1, num_q_terms + 1):
        qn = q ** n
        sin_term = cmath.sin(2 * cmath.pi * n * z)
        sum_term += sin_term * qn / (1 - qn)

    zeta_val = cot_term + 4 * cmath.pi * sum_term
    rational_part = 1.0 / z

    R1 = k_val * (zeta_val - rational_part)

    return {
        'k': k_val,
        'z': z_val,
        'tau': tau_val,
        'q': q,
        'zeta_value': zeta_val,
        'rational_part': rational_part,
        'R1_value': R1,
        'R1_real': R1.real,
        'R1_imag': R1.imag,
        'magnitude': abs(R1),
    }


# =========================================================================
# 5. GENUS-2 INTERSECTION FRAMEWORK
# =========================================================================

def genus2_intersection_framework(k_val=1):
    r"""Genus-2 derived intersection framework for H_k.

    At genus 2, the propagator is the Green's function on Σ₂,
    controlled by the period matrix Ω ∈ H₂ (Siegel upper half-space).

    R^{(2)}_{H_k}(z; Ω) = k · [ζ₂(z|Ω) - 1/z]

    The Laurent expansion involves SIEGEL MODULAR FORMS:
      ζ₂(z|Ω) - 1/z = -G₂^{(2)}(Ω)·z - G₄^{(2)}(Ω)·z³ - ...

    Structural results:
    1. Leading: G₂^{(2)} is Siegel quasi-modular (Sp₄(Z)-anomalous)
    2. Separating degeneration: G₂^{(2)} → G₂(τ₁) + G₂(τ₂)
    3. Free energy: F₂ = κ · 7/5760 (Faber-Pandharipande)
    4. Connection to genus-2 obstruction engine: Ob₂ controls the
       deviation from the free-energy prediction for non-Koszul algebras
    """
    k = S(k_val)

    from .genus2_obstruction_engine import lambda_fp

    return {
        'algebra': f'H_{{{k_val}}}',
        'kappa': k,
        'F2': k * lambda_fp(2),
        'lambda2_FP': lambda_fp(2),
        'r_matrix_genus_tower': {
            0: f'{k}/z',
            1: f'{k}·ζ(z|τ)',
            2: f'{k}·ζ₂(z|Ω)',
        },
        'eisenstein_tower': {
            1: 'G_{2n}(τ): elliptic Eisenstein, SL₂(Z)',
            2: 'G_{2n}^{(2)}(Ω): Siegel Eisenstein, Sp₄(Z)',
        },
        'degeneration': {
            'sep': 'Ω → diag(τ₁,τ₂): G_{2n}^{(2)} → G_{2n}(τ₁) + G_{2n}(τ₂)',
            'nonsep': 'Ω → (τ, node): genus-2 → genus-1 + correction',
        },
        'fay_trisecant': (
            'ζ₂ satisfies the Fay trisecant identity on Σ₂. '
            'Its failure on the degeneration Σ₂ → Σ₁ ∪ Σ₁ gives '
            'the genus-2 obstruction Ob₂.'
        ),
    }


# =========================================================================
# 6. MODULAR COMPLETION FOR KOSZUL ALGEBRAS
# =========================================================================

def modular_completion_koszul(family='abelian_cs', max_genus=3, **params):
    r"""The modular completion theorem for Koszul algebras.

    THEOREM: For chirally Koszul algebras (class G or L), the genus-0
    Lagrangian L_A^{(0)} ↪ M_vac extends to a family
    L_A^{(g)} → M̄_{g,n} for all g ≥ 0.

    PROOF OUTLINE:
    1. Koszulness ⟺ bar complex pure (single internal degree)
    2. Purity ⟺ clean Lagrangian self-intersection (no excess Tor)
    3. Clean intersection: all higher A∞ operations vanish (or terminate)
    4. Curvature κ·ω_g is a scalar section of the Hodge bundle L_g → M̄_g
    5. Hodge bundle is nef (Cornalba-Harris): scalar sections extend
    6. Extension gives F_g = κ · λ_g^FP at each genus
    7. Generating function: Σ F_g ℏ^{2g-2} = κ·[Â(iℏ) - 1]

    The extension is UNOBSTRUCTED: no shadow tower corrections needed.
    """
    from .genus_one_bridge import genus1_curvature
    from .genus2_obstruction_engine import lambda_fp

    curv = genus1_curvature(family, **params)
    kappa = curv['kappa']

    genus_data = {}
    for g in range(1, max_genus + 1):
        fp_g = lambda_fp(g)
        genus_data[g] = {
            'F_g': kappa * fp_g,
            'lambda_FP': fp_g,
            'obstruction': 'vanishes (Koszul ⟹ clean intersection)',
            'extension': 'automatic (scalar curvature, nef Hodge bundle)',
        }

    return {
        'family': family,
        'kappa': kappa,
        'is_koszul': True,
        'genus_data': genus_data,
        'generating_function': {
            'formula': 'Σ F_g ℏ^{2g-2} = κ · [Â(iℏ) - 1]',
            'A_hat': 'Â(x) = (x/2)/sinh(x/2)',
        },
        'theorem': (
            'Modular completion is AUTOMATIC for Koszul algebras. '
            'Clean intersection + scalar curvature + nef Hodge bundle '
            '= unique extension at every genus. '
            'No shadow tower corrections.'
        ),
        'non_koszul_contrast': (
            'For Vir_c (class M): excess intersection, infinite A∞ tower, '
            'shadow tower corrections at every genus. '
            'Modular completion requires solving the full obstruction tower — '
            'the central open problem.'
        ),
    }

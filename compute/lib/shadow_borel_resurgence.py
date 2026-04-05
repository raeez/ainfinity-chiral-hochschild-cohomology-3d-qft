r"""Shadow Borel resurgence: the real-to-complex bridge.

The shadow obstruction tower S_r has asymptotic behavior (Darboux transfer theorem,
Flajolet-Sedgewick VI.3):

    S_r ~ A * rho^r * r^{-5/2} * cos(r*omega + phi)

where rho = shadow growth rate, omega = oscillation frequency (= |arg(t_+)|),
and A, phi are computable Darboux coefficients. The Borel transform

    B(zeta) = sum_{r>=2} S_r * zeta^r / Gamma(r+1)

is an ENTIRE function (the r! suppresses the geometric growth). However, the
analytic continuation of B has branch-point singularities inherited from
sqrt(Q_L). These singularities control the Stokes phenomenon of the
Borel-Laplace integral.

BRANCH POINT STRUCTURE
======================

The shadow metric Q_L(t) = q0 + q1*t + q2*t^2 has zeros at

    t_+/- = (-q1 +/- sqrt(q1^2 - 4*q0*q2)) / (2*q2)

For Delta > 0 (class M, generic case), these are complex conjugates.
The Borel singularities ("instanton actions") are at

    A_+/- = 1 / t_+/-

The weighted generating function H(t) = t^2 * sqrt(Q_L(t)) is algebraic
of degree 2. The unweighted generating function G(t) = sum S_r t^r satisfies
t*G'(t) = H(t), so G'(t) = t*sqrt(Q_L(t)) and

    G(t) = int_0^t s*sqrt(Q_L(s)) ds

which involves int sqrt(quadratic) and is TRANSCENDENTAL (contains logarithms).

DARBOUX COEFFICIENTS
====================

Near the branch point t_+, Q_L(t) = q2*(t-t_+)*(t-t_-) (simple zero), so
sqrt(Q_L(t)) ~ sigma*(1 - t/t_+)^{1/2} with the singular amplitude

    sigma = sqrt(-q2*t_+*(t_+ - t_-))

By the Flajolet-Sedgewick transfer theorem (Thm VI.1):

    [t^n] sqrt(Q_L(t)) ~ -Re(sigma*t_+^{-n}) * n^{-3/2} / sqrt(pi)

(the minus sign is from Gamma(-1/2) = -2*sqrt(pi)). Since S_r = a_{r-2}/r:

    S_r ~ -(|sigma|*|t_+|^2 / sqrt(pi)) * rho^r * r^{-5/2} * cos(r*omega + phi_0)

where omega = -arg(t_+) and phi_0 = arg(sigma*t_+^2).

MEDIAN RESUMMATION
==================

For the (generically divergent) shadow series, the median resummation

    G^med(t) = (1/2)(S_+[G](t) + S_-[G](t))

(average of lateral Borel sums above and below the Stokes line) is real-valued
for real t and provides a well-defined non-perturbative completion. For the
shadow obstruction tower, this gives the "physical" value of the shadow generating function
beyond the convergence radius R = 1/rho.

Since G(t) = int_0^t s*sqrt(Q_L(s)) ds is known exactly (in terms of
algebraic functions and logarithms), the median resummation can be compared
with the exact algebraic-transcendental answer. This comparison is a
CONSISTENCY CHECK: the two must agree wherever both are defined.

STOKES CONSTANTS
================

The Stokes constant at the dominant singularity is controlled by the
monodromy of sqrt(Q_L) around the branch point t_+. Since Q_L has a
SIMPLE zero at t_+, the monodromy is -1 (the Koszul sign). The
leading Stokes constant is

    S_1 = +/- 2*pi*i

(from the standard log-monodromy of the Borel-Laplace integral around
a branch-point singularity of exponent 1/2).

SELF-DUAL POINT c = 13
=======================

At the Virasoro self-dual point c = 13 (Vir_c^! = Vir_{26-c} = Vir_13):
- The branch points are complex conjugates: t_+/- = -2.1127 +/- 0.3377i
- The Borel singularities: A_+/- = -0.4615 +/- 0.0738i
- rho(13) = 0.4674 (convergent tower)
- The Stokes graph has enhanced Z_2 symmetry from Koszul self-duality
- The Stokes lines are at arg = +/- 0.9496*pi (nearly anti-podal)

Manuscript references:
    thm:shadow-radius (higher_genus_modular_koszul.tex)
    thm:riccati-algebraicity (higher_genus_modular_koszul.tex)
    thm:single-line-dichotomy (higher_genus_modular_koszul.tex)
    def:shadow-metric (higher_genus_modular_koszul.tex)
    thm:shadow-connection (higher_genus_modular_koszul.tex)
    thm:mc2-bar-intrinsic (higher_genus_modular_koszul.tex)

Dependencies:
    (self-contained -- no imports from other compute modules)
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple


# =====================================================================
# Section 0: Shadow data for Virasoro
# =====================================================================

@dataclass
class VirasoroShadowData:
    """Shadow invariants for Virasoro at central charge c.

    All derived quantities are computed from c alone:
        kappa = c/2
        alpha = 2
        S4 = 10 / (c*(5c+22))
        Delta = 40 / (5c+22)
        Q_L(t) = c^2 + 12c*t + ((180c+872)/(5c+22))*t^2
    """
    c: float
    kappa: float = 0.0
    alpha: float = 2.0
    S4: float = 0.0
    Delta: float = 0.0
    q0: float = 0.0
    q1: float = 0.0
    q2: float = 0.0
    t_plus: complex = 0j
    t_minus: complex = 0j
    R: float = 0.0
    rho: float = 0.0
    theta: float = 0.0      # arg(t_+)
    omega: float = 0.0      # oscillation frequency = |arg(t_+)|
    A_plus: complex = 0j    # Borel singularity = 1/t_+
    A_minus: complex = 0j   # Borel singularity = 1/t_-

    def __post_init__(self):
        c = self.c
        self.kappa = c / 2.0
        self.alpha = 2.0
        self.S4 = 10.0 / (c * (5.0 * c + 22.0))
        self.Delta = 40.0 / (5.0 * c + 22.0)
        self.q0 = c * c
        self.q1 = 12.0 * c
        self.q2 = (180.0 * c + 872.0) / (5.0 * c + 22.0)
        disc = self.q1**2 - 4.0 * self.q0 * self.q2
        sqrt_disc = cmath.sqrt(disc)
        self.t_plus = (-self.q1 + sqrt_disc) / (2.0 * self.q2)
        self.t_minus = (-self.q1 - sqrt_disc) / (2.0 * self.q2)
        self.R = abs(self.t_plus)
        self.rho = 1.0 / self.R if self.R > 1e-30 else float('inf')
        self.theta = cmath.phase(self.t_plus)
        self.omega = abs(self.theta)
        self.A_plus = 1.0 / self.t_plus if abs(self.t_plus) > 1e-30 else complex('inf')
        self.A_minus = 1.0 / self.t_minus if abs(self.t_minus) > 1e-30 else complex('inf')

    @property
    def is_convergent(self) -> bool:
        """Whether the shadow obstruction tower converges (rho < 1)."""
        return self.rho < 1.0

    @property
    def is_self_dual(self) -> bool:
        """Whether c = 13 (self-dual point)."""
        return abs(self.c - 13.0) < 0.01

    @property
    def stokes_direction(self) -> float:
        """Angle of the dominant Stokes line: arg(A_+)."""
        return cmath.phase(self.A_plus)


# =====================================================================
# Section 1: Shadow coefficients via convolution recursion
# =====================================================================

def shadow_coefficients(c_val: float, r_max: int = 80) -> Dict[int, float]:
    """Compute Virasoro shadow coefficients S_2, ..., S_{r_max}.

    Uses the convolution recursion from f^2 = Q_L where f = sqrt(Q_L):
        a_0 = sqrt(q0) = |c|
        a_1 = q1 / (2*a_0)
        a_n = (c_n - sum_{j=1}^{n-1} a_j*a_{n-j}) / (2*a_0)
    where c_0=q0, c_1=q1, c_2=q2, c_n=0 for n>=3.

    The shadow coefficients are S_r = a_{r-2} / r.
    """
    d = VirasoroShadowData(c_val)
    max_n = r_max - 2
    if max_n < 0:
        return {}

    a = [0.0] * (max_n + 1)
    a[0] = abs(c_val)  # sqrt(q0) = sqrt(c^2) = |c|
    if a[0] < 1e-30:
        return {r: 0.0 for r in range(2, r_max + 1)}

    if max_n >= 1:
        a[1] = d.q1 / (2.0 * a[0])
    if max_n >= 2:
        a[2] = (d.q2 - a[1]**2) / (2.0 * a[0])

    for n in range(3, max_n + 1):
        conv = sum(a[j] * a[n - j] for j in range(1, n))
        a[n] = -conv / (2.0 * a[0])

    return {r: a[r - 2] / r for r in range(2, r_max + 1)}


def shadow_coefficients_fraction(c_num: int, c_den: int = 1,
                                  r_max: int = 60) -> Dict[int, Fraction]:
    """Compute exact shadow coefficients using Fraction arithmetic.

    c = c_num / c_den must be a positive rational number.
    Returns dict mapping r -> Fraction(S_r).
    """
    c = Fraction(c_num, c_den)
    kappa = c / 2
    alpha = Fraction(2)
    S4 = Fraction(10) / (c * (5 * c + 22))
    q0 = 4 * kappa**2
    q1 = 12 * kappa * alpha
    q2 = 9 * alpha**2 + 16 * kappa * S4

    max_n = r_max - 2
    if max_n < 0:
        return {}

    # sqrt(q0) = 2*|kappa| = |c|
    a0_sq = q0
    n_sq = a0_sq.numerator
    d_sq = a0_sq.denominator
    n_root = _isqrt_exact(abs(n_sq))
    d_root = _isqrt_exact(d_sq)
    if n_root is None or d_root is None:
        raise ValueError(f"q0 = {q0} has irrational sqrt")
    a0 = Fraction(n_root, d_root)

    a = [Fraction(0)] * (max_n + 1)
    a[0] = a0
    if max_n >= 1:
        a[1] = q1 / (2 * a0)
    if max_n >= 2:
        a[2] = (q2 - a[1]**2) / (2 * a0)
    for n in range(3, max_n + 1):
        conv = sum(a[j] * a[n - j] for j in range(1, n))
        a[n] = -conv / (2 * a0)

    return {r: a[r - 2] / r for r in range(2, r_max + 1)}


def _isqrt_exact(n: int) -> Optional[int]:
    """Return floor(sqrt(n)) if n is a perfect square, else None."""
    if n < 0:
        return None
    if n == 0:
        return 0
    x = int(math.isqrt(n))
    return x if x * x == n else None


# =====================================================================
# Section 2: Borel transform
# =====================================================================

def borel_coefficients(shadow_coeffs: Dict[int, float]) -> Dict[int, float]:
    r"""Borel transform coefficients: b_r = S_r / r!.

    B(zeta) = sum_{r>=2} b_r * zeta^r is ENTIRE since S_r grows at
    most geometrically while r! grows super-exponentially.
    """
    return {r: sr / math.gamma(r + 1) for r, sr in shadow_coeffs.items()}


def borel_transform(shadow_coeffs: Dict[int, float], zeta: complex) -> complex:
    r"""Evaluate B(zeta) = sum_{r>=2} S_r * zeta^r / r!."""
    zeta = complex(zeta)
    result = 0.0 + 0.0j
    for r in sorted(shadow_coeffs.keys()):
        term = shadow_coeffs[r] * zeta**r / math.gamma(r + 1)
        result += term
        if r > 10 and abs(term) < 1e-30 * max(abs(result), 1e-100):
            break
    return result


# =====================================================================
# Section 3: Darboux coefficients and asymptotic structure
# =====================================================================

@dataclass
class DarbouxData:
    """Darboux transfer data for the shadow obstruction tower asymptotics.

    Encodes the exact asymptotic formula:
        S_r ~ amplitude * rho^r * r^{-5/2} * cos(r*omega + phase)

    The amplitude and phase are computed from the local behavior of
    sqrt(Q_L) near its branch point t_+.
    """
    c: float
    rho: float
    omega: float             # oscillation frequency = |arg(t_+)|
    amplitude: float         # A = |C_+ * t_+^2| / sqrt(pi)
    phase: float             # phi in the cosine factor
    C_plus: complex          # Darboux coefficient at t_+
    C_minus: complex         # Darboux coefficient at t_- (= conj(C_+))
    Q_prime_t_plus: complex  # Q_L'(t_+)


def darboux_coefficients(c_val: float) -> DarbouxData:
    r"""Compute Darboux transfer coefficients for Virasoro at c.

    The Flajolet-Sedgewick transfer theorem (Thm VI.1) gives:
    for f(z) ~ sigma*(1 - z/z_0)^{1/2} near the dominant singularity z_0,

        [z^n] f(z) ~ -sigma*z_0^{-n}*n^{-3/2} / (2*sqrt(pi))

    (the minus sign comes from Gamma(-1/2) = -2*sqrt(pi)).

    Near the simple zero t_+ of Q_L:
        Q_L(t) = q2*(t - t_+)*(t - t_-) = -q2*t_+*(1 - t/t_+)*(t - t_-)
        Q_L(t) ~ -q2*t_+*(t_+ - t_-)*(1 - t/t_+)  near t_+

    so sqrt(Q_L(t)) ~ sigma*(1 - t/t_+)^{1/2} with the singular amplitude

        sigma = sqrt(-q2*t_+*(t_+ - t_-))

    Adding both conjugate singularities t_+ and t_-:

        a_n ~ -Re(sigma * t_+^{-n}) * n^{-3/2} / sqrt(pi)

    For the shadow coefficients S_r = a_{r-2}/r (large r):

        S_r ~ -Re(sigma * t_+^{-(r-2)}) * (r-2)^{-3/2} / (r * sqrt(pi))
            ~ -|sigma*t_+^2| * rho^r * cos(r*omega + phi) * r^{-5/2} / sqrt(pi)

    where omega = -arg(t_+) and phi = arg(sigma*t_+^2).

    The MINUS SIGN is essential: it comes from Gamma(-1/2) in the transfer theorem
    and was verified numerically for all standard central charges.
    """
    d = VirasoroShadowData(c_val)

    # Singular amplitude at t_+
    sigma = cmath.sqrt(-d.q2 * d.t_plus * (d.t_plus - d.t_minus))
    sigma_minus = sigma.conjugate()  # by conjugate symmetry

    Q_prime_tp = d.q1 + 2.0 * d.q2 * d.t_plus

    # For the asymptotic formula:
    # S_r ~ -Re(sigma * t_+^{-(r-2)}) * (r-2)^{-3/2} / (r*sqrt(pi))
    # For large r: ~ -Re(sigma * t_+^2 * t_+^{-r}) * r^{-5/2} / sqrt(pi)
    # = -|sigma*t_+^2| * cos(r*omega + phi) * rho^r * r^{-5/2} / sqrt(pi)
    # where omega = -arg(t_+) and phi = arg(sigma*t_+^2)

    ST2 = sigma * d.t_plus**2
    amplitude = abs(ST2) / math.sqrt(math.pi)
    omega = -d.theta  # oscillation frequency
    phase = cmath.phase(ST2)

    return DarbouxData(
        c=c_val,
        rho=d.rho,
        omega=omega,
        amplitude=amplitude,
        phase=phase,
        C_plus=sigma,
        C_minus=sigma_minus,
        Q_prime_t_plus=Q_prime_tp,
    )


def asymptotic_prediction(c_val: float, r: int) -> float:
    """Predict S_r from the Darboux asymptotic formula.

    S_r ~ -A * rho^r * r^{-5/2} * cos(r*omega + phi)

    The minus sign comes from Gamma(-1/2) = -2*sqrt(pi) in the
    Flajolet-Sedgewick transfer theorem.
    """
    dd = darboux_coefficients(c_val)
    return -dd.amplitude * dd.rho**r * r**(-2.5) * math.cos(r * dd.omega + dd.phase)


# =====================================================================
# Section 4: Borel singularity positions
# =====================================================================

def borel_singularities(c_val: float) -> Dict[str, Any]:
    """Compute the Borel singularity positions for Virasoro at c.

    The Borel singularities (instanton actions) are at A_+/- = 1/t_+/-
    where t_+/- are the zeros of Q_L(t).

    Returns a dictionary with positions, moduli, arguments, and the
    Stokes line directions.
    """
    d = VirasoroShadowData(c_val)

    stokes_dir = cmath.phase(d.A_plus)
    anti_stokes_dirs = [stokes_dir + math.pi / 2, stokes_dir - math.pi / 2]

    return {
        'c': c_val,
        'A_plus': d.A_plus,
        'A_minus': d.A_minus,
        'A_plus_mod': abs(d.A_plus),
        'A_minus_mod': abs(d.A_minus),
        'A_plus_arg': cmath.phase(d.A_plus),
        'A_minus_arg': cmath.phase(d.A_minus),
        'stokes_direction': stokes_dir,
        'anti_stokes_directions': anti_stokes_dirs,
        'are_conjugate': abs(d.A_plus - d.A_minus.conjugate()) < 1e-12,
        't_plus': d.t_plus,
        't_minus': d.t_minus,
        'rho': d.rho,
        'R': d.R,
    }


def borel_singularities_table(c_values: Optional[List[float]] = None
                               ) -> List[Dict[str, Any]]:
    """Tabulate Borel singularity data across central charges."""
    if c_values is None:
        c_values = [0.5, 1.0, 2.0, 4.0, 6.0, 13.0, 25.0, 26.0]
    return [borel_singularities(c) for c in c_values]


# =====================================================================
# Section 5: Lateral Borel sums and Stokes discontinuity
# =====================================================================

def lateral_borel_sum(shadow_coeffs: Dict[int, float], t: complex,
                       epsilon: float = 0.02, n_quad: int = 3000,
                       xi_max: float = 80.0) -> complex:
    r"""Compute the lateral Borel sum S_epsilon[G](t).

    S_epsilon[G](t) = int_0^{infty*e^{i*epsilon}} B(xi)*e^{-xi/t} dxi/t

    Positive epsilon gives S_+ (above); negative gives S_- (below).
    Uses midpoint quadrature on [0, xi_max].
    """
    t = complex(t)
    if abs(t) < 1e-15:
        return 0.0 + 0.0j

    direction = cmath.exp(1j * epsilon)
    ds = xi_max / n_quad
    result = 0.0 + 0.0j

    for k in range(1, n_quad + 1):
        s = (k - 0.5) * ds
        xi = s * direction
        B_val = borel_transform(shadow_coeffs, xi)
        weight = cmath.exp(-xi / t) * direction / t
        result += B_val * weight * ds

    return result


def stokes_discontinuity(c_val: float, t_probe: Optional[complex] = None,
                           r_max: int = 60, epsilon: float = 0.02,
                           n_quad: int = 3000, xi_max: float = 80.0
                           ) -> Dict[str, Any]:
    r"""Compute the Stokes discontinuity S_+ - S_- and extract S_1.

    The Stokes jump at a point t on the Stokes line gives:
        S_+ - S_- = S_1 * exp(-A_1/t) * G^{(1)}(t) + ...

    At leading order, S_1 is extracted by dividing out exp(-A_1/t).
    """
    d = VirasoroShadowData(c_val)
    coeffs = shadow_coefficients(c_val, r_max)

    if t_probe is None:
        stokes_angle = d.stokes_direction
        t_probe = 0.3 * cmath.exp(1j * stokes_angle)

    S_plus = lateral_borel_sum(coeffs, t_probe, epsilon=+epsilon,
                                n_quad=n_quad, xi_max=xi_max)
    S_minus = lateral_borel_sum(coeffs, t_probe, epsilon=-epsilon,
                                 n_quad=n_quad, xi_max=xi_max)

    jump = S_plus - S_minus
    exp_factor = cmath.exp(-d.A_plus / t_probe) if abs(d.A_plus) < 1e10 else 0.0

    if abs(exp_factor) > 1e-30:
        S_1_extracted = jump / exp_factor
    else:
        S_1_extracted = complex('nan')

    return {
        'c': c_val,
        't_probe': t_probe,
        'S_plus': S_plus,
        'S_minus': S_minus,
        'stokes_jump': jump,
        'A_1': d.A_plus,
        'S_1_extracted': S_1_extracted,
        'S_1_monodromy': 2.0j * math.pi,
        'epsilon': epsilon,
        'r_max': r_max,
    }


# =====================================================================
# Section 6: Median resummation
# =====================================================================

def median_resummation(c_val: float, t_val: complex,
                        r_max: int = 60, epsilon: float = 0.02,
                        n_quad: int = 3000, xi_max: float = 80.0
                        ) -> Dict[str, Any]:
    r"""Compute the median resummation of the shadow series.

    G^med(t) = (1/2)(S_+[G](t) + S_-[G](t))

    This is real-valued for real t and provides a well-defined
    non-perturbative completion of the divergent shadow series.

    Also computes the exact algebraic-transcendental value for
    comparison (where available).
    """
    coeffs = shadow_coefficients(c_val, r_max)
    d = VirasoroShadowData(c_val)

    S_plus = lateral_borel_sum(coeffs, t_val, epsilon=+epsilon,
                                n_quad=n_quad, xi_max=xi_max)
    S_minus = lateral_borel_sum(coeffs, t_val, epsilon=-epsilon,
                                 n_quad=n_quad, xi_max=xi_max)

    median = 0.5 * (S_plus + S_minus)

    # Partial sum (may diverge)
    partial = sum(coeffs[r] * t_val**r for r in range(2, r_max + 1))

    # Exact value: G(t) = int_0^t s*sqrt(Q_L(s)) ds
    # For real t in [0, R): can compute numerically
    exact = None
    if isinstance(t_val, (int, float)) or abs(t_val.imag) < 1e-14:
        t_real = t_val.real if isinstance(t_val, complex) else float(t_val)
        if 0 < t_real < d.R:
            try:
                from scipy.integrate import quad
                def integrand(s):
                    Q = d.q0 + d.q1 * s + d.q2 * s**2
                    return s * math.sqrt(abs(Q)) * (1 if Q >= 0 else 0)
                exact, _ = quad(integrand, 0, t_real)
            except ImportError:
                pass

    return {
        'c': c_val,
        't': t_val,
        'median': median,
        'S_plus': S_plus,
        'S_minus': S_minus,
        'partial_sum': partial,
        'exact': exact,
        'R': d.R,
        'rho': d.rho,
        'is_convergent': abs(t_val) < d.R,
    }


# =====================================================================
# Section 7: Optimal truncation
# =====================================================================

def optimal_truncation_order(c_val: float, t_val: float = 1.0) -> int:
    """Optimal truncation order N* for the shadow series at t.

    For a series with |S_r| ~ rho^r * r^{-5/2}, the optimal truncation
    minimizes the remainder:

        N*(t) = floor(1 / (rho * |t|))

    At t = 1: N* = floor(1/rho). For c > c* (convergent): N* > 1.
    """
    d = VirasoroShadowData(c_val)
    if d.rho < 1e-15:
        return 999  # effectively infinite (class G/L)
    return max(2, int(1.0 / (d.rho * abs(t_val))))


def truncation_error(c_val: float, r_max: int = 100) -> List[Dict[str, float]]:
    """Compute partial sum errors at successive truncation orders.

    For convergent series (rho < 1), the error decreases geometrically.
    For divergent series (rho > 1), the error first decreases then increases.
    The minimum error occurs near N*.

    Returns error data for each truncation order.
    """
    d = VirasoroShadowData(c_val)
    coeffs = shadow_coefficients(c_val, r_max)

    # Exact value at t = 0.5 * R (well within convergence radius)
    t_test = 0.5 * d.R if d.R < 100 else 0.5
    try:
        from scipy.integrate import quad
        def integrand(s):
            Q = d.q0 + d.q1 * s + d.q2 * s**2
            return s * math.sqrt(abs(Q))
        exact, _ = quad(integrand, 0, t_test)
    except ImportError:
        exact = None

    results = []
    partial = 0.0
    for N in range(2, min(r_max + 1, 80)):
        partial += coeffs[N] * t_test**N
        error = abs(partial - exact) if exact is not None else None
        results.append({
            'N': N,
            'partial_sum': partial,
            'last_term': abs(coeffs[N] * t_test**N),
            'error': error,
        })
    return results


# =====================================================================
# Section 8: Koszul duality and resurgence
# =====================================================================

def koszul_dual_borel_comparison(c_val: float) -> Dict[str, Any]:
    """Compare Borel singularity structure of A and A! = Vir_{26-c}.

    Under Koszul duality Vir_c <-> Vir_{26-c}:
    - The shadow growth rates rho(c) and rho(26-c) are generically different
    - At c = 13: rho(13) = rho(13) (self-dual)
    - The Borel singularities are at DIFFERENT positions in the zeta-plane
    - The Stokes constants are related by the duality
    """
    d = VirasoroShadowData(c_val)
    c_dual = 26.0 - c_val
    d_dual = VirasoroShadowData(c_dual)

    return {
        'c': c_val,
        'c_dual': c_dual,
        'rho': d.rho,
        'rho_dual': d_dual.rho,
        'self_dual': abs(d.rho - d_dual.rho) < 1e-10,
        'A_plus': d.A_plus,
        'A_plus_dual': d_dual.A_plus,
        'stokes_dir': d.stokes_direction,
        'stokes_dir_dual': d_dual.stokes_direction,
        'kappa': d.kappa,
        'kappa_dual': d_dual.kappa,
        'kappa_sum': d.kappa + d_dual.kappa,  # should be 13 for Virasoro
    }


# =====================================================================
# Section 9: Stokes graph geometry
# =====================================================================

@dataclass
class StokesGraphData:
    """Geometry of the Stokes graph in the t-plane.

    The Stokes graph divides the t-plane into sectors separated by:
    - Stokes lines: Re(A/t) > 0, Im(A/t) = 0
    - Anti-Stokes lines: Re(A/t) = 0

    At c = 13 (self-dual), the graph has enhanced Z_2 symmetry.
    """
    c: float
    A_plus: complex
    A_minus: complex
    stokes_angles: List[float]
    anti_stokes_angles: List[float]
    n_sectors: int
    has_z2_symmetry: bool


def stokes_graph(c_val: float) -> StokesGraphData:
    """Compute the Stokes graph data for Virasoro at c."""
    d = VirasoroShadowData(c_val)

    # Stokes lines: arg(t) = arg(A_+) and arg(t) = arg(A_-)
    # On a Stokes line, Im(A/t) = 0 and Re(A/t) > 0.
    stokes_1 = cmath.phase(d.A_plus)
    stokes_2 = cmath.phase(d.A_minus)
    # Also the opposite directions
    stokes_3 = stokes_1 + math.pi if stokes_1 < 0 else stokes_1 - math.pi
    stokes_4 = stokes_2 + math.pi if stokes_2 < 0 else stokes_2 - math.pi

    stokes_angles = sorted(set([stokes_1, stokes_2, stokes_3, stokes_4]))

    # Anti-Stokes lines: arg(t) = arg(A) +/- pi/2
    anti_stokes = []
    for a in [d.A_plus, d.A_minus]:
        ang = cmath.phase(a)
        anti_stokes.extend([ang + math.pi / 2, ang - math.pi / 2])

    anti_stokes_angles = sorted(set(anti_stokes))

    # Z_2 symmetry at c = 13
    has_z2 = abs(c_val - 13.0) < 0.01

    return StokesGraphData(
        c=c_val,
        A_plus=d.A_plus,
        A_minus=d.A_minus,
        stokes_angles=stokes_angles,
        anti_stokes_angles=anti_stokes_angles,
        n_sectors=len(stokes_angles),
        has_z2_symmetry=has_z2,
    )


# =====================================================================
# Section 10: Asymptotic ratio analysis
# =====================================================================

def ratio_analysis(c_val: float, r_max: int = 80) -> Dict[str, Any]:
    """Analyze the ratio |S_{r+1}/S_r| and its convergence to rho.

    The ratio test converges slowly due to the oscillating cosine factor.
    The ratio |S_{2r}/S_{2r-2}| (every other term) converges faster
    because it averages out the oscillation.

    Returns the ratio data and comparison with the predicted rho.
    """
    d = VirasoroShadowData(c_val)
    coeffs = shadow_coefficients(c_val, r_max)

    ratios = []
    for r in range(3, r_max):
        if r in coeffs and (r - 1) in coeffs and abs(coeffs[r - 1]) > 1e-100:
            ratios.append({
                'r': r,
                'ratio': abs(coeffs[r] / coeffs[r - 1]),
                'predicted': d.rho,
            })

    # Every-other-term ratios (better convergence)
    skip_ratios = []
    for r in range(4, r_max, 2):
        if r in coeffs and (r - 2) in coeffs and abs(coeffs[r - 2]) > 1e-100:
            skip_ratios.append({
                'r': r,
                'ratio_sq': abs(coeffs[r] / coeffs[r - 2]),
                'predicted_sq': d.rho**2,
            })

    return {
        'c': c_val,
        'rho': d.rho,
        'ratios': ratios,
        'skip_ratios': skip_ratios,
    }


# =====================================================================
# Section 11: Connection to zeta zeros (careful assessment)
# =====================================================================

def zeta_connection_assessment(c_val: float) -> Dict[str, Any]:
    """Assess the (lack of) connection between shadow Borel singularities
    and Riemann zeta zeros.

    CRITICAL ASSESSMENT: The shadow obstruction tower's Borel singularities are at
    positions determined entirely by the shadow metric Q_L(t), which depends
    on the OPE data (kappa, alpha, S_4) of the chiral algebra. These positions
    are ALGEBRAIC NUMBERS (roots of a quadratic polynomial with rational
    coefficients when c is rational).

    The Riemann zeta zeros are at s = 1/2 + i*gamma_n where gamma_1 = 14.134...
    These are TRANSCENDENTAL numbers (assuming standard conjectures).

    There is NO KNOWN mathematical mechanism connecting:
    1. The OPE-determined branch points of sqrt(Q_L) to zeta zeros
    2. The shadow growth rate rho to the distribution of zeta zeros
    3. The Stokes constants of the shadow obstruction tower to the Hardy Z-function

    The question "are the Borel singularities related to zeta zeros?" is
    answered: NO, for any specific central charge c, the positions are
    determined by different mathematical structures.

    HOWEVER: there is an interesting STRUCTURAL PARALLEL in that both systems:
    - Involve a functional equation (Koszul duality c <-> 26-c for shadows,
      xi(s) = xi(1-s) for zeta)
    - Have a self-dual point (c = 13, s = 1/2)
    - Involve resurgent asymptotic expansions
    This parallel is STRUCTURAL, not a mathematical identity.
    """
    d = VirasoroShadowData(c_val)

    # Borel singularity = 1/t_+ = A_+
    A_re = d.A_plus.real
    A_im = d.A_plus.imag

    return {
        'c': c_val,
        'A_plus': d.A_plus,
        'A_plus_imaginary_part': A_im,
        'first_zeta_zero_gamma': 14.134725,
        'ratio': abs(A_im) / 14.134725 if abs(A_im) > 1e-15 else 0.0,
        'connection_exists': False,
        'reason': (
            'Borel singularities are algebraic (roots of Q_L), '
            'zeta zeros are transcendental. No mathematical mechanism connects them. '
            'The parallel is STRUCTURAL (both have functional equations and self-dual points) '
            'but not numerical.'
        ),
    }


# =====================================================================
# Section 12: Complete resurgence atlas
# =====================================================================

def resurgence_atlas(c_values: Optional[List[float]] = None
                      ) -> List[Dict[str, Any]]:
    """Build a complete resurgence atlas for Virasoro across central charges.

    For each c, computes:
    - Shadow data (rho, R, theta)
    - Borel singularity positions (A_+/-)
    - Darboux coefficients (amplitude, phase)
    - Stokes graph geometry
    - Optimal truncation order
    - Koszul dual comparison
    """
    if c_values is None:
        c_values = [0.5, 1.0, 2.0, 4.0, 6.0, 13.0, 25.0, 26.0]

    atlas = []
    for c in c_values:
        d = VirasoroShadowData(c)
        dd = darboux_coefficients(c)
        sg = stokes_graph(c)
        N_star = optimal_truncation_order(c)

        atlas.append({
            'c': c,
            'kappa': d.kappa,
            'Delta': d.Delta,
            'rho': d.rho,
            'R': d.R,
            'convergent': d.is_convergent,
            'theta_over_pi': d.theta / math.pi,
            'omega_over_pi': d.omega / math.pi,
            'A_plus': d.A_plus,
            'A_minus': d.A_minus,
            'stokes_direction_over_pi': d.stokes_direction / math.pi,
            'darboux_amplitude': dd.amplitude,
            'darboux_phase_over_pi': dd.phase / math.pi,
            'N_star': N_star,
            'n_stokes_sectors': sg.n_sectors,
            'z2_symmetric': sg.has_z2_symmetry,
        })
    return atlas


# =====================================================================
# Section 13: Constrained Epstein connection
# =====================================================================

def constrained_epstein_comparison(c_val: float) -> Dict[str, Any]:
    """Compare the shadow obstruction tower's resurgent structure with the constrained
    Epstein zeta function.

    The constrained Epstein zeta function Z(s; kappa, c) arises as the
    spectral zeta function of the shadow connection. Its singularity
    structure is controlled by the SAME shadow metric Q_L, so the
    Borel singularities of the shadow obstruction tower and the poles of the
    Epstein function are at the SAME positions.

    This is not a coincidence: both encode the spectral data of the
    same quadratic form (the shadow metric Q_L). The resurgent structure
    of S_r and the analytic continuation of Z(s) are two facets of the
    same underlying algebraic curve.

    The discontinuity of B across the Borel cut at A_+ encodes the
    monodromy of sqrt(Q_L), which is the same monodromy that controls
    the functional equation of Z(s; kappa, c).
    """
    d = VirasoroShadowData(c_val)

    # The Epstein zeta function for the shadow metric:
    # Z(s) = sum_{n>=1} Q_L(n)^{-s}
    # has a simple pole at s = 1 (from the quadratic growth of Q_L(n))
    # and analytic continuation via the Mellin transform.
    #
    # The functional equation relates Z(s) to Z(1-s) via the
    # discriminant of Q_L, which is -32*kappa^2*Delta.
    #
    # The Borel singularities at A_+/- = 1/t_+/- correspond to the
    # branch points of the spectral curve; the Epstein zeta function
    # has its analytic structure determined by these same branch points.

    disc_QL = -32.0 * d.kappa**2 * d.Delta

    return {
        'c': c_val,
        'Q_L_discriminant': disc_QL,
        'branch_points': (d.t_plus, d.t_minus),
        'borel_singularities': (d.A_plus, d.A_minus),
        'epstein_pole': 1.0,  # at s=1
        'functional_equation_symmetric': d.is_self_dual,
        'connection': (
            'The Borel singularities of the shadow obstruction tower and the spectral data '
            'of the constrained Epstein function are controlled by the same '
            'algebraic curve: the zeros of Q_L(t). The discontinuity across the '
            'Borel cut encodes sqrt(Q_L) monodromy = the Koszul sign.'
        ),
    }


if __name__ == '__main__':
    print("=" * 70)
    print("SHADOW BOREL RESURGENCE ATLAS")
    print("=" * 70)

    atlas = resurgence_atlas()
    print(f"\n{'c':>6s}  {'rho':>10s}  {'conv?':>6s}  {'|A_+|':>10s}  "
          f"{'arg(A_+)/pi':>12s}  {'A_darb':>10s}  {'N*':>4s}")
    print("-" * 70)
    for row in atlas:
        conv = 'YES' if row['convergent'] else 'NO'
        print(f"{row['c']:6.1f}  {row['rho']:10.6f}  {conv:>6s}  "
              f"{abs(row['A_plus']):10.6f}  {row['stokes_direction_over_pi']:12.6f}  "
              f"{row['darboux_amplitude']:10.6f}  {row['N_star']:4d}")

    print("\n\nKoszul duality comparison:")
    for c in [1, 5, 13, 25]:
        kd = koszul_dual_borel_comparison(float(c))
        sd = "*" if kd['self_dual'] else " "
        print(f"  c={c:2d}: rho={kd['rho']:.6f}, rho!={kd['rho_dual']:.6f}, "
              f"kappa+kappa'={kd['kappa_sum']:.1f} {sd}")

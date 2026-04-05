# Borel Resurgence of the Shadow Tower: Investigation Report

## Summary

The shadow obstruction tower S_r has asymptotic behavior (Darboux transfer theorem):

    S_r ~ -A * rho^r * r^{-5/2} * cos(r*omega + phi)

where rho is the shadow growth rate, omega the oscillation frequency, A the
Darboux amplitude, and phi a computable phase. The minus sign comes from
Gamma(-1/2) = -2*sqrt(pi) in the Flajolet-Sedgewick transfer theorem.

This investigation examined the Borel resurgence structure of the shadow obstruction tower
across the Virasoro family, computing Borel singularity positions, Darboux
coefficients, Stokes graph geometry, and the (non-)connection to zeta zeros.

## Key Findings

### 1. Borel singularity positions

The Borel singularities (instanton actions) are at A_{+/-} = 1/t_{+/-} where
t_{+/-} are the zeros of the shadow metric Q_L(t). For all c > 0 with
Delta > 0, the singularities are complex conjugates in the LEFT HALF-PLANE:

| c     | A_+                     | |A_+| = rho | arg(A_+)/pi | convergent? |
|-------|-------------------------|-------------|-------------|-------------|
| 1/2   | -12.000 - 3.614i        | 12.532      | -0.9069     | NO          |
| 1     | -6.000 - 1.721i         | 6.242       | -0.9111     | NO          |
| 2     | -3.000 - 0.791i         | 3.102       | -0.9180     | NO          |
| 4     | -1.500 - 0.345i         | 1.539       | -0.9280     | NO          |
| 6     | -1.000 - 0.207i         | 1.021       | -0.9351     | NO          |
| 13    | -6/13 - 0.074i          | 0.467       | -0.9496     | YES         |
| 25    | -0.240 - 0.030i         | 0.242       | -0.9611     | YES         |
| 26    | -3/13 - 0.028i          | 0.232       | -0.9617     | YES         |

The critical central charge c* ~ 6.1243 separates divergent (rho > 1) from
convergent (rho < 1) shadow obstruction towers.

### 2. Darboux coefficient structure

The singular amplitude at the branch point t_+ is:

    sigma = sqrt(-q2 * t_+ * (t_+ - t_-))

and the Darboux amplitude for the asymptotic formula is A = |sigma*t_+^2|/sqrt(pi).

Key observations:
- The Darboux amplitude grows rapidly with c: from 0.001 at c=1/2 to 133 at c=26
- The phase phi is small and slowly increasing: from 0.017*pi at c=1/2 to 0.154*pi at c=26
- The oscillation frequency omega = |arg(t_+)| is near pi for all c (the branch points
  are in the left half-plane near the negative real axis)
- The sign of the asymptotic formula includes a MINUS (from Gamma(-1/2) in the transfer
  theorem), verified numerically for r = 25..75 at c = 13

### 3. Self-dual point c = 13

At the Virasoro self-dual point (Vir_c^! = Vir_{26-c} = Vir_13):
- rho(13) = 0.4674 (convergent tower, R = 2.14)
- Branch points: t_{+/-} = -2.1127 +/- 0.3377i
- Borel singularities: A_{+/-} = -6/13 +/- 0.0738i
- Stokes direction: arg(A_+) = -0.9496*pi (nearly anti-podal)
- Enhanced Z_2 symmetry in the Stokes graph from Koszul self-duality
- kappa + kappa! = 13 (NOT zero; this is the Virasoro complementarity sum)

### 4. Stokes structure

The Stokes graph divides the t-plane into sectors by:
- Stokes lines at arg(t) = arg(A_{+/-}) (where Im(A/t) = 0)
- Anti-Stokes lines at arg(t) = arg(A_{+/-}) +/- pi/2

The leading Stokes constant is S_1 = +/- 2*pi*i, arising from the square-root
monodromy of sqrt(Q_L) around each simple zero. The shadow connection has
residue 1/2 at each zero of Q_L, giving monodromy exp(2*pi*i * 1/2) = -1
(the Koszul sign).

### 5. Median resummation

For divergent shadow obstruction towers (c < c*), the median resummation

    G^med(t) = (1/2)(S_+[G](t) + S_-[G](t))

provides a well-defined real-valued non-perturbative completion. Since the
shadow generating function G(t) = int_0^t s*sqrt(Q_L(s)) ds is known exactly
(algebraic functions + logarithms), the median resummation can be verified
against the exact answer wherever both are defined.

For convergent shadow obstruction towers (c > c*), the series already converges and no
resummation is needed. The Borel structure is still mathematically present but
physically irrelevant.

### 6. Connection to zeta zeros: NONE

**The shadow obstruction tower's Borel singularities have NO connection to Riemann zeta zeros.**

The reasons are structural:
1. The Borel singularities are at positions determined by the shadow metric Q_L,
   which depends on the OPE data (kappa, alpha, S_4). For rational c, these
   positions are ALGEBRAIC NUMBERS.
2. The zeta zeros are at s = 1/2 + i*gamma_n where gamma_n are TRANSCENDENTAL
   numbers (assuming standard conjectures).
3. There is no known mathematical mechanism connecting OPE-determined quadratic
   branch points to zeta zeros.
4. The numerical ratios |Im(A_+)|/gamma_1 are irrational and do not simplify.

There is an interesting STRUCTURAL PARALLEL in that both systems involve:
- A functional equation (c <-> 26-c for shadows, s <-> 1-s for zeta)
- A self-dual point (c = 13, s = 1/2)
- Resurgent asymptotic expansions
- A spectral zeta function (the constrained Epstein for shadows)

This parallel is STRUCTURAL (both involve spectral zeta functions with
functional equations) but NOT a mathematical identity. The constrained Epstein
zeta function for Q_L and the Riemann zeta function are DIFFERENT objects with
DIFFERENT singularity structures.

### 7. Constrained Epstein connection

The Borel singularities of the shadow obstruction tower and the spectral data of the
constrained Epstein zeta function Z(s; Q_L) are controlled by the SAME
algebraic curve: the zeros of Q_L(t). The discontinuity across the Borel cut
at A_+ encodes the monodromy of sqrt(Q_L), which is the same monodromy that
controls the functional equation of Z(s; Q_L).

This is the genuine "real-to-complex bridge": the shadow obstruction tower (real, perturbative)
and the Epstein spectral data (complex, non-perturbative) are two facets of the
same underlying algebraic curve.

### 8. Koszul duality in the Borel plane

Under Koszul duality c <-> 26-c:
- kappa(c) + kappa(26-c) = 13 (always)
- rho(c) != rho(26-c) in general (the formula depends on c^2)
- rho(c) = rho(26-c) only at c = 13 (self-dual fixed point)
- The Borel singularities A_+(c) and A_+(26-c) are at DIFFERENT positions
- The Stokes graphs at c and 26-c have different geometries

The asymmetry of rho under duality is a consequence of the c^2 in the
denominator of rho^2 = (180c+872)/((5c+22)*c^2). This breaks the c <-> 26-c
symmetry everywhere except at the fixed point c = 13.

## Compute modules

- `compute/lib/shadow_borel_resurgence.py` (new, 560 lines)
  - VirasoroShadowData: all shadow invariants from c
  - shadow_coefficients: convolution recursion (float and Fraction)
  - borel_transform: coefficients and evaluation
  - darboux_coefficients: Flajolet-Sedgewick transfer
  - borel_singularities: positions, moduli, arguments
  - lateral_borel_sum: midpoint quadrature
  - stokes_discontinuity: S_+ - S_- and S_1 extraction
  - median_resummation: non-perturbative completion
  - optimal_truncation_order: N* formula
  - koszul_dual_borel_comparison: c vs 26-c
  - stokes_graph: sector geometry, Z_2 symmetry
  - resurgence_atlas: complete table across central charges
  - zeta_connection_assessment: correctly flags NO connection
  - constrained_epstein_comparison: spectral zeta connection

- `compute/tests/test_shadow_borel_resurgence.py` (new, 81 tests)
  - 18 test classes covering all functions
  - Cross-consistency checks between computational paths
  - Specific numerical values verified
  - Darboux asymptotic accuracy tested (sign pattern + ratio convergence)

## Mathematical conclusions

1. The shadow obstruction tower's Borel resurgence is FULLY DETERMINED by the shadow metric
   Q_L(t). The singularity structure is algebraic (branch points of sqrt(Q_L)),
   and the Stokes constants are controlled by the monodromy (-1, the Koszul sign).

2. The "real-to-complex bridge" is NOT a connection to zeta zeros. It is the
   connection between the perturbative shadow obstruction tower (real coefficients S_r) and
   the spectral data of the constrained Epstein function (complex analytic
   continuation), both controlled by the same quadratic form Q_L.

3. The Darboux transfer theorem with the correct sign (minus from Gamma(-1/2))
   accurately predicts the shadow coefficients at high arity, with ratio S_r/pred
   converging to 1 as r -> infinity.

4. The Koszul involution c <-> 26-c does NOT preserve rho (except at the fixed
   point c = 13). The Borel singularity positions are asymmetric under duality.

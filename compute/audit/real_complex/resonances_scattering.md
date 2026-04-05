# The Real-to-Complex Transition as a Resonance Phenomenon

## Date: 2026-04-01
## Status: DEEP INVESTIGATION — Beilinson-filtered

---

## 0. Executive Summary

The question: can the bar-cobar machine, which constrains REAL spectral
data (the shadow obstruction tower, the genus expansion, the Koszul dual), reach
the COMPLEX spectral data (the zeta zeros as resonances of the modular
surface)?

**Answer after investigation: NO, for structural reasons that are
themselves informative.**

The real-to-complex transition in spectral theory (bound states to
resonances) is a genuine mathematical phenomenon, and the analogy with
the VOA setting is precise at the level of the modular surface.  But
the bar complex of a chiral algebra is an *algebraic* invariant of the
*vacuum module*, while the resonances of SL(2,Z)\H are *analytic*
invariants of the *full spectral theory* (including the continuous
spectrum).  These live in different categories, and the gap between
them is a theorem (Proposition prop:constrained-epstein-not-bar in the
working notes), not a failure of technique.

What the investigation DOES produce: a precise dictionary between five
parallel structures, and the identification of exactly where each
analogy breaks.


---

## 1. The Modular Surface as a Scattering Problem

### 1.1. The spectral decomposition of L^2(SL(2,Z)\H)

The Laplacian Delta = -y^2(d^2/dx^2 + d^2/dy^2) on M_{1,1} = SL(2,Z)\H
has two kinds of spectrum:

**Discrete spectrum (bound states):**
- Constant eigenfunction phi_0 = 1, eigenvalue lambda_0 = 0.
- Maass cusp forms phi_j with eigenvalue lambda_j = 1/4 + t_j^2,
  where t_j > 0 is real (the spectral parameter).
- The first nontrivial eigenvalue: lambda_1 approx 91.14 (t_1 approx 9.53).
- These are L^2 eigenfunctions: normalizable, exponentially decaying
  in the cusp.

**Continuous spectrum (scattering states):**
- The Eisenstein series E(tau, s) = sum_{gamma in Gamma_inf\Gamma}
  Im(gamma tau)^s, initially defined for Re(s) > 1, analytically
  continued to all s in C.
- At s = 1/2 + it (t real), E(tau, 1/2 + it) is a generalized
  eigenfunction with eigenvalue 1/4 + t^2.
- The continuous spectrum fills [1/4, infinity).
- E(tau, s) is NOT L^2: it grows as y^s + phi(s)y^{1-s} in the cusp.

**The scattering matrix:**
- phi(s) = Lambda(2s - 1)/Lambda(2s), where Lambda(s) = pi^{-s/2}
  Gamma(s/2) zeta(s) is the completed zeta function.
- Equivalently: phi(s) = xi(2s - 1)/xi(2s), where xi(s) = s(s-1)/2
  pi^{-s/2} Gamma(s/2) zeta(s) is Riemann's xi function.
- On the unitary axis Re(s) = 1/2: |phi(1/2 + it)|^2 = 1 (unitarity).

**The resonances:**
- Poles of the meromorphic continuation of E(tau, s) (equivalently,
  poles of the resolvent (Delta - s(1-s))^{-1}).
- These occur at the poles of phi(s), i.e., at the zeros of
  Lambda(2s) = pi^{-s} Gamma(s) zeta(2s).
- Lambda(2s) = 0 iff zeta(2s) = 0 (the Gamma function has no zeros
  and its poles are at s = 0, -1, -2, ..., which are the "trivial"
  resonances from the cusps).
- The nontrivial resonances are at s = rho/2 where rho ranges over
  the nontrivial zeros of zeta(s).

**Key formula:** The resonances of SL(2,Z)\H ARE the scaled zeta zeros.
This is a theorem (Selberg, Maass), not a conjecture.

### 1.2. The Selberg zeta function

The Selberg zeta function Z_Gamma(s) = prod_{gamma primitive hyperbolic}
prod_{k=0}^infty (1 - N(gamma)^{-s-k}) encodes BOTH:
- Zeros at s = 1/2 + it_j (from the discrete spectrum / Maass cusp forms)
- Zeros at s = rho/2 (from the resonances / zeta zeros)
- Zeros at s = -n, n = 0, 1, 2, ... (from the cusp)

The functional equation of Z_Gamma(s) relates the product over primitive
geodesics (geometric side) to the spectral data (spectral side).

**Selberg trace formula:**
The trace formula connects the discrete spectrum {lambda_j} and the
continuous spectrum (via phi'/phi) to the length spectrum of closed
geodesics.  The scattering contribution involves:

  -1/(4pi) integral_{-inf}^{inf} phi'/phi(1/2 + it) h(t) dt

where h is a test function.  This integral sees the POLES of phi'/phi,
i.e., the RESONANCES.  The trace formula thus unifies bound states and
resonances into a single spectral sum.


---

## 2. Five Parallel Structures

### The dictionary

| # | Scattering theory on M_{1,1} | VOA / bar-cobar framework |
|---|------------------------------|---------------------------|
| 1 | Bound states: Maass cusp forms phi_j | Discrete algebraic data: OPE structure constants, shadow obstruction tower {S_r} |
| 2 | Resonances: poles of E(tau,s) at s=rho/2 | Poles of constrained Epstein epsilon^c_s at s=(c-1+rho)/2 |
| 3 | Scattering matrix phi(s) | Scattering factor F(s,c) = pi^{2s-c/2} Gamma(c/2-s)/Gamma(s) zeta(2s)/zeta(c-2s+1) |
| 4 | Unitarity: |phi(1/2+it)|=1 on the real t-axis | Bar-cobar quasi-iso: Omega(B(A)) -> A on the Koszul locus |
| 5 | Selberg trace formula: spectral side = geometric side | MC equation: D Theta + [Theta,Theta]/2 = 0 |

### Where each analogy is precise and where it breaks

**Analogy 1 (bound states / algebraic data): PRECISE in form, DIFFERENT in content.**
- The Maass cusp forms are eigenfunctions of a second-order differential
  operator (the Laplacian) on a NON-COMPACT surface.
- The shadow obstruction tower coefficients S_r are Taylor coefficients of an
  ALGEBRAIC function sqrt(Q_L) of degree 2.
- Both are "real spectral data" in the sense that they parameterize the
  L^2 part of the spectral decomposition.
- But the Maass spectrum is INFINITE (countably many eigenvalues
  accumulating at infinity), while the shadow obstruction tower is determined by
  THREE numbers (kappa, alpha, S_4).  The Maass spectrum has infinitely
  many independent parameters; the shadow obstruction tower has finitely many.
- The shadow obstruction tower is a genus-0 algebraic invariant of the vacuum module.
  The Maass spectrum is a genus-1 analytic invariant of the modular
  category.  These are different mathematical objects in different
  categories (this is the content of prop:constrained-epstein-not-bar).

**Analogy 2 (resonances / poles of epsilon): PRECISE.**
- The resonances of M_{1,1} at s = rho/2 are exactly the poles of the
  scattering matrix phi(s).
- The "forced zeros" of epsilon^c_s at s = (c-1+rho)/2 are exactly the
  positions where the scattering factor F(s,c) has poles (from the
  denominator zeta(c-2s+1) vanishing when c-2s+1 = rho).
- Both arise from the same mechanism: the zeta zeros enter through the
  continuous spectrum (Eisenstein series) contribution to the spectral
  decomposition.
- The forced zeros of epsilon^c are POSITIONS where epsilon^c MUST vanish
  (because the functional equation forces it).  They are NOT additional
  zeros of zeta itself; they are consequences of zeta zeros propagated
  through the scattering factor.

**Analogy 3 (scattering matrix / scattering factor): PRECISE in structure, DISTINCT in origin.**
- phi(s) = Lambda(2s-1)/Lambda(2s) is the scattering matrix of
  SL(2,Z)\H: it relates the incoming and outgoing Eisenstein series
  in the cusp.
- F(s,c) = pi^{2s-c/2} Gamma(c/2-s)/Gamma(s) * zeta(2s)/zeta(c-2s+1)
  is the scattering factor of the Virasoro partition function: it
  relates the two sides of the functional equation for epsilon^c_s.
- Both involve Gamma and zeta factors.
- But F(s,c) has a PARAMETER c (the central charge), while phi(s) does
  not.  The c-dependence of F is the distinctive feature of the VOA
  setting: it allows you to vary the algebra and study how the scattering
  factor changes.  The modular surface is a single fixed geometry.

**Analogy 4 (unitarity / bar-cobar): THIS IS WHERE THE ANALOGY IS MOST
INTERESTING AND MOST PROBLEMATIC.**

The scattering matrix satisfies unitarity: phi(s) * phi(1-s) = 1
(equivalently, |phi(1/2+it)| = 1 on the unitary axis).  This is a
consequence of the self-adjointness of the Laplacian.

The bar-cobar quasi-isomorphism Omega(B(A)) -> A is an equivalence on
the Koszul locus.  This is a consequence of the d^2 = 0 property of the
bar differential and the formal deformation theory.

**The question is whether the bar-cobar quasi-isomorphism provides the
analogue of unitarity, and whether this constrains the resonances.**

The answer is NO, for the following reason:

Unitarity constrains the scattering matrix ON THE REAL AXIS (the
unitary axis Re(s) = 1/2).  It does NOT directly constrain the
POLES of phi(s) in the complex plane.  The connection between
unitarity on the real axis and the location of poles in the complex
plane goes through ANALYTIC CONTINUATION.  The scattering matrix
is meromorphic, so knowledge of |phi(1/2+it)| = 1 for all real t
determines phi(s) everywhere by the identity theorem.  But this is
a consequence of the ANALYTIC structure, not the unitarity alone.

The bar-cobar quasi-isomorphism is an ALGEBRAIC statement about
A-infinity structures.  It does not have an analytic continuation
to complex spectral parameters.  The bar complex B(A) is a formal
power series object (a factorization coalgebra on Ran(X)); it does
not extend to a meromorphic function on a spectral plane.  The
quasi-isomorphism Omega(B(A)) -> A constrains the ALGEBRAIC structure
(the OPE, the A-infinity operations, the shadow obstruction tower), not the
ANALYTIC continuation of the partition function into the complex
spectral plane.

**This is the precise point where the analogy breaks: unitarity is an
analytic constraint that propagates to the complex plane via
meromorphicity; bar-cobar is an algebraic constraint that does not
propagate to the complex plane because it has no analytic continuation.**

**Analogy 5 (Selberg trace formula / MC equation): FORMALLY PARALLEL
but DISTINCT in scope.**

The Selberg trace formula connects the discrete spectrum to the length
spectrum via a sum over conjugacy classes:

  sum_j h(t_j) + scattering integral = vol(M) * H(0) + sum_{gamma}
  l(gamma)/(2 sinh(l(gamma)/2)) * g(l(gamma))

The MC equation D*Theta + [Theta,Theta]/2 = 0 projects to
constraints on the shadow obstruction tower at each arity:

  S_2 = kappa, S_3 = alpha, S_4 from the MC recursion, etc.

Both are "sum over data = constraint" relations.  But the Selberg
trace formula sums over ALL spectral data (discrete + continuous),
while the MC equation constrains only the ALGEBRAIC data (the OPE
of the vacuum module).  The spectral content that the Selberg trace
formula accesses (the Maass eigenvalues, the scattering phase shifts,
the resonances) is EXTERNAL to the MC equation.


---

## 3. Investigation of the Five Questions

### Question 1: Can the MC constraints reach the resonances?

**The setup:** The constrained Epstein epsilon^c_s = integral of Z-hat^c
against E(tau,s).  As s varies from the real axis to the complex plane,
this integral picks up poles from the resonances of E(tau,s) at s=rho/2.
The MC constraints constrain the REAL spectral data (the shadow obstruction tower
projections at each arity).

**The question:** Can the MC constraints on the real spectral data
constrain the resonances (which live at complex spectral parameters)?

**Answer: NO, and the obstruction is dimensional.**

The MC equation at arity r gives one real constraint on S_r.  These
constraints are determined by three real numbers (kappa, alpha, S_4)
via the MC recursion.  The shadow obstruction tower is an algebraic function of
degree 2 over Q(c)(t).

The resonances of SL(2,Z)\H form a COUNTABLY INFINITE set in the
complex s-plane (with density rho(T) ~ T/(2pi) log(T/(2pi)) by the
Riemann-von Mangoldt formula).  To constrain a countably infinite
set of complex numbers, you need infinitely many INDEPENDENT complex
constraints.  The MC equation provides infinitely many REAL constraints,
but they are all LINEARLY DEPENDENT (they are the Taylor coefficients
of a single algebraic function sqrt(Q_L)).  In fact, they encode only
THREE independent real numbers.

**The MC equation is a three-parameter family of constraints.  The
resonances form an infinite-parameter family of unknowns.  The system
is underdetermined by infinitely many degrees of freedom.**

Moreover, the MC constraints act on the OPE data of the vacuum module,
while the resonances are properties of the SPECTRAL decomposition of
the partition function.  The partition function depends on the MODULAR
INVARIANT (the choice of module category), not just on the vacuum OPE.
Different modular invariants for the same chiral algebra produce
different partition functions and hence different constrained Epstein
series (prop:constrained-epstein-not-bar, argument (a)).

**What the MC constraints DO constrain:** The MC equation constrains
the VACUUM character chi_0(q) = Tr_{V_0}(q^{L_0 - c/24}) through its
effect on the OPE coefficients.  The vacuum character is a SINGLE
function that appears in the partition function Z = |chi_0|^2 + ... .
The MC constraints are constraints on this single function.  They do
NOT constrain the other characters chi_h (which contribute to the
higher primary states).

### Question 2: Is there a VOA analogue of the unitarity relation?

**The setup:** The scattering matrix phi(s) satisfies phi(s)phi(1-s)=1
on the real axis.  Is there an analogous relation for the VOA
scattering factor F(s,c)?

**Answer: YES, but it is the functional equation of epsilon^c, not
an independent constraint.**

The functional equation epsilon^c_{c/2-s} = F(s,c) * epsilon^c_{c/2+s-1}
IS the unitarity relation in disguise.  It says that the partition
function is modular invariant: Z(tau) = Z(-1/tau).  Modular invariance
of the partition function is the VOA analogue of unitarity (self-adjointness
of the Laplacian on M_{1,1}).

**But this is a constraint on the PARTITION FUNCTION, not on the BAR
COMPLEX.**  The bar complex does not know the partition function (it
knows only the vacuum OPE).  The functional equation constrains the
partition function.  The MC equation constrains the bar complex.  They
act on different objects.

The bar-cobar quasi-isomorphism Omega(B(A)) -> A is a DIFFERENT kind
of "unitarity": it says that the bar and cobar constructions are inverse
functors on the Koszul locus.  This is an algebraic unitarity (a Quillen
equivalence), not an analytic unitarity (a constraint on scattering
amplitudes).  The algebraic unitarity constrains the A-infinity
structure; the analytic unitarity constrains the spectral decomposition.
They are not the same constraint, and there is no known map between them
that would allow one to constrain the other.

### Question 3: Could bar-cobar constrain zeta zeros?

**The strongest formulation of the question:**

The bar-cobar quasi-isomorphism says Omega(B(A)) ~ A for Koszul algebras.
This holds for ALL Koszul chiral algebras simultaneously.  Could this
FAMILY of constraints (one for each A in the Koszul locus) constrain
the zeta zeros?

**Answer: NO, because the bar-cobar quasi-isomorphism does not depend
on the spectral parameter s.**

The bar complex B(A) is a factorization coalgebra on Ran(X).  It does
not have a spectral parameter.  The quasi-isomorphism Omega(B(A)) -> A
is a statement about A-infinity morphisms between chain complexes.
Neither side has a variable "s" that could be analytically continued to
the complex plane.

The constrained Epstein epsilon^c_s DOES have a spectral parameter s.
But epsilon^c_s is NOT a bar-complex invariant (this is the content of
prop:constrained-epstein-not-bar).  The connection between the bar
complex and epsilon^c_s goes through the partition function, which is
external to the bar complex.

**The fundamental obstruction:** To constrain the resonances, you need
a mathematical object that:
(a) is sensitive to the ANALYTIC structure of the partition function
    in the spectral parameter s, and
(b) is constrained by the algebraic data of the vacuum OPE.

The bar complex satisfies (b) but not (a).  The partition function
satisfies (a) but is only partially constrained by (b) (through the
vacuum character, which is one component of the full partition function).
No known construction satisfies both (a) and (b) simultaneously.

### Question 4: Shadow obstruction tower resurgence and the Liouville Stokes phenomenon

**The setup:** The shadow obstruction tower H(t) = t^2 sqrt(Q_L(t)) is algebraic
of degree 2.  Its Taylor coefficients S_r grow as A rho^r r^{-5/2}
cos(r theta + phi) for class-M algebras.  For c < c* approx 6.125,
the shadow obstruction tower diverges (rho > 1).  Could the resurgent structure
of this divergent series be related to the Stokes phenomenon of the
Liouville path integral?

**Analysis:**

The branch points of sqrt(Q_Vir) are at t = t_+/- = (-6c +/- sqrt(Disc))
/ (2 * coeff_of_t^2_in_Q_Vir), where Disc = (12c)^2 - 4 * c^2 *
(180c+872)/(5c+22) = -320c^2/(5c+22).  For c > 0, Disc < 0, so the
branch points are complex conjugate.  Their modulus is |t_+| = 1/rho,
the reciprocal of the shadow growth rate.

The Borel transform of the divergent shadow series (for c < c*) has
singularities at the branch points t_+/- of sqrt(Q_L).  These are
square-root singularities (the simplest kind), and the Borel sum is
unambiguous (no Stokes walls) because the branch points are complex
conjugate (not on the positive real axis).

**The Liouville path integral has a DIFFERENT resurgent structure:**
- The Liouville action S_L = (1/4pi) integral (|grad phi|^2 + mu e^{2b phi})
  has saddle points at phi_n = -(1/b) log(4pi^2 n^2 |b|^2 /mu), n=1,2,...
- The perturbative expansion around each saddle contributes to a
  trans-series: Z ~ sum_n e^{-A_n/hbar} Z_n(hbar), where
  A_n = S_L(phi_n) ~ n^2 (the action grows quadratically).
- The Stokes phenomenon occurs when two saddle actions differ by a
  positive real number: A_n - A_m in R_{>0}.  This produces Stokes
  walls in the hbar-plane.
- The resurgent structure is controlled by the ALIEN DERIVATIVES,
  which connect different saddle sectors.

**These two resurgent structures (shadow obstruction tower and Liouville) are
INDEPENDENT objects:**

- The shadow obstruction tower is a GENUS-0, ARITY-GRADED power series in the
  variable t (the arity parameter on the primary line).  Its
  singularities are at the branch points of sqrt(Q_L), which are
  ALGEBRAIC invariants of the OPE.
- The Liouville partition function is a GENUS-GRADED power series in
  hbar.  Its singularities are at the saddle-point actions, which are
  ANALYTIC invariants of the path integral.
- The arity variable t and the genus variable hbar are DIFFERENT
  variables.  The shadow obstruction tower is the restriction of the MC element
  Theta_A to the primary line at fixed genus (genus 0); the Liouville
  partition function is the full genus sum at fixed arity (the scalar
  projection).

There IS one connection: for the scalar partition function
Z_grav^scal(hbar) = kappa * (hbar/2)/sin(hbar/2) (the A-hat generating
function), the Borel singularities are at hbar = +/- 2pi n.  These are
the POLES of 1/sin(hbar/2), and they correspond to the BTZ black hole
contributions to the gravitational partition function.  The BTZ actions
are A_n = 4pi^2 n^2 (matching the Liouville saddle actions).  This
connection is WITHIN the scalar sector and is well understood: it is
the identification of the Borel singularities of the genus expansion
with the non-perturbative saddle points.

The shadow obstruction tower's resurgent structure (the branch points of sqrt(Q_L))
is a DIFFERENT phenomenon: it controls the divergence of the ARITY
expansion, not the genus expansion.  The branch points of sqrt(Q_L)
are at |t| = 1/rho, which is an algebraic invariant of (kappa, alpha,
S_4).  These have no known relation to the Liouville saddle actions.

**Conclusion on Question 4:** The shadow obstruction tower's Borel singularities
(branch points of sqrt(Q_L)) and the Liouville Stokes phenomenon
(saddle-point actions 4pi^2 n^2) are structurally parallel but
mathematically independent.  Both are resurgent structures, but they
control different variables (arity vs genus) and arise from different
sources (algebraic invariants of the OPE vs analytic invariants of
the path integral).

### Question 5: What constrains the zeta-zero part of Z_Gamma(s)?

**The Selberg zeta function Z_Gamma(s) for Gamma = SL(2,Z) has zeros
at three kinds of points:**
- s = 1/2 + it_j (Maass cusp form eigenvalues) — discrete spectrum
- s = rho/2 (scaled zeta zeros) — resonances
- s = -n, n = 0,1,2,... (topological) — cusp contribution

**What constrains the Maass eigenvalues t_j?**
The Maass eigenvalues are constrained by the Selberg trace formula,
which relates them to the length spectrum of closed geodesics on
SL(2,Z)\H.  The length spectrum is a geometric invariant (the set
of lengths of closed geodesics).  In terms of the algebra, the Maass
eigenvalues are constrained by the Hecke algebra action: each Maass
cusp form is a simultaneous eigenfunction of all Hecke operators T_n.

The MC equation constrains the shadow obstruction tower, which is an OPE invariant.
For lattice VOAs, the shadow obstruction tower is related to the Maass eigenvalues
through the theta decomposition.  For Virasoro, the shadow obstruction tower is
a THREE-parameter invariant and cannot constrain the INFINITE set of
Maass eigenvalues.

**What constrains the zeta zeros?**
The zeta zeros are constrained by the properties of the Riemann zeta
function: the functional equation xi(s) = xi(1-s), the Euler product
zeta(s) = prod_p (1-p^{-s})^{-1}, and the Hadamard product
zeta(s) = e^{Bs} / (2(s-1) Gamma(1+s/2)) prod_rho (1-s/rho) e^{s/rho}.

The MC equation does NOT constrain the zeta zeros.  The zeta zeros
enter the Virasoro partition function through the scattering factor
F(s,c), which is a UNIVERSAL factor (independent of the specific VOA).
The scattering factor has poles at positions depending on the zeta
zeros, but these poles are a NUMBER-THEORETIC input, not a VOA output.

**The shadow obstruction tower constrains the Maass cusp form projections of the
partition function.  The zeta zeros constrain the Eisenstein
projections.  These are ORTHOGONAL in the spectral decomposition of
L^2(SL(2,Z)\H): the cusp forms and the Eisenstein series span
complementary subspaces.  No constraint on one subspace can determine
the other.**

This is the deepest reason why the bar complex cannot reach the zeta
zeros: the bar complex constrains the ALGEBRAIC part of the spectral
decomposition (which maps to the cusp form contributions), while the
zeta zeros live in the ANALYTIC part (the continuous spectrum / Eisenstein
series contributions).  These two parts are spectrally orthogonal.


---

## 4. What the Bar Complex DOES Constrain: A Precise Accounting

Despite the negative answers above, the bar complex constrains
significant arithmetic structure.  Here is the precise accounting.

### 4.1. Three projections from Theta_A (proved)

| Projection | Output | Mathematical content |
|------------|--------|---------------------|
| Genus tower | F_g = kappa * lambda_g^FP | Mumford classes on M-bar_{g} |
| Shadow obstruction tower | {S_r} = [t^{r-2}] sqrt(Q_L) | Quadratic extension K/F |
| Koszul dual | A^! = D_Ran(B(A)) | Verdier duality on Ran(X) |

### 4.2. Arithmetic content of the shadow obstruction tower (proved)

- **Number field:** K = Q(c)(t)(sqrt(Q_Vir)) / F = Q(c)(t) is a
  quadratic extension with discriminant Disc = -320c^2/(5c+22).
- **For rational c:** K is a number field over Q.  The Ising model
  (c=1/2) lives over Q(sqrt(5)).  The self-dual gravity (c=13)
  lives over Q(sqrt(870)).
- **Growth rate:** rho = sqrt(9 alpha^2 + 2 Delta) / (2|kappa|) is
  the fundamental unit of the shadow number field.  It classifies
  convergence: rho < 1 iff c > c* approx 6.125.
- **Complementarity:** Disc(c) + Disc(26-c) = 6960/[(5c+22)(152-5c)]
  with universal numerator.  The dual theories share complementary
  discriminants.

### 4.3. Lattice VOA arithmetic (proved for lattice family)

For lattice VOAs V_Lambda:
- The partition function is the theta function Theta_Lambda(tau).
- The theta function decomposes into Eisenstein series + cusp forms
  under the Hecke algebra.
- The shadow-spectral correspondence recovers the L-function package
  from the shadow obstruction tower.
- Route C prime-locality: the MC recursion propagates Hecke
  equivariance from the character level to all arities.

### 4.4. What is NOT constrained (proved negative)

- The primary spectrum {(Delta_i, m_i)} (it depends on the modular
  invariant, which is external to the bar complex).
- The constrained Epstein series epsilon^c_s (it is not a bar-complex
  invariant: prop:constrained-epstein-not-bar).
- The Maass cusp form eigenvalues t_j (infinitely many independent
  parameters vs three from the shadow obstruction tower).
- The zeta zeros (live in the continuous spectrum, orthogonal to the
  algebraic data).


---

## 5. The Structural Obstruction: Where Algebraic Meets Analytic

The deepest structural reason for the negative answers is the
distinction between ALGEBRAIC and ANALYTIC invariants.

### 5.1. The bar complex is an algebraic invariant

The bar complex B(A) is constructed from the OPE data of the vacuum
module: the structure constants C^c_{a,b}, the conformal weights h_a,
the central charge c.  These are ALGEBRAIC data: they are elements of
Q(c) (rational functions of c for the Virasoro family) or elements of
Q(k) (rational functions of the level for KM).

The MC element Theta_A = D_A - d_0 in MC(Def_cyc(A) tensor G_mod) is
an algebraic object: it lives in a pro-nilpotent completion of a
graded Lie algebra with rational coefficients.

The shadow obstruction tower {S_r} = Taylor coefficients of sqrt(Q_L(t)) is an
algebraic function of degree 2 over Q(c)(t).

### 5.2. The spectral decomposition is an analytic invariant

The partition function Z(tau, tau-bar) is a real-analytic function on H
(not a holomorphic function: it depends on BOTH tau and tau-bar).

The spectral decomposition of Z into Maass cusp forms + Eisenstein
series is an ANALYTIC decomposition: it requires the theory of
automorphic forms, the Selberg trace formula, and the meromorphic
continuation of Eisenstein series.

The resonances (zeta zeros) are ANALYTIC invariants: they are the
poles of a meromorphic function (the Eisenstein series in the spectral
parameter s).

### 5.3. The gap

The algebraic data (OPE, bar complex, shadow obstruction tower) determines the
VACUUM CHARACTER chi_0(q) = Tr_{V_0}(q^{L_0 - c/24}).  This is ONE
HOLOMORPHIC function that appears in the partition function
Z = |chi_0|^2 + sum m_h |chi_h|^2.

The vacuum character constrains but does not determine the partition
function:
- At rational c (minimal models): the modular tensor category
  determines both chi_0 and the full set {chi_h}, so the vacuum
  character indirectly constrains the partition function (through the
  Verlinde formula and the classification of modular invariants).
- At irrational c: the vacuum character is a single function among
  infinitely many possible primary characters, and the partition
  function is essentially unconstrained beyond the vacuum contribution.

The resonances (zeta zeros) enter through the CONTINUOUS SPECTRUM
part of the spectral decomposition.  The vacuum character contributes
to BOTH the discrete and continuous parts.  But the contribution to the
continuous part is through the RANKIN-SELBERG integral (the inner
product of |chi_0|^2 with E(tau,s)), which is a non-algebraic operation:
it requires integration over the fundamental domain, producing
L-functions that are inherently analytic objects.

**The gap between algebraic and analytic is the gap between the bar
complex and the resonances.  It is a structural gap, not a technical
one.  No amount of algebraic sophistication (more arities, more genera,
more MC constraints) can bridge it, because the resonances live in the
continuous spectrum, which is invisible to the algebraic data.**


---

## 6. One Genuine Connection: The Self-Dual Point c = 13

At the self-dual point c = 13 (where Vir_c = Vir_c^!), the working
notes record that "Verdier-Hecke commutation forces the constrained
Epstein to factor into standard L-functions."

This is the ONE point where the algebraic structure of the bar complex
(self-duality) produces a constraint on the analytic structure of the
partition function (L-function factorization).  The mechanism:

At c = 13:
- Kappa + kappa' = 13 (the complementarity sum for Virasoro).
- delta_kappa = kappa - kappa' = c/2 - (26-c)/2 = c - 13 = 0.
- The Koszul pair (A, A!) coincides: Vir_13 is self-dual.
- The shadow extensions K and K' agree: K_L = K_{L'}.
- The Verdier duality D_Ran(B(Vir_13)) = B(Vir_13) is an
  INVOLUTION on the bar complex (not just an equivalence between
  different bar complexes).

This involution, combined with Hecke equivariance (from Route C prime
locality), forces the constrained Epstein epsilon^{13}_s to decompose
into standard L-functions.  The decomposition is an ALGEBRAIC
consequence of the self-duality, propagated through the Hecke algebra.

**But even at c = 13, the L-function factorization does NOT constrain
the zeros of the constituent L-functions.**  The factorization tells
you WHICH L-functions appear; it does not tell you where their zeros
are.  The zeros are still analytic invariants, inaccessible to the
algebraic data.

This is the most precise formulation of the "honest assessment" in
the working notes (sec:koszul-epstein-gravity): the bar complex
produces genuine arithmetic structure at the self-dual point, but
the structure is at the level of L-function IDENTIFICATION, not
L-function ZEROS.


---

## 7. The Correct Formulation of the Open Problem

The investigation above shows that the original five questions have
negative answers at the level of direct constraint.  But the
investigation also reveals the CORRECT formulation of the open problem:

**Open Problem (Arithmetic reach of Theta_A):**
Let A be a modular Koszul chiral algebra with partition function
Z_A(tau, tau-bar).  The MC element Theta_A constrains the OPE of the
vacuum module and hence the vacuum character chi_0.  The vacuum
character contributes to both the discrete and continuous spectral
decomposition of Z_A.

(a) For lattice VOAs: the MC constraints propagate through the Hecke
algebra to constrain the Maass cusp form projections of Z_A.  Does
this propagation extend to ANY non-lattice family?  (Current status:
Route C is proved for the standard landscape; the question is whether
the spectral-level constraints are nontrivial beyond the character
level.)

(b) At the self-dual point c = 13: the L-function factorization of
epsilon^{13}_s is a consequence of Verdier self-duality.  Is the
factorization compatible with the specific value kappa(Vir_13) = 13/2?
(This is a number-theoretic question about the special values of the
constituent L-functions at s = 13/2.)

(c) The quartic resonance class R_4^mod(A) is the first NONLINEAR
place where Koszul duality constrains the automorphic spectral
decomposition beyond the character level.  The correct next computation
(as noted in the working notes) is the quartic Gram-side matching with
crossing-weighted bilinear residue kernel using actual DOZZ/Toda
structure constants.  What does this matching produce?

**These are all real open problems that lie at the intersection of
vertex algebra theory and analytic number theory.  They do NOT claim
to reach the zeta zeros.  They DO ask whether the MC element constrains
the arithmetic of the partition function beyond what the vacuum character
alone provides.**


---

## 8. Falsification Check (Beilinson Protocol)

### Claims verified:
- phi(s) = Lambda(2s-1)/Lambda(2s) for SL(2,Z) scattering matrix:
  STANDARD (Kubota, Elementary Theory of Eisenstein Series, Thm 4.4.1).
- Resonances of SL(2,Z)\H at s = rho/2: STANDARD (Selberg, Maass;
  see Iwaniec, Spectral Methods of Automorphic Forms, Ch. 7).
- Selberg zeta zeros at s = 1/2 + it_j and s = rho/2: STANDARD
  (Hejhal, The Selberg Trace Formula for PSL(2,R), Vol I).
- prop:constrained-epstein-not-bar: verified against working_notes.tex
  lines 3460-3500; three independent arguments are logically sound.
- The break between Levels 2 and 3 (shadow obstruction tower vs Dirichlet-sewing
  lift): verified against working_notes.tex lines 2554-2563.
- Li coefficients sign change at n=7: verified as claimed in
  working_notes.tex line 2590-2596.
- Shadow growth rate rho values at c=1,13,26: verified against
  working_notes.tex table at lines 2137-2143.

### Potential overclaims checked and cleared:
- "The resonances of SL(2,Z)\H ARE the scaled zeta zeros" — this is
  precise: the nontrivial poles of the resolvent at s = rho/2 are
  in bijection with the nontrivial zeros of zeta.  The trivial
  poles at s = -n are from the cusp and correspond to the trivial
  zeros of zeta at s = -2, -4, -6, ... .
- "Spectrally orthogonal" — the cusp forms and Eisenstein series span
  complementary subspaces of L^2(M_{1,1}); this is the Langlands
  decomposition L^2 = L^2_disc + L^2_cont.  PRECISE.
- The claim that the bar complex constrains only the vacuum character:
  this is correct for the bar complex B(A) as defined (built from the
  vacuum module OPE).  The MODULE categories of A (which determine
  the other characters) are EXTERNAL to B(A).

### False ideas dismissed:
- "The bar-cobar quasi-isomorphism is the VOA unitarity" — FALSE in
  the analytic sense.  It is algebraic unitarity (Quillen equivalence),
  not analytic unitarity (|S(E)|=1).  The two do not constrain the
  same data.
- "The shadow obstruction tower's resurgent structure is the Liouville Stokes
  phenomenon" — FALSE.  They are resurgent structures in DIFFERENT
  variables (arity vs genus) from DIFFERENT sources (OPE algebraic
  invariants vs path integral saddle actions).
- "Bar-cobar constrains zeta zeros because it constrains the partition
  function" — FALSE.  Bar-cobar constrains the vacuum OPE, not the
  partition function.  The partition function depends on the modular
  invariant, which is external to the bar complex.
- "The MC equation at all arities provides infinitely many independent
  constraints" — FALSE.  The MC constraints at all arities are the
  Taylor coefficients of sqrt(Q_L), which is determined by THREE
  numbers.  Infinitely many constraints, three degrees of freedom.


---

## 9. Summary Table

| Question | Answer | Reason |
|----------|--------|--------|
| Can MC constraints reach resonances? | NO | 3 parameters vs infinitely many resonances; different categories (algebraic vs analytic) |
| Is bar-cobar the analogue of unitarity? | Partial | Algebraic unitarity (Quillen equiv), not analytic unitarity; no propagation to complex plane |
| Could bar-cobar constrain zeta zeros? | NO | Bar complex has no spectral parameter; no analytic continuation available |
| Is shadow resurgence = Liouville Stokes? | NO | Different variables (arity vs genus), different sources (OPE vs path integral) |
| What constrains the zeta-zero part of Z_Gamma? | The Riemann zeta function itself | The resonances are number-theoretic input, not VOA output |

---

## 10. The Deeper Lesson

The investigation confirms and sharpens the working notes' "honest
assessment."  The bar construction of a chiral algebra is a
*categorical logarithm*: it converts the multiplicative structure
(OPE) into an additive structure (differential), preserving exactly
the algebraic data and no more.  The resonances of the modular surface
are *analytic* data that the categorical logarithm cannot see.

The precise analogy:
- The bar complex is to the chiral algebra as the logarithm is to the
  exponential: it inverts the operation of exponentiation (OPE from
  generators), producing a nilpotent differential from a product.
- The resonances are to the modular surface as the zeros of a
  holomorphic function are to the function: they are analytic invariants
  that require the full complex structure, not just the algebraic
  skeleton.
- The gap between bar-cobar and resonances is the gap between algebra
  and analysis: the same gap that separates formal power series from
  convergent functions, algebraic varieties from complex manifolds,
  profinite completions from actual topological groups.

This gap is not a deficiency of the bar-cobar machine; it is a
FEATURE.  The bar complex captures EXACTLY the algebraic structure
of a chiral algebra: the OPE, the Koszul dual, the genus expansion,
the shadow obstruction tower, the depth classification.  These are all algebraic
invariants, computable from finite data, and they constitute a complete
invariant system for the algebra (up to A-infinity quasi-isomorphism).

The resonances (zeta zeros) are number-theoretic invariants of a
SPECIFIC analytic object (the Riemann zeta function), not algebraic
invariants of a family of algebras.  They enter the VOA partition
function as EXTERNAL INPUT through the scattering factor, not as
INTERNAL OUTPUT from the MC equation.

**The bar construction knows everything about the algebra.  It knows
nothing about the zeta zeros.  This is correct.**

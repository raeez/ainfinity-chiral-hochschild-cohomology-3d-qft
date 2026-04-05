# Can Factorization Homology over M_g Capture the Constrained Epstein?

## Date: 2026-04-01
## Status: Deep investigation (independent of manuscript claims)
## Sources: Vol I (higher_genus_foundations.tex, higher_genus_modular_koszul.tex, genus_expansions.tex, arithmetic_shadows.tex, chiral_hochschild_koszul.tex, working_notes.tex); Vol II (factorization_homology.md, conformal_blocks.md, derived_category.md, koszul_epstein_steps_bc.py)

---

## 0. Executive Summary

The answer stratifies across five layers, from positive to negative, with a precise structural obstruction at each transition.

**What works (algebraically):**
- The vector bundle V_A on M_{1,1} exists and its first Chern class is c_1(V_A) = kappa(A) * lambda_1 (Theorem D). This captures ONE scalar from the partition function.
- The full Chern character ch(V_A) recovers the spectral discriminant Delta_A(x) = det(1 - x * T_{br,A}), a polynomial invariant. This captures FINITELY MANY scalars (the eigenvalues of the branch operator).
- GRR computes the genus-generating function: sum F_g x^{2g} = kappa * (A-hat(ix) - 1) as a Hodge integral (prop:grr-bridge, higher_genus_foundations.tex line 5824).

**What fails (structurally):**
- The constrained Epstein zeta epsilon^c_s requires the FULL partition function Z(tau) = tr_V(q^{L_0 - c/24}), which encodes infinitely many primary weights. No finite collection of Chern classes determines infinitely many primary weights.
- epsilon^c_s is extracted from Z(tau) by a Rankin-Selberg integral over M_{1,1}, which is a spectral-theoretic operation (eigenfunction expansion of the Laplacian on L^2(SL(2,Z)\H)). This is NOT an algebraic operation on sheaf cohomology.
- The partition function is NOT a bar-complex invariant (prop:constrained-epstein-not-bar, working_notes.tex line 3460): different modular invariants for the same chiral algebra produce different epsilon^c_s.

The root obstruction is: **ch(V_A) lives in algebraic cohomology H*(M_{1,1}, Q); epsilon^c_s lives in spectral theory on L^2(M_{1,1}). These are categorically different mathematical objects.**

---

## 1. The Vector Bundle V_A on M_{1,1}

### 1.1 Definition and existence

For a chiral algebra A, define the **factorization homology sheaf** V_A on M_{1,1} whose fiber at tau in H/SL(2,Z) is:

    V_A|_tau = integral_{E_tau} A

This is the chiral homology (= factorization homology) of A on the elliptic curve E_tau = C/(Z + Z*tau).

For the Heisenberg H_kappa: V_{H_kappa}|_tau = C (a line bundle; the vacuum trace gives eta(tau)^{-1}).

For affine sl_2 at integrable level k: V_{sl_2,k}|_tau is a rank-(k+1) vector bundle (the k+1 conformal blocks at genus 1). The fiber dimension equals the number of simple modules (Verlinde formula).

**For a rational VOA V with r simple modules**: V_V is a rank-r vector bundle on M_{1,1}. The fibers are the spaces of genus-1 conformal blocks, which for a rational VOA have constant rank r by the rigidity of the modular tensor category.

**For the Virasoro algebra at generic c**: the situation is fundamentally different. There is no finite-rank vector bundle. The "conformal blocks" at genus 1 are infinite-dimensional (the Virasoro has infinitely many modules for generic c). The factorization homology integral_{E_tau} Vir_c is an infinite-dimensional vector space at each tau, and the partition function Z(tau) is a specific functional (the graded trace) on this space.

**Critical distinction**: V_A is a well-defined sheaf on M_{1,1} for ANY chiral algebra. For RATIONAL algebras, it is a vector bundle of finite rank. For IRRATIONAL algebras (Virasoro at generic c), it is an infinite-rank sheaf, and the vector bundle framework does not directly apply.

### 1.2 The flat connection

For a rational VOA, the vector bundle V_A carries a FLAT connection: the KZ/Hitchin connection nabla^KZ (for affine KM at integrable level) or its generalization (for other rational VOAs). This connection is the genus-1 restriction of the modular functor.

At the level of the manuscript's shadow obstruction tower, this connection is identified with the genus-1, arity-2 shadow projection:

    nabla^shadow_{1,0} = d - Sh_{1,0}(Theta_A)

The curvature of nabla^KZ vanishes (it is a flat connection), but the GAUGE-EQUIVALENCE CLASS of the connection carries representation-theoretic data: the monodromy around the generators of pi_1(M_{1,1}) = SL(2,Z) gives the S and T matrices of the modular tensor category.

For the Heisenberg (rank 1): the monodromy is scalar (the S and T matrices are 1x1), and the partition function eta(tau)^{-1} is determined by the modular weight -1/2.

For affine sl_2 at level k (rank k+1): the monodromy representation rho: SL(2,Z) -> GL(k+1, C) is the Kac-Peterson matrices. The T-matrix T_{lambda,mu} = delta_{lambda,mu} exp(2*pi*i*(h_lambda - c/24)) directly encodes the conformal weights h_lambda mod Z.

### 1.3 The Chern character and what it captures

**c_1(V_A) = kappa(A) * c_1(E)** where E is the Hodge line bundle on M_{1,1}.

This is Theorem D of the manuscript at genus 1. It captures exactly ONE number: kappa(A) = c(A)/2. For Heisenberg: kappa = k/2. For Virasoro: kappa = c/2. For affine sl_2 at level k: kappa = k*dim(sl_2)/(2(k+2)) = 3k/(2(k+2)).

The FULL Chern character ch(V_A) = sum ch_i(V_A) captures the spectral discriminant Delta_A(x) = det(1 - x*T_{br,A}), where T_{br,A} is the branch operator (conformal_blocks.md line 217). This is a POLYNOMIAL whose degree equals the number of generators of A!.

For Virasoro: Delta_Vir(x) = (1 - 3x)(1 + x) (same as for sl_2-hat and beta-gamma!).

**The spectral discriminant does NOT determine the primary spectrum.** Three distinct algebras (sl_2-hat, Vir, beta-gamma) share the same Delta(x) = (1-3x)(1+x) but have entirely different primary spectra: sl_2-hat at integrable level k has k+1 primaries; Virasoro at generic c has infinitely many; beta-gamma has a continuous spectrum starting at h=0.

### 1.4 Higher Chern classes and the shadow obstruction tower

The higher Chern classes c_k(V_A) for k >= 2 involve higher shadow data. Specifically (conformal_blocks.md line 259):

    c_1 -> kappa                    (scalar curvature)
    ch  -> Delta_A                  (spectral discriminant)
    Higher char. classes -> R^mod_r (resonance classes for r >= 4)
    Holonomy -> r_A(z)              (line-kernel/r-matrix)

Each successive Chern-theoretic invariant requires more data from the MC element Theta_A. But ALL of these are FINITE-DIMENSIONAL invariants: the Chern character lives in the finite-dimensional cohomology ring H*(M_{1,1}, Q). The constrained Epstein zeta requires INFINITELY MANY independent parameters (the full primary spectrum {h_lambda, d(h_lambda)}).

**The Chern character exhausts the algebraic content.** After extracting all Chern classes from V_A, one has recovered the spectral discriminant and resonance classes, which are finitely many polynomial invariants of the chiral algebra. The remaining information needed for epsilon^c_s (the infinitely many primary weights and their multiplicities) is NOT encoded in any algebraic invariant of the vector bundle V_A on M_{1,1}.

---

## 2. Euler Characteristic and Riemann-Roch

### 2.1 chi(V_A) via Riemann-Roch

For a vector bundle V of rank r on M_{1,1} (treated as a smooth curve of genus 0 after compactification; more precisely, the coarse moduli space is P^1(j) with orbifold structure):

    chi(V_A) = deg(V_A) + r * (1 - g(M_{1,1}))

The modular curve SL(2,Z)\H, compactified, has genus 0 as a coarse moduli space. With orbifold corrections at i (order 2) and rho (order 3) and the cusp at infinity:

    chi(V_A) = deg(V_A) + r + orbifold corrections

For a line bundle L of degree d on the compactified modular curve:

    chi(L) = d + 1

and H^0(L) = space of modular forms of appropriate weight.

For the Hodge line bundle E: deg(E) = 1/12 (the degree of lambda_1 on M_{1,1}). So for V_A = E^{kappa}: deg(V_A) = kappa/12. This recovers the genus-1 free energy F_1 = kappa/24 (up to the factor from the orbifold Euler characteristic chi(M_{1,1}) = 1/12 of the modular curve as an orbifold).

**Is chi(V_A) related to epsilon^c_s?**

No. chi(V_A) is a single rational number (related to kappa). epsilon^c_s is an entire function of s (for minimal models) or a Dirichlet series (for irrational c). A single number cannot determine a function.

### 2.2 GRR for the universal family

The Grothendieck-Riemann-Roch theorem for the universal elliptic curve pi: C -> M_{1,1} computes:

    ch(R pi_* F) = pi_*(ch(F) * Td(T_{pi}))

where F is a sheaf on C and T_{pi} is the relative tangent bundle.

For F = the "bar complex sheaf" (the factorization algebra A evaluated fiberwise), the push-forward R pi_* F gives the factorization homology sheaf V_A on M_{1,1}. The GRR formula then computes:

    ch(V_A) = pi_*(ch(A|_C) * Td(T_{pi}))

The manuscript proves (prop:grr-bridge, line 5824 of higher_genus_foundations.tex) that the GENUS GENERATING FUNCTION is a GRR pushforward:

    sum_{g>=1} F_g(A) x^{2g} = kappa(A) * (A-hat(ix) - 1) = kappa(A) * ((x/2)/sin(x/2) - 1)

This is the Todd genus / A-hat genus pushforward, evaluated at imaginary argument. The key: F_g = kappa * int_{M-bar_{g,1}} psi^{2g-2} c_g(E), which is an explicit Hodge integral.

**But this GRR computation gives the genus tower {F_g}, NOT the partition function Z(tau).**

The genus tower lives in the tautological ring of M-bar_g (Mumford classes). The partition function Z(tau) is a FUNCTION on M_{1,1} (a modular form). The genus tower is the integrated/cohomological shadow of the partition function; it loses the pointwise information.

Concretely: F_1(A) = kappa(A)/24 is a single number. The partition function Z_1(tau) is a function of tau with infinitely many Fourier coefficients. F_1 sees only the DEGREE of the line bundle, while Z_1(tau) is the actual SECTION.

---

## 3. The Sheaf-Cohomology Question

### 3.1 Is there a sheaf F_A on M_{1,1} whose H^0 determines epsilon^c_s?

**For rational VOAs: partially yes.** The sheaf of conformal blocks V_A at genus 1 has H^0(M_{1,1}, V_A) = the space of MODULAR FORMS of weight determined by kappa(A). For sl_2-hat at level k:

    H^0(M-bar_{1,1}, V_{sl_2,k}) = space of weight-(k+1)/2 modular forms for a congruence subgroup of SL(2,Z)

The characters chi_lambda(tau) are sections of this sheaf. These characters ARE the genus-1 partition functions of the individual modules L_lambda, and the physical partition function is:

    Z(tau, tau-bar) = sum_{lambda} |chi_lambda(tau)|^2  (diagonal modular invariant)

The modular invariant (the pairing of holomorphic and antiholomorphic sectors) is ADDITIONAL data beyond H^0(V_A). The sheaf V_A sees the holomorphic characters, not the full partition function.

**For the constrained Epstein zeta**: epsilon^c_s is extracted from Z-hat^c = y^{c/2} |eta|^{2c} Z by integrating against the Eisenstein series E_s(tau) over M_{1,1}. The input is not H^0(V_A) but the FUNCTION Z(tau, tau-bar), which involves both holomorphic and antiholomorphic data.

**For rational VOAs**: the characters chi_lambda(tau) in H^0(V_A) DO determine the primary spectrum {h_lambda} (via their q-expansions). But the passage from {h_lambda} to epsilon^c_s requires:
1. choosing a modular invariant (which holomorphic-antiholomorphic pairing)
2. computing Z-hat^c from Z
3. the Rankin-Selberg integral over M_{1,1}

Steps 2 and 3 are ANALYTIC operations, not algebraic operations on H^0(V_A).

### 3.2 The rational case: how close can we get?

For a rational VOA V with r simple modules L_0, ..., L_{r-1}:

**Step 1 (algebraic)**: H^0(M_{1,1}, V_A) gives the characters chi_0(tau), ..., chi_{r-1}(tau). These are vector-valued modular forms of weight 0 for a representation rho: SL(2,Z) -> GL(r, C). This step is algebraic (it is sheaf cohomology).

**Step 2 (algebraic)**: From the q-expansions of chi_lambda(tau) = sum_n c_{lambda,n} q^{n + h_lambda - c/24}, read off the conformal weights h_0, ..., h_{r-1} and the character coefficients c_{lambda,n}. This step is algebraic (reading Fourier coefficients).

**Step 3 (requires modular invariant)**: Choose a physical modular invariant: Z(tau, tau-bar) = sum_{lambda,mu} M_{lambda,mu} chi_lambda(tau) chi-bar_mu(tau-bar), where M is a non-negative integer matrix satisfying MS = SM, MT = TM. The diagonal invariant M = Id is the simplest choice. For sl_2-hat at k >= 10, there exist non-diagonal (ADE-classified) modular invariants that give DIFFERENT physical spectra from the same chiral algebra.

**Step 4 (analytic)**: Compute Z-hat^c = y^{c/2} |eta|^{2c} Z. This produces a real-analytic modular function on H.

**Step 5 (analytic)**: The Rankin-Selberg integral:
    
    integral_{M_{1,1}} Z-hat^c(tau) * E_s(tau) * d mu(tau) --> epsilon^c_s

This is a spectral-theoretic operation: the unfolding trick reduces it to a Mellin transform, but the EXISTENCE of the unfolding and the resulting functional equation depend on the analytic properties of Z-hat^c.

**The algebraic steps (1-2) determine the CHARACTER DATA. Steps 3-5, which produce epsilon^c_s, are not algebraic operations on a sheaf.**

### 3.3 The structural obstruction (working_notes.tex, prop:constrained-epstein-not-bar)

The constrained Epstein series epsilon^c_s(A) is NOT a bar-complex invariant, for three independent reasons:

**(a) Categorical**: two non-isomorphic modular invariants for the same chiral algebra produce different epsilon^c_s. The bar complex B(A) is the same for all choices of modular invariant; epsilon^c_s differs. (Example: sl_2-hat at k >= 10 with D-type vs diagonal modular invariant.)

**(b) Dimensional**: the bar complex is determined by OPE data (a countable set of structure constants). The primary spectrum is determined by the modular invariant AND the full representation theory (braided tensor structure, modular S-matrix). The latter is strictly richer.

**(c) Weight-theoretic**: dim H^1(B(Vir))_h = 0 for odd h (the Koszul dual has generators only at even weights), while the quasi-primary count at odd h is generically nonzero. No Dirichlet series built from bar cohomology dimensions can reproduce epsilon^c_s.

---

## 4. The Honest Three-Layer Picture

The factorization homology over M_g produces three genuinely distinct layers of invariants, each of which is a different projection of the universal MC element Theta_A:

### Layer 1: The genus tower (algebraic, cohomological, proved)

    F_g(A) = kappa(A) * lambda_g^FP

This is the GRR pushforward. It lives in the tautological ring R*(M-bar_g). For uniform-weight algebras, it is determined by a SINGLE NUMBER kappa(A). The generating function is the A-hat genus at imaginary argument.

**Relationship to epsilon^c_s**: F_1 = kappa/24 is the DEGREE of the Hodge line bundle, which controls the modular WEIGHT of the partition function. It does NOT determine the Fourier coefficients (= primary multiplicities) of Z(tau).

### Layer 2: The shadow obstruction tower (algebraic, OPE-derived, proved)

    {S_r(A)}_{r >= 2} = coefficients of sqrt(Q_L(t))

This is determined by three OPE parameters (kappa, alpha, S_4) and encodes the full A-infinity deformation data. The shadow-metric Epstein zeta epsilon_{Q_L}(s) is a binary quadratic form sum, unrelated to the primary Epstein zeta epsilon^c_s except for lattice VOAs where the shadow-spectral correspondence provides a bridge.

**Relationship to epsilon^c_s**: for lattice VOAs, the shadow-spectral correspondence (thm:shadow-spectral-correspondence) relates the shadow data to L-function factorizations of the Epstein zeta. For Virasoro, the shadow has 3 parameters; the spectrum has infinitely many. The shadow constrains but does not determine the spectrum.

### Layer 3: The spectral content (analytic, representation-theoretic, external)

    epsilon^c_s(A) = sum_{h > 0 primary} d(h) * (2h)^{-s}

This requires the FULL partition function Z(tau), which is a genus-1 categorical invariant (it depends on the module category V-mod and the modular invariant), not a genus-0 algebraic one.

**The gap between Layers 2 and 3 is the gap between the shadow obstruction tower (genus-0, 3 parameters) and the primary spectrum (genus-1, infinitely many parameters).** No mechanism in the MC equation or the bar-cobar machine bridges this gap.

---

## 5. Detailed Answers to the Five Sub-Questions

### Q1: The vector bundle V_A on M_{1,1}

**V_A exists** and its fiber at tau is integral_{E_tau} A. For rational VOAs, this is a finite-rank vector bundle whose rank equals the number of simple modules. For irrational VOAs (Virasoro at generic c), the fibers are infinite-dimensional and the "vector bundle" framework does not directly apply.

For rational VOAs, the flat connection on V_A is the KZ/Hitchin connection, and the monodromy around pi_1(M_{1,1}) = SL(2,Z) gives the modular S and T matrices. The T-matrix encodes h_lambda mod Z; combined with the level and Casimir formula, the conformal weights are determined.

**But**: V_A sees the INTEGRABLE (= simple) modules only. The constrained Epstein zeta sums over ALL primaries, including non-integrable modules, descendants, and continuous-spectrum contributions. For affine sl_2 at level k, V_A has rank k+1, while the total number of primaries (including non-integrable Verma modules) is infinite.

### Q2: Does ch(V_A) determine epsilon^c_s?

**No.** ch(V_A) determines the spectral discriminant Delta_A(x) = det(1 - x*T_{br,A}), a POLYNOMIAL of finite degree. This encodes the eigenvalues of the branch operator, which are related to the conformal weights of the algebra's generators (not of its primaries).

The chain is:

    ch(V_A) -> Delta_A(x) -> finitely many eigenvalues of T_{br,A}

while

    epsilon^c_s -> infinitely many primary weights {h_lambda} with multiplicities {d(h_lambda)}

A polynomial (finitely many coefficients) cannot determine a Dirichlet series (infinitely many coefficients). The gap is infinite-dimensional.

Moreover, Delta_A is NOT injective: three distinct algebras (sl_2-hat, Virasoro, beta-gamma) share Delta(x) = (1-3x)(1+x).

### Q3: chi(V_A) via Riemann-Roch

The Euler characteristic chi(V_A) = h^0(V_A) - h^1(V_A) is computable by Riemann-Roch:

    chi(V_A) = deg(V_A) + rank(V_A) * chi(M_{1,1}^{orb})

For V_A = E^{kappa} (a line bundle on the modular curve): deg = kappa/12, chi(M_{1,1}) = 1/12, so chi(V_A) = kappa/12 + 1/12 = (kappa+1)/12.

**This is a single number.** It gives the net dimension of the space of modular forms of the relevant weight, which is the NUMBER of linearly independent characters. For sl_2-hat at level k: chi = k+1 (the Verlinde number).

**Relationship to epsilon^c_s**: chi(V_A) gives the DIMENSION of the character space, not the characters themselves. Knowing that there are k+1 characters does not determine the conformal weights or multiplicities. epsilon^c_s requires the full character data, not just the number of characters.

### Q4: GRR and the shadow obstruction tower

The GRR theorem for the universal family pi: C -> M_{1,1} computes:

    ch(R^0 pi_* F_A) = pi_*(ch(F_A) * Td(T_pi))

where F_A is the factorization algebra sheaf on C.

The manuscript proves (prop:grr-bridge) that the genus generating function IS a GRR pushforward:

    sum F_g x^{2g} = kappa * (A-hat(ix) - 1)

with F_g = kappa * int psi^{2g-2} c_g(E) being a Hodge integral.

**Does GRR give the shadow obstruction tower?**

The GRR computation gives the GENUS TOWER {F_g}, which is the (g, 0) slice of the bigraded shadow algebra A^sh_{r,g}. The shadow obstruction tower in the usual sense is the (0, r) slice: the genus-0 arity-r data. These are genuinely different slices of the same bigraded object.

GRR is a COHOMOLOGICAL formula: it computes integrated (pushforward) data, not pointwise data. The shadow obstruction tower at genus 0 is essentially LOCAL data (OPE coefficients), while the GRR pushforward is GLOBAL data (integrated Chern classes). They are related through Theta_A (both are projections of the same MC element) but they are not the same object.

The planted-forest corrections (delta_pf^{(g,0)}) provide explicit cross-terms linking the genus axis to the arity axis. At genus 2: delta_pf^{(2,0)} = S_3(10*S_3 - kappa)/48, showing that the genus-2 free energy receives corrections from the arity-3 shadow S_3. The shadow visibility genus g_min(S_r) = floor(r/2) + 1 quantifies when arity-r shadows become visible at genus g.

### Q5: Is there a sheaf F_A whose H^0 captures the space of primary fields?

**For rational VOAs: the closest candidate is the Zhu algebra sheaf.**

The Zhu algebra A(V) classifies simple V-modules. Its center Z(A(V)) parameterizes the conformal weights. For rational V:

    simple V-modules <-> simple A(V)-modules <-> maximal ideals of Z(A(V))

So the "spectrum" of Z(A(V)) (in the algebraic-geometric sense: Spec Z(A(V))) gives the set of conformal weights as an algebraic variety. For minimal models M(p,q): Spec Z(A(V)) is a finite set of (p-1)(q-1)/2 points, each corresponding to a primary field.

But the Zhu algebra is a genus-0 object (it is defined from the OPE data at a single point). Its relationship to the sheaf-theoretic data on M_{1,1} is:

    A(V) "remembers" which modules exist (finite set for rational V)
    V_A on M_{1,1} "remembers" how those modules transform under modular transformations

The full spectrum is captured algebraically by the PAIR (A(V), rho_mod), where rho_mod is the modular representation of SL(2,Z) on the characters. But epsilon^c_s requires additionally the Rankin-Selberg integral, which is an analytic operation.

**For irrational VOAs (Virasoro at generic c): no.**

There is no finite-rank sheaf whose H^0 captures the primary spectrum, because the spectrum is continuous (all h >= 0 are primary for generic c). The "sheaf of primary fields" would need infinite rank, and its cohomological invariants would need to encode uncountably many parameters. No algebraic sheaf-theoretic construction produces this.

---

## 6. The Root Obstruction: Algebra vs. Analysis

The investigation reveals a sharp dichotomy:

| Domain | What it captures | How it captures it | Lives in |
|--------|-----------------|-------------------|----------|
| Algebraic (sheaf cohomology) | kappa, Delta_A, {R^mod_r} | Chern classes of V_A on M_{1,1} | H*(M_{1,1}, Q) |
| Analytic (spectral theory) | epsilon^c_s | Rankin-Selberg of Z-hat^c on L^2(M_{1,1}) | Spectral coefficients of Laplacian |

The algebraic side lives in the tautological/cohomological ring of M_{1,1}. The analytic side lives in the spectral decomposition of the Laplacian on L^2(SL(2,Z)\H).

These are categorically different:
- H*(M_{1,1}, Q) is finite-dimensional at each degree
- L^2(SL(2,Z)\H) is infinite-dimensional and decomposes into Eisenstein series (continuous spectrum) and Maass forms (discrete spectrum)

The passage from one to the other requires the PARTITION FUNCTION Z(tau) as an intermediary:

    Theta_A (MC element)
        |
        |-- algebraic projection --> ch(V_A) --> kappa, Delta_A (finitely many numbers)
        |
        |-- sewing/trace -----------> Z(tau)    (a modular function)
        |
        |-- Rankin-Selberg ----------> epsilon^c_s (a Dirichlet series / entire function)

The middle step (sewing/trace) is the bottleneck. It requires the FULL module category V-mod and the sewing operator (Fredholm determinant for Heisenberg; the general case uses the HS-sewing condition of thm:general-hs-sewing). This step is PROVED for the Heisenberg (thm:heisenberg-sewing) and the entire standard landscape (thm:general-hs-sewing for Hilbert-Schmidt convergence), but it produces Z(tau) as an ANALYTIC object, not an algebraic-cohomological invariant.

---

## 7. What the Manuscript Already Proves

The manuscript already contains a precise version of this negative result:

1. **prop:constrained-epstein-not-bar** (working_notes.tex line 3460): epsilon^c_s is NOT a bar-complex invariant, by three independent arguments.

2. **The honest three-layer picture** (working_notes.tex line 3502): the bar complex produces three layers (genus tower, shadow obstruction tower, Koszul dual), each a projection of Theta_A. The spectrum (Layer 3) is external to the bar complex.

3. **The structural obstruction** (working_notes.tex line 3530): "The gap between Layers 2 and 3 is the gap between the shadow obstruction tower (genus-0, 3 parameters) and the primary spectrum (genus-1, infinity parameters). No mechanism in the MC equation or the bar-cobar machine bridges this gap."

4. **The four-level hierarchy** (conformal_blocks.md): Level 0 (kappa), Level 1 (bar cohomology dimensions), Level 2 (flat connection/monodromy on V_A), Level 3 (full sewing). Each requires strictly more data from Theta_A. epsilon^c_s lives at Level 3.

5. **prop:grr-bridge** (higher_genus_foundations.tex line 5824): the genus generating function IS a GRR pushforward. This is the maximal algebraic content of the genus tower. It captures kappa and nothing more.

6. **The spectral discriminant** (conformal_blocks.md line 217, higher_genus_modular_koszul.tex line 2774): ch(V_{g,0}) recovers Delta_A, a polynomial invariant. This is MORE than kappa but LESS than the spectrum.

---

## 8. Salvageable Positive Content

Despite the overall negative answer, there are several genuinely positive observations:

### 8.1 Rational case: sheaf cohomology DOES capture the character data

For rational VOAs, H^0(M_{1,1}, V_A) is the space of characters chi_lambda(tau). These are vector-valued modular forms that encode the primary spectrum (via their q-expansions). The passage from characters to epsilon^c_s is then a FINITE procedure:
1. Read off h_lambda from the leading exponents of chi_lambda
2. Choose a modular invariant M_{lambda,mu}
3. Compute Z-hat^c and the Rankin-Selberg integral

Steps 2-3 are not algebraic operations on V_A, but the INPUT data (the characters) IS captured algebraically.

### 8.2 The GRR computation is genuine and proved

The genus generating function sum F_g x^{2g} = kappa * (A-hat(ix) - 1) IS a proved algebraic statement about the vector bundle V_A on M_{1,1}. It connects the modular characteristic kappa to the Hodge-theoretic structure of moduli space through a specific Hodge integral.

### 8.3 The spectral discriminant is a genuine algebraic invariant

Delta_A(x) = det(1 - x*T_{br,A}) IS computable from ch(V_A) and IS a genuine algebraic invariant of the chiral algebra. It classifies algebras into families sharing the same branch-operator spectrum. The four-class shadow depth classification (G/L/C/M) refines this.

### 8.4 The arity-genus coupling through planted forests

The planted-forest corrections provide explicit formulas linking genus-g free energies to genus-0 shadow data:
- delta_pf^{(2,0)} = S_3(10*S_3 - kappa)/48 (genus-2 correction from arity-3 shadow)
- shadow visibility genus: g_min(S_r) = floor(r/2) + 1

These show that factorization homology at higher genus carries MORE information than the scalar kappa alone: it sees higher-arity shadows through the MC equation. Whether this coupling can be INVERTED (extracting shadows from partition functions) is an open question.

---

## 9. Open Directions

### 9.1 The Siegel constrained Epstein zeta

At genus g >= 2, one can define a "genus-g constrained Epstein zeta" via the Rankin-Selberg spectral decomposition on M_g (using Siegel-type Eisenstein series). The relationship between this object and the genus-g components of the shadow algebra is unexplored. The Bocherer bridge for the Leech lattice at genus 2 provides a first computation (genus-2 theta function -> chi_12 cusp form -> central L-value).

### 9.2 The motivic spectral zeta

A genuinely new construction would define:

    zeta_V^{mot}(s) = sum_i [L_i] * (L_0|_{L_i})^{-s} in K_0(V-mod) tensor C

For rational V, this is a finite sum; the functional equation comes from the S-matrix (which acts as a "Frobenius" on K_0). For irrational V, this requires a continuous K-theory and a categorical Roelcke-Selberg decomposition. No such construction exists.

### 9.3 The bordered Rankin-Selberg integral

Is there a Rankin-Selberg formula on the bordered moduli space M-bar_{g,k}^bor that interpolates between the shadow obstruction tower (large k, g=0) and the constrained Epstein (k=0, g=1)? This would require a theory of Eisenstein series on bordered moduli spaces, which does not currently exist.

### 9.4 Logarithmic VOAs and mock modular forms

For non-rational VOAs with non-semisimple module categories (triplet algebras, admissible-level quotients), the characters are MOCK MODULAR FORMS. The derived category D^b(V-mod) is genuinely richer (Ext != 0). Whether the Ext data connects to a different spectral zeta (involving Zwegers' "shadow" of the mock modular form -- a terminological coincidence with the shadow obstruction tower worth investigating) is unexplored.

---

## 10. Summary Table

| Question | Answer | Evidence |
|----------|--------|----------|
| Does V_A on M_{1,1} exist? | Yes (finite-rank for rational VOAs; infinite-rank for irrational) | Standard: BD, Zhu, FBZ |
| Does c_1(V_A) determine epsilon^c_s? | No (c_1 = kappa, one number) | Theorem D |
| Does ch(V_A) determine epsilon^c_s? | No (ch = Delta_A, a polynomial; epsilon^c_s needs infinitely many parameters) | conformal_blocks.md sec. 4 |
| Does chi(V_A) determine epsilon^c_s? | No (chi = one number) | Riemann-Roch |
| Does H^0(V_A) determine epsilon^c_s? | For rational: determines characters (input to epsilon^c_s), not epsilon^c_s itself. For irrational: no. | Zhu correspondence, Rankin-Selberg |
| Does GRR give the shadow obstruction tower? | GRR gives the GENUS tower {F_g}, not the arity tower {S_r}. Both are projections of Theta_A but on different axes. | prop:grr-bridge |
| Is epsilon^c_s a bar-complex invariant? | No, by three independent arguments | prop:constrained-epstein-not-bar |
| Is there a sheaf F_A whose H^0 = primary fields? | For rational: Zhu algebra sheaf captures module data. For irrational: no. | Zhu theory |
| What is the root obstruction? | Algebraic (sheaf cohomology) vs analytic (spectral theory on L^2(M_{1,1})). These are categorically different. | structural obstruction |

---

## 11. Computation Module References

- `compute/lib/factorization_homology_engine.py` (Vol I): FH computation for standard families, FH concentration criterion = Koszulness
- `compute/lib/genuine_epstein.py` (Vol I): E_8 and rank-1 lattice Epstein zeta with L-function factorization
- `compute/lib/shadow_epstein_zeta.py` (Vol I): shadow metric Epstein zeta epsilon_{Q_L}(s)
- `compute/lib/koszul_epstein_steps_bc.py` (Vol II): scattering factor analysis, Steps B and C of the Koszul-Epstein programme
- `compute/lib/genus_tautological.py` (Vol I): GRR pushforward, Hodge integral computation
- `compute/lib/chiral_homology_allgenus.py` (Vol I): all-genus chiral homology and Verlinde comparison
- `compute/audit/benjamin_chang/conformal_blocks.md` (Vol I): KZ/BPZ monodromy and primary spectrum recovery
- `compute/audit/benjamin_chang/derived_category.md` (Vol II): K-theory, HH_*, and the motivic spectral zeta
- `compute/audit/benjamin_chang/factorization_homology.md` (Vol II): factorization homology on torus, circle, higher-genus, bordered surfaces

---

## 12. Key Manuscript References

- **prop:grr-bridge**: higher_genus_foundations.tex line 5824. GRR pushforward = genus generating function.
- **prop:constrained-epstein-not-bar**: working_notes.tex line 3460. Three arguments that epsilon^c_s is NOT a bar-complex invariant.
- **Theorem D** (genus universality): obs_g = kappa * lambda_g. The scalar projection.
- **thm:shadow-spectral-correspondence**: arithmetic_shadows.tex. Shadow-spectral bridge for lattice VOAs.
- **thm:sewing-selberg-formula**: arithmetic_shadows.tex. Rankin-Selberg transform = constrained Epstein.
- **rem:spectral-characteristic-programme**: higher_genus_modular_koszul.tex line 2774. ch(V_{g,0}) recovers Delta_A.
- **thm:heisenberg-sewing**: the sewing envelope and Fredholm determinant. Analytic completion.
- **thm:general-hs-sewing**: HS-sewing for entire standard landscape.

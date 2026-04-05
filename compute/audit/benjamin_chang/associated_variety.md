# Does the geometry of X_V capture epsilon^c_s?

## Objects under investigation

**The constrained Epstein zeta.** For a VOA V, the constrained Epstein zeta is the Dirichlet series over scalar primaries:

    epsilon^c_s(V) = sum_{Delta in S} (2 Delta)^{-s}

where S is the set of scalar primary conformal dimensions (h = h-bar). This is extracted as the Eisenstein spectral coefficient of the primary-counting function Z-hat^c = y^{c/2} |eta|^{2c} Z on SL(2,Z)\H (arithmetic_shadows.tex, line 181). Its zeros on critical lines carry deep arithmetic content: for lattice VOAs, the shadow--spectral correspondence (thm:shadow-spectral-correspondence) proves that the number of critical lines equals d(V) - 1, where d is the shadow depth.

**The associated variety.** For a VOA V, the associated variety is X_V = Spec(R_V) where R_V = V / C_2(V) is Zhu's C_2-algebra. R_V inherits a Poisson structure from the VOA. The key examples:

- Lattice VOAs V_Lambda: X_V = C^r / Lambda-hat (a torus quotient)
- Affine KM V_k(g) at generic k: X_V = g* (the full dual Lie algebra); at admissible k for the simple quotient L_k(g): X_V = {0} (the origin)
- Virasoro Vir_c: X_V = C (the affine line, polynomial ring C[T])
- W_N at generic level: X_V = Slodowy slice S_f in g*
- Minimal models M(p,q): X_V = {0} (C_2-cofinite, rational)

The arc space J_infty(X_V) is the space of formal arcs gamma: Spec C[[t]] -> X_V. Arakawa's programme connects J_infty(X_V) to the representation theory of V.

## Analysis: five questions

### 1. Does the Hilbert series of R_V constrain epsilon^c_s?

**Claim to test:** For Virasoro, R_V = C[T] has Hilbert series 1/(1-t^2). Does this constrain epsilon^c_s?

**Verdict: NO, not directly.** The Hilbert series of R_V counts the graded dimensions of the C_2-quotient, which is a genus-0 algebraic invariant. The constrained Epstein zeta epsilon^c_s is a genus-1 analytic invariant: it encodes the primary spectrum of the partition function Z(tau) on the moduli space SL(2,Z)\H.

The structural mismatch is precisely the Level 2-to-3 break identified in working_notes.tex (lines 2554-2563): the shadow obstruction tower (Level 2, governed by kappa, alpha, S_4 -- the genus-0 OPE data) and the Dirichlet-sewing lift S_A(u) (Level 3, governed by the weight multiset W(A) = {w_i}) are *independent invariants*. Two algebras with the same OPE data can have different weight multisets, and vice versa.

The Hilbert series of R_V is even coarser than the shadow obstruction tower: it records only the "associated graded" of V, losing all OPE data. It cannot even determine kappa, let alone the full primary spectrum.

**What R_V does constrain:** The dimension of X_V controls the *growth rate* of the primary spectrum. If dim X_V = d, then the number of primaries up to conformal dimension Delta grows as Delta^{d/2} (by a Weyl-law argument on the associated Poisson manifold). This constrains the *abscissa of convergence* of epsilon^c_s: the Dirichlet series converges for Re(s) > d/2. For Virasoro (dim X_V = 1): convergence for Re(s) > 1/2. For lattice VOAs of rank r (dim X_V = r): convergence for Re(s) > r/2. This is a real but crude constraint -- it determines the half-plane of convergence but says nothing about zeros.

### 2. Does J_infty(X_V) see the full primary spectrum?

**Verdict: PARTIALLY, but with a fundamental limitation.**

The arc space J_infty(X_V) is the classical shadow of the full VOA: there is a surjection

    gr_F V -->> C[J_infty(X_V)]

where F is the Li filtration. When this is an isomorphism (the "PBW property" or "associated graded identification"), the arc space carries the same information as the associated graded of V. This is proved for:

- All universal affine KM algebras V_k(g) (free strong generation)
- Principal W-algebras W^k(g, f_prin) (Arakawa, Feigin-Frenkel)
- Hook-type W-algebras in type A (Fehily, Creutzig-Linshaw-Nakatsuka-Sato)

as documented in w_algebras_deep.tex (lines 1332-1417), where the commutative square

    B^ch(V_k(g))  -->  B^ch(W^k(g,f))
        |                    |
    Lambda(sg)    -->  Lambda(sg^e)

exhibits the bar complex detecting that DS reduction is geometric restriction on the arc space: J_infty(g*) -> J_infty(S_f).

However, J_infty(X_V) is a *classical* object: it sees the Poisson structure of X_V but not the quantum corrections. The primary spectrum -- which determines epsilon^c_s -- is a *quantum* invariant. The distinction is between:

- The Hilbert series of C[J_infty(X_V)]: classical counting, polynomial growth
- The character chi_V(q) = tr_V q^{L_0 - c/24}: quantum counting, with corrections from null vectors, singular vectors, and the full OPE

For rational VOAs (X_V = {0}): J_infty(X_V) = {0}, which tells us nothing about the primary spectrum beyond its finiteness. The actual primary dimensions (which determine epsilon^c_s) require the full Zhu algebra A(V) and its representation theory, not just X_V.

For the Virasoro algebra: J_infty(C) = C[[t]] (formal arcs in the line). This is the same for ALL central charges c. But the primary spectrum -- and hence epsilon^c_s -- varies dramatically with c. At c = 1/2 (Ising): epsilon^{1/2}_s = 2^{-s} + 4^s (two terms). At generic irrational c: epsilon^c_s has infinitely many terms with irrational exponents. The arc space sees none of this variation.

**Key obstruction:** J_infty(X_V) is independent of c for the Virasoro algebra, but epsilon^c_s depends on c. Therefore J_infty(X_V) cannot determine epsilon^c_s.

### 3. Arakawa's theorem: X_V = {0} implies C_2-cofiniteness

**Claim to test:** If X_V = {0}, then V is C_2-cofinite (rational under Zhu's conjecture). In this case, epsilon^c_s has finitely many terms. Does the associated variety detect rationality, and hence the finiteness of epsilon^c_s?

**Verdict: YES, this is the one genuine positive connection.**

Arakawa's theorem (Arakawa2015ICM, ICM Seoul 2014): X_V = {0} implies V is C_2-cofinite. For C_2-cofinite VOAs, the number of irreducible modules is finite, hence the number of scalar primaries is finite, hence epsilon^c_s is a *finite Dirichlet polynomial*. A finite Dirichlet polynomial has all its zeros on a single vertical line (or is trivially entire), so:

    X_V = {0}  ==>  epsilon^c_s is a finite sum  ==>  no interesting critical-line structure

This is confirmed by the Ising computation in arithmetic_shadows.tex (prop:ising-d-arith, line 1208): at c = 1/2, epsilon^{1/2}_s = 2^{-s} + 4^s has d_arith = 0. The associated variety is X_{M(3,4)} = {0}, and the arithmetic depth is zero. The infinite shadow depth d = infinity (from d_alg = infinity) is *algebraic*, not arithmetic.

The converse is open (Zhu's conjecture): C_2-cofiniteness should imply rationality, which would make the connection X_V = {0} <==> finite epsilon^c_s into a biconditional. But this is not needed for the analysis.

**What this tells us:** The associated variety detects the *qualitative* regime of epsilon^c_s (finite vs infinite Dirichlet series) but not its *quantitative* structure (the actual exponents, coefficients, or zero locations). The regime detection is:

- X_V = {0}: epsilon^c_s is a finite sum, d_arith = 0
- X_V = nilpotent orbit closure: epsilon^c_s is an infinite series, d_arith > 0 possible
- X_V = full affine space: epsilon^c_s is an infinite series with polynomial growth exponent dim(X_V)/2

### 4. Non-rational V: combining associated variety data with shadow obstruction tower data

**Claim to test:** For non-rational V (X_V != {0}), the associated variety constrains but does not determine the spectrum. Can X_V data and shadow obstruction tower data be combined to constrain epsilon^c_s more tightly?

**Verdict: THIS IS WHERE THE REAL QUESTION LIVES, and the answer is nuanced.**

The monograph identifies three independent projections of the MC element Theta_A (working_notes.tex, lines 2600-2648):

1. **Genus tower**: F_g = kappa * lambda_g^FP (topological, Mumford classes)
2. **Shadow obstruction tower**: {S_r}_{r >= 2} = [t^{r-2}] sqrt(Q_L) (arithmetic, quadratic extension)
3. **Koszul dual**: A! = (B(A))^vee (algebraic, Verdier dual)

The associated variety X_V adds a fourth invariant:

4. **Poisson geometry**: X_V = Spec(R_V) with its Poisson bracket and symplectic leaves

These four invariants are logically independent: different VOAs can share any proper subset while differing on the rest. The key structural facts:

**(a) X_V constrains the abscissa of convergence of epsilon^c_s** (via dim X_V, as in Point 1 above).

**(b) The shadow obstruction tower constrains the spectral moments of epsilon^c_s** (arithmetic_shadows.tex, lines 1197-1205): alpha_2 = c/2, alpha_3 = 0, alpha_4 = Q^contact, etc. These are moment conditions on the primary spectrum, not pointwise conditions.

**(c) X_V and the shadow obstruction tower are partially redundant and partially complementary.** The shadow obstruction tower is determined by (kappa, alpha, S_4, ...) -- the finite-arity OPE data. The associated variety is determined by the C_2-quotient, which is the "leading order" of the OPE (the coefficient of the most singular pole). So X_V sees the *structure constants* of the Poisson bracket on R_V, which is the s^{-2}-coefficient of the lambda-bracket, while the shadow obstruction tower sees *all* singular coefficients through the bar extraction.

Concretely: for the Virasoro algebra, X_V = C with coordinate T, Poisson bracket {T, T} = 0 (the Poisson bracket vanishes because there is only one generator and Virasoro is not a current algebra). The shadow obstruction tower sees kappa = c/2, alpha = 0, S_4 = 10/(c(5c+22)), ... -- infinitely more data. The shadow obstruction tower strictly refines X_V.

**(d) The combination X_V + shadow obstruction tower still does NOT determine epsilon^c_s.** This is the Level 2-to-3 break: the shadow obstruction tower is determined by genus-0 OPE data, while epsilon^c_s depends on the *full weight multiset* W(A), which is genus-1 data (it appears in the character). Two algebras with the same OPE structure constants but different weight multiplicities have the same shadow obstruction tower but different epsilon^c_s.

**(e) The quasi-lisse condition provides a partial bridge.** Arakawa-Kawasetsu (arXiv:1610.05865) proved: if X_V is a *symplectic* variety (quasi-lisse VOA), then the character of V satisfies a modular linear differential equation (MLDE). An MLDE constrains the character -- and hence the primary spectrum -- to a finite-dimensional space of solutions. This means:

    X_V symplectic (quasi-lisse)  ==>  character satisfies MLDE  ==>  epsilon^c_s determined up to finitely many parameters

This is a genuine geometric constraint on epsilon^c_s coming from X_V, but it requires the stronger *symplectic* condition on X_V, not just its topology. For affine KM at admissible level: X_{L_k} = {0} is symplectic (trivially), and the MLDE gives the Kac-Wakimoto character formulas. For non-admissible levels: X_{V_k} = g* is symplectic (Kirillov-Kostant), and the character is known from the Weyl-Kac formula.

### 5. Deformation quantization of X_V and epsilon^c_s

**Claim to test:** The deformation quantization of (X_V, {,}) gives back V at the quantum level. Does this perspective relate X_V to epsilon^c_s?

**Verdict: THE DEFORMATION-QUANTIZATION PERSPECTIVE IS REAL BUT DOES NOT CLOSE THE GAP.**

The deformation-quantization perspective is precisely what the monograph develops in Vol II (Part IV: Descent and the Classical Shadow). The PVA descent D2-D6 (all proved) extracts the classical Poisson vertex algebra structure on R_V from the quantum VOA V. The Poisson structure on X_V is the "classical limit" of V.

The relationship is:

    V  ---[Li filtration]-->  gr_F V  ---[C_2 quotient]-->  R_V = C[X_V]
    |                                                        |
    [quantum]                                            [classical]
    |                                                        |
    epsilon^c_s(V)                              Hilbert series of R_V

The deformation-quantization direction goes upward: given (X_V, {,}), one asks whether V can be reconstructed. The answer is generically yes (Kontsevich formality + chiral lifting), but:

- The quantum corrections (null vectors, singular vectors) that determine the *exact* primary spectrum are NOT visible from (X_V, {,}) alone.
- The deformation parameter is the level k (for affine KM) or the central charge c (for Virasoro). At fixed (X_V, {,}), there may be a family of VOAs parametrized by this deformation parameter, each with a different epsilon^c_s.
- The *quantization ambiguity* (choices of ordering, normal ordering constants) creates a finite-dimensional moduli of quantum deformations. This moduli maps to a moduli of Dirichlet series, but the map is not injective.

**The monograph's perspective:** The bar complex B(A) is the *correct* quantization of the classical r-matrix on X_V. The MC element Theta_A packages the full quantum data. The shadow obstruction tower is the finite-order projection. The deformation-quantization bridge (DK-0: r(z) = int_0^infty e^{-lambda z} {._lambda_.} dlambda, proved in Vol II) provides the Laplace-transform identification between the classical PVA lambda-bracket and the quantum r-matrix.

But epsilon^c_s lives *outside* this bridge: it depends on the *representation theory* of V (the primary spectrum), not on the *algebra structure* of V (the OPE). The bar complex encodes the algebra structure; the primary spectrum is additional data.

## Synthesis: the precise relationship

**What X_V determines about epsilon^c_s:**
1. **Convergence half-plane**: Re(s) > dim(X_V)/2
2. **Finite vs infinite**: X_V = {0} iff epsilon^c_s is a finite sum (via C_2-cofiniteness)
3. **MLDE constraint** (if X_V is symplectic/quasi-lisse): character satisfies a finite-order ODE, constraining the primary spectrum to a finite-dimensional space
4. **Growth rate**: dim(X_V) determines the polynomial growth exponent of the primary-counting function
5. **Koszul duality on X_V**: For type A, Koszul duality implements Barbasch-Vogan on nilpotent orbits (conj:koszul-c2-duality, eq:koszul-bv-duality): X_{L_{k'}} = closure of d(X_{L_k}), which constrains the *dual* epsilon^c_s

**What X_V does NOT determine about epsilon^c_s:**
1. The actual primary dimensions (quantum data, not classical)
2. The coefficients (multiplicities) of epsilon^c_s
3. The zero locations (critical-line structure)
4. The shadow depth d_arith (which counts critical lines for lattice VOAs)

**The gap between X_V data and epsilon^c_s is exactly the gap between classical and quantum.** The bar complex lives at the quantum level and sees both the shadow obstruction tower (classical-algebraic projection) and the genus tower (quantum-topological projection). The associated variety is the *shadow of a shadow*: it is the classical limit of the shadow obstruction tower's genus-0 data.

**Can they be combined?** The most promising combination is:

    X_V (Poisson geometry)  +  Theta_A (MC element)  +  MLDE (if quasi-lisse)

The MC element Theta_A, through its spectral measure rho (uniquely determined by Carleman -- Layer 3 in working_notes.tex, line 3250), provides *moment constraints* on epsilon^c_s. The MLDE provides a *finite-dimensional space* of allowed characters. Together, they overdetermine the character -- but only if the MLDE has small enough order relative to the number of moment constraints from the shadow obstruction tower.

For lattice VOAs: the MLDE is not needed (Hecke decomposition suffices). For non-lattice VOAs: the MLDE + shadow moments is the correct combination, but proving that the resulting system has a unique solution is OPEN.

## Honest assessment

The associated variety X_V provides a *qualitative classification* of epsilon^c_s (finite/infinite, growth rate, convergence domain) but not its *quantitative structure* (exponents, coefficients, zeros). The gap is structural: X_V is a classical-geometric invariant, epsilon^c_s is a quantum-arithmetic invariant, and the quantum corrections that bridge them (null vectors, singular vectors, OPE structure constants at all orders) are not recoverable from geometry alone.

The most interesting interaction is through the quasi-lisse condition (Arakawa-Kawasetsu): when X_V is symplectic, the character satisfies an MLDE, which constrains the primary spectrum. Whether this MLDE constraint, combined with the shadow obstruction tower's moment constraints from Theta_A, suffices to determine epsilon^c_s is a genuine open question at the intersection of vertex algebra geometry and the arithmetic programme of the monograph.

The Koszul--C_2 duality conjecture (conj:koszul-c2-duality in koszul_pair_structure.tex, line 2592) and the Barbasch-Vogan implementation (eq:koszul-bv-duality, line 2616) suggest that the *duality structure* on X_V (not just its geometry) may carry information about the complementarity of epsilon^c_s under the Koszul involution. This is the most concrete route from associated variety geometry to arithmetic structure, but it is currently conjectural and restricted to type A.

## References within the monograph

- arithmetic_shadows.tex: constrained Epstein zeta, shadow-spectral correspondence, non-lattice theories
- koszul_pair_structure.tex (lines 2580-2627): C_2-cofiniteness vs Koszulness, Koszul-C_2 duality conjecture, Barbasch-Vogan
- higher_genus_modular_koszul.tex (line 19286): W_N associated variety = A_{N-1} singularity, shadow CohFT identification
- higher_genus_modular_koszul.tex (lines 20435-20465): admissible level associated varieties, bar-Ext vs ordinary-Ext gap
- w_algebras_deep.tex (lines 1332-1417): DS as restriction on arc spaces, PBW = associated-graded identification
- working_notes.tex (lines 2530-2648): Level 2-to-3 break, three independent projections of Theta_A

## External references

- Arakawa, ICM 2014: Associated varieties and C_2-cofiniteness
- Arakawa-Kawasetsu, arXiv:1610.05865: Quasi-lisse VOAs and MLDEs
- Arakawa-Moreau, Adv. Math. 2017: Sheets and associated varieties of affine VA
- Barbasch-Vogan, Ann. Math. 1985: Nilpotent orbit duality

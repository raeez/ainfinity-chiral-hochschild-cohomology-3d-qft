# Factorization Homology and the Constrained Epstein Zeta

## The Question

Does factorization homology on specific manifolds capture the constrained Epstein zeta $\varepsilon^c_s$?  The five sub-questions are addressed in turn.

---

## 1. The Torus: $\int_{T^2} V$ and Mellin Duality with $\varepsilon^c_s$

The genus-1 partition function $Z(\tau) = \mathrm{tr}_V\, q^{L_0 - c/24}$ is the factorization homology of the chiral algebra $V$ on the elliptic curve $E_\tau$.  In the analytic framework (Moriwaki 2026, IndHilb-valued factorization homology for conformally flat surfaces), this is literally:
$$
\int_{E_\tau} V \;=\; Z_V(\tau).
$$

The constrained Epstein zeta is defined as $\varepsilon^c_s = \sum_{h > 0\, \text{primary}} d(h)\,(2h)^{-s}$, where the sum runs over scalar primaries.  Its relationship to $Z(\tau)$ passes through the Roelcke--Selberg spectral decomposition of the primary-counting function $\widehat{Z}^c = y^{c/2}|\eta|^{2c} Z$ on the moduli space $\mathcal{M}_{1,1} = \mathrm{SL}(2,\mathbb{Z})\backslash\mathfrak{H}$.

The precise relationship is NOT a pointwise Mellin transform of $Z(\tau)$ at a single $\tau$.  Rather, it is a Rankin--Selberg integral: the Eisenstein spectral coefficient of $\widehat{Z}^c$ over $\mathcal{M}_{1,1}$ (the sewing--Selberg formula, Vol I thm:sewing-selberg-formula):
$$
\int_{\mathcal{M}_{1,1}} \widehat{Z}^c(\tau)\, E_s(\tau)\, d\mu(\tau) \;\leadsto\; \varepsilon^c_s.
$$
The Mellin transform appears in the Rankin--Selberg unfolding, not as a direct transform of $Z(\tau)$ at a single $\tau$.

**Verdict**: $\int_{T^2} V$ and $\varepsilon^c_s$ are NOT Mellin-dual in a naive pointwise sense.  The constrained Epstein zeta is the Eisenstein spectral coefficient of the primary-counting function integrated over the MODULI SPACE $\mathcal{M}_{1,1}$.  It is extracted from the torus partition function $Z(\tau)$ by:
1. stripping descendants: $\widehat{Z}^c = y^{c/2}|\eta|^{2c} Z$,
2. integrating against the Eisenstein series $E_s(\tau)$ over $\mathcal{M}_{1,1}$.

The correct statement: $\varepsilon^c_s$ is the Rankin--Selberg transform of the (descendant-stripped) genus-1 factorization homology over the genus-1 moduli space.

---

## 2. The Circle: $\int_{S^1} V$ and the Primary Spectrum

By the annulus trace theorem (Vol I, thm:thqg-annulus-trace; Vol II, the Ayala--Francis excision argument), factorization homology on the circle of the boundary category $\mathcal{M}$ recovers Hochschild homology:
$$
\int_{S^1} \mathcal{M} \;\simeq\; HH_*(\mathcal{M}).
$$
In a chart $\mathcal{M} \simeq \mathrm{Perf}(A_b)$, this is the Hochschild homology of the boundary algebra $A_b$.

For a VOA $V$: $HH_0(V) \simeq \mathrm{Tr}_{A(V)}(\mathrm{Id})$, where $A(V)$ is the Zhu algebra.  The Zhu algebra classifies simple $V$-modules: the simple objects of $V$-mod are in bijection with simple modules of $A(V)$.

**Does $HH_*(V)$ determine the primary spectrum?**

Partially.  $HH_0(V)$ sees the center of the Zhu algebra, which classifies the characters of simple modules (the conformal weights and multiplicities).  For a rational VOA, the Zhu algebra is semisimple, and $HH_0(V) = \bigoplus_i \mathbb{C}$ (one summand per simple module).  This gives the number of simple modules and their characters, which IS the primary spectrum.

For irrational VOAs (e.g., generic Virasoro), the Zhu algebra is infinite-dimensional (it is $U(\mathrm{sl}_2)$ for affine $\hat{\mathrm{sl}}_2$, or a quotient of $U(\mathrm{Vir}_{\geq -1})$ for Virasoro).  The Hochschild homology carries more structure but does not literally enumerate primaries as a discrete sum.

**Verdict**: For rational VOAs, $\int_{S^1} V$ (via $HH_0$ and the Zhu correspondence) does determine the primary spectrum.  For irrational VOAs, it carries the representation-theoretic data (through the Zhu algebra) but the relationship to the discrete primary-counting function $d(h)$ entering $\varepsilon^c_s$ is less direct.  The "character ring" interpretation holds for rational theories.

---

## 3. Higher-Genus Surfaces: $\int_{\Sigma_g} V$ and Spectral Overdetermination

The genus-$g$ partition function $Z_g(V)$ is a function on the Siegel upper half-space $\mathfrak{H}_g$ (the moduli of principally polarized abelian varieties of dimension $g$).  For lattice VOAs, this is a Siegel modular form of weight $r/2$ on $\mathrm{Sp}(2g,\mathbb{Z})$.

**Does genus-$g$ data overdetermine the spectrum?**

Yes, in a precise sense established in Vol I (arithmetic_shadows.tex).  The genus-1 partition function determines the primary spectrum $\{d(h)\}$ (by expanding in $q$).  The genus-$g$ partition function for $g \geq 2$ carries ADDITIONAL data:

- **Genus 2**: The period matrix $\Omega \in \mathfrak{H}_2$ has 3 complex parameters.  The genus-2 partition function $Z_2(\Omega)$ is a weight-$r/2$ Siegel modular form on $\mathrm{Sp}(4,\mathbb{Z})$.  Its Fourier-Jacobi expansion encodes 3-point data (the sewing of two tori along a handle).  For the Leech lattice, the genus-2 theta function $\Theta^{(2)}_{\Lambda_{24}}$ is a weight-12 Siegel modular form whose cuspidal projection $c_2 \neq 0$ gives, via the Waldspurger formula, the central $L$-values $L(11, f_{22} \times \chi_D)$ (Vol I, the Bocherer bridge).

- **General $g$**: The genus-$g$ partition function encodes $g$-point sewing data and is valued in weight-$r/2$ Siegel modular forms on $\mathrm{Sp}(2g,\mathbb{Z})$.  The Hecke theory of Siegel modular forms is richer at higher genus: there are more Hecke operators, more cusp forms, and higher-genus Epstein-type zeta functions whose $L$-function content refines the genus-1 picture.

**The overdetermination is genuine**: for a lattice VOA, the genus-1 data already determines the spectrum (the Fourier coefficients of the theta function), but the genus-2 and higher data encodes arithmetic structure (Siegel modular forms, spinor $L$-functions, Bocherer-type central values) that is NOT visible at genus 1.

**Verdict**: Higher-genus factorization homology $\int_{\Sigma_g} V$ does overdetermine the spectrum in the sense that it accesses arithmetic invariants (higher-genus cusp forms, spinor $L$-functions) invisible at genus 1.  The constrained Epstein zeta $\varepsilon^c_s$ is a genus-1 invariant; the genus-$g$ analogue would be a "Siegel constrained Epstein zeta" (not yet defined in the manuscript) encoding the Rankin--Selberg spectral decomposition on $\mathcal{M}_g$.

---

## 4. The Key Structural Observation: Arity vs. Genus

The shadow obstruction tower lives on the **arity axis**: it consists of the projections $\mathrm{Sh}_r(\Theta_\mathcal{A}) = \pi_r(\Theta_\mathcal{A}) \in \mathcal{A}^{\mathrm{sh}}_{r,\bullet}$ at genus 0 and varying arity $r = 2, 3, 4, \ldots$.  These are the OPE-derived invariants: $\kappa$ (arity 2), the cubic shadow $C$ (arity 3), the quartic contact invariant $Q^{\mathrm{contact}}$ (arity 4), etc.

The constrained Epstein zeta lives on the **genus axis**: it is extracted from the genus-1 partition function $Z_1(\tau)$ by Rankin--Selberg transform over $\mathcal{M}_{1,1}$, with no marked points.

These are genuinely different axes of the bigraded shadow algebra $\mathcal{A}^{\mathrm{sh}}_{r,g}$ (Definition def:shadow-algebra, Vol I higher_genus_modular_koszul.tex):

| Axis | Data | Lives at | Invariant |
|------|------|----------|-----------|
| **Arity** (shadow obstruction tower) | $\mathcal{A}^{\mathrm{sh}}_{r,0}$ | genus 0, arity $r$ | $\kappa, C, Q, \mathrm{Sh}_r$ |
| **Genus** (partition functions) | $\mathcal{A}^{\mathrm{sh}}_{0,g}$ | genus $g$, arity 0 | $F_g = \kappa \lambda_g^{\mathrm{FP}}$ |
| **Mixed** | $\mathcal{A}^{\mathrm{sh}}_{r,g}$ | genus $g$, arity $r$ | Higher sewing-shadow data |

The **tridegree** (def:vol1-rigid-planted-forest-depth-filtration) assigns each component of $\Theta_\mathcal{A}$ a triple $(g, n, d)$ = (loop genus, arity, log boundary depth).  The shadow obstruction tower is the $(0, r, \bullet)$ slice.  The genus tower is the $(g, 0, \bullet)$ slice.  The constrained Epstein zeta is extracted from the $(1, 0, \bullet)$ component.

**Can they be unified via factorization homology on bordered surfaces?**

In principle, yes.  The open/closed MC element $\Theta^{\mathrm{oc}} = \Theta_\mathcal{A} + \sum \mu^{M_j}$ (Vol I, constr:thqg-oc-mc-element) packages both the closed-sector genus tower AND the open-sector module structures.  The bordered moduli spaces $\overline{\mathcal{M}}_{g,n}^{\mathrm{bor}}$ (with $k$ boundary components) interpolate between the arity axis (many boundary components at genus 0) and the genus axis (no boundary components at higher genus).

However, this unification is FORMAL, not analytic.  The bordered FM compactification provides the combinatorial framework (planted-forest corrections on bordered surfaces), but the analytic step -- relating the bordered factorization homology to the Rankin--Selberg spectral decomposition -- has not been established.  The reason:

- The shadow obstruction tower at genus 0 is ALGEBRAIC (determined by finitely many OPE coefficients $\kappa, \alpha, S_4$).
- The constrained Epstein zeta at genus 1 is ANALYTIC (its $L$-function content comes from the Hecke decomposition of the theta function, which involves infinite-dimensional spectral theory on $\mathcal{M}_{1,1}$).

The structural obstruction (Vol I, rem:structural-obstruction) is precisely this: algebraic constraints from the MC equation act on the real spectral line, while the zeta zeros live at complex spectral parameters.  No bordered surface interpolation changes this.

**Verdict**: The bigraded MC element $\Theta_\mathcal{A}$ does formally unify the arity and genus axes via its $(g, r)$-components.  Bordered surfaces provide the geometric substrate.  But the analytic content (L-function factorization, critical lines) is visible only after integrating over moduli, which is a different operation from computing factorization homology at a single surface.

---

## 5. Heisenberg on Bordered Surfaces: Interpolation?

For the Heisenberg algebra $\mathcal{H}_\kappa$ on a bordered surface $\Sigma_{g,k}$ (genus $g$, $k$ boundary components):

- **$k$ large, $g = 0$**: A genus-0 surface with $k$ punctures is $\mathbb{P}^1$ minus $k$ points.  The factorization homology $\int_{\mathbb{P}^1 \setminus \{p_1, \ldots, p_k\}} \mathcal{H}$ encodes $k$-point correlation functions of the Heisenberg on $\mathbb{P}^1$.  These are the genus-0, arity-$k$ shadow components: the Wick contractions with propagator $1/(z_i - z_j)$.  For $\mathcal{H}$, the shadow obstruction tower terminates at arity 2 (class G, shadow depth 2), so the $k$-point functions at $k \geq 3$ are determined by Wick's theorem from the 2-point function.  The shadow obstruction tower at genus 0 IS the collection of $k$-point functions.

- **$k = 0$, $g = 1$**: The torus partition function $Z_1(\tau) = \eta(\tau)^{-1}$ (up to normalization).  The constrained Epstein zeta $\varepsilon^1_s = 2\zeta(2s)$ (rank-1 lattice, one critical line).

- **$k = 1$, $g = 0$**: The disk with one puncture.  Factorization homology gives the Fock space $\mathcal{F}_\kappa$.

- **$k = 1$, $g \geq 1$**: A genus-$g$ surface with one boundary component.  The factorization homology $\int_{\Sigma_{g,1}} \mathcal{H}$ is a vector in $\mathcal{F}_\kappa$ (the state on the boundary circle).  This is the genus-$g$ vacuum amplitude with one outgoing state.  For the Heisenberg, this is essentially $|\eta(\Omega)|^{-2}$ evaluated on the genus-$g$ period matrix $\Omega$.

- **$k = 2$, $g = 0$**: The annulus.  $\int_{\mathrm{Ann}} \mathcal{H} = HH_*(\mathcal{H})$ by the annulus trace theorem.  The annulus partition function is $\mathrm{tr}_{\mathcal{F}} q^{L_0 - \kappa/24} = \eta(q)^{-1}$.

**Does this interpolate between the shadow obstruction tower and the constrained Epstein?**

Not in a useful sense.  The interpolation is:
- At $(g, k) = (0, r)$: the $r$-point correlation function, which is the arity-$r$ shadow at genus 0.
- At $(g, k) = (1, 0)$: the torus partition function, from which $\varepsilon^c_s$ is extracted.

But these live at different stages of the Chriss--Ginzburg architecture:
1. The shadow obstruction tower is extracted from $\Theta_\mathcal{A}$ by projecting to arity $r$ at genus 0.  This is a purely algebraic operation.
2. The constrained Epstein zeta is extracted from $Z_1(\tau)$ by a Rankin--Selberg integral over $\mathcal{M}_{1,1}$.  This is a spectral-theoretic operation.

The bordered surface $\Sigma_{g,k}$ provides a geometric interpolation in the sense that the compactified moduli space $\overline{\mathcal{M}}_{g,k}^{\mathrm{bor}}$ degenerates, as $k \to 0$, to the closed moduli space $\overline{\mathcal{M}}_g$, and as $g \to 0$, to the configuration space $\mathrm{FM}_k(\mathbb{C})$.  The MC equation on the full bordered modular operad connects these degenerations.

For the Heisenberg specifically: the bordered partition function $Z_{g,k}(\Omega; z_1, \ldots, z_k)$ is a section of a determinant line bundle on $\overline{\mathcal{M}}_{g,k}^{\mathrm{bor}}$.  At $k = 0$, it specializes to $Z_g(\Omega)$.  At $g = 0$, it specializes to the $k$-point function on $\mathbb{P}^1$.  But the arithmetic content (the $L$-function factorization) appears only after integrating over moduli, which is not part of the factorization homology at a single surface.

**Verdict**: The bordered surface picture provides a geometric interpolation between the genus-0 arity-$k$ data (shadow obstruction tower) and the genus-$g$ arity-0 data (partition functions), but this interpolation does NOT directly capture $\varepsilon^c_s$.  The constrained Epstein zeta requires the additional step of spectral decomposition on $\mathcal{M}_{1,1}$, which is extrinsic to the factorization homology on any single surface.

---

## Summary Assessment

| Question | Answer | Status |
|----------|--------|--------|
| Does $\int_{T^2} V$ capture $\varepsilon^c_s$? | No directly; $\varepsilon^c_s$ is the Rankin--Selberg transform of $\widehat{Z}^c$ over $\mathcal{M}_{1,1}$ | Proved (Vol I, thm:sewing-selberg-formula) |
| Does $\int_{S^1} V$ determine the primary spectrum? | For rational VOAs, yes (via Zhu); for irrational, partially | Proved for rational (Zhu correspondence) |
| Does $\int_{\Sigma_g} V$ overdetermine the spectrum? | Yes: genus-$g$ data accesses higher Siegel modular forms and spinor $L$-functions invisible at genus 1 | Proved for lattice VOAs (Bocherer bridge) |
| Can arity and genus axes be unified? | Formally yes, via the bigraded MC element $\Theta_\mathcal{A} \in \mathcal{A}^{\mathrm{sh}}_{r,g}$ | Structure exists; analytic content differs |
| Does $\int_{\Sigma_{g,k}} \mathcal{H}$ interpolate? | Geometrically yes; arithmetically no | Structural obstruction (rem:structural-obstruction) |

### The Root Distinction

Factorization homology $\int_\Sigma V$ produces a partition function at a SINGLE surface $\Sigma$.  The constrained Epstein zeta $\varepsilon^c_s$ is extracted by integrating the partition function over the MODULI SPACE of surfaces.  These are categorically different operations:

$$
\int_\Sigma V \;\;\text{(factorization homology at one surface)} \quad\neq\quad \int_{\mathcal{M}_{1,1}} \left(\int_{E_\tau} V\right) E_s(\tau)\, d\mu(\tau) \;\;\text{(Rankin--Selberg transform)}
$$

The shadow obstruction tower lives in the fiber of the universal family (it is the local deformation-theoretic data at a single point of the curve).  The constrained Epstein zeta lives on the base (it is the global spectral-theoretic data on the moduli space).  The MC element $\Theta_\mathcal{A}$ unifies the fiber direction (arity) with the base direction (genus) through its bigraded structure, but extracting $\varepsilon^c_s$ requires the additional passage from fiber to spectral decomposition on the base.

### Honest Open Questions

1. **The Siegel constrained Epstein zeta**: at genus $g \geq 2$, one can define a "genus-$g$ constrained Epstein zeta" via the Rankin--Selberg spectral decomposition on $\mathcal{M}_g$.  The relationship between this object and the genus-$g$ components $\mathcal{A}^{\mathrm{sh}}_{0,g}$ of the shadow algebra is unexplored.  For the Leech lattice at genus 2, the Bocherer bridge (Vol I) provides a first computation.

2. **Bordered Rankin--Selberg**: Is there a Rankin--Selberg formula on the bordered moduli space $\overline{\mathcal{M}}_{g,k}^{\mathrm{bor}}$ that interpolates between the shadow obstruction tower (large $k$, $g = 0$) and the constrained Epstein (k = 0, $g = 1$)?  This would require a theory of Eisenstein series on $\overline{\mathcal{M}}_{g,k}^{\mathrm{bor}}$, which does not currently exist.

3. **Arity-genus duality**: The MC equation links the genus axis to the arity axis through the planted-forest corrections: $\delta_{\mathrm{pf}}^{(g,0)} = S_3(10S_3 - \kappa)/48$ at genus 2 involves the arity-3 shadow $S_3$.  This suggests that the genus-$g$ partition function knows about the arity-$r$ shadows for $r \leq 2g$ (the shadow visibility genus: $g_{\min}(S_r) = \lfloor r/2 \rfloor + 1$, Vol I cor:shadow-visibility-genus).  Whether this coupling can be inverted -- extracting shadows from partition functions, or partition functions from shadows -- is open.

---

## Key References in the Manuscript

- **Vol I, arithmetic_shadows.tex**: The constrained Epstein zeta, sewing--Rankin--Selberg bridge, shadow--spectral correspondence, structural obstruction
- **Vol I, thqg_open_closed_realization.tex**: Annulus trace theorem, bordered FM construction, open/closed MC element
- **Vol I, higher_genus_modular_koszul.tex**: Shadow algebra bigrading, shadow extraction, genus spectral sequence, planted-forest corrections
- **Vol II, examples-worked.tex**: SL(2) CS annulus partition function, four-stage architecture instantiation
- **Vol II, modular_pva_quantization_core.tex**: IndHilb-valued factorization homology (Moriwaki), analytic vs algebraic levels
- **Vol II, factorization_swiss_cheese.tex**: Factorization homology on P^1 giving Rep_q(g) (Latyntsev)
- **Vol I, working_notes.tex (archive)**: Honest assessment of the descent chain, structural obstruction, the "fan not chain" picture

## Computation Modules

- `compute/lib/genuine_epstein.py` (Vol I): E_8 and rank-1 lattice Epstein zeta functions with L-function factorization
- `compute/lib/koszul_epstein_steps_bc.py` (Vol II): Scattering factor analysis, zero-forcing mechanism assessment
- `compute/lib/virasoro_epstein_attack.py` (Vol I): Virasoro constrained Epstein, Gap A/B/C framework
- `compute/lib/shadow_epstein_zeta.py` (Vol I): Shadow metric Epstein zeta

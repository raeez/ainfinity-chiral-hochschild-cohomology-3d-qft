# MTC Perspective on the Constrained Epstein Zeta

## Status: Investigation for Benjamin Chang

This note examines whether modular tensor category (MTC) data determines
the constrained Epstein zeta function, and what happens beyond the rational setting.

---

## 1. The Two Zeta Functions: Disambiguation

The monograph uses two distinct objects both called "Epstein zeta":

**(a) The shadow-metric Epstein zeta** (Vol I, thm:shadow-epstein-zeta in
`higher_genus_modular_koszul.tex`, line 15152):
$$
\varepsilon_{Q_L}(s) := \sideset{}{'}\sum_{(m,n) \in \mathbb{Z}^2} Q_L(m,n)^{-s},
\qquad \mathrm{Re}(s) > 1,
$$
where $Q_L(m,n) = 4\kappa^2 m^2 + 12\kappa\alpha\,mn + (9\alpha^2 + 2\Delta)n^2$
is the binary quadratic form associated with the shadow metric on a primary line $L$.
This is defined only for class **M** algebras ($\Delta \neq 0$) and has a functional
equation $\Lambda_{Q_L}(s) = \Lambda_{Q_L}(1-s)$ proved by Poisson summation.
**This is an invariant of the bar complex's nonlinear shadow obstruction tower, not of the module category.**

**(b) The constrained Epstein zeta / primary-counting zeta** (Vol I,
`genus_complete.tex`, line 2472, following Benjamin-Chang):
$$
\varepsilon^c_s = \sum_{\Delta \in S} (2\Delta)^{-s},
$$
where $S = \{h_i\}$ is the set of conformal dimensions of primary fields in
the VOA. This is the object in Benjamin-Chang's paper. It is determined by
the *module category* (the primary spectrum), not by the shadow obstruction tower.

**(c) The shadow-Epstein zeta** (Vol I, `higher_genus_modular_koszul.tex`,
line 17134):
$$
\varepsilon^c_s(\mathcal{A}) = \sum_{r \geq 2} S_r\, r^{-s},
$$
the Dirichlet series whose coefficients are the shadow obstruction tower moments $S_r$.
Its convergence is governed by the shadow growth rate $\rho(\mathcal{A})$.

These are **three distinct Dirichlet series**. The question "does MTC determine
$\varepsilon^c_s$?" makes sense for object (b), which depends on primary dimensions.
Objects (a) and (c) are invariants of the bar complex's shadow structure and live
in the OPE/deformation layer, not the module-category layer.

The Dirichlet-sewing lift $S_{\mathcal{A}}(u)$ (Vol I, `genus_complete.tex`,
line 1863, thm:dirichlet-weight-formula) is yet a fourth object:
$$
S_{\mathcal{A}}(u) = \zeta(u+1) \sum_{i=1}^r \bigl(\zeta(u) - H_{w_i-1}(u)\bigr),
$$
determined entirely by the weight multiset $W(\mathcal{A}) = \{w_1, \ldots, w_r\}$
(conformal weights of strong generators, Definition def:weight-multiset). This
is also distinct from $\varepsilon^c_s$ (which sums over *primaries*, not generators).

**Key structural point** (Vol I, `working_notes.tex`, line 2554--2563):
The shadow obstruction tower depends on $(\kappa, \alpha, S_4)$ (arity 2, 3, 4 OPE data).
The Dirichlet-sewing lift depends on the weight multiset $W(\mathcal{A})$.
These are *independent invariants*: two algebras with the same $W$ but different
OPE structure constants have the same $S_{\mathcal{A}}(u)$ but different shadow obstruction towers.

---

## 2. Rational VOAs and MTC

### Setup
For a rational VOA $V$ of central charge $c$, the module category
$\mathrm{Rep}(V)$ is a modular tensor category (MTC) $\mathcal{C}$ with:
- Finitely many simple objects $\{L_0, L_1, \ldots, L_N\}$ (primaries)
- Conformal dimensions $\{h_0 = 0, h_1, \ldots, h_N\}$
- S-matrix $S_{ij}$ (unitary, symmetric)
- T-matrix $T_{ij} = \delta_{ij}\, e^{2\pi i(h_i - c/24)}$
- Fusion coefficients $N^k_{ij} = \sum_l S_{il}S_{jl}\overline{S}_{kl}/S_{0l}$ (Verlinde formula)

### Does MTC determine $\varepsilon^c_s$?

**Yes, completely.** For a rational VOA with MTC $\mathcal{C}$:

$$
\varepsilon^c_s(V) = \sum_{i=0}^N (2h_i)^{-s}
$$

is a **finite Dirichlet polynomial** (excluding $h_0 = 0$ if the vacuum is removed).
The conformal dimensions $\{h_i\}$ are read off from the T-matrix:
$h_i = (1/2\pi i)\log T_{ii} + c/24 \pmod{\mathbb{Z}}$.

Since the S-matrix and T-matrix together generate a projective representation
of $\mathrm{SL}(2,\mathbb{Z})$, and the T-matrix encodes $\{h_i \bmod 1\}$, the
MTC data determines $\{h_i \bmod 1\}$ exactly. **But not $\{h_i\}$ itself.**

The T-matrix gives $e^{2\pi i h_i}$, hence $h_i \bmod 1$. To recover the
actual conformal dimension $h_i \in \mathbb{Q}_{\geq 0}$, one needs additional
information (e.g., the Zhu algebra, character formulas, or the positivity
constraint $h_i \geq 0$). For $h_i < 1$, the MTC determines $h_i$ unambiguously
(since $h_i \bmod 1 = h_i$ when $0 \leq h_i < 1$). For $h_i \geq 1$, there is
an integer ambiguity.

**In practice**, this ambiguity is resolved for all standard families by the
representation theory of the VOA. But abstractly, the MTC alone determines
$\varepsilon^c_s$ only up to integer shifts of conformal dimensions.

**Stronger statement for integrable affine VOAs:** For $V_k(\mathfrak{sl}_2)$
at integrable level $k \in \mathbb{Z}_{\geq 0}$, the conformal dimensions of the
$k+1$ simple modules $L_{k,j}$ ($0 \leq j \leq k$) are
$$
h_j = \frac{j(j+2)}{4(k+2)},
$$
all of which satisfy $0 \leq h_j < 1$ for $j \leq k$. So $\varepsilon^c_s$ is
fully determined by the MTC. This is verified at $k=1$ in the benchmark
(examples-worked.tex, line 2822): two simple modules with $h_0 = 0$,
$h_1 = 3/4(k+2) = 1/4$, giving
$$
\varepsilon^c_s(V_1(\mathfrak{sl}_2)) = (2 \cdot 0)^{-s} + (2 \cdot 1/4)^{-s}
= 2^s.
$$
(The $h_0 = 0$ term diverges and is typically regularized or excluded.)

### Does Verlinde constrain primary dimensions beyond characters?

Not directly. The Verlinde formula
$$
N^k_{ij} = \sum_l \frac{S_{il}\,S_{jl}\,\overline{S}_{kl}}{S_{0l}}
$$
constrains the S-matrix entries via $N^k_{ij} \in \mathbb{Z}_{\geq 0}$
(integrality of fusion coefficients). This in turn constrains possible S-matrices,
and hence (via the modular representation) constrains the T-matrix and
therefore $\{h_i \bmod 1\}$.

The deep constraint is the **Vafa theorem**: for an MTC with $N+1$ simple objects,
$h_i \in \mathbb{Q}$ and the order of the T-matrix divides a computable function
of $N$. So the set of possible $\{h_i\}$ is a priori constrained to a finite set
(given $c$ and $N$). The Verlinde integrality adds further restrictions.

But these are constraints on the *possible MTCs*, not new information about
$\varepsilon^c_s$ given an MTC. Once the MTC is specified, $\varepsilon^c_s$ is
determined (up to the integer ambiguity above).

---

## 3. The P^1 Benchmark

The benchmark (Vol II, `examples-worked.tex`, line 1785--2933) computes the
complete modular Koszul datum for $V_k(\mathfrak{sl}_2)$ on
$(\mathbb{P}^1, \{0, \infty\})$.

At integrable level $k$, the S-matrix is (eq:benchmark-s-matrix, line 2579):
$$
S_{ij} = \sqrt{\frac{2}{k+2}}\,\sin\!\biggl(\frac{\pi(i+1)(j+1)}{k+2}\biggr),
\qquad 0 \leq i, j \leq k.
$$

At $k=1$ (eq:benchmark-s-matrix-k1, line 2850):
$$
S = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}
= \frac{1}{\sqrt{2}} H,
$$
the Hadamard matrix (up to normalization). The Verlinde algebra is
$\mathbb{C} \times \mathbb{C}$ (two simple idempotents).

The S-matrix determines:
- **Fusion rules**: $L_1 \otimes L_1 \cong L_0$ (benchmark, line 2650).
- **Modular properties of characters**: $\chi_j(-1/\tau) = \sum_l S_{jl}\,\chi_l(\tau)$.
- **WZW partition function**: $Z(\tau) = \sum_j |\chi_j(\tau)|^2$ (diagonal modular
  invariant, eq:benchmark-wzw-partition, line 2248).

But the S-matrix does **not** determine:
- The shadow obstruction tower invariants $(\kappa, \alpha, S_4, \ldots)$, which come from
  the OPE structure constants of the *vacuum module*, not from the module category.
- The Dirichlet-sewing lift $S_{\mathcal{A}}(u)$, which depends on the weight
  multiset of *generators*.
- The shadow-metric Epstein zeta $\varepsilon_{Q_L}(s)$, which depends on the
  nonlinear shadow invariants.

The S-matrix does determine $\kappa$ indirectly: the central charge $c$ is
encoded in the T-matrix, and $\kappa = \dim(\mathfrak{g})(k + h^\vee)/(2h^\vee)$
depends on the level $k$, which can be recovered from $c = 3k/(k+2)$. But this
uses the *identification* of the VOA as $V_k(\mathfrak{sl}_2)$, not just the
abstract MTC data.

**Important distinction** (rem:benchmark-generic-vs-integrable, line 2596):
At generic (non-integrable) $k$, the algebraic derived center has $H^0 \cong \mathbb{C}$.
At integrable $k$, the *categorical* center has dimension $k+1$ (the Verlinde algebra).
"The bar complex $\bar{B}(V_k(\mathfrak{sl}_2))$ classifies twisting morphisms
(universal couplings between $V_k$ and $V_{-k-4}$); it does not see the Verlinde
algebra directly." The MTC is an invariant of the categorical center, not of the
bar complex.

---

## 4. Non-Rational VOAs: Vertex Tensor Categories

### What happens beyond rationality?

For a non-rational VOA (e.g., generic-level affine KM, Virasoro at irrational $c$,
free fields), the module category is *not* an MTC. The theory of
**vertex tensor categories** (Huang-Lepowsky-Zhang, arXiv:1012.4193) provides the
correct framework.

Key differences from the MTC setting:

| Property | Rational (MTC) | Non-rational (VTC) |
|----------|---------------|-------------------|
| Simple objects | Finitely many | Infinitely many (or continuum) |
| Semisimplicity | Yes | No |
| S-matrix | Finite, unitary | Not defined as a matrix |
| Tensor product | Exact, bifunctorial | May fail exactness |
| Braiding | Non-degenerate | May be degenerate |
| $\varepsilon^c_s$ | Finite Dirichlet polynomial | Dirichlet series (possibly divergent) |

For non-rational VOAs, the primary spectrum $\{h_i\}$ is typically infinite
(and may be a continuum for non-C_2-cofinite theories). The constrained Epstein
zeta $\varepsilon^c_s = \sum_i (2h_i)^{-s}$ becomes a genuine Dirichlet series
with infinitely many terms. Its analytic properties (meromorphic continuation,
functional equation, Euler product decomposition) are nontrivial.

### The Virasoro case

For the Virasoro algebra $\mathrm{Vir}_c$ at generic $c$, the module category is
not semisimple and has infinitely many simple objects (the Verma modules $M_{c,h}$
for $h \geq 0$). The primary spectrum is continuous: $S = [0, \infty)$. The
constrained Epstein zeta becomes an integral:
$$
\varepsilon^c_s(\mathrm{Vir}_c) = \int_0^\infty \rho(h)\,(2h)^{-s}\,dh,
$$
where $\rho(h)$ is the spectral density. This is the Mellin transform of the
spectral density, not a sum over finitely many primaries.

The Dirichlet-sewing lift (Vol I, genus_complete.tex, line 1883) for Virasoro is:
$$
S_{\mathrm{Vir}}(u) = \zeta(u+1)(\zeta(u) - 1),
$$
which depends only on the weight multiset $W = \{2\}$ (single generator of weight 2),
not on $c$. This is a manifestation of the structural obstruction
(working_notes.tex, line 2565--2585): the Virasoro character is not a modular form
of finite weight, the Hecke decomposition does not apply, and the shadow spectral
measure is a different object from the automorphic spectral measure.

### The structural obstruction (working_notes.tex, line 2530--2596)

The descent chain has five levels:

| Level | Object | Status |
|-------|--------|--------|
| 0 | Koszul duality $\mathcal{A} \mapsto \mathcal{A}^!$ | Proved |
| 1 | Bar as Fourier: $\bar{B}_X$, inversion, Plancherel | Proved |
| 2 | Shadow obstruction tower as quadratic extension | Proved |
| 3 | Dirichlet-sewing lift $S_{\mathcal{A}}(u)$ | Defined, but independent of Level 2 |
| 4 | Zeros of $\varepsilon^c_s$ on critical lines | Lattice only; structural obstruction |

The break between Levels 2 and 3 is the key: the shadow obstruction tower depends on OPE data
$(\kappa, \alpha, S_4)$, while the Dirichlet-sewing lift depends on the weight
multiset $W(\mathcal{A})$. These are independent invariants.

---

## 5. Logarithmic VOAs

For logarithmic VOAs (e.g., symplectic fermions at $c = -2$, triplet algebras $W(p)$),
the module category has:
- **Indecomposable but non-simple objects**: the category is not semisimple.
- **Jordan blocks for $L_0$**: on indecomposable modules, $L_0$ is not diagonalizable.
  The partition function involves $\log q$ terms: $\mathrm{tr}(q^{L_0}) = q^h(P(\log q))$
  where $P$ is a polynomial.
- **Non-semisimple braided tensor category**: Huang-Lepowsky-Zhang theory still
  provides a vertex tensor category, but it is not an MTC.

For logarithmic spectra, $\varepsilon^c_s$ requires modification:
- The sum $\sum_i (2h_i)^{-s}$ runs over the eigenvalues of $L_0$ on
  *generalized eigenspaces*, not eigenspaces.
- Multiplicities must account for Jordan block sizes.
- The "primary spectrum" now includes values where $L_0$ has nontrivial Jordan form.

The appropriate modification is:
$$
\varepsilon^{c,\log}_s = \sum_{\Delta} d(\Delta)\,(2\Delta)^{-s},
$$
where $d(\Delta)$ counts the dimension of the generalized eigenspace of $L_0$ with
eigenvalue $\Delta$ in the space of primary fields (i.e., fields killed by all
positive modes of the current algebra, but allowing $L_0$ to have Jordan blocks).

For the triplet algebra $W(p)$ at $c = 1 - 6(p-1)^2/p$:
- There are $2p$ simple modules, but the category is non-semisimple.
- The "S-matrix" is a $2p \times 2p$ matrix, but it does not satisfy unitarity
  (it satisfies $S^2 = C$ where $C$ is the charge conjugation matrix).
- The fusion rules are computed by a modified Verlinde formula
  (Fuchs-Hwang-Semikhatov-Tipunin, Gaberdiel-Runkel).

The vertex tensor category structure constrains $\varepsilon^{c,\log}_s$ in the same
way as the MTC constrains $\varepsilon^c_s$ in the rational case: the T-matrix
(suitably generalized) encodes $\{h_i\}$, and the braiding/fusion data provides
additional constraints. But the non-semisimplicity means that:

1. The "S-matrix" is not unitary, so the Verlinde formula gives non-integer
   coefficients (which are then interpreted as dimensions of Hom spaces in
   the non-semisimple category).
2. The modular representation is typically indecomposable but not irreducible.
3. The functional equation of $\varepsilon^{c,\log}_s$ (if it exists) would
   reflect the non-unitary modular representation, not the standard Poisson
   summation of the MTC case.

---

## 6. Summary of Findings

### Does MTC determine $\varepsilon^c_s$?

**For rational VOAs: YES** (up to integer ambiguity in $h_i$). The T-matrix
of the MTC encodes $\{h_i \bmod 1\}$, which for standard families determines
$\{h_i\}$ exactly and hence $\varepsilon^c_s$ as a finite Dirichlet polynomial.
The Verlinde formula provides additional integrality constraints on the S-matrix,
which constrain the possible MTCs (and hence possible $\varepsilon^c_s$), but
does not add information beyond what the MTC already determines.

**The MTC does NOT determine** the shadow obstruction tower invariants $(\kappa, \alpha, S_4)$,
the Dirichlet-sewing lift $S_{\mathcal{A}}(u)$, or the shadow-metric Epstein zeta
$\varepsilon_{Q_L}(s)$. These are bar-complex invariants, not module-category
invariants. The benchmark at line 2596--2608 makes this precise: "the bar complex
classifies twisting morphisms; it does not see the Verlinde algebra directly."

### Can VTC extend to non-rational VOAs?

**Partially.** The Huang-Lepowsky-Zhang vertex tensor category provides the
correct categorical framework, but:
- For non-C_2-cofinite theories: the tensor product may not exist on all modules.
- The primary spectrum may be continuous, making $\varepsilon^c_s$ an integral
  (Mellin transform) rather than a sum.
- The structural obstruction (Vol I, working_notes.tex) is that the shadow spectral
  measure and the automorphic spectral measure are different objects, and the map
  between them (Rankin-Selberg) is not sign-preserving. This obstruction is
  independent of the categorical framework.

### The honest bottom line

The hierarchy of spectral invariants is:

```
Module category (MTC/VTC)
    |
    v  (encodes h_i via T-matrix)
    |
Primary-counting zeta  eps^c_s = sum (2h_i)^{-s}
    |
    |  INDEPENDENT (genus_complete.tex, line 2554-2563)
    |
Bar complex / shadow obstruction tower
    |
    v  (encodes kappa, alpha, S_4, ...)
    |
Shadow-Epstein zeta  eps_{Q_L}(s)
    |
    v  (encodes S_r)
    |
Shadow-Dirichlet series  sum S_r r^{-s}
```

The module category and the bar complex are DIFFERENT layers of the VOA.
The module category classifies representations; the bar complex classifies
deformations (twisting morphisms). For lattice VOAs, these layers are
connected via the theta function (the character IS a modular form, and
the Hecke decomposition works). For Virasoro and W-algebras, the layers
are structurally independent (the character is NOT a modular form of finite weight).

---

## Source Files Referenced

- `/Users/raeez/chiral-bar-cobar-vol2/chapters/examples/examples-worked.tex`
  (P^1 benchmark, lines 1785--2933; S-matrix at line 2579; k=1 data at line 2822)
- `/Users/raeez/chiral-bar-cobar/chapters/theory/higher_genus_modular_koszul.tex`
  (shadow-Epstein zeta, thm:shadow-epstein-zeta at line 15152; shadow-Dirichlet at line 17134)
- `/Users/raeez/chiral-bar-cobar/chapters/connections/genus_complete.tex`
  (constrained Epstein zeta at line 2472; Dirichlet-sewing lift at line 1863;
  weight multiset at line 1838; structural independence at line 2554)
- `/Users/raeez/chiral-bar-cobar-vol2/archive/source_tex/working_notes.tex`
  (descent chain at line 2530; structural obstruction at line 2565; honest
  assessment at line 3282)

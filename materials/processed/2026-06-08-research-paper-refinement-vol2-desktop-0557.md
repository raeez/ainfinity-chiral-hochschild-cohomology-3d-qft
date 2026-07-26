# Research Paper Refinement Vol II: Desktop 05:57 Extraction

Source: `/Users/raeez/Desktop/Research Paper Refinement vol2.pdf`.
Copied to
`materials/raw/2026-06-08-research-paper-refinement-vol2-desktop-0557.pdf`.

SHA-256:
`bba071480ffbfbaa1dbe6072c11be291e54cf940e0d17d2aae2821836c50f0a3`.

The Desktop file differs from the earlier
`materials/raw/2026-06-08-research-paper-refinement-vol2.pdf`
(`9311a84b2a4d75760225894a614825f59dc43f55adebc17607d13148136f3f36`).
The extracted text is at `/tmp/vol2_refinement_20260608_desktop0557.txt`.

## Immediate Formula Obligations

The new source again prioritizes formula-level repairs over
architectural prose.

1. The two-coloured chain object must be written as a dg coloured
   dioperad, with the open-to-closed component equal to zero.
2. The \(A_\infty\)-chiral algebra must be a degree-\(1\)
   coderivation of \(T^c(s^{-1}\bar A)\), with the suspended sign
   convention fixed.
3. Sesquilinearity must be derived from translation equivariance.
   In consecutive variables the source singles out
   \[
   m_k(\ldots,\partial a_i,\ldots)
   =-(\lambda_i+\cdots+\lambda_{k-1})m_k(\ldots,a_i,\ldots),
   \quad i<k,
   \]
   and
   \[
   m_k(a_1,\ldots,\partial a_k)
   =
   \left(\partial_A+\sum_{j=1}^{k-1}j\lambda_j\right)
   m_k(a_1,\ldots,a_k).
   \]
4. Fulton--MacPherson residue signs must be governed by the orientation
   convention
   \(\operatorname{or}(\FM_k(\C))=(-d\varepsilon_S)\wedge
   \operatorname{or}(D_S)\), with the regularized residue
   \(\operatorname{Res}_{D_S}\omega=(\varepsilon_S\iota_{\partial_{\varepsilon_S}}\omega)|_{\varepsilon_S=0}\).
5. The Arnold relation must be the printed Jacobi proof; \(m_4\)
   claims belong to arity-five obstruction computations when the
   source says the arity-four boundary is only Jacobi.
6. The chiral Hochschild complex must be
   \(C^\bullet_{\mathrm{ch}}(A,A)\simeq \Coder(B^{\mathrm{ch}}(A))\),
   with \(d_{\mathrm{Hoch}}=[D_A,-]\).
7. The missing higher homotopy coherence must be a signed
   Boardman--Vogt tree homotopy formula.
8. The \(E_3\)-PBW route is not a slogan; it requires an associated
   graded theorem and a growth/vanishing bound.
9. The 3d HT bridge must be written from a BV datum
   \((E,Q,\omega_{\mathrm{BV}},I,B)\), including propagators, kernels,
   and the QME.
10. Heisenberg, affine, Virasoro, \(W_3\), class-\(\mathsf M\), DS, and
    gravity claims each require their own normalization and status
    gates before serving as theorems.

## Current Pass

A175 repairs item 3 on the active Vol II surface:
`chapters/theory/axioms.tex` now derives sesquilinearity from
translation-equivariant divided-power mode coefficients, and
`chapters/theory/raviolo.tex` uses the same coefficient extraction.

A176 repairs item 1 on the active Vol II surface:
`chapters/theory/sc_chtop_heptagon.tex` now writes
\(C_\bullet^{\log}\SCchtop\) as a dg coloured dioperad with the
open-to-closed component equal to the zero chain complex, and
`chapters/theory/foundations.tex` now writes the oriented
Boardman--Vogt tree complex with its edge-ordering line.

A177 repairs item 4 on the active Vol II surface:
`chapters/theory/fm-calculus.tex`, `chapters/theory/fm-proofs.tex`,
and `chapters/theory/orientations.tex` now use the inward radial
parameter convention
\(\operatorname{or}(\FM_k(\C))=(-d\varepsilon_S)\wedge
\operatorname{or}(D_S)\), the regularized residue
\((\varepsilon_S\iota_{\partial_{\varepsilon_S}}\omega)|_{\varepsilon_S=0}\),
and incidence signs \(\epsilon_{D_S}\) rather than exponentiated
sign symbols.

A178 repairs item 5 on the active Vol II surface:
`chapters/theory/pva-descent-repaired.tex` prints the Arnold
three-term relation as the Jacobi identity, while
`chapters/theory/fm-proofs.tex`,
`chapters/theory/modular_swiss_cheese_operad.tex`,
`chapters/examples/w-algebras-w3.tex`,
`chapters/connections/3d_gravity.tex`, and the Virasoro wheel compute
surface now distinguish the arity-four Stasheff boundary source from a
chain-level \(m_4\) primitive. With \(m_1=0\) the arity-four minimal
identity is \(m_2\circ m_3+m_3\circ m_2=0\); \(m_4\) first appears in
the arity-five minimal identity
\(m_2\circ m_4+m_3\circ m_3+m_4\circ m_2=0\).

A179 repairs item 8 on the active Vol II surface:
`chapters/theory/chiral_higher_deligne.tex` now states the missing
chiral-\(E_3\)-PBW route as a filtered derived-centre package rather
than a slogan: filtered \(E_3\)-envelope
\(R_x\simeq U_{E_3}(W_x)\), associated graded
\(\operatorname{gr}R_x\simeq \Free_{E_3}(H^\bullet(W_x[-2]))\),
convergent PBW spectral sequence, polynomial-growth/amplitude bounds,
and an \(E_1\)-page support bound in total degrees \(\le2\). The
chapter keeps this package conjectural and keeps Theorem H on the
established ordered-bar / Orlik--Solomon route.

A180 repairs item 7 on the active Vol II surface:
`chapters/theory/chiral_higher_deligne.tex` now defines the signed
two-coloured cobar tree operator
\[
h_{\mathrm{oc}}^{\mathrm{tree}}(T)=
\sum_{e\in E_{\mathrm{int}}(T)}(-1)^{\sigma(e,T)}T_e,\qquad
\sigma(e_j,T)=j-1+\sum_{v\prec e_j}|x_v|\pmod 2,
\]
where \(T_e\) inserts the Boardman--Vogt edge-length parameter and the
corresponding suspended cobar generator. The chapter states the exact
contraction identity \(dh+hd=\mathrm{id}-ip\), the mixed-cooperadic
compatibility identity, and the obstruction class
\([o_{\mathrm{oc}}]\in H^1\operatorname{Cone}(\Omega((\SCchtop)^!)
\to\End(Z^{\mathrm{der}}_{\mathrm{ch}}(\cA),\cA))\). The strict
chain-level \(\SCchtop\) lift remains conditional on vanishing of that
class; the \(E_2\)-chiral brace action and cohomological
two-coloured statement are not inflated into a strict Swiss-cheese
theorem.

A180 verification: `git diff --check` is clean on the edited surfaces,
`make fast` converged after two passes with zero undefined citations and
zero undefined references, and `make verify-licensing` reported zero
blocking violations.

A181 repairs item 9 on the active Vol II surface:
`chapters/theory/bv-construction.tex` now starts the 3d HT bridge from
the explicit BV datum
\[
(\mathcal E,Q,\omega_{\mathrm{BV}},I,B),\qquad
\mathcal E=\Omega^\bullet(\R_t)\widehat\otimes
\Omega^{0,\bullet}(\C_z)\otimes\mathfrak a[1],\quad
Q=d_t+\bar\partial_z+d_{\mathfrak a}.
\]
It prints the BV pairing, cyclic classical action, classical master
equation, heat-kernel propagator
\[
P_{\epsilon<L}=\int_\epsilon^L Q^*K_s\,ds,\qquad
QP_{\epsilon<L}=K_L-K_\epsilon,
\]
the Costello graph expansion for \(I[L]\), the graph-weight integral,
the scale-\(L\) QME
\[
QI[L]+\hbar\Delta_LI[L]+\tfrac12\{I[L],I[L]\}_L=0,
\]
and the boundary operation extraction
\[
m_k^{\mathrm{BV}}
=\sum_\Gamma
\frac{\hbar^{b_1(\Gamma)}}{|\operatorname{Aut}\Gamma|}
\operatorname{Res}_{\partial\FM_\Gamma}W_\Gamma(P_{0<L},I),
\]
to be compared with
\(\pi_1D_A|_{(s^{-1}A)^{\otimes k}}\) on the boundary algebra.

A181 verification: `git diff --check` is clean on the edited surfaces,
the focused stale-wording grep found no old generic BV placeholder,
`make fast` converged after two passes with zero undefined citations and
zero undefined references, and `make verify-licensing` reported zero
blocking violations.

A182 begins the item 10 normalization gate on the Heisenberg benchmark.
`compute/lib/examples/abelian_cs.py` now stores the evaluated collision
residue as
\[
r^{\mathrm{coll}}_{H_k}(z;q_1,q_2)=\frac{kq_1q_2}{z}
\]
and the quantum line braiding as
\[
R_{q_1,q_2}(z)=\exp\!\left(\hbar\,\frac{kq_1q_2}{z}\right).
\]
The code and prose now distinguish the pre-\(d\log\) Laplace/OPE
coefficient \(k/z^2\) from the post-bar collision kernel
\(k\,\Omega_{\cH}/z\), and `chapters/examples/examples-worked.tex`
prints the same tensor kernel and charge evaluation. The focused test
surface now checks both the general-level formula and the \(k=0\)
collapse.

A182 verification: `git diff --check` is clean on the edited surfaces,
`compute/.venv/bin/python -m py_compile` passes on the touched compute
files, the virtualenv symbolic check verifies
`abelian_r_matrix(z,hbar,k)=hbar*k*q1*q2/z` and
`abelian_quantum_line_braiding(z,hbar,k)=exp(hbar*k*q1*q2/z)`, focused
pytest reports `120 passed`, `make fast` converged after two passes with
zero undefined citations and zero undefined references, and
`make verify-licensing` reported zero blocking violations.

A183 repairs item 11 on the active Vol II surface. The affine
normalization gate now separates the three representatives:
\[
r_{\mathrm{tr}}(z)=\frac{k\Omega}{z},\qquad
r_{\mathrm{KZ}}(z)=\frac{\Omega}{(k+h^\vee)z},\qquad
r_{\mathrm{KZ,lev}}(z)=\frac{k\Omega}{(k+h^\vee)z}.
\]
At \(k=0\), only the singular trace-form and level-prefixed residues
vanish; the KZ representative is \(\Omega/(h^\vee z)\), and the
non-abelian simple-pole bracket \(f^{ab}_{c}J^c/(z-w)\) survives as the
regular \(z^0\) Lie part after \(d\log\)-absorption. Thus \(k=0\) is a
Heisenberg symmetric-collapse test, not a non-abelian affine
line-category collapse. The live passages in
`chapters/connections/spectral-braiding-core.tex`,
`chapters/connections/line-operators.tex`,
`chapters/connections/bar-cobar-review.tex`, and
`chapters/connections/ordered_associative_chiral_kd_core.tex` now state
the Sugawara/KZ relation as a comparison on the non-critical surface,
not as the scalar identity \(k=1/(k+h^\vee)\).

A183 verification: `git diff --check` is clean on the edited surfaces,
`compute/.venv/bin/python -m py_compile` passes on the touched compute
files, focused pytest reports `50 passed` for
`test_collision_residue_rmatrix.py` and `4 passed` for the KZ flatness
class, fixed-string stale scans find no remaining false affine
\(k=0\) collapse or scalar bridge identity, `make fast` converged after
two passes with zero undefined citations and zero undefined references,
and `make verify-licensing` reported zero blocking violations.

A184 repairs item 12's first Virasoro normalization gate. The gravity
engine no longer identifies the unabsorbed Virasoro OPE/Laplace kernel
\[
r^L(z)=\frac{c/2}{z^4}+\frac{2T}{z^2}+\frac{\partial T}{z}
\]
with the bar collision \(r\)-matrix. `gravity_laplace_kernel_poles()`
now stores the pre-\(d\log\) pole set \(\{4,2,1\}\), while
`gravity_r_matrix_poles()` stores the collision residue
\[
r^{\mathrm{coll}}(z)=\frac{c/2}{z^3}+\frac{2T}{z},
\]
with \(\partial T\) regular after the one-pole absorption. The active
gravity chapter and introduction now attach the CYBE computation to the
collision residue, not to the raw Laplace kernel.

A184 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched gravity engine and test files, focused pytest reports
`103 passed, 1 skipped` across `test_gravity_3d_engine.py` and
`test_collision_residue_rmatrix.py`, stale-claim scans are clean on the
active Vol II surface and find no live Vol I/Vol III contradiction,
`make fast` converged after two passes with zero undefined citations
and zero undefined references, and `make verify-licensing` reported
zero blocking violations.

A185 repairs item 12's Virasoro line-sector scalar formula. The closed
expression
\[
z^{2h}\exp\!\left(-\frac{c}{4z^2}\right)
\]
now appears as the primary-state Ward factor
\[
R^{\mathrm{prim}}_{h,c}(z),
\]
not as the full Virasoro spectral \(R\)-matrix. The proof starts from
the collision residue
\[
r_{\mathrm{Vir}}^{\mathrm{coll}}(z)=\frac{c/2}{z^3}+\frac{2T}{z},
\]
projects \(T\) to \(h\) on the highest-weight line, obtains the scalar
connection
\[
\nabla_h=d-\left(\frac{2h}{z}+\frac{c}{2z^3}\right)\,dz,
\]
and integrates it to
\[
R^{\mathrm{prim}}_{h,c}(z)
=z^{2h}\exp\!\left(-\frac{c}{4z^2}\right).
\]
The no-odd-powers claim is now proved by expanding the single monomial
\(-c/(4z^2)\):
\[
\exp\!\left(-\frac{c}{4z^2}\right)
=\sum_{k\ge0}\frac{(-c/4)^k}{k!}z^{-2k}.
\]

A185 also propagated the correction to live Vol I and the standalone
Virasoro note: \(S_3(\mathrm{Vir}_c)=2\) is now treated as the
weight-three scalar shadow, not as a finite proof of class
\(\mathsf M\); class \(\mathsf M\) remains a quartic/contact completed
tower statement.

A185 verification: Vol II `py_compile` passes, Vol II focused pytest
reports `57 passed, 1 skipped`, Vol I `py_compile` passes on the
touched compute files, Vol I focused pytest reports `151 passed`,
fixed-string stale scans find no remaining live scalar-as-full-\(R\)
or \(S_3\)-as-class-\(\mathsf M\)-proof advertisement on the checked
surfaces, `git diff --check` is clean, Vol II `make fast` converged
after two passes, Vol II `make verify-licensing` reported zero blocking
violations, and Vol I `make fast` exited 0 with the inherited unrelated
undefined reference `def:kappa-scalar-vs-mumford`.

A186 repairs item 12's HPL-transfer and class-\(\mathsf M\) ambient
gate. The Virasoro recursion theorem now states the actual transferred
minimal operation
\[
m_n^{\min}(x_1,\ldots,x_n)
=\sum_{\tau\in\mathrm{PRT}_n}(-1)^{\epsilon(\tau;x)}
p\,\mu_\tau(i x_1,\ldots,i x_n;h,\ldots,h),
\qquad
Qh+hQ=1-ip,
\]
with side conditions \(pi=1_H\), \(ph=hi=h^2=0\). On a positive
conformal-weight complement with antighost \(G_0\) satisfying
\([Q,G_0]=L_0\), the homotopy is normalized as
\[
h|_{A_w}=w^{-1}G_0,\qquad w>0.
\]
The BV/MC equation is now the packaged Stasheff identity after this
transfer, not a substitute for the transfer formula.

A186 also corrects the gravity tower language. Nonzero
\[
S_4(\mathrm{Vir}_c)=\frac{10}{c(5c+22)}
\]
is recorded as the first bar-weight-four obstruction to raw
direct-sum class-\(\mathsf M\) termination. The all-arity statement
\(m_k\ne0\) for generic \(c\) is the wheel-residue survival conclusion
under the KZ analytic SDR package, not an induction from \(m_3\) alone.
The chain-level class-\(\mathsf M\) statement is therefore placed in
the weight-completed/pro ambient \(\hypAmbientWtCpl\). The same
correction was propagated to the live Vol I ordered-bar paragraph.

A186 verification: Vol II `py_compile` passes on the touched gravity
engine and test file. Focused pytest reports `112 passed, 1 skipped`
across `test_gravity_3d_engine.py`,
`test_f1_bv_bar_classm_higher_genus.py`, and
`test_class_m_ambient_iv.py`. Fixed-string propagation scans find no
remaining live tower-overclaim pattern on active Vol II or Vol I
surfaces. `git diff --check` is clean on touched Vol II and Vol I
files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations. Vol I `make fast` exited 0 after
four passes with zero undefined citations and the inherited unrelated
undefined reference `def:kappa-scalar-vs-mumford`.

A187 repairs item 13's \(\mathcal W_3\) \(\Lambda\)-channel gap. The
W3 chapter now gives the missing projection
\[
P_\Lambda(X)=
\frac{\langle X,\Lambda\rangle}{\langle\Lambda,\Lambda\rangle}\Lambda,
\qquad
P_\Lambda m_2(W,W;\lambda)=\frac{32}{5c+22}\lambda\Lambda,
\]
on the generic non-isotropic \(\Lambda\)-line, defines the
coefficient-one channel
\[
\mu_\Lambda^\circ(W_a,W_b;\lambda_{ab})=\lambda_{ab}\Lambda,
\]
and inscribes the pair-corner FM kernel
\[
K_{ab|cd}(\lambda)
=\operatorname{Res}_{D_{ab|cd}}(\omega_{ab}\wedge\omega_{cd}).
\]
The coefficient-one profile of the arity-four summand is now
\[
(m_{4,\Lambda}^{\mathrm{tree}})^\circ
=
\frac{32}{5c+22}
\sum_{(ab|cd)\in\mathsf{Pair}^{\mathrm{pl}}_4}
\varepsilon_{ab|cd}\,
p\mu(h\mu_\Lambda^\circ(W_a,W_b),
     h\mu_\Lambda^\circ(W_c,W_d))K_{ab|cd}.
\]
The literal HPL projector formula is displayed separately with
\(P_\Lambda m_2\) inside the two internal brackets, so the
\(\beta\)-powers are supplied by the projections themselves and are
not double-counted. The primary-line normalization
\[
\Lambda_0|h\rangle=(h^2-3h/5)|h\rangle
\]
is also recorded, hence the \(\Lambda\)-channel vanishes at
\(h=0,3/5\).

A187 also corrects the W3 recursive master equation to the vanishing
Stasheff/Maurer--Cartan form and tightens the DNP W3 Hamiltonian
paragraph: local PVA/Jacobi coefficient data are not a commutativity
theorem; commutativity requires a supplied \(\mathcal W_3\)
Gaudin/RLL algebra. The active preface, W-algebra frontier comparison,
and standalone preface survey now refer to the normalized
\(WW\Lambda\) projection and the codimension-\(2\) boundary residue
rather than only the scalar coefficient.

A187 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched W3 engine and test file. Focused W3 pytest reports
`63 passed`. `git diff --check` is clean on touched files. Vol II
`make fast` converged after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reported zero
blocking violations and zero warnings.

A188 repairs item 14's class-\(\mathsf M\) Banach-completion gate. The
class-\(\mathsf M\) topologization chapter now defines the analytic
body
\[
\widehat A_\rho
=
\left\{a=\sum_{w\ge0}a_w:
\|a\|_\rho=\sum_{w\ge0}\|a_w\|_w\rho^w<\infty\right\},
\]
and separates it from the modular cochain norm. It also defines the
operation constants
\[
C_k(\rho)=
\sup_{\|a_i\|_\rho\le1}
\|m_k^\rho(a_1,\ldots,a_k)\|_\rho
\]
so that
\[
\|m_k^\rho(a_1,\ldots,a_k)\|_\rho
\le C_k(\rho)\prod_i\|a_i\|_\rho.
\]
The new Banach MC criterion states that
\[
\sum_{k\ge2}C_k(\rho)<\infty
\]
implies absolute convergence of the transferred coderivation and the
Maurer--Cartan/Stasheff expression on the unit ball of
\(\widehat A_\rho\), by the Weierstrass \(M\)-test and Cauchy-product
majorant. This remains a Banach statement and does not imply raw
direct-sum finite propagation.

A188 also records the status-labelled radii:
\[
\rho_*^{\mathrm{Vir}}(c)=|c|/6,\qquad
\rho_*^{\cW_3}(c)=|c|/10,\qquad
\rho_*^{\cW_N}(c)=|c|/[12(H_N-1)]
\]
with the \(N\ge4\) row conditional on the finite Riccati envelope and
harmonic \(\beta_N\)-law. The compute layer now has finite-truncation
guards for the Banach norm, multilinear operation inequality, MC
majorants, and a radius certificate marking \(N=2,3\) as proved
low-rank rows and \(N\ge4\) as conditional.

A188 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched WN/topologization compute files. Focused pytest reports
`24 passed` across `test_wn_tempered_closure.py` and
`test_topologization_class_m_original_complex.py`. `git diff --check`
is clean on touched files. Vol II `make fast` converged after two
passes with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A189 repairs item 15's DS--Hochschild HPL-transfer gate. The
DS-Hochschild theorem now starts from the BRST realization
\[
\cW_k(\fg,f)=
H^0(V_k(\fg)\widehat\otimes\cF_{\mathrm{gh}},Q_{\mathrm{DS}})
\]
and prints the KRW charge
\[
Q_{\mathrm{DS}}
=
\operatorname*{Res}_{z=0}\left[
\sum_\alpha(J^{u_\alpha}-\chi(u_\alpha))\psi_\alpha^*
-\frac12\sum_{\alpha,\beta,\gamma}f_{\alpha\beta}^{\gamma}
:\psi_\alpha^*\psi_\beta^*\psi_\gamma:
+Q_{\mathrm{neut}}
\right]dz .
\]
It now states the Hochschild SDR
\[
(\ChirHoch^\bullet(V_k,V_k)\widehat\otimes\cF_{\mathrm{gh}},Q_{\mathrm{DS}})
\rightleftarrows
\ChirHoch^\bullet(\cW_k,\cW_k),
\qquad
Q_{\mathrm{DS}}h+hQ_{\mathrm{DS}}=1-ip,
\]
and the actual transferred brace formula
\[
\mu_{\cW}^{n}
=
\sum_{T\in\mathrm{PRT}_n}\pm
p\,\mu_T^V(i,\ldots,i;h,\ldots,h).
\]
The line to use is now printed as
\[
\ChirHoch^\bullet(\cW_k(\fg,f),\cW_k(\fg,f))
\simeq
H^\bullet_{\mathrm{DS}}
\ChirHoch^\bullet(V_k(\fg),V_k(\fg))
\]
as \(E_2^{\mathrm{ch}}\)-brace algebras. The theorem now says that
the \(\SCchtop\)-pair lift requires the heptagon hypotheses, the
two-coloured cobar homotopy, and completed mixed estimates.

A189 also records the fractional-good-grading mechanism: for a
\(\frac1q\mathbb Z\)-good grading, pass to the cover \(z=s^q\), so
Kazhdan degree \(j/q\) becomes integral degree \(j\), construct
\(G'_f\) upstairs with \(T_{\mathrm{DS}}=[Q_{\mathrm{DS}},G'_f]\),
check \(\mu_q\)-invariance, and descend. The downstream
`hochschild.tex` HKR-scope remark no longer says class-\(\mathsf M\)
bulk is unconditional at chain level; it is cohomological and
weight-completed/pro unless bounded-shift HPL or Banach estimates are
supplied.

A189 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched DS-Hochschild IV test. Focused pytest reports `7 passed`
across `test_chd_ds_hochschild_iv.py` and
`test_chiral_higher_deligne.py`. `git diff --check` is clean on
touched files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A190 repairs item 16's explicit modular graph-complex gate. The modular
PVA core now states the relative modular graph complex as a theorem. It
separates
\[
\mathfrak{o}_{\mathrm{GK}}(\Gamma)
=
\det(kE(\Gamma))\otimes\det(H_1(\Gamma;k))^{-1}
\]
from the stable-curve Hodge convention
\[
\mathfrak{o}_{\mathrm{Hdg}}(\Gamma)
=
\det(kE(\Gamma))\otimes
\bigotimes_{v\in V(\Gamma)}
\det(H^1(\Sigma_{g(v)};k))^{-1},
\]
prints the relative Feynman transform
\[
\FT_{\mathrm{rel}}(C)
=
\prod_{[\Gamma]}
\left(\mathfrak{o}(\Gamma)\otimes
\bigotimes_{v}sC(g(v),\mathrm{Fl}(v))\right)_{\Aut(\Gamma)},
\]
and records the dual contraction convention
\[
D_{\mathrm{sep}}=\sum_{e\ \mathrm{sep}}\pm\operatorname{contract}_e,
\qquad
D_{\mathrm{nsep}}=\sum_{e\ \mathrm{nsep}}\pm\xi_e.
\]
The theorem now states
\[
D_0^2=0,\qquad D_1^2=0,\qquad D_0D_1+D_1D_0=0
\]
and the genus-expanded MC recursion
\[
D_0\Theta^{(0)}+\tfrac12[\Theta^{(0)},\Theta^{(0)}]=0,
\]
\[
D_0\Theta^{(1)}+D_1\Theta^{(0)}
+[\Theta^{(0)},\Theta^{(1)}]=0,
\]
\[
D_0\Theta^{(g)}+D_1\Theta^{(g-1)}
+\tfrac12\sum_{a+b=g}[\Theta^{(a)},\Theta^{(b)}]=0.
\]
The obstruction class is now
\[
[o_g]=
\left[
D_1\Theta^{(g-1)}
+\tfrac12\sum_{\substack{a+b=g\\a,b<g}}
[\Theta^{(a)},\Theta^{(b)}]
\right]
\in
H^2(\gr_F^gL,D_0+[\Theta^{(0)},-]).
\]
`relative_feynman_transform.tex` points its skeleton definition to this
theorem, and the compute engine now guards the orientation conventions,
\(D_0/D_1\) identities, recursion terms, and lower-genus obstruction
pairs.

A190 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched relative modular engine and tests. Focused pytest reports
`72 passed` across `test_factorization_modular_engine.py` and
`test_modular_bar_d_squared.py`. `git diff --check` is clean on touched
files. Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.

A191 repairs item 17's genus-one Arakelov kernel gate. The
factorization chapter now starts the genus-one discussion from the
normalized elliptic Green function
\[
G_\tau(z)
=
-\log\left|
\frac{\vartheta_1(z|\tau)}{\eta(\tau)}
\right|^2
+
\frac{2\pi(\operatorname{Im}z)^2}{\operatorname{Im}\tau},
\qquad
\partial\bar\partial G_\tau
=
\pi i(\delta_0-\mu_\tau).
\]
It fixes the residue-one chiral propagator convention
\[
P_\tau(z):=-\partial_zG_\tau(z)\,dz
=
\left(
\partial_z\log\vartheta_1(z|\tau)
+
\frac{2\pi i\,\operatorname{Im}z}{\operatorname{Im}\tau}
\right)dz,
\]
and states that the opposite \(+\partial_zG_\tau\,dz\) convention has
residue \(-1\) and reverses edge-orientation signs.

A191 also adds the scalar genus-tower theorem:
\[
F_g^{\mathrm{scal}}(\cA)
=
\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\qquad
\lambda_g^{\mathrm{FP}}
=
\int_{\overline{\cM}_{g,1}}\psi_1^{2g-2}\lambda_g
=
\frac{2^{2g-1}-1}{2^{2g-1}}\frac{|B_{2g}|}{(2g)!}.
\]
For multi-weight algebras the theorem now prints
\[
F_g(\cA)=
\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}
+\delta F_g^{\mathrm{cross}}(\cA),
\qquad
\delta F_2^{\mathrm{cross}}(\cW_3)=\frac{c+204}{16c}.
\]
The compute layer now guards the Green-kernel convention,
\(\lambda_g^{\mathrm{FP}}\), scalar \(F_g\), the \(\cW_3\) genus-two
cross term, and the scalar-plus-cross split.

A191 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched genus-one bridge engine and test. Focused pytest reports
`101 passed` across `test_genus_one_bridge.py` and
`test_f3_givental_stokes_a2_engine.py`. `git diff --check` is clean on
touched files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A192 repairs item 18's tensor-Arakelov component gate. Theorem D's
native Vol II tensor form now prints the basis expansion
\[
K(\cA)=\sum_{\alpha,\beta}K_{\alpha\beta}(\cA)e^\alpha e^\beta
\]
and the stable-graph component formula
\[
K_{\alpha\beta}(\cA)
=
\sum_{\Gamma\in\mathsf{G}^{\mathrm{st}}_{g,2}(\alpha,\beta)}
\frac{1}{|\Aut\Gamma|}
\int_{[\overline{\cM}_\Gamma]}
\left(\prod_{e\in E(\Gamma)}P_{\alpha_e\beta_e}\right)
\left(\prod_{v\in V(\Gamma)}C_v\right).
\]
The scalar curvature is now explicitly the trace
\[
\operatorname{tr}_{\cF}K(\cA)
=
\sum_\alpha K_{\alpha\alpha}(\cA)
=
\kappaChHodge(\cA)\omega_g.
\]
For \(\cW_3\), the theorem now records the pre-trace channel matrix
\[
K(\cW_3)=
\begin{pmatrix}K_{TT}&K_{TW}\\K_{WT}&K_{WW}\end{pmatrix},
\qquad
K_{TT}+K_{WW}=\frac{5c}{6}\omega_g,
\]
and states that \(K_{TW},K_{WT}\) are invisible to the scalar trace and
source
\[
\delta F_2^{\mathrm{cross}}(\cW_3)=\frac{c+204}{16c}.
\]

A192 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched tensor-Arakelov test. Focused pytest reports `43 passed`
across `test_wN_tensor_arakelov_weight_distribution.py` and
`test_f3_givental_stokes_a2_engine.py`. `git diff --check` is clean on
touched files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A193 repairs item 19's curved-MC anomaly-completion gate. The anomaly
completion core now distinguishes the strict central-shadow
transgression algebra from the associative curved twist. The new
definition prints
\[
B_{\Theta}^{\mathrm{curv}}
=
B\langle\eta,u\rangle/(u-\eta^2),
\qquad
d_\eta=d+[\eta,-],
\]
and imposes the curved Maurer--Cartan flatness equation
\[
d\eta+\frac12[\eta,\eta]+\Theta=0,
\qquad
d\eta+\eta^2+\Theta=0
\]
in the associative convention. The text records that the older
\(d\eta=\Theta\) formula belongs only to the strict central-shadow Ore
extension \(B_\Theta\), after \(\Theta\) is a closed central cocycle,
and forms no twisted differential \(d+[\eta,-]\).

The Virasoro class-\(\mathbf M\) example now gives
\[
B_{\Theta}^{\mathrm{Vir}}
=
\Bch(\mathrm{Vir}_c)_{\rho}^{\wedge}
\otimes k\langle\eta,u\rangle/(u-\eta^2),
\qquad
d\eta+\eta^2+\Theta_{\mathrm{Vir}}=0.
\]
The compute layer mirrors this split: `transgression_algebra` is
labelled as the strict central-shadow Ore extension,
`curved_mc_twist_data` prints
\(d_\eta^2=[\Theta+d\eta+\tfrac12[\eta,\eta],-]\), and
`virasoro_curved_class_m_completion()` guards the class-\(\mathbf M\)
completion.

A193 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched anomaly engine and test. Focused pytest reports `36 passed`
in `test_anomaly_completed_engine.py`. `git diff --check` is clean on
touched files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A194 repairs item 20's K3/Borcherds operator gate. The scalar theorem
still proves only
\[
Z_{\mathrm{BPS}}^{K3\times E}
=
(\Phi_{10}^{\mathrm{un}})^{-1}
=
\Delta_5^{-2},
\]
but the live surface now states the conditional operator theorem next
to that scalar bridge. The new theorem
`thm:k3-borcherds-operator-square` places
\[
\operatorPrim{X}\in
Z^0CC^{\mathrm{ch,cyc}}_\bullet(A_X)^\wedge_{\Lambda^{2,1}_{\mathrm{II}}}
\otimes \mathcal L^5\otimes\nu_{\Delta_5},
\qquad
A_X={\SpCh}_{E,C}\PhiFA_3(D^b\mathrm{Coh}(K3\times E)),
\]
and prints the chain-level protected-Pfaffian identity
\[
\protectedPfaff{\operatorPrim{X}}
=
\Delta_5(Z)
=
e^{2\pi i(\rho,Z)}
\prod_{\alpha\in L_+}
(1-e^{2\pi i(\alpha,Z)})^{c(-\alpha^2/2)}.
\]

A194 also prints the Hall side
\[
\CoHA(K3)\to D_{\mathrm{Hall}}(K3)\to\mathfrak g_{\Delta_5}
\]
and the chiral side
\[
D^b\mathrm{Coh}(K3)\xrightarrow{\PhiFA}E_d\text{-}\mathrm{FactAlg}
\xrightarrow{\SpCh}\mathrm{ChirAlg}_C
\xrightarrow{\Zderch}\SCchtop\text{-}\mathrm{Alg}.
\]
The theorem is conditional on the P1 operator package and carries
licensing tags \(\alpha+\beta+\gamma+\varepsilon\); the scalar
reciprocal is now explicitly only the character after the operator
identity and the square have been installed.

The main overview now names \(\operatorPrim{K3\times E}\), its cyclic
chiral Hochschild home, and the Hall/chiral square. The Construction
Problem paragraph now separates the BPS scalar character
\(\Delta_5^{-2}\) from the Oberdieck--Pixton reduced-DT scalar branch
\(-4096\Delta_5^{-2}\), and states that neither scalar branch is the
operator.

The P1 compute engine now exposes `borcherds_root_product_profile()`,
`k3_borcherds_operator_profile()`, and
`k3_borcherds_hall_chiral_square()`. The P1 tests guard the root-product
formula, cyclic-Hochschild placement, scalar non-promotion, both paths
of the square, the four required comparison blocks, and the manuscript
theorem source.

A194 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched P1 engine and test. Focused pytest reports `23 passed` in
`test_p1_protected_pfaffian_chain_level.py`. `git diff --check` is clean
on touched files. Targeted stale-phrase grep finds no remaining scalar
promotion phrase. Vol II `make fast` converged after two passes with
zero undefined citations and references, and Vol II `make
verify-licensing` reported zero blocking violations and zero warnings.

A195 repairs item 21's celestial residue-theorem gate. The bridge now
states the operator-level Mellin residue theorem rather than only the
soft hierarchy:
\[
\operatorname{Res}_{\Delta=1-p}\mathcal M[A_{n+1}]
=
\sum_{i=1}^{n}\mathcal D^A_{p+2,i}\mathcal M[A_n],
\qquad
\mathcal D^A_{a,i}
=
\rho_i(\pi_1D_A|(s^{-1}A)^{\otimes a}).
\]
Equivalently, with arity \(a=p+2\),
\[
\operatorname{Res}_{\Delta=3-a}\mathcal M[A_{n+1}]
=
\sum_i\rho_i(\pi_1D_A|(s^{-1}A)^{\otimes a})\mathcal M[A_n].
\]
The gravitational stress-tensor convention is recorded as the shifted
pole \(\Delta^{\mathrm{st}}=4-a\).

The Virasoro low-degree residues are now separated by quotient:
\[
[S_2]_{\mathrm{soft}}=c/2,\qquad
[S_3]_{\mathrm{soft}}=0,\qquad
[S_4]_{\mathrm{soft}}=\frac{10}{c(5c+22)},\qquad
[S_5]_{\mathrm{soft}}=-\frac{48}{c^2(5c+22)}.
\]
The raw ordered-bar representative remains \(S_3(\mathrm{Vir}_c)=2\);
it is not the reduced soft-residue class.

Propagation fixed the convention bridge in
`universal_celestial_holography.tex`, corrected the older
soft-graviton theorem from \(\Delta_r=2-r\) to the reduced/stress pair
\(\Delta_r^{\mathrm{red}}=3-r\),
\(\Delta_r^{\mathrm{st}}=4-r\), and replaced the stale celestial
transfer limit at \(\Delta\to1\) by
\(\operatorname*{Res}_{\Delta=1-p}\) with the same
\(\rho_i(\pi_1D_\mathcal A)\) operator.

The celestial OPE compute engine now exposes
`celestial_mellin_transform_profile()`, `soft_residue_dimension()`,
`residue_operator_profile()`, `mellin_residue_identity_profile()`, and
`virasoro_residue_coefficients()`. Tests guard the Mellin transform
weight, reduced/arity/stress pole conventions, the coderivation
operator formula, the residue identity, the \(S_3\) gauge-trivial soft
class, and the manuscript theorem source.

A195 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched celestial engine/tests. Focused pytest reports `92 passed`
across `test_celestial_ope_from_shadow.py` and
`test_universal_celestial.py`. `git diff --check` is clean on touched
files. Targeted stale-formula greps find no remaining
`Delta_r = 2 - r` or unqualified `Delta_s = 4 - r corresponds` witness.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.

A196 repairs item 22's gravity-scope gate. The Virasoro rung now
prints the exact algebraic statement
\[
\widehat{Z^{\mathrm{der}}_{\mathrm{ch}}(\mathrm{Vir}_c)}_\rho
\in
E_3^{\mathrm{top}}\text{-}\mathsf{Alg},
\qquad
0<\rho<|c|/6.
\]
It also defines the scalar completed-bar trace
\[
Z_{\bar B}^{\mathrm{grav}}(\hbar;\rho)
=
\operatorname{Tr}_{(\mathrm B^{\mathrm{ch}}(\mathrm{Vir}_c))^\wedge_\rho}
\exp\!\left(\sum_{g\ge0}\hbar^g\Theta^{(g)}_{\mathrm{Vir}_c}\right),
\]
and states that this trace is not automatically the
dynamical-metric gravitational path integral.

The Universal Holography functor surface now says that interpreting
\(\Phi_{\mathrm{hol}}(\mathrm{Vir}_c)\) as pure \(3d\) gravity uses the
Brown--Henneaux dictionary together with separate modular invariance,
vacuum dominance, saddle dominance, and Borel/Stokes hypotheses; the
functor alone supplies the completed algebraic HT object.

The gravity compute engine now exposes
`virasoro_exact_gravity_scope_profile()` and
`virasoro_scalar_bar_trace_profile()`. Tests guard the
\(0<\rho<|c|/6\) bound, the Brown--Henneaux dictionary, the required
physical hypotheses, and the fact that \(Z_{\bar B}^{\mathrm{grav}}\)
is a scalar boundary trace rather than a metric path integral.

A196 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched gravity engine and test. Focused pytest reports
`65 passed, 1 skipped` in `test_gravity_3d_engine.py`.
`git diff --check` is clean on touched files. Targeted greps find no
remaining `path-integral bridge` phrase; remaining
`dynamical-metric path integral` hits are negative scope qualifiers.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.

A197 repairs item 23's ordered chiral bar skeleton.  The foundations
chapter now contains
`thm:ordered-chiral-bar-residue-skeleton`, which states
\[
\mathrm B^{\mathrm{ch}}(A)=T^c(s^{-1}\bar A),\qquad
D_A=\sum_{k\ge1}D_{m_k},\qquad
\pi_1D_A|_{(s^{-1}A)^{\otimes k}}=s^{-1}m_ks^{\otimes k}.
\]
For a collision divisor \(D_S\subset\partial\FM_k(C)\), it records
\[
m_S=\operatorname{Res}_{D_S}
\left(
\prod_{e\in E(S)}d\log(z_{s(e)}-z_{t(e)})
\otimes\mathrm{OPE}_S
\right),
\]
with
\[
\operatorname{Res}_{D_S}(d\log\varepsilon_S\wedge\beta+\alpha)
=\beta|_{\varepsilon_S=0},\qquad
\operatorname{or}(\FM_k(C))=(-d\varepsilon_S)\wedge\operatorname{or}(D_S).
\]
The theorem also states that \(D_A^2=0\) is the residue image of
\(\partial^2[\FM_k(C)]=0\) together with the Arnold--Orlik--Solomon
relation, that
\[
\Zderch{A}=\operatorname{Coder}(\mathrm B^{\mathrm{ch}}(A)),
\qquad d=[D_A,-],
\]
and that the brace coordinates are ordered insertions.  The strict
\(W(\SCchtop)\) clause is conditional on a Drinfeld associator and
two-coloured contracting homotopy \(h_{\mathrm{oc}}\); the
\(E_3^{\mathrm{top}}\) clause is conditional on \(T=[Q,G]\) and
\(Q\)-exact translations.  Class \(\mathsf M\) identities are placed in
\(\widehat A_\rho\) with \(\|m_k\|_\rho\le C_k(\rho)\) and
\(\sum_kC_k(\rho)<\infty\).

The FM-boundary compute module now exposes
`local_residue_convention()` and
`ordered_chiral_bar_residue_skeleton()`.  Tests guard the
residue/orientation convention, coderivation projection, collision
residue, \(D_A^2\) sources, derived-centre differential, ordered
braces, mixed operations, strictification hypotheses, topologisation
condition, class \(\mathsf M\) completion, positive arity, and the
active theorem source.

A197 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched FM-boundary module and ordered KD tests. Focused pytest
reports `41 passed` in `test_ordered_chiral_kd_engine.py` and
`17 passed` across `test_infrastructure.py::TestFMBoundary` plus
`TestOrderedChiralBarResidueSkeleton`. `git diff --check` is clean on
touched files. Label-count greps find exactly one occurrence of each new
theorem/equation label. Vol II `make fast` converged after two passes
with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A198 repairs the boundary-linear local sector.  The active
`ht_bulk_boundary_line_core.tex` Kuranishi section now fixes the
divided-power Taylor convention
\[
F_I=\sum_{q\ge1}\frac1{q!}F_q^I,\qquad
F_C=\sum_{q\ge1}\frac1{q!}F_q^C,
\]
and prints the first massive-field recursions
\[
\varphi_2(t_1,t_2)=-A^{-1}F_2^I(t_1,t_2),
\]
\[
\varphi_3(t_1,t_2,t_3)
=
-A^{-1}\left[
F_3^I(t_1,t_2,t_3)+
\sum_{\mathrm{cyc}}F_2^I(t_1,\varphi_2(t_2,t_3))
\right].
\]
The effective cubic coupling remains
\[
\kappa_3(t_1,t_2,t_3)
=
F_3^C(t_1,t_2,t_3)
-
\sum_{\mathrm{cyc}}
F_2^C(t_1,A^{-1}F_2^I(t_2,t_3)).
\]

The local bulk/line theorem now has the coordinate derived-critical
model
\[
\cO(\dCrit(W_{\eff}))
=
(\kk[[u,\gamma]]\ot\Lambda(\xi,\chi),d_{W_{\eff}}),
\]
with
\[
d_{W_{\eff}}
=
\sum_i\frac{\partial W_{\eff}}{\partial u_i}\frac{\partial}{\partial\xi_i}
+
\sum_a\frac{\partial W_{\eff}}{\partial\gamma_a}
\frac{\partial}{\partial\chi_a}.
\]
The superseded `ht_bulk_boundary_line.tex` copy carries the same
unlabelled formulas to prevent drift.

The new Kuranishi compute oracle exposes `massive_i2()`,
`massive_i3()`, `kappa2()`, `kappa3()`,
`minimal_line_operation_profile()`, and
`derived_critical_dga_profile()`.  Tests guard the negative Schur
complement, cyclic cubic exchange, direct-minus-exchange \(\kappa_3\),
\(1/n!\) line-operation normalization, derived-critical DGA model, and
the active manuscript source.  Propagation greps found no
line-algebra-equals-bulk collapse.

A198 verification: `compute/.venv/bin/python -m py_compile` passes on
the new module/test, focused pytest reports `7 passed` in
`test_boundary_linear_kuranishi.py`, `git diff --check` is clean on
touched files, and label-count greps find exactly one occurrence of each
new active equation label. Vol II `make fast` converged after two
passes with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A199 repairs the bordered configuration-space diagonal model.  The
active ordered KD core now states
`thm:bordered-diagonal-bicomodule-chain-model`, with licensing data
\(\alpha+\beta+\gamma\).  It defines
\[
D_A^{p,q}
=
C_\bullet^{\log}\bigl(\overline{\mathrm{Conf}}_{p,q}(\mathbb H)\bigr)
\otimes A^{\otimes p}\otimes A^{\otimes q},
\]
with divisor types \(z_i\to z_j\), \(x_a\to x_b\), and
\(\Im z_i\to0,\Re z_i\to x_a\).  The diagonal boundary condition is
\(\Delta_A={}_AA_A\), with left action from the negative boundary
approach and right action from the positive boundary approach.

The residue differential is now printed as
\[
d_D=d_{\mathrm{int}}+d_\partial+d_{\mathrm{mix}},
\qquad
d_{\mathrm{mix}}
=
\sum_{S,I}\operatorname{Res}_{D_{S,I}^{\mathrm{mix}}}\otimes\nu_{S,I}.
\]
The annular closure is
\[
D_A^{\mathrm{ann}}
=
\bigoplus_{n\ge0}
C_\bullet^{\log}\bigl(
\overline{\mathrm{Conf}}_n(S^1\times\mathbb R_{\ge0})
\bigr)\otimes_{\mathbb Z/n}A^{\otimes n},
\]
and the cyclic coinvariants use the desuspended sign
\[
[s^{-1}a_1|\cdots|s^{-1}a_n]
\mapsto
(-1)^{\bar a_n(\bar a_1+\cdots+\bar a_{n-1})}
[s^{-1}a_n|s^{-1}a_1|\cdots|s^{-1}a_{n-1}],
\qquad
\bar a_i=|a_i|-1.
\]

The superseded ordered KD copy no longer marks this bordered geometric
realization as conjectural and no longer says annular closure is
pending.  The ordered KD compute engine now exposes
`bordered_diagonal_bicomodule_profile()`,
`annular_cyclic_rotation_sign()`, and
`annular_diagonal_closure_profile()`.  Tests guard the three stratum
types, boundary-approach convention, residue decomposition,
desuspended cyclic signs, annular closure, and active manuscript source.

A199 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched ordered KD engine/test, focused pytest reports `50 passed`
in `test_ordered_chiral_kd_engine.py`, `git diff --check` is clean on
touched files, and label-count greps find exactly one occurrence of each
new active label. Vol II `make fast` converged after two passes with
zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A200 repairs the punctured-disk versus circle scope.  The Hochschild
chapter now states `prop:punctured-disk-circle-rotation-gap`, including
\[
B^{D^*,\Sigma}_\bullet(\cA)
=
\bigoplus_{n\ge0}
C_\bullet(\Conf_n(D^*))\otimes_{\Sigma_n}\cA^{\otimes n},
\]
and the non-equivalence
\[
H_\bullet(\Conf_k(D^*))\not\simeq H_\bullet(\Conf_k(S^1))
\]
as operadic modules.  The gap classes are
\[
\theta_i=d\arg z_i,\qquad \omega_{ij}=d\log(z_i-z_j),
\]
and the rotation BV operator is
\[
\Delta_{\mathrm{rot}}
=
\sum_i\iota_{\partial_{\arg z_i}}.
\]

The former full \(D^*\to S^1\) quasi-isomorphism has been narrowed:
`prop:punctured-disk-S1-qiso` now has source
\(B^{D^*,\mathrm{ang}}_\bullet(\cA)\), the angular \(E_1\)-quotient.
The raviolo restriction summary now uses
\(C^{\mathrm{ch,ang}}(D^\times,\cA)\to C^{\mathrm{top}}(S^1,\cA)\)
and explicitly retains the full relative \(E_2/E_1\) complex.

The new gap oracle exposes `rotation_bv_operator()`,
`punctured_disk_circle_gap_profile()`, and
`angular_projection_effect()`.  Tests guard the BV formula, angular
quotient scope, non-equivalence of full operadic modules, annular trace
versus chiral centre with rotation, projection preserves/forgets data,
and manuscript source formulas.

A200 verification: `compute/.venv/bin/python -m py_compile` passes on
the new module/test, focused pytest reports `7 passed` in
`test_punctured_disk_circle_gap.py`, and `git diff --check` is clean on
touched files. Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A201 repairs the Type-A Baxter--Rees RTT compactification surface.  The
appendix now prints the component RTT coordinates
\[
T_{ij}(u)=\delta_{ij}+\sum_{r\ge1}T_{ij}^{(r)}u^{-r},
\qquad
T(u)=\sum_{i,j}E_{ij}\otimes T_{ij}(u),
\]
and the component relation
\[
(u-v)[T_{ij}(u),T_{kl}(v)]
=
\hbar\bigl(T_{kj}(u)T_{il}(v)-T_{kj}(v)T_{il}(u)\bigr).
\]

The former boundary-deformation remark no longer leaves the deformation
complex to later work.  It has been replaced by
`thm:typeA-rtt-baxter-rees-obstruction-tower`, with licensing tags
\(\alpha+\beta+\gamma+\delta\), the Rees algebra
\[
R_\beta Y_\hbar(\mathfrak{gl}_N)
=
\bigoplus_{d\ge0}\beta^dF_dY_\hbar(\mathfrak{gl}_N),
\]
the formal family
\[
Y_{\mathrm{Bax}}^{\mathrm{RTT}}
=
\operatorname{Spf}\widehat{R_\beta Y_\hbar(\mathfrak{gl}_N)}
\to
\operatorname{Spf}\mathbb C[[\beta]],
\]
finite weight-window continuations
\[
R_{V,W}^{(\le w)}(u,\beta)
=
\sum_{d=0}^{D(w)}\beta^dR_d^{(\le w)}(u),
\qquad
\Theta_\Phi^{(\le w)}(\beta)
=
\sum_{d=0}^{D'(w)}\beta^d\Theta_d^{(\le w)},
\]
and the RTT deformation-complex equations
\[
d_{\Theta_0}\dot\Theta=0,
\qquad
\operatorname{Ob}_2(\dot\Theta;\Theta_2)
=
\frac12[\dot\Theta,\dot\Theta]+d_{\Theta_0}\Theta_2.
\]
Thus "boundary tensor geometry" now means the explicit obstruction tower
in \(C^\bullet_{\mathrm{RTT}}\), not an additional primitive structure.

The new oracle `compute/lib/typeA_baxter_rees_obstruction.py` exposes
`rtt_component_relation()`, `baxter_rees_family_profile()`,
`weightwise_continuation_profile()`,
`rtt_obstruction_tower_profile()`, and `beta_polynomial_terms()`.
Tests guard the RTT component identity, the Rees generic/special
fibres, finite weight-window \(\beta\)-polynomials, tangent and second
obstruction equations, source labels, and absence of the old deferred
deformation-complex sentence.  A wave-13 IV decorator now covers
`thm:typeA-rtt-baxter-rees-obstruction-tower`.

A201 verification: `compute/.venv/bin/python -m py_compile` passes on
the new module/test and IV file, focused pytest reports `8 passed`
across the new oracle and the two relevant wave-13 IV tests, and
`git diff --check` is clean on touched files. Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.

A202 repairs the shifted RTT rank-one test.  The shifted-RTT appendix
now states `thm:shifted-rtt-pairing-annihilator-kleinian-test`, with
the pairing-annihilator ideal
\[
I_\mu
=
\{x\in \widehat Y_{\hbar}^{\mathrm{RTT}}:
\langle x,y\rangle=0\ \forall y\in
Y_{\hbar}^{\mathrm{RTT},\vee}_{\ge\mu}\},
\]
and the exact coideal criterion
\[
\Delta(I_\mu)\subset
I_\mu\widehat\otimes\widehat Y_{\hbar}^{\mathrm{RTT}}
+
\widehat Y_{\hbar}^{\mathrm{RTT}}\widehat\otimes I_\mu .
\]
The chamber-shifted currents are
\[
T_{ij}^{(\mu)}(u)
=
u^{-\mu_i+\mu_j}
\left(\delta_{ij}+\sum_{r\ge1}T_{ij}^{(r)}u^{-r}\right),
\]
and boundary data enter through
\[
\operatorname{qdet}T(u)
=
\sum_{\sigma\in S_N}(-1)^\sigma
T_{\sigma(1),1}(u)\cdots T_{\sigma(N),N}(u-(N-1)\hbar)
=P_\mu(u).
\]

The rank-one determinantal/Casimir test is now the general
\(A_{m-1}\) quotient
\[
A_{\hbar}^{(m)}
=
\mathbb C\langle x,y,z\rangle/
([z,x]-\hbar x,\ [z,y]+\hbar y,\
yx-\prod_{a=0}^{m-1}(z-a\hbar),\
xy-\prod_{a=1}^{m}(z+a\hbar)),
\]
with
\[
\operatorname{gr}A_{\hbar}^{(m)}
\cong
\mathbb C[x,y,z]/(xy-z^m).
\]
Thus a proposed shifted line algebra is rejected if it fails the coideal
condition, kills the root sector, or fails this Kleinian associated
graded.

The new oracle `compute/lib/shifted_rtt_kleinian_test.py` exposes
`shifted_generator_profile()`, `pairing_annihilator_profile()`,
`quantum_determinant_profile()`, `kleinian_boundary_relations()`, and
`shifted_rtt_candidate_passes_rank_one_test()`.  Tests guard the shifted
exponent, pairing-annihilator/coideal criterion, quantum determinant
term count and shifts, \(A_{m-1}\) product factors, associated graded,
candidate rejection, invalid inputs, and manuscript source labels.  A
wave-16 IV decorator now covers
`thm:shifted-rtt-pairing-annihilator-kleinian-test`.

A202 verification: `compute/.venv/bin/python -m py_compile` passes on
the new module/test and IV file, focused pytest reports `8 passed`
across the new oracle and wave-16 IV test, and `git diff --check` is
clean on touched files. Vol II `make fast` converged after two passes
with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A203 repairs the strong-versus-weak chiral Yangian comparison.  The
line-operator chapter now states
`thm:strong-weak-chiral-yangian-comparison-obstruction`.  The weak datum
is
\[
Y_{\mathrm{weak}}(A)=(A,R(z),\otimes_z),
\]
with the Yang--Baxter equation.  The strong datum is
\[
Y_{\mathrm{str}}^{\mathrm{ch}}(A)
=
(A,\Delta_z,R(z),S(z),\epsilon,\{m_k\}_{k\ge2},
\Theta_{\mathrm{mod}}),
\]
with
\[
d\Theta_{\mathrm{mod}}
+\frac12[\Theta_{\mathrm{mod}},\Theta_{\mathrm{mod}}]=0.
\]

The theorem now prints the required comparison map
\[
\Psi_{\mathrm{RTT}\to\mathrm{mod}}\colon
\operatorname{RTTDef}(A,R)\to
\operatorname{MC}_{\mathrm{mod}}(A),
\]
its tangent condition
\[
d\Psi_{\mathrm{RTT}\to\mathrm{mod}}(\dot R(z))=\Theta_1,
\]
and the second obstruction
\[
\operatorname{ob}_2(\dot R)
=
\left[
d_{\Theta_0}\Theta_2+\frac12[\Theta_1,\Theta_1]
\right]
\in H^2_{\mathrm{mod}}(A).
\]
Thus "monodromy equals the \(R\)-matrix" is a strong modular-MC theorem
only when the KZ/KZB monodromy functor lands through this comparison and
the obstruction classes vanish; otherwise it remains a weak
RTT/factorization statement.

The existing two-definition remarks in the ordered KD core,
dg-shifted factorization bridge, and HT physical-origins chapter now
point to this theorem.  The new oracle
`compute/lib/strong_weak_chiral_yangian.py` exposes
`weak_chiral_yangian_profile()`, `strong_chiral_yangian_profile()`,
`comparison_map_profile()`, and
`monodromy_rmatrix_statement_status()`.  Tests guard the weak and strong
data, comparison map, tangent image, obstruction formula, monodromy
scope, source labels, and propagation references.  A wave-16 IV
decorator now covers
`thm:strong-weak-chiral-yangian-comparison-obstruction`.

A203 verification: `compute/.venv/bin/python -m py_compile` passes on
the new module/test and IV file, focused pytest reports `7 passed`
across the new oracle and wave-16 IV test, and `git diff --check` is
clean on touched files. Vol II `make fast` converged after two passes
with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.

A204 repairs the super-Yangian Berezinian bridge surface.  The
super-Yangian chapter now states
`def:canonical-berezinian-bridge-datum`, including the super-permutation
and rational super \(R\)-matrix
\[
P_s(v\otimes w)=(-1)^{\parity v\,\parity w}w\otimes v,
\qquad
R(u)=1-\frac{\hbar P_s}{u},
\]
the super-RTT relation, and the graded matrix multiplication law
\[
(E_{ij}\otimes a)(E_{kl}\otimes b)
=
(-1)^{(\parity a+\parity i+\parity j)(\parity k+\parity l)}
\delta_{jk}E_{il}\otimes ab.
\]

The missing bridge is now explicit as a chain-level pairing
\[
\langle-,-\rangle_{\mathrm{can}}\colon A\otimes A^!\to\mathbb C
\]
with
\[
\operatorname{Tr}_{A/A_{\nil}}\langle a,a^!\rangle_{\mathrm{can}}
=\str(aa^!),\qquad
\operatorname{Ber}\bigl(\langle e_i,e^!_j\rangle_{\mathrm{can}}\bigr)
=
\operatorname{Ber}(A),
\]
and
\[
\langle m_k(a_1,\dots,a_k),b\rangle_{\mathrm{can}}
=
(-1)^\epsilon
\langle a_1\otimes\cdots\otimes a_k,\Delta_k b\rangle_{\mathrm{can}},
\]
where
\[
\epsilon
=
\parity b(\bar a_1+\cdots+\bar a_k)
+
\sum_{i<j}\bar a_i\bar a_j
\quad\in\mathbb Z/2.
\]
The new `prop:berezinian-bridge-criterion` states the formal result:
if this datum exists, the supertrace and Berezinian normalisations are
two shadows of one chain-level duality; without it they remain parallel
scalar identities.

The F13 engine now exposes `super_permutation_sign()`,
`graded_matrix_product_sign()`, `mk_delta_compatibility_sign()`, and
`canonical_berezinian_bridge_profile()`.  Tests guard the odd-odd
super-permutation sign, graded matrix product sign, \(m_k/\Delta_k\)
compatibility sign, bridge profile, source labels, and an independent IV
sign check for `prop:berezinian-bridge-criterion`.

A204 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched F13 engine/test and super IV file, focused pytest reports
`25 passed`, and `git diff --check` is clean on touched files. Vol II
`make fast` converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.

A205 repairs the half-space BV reflected-weight mechanism.  The affine
half-space BV chapter now states
`thm:image-charge-reflected-arnold-cancellation`, with
\[
P_H(x,y)
=
P_{\mathbb R\times\mathbb C}(x,y)
-
\sigma P_{\mathbb R\times\mathbb C}(x,\bar y),
\qquad
\sigma(\phi)=(-1)^{w_\phi}\phi.
\]
The boundary OPE coefficient is
\[
\operatorname{OPE}_{\partial}(\phi_i,\phi_j)
=
\operatorname{OPE}_{\mathrm{bulk}}(\phi_i,\phi_j)
-
(-1)^{w_j}\operatorname{OPE}_{\mathrm{bulk}}(\phi_i,\sigma\phi_j),
\]
and the reflected weights obey
\[
w_{\mu\circ\nu}
\equiv
w_\mu+w_\nu+\operatorname{codim}D_{\mu,\nu}
\pmod2.
\]

The first non-affine obstruction is now printed as
\[
\mathcal O_2
=
\operatorname{Res}_{D_{12}}
[P(x_1,x_2)-\sigma P(x_1,\bar x_2)]\otimes\pi_2,
\]
and the cubic obstruction as
\[
\mathcal O_3
=
\sum_{\mathrm{cyc}}\operatorname{Res}_{D_{ij}D_{jk}}
P_H(x_i,x_j)P_H(x_j,x_k)\otimes[\pi_2,\pi_2].
\]
Under the doubling package,
\[
\mathcal O_3=0
\]
by the reflected Arnold relation
\[
\omega^H_{12}\wedge\omega^H_{23}
+\omega^H_{23}\wedge\omega^H_{31}
+\omega^H_{31}\wedge\omega^H_{12}=0.
\]

The affine half-space BV engine now exposes
`image_charge_propagator_profile()`,
`reflected_weight_composition()`, `reflected_obstruction_profile()`,
and `reflected_arnold_relation_profile()`.  Tests guard the
image-charge profile, reflected-weight composition, obstruction
formulas, Arnold cancellation, source labels, and a wave-14 IV witness
for `thm:image-charge-reflected-arnold-cancellation`.

A205 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched affine engine/test and IV file, focused pytest reports
`6 passed`, and `git diff --check` is clean on touched files. Vol II
`make fast` converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.

A206 repairs the irregular KZB boundary formal-type surface.  The
modular Swiss-cheese chapter now contains
`prop:irregular-kzb-q-formal-type-stokes-cocycle`, which writes the
ramified \(q\)-coordinate KZB connection as
\[
\nabla
=d-\left(
\frac{A_m}{q^m}+\frac{A_{m-1}}{q^{m-1}}+\cdots+
\frac{A_1}{q}+A_0
\right)dq
-\sum_\alpha B_\alpha(q)\,dt_\alpha.
\]
The formal type is
\[
\Theta_\partial(q)
=
\sum_{r=1}^{m-1}\frac{A_{r+1}}{r q^r}+A_1\log q,
\]
and the JMU/Levelt--Turrittin normal form is
\[
\nabla_q\simeq d_q-d_q\Theta_\partial-\Lambda\,\frac{dq}{q}
\]
after \(G(q)\in\operatorname{Aut}(V)[[q^{1/N}]]\).

The sectorial solution and Stokes transition are now printed as
\[
Y_\ell(q)=H_\ell(q)e^{\Theta_\partial(q)}q^\Lambda,\qquad
S_{\ell,\ell'}=Y_\ell^{-1}Y_{\ell'}.
\]
At a clutching corner,
\[
\Theta_{\partial,\Gamma_1\#\Gamma_2}
=
\Theta_{\partial,\Gamma_1}\oplus
\Theta_{\partial,\Gamma_2}\oplus
\Theta_{\mathrm{node}},
\qquad
S_{\ell,\ell'}^{\Gamma_1\#\Gamma_2}
=
S_{\ell,\ell'}^{\Gamma_1}
S_{\ell,\ell'}^{\Gamma_2}
S_{\ell,\ell'}^{\mathrm{node}},
\]
with order fixed by tangential basepoints.  Generic nonrational
associativity is the Stokes cocycle
\[
S_{\ell_1,\ell_2}S_{\ell_2,\ell_3}S_{\ell_3,\ell_1}=1.
\]

The generic-level theorem now points to this \(q\)-normal form.  The
curved-Dunn chapter now states explicitly that this KZB input supplies
boundary composition only after formal type and mixed Stokes data; it is
not the curved-Dunn \(H^2\)-acyclicity mechanism.

The new irregular KZB engine exposes `boundary_connection_profile()`,
`formal_type_profile()`, `formal_gauge_profile()`,
`stokes_sector_profile()`, `clutching_profile()`,
`covered_level_locus()`, `uncovered_level_locus()`, and
`kzb_vs_curved_dunn_profile()`.  Tests guard the \(q\)-normal form, JMU
gauge, Stokes cocycle, clutching identities, exact level-locus rows, the
uncovered rational nonintegral noncritical row, manuscript labels, and a
new IV witness for `prop:irregular-kzb-q-formal-type-stokes-cocycle`.

A206 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched irregular KZB engine/test and IV file, focused pytest
reports `7 passed` for `test_irregular_kzb_stokes.py` and `3 passed, 11
deselected` for `test_climax_theorems_iv.py -k kzb`, and `git diff
--check` is clean on touched files.  Full worktree whitespace scanning
is polluted by generated `out/main.log`.  Vol II `make fast` converged
after two passes with zero undefined citations and references, and Vol
II `make verify-licensing` reported zero blocking violations and zero
warnings.

A207 repairs the anomaly-completed Steinberg transgression relation.
The core anomaly-completion chapter now distinguishes the full
noncommutative Ore transgression from its central-shadow specialization:
\[
B_{\Theta}^{\mathrm{Ore}}
=
\frac{B\langle \eta\rangle}
{\left\langle
\eta b-(-1)^{|b|}b\eta-\iota_{\Theta}(b)
\mid b\in B
\right\rangle},
\qquad |\eta|=1,\qquad d\eta=\Theta.
\]
When \(\iota_\Theta=0\) and
\(\Theta\in Z^2(B)\cap Z(B)\), this reduces to the central-shadow
transgression used in the abelian computations.

The \(SU(2)\) anomaly-completed Steinberg construction now uses this
Ore-contraction relation in Step 4 and keeps the central-shadow
specialization explicit.  The summary table records the boundary
algebra \(V_k(\mathfrak{sl}_2)\), dual level \(-k-4\), modular
characteristic
\[
\kappa=3(k+2)/4,
\]
the anomaly class \(\Theta=\kappa\omega_1\), geometric source
\[
H^3(SU(2);\mathbb Z)\cong\mathbb Z,\qquad H_3=k c_3,
\]
secondary anomaly \(u=\eta^2\), and the genus-Clifford dichotomy
\[
\mathfrak G_g(B_\Theta)[u^{-1}]
\cong
\operatorname{Mat}_{2^g}(B_\Theta[u^{-1}]),
\qquad
\mathfrak G_g(B_\Theta)/(u)
\cong
(B_\Theta/(u))\otimes
\Lambda(\alpha_1,\beta_1,\dots,\alpha_g,\beta_g).
\]
The nonabelian gerbe/string-structure remark already states that
\(d\eta=\Theta\) is the cochain-level witness of the string
trivialization and \(u=\eta^2\) is the quadratic class of that
trivialization.

The anomaly engine now exposes `ore_transgression_relation()` and
`su2_anomalous_steinberg_profile()`.  Tests guard the \(\iota_\Theta\)
term, the \(\iota_\Theta=0\) specialization, the source labels,
\(\kappa=3(k+2)/4\), \(k c_3\), \(u=\eta^2\), and the
Morita/exterior dichotomy.

A207 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched anomaly engine/test files, focused pytest reports
`19 passed, 21 deselected` for
`test_anomaly_completed_engine.py -k 'ore or su2 or Transgression'`, and
`git diff --check` is clean on touched files.  Vol II `make fast`
converged after two passes with zero undefined citations and references,
and Vol II `make verify-licensing` reported zero blocking violations and
zero warnings.

A208 repairs the original-complex class-\(\mathsf M\)
Kac--Shapovalov criterion.  The class-\(\mathsf M\) topologisation
chapter now defines the Virasoro KS-original tempered locus by the
Verma basis
\[
L_{-\lambda}|h\rangle
=
L_{-\lambda_1}\cdots L_{-\lambda_\ell}|h\rangle,\qquad |\lambda|=n,
\]
the Gram matrix
\[
G_n(c,h)_{\lambda,\mu}
=
\langle h|
L_{\mu_m}\cdots L_{\mu_1}
L_{-\lambda_1}\cdots L_{-\lambda_\ell}
|h\rangle,
\]
summability
\[
\sum_{k\ge2}\frac{\|m_k\|_{\mathrm{KS},n}}{k!}<\infty
\quad\text{for every }n,
\]
and finite propagation
\[
\forall n\ \exists K(n):
\pi_{\le n}m_k|_{A_{\le n}^{\otimes k}}=0\quad(k>K(n)).
\]

The new theorem `thm:virasoro-ks-original-complex-criterion` records
the Kac determinant
\[
\det G_n(c,h)
=
C_n\prod_{\substack{r,s\ge1\\rs\le n}}
(h-h_{r,s}(c))^{p(n-rs)}
\]
and the finite-window transfer estimate
\[
\|m_k\|_{\mathrm{KS},n}
\le
C^k
\max_{|\lambda|\le n}\|G_{|\lambda|}(c,h)^{-1}\|\,
P_k(n,c,h).
\]
The theorem identifies the raw original-complex locus with the
summability-plus-finite-propagation condition.  Outside that locus the
class-\(\mathsf M\) statement remains in the KS-\(\rho\) Banach or
weight-completed ambient; generic central charge alone does not imply
raw direct-sum closure because the quartic obstruction survives and
higher transferred operations can feed lower finite weights through
inverse Shapovalov propagators.

The new KS oracle `virasoro_ks_original_complex.py` exposes
`verma_basis_vector()`, `gram_matrix_profile()`,
`kac_determinant_profile()`, `mk_norm_bound_profile()`,
`finite_propagation_profile()`, and
`original_tempered_locus_profile()`.  Tests guard the Verma basis,
Gram formula, Kac determinant shape, inverse-Gram norm source,
finite-propagation condition, raw/completed ambient split, and source
labels for the new theorem.

A208 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched KS engine/test files, focused pytest reports `3 passed, 4
deselected` for
`test_topologization_class_m_original_complex.py -k 'ks or Kac or source or criterion'`,
and `git diff --check` is clean on touched files.  Vol II `make fast`
converged after two passes with zero undefined citations and references,
and Vol II `make verify-licensing` reported zero blocking violations and
zero warnings.

A209 repairs the logarithmic triplet \(W(p)\) amplitude estimate.  The
triplet chapter now separates the semisimple regular channel split from
the full logarithmic shadow:
\[
S_r(\cW(p))
=
S_r^{\mathrm{reg}}(\cW(p))
+S_r^{\mathrm{nil}}(\cW(p))
+S_r^{\log}(\cW(p)),
\qquad
S_r^{\mathrm{reg}}=S_r^{TT}+S_r^{TW}+S_r^{WW}.
\]
The new conjecture
`conj:wp-logarithmic-boundary-amplitude-bound` records the missing
boundary-changing estimate
\[
\left|\left\langle
\phi_{0,1}(z_1)\cdots\phi_{0,1}(z_r)
\right\rangle\right|
\le
C^r(r!)^{1-\varepsilon}
\prod_{i<j}|z_i-z_j|^{-N_{ij}}
\left(1+\sum_{i<j}|\log|z_i-z_j||\right)^M,
\]
and its shadow consequence
\[
|S_r^{\log}(\cW(p))|\le C^r(r!)^{1-\varepsilon}.
\]
Finite-dimensional Zhu algebra and \(C_2\)-cofiniteness are explicitly
not sufficient for this estimate; they control finite graded fibres, not
arity growth through powers of \(\log(z_i-z_j)\).  The full
`conj:tempered-stratum-contains-wp` now depends on both the regular
\(TW/WW\) amplitude bound and this logarithmic boundary estimate.

The new oracle `wp_logarithmic_boundary_amplitude.py` records the
reg/nil/log decomposition, the \(\phi_{0,1}\) correlator profile, and
the Stirling consequence of a subfactorial logarithmic bound:
\[
\log\left(\left(C^r(r!)^{1-\varepsilon}/r!\right)^{1/r}\right)
=
\log C-\varepsilon\log(r!)/r\to-\infty.
\]
The F14 Adamovi\'c--Milas compute surface was downgraded to regular
semisimple channel constants and finite-\(r\) checks; it no longer reads
as a proof of full \(W(p)\) tempering.

A209 verification: `python3 -m py_compile` passes on the touched W(p)
compute/test files, focused pytest reports `17 passed` for
`test_logarithmic_wp_tempered.py -k 'logarithmic or manuscript or pole_envelope'`
and `5 passed, 5 deselected` for
`test_f14_w_p_triplet_amplitude.py -k 'five_class or amplitude_constants or triangle or stirling'`.
Vol II `make fast` converged after two passes with zero undefined
citations and references, Vol II `make verify-licensing` reported zero
blocking violations and zero warnings, and `git diff --check` is clean
on touched files.

A210 repairs the irrational-coset criterion.  The phrase VSKR + BGG no
longer means "bounded generator-to-generator OPE pole order."  It now
means Virasoro sub-channel in the Kac-regular locus plus the actual
BGG/Verma-shadow resolution
\[
\cdots\to\bigoplus_{\ell(w)=2}M(w\!\cdot\!\lambda)
\to\bigoplus_{\ell(w)=1}M(w\!\cdot\!\lambda)
\to M(\lambda)\to L(\lambda)\to0
\]
and the alternating coset shadow coefficient
\[
S_r^{\mathrm{coset}}(\lambda)=
\sum_{w\in W_{\mathrm{aff}}}(-1)^{\ell(w)}
S_r(M(w\!\cdot\!\lambda))
q^{\Delta(w\!\cdot\!\lambda)-\Delta(\lambda)}.
\]
The required analytic input is the Kazhdan-ramification lower bound
\[
\operatorname{Re}(\Delta(w\!\cdot\!\lambda)-\Delta(\lambda))
\ge a\ell(w)^2-b,
\]
which turns the affine Weyl length sum into a Gaussian-weighted
polynomial envelope.  Finite local OPE pole order remains a local
input to the estimate, but it is not the BGG criterion and does not
replace the ramification inequality.

The theorem `thm:irrational-coset-tempered` now states the proved
conditional implication: VSKR plus the BGG ramification clause implies
\((H_S)\).  The parafermion radius, Virasoro-in-affine radius,
infinite-Zhu witness, and non-rational affine-minimal corollary are
marked conditional where their tempering conclusion uses that clause.
The summary surfaces no longer claim that every irrational coset is
tempered unconditionally or that finite Zhu/C_2-cofiniteness is the
sharp criterion.

The oracle `irrational_coset_tempered_engine.py` now records the BGG
resolution profile, ramification lower bound, Gaussian exponent, and a
sufficient comparison showing Gaussian length decay dominates
polynomial affine-Weyl growth.  `TemperingCertificate` includes the
field `bgg_ramification_assumed`, and the tests assert that the
criterion uses the BGG alternating sum rather than Zhu dimension,
central charge, or bounded OPE pole order as a proxy.

A210 verification: `python3 -m py_compile` passes on
`compute/lib/irrational_coset_tempered_engine.py` and
`compute/tests/test_irrational_cosets_tempered.py`; focused pytest
reports `5 passed, 5 deselected` for
`test_irrational_cosets_tempered.py -k 'bgg or ramification or refined or certificate or rho'`.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.

A211 repairs the finite-orbifold descent datum for Monster and
Schellekens rows.  The Monster/Schellekens chapters already separated
VOA presentation, level matching, Dijkgraaf--Witten/BV anomaly, and
finite-orbifold BV descent, but the actual descent object was not
displayed.  The live surface now records the homotopy fixed point
totalisation
\[
V^{hG}
=
\Tot\!\left[
V\rightrightarrows\Map(G,V)\to\Map(G^2,V)\to\cdots
\right],
\]
together with the obstruction class
\([\alpha_{\mathrm{orb}}]\in H^3(BG;\mathrm U(1))\).  Chain-level
\(E_3\) descent requires this class to vanish.  If
\(\alpha_{\mathrm{orb}}=d\beta\), the associator is corrected by
\[
\Phi^\beta_{g,h,k}
=
\beta(g,h)\beta(gh,k)\beta(h,k)^{-1}\beta(g,hk)^{-1}
\Phi_{g,h,k}.
\]
The Leech \(\mathbb Z/2\) Monster row remains conditional on the local
sign \(\alpha_{\mathrm{orb}}(\sigma,\sigma,\sigma)=+1\), and the
Schellekens rows remain stratified as
\[
71=24_{\mathrm{Niemeier}}+1_{V_1=0}+46_{\mathrm{nonlat}},
\]
with cyclic-orbifold level matching furnishing the VOA presentation but
not the BV trivialization.

The new oracle `finite_orbifold_descent.py` records the homotopy fixed
point profile, the DW obstruction group, the \(\beta\)-corrected
associator formula/value, and the Schellekens row partition.  The
Monster and Schellekens tests now assert that non-trivial finite
orbifold strata use \(V^{hG}\), require a zero DW class, and keep the
associator correction separate from determinant signs and level
matching.

A211 verification: `python3 -m py_compile` passes on
`compute/lib/finite_orbifold_descent.py`, `compute/lib/z2_group_cohomology.py`,
`compute/tests/test_monster_chain_level_e3_top.py`, and
`compute/tests/test_schellekens_71_alpha_classification.py`; focused
pytest reports `14 passed` for the two Monster/Schellekens test files.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.

A212 repairs the six-dimensional holomorphic Chern--Simons anomaly
surface.  The live chapters already contained the 6d hCS anomaly
package, but the chiral-avatar theorem stated the QME before displaying
the obstruction, and one shadow-tower remark could be read as if the
algebraic shadow constructed the Green--Schwarz/axion cancellation.  The
theorem now displays the degree-one local BV anomaly
\[
\mathcal A^{(1)}
=
\int_Y \Omega_Y\wedge
\Tr_{\mathfrak g}(A\,\partial A\,\partial A\,\partial A)
\in
H^1_{\mathrm{loc}}(\Obs^{\mathrm{cl}}(Y),Q+\{S,-\}),
\]
equivalently the invariant quartic class
\[
[\mathrm{anom}_1]=\int_Y\Tr_{\mathrm{ad}}A(F_A)^3.
\]
The QME holds only on the cancellation locus.  The cancellation theorem
now states the scalar quartic trace identity
\[
\Tr_{\mathfrak g}(X^4)=
\lambda_{\mathfrak g}(\Tr_{\mathfrak g}(X^2))^2,
\]
with Green--Schwarz/axion counterterm data where the theorem requires
them.  The \(K3\) factor remains the Euler-class integration
\[
\frac1{24}\int_{K3}c_2(TK3)=\frac{24}{24}=1.
\]
The main summary and preface now say explicitly that the shadow tower
detects the residual quartic class but does not itself construct the
BV cancellation.

The new oracle `six_d_hcs_anomaly.py` records the local anomaly
profile, the \(K3\) Euler factor, the Deligne quartic-identity locus,
and the shadow-versus-gauge-cancellation distinction.  Tests guard the
exact local formula, the \(1/24\cdot24\) factor, the Deligne rigid
locus versus \(E_6\) and \(A_2\), and the non-construction role of the
shadow tower.

A212 verification: `python3 -m py_compile` passes on
`compute/lib/six_d_hcs_anomaly.py` and
`compute/tests/test_six_d_hcs_anomaly.py`; focused pytest reports
`4 passed` for `test_six_d_hcs_anomaly.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.

A213 repairs the K3 protected-Pfaffian operator surface.  The
manuscript already treated CP1 as conditional and already rejected
promotion from the scalar \(\Delta_5\) identity to an acting boundary
object, but several high-level surfaces still named only
\(\operatorPrim{X}\).  The live theorem now introduces the actual
Dirac--Igusa operator
\[
\mathfrak D_X\in \operatorname{End}_{\mathrm{ChirHoch}}(H_X)
\]
and the protected Pfaffian section
\[
\protectedPfaff{\mathfrak D_X}\in
H^0(\Lambda^{2,1}_{\mathrm{II}},\mathcal L^5\otimes\nu_{\Delta_5})
\]
with automorphic comparison
\[
\iota_{\mathrm{aut}}(\protectedPfaff{\mathfrak D_X})=\Delta_5.
\]
The Hall side is now labelled by
\[
I_{\mathrm{Hall}}\colon \CoHA(K3)\to D_{\mathrm{Hall}}(K3)
\to\mathfrak g_{\Delta_5},
\]
and the theorem includes the Pfaffian/Borcherds square
\[
\begin{tikzcd}
\ChirHoch^\bullet(\SpCh\PhiFA(K3)) \arrow[r,"\protectedPfaff"]
\arrow[d,"I_{\mathrm{Hall}}"'] &
H^0(\mathcal L^5) \arrow[d,"\mathrm{Borcherds}"]\\
\mathfrak g_{\Delta_5} \arrow[r,"\mathrm{den}"'] &
\mathbb C\cdot\Delta_5 .
\end{tikzcd}
\]
The scalar identity \(\Delta_5\) or \(\Delta_5^{-2}\) remains only the
character shadow until this square is present.

The P1 oracle now records the operator membership, Pfaffian section,
automorphic identity, \(I_{\mathrm{Hall}}\), and the
`Borcherds o Pf_prot = den o I_Hall` commutativity condition.  The P1
tests guard both the profile and the manuscript source so the theorem
cannot silently regress to a scalar statement.

A213 verification: `python3 -m py_compile` passes on the P1 and
Hall-residual helper/test files; focused pytest reports `5 passed, 19
deselected` for
`test_p1_protected_pfaffian_chain_level.py -k 'operator or hall_chiral or pfaffian_square or manuscript or scalar_no_promotion'`
and `5 passed` for `test_hall_borcherds_gravity_residual.py`.  Vol II
`make fast` converged after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reported zero
blocking violations and zero warnings.

A214 repairs the chiral Springer antipode package.  The frontier already
had a Steinberg involution and a Drinfeld-double antipode section, but it
did not isolate the Hopf-side complementary Steinberg carrier.  The live
surface now states
\[
\Stch(\cA)=\cL_{\cA}\times^h_{\Mvac}\cL_{\cA}^{\vee},
\]
with \(\cL_{\cA}^{\vee}\) the Verdier--Koszul complementary
Lagrangian.  The convolution is explicitly
\[
F\star G=p_{13,*}(p_{12}^{*}F\otimes p_{23}^{*}G)
\]
on
\[
\cL_{\cA}\times^h_{\Mvac}\cL_{\cA}^{\vee}
\times^h_{\Mvac}\cL_{\cA}.
\]
Orientation reversal is the transposition
\[
\iota(x,y)=(y,x)
\]
followed by Verdier duality:
\[
S_{\cA}(F)=\mathbb D_{\mathrm{Verdier}}(\iota^*F).
\]

The Hopf identities are now displayed in the Springer-antipode
conjecture:
\[
m(S_{\cA}\otimes 1)\Delta_{\cA}=\eta_{\cA}\varepsilon_{\cA},
\qquad
m(1\otimes S_{\cA})\Delta_{\cA}=\eta_{\cA}\varepsilon_{\cA},
\]
\[
\Delta_{\cA}S_{\cA}
=(S_{\cA}\otimes S_{\cA})\Delta_{\cA}^{\op},
\qquad
\varepsilon_{\cA}S_{\cA}=\varepsilon_{\cA}.
\]
The higher-class obstruction is the full-MC orientation-reversal
homotopy
\[
S_{\cA}(\Theta_{\cA})
=-\Theta_{\cA}+d_{\Theta_{\cA}}(\Xi).
\]
Without such a homotopy, the construction gives only an
associated-graded anti-equivalence of the Springer convolution
category, not the antipode of the full chiral Drinfeld double.  The
ordered-associative Drinfeld-double programme summary now carries the
same caveat.

The new oracle `chiral_springer_antipode.py` records the complementary
carrier, convolution correspondence, Verdier antipode, four Hopf
identities, and the class-\(\mathbf G\) versus higher-MC obstruction
split.  Tests guard both the finite profile and the manuscript source
strings.

A214 verification: `python3 -m py_compile` passes on
`compute/lib/chiral_springer_antipode.py` and
`compute/tests/test_chiral_springer_antipode.py`; focused pytest reports
`6 passed` for `test_chiral_springer_antipode.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean.

A215 repairs the Maloney--Witten/Page-Stokes overidentification at the
top-level synopsis.  The detailed gravity theorem already said the
Virasoro scalar bar trace is only the chain-level scalar seed and that
the Maloney--Witten modular orbit sum requires extra saddle data.  The
remaining stale surface was the Part VI synopsis in `main.tex`, which
still compressed the Page/de~Sitter statements too close to scalar
complementarity and Borel resummation.

The repaired statement is:
\[
\Tr_{\mathrm{Bord}}(\Vir_c)\exp(\Theta)
\]
is the perturbative thermal-AdS\(_3\)/BTZ seed after the
Brown--Henneaux dictionary and saddle labelling are supplied.  It is
not
\[
Z_{\mathrm{MW}}(\tau,\bar\tau)=
\sum_{\gamma\in\mathrm{SL}_2(\mathbb Z)/\Gamma_\infty}
Z_{\mathrm{thermal}}(\gamma\tau,\gamma\bar\tau).
\]
The scalar genus series has radius \(4\pi^2\) in \(\hbar^2\), and no
Page or BTZ Stokes jump follows from that scalar tower alone.  Page
time and de~Sitter entropy now live in the separate raw-transseries
lane under sectorial Borel summability, Stokes, modular invariance,
vacuum dominance, and saddle-extraction hypotheses.

The gravity oracle now includes `maloney_witten_scope_profile()`,
which records the scalar trace, the MW orbit sum, the false
\(\Phi_{\mathrm{hol}}\Rightarrow Z_{\mathrm{MW}}\) implication, the
\(4\pi^2\) scalar radius, the entire ordinary Borel transform, and the
extra Page/Stokes hypotheses.  Tests guard both the profile and the
`main.tex` synopsis.

A215 verification: `python3 -m py_compile` passes on
`compute/lib/gravity_3d_engine.py` and
`compute/tests/test_gravity_3d_engine.py`; focused pytest reports
`6 passed, 63 deselected` for
`test_gravity_3d_engine.py -k "ExactGravityScope or maloney or scalar_genus or main_synopsis"`.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.  Touched-file
`git diff --check` was clean.

A216 repairs the wild-boundary classification.  The Koszul atlas is now
stated as the scoped sequent
\[
A\in\operatorname{Kosz}_{\mathrm{ch}}
\Longrightarrow
r_{\mathrm{sh}}(A)\in\{2,3,4,\infty\}.
\]
The shifted algebraic coordinate is no longer the primitive invariant:
\[
d_{\mathrm{alg}}=r_{\mathrm{sh}}-2
\]
on the finite Koszul classes, while class \(\mathbf M\) has
\[
r_{\mathrm{sh}}=d_{\mathrm{alg}}=\infty .
\]
Class \(\mathbf W\) is not a refinement of class \(\mathbf M\).  It is
the outside-Koszul boundary, modelled algebraically by Kronecker
\(K_m\), \(m\ge3\), where
\[
A\notin\operatorname{Kosz}_{\mathrm{ch}},\qquad
E_r(B^{\mathrm{ord}}(A))\not\Rightarrow E_\infty
\quad\text{for every finite }r,
\]
and \(r_{\mathrm{sh}}(A)\) is undefined.

The introduction now states the actual open problem: Virasoro and
\(\cW_N\) already furnish class \(\mathbf M\), so the problem is not
the existence of infinite depth.  The open problem is the existence of
genuine logarithmic \(\SCchtop\)-algebras whose ordered bar spectral
sequence refuses the Koszul shadow altogether.  The later introduction
table was also corrected from the ambiguous notation \(d(\cA)\) to
the Koszul shadow-depth notation \(r_{\mathrm{sh}}(\cA)\), leaving
total depth \(d=1+d_{\mathrm{arith}}+d_{\mathrm{alg}}\) as a separate
derived quantity.

The new oracle `shadow_depth_atlas.py` records the four Koszul atlas
rows G/L/C/M, their \(r_{\mathrm{sh}}\) values \(2,3,4,\infty\), their
shifted \(d_{\mathrm{alg}}\) values \(0,1,2,\infty\), and the
non-Koszul wild boundary row with undefined shadow depth.  Tests guard
the exact values, the non-refinement of class \(\mathbf W\), the real
open problem, and the live source strings.

A216 verification: `python3 -m py_compile` passes on
`compute/lib/shadow_depth_atlas.py` and
`compute/tests/test_shadow_depth_atlas.py`; focused pytest reports
`4 passed` for `test_shadow_depth_atlas.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean, and the new untracked files plus this processed ledger had no
trailing whitespace.

A217 repairs the K3/Class-\(\mathcal S\) closure statement.  The
comparison at the K3 base point is not the ungated equivalence
\[
\operatorname{Schur}\!\left(T[A_1,\Sigma_{0,24}]\right)
\cong \mathbf H_{\Delta_5}.
\]
The theorem-level statement is the finite Hall/chiral gate:
\[
\operatorname{rad}_{\mathrm{Hall},N}/\operatorname{rad}_{N}=0,\qquad
D^{\mathrm{fin}}_{\mathrm{Hall}}\ \text{exists},
\]
\[
\operatorname{Borch}\circ\operatorname{Hall}
\ \text{is height-compatible},\qquad
\operatorname{Schur}\to\Zderch
\ \text{is compatible with the }\SCchtop\text{ trace}.
\]
The reduced compact Hall source, finite radical-quotient
Hall--Drinfeld doubles, finite Hall--Borcherds recognition maps, and
\(\SCchtop\) trace compatibility must be installed compatibly in
height.  Without these gates, the K3/Class-\(\mathcal S\) comparison
is a shadow comparison, not an object equivalence.

The top-level K3 bridge in `main.tex` now displays this gate before the
Dirac--Igusa protected-Pfaffian operator.  The class-\(\mathcal S\)
K3 remark and the K3/Borcherds protected-Pfaffian operator-square
theorem in `3d_gravity.tex` now state the same finite-height
conditions, so the scalar \(\Delta_5\) lane cannot be read as a global
construction of \(\mathbf H_{\Delta_5}\).

The P1 oracle now includes `k3_class_s_closure_gate_profile()`, which
records the false equivalence, the correct conditional status, the
finite theorem gates, and the four comparison blocks.  The existing
Hall/chiral-square profile now uses the A217 gate package as its
`requires` tuple.  Tests guard both the oracle and the manuscript
source strings.

A217 verification: `python3 -m py_compile` passes on
`compute/lib/p1_protected_pfaffian_chain_level.py` and
`compute/tests/test_p1_protected_pfaffian_chain_level.py`; full
focused pytest reports `25 passed` for
`test_p1_protected_pfaffian_chain_level.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean, and this processed ledger had no trailing whitespace.

A218 installs the deletion ledger demanded by the PDF.  The Universal
Holography master theorem now carries
`tab:universal-holography-deletion-ledger` immediately after the
theorem statement.  The table makes nine corrections explicit:
\(\SCchtop\) recognition is conditional on the H4 product-formal
local-shadow theorem; arbitrary-logarithmic bulk observables require a
chosen HT realisation plus H1--H2 and exact-sector hypotheses; class
\(\mathbf M\) raw direct-sum chain-level form is false and is replaced
by weight-completed/pro or finite-propagation ambients;
\(E_3\)-PBW concentration remains needs-verification with the
Orlik--Solomon/FM proof as the established local mechanism;
DS--Hochschild is not proved for all nilpotents; chain-level
associator-free mixed structure is false without a
\(\Phi\)-dependent representative; the bar scalar trace is not the
Maloney--Witten sum; the scalar genus tower does not determine the
full tensor channel; and the K3/Class-\(\mathcal S\) comparison is
conditional on the finite Hall--Borcherds gate.

The new oracle `deletion_ledger.py` records the nine claims, corrected
statuses, replacement forms, false-claim subset, and allowed maximal
theorem status alphabet:
\[
\{\mathrm{ProvedHere},\mathrm{ProvedElsewhere},
\mathrm{Conditional},\mathrm{Conjectured}\}.
\]
The source guard checks that the theorem-facing table stays present
and contains the required rows.

A218 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`4 passed` for `test_deletion_ledger.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean, and this processed ledger had no trailing whitespace.

A219 installs the corrected maximal theorem form.  The Universal
Holography master theorem now carries the sequent
\[
\Xi(A)\vdash_{\mathcal A(A)}
\bigl(\Phihol(A)=T_A,
\mathrm{Obs}^{\partial}(T_A)\simeq A,
\mathrm{Obs}^{\mathrm{bulk}}(T_A)\simeq
Z^{\mathrm{der}}_{\mathrm{ch}}(A)\bigr),
\]
with
\[
\sigma(A)\in
\{\ClaimStatusProvedHere,\ClaimStatusProvedElsewhere,
\ClaimStatusConditional,\ClaimStatusConjectured\}.
\]
The ambient is now explicit:
\[
\mathcal A(A)=
\begin{cases}
\mathrm{Ch}(\mathrm{Vect})&
A\in\{\mathbf G,\mathbf L,\mathbf C\},\\
\widehat{\mathrm{Ch}}_{\mathrm{wt},\rho}\ \text{or}
\mathrm{pro}\text{-}\mathrm{Ch}(\mathrm{Vect})&
A\in\mathbf M,\\
\varprojlim_N\mathcal A_{\le N}^{\mathrm{wt}}&
A=W_{\infty}[\lambda].
\end{cases}
\]
The pointwise theorem now says class \(\mathbf C\) is original-complex
on its stated non-critical contact loci; class \(\mathbf M\) is
weight-completed/pro; \(W_{\infty}[\lambda]\) is a bounded-weight
pro-window tower, not a uniform Banach algebra; and the raw direct-sum
class-\(\mathbf M\) statement is not a theorem.  The sequent remark
also states that arbitrary logarithmic \(\SCchtop\)-algebras have the
algebraic package but not the physical observable identification until
H1, H2, and exact-sector comparison are supplied; that \(k=-h^\vee\)
is outside the non-critical affine source; and that K3/Borcherds is
only scalar until the Hall--Borcherds operator gates are built.

The deletion-ledger oracle now includes
`corrected_maximal_theorem_form()`, which records the sequent, status
alphabet, ambient map, critical-affine exclusion, logarithmic
H1/H2/exact-sector gate, raw-direct-sum rejection, and K3/Borcherds
scalar-only caveat.  Tests guard both the oracle and the live source.

A219 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`6 passed` for `test_deletion_ledger.py`.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean, and this processed ledger had no trailing whitespace.

A220 propagates the deletion-ledger correction that good grading is
not DS--Hochschild transport.  The good-grading scope remark in
`axioms.tex` now says that a good grading is the input for forming the
Kazhdan DS-BRST complex, not a DS--Hochschild transport theorem, not a
chain-level \(\SCchtop\) lift, and not a Koszul-preservation
statement.  The deletion-ledger replacement now names principal, hook,
and named cover-descent fibres such as Bershadsky--Polyakov as
verified transport cases; subregular, minimal beyond named fibres,
exceptional good-graded, and exotic nilpotents require case-by-case
BRST admissibility, invariant transferred braces, and Hochschild
comparison.

The 3d-gravity DS theorem now marks its good-graded statement as
cohomological topologization; the stronger DS--Hochschild chain-level
transport and chain-level \(\SCchtop\) lift are separate assertions
requiring the principal, hook, or cover-descent hypotheses of
`thm:chd-ds-hochschild`.  The THQG boundary theorem now requires a
good DS grading and an applicable Costello--Gaiotto DS boundary
condition.  The two Bershadsky--Polyakov Kazhdan-denominator slips in
`e_infinity_topologization.tex` were corrected from \(d_f=1\) to
\(d_f=2\).

The deletion-ledger oracle now includes `ds_nilpotent_scope_profile()`,
which records the separation between DS-BRST existence and
DS--Hochschild transport, the verified transport cases, the
case-by-case nilpotent families, and the transport gates
non-critical level, DS HPL special deformation retract,
\(\mu_q\)-invariant transferred braces, \(\mu_q\)-invariant
transported antighost, and Hochschild comparison.  The DS-Hochschild,
E3-topological DS, climax wave-4, and FM81 tests now prevent a generic
good grading from implying DS--Hochschild transport or Koszul
preservation.

A220 verification: `python3 -m py_compile` passes on the touched
helper and tests; focused pytest reports `29 passed` for
`test_deletion_ledger.py`, `test_chd_ds_hochschild_iv.py`,
`test_e3_topological_ds_general.py`,
`test_e3_topological_ds_general_iv.py`,
`test_climax_theorems_wave4_iv.py`, and
`test_fm81_fractional_ghost.py`.  Residual broad-phrase greps are
clean; the remaining `DS--Hochschild all nilpotents` occurrences are
the intentional deletion-ledger row and guards.  Vol II `make fast`
converged after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reported zero blocking
violations and zero warnings.  Touched-file `git diff --check` was
clean.

A221 propagates the deletion-ledger correction that a chain-level
associator-free mixed structure is false, not a remaining open form.
The Part VIII synthesis chapter no longer calls the
associator-independent chain-level chiral Deligne--Tamarkin action a
formal question.  It now says the chain-level associator-independent
formulation is not a theorem, because the Drinfeld-associator choices
form a non-trivial \(\mathrm{GRT}_{1}(\mathbb Q)\)-torsor.  The
positive statement is separated: cohomology and the bar-side
invariants \(\kappa\), the shadow tower, and Koszul depth are
associator-free; the chain-level mixed object is a
\(\Phi\)-dependent representative.  The atlas summary now lists the
\(\Phi\)-dependent chain-level Deligne--Tamarkin torsor as the outside
obligation and explicitly says the associator-independent chain-level
collapse is false.

The deletion-ledger oracle now includes
`associator_chain_scope_profile()`, recording the split between
cohomological associator-independence and chain-level
\(\Phi\)-dependence.  The tests guard the oracle and the live Part
VIII and Higher-Deligne source strings.

A221 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`10 passed` for `test_deletion_ledger.py`.  Residual associator greps
leave only the intentional deletion-ledger row and the Koszulness
sentence that the Higher Massey tower obstructs associator-independent
chain-level chiral Deligne--Tamarkin.  Vol II `make fast` converged
after two passes with zero undefined citations and references, and
Vol II `make verify-licensing` reported zero blocking violations and
zero warnings.  Touched-file `git diff --check` was clean.

A222 propagates the deletion-ledger correction that the raw bar trace
is not the Maloney--Witten sum.  The off-Koszul bridge no longer
labels its third clause "Orbit sum = Maloney--Witten partition
function"; it now says "Conditional orbit-sum comparison with the
Maloney--Witten partition function."  The text separates the raw
chain-level Virasoro bar trace from its Borel--Zwegers completed seed
and from the physical modular saddle sum.  The Maloney--Witten object
requires Borel summability, Zwegers completion, Brown--Henneaux data,
saddle labelling, and an ensemble prescription.  The class-\(\mathsf G\)
Heisenberg case now states that its \(Z^{\mathrm{MW}}_{\mathrm{grav},
\mathcal H_k}\) notation names a completed Rademacher orbit object,
not a pure-gravity path integral without an additional saddle
interpretation.

The deletion-ledger oracle now includes
`maloney_witten_bridge_scope_profile()`, recording that the raw bar
trace is only a chain-level perturbative seed and listing the
additional data required to assemble a Maloney--Witten orbit sum.  The
tests guard both the oracle and the live `3d_gravity.tex` / climax
source strings.

A222 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`12 passed` for `test_deletion_ledger.py`.  Residual equality greps
leave only the intentional deletion-ledger row and guards.  Vol II
`make fast` converged after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reported zero
blocking violations and zero warnings.  Touched-file
`git diff --check` was clean.

A223 propagates the deletion-ledger correction that an abstract
logarithmic \(\SCchtop\)-algebra does not determine physical bulk
observables.  The Hochschild package-content remark now carries an
explicit scope: the datum exists only after a chosen HT
prefactorization realization supplies the product-formal local shadow
and the boundary-linear exact-sector hypotheses.  Without that
realization, the abstract logarithmic \(\SCchtop\)-algebra has only
the algebraic local-shadow package and does not determine a physical
bulk observable complex.

The deletion-ledger oracle now includes
`arbitrary_logarithmic_bulk_scope_profile()`, recording that the
bulk/Hochschild comparison requires chosen HT realization,
product-formal local shadow, H1--H2 physics bridge, and
boundary-linear exact-sector comparison.  Tests guard the oracle and
the live Hochschild/brace source strings.

A223 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`14 passed` for `test_deletion_ledger.py`.  Residual arbitrary-log
bulk greps leave only explicit open/scope warnings and the new guard.
Vol II `make fast` converged after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reported
zero blocking violations and zero warnings.  Touched-file
`git diff --check` was clean.

A224 propagates the deletion-ledger correction that \(\SCchtop\)
recognition is not global.  The example and bar-cobar surfaces no
longer refer merely to "the recognition theorem" or "recognition
conditions" where the theorem being used is the product-formal
local-shadow theorem.  They now say "product-formal local-shadow
recognition theorem" or "product-formal local-shadow recognition
conditions."  The Wave 17 independent-verification decorator now
derives from the product-formal local-shadow recognition theorem, not
from a bare programme recognition theorem for logarithmic
\(\SCchtop\)-algebras.

The deletion-ledger oracle now includes `recognition_scope_profile()`,
recording that there is no global recognition theorem: the recognized
surface is product-formal local-shadow rectangle data, not arbitrary
global Ran-space factorization data or a physical HT bulk without a
chosen prefactorization realization.  Tests guard the oracle and live
source strings.

A224 verification: `python3 -m py_compile` passes on the touched
helper and tests; focused pytest reports `36 passed` for
`test_deletion_ledger.py` and `test_climax_theorems_wave17_iv.py`.
Residual global-recognition greps leave only the intentional
deletion-ledger row and guards.  Vol II `make fast` converged after
two passes with zero undefined citations and references, and Vol II
`make verify-licensing` reported zero blocking violations and zero
warnings.  Touched-file `git diff --check` was clean.

A225 propagates the deletion-ledger correction that \(E_3\)-PBW is not
the proved Hochschild concentration mechanism.  The Vol II live
surfaces now separate the two routes.  Chiral PBW/Koszul collapse is
the filtration or associated-graded input.  The established vanishing
mechanism for Theorem H is the ordered-bar
Arnold--Orlik--Solomon/FM concentration proof.  The possible
chiral-\(E_3\)-PBW proof remains a conjectural second route requiring a
filtered \(E_3\) envelope, free \(E_3\) associated graded, convergent
PBW spectral sequence, polynomial-growth/amplitude bounds, and
\(E_1\)-page support in total degrees \(\leq2\).

The MC1 row in `chapters/connections/concordance.tex` now says
PBW/Koszul collapse is the input and ordered-bar
Arnold--Orlik--Solomon/FM concentration is the vanishing mechanism.
The line-operator hierarchy, anomaly perfectness proof, and Hochschild
HCA remark carry the same separation.  The active line-operator proof
now states that the PBW/Koszul concentration input alone does not force
all transferred \(A_\infty\)-operations to vanish.

The deletion-ledger oracle now includes
`pbw_concentration_scope_profile()`, recording the proved ordered-bar
route and the missing spectral-sequence clauses for the conjectural
chiral-\(E_3\)-PBW route.  Tests guard the oracle and the affected live
and superseded manuscript strings.

A225 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`21 passed` for `test_deletion_ledger.py` and
`test_chiral_higher_deligne.py`.  Vol II residual greps leave only the
intentional deletion-ledger row and guards.  Cross-volume grep of
Vol I and Vol III found no `E3-PBW` / `PBW proves concentration`
advertisement; remaining Vol I hits are named MC1/PBW-concentration
theorem surfaces, not the deleted \(E_3\)-PBW route.  Vol II
`make fast` converged after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reported zero
blocking violations and zero warnings.

A226 propagates the deletion-ledger correction that the scalar genus
tower does not determine the full tensor channel.  The scalar tower is
the uniform-weight closed-sector trace of the modular bar curvature.
It reconstructs the whole package only in the rank-one abelian
Heisenberg/Gaussian exception, where there are no mixed tensor entries.
Outside that lane the full tensor channel requires a chosen channel
splitting, scalar diagonalisation, a conformal-block basis,
off-diagonal stable-graph component integrals, and non-scalar ordered
or field-valued coefficients.

The Heisenberg Rosetta Stone now says the scalar \(k\) determines the
package only because the tensor channel is rank-one abelian, and that
this does not extend to multi-channel tensor claims.  The 3d-gravity
uniform-weight paragraph now says only closed scalar observables are
\(\kappa\)-polynomials and explicitly denies reconstruction of the full
tensor channel.  The movements table says "complete scalar genus
tower"; the holographic-reconstruction remark says \(\kappa\)
determines the rank-one Gaussian shadow, not the full tensor-channel
package; and the programme climax scalar partition seed is explicitly
not the full tensor-channel object.

The deletion-ledger oracle now includes
`scalar_tensor_channel_scope_profile()`, recording the rank-one
exception, the full tensor-channel requirements, and the data forgotten
by scalar trace.  Tests guard the oracle, the corrected manuscript
strings, and the existing W3 tensor-Arakelov source where the mixed
entries \(K_{TW}\) and \(K_{WT}\) are invisible after scalar trace.

A226 verification: `python3 -m py_compile` passes on
`compute/lib/deletion_ledger.py` and
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`26 passed` for `test_deletion_ledger.py` and
`test_wN_tensor_arakelov_weight_distribution.py`.  Vol II residual
greps leave only the intentional deletion-ledger row, guards, and
corrected negative statements.  Cross-volume grep of Vol I and Vol III
found no active positive scalar-to-full-tensor slogan; the active
Vol I hit is already a negative statement saying the scalar shadow
tower is not the full tensor bar object.  Vol II `make fast` converged
after two passes with zero undefined citations and references, and
Vol II `make verify-licensing` reported zero blocking violations and
zero warnings.

A227 tightens the K3 base-point introduction.  The introduction no
longer says that the \(\SCchtop\) construction lands directly on
class-\(\mathcal S\) \(A_1\) on \(\Sigma_{0,24}\).  It now says the K3
base point uses the class-\(\mathcal S\) \(A_1\) Schur-sector candidate
with central charges \((c_{4d},c_{2d})=(107/6,-214)\), Coulomb rank
\(21\), and flavour \(\mathfrak{su}(2)^{24}\).  The \(\SCchtop\)
landing statement is explicitly conditional after Hall--Borcherds
source recognition and the Schur comparison theorem.

The existing finite-gate oracle in
`p1_protected_pfaffian_chain_level.py` remains the operative K3/Class-S
comparison profile.  A new source guard in
`test_p1_protected_pfaffian_chain_level.py` requires the conditional
landing wording and forbids the old ungated "construction lands on"
sentence.

A227 verification: `python3 -m py_compile` passes on
`compute/lib/p1_protected_pfaffian_chain_level.py` and
`compute/tests/test_p1_protected_pfaffian_chain_level.py`; focused
pytest reports `26 passed` for
`test_p1_protected_pfaffian_chain_level.py`.  Vol II residual grep
leaves only the new conditional wording and the guard.  Cross-volume
grep found no active instance of the old \(\SCchtop\)-lands-on-class-S
slogan; other "automatic" hits in Vol III/notes are either explicit
antipattern rows or unrelated K3/CY\(2\) structural facts.  Vol II
`make fast` converged after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reported zero
blocking violations and zero warnings.

A228 corrects the K3 base-point Theorem H summary.  The introduction no
longer says "Hochschild concentration in degree
\(2=\dim_{\mathbb C}(K3)\)."  It now says that the chiral Hochschild
support lies in degrees \(\{0,1,2\}\), with top allowed chiral degree
\(2=\dim_{\mathbb C}(K3)\) after the chiral-specialisation comparison.
This keeps the curve-level chiral Theorem H window separate from the
categorical HKR amplitude of \(D^b\mathrm{Coh}(K3)\), which the
Hochschild supplement records as \(\{0,2,4\}\), not concentration only
in degree \(2\).

A source guard in `test_p1_protected_pfaffian_chain_level.py` now
requires the \(\{0,1,2\}\) chiral window and forbids the old
single-degree wording.

A228 verification: `python3 -m py_compile` passes on
`compute/tests/test_p1_protected_pfaffian_chain_level.py`; focused
pytest reports `27 passed` for
`test_p1_protected_pfaffian_chain_level.py`.  Vol II residual grep
leaves the corrected K3 introduction, the guard, the Hochschild
supplement's explicit negative statement that categorical K3 Hochschild
is not only degree \(2\), and unrelated degree-zero bar-complex uses.
Cross-volume grep of Vol I and Vol III found no active K3
single-degree Theorem H slogan; active Theorem H surfaces use degrees
\(\{0,1,2\}\).  Vol II `make fast` converged after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reported zero blocking violations and zero warnings.

A229 corrects the Phi two-stage IV witness for the K3 rank-24 source.
The witness no longer identifies \(24\) with
\(HH^0(D^b(Coh(K3)))\).  HKR gives
\(HH^0(D^b(Coh(K3)))=H^0(K3,\mathcal O_{K3})=\mathbb C\), hence
dimension \(1\).  The rank-24 number is now recorded as the total even
Mukai/HKR rank \(1+22+1=24\), matching \(\chi(K3)=24\).

`compute/tests/test_phi_two_stage_iv.py` now states this in the
stage-1 docstring, uses the total Mukai/HKR rank in the independent
verification rationale, and adds `categorical_hkr_k3_profile()` with a
guard test requiring `hh0_dim == 1`, even HKR dimensions `(1, 22, 1)`,
and the explicit source label "total even Mukai/HKR rank, not HH^0".

A229 verification: `python3 -m py_compile` passes on
`compute/tests/test_phi_two_stage_iv.py`; focused pytest reports
`5 passed` for `test_phi_two_stage_iv.py`.  Exact residual greps over
Vol II chapters, compute layers, Vol I, and Vol III find no active
"\(HH^0(D^b(Coh(K3))) = \chi(K3) = 24\)" or "formula tied to
Hochschild zeroth cohomology" claim.  No live TeX changed, so no
`make fast` or `make verify-licensing` run was required for A229.

A230 corrects the same Phi two-stage IV witness at the \(S^1\)
specialisation step.  The old rationale said factorisation rank was
preserved by \(S^1\)-integration "up to" the universal
Euler-characteristic factor \(\chi(S^1)=0\), while still returning
rank \(24\).  That is not a typed factorisation-homology statement:
\(\int_{S^1}\) is Hochschild/factorisation homology for an
\(E_2\)-factorisation algebra, not ordinary Euler-characteristic
multiplication.

`compute/tests/test_phi_two_stage_iv.py` now says that
\(\SpCh_{S^1,C}\) applies the typed Hochschild trace from an
\(E_2\) holomorphic factorisation algebra to an \(E_1\) chiral algebra
on \(C\).  The preserved numerical witness is the primitive
Mukai/Heisenberg--Fock current row on the abelian harmonic branch, not
all Hochschild dimensions and not any product by \(\chi(S^1)\).
`stage2_specialisation_scope_profile()` and
`test_s1_specialisation_is_not_euler_characteristic_multiplication()`
guard the source/target types, the non-use of \(\chi(S^1)\), and the
branch-level rank scope.

A230 verification: `python3 -m py_compile` passes on
`compute/tests/test_phi_two_stage_iv.py`; focused pytest reports
`6 passed` for `test_phi_two_stage_iv.py`.  Exact residual greps over
Vol II, Vol I, and Vol III find no active stale phrase "universal
Euler characteristic factor", "rank is preserved by S^1-integration",
or "chi-graded factorisation"; remaining \(\chi(S^1)=0\) hits are
generic sphere/descent examples outside the Phi stage-2 rank witness.
No live TeX changed, so no `make fast` or `make verify-licensing` run
was required for A230.

A231 corrects the \(K3\) topology input in the six-dimensional hCS
anomaly normalisation.  The proof of
Proposition~\(\ref{prop:6dhcs-chi-K3-24}\) no longer says that
\(\chi(K3)=24\) follows from the Hirzebruch signature formula.  The
anomaly cancellation uses the Euler class:
\[
 \frac1{24}\int_{K3}e(T_{K3})
 =\frac1{24}\int_{K3}c_2(T_{K3})
 =1.
\]
The signature formula is now stated only as the separate check
\(\sigma(K3)=-16=\frac13\int_{K3}p_1(T_{K3})\), with
\(p_1=c_1^2-2c_2=-2c_2\).

`compute/lib/six_d_hcs_anomaly.py` now records the Euler source, the
non-source role of the signature formula, \(\int p_1=-48\), and
\(\sigma=-16\).  `test_six_d_hcs_anomaly.py` guards both the oracle and
the live manuscript sentence.

A231 verification: `python3 -m py_compile` passes on
`compute/lib/six_d_hcs_anomaly.py` and
`compute/tests/test_six_d_hcs_anomaly.py`; focused pytest reports
`5 passed` for `test_six_d_hcs_anomaly.py`.  Vol II `make fast`
converges after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reports zero blocking
violations and zero warnings.  Exact greps over Vol II, Vol I, and
Vol III find no active sentence deriving
\(\chi(K3)=24\) from the Hirzebruch signature formula; the remaining
signature-formula paragraph computes \(\int p_1=-48\) and separately
records \(\chi(K3)=24\).

A232 corrects the last source conflation in the same hCS anomaly proof.
The \(1/24\) in the cancellation is not produced by
Hirzebruch--Riemann--Roch.  It is the analytic
heat-kernel/Bernoulli coefficient.  The \(24\) is the topological
Euler integral over \(K3\).  HRR remains only as the compatible K3
index check \(\chi(\mathcal O_{K3})=2\) and \(\chi(T_{K3})=-20\).

The proof now uses the product identity
\[
 \left(\frac{|B_2|}{4}\right)\int_{K3}e(T_{K3})
 =
 \frac1{24}\cdot24
 =
 1.
\]
`compute/lib/six_d_hcs_anomaly.py` records the analytic
`heat_kernel_source` and the limited `hrr_role`; the source guard in
`test_six_d_hcs_anomaly.py` requires the product formula and forbids
the old "follows from the Euler-class term in Riemann--Roch" phrase.

A232 verification: `python3 -m py_compile` passes on
`compute/lib/six_d_hcs_anomaly.py` and
`compute/tests/test_six_d_hcs_anomaly.py`; focused pytest reports
`5 passed` for `test_six_d_hcs_anomaly.py`.  Vol II `make fast`
converges after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reports zero blocking
violations and zero warnings.  Exact greps over Vol II, Vol I, and
Vol III find no active stale source phrase; the only remaining
occurrence is the negative guard in the focused test.

A233 propagates the K3 Euler-source correction into the remaining
compute witnesses.  The P2 gravity-line quartic-effectiveness test no
longer treats \(\chi_{\mathrm{top}}(K3)=24\) as an HRR output.  The
P2 witness records
`chi_top_K3_source = "Euler class / Chern-Gauss-Bonnet"` and an
`hrr_role_for_K3` clause saying HRR checks holomorphic indices such as
\(\chi(\mathcal O_{K3})=2\), not the topological Euler number.

The cross-volume Vol I M24 bridge test now names its path "Milnor 1958
+ Euler topology" and combines \(\chi(K3)=24\) with the CY\(_2\)
shadow-tower relation, not with an HRR source statement.

A233 verification: `python3 -m py_compile` passes on
`compute/lib/p2_gravity_line_pentagon_trace.py` and
`compute/tests/test_p2_gravity_line_pentagon_trace.py`; focused Vol II
pytest reports `27 passed` for
`test_p2_gravity_line_pentagon_trace.py`.  Cross-volume Vol I
`python3 -m py_compile` passes on
`compute/tests/test_cy_m24_bar_bridge_engine.py`; focused Vol I pytest
reports `110 passed` for that file.  Exact greps over Vol II, Vol I,
and Vol III find no active phrase attaching HRR or
Hirzebruch--Riemann--Roch to \(\chi_{\mathrm{top}}(K3)=24\) or
\(\chi(K3)=24\).  No live TeX changed, so no `make fast` or
`make verify-licensing` run was required for A233.

A234 repairs a residual bulk/Hochschild scope leak in the standard
family computations.  The Hochschild chapter no longer says that, for
each standard family, \(C^*(\cA,\cA)\) is simply "computed and
identified with the bulk observables."  It now says that the algebraic
object computed below is the chiral Hochschild \(E_2\)-complex, and
that its interpretation as physical bulk observables requires a chosen
3d HT prefactorization model satisfying the physics bridge and the
boundary-linear exact-sector hypotheses; the identification is the
scoped comparison of Theorem~\(\ref{thm:bulk_hochschild}\).

`compute/tests/test_deletion_ledger.py` now includes
`test_standard_family_bulk_computation_setup_is_scoped()`, which
requires the scoped replacement and forbids the old unqualified
sentence.

A234 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`21 passed` for `test_deletion_ledger.py`.  Fixed-string greps over
Vol II, Vol I, and Vol III find no live copy of the old "computed and
identified with the bulk observables" sentence; the only exact old
phrase is the negative assertion in the new test.  Vol II `make fast`
converges after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reports zero blocking
violations and zero warnings.

A235 propagates the A234 scope correction into live summary and example
entry-points.  The conclusion, introduction, raviolo roadmap, abelian
Chern--Simons computation, and worked affine example no longer state
that bulk observables simply identify with chiral Hochschild cochains.
They now distinguish the algebraic closed-sector Hochschild complex
from physical bulk observables, which appear only after a chosen HT
prefactorisation/realization and the boundary-linear exact-sector
comparison.

`compute/tests/test_deletion_ledger.py` now includes
`test_bulk_hochschild_summary_echoes_are_scoped()`, requiring the
scoped replacement phrases in `conclusion.tex`,
`examples-computing.tex`, `examples-worked.tex`, `introduction.tex`,
and `raviolo.tex`, while forbidding the retired unqualified formulas.

A235 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`22 passed` for `test_deletion_ledger.py`.  Active-surface greps over
Vol II, Vol I, and Vol III chapters/compute find no live unscoped
"bulk observables identify with chiral Hochschild cochains" formulas;
only the negative assertions in the new test guard remain.  Vol II
`make fast` converges after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reports zero
blocking violations and zero warnings.

A236 removes the same bulk/Hochschild overidentification from headings,
tables, and source breadcrumbs.  The Hochschild standard-family
computations are now titled "Chiral Hochschild closed sector for ...";
the local body text speaks of the Hochschild closed-sector state
algebra, closed-sector-to-boundary comparison map, and closed-sector
deformation directions.  The W-algebra copies now use the same visible
heading.  The Hochschild bridge table, the introduction table, the
concordance row, and the `main.tex` input comment no longer say
"Bulk \(\simeq\) chiral Hochschild" or "chiral Hochschild = bulk";
they name the algebraic closed sector and the scoped physical-bulk
comparison instead.

`compute/tests/test_deletion_ledger.py` now extends
`test_bulk_algebra_computation_headings_are_closed_sector()` to guard
the repaired computation band, W-algebra copies, table rows, and source
comment, while preserving reference labels such as `comp:bulk-*`.

A236 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`23 passed` for `test_deletion_ledger.py`.  Active-surface greps over
Vol II, Vol I, and Vol III chapters/compute find no live "Bulk algebra
for", "fixed-level bulk state algebra", "Bulk \(\simeq\) chiral
Hochschild", or "chiral Hochschild = bulk" occurrence outside the
negative guard.  Vol II `make fast` converges after two passes with
zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.

A237 repairs the Lagrangian HKR restatement immediately after the
bulk--Hochschild theorem.  The Hochschild chapter no longer identifies
\(\cO(T^*[-1]\cL_b)\) with the ordinary restriction
\(\cO(\cM_{\mathrm{vac}})|_{\cL_b}\), and no longer says that the bulk
algebra is the ambient shifted-symplectic stack seen from the
Lagrangian.  It now uses the formal Darboux-neighbourhood algebra
\(\cO(\widehat{\cM_{\mathrm{vac}}}_{\cL_b})_{\mathrm{Darboux}}\),
states that ordinary restriction is the boundary chart
\(\cO(\cL_b)\), and keeps the bulk cochain algebra on the shifted
cotangent fibre/formal self-intersection in the boundary-linear exact
sector.

`compute/tests/test_deletion_ledger.py` now includes
`test_lagrangian_hkr_bulk_is_not_ambient_restriction()`, requiring the
Darboux replacement and forbidding the retired ambient-stack slogans.

A237 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`24 passed` for `test_deletion_ledger.py`.  Fixed-string greps over
live Vol II, Vol I, Vol III, and compute source surfaces find the old
ambient-stack slogans only in the negative guard; the audit log keeps
the phrase as the recorded defect.  Vol II `make fast` converges after
two passes with zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A238 repairs the affine Chern--Simons boundary/bulk collapse in the
spectral-braiding and complete-example surfaces.  The abelian
spectral-braiding example no longer says that the bulk algebra is
\(\widehat{\mathfrak u(1)}_k\); it says this is the boundary current
chart, that the physical bulk object is the scoped chiral derived
centre \(C^\bullet_{\mathrm{ch}}(\widehat{\mathfrak u(1)}_k,
\widehat{\mathfrak u(1)}_k)\), and that the Heisenberg Yangian is the
open-colour Koszul dual line algebra.  The complete-example copies now
write \(\mathcal A_{\partial}\), not
\(\mathcal A_{\mathrm{bulk}}\), for the affine Kac--Moody PVA, and
add the scoped derived-centre bulk clause.

`compute/tests/test_deletion_ledger.py` now includes
`test_affine_cs_current_algebra_is_boundary_not_bulk()`, guarding the
live and split files against the retired affine-current-as-bulk
formulas.

A238 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`25 passed` for `test_deletion_ledger.py`.  Exact fixed-string greps
over Vol II, Vol I, Vol III, and compute surfaces find the old formulas
only in the negative guard.  Vol II `make fast` converges after two
passes with zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A239 repairs the \(\mathcal N=4\) celestial symmetry-matching surface.
The frontier chapter no longer states that global symmetry matching
already is the equivalence
\(\cA_{\mathrm{bulk}}\simeq Z_{\mathrm{der}}(\cA_\partial)\).  It now
names the comparison map
\(\gamma_{\cN=4}\colon\cA_{\mathrm{bulk}}\to
Z_{\mathrm{der}}(\cA_\partial)\) and says that this map becomes a
quasi-isomorphism precisely when the mixed obstruction class
\(\Obs_{\cN=4}^{\mathrm{mix}}\in H^2(\Gamma_{\cN=4})\) vanishes.  The
bar/BV graph coalgebra matching is likewise only perturbative until
that obstruction gate is closed.

`compute/tests/test_deletion_ledger.py` now includes
`test_n4_symmetry_matching_is_obstruction_gated()`, requiring the
comparison-map formulation in the active frontier and split celestial
files and forbidding the retired unconditional derived-centre
equivalence.

A239 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`26 passed` for `test_deletion_ledger.py`.  Exact fixed-string greps
over Vol II, Vol I, Vol III, and compute surfaces find the old
\(\mathcal N=4\) unconditional derived-centre slogans only in the
negative guard.  Vol II `make fast` converges after two passes with
zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A240 repairs the boundary-linear bulk/line exact-sector wording.  The
boundary-linear theorem no longer says that the derived centre "equals"
the bulk algebra after displaying quasi-isomorphisms.  It now says
\(\cO(\dCrit(W))\) is canonically quasi-isomorphic to
\(\Zder(B_{L,W})\), i.e. equal in the derived Morita category and not
as a literal dg presentation.  The local bulk/line remark no longer
says that bulk is the derived centre of the line algebra; it says a
pointed line algebra \(K_\kappa\) is a Morita presentation of the local
line category, and the completed local bulk algebra is computed by
\(\HH^\bullet(K_\kappa)\), equivalently by the derived centre of that
local line category.

`compute/tests/test_deletion_ledger.py` now includes
`test_boundary_linear_bulk_line_is_morita_scoped()`, guarding the active
core and split file against the retired equality and line-algebra-centre
slogans.

A240 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`27 passed` for `test_deletion_ledger.py`.  Exact fixed-string greps
over live Vol II, Vol I, Vol III, and compute source surfaces find the
old phrases only in the negative guard; the audit log keeps one phrase
as the recorded defect.  Vol II `make fast` converges after two passes
with zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A241 repairs the worked-example derived-centre bulk language.  The
universal defect paragraph no longer says that
\(C^\bullet_{\mathrm{ch}}(\cA_{\mathrm{open}},\cA_{\mathrm{open}})\)
computes the universal bulk algebra.  It says this complex is the
closed-sector derived-centre vertex, and that identification with the
universal physical bulk algebra uses the Kodaira--Spencer HT
realisation and the scoped bulk--Hochschild comparison, not Verdier
duality or bar-cobar inversion.  The 5d HT structure theorem now says
that, in that scoped 5d HT realisation, bulk observables are represented
by the chiral derived centre.  The affine benchmark now says that, on
the bulk--Hochschild comparison surface, the benchmark bulk object is
represented by the chiral derived centre.

`compute/tests/test_deletion_ledger.py` now includes
`test_examples_worked_derived_centre_bulk_is_scoped()`, requiring the
scoped formulations and forbidding the retired automatic universal-bulk
phrases.

A241 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`28 passed` for `test_deletion_ledger.py`.  Exact fixed-string greps
over active Vol II, Vol I, Vol III, and compute source surfaces find no
surviving stale occurrence; the raw full-surface grep sees only the
negative guard and the audit record of the retired phrases.  Vol II
`make fast` converges after two passes with zero undefined citations
and references, and Vol II `make verify-licensing` reports zero
blocking violations and zero warnings.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched files.

A242 repairs the gravity-sector derived-centre bulk reading.  The
gravitational \(\Sigma_n\)-descent remark no longer says that
\(Z^{\mathrm{der}}_{\mathrm{ch}}(\mathrm{Vir}_c)\) is the universal
bulk algebra of the gravitational theory.  It says this object is the
algebraic closed-colour vertex attached to the Virasoro boundary chart,
computed by chiral Hochschild cochains and not by the bar complex.  Its
physical gravitational reading now requires the Virasoro HT
prefactorization realisation, the Brown--Henneaux boundary chart
\(c=3\ell/(2G_N)\), and the \(\Sigma_n\)-closed-sector descent; pure
3d gravity further passes to the central-extension quotient
\(\C[\![c]\!]\).

`compute/tests/test_deletion_ledger.py` now includes
`test_gravity_derived_centre_bulk_is_brown_henneaux_scoped()`, requiring
the scoped gravity formulation and forbidding the retired
derived-centre-as-whole-physical-bulk phrases.

A242 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`29 passed` for `test_deletion_ledger.py`.  Exact fixed-string greps
over Vol II, Vol I, Vol III, and compute source surfaces find no
surviving stale occurrence when the negative guard is excluded; the raw
non-audit grep sees only the negative guard.  Vol II `make fast`
converges after two passes with zero undefined citations and
references, and Vol II `make verify-licensing` reports zero blocking
violations and zero warnings.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched files.

A243 repairs the concordance and Hochschild summary versions of the
derived-centre bulk slogan, then propagates the same correction across
the Vol I and Vol III compute surfaces that still encoded
"derived centre = universal bulk" as oracle text.  The active
concordance architecture no longer says that
\(\cZ^{\mathrm{der}}_{\mathrm{ch}}(\cA)\) serves as the universal bulk
algebra.  It says this object is the algebraic closed-sector vertex of
the boundary chart and becomes a physical bulk observable algebra only
after a chosen HT prefactorization realisation satisfies the
bulk--Hochschild comparison hypotheses.  The active Hochschild package
bullet now says the boundary-linear exact-sector comparison identifies
physical bulk with the Morita-invariant derived centre of the boundary
category.  The non-input THQG extension, trimmed preface, draft
foundation echo, Vol III chiral/CY bridge strings, and Vol I
derived-centre/CSFT/K3E/Costello-BV compute oracles now carry the same
stage distinction.

`compute/tests/test_deletion_ledger.py` now includes
`test_concordance_derived_centre_bulk_summaries_are_scoped()`, requiring
the scoped formulations and forbidding the retired universal-bulk and
bulk-is-derived-centre slogans on the repaired surfaces.

A243 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`30 passed` for `test_deletion_ledger.py`.  Cross-volume Python compile
checks pass for the touched Vol I and Vol III compute modules.  Focused
Vol I pytest reports `932 passed` across the touched derived-centre,
CSFT, K3E, Costello-BV, and theorem guards.  Focused Vol III pytest
reports `94 passed` for `test_coha_drinfeld_bulk.py`.  Broad exact grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving non-audit occurrence of the retired universal-bulk and
derived-centre-equals-bulk forms.  Vol II `make fast` converges after
two passes with zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched Vol II, Vol I, and Vol III files.

A244 repairs active reader-facing bulk shorthand after the
universal-bulk slogans were removed.  The preface no longer writes
\(\Abulk=\cH_k\) or says that the bulk is the abelian chiral algebra.
It says \(\Zderch(\cH_k)=\bulkChirHoch{\cH_k}\) is the algebraic
closed-sector vertex, whose physical bulk reading is the rank-one exact
HT comparison and not literal equality with the boundary chart.  The
introduction now says the level-\(\mathsf Z\) object is the algebraic
derived chiral centre and physical bulk requires the scoped HT
open/closed comparison.  Its Virasoro gravity paragraph now says the
chiral Hochschild object is the algebraic closed-sector vertex and the
physical 3d-gravity reading passes through Brown--Henneaux.  The
super-Yangian universal-holography sentence now says the super derived
chiral centre is the algebraic closed sector, with physical bulk in the
same HT comparison package.

`compute/tests/test_deletion_ledger.py` now includes
`test_active_bulk_shorthand_is_closed_sector_scoped()`, requiring these
scoped active formulations and forbidding the four retired shorthand
phrases.

A244 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`31 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the four retired active phrases outside the
negative guard.  Vol II `make fast` converges after two passes with
zero undefined citations and references, and Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A245 repairs the physical-origin opening.  The active HT physical
origins chapter no longer says that every chiral algebra in the
\(E_1\) core arises as a Costello--Li holomorphic-topological boundary
theory.  It says the physical-origin lane contains those chiral
algebras for which a Costello--Li HT field theory, or an explicit
open/closed comparison datum, has been constructed; the standard
families enter through named realisations, and an abstract
\(E_1\)-core chiral algebra does not thereby acquire a physical HT
origin.  The same paragraph now calls
\(\mathcal Z^{\mathrm{der}}_{\mathrm{ch}}(\cA)\) the algebraic
closed-sector vertex, with physical bulk reading requiring the chosen
HT open/closed comparison.

`compute/tests/test_deletion_ledger.py` now includes
`test_ht_physical_origin_requires_constructed_realisation()`, requiring
the constructed-realisation gate and forbidding the retired universal
HT-origin and bulk-algebra phrases.

A245 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`32 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired phrases outside the negative guard.
Vol II `make fast` converges after two passes with zero undefined
citations and references, and Vol II `make verify-licensing` reports
zero blocking violations and zero warnings.  Scoped `git diff --check`
and trailing-whitespace scans pass on the touched files.

A246 repairs bare HT path-integral derivations in the physical-origins
chapter.  The localization theorem is now conditional perturbative
localization: it fixes a Costello--Li holomorphic twist, a BV gauge
fixing, and a renormalisation scheme, identifies the perturbative BV
theory near the \(Q\)-fixed locus with a derived holomorphic moduli
problem, and says the ordinary moduli space is the underived
truncation.  The boundary OPE is now the renormalised boundary
factorization product of the gauge-fixed bulk theory, and the HCS
Kac--Moody example is a computation from gauge-fixed HCS BV Feynman
rules.  The BV/bar-cobar pairing is now the algebraic
configuration-space model for the perturbative BV integral, not a bare
measure-theoretic path-integral identity; the same correction is
propagated to `bv_ht_physics.tex`.

`compute/tests/test_deletion_ledger.py` now includes
`test_ht_path_integral_claims_are_perturbative_bv_scoped()`, requiring
the gauge-fixing/renormalisation language and forbidding the retired
bare path-integral derivation phrases.

A246 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`33 passed` for `test_deletion_ledger.py`.  Broad exact fixed-string
grep over Vol II, Vol I, and Vol III manuscript/compute surfaces finds
no surviving occurrence of the retired phrases outside the negative
guard.  Vol II `make fast` converges after two passes with zero
undefined citations and references, and Vol II `make verify-licensing`
reports zero blocking violations and zero warnings.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
files.

A247 repairs the remaining AGT and class-S localization overclaims in
the HT physical-origins lane.  The AGT dictionary now says equivariant
localization computes the Nekrasov partition function by fixed-point
data on the instanton/Hitchin moduli problem, and that under the AGT
comparison hypotheses and the Nekrasov--Shatashvili limit this
equivariant partition function matches the \(\mathcal W(G)\) conformal
block.  It explicitly denies a literal unregularised 4d path integral
over \(\mathcal M_{\mathrm{Hit}}\).  The class-S/N=4 theorem now says
the Costello holomorphic twist, BV gauge fixing, and dimensional
reduction produce perturbative holomorphic BF theory whose derived
classical moduli problem is the \(\bar\partial\)-connection stack.  Its
localization step is perturbative BV localization in the formal
neighbourhood of the \(Q\)-fixed locus.  The HCS OPE is now computed
from renormalised HCS/BF BV Feynman rules.  The compact
`bv_ht_physics.tex` echo carries the same corrections.

`compute/tests/test_deletion_ledger.py` now includes
`test_ht_agt_and_class_s_localization_are_scoped()`, requiring the
equivariant localization / AGT / perturbative BV language and
forbidding the retired localization/path-integral slogans.

A247 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`34 passed` for `test_deletion_ledger.py`.  Broad exact fixed-string
grep over Vol II, Vol I, and Vol III manuscript/compute surfaces finds
no surviving occurrence of the retired phrases outside the negative
guard.  Vol II `make fast` converges after two passes with zero
undefined citations and references.  `make verify-licensing` reports
zero blocking violations and zero warnings.  Scoped `git diff --check`
and trailing-whitespace scans pass on the touched files.

A248 repairs an active affine-example promotion of algebraic centre
data to a full Chern--Simons path-integral statement.  The generic
level equality \(Z(\widehat{\fg}_k)=\C\) now states exactly that the
chiral Hochschild closed-sector complex has no central local class
beyond the unit.  Its physical reading is conditional on a chosen
perturbative Chern--Simons/BV boundary comparison, where the unit is the
vacuum summand of local closed observables.  The paragraph now denies
that this proves a theorem about the full non-perturbative
Chern--Simons functional integral, whose global flat-connection sectors,
framing dependence, and finite-level modular data require separate
\(3d\) input.

`compute/tests/test_deletion_ledger.py` now includes
`test_affine_generic_center_not_promoted_to_cs_path_integral()`,
requiring the local Hochschild / perturbative BV comparison language
and forbidding the retired path-integral promotion.

A248 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`35 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired CS path-integral promotion outside
the negative guard.  Vol II `make fast` converges after two passes with
zero undefined citations and references.  `make verify-licensing`
reports zero blocking violations and zero warnings.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
files.

A249 repairs the hCS finiteness bridge in `bv_brst.tex`.  The bridge no
longer promotes hCS finiteness to a literal six-dimensional
path-integral theorem on \(\bR\times K3\times E\).  It now says the
renormalised perturbative 6d hCS/BV amplitudes are finite at every
genus in the Fulton--MacPherson compactified configuration-space model.
The paragraph explicitly states that this is about gauge-fixed
counterterms and factorisation-algebra amplitudes after installing a
regulator, not an ordinary non-perturbative analytic hCS measure
integral.

`compute/tests/test_deletion_ledger.py` now includes
`test_bv_brst_hcs_finiteness_is_perturbative_regulator_statement()`,
requiring the perturbative BV / regulator language and forbidding the
retired hCS path-integral formulation.

A249 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`36 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired hCS path-integral formulation
outside the negative guard.  Vol II `make fast` converges after two
passes with zero undefined citations and references.  `make
verify-licensing` reports zero blocking violations and zero warnings.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A250 repairs the \(3d\) gravity comparison with Chern--Simons torsion.
The text no longer promotes Reidemeister torsion with a bare
\(|\eta(\tau)|^{-2}\) answer to the full Chern--Simons functional
integral.  It now states the semiclassical one-loop datum: after
fixing a flat connection and boundary polarization, the quadratic BV
determinant is the Ray--Singer/Reidemeister torsion.  For the thermal
AdS\(_3\) vacuum this torsion is the eta/vacuum-character determinant
factor after the chosen boundary-mode normalization, not the full path
integral over all flat connections or the complete modular sum.

`compute/tests/test_deletion_ledger.py` now includes
`test_gravity_cs_torsion_is_semiclassical_not_full_path_integral()`,
requiring the flat-connection / boundary-polarization wording and
forbidding the retired torsion-as-full-path-integral slogan.

A250 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`37 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired torsion-as-full-path-integral
slogan outside the negative guard.  Vol II `make fast` converges after
two passes with zero undefined citations and references.  `make
verify-licensing` reports zero blocking violations and zero warnings.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A251 repairs the Gopakumar--Vafa scalar-lane comparison in
`modular_pva_quantization_core.tex`.  The genus-one scalar term no
longer matches a CS path-integral determinant as such; it matches the
gauge-fixed semiclassical one-loop determinant of the Chern--Simons BV
complex at the trivial flat connection, after the same boundary-mode
normalization.  The text explicitly denies that this is the full
matrix-model functional integral.

`compute/tests/test_deletion_ledger.py` now includes
`test_modular_pva_gv_one_loop_is_gauge_fixed_seed()`, requiring the
gauge-fixed BV seed language and forbidding the retired CS
path-integral formulation.

A251 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`38 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired GV/CS path-integral determinant
phrase outside the negative guard.  Vol II `make fast` converges after
two passes with zero undefined citations and references.  `make
verify-licensing` reports zero blocking violations and zero warnings.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A252 repairs the Heisenberg recursion proof in
`thqg_modular_bootstrap.tex`.  The genus-\(g\) object is now the
Gaussian determinant-line amplitude of the gauge-fixed free-boson
factorization algebra, equal to the \(k\)-th determinant-line power.
The proof explicitly states that this is the perturbative determinant
model for the Heisenberg sector, not an unqualified
measure-theoretic path integral.  The linearity in \(k\), the connected
logarithm argument, and the conclusion \(k\omega_g\) are preserved.

`compute/tests/test_deletion_ledger.py` now includes
`test_heisenberg_recursion_uses_gaussian_determinant_not_path_integral()`,
requiring the determinant-line wording and forbidding the retired
genus-\(g\) Heisenberg path-integral phrase.

A252 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`39 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired genus-\(g\) Heisenberg
path-integral phrase outside the negative guard.  Vol II `make fast`
converges after two passes with zero undefined citations and
references.  `make verify-licensing` reports zero blocking violations
and zero warnings.  Scoped `git diff --check` and trailing-whitespace
scans pass on the touched files.

A253 repairs the BV integration sections in `ht_physical_origins.tex`
and the compact split echo `bv_ht_physics.tex`.  The displayed BV
expression is now a gauge-fixed perturbative BV functional formally
denoted by the integral sign.  The symbol \([D\phi]\) is identified as
the induced formal BV Berezinian density after the renormalisation
prescription has been fixed; it is not an ordinary measure on a space
of fields.  The configuration-space object is now a density, and the
comparison table says "Perturbative BV functional" rather than "Path
integral."

`compute/tests/test_deletion_ledger.py` now strengthens
`test_ht_path_integral_claims_are_perturbative_bv_scoped()`, requiring
the formal Berezinian / perturbative-functional language and forbidding
the retired measure/path-integral formulations.

A253 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`39 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired BV measure/path-integral
formulations outside the negative guard.  Vol II `make fast` converges
after two passes with zero undefined citations and references.  `make
verify-licensing` reports zero blocking violations and zero warnings.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A254 repairs the perturbative-finiteness scalar-shadow dictionary in
`thqg_perturbative_finiteness.tex`.  The shadow free energy
\(F_g(\cA)\) is now the genus-\(g\) scalar vacuum-energy shadow: the
connected no-insertion coefficient of the algebraic shadow determinant
expansion.  Its cosmological-constant reading requires a chosen bulk
comparison datum and is not a contribution to an independently
constructed path integral.  The pure-gravity formula is now a
scalar-shadow partition series, not the full gravitational partition
function over bulk geometries.  The holographic dictionary row now says
"genus-\(g\) scalar vacuum amplitude"; the comparison with known
results now cites the metric one-loop determinant rather than a gravity
path-integral equality.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_finiteness_shadow_free_energy_not_path_integral()`, requiring
the scalar-shadow wording and forbidding the retired path-integral
formulations.

A254 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`40 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired scalar-shadow/path-integral
formulations outside the negative guard.  Vol II `make fast` converges
after two passes with zero undefined citations and references.  `make
verify-licensing` reports zero blocking violations and zero warnings.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A255 repairs the theorem-level partition-function scope in
`thqg_perturbative_finiteness.tex` and the split Fredholm echo
`thqg_fredholm_partition_functions.tex`.  The local definition is now
the scalar shadow genus series
`Z_{\mathrm{sh}}^{\mathrm{scal}}`; the gravitational notation
`Z_{\mathrm{grav}}^{\mathrm{scal}}` is introduced only after a bulk
comparison datum has been fixed.  The `\ClaimStatusProvedHere` theorem
now proves convergence and finite-depth properties for the completed
degree-summed shadow partition series
`Z^{\mathrm{sh,deg}}`, not for a full gravitational partition
function over bulk metrics.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_degree_summed_series_not_full_gravity_partition_function()`,
requiring scalar-shadow and degree-summed-shadow wording in the live
chapter and the Fredholm split file, and forbidding the retired
full-gravitational-partition-function forms.

A255 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`41 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving live occurrence of the retired theorem-level forms outside
the negative guard; remaining hits are Vol I archived audit/orphan
markdown snapshots, not live manuscript or compute inputs.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` reports zero blocking
violations and zero warnings.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched files.

A256 repairs the one-loop motivic-\(\zeta(2)\) OPE proof in
`chiral_ce_factalg_gen_rel.tex`.  The coefficient is now identified as
the renormalised BV Feynman-rule residue after fixing the Costello
renormalisation scheme and Bergman propagator.  The displayed
\(\P^1\times\P^1\) integral is explicitly a configuration-space
integral, not an unqualified measure-theoretic path integral over
fields.

`compute/tests/test_deletion_ledger.py` now includes
`test_chiral_ce_one_loop_uses_renormalised_feynman_residue()`, requiring
the renormalised-Feynman-residue wording and forbidding the retired
path-integral phrase.

A256 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`42 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired phrase outside the negative guard.
Vol II `make fast` converges after two passes with zero undefined
citations and references.  `make verify-licensing` reports zero
blocking violations and zero warnings.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched files.

A259 repairs the BV integration naming surface in
`ht_physical_origins.tex` and the split echo `bv_ht_physics.tex`.
The displayed BV expression is now consistently introduced as a
gauge-fixed perturbative BV functional using an induced formal BV
Berezinian density after a renormalisation prescription has been fixed.
The section headings no longer call this a BV path integral, and the
definition titles no longer call it a BV partition function.

The final comparison remarks now state that the bar-cobar pairing
requires a configuration-space density model for BV gauge-fixing data.
The gauge-fixing Lagrangian fixes the Fulton--MacPherson
regularisation, and the resulting density pairs the chiral bar and
cobar objects by a compactified configuration-space integral.  This is
the construction present in the manuscript; it is not a
measure-theoretic identification of a BV path integral with a
bar-cobar pairing.

`compute/tests/test_deletion_ledger.py` now strengthens
`test_ht_path_integral_claims_are_perturbative_bv_scoped()` so it
requires the perturbative-functional heading, the gauge-fixed
definition title, the formal-Berezinian-density step, and the
configuration-space density comparison.  The same guard forbids the
retired BV-path-integral, BV-partition-function, and BV-measure labels.

A259 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`44 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired BV-path-integral identification or
partition-function comparison outside the negative guard.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A260 repairs the BV--HT extension supplement
`thqg_bv_ht_extensions.tex`.  The section that had developed a
conjectural equality between a BV path integral and the bar-cobar
pairing is now explicitly a construction of perturbative BV graph
functionals.  The constructed object is a graph-weight BV density on
the genus-\(g\) chiral bar complex, and its integrated output is the
scalar graph shadow \(Z_g^{\mathrm{BV,gr}}\).

The Heisenberg computation is kept as a Gaussian determinant shadow:
the determinant formula remains, but the prose no longer constructs an
ordinary field-space measure on maps \(\Sigma_g\to\mathbb R^N\).  The
Mumford and anomaly-cancellation paragraphs now speak of scalar graph
shadows rather than BV partition functions.  The Virasoro gravity
paragraph now says the graph formula is a scalar comparison surface for
the Maloney--Witten problem, not a construction of the gravitational
path integral over hyperbolic \(3\)-manifolds; the tower terms are
perturbative shadow corrections in the bar complex.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_bv_extensions_use_graph_shadow_not_bv_partition_function()`,
requiring the graph-functional, graph-density, scalar-shadow, and
Virasoro-scalar-shadow wording and forbidding the retired
BV-measure/BV-partition-function/path-integral phrases.

A260 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`45 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired BV-path-integral equality,
BV-partition-function definition, Maloney--Witten "precise modular
sum" phrase, gravitational-instanton phrase, or field-space
path-integral-over-maps phrase outside the negative guard.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A257 repairs the free-boson one-loop amplitude in
`feynman_connection.tex`, with propagation to the parallel Vol I file
`../chiral-bar-cobar/chapters/connections/feynman_connection.tex`.
The elliptic formula
\((\operatorname{Im}\tau)^{-1/2}|\eta(\tau)|^{-2}\) is now the
zero-mode-normalised, zeta-regularised Gaussian determinant of the free
boson on \(E_\tau\), i.e. the analytic determinant-line amplitude of
the free theory.  The bar complex still supplies the zero-mode-reduced
vacuum line and scalar anomaly; it does not construct the analytic
function as a field-space measure integral.

`compute/tests/test_deletion_ledger.py` now includes
`test_feynman_connection_one_loop_uses_zeta_regularised_determinant()`,
requiring the determinant-line wording and forbidding the retired
\(E_\tau\) Gaussian-path-integral phrase.

A257 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`43 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired \(E_\tau\) phrase.  Broader
"Gaussian path integral" hits remain in Vol I examples/compute/
standalone surfaces and require separate local audits rather than this
specific elliptic determinant patch.  Vol II `make fast` converges
after two passes with zero undefined citations and references.  Vol II
`make verify-licensing` reports zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched Vol II files and the propagated Vol I file.

A258 repairs the Costello--Francis--Gwilliam comparison in
`spectral-braiding-core.tex`.  The filtered structure of the CFG
\(E_3\)-algebra now records the loop expansion of the renormalised
perturbative BV Feynman rules; the remark explicitly does not assert a
non-perturbative BV path integral.  The factorisation-homology trace
comparison and the annular-bar shadow truncation comparison are
unchanged.

`compute/tests/test_deletion_ledger.py` now includes
`test_spectral_braiding_cfg_uses_perturbative_bv_feynman_expansion()`,
requiring the perturbative-Feynman wording and forbidding the retired
BV-path-integral formulation.

A258 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`44 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving live occurrence of the retired CFG phrase outside the
negative guard; remaining matches are Vol I archived repair/rescue
markdown notes, not live manuscript or compute inputs.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` reports zero blocking
violations and zero warnings.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched files.

A261 repairs the active Heisenberg all-genera BV/bar proof in
`bv_brst.tex`, with propagation to the Vol I echo
`../chiral-bar-cobar/chapters/connections/bv_brst.tex` and the Vol I
compute proof `../chiral-bar-cobar/compute/lib/heisenberg_bv_bar_proof.py`.
The scalar equality
\[
F_g^{\mathrm{BV}}(\cH_\kappa)=F_g^{\mathrm{bar}}(\cH_\kappa)
=\kappa\,\lambda_g^{\mathrm{FP}}
\]
is unchanged.  What changed is the analytic object feeding the proof:
the free-boson side is now the zero-mode-reduced gauge-fixed
determinant-line scalar
\[
Z_g^{\mathrm{det}}(\cH_\kappa;\Sigma_g)
=(\det'_\zeta \bar\partial_{\Sigma_g})^{-\kappa},
\]
not an asserted BV partition function or Gaussian field-space
integral.  The Selberg factorisation and the Quillen/GRR argument now
refer to determinant-line scalars, and the higher-genus BRST remark
requires moduli-space density, gauge-fixing data, and Costello
counterterms rather than a path-integral measure.

`compute/tests/test_deletion_ledger.py` now includes
`test_bv_brst_heisenberg_uses_determinant_line_scalar_not_bv_partition()`,
checking the Vol II manuscript surface, the Vol I manuscript echo, and
the Vol I compute proof commentary.  The guard forbids the retired
BV-partition-function, Gaussian-functional-integral,
path-integral-measure, and WZW-path-integral phrases.

A261 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py` and on
`../chiral-bar-cobar/compute/lib/heisenberg_bv_bar_proof.py`; focused
pytest reports `46 passed` for `test_deletion_ledger.py`.  Exact
fixed-string grep over Vol II and Vol I manuscript/compute surfaces
finds no surviving live occurrence of the retired Heisenberg
BV-partition-function sentence, Gaussian-functional-integral phrase,
moduli path-integral-measure phrase, WZW path-integral phrase, or old
\(Z_g^{\mathrm{BV}}\) determinant symbol.  Remaining hits are archived
audit/resume markdown snapshots or negative guard strings.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A262 repairs the active K3 Swiss-cheese remark in
`factorization_swiss_cheese.tex`.  The paragraph no longer calls the
object a Swiss-cheese BV partition function at \(K3\times T^2\).  It
now states the constructed object: a perturbative Swiss-cheese anomaly
line with scalar Deligne characteristic in
\(H^2(\overline{\mathcal A_2},\mathcal O^*)\).  The \(\Delta_5\)
statement is gated by Pfaffian orientation, Hall--Borcherds comparison,
and K3 elliptic-genus normalisation.  The Borcherds weight formula
\(\kappa_{\mathrm{BKM}}=\mathrm{wt}(\Delta_5)
=c_{\phi_{0,1}^{K3}}(0)/2=5\) is unchanged.

`compute/tests/test_deletion_ledger.py` now includes
`test_factorization_swiss_cheese_k3_delta5_is_scalar_anomaly_characteristic()`,
requiring the scalar-Deligne-characteristic wording and forbidding the
retired BV-partition-function / unique-Deligne-trivialiser phrasing.

A262 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`47 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired K3 Swiss-cheese BV-partition-function
sentence.  The only remaining "unique trivialiser in Deligne cohomology"
hit is the negative guard.  Vol II `make fast` converges after two
passes with zero undefined citations and references.  `make
verify-licensing` exits 0 with zero blocking violations; its printed
missing-tag inventory is the pre-existing nonblocking audit surface.

A263 repairs the ordered associative Drinfeld programme meta-remark in
`ordered_associative_chiral_kd_frontier.tex`.  The boundary comparison
subproblem no longer says that the hemisphere partition function equals
the cyclic pairing on the ordered bar complex.  It now states the
precise comparison target: after boundary polarization, localization
density, and the hemisphere comparison map are fixed, the Dimofte
hemisphere scalar amplitude maps to the cyclic form on
\(B^{\mathrm{ord}}(\cA)\).  The paragraph explicitly says this is not
an identification of a physical partition function with the raw ordered
bar complex.  The universal reconstructor conclusion is now conditional:
it begins only when all four data are supplied.

`compute/tests/test_deletion_ledger.py` now includes
`test_ordered_drinfeld_boundary_pairing_is_gated_scalar_shadow()`,
requiring the gated scalar boundary-pairing language and forbidding the
retired hemisphere-partition-function equality and unconditional
programme conclusion.

A263 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`48 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired hemisphere-partition-function
equals ordered-bar cyclic pairing sentence or the retired
Dimofte-hemisphere-equals-algebraic-pairing sentence.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A264 repairs the active Heisenberg base case in
`modular_pva_quantization_frontier.tex`.  The hemisphere comparison no
longer identifies a \(D^3\) abelian-Chern--Simons value with a
perturbative path integral around a Gaussian saddle, nor calls the
displayed value a direct perturbative Gaussian integral.  The section
now names the actual object: a localised hemisphere scalar amplitude,
defined only after boundary polarisation, Omega-background
localisation density, gluing density, and one-loop determinant
regularisation are fixed.  The displayed abelian formula is now
\[
Z_{D^3}^{\mathrm{scal}}\bigl[\mathrm{U}(1)_k\bigr]
=\mathfrak c\,k^{-1/2},
\]
and is stated as the formal Gaussian determinant-line prediction.  The
DFPY16 citation remains only as the \(S^3\) sphere-localisation
analogue; the manuscript still says that a primary \(D^3\) hemisphere
source would be needed to fix the analytic normalisation.

`compute/tests/test_deletion_ledger.py` now includes
`test_modular_pva_heisenberg_hemisphere_uses_scalar_amplitude()`,
requiring the scalar-amplitude heading, label, indices,
\(Z_{D^3}^{\mathrm{scal}}\) symbol, and determinant-line wording, and
forbidding the retired path-integral, direct-Gaussian-integral, old
index, old label, and old \(Z_{D^3}\) strings.

A264 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`49 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving occurrence of the retired Gaussian-saddle path-integral
phrase, the retired direct-perturbative-Gaussian-integral phrase, or
the retired \(D^3\) abelian-CS hemisphere-partition-function sentence
outside the negative guard.  The only remaining retired index string is
the negative guard.  Vol II `make fast` converges after two passes with
zero undefined citations and references.  `make verify-licensing`
exits 0 with zero blocking violations; its printed missing-tag
inventory is the pre-existing nonblocking audit surface.

A265 repairs the rank-\(d\) Heisenberg corollary in
`thqg_modular_bootstrap.tex`.  The corollary no longer calls
\[
(\det\operatorname{Im}\Omega_g)^{-dk/2}
\]
the genus-\(g\) partition function, and no longer says it is the
exponential of the Faber--Pandharipande free energy.  The manuscript
now separates the two scalar objects:
\[
F_g^{(d)}=dk\,\lambda_g^{\mathrm{FP}}
\]
is the integrated connected Hodge trace, while
\[
\chi_{g,\mathrm{det}}^{(d)}(\Omega_g)
=
(\det\operatorname{Im}\Omega_g)^{-dk/2}
\]
is the zero-mode Gaussian determinant-line scalar representative over
the period matrix locus.  The proof now says the Fredholm computation
gives the determinant-line representative, whose logarithmic curvature
is \(dk\) times the Hodge form; applying the scalar trace gives
\(F_g^{(d)}\).  The result is explicitly not a physical partition
function.

`compute/tests/test_deletion_ledger.py` now strengthens
`test_heisenberg_recursion_uses_gaussian_determinant_not_path_integral()`
to require the determinant-line scalar wording,
\(\chi_{g,\mathrm{det}}^{(d)}\), and the Hodge determinant-line
representative statement, while forbidding the retired
partition-function heading, old rank-\(d\) \(Z\)-symbol, old label
comment, and free-energy exponentiation sentence.

A265 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`49 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving live occurrence of the retired exponentiation sentence,
genus-\(g\) partition-function heading, old rank-\(d\) \(Z\)-symbol,
or old rank-\(d\) label comment outside the negative guard.  Remaining
external hits are archived Vol I audit/orphan markdown snapshots, not
live manuscript or compute inputs.  Vol II `make fast` converges after
two passes with zero undefined citations and references.  `make
verify-licensing` exits 0 with zero blocking violations; its printed
missing-tag inventory is the pre-existing nonblocking audit surface.

A266 repairs the effective scalar Virasoro object in
`thqg_3d_gravity_movements_vi_x.tex`.  The \(u\)-eigenvalue
decomposition no longer names \(Z_g^{\mathrm{scal}}\) as a partition
function.  It is now the genus-\(g\) effective scalar shadow series,
and the all-genera object
\[
Z_{\mathrm{grav}}^{\mathrm{scal}}
=\exp(\mathcal F_{\mathrm{grav}}^{\mathrm{scal}})
\]
is stated as a scalar generating series.  The \(c=13\) specialization
no longer says self-duality means inversion of the scalar series.
Self-duality is now the reflection invariance
\(\mathcal S\colon c\mapsto 26-c\), with
\(\mathcal S(Z_{13}^{\mathrm{scal}})=Z_{13}^{\mathrm{scal}}\).

`compute/tests/test_deletion_ledger.py` now includes
`test_gravity_effective_scalar_shadow_series_not_partition_function()`,
requiring the scalar-shadow-series title, index, decomposition,
scalar-generating-series proof sentence, normalized-\(\hat A\)
exponential-series wording, and reflection-invariance statement, while
forbidding the retired partition-object title/index/body forms, old
exponentiation opener, and self-duality-as-inversion sentence.

A266 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`50 passed` for `test_deletion_ledger.py`.  Live TeX grep over
`main.tex`, `chapters/`, and `appendices/` finds no surviving exact
occurrence of the retired effective-scalar partition-object phrase, the
retired self-duality-as-inversion sentence, or the retired
\(\hat A\)-exponential wording.  Broader repository hits are old audit
notes or negative guard strings, not active manuscript input.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A267 repairs the modular-bootstrap output statement in
`thqg_modular_bootstrap.tex`.  The MC recursion no longer claims to
determine a complete gravitational partition function from boundary
data alone.  It is now algebraically complete for the fixed
bar--modular shadow: it determines connected scalar shadow classes.
The one-loop theorem now writes \(Z_1^{\mathrm{scal}}\).  The Gaussian
non-renormalization theorem now uses the connected scalar amplitude
\[
\mathcal A_g^{\mathrm{scal,conn}}(\cA)
=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\]
and the disconnected all-genus object is
\[
Z_{\mathrm{sh}}^{\mathrm{scal}}(\cA;\hbar)
=\exp\Bigl(\sum_{g\ge 1}\hbar^g
\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}\Bigr).
\]
The gravitational reading is conditional on a separate comparison
between \(\Theta_\cA\), a bulk saddle expansion, and the completed
physical partition function.

`compute/tests/test_deletion_ledger.py` now includes
`test_modular_bootstrap_outputs_scalar_shadow_not_gravity_partition()`,
requiring the fixed-shadow completeness statement, scalar-shadow
homotopy consequence, one-loop scalar-shadow theorem, connected scalar
amplitude, disconnected scalar-shadow generating series, and
separate-comparison sentence, while forbidding the retired genus
partition-function, one-loop partition-function, genus-power Gaussian,
and complete-gravity-output forms.

A267 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`51 passed` for `test_deletion_ledger.py`.  Exact fixed-string grep
over Vol II, Vol I, and Vol III manuscript/compute surfaces finds no
surviving live occurrence of the retired complete-gravity-output
sentence, all-genera partition-function item, genus-\(g\)
partition-determined sentence, partition-uniqueness sentence,
one-loop gravity title, or old genus-power Gaussian formula. Remaining
external hits are archived Vol I orphan/audit markdown snapshots or
negative guard strings, not live manuscript or compute inputs.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A268 repairs the holographic comparison paragraph in
`thqg_symplectic_polarization.tex`.  Complementarity no longer appears
to compare a boundary genus-\(g\) partition function directly with an
unqualified gravitational path integral.  The paragraph now requires a
physical AdS\(_3\)/CFT\(_2\) comparison datum: a completed
boundary-CFT genus-\(g\) partition function together with a chosen bulk
saddle expansion, or a gravitational functional integral only after
the bulk measure and contour prescription have been specified.
Complementarity supplies only the algebraic bulk--boundary
polarization entering that comparison; it does not construct the bulk
measure, prove modular invariance, select saddles, or imply a physical
entropy formula.

The two BTZ headings in `thqg_symplectic_polarization.tex` now use
"comparison package" rather than "partition function."  Propagation
also corrected the active `thqg_perturbative_finiteness.tex` index
entry from `BTZ black hole!partition function` to
`BTZ black hole!scalar shadow comparison`, matching the local theorem
which proves only the BTZ-normalized scalar shadow series and
explicitly stops before the analytic BTZ partition function.

`compute/tests/test_deletion_ledger.py` now includes
`test_symplectic_polarization_holography_uses_comparison_datum()` and
strengthens
`test_thqg_finiteness_shadow_free_energy_not_path_integral()`.  The
guards require comparison-datum wording, BTZ comparison-package
headings, the scalar-shadow comparison index, and the explicit
negation of bulk-measure/modular-invariance/saddle-selection/entropy
production, while forbidding the retired dictionary/path-integral
paragraph and stale BTZ partition-function headings/index.

A268 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`52 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript paths, live Vol I and Vol III
manuscript/compute paths, and Vol II compute paths finds no surviving
occurrence of the retired dictionary/path-integral sentence, old
complementarity disclaimer, BTZ partition-function headings, or stale
BTZ partition-function index outside negative guard strings.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A269 repairs the genus-one entanglement paragraph in
`thqg_symplectic_polarization.tex`.  The symplectic-polarization
surface no longer claims that the entanglement entropy of the
bulk--boundary system is computed by
\(-\mathrm{Tr}(H_\cA^{(1)}\log H_\cA^{(1)})\), and no longer obtains
a reduced density matrix by tracing over \(Q_1(\cA^!)\).  The
manuscript now says what complementarity proves:
\[
C_1(\cA)=\mathrm{rank}\,H_\cA^{(1)}=\dim Q_1(\cA)=1.
\]
This is a genus-one Lagrangian capacity.  A von Neumann entropy
statement is conditional on a physical comparison datum realizing
\(Q_1(\cA)\) and \(Q_1(\cA^!)\) as trace-class Hilbert factors and
supplying a normalized positive state.  Complementarity itself does
not construct a density matrix, a partial trace, or a Hilbert tensor
product.

A269 also repairs the Ryu--Takayanagi paragraph: Hessian rank is now
capacity, not area.  A physical RT or quantum-extremal-surface reading
requires a boundary subregion, bulk metric, homology constraint, and
variational minimal-surface problem.  The genus-one theorem was also
tightened: it records a single scalar degree-two channel and vanishing
nonlinear shadow projections, rather than calling the potential
quadratic while writing a linear formula.  The stale metadata claim in
`metadata/claims.jsonl` and `metadata/dependency_graph.dot` now reads
"Genus-\(1\) Lagrangian capacity bound" with conditional status.

`compute/tests/test_deletion_ledger.py` now includes
`test_symplectic_polarization_entropy_is_capacity_bound_not_entropy()`,
and the A268 comparison-datum guard was updated for the new section
title.  The guard requires the capacity computation, Hilbert-completion
hypothesis, direct-sum-not-tensor-product scope, RT comparison-datum
language, and forbids the retired entropy-from-Hessian, reduced-density,
Hessian-as-area, and BTZ-entanglement summary forms.

A269 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`53 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript and metadata surfaces, Vol II compute outside
negative guards, and live Vol I/Vol III manuscript/compute/metadata
surfaces finds no surviving occurrence of the retired local forms.
Vol II `make fast` converges after two passes with zero undefined
citations and references.  `make verify-licensing` exits 0 with zero
blocking violations; its printed missing-tag inventory is the
pre-existing nonblocking audit surface.

A270 repairs the complementarity-potential gravitational-action
remark in `thqg_symplectic_polarization.tex`.  The manuscript no
longer says that \(S_\cA\) plays the role of the gravitational action,
that its Taylor jets are gravitational coupling constants, or that the
Virasoro shadow tower produces gravitational counterterms.  It now
states the mathematical object: \(S_\cA\) is the formal function whose
graph is the dual Lagrangian, and whose Taylor coefficients are the
shadow jets of the bar-intrinsic Maurer--Cartan element.  A physical
action/coupling/counterterm reading requires a separate comparison
datum: bulk field space, BV/gauge-fixing data, boundary polarization,
observable normalization, and renormalization scheme.

The shadow-depth table now records algebraic roles and possible
post-comparison channels: quadratic shadow channel, cubic shadow
channel, quartic/contact channel, and higher vertex or composite
channels after comparison.  The G/L/C/M classes now classify algebraic
obstruction length, not perturbative complexity of a metric
gravitational dual.  The Virasoro/W\(_N\) statement is
class-\(\mathsf M\) algebra, not perturbative non-renormalizability of
3d gravity.

Propagation repaired two live echoes.  In `3d_gravity.tex`, the four
terms of \(m_3(T,T,T)\) are now algebraic collision-channel roles;
the central cubic channel is proportional to \(1/G\) only after the
Brown--Henneaux normalization is installed.  In
`thqg_fm_calculus_extensions.tex`, the central cubic term is a central
cubic comparison channel arising from triple collision in
\(\FM_3(\C)\), comparable with a cubic gravitational vertex only after
a metric comparison datum is fixed.

`compute/tests/test_deletion_ledger.py` now includes
`test_symplectic_polarization_potential_is_comparison_data_not_action()`.
The guard requires the formal-function statement, physical comparison
datum, post-comparison channel table, algebraic obstruction-length
classification, and propagated central-cubic comparison wording, while
forbidding the retired gravitational-action, coupling-constant,
three-graviton-coupling, perturbative-complexity,
non-renormalizability, and counterterm slogans.

A270 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`54 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript/metadata surfaces, Vol II compute outside
tests/audit, and live Vol I/Vol III manuscript/compute/metadata
surfaces finds no surviving occurrence of the retired local forms
outside negative guards and historical audit notes.  Vol II
`make fast` converges after two passes with zero undefined citations
and references.  `make verify-licensing` exits 0 with zero blocking
violations; its printed missing-tag inventory is the pre-existing
nonblocking audit surface.

A271 repairs the selected-3-manifold examples in
`thqg_perturbative_finiteness.tex`.  The subsection is no longer "The
partition function on specific 3-manifolds"; it is now
"Scalar-shadow comparison with selected 3-manifold saddles."  The
solid-torus example names the boundary character
\(Z_{\mathrm{char}}(\cA;\tau)\), and the Virasoro vacuum-character
factor is compared with the holomorphic one-loop graviton determinant
only after thermal-AdS boundary-mode normalization is fixed.

The genus-two handlebody example no longer states that
\(Z_{H_2}=\sum_\alpha |\mathcal F_\alpha(\Omega)|^2\) is the physical
handlebody partition function, no longer assumes vacuum-block
dominance, and no longer identifies \(F_2=7c/11520\) with a two-loop
handlebody correction.  It now names the missing physical data:
conformal-block pairing, measure or sum over intermediate labels,
handlebody saddle prescription, and metric action normalization.  The
block norm \(Z_{H_2}^{\mathrm{block}}\) and large-\(c\) metric ansatz
are comparison data.  The scalar coefficient
\[
F_2^{\mathrm{scal}}(\mathrm{Vir}_c)=7c/11520
\]
is only the connected genus-two Hodge-trace shadow coefficient on the
uniform-weight scalar lane, not a computed two-loop handlebody
determinant and not a physical handlebody partition function.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_handlebody_example_is_comparison_not_partition_function()`.
The guard requires the comparison subsection title, boundary-character
notation, thermal-AdS normalization sentence, handlebody comparison
inputs, \(Z_{H_2}^{\mathrm{block}}\), \(S_{H_2}^{\mathrm{met}}\),
\(F_2^{\mathrm{scal}}\), and the explicit negative sentence, while
forbidding the retired partition-function, vacuum-dominance,
\(S_{\mathrm{grav}}\), and two-loop-correction formulations.

A271 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`55 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript/metadata surfaces, Vol II compute outside
tests/audit, and live Vol I/Vol III manuscript/compute/metadata
surfaces finds no surviving occurrence of the retired local forms
outside negative guards.  Vol II `make fast` converges after two
passes with zero undefined citations and references.  `make
verify-licensing` exits 0 with zero blocking violations; its printed
missing-tag inventory is the pre-existing nonblocking audit surface.

A272 repairs the two-parameter convergence theorem in
`thqg_perturbative_finiteness.tex`.  The theorem no longer names the
genus-summed object \(Z_{\mathrm{grav}}(T;\hbar,t)\) or indexes the
statement as a gravitational partition-function convergence theorem.
It is now "Completed degree-summed shadow two-parameter convergence,"
with the index `shadow partition series!two-parameter convergence`.

The fixed-genus summand is now
\[
Z_g^{\mathrm{sh,deg}}(T;t)
=\sum_{r\geq0}t^r
\int_{\overline{\mathcal M}_{g,n(g,r)}}
\operatorname{Sh}_{g,r}(\Theta_\cA),
\]
and the Gaussian genus-summed object is
\[
Z_{\mathrm{sh}}^{\mathrm{deg}}(T;\hbar,t)
=\sum_{g\geq0}\hbar^g
Z_g^{\mathrm{sh,deg}}(T;t).
\]
The negative conclusion now says that HS-sewing and recursive
existence do not provide a bidisk for the completed degree-summed
shadow genus expansion outside the explicit degree estimates and the
Gaussian scalar lane.  The proof sentence was repaired accordingly:
the estimates make the degree-weighted shadow trace analytic in \(t\),
not a physical scalar trace.

A272 also updates the stale metadata rows and dependency-graph nodes
for `thm:thqg-I-full-convergence` and
`thm:thqg-I-2d-convergence`, so the metadata no longer advertises a
"completed full partition function" or the old two-parameter
gravitational title.

`compute/tests/test_deletion_ledger.py` now strengthens
`test_thqg_degree_summed_series_not_full_gravity_partition_function()`
to require the two-parameter shadow title, shadow index,
\(Z_g^{\mathrm{sh,deg}}(T;t)\),
\(Z_{\mathrm{sh}}^{\mathrm{deg}}(T;\hbar,t)\), the degree-weighted
shadow-trace proof sentence, and the completed shadow-genus bidisk
wording, while forbidding the retired gravitational title, index,
\(Z_{\mathrm{grav}}(T;\hbar,t)\), full-degree-summed bidisk wording,
and scalar-trace-in-\(t\) sentence.

A272 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`55 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript/metadata surfaces, Vol II compute outside
tests/audit, and live Vol I/Vol III manuscript/compute/metadata
surfaces finds no surviving occurrence of the retired two-parameter
forms outside negative guards.  A broader Vol I scan sees the old
strings only in the historical
`fix_wave_B_20260413_171623/B11_orphaned_chapters.md` archive, not in
the live Vol I surface.  Vol II `make fast` converges after two passes
with zero undefined citations and references.  `make verify-licensing`
exits 0 with zero blocking violations; its printed missing-tag
inventory is the pre-existing nonblocking audit surface.

A273 repairs the scalar closed-form cluster in
`thqg_perturbative_finiteness.tex`.  The repaired principle is:
before a physical comparison datum is supplied, the closed form
\[
Z_{\mathrm{sh}}^{\mathrm{scal}}(\cA;\hbar)
=\kappaChHodge(\cA)\frac{\sqrt{\hbar}/2}{\sin(\sqrt{\hbar}/2)}
\]
is a scalar Hodge-trace shadow series, not a gravitational partition
function.

The convergence-radius proposition now reads "Convergence radius of the
scalar shadow closed form" and its index is
`convergence radius!scalar shadow closed form`.  The main finiteness
theorem now uses \(Z_{\mathrm{sh}}^{\mathrm{scal}}\) for the completed
scalar series.  The pole theorem, partial fraction formula, shadow
entropy definition, quantitative table, tensor-product formula,
large-\(N\) \(\mathcal W_N\) limit, and analytic-extension loci now use
\(\kappaChHodge\), \(Z_{\mathrm{sh}}^{\mathrm{scal}}\), and
\(R_{\mathrm{sh}}^{\mathrm{scal}}\) rather than bare \(\kappa\),
\(Z_{\mathrm{grav}}\), or \(R_{\mathrm{grav}}\).

The large-\(c\) and Newton statements were retitled as
Brown--Henneaux-normalized scalar-shadow statements.  The
Brown--Henneaux substitution \(c=3\ell/(2G)\), \(\hbar=2G/(3\ell)\)
now gives the scalar shadow closed form, not a constructed
gravitational partition function.

The non-perturbative conjecture was retitled "Non-perturbative
comparison datum."  It now states that a physical completion requires a
bulk completion, contour, measure, unitarity/positivity structure, and
saddle sector.  Positivity for real \(\hbar>0\) belongs to that
completed physical datum; it is not a consequence of the meromorphic
scalar shadow closed form alone.  The dictionary table now calls
\(Q^{\mathrm{contact}}\) a quartic/contact comparison channel rather
than a two-loop graviton self-energy.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_scalar_closed_forms_are_shadow_data_before_comparison()`.
The guard requires the scalar-shadow titles, Brown--Henneaux
scalar-shadow normalization, non-perturbative comparison-datum wording,
measure/unitarity/saddle inputs, positivity-datum sentence, branch
choice for shadow entropy, \(Z_{\mathrm{sh}}^{\mathcal W_N,\mathrm{scal}}\),
\(R_{\mathrm{sh}}^{\mathrm{scal}}\), and repaired metadata titles.  It
forbids the retired scalar-gravitational partition-function,
gravitational-coupling/Newton, unlicensed completion/positivity,
\(Z_{\mathrm{grav,scal}}\), \(R_{\mathrm{grav}}\), and
two-loop-graviton-self-energy forms.

A273 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`56 passed` for `test_deletion_ledger.py`.  Fixed-string grep over
live Vol II manuscript/metadata surfaces, Vol II compute outside
tests/audit, and live Vol I/Vol III manuscript/compute/metadata
surfaces finds no surviving occurrence of the retired
scalar-gravitational forms.  The only Vol II compute hits are the
intentional negative guard strings.  Vol II `make fast` converges after
two passes with zero undefined citations and references.  `make
verify-licensing` exits 0 with zero blocking violations; its printed
missing-tag inventory is the pre-existing nonblocking audit surface.

A274 repairs the Faber--Pandharipande scalar-lane notation in
`thqg_perturbative_finiteness.tex`.  The repaired principle is:
\[
F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}
\]
on the proved uniform-weight scalar lane; the scalar is the Hodge-row
component of the \(\kappa\)-tuple, not a bare untyped \(\kappa\).

The proof of the closed-form free energy, the universality remark, the
\(\hat A\)-genus generating-function derivation, the two generating
functions remark, the absolute-convergence theorem, its effective
two-sided bounds, the pole residues, the genus \(1\) through \(10\)
explicit values, and the specialization table now use
\(\kappaChHodge(\cA)\) or the relevant specialization
\(\kappaChHodge(\mathcal H_k)\), \(\kappaChHodge(\widehat{\mathfrak{sl}}_{2,k})\),
or \(\kappaChHodge(\mathrm{Vir}_c)\).

The same repair was propagated inside the target file to the algebraic
shadow finiteness proof, the genus-\(g\) amplitude decomposition, the
Brown--Henneaux semiclassical paragraph, the BTZ-normalized scalar
shadow proof, the mechanism summary, the effective genus bound and tail
bound, the low-degree primitive bound, the tensor-product linearity
proof, the critical-level, lattice, free-boson, Virasoro self-dual, and
multi-weight diagonal discussions.  Sectoral symbols
\(\kappa_i,\kappa_T,\kappa_W,\kappa_{\mathrm{BKM}}\) remain, because
they name components rather than the Hodge-row scalar.

`compute/tests/test_deletion_ledger.py` now includes
`test_thqg_faber_pandharipande_scalar_lane_uses_kappa_ch_hodge()`.  The
guard slices the FP scalar-lane cluster and requires the repaired
\(\kappaChHodge(\cA)\) formula, generating-function derivation,
absolute-value bound, residue, explicit \(F_1\) value, and table
caption.  It forbids the retired bare-\(\kappa\) formulas in the cluster
and exact `\kappa \cdot \lambda_g^{\mathrm{FP}}`, `\kappa/24`, and
`|\kappa|` stale strings in the target file.

A274 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`57 passed` for `test_deletion_ledger.py`.  Fixed-string grep confirms
`thqg_perturbative_finiteness.tex` has no surviving exact retired
bare-\(\kappa\) FP-lane forms; remaining \(\kappa\) matches in that file
are sectoral components.  A broader Vol II/Vol I/Vol III scan found
many still-live bare-\(\kappa\) FP-lane echoes outside this target file,
so global propagation remains a follow-on obligation.  Vol II `make
fast` converges after two passes with zero undefined citations and
references.  `make verify-licensing` exits 0 with zero blocking
violations and zero warnings.

A275 propagates the A274 repair across the live Vol II FP scalar-lane
echo surface.  The repaired form is again
\[
F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\]
with genus-one and coefficient specializations using
\(\kappaChHodge(\cA)\) or the relevant family specialization rather
than an untyped scalar \(\kappa\).

The propagation touches the comparison tables and examples in
`holomorphic_topological.tex`, `thqg_holographic_reconstruction.tex`,
`thqg_fredholm_partition_functions.tex`, `3d_gravity.tex`,
`thqg_3d_gravity_movements_vi_x.tex`, `thqg_bv_ht_extensions.tex`,
`thqg_gravitational_complexity.tex`,
`thqg_gravitational_s_duality.tex`,
`thqg_bv_construction_extensions.tex`, `rosetta_stone.tex`,
`w-algebras-stable.tex`, `modular_swiss_cheese_operad.tex`,
`modular_pva_quantization_core.tex`, `bv_brst.tex`,
`foundations_recast_draft.tex`, `ht_bulk_boundary_line.tex`,
`ht_bulk_boundary_line_frontier.tex`,
`thqg_celestial_holography_extensions.tex`,
`thqg_modular_bootstrap.tex`,
`thqg_spectral_braiding_extensions.tex`, `examples-worked.tex`,
`relative_feynman_transform.tex`, and the compute docstring surface in
`compute/lib/f10_resurgence_cross_channel.py`.

`compute/tests/test_deletion_ledger.py` now includes
`test_live_vol2_fp_lane_echoes_use_kappa_ch_hodge()`.  The guard scans
live Vol II chapters, appendices, and compute libraries for exact and
regex variants of the retired bare-\(\kappa\) FP scalar-lane formulas,
including \(F_g=\kappa\lambda_g^{\mathrm{FP}}\),
\(\kappa\cdot\lambda_g^{\mathrm{FP}}\), \(\kappa/24\), and
\(\kappa c_g\).  It intentionally leaves subscripted effective
parameters such as \(\kappa_{\mathrm{eff}}\) for a later semantic audit.

A275 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py` and
`compute/lib/f10_resurgence_cross_channel.py`; focused pytest reports
`58 passed` for `test_deletion_ledger.py`.  Fixed-string and regex grep
over live Vol II `chapters`, `appendices`, and `compute/lib` finds no
surviving exact bare-\(\kappa\) FP scalar-lane formulas outside
intentional negative test strings.  Vol II `make fast` converges after
two passes with zero undefined citations and references.  `make
verify-licensing` exits 0 with zero blocking violations and zero
warnings.  Broader cross-volume grep still finds stale Vol I echoes and
Vol II subscripted effective scalar forms; those remain follow-on
obligations, not solved by A275.

A276 resolves the Vol II subscripted effective-scalar part of that
follow-on obligation in the gravity lane.  The key point is not to
delete \(\kappa_{\mathrm{eff}}\), but to type it:
\[
\kappa_{\mathrm{eff}}
:=\kappaChHodge(\mathrm{Vir}_c\otimes bc_2)
=\kappaChHodge(\mathrm{Vir}_c)+\kappaChHodge(bc_2)
=(c-26)/2.
\]
Thus \(\kappa_{\mathrm{eff}}\) is the Hodge scalar of the coupled
matter--ghost package, not a new primitive scalar and not the intrinsic
Virasoro invariant \(\kappaChHodge(\mathrm{Vir}_c)=c/2\).

The exposed effective FP formulas in `3d_gravity.tex` and
`thqg_3d_gravity_movements_vi_x.tex` now use
\[
F_g^{\mathrm{eff}}(\mathrm{Vir}_c\otimes bc_2)
=\kappaChHodge(\mathrm{Vir}_c\otimes bc_2)
\lambda_g^{\mathrm{FP}},
\]
with \(\kappa_{\mathrm{eff}}\) only as the following abbreviation.  The
de Sitter proposition now writes
\(\kappa_{\mathrm{dS}}
=\kappaChHodge(\mathrm{Vir}_{c_{\mathrm{dS}}})\) and
\[
F_g^{\mathrm{dS}}(\mathrm{Vir}_{c_{\mathrm{dS}}})
=\kappaChHodge(\mathrm{Vir}_{c_{\mathrm{dS}}})
\lambda_g^{\mathrm{FP}}.
\]
This separates intrinsic Virasoro curvature, matter--ghost effective
curvature, and de Sitter analytic-continuation curvature.

`compute/tests/test_deletion_ledger.py` now includes
`test_gravity_effective_fp_lane_names_hodge_scalar_source()`.  The guard
requires the matter--ghost Hodge scalar definition, the typed effective
FP formula, the typed de Sitter FP formula, and the typed all-genera
generating function.  It forbids the retired standalone effective and
de Sitter FP formulas.

A276 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`; focused pytest reports
`59 passed` for `test_deletion_ledger.py`.  Exact and regex grep over
live Vol II `chapters`, `appendices`, `compute/lib`, and tests finds no
surviving standalone
`F_g = \kappa_{\mathrm{eff}}\lambda_g^{\mathrm{FP}}`,
`F_g^{\mathrm{eff}} = \kappa_{\mathrm{eff}}\lambda_g^{\mathrm{FP}}`,
or `F_g = \kappa_{\mathrm{dS}}\lambda_g^{\mathrm{FP}}` formula.  The
remaining \(\kappa_{\mathrm{eff}}\lambda_g^{\mathrm{FP}}\) matches in
the touched files are trailing abbreviation equalities immediately
after a displayed \(\kappaChHodge(\cdots)\) source formula.  Vol II
`make fast` converges after two passes with zero undefined citations and
references.  `make verify-licensing` exits 0 with zero blocking
violations and zero warnings.  Vol I still has a broad bare-\(\kappa\)
FP backlog from A275 and archived effective-notation echoes; A276 does
not claim those are solved.

A277 begins that Vol I propagation backlog with the live Heisenberg
example cluster.  The target file already states
\(\kappa(\cH_\kappa)=\kappa\), so the repair is not to rename the level
parameter.  The repair is to make the invariant map visible in every
FP scalar-lane formula:
\[
F_g(\cH_\kappa)
=\kappa(\cH_\kappa)\lambda_g^{\mathrm{FP}}.
\]

In `/Users/raeez/chiral-bar-cobar/chapters/examples/heisenberg_eisenstein.tex`,
the genus-expansion table, complementarity remark, genus-2 Bernoulli
computation, scalar projection paragraph, tautological genus-\(g\)
display, open/closed MC scalar projection, and clutching identity now
use \(\kappa(\cH_\kappa)\) or
\(\kappa(\mathcal H_\kappa)\).  The OPE parameter remains
\(\kappa\), and the final equality
\(\kappa(\cH_\kappa)=\kappa\) supplies the numerical level value.

Vol I now has the guard
`compute/tests/test_heisenberg_fp_scalar_typing.py`, with
`test_heisenberg_fp_lane_names_modular_characteristic_source()`.  It
requires the typed Heisenberg FP and scalar-MC formulas and forbids the
retired compact and dotted bare-\(\kappa\) FP spellings in
`heisenberg_eisenstein.tex`.

A277 verification: `python3 -m py_compile` passes on the new Vol I
guard; focused pytest reports `1 passed`.  Targeted grep confirms
`heisenberg_eisenstein.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no
`\kappa\lambda_g^{\mathrm{FP}}`, and no compact
`F_g(...)=\kappa\lambda_g^{\mathrm{FP}}` FP scalar formula.  Vol I
`make fast` reaches a stable warning state after three passes with zero
undefined citations and one persistent undefined reference from the
pre-existing manuscript surface, not the touched Heisenberg lines.
Vol I `make verify` exits 0 with all checks passed.  The broader Vol I
bare-\(\kappa\) FP backlog remains open outside the Heisenberg cluster.

A278 extends the Vol I propagation from the Heisenberg example to the
central `genus_expansions.tex` surface.  This chapter is the local
source for the scalar-lane formula, the genus-one all-family clause, the
standard-landscape table, and the warning that interacting multi-weight
families acquire cross-channel corrections.  The repair therefore had
to preserve both statements:
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
\]
on the proved uniform-weight/free-field scalar lane, and failure of the
same scalar formula for interacting multi-weight families at
\(g\ge2\).

In `/Users/raeez/chiral-bar-cobar/chapters/examples/genus_expansions.tex`,
the lattice proof, Verlinde comparison, bivariate generating function,
universal generating-function proof, free-field complementarity proof,
standard-landscape table header, planted-forest warning, growth-rate
remark, scalar-duality theorem, and affine loop-expansion proposition
now use typed forms such as \(F_g(\cA)=\kappa(\cA)
\lambda_g^{\mathrm{FP}}\),
\(F_1(\cA)=\kappa(\cA)/24\), and
\(F_g(\widehat{\fg}_k)=\kappa(\widehat{\fg}_k)
\lambda_g^{\mathrm{FP}}\).  The \(F_2,F_3\) columns in the landscape
table are now labelled scalar entries
\(F_2^{\mathrm{sc}},F_3^{\mathrm{sc}}\), matching the surrounding
multi-weight warning.

Vol I now has the guard
`compute/tests/test_genus_expansions_scalar_typing.py`, with
`test_genus_expansions_fp_lane_uses_typed_kappa()`.  It requires the
typed generic, Heisenberg, affine, and genus-one formulas and forbids
the retired bare-\(\kappa\) FP spellings in `genus_expansions.tex`.

A278 verification: `python3 -m py_compile` passes on the two Vol I
scalar-typing guards; focused pytest reports `2 passed` for
`test_heisenberg_fp_scalar_typing.py` and
`test_genus_expansions_scalar_typing.py`.  Targeted grep confirms
`genus_expansions.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no `\kappa/24`, no
`F_g = \kappa...`, and no `F_1 = \kappa...` scalar formula.  Vol I
`make fast` reaches the known stable warning state after three passes
with zero undefined citations and one persistent undefined reference
`def:kappa-scalar-vs-mumford`, outside the edited lines.  Vol I
`make verify` exits 0 with all checks passed.  Broader Vol I live files
still contain bare-\(\kappa\) FP echoes outside
`genus_expansions.tex`.

A279 continues the same Vol I propagation into the standard-family
landscape census.  The target is not the numerical table values
themselves: rows such as the Heisenberg row may still display the
evaluated value \(\kappa/24\).  The target is the theorem-level and
caption-level formula in which the modular characteristic had been
written as an untyped scalar:
\[
F_g=\kappa\lambda_g^{\mathrm{FP}}.
\]

In `/Users/raeez/chiral-bar-cobar/chapters/examples/landscape_census.tex`,
the free-energy landscape caption and headers now state
\[
F_1(\cA)=\kappa(\cA)/24,\qquad
F_2^{\mathrm{scalar}}(\cA)=7\kappa(\cA)/5760.
\]
The scalar lane now reads
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
and the multi-weight lane now includes the typed correction
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
        +\delta F_g^{\mathrm{cross}}(\cA).
\]
The genus-one verification paragraph, \(\theta_1\) formula, Polyakov
comparison, anomaly-ratio table, DS anomaly-ratio proof, and
genus-obstruction table caption now use \(\kappa(\cA)\) and \(c(\cA)\)
at the statement level.  The DS corollary now writes
\[
F_1(\mathcal W^k(\mathfrak g))
  =\kappa(\mathcal W^k(\mathfrak g))/24,\qquad
\kappa(\mathcal W^k(\mathfrak g))
  =c(\mathcal W^k(\mathfrak g))\varrho(\mathfrak g),
\]
so the central-charge denominator and modular characteristic source
remain typed.

Vol I now has the guard
`compute/tests/test_landscape_census_scalar_typing.py`, with
`test_landscape_census_scalar_lane_types_modular_characteristic()`.  It
requires the typed landscape caption, scalar lane, multi-weight
correction, \(\theta_1(\cA)\) source, and
\(\mathcal W^k(\mathfrak g)\) anomaly-ratio formulas; it forbids the
retired caption/prose spellings with bare \(\kappa\).

A279 verification: `python3 -m py_compile` passes on all three Vol I
scalar-typing guards.  Focused pytest reports `3 passed` for
`test_landscape_census_scalar_typing.py`,
`test_genus_expansions_scalar_typing.py`, and
`test_heisenberg_fp_scalar_typing.py`.  Targeted grep confirms
`landscape_census.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no `F_g = \kappa...`, no
`F_1 = \kappa...`, and no `\theta_1 = \kappa \cdot \mu` scalar-source
formula.  Vol I `make fast` reaches the known stable warning state
after three passes with zero undefined citations and one persistent
undefined reference outside the edited lines.  Vol I `make verify`
exits 0 with all checks passed.  Broader Vol I live files still contain
bare-\(\kappa\) FP echoes outside the Heisenberg, genus-expansion, and
landscape-census surfaces.

A280 moves the same scalar-characteristic discipline into the active
Vol I \(Y\)-algebra chapter.  The local theorem already proves
\[
\kappa(Y_{1,1,1}[\Psi])=\Psi,
\]
so the repair is not to remove the coupling \(\Psi\).  The repair is to
show the invariant map before evaluating it at \(\Psi\), and to correct
the Faber--Pandharipande coefficient from an unmarked
\(\overline{\cM}_g\) phrase to the marked scalar coefficient.

In `/Users/raeez/chiral-bar-cobar/chapters/examples/y_algebras.tex`,
the class-\(\mathsf G\) remark, genus-one proof, higher-genus tower
proposition, multi-weight proof, and HT-origin datum now write
\[
F_g(Y_{1,1,1}[\Psi])
  =\kappa(Y_{1,1,1}[\Psi])\lambda_g^{\mathrm{FP}}
  =\Psi\lambda_g^{\mathrm{FP}}.
\]
The multi-weight proof now displays
\[
F_g(Y_{1,1,1}[\Psi])
 =\kappa(Y_{1,1,1}[\Psi])\lambda_g^{\mathrm{FP}}
 +\delta F_g^{\mathrm{cross}}(Y_{1,1,1}[\Psi]),
\]
then kills the cross-channel term by the \(c=0\) \(T\)-channel and
\(\kappa_T=0\) argument.  The coefficient is now defined as
\[
\lambda_g^{\mathrm{FP}}
 =\int_{\overline{\cM}_{g,1}}\psi_1^{2g-2}\lambda_g,
\]
not as an unmarked integral on \(\overline{\cM}_g\).  The HT-origin
paragraph now names
\(\kappa(Y_{1,1,1}[\Psi])=\Psi\) for the trivalent case and
\(\kappa(A)\omega_g\) for the general genus-\(g\) curvature.

Vol I now has the guard
`compute/tests/test_y_algebras_scalar_typing.py`, with
`test_y111_genus_tower_types_kappa_source_and_fp_coefficient()`.  It
requires the typed \(Y_{1,1,1}\) genus-one formula, the
\(\kappa(Y_{1,1,1}[\Psi])=\Psi\) source, the typed genus-\(g\) scalar
formula, the typed cross-channel correction, the marked
Faber--Pandharipande coefficient, and the HT-origin curvature wording;
it forbids the retired bare-\(\kappa\) and unmarked-FP descriptions.

A280 verification: `python3 -m py_compile` passes on all four Vol I
scalar-typing guards.  Focused pytest reports `4 passed` for
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`,
`test_genus_expansions_scalar_typing.py`, and
`test_heisenberg_fp_scalar_typing.py`.  Targeted grep confirms
`y_algebras.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no `F_g = \kappa...`, no
`F_1 = \kappa...`, and no
`intersection numbers on~$\overline{\cM}_g$` FP-description formula.
Vol I `make fast` reaches the known stable warning state after three
passes with zero undefined citations and one persistent undefined
reference outside the edited lines.  Vol I `make verify` exits 0 with
all checks passed.  Broader Vol I live files still contain
bare-\(\kappa\) FP echoes outside the Heisenberg, genus-expansion,
landscape-census, and \(Y_{1,1,1}\) surfaces.

A281 returns to the active Vol I Heisenberg chapter to remove the last
marked-integral scalar-source echoes that escaped A277.  The level
parameter remains \(\kappa\), but the class source must be the invariant
map \(\kappa(\mathcal H_\kappa)\) until the final numerical evaluation.

In `/Users/raeez/chiral-bar-cobar/chapters/examples/heisenberg_eisenstein.tex`,
the genus-2 Hodge obstruction now reads
\[
m_0^{(2)}=\kappa(\mathcal H_\kappa)\lambda_2.
\]
The genus-2 free energy now factors through the marked scalar
coefficient as
\[
F_2(\mathcal H_\kappa)
 =\kappa(\mathcal H_\kappa)
  \int_{\overline{\mathcal M}_{2,1}}\psi_1^2\lambda_2
 =7\kappa(\mathcal H_\kappa)/5760
 =7\kappa/5760.
\]
The determinant-line comparison remark now writes the general scalar
projection as
\[
F_g(\mathcal H_\kappa)
 =\kappa(\mathcal H_\kappa)
  \int_{\overline{\mathcal M}_{g,1}}\psi_1^{2g-2}\lambda_g.
\]
This preserves the numerical Heisenberg evaluation while preventing
the level parameter from masquerading as the invariant map.

The Vol I guard
`compute/tests/test_heisenberg_fp_scalar_typing.py` now also requires
the typed genus-2 Hodge obstruction, the typed genus-2 marked integral,
and the typed genus-\(g\) marked projection.  It forbids the retired
`\kappa \cdot \int_{\overline{\mathcal M}...}` and
\(m_0^{(2)}=\kappa\lambda_2\) forms.

A281 verification: `python3 -m py_compile` passes on all four Vol I
scalar-typing guards.  Focused pytest reports `4 passed` for
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Targeted grep confirms
`heisenberg_eisenstein.tex` has no surviving exact
`\kappa \cdot \int`, no `m_0^{(2)}=\kappa\,\lambda_2`, no
`F_g = \kappa...`, and no `\kappa \cdot \lambda_g^{\mathrm{FP}}`
scalar formula.  Vol I `make fast` reaches the known stable warning
state after three passes with zero undefined citations and one
persistent undefined reference outside the edited lines.  Vol I
`make verify` exits 0 with all checks passed.  Broader Vol I live files
still contain bare-\(\kappa\) FP echoes outside the four converged
scalar-typing surfaces.

A282 moves the scalar-characteristic typing into the active Vol I
free-fields chapter, concentrating on the free fermion.  The free
fermion has fixed scalar modular characteristic
\[
\kappa(\cF)=\frac14,
\]
so the repair is to show the invariant map before evaluating it, not to
erase the numerical value.

In `/Users/raeez/chiral-bar-cobar/chapters/examples/free_fields.tex`,
the free-fermion shadow metric now reads
\[
Q_{\mathcal F}(t)=(2\kappa(\mathcal F))^2=\frac14,\qquad
S_2=\kappa(\mathcal F)=\frac14.
\]
The genus proof now writes
\[
\mathrm{obs}_g(\mathcal F)
 =\kappa(\mathcal F)\lambda_g,\qquad
F_g(\mathcal F)
 =\kappa(\mathcal F)\lambda_g^{\mathrm{FP}},
\]
and defines
\[
\lambda_g^{\mathrm{FP}}
 =\int_{\overline{\mathcal M}_{g,1}}\psi_1^{2g-2}\lambda_g
 =\frac{2^{2g-1}-1}{2^{2g-1}}\frac{|B_{2g}|}{(2g)!}.
\]
The spin-structure discussion now separates the spin-dependent
Pfaffian partition function from the scalar shadow:
\[
F_g(\cF)=\kappa(\cF)\lambda_g^{\mathrm{FP}}
       =\frac14\lambda_g^{\mathrm{FP}}.
\]
The genus-one theorem now reads
\[
F_1(\cF)=\kappa(\cF)\lambda_1^{\mathrm{FP}}
        =\kappa(\cF)/24.
\]
The Heisenberg comparison table still contains the evaluated
Heisenberg row value \(\kappa/24\); this is not a theorem-level
free-fermion invariant formula.

Vol I now has the guard
`compute/tests/test_free_fields_scalar_typing.py`, with
`test_free_fermion_scalar_lane_types_modular_characteristic()`.  It
requires the typed shadow metric, typed \(S_2\), typed obstruction and
free-energy formulas, the marked FP coefficient, the spin-structure
scalar formula, and the typed genus-one formula.  It forbids the
retired bare-\(\kappa\) scalar-source forms and the free-fermion
\(\kappa=1/4\) shorthand.

A282 verification: `python3 -m py_compile` passes on all five Vol I
scalar-typing guards.  Focused pytest reports `5 passed` for
`test_free_fields_scalar_typing.py`,
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Targeted grep confirms
`free_fields.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no `F_g = \kappa...`, no
`F_1 = \kappa...`, and no `\kappa = 1/4` scalar-source formula.  Vol I
`make fast` reaches the known stable warning state after three passes
with zero undefined citations and one persistent undefined reference
outside the edited lines.  Vol I `make verify` exits 0 with all checks
passed.  Broader Vol I live files still contain bare-\(\kappa\) FP
echoes outside the five converged scalar-typing surfaces.

A283 moves the same scalar-characteristic discipline into the active
Vol I W-algebra chapter, at the AGT shadow correspondence.  The old
sentence multiplied the scalar source by an all-weight/cross-channel
package and then called only the scalar factor the shadow amplitude.
This made three different objects occupy one formula: the scalar
constant-map shadow, the multi-weight cross-channel shadow, and the
representation/instanton-dependent Nekrasov terms.

In `/Users/raeez/chiral-bar-cobar/chapters/examples/w_algebras.tex`,
the theorem now writes the symmetric-limit genus coefficient as
\[
F_g^{\mathrm{Nek}}(\vec a,q)
=F_g^{\mathrm{sc}}(\mathcal W_N)
 +\delta F_g^{\mathrm{cross}}(\mathcal W_N)
 +F_g^{\mathrm{rep}}(\vec a,q),
\qquad
F_g^{\mathrm{sc}}(\mathcal W_N)
=\kappa(\mathcal W_N)\lambda_g^{\mathrm{FP}}.
\]
The following prose names the first summand as the scalar shadow
amplitude, the second as the multi-weight cross-channel shadow, and
the third as representation/instanton data.  Proof part(ii) now says
only the scalar part of the Nekrasov genus expansion reduces to the
constant-map contribution and its \(\hat A\)-genus coefficient; the
non-scalar all-weight contribution is the cross-channel term.

Vol I now has the guard
`compute/tests/test_w_algebras_scalar_typing.py`, with
`test_w_algebras_agt_shadow_decomposition_types_scalar_source()`.  It
requires the typed W-algebra scalar source
\(\kappa(\mathcal W_N)=c\cdot(H_N-1)\), the three-term Nekrasov
decomposition, the scalar/cross/representation split, and the proof
wording.  It forbids the retired bare-\(\kappa\) AGT product, the old
`ALL-WEIGHT + \delta` package, and the "leading universal term"
phrase.

A283 verification: `python3 -m py_compile` passes on all six Vol I
scalar-typing guards.  Focused pytest reports `6 passed` for
`test_w_algebras_scalar_typing.py`,
`test_free_fields_scalar_typing.py`,
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Targeted grep confirms
`w_algebras.tex` has no surviving exact
`\kappa \cdot \lambda_g^{\mathrm{FP}}`, no `F_g = \kappa...`, no
retired `ALL-WEIGHT + \delta` product, and no "leading universal
term" phrase.  Vol I `make fast` reaches the known stable warning
state after three passes with zero undefined citations and one
persistent undefined reference outside the edited lines.  Vol I
`make verify` exits 0 with all checks passed.

A284 moves the scalar-characteristic typing into the active Vol I
minimal-model chapter, at the Ising scalar tower.  The old passage
wrote the scalar-level free energy as
\(F_g(\mathrm{Ising})=\kappa\lambda_g^{\mathrm{FP}}\) and then
evaluated it at \(\kappa=1/4\).  That notation hid three distinctions:
the typed Virasoro characteristic
\(\kappa^{\mathrm{Vir}}(\mathrm{Vir}_{1/2})\), the scalar
uniform-weight tower, and the full genus amplitude corrected by
higher-degree shadow terms for \(g\geq2\).

In `/Users/raeez/chiral-bar-cobar/chapters/examples/minimal_model_examples.tex`,
the divergence-paradox remark now writes
\[
F_g^{\mathrm{sc}}(\mathrm{Ising})
=\kappa^{\mathrm{Vir}}(\mathrm{Vir}_{1/2})
 \lambda_g^{\mathrm{FP}}.
\]
The scalar-free-energy proposition now writes
\[
F_g^{\mathrm{sc}}(\mathrm{Ising})
=\kappa^{\mathrm{Vir}}(\mathrm{Vir}_{1/2})
 \lambda_g^{\mathrm{FP}}
=\frac14\,\lambda_g^{\mathrm{FP}},
\]
and the table column is \(F_g^{\mathrm{sc}}(\mathrm{Ising})\).  The
proof now derives the coefficient from
\(\kappa^{\mathrm{Vir}}(\mathrm{Vir}_{1/2})=1/4\), rather than from a
bare scalar parameter.

Vol I now has the guard
`compute/tests/test_minimal_model_scalar_typing.py`, with
`test_ising_scalar_free_energy_uses_typed_virasoro_characteristic()`.
It requires the typed Ising scalar source, the scalar-level formula,
the scalar table heading, and the full-amplitude correction sentence.
It forbids the retired bare-\(\kappa\) Ising formula, the
`at \(\kappa=1/4\)` phrasing, the old table heading, and the
"individual genus-\(g\) amplitudes" wording.

A284 verification: `python3 -m py_compile` passes on all seven Vol I
scalar-typing guards.  Focused pytest reports `7 passed` for
`test_minimal_model_scalar_typing.py`,
`test_w_algebras_scalar_typing.py`,
`test_free_fields_scalar_typing.py`,
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Targeted grep confirms no
live occurrence of the retired Ising bare-\(\kappa\) scalar formula,
no active `F_g=(1/4)\lambda_g^{\mathrm{FP}}` proof sentence, and no
active "individual genus-\(g\) amplitudes" phrase outside the negative
guard.  The only `at \(\kappa=1/4\)` hits are historical
`final_gaps_20260413_213946` snapshots, not live manuscript or compute
code.  Vol I `make fast` reaches the known stable warning state after
three passes with zero undefined citations and one persistent
undefined reference outside the edited lines.  Vol I `make verify`
exits 0 with all checks passed.

A285 moves the scalar-characteristic discipline into two active Vol I
connection surfaces.  In
`/Users/raeez/chiral-bar-cobar/chapters/connections/genus_complete.tex`,
the Heisenberg Chern--Weil example no longer lets the level parameter
\(\kappa\) stand for the scalar projection.  The trace now reads
\[
\operatorname{tr}(\Theta_{\mathcal H})
=\sum_{g\geq1}\kappa(\mathcal H_\kappa)\lambda_g^{\mathrm{FP}},
\]
and the text states that \(\kappa(\mathcal H_\kappa)\) evaluates to
the level \(\kappa\) for the rank-one Heisenberg algebra.  The scalar
free-energy line now reads
\[
F_g^{\mathrm{sc}}(\mathcal H_\kappa)
=\kappa(\mathcal H_\kappa)\lambda_g^{\mathrm{FP}},
\]
while the level-dependent spectral and line-kernel channels remain
\(\Delta(x)=1-\kappa x\) and \(r_{\mathcal H}(z)=\kappa/z\).

Propagation found a live theorem-summary survivor in
`/Users/raeez/chiral-bar-cobar/chapters/connections/grand_unification_platonic.tex`.
The obstruction-tower clause now writes
\[
\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g
\]
on the uniform-weight stratum, and the multi-weight scalar lane as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
The corollary, proof-verification paragraph, and summary table carry
the same typed \(\mathrm{obs}_g(\cA)\) statement.

Vol I now has two additional guards:
`compute/tests/test_genus_complete_scalar_typing.py`, with
`test_genus_complete_heisenberg_chern_weil_types_scalar_projection()`,
and `compute/tests/test_grand_unification_scalar_typing.py`, with
`test_grand_unification_types_obstruction_and_cross_channel_scalar_lane()`.
They require the typed Heisenberg scalar projection, the typed
theorem-D obstruction readout, and the scalar/cross-channel split, and
forbid the retired bare-\(\kappa\) FP and theorem-summary forms.

A285 verification: `python3 -m py_compile` passes on all nine Vol I
scalar-typing guards.  Focused pytest reports `9 passed` for
`test_grand_unification_scalar_typing.py`,
`test_genus_complete_scalar_typing.py`,
`test_minimal_model_scalar_typing.py`,
`test_w_algebras_scalar_typing.py`,
`test_free_fields_scalar_typing.py`,
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Targeted propagation grep
finds no live source occurrence of the repaired Heisenberg
Chern--Weil FP forms or the repaired grand-unification
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) formula; historical
`final_gaps_20260413_213946` snapshots still contain old recorded
line excerpts.  A broader \(\mathrm{obs}_g=\kappa\lambda_g\) pattern
still has many live Vol I echoes outside this target and remains the
next propagation backlog.  Vol I `make fast` reaches the known stable
warning state after three passes with zero undefined citations and one
persistent undefined reference outside the edited lines.  Vol I
`make verify` exits 0 with all checks passed.

A286 moves the obstruction-class scalar-characteristic repair into the
active Vol I examples surface.  The old examples prose still contained
theorem-level identities of the form
\(\mathrm{obs}_g=\kappa\lambda_g\) or
\(\mathrm{obs}_1=\kappa\lambda_1\), even where the surrounding
passages had already distinguished level parameters, central charge,
and typed invariants.  The repair rewrites the generic genus expansion
as
\[
\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g
\]
on the uniform-weight scalar lane, and rewrites the Virasoro instance
as
\[
\mathrm{obs}_g(\mathrm{Vir}_c)
=\kappa(\mathrm{Vir}_c)\lambda_g.
\]
The free-fermion spin-structure remark now reads
\(\mathrm{obs}_g(\cF)=\kappa(\cF)\lambda_g=\frac14\lambda_g\),
and the affine Kac--Moody remark now reads
\[
\mathrm{obs}_g(\widehat{\fg}_k)
=\kappa(\widehat{\fg}_k)\lambda_g.
\]

The W-algebra surface needed the strongest correction.  Its genus-one
theorem/table now types
\[
\mathrm{obs}_1(\mathcal W_N^k)
=\kappa(\mathcal W_N^k)\lambda_1
\]
and separates this unconditional genus-one scalar statement from the
multi-weight genus-\(\geq2\) cross-channel corrections.  The W3 scalar
conductor no longer appears as \(\kappa+\kappa'\); it is written as
\[
\kappa(\mathcal W_3^k)+\kappa((\mathcal W_3^k)^!)=250/3,
\]
with the fixed-companion value at \(c=50\),
\[
\kappa(\mathcal W_3^k)=\kappa((\mathcal W_3^k)^!)=125/3.
\]
The same pass removes the remaining `\kappa = c` and `\kappa/c`
scalar-characteristic shorthands from `w_algebras.tex`, while leaving
local channel components and OPE parameters intact.

Vol I now has the additional guard
`compute/tests/test_kac_moody_scalar_typing.py`, with
`test_kac_moody_higher_genus_obstruction_uses_typed_characteristic()`.
The existing genus-expansion, free-field, and W-algebra scalar-typing
guards were strengthened to require the typed obstruction classes and
W3 conductor/fixed-point values, and to forbid the retired bare
obstruction and conductor formulas.

A286 verification: `python3 -m py_compile` passes on all ten Vol I
scalar-typing guards.  Focused pytest reports `10 passed` for
`test_grand_unification_scalar_typing.py`,
`test_genus_complete_scalar_typing.py`,
`test_kac_moody_scalar_typing.py`,
`test_minimal_model_scalar_typing.py`,
`test_w_algebras_scalar_typing.py`,
`test_free_fields_scalar_typing.py`,
`test_heisenberg_fp_scalar_typing.py`,
`test_y_algebras_scalar_typing.py`,
`test_landscape_census_scalar_typing.py`, and
`test_genus_expansions_scalar_typing.py`.  Exact grep over
`chapters/examples` finds no remaining live
\(\mathrm{obs}_g=\kappa\cdot\lambda_g\) or
\(\mathrm{obs}_g=\kappa\lambda_g\) fixed-string form; focused W-algebra
grep finds no remaining `\kappa + \kappa'`, `\kappa = c`, or
`\kappa/c` scalar shorthand.  Vol I `make fast` reaches the known
stable warning state after three passes with zero undefined citations
and one persistent undefined reference.  Vol I `make verify` exits 0
with all checks passed.  Broader non-examples Vol I surfaces still
contain bare obstruction-class echoes and remain the next propagation
backlog.

A287 moves the same obstruction-characteristic discipline into the
active Vol I theorem and introduction layer.  The target set was
`introduction.tex`, `part_ii_platonic_introduction.tex`,
`higher_genus_foundations.tex`, `higher_genus_modular_koszul.tex`,
`higher_genus_complementarity.tex`, and
`thqg_introduction_supplement_body.tex`.  These files still contained
the reader-facing forms
\(\mathrm{obs}_g=\kappa\lambda_g\),
\(\mathrm{obs}_1=\kappa\lambda_1\), or a typed coefficient
\(\kappa(\cA)\) with an untyped left-hand side
\(\mathrm{obs}_g\).  The repair separates the abstract obstruction
cochain from the scalar-characteristic theorem:
\[
\mathrm{obs}_g=\sum_{g_1+g_2=g}d_{g_1}d_{g_2}
\quad\hbox{is the obstruction cochain,}
\]
whereas the scalar theorem reads
\[
\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g
\]
on the proved scalar lane, with the genus-one specialization
\[
\mathrm{obs}_1(\cA)=\kappa(\cA)\lambda_1
\]
unconditional for all-weight modular Koszul algebras.

The introduction and Part II introduction now state Theorem D with
\(\mathrm{obs}_g(\cA)\), not a bare obstruction class.  Where the text
mentions free energies, the scalar lane is written as
\[
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
and the all-weight formula is separated as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA).
\]
`higher_genus_foundations.tex` now types the Heisenberg, affine
Kac--Moody, W3, Hodge-filtration, additivity, Chern--Simons benchmark,
and obstruction-lifting formulas.  The family-evaluated forms now read,
for example,
\[
\mathrm{obs}_g(V_k(\fg))
=\kappa(V_k(\fg))\lambda_g
=\frac{(k+h^\vee)\dim(\fg)}{2h^\vee}\lambda_g
\]
and
\[
\mathrm{obs}_1(\mathcal W_3^k)
=\kappa(\mathcal W_3^k)\lambda_1
=(5c/6)\lambda_1.
\]
`higher_genus_modular_koszul.tex` now carries the typed formula in
the MK axiom list, the proof-routing diagram, and the four-test
interface.

Vol I now has
`compute/tests/test_theorem_layer_obstruction_scalar_typing.py`, with
two guards: one requiring the typed theorem-layer obstruction formulas
in each touched file, and one forbidding the retired bare
\(\mathrm{obs}_g=\kappa\lambda_g\) and
\(\mathrm{obs}_1=\kappa\lambda_1\) slogans.

A287 verification: `python3 -m py_compile` passes on all scalar-typing
guards.  Focused pytest reports `12 passed` for the scalar-typing
guard suite, including the new theorem-layer guard.  Exact target-set
grep finds no remaining `\mathrm{obs}_g=\kappa`,
`\mathrm{obs}_g = \kappa`, `\mathrm{obs}_1=\kappa`, or
`\mathrm{obs}_1 = \kappa` scalar slogan in the six active
theorem-layer files.  The only remaining `\mathrm{obs}_g =` in that
target set is the abstract obstruction-cochain definition in
`higher_genus_foundations.tex`, not a scalar-characteristic formula.
Vol I `make fast` reaches the known stable warning state after three
passes with zero undefined citations and one persistent undefined
reference.  Vol I `make verify` exits 0 with all checks passed.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched Vol I files and the new guard.  The remaining backlog is now
localized: standalone documents, `concordance.tex`, active
example-only genus-one echoes, and the FP free-energy shorthand
population in `higher_genus_modular_koszul.tex`.

A288 removes the remaining active example-only genus-one obstruction
slogans.  The target set was
`w_algebras_deep.tex`, `y_algebras.tex`,
`w3_holographic_datum.tex`, and `n2_superconformal.tex`.  Each file
already named the relevant typed scalar characteristic, but still
used the unbased genus-one statement
\(\mathrm{obs}_1=\kappa\lambda_1\) or
\(\mathrm{obs}_1=\kappa\cdot\lambda_1\).

The W-rank comparison now reads
\[
\mathrm{obs}_1(\Walg_N)=\kappa(\Walg_N)\lambda_1.
\]
The \(Y_{1,1,1}\) computation now reads
\[
\mathrm{obs}_1(Y_{1,1,1}[\Psi])
=\kappa(Y_{1,1,1}[\Psi])\lambda_1
=\Psi\lambda_1.
\]
The W3 holographic datum now uses
\[
\mathrm{obs}_1(\Walg_3)=\kappa(\Walg_3)\lambda_1
\]
both in its convention paragraph and in the uniform-weight/cross-channel
remark.  The \(N=2\) SCA shadow-tower remark now reads
\[
\mathrm{obs}_1(\mathrm{SCA}_c)
=\kappa(\mathrm{SCA}_c)\lambda_1.
\]

Vol I now has
`compute/tests/test_example_genus_one_obstruction_scalar_typing.py`,
with guards requiring these four typed formulas and forbidding the
retired bare genus-one obstruction slogans in the active example
target set.

A288 verification: `python3 -m py_compile` passes on all scalar-typing
guards.  Focused pytest reports `14 passed` for the scalar-typing
guard suite, including the new active-example guard.  Exact grep over
active `chapters/examples` finds no remaining
`\mathrm{obs}_1 = \kappa` or `\mathrm{obs}_1=\kappa`
fixed-string occurrence.  Vol I `make fast` reaches the known stable
warning state after three passes with zero undefined citations and one
persistent undefined reference.  Vol I `make verify` exits 0 with all
checks passed.  Scoped `git diff --check` and trailing-whitespace
scans pass on the touched Vol I files and the new guard.  Remaining
stale genus-one copies are now restricted to standalone documents and
`concordance.tex`; broader \(\mathrm{obs}_g\) and FP free-energy
shorthands remain in standalone, concordance, type-system, and selected
active theory-summary surfaces.

A289 removes the next active theory/appendix summary layer of
obstruction slogans.  The target set was `appendices/type_system.tex`,
`master_reconstruction.tex`, `theorem_A_infinity_2.tex`,
`chiral_climax_platonic.tex`, and
`mc5_class_m_chain_level_platonic.tex`.  These active inputs still
contained compact forms such as
\(\mathrm{obs}_g=\kappa\lambda_g\) or
\(\mathrm{obs}_g=\kappa\cdot\lambda_g\), and one climax summary still
wrote the all-weight package as
\(F_g=\kappa\lambda_g^{\mathrm{FP}}+\delta F_g^{\mathrm{cross}}\).

The type-system appendix now writes
\[
\mathrm{obs}_g(\Ab)=\kappa(\Ab)\lambda_g.
\]
The master reconstruction and Theorem A bridge now write
\[
\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g.
\]
The climax theorem now separates the all-weight free energy from its
scalar summand:
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
with the exact obstruction identity
\(\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g\) only on the
uniform-weight surface.  The K3 MC5 paragraph now types the \(g=2\)
formality obstruction as
\[
\mathrm{obs}_2(\mathbf H_{\Delta_5})
=\kappa_{\mathrm{ch}}(\mathbf H_{\Delta_5})\lambda_2
\]
and its Heegner restriction as
\[
\mathrm{ob}_{\mathrm{form}}|_{H_n}
=c_n=\kappa_{\mathrm{ch}}(\mathbf H_{\Delta_5})\lambda_2|_{H_n}.
\]

Vol I now has
`compute/tests/test_active_theory_summary_obstruction_scalar_typing.py`,
with guards requiring the five active summary repairs and forbidding
bare \(\mathrm{obs}_g=\kappa\lambda_g\) and
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) forms in the target set.

A289 verification: `python3 -m py_compile` passes on all scalar-typing
guards.  Focused pytest reports `16 passed` for the scalar-typing
guard suite, including the new active-theory-summary guard.  Exact
target-set grep finds no remaining `\mathrm{obs}_g=\kappa`,
`\mathrm{obs}_g = \kappa`, `F_g=\kappa`, or `F_g = \kappa`
occurrence in the five active summary files.  Vol I `make fast`
reaches the known stable warning state after three passes with zero
undefined citations and one persistent undefined reference.  Vol I
`make verify` exits 0 with all checks passed.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I files and the new guard.  Remaining stale copies are now
concentrated in standalone documents, `concordance.tex`, audit/cache
prose, and the larger active FP free-energy shorthand backlog in
`higher_genus_modular_koszul.tex`.

A290 removes the active FP free-energy shorthand backlog in
`higher_genus_modular_koszul.tex`.  The chapter still had theorem-level
forms such as
\[
F_g=\kappa\lambda_g^{\mathrm{FP}},\qquad
F_g=\kappa\cdot\lambda_g^{\mathrm{FP}},
\]
as well as UW/shadow/tropical and genus-specific variants like
\(F_2=\kappa\cdot7/5760\) and
\(F_3=\kappa\cdot31/967680\).  These forms erased the distinction
between the scalar summand and the full all-weight free energy.

The active chapter now writes the scalar lane as
\[
F_g^{\mathrm{sc}}(\cA)
=\kappa(\cA)\lambda_g^{\mathrm{FP}}
\]
and the full all-weight decomposition as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA).
\]
Family-specific passages now type the source by
\(\kappa(\cW_3)\), \(\kappa(\mathrm{Vir}_c)\),
\(\kappa(\cH_\kappa)\), \(\kappa(V_k(\mathfrak{sl}_2))\), or
\(\kappa(\mathcal A)\).  Cross-channel ratio denominators now use
\(F_g^{\mathrm{sc}}\), and the Riccati/Polyakov paragraph keeps the
local metric coordinate \(\kappa\) while avoiding an untyped
free-energy theorem.

Vol I now has
`compute/tests/test_higher_genus_modular_koszul_fp_scalar_typing.py`,
which requires the typed scalar-source forms and forbids the retired
bare FP free-energy, UW/shadow/tropical, ratio, and genus-specific
variants in the active chapter.

A290 verification: the focused new guard reports `2 passed`; the full
scalar-typing guard suite reports `18 passed`; `python3 -m py_compile`
passes on all scalar-typing guards.  Exact fixed-string and regex
greps over `higher_genus_modular_koszul.tex` find no remaining
`F_g = \kappa \cdot \lambda_g^{\mathrm{FP}}`,
`F_g=\kappa\cdot\lambda_g^{\mathrm{FP}}`, `F_g = \kappa`,
`F_g=\kappa`, or bare
\(\kappa\cdot\lambda_g^{\mathrm{FP}}\) FP-free-energy occurrence.
Vol I `make fast` reaches the known stable warning state after three
passes with zero undefined citations and one persistent undefined
reference.  Vol I `make verify` exits 0 with all checks passed.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched Vol I chapter, the new guard, and these Vol II ledgers.

A291 moves the same scalar/free-energy discipline into the active
Vol II-native Theorem D surface
`chapters/theory/theorems_C_D_native_vol2_platonic.tex`.  The chapter
still presented the uniform-weight coefficient as
\[
F_g=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}
\]
and described the Vol I theorem as if this were the all-weight tensor
object.  It also identified off-diagonal tensor entries
\(K_{w_1,w_2}\) directly with
\(\delta F_g^{\mathrm{cross}}\), collapsing a tensor-valued two-form
with its Faber--Pandharipande numerical trace.

The native theorem surface now separates the three objects:
\[
K(\cA)\in
\operatorname{Sym}^2(\cF^\vee)\otimes\Omega^2(\overline{\cM}_{g,n}),
\qquad
F_g^{\mathrm{sc}}(\cA)
=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\]
and
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA).
\]
The off-diagonal \(K_{w_1,w_2}\) remains a tensorial \(2\)-form; its
Faber--Pandharipande numerical trace contributes
\(\delta F_g^{\mathrm{cross}}|_{w_1,w_2}\).  The bare curvature
slogan \(d^2=\kappa\omega_g\) is now typed as
\[
d^2=\kappaChHodge(\mathcal A)\omega_g.
\]

The existing Vol II guard
`compute/tests/test_wN_tensor_arakelov_weight_distribution.py` now
contains
`test_tensor_arakelov_source_separates_scalar_trace_from_full_energy`,
which requires the scalar/free-energy split and the trace wording, and
forbids the retired \(F_g=\kappaChHodge\lambda_g^{\mathrm{FP}}\),
\(d^2=\kappa\omega_g\), and
\(K_{w_1,w_2}=\delta F_g^{\mathrm{cross}}\) collapses.

A291 verification: focused pytest for
`test_wN_tensor_arakelov_weight_distribution.py` reports `7 passed`;
the existing deletion-ledger reader for the same file reports
`1 passed`; `python3 -m py_compile` passes on the modified guard.
Targeted regex grep over `theorems_C_D_native_vol2_platonic.tex`
finds no surviving `F_g = \kappaChHodge`, `F_g=\kappaChHodge`,
`F_g(\cA)=\kappaChHodge`, bare
\(\kappaChHodge(\cA)\cdot\lambda_g^{\mathrm{FP}}\),
\(d^2=\kappa\cdot\omega_g\), or direct
\(K_{w_1,w_2}=\delta F_g^{\mathrm{cross}}\) occurrence.  Vol II
`make fast` converges after two passes with zero undefined citations
and zero undefined references.  Vol II `make verify-licensing` exits
0 with zero blocking violations and zero warnings.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
theorem, guard, and ledgers.

A292 carries the same correction into the active foundations chapter.
Two passages in `chapters/theory/foundations.tex` still summarized the
genus expansion as
\[
F_g=\kappaChHodge(\cA)\cdot\lambda_g^{\mathrm{FP}}.
\]
In the bar-complex overview this made the scalar FP coefficient look
like the full theorem-level genus expansion.  In the open/closed MC
consistency paragraph it made the conditional higher-genus scalar
continuation look like the all-weight free energy.

The foundations overview now names the scalar coefficient and the
all-weight decomposition separately:
\[
F_g^{\mathrm{sc}}(\cA)
=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\qquad
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA).
\]
The open/closed MC paragraph now writes
\[
F_1^{\mathrm{sc}}(\cA)=\kappaChHodge(\cA)/24
\]
at genus one, keeps the higher-genus scalar continuation conditional,
and adds \(\delta F_g^{\mathrm{cross}}(\cA)\) for the all-weight
numerical free energy.

Vol II now has
`compute/tests/test_foundations_theorem_d_scalar_trace_typing.py`,
which requires these foundations summaries and forbids the retired
bare \(F_g=\kappaChHodge\lambda_g^{\mathrm{FP}}\),
\(F_g(\cA)=\kappaChHodge\lambda_g^{\mathrm{FP}}\), and
\(F_1=\kappaChHodge(\cA)/24\) forms in the active file.

A292 verification: the focused new foundations guard reports
`2 passed`; the paired foundations plus tensor-Arakelov guard run
reports `9 passed`; `python3 -m py_compile` passes on both guards.
Targeted regex grep over `foundations.tex` finds no surviving
`F_g = \kappaChHodge`, `F_g=\kappaChHodge`,
`F_g(\cA)=\kappaChHodge`, `F_1 = \kappaChHodge(\cA)/24`, or bare
\(\kappaChHodge(\cA)\cdot\lambda_g^{\mathrm{FP}}\) occurrence.  Vol
II `make fast` converges after two passes with zero undefined
citations and zero undefined references.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans
pass on the touched foundations file, the new guard, and these ledgers.

A293 carries the correction into the active BV/BRST connection chapter
`chapters/connections/bv_brst.tex`.  Three live statements still wrote
the general scalar genus tower as
\[
F_g=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}
\]
or as the class-\(\mathsf C\) pattern
\[
F_g=\kappaChHodge(\cA)\cdot\lambda_g^{\mathrm{FP}}.
\]
This notation made the scalar FP trace look like the full all-weight
free energy, although the same paragraph already identifies
\(\delta F_g^{\mathrm{cross}}\) as the multi-weight obstruction to
chain-level BV/bar identification.

The active BV/BRST chapter now writes the general scalar lane as
\[
F_g^{\mathrm{sc}}(\cA)
=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\]
keeps the Heisenberg all-genera equality as a family-specific scalar
determinant-line statement, and states the non-CY predictions as a
class-\(\mathsf C\) scalar pattern rather than a total free-energy
formula.

Vol II now has
`compute/tests/test_bv_brst_theorem_d_scalar_trace_typing.py`, which
requires the BV/BRST scalar \(F_g^{\mathrm{sc}}\) forms and the
explicit multi-weight correction sentence, while forbidding the retired
bare \(F_g=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}\) and
\(F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}\) forms in the
active file.

A293 verification: the focused new BV/BRST guard reports `2 passed`;
the paired BV/BRST, foundations, and tensor-Arakelov guard run reports
`11 passed`; `python3 -m py_compile` passes on all three guards.
Targeted regex grep over `bv_brst.tex` finds no surviving
`F_g = \kappaChHodge`, `F_g=\kappaChHodge`,
`F_g(\cA)=\kappaChHodge`, or bare
\(\kappaChHodge(\cA)\cdot\lambda_g^{\mathrm{FP}}\) occurrence.  Vol
II `make fast` converges after two passes with zero undefined
citations and zero undefined references.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans
pass on the touched BV/BRST file, the new guard, and these ledgers.

A294 propagates the same scalar/full-energy separation across the live
Vol II chapter surface.  The stale compact formulas
\[
F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}},
\qquad
F_g=\kappaChHodge(\mathrm{Vir}_c)\lambda_g^{\mathrm{FP}},
\]
still occurred in holomorphic-topological comparison material, the
THQG perturbative/fredholm/gravity/s-duality chapters, the HT
bulk-boundary line copies, the rosetta atlas, the \(W\)-algebra stable
examples, and the recast foundations draft.  Those forms erased the
distinction between the protected scalar FP trace and the full
all-weight free energy.

The repaired surfaces now use
\[
F_g^{\mathrm{sc}}(\cA)
=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}
\]
for scalar-lane statements and
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
 +\delta F_g^{\mathrm{cross}}(\cA)
\]
for the all-weight numerical free energy.  Uniform-weight equalities
are kept only after saying that no mixed-weight channel contributes.
The THQG multi-weight prose now says that the full free energy receives
the cross-channel correction; it no longer says that the scalar formula
itself receives or fails by that correction.

Vol II now has
`compute/tests/test_scalar_trace_full_energy_surface_guard.py`, which
scans `chapters/connections`, `chapters/theory`, and
`chapters/examples` for the retired bare \(F_g=\kappaChHodge\) and
\(F_g(\cA)=\kappaChHodge\) patterns and requires repaired scalar/full
decomposition anchors.  The BV/BRST guard now also catches the spaced
\(F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}\) form.

A294 verification: the scalar/full-energy guard suite reports
`17 passed`; `python3 -m py_compile` passes on all five scalar/tensor
guards.  Fixed-string greps over `chapters/connections`,
`chapters/theory`, and `chapters/examples` find no surviving
`F_g(\cA)=\kappaChHodge`, `F_g(\cA) = \kappaChHodge`,
`F_g=\kappaChHodge`, `F_g = \kappaChHodge`, or `scalar formula
receives` occurrence.  Vol II `make fast` converges after two passes
with zero undefined citations and zero undefined references.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and zero
warnings.

A295 carries the same invariant into the executable genus-2 obstruction
engine.  The engine already computed the correct coefficient
\(\kappa\lambda_g^{\mathrm{FP}}\), but its public documentation and
formula metadata still called that value the genus-\(g\) free energy.
That made the compute oracle advertise the scalar FP trace as the full
all-weight numerical energy.

The engine now exposes
`F_g_scalar_free_energy(kappa, g)` as the canonical scalar trace and
keeps `F_g_free_energy(kappa, g)` only as a documented compatibility
alias.  It also adds
`full_free_energy_decomposition(kappa, g, delta_cross)`, returning
\[
F_g^{\mathrm{sc}},\qquad
\delta F_g^{\mathrm{cross}},\qquad
F_g=F_g^{\mathrm{sc}}+\delta F_g^{\mathrm{cross}}.
\]
The module docstring, graph-amplitude descriptions, family checks,
ratio statements, and genus-2 complementarity comments now state
scalar traces rather than unqualified free energies.

`compute/tests/test_genus2_obstruction_engine.py` now asserts that the
legacy function is a scalar-trace alias and that a nonzero
\(\delta F_g^{\mathrm{cross}}\) changes the full free energy without
changing \(F_g^{\mathrm{sc}}\).

A295 verification: focused pytest for
`test_genus2_obstruction_engine.py` reports `58 passed`; the paired
scalar/tensor manuscript+compute guard run reports `75 passed`;
`python3 -m py_compile` passes on the touched engine and all
scalar/tensor guards.  Target-file grep finds no surviving
`Genus-g free energy`, `F_g(A) = kappa`, `F_g = kappa`, `scalar formula
receives`, or `free energy F_g =` form in the engine and its tests.
Vol II `make fast` converges after two passes with zero undefined
citations and zero undefined references.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and zero
warnings.  The cross-volume scan still exposes live Vol I/Vol III
manuscript debt for the next propagation pass.

A296 propagates the scalar/full-coefficient discipline into the live
Vol III manuscript.  The remaining Vol III debt was not only the exact
"scalar formula receives/acquires" wording: active examples also wrote
\(F_g=\kappa_{\mathrm{ch}}\lambda_g^{\mathrm{FP}}\) and then appended
cross-channel corrections, which made the scalar trace and the full
all-weight genus coefficient share the same symbol.

The repaired Vol III surface now writes
\[
F_g^{\mathrm{sc}}
 = \kappa_{\mathrm{ch}}\lambda_g^{\mathrm{FP}}
\]
or the categorical analogue on the uniform-weight lane, and reserves
\[
F_g=F_g^{\mathrm{sc}}+\delta F_g^{\mathrm{cross}}
\]
for the full coefficient.  The repair touched the modular Koszul
bridge, braided factorisation, modular trace, CY-to-chiral DT
comparison, K3 \(\times\) E programme, K3--BKM comparison, and
bar-cobar bridge surfaces.  Class-G and single-weight loci now state
explicitly why the scalar trace equals the full coefficient there.

Vol III now has
`compute/tests/test_scalar_trace_full_energy_surface.py`, which scans
the seven live source files for the retired scalar-correction slogans
and requires repaired scalar/full-coefficient anchors.

A296 verification: the new Vol III guard plus the modular bridge suite
reports `55 passed`; `python3 -m py_compile` passes on the new guard.
Fixed-string scans find no surviving `scalar formula receives`, `scalar
formula acquires`, `scalar formula fails and requires`, or `all-weight,
with cross-channel correction` occurrence in Vol III `chapters` and
`compute`, and no surviving TeX `F_g = \kappa_{\mathrm{ch}}` formula.
Vol III `make fast` exits 0 after four passes with zero undefined
citations, one stable undefined-reference class, and no rerun request.
The unresolved labels
`sec:bar-cobar-bi-based-k3e` and
`def:bar-cobar-bi-based-ran-datum` are residual structural label debt in
the dirty Vol III tree, not introduced by this pass.  Scoped
`git diff --check` and trailing-whitespace scans pass on all touched
files.

A297 closes that residual Vol III label debt.  The labels themselves
were present in `chapters/connections/bar_cobar_bridge.tex`, but they
were absent from the final aux file because an earlier malformed
`\textup` parenthetical,
`conifold\textup{), the toric`, left TeX scanning the argument until
the input file ended.  The result was a structurally false build
surface: source-level labels existed, while the compiled monograph could
not resolve references to the bi-based Ran-space section or datum.

The repair changes the parenthetical to
`conifold\textup{)}, the toric`, so the bi-based section and the
definition are actually processed by LaTeX.  Vol III now has
`compute/tests/test_bar_cobar_bridge_label_surface.py`, which requires
the two source labels and forbids the malformed local-CY parenthetical.

A297 verification: focused pytest for the new guard reports `2
passed`; paired with the scalar/full-energy guard it reports `7
passed`; `python3 -m py_compile` passes on both guards.  Vol III
`make fast` converges after three passes with zero undefined citations,
zero undefined references, and no rerun request.  The final `main.aux`
contains both `\newlabel{sec:bar-cobar-bi-based-k3e}` and
`\newlabel{def:bar-cobar-bi-based-ran-datum}`.  Targeted log/source
scans find no surviving runaway `\textup` signature or undefined
reference warning on the repaired surface.  Scoped `git diff --check`
and trailing-whitespace scans pass on the touched files.

A298 propagates the scalar/full-coefficient discipline into the live
Vol I manuscript and compute surface.  The remaining Vol I debt
included exact retired wording that the scalar formula receives a
cross-channel correction, plus theorem-summary surfaces where the
Faber--Pandharipande trace still used bare \(F_g\).

The repaired Vol I surface now reserves
\[
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
\]
for the scalar trace and writes
\[
F_g=F_g^{\mathrm{sc}}+\delta F_g^{\mathrm{cross}}
\]
for the full all-weight coefficient in the multi-weight genus
\(g\ge2\) regime.  The repair touched the N=2 SCA example,
concordance, \(E_1\) modular Koszul chapter, higher-genus
complementarity, higher-genus foundations, higher-genus modular
Koszul chapter, and the CSFT action docstring.

Vol I now has
`compute/tests/test_scalar_full_coefficient_typing.py`, which requires
the scalar/full anchors across those touched surfaces and forbids the
retired scalar-correction slogans.

A298 verification: `python3 -m py_compile` passes on
`compute/lib/csft_from_bar.py` and the new guard.  Focused pytest for
the new guard reports `3 passed`; the scalar-typing guard suite reports
`21 passed`.  Fixed-string scans find no surviving `scalar formula
receives`, `scalar formula acquires`, `scalar formula fails and
requires`, or `all-weight, with cross-channel correction` occurrence in
Vol I `chapters` and `compute`; the targeted touched-file regex scan
finds no surviving bare full-coefficient
`F_g=\kappa\lambda_g^{\mathrm{FP}}` surface.  Vol I `make fast`
converges after three passes with zero undefined citations, one stable
undefined reference, and no rerun request.  Vol I `make verify` exits
0 with `verify_edit: all checks passed`.  Scoped `git diff --check`
and trailing-whitespace scans pass on the touched Vol I files.

A299 closes the flexible/displayed Vol I variants left after A298.
The strengthened guard exposed theorem and summary surfaces that no
longer used the exact retired phrase but still wrote displayed
decompositions in the collapsed form
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
       +\delta F_g^{\mathrm{cross}}(\cA),
\]
or used bare \(F_g(\cA)\) for the scalar \(\widehat A\)-series.

The repaired surfaces now write the full coefficient as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
       +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
and write scalar generating functions as
\[
\sum_g F_g^{\mathrm{sc}}(\cA)x^{2g}.
\]
The pass also removed the misleading `ALL-WEIGHT + ...` parenthetical
labels from the N=2 SCA and higher-genus theorem surfaces.

A299 strengthened
`compute/tests/test_scalar_full_coefficient_typing.py` to include the
preface and to use regex checks for spaced/displayed collapsed
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) variants.

A299 verification: `python3 -m py_compile` passes on
`compute/lib/csft_from_bar.py` and the strengthened guard.  The focused
guard reports `3 passed`; the scalar-typing guard suite reports
`21 passed`.  Targeted scans over the guarded Vol I files find no
surviving `ALL-WEIGHT +`, `ALL-WEIGHT;`, `free energy receives a
cross-channel`, `scalar formula receives`, or
`all-weight, with cross-channel correction` occurrence, and no surviving
bare full-coefficient \(F_g=\kappa\lambda_g^{\mathrm{FP}}\) variant.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I files.

A300 propagates the scalar/full-coefficient discipline into the Vol I
summary layer: `guide_to_main_results.tex`, `outlook.tex`,
`master_concordance.tex`, `feynman_diagrams.tex`, and
`genus_complete.tex`.  These surfaces still let the scalar
Faber--Pandharipande trace appear as the full all-weight coefficient,
usually through formulas of the form
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
       +\delta F_g^{\mathrm{cross}}(\cA).
\]

The repaired summary surface now writes the all-weight coefficient as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
       +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
or, in the integrable-hierarchy discussion where the diagonal channel
is the actual object,
\[
F_g(\cA)=F_g^{\mathrm{diag}}(\cA)
       +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{diag}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
The uniform-weight obstruction statement
\(\operatorname{obs}_g(\cA)=\kappa(\cA)\lambda_g\) remains a scalar-lane
statement and was not promoted to an all-weight formula.

A300 adds
`compute/tests/test_scalar_full_summary_surfaces.py`, which requires
the repaired scalar/full or diagonal/full anchors across the five Vol I
summary files and forbids the retired summary-level collapse.

A300 verification: `python3 -m py_compile` passes on the new guard.
The focused summary guard reports `3 passed`; paired with
`compute/tests/test_scalar_full_coefficient_typing.py` it reports
`6 passed`.  Targeted scans over the A300 cluster find no surviving
bare full-coefficient \(F_g=\kappa\lambda_g^{\mathrm{FP}}\) variant and
no surviving retired all-weight/scalar-correction wording.  Vol I
`make fast` converges after three passes with zero undefined citations,
one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I files and guard.

A301 repairs the live Vol I `genus_expansions.tex` chapter.  The
chapter still exported the scalar Faber--Pandharipande contraction as
bare \(F_g\) in its opening summary, tautological integral, Verlinde
comparison, universal generating-function theorem, complementarity
discussion, master landscape table introduction, growth-rate argument,
and Chern--Simons loop-expansion interpretation.  It also sometimes
said the scalar formula "fails" or "receives" a correction, rather than
separating the scalar/diagonal coefficient from the full all-weight
coefficient.

The repaired chapter now writes the scalar lane as
\[
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
the \(\mathcal W_3\) diagonal lane as
\[
F_g^{\mathrm{diag}}(\mathcal W_3^k)
  =\kappa(\mathcal W_3^k)\lambda_g^{\mathrm{FP}}
  ={5c(k)\over 6}\lambda_g^{\mathrm{FP}},
\]
and the full multi-weight coefficient as
\[
F_g=F_g^{\mathrm{sc}}+\delta F_g^{\mathrm{cross}}
\quad\hbox{or}\quad
F_g=F_g^{\mathrm{diag}}+\delta F_g^{\mathrm{cross}},
\]
according to the local lane.
The universal generating function, universal ratios, growth-rate
statement, and loop-expansion comparison now use
\(F_g^{\mathrm{sc}}\) when they are using the
Faber--Pandharipande scalar trace.

A301 replaces Vol I
`compute/tests/test_genus_expansions_scalar_typing.py` with a guard
that requires the repaired scalar/full and diagonal/full anchors in
`genus_expansions.tex` and forbids the old bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) forms for \(\cA\),
\(\mathcal A\), \(\hat{\mathfrak g}_k\), and generic \(\kappa\).

A301 verification: `python3 -m py_compile` passes on the repaired
guard.  The focused guard reports `3 passed`; paired with
`compute/tests/test_scalar_full_coefficient_typing.py` and
`compute/tests/test_scalar_full_summary_surfaces.py` it reports
`9 passed`; the broader scalar-typing guard set reports `26 passed`.
Targeted scans over `genus_expansions.tex` find no surviving bare
full-coefficient Faber--Pandharipande formula and no retired
scalar/full phrasing in that chapter.  Vol I `make fast` converges
after three passes with zero undefined citations, one stable undefined
reference, and no rerun request.  Vol I `make verify` exits 0 with
`verify_edit: all checks passed`.  Scoped `git diff --check` and
trailing-whitespace scans pass on the touched Vol I file and guard.

A302 repairs the live Vol I Heisenberg exact Gaussian lane in
`heisenberg_eisenstein.tex`.  The Heisenberg chapter still wrote the
Faber--Pandharipande tower as bare
\[
F_g(\mathcal H_\kappa)
  =\kappa(\mathcal H_\kappa)\lambda_g^{\mathrm{FP}}
\]
at its summary table, scalar-genus introduction, complementarity
paragraph, genus-2 computation, Faber--Pandharipande theorem,
shadow-calculus Bernoulli tower, and closed-sector compatibility
proof.  On this lane the equality with full \(F_g\) is true, but only
because the Heisenberg is exactly Gaussian and the mixed-channel
summand vanishes.

The repaired chapter now writes
\[
F_g(\mathcal H_\kappa)
  =F_g^{\mathrm{sc}}(\mathcal H_\kappa)
  =\kappa(\mathcal H_\kappa)\lambda_g^{\mathrm{FP}},
\qquad
\delta F_g^{\mathrm{cross}}(\mathcal H_\kappa)=0.
\]
The same scalar/full equality is propagated through the dual
Heisenberg contribution, genus-2 ratio, \(\kappa=1\) table caption,
Bernoulli tower, \(\hat A\)-generating function, and closed-sector
projection.

A302 replaces Vol I
`compute/tests/test_heisenberg_fp_scalar_typing.py` with a guard that
requires the exact Gaussian anchors \(F_g=F_g^{\mathrm{sc}}\), the
vanishing \(\delta F_g^{\mathrm{cross}}=0\), the dual contribution,
the scalar integral theorem, the \(\hat A\)-series equality, and the
closed-sector compatibility equality.  The guard forbids the old bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) and bare
\(F_2=\kappa(\mathcal H_\kappa)\cdots\) forms in the chapter.

A302 verification: `python3 -m py_compile` passes on the repaired
guard.  The focused Heisenberg guard reports `3 passed`; paired with
the A300--A301 scalar/full guards it reports `12 passed`; the broader
scalar-typing guard set reports `28 passed`.  Targeted scans over
`heisenberg_eisenstein.tex` find no surviving bare full-coefficient
Faber--Pandharipande formula in the Heisenberg lane.  Vol I
`make fast` converges after three passes with zero undefined citations,
one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A303 extends the exact-lane scalar/full typing to the live Vol I
free-field chapter `free_fields.tex`.  The free fermion genus theorem,
proof, Heisenberg comparison, Pfaffian/spin discussion, genus table,
and \(\hat A\)-universality remark still exported the
Faber--Pandharipande tower as bare
\[
F_g(\mathcal F)=\kappa(\mathcal F)\lambda_g^{\mathrm{FP}}
\quad\hbox{or}\quad
F_g(\mathcal F)={1\over4}\lambda_g^{\mathrm{FP}}.
\]
That formula is true on the class~G lane only after naming the reason:
the free fermion has no mixed-channel genus correction on this shadow
surface.

The repaired chapter now writes
\[
F_g(\mathcal F)
  =F_g^{\mathrm{sc}}(\mathcal F)
  =\kappa(\mathcal F)\lambda_g^{\mathrm{FP}},
\qquad
\delta F_g^{\mathrm{cross}}(\mathcal F)=0.
\]
The Heisenberg ratio and \(\hat A\)-generating functions are typed with
\(F_g^{\mathrm{sc}}\).  The Pfaffian and spin-structure data remain in
the partition function \(Z_g\), not in this scalar/full class~G
coefficient.

A303 replaces Vol I
`compute/tests/test_free_fields_scalar_typing.py` with a guard that
requires \(F_g=F_g^{\mathrm{sc}}\), the vanishing
\(\delta F_g^{\mathrm{cross}}=0\), the theorem generating-function
equality, the \(\cF\) spin/Pfaffian equality, the \(\hat A\)-series
equality, and the first typed value.  The guard forbids bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\),
\(F_g(\cF)=\kappa(\cF)\lambda_g^{\mathrm{FP}}\), and
\(F_g(\cF)=\frac14\cdots\) variants.

A303 verification: `python3 -m py_compile` passes on the repaired
guard.  The focused free-field guard reports `3 passed`; paired with
the A300--A302 scalar/full guards it reports `15 passed`; the broader
scalar-typing guard set reports `30 passed`.  Targeted scans over
`free_fields.tex` find no surviving bare full-coefficient
Faber--Pandharipande formula in the free-field class~G lane.  Vol I
`make fast` converges after three passes with zero undefined citations,
one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A304 repairs the live Vol I `landscape_census.tex` export surface for
the scalar/full genus-coefficient distinction.  The census lane column,
genus-one anomaly remark, genus-obstruction introduction, and
Faber--Pandharipande table/caption still stated the scalar
Faber--Pandharipande term as the unqualified full coefficient
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
For multi-weight rows this is only the scalar term.  The full
coefficient is
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]

The repaired census now states \(F_g=F_g^{\mathrm{sc}}\) only on
scalar/exact rows, writes the multi-weight decomposition explicitly,
changes the universal ratio to
\(F_2^{\mathrm{sc}}/F_1^{\mathrm{sc}}=7/240\), and retags the
Faber--Pandharipande table as scalar.  The genus-one remark now writes
the all-genus FP formula as \(F_g^{\mathrm{sc}}\) and reserves the
unqualified equality for \(F_1=F_1^{\mathrm{sc}}=\kappa/24\).

A304 adds Vol I
`compute/tests/test_landscape_census_scalar_typing.py` with positive
anchors for the scalar coefficient, exact-lane equality, multi-weight
decomposition, genus-one equality, scalar ratio, obstruction class, and
W-algebra anomaly-ratio formulas.  The guard forbids the old bare
\(F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}\), bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\), scalar-formula-fails wording,
and untyped FP table caption.

A304 verification: `python3 -m py_compile` passes on the guard.  The
focused landscape-census guard reports `3 passed`; paired with the
A300--A303 scalar/full guards it reports `18 passed`; the broader
scalar-typing guard set reports `32 passed`.  Targeted stale-form
scans over `landscape_census.tex` and the guard are clean.  The broad
Vol I scalar/full scan no longer reports `landscape_census.tex`; the
remaining debt is in `bar_construction.tex`,
`koszulness_vii_multiweight_platonic.tex`,
`frontier_modular_holography_platonic.tex`,
`entanglement_modular_koszul.tex`, `fourier_seed.tex`,
`arithmetic_shadows.tex`, and `shadow_L_function_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A319 strengthens item 8 of `Research Paper Refinement vol2.pdf`: the
\(E_3\)-PBW route now requires not only a filtered associated-graded
theorem, growth, amplitude, and page support, but also a Rees-flat
no-hidden-extension lift.  The conjecture header in
`chapters/theory/chiral_higher_deligne.tex` now carries
`\ClaimStatusConjectured` and licensing tags \(\gamma+\delta\).  The
package requires a complete, separated, exhaustive arity/collision
filtration with filtered \(E_3\)-operations, an associated-graded
isomorphism
\[
\theta_x\colon
\Free_{E_3}(H^\bullet(W_x[-2]))
\xrightarrow{\sim}\operatorname{gr}_F R_x,
\qquad
\ker\theta_x=\operatorname{coker}\theta_x=0,
\]
and the Rees condition
\[
\mathcal R_F R_x/(u)\cong\operatorname{gr}_F R_x,\qquad
\mathcal R_F R_x[u^{-1}]_{u=1}\simeq R_x,\qquad
H^\bullet(\mathcal R_F R_x)\ \text{is }u\text{-torsion-free}.
\]
The PBW spectral sequence is now written through the associated graded:
\[
E_1^{p,q}
=H^{p+q}(\operatorname{gr}^p_F R_x)
\cong H^{p+q}(\Free^p_{E_3}(W_x[-2]))
\Longrightarrow H^{p+q}(R_x).
\]

A319 also repairs the formal consequence proof: the differential
\(d_r\) has bidegree \((r,1-r)\) and therefore total degree \(+1\), not
total degree \(0\).  Vanishing in total degree \(>2\) persists because
each later \(E_{r+1}^{p,q}\) is a subquotient of the already-zero group
\(E_r^{p,q}\).  Polynomial growth remains only a finiteness/convergence
input; it does not kill cohomological degrees by itself.

A319 propagates the Rees-flat requirement to
`compute/lib/deletion_ledger.py` and `compute/tests/test_deletion_ledger.py`.
The strengthened Chiral Higher Deligne guard now requires the
conjecture status/licensing tag, Rees equation, \(u\)-torsion-free
condition, \(\theta_x\), kernel/cokernel vanishing, the
\(\operatorname{gr}^p_F R_x\) \(E_1\)-page formula, total-degree \(+1\)
differential, and the absence of the retired "preserves total degree"
sentence.  During the widened deletion-ledger run, A319 also repaired a
local scalar-lane drift in
`chapters/connections/thqg_perturbative_finiteness.tex`: the closed
generating function and adjacent remark now use \(F_g^{\mathrm{sc}}\)
rather than bare \(F_g\) on the proved scalar lane.

A319 verification: `python3 -m py_compile` passes on the touched
guards.  The paired Chiral Higher Deligne and deletion-ledger test set
reports `63 passed`; the widened Chiral Higher Deligne /
Boardman--Vogt / ordered-bar / deletion-ledger set reports `119
passed`.  Vol II `make fast` converges after two passes in
`/tmp/mkd-chiral-bar-cobar-vol2-20260616000232-70765`, producing
`2479pp` with zero undefined citations, zero undefined references, and
zero rerun requests.  Direct final-log scans find no fatal TeX errors,
undefined controls, LaTeX error headers, or unresolved
reference/citation warnings; the visible tail warnings are pdfTeX
named-destination warnings.  Vol II `make verify-licensing` exits 0
with zero blocking violations and 127 warning-class missing-tag lines
outside the touched surface.  Scoped `git diff --check` passes on the
touched files.

A318 repairs item 7 of `Research Paper Refinement vol2.pdf`: the
two-coloured Chiral Higher Deligne lift now names the signed
Boardman--Vogt tree homotopy rather than appealing to unnamed higher
coherence.  The main theorem
`thm:chiral-higher-deligne` is retitled as a conditional
two-coloured theorem with `\ClaimStatusConditional` on clauses
(2)--(3) and licensing tags \(\alpha+\beta+\gamma+\delta\).  The
new theorem-status paragraph separates the proved associator-fixed
\(E_2\) brace action from the conditional two-coloured cobar lift and
the conditional topologised \(E_3\) refinement.

The repair makes the proof obligation explicit in
`chapters/theory/chiral_higher_deligne.tex`.  The tree homotopy is the
degree \(-1\) operator
\[
h_{\mathrm{oc}}^{\mathrm{tree}}(T)
=
\sum_{e\in E_{\mathrm{int}}(T)}
(-1)^{\sigma(e,T)}T_e,\qquad
\sigma(e_j,T)=(j-1)+\sum_{v\prec e_j}|x_v|\pmod 2,
\]
where \(T_e\) inserts the Boardman--Vogt edge-length parameter and
\(\mathfrak{o}(T)=\det(\Bbbk^{E_{\mathrm{int}}(T)})\) fixes the
orientation-line sign.  It is zero on
\(\mathcal C_{\mathrm{oc}}(\mathsf{cl}^{p},\mathsf{op}^{q};\mathsf{cl})\)
for \(q>0\).  The contraction datum is now the normalized SDR
\[
d h_{\mathrm{oc}}^{\mathrm{tree}}
+h_{\mathrm{oc}}^{\mathrm{tree}}d
=\mathrm{id}-i\circ p,\qquad
p h_{\mathrm{oc}}^{\mathrm{tree}}
=h_{\mathrm{oc}}^{\mathrm{tree}}i
=(h_{\mathrm{oc}}^{\mathrm{tree}})^2=0,
\]
together with the mixed-cooperadic derivation identity and obstruction
class
\[
[o_{\mathrm{oc}}]\in
H^1\operatorname{Cone}\!\left(
\Omega((\SCchtop)^!)\to
\End(Z^{\mathrm{der}}_{\mathrm{ch}}(\cA),\cA)
\right).
\]
The proof now uses the vanishing of this class and the normalized tree
homotopy to contract the convolution complex of
\(W(\SCchtop)\)-bulk-boundary maps.  The DS-Hochschild SDR notation in
the same chapter was normalized from \(1-ip\) to
\(\mathrm{id}-i\circ p\).  The foundations summary and compute oracle
were propagated to the same formula, and the ordered-bar compute guard
now uses the reduced bar \((s^{-1}\bar A)^{\otimes k}\).

A318 adds Vol II
`compute/tests/test_chd_boardman_vogt_tree_homotopy.py`, requiring the
conditional theorem header, the signed tree formula, the open-to-closed
void, the SDR equations, the mixed derivation identity, the obstruction
class, and the shared foundations / compute summary.  It also
strengthens `compute/tests/test_ordered_chiral_kd_engine.py` to require
the `id - i circ p` formula, side conditions, and reduced-bar source
formula.

A318 verification: `python3 -m py_compile` passes on the new guard.
The focused Boardman--Vogt guard reports `6 passed`.  The widened
Chiral Higher Deligne / ordered-bar / \(\SCchtop\) chain-dioperad test
set reports `63 passed`.  Vol II `make fast` converges after two passes
in `/tmp/mkd-chiral-bar-cobar-vol2-20260615234803-76685`, producing
`2477pp` with zero undefined citations, zero undefined references, and
zero rerun requests.  Direct final-log scans find no fatal TeX errors,
undefined controls, LaTeX error headers, or unresolved
reference/citation warnings; the visible tail warnings are pdfTeX
named-destination warnings.  Vol II `make verify-licensing` exits 0
with zero blocking violations and 127 warning-class missing-tag lines
outside the touched surface.  Scoped `git diff --check` passes on the
touched files.

A317 repairs item 6 of `Research Paper Refinement vol2.pdf`: the active
Vol II Hochschild surface now defines the chiral Hochschild complex as
the coderivation complex of the chiral bar coalgebra, not as a product
of vertex-algebra homomorphism spaces.  In
`chapters/connections/hochschild.tex`,
\[
C^\bullet_{\mathrm{ch}}(A,A)
:=\Coder_0(B^{\mathrm{ch}}(A))
\simeq \Coder(\barBch(A)),\qquad
d_{\mathrm{Hoch}}=[D_A,-],
\]
where \(B^{\mathrm{ch}}(A)=k\oplus\barBch(A)\) and
\(D_A=\sum_{r\ge1}D_{m_r}\).  The later spectral-coordinate
definition is now explicitly a coordinate description of coderivation
corestrictions; it no longer identifies cochains with
\(\Hom_{\mathrm{VA}}\) homomorphisms.

A317 propagates the same correction to
`chapters/theory/foundations.tex` and `chapters/theory/axioms.tex`.
The foundations bridge now writes
\(\operatorname{Coder}_0(B^{\mathrm{ch}}(A_b))\) with differential
\([D_{A_b},-]\), and the axioms bridge replaces the old
\(\Coder(\BarTwc{A})\), \([d_B,-]\) display by the
\(B^{\mathrm{ch}}\)-coderivation model.  The strict PVA binary
formula is stated only as the \(m_1,m_2\) truncation; the full
\(A_\infty\)-chiral differential uses all \(D_{m_r}\).

A317 also repairs build-surface errors found by direct log inspection:
`\Disk_3` in `hochschild.tex` is now `\mathrm{Disk}_3`,
`A_{\nil}` in `super_chiral_yangian.tex` is now
`A_{\mathrm{nil}}`, `\effKoszul` in
`koszulness_moduli_M_kosz.tex` is in math mode, and the
`itemize`/`enumerate` closers in
`beta_N_closed_form_all_platonic.tex` are balanced.

A317 adds
`compute/tests/test_chiral_hochschild_coderivation_model.py`, requiring
the \(B^{\mathrm{ch}}\)-coderivation model and commutator differential
across the Hochschild, foundations, and axioms surfaces.  The guard
forbids the old \(\Hom_{\mathrm{VA}}\) formula, the binary-only
\(\delta_Q+\delta_{\mathrm{Hoch}}\) differential, and the
\([d_B,-]\) echo.

A317 verification: `python3 -m py_compile` passes on the guard.  The
focused guard reports `2 passed`; the widened
Hochschild/Deligne/axiom cluster reports `77 passed`.  A read-only
stack check verifies the beta-N list environments.  Vol II `make fast`
converges after two passes with `2477pp`, zero undefined citations,
zero undefined references, and zero rerun request; direct final-log
inspection finds no TeX error headers, no undefined controls, and no
undefined reference/citation warnings.  Vol II `make verify-licensing`
exits 0 with zero blocking violations and 127 warning-class missing-tag
lines outside this touched surface.  Scoped `git diff --check` and
trailing-whitespace scans pass.

A316 repairs item 5 of `Research Paper Refinement vol2.pdf`: the active
Fulton--MacPherson proof now prints the Arnold/Jacobi mechanism and
the arity-four scope correctly.  In `chapters/theory/fm-proofs.tex`,
the \(\FM_4(\C)\) assembly no longer says that "nine terms correspond
to nine non-vanishing boundary strata."  It now says that the ten
full-chain \(n=4\) terms are supported on six non-vanishing consecutive
boundary divisors:
\[
D_{\{1,2,3,4\}},\quad
D_{\{1,2,3\}},D_{\{2,3,4\}},\quad
D_{\{1,2\}},D_{\{2,3\}},D_{\{3,4\}}.
\]
When \(m_1=0\), the minimal arity-four identity is the five-term
relation
\[
m_2\circ m_3+m_3\circ m_2=0,
\]
and the first minimal-model identity involving \(m_4\) is arity five:
\[
m_2\circ m_4+m_3\circ m_3+m_4\circ m_2=0.
\]
The repaired proof also uses the reduced inner FM factors
\(\FM_2^{\mathrm{red}}(\C)\) and
\(\FM_3^{\mathrm{red}}(\C)\), and removes the contradictory
\(\epsilon(1,3)=(3-1)|a_1|\) sign paragraph.

A316 adds Corollary~\(\ref{cor:arnold-pva-jacobi-residue}\), stating
that the Arnold form identity is exactly the PVA Jacobi residue
identity at arity three.  The three residues of
\(\Omega_3^{\mathrm{sing}\text{-}\mathrm{sing}}\) on
\(D_{\{1,2\}},D_{\{2,3\}},D_{\{1,3\}}\), with signs \(+,+,-\), give
\[
\{\{a_\lambda b\}_{\lambda+\mu}c\},\qquad
\{a_\lambda\{b_\mu c\}\},\qquad
(-1)^{(\degree{a}+1)(\degree{b}+1)}
\{b_\mu\{a_\lambda c\}\}.
\]

A316 propagates the same correction to
`chapters/examples/rosetta_stone.tex`: class-\(\mathbf L\) examples no
longer say \(m_4=0\) "by the Jacobi identity."  They now say the
arity-four Lie-Jacobi source vanishes on cohomology, so the minimal
model has \(m_4^H=0\); this is not a chain-level determination of
\(m_4\) by the arity-four identity.

A316 adds Vol II `compute/tests/test_arnold_jacobi_arity_four.py`.
The focused guard reports `3 passed`; the widened
FM/\(A_\infty\)/PVA-descent/example suite reports `252 passed`.  Stale
phrase scans find the retired nine-strata and \(m_4\)-by-Jacobi
formulations only inside the negative guard.  Vol II `make fast`
converges after two passes with zero undefined citations, zero
undefined references, and zero rerun request; direct final-log
inspection finds no concrete citation or reference warning.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and 127
warning-class missing-tag lines outside the touched theorem/example
surface.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched files.

A315 repairs item 4 of `Research Paper Refinement vol2.pdf`: the
cluster-factorization proof in Vol II now uses the
Fulton--MacPherson residue calculus itself to determine the signs.  The
active proof no longer ends with an appeal that the orientation and sign
conventions live in Appendix~\(\ref{app:FM_Stokes}\).  It prints the
local normal form
\[
\omega_k=d\log\varepsilon_S\wedge\beta_S+\alpha_S,
\qquad
\Res_{D_S}(\omega_k)
=\beta_S|_{\varepsilon_S=0}
=
\lim_{\varepsilon_S\to0}
(\varepsilon_S\iota_{\partial_{\varepsilon_S}}\omega_k)|_{D_S},
\]
and explicitly rejects the naive contraction
\(\iota_{\partial_{\varepsilon_S}}\omega_k|_{D_S}\), which diverges for
logarithmic forms.

The same proof now records the boundary orientation
\[
\operatorname{or}(\FM_k(\C))
=
(-d\varepsilon_S)\wedge\operatorname{or}(D_S),
\]
uses the reduced inner cluster factors
\[
D_\pi\cong
\FM_r(\C)\times\prod_{j=1}^r\FM_{k_j}^{\mathrm{red}}(\C),
\]
and defines the collision part of the bar coderivation as
\[
D_{\mathrm{coll}}
=
\sum_{\substack{S\subset\{1,\ldots,k\}\\ |S|\ge2}}
\epsilon_{D_S}\Res_{D_S}\otimes m_S.
\]
The square-zero statement is tied to
\(\partial^2[\FM_k(\C)]=0\): at a codimension-two corner the two
iterated residues cancel because
\[
d\varepsilon_T\wedge d\varepsilon_S
=
-d\varepsilon_S\wedge d\varepsilon_T.
\]

A315 adds Vol II `compute/tests/test_fm_residue_signs_axioms.py` to
require the residue formula, the regularized contraction, the
orientation formula, the collision differential, the reduced inner FM
strata, the boundary-of-boundary source, and the normal-covector sign.
The guard forbids the retired verbal sign appeal and the phrase "up to
shuffle signs" on the active cluster theorem.

A315 verification: `python3 -m py_compile` passes on the new guard.
The focused guard reports `3 passed`; the widened
FM/\(A_\infty\)/PVA-descent/operadic suite reports `206 passed`.  Vol
II `make fast` converges after two passes with zero undefined
citations, zero undefined references, and zero rerun request; direct
final-log inspection finds no concrete citation or reference warning.
Vol II `make verify-licensing` exits 0 with zero blocking violations and
127 warning-class missing-tag lines outside the touched theorem surface.
Scoped `git diff --check` and trailing-whitespace scans pass on the
touched files.

A314 repairs item 3 of `Research Paper Refinement vol2.pdf`: the
sesquilinearity proof in Vol II now derives the formulas from
translation equivariance of chiral distribution modes in the
divided-power \(\lambda\)-bracket convention.  The proof no longer
passes through a formal \(\partial_{\lambda_i}m_k\) computation.

The active formulas are
\[
m_k(\ldots,\partial a_i,\ldots)
=-(\lambda_i+\cdots+\lambda_{k-1})m_k(\ldots,a_i,\ldots),
\qquad 1\le i<k,
\]
and
\[
m_k(a_1,\ldots,\partial a_k)
=\left(\partial_A+\sum_{j=1}^{k-1}j\lambda_j\right)
m_k(a_1,\ldots,a_k)
\]
in consecutive variables.  The proof first uses independent variables
\(\mu_i=z_i-z_k\), expands
\[
m_k(a_1,\ldots,a_k;\mu)
=
\sum_{\mathbf n}m_{k,\mathbf n}(a_1,\ldots,a_k)
\prod_{r=1}^{k-1}\frac{\mu_r^{n_r}}{n_r!},
\]
applies
\[
m_{k,\mathbf n}(\ldots,\partial a_i,\ldots)
=-n_i\,m_{k,\mathbf n-\mathbf e_i}(\ldots),
\]
and then converts
\[
\sum_{i=1}^{k-1}\mu_i
=
\sum_{i=1}^{k-1}\sum_{j=i}^{k-1}\lambda_j
=
\sum_{j=1}^{k-1}j\lambda_j.
\]
This is exactly the factor which the PDF flags as easy to lose.

A314 extends Vol II
`compute/tests/test_ainfty_chiral_coderivation_axioms.py` to require
the \(j\lambda_j\) right-input formula, the divided-power coefficient
identity, the vertex-algebra covariance \((Ta)_{(n)}=-n a_{(n-1)}\),
the independent-to-consecutive conversion, and the explicit exclusion of
formal \(\partial_{\lambda_i}\)-derivatives.

A314 verification: `python3 -m py_compile` passes on the guard.  The
focused guard reports `7 passed`.  The widened \(A_\infty\), PVA,
PVA-descent, and infrastructure suite reports `242 passed`.  Vol II
`make fast` converges after two passes with zero undefined citations,
zero undefined references, and zero rerun request; direct final-log
inspection finds no concrete citation or reference warning.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and 127
warning-class missing-tag lines outside the touched axioms surface.
Scoped `git diff --check`, trailing-whitespace scans, and stale
sesquilinearity-route scans pass on the touched files.

A313 repairs the Vol II concrete \(A_\infty\)-chiral axiom surface
against item 2 of `Research Paper Refinement vol2.pdf`.  The chapter
now begins the definition with the actual coalgebra:
\[
\mathrm B^{\mathrm{ch}}(A)=T^c(s^{-1}\bar A)
=\bigoplus_{n\ge0}(s^{-1}\bar A)^{\otimes n},
\qquad
\barBch(A)=\bigoplus_{n\ge1}(s^{-1}\bar A)^{\otimes n}.
\]
The \(n=0\) summand is the coaugmentation \(k\), while
\(\barBch(A)\) is the reduced coideal carrying the bar codifferential.

The definition now states \(m_1=Q\), the suspended degree
\(\bar a_i=|a_i|-1=|s^{-1}a_i|\), the cogenerator projection
\[
\pi D_A[s^{-1}a_1|\cdots|s^{-1}a_k]
=s^{-1}m_k(a_1,\ldots,a_k),
\]
the full coderivation formula with
\(\epsilon(i,k)=\sum_{r<i}\bar a_r\), the deconcatenation coalgebra
condition
\[
\Delta D_A=(D_A\otimes\id+\id\otimes D_A)\Delta,
\]
and the single square-zero equation \(D_A^2=0\).  The later
rectification proof now uses the completed reduced cofree coalgebra
\(\widehat{\barBch}(A)=\prod_{n\ge1}(s^{-1}\bar A)^{\otimes n}\).

A313 propagates the same reduced/coaugmented split to the main theorem
summary, the ordered-chiral-bar residue skeleton in `foundations.tex`,
and the active Heisenberg Rosetta Stone ordered-bar cofibration.  It
also repairs the stale `thm:SC-self-duality` reference to the active
nonselfdual Swiss-cheese Koszul-duality proposition and adds four
verified prefixed Vol I fallback anchors exposed by direct log
inspection.

A313 adds Vol II
`compute/tests/test_ainfty_chiral_coderivation_axioms.py`, requiring
the first coalgebra display, the \(m_k\)-type and degree, \(m_1=Q\),
the suspended sign, coderivation corestriction, the coalgebra identity,
\(D_A^2=0\), the completed reduced cofree coalgebra, spectral
substitution in the main Stasheff identity, the propagated
reduced/coaugmented split, the nonselfdual Swiss-cheese reference, and
the four Vol I anchors.

A313 verification: `python3 -m py_compile` passes on the new guard.
The focused guard reports `6 passed`.  The widened \(A_\infty\), PVA,
\(\SCchtop\), and chiral-Higher-Deligne test set reports `151 passed`.
Vol II `make fast` converges after two passes with zero undefined
citations, zero undefined references, and zero rerun request; direct
final-log inspection finds no concrete citation or reference warning.
Vol II `make verify-licensing` exits 0 with zero blocking violations
and zero warnings.  Scoped `git diff --check`, trailing-whitespace
scans, stale bar-collapse scans, and stale self-duality-label scans
pass on the touched files.

A312 repairs the Vol II local \(\SCchtop\) heptagon definition against
the first explicit demand of `Research Paper Refinement vol2.pdf`: the
two-coloured object must be a dg coloured dioperad, not a named
metaphor.  The active definition already printed the logarithmic chain
components, the zero open-to-closed chain complex, the local
Fulton--MacPherson and ordered-interval composition maps, and a
Boardman--Vogt arity formula.  The remaining defect was the sign line:
the tree orientation only recorded \(\det(k^{E_{\mathrm{int}}(T)})\).

The repaired definition now uses the full orientation convention
\[
\mathfrak{o}(T)
=\det(k^{E_{\mathrm{int}}(T)})
\otimes\bigotimes_{v\in V(T)}\mathfrak{o}(v),
\]
and states that the \(d_{\mathrm{contract}}\) and
\(d_{\mathrm{compose}}\) signs move the edge through both the edge
determinant and the vertex-orientation tensor factors.  This is the
chain-level datum that turns two-coloured coherence into a signed
Boardman--Vogt statement.

A312 adds Vol II `compute/tests/test_sc_chtop_chain_dioperad.py`
requiring the logarithmic chain coloured dioperad components, the zero
open-to-closed chain complex, both local composition maps, the
\(W(\SCchtop)\) arity formula, the
\(d_{\mathrm{chain}}+d_{[0,1]}+d_{\mathrm{contract}}+
d_{\mathrm{compose}}\) differential, and the vertex-orientation tensor
factors.

A312 verification: `python3 -m py_compile` passes on the new guard.
The focused guard reports `4 passed`.  The widened factorization,
heptagon, chiral-Higher-Deligne, and Swiss-cheese climax-wave test set
reports `93 passed`.  Vol II `make fast` converges after two passes
with zero undefined citations and references.  Vol II
`make verify-licensing` exits 0 with zero blocking violations and zero
warnings.  Scoped `git diff --check` and trailing-whitespace scans pass
on the touched Vol II file and guard.

A311 repairs the active Vol I frontier modular-holography chapter
`frontier_modular_holography_platonic.tex`.  Four frontier formulas
still wrote the protected scalar Faber--Pandharipande coefficient as
the bare full coefficient
\[
F_g=\kappa\lambda_g^{\mathrm{FP}}:
\]
the independent-sum additivity proof, the Nekrasov/gauge-origami
comparison, the twisted D3 amplitude statement, and the Burns-space
genus tower.  Nearby height-interpolation prose also said that
\(\kappa\) controls \(F_g\) on the uniform-weight lane.  This erased
the scalar protected sector versus all-weight graph/cross-channel
distinction.

The repaired frontier surface now writes additivity for
\(F_g^{\mathrm{sc}}\), interprets gauge origami as the scalar vacuum
sector, makes the D3 statement a protected scalar amplitude on the
uniform-weight scalar lane, types the Burns values as
\(F_i^{\mathrm{sc}}(\cA_{\mathrm{Burns}})\), and says that
\(\kappa\) controls \(F_g^{\mathrm{sc}}\) on the uniform-weight scalar
lane.  The Burns \(T\)-line planted-forest correction remains visible
as non-scalar/full data.

A311 adds Vol I
`compute/tests/test_frontier_modular_holography_scalar_typing.py`
requiring the typed additivity/origami, D3, Burns, and height
statements.  The guard forbids bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\), the old `at the scalar level`
variant, `F_g on the uniform-weight lane`, and the old `shadow genus
expansion $F_g=...` phrase.

A311 verification: `python3 -m py_compile` passes on the new guard.
The focused frontier guard reports `3 passed`; paired with the
A300--A310 scalar/full guards it reports `39 passed`; the broader
scalar-typing guard set reports `53 passed`.  The frontier and twisted
holography computation tests report `261 passed`.  Targeted stale-form
scans over the frontier file are clean.  The broad Vol I scalar/full
scan over `chapters/connections`, `chapters/theory`, and
`compute/tests` now returns no hits for the
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\) collapse class repaired in
A300--A311.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A310 repairs the active Vol I entanglement/modular-Koszul chapter
`entanglement_modular_koszul.tex`.  The BTZ theorem proof and the JT
evidence paragraph still wrote the scalar Faber--Pandharipande input
as the bare full free-energy coefficient
\[
F_g=\kappa\lambda_g^{\mathrm{FP}}.
\]
This is false as an all-weight Virasoro class~M statement: the
planted-forest/cross-channel terms contribute beyond the scalar lane at
\(g\ge2\).  The same file also retained three untyped
`uniform-weight lane` phrases attached to scalar formulas.

The repaired BTZ proof now states
\[
F_g^{\mathrm{sc}}=\kappa\lambda_g^{\mathrm{FP}},\qquad
F_g=F_g^{\mathrm{sc}}+\delta F_g^{\mathrm{cross}},
\]
and then uses only scalar closed-sector inputs for the saddle
comparison.  The JT evidence paragraph now compares the leading-order
JT amplitude with the scalar shadow free energy on the uniform-weight
scalar lane.  The earlier same-file references now say
`uniform-weight scalar lane`.

A310 adds Vol I
`compute/tests/test_entanglement_modular_koszul_scalar_typing.py`
requiring the typed BTZ/JT scalar formulas, the full/cross
decomposition, the scalar closed-sector input phrase, the
uniform-weight scalar lane, and the Virasoro class~M planted-forest
warning.  The guard forbids bare
\(F_g=\kappa\lambda_g^{\mathrm{FP}}\), the old `shadow free energy
$F_g=...` phrase, the untyped `uniform-weight lane`, and
\(F_g^{\mathrm{scal}}\).

A310 verification: `python3 -m py_compile` passes on the new guard.
The focused entanglement guard reports `3 passed`; paired with the
A300--A309 scalar/full guards it reports `36 passed`; the broader
scalar-typing guard set reports `50 passed`.  Targeted stale-form
scans over `entanglement_modular_koszul.tex` are clean.  The broad
Vol I scalar/full scan no longer reports
`entanglement_modular_koszul.tex`; the remaining debt is four hits in
`frontier_modular_holography_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A309 repairs the active Vol I arithmetic-shadow chapter
`arithmetic_shadows.tex`.  The Chern--Simons/Bernoulli remark compared
the shadow \(\hat A\)-genus coefficient with perturbative
Chern--Simons Bernoulli factors, but still wrote the comparison as the
bare full coefficient
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
The same file also retained two same-convention leaks: a lattice
sentence using bare \(F_g(V_\Lambda)=r\lambda_g^{\mathrm{FP}}\) with
an `UNIFORM-WEIGHT` marker, and an oscillator-factorization sentence
using \(F_g^{\mathrm{scal}}\).

The repaired arithmetic-shadow surface now writes the lattice and
oscillator coefficients as
\[
F_g^{\mathrm{sc}}(V_\Lambda)=r\lambda_g^{\mathrm{FP}},
\]
and the Bernoulli comparison as
\[
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
The closing warning now says that \(F_g^{\mathrm{sc}}\), not bare
\(F_g\), is the genus-\(g\) scalar projection of the shadow obstruction
tower, distinct from the degree-\(r\) shadow coefficient \(S_r\).

A309 adds Vol I
`compute/tests/test_arithmetic_shadows_scalar_typing.py` requiring the
typed lattice and Bernoulli scalar FP formulas and the
\(F_g^{\mathrm{sc}}\) versus \(S_r\) warning.  The guard forbids bare
\(F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}\), bare
\(F_g(V_\Lambda)=r\lambda_g^{\mathrm{FP}}\), `UNIFORM-WEIGHT`,
\(F_g^{\mathrm{scal}}\), and the old `Note: $F_g$ is...` phrase.

A309 verification: `python3 -m py_compile` passes on the new guard.
The focused arithmetic-shadows guard reports `3 passed`; paired with
the A300--A308 scalar/full guards it reports `33 passed`; the broader
scalar-typing guard set reports `47 passed`.  Targeted stale-form
scans over `arithmetic_shadows.tex` are clean.  The broad Vol I
scalar/full scan no longer reports `arithmetic_shadows.tex`; the
remaining debt is six hits in `entanglement_modular_koszul.tex` and
`frontier_modular_holography_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A305 repairs the core Vol I `bar_construction.tex` summary of
Theorem~D inside the nilpotence-periodicity correspondence.  The
chapter still wrote the Theorem~D genus expansion as the bare full
coefficient
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
\]
with a parenthetical uniform-weight qualifier and a separate
multi-weight correction sentence.  That compressed three distinct
objects: the scalar Faber--Pandharipande coefficient, the full
coefficient on exact uniform-weight lanes, and the full coefficient on
multi-weight lanes.

The repaired paragraph now states
\[
F_g^{\mathrm{sc}}(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
then names the exact-lane equality
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
\]
only for uniform-weight exact lanes, and the multi-weight full
coefficient
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  +\delta F_g^{\mathrm{cross}}(\cA)
\]
for \(g\ge2\).

A305 adds Vol I
`compute/tests/test_bar_construction_scalar_typing.py` requiring the
scalar FP anchor, exact-lane equality, multi-weight decomposition, and
the existing fiberwise-curvature anchor
\(\dfib^2=\kappa(\cA)\omega_g\).  The guard forbids the retired bare
full-coefficient FP formula, the old genus-expansion phrase, the
`UNIFORM-WEIGHT` parenthetical, and the old "multi-weight extension
receives" wording.

A305 verification: `python3 -m py_compile` passes on the new guard.
The focused bar-construction guard reports `3 passed`; paired with the
A300--A304 scalar/full guards it reports `21 passed`; the broader
scalar-typing guard set reports `35 passed`.  Targeted stale-form scans
over `bar_construction.tex` and the guard are clean.  The broad Vol I
scalar/full scan no longer reports `bar_construction.tex`; the
remaining debt is in `koszulness_vii_multiweight_platonic.tex`,
`frontier_modular_holography_platonic.tex`,
`entanglement_modular_koszul.tex`, `fourier_seed.tex`,
`arithmetic_shadows.tex`, and `shadow_L_function_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A306 repairs the active Vol I characterization~(vii) surface
`koszulness_vii_multiweight_platonic.tex`.  The chapter already
contained the multi-weight graph-sum correction, but its displayed
formulas still wrote the diagonal Faber--Pandharipande term as a bare
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
That conflated the diagonal scalar term, the scalar trace after the
cross-channel graph sum, and the exact equality on uniform-weight
lanes.

The repaired chapter now writes the multi-weight scalar trace as
\[
F_g(\cA)
  =F_g^{\mathrm{sc}}(\cA)
  +\delta F_g^{\mathrm{cross}}(\cA),
\qquad
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]
The proof identifies
\[
F_g^{\mathrm{sc}}(\cA)
  =\sum_i\kappa_i\lambda_g^{\mathrm{FP}}
  =\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
and the uniform-weight corollary now states
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]

A306 adds Vol I
`compute/tests/test_koszulness_vii_multiweight_scalar_typing.py`
requiring the multi-weight decomposition, diagonal scalar FP term,
channel-sum computation of \(F_g^{\mathrm{sc}}\), exact
uniform-weight equality, \(\delta F_1^{\mathrm{cross}}=0\), and the
recovered equality \(F_g=F_g^{\mathrm{sc}}\).  The guard forbids bare
\(F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}\) and the retired
scalar-formula-without-correction wording.

A306 verification: `python3 -m py_compile` passes on the new guard.
The focused characterization~(vii) guard reports `3 passed`; paired
with the A300--A305 scalar/full guards it reports `24 passed`; the
broader scalar-typing guard set reports `38 passed`.  Targeted
stale-form scans over `koszulness_vii_multiweight_platonic.tex` and
the guard are clean.  The broad Vol I scalar/full scan no longer
reports `koszulness_vii_multiweight_platonic.tex`; the remaining debt
is in `frontier_modular_holography_platonic.tex`,
`entanglement_modular_koszul.tex`, `fourier_seed.tex`,
`arithmetic_shadows.tex`, and `shadow_L_function_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A307 repairs the Vol I shadow \(L\)-function genus-slot bridge
`shadow_L_function_platonic.tex`.  The theorem already rejected the
false evaluation
\[
F_g\leftrightarrow L^{\mathrm{sh}}(1-2g)
\]
and used the correct coefficient projector at \(r=2g-2\).  Its genus
formulas, however, still wrote the uniform-weight and all-weight
Faber--Pandharipande term as a bare
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}.
\]

The repaired theorem now writes the uniform-weight lane as
\[
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}},
\]
and the all-weight lane as
\[
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}},\qquad
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  +\delta F_g^{\mathrm{cross}}(\cA).
\]
The coefficient-extraction statement is unchanged.

A307 adds Vol I `compute/tests/test_shadow_l_function_scalar_typing.py`
requiring the typed uniform-weight equality, scalar FP anchor,
all-weight decomposition, coefficient projector, and shadow-slot
projector.  The guard forbids bare
\(F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}\), the old all-weight
bare-plus-cross formula, and the retired "scalar formula is replaced
by" wording.

A307 verification: `python3 -m py_compile` passes on the new guard.
The focused shadow-\(L\) guard reports `3 passed`; paired with the
A300--A306 scalar/full guards it reports `27 passed`; the broader
scalar-typing guard set reports `41 passed`.  Targeted stale-form scans
over `shadow_L_function_platonic.tex` and the guard are clean.  The
broad Vol I scalar/full scan no longer reports
`shadow_L_function_platonic.tex`; the remaining debt is in
`fourier_seed.tex`, `arithmetic_shadows.tex`,
`entanglement_modular_koszul.tex`, and
`frontier_modular_holography_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A308 repairs the active Vol I Fourier seed theorem surface
`fourier_seed.tex`.  In Theorem~\ref{thm:fourier-four-properties}, the
characteristic-function clause already had the correct cohomological
obstruction class
\[
\mathrm{obs}_g(\cA)=\kappa(\cA)\lambda_g,
\]
but the following all-weight numerical line still compressed the scalar
Faber--Pandharipande coefficient and the full all-weight coefficient as
\[
F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}
  +\delta F_g^{\mathrm{cross}}(\cA).
\]

The repaired sentence now separates them:
\[
F_g^{\mathrm{sc}}(\cA)
  =\kappa(\cA)\lambda_g^{\mathrm{FP}},\qquad
F_g(\cA)=F_g^{\mathrm{sc}}(\cA)
  +\delta F_g^{\mathrm{cross}}(\cA).
\]
The \(\hat A\)-series remains the scalar-lane generating function.

A308 adds Vol I `compute/tests/test_fourier_seed_scalar_typing.py`
requiring the obstruction-class formula, scalar FP coefficient,
all-weight decomposition, \(\delta F_1^{\mathrm{cross}}=0\), and the
scalar-lane generating-function phrase.  The guard forbids bare
\(F_g(\cA)=\kappa(\cA)\lambda_g^{\mathrm{FP}}\) and the old bare
all-weight line.

A308 verification: `python3 -m py_compile` passes on the new guard.
The focused Fourier-seed guard reports `3 passed`; paired with the
A300--A307 scalar/full guards it reports `30 passed`; the broader
scalar-typing guard set reports `44 passed`.  Targeted stale-form scans
over `fourier_seed.tex` and the guard are clean.  The broad Vol I
scalar/full scan no longer reports `fourier_seed.tex`; the remaining
debt is in `arithmetic_shadows.tex`, `entanglement_modular_koszul.tex`,
and `frontier_modular_holography_platonic.tex`.
Vol I `make fast` converges after three passes with zero undefined
citations, one stable undefined reference, and no rerun request.  Vol I
`make verify` exits 0 with `verify_edit: all checks passed`.  Scoped
`git diff --check` and trailing-whitespace scans pass on the touched
Vol I file and guard.

A320 repairs item 9's 3d HT bridge.  The live theorem
`thm:physics-bridge` now starts from the explicit BV datum
\[
(\mathcal E,Q,\omega_{\mathrm{BV}},I,B),\qquad
\mathcal E\simeq
\Omega^\bullet(\R_t)\widehat\otimes\Omega^{0,\bullet}(\C_z)
\otimes\mathfrak a[1],\qquad
Q=d_t+\dbar_z+d_{\mathfrak a},
\]
with \(B\) a \(Q\)-compatible Lagrangian boundary condition and
\(\{S_{\mathrm{cl}},S_{\mathrm{cl}}\}_{\mathrm{BV}}=0\).
The factorized-parametrix hypothesis now requires the kernel homotopy
\[
Q_{\boxtimes}P_{\epsilon<L}=K_L-K_\epsilon
\]
and the finite product HT expansion of \(P_{\mathrm{sing}}\), rather
than the retired \(QG=\delta_{\C}\boxtimes\delta_{\R}\) shorthand.

The construction chapter `bv-construction.tex` now displays the same
kernel identity, defines \(Q_{\boxtimes}\) on distributional kernels,
adds `eq:3dht-propagator-kernel-homotopy`,
\[
[Q,P_{\epsilon<L}]=K_L-K_\epsilon,
\]
and keeps the scale-QME
\[
QI[L]+\hbar\Delta_LI[L]+\frac12\{I[L],I[L]\}_L=0
\]
before boundary \(A_\infty\)-chiral operations are extracted.  The
concordance table, the theorem's hypothesis summary, the affine
half-space checklist, and the IV metadata now all name the BV datum and
kernel homotopy.

A320 adds `compute/tests/test_bv_ht_bridge_bv_datum.py`, which guards
the tuple, product HT field complex, CME, propagator kernel homotopy,
factorized singular part, scale-QME, and \(Q_{\mathrm{BRST}}^L\)
formula, and forbids the retired `QG=\delta_{\C}\boxtimes\delta_{\R}`
form in the theorem.

A320 verification: `python3 -m py_compile` passes on the new and IV
guards.  Focused pytest reports `24 passed`; widened bridge/KZ/Obs
tests report `36 passed`; affine and nearby IV tests report `82
passed`.  `make fast` converges after two passes with 2479 pages, zero
undefined citations, zero undefined references, and zero rerun
requests.  Direct final-log scans find no fatal TeX errors, undefined
controls, or unresolved citation/reference warnings.  `make
verify-licensing` exits 0 with zero blockers and the existing 127
warning-class missing-tag lines outside this touched surface.  Scoped
`git diff --check` passes on the touched files.

A321 strengthens item 10 on the DK-0 standard-family theorem.  The
former `DK-0 Laplace bridge for five families; \ClaimStatusProvedHere`
surface in `chapters/examples/examples-worked.tex` is now
`Family-gated DK-0 Laplace bridge; \ClaimStatusProvedHereConditional`
with licensing tags \(\alpha+\beta+\gamma\).  The theorem still proves
the coefficient identity
\[
r^{L,ab}(z)=\sum_{n\ge0}c_n^{ab}\frac{n!}{z^{n+1}},
\qquad
\{a_\lambda b\}=\sum_{n\ge0}c_n^{ab}\lambda^n,
\]
but only under printed family gates.

The new gate table separates:
Heisenberg \(\{J_\lambda J\}=k\lambda\) and
\(r^{\mathrm{coll}}_{\cH}=k\Omega_{\cH}/z\);
affine trace-form normalization from the KZ representative
\(\Omega/((k+h^\vee)z)\), with the critical KZ surface excluded;
Virasoro pre-\(d\log\) kernel from collision residue
\((c/2)z^{-3}+2Tz^{-1}\), with no full fusion or
descendant-sensitive spectral \(R\)-operator asserted;
the free-field \(\beta\gamma/bc\) parity/statistics-exchange gate; and
the standard Zamolodchikov \(\mathcal W_3\) normalization away from the
singular denominator \(5c+22=0\).

A321 adds `compute/tests/test_dk0_family_normalization_gates.py`, which
requires the family-gated theorem title, the gate table, all five
family-specific gates, and a toy admissibility check that fails if any
family gate is absent.  Fixed-string scans find no remaining copy of
the retired ungated theorem title; the phrase "For each of the five
standard families" survives only inside the gated theorem surface.

A321 verification: `python3 -m py_compile` passes on the new guard and
nearby example/spectral tests.  Focused pytest reports `123 passed`
across `test_dk0_family_normalization_gates.py`, `test_examples.py`,
`test_collision_residue_rmatrix.py`, and `test_spectral_braiding.py`.
`make fast` converges after two passes with 2479 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings.  `make verify-licensing` exits
0 with zero blockers and the existing 127 warning-class missing-tag
lines outside this touched surface.  Scoped `git diff --check` passes
on the touched files.

A322 continues item 10 on the \(\mathcal W_N\) self-dual-halving
surface.  The live `w-algebras-frontier.tex` proposition and the
parallel `w-algebras-stable.tex` split copy no longer write the
curvature halving as bare full-genus halving for multi-weight
\(\mathcal W_N\).  The genus-one statement is now scalar:
\[
F_1^{\mathrm{sc}}(c^*)
  =\tfrac12F_1^{\mathrm{sc}}(\alpha_N).
\]
For \(N=2\) (Virasoro, uniform-weight), the scalar lane is the full
genus tower, so \(F_g(c^*)=\tfrac12F_g(\alpha_N)\) remains a full
all-genus statement.

For \(N\ge3\), the text now displays the full multi-weight
decomposition
\[
F_g(\cW_{N,c})
 =
F_g^{\mathrm{sc}}(\cW_{N,c})
 +\delta F_g^{\mathrm{cross}}(\cW_{N,c}),
\]
and states that full halving would require the separate condition
\[
\delta F_g^{\mathrm{cross}}(\cW_{N,c^*})
=\tfrac12\,\delta F_g^{\mathrm{cross}}(\cW_{N,\alpha_N}),
\]
which is not asserted.

A322 adds `compute/tests/test_w_algebra_self_dual_halving_scope.py`,
requiring the scalar-lane statement, the full multi-weight
decomposition, the displayed cross-channel equality as a non-asserted
extra requirement, and a toy model where scalar halving does not imply
full halving.  Exact stale scans find no surviving live/split occurrence
of the retired bare \(F_1(c^*)\) halving phrases or the old
"unconditionally for all families" wording on the checked W-algebra
surface; the only remaining `F_1(c^*)` hits are negative assertions in
the guard.

A322 verification: `python3 -m py_compile` passes on the new guard and
nearby W-algebra tests.  Focused guard pytest reports `3 passed`; the
widened W-algebra/example bundle reports `91 passed`.  `make fast`
converges after two passes with 2479 pages, zero undefined citations,
zero undefined references, and zero rerun requests.  Direct final-log
scans find no fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.  `make verify-licensing` exits 0 with zero
blockers and the existing 127 warning-class missing-tag lines outside
this touched surface.  Scoped `git diff --check` passes on the touched
files.

A323 continues item 10 on the Pixton/MC-tower surface.  The live
theorem `thm:pixton-semisimple` in
`chapters/connections/ht_bulk_boundary_line.tex` no longer claims that
one semisimple MC-tower argument covers all standard families at generic
level.  The theorem is now:

`Pixton ideal from the scalar MC tower on the semisimple r-spin locus;
\ClaimStatusProvedHereConditional; licensing tags
\alpha+\beta+\gamma+\delta`.

The statement now begins with a one-parameter chirally Koszul scalar
family \(\{\mathcal A_\kappa\}_{\kappa\in U}\) whose domain contains all
\[
\kappa_r=\frac{r-2}{2r},\qquad r\ge3.
\]
It assumes:
1. the scalar shadow of \(\mathcal A_{\kappa_r}\) is a semisimple
   rank-one CohFT with flat unit;
2. the Givental--Teleman \(R\)-matrix is the
   Faber--Shadrin--Zvonkine \(r\)-spin \(R\)-matrix;
3. the scalar modular MC map sends boundary strata to the tautological
   relations of that CohFT.

The conclusion is the union over all \(r\ge3\) of scalar MC-descended
relations, which by FSZ and PPZ generates the Pixton ideal.  The text
now explicitly says that this is only a scalar-lane statement: it does
not assert a full tensor-channel theorem for Virasoro or
\(\mathcal W_N\), and admissible-level, resonant, and logarithmic
non-semisimple specialisations need a separate replacement for
Givental--Teleman semisimplicity.

A323 adds `compute/tests/test_pixton_semisimple_gate.py`.  The guard
requires the conditional theorem title and licensing tags, the
\(\kappa_r\) family condition, the FSZ/PPZ/Givental--Teleman gates, and
the scalar/full-channel exclusions.  It also includes a toy rational
check: a fixed \(\kappa=\kappa_3\) is a strict subset of the values
\(\{\kappa_r\}_{r=3}^{9}\), so a single fixed algebra cannot supply the
all-\(r\) PPZ route.

A323 verification: `python3 -m py_compile` passes on the new guard.
Focused guard pytest reports `4 passed`; the widened related bundle
reports `170 passed` across `test_pixton_semisimple_gate.py`,
`test_deletion_ledger.py`, `test_modular_pva_quantization.py`, and
`test_semisimple_purity.py`.  Fixed-string stale scans find no remaining
live copy of the retired "all standard families at generic level"
sentence or the old single-algebra theorem opening outside the new
negative guard assertions.  `make fast` converges after two passes with
2479 pages, zero undefined citations, zero undefined references, and
zero rerun requests.  Direct final-log scans find no fatal TeX errors,
undefined controls, or unresolved citation/reference warnings.
`make verify-licensing` exits 0 with zero blockers and the existing 127
warning-class missing-tag lines outside this touched surface.  Scoped
`git diff --check` passes on the touched files.

A324 continues item 10 on the split modular-PVA extension surface.  The
file `chapters/connections/thqg_modular_pva_extensions.tex` is not in
the live `main.tex` input graph, but it is a merge source and previously
advertised two global statements that the live chapter does not justify:

1. `HS-sewing for the standard landscape; \ClaimStatusProvedHere`;
2. `Extended quantization dictionary; \ClaimStatusProvedHere`, with a
   table saying every standard family has a one-dimensional relevant
   \(H^1\) and an existing full modular lift.

The HS statement is now:

`HS-sewing under polynomial/subexponential analytic gates;
\ClaimStatusConditional; licensing tags \gamma+\delta`.

It requires a positive-energy Hilbert completion, subexponential sector
growth, polynomial Hilbert--Schmidt OPE growth, and exclusion of
resonant/logarithmic/non-positive-energy/non-principal DS
specialisations unless separate analytic norm and character bounds are
supplied.  The \(W_N\) row is conditional on the DS character formula
for the chosen principal generic-level completion and on polynomial
bounds for normalized \(W_N\) structure constants.

The dictionary is now:

`Gated modular-PVA quantization dictionary;
\ClaimStatusConditional; licensing tags \alpha+\gamma+\delta`.

Its table no longer says that every row has a constructed full
all-loop lift.  It records sector-specific gates: fixed lattice has no
universal \(\partial_k\), affine requires KZ analytic SDR, Virasoro is a
scalar-lane genus-one torsor, \(W_3\) is conditional on normal form plus
the scalar one-loop calculation, \(W_4,W_5\) are finite-window checked
sectors, and \(W_N\) remains a conjectural inductive DS/KZ gate rather
than a constructed full lift.

A324 adds `compute/tests/test_modular_pva_extension_gates.py`, which
guards the conditional HS theorem, the analytic hypotheses, the
non-principal exclusions, the gated dictionary title, the fixed-lattice
no-\(\partial_k\) row, the \(W_N\) not-constructed row, and a toy check
showing that a fixed lattice row cannot satisfy the old universal
single-parameter claim.

A324 verification: `python3 -m py_compile` passes on the new guard.
Focused guard pytest reports `4 passed`; the widened modular-PVA bundle
reports `195 passed` across `test_modular_pva_extension_gates.py`,
`test_modular_pva_quantization.py`, `test_modular_obstruction_engine.py`,
and `test_p3_pva_quantum_lift.py`.  Fixed-string stale scans find no
surviving copy in the split file of the retired uniform standard-family
phrases.  `make verify-licensing` exits 0 with zero blockers and the
existing 127 warning-class missing-tag lines outside this touched
surface.  `make fast` converges after two passes with 2479 pages, zero
undefined citations, zero undefined references, and zero rerun requests.
Direct final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings.  Scoped `git diff --check`
passes on the touched files.

A325 repairs the live class-\(\mathsf M\) holography corollary in
`chapters/theory/chiral_higher_deligne.tex`.  The old corollary stated
class-M universal holography as if the derived chiral centre were
already the 3d HT bulk and then said the Universal Holography functor
was chain-level across all four classes.  The proof silently used the
weight-completed/pro ambient, DS--Hochschild transport, and an HT
bulk--Hochschild comparison map.

The corollary is now:

`Class-M bulk-centre comparison in the weight-completed ambient;
\ClaimStatusConditional; licensing tags \beta+\gamma+\delta`.

The statement assumes:
1. strict Mittag--Leffler weight-completed/pro/\(J\)-adic ambient
   \(\hypAmbientWtCpl\);
2. DS--Hochschild transport for the chosen principal or hook-type
   reduction, including cover descent for fractional good gradings;
3. a constructed 3d HT BV/factorization model and a bulk--Hochschild
   comparison map
   \[
   \chi_{\mathrm{HT}}\colon
   Z^{\mathrm{der}}_{\mathrm{ch}}(\mathcal A)
   \to
   \mathcal A_{\mathrm{bulk}}^{\mathrm{3d\,HT}}(\mathcal A);
   \]
4. the heptagon hypotheses and the two-coloured cobar homotopy for the
   chain-level mixed refinement.

The conclusion is now split correctly.  \(\chi_{\mathrm{HT}}\) is a
cohomological \(E_2\)-chiral/Gerstenhaber equivalence.  Only under the
fourth input does it refine to a mixed chiral-topological chain-level
equivalence.  The statement explicitly records that the bounded
direct-sum class-M chain-level identification is false because
\[
S_4(\mathrm{Vir}_c)=\frac{10}{c(5c+22)}\ne0
\]
populates bar-weight \(4\), and that G/L/C original-complex
identifications are separate inputs rather than consequences of the
class-M comparison.

A325 adds `compute/tests/test_chd_class_m_holography_gates.py`, which
guards the conditional title and tags, the four input gates, the
cohomology/chain-level split, the direct-sum falsity sentence, and the
absence of the retired unconditional phrases.

A325 verification: `python3 -m py_compile` passes on the new guard and
nearby CHD/class-M tests.  Focused pytest reports `73 passed` across
`test_chd_class_m_holography_gates.py`,
`test_chiral_higher_deligne.py`, `test_deletion_ledger.py`, and
`test_class_m_ambient_iv.py`.  The widened CHD/holography bundle reports
`75 passed` after adding `test_universal_holography_functor.py`.
Fixed-string stale scans find no surviving live copy of the retired
unconditional universal-holography phrases outside the negative guard.
`make verify-licensing` exits 0 with zero blockers and the existing 127
warning-class missing-tag lines outside this touched surface.  `make
fast` converges after two passes with 2479 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings.  Scoped `git diff --check`
passes on the touched files.

A326 repairs the live Universal Holography climax theorem in
`chapters/connections/programme_climax_platonic.tex`.  The old live
statement still treated \(A\) in the non-logarithmic standard landscape
as enough data to produce a canonical three-dimensional HT gauge theory
\(T_A\).  That collapsed the algebraic stage-chain candidate with an
actual Costello--Gwilliam HT BV/factorization realization.

The master theorem is now \(\Xi\)-relative:
\[
\Phihol\colon
\ChirAlgclimax^{\omega,\mathrm{BL,adm},\Xi}_{X}
\to
\HTQFTclimax^{\Xi}_{X\times\mathbb R}.
\]
The datum \(\Xi(A,\omega,\mathrm{BL})\) consists of a
Costello--Gwilliam HT BV/factorization model \(\mathcal F_A\), a
boundary comparison \(\eta^\partial_A\), a bulk--Hochschild comparison
\(\chi_{\mathrm{HT},A}\), the ambient \(\mathcal A(A)\), the
class-\(\mathbf M\) gate \(\hypAmbientWtCpl\), DS--Hochschild and
cover-descent data when needed, and the endpoint gates
\(\hypProchazka,\hypCKL,\hypPRSh,\hypYamada\) for
\(W_\infty[\lambda]\).

The pointwise theorem `thm:programme-climax` now says that, assuming
\(\Xi(A)\), the supplied model \(T_A:=\mathcal F_A\) is the value of
\(\Phihol\).  The theorem no longer asserts a bare canonical physical
theory from \(A\) alone.  Boundary recovery is the map
\(\eta^\partial_A\), and bulk recovery is the map
\(\chi_{\mathrm{HT},A}\) in the declared ambient.

The opening display and bicoloured-primitive factorization corollary
were propagated to the same \(\Xi\)-typed source/target.  Functoriality
now preserves Drinfeld--Sokolov reductions and Heisenberg-coset
extractions only when the \(\Xi\)-data commute with those operations.

A326 adds `compute/tests/test_programme_climax_holography_gates.py`,
which guards the \(\Xi\)-relative source/target, the
Costello--Gwilliam HT realization, the bulk--Hochschild comparison map,
the class-\(\mathbf M\) and endpoint gates, and absence of the retired
canonical-existence phrase.

A326 verification: `python3 -m py_compile` passes on the new guard.
The virtualenv pytest bundle reports `68 passed` across
`test_programme_climax_holography_gates.py`, `test_deletion_ledger.py`,
`test_chd_class_m_holography_gates.py`, and
`test_universal_holography_functor.py`.  `make verify-licensing` exits
0 with zero blockers and 124 warning-class missing-tag lines elsewhere;
the touched theorem/corollary are no longer among the warnings.  `make
fast` converges after two passes with 2479 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings; existing pdfTeX destination
warnings remain in the log tail.

A327 propagates the \(\Xi\)-relative theorem form through the live rung
corollaries and downstream trace statements in
`chapters/connections/programme_climax_platonic.tex`.  The old rung
statements still wrote \(\Phihol(H_k)\), \(\Phihol(V_k(\mathfrak g))\),
\(\Phihol(\Vir_c)\), \(\Phihol(W_{N,c})\), and
\(\Phihol(W_\infty[\lambda])\), and their proofs still specialised the
master theorem at \(A\) alone.  That notation silently erased the
realization datum installed in A326.

The sequent display is now:
\[
\Xi(A)\vdash_{\mathcal A(A)}
\bigl(\Phihol(A,\Xi_A)=T_{A,\Xi},\;
\mathrm{Obs}^{\partial}(T_{A,\Xi})\simeq A,\;
\mathrm{Obs}^{\mathrm{bulk}}(T_{A,\Xi})\simeq
Z^{\mathrm{der}}_{\mathrm{ch}}(A)\bigr).
\]
The rung corollaries now carry \(\Xi_{H_k}\),
\(\Xi_{V_k(\mathfrak g)}\), \(\Xi_{\Vir_c}^{\mathrm{BH}}\),
\(\Xi_{W_N}\), and \(\Xi_{W_\infty,\lambda}\).  Each proof applies the
master theorem to the corresponding quadruple
\((A,\omega,\mathrm{BL},\Xi)\).  The F3 frontier, gravity trace theorem,
Conway superselection remark, Vol I/II/III comparison paragraph, and
closing bar-complex summary were propagated to the same discipline.

A327 extends `compute/tests/test_programme_climax_holography_gates.py`
and updates `compute/tests/test_deletion_ledger.py` so the guards now
require the typed sequent and typed rung invocations, and reject the old
`\Phihol(A)=T_A` and `Specialise ... at` patterns.

A327 verification: `python3 -m py_compile` passes on the touched guard
tests.  The focused virtualenv pytest bundle reports `69 passed` across
`test_programme_climax_holography_gates.py`, `test_deletion_ledger.py`,
`test_chd_class_m_holography_gates.py`, and
`test_universal_holography_functor.py`.  Fixed-string scans find no
surviving live copy of `\Phihol(A)=T_A` or `Specialise
\cref{thm:universal-holography-master} at`.  `make verify-licensing`
exits 0 with zero blockers and 123 warning-class missing-tag lines
elsewhere.  `make fast` converges after two passes with 2481 pages, zero
undefined citations, zero undefined references, and zero rerun requests.
Direct final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings; existing pdfTeX destination
warnings remain in the log tail.  Scoped `git diff --check` passes.

A349 repairs `prop:gravity-m4`, the next gravity-climax theorem surface
after A348.  The displayed quaternary Virasoro formula is no longer
presented as a canonical cohomology operation determined by the
arity-four identity.  It is now the Stasheff-gauge chain representative
selected after fixing consecutive collision coordinates, the
divided-power \(\lambda\)-bracket convention, the completed ambient
\(\hypAmbientWtCpl\), and a Virasoro BRST contraction
\(h_{\mathrm{Vir}}\).

The proposition is now titled
`Quaternary Virasoro A-infinity operation in the Stasheff gauge;
\ClaimStatusProvedHere; licensing
\(\alpha+\gamma+\varepsilon\) via consecutive collision coordinates,
\(\hypAmbientWtCpl\), and the Virasoro BRST contraction`.  Its formula
is written as
\[
m_4^{h_{\mathrm{Vir}}}(T,T,T,T)
=-h_{\mathrm{Vir}}\mathrm{Obs}_4(T,T,T,T),
\]
with
\[
m_4^{h_{\mathrm{Vir}}}
=
(\lambda_{12}-\lambda_{34})(\lambda_{12}-\lambda_{23}+\lambda_{34})
\left[
2\partial T
4(\lambda_{12}+\lambda_{23}+\lambda_{34})T
\right.
\]
\[
\left.
\hspace{7em}
\frac{c}{6}(\lambda_{12}+\lambda_{23}+\lambda_{34})^3
\right].
\]
The proof now uses the five arity-four Stasheff source terms, not the
associator equation as a false stand-in for the Stasheff relation.

A349 also repairs `compute/lib/symbolic_stasheff.py`: its docstrings
now say the output is a Stasheff-gauge representative after choosing
the Virasoro SDR/BRST contraction.  The stale `no homotopy needed`
language is removed.  The new guard
`compute/tests/test_gravity_m4_stasheff_gauge.py` checks the theorem
title, status, licensing, gauge language, source-tree citations, and
the symbolic coefficients
\[
[\partial T]=2C,\qquad [T]=4C\Sigma,\qquad [1]=cC\Sigma^3/6,
\]
where \(C=(\ell_1-\ell_3)(\ell_1-\ell_2+\ell_3)\) and
\(\Sigma=\ell_1+\ell_2+\ell_3\).

A349 verification: Vol II py-compile passes on the new guard and
`compute/lib/symbolic_stasheff.py`.  The focused guard reports
`2 passed`; the widened Virasoro formula bundle reports
`117 passed, 22 subtests passed` across
`compute/tests/test_gravity_m4_stasheff_gauge.py`,
`compute/tests/test_field_sector_generating_function.py`, and
`compute/tests/test_swiss_cheese_virasoro_wheels.py`.  Propagation
scans find no live copy of the retired `no homotopy needed`, bare
proposition title, or associator-proof phrase outside negative guard
assertions.  `git diff --check` passes.  `make verify-licensing` exits
0 with zero blockers and 119 warning-class missing-tag lines elsewhere;
the first gravity-climax warning is now `prop:gravity-m5`.  Vol II
`make fast` converges after two passes with 2491 pages, zero undefined
citations, zero undefined references, and zero rerun requests.

A348 repairs `prop:formality-depth-discriminant` and the surrounding
gravity-climax discriminant prose.  The theorem is now
`SC-formality, support depth, and scalar mixed-shell discriminant;
licensing \(\gamma+\varepsilon\) via
\(\hypAmbientWtCpl+\effKoszul\)`.  It defines the scalar-lane
projection
\[
\Delta_{\mathrm{sc}}(\cA)=
8\,\kappaChHodge(\cA)\,S^{\mathrm{sc}}_4(\cA),\qquad
S^{\mathrm{sc}}_4(\cA)=\pi_{\mathrm{sc}}(\Theta_{\cA,4})/x^4.
\]
The support-depth invariant is \(r^\Theta_{\max}\) read from
\(\Theta_\cA\).  The scalar discriminant may vanish or be nonzero on a
finite contact class, so it is not a finite-versus-infinite classifier.

A348 computes the non-degenerate Virasoro scalar lane as
\[
\Delta_{\mathrm{sc}}(\mathrm{Vir}_c)
=8\cdot(c/2)\cdot 10/[c(5c+22)]
=40/(5c+22),\qquad c\ne0,-22/5.
\]
The \(c=0\) scalar Hessian degeneration and the \(c=-22/5\) Kac
denominator degeneration are now named explicitly.  Infinite support
depth remains the non-terminal Virasoro wheel in
\(\hypAmbientWtCpl\), not the scalar fraction itself.

A348 propagates this correction through the `3d_gravity.tex` chapter
opening, physics interpretation, and scrambling dichotomy; through Vol
II `examples-worked.tex` and `ht_bulk_boundary_line.tex`; through Vol I
`arithmetic_shadows.tex`, `introduction.tex`, the preface surfaces, and
the high-level single-line Riccati summaries in
`higher_genus_modular_koszul.tex`; and through Vol III
`e2_chiral_algebras.tex` and `quantum_chiral_algebras.tex`.  The Vol I
Hitchin comparison now uses the corrected normalization
\(\Delta_{\mathrm{sc}}=40/(5c+22)\) and
\(\Delta_H=-8c^2\Delta_{\mathrm{sc}}\).

A348 extends `compute/tests/test_gravity_support_depth_classification.py`
to guard the chapter-opening scalar-lane language, the licensed
proposition, the \(40/(5c+22)\) computation, the degeneracy exclusions,
and the finite-contact versus mixed-shell scrambling split.  The guard
rejects the retired `Formality, contact depth, and scalar
discriminant`, `single number attached to any chiral algebra`, and
`\Delta`-classifier language.

A348 verification: Vol II py-compile passes on the touched guard.  The
focused guard reports `5 passed`; the widened gravity bundle reports
`34 passed` across
`compute/tests/test_gravity_support_depth_classification.py`,
`compute/tests/test_climax_theorems_wave18_iv.py`, and
`compute/tests/test_heisenberg_scalar_linearity.py`.  Active Vol II
classifier scans are clean outside negative guard assertions; exact
Vol I/Vol III scans find no surviving copies of the patched
`shadow discriminant controls whether` / `forces infinite tower`
claims.  `git diff --check` passes on the touched Vol II, Vol I, and
Vol III files.  `make verify-licensing` exits 0 with zero blockers and
120 warning-class missing-tag lines elsewhere.  Vol II `make fast`
converges after two passes with 2491 pages, zero undefined citations,
zero undefined references, and zero rerun requests.

A347 repairs the next gravity-climax theorem surface:
`cor:gauge-gravity-dichotomy` in `chapters/connections/3d_gravity.tex`.
The corollary no longer claims that the whole standard HT landscape is
organised by \((\Ainf,\Delta)\)-complexity data.  It is now the
gauge-gravity-matter support-depth trichotomy, licensed by
\(\gamma+\varepsilon\) via \(\hypAmbientWtCpl+\effKoszul\), and scoped
to the effective Koszul support branch plus the completed ambient for
the class-\(\mathbf M\) wheel.

A347 states the invariant triple explicitly:
\[
\bigl(r^\Theta_{\max},\
\Delta_z|_{\mathrm{free\ generators}},\
\Delta_z|_{\mathrm{closed\ composites}}\bigr).
\]
The table now uses support depth \(r^\Theta_{\max}\).  Gauge theory has
finite Lie transfer; pure gravity and BP matter-coupled gravity share
the non-terminal Virasoro support wheel; BP differs by the
non-primitive composite coproduct coming from the surviving weight-one
current.

A347 also removes residual `The shadow depth is \(r_{\max}\)` claims
from the gravity primitive-package surface, the Vol II FM-calculus
extension file, and the Vol III local \(\mathbb P^2\) theorem row.  The
replacement is \(r^\Theta_{\max}=\infty\), with the composite-closed
support packet named as the non-terminating object.

A347 extends `compute/tests/test_gravity_support_depth_classification.py`
to guard the licensed corollary, the invariant triple, the
non-terminal support wheel, and the Virasoro \(r^\Theta_{\max}\)
wording.  It rejects the retired whole-landscape trichotomy and
`shadow depth \(r_{\max}\)` phrases.

A347 verification: Vol II py-compile passes on the touched guard.  The
focused Vol II bundle
`compute/tests/test_gravity_support_depth_classification.py
compute/tests/test_climax_theorems_wave18_iv.py
compute/tests/test_heisenberg_scalar_linearity.py` reports `33 passed`.
Fixed-string scans of Vol II, Vol I, and Vol III surfaces find no
surviving copy of the retired corollary/shadow-depth phrases.  `make
verify-licensing` exits 0 with zero blockers and 121 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2491 pages, zero undefined citations, zero undefined
references, and zero rerun requests; final-log scan finds no fatal TeX
errors, undefined controls, or unresolved citation/reference warnings.

A347 residual: the next gravity-climax warning is
`prop:formality-depth-discriminant`.  It already has partial
support-depth wording but still lacks statement-level licensing.

A346 propagates the A345 support-packet repair into the gravity climax.
The introduction no longer says that the \(E_1\) chiral coalgebra is
classified by maximal OPE pole order.  Its section is now the
support-packet principle:
\[
\Theta_\cA=\sum_{r\ge2}\Theta_{\cA,r},\qquad
r^\Theta_{\max}(\cA)=\sup\{r\ge2\mid\Theta_{\cA,r}\ne0\},
\]
after normal-ordered composite closure.  OPE pole order is input data,
not the invariant.

A346 repairs `chapters/connections/3d_gravity.tex` at
`prop:pole-order-classification`, preserving the label for reference
stability.  The proposition is now `Composite-closed support-depth
classification`, with licensing
\(\gamma+\varepsilon\) via \(\hypAmbientWtCpl+\effKoszul\).  It
classifies standard-atlas chirally Koszul boundary algebras by
\(r^\Theta_{\max}\): binary central class \(\mathbf G\),
Lie--Jacobi cubic class \(\mathbf L\), finite contact class
\(\mathbf C\), and the non-terminal Virasoro wheel in class
\(\mathbf M\).

A346 also repairs downstream gravity prose.  The quartic OPE channel is
now only the first visible Virasoro witness; the classified datum is the
non-terminating composite-closed support packet.  The independent
verification metadata now says FBZ/Kac/DS--Kac verify the OPE inputs,
not the support-packet classification itself.

A346 propagates the same correction into Vol I:
`../chiral-bar-cobar/appendices/ordered_associative_chiral_kd.tex`,
`../chiral-bar-cobar/chapters/theory/ordered_associative_chiral_kd.tex`,
and the O3 helper test.  The Vol I remarks now distinguish the binary
Heisenberg support packet from the finite \(\beta\gamma\) contact
envelope; the O3 test advertises only a pole-threshold input check.

A346 adds `compute/tests/test_gravity_support_depth_classification.py`.
The guard requires the support-packet principle, \(\Theta_\cA\),
\(r^\Theta_{\max}\), \(\effKoszul\), and \(\hypAmbientWtCpl\), and
rejects the old raw-pole theorem/principle wording.

A346 verification: Vol II py-compile passes on the touched tests.  The
focused Vol II bundle
`compute/tests/test_gravity_support_depth_classification.py
compute/tests/test_climax_theorems_wave18_iv.py
compute/tests/test_heisenberg_scalar_linearity.py` reports `32 passed`.
The propagated Vol I O3 suite reports `61 passed`.  Fixed-string scans
of active Vol II, Vol I, and Vol III surfaces find no surviving live
copy of the retired pole-order-classification/principle phrases; the
remaining hit is a historical audit note.  `make verify-licensing`
exits 0 with zero blockers and 122 warning-class missing-tag lines
elsewhere.  Vol II `make fast` converges after two passes with 2489
pages, zero undefined citations, zero undefined references, and zero
rerun requests; final-log scan finds no fatal TeX errors, undefined
controls, or unresolved citation/reference warnings.

A346 residual: older standard-family surfaces still sometimes use
`shadow depth` as a shorthand for the Koszul support-depth convention.
Future theorem statements should prefer \(r^\Theta_{\max}\) whenever
there is any risk of confusing support depth with generator pole order.

A345 repairs the Rosetta pole-structure theorem in
`chapters/examples/rosetta_stone.tex`.  The old theorem classified
standard families by the maximal pole order of a chosen OPE generator
span.  That was not an invariant of the closed chiral algebra: the
simple-pole systems repaired in A344 have
\[
m^{\mathrm{gen}}_j=0\qquad (j\ge3)
\]
on the chosen generators, while their normal-ordered composite closure
carries a quartic class-\(\mathbf C\) support component.

A345 replaces the theorem by
`Composite-closed support-depth stratification`, preserving the old
label `thm:pole-structure-dichotomy` for reference stability.  The
classification invariant is now
\[
r^\Theta_{\max}(\cA)
=\sup\{r\ge2\mid\Theta_{\cA,r}\ne0\},
\]
computed from the canonical bar-intrinsic support packet after
composite closure, not from raw generator pole order.

A345 also repairs the atlas convention.  The local section is now the
support-packet atlas, and it explicitly records
\[
\operatorname{Sh}_r(\cA)=\pi_{\mathrm{sc}}(\Theta_{\cA,r}),
\]
so scalar shadows are projections of the composite-closed packet.
Adjacent W-algebra atlas references, the lattice VOA compute convention
note, and the degenerate Virasoro example were updated to the same
language.

A345 extends `compute/tests/test_heisenberg_scalar_linearity.py` to
guard the repaired theorem and atlas bridge.  The guard requires the
composite-closed support package, \(r^\Theta_{\max}\), the
simple-pole/composite distinction, and the support-packet atlas
definition; it rejects the retired raw-pole theorem title, maximal-pole
hypothesis, and old shadow-obstruction atlas wording.

A345 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched compute/test files.  The focused pytest bundle
`compute/tests/test_heisenberg_scalar_linearity.py
compute/tests/test_shadow_depth_atlas.py` reports `13 passed`.
Fixed-string scans of active Vol II and compute surfaces find no
surviving copy of the retired raw-pole theorem phrases or old
`shadow obstruction tower atlas` wording outside the negative guard.
`make verify-licensing` exits 0 with zero blockers and 123
warning-class missing-tag lines elsewhere.  Vol II `make fast`
converges after two passes with 2489 pages, zero undefined citations,
zero undefined references, and zero rerun requests; direct final-log
scan finds no fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.

A345 residual: the broader \(r_{\mathrm{sh}}\) shorthand remains in the
Vol I/Vol II introduction and relative-Feynman-transform surfaces as
the Koszul shadow-depth convention.  The raw-pole classification and
old atlas-name drift are removed from the active Rosetta support
surface.

A344 repairs the simple-pole Rosetta examples in
`chapters/examples/rosetta_stone.tex`.  The \(\beta\gamma\),
free-fermion, \(bc\), and symplectic-boson sections no longer treat a
generator-span calculation as the full closed chiral algebra.

A344 fixes the mathematical content.  On the chosen generator span, the
simple-pole OPE gives only the binary residue and the transferred
support packet has no primitive operation above arity \(2\):
\[
m^{\mathrm{gen}}_j=0\qquad (j\ge3).
\]
This statement is explicitly not a classification of the algebra after
normal-ordered composite closure.

A344 scopes the composite algebra.  The \(\beta\gamma\) ghost current
\({:}\beta\gamma{:}\), the free-fermion current
\({:}\psi\psi{:}\), the \(bc\) ghost current \({:}bc{:}\), and the
symplectic-boson current \({:}\chi^+\chi^-{:}\) produce the quartic
composite support layer.  The closed composite algebra is therefore
class \(\mathbf C\) with
\[
r^\Theta_{\max}=4,
\]
while the generator collision-depth spectrum remains \(\{0\}\).

A344 extends `compute/tests/test_heisenberg_scalar_linearity.py` to
guard the whole \(\beta\gamma\)-through-symplectic-boson block.  The
guard requires \(m^{\mathrm{gen}}_j=0\), composite
\(r^\Theta_{\max}=4\), and the explicit warning that generator-span
vanishing does not classify the closed composite algebra.  It rejects
the retired `All m_k=0`, `depth 4`, bare `r_max=4`, and old
generator-depth headings.

A344 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched guard.  The focused Heisenberg/Rosetta scalar-linearity
suite reports `7 passed`.  Fixed-string scans of `rosetta_stone.tex`
find no surviving exact copy of the repaired stale forms.  `make
verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2487 pages, zero undefined citations, zero undefined
references, and zero rerun requests; direct final-log scan finds no
fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.

A344 residual: the older pole-structure atlas still uses global
shadow-depth terminology and needs a separate theorem-level
support-depth pass.

A343 repairs the later Heisenberg Rosetta computations in
`chapters/examples/rosetta_stone.tex`.  After A342 fixed the opening,
the `heisenberg-hydrogen-atom` and `heisenberg-e1-ordered-shadow`
computations still collapsed support degree, scalar Hodge data, and
transport into a single scalar \(k\).  They also contained a formula
contradiction: the text had been corrected to \(\mathcal R(z)=z^k\),
but the same computation still expanded \(e^{k\hbar/z}\) as a Taylor
series in \(1/z\).

A343 fixes the support statement.  The level \(k\) now controls only
the rank-one scalar/OPE data.  The canonical support packet has no
primitive operation above arity \(2\); all higher support operations
vanish as \(m^H_j=0\) for \(j\ge3\); and
\[
r^\Theta_{\max}(\cH_k)=2
\]
is a support degree.  This is kept separate from the coalgebra
statement that \(\barBch(\cH_k)\) is cogenerated in degree \(1\).

A343 fixes the transport formula.  The ordered-shadow computation now
uses logarithmic transport on a chosen branch:
\[
\mathcal R_\hbar(z;z_0)
=\exp(k\hbar\,\xi)
=\sum_{n=0}^{50}\frac{(k\hbar)^n}{n!}\xi^n+O(\xi^{51}),
\qquad
\xi=\log(z/z_0).
\]
The expansion is entire in the logarithmic coordinate and multivalued
as a function of \(z\), with monodromy \(\exp(2\pi i k\hbar)\).  It is
not a Taylor series in \(1/z\).

A343 also scopes the genus and partition claims.  The genus coefficient
is now the scalar Hodge coefficient
\[
F_g^{\mathrm{sc}}(\cH_k)=k\lambda_g^{\mathrm{FP}},
\]
and the eta/partition statement is a determinant-line shadow of the
Heisenberg bar complex, not an unqualified assertion that the full bar
complex determines all partition data.

A343 extends `compute/tests/test_heisenberg_scalar_linearity.py` to
isolate the two Rosetta computations.  The guard requires
\(r^\Theta_{\max}\), \(m^H_j=0\), the degree-one cogenerator statement,
\(F_g^{\mathrm{sc}}\), determinant-line shadow language, and the
logarithmic transport expansion.  It rejects the retired single-scalar,
shadow-tower, degree-two-cogenerator, \(e^{1/z}\), \(1/z\)-Taylor, and
unscoped \(F_g\) phrases.

A343 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched guard.  The focused Heisenberg scalar-linearity suite
reports `6 passed`.  Fixed-string scans of `rosetta_stone.tex` find no
surviving exact copy of the repaired local stale phrases.  `make
verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2487 pages, zero undefined citations, zero undefined
references, and zero rerun requests; direct final-log scan finds no
fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.

A343 residual: later Rosetta examples still contain generator-level
phrases of the form `All $m_k=0$ for $k\ge3$ on generators`, and the
older pole-structure atlas still uses shadow-depth vocabulary.  Those
require a separate pass.

A342 repairs the live Rosetta opening in
`chapters/examples/rosetta_stone.tex`.  The opening no longer says
that Heisenberg simply "has shadow depth \(r_{\max}=2\)", that
\(\Theta^{\mathrm{oc}}\) terminates as a primitive fact, that
\(F_g=k\lambda_g^{\mathrm{FP}}\) is the full genus tower, that the
line category is literally \(\cH_{-k}\)-modules, or that the derived
centre is the physical free-boson bulk.

A342 fixes the mathematical content.  The opening now treats
\(\cH_k\) as the Gaussian support-depth test case:
\[
r^\Theta_{\max}(\cH_k)=2
\]
on the canonical bar-intrinsic branch, so the support packet of
\(\Theta^{\mathrm{oc}}\) has no primitive component above degree \(2\).
The closed-form data are explicitly rank-one scalar data:
\[
\kappaChHodge(\cH_k)=k,\qquad
r_{\cH}(z)=k\Omega_{\cH}/z,\qquad
F_g^{\mathrm{sc}}(\cH_k)=k\lambda_g^{\mathrm{FP}}\quad(g\ge1).
\]

A342 also separates the comparison data.  The completed abelian line
comparison models the relevant open-colour line sector by
\(\cH_{-k}\)-modules after the chosen comparison; this is not the
chiral Koszul dual \(\cH_k^!=\Sym^{\mathrm{ch}}(V^*)\).  The derived
centre is the algebraic bulk candidate, and comparison with the
physical free-boson bulk requires the HT bulk--Hochschild comparison
map.  The scalar \(k\) determines only the rank-one abelian scalar/OPE
data; line and bulk readings require their comparison maps.

A342 extends `compute/tests/test_heisenberg_scalar_linearity.py` to
guard the opening.  The test requires the repaired support-depth,
scalar-lane, line-comparison, and bulk-comparison language, and rejects
the old primitive-equality phrases `has shadow depth`, `terminates at
degree 2`, `genus tower $F_g = k`, `The line category is`, `the
derived center is the free boson bulk`, and `The single scalar~$k$
determines the Heisenberg package`.

A342 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched guard.  The focused Heisenberg scalar-linearity suite
reports `4 passed`.  Fixed-string scans of `rosetta_stone.tex` find no
surviving exact copy of the stale opening phrases.  `make
verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2487 pages, zero undefined citations, zero undefined
references, and zero rerun requests; direct final-log scan finds no
fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.

A342 residual: this repairs the opening only.  Later Rosetta atlas
passages still contain older support-depth/shadow-tower vocabulary and
need a separate theorem-by-theorem propagation pass.

A341 repairs the celestial soft hierarchy surfaces.  The active
`chapters/connections/thqg_soft_graviton_theorems.tex`, the extension
`chapters/connections/thqg_celestial_holography_extensions.tex`, and
the engine `compute/lib/celestial_ope_from_shadow.py` no longer treat
support depth as a count of independent celestial soft factors, Ward
identities, or OPE coefficients.

A341 fixes the executable semantics.  `ShadowDepthClass` now carries
`support_channel_bound`, the algebraic arity cutoff
\(r^\Theta_{\max}\), not the old \(r_{\max}-1\) count.
`ward_identity_support_profile` reports that the physical Ward count
requires a celestial comparison datum.  The support-channel table now
records arity bounds and marks the physical count as comparison
dependent.

A341 fixes the active theorem text.  The degree-truncation table in
`thqg_soft_graviton_theorems.tex` now lists support channels:
arity \(2\), arities \(2,3\), arities \(2,3,4\), and no finite
cutoff.  The text states that a comparison map sends algebraic channels
to soft operators, and that this map, not support depth alone, decides
which images are nonzero and independent.

A341 also repairs the class-M statement in
`thqg_celestial_holography_extensions.tex`.  Class M gives no finite
algebraic celestial support cutoff.  An infinite physical
soft-operator package follows only after a chosen celestial comparison
map detects infinitely many support channels.  Nonzero support or
obstruction data no longer directly imply nonzero physical soft
factors.

A341 extends `compute/tests/test_celestial_ope_from_shadow.py` to
require the support-channel bound \(r^\Theta_{\max}\), the
comparison-required Ward profile, the class-M no-finite-cutoff
statement, and the active soft-graviton support-channel table.  The
guard rejects the retired `r_max - 1` rule, `n_ward_identities`,
`n_celestial_factors`, `n_soft_factors`, `soft theorem tower is
infinite`, `contributes a nonzero soft factor`, `number of independent
OPE coefficients`, and the old `# cel. factors` table heading.

A341 verification: `compute/.venv/bin/python -m py_compile` passes on
the touched celestial engine and test.  The full celestial OPE test
file reports `83 passed`; the deletion-ledger suite reports
`62 passed`.  Fixed-string propagation scans over the touched
celestial theorem files and engine find the retired exact count and
physical-infinitude phrases only in negative guard assertions.
`make verify-licensing` exits 0 with zero blockers and 123
warning-class missing-tag lines elsewhere.  Vol II `make fast`
converges after two passes with 2487 pages, zero undefined citations,
zero undefined references, and zero rerun requests; direct final-log
scans find no fatal TeX errors, undefined controls, or unresolved
citation/reference warnings, with existing pdfTeX destination warnings
still present.

A341 residual: broader support-depth/shadow-tower terminology remains
in `examples-worked.tex`, `rosetta_stone.tex`, and older W-algebra
example surfaces.  The active celestial soft-factor count claim is
repaired here.

A340 repairs the live holomorphic-topological chapter in
`chapters/connections/holomorphic_topological.tex`.  The chapter no
longer uses bare \(r_{\max}\) to classify HT systems or boundary VOAs,
and it no longer claims that support depth counts protected BPS
multiplet types.

A340 fixes the mathematical content.  The MC theorem now says finite
projections \(\Theta_{\cA_T}^{\le r}\) form the support filtration of
the bar-intrinsic MC packet.  Obstruction classes test extension of
those truncations, but they do not determine the \(H^1\)-lift
coordinates.  The HT archetype corollary uses
\(r^\Theta_{\max}(\cA_T)\), the support depth of the canonical
bar-intrinsic MC packet.

A340 also rewrites the BPS paragraph.  Support depth admits a physical
reading only after a BPS comparison datum is fixed: a stability
condition, a central-charge map, and a comparison from support
coordinates of \(\Theta_{\cA_T}\) to protected sectors.  Support depth
bounds algebraic support channels available to such a comparison; it
does not count protected BPS multiplet types by itself, and the active
Kontsevich--Soibelman factors do not depend on support depth alone.

A340 propagates the correction through the standard HT examples,
boundary VOA table, twistor-anomaly remarks, celestial convolution
remark, Costello comparison, Gaiotto formality comparison, and Gaiotto
triality paragraph.  These surfaces now state finite or unbounded
support packets, not obstruction-tower termination or bare shadow-depth
classes.

A340 adds
`test_ht_support_depth_is_not_a_bps_multiplet_counter` in
`compute/tests/test_holomorphic_topological_scalar_trace_typing.py`.
The guard requires the support-depth title, \(r^\Theta_{\max}\)
coordinate, BPS comparison datum, standard-example support classes,
boundary-VOA support-depth table, formality support-depth paragraph,
residual-anomaly support-packet wording, and genus-zero support-packet
projections.  It rejects `r_{\max}`, `r_max`, `\rmax`, shadow-depth and
shadow-class wording, shadow obstruction tower, tower-termination
language, and the old BPS-multiplet counter claim.

A340 verification: `compute/.venv/bin/python -m py_compile` passes on
`compute/tests/test_holomorphic_topological_scalar_trace_typing.py`.
The focused A340 pytest reports `1 passed`; the full HT scalar-trace
guard file reports `5 passed`; the deletion-ledger suite reports
`62 passed`.  Targeted stale-phrase grep finds no surviving
\(r_{\max}\), shadow-depth/class, shadow-obstruction-tower,
tower-termination, or BPS-counter wording in
`holomorphic_topological.tex`.  `make verify-licensing` exits 0 with
zero blockers and 123 warning-class missing-tag lines elsewhere.  Vol
II `make fast` converges after two passes with 2487 pages, zero
undefined citations, zero undefined references, and zero rerun
requests; direct final-log scans find no fatal TeX errors, undefined
controls, or unresolved citation/reference warnings, with existing
pdfTeX destination warnings still present.

A340 residual: broad propagation search still finds older support-depth
language outside the HT target, especially in `examples-worked.tex`,
`rosetta_stone.tex`, older W-algebra/example files,
`thqg_celestial_holography_extensions.tex`,
`unified_chiral_quantum_group.tex`, and compute docstrings.

A339 repairs the proved-examples fingerprint surface in
`chapters/examples/examples-complete-proved.tex`.  The Chern--Simons
package and the fingerprint theorem no longer use \(r_{\max}\) or the
claim that the support-depth coordinate fixes the truncation of a
shadow obstruction tower and hence the homological length of \(B(A)\).

A339 fixes the mathematical content.  The CS package now has support
shadow depth \(r^\Theta_{\max}=3\): the bar-intrinsic support packet
has the scalar quadratic component and the cubic Lie--Jacobi component,
with no primitive support coordinate above degree \(3\).  This does not
say obstruction-class vanishing determines the packet, and it does not
delete the scalar Hodge genus tower.

A339 also rewrites the fingerprint coordinate as
\(r^\Theta_{\max}\), the support shadow depth of the canonical
bar-intrinsic MC packet on the standard branch.  In the fingerprint
proof this coordinate fixes only the support envelope
\(\Theta_{A,r}=0\) for \(r>r^\Theta_{\max}(A)\).  It is not an
obstruction-depth invariant and does not by itself determine the
homological length of \(B(A)\); the bar weight profile is fixed only
after the OPE-residue and strong-generator data are included.

A339 extends
`test_affine_generic_center_not_promoted_to_cs_path_integral` so the
guard requires the support-packet warning, the \(r^\Theta_{\max}\)
fingerprint coordinate, the obstruction-depth/homological-length
separation, and the OPE-residue/strong-generator qualification.  It
rejects `\rmax`, `r_{\max}`, `r_max`, `rmax`, the old
truncation/homological-length sentence, the old CS obstruction-tower
termination sentence, the old finite-degree MC termination wording, and
`infinite shadow depth`.

A339 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused examples and
symplectic guard pair reports `2 passed`; the full deletion-ledger
suite reports `62 passed`.  Fixed-string and regex scans find no
targeted legacy \(r_{\max}\), obstruction-length, or finite-degree
termination forms in `examples-complete-proved.tex`.
`make verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2487 pages, zero undefined citations, zero undefined
references, and zero rerun requests; final-log scans find no fatal TeX
errors, undefined controls, or unresolved citation/reference warnings.
Environment-balance counts match for proofs, theorems, propositions,
remarks, definitions, enumerates, and itemizes in the edited chapter.
Tracked scoped `git diff --check` passes on the touched chapter and
audit ledger; a direct trailing-whitespace scan passes on all touched
paths.

A339 residual: live examples and theory files still contain older
support-depth terminology outside the proved-examples target, especially
`examples-worked.tex`, `rosetta_stone.tex`,
`holomorphic_topological.tex`, and
`unified_chiral_quantum_group.tex`.

A338 repairs the live modular-bootstrap spectral-sequence surface in
`chapters/connections/thqg_modular_bootstrap.tex`.  The chapter no
longer uses shadow depth \(r_{\max}\) to determine collapse pages,
nonzero differentials, graph-degree bounds, or gravitational
complexity.  It now uses support shadow depth
\(r^\Theta_{\max}(\cA)\) as a support-degree envelope for the
bar-intrinsic support packet.

A338 fixes the key implication.  Infinite support depth says
\(\Theta_{\cA}^{\leq N}\ne\Theta_{\cA}\) for every finite \(N\); it
does not say \(o_{N+1}(\cA)\ne0\), and it does not imply \(d_N\ne0\)
after passage to the \(E_N\)-page.  Finite support depth gives no
support component above \(r^\Theta_{\max}\); obstruction vanishing
controls extendability of a truncation and does not erase a chosen
\(H^1\)-lift coordinate.

A338 rewrites the mixed-type proposition as a no-finite-support-cutoff
statement, rewrites the spectral-sequence theorem as a support-projection
theorem, and propagates the same discipline through the Gaussian scalar
degeneration corollary, Koszulness/bootstrap remark, gravitational
dictionary, non-renormalization proof, support-degree bound, and final
modular-bootstrap classification.  The gravitational comparison language
now says support-degree complexity, not metric-gravity graph-degree
classification.

A338 extends
`test_modular_bootstrap_outputs_scalar_shadow_not_gravity_partition` so
the guard requires the no-finite-support-cutoff proposition, the
support-projection theorem, the explicit \(o_{N+1}\)/\(d_N\) separation,
the support-depth dictionary row, the support-degree bound, and the
finite-page-degeneration refusal.  It rejects `r_{\max}`, `r_max`,
`rmax`, the old obstruction-tower and collapse-page slogans, the old
nonzero-every-\(d_r\) claim, and the old graph-sum/vertex-degree proof.

A338 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused modular-bootstrap
guard reports `1 passed`; the full deletion-ledger suite reports
`62 passed`.  Fixed-string and regex scans find no targeted legacy
depth/collapse/graph-degree forms in `thqg_modular_bootstrap.tex`.
`make verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  Vol II `make fast` converges after two
passes with 2487 pages, zero undefined citations, zero undefined
references, and zero rerun requests; final-log scans find no fatal TeX
errors, undefined controls, or unresolved citation/reference warnings.
Environment-balance counts match for proofs, theorems, propositions,
remarks, corollaries, enumerates, and itemizes.  Tracked scoped
`git diff --check` passes on the touched chapter and audit ledger; a
direct trailing-whitespace scan passes on all touched paths.

A338 residual: support-depth terminology and older \(r_{\max}\)-style
fingerprints remain outside the modular-bootstrap target, especially in
the live example/theory layers `examples-complete-proved.tex`,
`examples-worked.tex`, `rosetta_stone.tex`,
`holomorphic_topological.tex`, and `unified_chiral_quantum_group.tex`.

A337 extends the support-depth repair into
`chapters/connections/thqg_perturbative_finiteness.tex`.  The chapter
no longer treats finite shadow depth or obstruction stabilization as a
bound on degree contributions.  Its degree-summed shadow series,
genus-\(g\) amplitude decomposition, examples, tables, analytic
extension loci, convergence theorems, and holographic dictionary now
use support shadow depth \(r^\Theta_{\max}(\cA)\).

A337 fixes the mathematical implication rather than only its notation.
Finite support depth means the selected bar-intrinsic support packet has
no nonzero component above \(r^\Theta_{\max}\).  Eventual vanishing of
obstruction classes records \(H^2\)-extendability; it does not kill the
corresponding \(H^1\)-lift coordinate.  The infinite case is now a
completed degree product, while scalar or normed numerical convergence
is kept as a separate analytic hypothesis.

A337 extends
`test_thqg_degree_summed_series_not_full_gravity_partition_function` so
the guard requires finite/infinite support shadow depth, the
bar-intrinsic support packet, the explicit \(H^1\)-lift-coordinate
warning, and the corrected holographic dictionary row.  It rejects the
old finite/infinite shadow-depth phrases, the graph-bound slogan, the
obstruction-stabilization slogan, and any remaining `r_{\max}`,
`r_max`, `rmax`, \(d_\infty\), or \(f_\infty\) in the chapter.

A337 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused guard reports
`1 passed`; the widened regression bundle reports `68 passed` across
the deletion-ledger, universal-holography, and scalar-trace/full-energy
guards.  Fixed-string and regex scans find no remaining legacy depth
notation or old graph-bound slogan in
`thqg_perturbative_finiteness.tex`.  `make verify-licensing` exits 0
with zero blockers and 123 warning-class missing-tag lines elsewhere.
Vol II `make fast` converges after two passes with 2487 pages, zero
undefined citations, zero undefined references, and zero rerun
requests; final-log scans find no fatal TeX errors, undefined controls,
or unresolved citation/reference warnings.  Tracked scoped
`git diff --check` passes on the touched chapter and audit ledger; a
direct trailing-whitespace scan passes on all four touched paths.

A337 residual: broad searches still find other Vol II and cross-volume
"shadow depth" and "obstruction tower terminates" variants outside the
perturbative-finiteness target, including gravitational-Yangian,
FM-calculus extension, examples, 3d-gravity, and compute-documentation
surfaces.

A336 clears the later two-chapter legacy \(r_{\max}\) surface that A335
left as residual.  In `thqg_gravitational_complexity.tex` and
`thqg_soft_graviton_theorems.tex`, every literal \(r_{\max}\),
`r_max`, `rmax`, \(d_\infty\), and \(f_\infty\) occurrence has been
removed.  The later tables, examples, soft algebra, soft OPE,
open-problem, and categorical-conjecture surfaces now use
\(r^\Theta_{\max}\) or the local abbreviation
\(r_\Theta=r^\Theta_{\max}(\cA)\).

A336 fixes two mathematical overclaims, not only notation.  First, the
subleading and sub-subleading soft statements now assume the actual
nonzero support components
\[
\Theta_{\cA,3}=\mathfrak C(\cA)\ne0,\qquad
\Theta_{\cA,4}=\mathfrak Q(\cA),
\]
rather than treating \(r^\Theta_{\max}\ge3,4\) as a proxy.  Second, the
deformation-stability theorem is now a finite-window algebraic
constancy statement for support coordinates on a fixed canonical
standard-family gauge slice.  It explicitly makes no monotonicity claim
without control of the \(H^1\)-lift coordinates.

A336 also scopes the G2 main theorem to the canonical bar-intrinsic
branches of the standard gravitational landscape.  The statement now
says that no classification of arbitrary components of the filtered MC
moduli problem is asserted.  The open gap question, Steinberg
stratification conjecture, categorical detection problem, and
\(t\)-structure amplitude conjecture now lift \(r^\Theta_{\max}\), not
the obsolete obstruction-depth proxy.

A336 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused support-depth
guard reports `1 passed`.  The widened regression bundle reports
`143 passed` across the deletion-ledger, shadow-depth, holography, and
climax IV gates.  Fixed-string/regex scans find no remaining
`r_{\max}`, `r_max`, `rmax`, \(d_\infty\), or \(f_\infty\) in the two
target chapters.  `make verify-licensing` exits 0 with zero blockers
and 123 warning-class missing-tag lines elsewhere.  Vol II `make fast`
converges after two passes with 2487 pages, zero undefined citations,
zero undefined references, and zero rerun requests; final-log scans
find no fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.  Scoped `git diff --check` passes.

A336 residual: this clears the two-chapter legacy \(r_{\max}\) surface
identified by A335.  Broader cross-volume searches for conceptually
equivalent but differently notated obstruction-depth/support-depth
conflations remain future work.

A335 repairs the central \(r_{\max}\) calculus in the gravitational
complexity and soft-graviton surfaces.  The gravitational complexity is
now
\[
\rho_{\mathrm{grav}}(\cA)=r^\Theta_{\max}(\cA)
=\inf\{r\ge2:\Theta_{\cA,n}=0\ \text{for all }n>r\},
\]
while
\[
d_{\mathrm{obs}}(\cA)
=\inf\{r\ge2:o_{r'+1}(\cA)=0\ \text{for every }r'\ge r\}
\]
records only eventual \(H^2\)-vanishing.  The text now says explicitly
that \(d_{\mathrm{obs}}\) does not set the higher
\(H^1(J^n(\cA),d_2)\)-lift coordinates to zero.

A335 rewrites the operadic-complexity comparison as
\[
r^\Theta_{\max}(\cA)=d^\Theta_\infty(\cA)=f^\Theta_\infty(\cA),
\]
with \(d_{\mathrm{obs}}\) a separate invariant.  The quasi-isomorphism
statement is correspondingly filtered and gauge-component scoped.

A335 replaces the obstruction-only four-class theorem by a
standard-family support-depth classification.  The theorem now applies
to the canonical bar-intrinsic branches of the standard landscape:
Gaussian \(r^\Theta_{\max}=2\), Lie \(r^\Theta_{\max}=3\), contact
\(r^\Theta_{\max}=4\), and mixed \(r^\Theta_{\max}=\infty\).  The
finite rows include the vanishing of the higher \(H^1\)-lift
coordinates; the theorem explicitly does not classify arbitrary
components of the full filtered MC moduli problem.

A335 also repairs the soft hierarchy.  The chapter now defines
\(r_\Theta(\cA)=r^\Theta_{\max}(\cA)\) as the support soft depth and
uses it in the degree decomposition and Ward theorem.  Finite Ward
systems are determined by the finite support packet
\(\{\Sh_2,\ldots,\Sh_{r_\Theta}\}\) together with vanishing higher lift
coordinates, not by obstruction vanishing alone.  The later
\(L_\infty\) remarks now say support depth, not formality-obstruction
depth.

A335 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused support-depth
guard reports `1 passed`.  The widened regression bundle reports
`143 passed` across the deletion-ledger, shadow-depth, holography, and
climax IV gates.  Fixed-string scans find no surviving exact copy of
the retired local definitions/slogans in the two edited chapters, and
no duplicate adjacent proof environment.  Scoped `git diff --check`
passes.  `make verify-licensing` exits 0 with zero blockers and 123
warning-class missing-tag lines elsewhere.  Vol II `make fast`
converges after two passes with 2487 pages, zero undefined citations,
zero undefined references, and zero rerun requests; final-log scans
find no fatal TeX errors, undefined controls, or unresolved
citation/reference warnings.

A335 residual: a broad fixed-string scan still finds 55 occurrences of
legacy \(r_{\max}\) notation in later soft-graviton and gravitational
complexity sections.  Those later tables, examples, deformation loci,
and categorical detection paragraphs remain the next propagation
surface.

A334 repairs the standard-family reconstruction rows after the A333
support-depth split.  The primary reconstruction chapter and the
celestial copy still used the old \(r_{\max}\) notation on the
Heisenberg, affine, \(\beta\gamma\), and Virasoro rows.  The repaired
statements now use
\[
r^\Theta_{\max}(\mathcal H_k)=2,\qquad
r^\Theta_{\max}(\widehat{\mathfrak g}_k)=3,\qquad
r^\Theta_{\max}(V_{\beta\gamma})=4,\qquad
r^\Theta_{\max}(\mathrm{Vir}_c)=\infty.
\]
For the finite rows, the statements include the canonical gauge-slice
condition that the higher \(H^1\)-lift coordinates vanish.  Thus the
finite-depth claim is a support-packet statement, not an
obstruction-depth statement.

A334 keeps the affine normalization fixed.  The local compute
convention and tests say
\[
\kappaChHodge(\widehat{\mathfrak g}_k)
=\frac{\dim(\mathfrak g)(k+h^\vee)}{2h^\vee},
\]
the shifted modular characteristic, not the Sugawara central charge
\(c/2\).  The celestial affine corollary now says the scalar shadow
series is determined by the level, \(\dim(\mathfrak g)\), and the
structure constants; it no longer says `rank \(\dim(\mathfrak g)\)'.

A334 also replaces scalar modular partition-function language on this
surface by scalar modular shadow-series language.  The \(\beta\gamma\)
applied row no longer asserts an unjustified four-number full
partition-function theorem; it uses the finite contact packet.

A334 propagation replaces exact Virasoro echoes
\(r_{\max}(\mathrm{Vir}_c)=\infty\) by
\(r^\Theta_{\max}(\mathrm{Vir}_c)=\infty\) in downstream Vol II
summaries and in the Vol I echo
`~/chiral-bar-cobar/chapters/theory/higher_genus_modular_koszul.tex`.
A broader scan still finds a large legacy \(r_{\max}\) calculus in
soft-graviton and gravitational-complexity chapters; that is the next
support-depth propagation surface, not completed by A334.

A334 verification: `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused guard reports
`1 passed`.  The widened regression bundle reports `242 passed` across
the deletion ledger, celestial, shadow-depth, holography, and class-M
gate tests.  Fixed-string scans find the retired scalar-partition and
class-row formulas only inside negative guard assertions.  Scoped
`git diff --check` passes on all touched Vol II files and separately on
the touched Vol I file.  `make verify-licensing` exits 0 with zero
blockers and 123 warning-class missing-tag lines elsewhere.  Vol II
`make fast` converges after two passes with 2487 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Vol I
`make fast` reaches a stable warning state after three passes with zero
undefined citations and 1201 persistent undefined references on that
dirty tree; its log still contains pre-existing TeX errors/undefined
controls away from the touched propagation line.

A333 repairs the finite-shadow-depth criterion exposed by A332.  The
old shadow-depth theorem still treated eventual vanishing of the
\(H^2\)-obstruction tower as if it forced actual truncation of the
filtered MC element.  That is false in the obstruction calculus:
\(H^2\)-classes decide whether the next lift exists, while the set of
lifts is an \(H^1(J^n(\cA),d_2)\)-torsor.

The live theorem is now
`Obstruction depth and support depth of the shadow tower;
\ClaimStatusConditional; licensing tags \(\gamma+\epsilon\)`.  It
defines the obstruction depth
\[
d_{\mathrm{obs}}(\cA)=
\inf\{r\ge2: o_{r'+1}(\cA)=0\ \text{for all } r'\ge r\}
\]
and the support depth
\[
r^\Theta_{\max}(\cA)=
\inf\{r\ge2:\Theta_{\cA,n}=0\ \text{for all } n>r\}.
\]
The theorem proves
\[
d_{\mathrm{obs}}(\cA)\le r^\Theta_{\max}(\cA),
\]
and states the corrected consequence: finite reconstruction is
governed by \(r^\Theta_{\max}\), not by obstruction depth alone.
Vanishing of \(o_n(\cA)\) above a degree says that later extensions are
unobstructed; it does not set the later \(H^1\)-lift coordinates to
zero.

A333 propagates this distinction through the support-shadow-depth
definition, the landscape table, the Koszulness comparison, the
operadic complexity theorem, the Postnikov \(k\)-invariant discussion,
and the deformation-theoretic interpretation.  The \(k\)-invariants
\(k_r=o_{r+1}\) now determine the obstruction-depth profile only.  The
bar-intrinsic MC gauge class still needs the \(H^1\)-lift coordinates
recorded in the shadow packet.

A333 also repairs the celestial finite-jet surface.  The finite
\(q=0\) reconstruction theorem now assumes finite support shadow depth:
above \(r^\Theta_{\max}\), both obstruction classes \(o_n(\cA)\) and
lift coordinates \(\lambda_n\) vanish in the chosen gauge slice.

A333 verification: fixed-string propagation scans over active Vol II
chapters, appendices, compute tests, Vol I chapters, and Vol III
chapters find the retired exact finite-depth formulas only in negative
guard assertions.  `python3 -m py_compile` passes on
`compute/tests/test_deletion_ledger.py`.  The focused shadow-depth
guard reports `1 passed`.  The widened holography/shadow-depth bundle
reports `265 passed` across `test_deletion_ledger.py`,
`test_climax_theorems_wave12_iv.py`,
`test_universal_holography_functor.py`,
`test_programme_climax_holography_gates.py`,
`test_part_vi_platonic_introduction.py`,
`test_chd_class_m_holography_gates.py`,
`test_holographic_ht_engine.py`, `test_shadow_depth_atlas.py`, and
`test_swiss_cheese_virasoro_wheels.py`.

A328 repairs the live Universal Holography functor chapter itself.  The
chapter still stated a bare functor
\[
\Phi_{\mathrm{hol}}\colon
\ChirAlg^{\omega,\mathrm{BL,adm}}_X\to
\HTQFT_{X\times\mathbb R}
\]
and its proof still constructed a physical HT bulk from
\(\mathcal A\) alone by Costello--Gwilliam plus Hochschild.  Move 4
still claimed that chain-level equalities compose to equality of
class-\(\mathbf M\) centres and bulk.  This contradicted A326--A327,
where Universal Holography had already become a theorem on supplied
HT realization data.

The live functor theorem now has the \(\Xi\)-decorated arity
\[
\Phi_{\mathrm{hol}}\colon
\ChirAlg^{\omega,\mathrm{BL,adm},\Xi}_X
\longrightarrow
\HTQFT^\Xi_{X\times\mathbb R}.
\]
An object is
\((\mathcal A,\omega,T_\omega,\Xi_\mathcal A)\).  The datum
\(\Xi_\mathcal A\) consists of a Costello--Gwilliam HT
BV/factorization model \(\mathcal F_{\mathcal A,\Xi}\), a boundary
comparison
\(\eta^\partial_\mathcal A\colon
\mathrm{Obs}^{\partial}(\mathcal F_{\mathcal A,\Xi})\simeq
\mathcal A\), a bulk--Hochschild comparison
\[
\chi_{\mathrm{HT},\mathcal A}\colon
Z^{\mathrm{der}}_{\mathrm{ch}}(\mathcal A)\to
\mathrm{Obs}^{\mathrm{bulk}}(\mathcal F_{\mathcal A,\Xi}),
\]
and the declared ambient.  The theorem is proved on supplied
realization data and conditional on the existence of \(\Xi\).

The proof now separates the algebraic and physical faces.  The
Costello--Gwilliam envelope verifies the holomorphic boundary face, and
chiral Hochschild computes the algebraic centre
\(\mathrm{ChirHoch}^{\bullet}(\mathcal A,\mathcal A)=
Z^{\mathrm{der}}_{\mathrm{ch}}(\mathcal A)\).  The assertion that this
centre is the HT bulk is exactly the map
\(\chi_{\mathrm{HT},\mathcal A}\), not a consequence of Dunn
additivity alone.  For class \(\mathbf M\), DS--Hochschild controls the
derived-centre side in the weight-completed ambient; HT bulk comparison
still requires \(\Xi\)-compatibility and \(\chi_{\mathrm{HT}}\).

A328 propagates the same arity to active downstream surfaces.  The
heptagon bulk-to-boundary theorem now works with a chosen
\(\Xi\)-representative and \(\eta^\partial\).  The physical UV
finiteness remark now refers to
\(\Phi_{\mathrm{hol}}(\mathcal A,\omega,T_\omega,\Xi_\mathcal A)\) and
to \(\chi_{\mathrm{HT},\mathcal A}\).  The Part VI Monster theorem is
conditional on finite-orbifold BV descent and
\(\Xi_{V^\natural}^{\mathrm{orb}}\).  The Vol I holographic datum
master now states the HT lane as a \(\Xi\)-realized map with
\(\eta^\partial\).

A328 updates `compute/tests/test_universal_holography_functor.py` to
guard the \(\Xi\)-decorated source/target, \(\Xi_\mathcal A\),
\(\eta^\partial_\mathcal A\), \(\chi_{\mathrm{HT},\mathcal A}\), and
absence of the old bare display, `derived centre $=$ bulk`, and
`chain-level equalities at each step compose`.  The same guard checks
the active heptagon, perturbative-finiteness, and Part VI propagation
surfaces.

A328 verification: `python3 -m py_compile` passes on the touched guard
tests.  The focused virtualenv pytest bundle reports `29 passed` across
`test_universal_holography_functor.py`,
`test_programme_climax_holography_gates.py`,
`test_chd_class_m_holography_gates.py`,
`test_part_vi_platonic_introduction.py`,
`test_sc_chtop_heptagon_iv.py`, and
`test_sc_chtop_chain_dioperad.py`.  Fixed-string propagation scans find
no active surviving bare
\(\Phi_{\mathrm{hol}}(\mathcal A)\),
\(\Phi_{\mathrm{hol}}(V^\natural)\), or
\(\Phi_{\mathrm{hol}}(\mathrm{Vir}_c)\) over the checked Vol II,
Vol I, Vol III, and compute surfaces outside negative guard assertions.
`make verify-licensing` exits 0 with zero blockers and 123 warning-class
missing-tag lines elsewhere.  `make fast` converges after two passes
with 2483 pages, zero undefined citations, zero undefined references,
and zero rerun requests.  Direct final-log scans find no fatal TeX
errors, undefined controls, or unresolved citation/reference warnings;
existing pdfTeX destination warnings remain in the log tail.  Scoped
Vol II and Vol I `git diff --check` passes.

A329 repairs the Hochschild chapter's bulk--Hochschild theorem so it
matches the \(\Xi\)-relative Universal Holography discipline.  The live
theorem `thm:bulk_hochschild` still carried the title
`Bulk = chiral Hochschild cochains` and opened by saying that bulk
local operators form the chiral derived centre.  The surrounding prose
was conditional, but the theorem surface still made the physical bulk
look definitional from \(A_\partial\) alone.

The theorem is now titled `Bulk--Hochschild comparison for a chosen
physical prefactorization model; \ClaimStatusConditional; licensing
tags \(\beta+\gamma+\delta\)`.  It first names the algebraic
closed-sector object
\[
Z^{\mathrm{der}}_{\mathrm{ch}}(A_\partial)
\simeq C^\bullet_{\mathrm{ch}}(A_\partial,A_\partial)
\simeq \mathrm{Coder}_0(B^{\mathrm{ch}}(A_\partial)).
\]
For a chosen physical HT prefactorization model \(\mathsf{Obs}\), the
proved comparison is the filtered Ran quasi-isomorphism
\[
\Psi_{\mathrm{Ran}}\colon
\mathcal O_{\mathrm{bulk}}\xrightarrow{\simeq}
C^\bullet_{\mathrm{ch}}(A_\partial,A_\partial).
\]
Equivalently, its inverse homotopy class is the bulk--Hochschild map
\[
\chi_{\mathrm{HT},A_\partial}\colon
Z^{\mathrm{der}}_{\mathrm{ch}}(A_\partial)\xrightarrow{\simeq}
\mathcal O_{\mathrm{bulk}}.
\]
The theorem now explicitly says this comparison is not a definition of
a physical bulk from \(A_\partial\) alone.

A329 propagates the notation through the proof and neighbouring scope
statements.  The associated-graded comparison is now
\(\mathrm{gr}(\Psi_{\mathrm{Ran}})\).  The categorical Morita remark
says the chosen physical bulk is identified by the comparison with the
derived centre.  The bulk--boundary--line factorization theorem now
labels its third clause `Bulk--centre comparison` and says the first
equivalence is \(\Psi_{\mathrm{Ran}}\), equivalently
\(\chi_{\mathrm{HT},A_\partial}\), not a physical bulk construction
from \(A_\partial\) alone.  The open-problem remark now states that
bulk reconstruction remains open; the theorem compares a chosen model
with boundary chiral Hochschild cochains and makes no uniqueness claim.

A329 extends `compute/tests/test_chiral_hochschild_coderivation_model.py`
to require the new title, \(\Psi_{\mathrm{Ran}}\),
\(\chi_{\mathrm{HT},A_\partial}\), and the negative assertions against
the retired title, retired first sentence, `Bulk = derived center`, and
the old bulk-from-boundary wording.

A329 verification: `python3 -m py_compile` passes on the touched guard.
The immediate guard bundle reports `6 passed`; the widened Hochschild IV
bundle reports `31 passed` across `test_climax_theorems_wave11_iv.py`,
`test_climax_theorems_wave17_iv.py`,
`test_chiral_hochschild_coderivation_model.py`, and
`test_global_triangle_boundary_linear.py`.  Fixed-string propagation
scans over active Vol II chapters, appendices, compute tests, Vol I
chapters, and Vol III chapters find no surviving copy of `Bulk =
chiral Hochschild cochains`, `bulk local operators form the chiral
derived centre`, `Bulk = derived center`, or `the bulk is identified
from the boundary chiral Hochschild cochains` outside negative guard
assertions.  `make verify-licensing` exits 0 with zero blockers and 123
warning-class missing-tag lines elsewhere.  `make fast` converges after
two passes with 2483 pages, zero undefined citations, zero undefined
references, and zero rerun requests.  Direct final-log scans find no
fatal TeX errors, undefined controls, or unresolved citation/reference
warnings; existing pdfTeX destination warnings remain in the log tail.
Scoped `git diff --check` passes.

A330 repairs the boundary-linear Landau--Ginzburg bulk/line surface.
The local HKR/Morita theorem is correct, but its surrounding prose still
advertised it as `bulk = shifted cotangent` and said the bulk is
recovered from the boundary algebra or from line Hochschild cochains.
That wording collapsed a model-specific exact-sector comparison into a
general physical reconstruction theorem.

The active core table now says: `model bulk compared with the shifted
cotangent of the derived boundary zero locus`.  The first-principles
prose now says the local operator algebra is the chosen bulk observable
algebra, and that the bridge from open data to bulk data, when it
exists, is a derived-centre comparison in the relevant Morita category.

The prelude to `thm:boundary-linear-bulk-boundary` now states the exact
scope.  In the boundary-linear LG model, the bulk critical locus is the
shifted cotangent of the derived boundary zero locus, hence the model
bulk algebra is compared with the boundary derived centre in the
derived Morita category.  This is explicitly a local exact-sector
comparison, not reconstruction of an arbitrary physical bulk from a
boundary algebra alone.

A330 propagates the same wording to the split copy and active frontier
summary.  The local model bulk is compared with
\(\HH^\bullet(K_\kappa)\) by the derived-centre/Morita comparison, and
the completed model bulk algebra is identified by Morita comparison with
\(\HH^\bullet(K_\kappa)\), equivalently with the derived centre of the
local line category.

A330 extends
`compute/tests/test_deletion_ledger.py::test_boundary_linear_bulk_line_is_morita_scoped`
to require the model-bulk/Morita-comparison wording across the core,
split, and frontier surfaces, and to reject the retired `bulk =
shifted cotangent`, `bulk is recovered from the boundary algebra`,
`bulk is recovered from lines by Hochschild cochains`, and completed
local-bulk-is-computed-by-\(\HH^\bullet(K_\kappa)\) phrases.

A330 verification: `python3 -m py_compile` passes on the touched guard
tests.  The narrow boundary-linear bundle reports `8 passed`; the
widened local/IV bundle reports `80 passed` across
`test_deletion_ledger.py`, `test_boundary_linear_kuranishi.py`,
`test_climax_theorems_wave10_iv.py`, and
`test_climax_theorems_wave15_iv.py`.  Fixed-string propagation scans
over active Vol II chapters, appendices, compute tests, Vol I chapters,
and Vol III chapters find no surviving copy of the retired phrases
outside negative guard assertions.  `make verify-licensing` exits 0
with zero blockers and 123 warning-class missing-tag lines elsewhere.
`make fast` converges after two passes with 2483 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings; existing pdfTeX destination
warnings remain in the log tail.  Scoped `git diff --check` passes.

A331 repairs the active holographic reconstruction theorem surface.
The earlier theorem `thm:holographic-recon-shadow`, the later
class-\(\mathsf M\) corollary, and the Part VI theorem schema still
said that an admissible boundary algebra determines a canonical
relative 3d HT theory whose bulk is the derived chiral centre.  That
was the same false primitive equality repaired in A326--A330: the
physical HT model, boundary comparison, bulk--Hochschild map, and
ambient must be supplied as realization data.

The theorem `thm:holographic-recon-shadow` is now titled `Universal
Holography, admissible shadow projection with realization datum;
\ClaimStatusConditional; licensing tags \(\beta+\gamma+\delta\)`.
Its hypotheses include a chosen HT BV/factorization model
\(\mathcal F_{\cA,\Xi}\), a boundary comparison
\[
\eta^\partial_\cA\colon
\Obs^\partial(\mathcal F_{\cA,\Xi})\xrightarrow{\simeq}\cA,
\]
a bulk--Hochschild comparison
\[
\chi_{\mathrm{HT},\cA}\colon
\Zderch(\cA)\to
\Obs^{\mathrm{bulk}}(\mathcal F_{\cA,\Xi}),
\]
and the ambient \(\mathcal A(\cA)\).  The conclusion is the supplied
model \(\mathcal T_{\cA,\Xi}:=\mathcal F_{\cA,\Xi}\), and the functor
has typed source and target
\[
\mathsf{ChirAlg}^{\omega,\mathrm{adm},\Xi}_X
\longrightarrow
\mathsf{HT\text{-}QFT}^{\Xi}_{X\times\mathbb R}.
\]

A331 also retitles `thm:class-M-chain-bulk` as a class-\(\mathsf M\)
bulk--centre comparison.  The statement now assumes
\(\Xi_\cA^{\mathrm{DS}}\) and says the class-\(\mathsf M\) global
triangle is a triangle of comparison maps in the declared ambient, not
reconstruction of physical bulk from \(\cA\) alone.

The Part VI theorem schema now uses \(\cT_{\cA,\Xi}\),
\(\eta^\partial_\cA\), and \(\chi_{\mathrm{HT},\cA}\).  Its proof no
longer invokes \(\Phi_{\mathrm{hol}}^{-1}\) on a derived-centre image
or bare \(\Phi_{\mathrm{hol}}(\cdots)\) images without the
\(\Xi\)-realized comparison surface.

A331 extends
`compute/tests/test_deletion_ledger.py::test_holographic_reconstruction_is_xi_realized_comparison`
to require the \(\Xi\)-realized theorem title, typed source/target,
\(\mathcal T_{\cA,\Xi}:=\mathcal F_{\cA,\Xi}\),
\(\eta^\partial_\cA\), \(\chi_{\mathrm{HT},\cA}\), the class-M
bulk--centre comparison title, and the Part VI compared-bulk wording.
It rejects the retired canonical-relative-theory sentence,
derived-centre-as-bulk wording, universal bulk identification wording,
and the unlicensed \(\Phi_{\mathrm{hol}}^{-1}\) reconstruction phrase.
`compute/tests/test_climax_theorems_wave12_iv.py::test_class_m_chain_bulk`
now requires the DS--Hochschild bridge, a \(\Xi\)-realization, a
\(\chi_{\mathrm{HT}}\) comparison, and the weight-completed/pro
ambient.

A331 verification: `python3 -m py_compile` passes on the touched
guard tests.  The focused guard bundle reports `2 passed`.  The widened
holography/Part VI bundle reports `84 passed` across
`test_deletion_ledger.py`, `test_climax_theorems_wave12_iv.py`,
`test_universal_holography_functor.py`,
`test_programme_climax_holography_gates.py`,
`test_part_vi_platonic_introduction.py`, and
`test_chd_class_m_holography_gates.py`.  Fixed-string propagation scans
over active Vol II chapters, appendices, compute tests, Vol I chapters,
and Vol III chapters find no surviving copy of the exact retired
phrases outside negative guard assertions; a second scan finds no
untyped local source/target display or bare \(\mathcal T_\cA\) notation
on the touched theorem surfaces.  `make verify-licensing` exits 0 with
zero blockers and 123 warning-class missing-tag lines elsewhere.  `make
fast` converges after two passes with 2487 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings; existing pdfTeX destination
warnings remain in the log tail.  Scoped `git diff --check` passes.

A332 repairs the shadow-tower reconstruction theorem adjacent to the
\(\Xi\)-realized holographic reconstruction surface.  The theorem had
said that the full datum \(\Theta_\cA\) is determined by Hessian,
cubic, quartic, and the obstruction tower.  The obstruction-extension
sequence shows the missing datum: \(H^2\)-obstruction classes decide
liftability, but when the obstruction vanishes the set of lifts is an
\(H^1(J^{r+1}(\cA),d_2)\)-torsor.  Obstruction classes alone do not
choose the filtered MC lift.

The theorem is now titled `Filtered shadow-tower reconstruction of the
bar-intrinsic MC gauge class; \ClaimStatusConditional; licensing tags
\(\alpha+\gamma+\epsilon\)`.  Its input is a shadow reconstruction
packet
\[
\mathfrak S(\cA)=
(S_2,S_3,S_4,\{(o_n(\cA),\lambda_n)\}_{n\ge5}),
\]
where \(o_n(\cA)\in H^2(J^n(\cA),d_2)\) is the obstruction class and
\(\lambda_n\in H^1(J^n(\cA),d_2)\) is the lift coordinate in the
chosen bar-intrinsic gauge slice.

A332 also repairs the local Mittag--Leffler corollary.  It is now
`Mittag--Leffler for the bar-intrinsic shadow branch;
\ClaimStatusConditional; licensing tags \(\gamma+\epsilon\)`.  It
defines \(E^\Theta_\cA(r)\subset E_\cA(r)\), proves
\(\varprojlim^1 E^\Theta_\cA(r)=0\), and explicitly makes no
surjectivity claim for arbitrary components of \(E_\cA(r)\).

The strict \(E_3\) projection corollary is now \(\Xi\)-realized:
\[
Z^{\mathrm{der}}_{\mathrm{ch}}(\mathcal A)
\xrightarrow{\chi_{\mathrm{HT},\mathcal A}\simeq}
\mathrm{Obs}^{\mathrm{bulk}}(\mathcal T_{\mathcal A,\Xi}).
\]

A332 extends
`compute/tests/test_deletion_ledger.py::test_shadow_tower_reconstruction_requires_lift_coordinates`
to require the reconstruction packet, \(\lambda_n\)-coordinates,
bar-intrinsic gauge slice, branch-level ML statement, and
\(\Xi\)-realized strict \(E_3\) comparison.  It rejects the retired
obstruction-only reconstruction theorem, bare strict-\(E_3\) projection,
global \(\varprojlim^1 E_\cA(r)=0\) claim, and `Surjectivity for
general ...` proof sentence.

A332 verification: `python3 -m py_compile` passes on the touched
guard tests.  The focused guard bundle reports `2 passed`.  The widened
holography bundle reports `85 passed` across
`test_deletion_ledger.py`, `test_climax_theorems_wave12_iv.py`,
`test_universal_holography_functor.py`,
`test_programme_climax_holography_gates.py`,
`test_part_vi_platonic_introduction.py`, and
`test_chd_class_m_holography_gates.py`.  Fixed-string propagation scans
over active Vol II chapters, appendices, compute tests, Vol I chapters,
and Vol III chapters find the retired exact formulas only in negative
guard assertions; a broader scan finds no surviving theorem title or
local ML statement advertising obstruction-only reconstruction or
global extension-moduli surjectivity.  `make verify-licensing` exits 0
with zero blockers and 123 warning-class missing-tag lines elsewhere.
`make fast` converges after two passes with 2487 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  Direct
final-log scans find no fatal TeX errors, undefined controls, or
unresolved citation/reference warnings; existing pdfTeX destination
warnings remain in the log tail.  Scoped `git diff --check` passes.

A350 repairs `prop:gravity-m5`, the next gravity-climax theorem surface
after A349.  The proposition no longer advertises a quinary Virasoro
operation with hidden polynomials.  It is now a Stasheff-gauge
representative \(m_5^{h_{\mathrm{Vir}}}\) with
`\ClaimStatusProvedHere` and licensing
\(\alpha+\gamma+\varepsilon\) via consecutive collision coordinates,
`\hypAmbientWtCpl`, and the Virasoro BRST contraction.

A350 also repairs the executable witness.  `compute/lib/symbolic_stasheff.py`
now treats field labels by derivative order and applies
\((\lambda+\partial)^n\) generically in the final slot.  This removes
the old truncation where `_apply_partial` stopped at `d3T` and the
arity-five witness could not produce the leading \(d^4T\)-coefficient.
`mk_exact_numerical(5,...)` now exposes the full six-component field
vector \(d^4T,d^3T,d^2T,dT,T,1\).

The repaired symbolic recursion gives the symmetric field/scalar vector
\[
(-1,\,-14,\,-71,\,-154,\,-120,\,-5c),
\]
so
\[
m_5^{h_{\mathrm{Vir}}}(T^5;\lambda,\lambda,\lambda,\lambda)
=-\partial^4T-14\lambda\partial^3T-71\lambda^2\partial^2T
-154\lambda^3\partial T-120\lambda^4T-5c\lambda^6 .
\]
The theorem now displays the full \(P_1(\ell)\), \(P_0(\ell)\), and
\(Q_5(\ell)\) polynomials, with
\(Q_5(\lambda,\ldots,\lambda)=-60\lambda^6\).  Its proof gives the
explicit nine-term arity-five source
\[
m_2\circ m_4+m_4\circ m_2+m_3\circ m_3
\]
with collapsed consecutive spectral parameters and the slot
sesquilinearity rules.  The source-tree witness is
`stasheff_rhs_arity5` / `m5_virasoro_symbolic`.

A350 repairs the downstream quintic-shadow paragraph as well.  The old
text claimed the raw symmetric scalar part was
\(-(17c/6)\lambda^6\) and used a Shapovalov division to recover
\(S_5=-48/[c^2(5c+22)]\).  That was incompatible with the live
proposition and the computation.  The paragraph now separates the two
scalar-shadow verifications of \(S_5\) from the raw arity-five operator
check.  The raw scalar contact term is \(-5c\lambda^6\); it is not the
shadow coefficient.  A genuine third route would require constructing
the projection
\(\Pi^{(5)}_{\mathrm{sh}}\colon m_5^{h_{\mathrm{Vir}}}\mapsto S_5\),
which is named as a missing map rather than silently assumed.

A350 adds `compute/tests/test_gravity_m5_stasheff_gauge.py`.  The guard
requires the theorem's status/licensing/source witness, rejects the old
placeholder prose, checks all six exact symbolic coefficients, checks
the symmetric specialization, checks agreement with the numerical
`StasheffEngine`, and guards the separation between raw \(m_5\)
arithmetic and the normalized shadow coefficient.

A350 verification: `compute/.venv/bin/python -m py_compile` passes on
`compute/lib/symbolic_stasheff.py` and the new guard.  The focused
pytest bundle reports `189 passed, 1 skipped, 22 subtests passed`
across `test_gravity_m5_stasheff_gauge.py`,
`test_gravity_m4_stasheff_gauge.py`,
`test_field_sector_generating_function.py`,
`test_swiss_cheese_virasoro_wheels.py`, and
`test_gravity_3d_engine.py`.  Fixed-string scans find no live active
copy of the retired placeholder phrases, the `d3T in right slot` skip,
or the false `17c/6` scalar normalization outside negative guard
assertions and historical/quarantined notes.  Scoped `git diff --check`
passes.  `make verify-licensing` exits 0 with zero blockers and 118
warning-class missing-tag lines elsewhere; the next visible gravity
warning is `thm:gap-migration`.  `make fast` converges after two passes
with 2493 pages, zero undefined citations, zero undefined references,
and zero rerun requests.

A351 repairs `thm:gap-migration` and the adjacent \(m_6\) computation.
The old `comp:m6-depth-spectrum` claimed
\(\operatorname{Spec}(m_6)=\{1,2,3,4,5,7\}\).  Recomputing with
`compute/m7_m10_depth_frontier.py::StasheffEngine` at generic spectral
parameters gives field depths \(\{2,3,4,5\}\) and scalar depth~\(7\):
the \(d^5T\) and \(d^4T\) coefficients vanish, while
\(d^3T,d^2T,dT,T\) and the scalar term are nonzero.  Thus the correct
generic spectrum is
\[
\operatorname{Spec}(m_6)=\{2,3,4,5,7\}.
\]

A351 retitles `thm:gap-migration` as a conditional theorem with
licensing \(\gamma+\varepsilon\) via `\hypAmbientWtCpl+\effKoszul`.
It now separates: the Virasoro weight-depth identity, the structural
gap \(d=n\), the even-arity secondary cancellation
\(\operatorname{Spec}(m_4|_T)=\{2,3\}\) and
\(\operatorname{Spec}(m_6|_T)=\{2,3,4,5\}\), and the conditional
principal \(\mathcal W_N\) T-propagated binding-sector formula
\[
d_{\mathrm{gap}}^{\mathrm{bind}}(\mathcal W_N,n)=2N+n-4.
\]
The formula is no longer advertised as a total-spectrum theorem, nor as
a raw all-\(W_N\) sector statement.

A351 rewrites the \(m_6\) table.  Depths \(0\) and \(1\) are secondary
Stasheff cancellations; depth \(6\) is the structural weight-\(1\) gap;
depth \(7\) is the scalar contact term.  The old free-monomial counts
are gone because they were not computed support counts.  The
period-\(2\) verification remark now treats degree~\(2\) as the base
case rather than an odd case.

A351 propagates the sector scope to `chapters/theory/introduction.tex`,
`chapters/examples/w-algebras-w3.tex`,
`chapters/examples/w-algebras-stable.tex`,
`chapters/examples/w-algebras-frontier.tex`,
`chapters/examples/rosetta_stone.tex`, and
`chapters/connections/conclusion.tex`.  In particular, the W4 frontier
now distinguishes the raw \((W_4,W_4,W_4)\) gap from the theorem's
T-propagated \((W_4,W_4,T)\) value.

A351 also updates `compute/m6_depth_spectrum.py` to track `d5T` and
use derivative-order-aware partial shifts, so the exploratory exact
script can no longer hide a possible depth-\(0\) term.  The new guard
`compute/tests/test_gap_migration_m6.py` checks the generic \(m_6\)
field depths, the corrected theorem/table language, and the propagated
T-propagated binding-sector scope.

A351 verification: `compute/.venv/bin/python -m py_compile` passes on
`compute/m6_depth_spectrum.py` and the new guard.  The focused pytest
bundle reports `192 passed, 1 skipped, 22 subtests passed` across
`test_gap_migration_m6.py`, the \(m_4/m_5\) Stasheff guards,
field-sector, Virasoro-wheel, and 3d gravity engine tests.  Fixed scans
over active chapters and compute files find no live stale
`Spec(m_6)={1,...}` spectrum, no free-monomial-count sentence, and no
unscoped `was proved in general` / total-spectrum \(2N+n-4\) claim
outside negative guard assertions.  Scoped `git diff --check` passes.
`make verify-licensing` exits 0 with zero blockers and 117
warning-class lines elsewhere; the next gravity warning is
`thm:period-2-parity`.  `make fast` converges after two passes with
2493 pages, zero undefined citations, zero undefined references, and
zero rerun requests.

A352 repairs `thm:period-2-parity` and the dependent
`thm:gravity-c-linearity` lane.  The old period-\(2\) theorem had two
mathematical defects: it said odd symmetric field depths stopped at
\(\{0,\ldots,k-2\}\), although
\(\varphi_k(x)=(-1)^nC_n\prod_{m=2}^k(x+m)\) has degree \(k-1\);
and it identified the field signed sum
\(\varphi_k(1)\) with the scalar shadow through the false relation
\(P_k=12S_k/c\).

A352 retitles the theorem as conditional with licensing
\(\gamma+\varepsilon\) via `\hypAmbientWtCpl+\effKoszul`.  The
rightmost-reduction recurrence is now the load-bearing all-arity
hypothesis, and the scalar clause separately assumes the symmetric
scalar--\(T\) proportionality.  The field depth range is now
\(d=0,\ldots,k-1\).  The field signed sum is renamed
\(\Sigma_k^{\mathrm{fld}}=\varphi_k(1)\), while the scalar formula is
proved from
\[
[1]m_k=(c/12)P_k=(c/24)[T]m_k,
\qquad
P_k(1,\ldots,1)=\frac12\varphi_k(0)
=(-1)^nC_n\,k!/2 .
\]

A352 also makes `thm:gravity-c-linearity` conditional/licensed and
corrects its field-derivative range to \(0\le j\le k-1\).  The
two-colour, partition-independence, and ordered-shadow refinement
paragraphs now distinguish the scalar shadow from the \(x=1\) field
signed sum and state the possible field-depth range \(0,\ldots,k-1\)
with structural gap \(d=k\).

A352 extends `compute/tests/test_catalan_factorisation.py` through even
\(k=12\) and odd \(k=13\), adds scalar Catalan and scalar--\(T\)
proportionality checks through odd \(k=13\), and adds a manuscript
guard rejecting the old scalar/field conflation.  The adjacent
`test_gap_migration_m6.py` delimiter was updated to the new theorem
title.

A352 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_catalan_factorisation.py -q` passes with
`53 passed`.  Adjacent focused tests pass:
`test_gap_migration_m6.py` (`3 passed`),
`test_gravity_m5_stasheff_gauge.py` (`4 passed`),
`test_gravity_m4_stasheff_gauge.py` (`2 passed`), and
`test_field_sector_generating_function.py`
(`9 passed, 22 subtests passed`).  `make verify-licensing` exits 0
with zero blockers and 115 warning-class lines elsewhere; the next
gravity warning is `thm:gravity-finite-presentation`.  `make fast`
converges after two passes with zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped
`git diff --check` passes.  A full `compute/tests` run was interrupted
after one failure marker because it ran long; rerunning with `-x`
stops at unrelated guard
`test_chd_ds_hochschild_iv.py::test_chd_source_contains_ds_hpl_transfer_formulas`
after `713 passed`, because `chapters/theory/chiral_higher_deligne.tex`
lacks the exact source string
`Q_{\mathrm{DS}}h+hQ_{\mathrm{DS}}=1-ip`.

A353 repairs the DS-Hochschild compatibility theorem exposed by the
A352 full-suite stop.  The old `thm:chd-ds-hochschild` surface treated
principal/hook DS-Hochschild compatibility at all non-critical levels as
unconditional and used a false generic `C_2`-cofiniteness reading of
Arakawa.  Generic non-critical \(W\)-algebras are not lisse in that
sense; the theorem needs an explicit DS deformation retract, a
functorial chiral HKR/brace model, and weightwise completed HPL
convergence.

A353 retitles `thm:chd-ds-hochschild` as
`\ClaimStatusConditional`, licensed by
\(\alpha+\beta+\gamma+\delta\).  The theorem now assumes the good
grading and fractional cover/descent datum when required, the completed
DS-Hochschild SDR
\[
Q_{\mathrm{DS}}h+hQ_{\mathrm{DS}}=\mathrm{id}-i\circ p,
\qquad
p i=\mathrm{id},\quad ph=hi=h^2=0,
\]
the chiral HKR/brace model, completed HPL tree convergence in
`\hypAmbientWtCpl`, and the separate mixed-lift package.  The proof now
uses KRW for the BRST complex and Arakawa only on the stated
finite/lisse/admissible surfaces; it no longer asserts generic
non-critical `C_2`-cofiniteness.

A353 propagates the same package to `chapters/frame/part_viii_synthesis.tex`,
`notes/part_viii_vol_iv_roadmap_synthesis.md`, and
`notes/first_principles_cache_comprehensive.md`, so gravity-chain and
Part VIII summaries inherit the DS SDR/HKR/completed-HPL hypotheses
instead of advertising an unconditional Arakawa-plus-HPL closure.  The
DS guard now requires the unambiguous
`\mathrm{id}-i\circ p` formula and rejects the ambiguous `1-ip`
shorthand.  The deletion-ledger scalar-lane guard was also synced to
the live THQG proposition title, `Shadow scalar coefficient in closed
form`.

A353 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_chd_ds_hochschild_iv.py
compute/tests/test_chd_boardman_vogt_tree_homotopy.py -q` passes with
`10 passed`.  The DS-Hochschild plus universal-holography bundle passes
with `8 passed`; the climax DS-reduction plus CHD bundle passes with
`5 passed`; the focused THQG scalar-lane guard passes with `1 passed`;
the whole deletion ledger passes with `62 passed`.  Fixed-string scans
find no remaining `1-ip` on the repaired DS/HPL surfaces.  Scoped
`git diff --check` passes.  A full `compute/tests -q -x` rerun advanced
past the previous CHD and THQG failures but the tool session exited
abnormally near 23% with code `-1` and no pytest failure summary, so it
is not counted as a full-suite pass.  `make verify-licensing` exits 0
with zero blockers and 115 warning-class lines elsewhere.  `make fast`
converges after two passes at 2493 pages with zero undefined citations,
zero undefined references, and zero rerun requests.

A355 repairs `prop:graviton-tracelessness`, the next gravity licensing
warning.  The old proposition defined the derivative-shape polynomial
with range \(0,\ldots,k-2\), contradicting the repaired period-\(2\)
field polynomial where odd arities have degree \(k-1\).  It also
attempted an independent Stasheff induction and asserted an automatic
general-generator formula \(G_k^{(a)}(-h_a)=0\).

A355 retitles the proposition as `Symmetric-point graviton
tracelessness`, with `\ClaimStatusConditional` and licensing
\(\gamma+\varepsilon\) via `\hypAmbientWtCpl+\effKoszul`.  The field
polynomial is now
\[
G_k(x)=\sum_{j=0}^{k-1}a_j^{(k)}x^j,\qquad
a_j^{(k)}=[\partial^jT]\,m_k(T^{\otimes k};1,\ldots,1).
\]
The proof is a direct consequence of `thm:period-2-parity`: \(G_2=x+2\),
even \(k\ge4\) gives \(G_k=0\), and odd \(k=2n+3\) gives
\[
G_k(x)=(-1)^nC_n\prod_{m=2}^{k}(x+m),
\]
so \(x+2\) is an explicit factor.  The general-generator claim is now
only a Ward-factor criterion: for a generator \(a\), the analogous
statement requires its own field polynomial to have factor \(x+h_a\).

A355 adds a Catalan manuscript guard requiring the proposition
status/licensing, the \(0,\ldots,k-1\) range, the period-\(2\) theorem
reference, the product \(\prod_{m=2}^{k}(x+m)\), and the Ward-factor
language.  It rejects the old \(k-2\) range, the independent induction
phrase, and the automatic general-generator formula.  The A354
gravity-engine delimiter was also synced to the strengthened
proposition title.

A355 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_catalan_factorisation.py -q` passes with
`54 passed`.  `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py -q` passes with
`69 passed, 1 skipped`.  `make verify-licensing` exits 0 with zero
blockers and 113 warning-class lines elsewhere; the next gravity
warning is now `thm:graviton-resolvent-closed`.  `make fast` converges
after two passes at 2493 pages with zero undefined citations, zero
undefined references, and zero rerun requests.

A354 repairs `thm:gravity-finite-presentation`, the next gravity
licensing warning.  The old theorem said the Virasoro tower was
uniquely determined by two data, \(m_2\) and \(h\), but its proof used
the DS perturbation, inclusion, and projection.  The upstream general
finite-presentation proposition carried the same compression, and the
gravity compute oracle still wrote the ambiguous SDR identity
`Qh + hQ = 1 - ip`.

A354 upgrades `prop:finite-presentation-general` to finite
HPL-datum presentation:
\[
\mathfrak H_{\cA}=(V,d_0,m_2^0,\delta,i,p,h),
\]
with `\ClaimStatusProvedHere` and licensing \(\gamma+\delta\) via
`\hypAmbientWtCpl` and the complete HPL filtration.  The proof now uses
the planar-tree formula with \(p\) at the root, \(i\) at leaves, \(h\)
on internal edges, and \(m_2\) or the reduced perturbation on binary
vertices.

A354 retitles `thm:gravity-finite-presentation` as conditional,
licensed by \(\gamma+\delta+\varepsilon\) via
`\hypAmbientWtCpl`, `\hypKZSDR`, and `\effKoszul`.  The theorem now
fixes the affine level \(k_{\mathrm{aff}}\ne -2\), writes
\[
c=1-6(k_{\mathrm{aff}}+1)^2/(k_{\mathrm{aff}}+2),
\]
and assumes the DS/KZ contraction
\[
(V_{\mathrm{DS}},Q_{\mathrm{DS}},m_2^{\mathrm{aff}},
\delta_{\mathrm{DS}},i,p,h_{\mathrm{DS}}),
\qquad
Q_{\mathrm{DS}}h_{\mathrm{DS}}+h_{\mathrm{DS}}Q_{\mathrm{DS}}
=\mathrm{id}-i\circ p .
\]
Finite presentation now means that this finite DS/KZ HPL datum
produces the infinite tower; non-truncation remains conditional on the
wheel-residue survival surface of `prop:vir-truncation`.

A354 also propagates the finite-datum wording to the dg-shifted Yangian
bridge and replaces stale `1 - ip` notation in the gravity engine,
first-principles note, and attack-finding note with `id - i∘p`.  The
gravity engine guard now requires the full finite HPL datum tuple and
rejects the old “two data” theorem block.

A354 verification: `compute/.venv/bin/python -m py_compile
compute/lib/gravity_3d_engine.py` passes.  The focused HPL transfer
tests pass with `6 passed`; the full gravity engine file passes with
`69 passed, 1 skipped`.  Fixed-string scans find no live `1 - ip`, no
live `presented by $(m_2,h)`, and no live “one binary operation and one
homotopy” outside negative test guards.  `make verify-licensing` exits
0 with zero blockers and 114 warning-class lines elsewhere; the next
gravity warning is now `prop:graviton-tracelessness`.  `make fast`
converges after two passes at 2493 pages with zero undefined citations,
zero undefined references, and zero rerun requests.

A356 repairs `thm:graviton-resolvent-closed`, the next exposed gravity
resolvent theorem after A355.  The old theorem called
\[
G_{\mathrm{scal}}(t)=\int_0^t s\sqrt{Q_{\mathrm{Vir}}(s)}\,ds
\]
an elliptic integral, although \(Q_{\mathrm{Vir}}\) is quadratic and
the double cover \(u^2=Q_{\mathrm{Vir}}(t)\) has genus zero.  It also
omitted the square-root branch \(\sqrt{Q_{\mathrm{Vir}}(0)}=c\),
displayed a branch-point formula off by a factor of two, treated scalar
branch points as graviton bound-state spectrum without a spectral
comparison datum, and left the later \(S_9\) table with stale raw rows
\(S_2=c/12\), \(S_3=-c\).

A356 retitles the theorem as `Scalar shadow resolvent on the Virasoro
metric branch`, with `\ClaimStatusConditional` and licensing
\(\gamma+\delta\) via `\hypAmbientWtCpl+\hypStokes`.  The theorem now
chooses the branch \(\sqrt{Q_{\mathrm{Vir}}(0)}=c\), identifies the
integral as an Abelian integral on the rational double cover, extracts
\(S_r\) for all \(r\ge2\), and gives the corrected branch points
\[
t_\pm=\frac{c(5c+22)}{180c+872}
\left(-6\pm4i\sqrt{\frac{5}{5c+22}}\right).
\]
The normalized low-degree rows are now \(S_2=c/2\) and \(S_3=2\); the
raw ternary scalar proportional to \(-c\) is explicitly not denoted
\(S_3\).

A356 adds manuscript guards to
`compute/tests/test_shadow_borel_resurgence.py` requiring the theorem
status/licensing, branch choice, genus-zero Abelian-integral wording,
corrected branch-point formula, normalized \(S_2,S_3\) table rows, and
exact fraction values through \(S_5\) at \(c=13\).  The guard rejects
the retired `elliptic integral`, `No real pole`, `bound-state`, stale
half-factor branch formula, stale \(c/12,-c\) table rows, and `not its
Taylor coefficients` wording.  The exploratory compute-script narrative
was also corrected to an algebraic-logarithmic Abelian integral on a
genus-zero double cover.

A356 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_shadow_borel_resurgence.py -q` passes with
`84 passed`.  `compute/.venv/bin/python -m pytest
compute/tests/test_celestial_ope_from_shadow.py -q` passes with
`83 passed`.  `make verify-licensing` exits 0 with zero blockers and
112 warning-class lines elsewhere.  `make fast` converges after two
passes at 2493 pages with zero undefined citations, zero undefined
references, and zero rerun requests.  Scoped `git diff --check` passes.

A357 repairs `thm:gravity-koszul-triangle`, the next gravity licensing
warning after A356.  The old theorem had no status/licensing and still
allowed the central-extension line \(\C[\![c]\!]\) to be read as the
whole bulk \(\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)\).  Its introduction
and Steinberg remark also said that the bar complex extracts boundary,
lines, and bulk by three functors, erasing the distinction between
bar-cobar inversion, line-side Koszul duality, and the chiral
Hochschild derived centre.

A357 retitles the theorem as `Central-sector gravitational Koszul
triangle`, with `\ClaimStatusConditional` and licensing
\(\alpha+\beta+\gamma+\varepsilon\) via
`\hypBHdict+\hypAmbientWtCpl+\effKoszul`.  The theorem now works in the
Brown--Henneaux boundary chart \(b_{\mathrm{BH}}\), in the completed
ambient, and on the effective Koszul locus.  Its vertices are:
\[
A_{b_{\mathrm{BH}}}=\mathrm{Vir}_c,\qquad
\widehat{\mathcal C}^{b_{\mathrm{BH}},\mathrm{pert}}_{\mathrm{line}}
\simeq \mathcal A^!_{\mathrm{line}}\text{-}\mathbf{mod}^{A_\infty},
\qquad
\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)
\xrightarrow{\pi_{\mathrm{cent}}}
\HH^0_{\mathrm{GF}}\oplus\HH^2_{\mathrm{GF}}[-2].
\]
The connective classical shadow of the last projection is
\(\C[\![c]\!]\), read as the cosmological-constant line; it is
explicitly not the full derived chiral centre.  The proof now treats
Brown--Henneaux as an external chart datum and keeps the
\(\mathrm{Vir}_{26-c}\)-module realization of lines conjectural.

A357 adds manuscript guards to `compute/tests/test_gravity_3d_engine.py`
requiring the central-sector theorem title, status/licensing,
\(b_{\mathrm{BH}}\), \(\pi_{\mathrm{cent}}\), the GF two-term target,
the completed line category, and the explicit non-identification of
\(\C[\![c]\!]\) with the full derived centre.  The guard rejects the
old theorem title, old `Bulk central-extension projection` vertex,
three-functor bar-extraction phrase, and homotopy-Koszulity shortcut.
The stale Catalan delimiter was updated to the A356 scalar-resolvent
title.  The FM125 line in
`notes/part_VI_climax_platonic_reconstitution.md` now records the
correct projection status instead of claiming a full derived
equivalence.

A357 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py -q` passes with
`71 passed, 1 skipped`.  `compute/.venv/bin/python -m pytest
compute/tests/test_catalan_factorisation.py -q` passes with
`54 passed`.  `compute/.venv/bin/python -m pytest
compute/tests/test_climax_theorems_wave17_iv.py -q` passes with
`20 passed`.  `make verify-licensing` exits 0 with zero blockers and
111 warning-class lines elsewhere; the next gravity warning is
`prop:gravity-koszul-dual`.  `make fast` converges after two passes at
2495 pages with zero undefined citations, zero undefined references,
and zero rerun requests.  Scoped `git diff --check` passes.

A358 repairs `prop:gravity-koszul-dual`, the next exposed gravity
licensing warning after A357.  The old proposition and its compute
guards still encoded the false primitive slogan
\(\mathrm{Vir}_c^!=\mathrm{Vir}_{26-c}\).  That collapsed the
ordered-bar Koszul object
\[
A^{i,\mathrm{ord}}_{\mathrm{Vir}}
  =H^\bullet(\barB^{\mathrm{ord}}(\mathrm{Vir}_c))^\vee,
\]
the Verdier dual \(A^!\), the abstract homotopy line-side dual
\((\mathrm{Vir}_c)^{!_{\mathrm{line}},\infty}\), and the strict
same-family Virasoro representative \(\mathrm{Vir}_{26-c}\).

A358 retitles the proposition as `Same-family Virasoro representative
of the line-side dual; \ClaimStatusConditional`, with licensing
\(\beta+\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effKoszul`.  The representative now enters only
after the comparison datum
\[
\chi_{\mathrm{Vir}}^{\mathrm{line}}\colon
(\mathrm{Vir}_c)^{!_{\mathrm{line}},\infty}
\xrightarrow{\sim}\mathrm{Vir}_{26-c}.
\]
The scalar identity is only
\[
\kappaChHodge(\mathrm{Vir}_c)+\kappaChHodge(\mathrm{Vir}_{26-c})
=c/2+(26-c)/2=13,
\]
with fixed point \(c=13\).  The proposition explicitly refuses to
identify the Verdier dual, the ordered-bar Koszul object, or the
abstract line package with \(\mathrm{Vir}_{26-c}\) away from this
comparison surface.  The concrete equivalence of line operators with
highest-weight \(\mathrm{Vir}_{26-c}\)-modules remains
Conjecture~`conj:gravity-line-identification`.

A358 also propagates the correction through active compute and memory
surfaces: gravity, P2 line OPE, Pentagon trace, bulk-boundary duality,
Hochschild, cross-volume, shadow-resurgence, holographic HT,
modular-PVA, genus-one/genus-two, gauge, SC bar-cobar, YM synthesis,
first-principles-cache, and current swarm/audit notes.  The helper
`verify_theorem_j_virasoro` now marks the concrete
\(\mathrm{Vir}_{26-c}\)-module model as conjectural rather than
proved.

A358 verification: fixed-string scans over active
`chapters/connections/3d_gravity.tex`, `compute/lib`, and
`compute/tests` find no surviving raw Virasoro bang-dual equality and no live
`Koszul dual ... Vir` / `A! ... Vir_{26-c}` slogan.  Focused pytest
passed across the touched surface: gravity (`72 passed, 1 skipped`),
P2 OPE (`21 passed`), P2 Pentagon trace (`27 passed`),
bulk-boundary (`69 passed`), Hochschild (`49 passed`), cross-volume
(`68 passed`), shadow resurgence (`84 passed`), holographic HT
(`70 passed`), modular PVA (`86 passed`), DS-BRST (`72 passed`),
genus-one kappa (`67 passed`), celestial shadow (`83 passed`), YM
synthesis (`41 passed`), genus-one bridge (`64 passed`), genus-two
obstruction (`59 passed`), gauge orbit (`40 passed`), SC bar-cobar
(`60 passed`), line operators (`26 passed`), and cross-engine
consistency (`45 passed, 2 xfailed`).  `make verify-licensing` exits
0 with zero blockers and 110 warning-class lines elsewhere; the next
gravity warning is `prop:gravity-ybe`.  `make fast` converges after
two passes with zero undefined citations, zero undefined references,
and zero rerun requests.  Scoped `git diff --check` passes.

A359 repairs the remaining Virasoro-specific self-duality language
exposed by A358.  The correct statement is not that
\(\mathrm{Vir}_c\) is primitively self-dual.  The strict same-family
line-side comparison sends \(c\mapsto 26-c\), so its fixed-point
equation is \(c=26-c\), hence \(c=13\).  The value \(c=26\) is instead
the critical-string/effective-anomaly-cancellation value and is not
fixed by the comparison.

A359 rewrites the active gravity chapter so the scalar resolvent,
Cardy capacity, Page/entropy formulas, determinant normalization, and
Nariai continuation refer to the same-family comparison fixed point
\(c=13\).  The phrase "self-dual 3d gravity" is removed from the live
gravity text.

A359 also renames Virasoro-specific executable oracles from
self-duality language to fixed-point language:
`virasoro_line_comparison_fixed_point`,
`line_comparison_fixed_point`, `line_comparison_fixed`,
`comparison_fixed`, `is_comparison_fixed`, and
`comparison_fixed_analysis`.  The shadow-resurgence and exact-WKB
return dictionaries now record equality of the \(c\) and \(26-c\)
comparison data, not object self-duality.  The cross-volume Virasoro
table now uses `comparison_fixed_at` and `NOT_fixed_at_26`.

A359 verification: focused pytest passed on the touched suites:
gravity (`72 passed, 1 skipped`), bulk-boundary (`69 passed`),
Hochschild (`49 passed`), shadow resurgence (`84 passed`), exact WKB
(`45 passed`), cross-volume (`68 passed`), DS-BRST (`72 passed`),
holographic HT (`70 passed`), genus-one kappa (`67 passed`),
celestial shadow (`83 passed`), universal celestial (`11 passed`),
modular PVA (`86 passed`), genus-two obstruction (`59 passed`),
genus-two graph sum (`84 passed`), raw direct-sum obstruction
(`4 passed`), class-M direct-sum obstruction (`7 passed`), Swiss
cheese Virasoro wheels (`106 passed`), and the combined support bundle
`test_ym_synthesis_engine.py test_koszul_epstein_steps_bc.py
test_affine_kac_moody_pva.py test_adversarial_verification.py`
(`237 passed`).  Fixed-string scans over active `3d_gravity.tex`,
`compute/lib`, and `compute/tests` find no Virasoro-specific
self-duality slogan and no \(c=26\) self-dual mislabel.  Remaining
self-duality hits concern non-Virasoro \(E_\infty\), Yangian,
Heisenberg, affine, or \(W\)-algebra statements.  `make
verify-licensing` exits 0 with zero blockers and 110 warnings
elsewhere; `make fast` converges after two passes at 2495 pages with
zero undefined citations, zero undefined references, and zero rerun
requests.  Scoped `git diff --check` passes.

A360 repairs `prop:gravity-ybe`, the next exposed gravity licensing
warning after A359.  The old proposition stated "The CYBE" without
claim status, licensing tags, ambient completion, or the distinction
between the pre-\(d\log\) Laplace/OPE kernel
\[
r^L(z)=\partial T/z+2T/z^2+(c/2)/z^4
\]
and the bar collision kernel
\[
r^{\mathrm{coll}}_{\mathrm{Vir}}(z)=(c/2)/z^3+2T/z.
\]
That wording risked reading the Virasoro cubic-pole, operator-valued
kernel as an ordinary finite-dimensional Casimir solution
\(\Omega/z\), and risked conflating the classical collision-residue
CYBE with the strict quantum YBE for the fusion kernel.

A360 retitles the proposition as `Collision-residue gravitational
CYBE; \ClaimStatusConditional`, with licensing
\(\alpha+\beta+\gamma\) via `\hypBHdict+\hypAmbientWtCpl`.  The
statement now fixes the Brown--Henneaux boundary chart, defines
\(r^{\mathrm{coll}}_{\mathrm{Vir}}\), invokes
Theorem~`thm:thqg-V-cybe-from-arnold`, and states the displayed
three-channel CYBE only in the completed line-side endomorphism
algebra.  It explicitly says that before the line comparison is
installed the equality is the arity-three MC equation for
\(\Theta_{\mathrm{Vir}_c}\), and that it is neither the ordinary
finite-dimensional Casimir identity nor the strict quantum fusion-kernel
YBE before associator/normalisation data are chosen.

A360 propagates the same convention through the active surface:
`spectral-braiding-core.tex` now treats the explicit computation as
Virasoro PVA Jacobi in Laplace variables, with the collision-residue
CYBE as the completed line-side shadow.  The non-simply-laced B/C
paragraph now cites the general collision-residue theorem rather than
the Virasoro proposition.  The worked Virasoro example no longer denies
the collision-residue CYBE; it rejects only the wrong
finite-dimensional Casimir reading.  The introduction summary table now
distinguishes "Laplace Jacobi / collision-residue CYBE".  The
Virasoro r-matrix descriptors in `holographic_ht_engine.py` and
`bulk_boundary_duality_engine.py` now store
\((c/2)/z^3+2T/z\), not the Laplace/OPE kernel.  Stale
`self-dual at c=13` wording encountered in the touched compute surfaces
was corrected to comparison-fixed language.

A360 adds `virasoro_collision_cybe_scope()` and source-level guards in
`test_gravity_3d_engine.py`, plus Virasoro collision-expression checks
in `test_holographic_ht_engine.py` and
`test_bulk_boundary_duality_engine.py`.

A360 verification: focused pytest passed across
`test_gravity_3d_engine.py`, `test_holographic_ht_engine.py`,
`test_bulk_boundary_duality_engine.py`, and
`test_spectral_braiding.py` with `236 passed, 1 skipped`.  Stale scans
over the touched active surfaces found no surviving Virasoro-specific
`Gravitational YBE` theorem title, no
`Virasoro CYBE (quadratic, non-constant)` compute descriptor, no
\(z^{-4}\) Laplace kernel stored as collision residue, and no
`self-dual at c=13` wording.  `make verify-licensing` exits 0 with
zero blockers and 109 warnings elsewhere; the next gravity warning is
`prop:ds-sdr`.  `make fast` converges after two passes at 2495 pages
with zero undefined citations, zero undefined references, and zero
rerun requests.

A361 repairs `prop:ds-sdr`, the next exposed gravity licensing warning
after A360.  The old proposition had no claim status or licensing tag
and stated an SDR for \(C^{\mathrm{DS}}\) even though the displayed
maps only contract the associated-graded linear DS complex with
differential \(d_0\).  Since the following paragraph passes to
\(d_{\mathrm{BRST}}=d_0+\delta\), the old wording risked asserting a
full BRST deformation retract before the HPL perturbation datum had
been installed.  The formulas also used \(W=(k+2)T\) without naming the
critical-level exclusion \(k\ne -2\).

A361 introduces
\[
C^{\mathrm{DS}}_{\mathrm{lin}}
\]
as the \(\partial\)-stable differential-polynomial
\(\C[\partial]\)-module generated by the acyclic pairs
\((V,b)\), \((U,c)\) and the cohomology summand \(W\), with
\[
d_0(b)=V,\qquad d_0(U)=2c,\qquad d_0(W)=d_0(V)=d_0(c)=0.
\]
The proposition is retitled `Linear DS deformation retract;
\ClaimStatusProvedHere`, with licensing \(\alpha+\gamma\) via the
principal DS chart at \(k\ne -2\) and \(\hypAmbientWtCpl\).  Its
statement now proves the \(\C[\partial]\)-linear SDR identities
\[
pi=\operatorname{id},\qquad h^2=0,\qquad ph=0,\qquad hi=0,\qquad
\operatorname{id}-ip=d_0h+hd_0
\]
on \((C^{\mathrm{DS}}_{\mathrm{lin}},d_0)\), and explicitly says that
the full \(d_{\mathrm{BRST}}\) retract is the perturbed HPL object, not
this elementary linear retract.

A361 adds `virasoro_ds_linear_sdr_identities()`, an executable check
of \(h^2=0\), \(ph=0\), \(hi=0\), and
\(\operatorname{id}-ip=d_0h+hd_0\) on the five generators
\(W,V,b,U,c\).  It rejects \(k=-2\).  The source-level regression guard
in `test_gravity_3d_engine.py` now requires the linear theorem title,
\(C^{\mathrm{DS}}_{\mathrm{lin}}\), \(\C[\partial]\), \(k\ne -2\),
\(\hypAmbientWtCpl\), and absence of the retired full-complex wording.

A361 verification: `py_compile` passed on the touched Python files.
Focused pytest passed with `77 passed, 1 skipped` on
`compute/tests/test_gravity_3d_engine.py`.  Fixed-string scans over
active Vol II and adjacent Vol I/III surfaces found no live copy of the
old `SDR for the DS complex` theorem title outside the regression
guard.  `make verify-licensing` exits 0 with zero blockers and 108
warnings elsewhere; the next gravity warning is
`thm:ds-hpl-transfer`.  `make fast` converges after two passes at 2497
pages with zero undefined citations, zero undefined references, and
zero rerun requests.

A362 repairs `thm:ds-hpl-transfer`, the next exposed gravity licensing
warning after A361.  The old theorem asserted that the transferred
data \((H_{\mathrm{Vir}},\{m_n^{\mathrm{Vir}}\},
\Delta_z^{\mathrm{Vir}},r^{\mathrm{Vir}}(z))\) "is a dg-shifted
Yangian on the Virasoro algebra", while the following honest-gaps
remark said that sign verification and the specific dg-shifted-Yangian
axioms still remained.  The theorem therefore claimed more than its
proof supplied.  The same local surface still carried bare
\(r^{\mathrm{Vir}}(z)=(c/2)/z^3+2T/z\) copies after the A360
distinction between the pre-\(d\log\) Laplace/OPE kernel and the
collision residue.

A362 retitles the theorem as `DS-HPL Virasoro spectral transfer;
\ClaimStatusConditional`, with licensing
\(\alpha+\beta+\gamma+\delta\) via the principal DS chart at
\(k\ne -2\), the source affine dg-shifted-Yangian comparison,
\(\hypAmbientWtCpl\), and \(\hypKZSDR\).  The theorem now constructs
the homotopy-coherent Virasoro spectral package
\[
\mathsf Y_{\mathrm{Vir}}^{\mathrm{HPL}}
=
\bigl(H_{\mathrm{Vir}},\{m_n^{\mathrm{Vir}}\}_{n\ge2},
\Delta_z^{\mathrm{Vir}},
r_{\mathrm{Vir}}^{\mathrm{coll}}(z)\bigr),
\qquad
r_{\mathrm{Vir}}^{\mathrm{coll}}(z)=\frac{c/2}{z^3}+\frac{2T}{z}.
\]
It states the \(\Ainf\) Yang--Baxter relation as the
homotopy-coherent relation obtained by transferring coassociativity of
\(\Delta_z^{\mathrm{aff}}\), and explicitly says that a strict or
dg-shifted-Yangian presentation of \(H_{\mathrm{Vir}}\) requires the
additional rationality, translation-compatibility, and sign-normalised
higher-product verifications.

A362 rewrites the honest-gaps remark so it no longer contradicts the
theorem.  The remaining obligations are: closed sign normalisation of
the transferred products in degree \(\ge3\), translation compatibility
of the transferred spectral parameter action, and rational
strictification of \(r_{\mathrm{Vir}}^{\mathrm{coll}}\) as the
\(r\)-matrix of a dg-shifted Yangian rather than only the target twist
of a homotopy-coherent \(\Ainf\) package.  Nearby language in
`3d_gravity.tex` now says "transferred Virasoro spectral package" and
"homotopy-coherent Virasoro spectral package before
dg-shifted-Yangian strictification" instead of treating the HPL output
as an already proved gravitational Yangian.  Active Vol II copies in
`spectral-braiding-core.tex`, `bar-cobar-review.tex`,
`line-operators.tex`, and `thqg_holographic_reconstruction.tex` now use
\(r_{\mathrm{Vir}}^{\mathrm{coll}}\) and identify the cubic+simple
kernel as the post-\(d\log\)-absorption collision residue.

A362 adds `virasoro_ds_hpl_transfer_scope()` to
`compute/lib/gravity_3d_engine.py`, recording the conditional claim
status, the \(\alpha+\beta+\gamma+\delta\) licensing tags, the HPL
proved core, the collision kernel, and the unproved strict
dg-shifted-Yangian promotion obligations.  New source guards in
`test_gravity_3d_engine.py` require the theorem title, licenses,
\(\hypAmbientWtCpl\), \(\hypKZSDR\),
\(\mathsf Y_{\mathrm{Vir}}^{\mathrm{HPL}}\),
\(r_{\mathrm{Vir}}^{\mathrm{coll}}\), the negative "not a strict or
dg-shifted Yangian" clause, and the honest-gaps obligations; they
reject the retired theorem title and strict-Yangian conclusion.

A362 verification: `python3 -m py_compile` passed on the touched
Python files.  The configured pytest binary passed
`compute/tests/test_gravity_3d_engine.py` with `80 passed, 1 skipped`.
Fixed-string scans over active Vol II found no live copy of `is a
dg-shifted Yangian on the Virasoro algebra`, no live `transferred
dg-shifted Yangian`, no live `entire dg-shifted Yangian of gravity`,
and no active prose copy of
\(r^{\mathrm{Vir}}(z)=(c/2)/z^3+2T/z\) outside negative regression
guards.  `make verify-licensing` exits 0 with zero blockers and 107
warnings elsewhere; `thm:ds-hpl-transfer` is no longer a warning.
`make fast` converges after two passes at 2497 pages with zero
undefined citations, zero undefined references, and zero rerun
requests.

A363 repairs the coproduct-primitivity cluster following A362:
`prop:coproduct-degree2-vanishing` and
`thm:gravitational-primitivity`.  The old surface had no claim-status
or licensing tags and stated the all-degree theorem too strongly.  The
degree-two proof said the source tree acquires ghost degree \(-1\) from
an \(h\)-insertion, while the all-degree proof later counted source
trees as ghost degree \(+1\).  The theorem also asserted ordinary
primitive behavior for every element of \(W_k(\mathfrak g)\), thereby
collapsing generator-level transferred coproduct primitivity with a
Hopf-primitivity statement for arbitrary composites.

A363 retitles the base proposition as `Degree-two Virasoro generator
coproduct; \ClaimStatusConditional`, with licensing
\(\alpha+\gamma+\delta\) via the principal
\(\widehat{\mathfrak{sl}}_2\) DS chart at \(k\ne-2\),
\(\hypAmbientWtCpl\), and \(\hypKZSDR\).  It now states exactly
\[
\Delta_{z,2}^{\mathrm{Vir}}(T,T)=0,\qquad
\Delta_{z,1}^{\mathrm{Vir}}(T)=\tau_z(T)\otimes1+1\otimes T
\]
on the ghost-number-zero BRST projection surface, and explicitly says
this is a generator statement; composites are governed by transferred
\(\Ainf\)-morphism identities.

A363 retitles the all-degree theorem as `Principal-DS generator
coproduct primitivity; \ClaimStatusConditional`, with licensing
\(\alpha+\gamma+\delta\) via the principal DS chart, ghost-zero BRST
projection, \(\hypAmbientWtCpl\), and \(\hypKZSDR\).  It assumes the
signed ghost-defect condition and concludes only
\[
\Delta_{z,n}^{W}(W_s)=0\qquad(n\ge2)
\]
for chosen homogeneous principal DS cohomology generators \(W_s\).  It
does not assert ordinary Hopf-primitivity for arbitrary composites.
The new remark `rem:principal-ds-ghost-defect-obligation` isolates the
missing technical detail needed to remove the conditional status:
prove, with signs and degrees in the HPL morphism formula, that every
nonlinear source-side and target-side tree for a principal DS generator
has pre-projection tensor ghost bidegree different from \((0,0)\).

A363 propagates this corrected scope through the gauge-gravity
comparison, the \(\mathcal W_N\) extension, the non-principal/BP
discussion, the scope remark, the primitive-package worked example,
the introduction, and the stable \(W\)-algebra example.  The adjacent
Vol I bridge principles were also corrected: they now advertise the
conditional generator-level statement rather than unconditional
all-principal/all-degree primitivity, and the Vol I ordered-bar theorem
no longer claims `\ClaimStatusProvedHere` for the all-degree result.

A363 adds `principal_ds_coproduct_primitivity_scope()` to
`compute/lib/gravity_3d_engine.py`.  The new guard records the
conditional status, the degree-two Virasoro base case, the
\(\alpha+\gamma+\delta\) licenses, the required signed ghost-defect
lemma, and the non-assertion of ordinary Hopf-primitivity for
composites.  New tests in `test_gravity_3d_engine.py` require the
scoped theorem titles, licenses, \(\hypAmbientWtCpl\), \(\hypKZSDR\),
generator language, and signed ghost-defect condition; they reject the
old unindexed \(\Delta_{z,n}^W=0\), the `for every x` formula, the
retired theorem titles, and the "source trees have ghost \(+1\)" proof
slogan.

A363 verification: `python3 -m py_compile` passed on the touched
Python files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`82 passed, 1 skipped`.  Fixed-string scans over active Vol II and
adjacent Vol I/III surfaces found no remaining live copy of `kills
every HPL tree`, `annihilates every HPL tree`, `strictly primitive at
all degrees`, or `all principal DS reductions` outside negative guards.
`make verify-licensing` exits 0 with zero blockers and 105 warnings
elsewhere; the next gravity warning is `thm:brown-henneaux`.
`make fast` converges after two passes at 2497 pages with zero
undefined citations, zero undefined references, and zero rerun
requests.

A364 repairs `thm:brown-henneaux`, the next exposed gravity licensing
warning after A363.  The old theorem had no claim-status or licensing
tag and stated the Virasoro asymptotic-symmetry algebra without
separating the imported Brown--Henneaux/Witten physics theorem from the
monograph's chiral HT boundary chart.  Nearby remarks already said the
bar complex does not derive \(G_N\), but the theorem itself did not
carry that epistemic split.  Several active Vol II and adjacent Vol I
surfaces also still used bare \(G\) in the Brown--Henneaux
normalization.

A364 retitles the theorem as `Brown--Henneaux chiral boundary chart;
\ClaimStatusProvedElsewhere`, with licensing
\(\alpha+\beta+\gamma\) via \(\hypBHdict\).  The theorem now states
one oriented \(SL(2,\mathbb R)\) Chern--Simons factor at
\(k=\ell/(4G_N)\), the chiral central charge
\(c_{\mathrm{ch}}=6k=3\ell/(2G_N)\), and two commuting Virasoro copies
for full AdS\(_3\).  It cites `BrownHenneaux,Witten88` and explicitly
says this is an external asymptotic-symmetry/Chern--Simons comparison
theorem, not a theorem of the bar complex.

A364 propagates the same convention through the gravity chapter's
BTZ/Cardy derivation, conical-defect formulae, large-\(c\) scaling
remarks, the introduction, preface, universal-holography and landscape
surfaces, `thqg_perturbative_finiteness.tex`, and Vol I
`entanglement_modular_koszul.tex`.  The active scanned surfaces now use
\(G_N\) consistently for the Brown--Henneaux dictionary; the old
`3\ell/(2G)`, `\frac{3\ell}{2G}`, and `\ell/(4G)` forms survive only
as an intentional retired-string regression guard.

A364 adds `brown_henneaux_chiral_chart_scope()` to
`compute/lib/gravity_3d_engine.py`.  The new guard records
`ProvedElsewhere`, licensing \(\alpha+\beta+\gamma\), the
\(k=\ell/(4G_N)\) and \(c_{\mathrm{ch}}=6k=3\ell/(2G_N)\) dictionary,
two chiral Virasoro copies, and the non-derivation of both \(G_N\) and
a dynamical metric path integral from the bar complex.  New tests in
`test_gravity_3d_engine.py` require the theorem's status, license,
citations, chiral central charge, and `not a theorem of the bar
complex`; a source guard rejects bare-\(G\) Brown--Henneaux notation
across the active Vol II surfaces.

A364 verification: `python3 -m py_compile` passed on the touched Python
files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`85 passed, 1 skipped`.  Fixed-string scans over active Vol II and
adjacent Vol I surfaces found no live `3\ell/(2G)`,
`\frac{3\ell}{2G}`, `\ell/(4G)`, `G \to 0`, `G/\ell`,
`G^2/\ell^2`, or `8Gh` residue except the intentional retired-string
guard.  `make verify-licensing` exits 0 with zero blockers and 104
warnings elsewhere; the next gravity warning is `thm:gravity-mc`.
`make fast` converges after two passes at 2497 pages with zero
undefined citations, zero undefined references, and zero rerun
requests.  Scoped `git diff --check` passed in Vol II and in the
touched Vol I file.

A365 repairs `thm:gravity-mc`, the next exposed gravity licensing
warning after A364.  The old theorem had no claim-status or licensing
tag and still listed the first MC projection as `closed face =
\mathrm{Vir}_c`.  That conflicted with the repaired boundary/bulk
discipline: \(\mathrm{Vir}_c\) is the Brown--Henneaux boundary chart,
while the closed/bulk vertex is
\(\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)\) only after the
bulk--Hochschild comparison.  The theorem also said all five faces
were proved unconditionally at the abstract level, erasing the
physics-bridge BV/parametrix hypotheses, the completed ambient, and
the line-side Koszul-effectiveness condition.

A365 retitles the theorem as `Brown--Henneaux Virasoro MC bridge;
\ClaimStatusConditional`, with licensing
\(\alpha+\beta+\gamma+\delta+\varepsilon\) via
\(\hypBHdict+\hypAmbientWtCpl+\effKoszul\).  It now fixes
\(b_{\mathrm{BH}}\), \(k=\ell/(4G_N)\), and
\(c_{\mathrm{ch}}=6k=3\ell/(2G_N)\); assumes the exact
physics-bridge hypotheses; and states only the licensed projections:
boundary/chiral Virasoro face, completed open-colour line module
category, abstract spectral line braiding with
\[
r^{\mathrm{coll}}_{\mathrm{Vir}}(z)
= (c_{\mathrm{ch}}/2)/z^3+2T/z,
\]
bar-intrinsic Virasoro genus-shadow MC element, and PVA
\(\lambda\)-bracket shadow.

A365 adds the missing negative statement: the closed/bulk vertex
attached to this boundary chart is
\(\Zder^{\mathrm{ch}}(\mathrm{Vir}_{c_{\mathrm{ch}}})\) after the
bulk--Hochschild comparison, not the boundary Virasoro algebra.  The
theorem now explicitly excludes the concrete
\(\mathrm{Vir}_{26-c_{\mathrm{ch}}}\)-module model, the Virasoro
fusion-kernel realization, and the reading of BTZ geometries as MC
deformations from its proved content.

A365 rewrites `def:btz-deformation` as `Candidate BTZ deformation`.
A BTZ geometry is represented by a candidate MC deformation only after
the Brown--Henneaux, Cardy, and BTZ comparison data are supplied.  The
\(\kappa\)-variable BTZ remark and both preface overview copies now
use the same candidate/comparison-surface language.

A365 adds `gravitational_mc_bridge_scope()` to
`compute/lib/gravity_3d_engine.py`.  New tests require the theorem's
conditional status, five licensing tags, boundary/bulk separation,
physics-bridge hypotheses, collision kernel, and candidate BTZ
language; they reject the old `closed face = Vir_c`,
unconditional-items sentence, cubic-CS shortcut proof, and flat BTZ
deformation definition.

A365 verification: `python3 -m py_compile` passed on the touched Python
files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`88 passed, 1 skipped`.  Propagation scans over active Vol II,
compute, Vol I, and Vol III found no live `closed face =
\(\mathrm{Vir}_c\)`, no unconditional-items sentence, and no flat
`BTZ black holes as MC deformations` / `is a MC deformation` residue
outside negative guards and old audit notes.  `make verify-licensing`
exits 0 with zero blockers and 103 warnings elsewhere; the next
gravity warning is `prop:gravity-D-S-commute`.  `make fast` converges
after two passes at 2497 pages with zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped `git diff
--check` passed.

A366 repairs `prop:gravity-D-S-commute`, the next exposed gravity
licensing warning after A365.  The old proposition had no claim-status
or licensing tag and said that Verdier duality and the modular
\(S\)-transform commute on genus-one MC elements.  That was too broad:
the bare substitution \(\tau\mapsto-1/\tau\) is not a physical torus
partition-function theorem, and the Hodge line, Dedekind multiplier
line, and \(E_2\) anomaly data must be kept in the correct ambient.

A366 retitles the proposition as `Formal Verdier--$S$ square for the
genus-one shadow; \ClaimStatusConditional`, with licensing
\(\beta+\gamma\) via the Verdier level-swap comparison and
\(\hypAmbientWtCpl\).  It introduces \(S_{\mathrm{sh}}\) as pullback
on the completed genus-one shadow complex, with the Hodge line and
Dedekind multiplier line kept in the coefficient system.  The proved
display is the formal shadow square
\[
S_{\mathrm{sh}}\mathbb D
\Theta_{\mathrm{Vir}_c}^{(1),\mathrm{sh}}(\tau)
=
\mathbb D S_{\mathrm{sh}}
\Theta_{\mathrm{Vir}_c}^{(1),\mathrm{sh}}(\tau)
=
\Theta_{\mathrm{Vir}_{26-c}}^{(1),\mathrm{sh}}(-1/\tau),
\]
with degree-zero output
\(((26-c)/2)S_{\mathrm{sh}}(\omega_1)\).  The proof now explicitly
excludes the physical torus partition function and refers the
Dedekind multiplier and quasi-modular \(E_2\) anomaly to
`prop:genus1-modular-anomaly`.

A366 also rescopes the paired scalar Cardy expression.  It now appears
only after the physical modular-invariance and vacuum-dominance
hypotheses of `conj:gravity-cardy`, with \(0\le c\le26\) and shifted
energy \(\Delta\ge0\).  The later Vol II BTZ complementarity remark
and the adjacent Vol I `prop:ent-btz-complementarity` were propagated
to the same scalar/conditional/real-domain language, including the
comparison fixed point \(c=13\).

A366 adds `verdier_modular_s_shadow_scope()` to
`compute/lib/gravity_3d_engine.py`.  New guards in
`test_gravity_3d_engine.py` require the conditional status,
\(\beta+\gamma\) licensing, \(\hypAmbientWtCpl\),
\(S_{\mathrm{sh}}\), the Hodge and Dedekind multiplier coefficient
lines, the degree-zero output, the non-partition-function scope, and
the real-domain Cardy hypotheses; they reject the old proposition
title, old genus-one MC wording, and old flat paired Cardy formula.

A366 verification: `python3 -m py_compile` passed on the touched
Python files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`91 passed, 1 skipped`.  Fixed-string scans over active Vol II,
compute, and adjacent Vol I found no live old `commute on genus-$1$
MC elements`, no live `The commutativity produces a paired Cardy
formula`, and no live unscoped `After the Cardy/BTZ entropy dictionary
is imposed, the Koszul pair` wording outside negative guards.  `make
verify-licensing` exits 0 with zero blockers and 102 warnings
elsewhere; the next gravity warning is `prop:grav-open-six-requirements`.
`make fast` converges after two passes at 2497 pages with zero
undefined citations, zero undefined references, and zero rerun
requests.

A367 repairs `prop:grav-open-six-requirements`, the next exposed
gravity licensing warning after A366.  The old proposition had no
claim-status or licensing tag and treated the algebraic
`Perf-Vir_c` package with a primitive coproduct and \(r\)-matrix as
though it were the gravitational open-sector category itself.  It also
identified the central-extension projection \(\C[\![c]\!]\) too
closely with the bulk, asserted primitive coproducts by construction,
treated the \(\mathrm{Vir}_{26-c}\)-module line model as a requirement
rather than Conjecture~\ref{conj:gravity-line-identification}, and
read the vacuum character as a genus-one closed-string partition
function.

A367 renames the algebraic definition as `Brown--Henneaux algebraic
line test package`.  It now defines
\[
\cC_{\mathrm{grav}}^{\mathrm{test}}
=
\bigl(
\widehat{\mathrm{Perf}}_{\hypAmbientWtCpl}(\mathrm{Vir}_c),
\Delta_z^{\mathrm{Vir}},
r_{\mathrm{Vir}}^{\mathrm{coll}}(z)
\bigr)
\]
as a completed comparison model for the Brown--Henneaux boundary
chart, not as the full gravitational open-sector factorization
dg-category.  The preceding physical definition now only denotes the
boundary-condition perturbation category before the Virasoro-module
comparison is imposed.

A367 retitles the proposition as `Brown--Henneaux algebraic line test
package; \ClaimStatusConditional`, with licensing
\(\alpha+\beta+\gamma+\delta+\varepsilon\) via
\(\hypBHdict+\hypAmbientWtCpl+\effKoszul\).  The six components now
state exactly: the Brown--Henneaux boundary chart; the full bulk
\(\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)\) with \(\C[\![c]\!]\) only as
central-extension projection; generator-level primitive coproduct for
\(T\), with all-degree principal-DS primitivity conditional on the
signed ghost-defect lemma; the \(\mathrm{Vir}_{26-c}\)-module model as
conjectural; the annulus trace as the scalar one-boundary vacuum
character; and the genus tower as the algebraic scalar shadow whose
physical reading requires modular, vacuum-dominance, and saddle data.

A367 adds `brown_henneaux_line_test_package_scope()` to
`compute/lib/gravity_3d_engine.py`.  New guards in
`test_gravity_3d_engine.py` require the conditional status, all five
licensing tags, central-projection separation, generator-level
coproduct scope, conjectural line model, one-boundary trace scope, and
non-construction of the level-\(\mathsf A\) gravity-line operator
algebra.  The source guard rejects the old proposition title, the old
`Definition ... satisfies` sentence, `Primitive coproduct by
construction`, `the tensor product is trivial`, the full-open-sector
definition sentence, and the genus-one closed-string
partition-function wording.

A367 verification: `python3 -m py_compile` passed on the touched
Python files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`93 passed, 1 skipped`.  Fixed-string scans over active Vol II,
compute, and adjacent Vol I/III surfaces found no live copy of the old
proposition title, `Primitive coproduct by construction`, `the tensor
product is trivial`, or the genus-one closed-string partition-function
wording outside negative guards.  `make verify-licensing` exits 0
with zero blockers and 101 warnings elsewhere; the next gravity
warning is `thm:gravity-mc-package`.  `make fast` converges after two
passes at 2499 pages with zero undefined citations, zero undefined
references, and zero rerun requests.  Scoped `git diff --check`
passed.

A368 repairs `thm:gravity-mc-package`, the next exposed gravity
licensing warning after A367.  The old theorem had no claim-status or
licensing tag and identified \(\Theta_{\mathrm{grav}}\) with a sum
containing both genus-zero \(A_\infty\) operations
\(\sum_{k\ge2}\alpha_k\) and positive-genus terms.  This conflicted
with the Vol I bar-intrinsic theorem, where
\(\Theta_\cA=D_\cA-d_0\) is the positive-genus correction to the
genus-completed bar differential.  The same local surface also still
wrote the full chiral derived centre as \(\C[\![c]\!]\), collapsing
the full bulk object to its central-extension projection.

A368 repairs the adjacent bulk paragraph first: the full object is now
\(\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)\), with
\[
\pi_{\mathrm{cent}}\colon
\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)\to\C[\![c]\!]
\]
only the scalar central-extension projection.  The text explicitly
says this map is not an equivalence with the full derived chiral
centre.

A368 retitles the theorem as `Bar-intrinsic Virasoro MC shadow;
\ClaimStatusConditional`, with licensing
\(\alpha+\beta+\gamma+\varepsilon\) via
\(\hypBHdict+\hypAmbientWtCpl+\effKoszul\).  The theorem now defines
\[
D_{\mathrm{Vir}_c}
=d_0+\sum_{g\ge1}\hbar^g d_{\mathrm{Vir}_c}^{(g)},\qquad
\Theta_{\mathrm{grav}}
=D_{\mathrm{Vir}_c}-d_0
=\sum_{g\ge1}\hbar^g d_{\mathrm{Vir}_c}^{(g)},
\]
and states the MC equation
\([d_0,\Theta_{\mathrm{grav}}]+\frac12[\Theta_{\mathrm{grav}},
\Theta_{\mathrm{grav}}]=0\) in the completed coinvariant modular
convolution algebra.

A368 separates the projections.  Genus-zero operations \(m_k\) are
part of \(d_0\), not part of \(\Theta_{\mathrm{grav}}\).  The formula
\(F_g(\mathrm{Vir}_c)=\kappaChHodge(\mathrm{Vir}_c)
\lambda_g^{\mathrm{FP}}\) is only the proved uniform-weight scalar
trace; the full non-scalar genus component is the stable-graph
coderivation \(d_{\mathrm{Vir}_c}^{(g)}\).  The finite values
\(\kappaChHodge=c/2\), \(\mathfrak C_{\mathrm{Vir}}=-c\),
\(\mathfrak Q^{\mathrm{contact}}_{\mathrm{Vir}}=10/[c(5c+22)]\), and
generic \(o_5\ne0\) are finite projections of the MC element, not
substitutes for it.

A368 adds `virasoro_bar_intrinsic_mc_shadow_scope()` to
`compute/lib/gravity_3d_engine.py`.  New guards in
`test_gravity_3d_engine.py` require the conditional status,
\(\alpha+\beta+\gamma+\varepsilon\) licensing, the exact
\(D_{\mathrm{Vir}_c}-d_0\) definition, the exclusion of genus-zero
operations from \(\Theta_{\mathrm{grav}}\), the uniform-weight
scalar-lane scope, the non-scalar stable-graph component, the bulk
central-projection separation, and the non-physical-partition-function
scope.  Source guards reject the old theorem title, old `full modular
MC element` sentence, old \(\sum_{k\ge2}\alpha_k\) formula, old
\(D_A-d_0\) string, and full derived-centre equivalence.

A368 verification: `python3 -m py_compile` passed on the touched
Python files.  `/opt/homebrew/bin/pytest
compute/tests/test_gravity_3d_engine.py -q` passed with
`95 passed, 1 skipped`.  Fixed-string scans over active Vol II,
compute, and adjacent Vol I/III surfaces found no live old theorem
title, no live `full modular MC element`, no live
\(\sum_{k\ge2}\alpha_k\) MC formula, and no live full-derived-centre
equivalence outside negative guards.  `make verify-licensing` exits 0
with zero blockers and 100 warnings elsewhere; the next gravity
warning is `prop:genus0-product-decomposition`.  `make fast` converges
after two passes at 2499 pages with zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped
`git diff --check` passed.

A369 repairs `prop:genus0-product-decomposition`, the next exposed
gravity licensing warning after A368.  The old proposition had no
claim-status or licensing tag and stated the genus-zero operad as the
external product \(\SCchtop\times\Eone^{\mathsf{tr}}\).  That phrase
was too strong in the wrong direction: the construction permits mixed
\(\mathsf{ch}\)- and \(\mathsf{bdy}\)-inputs into
\(\mathsf{tr}\)-output, so the object is not the ordinary product of
two coloured operads on disjoint colour sets.

A369 retitles the proposition as `Directed genus-$0$ product
decomposition; \ClaimStatusProvedHere`, with licensing
\(\alpha+\gamma\) via the genus-zero HT slab chart and
\(\hypAmbientWtCpl\).  The statement now defines the strict directed
colour-filtered product
\[
\mathcal O^{\Ainf\text{-ch}}\big|_{g=0}
\cong
\SCchtop\boxtimes_{\mathrm{dir}}\Eone^{\mathsf{tr}},
\qquad
\mathsf{ch}\le\mathsf{bdy}\le\mathsf{tr},
\]
and lists the non-empty operation spaces
\(\FM_k(\C)\), \(\FM_k(\C)\times\Eone(m)\), and
\(\FM_k(\C)\times\Eone(m)\times\Eone(p)\).  It explicitly says that
mixed \(\mathsf{ch}/\mathsf{bdy}\) inputs into \(\mathsf{tr}\) remain
present.

A369 rewrites the proof through operation spaces, the chosen
genus-zero slab chart, compactified configuration-space splitting, and
componentwise FM/Stasheff insertion.  The Cauchy--Heaviside propagator
is retained as an analytic check, not as the definition of the
operad.  The algebra-over-operad remark now identifies the output as
an \(A_\infty\)-algebra object in the module category over the
\(\SCchtop\)-algebra \((A_{\mathsf{ch}},A_{\mathsf{bdy}})\), rather
than claiming that the transverse object simply carries an
\(\SCchtop\)-algebra structure.

A369 propagates
\(\SCchtop\boxtimes_{\mathrm{dir}}\Eone^{\mathsf{tr}}\) through the
genus-one construction, the curved-Dunn Kunneth paragraph, the
modular-operad definition, and the unitality proof.  The conclusion
table now names the genus-zero result as a directed product
decomposition.  Cross-volume scans found no live Vol I/III copy of the
old formula requiring mathematical propagation.

A369 adds `genus0_directed_product_decomposition_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and creates
`compute/tests/test_genus0_directed_product_iv.py`.  The IV witness
enumerates colour profiles from the order
\(\mathsf{ch}\le\mathsf{bdy}\le\mathsf{tr}\) and checks that mixed
\((\mathsf{ch},\mathsf{bdy},\mathsf{tr})\to\mathsf{tr}\) operations
are non-empty while the ordinary disjoint-colour product would make
them empty.

A369 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py
compute/tests/test_genus0_directed_product_iv.py -q -ra` passed with
`99 passed, 1 skipped`.  The local IV registry reports one
non-tautological entry for `prop:genus0-product-decomposition`.
`PATH=compute/.venv/bin:$PATH make verify-independence` imports
cleanly and reports zero tautological entries, but still exits 2 on
127 pre-existing orphan registry entries unrelated to A369.  `make
verify-licensing` exits 0 with zero blockers and 99 warnings
elsewhere.  `make fast` converges after two passes with zero undefined
citations, zero undefined references, and zero rerun requests.  Scoped
`git diff --check` passed.

A370 repairs `prop:qt-equivariance`, the next gravity equivariance
gap.  The old proposition had no claim-status or licensing tag and
asserted that quasi-triangularity plus YBE proves strict
\(\Aut(\Gamma)\)-equivariance for \(E_1\)-chiral quantum groups.  The
proof made two false moves: it treated a braided edge permutation as
literal commutation of monodromy products, and it derived
\(\mathrm{Mon}(R)_{21}=\mathrm{Mon}(R)^{-1}\) from
quasi-triangularity with \(\Delta=\id\).  The spectral-braiding core
already shows the correct relation:
\(R_{21}(-z)R_{12}(z)=\id\) is an inverse-braiding/orientation datum,
not a consequence of quasi-triangularity alone.

A370 rewrites the proposition as `Braided annular clutching is
graph-equivariant`, with `\ClaimStatusProvedHere` and licensing
\(\alpha+\beta+\gamma\) via \(\hypAmbientWtCpl\), the stable-graph
sewing orientation, and the quasi-triangular/braid-descent comparison
datum.  The theorem now assumes exactly four annular structures:
seamwise quasi-triangular comparison, YBE braid coherence, pure-braid
descent from braid lifts to the ordinary \(\Aut(\Gamma)\)-action after
sewing, and inverse orientation
\(\mathrm{Mon}(R)_{e^{\mathrm{op}}}=\mathrm{Mon}(R)_e^{-1}\).

A370 replaces the proof by the correct three-generator argument.
Endpoint exchange uses quasi-triangularity.  Parallel-edge
permutation uses positive braid lifts; YBE and distant-collar
commutativity give reduced-word coherence by Matsumoto's theorem, but
pure-braid descent is the additional datum needed for an honest graph
automorphism action.  Self-edge flips use inverse orientation and
explicitly do not follow from setting \(\Delta=\id\).

A370 propagates the repaired scope through the proved-loci remark,
equivariance-scope remark, Stokes--\(\FM_n\) hierarchy,
axiom-separation remark, modular-bar consequence paragraph,
`chapters/theory/introduction.tex`,
`chapters/theory/modular_swiss_cheese_operad.tex`, and
`standalone/stokes_gap_kzb_regularity.tex`.  Active-source scans now
find no live copy of the old proposition title, the old
unconditional-\(E_1\) sentence, the "QT+YBE is sufficient" shortcut,
or the false \(\Delta=\id\) derivation.

A370 adds `quantum_group_clutching_equivariance_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and creates
`compute/tests/test_qt_equivariance_annular_iv.py`.  The IV witness
checks the Artin braid presentation: YBE gives braid coherence
\(s_1s_2s_1=s_2s_1s_2\), but it does not impose the symmetric
relation \(s_i^2=1\), and inverse orientation is a separate inverse
braid datum.

A370 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py
compute/tests/test_qt_equivariance_annular_iv.py -q -ra` passed with
`101 passed, 1 skipped`.  `python3
compute/scripts/audit_independent_verification.py --show-orphans`
now recognizes `prop:qt-equivariance` as independently covered,
raises coverage to 135 entries, and reports zero tautological entries;
the global audit still exits 2 on 126 pre-existing orphan registry
entries unrelated to A370.  `make verify-licensing` exits 0 with zero
blockers and 98 warnings elsewhere.  `make fast` converges after two
passes with 2499 pages, zero undefined citations, zero undefined
references, and zero rerun requests.  Scoped `git diff --check`
passed.

A371 repairs `prop:modular-operad-unitality`, the next gravity
licensing warning after A370.  The old proposition had no
claim-status or licensing tag and proved all-genera unitality by
claiming that the unit vertex is a punctured sphere whose sewing seam
bounds a simply connected disk, hence has trivial \(R\)-monodromy.
That is the wrong mechanism.  The modular unit is an unstable
genus-\(0\) two-flag identity component adjoined to the stable-graph
category, and the seam is trivial because the annular sewing datum is
counit-normalised.

A371 rewrites the proposition as `Unit-normalised annular sewing is
unital at every genus`, with `\ClaimStatusProvedHere` and licensing
\(\alpha+\beta+\gamma\) via \(\hypAmbientWtCpl\), the unit-extended
stable-graph chart, and the counit-normalised annular sewing
comparison.  The theorem now assumes a strictly unital strongly
admissible \(\Eone\)-chiral algebra and a supplied annular clutching
datum normalised at the unit.  It names the three colour identities,
the unstable identity vertex, the diagonal bicomodule \(C_\Delta\),
the counit laws
\((\varepsilon\otimes\id)\Delta_z=\id\) and
\((\id\otimes\varepsilon)\Delta_z\simeq\tau_z\), and the
\(R\)-matrix unit normalisation
\((\varepsilon\otimes\id)R=1=(\id\otimes\varepsilon)R\).

A371 replaces the proof by a local cotensor calculation:
the exceptional edge sews the adjacent vertex comodule through the
diagonal bicomodule, and
\[
M\square_C C_\Delta\cong M\cong C_\Delta\square_C M.
\]
The \(R\)-matrix insertion is the identity by counit normalisation,
not by a simply-connected-disk argument.  Strict \(A_\infty\) units
kill higher trees at the identity component.  The infinite
class-\(\mathbf M\) tower can still act at the ordinary vertex, but
the exceptional edge adds no operation.  The theorem now explicitly
proves unitality of a supplied normalised sewing datum; it does not
construct or prove associativity of the general positive-genus
clutching maps.

A371 propagates this scope through the axiom-separation remark,
`chapters/theory/introduction.tex`,
`chapters/theory/modular_swiss_cheese_operad.tex`, and
`standalone/stokes_gap_kzb_regularity.tex`.  Active-source scans find
no live copy of the old simply-connected unit argument, the
contractible-loop monodromy phrase, or `unitality (vacuum axiom)` on
the touched surfaces.

A371 adds `modular_operad_unitality_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and creates
`compute/tests/test_modular_unitality_annular_iv.py`.  The IV witness
verifies the elementary coalgebra counit identities, that
augmentation-ideal \(R\)-terms vanish under either counit, and that
cotensoring with the diagonal bicomodule preserves the comodule basis.

A371 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py
compute/tests/test_genus0_directed_product_iv.py
compute/tests/test_qt_equivariance_annular_iv.py
compute/tests/test_modular_unitality_annular_iv.py -q -ra` passed with
`107 passed, 1 skipped`.  `make verify-independence` sees 967
ProvedHere labels and 136 independently covered claims, with zero
tautological entries; it still exits 2 on 126 pre-existing orphan
registry entries unrelated to A371.  `make verify-licensing` exits 0
with zero blockers and 97 warnings elsewhere.  `make fast` converges
after two passes with 2501 pages, zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped
`git diff --check` passed.

A372 repairs `prop:ainf-chiral-modular-bar-reduction`, the next
gravity licensing warning after A371.  The old proposition had no
claim-status or licensing tag and mixed the uncurved abstract
stable-graph theorem with the curved positive-genus nodal bar
calculation.  The abstract theorem `thm:modular-bar` proves
\(D^2=0\) for a modular bar datum whose internal differential is
square-zero.  The gravity proposition instead listed
\(d^2=\kappaChHodge(\cA)\omega_{g(v)}\) at positive genus as an
internal-differential axiom; that curvature is an obstruction term
until an \(S\)-tail or twisting datum absorbs it.

A372 rewrites the proposition as `Modular-bar criterion for annular
clutching data`, with `\ClaimStatusProvedHere` and licensing
\(\alpha+\beta+\gamma+\varepsilon\) via \(\hypAmbientWtCpl\), the
nodal chart, annular comparison maps, completed ambient, and the
\(S\)-tail curvature datum.  The theorem now states a criterion: a
supplied annular one-edge system gives a square-zero modular bar
differential only after four data are verified: an uncurved
square-zero internal differential or \(S\)-tail/twisting correction
\(d_S=d+D_S\) with \(d_S^2=0\); degree-\(+1\) one-edge expansion
maps; anticommutation with \(d_S\); and codimension-two cancellation.

A372 replaces the proof by the formal square expansion
\[
(D_{\mathrm{int}}+D_{\mathrm{exp}})^2
=D_{\mathrm{int}}^2+
[D_{\mathrm{int}},D_{\mathrm{exp}}]+D_{\mathrm{exp}}^2.
\]
At genus \(0\), the internal differential is flat.  At positive
genus, \(d^2=\kappaChHodge(\cA)\omega_g\) is a curvature computation,
not a modular-bar datum, until an \(S\)-tail/twisting term supplies
\(d_S^2=0\).  The proof now explicitly says the criterion does not
construct general positive-genus clutching maps or prove concrete
operadic associativity.

A372 also repairs the dependent proof of
`cor:affine-modular-bar-datum`: it now applies the criterion on the
covered affine loci.  Genus \(0\) is uncurved; at integrable level the
KZB regular-singular extension supplies the compensated sewing datum
used by the affine composition theorem.  The proof no longer says the
reduction proposition itself verifies all four axioms.

A372 adds `modular_bar_reduction_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and creates
`compute/tests/test_modular_bar_reduction_iv.py`.  The IV witness
checks the square expansion, shows raw curvature leaves a nonzero
obstruction, shows the compensating term restores square-zero, and
checks signed codimension-two cancellation.

A372 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py
compute/tests/test_genus0_directed_product_iv.py
compute/tests/test_qt_equivariance_annular_iv.py
compute/tests/test_modular_unitality_annular_iv.py
compute/tests/test_modular_bar_reduction_iv.py -q -ra` passed with
`110 passed, 1 skipped`.  `make verify-independence` sees 968
ProvedHere labels and 137 independently covered claims, with zero
tautological entries; it still exits 2 on 126 pre-existing orphan
registry entries unrelated to A372.  `make verify-licensing` exits 0
with zero blockers and 96 warnings elsewhere.  `make fast` converges
after two passes with 2501 pages, zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped
`git diff --check` passed.

A373 repairs `cor:affine-modular-bar-datum`, the next gravity
licensing warning after A372.  The corollary now has
`\ClaimStatusProvedHere` and licensing
\(\alpha+\beta+\gamma+\varepsilon\) via \(\hypAmbientWtCpl\).  Its
statement names the exact KZ/KZB covered loci: genus \(0\) at every
non-critical level, and arbitrary genus at integrable level
\(k\in\mathbb Z_{\ge0}\) in the semisimple integrable
positive-energy module category, with regular-singular KZB extension
and curvature compensation datum.  It explicitly asserts no modular
bar datum at generic non-integral level and genus \(g\ge1\).

A373 rewrites the proof through the A372 criterion.  At genus \(0\),
the internal differential is uncurved and KZ logarithmic flatness gives
the Drinfeld pentagon and hexagon.  At integrable level, the KZB
regular-singular extension supplies the \(S\)-tail/twisting correction
so \(d_S^2=0\), and the absence of Stokes-sector ambiguity gives the
codimension-two cancellation.  The following remark now says that the
corollary alone supplies only the modular bar datum and the
square-zero coalgebra; the full modular-operad axioms additionally
need braided-annular equivariance, unitality, pure-braid descent, and
inverse-orientation data.

A373 adds `affine_modular_bar_datum_scope()` to
`compute/lib/gravity_3d_engine.py`, source guards to
`compute/tests/test_gravity_3d_engine.py`, and creates
`compute/tests/test_affine_modular_bar_datum_iv.py`.  The IV witness
checks the affine level/genus predicate and the four modular-bar
criterion inputs, especially that generic non-integral positive genus
fails without the codimension-two/Stokes input.

A373 verification: `compute/.venv/bin/python -m pytest
compute/tests/test_gravity_3d_engine.py
compute/tests/test_modular_bar_reduction_iv.py
compute/tests/test_affine_modular_bar_datum_iv.py -q -ra` passed with
`109 passed, 1 skipped`.  `make verify-independence` sees 969
ProvedHere labels and 138 independently covered claims, with zero
tautological entries; it still exits 2 on 126 pre-existing orphan
registry entries unrelated to A373.  `make verify-licensing` exits 0
with zero blockers and 95 warnings elsewhere.  `make fast` converges
after two passes with 2501 pages, zero undefined citations, zero
undefined references, and zero rerun requests.  Scoped
`git diff --check` passed.

A374 repairs `thm:E3-topological-km`, the next gravity licensing
warning after A373.  The old theorem had no claim-status or
scanner-visible licensing tag and read as a strict theorem on
\(\Zder^{\mathrm{ch}}(V_k(\fg))\) itself.  The proof, however, only
supports the topologisation after passing to \(Q_{\mathrm{CS}}\)-
cohomology: the Sugawara antighost identity kills holomorphic
translations on cohomology, not as a strict raw-chain-level
operator identity on the original BV complex.

A374 retitles the theorem as a cohomological affine Kac--Moody
\(E_3\)-topological structure, adds `\ClaimStatusProvedHere`, and
introduces the scanner-visible hypothesis macro `\hypKMHTBV` with
licensing \(\alpha+\beta+\gamma+\varepsilon\).  The theorem now
targets
\[
H^\bullet_{Q_{\mathrm{CS}}}\Zder^{\mathrm{ch}}(V_k(\fg))
\]
and explicitly states that it does not assert a strict
raw-chain-level \(E_3^{\mathrm{top}}\)-structure before cohomology.
The proof now contains the cohomological translation calculation
\[
[\partial_z\cO]
=[T_{\mathrm{Sug},(0)}\cO]
=[[Q_{\mathrm{CS}},G_{\mathrm{Sug}}]_{(0)}\cO]
=0
\]
in \(Q_{\mathrm{CS}}\)-cohomology, followed by local constancy, the
single-colour external product, and Dunn additivity.  The strict
identity
\([Q_{\mathrm{CS}},\widetilde G_{\mathrm{Sug}}]=T_{\mathrm{Sug}}\)
on the original BV complex is deferred to
`rem:frontier-class-L-strict-chain-level`.

A374 propagates this cohomological scope through the introduction,
foundations, equivalence and locality summaries, and the heptagon collapse theorem.
The heptagon proof now says the \(E_3\)-topological faces are
equivalent on \(Q_{\mathrm{tot}}\)-cohomology and that Step (i)
supplies a cohomological Sugawara homotopy, not the excluded strict
operator-chain identity.  The PVA quantum-lift proof anchor now points
to `chapters/connections/3d_gravity.tex:thm:E3-topological-km`.

A374 adds `affine_e3_topological_km_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and tightens
`compute/tests/test_e3_topological_km.py` plus
`compute/tests/test_e3_topological_km_iv.py` so the independent
witnesses require non-criticality, \(Q_{\mathrm{CS}}\)-cohomology,
and exclusion of the strict raw-chain-level strengthening.

A374 verification: Python compilation of the edited compute files
passed.  Focused pytest on `compute/tests/test_gravity_3d_engine.py`,
`compute/tests/test_e3_topological_km.py`,
`compute/tests/test_e3_topological_km_iv.py`, and
`compute/tests/test_p3_pva_quantum_lift.py` passed with
`140 passed, 1 skipped`.  `make verify-independence` sees 969
ProvedHere labels and 121 independently covered claims, with zero
tautological entries; it still exits 2 because 104 test modules fail
to import in the make target and 86 pre-existing orphan registry
entries remain unrelated to A374.  The target label is decorated in
both `test_e3_topological_km.py` and
`test_e3_topological_km_iv.py`.  `make verify-licensing` exits 0
with zero blockers and 94 warnings elsewhere, and
`thm:E3-topological-km` is no longer warned.  `make fast` converges
after two passes with 2503 pages, zero undefined citations, zero
undefined references, and zero rerun requests.

A375 repairs the principal and good-graded Drinfeld--Sokolov
\(E_3\)-topologisation theorems.  The old proof treated Cartan
currents as \(Q_{\mathrm{CS}}\)-exact in the unreduced holomorphic
Chern--Simons bulk and used unverified standalone antighost
coefficients.  That was stronger than the mathematics supports.  The
correct object is the total DS differential
\[
Q_{\mathrm{DS}}=Q_{\mathrm{CS}}+Q_{\mathrm{red}},
\qquad
Q_{\mathrm{DS},f}=Q_{\mathrm{CS}}+Q_{\mathrm{red},f}
\]
on the DS bulk--boundary reduction complex.

A375 introduces the scanner-visible hypothesis macro `\hypDSBRST`.
The principal theorem `thm:E3-topological-DS` is now
`\ClaimStatusConditional`, licensed
\(\alpha+\beta+\gamma+\varepsilon\) via `\hypDSBRST`, and targets
\[
H^\bullet_{Q_{\mathrm{DS}}}\Zder^{\mathrm{ch}}(\cW).
\]
Its only non-formal input is the normalized cohomological identity
\[
T_{\mathrm{DS}}=[Q_{\mathrm{DS}},G'_{\mathrm{DS}}].
\]
The proof now derives topologisation by the cohomological translation
calculation and Dunn additivity.  It explicitly excludes unreduced
Cartan-current exactness and strict raw-chain-level topologisation.
The \(\mathfrak{sl}_2\) remark now treats the Cartan-antighost
coefficient as part of DS BRST normalization, not an invariant
standalone \(1/4\)-formula.

A375 rewrites `thm:E3-topological-DS-general` in the same conditional
form for good-graded \(f\), with target
\[
H^\bullet_{Q_{\mathrm{DS},f}}\Zder^{\mathrm{ch}}(\cW^k(\fg,f)).
\]
The non-principal obstruction remark now says the ghost filtration,
branched-cover descent, and primitive normalization are exactly the
content of `\hypDSBRST`, not automatic consequences.  The
Bershadsky--Polyakov remark now names the total differential
\(Q_{\mathrm{DS},f_{\min}}\) and removes the unverified explicit
antighost coefficient.

A375 propagates the total-DS-differential convention through the
introduction, preface, Virasoro and \(W_3\) examples, worked examples,
chiral-Higher-Deligne comparison, \(E_\infty\)-topologisation ladder,
fractional-ghost descent, and BP strict-chain comparison.  Stale
references to deleted labels such as `eq:current-Q-exact`,
`eq:T-imp-formula`, and `eq:G-prime-general` were removed.

A375 adds `principal_ds_e3_topological_scope()` and
`good_graded_ds_e3_topological_scope()` to
`compute/lib/gravity_3d_engine.py`, adds source guards to
`compute/tests/test_gravity_3d_engine.py`, and tightens the DS IV
predicates so they require non-criticality, the total DS differential,
the DS primitive package, cohomology, and exclusion of raw-chain-level
assertions.

A375 verification: Python compilation of the edited compute files
passed.  Focused pytest on `compute/tests/test_gravity_3d_engine.py`,
`compute/tests/test_climax_theorems_wave3_iv.py`,
`compute/tests/test_e3_topological_km.py`,
`compute/tests/test_e3_topological_km_iv.py`,
`compute/tests/test_e3_topological_ds_general.py`, and
`compute/tests/test_e3_topological_ds_general_iv.py` passed with
`124 passed, 1 skipped`.  `make verify-licensing` exits 0 with zero
blockers and 92 warnings elsewhere; both DS topologisation theorem
labels are no longer warned.  `make verify-independence` sees 970
ProvedHere labels and 122 independently covered claims, with zero
tautological entries; it still exits 2 because 104 test modules fail
to import in the make target and 81 pre-existing orphan registry
entries remain unrelated to A375.  `make fast` converges after two
passes with 2501 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A376 repairs `thm:genus1-mc-virasoro`, the next gravity licensing
warning after A375.  The old theorem conflated the scalar genus-one
Hodge trace, the vacuum-character Casimir term, and the genuine torus
Virasoro Ward blocks.  The scalar self-sewing is
\[
\Theta^{(1)}_{n=0}
=\kappaChHodge(\mathrm{Vir}_c)\omega_1
=(c/2)\omega_1,
\]
because \(T_{(3)}T=c/2\), and its Hodge trace is
\[
\kappaChHodge(\mathrm{Vir}_c)
\int_{\overline{\mathcal M}_{1,1}}\lambda_1=c/48.
\]
The vacuum-character one-point block instead has constant term
\[
\operatorname{CT}_q A^{(1)}_1(T;\tau)=-c/24,
\qquad
A^{(1)}_1(T;\tau)=
q\partial_q\log
\bigl(q^{-c/24}\prod_{m\ge2}(1-q^m)^{-1}\bigr).
\]

A376 introduces the scanner-visible hypothesis macro
`\hypVirTorusBlock`.  The theorem is now `\ClaimStatusConditional`,
licensed \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypVirTorusBlock`.  It states that the displayed
stable-graph expansion is a completed modular-bar expression: the
scalar self-sewing and universal Ward singularities are fixed, while
regular torus blocks and convergence of the full all-\(n\) series are
part of the torus-block hypothesis.

A376 also corrects the two-point central singularity.  The kernel is
not \((c/2)\wp(z,\tau)\).  With
\(P_2(z,\tau)=\wp(z;\tau)+E_2(\tau)/12\) and
\(P_4=(1/6)\partial_z^2P_2=z^{-4}+O(1)\), the repaired statement is
\[
\operatorname{Sing}_{z_1=z_2}A^{(1)}_2(T,T;z_{12},\tau)
=\frac{c}{2}P_4(z_{12},\tau)
+2P_2(z_{12},\tau)A^{(1)}_1(T;\tau).
\]
The regular term is now explicitly module-dependent torus-block data,
not a consequence of \(\kappaChHodge\).

A376 adds `genus1_virasoro_mc_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards to
`compute/tests/test_gravity_3d_engine.py`.  The new oracle records the
conditional status, \((\beta,\gamma,\delta)\) licensing, hypotheses
`hypAmbientWtCpl` and `hypVirTorusBlock`, the \(c=24\) checks
\(\kappa=12\), \(F_1=1/2\), and Ward constant \(-1\), plus negative
predicates excluding a \(\wp\)-central kernel, scalar determination of
regular torus blocks, and equality of the Casimir term with the Hodge
trace.

A376 verification: Python compilation of
`compute/lib/gravity_3d_engine.py` and
`compute/tests/test_gravity_3d_engine.py` passed.  Focused pytest on
`compute/tests/test_gravity_3d_engine.py`,
`compute/tests/test_genus_one_bridge.py`, and
`compute/tests/test_genus1_kappa_verification.py` passed with
`246 passed, 1 skipped`.  Stale-formula grep found no active
manuscript or compute hit for the retired \((c/2)\wp(z_{12};\tau)\),
\(E_2(\tau)T/24\), partition-counting, or genus-one free-energy
density formulas, except inside the new regression tests that forbid
them.  `make verify-licensing` exits 0 with zero blockers and 91
warnings elsewhere; `thm:genus1-mc-virasoro` is no longer warned and
the next gravity warning is `thm:genus1-amplitudes`.  `make
verify-independence` still exits 2 globally with 970 ProvedHere labels,
122 independently covered claims, zero tautological entries, 104 import
failures in the make target, and 81 pre-existing orphan registry
entries; A376 adds no ProvedHere claim.  `make fast` converges after
two passes with zero undefined citations, zero undefined references,
and zero rerun requests.

A377 repairs `thm:genus1-amplitudes`, the next gravity licensing
warning after A376.  The old theorem had no claim-status or
scanner-visible hypothesis tag.  It also omitted the simple-pole
derivative contribution in the three-point Virasoro Ward recursion.
The Virasoro OPE is
\[
T(z)T(w)\sim
\frac{c/2}{(z-w)^4}
+\frac{2T(w)}{(z-w)^2}
+\frac{\partial T(w)}{z-w}.
\]
The two-point simple-pole term vanishes because the one-point block is
insertion-position independent, but the three-point recursion must
contain
\[
P_1(z_{ij},\tau)\partial_{z_j}
A^{(1)}_2(T,T;z_j,z_k;\tau).
\]

A377 retitles the theorem as `Genus-$1$ Virasoro Ward amplitudes`, adds
`\ClaimStatusConditional`, and licenses it \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypVirTorusBlock`.  The theorem now works in the
completed genus-one ambient, fixes the torus Ward-block normalization,
and separates the scalar Hodge trace, determinant-line shadow, and
chiral Virasoro stress-tensor block.  The determinant-line derivative
is
\[
q\partial_q\log\eta^{-c/2}=-(c/48)E_2(\tau),
\]
while the generic Virasoro vacuum one-point block is
\[
q\partial_q\log
\bigl(q^{-c/24}\prod_{n\ge2}(1-q^n)^{-1}\bigr)
=-\frac{c}{24}+\sum_{n\ge2}\frac{nq^n}{1-q^n}.
\]

A377 adds \(P_1(z,\tau)=z^{-1}+O(z)\) to the kernel conventions and
repairs the three-point singular formula to
\[
\operatorname{Sing}_{z_i=z_j} A^{(1)}_3
=
\operatorname{Sing}_{z_i=z_j}
\left[
\frac{c}{2}P_4(z_{ij},\tau)A^{(1)}_1
+2P_2(z_{ij},\tau)A^{(1)}_2(z_{jk})
+P_1(z_{ij},\tau)\partial_{z_j}A^{(1)}_2(z_j,z_k)
\right].
\]
The regular two-point remainder is now explicitly torus-block data
supplied by `\hypVirTorusBlock`, regular in \(z_{12}\), and not
determined by \(\kappaChHodge\).  The theorem also excludes extra
minimal-model singular-vector quotients from the generic character
formula and says that the statement is chiral Ward-block data, not a
full non-chiral physical torus partition function.

A377 adds `genus1_virasoro_amplitudes_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards to
`compute/tests/test_gravity_3d_engine.py`.  The oracle records the
conditional status, \((\beta,\gamma,\delta)\) licensing,
determinant-line coefficient \(-\kappa/24\), generic vacuum Casimir
constant \(-c/24\), product start \(n=2\), the required
\(P_1\partial A_2\) term, and negative predicates excluding scalar
determination of regular blocks, minimal-model quotients, and a full
physical torus partition-function claim.

A377 propagates the coefficient distinction into Vol I compute
surfaces.  `../chiral-bar-cobar/compute/lib/conformal_bootstrap_mc_engine.py`
now separates the generic vacuum one-point constant \(-c/24\) from the
determinant-line \(E_2\)-coefficient \(-\kappa/24\), records product
start \(n=2\), and marks the one-point block as not a pure \(E_2\)
multiple.  `../chiral-bar-cobar/compute/lib/moonshine_kappa_resolution_engine.py`
now states that scalar genus-one curvature extracts \(T_{(3)}T=c/2\),
the determinant shadow has coefficient \(-c/48\), and the generic
Virasoro Ward block has Casimir constant \(-c/24\) plus oscillator
terms.

A377 verification: Python compilation passed for the edited Vol II
compute files and the two touched Vol I compute files.  Focused pytest
on `compute/tests/test_gravity_3d_engine.py`,
`compute/tests/test_genus_one_bridge.py`, and
`compute/tests/test_genus1_kappa_verification.py` passed with
`248 passed, 1 skipped`.  Active manuscript/compute grep found no old
pure \(A^{(1)}_1=-(c/24)E_2\) claim except inside the new regression
test that forbids it; remaining hits live only in quarantine/audit
artifacts.  `make verify-licensing` exits 0 with zero blockers and 90
warnings elsewhere; `thm:genus1-amplitudes` is no longer warned and the
next gravity warning is `prop:genus1-modular-anomaly`.  `make
verify-independence` still exits 2 globally with 970 ProvedHere labels,
122 independently covered claims, zero tautological entries, 104 import
failures in the make target, and 81 pre-existing orphan registry
entries; A377 adds no ProvedHere claim.  `make fast` converges after
two passes with 2501 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A378 repairs `prop:genus1-modular-anomaly`, the next gravity licensing
warning after A377.  The old proposition had no claim-status or
scanner-visible hypothesis tag.  It also made two false or unlicensed
identifications.  First, it said that the anomalous term
\[
12\tau/(2\pi i)
\]
in
\[
E_2(-1/\tau)=\tau^2E_2(\tau)+12\tau/(2\pi i)
\]
was \(\kappa^{-1}\) times the non-separating degeneration class.  The
coefficient is universal modular-form data; \(\kappa\) enters only
after multiplying by the determinant-line derivative
\[
q\partial_q\log\eta^{-\kappa}
=-\frac{\kappa}{24}E_2.
\]
Second, it identified the scalar non-holomorphic correction
\(-3/(\pi\operatorname{Im}\tau)\) with the bar-transgression generator
and asserted the unlicensed equality \(\eta^2=u=\kappa\omega_1\).

A378 introduces the scanner-visible hypothesis macro `\hypEisComp` for
the analytic comparison between the algebraic transgression class and
the non-holomorphic Eisenstein completion.  The proposition is now
`\ClaimStatusConditional`, licensed \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypEisComp`.  It states the standard modular-form
part unconditionally:
\[
\widehat E_2(\tau)=E_2(\tau)-3/(\pi\operatorname{Im}\tau)
\]
is modular of weight \(2\), while the comparison with the completed
bar-transgression ambient is exactly the extra datum `\hypEisComp`.

A378 rewrites the curvature comparison as
\[
d^2_{\bar B}=\kappaChHodge(\mathrm{Vir}_c)\omega_1,\qquad
q\partial_q\log\eta^{-\kappaChHodge(\mathrm{Vir}_c)}
=-\frac{\kappaChHodge(\mathrm{Vir}_c)}{24}E_2.
\]
The proposition now explicitly refuses to identify the scalar function
\(-3/(\pi\operatorname{Im}\tau)\) with the algebra generator itself and
refuses the unlicensed equality
\(\eta^2=\kappaChHodge(\mathrm{Vir}_c)\omega_1\).

A378 adds `genus1_modular_anomaly_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards to
`compute/tests/test_gravity_3d_engine.py`.  The oracle records the
universal \(E_2\) anomaly coefficient \(12\), `q_dq_log_eta=E2/24`,
determinant-line coefficient \(-\kappa/24\), bar curvature class
\(\kappaChHodge(\mathrm{Vir}_c)\omega_1\), and negative predicates
excluding a \(\kappa^{-1}\) anomaly, literal scalar-function =
transgression-generator identification, and \(\eta^2\)-equals-curvature
assertion.

A378 verification: Python compilation of the edited compute files
passed.  Focused pytest on `compute/tests/test_gravity_3d_engine.py`,
`compute/tests/test_genus_one_bridge.py`, and
`compute/tests/test_genus1_kappa_verification.py` passed with
`250 passed, 1 skipped`.  Stale-claim grep found the retired
\(\kappa^{-1}\) anomaly and \(\eta^2=u=\kappa\omega_1\) strings only
inside the new regression tests that forbid them.  `make
verify-licensing` exits 0 with zero blockers and 89 warnings elsewhere;
`prop:genus1-modular-anomaly` is no longer warned and the next gravity
warning is `thm:genus1-r-matrix`.  `make verify-independence` still
exits 2 globally with 970 ProvedHere labels, 122 independently covered
claims, zero tautological entries, 104 import failures in the make
target, and 81 pre-existing orphan registry entries; A378 adds no
ProvedHere claim.  `make fast` converges after two passes with 2501
pages, zero undefined citations, zero undefined references, and zero
rerun requests.

A379 repairs `thm:genus1-r-matrix`, the next gravity licensing warning
after A378.  The old theorem was already conditional on a KZB
heat-operator normalization, but it still stated a full modular
\(R\)-matrix expansion
\[
R^{\mathrm{mod}}(z;\hbar,\tau)
=R_0(z)+\hbar^2r_1(z;\tau)+O(\hbar^4)
\]
and wrote the contact term as the bare universal expression
\(2T\cdot E_2(\tau)\).  The proof then read as if Fay plus the heat
equation produced a \(\tau\)-dependent family of \(R\)-matrices.  The
available argument is narrower: it gives a connection-level KZB shadow
kernel after the theta/Kronecker heat-operator normalization is fixed.

A379 introduces the scanner-visible hypothesis macro `\hypKZBHeat`.
The theorem is now `Genus-$1$ Virasoro KZB shadow kernel`, with
`\ClaimStatusConditional` and licensing \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypKZBHeat`.  The full \(R^{\mathrm{mod}}\)
assertion is replaced by the shadow-connection expansion
\[
\nabla^{\mathrm{KZB},\mathrm{sh}}(z;\hbar,\tau)
=\nabla^{(0),\mathrm{sh}}(z)
+\hbar^2r^{\mathrm{KZB},\mathrm{sh}}_1(z;\tau)+O(\hbar^4).
\]
The repaired formula is
\[
r^{\mathrm{KZB},\mathrm{sh}}_1(z;\tau)
=\frac c2P_2(z,\tau)+\mathcal C_T^{\mathrm{KZB}}(\tau)T,
\qquad
P_2(z,\tau)=\wp(z;\tau)+E_2(\tau)/12,
\]
with \(\mathcal C_T^{\mathrm{KZB}}=2E_2\) only in the chosen
normalization.  The theorem now explicitly separates this scalar
elliptic \(P_2\)-part from the torus stress-tensor central singular
kernel \((c/2)P_4\), and it refuses full quantum \(R\)-matrix
construction, nonlinear quantum Yang--Baxter, and higher
\(\hbar^{2g}\) claims.

A379 also fixes two nearby summary leaks.  The later gravity summary
now names the genus-one KZB shadow kernel with scalar elliptic part
\((c/2)P_2\) plus the `\hypKZBHeat`-normalised contact term, instead of
advertising a genus-one \(R\)-matrix correction \(\sim\wp\).  The
Eisenstein-series summary now says that, under `\hypEisComp`,
\(\widehat E_2\) corresponds to the analytic image of the completed
transgression neutralisation of the bar curvature, not to a literal
transgression generator killing curvature.

A379 adds `genus1_virasoro_kzb_shadow_kernel_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards to
`compute/tests/test_gravity_3d_engine.py`.  It also repairs the older
KZB helper `compute/lib/genus1_intersection.py`, replacing the false
metadata phrase “heat equation for \(\wp\)” by the theta/Kronecker heat
equation producing \(\wp\)-terms, with a guard in
`compute/tests/test_genus1_intersection.py`.

A379 verification: Python compilation passed for the edited compute and
test files.  The new focused guards passed with `4 passed`; the broader
genus-one gravity test surface passed with `279 passed, 1 skipped`.
Cross-volume stale-claim grep over active Vol II plus Vol I/III
chapter/compute surfaces found no live old copies of the old theorem
title, \(2T\cdot E_2\), self-energy wording, bare-\(\wp\) heat equation,
or unlicensed transgression-element sentence.  `make verify-licensing`
exits 0 with zero blockers and 88 warnings elsewhere;
`thm:genus1-r-matrix` is no longer warned and the next gravity warning
is `thm:gravity-weinberg-from-ainf`.  `make verify-independence` still
exits 2 globally with 970 ProvedHere labels, 122 independently covered
claims, zero tautological entries, 104 import failures in the make
target, and 81 pre-existing orphan registry entries; A379 adds no
ProvedHere claim.  `make fast` converges after two passes with 2501
pages, zero undefined citations, zero undefined references, and zero
rerun requests.  `git diff --check` is clean on the touched Vol II
files.

A380 repairs `thm:gravity-weinberg-from-ainf`, the next gravity
licensing warning after A379.  The old soft-theorem band identified the
degree-\(2\) Virasoro shadow residue with the physical Weinberg theorem
without naming the celestial comparison datum.  Its dictionary also
used the wrong OPE-mode shift:
\[
S_i^{(p)}=T_{(p+1)}\mathcal O_i
=\operatorname{Res}(z-z_i)^{p+1}T(z)\mathcal O_i ,
\]
although the residue convention gives
\[
W_i^{(j)}=T_{(j)}\mathcal O_i
=\operatorname{Res}(z-z_i)^jT(z)\mathcal O_i\,dz .
\]
The proof then computed only the simple-pole residue
\(T_{(0)}\mathcal O_i=\partial_{z_i}\mathcal O_i\), claimed that
translation invariance produced Weinberg's factor, and displayed the
photon-like expression \((p_i\cdot\epsilon)/(p_i\cdot q)\).

A380 introduces `\hypCelSoft` and rewrites the theorem as
`Degree-$2$ Virasoro Ward residue and Weinberg comparison`, with
`\ClaimStatusConditional` and licensing \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypCelSoft`.  The algebraic content is now the
primary Ward kernel
\[
S_2(z_0)=\sum_i\left(\frac{h_i}{(z_0-z_i)^2}
+\frac{\partial_{z_i}}{z_0-z_i}\right),
\qquad
W_i^{(0)}=\partial_{z_i},\quad W_i^{(1)}=h_i.
\]
The physical statement is separated: under
\(\mathsf{CelSoft}_{\hypCelSoft}\), the BMS-normalised image of the
simple-pole package is the spin-two Weinberg factor
\[
S_{\mathrm{Wein}}^{(0)}
=\sum_i\frac{\varepsilon_{\mu\nu}p_i^\mu p_i^\nu}{q\cdot p_i},
\]
with the gravitational coupling absorbed into the comparison datum.
Without `\hypCelSoft`, the theorem explicitly asserts only the
Virasoro Ward-residue identity.  The companion proof in
`thqg_soft_graviton_theorems.tex` now uses the same spin-two comparison
phrase.

A380 adds `gravity_weinberg_ward_residue_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards in
`compute/tests/test_gravity_3d_engine.py`.  The guards prevent return
of the \(T_{(p+1)}\) shift, the photon-like factor, the translation-
invariance derivation, the old theorem title, and the unlicensed
infinite-soft-identity sentence.

A380 verification: Python compilation passed for the edited gravity
engine and test.  The focused guard passed with `3 passed`; the broader
gravity engine surface passed with `124 passed, 1 skipped`.
Cross-volume stale-claim grep over active Vol II plus Vol I/III
chapter/compute surfaces found no live old copies of the retired mode
shift, photon-like factor, old theorem title, unlicensed infinite-soft
sentence, or companion \(p_i\cdot\varepsilon_s/(p_i\cdot q_s)\) phrase.
`make verify-licensing` exits 0 with zero blockers and 87 warnings
elsewhere; the next gravity warning is `thm:gravity-cs-from-ainf`.
`make fast` converges after two passes with 2501 pages, zero undefined
citations, zero undefined references, and zero rerun requests.  `git
diff --check` is clean on the touched Vol II files.

A381 repairs `thm:gravity-cs-from-ainf`, the next gravity licensing
warning after A380.  The old theorem identified the double-pole residue
of the degree-\(2\) Virasoro Ward kernel,
\(h_i/(z_0-z_i)^2\), with the physical Cachazo--Strominger subleading
soft theorem.  Its proof computed only the conformal-weight term and
then equated the global conformal Ward identity with the physical
momentum-space theorem.  That omitted the simple-pole orbital term, the
anti-holomorphic completion, the Mellin soft-energy denominator, and
the degree-\(3\) shadow channel used by the celestial hierarchy.

A381 retitles the theorem as `Degree-$2$ conformal Ward package and the
Cachazo--Strominger comparison boundary`, with
`\ClaimStatusConditional` and licensing \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypCelSoft`.  The algebraic content is now
\[
\mathcal W_Y^{(2),\mathrm{hol}}
=\sum_i\bigl(Y(z_i)\partial_{z_i}+h_i\,\partial Y(z_i)\bigr),
\]
derived from \(W_i^{(0)}=\partial_{z_i}\) and \(W_i^{(1)}=h_i\).
The theorem states that the double-pole residue is only the weight
term.  Under \(\mathsf{CelSoft}_{\hypCelSoft}\), after
anti-holomorphic completion and the degree-\(3\) shadow channel are
included, the superrotation Ward package is read as
\[
S_{\mathrm{CS}}^{(1)}
=\sum_i\frac{\varepsilon_{\mu\nu}p_i^\mu q_\rho J_i^{\rho\nu}}
{q\cdot p_i},
\]
with coupling and phase conventions absorbed into the comparison
datum.

A381 also repairs the local hierarchy table so \(S^{(1)}\) is attached
to \(S_3\), while \(S_2|_{\mathrm{double}}\) is only the weight term.
The split movement summary in
`chapters/connections/thqg_3d_gravity_movements_vi_x.tex` now says the
Cachazo--Strominger channel comes from the degree-\(3\) shadow together
with the \(S_2\) weight term, and the companion theorem in
`thqg_soft_graviton_theorems.tex` now uses the spin-two contraction
\(\varepsilon_{\mu\nu}p_i^\mu q_{s,\rho}J_i^{\rho\nu}\).

A381 adds `gravity_cachazo_strominger_ward_package_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards in
`compute/tests/test_gravity_3d_engine.py`.  The guards prevent return
of the old theorem title, the double-pole-as-Cachazo--Strominger
slogan, the stale hierarchy row, and the vector-like physical
contraction.

A381 verification: Python compilation passed for the edited gravity
engine/test.  The focused guard passed with `3 passed`; the full
gravity engine surface passed with `127 passed, 1 skipped`.
Cross-volume stale-claim grep over active Vol II plus Vol I/III
chapter/compute surfaces found no live old copies of the retired
double-pole soft-theorem claim, old theorem title, stale hierarchy row,
global-Ward-to-physical theorem sentence, or vector-like
\(q_s^\mu J^i_{\mu\nu}\varepsilon_s^\nu\) contraction.  Remaining
formula echoes are legitimate BPZ/conformal-block equations or negative
source guards.  `make verify-licensing` exits 0 with zero blockers and
86 warnings elsewhere; the next gravity warning is
`thm:gravity-chy-from-ainf`.  `make fast` converges after two passes
with 2501 pages, zero undefined citations, zero undefined references,
and zero rerun requests.  `git diff --check` is clean on the touched
Vol II files.

A382 repairs `thm:gravity-chy-from-ainf`, the next gravity licensing
warning after A381.  The old theorem made the ternary \(\Ainf\)
operation \(m_3\) appear to be the algebraic source of the
Cachazo--He--Yuan sub-subleading soft theorem and described it as the
algebraic incarnation of the physical \(J^{\mu\nu}J^{\rho\sigma}\)
factor.  That conflated the cubic Stasheff precursor, the degree-\(4\)
quartic contact channel, and the celestial physical comparison.

A382 retitles the theorem as `Quartic contact sub-subleading shadow and
Cachazo--He--Yuan comparison`, with `\ClaimStatusConditional` and
licensing \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypCelSoft`.  The theorem now states that
\(T_{(2)}\mathcal O_i=0\) for primary insertions, so there is no binary
sub-subleading residue.  The cubic element contributes through the
degree-\(4\) self-interaction bracket
\([\Sh(\mathfrak C),\Sh(\mathfrak C)]\), while the first non-trivial
Virasoro graviton channel is the quartic contact coefficient
\[
\mathfrak Q^{\mathrm{contact}}_{\mathrm{Vir}}
=\frac{10}{c(5c+22)}\sim \frac{2}{c^2}.
\]
The algebraic sub-subleading shadow is now
\[
S^{(2)}
=
\Sh_{0,n}(\mathfrak Q(\mathrm{Vir}_c))\big|_{\mathrm{soft}}
+
[\Sh_{0,n}(\mathfrak C(\mathrm{Vir}_c)),
\Sh_{0,n}(\mathfrak C(\mathrm{Vir}_c))]\big|_{\mathrm{soft}},
\]
and it becomes the CHY/Hamada--Shiu--Li--Strominger package only after
the \(\mathsf{CelSoft}_{\hypCelSoft}\) comparison datum is applied.
The local table now records \(S^{(2)}\) as \(S_4\) plus
\([S_3,S_3]\), with \(S_3\) exact for Virasoro, and the operation side
as \(m_4\) plus \([m_3,m_3]\).

A382 adds `gravity_chy_quartic_contact_scope()` to
`compute/lib/gravity_3d_engine.py` and source guards in
`compute/tests/test_gravity_3d_engine.py`.  The guards prevent return
of the old theorem title, the raw-\(m_3\)-as-CHY slogan, the stale
hierarchy row, and the old physical \(J J\) incarnation phrase.

A382 verification: Python compilation passed for the edited gravity
engine/test.  The focused guard passed with `2 passed`; the full
gravity engine surface passed with `129 passed, 1 skipped`.  Stale
greps found no live old copies of the retired CHY/ternary phrases; the
surviving occurrences are negative source guards or an unrelated
spectral \(R\)-matrix proposition.  `make verify-licensing` exits 0
with zero blockers and 85 warnings elsewhere; the next gravity warning
is `thm:gravity-infinite-soft-tower`.  The first `make fast` attempt
was host-killed with exit 137 during pass 2, without a TeX diagnostic;
a clean rerun converged after two passes with 2503 pages, zero
undefined citations, zero undefined references, and zero rerun
requests.  `git diff --check` is clean on the touched Vol II files.

A383 repairs `thm:gravity-infinite-soft-tower`, the next gravity
licensing warning after A382.  The old theorem said soft order \(p\)
was controlled by shadow degree \(r=p+2\) via the raw operation
\(m_{p+2}\), and that class-\(\mathbf M\) produced infinitely many
independent algebraic shadow classes.  That collapsed normalized
shadow projections, lower Maurer--Cartan bracket terms, raw transferred
operations, and physical higher-soft interpretations.

A383 retitles the theorem as `Infinite Virasoro soft-shadow hierarchy
on the scalar metric branch`, with `\ClaimStatusConditional` and
licensing \(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypCelSoft`.  The theorem now fixes the completed
class-\(\mathbf M\) ambient and the scalar branch
\[
H(t)=t^2\sqrt{Q_{\mathrm{Vir}}(t)},\qquad
Q_{\mathrm{Vir}}(t)=(c+6t)^2+\frac{80t^2}{5c+22},
\quad \sqrt{Q_{\mathrm{Vir}}(0)}=c.
\]
The degree-\(r\) channel is the normalized projection
\(\Sh_{0,n}(\Theta_r)\) plus lower Maurer--Cartan brackets of total
degree \(r\), not raw \(m_r\) alone.  The scalar coefficients are
\(S_r=[t^r]H(t)/r\), and non-termination follows because
\(Q_{\mathrm{Vir}}(t)\) is not a square away from \(c=-22/5\).  The
statement explicitly says the coefficients are coupled by the MC
equation and are not independent primitive generators.

A383 also makes the quintic channel exact:
\[
S_5=-\frac{48}{c^2(5c+22)},\qquad
[\mathfrak C,\mathfrak Q]_H\text{ sub-channel}
=\frac{20}{c^2(5c+22)}.
\]
The full degree-\(5\) package is a conditional \(S^{(3)}\) soft-shadow
channel only after \(\mathsf{CelSoft}_{\hypCelSoft}\) is applied; the
\(O(G_N^3)\) reading now explicitly uses the Brown--Henneaux
substitution \(c=3\ell/(2G_N)\).  The local soft table now uses
normalized shadow projections plus lower brackets, and the split
soft-graviton file now treats \(r^\Theta_{\max}=\infty\) as algebraic
support-depth non-termination rather than a count of physical Ward
identities.

A383 adds `virasoro_shadow_metric_coefficients()` and
`gravity_infinite_soft_shadow_hierarchy_scope()` to
`compute/lib/gravity_3d_engine.py`, with source guards in
`compute/tests/test_gravity_3d_engine.py`.

A383 verification: Python compilation passed for the edited gravity
engine/test.  The focused guard passed with `4 passed`; the full
gravity engine surface passed with `133 passed, 1 skipped`.  The
celestial OPE bridge guard passed with `83 passed`.  Stale greps found
no live old copy of the raw-\(m_{p+2}\) control sentence, stale
soft-table rows, `Non-trivial for p even`, or the split-file physical
Ward-identity overclaim; surviving hits are negative guards or
historical audit notes.  `make verify-licensing` exits 0 with zero
blockers and 84 warnings elsewhere; the next gravity warning is
`thm:stasheff-cancellation`.  `make fast` converged after two passes
with 2503 pages, zero undefined citations, zero undefined references,
and zero rerun requests.  `git diff --check` is clean on the touched
Vol II files.

A384 repairs `thm:stasheff-cancellation` and
`prop:denominator-formula`, the next gravity licensing targets after
A383.  The old theorem claimed that Kac-table poles in Shapovalov
inverses cancel in the full \(\Ainf\) computation / total \(m_k\).
That contradicted the nearby \(S_5\) computation, where raw \(m_5\) is
not the normalized shadow coefficient and the all-arity projection
\(\Pi_{\mathrm{sh}}^{(r)}\) is not constructed.  The denominator
proposition also used literal `denom(S_r)`, although direct expansion
has rational scalar denominators \(3,7,11,\ldots\) in some
coefficients.

A384 introduces `\effScalarShadowProj` and retitles the theorem as
`Scalar shadow pole cancellation by Stasheff projection`, with
`\ClaimStatusConditional` and licensing \(\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effScalarShadowProj`.  The theorem now states only
the true pole-divisor result: after normalized scalar shadow projection,
\[
S_r=[t^r]\,t^2\sqrt{Q_{\mathrm{Vir}}(t)}/r,\qquad
Q_{\mathrm{Vir}}(t)=(c+6t)^2+\frac{80t^2}{5c+22},
\]
has nonconstant central-charge pole divisor
\[
c^{r-3}(5c+22)^{\lfloor(r-2)/2\rfloor}.
\]
Rational scalar denominators are excluded from that divisor.  The text
now says raw HPL summands and raw transferred operations may still
contain Shapovalov/Kac factors, and that the theorem does not construct
the all-arity operator-level projection to the scalar branch.

A384 also rewrites the planarity remark as a projected topological
recursion analogy, retitles the proposition as `Pole-divisor formula
for scalar shadow coefficients`, replaces `denom(S_r)` by
\(\operatorname{PoleDiv}_{\C(c)}(S_r)\), and scopes the
Shapovalov/bootstrap paragraph to scalar shadow projection rather than
unprojected \(m_k^H\).

A384 adds `virasoro_shadow_metric_pole_profile()` and
`stasheff_scalar_pole_cancellation_scope()` to
`compute/lib/gravity_3d_engine.py`, with source guards in
`compute/tests/test_gravity_3d_engine.py`.

A384 verification: Python compilation passed for the edited gravity
engine/test.  The focused guard passed with `3 passed`; the full
gravity engine surface passed with `136 passed, 1 skipped`.  Direct
metric-branch expansion through \(r=12\) verified the nonconstant pole
divisor and separated rational scalar denominators.  Stale greps found
no live old copy of the full-\(\Ainf\) cancellation claim, old
`total m_k` phrase, or literal Kac-denominator/non-planar tree
amplitude claim outside negative guards.  `make verify-licensing`
exits 0 with zero blockers and 82 warnings elsewhere; the next gravity
warning is `thm:shadow-closed-form`.  `make fast` converged after two
passes with 2503 pages, zero undefined citations, zero undefined
references, and zero rerun requests.  `git diff --check` is clean on
the touched Vol II files.

A385 repairs `thm:shadow-closed-form`, the next gravity licensing
warning after A384.  The old theorem stated the correct Catalan
closed-form expression but did not carry claim status or licensing and
did not identify the object as the normalized scalar \(T\)-line shadow
coefficient after `\effScalarShadowProj`.  Its proof also used the
false identity
\[
\binom{1/2}{m}=(-1)^{m-1}C_{m-1}/(4^m m),
\]
which fails already at \(m=1\).

A385 retitles the theorem as `Closed-form normalized scalar shadow
coefficient formula`, with `\ClaimStatusConditional` and licensing
\(\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effScalarShadowProj`.  The statement now fixes the
scalar branch
\[
Q_{\mathrm{Vir}}(t)=(c+6t)^2+\frac{80t^2}{5c+22},\qquad
H(t)=t^2\sqrt{Q_{\mathrm{Vir}}(t)},\qquad
\sqrt{Q_{\mathrm{Vir}}(0)}=c,
\]
excludes \(c=0,-22/5\), and states only
\[
S_r=[t^r]H(t)/r
=\frac{(-6)^{r-4}D}{2r\,c^{r-3}}F_r(D/144).
\]
It explicitly says this is not a raw \(m_r\) formula and not a full
ordered-bar invariant before scalar projection.  The proof now uses
the branch substitution \(u=-6t/c\), \(x=D/144\), so
\[
Q_{\mathrm{Vir}}(t)/c^2=(1-u)^2+4xu^2,
\]
and the coefficient follows from the shape-factor generating function.
The corrected half-binomial identity is
\[
\binom{1/2}{m}=(-1)^{m-1}C_{m-1}/2^{2m-1}.
\]

A385 adds `virasoro_catalan_shape_factor()`,
`virasoro_shadow_closed_form_coefficient()`, and
`shadow_closed_form_scope()` to `compute/lib/gravity_3d_engine.py`,
with `TestShadowClosedFormScope` in
`compute/tests/test_gravity_3d_engine.py`.  The companion script
`compute/ordered_e1_shadow_catalan.py` now describes normalized scalar
\(T\)-line coefficients and uses nonconstant pole-divisor language.
The WN Catalan script now separates the proved \(N=2,3\) closed forms
from the conditional \(N\ge4\) pattern.

A385 verification: Python compilation passed for the touched compute
files.  The focused guard passed with `6 passed`; the full gravity
engine surface passed with `142 passed, 1 skipped`.  Stale fixed-string
sweeps found no live old theorem title, false \(4^m m\)-denominator
identity, exact-denominator claim, or universal WN theorem overclaim in
the touched source files.  `make verify-licensing` exits 0 with zero
blockers and 81 warnings elsewhere; the next gravity warning is
`thm:w3-w-line-closed-form`.  `make fast` converged after two passes
with 2503 pages, zero undefined citations, zero undefined references,
and zero rerun requests.

A386 repairs `thm:w3-w-line-closed-form`, the next gravity licensing
warning after A385.  The theorem lacked claim-status/licensing data,
used the false half-binomial identity
\[
\binom{1/2}{j}=(-1)^{j-1}C_{j-1}/(4^j j),
\]
and asserted raw \(c\mapsto100-c\) complementarity for
\(S_{2r}^W\), which is false already for \(r=2\).  A stale
multichannel compute witness also used the wrong W-line square-root
branch \(c/3\), `gamma=61440/(c^2(5c+22)^3)`, and
\(S_2^{WW}=c/6\).

A386 retitles the theorem as the normalized scalar \(W\)-line theorem,
with `\ClaimStatusConditional` and licensing
\(\alpha+\gamma+\varepsilon\) via the \(W\)-line branch,
`\hypAmbientWtCpl+\effScalarShadowProj`.  The statement now fixes
\[
Q_W(t)=\frac{4c^2}{9}\left(1+\delta_3t^2\right),\qquad
\delta_3=\frac{122880}{c^2(5c+22)^3},\qquad
\sqrt{Q_W(0)}=2c/3,
\]
with \(c\ne0,-22/5\), and gives
\[
S_2^W=c/3,\qquad S_{2r+1}^W=0,\qquad
S_{2r}^W=
\frac{(-1)^rC_{r-1}30720^{r-1}}
     {3(2r-3)c^{2r-3}(5c+22)^{3(r-1)}}.
\]
The proof now derives this from
\(H_W(t)=(2c/3)t^2\sqrt{1+\delta_3t^2}\) and
\(\binom{1/2}{n}=(-1)^{n-1}C_{n-1}/2^{2n-1}\).
The theorem explicitly excludes the full two-variable \(\mathcal W_3\)
shadow tensor and unprojected transferred operations.

A386 replaces the false raw complementarity by the pole-cleared
constant
\[
N_{2r}^W=S_{2r}^Wc^{2r-3}(5c+22)^{3(r-1)}
=\frac{(-1)^rC_{r-1}30720^{r-1}}{3(2r-3)}.
\]
It adds `w3_w_line_shadow_coefficient()` and
`w3_w_line_closed_form_scope()` to `compute/lib/gravity_3d_engine.py`
with `TestW3WLineClosedFormScope`, and updates
`compute/w3_shadow_coefficients.py`, `compute/w3_shadow_closed_form.py`,
`compute/w3_multichannel_shadow.py`, and
`compute/w3_multichannel_shadow_results.json` to the \(2c/3\) branch
and \(\delta_3=122880/(c^2(5c+22)^3)\).

A386 verification: Python compilation passed for the touched compute
files.  The focused W-line guard passed with `6 passed`; the full
gravity engine surface passed with `148 passed, 1 skipped`.  Direct
comparison of the multichannel \(WW\) coefficients against
`w3_w_line_shadow_coefficient()` was zero through arity 14.  The W3
diagnostic scripts both ran cleanly.  Stale fixed-string sweeps found
no live old W-line title, false half-binomial denominator,
raw-complementarity slogan, \(Q^{WW}=(c/3)^2\) branch, or
\(S_2^{WW}=c/6\) line in the touched sources.  `make
verify-licensing` exits 0 with zero blockers and 80 warnings
elsewhere; the next gravity warning is `thm:catalan-dynkin-parity`.
`make fast` converged after two passes with 2505 pages, zero undefined
citations, zero undefined references, and zero rerun requests.

A387 repairs `thm:catalan-dynkin-parity`, the next gravity licensing
warning after A386.  The theorem had no claim-status/licensing tag and
asserted parity for every one-generator class-\(\mathbf M\) algebra
whose generator closes under \(m_2\).  The proof actually uses the
Virasoro-normalized root-killing point \(\varphi_j(-j)=0\), so an
arbitrary conformal weight \(s\) or a different \(L_0\)-string is not
covered.  It also promoted a diagonal symmetric-point statement to a
claim about the full ordered spectral polynomial.

A387 retitles the result as `Catalan--Dynkin parity for a
Virasoro-normalized symmetric branch`, with `\ClaimStatusConditional`
and licensing \(\alpha+\gamma+\varepsilon\) via the chosen generator
branch and `\hypAmbientWtCpl+\effKoszul`.  The theorem now fixes a
chosen one-generator branch \((\mathcal A,a)\), assumes
\[
\varphi_2^a(x)=x+2,\qquad
\varphi_3^a(x)=(x+2)(x+3),
\]
and assumes the symmetric-point rightmost-reduction identity.  Its
conclusion is exactly
\[
\varphi_{2q}^a(x)=0,\qquad
\varphi_{2n+3}^a(x)=(-1)^nC_n\prod_{m=2}^{2n+3}(x+m).
\]
The theorem now says explicitly that this is not a theorem about every
one-generator class-\(\mathbf M\) algebra and that reflection parity is
only a diagonal-quotient statement, not a statement about the full
ordered spectral polynomial.

A387 adds `catalan_dynkin_field_polynomial()` and
`catalan_dynkin_parity_scope()` to `compute/lib/gravity_3d_engine.py`,
with `TestCatalanDynkinParityScope` in
`compute/tests/test_gravity_3d_engine.py`.  It also repairs
`compute/tests/test_session_results.py`, which still carried the
retired W3 branch `gamma=61440/(c^2(5c+22)^3)`; the test now uses the
A386 pole-cleared coefficient
\(N_{2r}^W=S_{2r}^Wc^{2r-3}(5c+22)^{3(r-1)}\).

A387 verification: Python compilation passed for the touched compute
files.  The focused Catalan--Dynkin guard passed with `4 passed`;
`compute/tests/test_session_results.py` passed with `40 passed`;
`compute/tests/test_catalan_factorisation.py` passed with `54 passed`;
the full gravity engine surface passed with `152 passed, 1 skipped`.
Stale fixed-string sweeps found no live arbitrary-class-\(\mathbf M\)
Catalan--Dynkin claim and no live retired W3 `gamma=61440` branch in
the checked sources.  `make verify-licensing` exits 0 with zero
blockers and 79 warnings elsewhere; the next gravity warning is
`thm:crossing-stasheff`.  `make fast` converged after two passes with
2505 pages, zero undefined citations, zero undefined references, and
zero rerun requests.

A388 repairs `thm:crossing-stasheff`, the next gravity licensing
warning after A387.  The theorem had identified
\(d_{\bar B}^2=0\) with the four-point conformal bootstrap crossing
equation.  That was a stage collapse: the bar identity proves the
chain-level Stasheff/Borcherds compatibility of singular chiral OPE
channels after choosing a collision chart and OPE-channel comparison;
full conformal bootstrap crossing also requires analytic
conformal-block convergence, spectrum data, and positivity/unitarity.
The proof also incorrectly introduced a cohomological \(m_3^H\)
correction to four-point crossing, although minimal arity three gives
strict associativity of transferred \(m_2\).

A388 retitles the theorem as a conditional chiral OPE-channel
crossing statement, with `\ClaimStatusConditional` and licensing
\(\alpha+\beta+\gamma\) via an ordered \(\FM_3(\C)\) collision chart,
chiral OPE-channel comparison, and `\hypAmbientWtCpl`.  The theorem now
defines \(\mathcal A_3^\bullet(i,j,k)\) and states
\[
\mathcal A_3^\bullet(i,j,k)+\mathcal C_3^\bullet(i,j,k)=0
\]
as the singular-OPE channel identity after pairing the output with a
fourth state.  The proof identifies the two binary bar faces with the
singular \(s\)- and \(t\)-channel residues and identifies
\(d\,m_3+m_3\,d\) as the chain-level ternary contact boundary.  The
Virasoro clause is now the chosen Stasheff gauge \(m_3=-A_3\), not a
failure of strict minimal \(m_2\)-associativity.

A388 adds `crossing_stasheff_scope()` to
`compute/lib/gravity_3d_engine.py` and `TestCrossingStasheffScope` to
`compute/tests/test_gravity_3d_engine.py`.  The guard records the
\(\alpha+\beta+\gamma\) hypotheses, rejects the full-bootstrap reading,
checks \(m_3=-A_3\) on symbolic and rational Virasoro samples, and
guards the theorem source against the retired crossing slogan and full
bootstrap hierarchy title.

A388 verification: Python compilation passed for the touched compute
files.  The focused crossing/Stasheff guard passed with `3 passed`;
the full gravity engine surface passed with `155 passed, 1 skipped`.
Stale fixed-string sweeps found no live Vol II copy of the retired
theorem title, the full-bootstrap-crossing assertion, the
cohomological \(m_3^H\) crossing correction, the full-bootstrap
hierarchy remark title, or the `bar complex gives exact values` slogan.
The only exact cross-volume hits were archival Vol I resume markdown
files, not active manuscript inputs.  `make verify-licensing` exits 0
with zero blockers and 78 warnings elsewhere; the next gravity warning
is `thm:bootstrap-shapovalov`.  `make fast` converged after two passes
with 2505 pages, zero undefined citations, zero undefined references,
and zero rerun requests.

A389 repairs `thm:bootstrap-shapovalov`, the next gravity licensing
warning after A388.  The theorem had asserted a coefficient-wise
inequality
\[
|C^k_{ij,\mathrm{unitary}}|^2\le
|C^k_{ij,\mathrm{Verma}}|^2
\]
by projecting out null states.  That is not invariant without a fixed
orthonormal finite-level basis and normalization: individual OPE
coordinates can change under basis change.  The correct
Shapovalov/bootstrap positivity statement is a positive-semidefinite
Gram-norm contraction for the whole finite-level channel vector.

A389 retitles the theorem as `Shapovalov Gram positivity and channel
contraction`, with `\ClaimStatusConditional` and licensing
\(\alpha+\beta+\gamma+\delta\) via finite-level Shapovalov
normalisation, `\hypAmbientWtCpl`, and a unitary Hilbert quotient.  The
statement now fixes \(M_{c,h}[N]\), the Shapovalov Gram matrix
\(G_N(c,h)\), and a \(\beta\)-normalised channel vector
\(v_{ij}^{(N)}\).  The pointwise coefficient inequality is replaced by
\[
\|P_Nv_{ij}^{(N)}\|_{G_N}^2
\le
\|v_{ij}^{(N)}\|_{G_N}^2,
\]
where \(P_N\) is the orthogonal projection to
\((\operatorname{Rad}G_N)^\perp/\operatorname{Rad}G_N\).

A389 also corrects the minimal-model wording: unitary Virasoro minimal
models are the \(c=1-6/(m(m+1))\), \(m\ge2\), series, not arbitrary
\(c_{p,q}\).  The proof now separates inverse Shapovalov propagators
and Kac determinant denominators in raw HPL summands from quotient
positivity and from null-vector Ward/BPZ analytic-bootstrap data.  The
summary table now says `Channel positivity` and refers to
`eq:bootstrap-channel-contraction`.

A389 adds `shapovalov_channel_norm_squared()`,
`shapovalov_projected_channel_norm_squared()`, and
`shapovalov_bootstrap_scope()` to `compute/lib/gravity_3d_engine.py`,
with `TestShapovalovBootstrapScope` in
`compute/tests/test_gravity_3d_engine.py`.  The guard checks a diagonal
Gram contraction example, records the \(\alpha+\beta+\gamma+\delta\)
hypotheses, rejects a coordinate-wise OPE bound, and guards the theorem
source against the retired coefficient inequality and
generic-minimal-model unitarity wording.

A389 verification: Python compilation passed for the touched compute
files.  The focused Shapovalov guard passed with `3 passed`; the full
gravity engine surface passed with `158 passed, 1 skipped`.  Stale
fixed-string sweeps found no live Vol II, Vol I, or Vol III copy of the
retired theorem title, coefficient-wise inequality, or `coefficients
are bounded above by the Verma-module values` sentence.  The old
`eq:bootstrap-bound` label has no live source references.  `make
verify-licensing` exits 0 with zero blockers and 77 warnings elsewhere;
the next gravity warning is `prop:large-c-bootstrap`.  `make fast`
converged after two passes with 2505 pages, zero undefined citations,
zero undefined references, and zero rerun requests.

A390 repairs `prop:large-c-bootstrap`, the next gravity licensing
warning after A389.  The proposition had asserted that the shadow
obstruction tower matches the large-\(c\) conformal bootstrap and that
the \(n\)-th bootstrap correction is controlled by \(S_{n+2}\).  That
collapsed the scalar Virasoro shadow branch with the full analytic
large-\(c\) bootstrap, which also requires non-vacuum blocks, OPE
density, single-valued crossing, positivity, and analytic convergence.
The proof also contained the wrong exponent
\([t^r]H=O(c^{3-r})\), contradicting \(S_4\sim2/c^2\).

A390 retitles the proposition as `Scalar large-c Virasoro shadow
branch`, with `\ClaimStatusConditional` and licensing
\(\gamma+\varepsilon+\delta\) via
`\hypAmbientWtCpl+\effScalarShadowProj` and an identity-block contact
comparison.  The statement fixes
\[
H(t)=t^2\sqrt{Q_{\mathrm{Vir}}(t)},\qquad
Q_{\mathrm{Vir}}(t)=(c+6t)^2+\frac{80t^2}{5c+22},
\qquad S_r=[t^r]H(t)/r.
\]
It now states the explicit asymptotic
\[
S_r=\frac{8(-6)^{r-4}}{r}c^{2-r}+O(c^{1-r})
\qquad (r\ge4),
\]
derived from the closed-form coefficient theorem using
\(D=80/(5c+22)=16c^{-1}+O(c^{-2})\) and \(F_r(0)=1\).  The radius
statement is the scalar radius
\(R_{\mathrm{scal}}=c\sqrt{(5c+22)/(180c+872)}\sim c/6\).

A390 adds `virasoro_large_c_shadow_asymptotics()` and
`large_c_bootstrap_scope()` to `compute/lib/gravity_3d_engine.py`, with
`TestLargeCBootstrapScope` in `compute/tests/test_gravity_3d_engine.py`.
The guard computes the leading constants
\(2,-48/5,48,-1728/7,1296,-6912\) from the scalar coefficients, records
the \(\gamma+\varepsilon+\delta\) hypotheses, rejects the full-bootstrap
reading, and guards against the retired \(O(c^{3-r})\) exponent.

A390 verification: Python compilation passed for the touched compute
files.  The focused large-\(c\) guard passed with `3 passed`; the full
gravity engine surface passed with `161 passed, 1 skipped`.  Stale
fixed-string sweeps found no live Vol II, Vol I, or Vol III copy of the
retired theorem title, full-large-\(c\)-bootstrap assertion,
\(S_{n+2}\)-controls-bootstrap sentence, or wrong exponent.  The only
exact old-exponent hits were archival Vol I healing markdown files, not
active manuscript inputs.  `make verify-licensing` exits 0 with zero
blockers and 76 warnings elsewhere; the next gravity warning is
`thm:otoc-r-matrix`.  `make fast` converged after two passes with 2507
pages, zero undefined citations, zero undefined references, and zero
rerun requests.

A391 repairs `thm:otoc-r-matrix`, the next gravity licensing warning
after A390.  The theorem had written a scalar thermal sum
\(F_{\mathrm{OTOC}}^{\mathrm{mon}}\) with
\(|C_{VWp}|^2 e^{-\beta h_p}/Z(\beta)\) and a braiding phase
\(e^{-2\pi i(h_p-h_V-h_W)}\).  That collapsed the chiral
block-local-system monodromy supplied by the \(R\)-matrix with the full
normalized thermal OTOC, treated the braiding as scalar in every
channel, and disagreed in sign/orientation with `eq:r-matrix-diagonal`.

A391 retitles the theorem as `OTOC block monodromy from the
\(R\)-matrix`, with `\ClaimStatusConditional` and licensing
\(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypVirTorusBlock`, a conformal-block local system, an
OTO continuation path, and thermal trace data.  The theorem now states
\[
F_{\mathrm{OTOC}}^{\chi,\mathrm{mon}}(t)
=
\sum_{p,q}a_p^\beta\,
\rho_{\mathrm{blk}}(\gamma_{\mathrm{OTO}})_{p}^{\,q}\,
\mathcal F_q(z_{\mathrm{OTO}}(t)).
\]
The diagonal phase \(\exp(2\pi i(h_p-h_V-h_W))\) is asserted only after
choosing a multiplicity-free channel diagonalizing the braiding; the
opposite OTO orientation gives the inverse phase.  The theorem
explicitly does not determine the thermal coefficients, the
anti-holomorphic sector, the normalized full OTOC, or the Lyapunov
exponent.

A391 also repairs the adjacent MSS proof paragraph: monodromy supplies
only the braid action on the chiral block vector; Lorentzian growth
comes from analytic continuation of the Virasoro block plus thermal
coefficients and HHLL identity-block dominance.  The summary table now
says `OTOC block continuation` rather than full OTOC.

A391 adds `otoc_braiding_phase()` and `otoc_r_matrix_scope()` to
`compute/lib/gravity_3d_engine.py`, with `TestOTOCRMatrixScope` in
`compute/tests/test_gravity_3d_engine.py`.  The guard checks the
positive and inverse orientation phases, records the
\(\beta+\gamma+\delta\) hypotheses including
`\hypAmbientWtCpl+\hypVirTorusBlock`, rejects a full normalized OTOC
statement, and guards the theorem source against the retired scalar
thermal sum.

A391 verification: Python compilation passed for the touched compute
files.  The focused OTOC/R-matrix guard passed with `3 passed`; the
full gravity engine surface passed with `164 passed, 1 skipped`.  Stale
fixed-string sweeps found no live Vol II, Vol I, or Vol III copy of the
retired theorem title, scalar \(F_{\mathrm{OTOC}}\) sum,
Boltzmann-weight formula, inverse sign, or `phase becomes` paragraph.
`make verify-licensing` exits 0 with zero blockers and 75 warnings
elsewhere; the next gravity warning is `thm:mss-bound-bar`.  `make
fast` converged after two passes with 2507 pages, zero undefined
citations, zero undefined references, and zero rerun requests.

A392 repairs `thm:mss-bound-bar`, the next gravity licensing warning
after A391.  The old theorem title `MSS bound and annular bar
curvature` lacked a hypothesis/effectivity tag and risked treating the
annular bar curvature or wrap-around differential as the proof of the
Maldacena--Shenker--Stanford bound.  Its proof also contained the
overbroad strip-function slogan that analyticity on a strip and
growth \(e^{\lambda t}\) alone force \(\lambda\le2\pi/\beta\).

A392 retitles the theorem as an MSS analytic strip bound with
annular-bar input, with `\ClaimStatusConditional` and licensing
\(\beta+\gamma+\delta\) via
`\hypAmbientWtCpl+\hypVirTorusBlock+\hypModularCardy` and the
normalized thermal OTOC strip hypotheses.  The statement now separates
the analytic theorem from the algebraic datum: the Virasoro annular bar
complex supplies the genus-\(1\) curvature and wrap-around periodicity
datum, but it does not supply boundedness on the half-strip,
positivity, thermal trace coefficients, or the Lyapunov exponent, and
it does not identify the wrap-around mode with \(\lambda_L\).

A392 rewrites the proof so the inequality is explicitly the MSS strip
theorem applied to the normalized physical thermal OTOC.  The HHLL
identity block now contributes, in the pre-scrambling range,
\[
\mathcal F_{\mathrm{id}}^{\mathrm{HHLL}}(t)
=1+\frac{48h_W^2}{c}e^{2\pi t/\beta}
 +O(c^{-2}e^{4\pi t/\beta}),
\]
so saturation \(\lambda_L=2\pi/\beta\) is asserted only under
identity-block dominance and the thermal normalization package.  The
KMS paragraph now says \(d_{\mathrm{wrap}}\) records annular
periodicity after the trace/block normalization is supplied; it is not
the source of MSS boundedness or a Lyapunov eigenmode.

A392 also repairs the local chaos dictionary.  The row now reads
`MSS Lyapunov bound \(\lambda_L\le2\pi/\beta\)` sourced by the MSS
strip theorem plus annular periodicity datum, and the closing sentence
now says the bar complex records algebraic inputs rather than
containing the complete quantum chaotic content.  The downstream
scrambling-time proof no longer cites a non-existent `Argument~1` of
the theorem.

A392 adds `mss_bound_value()` and `mss_annular_bar_scope()` to
`compute/lib/gravity_3d_engine.py`, with `TestMSSAnnularBarScope` in
`compute/tests/test_gravity_3d_engine.py`.  The guard checks the
\(2\pi/\beta\) value, records the \(\beta+\gamma+\delta\) hypotheses,
rejects the annular-curvature-alone reading, and forbids the retired
theorem title, generic strip-function slogan, wrap-around-mode claim,
and old chaos-dictionary row.

A392 verification: Python compilation passed for the touched compute
files.  The focused MSS guard passed with `3 passed`; the full gravity
engine surface passed with `167 passed, 1 skipped`.  Stale fixed-string
sweeps found the retired formulations only as negative guard strings in
`compute/tests/test_gravity_3d_engine.py`, not as live Vol II, Vol I,
or Vol III claims.  `git diff --check` passed on the touched files.
`make verify-licensing` exits 0 with zero blockers and 74 warnings
elsewhere; the next gravity warning is `prop:scrambling-time`.  `make
fast` converged after two passes with 2507 pages, zero undefined
citations, zero undefined references, and zero rerun requests.

A393 repairs `prop:scrambling-time`, the next gravity licensing warning
after A392.  The old proposition stated
\(t_*=(\beta/2\pi)\log c\) as an exact large-\(c\) scrambling time for
the Virasoro system and treated the scalar shadow tower as if it were
the full thermal OTOC expansion.  Its interpretation remark also
claimed equal Fulton--MacPherson stratum weights at \(t_*\), interior
dominance for \(t>t_*\), and full scrambling from the ordered bar
complex alone.

A393 retitles the proposition as a conditional HHLL identity-block
scrambling scale, with `\ClaimStatusConditional` and licensing
\(\gamma+\varepsilon+\delta\) via
`\hypAmbientWtCpl+\effScalarShadowProj+\hypVirTorusBlock+\hypModularCardy`
and nonzero \(O(1)\) probe normalization.  The statement now introduces
\[
u(t)=A_{\mathrm{id}}c^{-1}e^{2\pi t/\beta},\qquad
A_{\mathrm{id}}\in\mathbb R_{>0},\quad A_{\mathrm{id}}=O(1),
\]
and defines \(t_*\) as the first scale where \(u(t)\) is order \(1\).
The result is
\[
t_*=(\beta/2\pi)\log(c/A_{\mathrm{id}})+O_\beta(1)
    =(\beta/2\pi)\log c+O_\beta(1),
\]
so the exact equality without the additive normalization constant is
not asserted.

A393 rewrites the proof so the physical expansion parameter comes from
the HHLL identity-block/thermal normalization package.  The scalar
large-\(c\) formula
\[
S_r=8(-6)^{r-4}c^{2-r}/r+O(c^{1-r}),\qquad r\ge4,
\]
is retained only as the scalar identity/contact branch after
`\effScalarShadowProj`; it is not identified with the full thermal
OTOC series.  The bar-complex reading is now the loss of uniform
fixed-degree approximation on that selected branch at \(u(t)\asymp1\),
not a proof of physical scrambling.

A393 updates the chaos dictionary row to
`Scrambling scale \(t_*=(\beta/2\pi)\log c+O_\beta(1)\)` sourced by the
HHLL identity block plus scalar shadow truncation scale.  The
bar-complex interpretation now explicitly says that thermal
Hilbert-space and OTOC hypotheses are required for the physical
scrambling interpretation.

A393 adds `scrambling_time_from_amplitude()` and
`scrambling_time_scope()` to `compute/lib/gravity_3d_engine.py`, with
`TestScramblingTimeScope` in `compute/tests/test_gravity_3d_engine.py`.
The guard checks the threshold equation, records the
\(\gamma+\varepsilon+\delta\) hypotheses, requires the
\(A_{\mathrm{id}}\) dependence, and forbids the retired exact-log
proposition, full-OTOC-from-shadow claim, equal-strata claim, and old
chaos-dictionary row.

A393 verification: Python compilation passed for the touched compute
files.  The focused scrambling-time guard passed with `3 passed`; the
full gravity engine surface passed with `170 passed, 1 skipped`.
Stale fixed-string sweeps found the retired formulations only as
negative guard strings in `compute/tests/test_gravity_3d_engine.py`,
not as live Vol II, Vol I, or Vol III claims.  `git diff --check`
passed on the touched files.  `make verify-licensing` exits 0 with zero
blockers and 73 warnings elsewhere; the next gravity warning is
`thm:ds-ordered-bar-intertwine`.  `make fast` converged after two
passes with 2507 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A394 repairs `thm:ds-ordered-bar-intertwine`, the next gravity
licensing warning after A393.  The old theorem stated that principal DS
reduction commutes with the ordered bar construction as a bare
dg-coalgebra equivalence and that \(R\)-descent commutes canonically.
Its proof only supplied fibrewise BRST compatibility, a PBW
spectral-sequence comparison, and informal \(R\)-descent functoriality.
The local block also had duplicate proof/remark/center artifacts.

A394 retitles the theorem as a filtered principal DS comparison for the
ordered bar, with `\ClaimStatusConditional` and licensing
\(\alpha+\beta+\gamma+\delta\) via the principal DS chart,
`\hypAmbientWtCpl`, `\hypKZSDR`, finite-weight BRST concentration, and
HPL-transferred \(R\)-descent.  The theorem now assumes an ordered
BRST-bar bicomplex, complete finite-weight PBW/ghost convergence,
ghost-degree-zero concentration, and supplied transferred \(R\)-descent.
The conclusion is a filtered quasi-isomorphism
\[
B^{\mathrm{ord}}(\mathcal W_k(\mathfrak g))
\xrightarrow{\simeq_{\mathrm{filt}}}
H^0_{\mathrm{DS}}(B^{\mathrm{ord}}(V_k(\mathfrak g)))
\]
in the homotopy category of ordered bar dg-coalgebras.

A394 rewrites the proof so \(d_{\mathrm{BRST}}\) acts over the ordered
configuration cover, anticommutes with
\(d_{\mathrm{bar}}^{\mathrm{ord}}\) by the chiral-derivation property
of the BRST current, and uses strong convergence plus ghost
concentration to justify the \(H^0_{\mathrm{DS}}\) target.  Without
concentration the right-hand side would be
\(H^\bullet_{\mathrm{DS}}\).  The \(R\)-descent paragraph now uses the
HPL-transferred \(R\)-matrix datum rather than saying DS leaves the
descent map unchanged.

A394 also scopes the gravitational Yangian remark: the ordered-bar
theorem gives the candidate line-side model, while strict
dg-shifted-Yangian promotion still requires the honest-gaps hypotheses.
The ordered trichotomy table is now an ordered-bar shadow on the
filtered principal surface, not a primitive classification before
\(R\)-descent is fixed.

A394 adds `ds_ordered_bar_intertwining_scope()` to
`compute/lib/gravity_3d_engine.py`, and extends
`TestVirasoroHPLTransfer` in `compute/tests/test_gravity_3d_engine.py`.
The guard requires the filtered comparison, `\hypAmbientWtCpl`,
`\hypKZSDR`, BRST concentration, HPL-transferred \(R\)-descent, and
negative assertions against unconditional ordered functoriality,
unchanged \(R\)-descent, strict Yangian promotion, and primitive
ordered classification.

A394 verification: Python compilation passed for the touched compute
files.  The focused HPL/DS guard class passed with `18 passed`; the full
gravity engine surface passed with `172 passed, 1 skipped`.  Stale
fixed-string sweeps found retired formulations only in negative guard
strings, except for archival Vol I healing markdown not included in
active manuscript inputs.  `git diff --check` passed on the touched
files.  `make verify-licensing` exits 0 with zero blockers and 72
warnings elsewhere; the next gravity warning is
`prop:gravity-page-curve`.  `make fast` converged after two passes with
2509 pages, zero undefined citations, zero undefined references, and
zero rerun requests.

A395 repairs `prop:gravity-page-curve`, the next gravity licensing
warning after A394.  The old proposition stated a `Gravitational Page
curve` from gravitational Koszul duality without a claim-status tag,
without the real comparison window \(0<c<26\), and without the
evaporation interval on which the decreasing branch remains
non-negative.  The following Borel conjecture also contradicted the
linear Page-time formula by asserting
\(t_P\sim\hbar^{-1}e^{S_{\mathrm{BH}}}\).

A395 retitles the result as a conditional Page profile from the
two-sector raw transseries, with `\ClaimStatusConditional` and
licensing \(\beta+\gamma+\delta+\varepsilon\) via
`\hypAmbientWtCpl+\effScalarShadowProj+\hypModularCardy`, the raw Page
transseries, and the real same-family window \(0<c<26\).  It defines
\[
t_{\mathrm{evap}}=6S_{\mathrm{BH}}/(26-c)
\]
and states
\[
S_{\mathrm{rad}}^{\mathrm{model}}(t)=
\min\left(\frac c6t,\,
S_{\mathrm{BH}}-\frac{26-c}{6}t\right)
\]
only on \(0\le t\le t_{\mathrm{evap}}\).  The formula is explicitly
not asserted outside the real window, beyond evaporation, or from
scalar complementarity alone.

A395 rewrites the evidence and island remarks so the decreasing branch
uses the same-family line--Verdier comparison sector
\(B(\cA_{\mathrm{line}})\), not a bare physical saddle \(B(\cA^!)\).
The Page time is the branch equality
\[
\frac c6t_P=S_{\mathrm{BH}}-\frac{26-c}{6}t_P,\qquad
t_P=3S_{\mathrm{BH}}/13.
\]
The Borel conjecture now states the Stokes wall as the equal-real-part
condition
\[
\Re\left((I_{\mathrm I}(t)-I_{\mathrm H}(t))/\hbar
-\log\mathfrak s(t)\right)=0,
\]
and says the Borel singularity is not a second formula for the Page
time.

A395 adds `page_curve_profile()` and `page_curve_scope()` to
`compute/lib/gravity_3d_engine.py`, with `TestPageCurveProfile` in
`compute/tests/test_gravity_3d_engine.py`.  The guard checks the branch
crossing, rejects numeric values outside \(0<c<26\) and
\(S_{\mathrm{BH}}>0\), records the conditional hypotheses, and forbids
the retired gravitational-Koszul title, \(B(\cA^!)\)-dominance wording,
old Borel-scale label, and exponential Page-time formula.

A395 verification: Python compilation passed for the touched compute
files.  The focused Page/Maloney guard set passed with `7 passed`; the
full gravity engine surface passed with `176 passed, 1 skipped`.  Stale
fixed-string sweeps found retired strings absent from the active Vol II
target and compute helper, except as negative guard strings in
`compute/tests/test_gravity_3d_engine.py`.  Cross-volume sweeps found
older Vol I Page-curve advertisements in standalone and compute layers;
these are quarantined as a separate Vol I rectification obligation
because the current pass repaired the Vol II target and the Vol I
standalone file is already dirty.  `git diff --check` passed on the
touched files.  `make verify-licensing` exits 0 with zero blockers and
71 warnings elsewhere; the next gravity warning is
`prop:gravity-desitter`.  `make fast` converged after two passes with
2509 pages, zero undefined citations, zero undefined references, and
zero rerun requests by the build counter; direct fresh-log scans found
no individual undefined citation/reference warnings or undefined-control
failures, while the generic `There were undefined references` line and
existing pdfTeX named-destination warnings remain.

A396 repairs `prop:gravity-desitter`, the next gravity licensing warning
after A395.  The old proposition stated a de Sitter shadow obstruction
tower under the bare Wick rotation \(\ell\mapsto i\ell\), used
\(\pi\ell_{\mathrm{dS}}/(2G)\) instead of the repository's
Brown--Henneaux \(G_N\) normalization, identified
\(c_{\mathrm{dS}}=13\) with a Nariai limit without constructing a
Nariai geometry, and included the Banks finite-Hilbert-space conjecture
as a proposition item.

A396 retitles the proposition as a de Sitter scalar shadow and
horizon-normalized entropy statement, with `\ClaimStatusConditional` and
licensing \(\alpha+\gamma+\varepsilon\) via the de Sitter real-section
metric normalization and `\hypAmbientWtCpl+\effScalarShadowProj`.  It
now fixes
\[
\ell_{\mathrm{dS}}>0,\qquad G_N>0,\qquad
c_{\mathrm{dS}}=3\ell_{\mathrm{dS}}/(2G_N)>0
\]
before making any scalar claim.

A396 replaces the entropy formula by the three-dimensional
Gibbons--Hawking chain
\[
S_{\mathrm{dS}}
=A_{\mathrm{hor}}/(4G_N)
=2\pi\ell_{\mathrm{dS}}/(4G_N)
=\pi\ell_{\mathrm{dS}}/(2G_N)
=\pi c_{\mathrm{dS}}/3.
\]
The scalar tower is now explicitly
\[
F_g^{\mathrm{dS},\mathrm{sc}}(\mathrm{Vir}_{c_{\mathrm{dS}}})
=\kappaChHodge(\mathrm{Vir}_{c_{\mathrm{dS}}})
\lambda_g^{\mathrm{FP}}.
\]
The \(c_{\mathrm{dS}}=13\) line is only the same-family fixed point,
with \(\kappa_{\mathrm{dS}}=13/2\) and
\(S_{\mathrm{dS}}=13\pi/3\); the proposition explicitly asserts no
Nariai geometry, no de Sitter Hilbert space, and no dS/CFT correlator
construction.

A396 moves the Banks-type formula to a separate
`\ClaimStatusHeuristic` remark:
\[
\dim(\mathcal H_{\mathrm{dS}})\stackrel{\mathrm{heur}}{=}
\exp(S_{\mathrm{dS}})=\exp(\pi c_{\mathrm{dS}}/3).
\]
The remark says this is compatible with the horizon entropy
normalization but is not a theorem of the Virasoro scalar shadow tower.

A396 adds `desitter_central_charge()`,
`desitter_horizon_entropy_from_radius()`,
`desitter_horizon_entropy_from_c()`, and `desitter_shadow_profile()` to
`compute/lib/gravity_3d_engine.py`, with `TestDesitterShadowScope` in
`compute/tests/test_gravity_3d_engine.py`.  The guard checks the entropy
formula three ways, the \(c=13\) fixed point, positive real-section
inputs, and negative assertions against old Wick-rotation, Nariai,
Banks, and bare-\(G\) text.  A396 also updates
`compute/tests/test_deletion_ledger.py` so Brown--Henneaux normalization
guards expect \(G_N\).

A396 verification: Python compilation passed for the touched compute
files.  The focused de Sitter gravity guard passed with `5 passed`; the
focused deletion-ledger guard passed with `1 passed`; the full gravity
engine surface passed with `180 passed, 1 skipped`; and the full
deletion-ledger surface passed with `62 passed`.  Active Vol II stale
sweeps found retired de Sitter strings only in negative guard
assertions.  Cross-volume sweeps found older Vol I de Sitter copies in
`entanglement_modular_koszul.tex`, `standalone/three_dimensional_quantum_gravity.tex`,
and `vol2_physical_consequences_transfer.tex`; these are quarantined as
a separate Vol I propagation obligation because the first two Vol I
files are already dirty.  `git diff --check` passed on the touched
files.  `make verify-licensing` exits 0 with zero blockers and 70
warnings elsewhere; the next gravity warning is `prop:gravity-jt`.
`make fast` converged after two passes with 2509 pages, zero undefined
citations, zero undefined references, and zero rerun requests by the
build counter; direct fresh-log scans found no individual undefined
citation/reference warnings or undefined-control failures, while the
generic `There were undefined references` line and existing pdfTeX
named-destination warnings remain.

A397 repairs `prop:gravity-jt`, the next gravity licensing warning
after A396.  The old proposition said that the \(c\to\infty\)
Schwarzian limit makes the shadow metric degenerate to
\(y=\sin(2\pi\sqrt{x})/(4\pi)\), that topological recursion on this
curve reproduces Weil--Petersson volumes, and that the Bernoulli decay
of the shadow partition function analytically completes the divergent
\((2g)!\) JT perturbative series.  The downstream class-\(\mathcal S\)
remark repeated the same JT curve/completion claim and added a
Schur-series UV-completion assertion.

A397 retitles the statement as a conditional Schwarzian/JT scalar
limit with `\ClaimStatusConditional` and licensing
\(\beta+\gamma+\delta+\varepsilon\) via the Schwarzian comparison
datum, `\hypAmbientWtCpl+\effScalarShadowProj`, and the JT
spectral-curve normalization.  The theorem now fixes
\[
\mathfrak s_{\mathrm{JT}}
=(\text{large-}c\text{ scaling},\,z,\,x=z^2,\,
y_{\mathrm{WP}},\,K_{\mathrm{EO}},\,\mathcal C_E,\,
\text{matrix/Stokes completion datum})
\]
before making a JT comparison.

A397 replaces the branch-ambiguous curve formula by
\[
x=z^2,\qquad
y_{\mathrm{WP}}(z)=\frac{\sin(2\pi z)}{4\pi},\qquad
\rho_0(E)=\frac{1}{i\pi}y_{\mathrm{WP}}(i\sqrt E)
=\frac{\sinh(2\pi\sqrt E)}{4\pi^2}.
\]
The WP/topological-recursion claim is now explicitly external
Mirzakhani--Eynard--Orantin--Saad--Shenker--Stanford comparison data,
not a theorem of the Vol II scalar shadow tower.  The Bernoulli scalar
series is stated to be convergent but not an analytic completion of the
JT \((2g)!\)-series; a matrix-integral or Stokes datum, for example
`\hypStokes` on the raw gravitational sector, is required.

A397 adds `jt_wp_spectral_curve_y()`, `jt_disk_density()`,
`jt_wp_to_density_balance()`, and `jt_schwarzian_scope()` to the Vol II
gravity engine, with `TestJTSchwarzianScope` guarding the Taylor
normalization, the exact \(z=i\sqrt E\) density conversion, nonnegative
energy scope, the conditional hypothesis package, and absence of the
retired JT correspondence/completion wording.

A397 also propagates the same repair into live Vol I surfaces:
`entanglement_modular_koszul.tex` now states a conditional
Schwarzian/JT comparison with contour data, `higher_genus_modular_koszul.tex`
uses \(x=z^2,\ y_{\mathrm{WP}}\) for the classical sine curve, and the
Vol I JT compute docstrings no longer say the scalar shadow tower
produces the double-scaled matrix model.

A397 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused Vol II JT/de Sitter guard passed
with `9 passed`; the full Vol II gravity engine surface passed with
`185 passed, 1 skipped`; `make verify-licensing` reports zero blockers
and 69 warnings, with the next gravity warning
`thm:3dg-borcherds-ensemble-bridge`; and scoped Vol II
`git diff --check` passed.  Vol I Python compilation passed for the
touched JT compute files; the focused Vol I entanglement scalar-typing
guard passed with `4 passed`; the broader Vol I JT compute surface
passed with `214 passed, 2 skipped`; and scoped Vol I
`git diff --check` passed.  Cross-volume fixed-string scans for the
old affirmative degeneration, topological-recursion reproduction, and
scalar-to-JT completion slogans are clean except for negative test
needles.

A398 repairs `thm:3dg-borcherds-ensemble-bridge`, the next gravity
licensing warning after A397.  The old scalar theorem already kept the
K3\(\times E\) Borcherds input separate from Maloney--Witten gravity,
but the theorem head lacked explicit claim status and licensing data.
The statement also needed executable guards separating the BPS scalar
convention from the Oberdieck--Pixton reduced-DT normalization, the
one-variable G\"ottsche coefficient \(p_{24}(5)=176256\), and the
Bruinier reduced obstruction class.

A398 retitles the theorem as a `\ClaimStatusProvedElsewhere` scalar
shadow with licensing \(\beta+\gamma+\varepsilon\) via the
Gritsenko--Nikulin half-K3 normalisation,
`\hypAmbientWtCpl+\effScalarShadowProj`, and the
DMVV--G\"ottsche--Bruinier scalar comparison.  The obstruction sentence
now names the half-K3 Gritsenko--Nikulin scalar convention before
stating \(c_2^\triangle=0\) and \(c_3^{\mathrm{Br}}=-64[H_3]\), keeping
it separate from the reduced \(\phi_{-2,1}\) convention whose
coefficient is \(-8\).

A398 strengthens the scalar-residual remark: the scalar theorem proves
neither the ordered Virasoro bar trace identity nor a Maloney--Witten
equality, and it does not prove the conditional operator square
`thm:k3-borcherds-operator-square`.

A398 adds `colored_partition_number()` and
`k3_borcherds_scalar_bridge_scope()` to the Vol II gravity engine, with
`TestK3BorcherdsScalarBridge` guarding the \(24\)-coloured partition
values through \(176256\), the OP/BPS scalar prefactors, the half-K3
\(-64[H_3]\) Bruinier convention, the separate reduced \(-8\) convention
flag, and the negative scalar-to-operator promotion flags.

A398 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused Borcherds scalar bridge guard
passed with `4 passed, 186 deselected`; the full Vol II gravity engine
surface passed with `189 passed, 1 skipped`; scoped Vol II
`git diff --check` passed; `make verify-licensing` reports zero
blockers and 68 warnings, with the next gravity warning
`thm:k3-borcherds-operator-square`; and `make fast` converged after two
passes with 2511 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A399 repairs `thm:k3-borcherds-operator-square`, the next gravity
licensing warning after A398.  The operator theorem now states a
conditional implication from an explicit P1 datum and finite Hall gates,
not an unconditional construction of the protected-Pfaffian operator
from the scalar reciprocal.

A399 retitles the theorem with `\ClaimStatusConditional` on the P1 datum
and finite Hall gates, with licensing
\(\alpha+\beta+\gamma+\varepsilon\) via the oriented Hall chart, the
Hall--Borcherds comparison, `\hypAmbientWtCpl`, and
`\effPfaffOrient+\effPBWnoExtra`.  The statement assumes
\[
\mathfrak p_1=(b_{\mathrm{Hall}},
\Theta_{\mathrm{Hall}\to\mathrm{Borch}},
\mathcal C^{\mathrm{ch,cyc}}_{X,\Lambda},
\operatorPrim{X},\mathfrak D_X,o_{\mathrm{Pf}},\iota_{\mathrm{aut}})
\]
and says explicitly that
\(\mathrm{Tr}_{\mathrm{cyc}}(\mathfrak D_X)=\operatorPrim{X}\) is a
cyclic-trace condition, not a scalar-character statement.

A399 makes the finite Hall gate system part of the theorem:
\(\operatorname{rad}_{\mathrm{Hall},N}/\operatorname{rad}_{N}=0\),
\(D^{\mathrm{fin}}_{\mathrm{Hall}}\) exists,
\(\operatorname{Borch}\circ\operatorname{Hall}\) is height-compatible,
and the class-\(\mathcal S\) Schur comparison is
\(\SCchtop\)-trace compatible.  Without those gates the diagram is a
shadow comparison, not an object equivalence.

A399 also upgrades the P1 executable profile: `k3_borcherds_operator_profile()`
and `k3_borcherds_hall_chiral_square()` now record conditional status,
the P1 datum, finite gates, orientation, PBW/no-extra-root
effectiveness, scope residuals, and negative flags against unconditional
operator construction.  The manuscript guard rejects the retired
existence wording.

A399 verification: Python compilation passed for the touched P1 engine
and guard.  The focused P1/operator-square guard passed with `27
passed`; the full P1 protected-Pfaffian guard passed with `27 passed`;
the focused Borcherds scalar bridge guard still passed with `4 passed,
186 deselected`; the full Vol II gravity engine surface passed with
`189 passed, 1 skipped`; scoped Vol II `git diff --check` passed; stale
existence/equivalence wording appears only in negative test needles;
`make verify-licensing` reports zero blockers and 67 warnings, with the
next gravity warning `prop:3dg-k3-scalar-no-promotion`; and `make fast`
converged after two passes with 2511 pages, zero undefined citations,
zero undefined references, and zero rerun requests.

A400 repairs `prop:3dg-k3-scalar-no-promotion`, the next gravity
licensing warning after A399.  The proposition is now a proved
scalar-character non-faithfulness theorem, not an unlicensed caveat.

A400 retitles the proposition with `\ClaimStatusProvedHere` and
licensing \(\beta+\gamma+\varepsilon\) via scalar non-faithfulness,
`\hypAmbientWtCpl`, and `\effScalarShadowProj`.  The statement begins in
the completed \(\SCchtop\)-ambient and under the scalar-character
projection.  It says that
\(\Phi_{10}^{\mathrm{un}}=\Delta_5^2\),
\(Z_{\mathrm{BPS}}^{K3\times E}=(\Phi_{10}^{\mathrm{un}})^{-1}\), and
\(\kappa_{\mathrm{BKM}}=5\) give only the Borcherds coordinate of the
K3\(\times E\) tuple \(\{0,3,5,24\}\).

A400 strengthens the theorem: these scalar data do not determine a
filtered \(\SCchtop\)-morphism to the gravitational boundary line, nor a
chain-level identification with the ordered Virasoro bar trace.  A
promotion requires four \(\beta\)-comparison data: positive-half
Hall--Borcherds bialgebra morphism, Drinfeld-double/current-envelope
extension, filtered gravity-line morphism, and derived-centre trace
compatibility.

A400 rewrites the proof through `lem:scalar-non-faithfulness`.  The
scalar character functor forgets the differential, Hall product,
Drinfeld pairing, current-envelope locality, mixed open--closed
operations, and derived-centre trace; adjoining a filtered acyclic
summand preserves scalar character while changing the chain object.

A400 adds `k3_scalar_no_promotion_scope()` to the Vol II gravity engine,
guarding the scalar tuple, the acyclic-summand witness, the four missing
comparison data, and negative flags against determining the filtered
\(\SCchtop\) morphism, ordered Virasoro bar trace, Maloney--Witten sum,
or object equivalence.

A400 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused K3 Borcherds scalar/no-promotion
guard passed with `6 passed, 186 deselected`; the full Vol II gravity
engine surface passed with `191 passed, 1 skipped`; stale positive
promotion wording appears only in negative test needles; scoped Vol II
`git diff --check` passed; `make verify-licensing` reports zero blockers
and 66 warnings, with the next gravity warning
`prop:3dg-heptagon-growth-bound`; and `make fast` converged after two
passes with zero undefined citations, zero undefined references, and
zero rerun requests.  `pdfinfo out/main.pdf` reports 2511 pages.

A401 repairs `prop:3dg-heptagon-growth-bound`, the next gravity
licensing warning after A400.  The old statement bounded ordered-bar
coefficients by \(C(W,M)^n n!\), depending only on generation weight and
OPE norm.  That misses the finite generator rank: a family can have
arbitrarily many generators of weight \(\le W\) with the same OPE norm.

A401 retitles the proposition with `\ClaimStatusProvedHere` and
licensing \(\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effKoszul` and a finite strong-generator profile.  It
now fixes
\[
G=\bigoplus_{1\le j\le W}G_j,\qquad R=\dim G<\infty,
\]
a Lyndon--PBW basis, and a finite-type norm with OPE constants bounded
by \(M\).

A401 changes the bound to
\[
|b_n(A)|\le C(W,R,M)^n n!,\qquad
C(W,R,M)=4R\max(1,M)e^{\pi\sqrt{2/3}}.
\]
The proof now displays the separate factors: Catalan tree shapes,
generator labels \(R^n\), partition decompositions, OPE norm, and PBW
symmetrisation.  It also states that \(R\)-dependence is essential.

A401 propagates the corrected constant to the factorisation-BV Gevrey
restatement and the master-bridge assembly remark.  The later
factorisation-BV proposition remains a separate licensing warning; this
pass only corrected the dependent constant and finite-type hypothesis.

A401 adds `heptagon_growth_bound_constant()` and
`heptagon_growth_bound_scope()` to the Vol II gravity engine, with
`TestHeptagonGrowthBound` guarding the \(R\)-factor, source text,
downstream constant propagation, local Borel radius, and negative flags
against sectorial continuation, Borel singularity location, and
Maloney--Witten interpretation.

A401 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused heptagon growth guard passed with
`4 passed, 192 deselected`; the full Vol II gravity engine surface
passed with `195 passed, 1 skipped`; stale `C(W,M)` / `C(W, M)` and old
two-parameter wording are absent from the searched surfaces; scoped Vol
II `git diff --check` passed; `make verify-licensing` reports zero
blockers and 65 warnings, with the next gravity warning
`prop:3dg-zwegers-heisenberg-shadow`; and `make fast` converged after
two passes with 2511 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A402 repairs `prop:3dg-zwegers-heisenberg-shadow`, the next gravity
licensing warning after A401.  The old statement asserted vanishing of
the Zwegers shadow for \(\mathcal H_k\) at arbitrary real level while
using the rank-one Fock character \(\eta^{-1}\).  It also placed the
determinant-line power \(\eta^{-k}\) next to the Fock trace, which
blurred oscillator rank, Heisenberg level, and scalar determinant
character.

A402 retitles the proposition with `\ClaimStatusProvedHere` and
licensing \(\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effScalarShadowProj` and the inverse eta multiplier.
It now assumes nonzero level \(k\in\mathbb R^\times\), the standard
rank-one Fock vacuum trace, and proves
\[
\mathcal F^{\mathrm{raw}}_{\mathcal H_k}(\tau)=\eta(\tau)^{-1}.
\]
The statement explicitly does not assert modularity for unbranched real
powers \(\eta(\tau)^{-r}\).

A402 rewrites the proof through the rescaling invariance of the
rank-one Heisenberg Fock spectrum, partition multiplicities \(p(n)\),
Dedekind's inverse-eta multiplier, and the ordinary Borel transform
being entire of exponential type zero.  The downstream off-Koszul
class-\(\mathsf G\) bridge and master assembly now use the rank-one Fock
trace and inverse eta multiplier; the modular PVA covariance paragraph
now attaches the weight \(-k/2\) of \(\eta^{-k}\) only to the chosen
determinant-character branch/multiplier lane.

A402 adds `heisenberg_zwegers_shadow_scope()` and
`TestHeisenbergZwegersShadow`, with negative guards against unbranched
real eta-power modularity, level-dependent oscillator multiplicities,
finite Borel singularities, mock completion, and Maloney--Witten
path-integral promotion.

A402 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused Heisenberg guard passed with
`3 passed, 196 deselected`; the full Vol II gravity engine surface
passed with `198 passed, 1 skipped`; scoped stale-phrase scans found
the retired proposition language only in negative test needles; scoped
Vol II `git diff --check` passed; `make verify-licensing` reports zero
blockers and 64 warnings, with the next gravity warning
`prop:3dg-hardy-ramanujan-cardy`; and `make fast` converged after two
passes with 2511 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A403 repairs `prop:3dg-hardy-ramanujan-cardy`, the next gravity
licensing warning after A402.  The old proposition displayed the correct
generic Virasoro vacuum coefficient formula \(p(n)-p(n-1)\), but it had
no licensing/status surface and no proof environment, and nearby prose
fused the PBW/Hardy--Ramanujan oscillator law with the physical Cardy
density.

A403 retitles the proposition with `\ClaimStatusProvedHere` for the
PBW--Hardy--Ramanujan formula and `\ClaimStatusProvedElsewhere` for
Cardy, with licensing \(\gamma+\delta\) via
`\hypAmbientWtCpl+\hypModularCardy`.  The statement now works in the
completed scalar character lane, assumes generic \(c>1\) with only the
translation null \(L_{-1}\mathbf1=0\), and states explicitly that the
Cardy theorem is not a Verma-module coefficient identity.

A403 adds the proof: the PBW basis is indexed by partitions with parts
\(\ge2\), so
\[
\sum_{n\ge0}\rho_{\mathrm{vac}}(n)q^n
=\prod_{m\ge2}(1-q^m)^{-1}
=(1-q)\prod_{m\ge1}(1-q^m)^{-1},
\]
hence \(\rho_{\mathrm{vac}}(n)=p(n)-p(n-1)\).  The
Hardy--Ramanujan asymptotic for \(p(n)\), after first difference, gives
the prefactor \(\pi/(12\sqrt2)\).

A403 propagates the split to the downstream Zwegers theorem and
conditionality remark: the Eichler derivative recovers the physical
Cardy density only under `\hypModularCardy`; the
PBW--Hardy--Ramanujan computation supplies only the bare oscillator
law.  A Vol I check found that the mirror surface already separates
universal Virasoro vacuum growth from physical Cardy growth, so no Vol I
edit was required.

A403 adds `partition_number()`,
`virasoro_vacuum_verma_coefficient()`,
`virasoro_vacuum_verma_asymptotic()`, and
`virasoro_hardy_ramanujan_cardy_scope()` to the gravity engine, plus
`TestVirasoroHardyRamanujanCardy` guarding low-degree coefficients, the
asymptotic prefactor, status/licensing, and negative flags against
treating Cardy as a Verma coefficient identity.

A403 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused Hardy--Ramanujan/Cardy guard
passed with `5 passed, 199 deselected`; the full Vol II gravity engine
surface passed with `203 passed, 1 skipped`; stale fused
Hardy--Ramanujan--Cardy wording appears only in negative test needles;
scoped Vol II `git diff --check` passed; `make verify-licensing`
reports zero blockers and 63 warnings, with the next gravity warning
`prop:3dg-fact-bv-l0-gevrey`; and `make fast` converged after two
passes with 2511 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A404 repairs `prop:3dg-fact-bv-l0-gevrey`, the next gravity licensing
warning after A403.  The old proposition claimed that the
factorisation-BV tensor-product coefficients obeyed the same constant
\(C(W,R,M)=4R\max(1,M)e^{\pi\sqrt{2/3}}\) as the ordered-bar heptagon
bound.  That was not justified: the factorisation-BV side has extra
BV generator labels, a regularised propagator/vertex norm, and a
weight-splitting convolution.

A404 retitles the proposition with `\ClaimStatusProvedHere` and
licensing \(\gamma+\varepsilon\) via
`\hypAmbientWtCpl+\effKoszul` plus a finite factorisation-BV graph
profile.  The statement now fixes the ordered-bar profile \((W,R)\),
the BV graph profile \((U,D)\), OPE norm \(M\), and regularised
propagator/vertex norm \(P\).

A404 replaces the constant by
\[
C_{\mathrm{fact}}(W,U,R,D,M,P)
=8RD\max(1,M,P)e^{\pi\sqrt{2/3}}.
\]
The proof now counts weight splittings \(n+1\le2^n\), heptagon
tree/PBW data on the ordered-bar factor, finite BV graph labels \(D^i\),
regularised kernel factors \(P^i\), and \(i!\,j!\le n!\).  It explicitly
does not construct sectorial Borel continuation, a Borel-plane saddle,
or a Maloney--Witten interpretation.

A404 propagates the distinction to the class-\(\mathsf M\) Borel
summability proof: the heptagon estimate and the finite-profile
factorisation-BV estimate are parallel local-convergence estimates under
their own finite-type norms, not equivalent statements with the same
constant.

A404 adds `factorization_l0_gevrey_constant()` and
`factorization_l0_gevrey_scope()` to the gravity engine, plus
`TestFactorizationL0Gevrey` guarding the \(D\)- and \(P\)-dependence,
local-only scope, source text, and negative flags against the retired
same-constant/equivalent-presentation wording.

A404 verification: Python compilation passed for the touched Vol II
gravity engine and guard.  The focused factorisation/heptagon guard
passed with `7 passed, 200 deselected`; the full Vol II gravity engine
surface passed with `206 passed, 1 skipped`; stale same-constant and
equivalence scans found no live manuscript copy; `make verify-licensing`
reports zero blockers and 62 warnings, with the next gravity warning
`prop:3dg-fact-bv-bcov-instanton`; and `make fast` converged after two
passes with 2513 pages, zero undefined citations, zero undefined
references, and zero rerun requests.

A405 repairs `prop:3dg-fact-bv-bcov-instanton`, the next gravity
licensing warning after A404.  The old proposition conflated the
renormalised factorisation-BV action, the Borel coordinate location, and
the Virasoro beta/residue coefficient by writing
\[
S^{\mathrm{eff}}[\Phi_\ast]=2\pi\,c(c-25)/24
\quad\Rightarrow\quad
\zeta_\ast=1/(2\pi).
\]
This was not a legitimate instanton-action statement: after the Borel
coordinate is chosen, the singularity location is the action difference,
whereas \(c(c-25)/24\) is the candidate one-loop residue.

A405 retitles the proposition with `\ClaimStatusConditional` and
licensing \(\alpha+\gamma+\delta\) via Borel-coordinate normalisation,
`\hypAmbientWtCpl+\hypStokes`.  The corrected datum is
\[
A_{\mathrm{BV}}(\Phi_\ast)
=S^{\mathrm{eff}}[\Phi_\ast]-S^{\mathrm{eff}}[\Phi_0]
=\frac{1}{2\pi},
\qquad
\rho_c=\frac{c(c-25)}{24}.
\]
The conditional Stokes datum is therefore
\((\zeta_\ast,\operatorname{Res}_{\zeta_\ast})
=(1/(2\pi),c(c-25)/24)\).  The text now says explicitly that the local
QME and the finite-profile Gevrey bound do not prove the existence of
\(\Phi_\ast\), sectorial continuation, saddle-to-ordered-bar
comparison, or a Maloney--Witten interpretation.

A405 propagates the distinction to the alien-derivative closure theorem:
it now invokes `prop:3dg-fact-bv-bcov-instanton` together with
`\hypStokes` before using the singularity at \(\zeta_\ast=1/(2\pi)\).

A405 adds `factorization_bcov_candidate_location()`,
`factorization_bcov_candidate_residue()`, and
`factorization_bcov_instanton_scope()` to the gravity engine, plus
`TestFactorizationBcovInstanton` guarding special residue values,
status/licensing, the action/residue separation, and negative phrases
against the retired action formula.

A405 verification: Python compilation passed for the touched gravity
engine and guard.  The focused BCOV guard passed with `4 passed, 207
deselected`; the full gravity engine surface passed with `210 passed, 1
skipped`; stale action/residue scans found old text only as negative
test needles and no Vol I / Vol III copy; `make verify-licensing`
reports zero blockers and 61 warnings, with
`prop:3dg-fact-bv-bcov-instanton` no longer listed.  Scoped
`git diff --check` passed.  `make fast` converged after two passes with
2513 pages, zero undefined citations, zero undefined references, and
zero rerun requests.  Direct final-log scans found no fatal TeX errors,
undefined controls, or unresolved citation/reference warnings; the tail
still contains existing pdfTeX destination warnings.

A406 repairs the active PVA descent warning
`prop:PVA-from-symplectic` in `pva-descent-repaired.tex`.  The
proposition identified the PVA bracket on \(H^\bullet(\A,Q)\) with the
Steinberg Poisson bracket, while its proof only gave that conclusion
after the comparison
\[
\chi\colon H^\bullet(\barB(\A))\to H^\bullet(\A,Q)
\]
is an isomorphism.  Without that isomorphism, the shifted-symplectic
Steinberg bracket lives on the bar/Koszul shadow; the direct
Fulton--MacPherson proof supplies the boundary PVA bracket.

A406 retitles the proposition with `\ClaimStatusConditional` and
licensing \(\alpha+\beta+\gamma+\varepsilon\) via boundary chart \(b\),
the comparison map \(\chi\), the derived shifted-symplectic ambient,
and `\effKoszul`.  The statement now assumes explicitly that the
bar-to-boundary comparison map is an isomorphism on the
Koszul-effective comparison locus, and the proof records what remains
true without that effectiveness.

A406 adds `TestPvaDescentSourceGuards` to
`compute/tests/test_pva_descent_chain_level.py`, requiring the
conditional status, `\effKoszul`, the comparison map \(\chi\), and the
absence of the old unconditional title.

A406 verification: Python compilation passed for the touched PVA test.
The focused source guard passed with `1 passed, 59 deselected`; the full
PVA chain-level test surface passed with `60 passed`; stale
unconditional Steinberg wording appears only as a negative test needle;
no Vol I / Vol III mirror was found; `make verify-licensing` reports
zero blockers and 60 warnings, with `prop:PVA-from-symplectic` no longer
listed.  Scoped `git diff --check` passed.  `make fast` converged after
two passes with 2513 pages, zero undefined citations, zero undefined
references, and zero rerun requests.  Direct final-log scans found no
fatal TeX errors, undefined controls, or unresolved citation/reference
warnings; the tail still contains existing pdfTeX destination warnings.

A407 repairs the degree-two PVA descent identity in
`prop:m1_m2`.  The active repaired file had
\[
Q(m_2(a,b))=m_2(Qa,b)+(-1)^{|a|}m_2(a,Qb),
\]
but the defining desuspended Stasheff convention in `axioms.tex`, the
BV-construction proof, and the FM-calculus proof give
\[
Q(m_2(a,b))+m_2(Qa,b)+(-1)^{|a|}m_2(a,Qb)=0.
\]
Thus the moved-form equality has two minus signs.  The old positive-RHS
formula was a sign error, not merely an unlicensed theorem title.

A407 retitles the active and legacy degree-two propositions with
`\ClaimStatusProvedHere` and licensing \(\gamma\) via
`\hypAmbientPro`.  It replaces the formula by the sum-zero identity,
adds the moved-form equality, and corrects the downstream descent proof:
regular and singular coefficient identities now carry the minus signs,
and representative changes are \(Q\)-exact by
\[
m_2(Qx,b)=-Qm_2(x,b)-(-1)^{|x|}m_2(x,Qb).
\]

A407 extends `TestPvaDescentSourceGuards` so the PVA source guard
requires the sum-zero formula and rejects the old positive-RHS formula.
It also updates `ainfty_sign_m1_m2()` in
`compute/lib/cross_volume_deep_bridge.py` to advertise the same
convention.

A407 verification: Python compilation passed for the touched PVA guard
and cross-volume bridge helper.  Focused PVA source guards passed with
`2 passed, 59 deselected`; the focused cross-volume bridge sign guard
passed with `7 passed, 61 deselected`; the full PVA chain-level test
surface passed with `61 passed`.  Stale positive-RHS formula scans found
only the negative test needle in this repo; the single external hit is
an old Vol I resume note.  `make verify-licensing` reports zero blockers
and 58 warnings, with both active and legacy `prop:m1_m2` warnings
removed.  `make fast` converged after two passes at 2513 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
controls, or unresolved citation/reference warnings.

A408 repairs the next binary PVA descent surface.  After A407 corrected
the \(m_1m_2\) sign, the active `lem:lambda_descends` still had
`\ClaimStatusProvedHere` without the theorem-level \(\gamma\) ambient
tag, although its proof uses the logarithmic chain/pro ambient, the
regular/singular projection of \(m_2\), and the Borel transform from
the singular OPE kernel to the polynomial \(\lambda\)-bracket.  A408
retitles the active lemma with licensing via `\hypAmbientPro` and makes
the statement name `thm:cohomology-PVA-main` and
`def:borel-transform-pva` as the constructional inputs.

A408 also repairs the legacy `pva-descent.tex` shadow: the old lemma
spoke about descent of an undifferentiated \(m_2\) class and still
referred to removed labels.  The legacy statement now says that the
regular product \(\mu(a,b)\) and singular polynomial bracket
\(\{a_\lambda b\}\) are \(Q\)-closed and representative-independent,
with the same \(\hypAmbientPro\) licensing.  Removed-label references
to `prop:m1_m2` and `lem:lambda_descends` were replaced by local prose.

A408 verification: Python compilation passed for the PVA guard file.
Focused PVA source guards passed with `4 passed, 59 deselected`; the
full PVA chain-level test surface passed with `63 passed`.  Fixed-string
scans found no surviving stale `prop:gravity-koszul-dual-branch`, no
legacy `Lemma \ref{lem:lambda_descends}`, and no legacy
`Proposition \ref{prop:m1_m2}` reference.  `make verify-licensing`
reports zero blockers and 56 warnings.  `make fast` converged after two
passes at 2513 pages, with zero undefined citations, zero undefined
references, and zero rerun requests; the fresh final log scan found no
named unresolved reference/citation warnings and no generic
undefined-reference summary.  Existing pdfTeX destination warnings
remain outside this pass.

A409 repairs PVA1 sesquilinearity.  The active `prop:PVA1_proof`
had proved status without a theorem-level \(\gamma\) ambient tag, and
the proof compressed the polynomial \(\lambda\)-bracket calculation
into “chain-level sesquilinearity descends coefficientwise.”  A409
retitles it with licensing via `\hypAmbientPro` and replaces the proof
by the explicit divided-power mode computation:
\[
(\partial a)_{(n)}b=-n\,a_{(n-1)}b,\qquad
a_{(n)}(\partial b)=\partial(a_{(n)}b)+n\,a_{(n-1)}b,
\]
with \(a_{(-1)}b=0\), followed by the Borel transform.  The second
calculation now visibly contains the \(+\lambda\)-term in
\(\{a_\lambda\partial b\}\).

A409 also repairs the legacy `pva-descent.tex` PVA1 shadows.  It adds
the same \(\hypAmbientPro\) licensing to the legacy PVA1 lemma and
proposition, removes stale references to removed
`eq:PVA1a/eq:PVA1b` and `prop:PVA1_proof` labels, replaces the old
\(\{a_\lambda b\}=[m_2(a,b)]\) compression by the singular polynomial
bracket, and fixes a malformed align row break hidden behind a
`% label removed` comment.

A409 verification: Python compilation passed for the PVA guard file.
Focused PVA source guards passed with `6 passed, 59 deselected`; the
full PVA chain-level test surface passed with `65 passed`.
Fixed-string scans over the active and legacy PVA descent files found
no surviving `eq:PVA1a`, `eq:PVA1b`,
`Definition \ref{def:ainfty_chiral}`, or
`\{a_\lambda b\} = [m_2(a,b)]` strings.  `make verify-licensing`
reports zero blockers and 53 warnings.  `make fast` converged after two
passes at 2513 pages, with zero undefined citations, zero undefined
references, and zero rerun requests; the fresh final log scan found no
named unresolved reference/citation warnings and no generic
undefined-reference summary.  Existing pdfTeX destination warnings
remain outside this pass.

A410 repairs the exchange-monodromy surface.  The active
`prop:exchange-lagrangian-monodromy` had proved status without a
theorem-level license and said the holomorphic component “winds once”
around the collision circle.  The displayed path gives instead
\[
z_1(s)-z_2(s)=2r e^{\pi i s},
\]
so the exchange is a half-turn from angle \(0\) to angle \(\pi\).  The
old proof also attributed the skew-symmetry sign to the logarithm branch
shift, but \(d\log(-\zeta)=d\log\zeta\); the sign is the oriented
Stokes boundary sign in the exchange-cylinder identity.

A410 retitles the active proposition as “Exchange half-monodromy as
Lagrangian monodromy” with licensing \(\alpha+\gamma\) via boundary
chart \(b\) and `\hypAmbientPro`.  The proof now says the square of the
exchange projects to the full generator of
\(\pi_1(\FM_2(\C))\cong\Z\), while no braid-group \(\Z\) is assigned to
the full ordered configuration space of two points in \(\C\times\R\).
The legacy `pva-descent.tex` skew-symmetry shadow was given the same
licensing and half-turn proof.  Nearby context in
`pva-expanded-repaired.tex` and `axioms.tex` now names the configuration
space of two points in the three-dimensional bulk rather than calling
the configuration space itself three-dimensional.

A410 verification: Python compilation passed for the PVA guard file.
Focused PVA source guards passed with `9 passed`; the full PVA
chain-level test surface passed with `68 passed`.  Retired-phrase scans
found no live copy of the old exchange title, no `carries the winding
class`, no logarithm-branch sign sentence, and no full-configuration
braid-\(\Z\) claim.  `make verify-licensing` reports zero blockers and
51 warnings.  `make fast` converged after two passes at 2513 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the fresh final log scan found no unresolved reference/citation
warnings, no undefined control sequence, and no fatal error.  Existing
pdfTeX destination warnings remain outside this pass.

A411 repairs the PVA2 exchange homotopy and its sign shadows.  The
active `lem:PVA2_proof` still lacked theorem-level licensing and still
described \(\tau^\ast\) as \(\zeta\mapsto-\zeta\), “which becomes
\(\lambda\mapsto-\lambda\) under Borel transform,” followed by a
separate translation transfer.  That intermediate operation is not a
PVA operation: after the divided-power Borel transform and translation
covariance, the vertex-algebra exchange is the single substitution
\(\lambda\mapsto-\lambda-\partial\).

A411 retitles `lem:PVA2_proof` with licensing \(\alpha+\gamma\) via the
oriented exchange half-monodromy and `\hypAmbientPro`.  Its proof now
computes \(\zeta(s)=2r e^{\pi i s}\), applies Cartan--Stokes to the
oriented interval family, records \(d\log(-\zeta)=d\log\zeta\), and
states explicitly that there is no separate PVA operation
\(\lambda\mapsto-\lambda\).  The same exchange license was added to
`prop:product-commutative` and `prop:PVA2_proof`.

A411 also removes stale shifted-PVA2 wording from
`pva-expanded-repaired.tex`, `appendices/pva-expanded.tex`, and
`pva-preview.tex`.  The visible PVA2 convention is now consistently
\[
\{a_\lambda b\}=-(-1)^{|a||b|}\{b_{-\lambda-\partial}a\},
\]
while shifted signs remain confined to nested Jacobi/Leibniz residue
orientation.

A411 verification: Python compilation passed for the PVA guard file.
Focused PVA source guards passed with `12 passed`; the full PVA
chain-level surface passed with `71 passed`; the independent PVA
coefficient engine passed with `11 passed`.  `make verify-licensing`
reports zero blockers and 48 warnings.  `make fast` converged after two
passes at 2513 pages, with zero undefined citations, zero undefined
references, and zero rerun requests; the final log scan found no
unresolved reference/citation warnings, no undefined control sequence,
and no fatal error.  Existing pdfTeX destination warnings remain outside
this pass.

A412 repairs product associativity in the PVA descent.  The active
`prop:product-associative` had proved status without theorem-level
\(\gamma\) licensing, and its proof said that the
\((-1)^{|a|}\) Stasheff sign becomes trivial on the symmetrized
product.  That mechanism is false: graded commutativity on cohomology
does not erase the shifted sign in the chain-level \(A_\infty\)
identity.  The correct mechanism is the desuspension already present in
the main theorem: \(\mu(a,b)=\mu_2^{\mathrm{reg}}(a,b;0)\), where
\(m_2^{\mathrm{reg}}\) is the suspended bar component.

A412 retitles `prop:product-associative` with licensing \(\gamma\) via
`\hypAmbientPro` and the unsuspended regular-product convention.  The
proof now applies \(D_A^2=0\) to \(Q\)-closed inputs, projects to the
regular--regular constant term, transports the binary vertices from
\(m_2^{\mathrm{reg}}\) to \(\mu_2^{\mathrm{reg}}\), and obtains the
ordinary associator
\[
\mu(\mu(a,b),c)-\mu(a,\mu(b,c))
\]
as an element of \(\operatorname{im}Q\), with
\(K_3^{\mathrm{reg}}(a,b,c)\) supplying the ternary homotopy.  No
vanishing of \(m_3\) is assumed.

A412 also repairs the nearby unsymmetrized/symmetrized-product remark,
the raviolo PVA proof, and the `pva-preview.tex` associativity shadow.
Those copies now say that the raw Stasheff sign is absorbed by the bar
desuspension convention, not discarded by graded commutativity.

A412 verification: Python compilation passed for the PVA guard file.
Focused PVA source guards passed with `14 passed`; the full PVA
chain-level surface passed with `73 passed`; `test_pva_axioms.py`
passed with `36 passed`.  Retired-phrase scans over chapters and
appendices found no surviving sign-trivialisation or symmetrization
mechanism.  `make verify-licensing` reports zero blockers and 47
warnings.  Scoped `git diff --check` passed.  `make fast` converged
after two passes at 2515 pages, with zero undefined citations, zero
undefined references, and zero rerun requests; the final log scan found
no fatal TeX errors, undefined control sequences, or unresolved
reference/citation warnings.  Existing pdfTeX destination warnings
remain outside this pass.

A413 repairs the Jacobi three-face proof.  The active
`lem:PVA3_proof` and `thm:Jacobi` had proved status without
theorem-level \(\gamma\) licensing, and the proof described the three
raw residues as if they were already the printed Jacobi identity.  That
was too imprecise: the printed identity fixes the signed Jacobiator
\[
J_{\lambda,\mu}(a,b,c)
=
\{a_\lambda\{b_\mu c\}\}
-
(-1)^{(|a|+1)(|b|+1)}
\{b_\mu\{a_\lambda c\}\}
-
\{\{a_\lambda b\}_{\lambda+\mu}c\}.
\]
The proof must show \(J_{\lambda,\mu}\in\operatorname{im}Q\).

A413 retitles the singular three-face lemma and Jacobi theorem with
\(\gamma\) licensing via `\hypAmbientPro`, the oriented \(\FM_3\)
incidence convention, the divided-power Borel transform, and
Arnold--Orlik--Solomon corner cancellation.  The lemma now identifies
the oriented contributions
\[
\mathfrak J_{23},\quad \mathfrak J_{13},\quad \mathfrak J_{12}
\]
as respectively
\[
\{a_\lambda\{b_\mu c\}\},\quad
-(-1)^{(|a|+1)(|b|+1)}\{b_\mu\{a_\lambda c\}\},\quad
-\{\{a_\lambda b\}_{\lambda+\mu}c\}.
\]
The theorem proof introduces \(K_3^{\mathrm{Jac}}\), the doubly
singular component of the interior \(m_3\)-homotopy, and proves the
signed Jacobiator is \(Q\)-exact.

A413 also propagates the same normalization to `fm-proofs.tex`,
`pva-expanded-repaired.tex`, `appendices/pva-expanded.tex`,
`pva-preview.tex`, and `raviolo.tex`.  These shadows now state that the
non-consecutive \(D_{\{1,3\}}\) channel is present for the dressed
doubly singular PVA form, and its sign is the incidence orientation
together with the shifted transposition past \(b\).

A413 verification: Python compilation passed for the touched guard
files.  Focused source guards plus the Arnold/Jacobi guard passed with
`19 passed`; the full PVA chain-level surface passed with `75 passed`;
`test_pva_axioms.py` plus the Arnold/Jacobi guard passed with
`39 passed`.  Retired-phrase scans over chapters and appendices found
no surviving unsigned-residue Jacobi phrasing.  `make verify-licensing`
reports zero blockers and 45 warnings.  Scoped `git diff --check`
passed.  `make fast` converged after two passes at 2515 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
control sequences, or unresolved reference/citation warnings.  Existing
pdfTeX destination warnings remain outside this pass.

A414 repairs the Steinberg/Lagrangian-convolution reading of Jacobi.
After A407 made the Steinberg Poisson comparison conditional on the
bar-to-boundary comparison map \(\chi\) and \(\effKoszul\),
`prop:jacobi-lagrangian-convolution` still claimed unconditionally that
the boundary PVA Jacobi identity is associativity of Lagrangian triple
convolution in \(\Steinb_b\).  That collapsed the Steinberg bar/Koszul
shadow into boundary cohomology.

A414 rewrites the preceding remark as a conditional Steinberg reading
and retitles the proposition with `\ClaimStatusConditional` and
licensing \(\alpha+\beta+\gamma+\varepsilon\) via boundary chart \(b\),
comparison map \(\chi\), the derived shifted-symplectic ambient, and
\(\effKoszul\).  The statement now assumes the hypotheses of
`prop:PVA-from-symplectic` and says the signed Jacobiator
\(J_{\lambda,\mu}\) is the \(\chi\)-image of the associator of
Lagrangian triple convolution.  Without the isomorphism hypothesis on
\(\chi\), it is only a Steinberg bar/Koszul-shadow interpretation, not
a proof of the boundary PVA identity.

A414 verification: Python compilation passed for the PVA guard file.
Focused source guards passed with `17 passed`; the full PVA chain-level
surface passed with `76 passed`; `test_pva_axioms.py` plus the
Arnold/Jacobi guard passed with `39 passed`.  `make verify-licensing`
reports zero blockers and 44 warnings.  Scoped `git diff --check`
passed.  `make fast` converged after two passes at 2515 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
control sequences, or unresolved reference/citation warnings.  Existing
pdfTeX destination warnings remain outside this pass.

A415 repairs the mixed regular--singular Leibniz surface.  The old
appendix still printed a quantum/Wick integral correction in the
Leibniz rule and explained it as an \(m_3\)-boundary contribution.  The
classical PVA descent does not have this term: the mixed sector has one
singular binary channel and one regular constant term, while the double
singular bracket belongs to the Jacobi sector.

A415 normalizes Leibniz by the signed defect
\[
\mathcal L_\lambda(a;b,c)
=
\{a_\lambda(bc)\}
-\{a_\lambda b\}c
-(-1)^{(|a|+1)|b|}b\{a_\lambda c\}.
\]
The oriented faces \(\mathfrak L_{23}\), \(\mathfrak L_{12}\), and
\(\mathfrak L_{13}\) give respectively the left term, the negative
first right-hand term, and the negative signed second right-hand term.
The total collision contributes the \(Q\)-exact mixed homotopy
\(K_3^{\mathrm{Leib}}\), so
\(\mathcal L_\lambda(a;b,c)\in\operatorname{im}Q\) and the defect
vanishes on cohomology.

A415 propagates the same normalization to
`pva-expanded-repaired.tex`, `appendices/pva-expanded.tex`,
`pva-preview.tex`, and `raviolo.tex`, and removes the stale raviolo
`thm:Leibniz-PVA` reference.  Cross-volume propagation also corrected
the live Vol I bridge `chapters/theory/en_koszul_duality.tex` and the
two non-backup Vol I staging copies that repeated the false integral
formula.

A415 verification: Python compilation passed for the PVA guard file.
Focused source guards passed with `19 passed`; the full PVA chain-level
surface passed with `78 passed`; `test_pva_axioms.py` plus the
Arnold/Jacobi guard passed with `39 passed`.  `make verify-licensing`
reports zero blockers and 42 warnings.  Scoped `git diff --check`
passed in Vol II and in the touched Vol I files.  `make fast`
converged after two passes at 2515 pages, with zero undefined
citations, zero undefined references, and zero rerun requests; the
final log scan found no fatal TeX errors, undefined control sequences,
or unresolved reference/citation warnings.  Existing pdfTeX destination
warnings remain outside this pass.

A416 repairs the vacuum/unit PVA axiom.  The active theorem had no
licensing package and displayed only \(\{\mathbf 1_\lambda a\}=0\);
the old appendix derived a spurious odd-degree obstruction from a
symmetrized product calculation before resolving it informally; the
executable D6 oracle checked only the left vacuum bracket.

A416 retitles `prop:PVA4_proof` with \(\gamma\) licensing via the
strict unital \(A_\infty\)-chiral structure in `\hypAmbientPro`, the
unsuspended regular product, regular/singular projection, and the
divided-power Borel transform.  The statement now proves the full
two-sided vacuum package:
\[
[\mathbf 1]\cdot[a]=[a]\cdot[\mathbf 1]=[a],\qquad
\{[\mathbf 1]_\lambda[a]\}=0,\qquad
\{[a]_\lambda[\mathbf 1]\}=0,\qquad
\partial[\mathbf 1]=0.
\]
The proof uses the desuspended strict unit identities and the
vanishing of both singular projections with a unit insertion.

A416 also licenses `cor:PVA-structure` with \(\alpha+\gamma\),
propagates the same two-sided unit convention to `pva-descent.tex`,
`appendices/pva-expanded.tex`, and `raviolo.tex`, and strengthens
`verify_d6_unit` so it checks both brackets, both product units, and
\(\partial1=0\).

A416 verification: Python compilation passed for the PVA library and
guard file.  Focused source guards plus D6 checks passed with
`26 passed`; the full PVA chain-level surface passed with `80 passed`;
`test_pva_axioms.py` plus the Arnold/Jacobi guard passed with
`39 passed`.  `make verify-licensing` reports zero blockers and
39 warnings.  Scoped `git diff --check` passed.  `make fast`
converged after two passes at 2515 pages, with zero undefined
citations, zero undefined references, and zero rerun requests; the
final log scan found no fatal TeX errors, undefined control sequences,
or unresolved reference/citation warnings.  Existing pdfTeX destination
warnings remain outside this pass.

A417 repairs the higher-operation/topological-nullhomotopy surface.
The active PVA chapter already denied full \(A_\infty\) formality, but
`lem:topological-contraction` and `prop:m3_vanish` still hid the
distinction between open ordered contractibility, a relative
compactified topological filling, and \(Q\)-exactness of direct
higher operations.  The old appendix also used "vanishing of higher
operations" inside the Jacobi proof, and the raviolo Landau-Ginzburg
example said \(m_3\) vanishes by the usual \(A_\infty\) argument.

A417 retitles the contraction lemma as an open ordered topological
contraction with \(\gamma\) licensing through the ordered topological
chain ambient in `\hypAmbientPro` and the translation quotient.  The
lemma now states only the open ordered factor; it explicitly does not
construct a relative bounding chain in the compactified operadic pair.

A417 retitles `prop:m3_vanish` as a conditional \(Q\)-exactness
criterion with \(\gamma+\varepsilon\) licensing through
`\hypAmbientPro`, the relative compactified ordered topological
factor, compatible topological nullhomotopies, and the factorisation
\[
\widetilde\omega_k
=
\omega_k^{\mathrm{hol}}\wedge\omega_k^{\mathrm{top}}.
\]
Under that datum \(m_k\) represents zero on \(H^\bullet(\A,Q)\);
without the supplied relative bounding chains no higher-operation
vanishing statement is asserted.

A417 propagates the same distinction to `pva-descent.tex`,
`appendices/pva-expanded.tex`, and `raviolo.tex`.  The legacy split now
uses the positive orthant \(\R_{>0}^{k-1}\), the relative boundary
\(\partial_{\mathrm{rel}}\Gamma_{k-1}\), and treats higher
\(Q\)-exactness as a separate conditional criterion.  The appendix
Jacobi proof now uses only the local arity-\(3\) total-collision
boundary \(m_1(m_3(a,b,c))=Qm_3(a,b,c)\), not global vanishing of
higher operations.  The raviolo example now permits an \(m_3\) term
unless the compatible topological nullhomotopy datum is supplied.

A417 verification: Python compilation passed for the PVA guard file.
Focused source guards passed with `23 passed`; the full PVA chain-level
surface passed with `82 passed`; `test_pva_axioms.py` plus the
Arnold/Jacobi guard passed with `39 passed`.  Retired-phrase scans
found the old higher-vanishing, usual-\(A_\infty\), and open-simplex
language only as negative guard needles in tests; active Vol I/III
chapter/staging scans found no live copy.  `make verify-licensing`
reports zero blockers and 36 warnings.  Scoped `git diff --check`
passed.  `make fast` converged after two passes at 2515 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
control sequences, unresolved citation/reference warnings, or rerun
requests.  Existing pdfTeX destination warnings remain outside this
pass.

A418 repairs PVA functoriality and the K3 recognized pro-limit.  The
active functor theorem had no licensing tag and did not state the
morphism class used by the proof: unit preservation was invoked "by
definition" although the definition did not include
\(f(\mathbf 1_A)=\mathbf 1_B\), and regular/singular projection plus
Borel coefficient extraction were used without naming continuity in
the logarithmic pro ambient.  The K3 windowwise proposition likewise
formed the completed pro-limit PVA without naming strict transitions
or the weight-completed ambient.

A418 replaces the morphism definition by strict unital
\(\hypAmbientPro\)-continuous morphisms.  Such a morphism is a cochain
map commuting with \(\partial\), preserving the unit, preserving all
logarithmic kernels coefficientwise, and commuting with completion,
Laurent coefficient extraction, regular/singular projection, and the
polynomial Borel transform.

A418 retitles `thm:pva-descent-functor` as strict functoriality with
\(\gamma\) licensing through \(\hypAmbientPro\) and strict unital
coefficientwise morphisms.  Product preservation is now proved from
regular projection and constant-term extraction; bracket preservation
from coefficientwise singular extraction and the Borel transform; unit
preservation from the explicit strict-unital clause.

A418 also retitles `prop:pvadr-windowwise-descent` as windowwise PVA
descent and recognized pro-limit with \(\alpha+\beta+\gamma\)
licensing through the height-\(N\) Hall--Borcherds recognition datum,
the maps \(\rho_N\), strict \(\hypAmbientPro\) morphisms, and
\(\hypAmbientWtCpl\).  The completed PVA
\(\widehat{\mathcal D}_{\Delta_5}\) is now asserted only after
compatible strict transitions carrying radicals into radicals and
commuting with recognition maps; the proof descends these transitions
by strict PVA functoriality before taking the relationwise inverse
limit in \(\hypAmbientWtCpl\).

A418 verification: Python compilation passed for the PVA guard file.
Focused source guards passed with `25 passed`; the full PVA chain-level
surface passed with `84 passed`; `test_pva_axioms.py` plus the
Arnold/Jacobi guard passed with `39 passed`.  Retired-phrase scans
found no active Vol II, Vol I, or Vol III chapter/appendix copy of the
old functor title, old window title, unqualified transition-map phrase,
or unsupported unit-preservation sentence.  `make verify-licensing`
reports zero blockers and 34 warnings.  Scoped `git diff --check`
passed.  `make fast` converged after two passes at 2517 pages, with
zero undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
control sequences, unresolved citation/reference warnings, or rerun
requests.  Existing pdfTeX destination warnings remain outside this
pass.

A419 repairs the legacy `pva-descent.tex` split warnings.  The file is
not the active input surface, but it is still scanned by the Beilinson
gate and still functions as a shadow copy for future audits.  Its
Jacobi, higher-operation, and Leibniz blocks were not merely missing
macros: Jacobi did not expose the signed Jacobiator/AOS mechanism in
the theorem title, higher-operation exactness did not state the
relative compactified topological filling package, and Leibniz still
printed the old unsigned rule.

A419 gives the legacy Jacobi block \(\gamma\) licensing through
\(\hypAmbientPro\), oriented \(\FM_3\) incidence, and
Arnold--Orlik--Solomon corner cancellation.  It now states
\[
J_{\lambda,\mu}
=
\{a_\lambda\{b_\mu c\}\}
-
(-1)^{(|a|+1)(|b|+1)}\{b_\mu\{a_\lambda c\}\}
-
\{\{a_\lambda b\}_{\lambda+\mu}c\}
\in\operatorname{im}Q,
\]
and the nearby three-face table uses the same
\(\mathfrak J_{23},\mathfrak J_{13},\mathfrak J_{12}\) signs as the
active theorem.

A419 gives the legacy higher-operation block
\(\gamma+\varepsilon\) licensing through \(\hypAmbientPro\), the
relative compactified ordered topological factor, compatible
topological nullhomotopies, and the factorisation
\(\widetilde\omega_k=\omega_k^{\mathrm{hol}}\wedge
\omega_k^{\mathrm{top}}\).  It now says ordinary contractibility of
the open ordered configuration space is not enough and that no
higher-operation vanishing statement is asserted without supplied
relative bounding chains.

A419 gives the legacy Leibniz block \(\gamma\) licensing through
\(\hypAmbientPro\), oriented \(\FM_3\) incidence, the unsuspended
regular product, and the divided-power Borel transform.  It now states
the signed defect
\[
\mathcal L_\lambda(a;b,c)
=
\{a_\lambda(bc)\}
-\{a_\lambda b\}c
-(-1)^{(|a|+1)|b|}b\{a_\lambda c\}
=QK_3^{\mathrm{Leib}}(a;b,c).
\]

A419 verification: Python compilation passed for the PVA guard file.
Focused source guards passed with `26 passed`; the full PVA chain-level
surface passed with `85 passed`; `test_pva_axioms.py` plus the
Arnold/Jacobi guard passed with `39 passed`.  Retired-phrase scans
found the old legacy warning titles only as negative guard needles,
found no unsigned Leibniz live phrase, and found no remaining
`pva-descent.tex` licensing warning.  `make verify-licensing` reports
zero blockers and 31 warnings.  Scoped `git diff --check` passed.
`make fast` converged after two passes at 2517 pages, with zero
undefined citations, zero undefined references, and zero rerun
requests; the final log scan found no fatal TeX errors, undefined
control sequences, unresolved citation/reference warnings, or rerun
requests.  Existing pdfTeX destination warnings remain outside this
pass.

A420 repairs the principal \(\cW_N\) tempering/radius surface.  The four
gate-scanned statements in `wn_tempered_closure_platonic.tex` lacked
explicit licensing, and the radius text conflated finite Riccati upper
envelopes with sharp ordinary-radius equality.  A finite envelope plus
Stirling proves
\[
\limsup_r (|S_r|/r!)^{1/r}=0,
\]
but the Cauchy--Hadamard equality
\[
\rho_*^{\cW_N}(c)=|c|/\beta_N
\]
requires leading-root sharpness
\[
\limsup_r |S_r(\cW_{N,c})|^{1/r}=\beta_N/|c|.
\]
Without sharpness, the finite envelope gives only the lower bound
\(\rho_*^{\cW_N}(c)\ge |c|/\beta_N\).

A420 retitles the Stirling proposition with \(\gamma\) licensing through
the finite shadow-envelope ambient and \(\hypAmbientWtCpl\); retitles
`thm:wn-tempered-all-N` with \(\alpha+\gamma+\varepsilon\) licensing
through fixed rank, generic central charge, finite Riccati envelope,
finite propagation, and \(\hypAmbientWtCpl\); and retitles the WN
radius proposition and tempered-stratum corollary with the corresponding
sharpness and finite-envelope hypotheses.  The theorem proof now says
the finite envelope is an explicit assumption and is not derived from
Fateev--Lukyanov.

A420 propagates the same sharp-radius distinction to
`beta_N_closed_form_all_platonic.tex`, including the infinite-rank
\(\cW_\infty\) consequence: the sharp-radius scale tends to zero as
\(\beta_N\to\infty\), so the finite-rank argument supplies no uniform
\(N\)-independent Banach radius.  The compute scope map now records
`leading_root_sharpness = assumed_for_radius_equality`.

A420 guard coverage: `test_wn_tempered_closure.py` now reads the
manuscript and requires the WN claim-status/licensing surface, the
finite-envelope assumption, the Fateev--Lukyanov non-derivation,
leading-root sharpness, and the root-test lower-bound fallback.
`test_w4_beta_direct.py` now guards the current W4 obstruction language:
the spin-lane witness is not a Riccati bridge, and the full-Miura
\(A_5^{\cW_4}\) computation remains absent.

A420 verification: Python compilation passed for the touched WN/W4
libraries and tests.  The focused WN/Beta/W4 bundle passed with
`43 passed`.  `make verify-licensing` reports zero blockers and
27 warnings, removing all four `wn_tempered_closure_platonic.tex`
warnings from the previous 31-warning surface.  Scoped
`git diff --check` passed.  `make fast` converged after two passes at
2517 pages, with zero undefined citations, zero undefined references,
and zero rerun requests; the final log scan found no fatal TeX errors,
undefined control sequences, unresolved citation/reference warnings, or
rerun requests.  Existing pdfTeX destination warnings remain outside
this pass.

A421 repairs the BV/HT physics theorem surface.  The split
`bv_ht_physics.tex` file still carried seven licensing warnings, but the
mathematical defect was not only missing labels: boundary chiral
algebras, observables/bar comparison, and 6d hCS \(E_3\) structures were
being stated too close to primitive physical slogans.  The dangerous
failure mode was the old \(\C^3\)-counts-as-\(E_3\) shortcut: Dunn
additivity applies to residual locally constant directions only after
the Dolbeault--Koszul reduction, not to the unreduced
six-real-dimensional Dolbeault complex.

A421 retitles the HT-from-4d theorem with \(\alpha+\gamma\) licensing
through the Costello holomorphic-twist datum, BV gauge fixing,
dimensional reduction, and \(\hypKMHTBV\).  It retitles the boundary
theorem with \(\alpha+\beta+\gamma\) licensing and makes the boundary
chiral algebra conditional on a chosen HCS boundary condition,
bulk-to-boundary OPE comparison, and level normalisation.  It retitles
the central-charge theorem with \(\alpha+\beta\) licensing through fixed
affine level, Drinfeld--Sokolov BRST reduction, and \(\hypDSBRST\).
It retitles the observables theorem with \(\beta+\gamma\) licensing and
makes the bar-cobar identification conditional on the supplied BRST
complex and BV/bar comparison package.

A421 repairs the 6d hCS \(E_3\)/Dunn/\(P_3\) surface in the split file
and propagates the same condition to the main-input
`six_d_hcs_e3_chiral_avatar_platonic.tex`.  The active theorem now says
the \(E_3\) structure is on the \(\bar\partial\)-Dolbeault-reduced,
translation-equivariant observables, licensed by \(\gamma\) through
Dolbeault--Koszul reduction, Costello--Gwilliam factorisation, and
\(\hypKMHTBV\).  It explicitly excludes an \(E_3\) assertion on the
unreduced six-real-dimensional Dolbeault complex.

A421 guard coverage: `test_bv_ht_bridge_bv_datum.py` now reads the
split BV/HT file, the active physical-origins copy, and the live 6d hCS
avatar chapter.  It requires the theorem-head licensing, boundary/OPE
and level data, DS-BRST status, observables BV/bar condition, and the
Dolbeault-reduced \(E_3\)/Dunn/\(P_3\) language.

A421 verification: Python compilation passed for the touched BV/HT and
6d hCS guard files.  The focused BV/HT and 6d hCS anomaly bundle passed
with `12 passed`.  The propagation scan found only repaired
Dolbeault--Koszul statements or unrelated \(E_1\) notation.  `make
verify-licensing` reports zero blockers and 20 warnings, removing the
seven `bv_ht_physics.tex` warnings from the previous 27-warning
surface.  Scoped `git diff --check` passed.  `make fast` converged after
two passes at 2517 pages, with zero undefined citations, zero undefined
references, and zero rerun requests; the final log scan found no fatal
TeX errors, undefined control sequences, unresolved citation/reference
warnings, or rerun requests.  Existing pdfTeX destination warnings
remain outside this pass.

A422 treats the iterated-Sugawara and Part VI endpoint surface.  The
initial warning cluster in `e_infinity_topologization.tex` and
`programme_climax_platonic.tex` exposed a real proof-scope defect:
finite-rung theorems and endpoint summaries consumed completed
BV--BRST ambient data, DS/BV comparison, translation lifts, antighost
coherence, Faber scalar projection, and \(W_\infty\) endpoint
hypotheses without putting those data in the theorem heads.

A422 adds named hypothesis macros for the missing chart/ambient lanes
and retitles the theorem heads accordingly.  The finite ladder now
carries \(\hypAmbientWtCpl+\hypTLift\); iterated Sugawara carries
\(\hypDSBRST+\hypTLift+\hypAmbientWtCpl\); Casimir commutativity
carries \(\hypDSBRST+\hypAmbientWtCpl+\hypTLift\); closed antighost
homotopies carry \(\hypDSBRST+\hypAmbientWtCpl\); the pure operadic
spiral carries \(\hypOperadicAmbient\).  The Heisenberg, affine
Kac--Moody, and finite \(W_N\) rungs now carry their actual
HT/BV or DS/T-lift packages.  The \(N=\infty\) statements are no longer
unconditional: they name
\(\hypProchazka,\hypCKL,\hypPRSh,\hypYamada\).

A422 also repairs the scalar lane in the partition-function and F6
blocks: the scalar is \(\kappaChHodge\), not bare \(\kappa\), and the
Faber scalar residue is explicitly under
\(\hypFaberScalar+\effScalarShadowProj\).  Source-level guards now read
the topologisation, Part VI, and programme-climax files and require
these packages to remain visible.

A422 verification: focused Python tests passed with `21 passed`; `make
verify-licensing` reports zero blockers and zero warnings; scoped
`git diff --check` passed; `make fast` converged after two passes at
2517 pages with zero undefined citations, zero undefined references,
and zero rerun requests; final log scan found no fatal TeX errors,
undefined controls, unresolved citation/reference warnings, or rerun
requests.  A user-run `make release` also completed during this pass,
but the current-tree verification is the later `make fast`.

A423 repairs the κ-primitivity proof body.  The Class A theorem had
claimed that Hodge numbers vary continuously and that
\(\kappaChHodge(X_t)\) varies over a Hodge moduli factor.  This was the
wrong invariant: the programme's Vol I Hodge row is the
Serre-cancelled supertrace \(\sum_q(-1)^q h^{0,q}\), which vanishes on
the K3-fibered CY3 lane.  The theorem now proves the true statement:
on the CHL window \(N\in\{1,2,3,4,6\}\), the Hodge and fibre rows are
constant while
\(\kappaBKM(\Phi_N)=c_N(0)/2\in\{5,2,1,1,1/2\}\).  Hence the
Borcherds row is not recovered from the Hodge/fibre rows, and the
additive collapse fails for the correct reason.

A423 also repairs the Class B BCOV theorem.  The BCOV row is not
independent from the compact-support Euler row: it satisfies the exact
relation \(\kappaCat^{\mathrm c}=24\kappa_{\mathrm{BCOV}}\).  The
theorem now states row separation only after quotienting by this
relation and proves the surviving witness matrix has determinant
\(3\).  The adjacent summary no longer claims five-row independence on
the full CY3 moduli stack.  A local false citation tying the
\(\kappaChHeis=3\) row to a Pope--Romans--Shen K3-lattice statement was
removed.

A423 guard coverage: `test_part_vi_platonic_introduction.py` now rejects
the old false Hodge-deformation phrases and requires Serre cancellation,
the CHL Borcherds value set, the BCOV/Euler quotient relation, and the
determinant computation.  Verification passed with `12 passed` across
the focused κ bundle; `make verify-licensing` reports zero blockers and
zero warnings; scoped `git diff --check` passed; `make fast` converged
after two passes at 2517 pages with zero undefined citations, zero
undefined references, and zero rerun requests; the final log scan found
no fatal TeX errors, undefined controls, unresolved citation/reference
warnings, or rerun requests.

A424 repairs the \(\beta_N\) denominator arithmetic.  The old corollary
claimed that the reduced denominator of \(H_N-1\) is
\(\mathrm{lcm}(2,\ldots,N)\), hence that
\(\operatorname{den}(\beta_N)=\mathrm{lcm}(1,\ldots,N)/12\) for all
\(N\ge5\).  This is false: \(H_6-1=29/20\) while
\(\mathrm{lcm}(1,\ldots,6)=60\), and
\(\beta_{18}=10190221/340340\) while
\(\mathrm{lcm}(1,\ldots,18)/12=1021020\).  The theorem now states the
true reduced-denominator formula
\[
d_N=\operatorname{den}(\beta_N)
   =\operatorname{den}(H_N-1)/\gcd(\operatorname{den}(H_N-1),12),
\]
with divisibility by the lcm quotient, not equality.

A424 also supplies the missing all-\(N\ge5\) non-integrality proof.  By
Bertrand's postulate choose \(p\) with \(N/2<p\le N\), \(p\ge5\).  In
\(L_N\sum_{j=2}^{N}1/j\), the \(j=p\) term is nonzero modulo \(p\) and
all other terms are divisible by \(p\); hence \(p\) survives in the
reduced denominator of \(H_N-1\), and since \(p\nmid12\), also in
\(\beta_N\).  The compute layer now guards the \(N=6\) and \(N=18\)
cancellations, the denominator formula, and the Bertrand-prime
obstruction for \(5\le N<100\).

A424 verification: Python compilation passed for the touched beta
surfaces; the focused beta/W4/WN bundle passed with `46 passed`; stale
exact-denominator greps over active manuscript, compute, and FRONTIER
surfaces were clean; `make verify-licensing` reports zero blockers and
zero warnings; scoped `git diff --check` passed; `make fast` converged
after two passes at 2517 pages with zero undefined citations, zero
undefined references, and zero rerun requests; the final log scan found
no fatal TeX errors, undefined controls, unresolved citation/reference
warnings, or rerun requests.

A425 repairs the \(\beta_N\) claim-status drift.  The harmonic law
\(\beta_N=12(H_N-1)=\sum_{s=2}^{N}12/s\) is no longer advertised as a
proved all-rank closure.  Its label is
`conj:beta-N-harmonic-closed-form`; the separate \(\kappa\)-ratio
scaling input is `conj:kappa-ratio-scaling-law`; and the WN tempering
closure remains the proved finite-rank theorem `thm:wn-tempered-all-N`,
which uses finite \(\beta_N\) plus Stirling dominance rather than the
all-rank harmonic identification.  FRONTIER, the WN closure chapter,
the Hochschild ladder summary, active Vol II notes, and active Vol I
mirrors now distinguish the proved \(N=3\) low-rank value
\(\beta_{W_3}=10\) from the conjectural all-rank harmonic law.

A425 verification: Python compilation passed for the touched beta/WN
tests and the Vol I W3 test.  The focused Vol II beta bundle passed
with `47 passed`; the Vol I W3 invariant test passed with `5 passed`.
Fixed-string stale-claim greps over live Vol II surfaces, active notes,
and active Vol I mirrors were clean.  `make verify-licensing` reports
zero blockers and zero warnings; scoped `git diff --check` passed; and
`make fast` converged after two passes at 2517 pages with zero
undefined citations, zero undefined references, and zero rerun
requests.  The final log scan found no fatal TeX errors, undefined
controls, unresolved citation/reference warnings, or rerun requests.

A426 removes the last live `\ClaimStatusNeedsVerification` marker from
the manuscript.  The deleted claim "\(E_3\)-PBW proves concentration"
is now classified as `\ClaimStatusConjectured`, because the precise
missing theorem is already named:
`conj:H-concentration-via-E3-rigidity`, conditional on the Rees-flat
chiral-\(E_3\)-PBW package for the derived centre.  The established
concentration route remains the ordered-bar Arnold--Orlik--Solomon/FM
mechanism.  The misleading `thm:H-concentration-via-E3-rigidity` alias
has been removed from the CHD conjecture; only the `conj:` label remains.

A426 verification: Python compilation passed for the deletion-ledger
and CHD guard surfaces.  The focused deletion-ledger/CHD bundle passed
with `67 passed`.  Live-surface greps found no `NeedsVerification`
marker and no `thm:H-concentration-via-E3-rigidity` alias in
`chapters/`, `compute/lib`, or `FRONTIER.md`.  `make verify-licensing`
reports zero blockers and zero warnings; scoped `git diff --check`
passed; and `make fast` converged after two passes at 2517 pages with
zero undefined citations, zero undefined references, and zero rerun
requests.  The final log scan found no fatal TeX errors, undefined
controls, unresolved citation/reference warnings, or rerun requests.

A427 repairs the higher-genus Hochschild bridge.  The genus-\(g\)
comparison is no longer described as a filtered quasi-isomorphism of
curved complexes.  The manuscript now defines a filtered
\(R_g\)-curved comparison, states
\[
\Psi_g\colon
\mathrm{C}^{\bullet}_{\mathrm{ch,top}}(\cA;\Sigma_g)
\longrightarrow
\ChirHoch^\bullet_g(\cA;\Sigma_g),
\]
and restricts ordinary quasi-isomorphism language to the
\(\omega_g=0\) special fibre or to the zero-curvature case.  The
curvature obstruction is typed as
\(\Theta_g\in H^2_{\mathrm{cyc}}(\cA,\cA)\), not as cohomology of a
curved Hochschild complex.

A427 also removes the surviving full-bridge overclaim.  Genus \(0\)
is an ordinary dg Hochschild bridge; genus \(g\ge1\) is a filtered
curved \(R_g\)-comparison unless the \(K=0\) or
\(\kappaChHodge(\cA)=0\) hypothesis removes the curvature.  The local
stalk proof now uses a holomorphic coordinate on a small analytic disc
and formal completion \(\operatorname{Spf}\C[[z]]\), not
uniformization of the disc by all of \(\C\).  The
anomaly-completed-topological-holography heuristic now invokes the
filtered curved comparison and explicitly excludes an ordinary
quasi-isomorphism at nonzero curvature.

A427 verification: Python compilation passed for the new curved-bridge
guard.  The focused Hochschild/genus bundle passed with `117 passed`.
Stale-phrase greps over active Vol II, Vol I, and Vol III sources found
no live copy of the retired curved-quasi-isomorphism or uniformization
claims; remaining Vol I hits are archived audit snapshots.  `make
verify-licensing` reports zero blockers and zero warnings; scoped
`git diff --check` passed; and `make fast` converged after two passes at
2517 pages with zero undefined citations, zero undefined references,
and zero rerun requests.  The final log scan found no fatal TeX errors,
undefined controls, unresolved citation/reference warnings, or rerun
requests.

A428 repairs the Zhu bridge.  The manuscript no longer treats
\(A(V)\) as a quotient of the full mode algebra and no longer claims
that rational \(C_2\)-cofinite \(V\) gives a quasi-isomorphism
\(\ChirHoch^{\le2}(V)\to\HH^{\le2}(A(V))\).  The live theorem is now
the degree-zero centre shadow
\[
\mathrm{Zhu}^0\colon Z_{\mathrm{ch}}(V)\to Z(A(V)).
\]
It states that \(A(V)\simeq \mathcal U(V)_0/I_{\mathrm{Zhu}}\), so
the zero-mode assignment is not an algebra homomorphism on the full
mode algebra.

A428 also installs the rational counterexample and the critical-level
direction.  For simple rational \(C_2\)-cofinite \(V\),
\[
A(V)\simeq\bigoplus_{i\in I(V)}\End(M_i(0)),\qquad
\HH^0(A(V))\simeq\C^{I(V)},\qquad \HH^j(A(V))=0\ (j>0).
\]
The simple chiral centre maps to the scalar diagonal; for the Ising
model this is the strict map \(\C\to\C^3\).  At critical affine level,
the Feigin--Frenkel centre maps onto \(Z(U(\fg))\) by the
Feigin--Frenkel--Reshetikhin one-point/Gaudin specialization followed
by Harish--Chandra projection; the kernel is infinite-dimensional.
The lost chiral information is kernel/higher Hochschild data, not an
infinite-codimension cokernel.

A428 verification: Python compilation passed for the new Zhu guard.
The focused Hochschild source bundle passed with `59 passed`.  Stale
greps over active Vol II, Vol I, and Vol III surfaces found no live
copy of the retired rational-quasi-isomorphism, full-mode quotient,
critical inclusion, or infinite-codimension-cokernel claims.  `make
verify-licensing` reports zero blockers and zero warnings; scoped
`git diff --check` passed; and `make fast` converged after two passes
at 2517 pages with zero undefined citations, zero undefined references,
and zero rerun requests.  The final log scan found no fatal TeX errors,
undefined controls, unresolved citation/reference warnings, or rerun
requests.

A429 repairs the negative-cyclic/K-theory stage claim.  The manuscript
no longer says that \(\mathrm{HC}^{-}(A)\) computes \(K\)-theory or
that a Bridgeland orientation identifies negative cyclic homology with
chiral \(K\)-theory.  The four-object Hochschild discipline now states
that \(\mathrm{HC}^{-}(A)\) is the negative-cyclic target of the Chern
character from \(K\)-theory; Bridgeland data orient the trace, choose a
chamber, and supply filtration data, but do not create an absolute
equivalence.

A429 also reverses the displayed comparison to the correct direction
\[
K^{\mathrm{ch}}(\mathrm{Mod}(A))
\xrightarrow{\operatorname{ch}^{-}_{\sigma}}
\mathrm{HC}^{-}(A)[\sigma].
\]
The theorem now separates the level-\(\mathsf A\) acting object
\(K^{\mathrm{ch}}(\mathrm{Mod}(A))\) from the level-\(\mathsf S\)
negative-cyclic trace shadow \(\mathrm{HC}^{-}(A)[\sigma]\).  It also
records the base-field obstruction: for \(A=\C\),
\(K_0(\C)=\Z\) whereas \(\mathrm{HC}^{-}_0(\C)\) is a complex vector
space containing \(\C\), and \(K_1(\C)=\C^\times\) whereas
\(\mathrm{HC}^{-}_1(\C)=0\).  The unsupported BZN/HKKP
K-theoretic-Hochschild-Bridgeland attribution has been removed from
the proof.  The architectural critique note now records the same
Chern-character direction \(K\to\mathrm{HC}^{-}\).

A429 verification: Python compilation passed for the new negative
cyclic guard.  The focused Hochschild source bundle passed with
`10 passed`.  Fixed-string stale scans over active Vol II, Vol I, and
Vol III surfaces found no live copy of the retired negative-cyclic
equivalence; remaining hits are the negative guard and the attack note.
`make verify-licensing` reports zero blockers and zero warnings; scoped
`git diff --check` passed; and `make fast` converged after two passes at
2517 pages with zero undefined citations, zero undefined references,
and zero rerun requests.  The final log scan found no fatal TeX errors,
undefined controls, unresolved citation/reference warnings, or rerun
requests.

A430 repairs the DS--Hochschild transfer route.  The manuscript no
longer treats the class-\(\mathsf M\) bridge as a consequence of a
generic finiteness slogan or of applying Hochschild cochains to a DS
algebra retract.  The live theorem path now requires the completed DS
bar-coalgebra SDR on reduced chiral bar coalgebras, the induced
transfer on \(\Coder_0(B^{\mathrm{ch}}(-))\), the smooth
chiral-HKR/derived-vacuum comparison, and the bounded-shift
weight-completed HPL convergence estimate.  The finite input is the
Kazhdan/PBW finite-weight profile; admissible/lisse \(C_2\)-finite
loci remain separate finite surfaces.

A430 propagation tightened the THQG reconstruction, Universal
Celestial, Universal Holography Functor, climax deletion/rung, bar-cobar
review, SC-heptagon, theorem-D ambient, and first-principles cache
surfaces.  The UHF DS bridge is now
\(\ClaimStatusProvedHereConditional\) with the hypothesis package in the
proposition.  Generic "Hochschild functoriality" language was replaced
by the Hochschild/coderivation construction plus an explicit comparison
datum; THQG now states that the closed-sector map belongs to the
\(\Xi\)-compatible morphism data and is not automatic functoriality of
cochains along an algebra map.

A430 verification: Python compilation passed for the DS/Hochschild guard
files.  The focused bundle
`test_ds_hochschild_scope.py`, `test_chd_ds_hochschild_iv.py`,
`test_chiral_higher_deligne.py`, and
`test_chd_class_m_holography_gates.py` passed with `17 passed`.
Fixed-string stale scans over active Vol II surfaces, `FRONTIER.md`,
and the first-principles cache found the retired route only inside
negative regression tests.  `make verify-licensing` reports zero
blockers and zero warnings; and `make fast` converged after two passes
at 2517 pages with zero undefined citations, zero undefined references,
and zero rerun requests.  The final log scan found no fatal TeX errors,
unresolved citation/reference warnings, or rerun requests.

A431 repairs the three-Hochschild unification scope attacked in F8 of
the 2026-06-10 Hochschild findings.  The theorem no longer says that
the mode-comparison map is a degree-zero isomorphism "always".  It is a
degree-zero isomorphism only on the non-critical scalar-centre locus; at
critical affine level \(k=-h^\vee\) the Feigin--Frenkel zero-mode
specialisation has infinite-dimensional kernel in the uncompleted mode
ambient.  The proof and scope remark now state that this is not an
ordinary mode-Hochschild isomorphism.

A431 also removes the false Weyl-algebra isomorphism claim for the
Gelfand--Fuks comparison.  The map \(\zeta\) is now a CE/HKR character
to the scalar Hochschild class coming from the deformation-quantisation
trace.  For a Weyl mode algebra,
\(HH^\bullet(A_{\mathrm{mode}},A_{\mathrm{mode}})=k\) in degree zero;
positive-degree Gelfand--Fuks classes lie in the kernel.  The
first-principles cache records the same critical kernel and distinguishes
ordinary uncompleted mode Hochschild from any completed critical
enveloping variant.

A431 verification: Python compilation passed for the Hochschild scope
guard bundle.  The focused tests
`test_three_hochschild_unification_scope.py`,
`test_zhu_bridge_scope.py`,
`test_chiral_hochschild_coderivation_model.py`, and
`test_hochschild_bulk_bridge.py` passed with `58 passed`.
Fixed-string stale scans found the retired claims only inside the
negative guard.  `make verify-licensing` reported zero blockers and
zero warnings; and `make fast` converged after two passes at 2517 pages
with zero undefined citations, zero undefined references, and zero
rerun requests.

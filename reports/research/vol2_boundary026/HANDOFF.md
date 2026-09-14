# Volume II: compact quantum current contraction

Candidate001 constructs the next supported chain comparison after the scalar bar calculation. It does not repeat that calculation. Independent acceptance remains required.

## Frozen question and boundary

The assignment asks for the first missing actual map from the higher-operation/current construction to its native boundary consumer. The exact composed002 higher-operation chapter ends by requiring a comparison with the full quantum boundary differential, collision operations, and bulk action. Its extracted-field interval contraction does not provide that comparison.

The first selected obligation is a chain map for the actual positive quantum differential, with compact primitives for the separated-product comparison. The construction supplies that map, a strong deformation retract of an actual supported current subcomplex, every polynomial Wick primitive, a split inclusion into the full smooth complex, and the native interval consumer. Compatible higher Taylor maps and collision operations remain open. This is a bounded construction, not completion of the full open–closed target.

## Source custody

- Worktree: `/Users/raeez/mathematics/worktrees/frontier-vol2-boundary-026-20260914`.
- Branch: `intake/frontier-vol2-boundary-026-20260914`.
- Base HEAD and native repository HEAD: `7c0a1c3203b38dd47281d2dffd607f3d880e6056`.
- Native root is read-only. No native checkout or central PDF was changed.
- The actual constructed BF/higher-operation closure is the composed002 source, located through the assigned old020 records. Its complete two chapters were read. Historical status prose was not used as proof.
- The pre-edit BF source hash is `cc69c1ff2bb37c9bab9f80c4b3f2f867de73033ac2923d6c0b1b73a855a87f91`.
- The pre-edit higher-operation source hash is `6b3ee1938d16a7a01e7d0b6c2aaf43fc27ee329dab95d7bf2fce0b5e1ed67dfe`.
- `source-freeze.json` and `baseline/` preserve the exact inputs before edits.
- No staging, commit, push, cleanup, or publication occurred. No child delegation occurred.
- The operational requirement is `gpt-6-astra` with `ultra`. Independently observed runtime controls were unavailable and remain unverified.

## Construction and deciding case

All cohomological degrees are stated in the source. The ground field is complex numbers. The parameters are polynomial in hbar initially. A finite list of distinct point/jet pairs may include several jets at one center. Different centers have disjoint small disks. One radial cutoff equals one on a centered disk containing all smaller disks.

The actual compact representatives are

\[
u_i=-\frac{\bar\partial\chi_{a_i}}{(z-a_i)^{n_i+1}},\quad
e_i=-\frac{\bar\partial\chi}{(z-a_i)^{n_i+1}},\quad
f_i=\frac{\chi-\chi_{a_i}}{(z-a_i)^{n_i+1}}.
\]

They satisfy `dbar f_i = u_i - e_i`. Forms have observable degree zero. Functions have degree minus one. The residue normalization is `1/(2 pi i)`. The bulk parameter remains kappa. The extracted current level is minus kappa.

The proof calculates `P(u_i,u_j)=p_ij`, `P(e_i,u_j)=P(e_i,e_j)=0`, and `mu(f_i,u_j)=p_ij`, `mu(f_i,e_j)=0`. Here

\[
p_{ij}=\frac{\kappa}{n_i!n_j!}
\partial_{a_i}^{n_i}\partial_{a_j}^{n_j}(a_i-a_j)^{-2}
\]

for distinct centers, and zero at a common center. The common-center value uses radial representatives and rotation weights, not a singular specialization of the rational function.

Introduce even variables `x_i,y_i` and odd variables `theta_i`, represented by `u_i,e_i,f_i`. The actual supported subcomplex is

\[
Q=\mathbb C[x,y]\otimes\Lambda(\theta)[\hbar],\qquad
D_Q=\sum_i(x_i-y_i)\partial_{\theta_i}
+\hbar\sum_{i<j}p_{ij}
(\partial_{x_j}\partial_{\theta_i}
+\partial_{x_i}\partial_{\theta_j}).
\]

Left odd derivatives fix the signs. Linear independence of the compact functions and forms makes the substitution map injective. Continuous dual coordinate functionals then establish injectivity in the smooth tensor convention. This constructs a subcomplex of the native observable complex, rather than changing the native carrier silently.

Set `delta=x-y`, `A=C[y][hbar]`, and

\[
d_0=\sum_i\delta_i\partial_{\theta_i},\quad
K_Q=\sum_{i<j}p_{ij}\partial_{x_i}\partial_{x_j},\quad
N_Q=e^{-\hbar K_Q}.
\]

On a monomial of positive total delta/theta degree q, set

\[
h_0=\frac1q\sum_i\theta_i\partial_{\delta_i}.
\]

Set it to zero in degree q=0. The supercommutator with d0 is the Euler derivation. Thus `d0 h0+h0 d0=1-i pi`. Also `D_Q=N_Q^{-1}d0N_Q`. Conjugation gives

\[
\Pi=\pi N_Q,\qquad H=N_Q^{-1}h_0N_Q,\qquad
D_QH+HD_Q=1-\iota\Pi.
\]

All side conditions are proved: `Pi i=1`, `Hi=Pi H=H^2=0`. Every exponential terminates on polynomials. No inverse of kappa or hbar occurs. The argument therefore includes zero level. It works in all polynomial degrees, all finite jet lists, and their fixed-center directed union. Coefficientwise hbar completion of this specified subcomplex preserves the identities.

The first deciding higher case is cubic:

\[
H(x_1x_2x_3)=
\sum_i\theta_i\left(y_jy_k+
\frac{\delta_jy_k+y_j\delta_k}{2}+
\frac{\delta_j\delta_k}{3}-\frac{2\hbar p_{jk}}3\right),
\]

where `{j,k}` is the complement of i. Omitting the last term leaves

\[
\frac{2\hbar}{3}(p_{12}\delta_3+p_{13}\delta_2+p_{23}\delta_1).
\]

This is nonzero whenever the covariance is nonzero. The binary primitive cannot detect this correction. The zero-level case was checked separately and has no quantum correction.

For every polynomial current product F, the full boundary chain identity is

\[
F(u)-\iota_U\left(\left.e^{-\hbar K_Q}F(x)\right|_{x=y}\right)
=D\rho HF(x).
\]

The expression in parentheses is the finite matching sum with contraction `-hbar p_ij`. This proves every polynomial separated-product comparison before passing to cohomology. The original compact collar chain map carries these exact primitives to the bulk–boundary pushforward.

## Full boundary split inclusion and interval consumer

For each finite point/jet list, finite Hermite interpolation supplies holomorphic polynomials F_i with divided jets `F_i^(n_j)(a_j)/n_j! = delta_ij`. The proof constructs them by injectivity and dimension of the full finite jet evaluation map. Moment integration defines a continuous classical chain map R on all compact smooth kernels. The quantum map is `L_U=R N_U`.

It satisfies `L_U i_U=1` on the annular polynomial algebra, and `L_U rho=Pi` on Q. Thus the supported polynomial states split into the cohomology of the full smooth boundary complex. The claim is injection and a chain retraction onto this finite sector. It is not a quasi-isomorphism of the whole smooth boundary complex.

The native higher-operation consumer now defines

\[
\Phi_a(\alpha\otimes F)=\epsilon_a(\alpha)\iota_U(F),\qquad
\Psi(b)=(e_0+e_1)\otimes L_U(b).
\]

Both maps intertwine the actual differentials. Their composites are calculated, including the old interval homotopy k_a. Tensoring the compact contraction with the interval uses

\[
\mathscr H(\alpha\otimes q)=(-1)^{|\alpha|}\alpha\otimes Hq.
\]

The mixed differential terms cancel. This gives the explicit first map from the finite current part of the interval model into the full boundary complex. Undivided jets are `n_i! y_i`, preserving the higher-operation chapter's convention. Its physical parameter is k, so the BF specialization is kappa=-k.

The consumer explicitly distinguishes this map of complexes from an ordered higher morphism. The latter still requires compatible Taylor maps on the same supported complexes. The target's ordinary multiplication fails the Leibniz identity on overlapping supports. No strict multiplicative comparison is asserted there.

## Exact source and dependencies

The source aggregate is `594ed4b6bd0008a48c5dcd5715fcb02b173b646ae2c84c40a5b89189b912cf17`.

| Path under `research-candidates/vol2_boundary026/` | SHA-256 | Role |
| --- | --- | --- |
| `chapters/compact-current-homotopies.tex` | `fbbd7da519ee2ceb36c39bc31076f7e7d67ad0d445dffc8cc22507364e64eda8` | Complete new proof, 334 lines |
| `chapters/zero-level-bf.tex` | `0ba02032127976cbca0bd1d98b8d3ff6ffb99b8da06b6ba33bfafc71477725b9` | Preserved BF chapter with one input insertion |
| `chapters/Volume_II_Free_Field_Higher_Operations.tex` | `82981ecc190f315a21a235a9959c7ad2fb110e06e42e65f1b51716757d1e214c` | Existing higher operations plus actual interval chain-map consumer |

The BF input insertion is immediately before `Fields, modes and regular products`. Its new labels are `sec:bf-compact-current-homotopies`, `lem:bf-compact-representative-pairings`, `prop:bf-compact-current-complex`, `thm:bf-compact-quantum-contraction`, and `prop:bf-compact-full-boundary-split`. The higher chapter's new consumer is `prop:interval-supported-current-chain-map`, at source line 654.

The source dependency chain is:

1. Native compact boundary differential and smooth-kernel carrier.
2. Native Cauchy pairing, continuous contraction, and finite quantum conjugation.
3. Native separated divided-jet contraction formula.
4. New representative pairings and actual finite supported subcomplex.
5. New Euler contraction, conjugation, all-polynomial primitives, and moment split inclusion.
6. Existing interval contraction and new native chain-map consumer.
7. Existing compact collar map carries the proved identities into its pushforward.

The complete two-chapter include closure is in `reader.tex`. It uses the actual shared template symlink. The root compositor should integrate the three mathematical source paths and recheck references in its own exact candidate. The source has no new external theorem dependency. As a primary-source cross-check, the regular-sequence Koszul vanishing agrees with [Stacks Project Lemma 15.31.2, tag 062F](https://stacks.math.columbia.edu/tag/062F). The explicit Euler homotopy proves the stronger chain statement internally. No novelty claim is made.

## Checks and limitations

`verify_quantum_contraction.py` uses Python 3.14.6 and exact rational arithmetic. It independently encodes the differential's odd signs and pair contractions, then checks the proposed contraction. It checks 672 monomials at each of three bulk levels, for a total of 2016. The delta degree is at most 6, with all eight odd masks. The tests cover `D^2=0`, `ND=d0N`, `DH+HD=1-iPi`, all side conditions, the displayed binary and cubic formulas, and the nonzero omission defect. The code checks these finite ranges only. The proof supplies all finite lists and degrees.

The internal reader builds to 35 pages with pdfLaTeX. The final two deterministic passes have identical PDF and AUX hashes. There are zero undefined references, undefined controls, duplicate labels, overfull boxes, or underfull boxes. One epstopdf warning says shell escape is disabled. No EPS conversion is required. All 324 final recorder inputs, including fonts, packages, template, and stable auxiliaries, are archived.

Physical pages 14–20 and 33–35 were rendered and inspected. These include the complete new proof, preceding and following transitions, and the full new consumer. Final pages 15–19 were inspected after the last hypothesis edit. Page 35 was inspected after the wording that separates strict field products from supported homotopies. No clipping, overlap, broken equations, or process prose was observed. This is constructor visual checking, not independent acceptance of the full inherited 35-page source.

The native `make platonic` command was not run because it would build the unmodified 495-page worktree entrypoint outside the assigned source subtree. The focused build includes the full changed dependency closure. No changed source uses a `ProvedHere` status macro. The legacy independence registry does not scan this candidate directory. It is not treated as theorem verification.

The manuscript and PDF text were scanned for operational vocabulary and closed-repository paths. There were no hits in the new proof or native consumer. The reader title and metadata contain mathematical content only. No standalone PDF was opened. Rendered PNGs were inspected internally. The central combined reader remains root-owned.

## Failed routes and execution failures

The unchanged classical homotopy h0 is refuted at the cubic input when covariance is nonzero. Its exact residual is displayed above and verified by the executable check. Conjugating both sides repairs that failure.

A strict ordinary differential graded algebra interpretation of Q is refuted by `D(theta_i x_j)-(D theta_i)x_j=hbar p_ij`. The construction preserves the native second-order differential and uses chain maps and compact homotopies. It does not rename a transported product as the native factorization product.

The first computation attempt used `/usr/bin/python3` version 3.9.6. It failed at `int.bit_count`, before any mathematical verdict. The successful command explicitly uses `/opt/homebrew/bin/python3` version 3.14.6. No test assertion was weakened.

The initial LuaLaTeX invocation failed because its output directory did not yet exist. The next LuaLaTeX build reached an inherited `llbracket` command and failed because that engine path does not define it. Both logs are retained. The native pdfLaTeX path succeeded without a mathematical-source notation workaround.

The live shared template hash is `4bd3e0dcea54b2209a2d02d5964e96ae2f11700dbe3b44fb0baf4e61d7352de3`. It differs from the composed002 template hash `07336dc2503195a619a5b566e8730e891d4c81b2465530873214d7087b83e4b5`. The template was not edited or replaced. The actual consumed bytes are archived. This candidate's render therefore does not inherit the prior template verdict.

## Residual obligations

- Construct compatible Taylor maps for an ordered higher morphism on the same full supported complexes.
- Define collision-chain operations and prove their compatibility with boundary strata.
- Construct the corresponding bulk action comparison at that higher level.
- Prove any proposed comparison after completion of the full smooth observable complex. Completion here applies only to the specified retract.
- Construct the relevant renormalized configuration integrals and interacting quantum master equation when claimed.
- Obtain fresh independent exact-byte review. This handoff is not an acceptance record.

## Handoff artifacts

`candidate001-manifest.json` has SHA-256 `c5f0f333c957cfcd65360add0c9471efac81f54c1865d907c5d6dd78e2d4fe0d`.
`native-delta.patch` has SHA-256 `8968f52893172911385fdc540980f6544a613bb524b834623f886ca04667cbae`.
The 35-page internal PDF has SHA-256 `1e20034c2d80bfc1a7452604e288bb24933f49ecf97c9db79c88a4cb85b5b8b6`.
The manifest records the source archive, compiler-input archive, executable checks, build commands, and hashes.

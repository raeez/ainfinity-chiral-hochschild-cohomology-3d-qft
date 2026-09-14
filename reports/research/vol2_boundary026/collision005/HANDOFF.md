# Three-current residue comparison

## Result and limit

The new body constructs an actual chain map for three elementary currents, including every elementary holomorphic jet. Its source contains compact Dolbeault test chains on the three-point insertion space and each pair diagonal. The residue traces enter the source differential with coefficient `-lambda`, where `lambda=-hbar*kappa`. The target is the original compact smooth boundary complex, or its polynomial interval tensor product. All outputs keep the fixed annular support.

For all three polynomial insertions, the body proves the maximal domain with the unchanged Dolbeault test differential and the prescribed Wick smearing. Its degree-minus-one part is exactly the kernel of the observable-valued contact map. This maximality permits arbitrary negative-degree comparison values and changes of the degree-zero comparison by exact observables.

The requested regular-product ternary primitive cannot exist on the stated boundary carrier at nonzero level. Endpoint evaluation and the `z^2` moment detect `2 lambda q3`. The regular product therefore cannot be the induced binary product of an uncurved A-infinity structure on these injected cohomology classes. This is a proved obstruction to that specified requirement, not a nonexistence theorem for collision resolutions with additional geometric source data.

The full composite collision comparison remains unresolved. In particular, the three quadratic insertions require a differential between contraction strata compatible with the triangle contact current. The elementary residue complex does not prove that differential or its comparison with iterated regular Laurent coefficients. It does not construct a nonzero higher ordered map into bare boundary observables.

## Frozen candidate

The source directory is `research-candidates/vol2_boundary026/collision005/`. The report directory is this directory. Both are inside the assigned boundary worktree.

The exact review entrypoint is `manifest.json`, SHA-256 `6a41dbd3f5988eb823d803a0a02a42eb675376f494062f6d4903af524d452cad`.

The new body is `chapters/collision-chains.tex`, SHA-256 `d6e160ede05e04b03e05e8e6f774f206b7bc5bcb98770932cfeb5e63af10d983`. The complete source aggregate is `b0149fa032e8cb51920a23ad9cf61ce26aed71771d51ff872401b2aab1a9cb05`.

The local PDF is `build-pdftex/reader.pdf`, with 60 pages and SHA-256 `84589938f68b5a125dbc4e35048e5f38e8f3d2e3a1f1c8c8486e1e355f4360aa`. The new body occupies physical pages 54–59. Bibliography is page 60. The complete input closure, including 337 compiler inputs, is recorded and archived locally.

All seven pure mathematical bodies in `integration-packet-v2.json` retain their supplied SHA-256 values. That packet and all three exact review reports were hash-checked. `input-freeze.json` records the source copies. No source from quantum001, collision002, higher-retract003, ordered026, or collision-domain026 was modified. A filename search restricted to Volume II returned no candidate004 or collision004 partial in the current worktree collection.

The reader adds only the new mathematical chapter and its integration entry. The three bibliography entries have primary-source links. `native-delta.patch` contains the new chapter as an addition. This is a candidate integration closure, not a central manuscript edit.

## Proof anchors

| Claim | New body lines | Decisive argument |
| --- | --- | --- |
| Typed compact Dolbeault test differential | 46–66 | Transpose convention and anticommuting exterior contractions |
| Residue trace and contact identity | 68–120 | `r=pi trace(partial_x contraction_dbarx)` and `delta_Z r+r delta_X=0` |
| Three-current source and target map | 124–205 | Off-diagonal differential `-lambda r`; exact cancellation in degree minus one |
| Arbitrary elementary jets | 207–230 | Explicit pole derivative and transpose factorial cancellation |
| Annular detection in the full boundary complex | 261–306 | Convergent annular series, Cauchy conjugation, and every finite moment projection |
| Maximal unchanged test subcomplex | 308–359 | Exactness of the contact output implies its literal vanishing |
| Nonzero-level ternary obstruction | 361–422 | Constant Laurent coefficient, endpoint evaluation, and nonzero `z^2` moment |

The source itself distinguishes the pair-diagonal source from the original quotient test complex. It does not conceal the added source differential. The maximal-domain theorem separately addresses the unchanged full polynomial source.

The elementary complex has no second contraction term because a contraction of two of three elementary fields leaves one field. That explanation does not apply to composite insertions. The copied distributional triangle and its contact formula remain available as exact inputs to that next construction.

## Verification

Run `/opt/homebrew/bin/python3 reports/research/vol2_boundary026/collision005/check_residue.py` from the assigned worktree. `calculation.json` records Python 3.14.6 and SymPy 1.14.0.

The exact calculation passes 32 test-differential square checks, 480 residue anticommutation checks, and 25 jet normalizations. It transforms all three pair-diagonal charts from the original insertion coordinates. It checks exterior degrees zero through three and five different jet pairs. It also recomputes the regular associator `2*lambda*q3` and the normal test germ with nonzero residue. The script does not import any prior verification helper.

These checks verify the finite polynomial jets and signs. The source proofs establish the statements for all compact smooth tests. No finite calculation certifies a composite collision-stratum algebra.

The reproducible LaTeX invocation and environment overrides are in `build-command-final.json`. Three final pdflatex passes exited zero with identical PDF and AUX bytes. The final log has no overfull, underfull, undefined-reference, multiply-defined-reference, or missing-character diagnostic. `verification-final.json` records the checks.

An initial integration wrapper unnecessarily loaded `stmaryrd`. Combined with the full seven-body closure, that exceeded pdfTeX's math alphabet limit. The new wrapper removed that unused package. The established wrapper configuration then built successfully. The old bodies and shared template were unchanged. An initial calculation under the system Python failed because SymPy was absent; the explicit Homebrew Python command above resolved that environment issue.

The new source, rendered pages 54–59, bibliography, and mathematical PDF metadata were inspected. No manuscript-firewall violation, clipping, or overlap was identified. The source scan and rendered closed-path scan returned zero hits.

## Reader reference gate

Reader acceptance remains blocked by the existing shared-template equation defect. All six new equation labels inherit section, lemma, or proposition numbers and destinations. The displays lack their intended equation tags. `verification-final.json` records the exact AUX entries. This is not repaired by the absence of LaTeX warnings.

The template symlink remains owned by the shared template. This task made no template changes. A corrected template requires a new exact build, manifest, and render check. The current PDF is a working artifact and was not copied to a central accepted location.

## Primary source and independence

Keller, *Introduction to A-infinity algebras and modules*, section 3.6, pages 11–12, was checked in the author's PDF. It supplies the relation between suspended coderivations and ordinary operations, including degree `2-n`. The local copy and extracted text are under `materials/raw/`. The arXiv record was also checked. The new construction has direct proofs; Keller does not supply the compact residue source or its support statements.

The underlying Brouder–Dang–Hélein and Brunetti–Fredenhagen source records were checked against their arXiv entries. Their precise theorem uses remain the frozen precursor's responsibility. No imported source is used to infer a residue-current multiplication absent from the written construction. No novelty claim is made.

No fresh independent review of collision005 has yet been received. The source self-check and calculations are diagnostic. Exact review must cover the same manifest, especially the maximal annular-detection argument and the scope of the ternary obstruction.

The programme requires gpt-6-astra with ultra effort for mathematical work. Confirming runtime metadata was unavailable. Compliance remains unverified; no mismatch was observed or inferred.

## Residual obligations and custody

The next mathematical obligation is a composite collision-stratum source carrying the actual E3 triangle contact. Its differential must square to zero and its map must preserve compact output support. A comparison with regular Laurent coefficients must retain the nonzero `2 lambda q3` class. The unchanged interval differential cannot remove that class.

The separated ordered Taylor construction retains fixed disjoint disk unions. Its nonzero higher terms have target `A_I tensor B(U_S)` and vanish after endpoint evaluation. Neither the present residue map nor the maximal-domain theorem changes that accepted boundary.

All writes stayed in the two assigned collision005 directories. No staging, commit, push, branch change, central PDF write, standalone PDF application opening, or descendant delegation occurred. The integration owner retains all release authority.

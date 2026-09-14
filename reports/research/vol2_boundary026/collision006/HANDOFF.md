# Pure-template integration for the three-current reader

## Result

The complete 60-page reader builds with the exact reviewed pure template. All 30 equation labels have distinct equation destinations and matching visible tags. All 34 equation-reference occurrences have the correct visible number and PDF link target. Local reference and render checks pass.

The eight mathematical bodies are byte-identical to collision005. Its scoped independent mathematical PASS therefore remains attached to the same source bytes. This integration still requires fresh review of its exact reader candidate. It is a working artifact, not an accepted release.

The reader wrapper marks all 30 existing equation labels with the supported `noeqref` command. This also gives unreferenced labels their own equation numbers and destinations. No label or proof was edited in a mathematical chapter. The template remains a symlink to the chosen frozen source. No global style file or existing consumer symlink was changed.

## Exact bindings

| Object | SHA-256 |
| --- | --- |
| `manifest.json` | `937af5c9f24e8548b8b397b220715fc0bd92640b33a15066bcad1b546d4ead26` |
| Complete source aggregate | `1c25afa9b2781b63c3303a9d8ab9a8e175acb997ceddc415d9e69788d898320b` |
| Working PDF, 60 pages | `243f0c0f4da7195ed52d9b4e2d33bc8ed128c0cf50c53c99ff39aaf8e9203ba4` |
| Chosen pure template | `3f435a5fc820e94363a299ffa47d884c1e158508b74e80ce478822b528894bc6` |
| Unchanged collision-chain body | `d6e160ede05e04b03e05e8e6f774f206b7bc5bcb98770932cfeb5e63af10d983` |
| Collision005 mathematical review | `e6613e37b80989a4594aa0b4ba5bae051e7c8dfb2329c0bddbb9f61e93fc376b` |
| Frozen independent template review | `63f348578cff7af12b6eb6d1b00308a22c3114023f49b2700610d212dcc64f83` |

The source entrypoint is `research-candidates/vol2_boundary026/collision006/reader.tex` in the assigned boundary worktree. The local PDF is `reports/research/vol2_boundary026/collision006/build/reader.pdf`.

The chosen template source is `candidate/source/raeez-math-template.sty` in the template owner's frozen equation-reference packet. Its author manifest hash is `4c9d643639048366a26adea3b2dbae42cf31109f17673f598113e986d48325fa`. The operational template with hash prefix `84fedc` was not used. The pure package has additional inherited API differences from the former operational `4bd3` package. The full consumer build and text comparison check their effect on this reader only.

`source-freeze.tar.gz` contains the complete source closure. `compiler-inputs.json` and its archive preserve all 337 non-output compiler inputs. `integration-delta.patch` compares the wrapper and selected template against collision005. The manifest records base commit, branch, each source hash, aggregate diff hash, build command, review inputs, and PDF identity.

## Reference and content checks

`check_reader.py` reads the source labels and reference occurrences. It checks AUX numbers against every actual PDF tag, named destination, destination page, and reference annotation. It also checks that no tag overlaps another extracted word.

The numbering is consecutive within each chapter: `1.1–1.3`, `2.1–2.2`, `3.1–3.8`, `4.1–4.11`, and `5.1–5.6`. All 30 destinations are distinct `equation.*` targets. The equation-reference count is 34, matching the literal source occurrences. All 45 non-equation labels retain their numbers, pages, and titles.

The whole-reader text comparison removes only whitespace, numeric internal-reference rectangles, equation tags, and the header/footer bands. The remaining difference consists of three inserted line-break hyphens: `correspond-ing` on page 53, `pre-cisely` on page 55, and `con-traction` on page 56. Raw PDF lines independently confirm each split. No other character difference remains. `reader-checks.json` records the exact scope and contexts. The source-byte comparison separately proves that all eight mathematical bodies remain unchanged.

All pages containing equation tags were visually inspected: 16, 17, 35, 36, 41–45, 47–52, 54–56, and 58. The remaining collision proof pages, its transition, and bibliography were also inspected: 53, 57, 59, and 60. Thus 23 page PNGs were inspected individually. Tags and proofs are legible. No clipping, overlap, missing mathematical glyph, or manuscript-firewall violation was identified.

## Build and reproduction

`build-command.json` records four successful pdflatex passes from the exact isolated source directory. The final two PDF hashes agree. The last three AUX hashes agree. The command uses no shell escape, recorder output, and fixed source-date variables.

Run the artifact checker from the assigned worktree:

    /opt/homebrew/bin/python3 reports/research/vol2_boundary026/collision006/check_reader.py

It uses PyMuPDF 1.27.2 and pypdf 6.14.2 with Python 3.14.6. The build uses pdfTeX 1.40.27 from TeX Live 2025. No mathematical test was rerun because no mathematical body changed.

The final build has no overfull, underfull, undefined-reference, undefined-citation, duplicate-label, missing-character, or duplicate-destination diagnostic. Two inherited warnings remain: disabled shell escape for epstopdf, and amsrefs' preferred citation syntax. Neither warning affects the checked reference behavior. The output recorder confines all generated TeX files to the owned build directory.

## Scope and custody

The mathematics remains the accepted elementary three-current residue complex with all elementary jets, the all-polynomial maximal domain for the unchanged test differential, and the nonzero regular-product ternary obstruction. Composite collision-stratum maps remain unconstructed. This reader integration does not add a collision algebra or a nonzero higher map into bare boundary observables.

The programme requirement remains gpt-6-astra with ultra effort for mathematical work. Runtime controls were not exposed for inspection in this task and remain unverified. No configuration mismatch was observed or inferred.

All writes stayed under the two collision006 roots. No staging, commit, push, global style write, central PDF write, standalone PDF application opening, or descendant delegation occurred. Collision005 and all older source candidates remain unchanged. The reserved Rees successor was not started, and no moving Rees source was imported.

Fresh integration review must use the exact manifest above. The integration owner retains release authority.

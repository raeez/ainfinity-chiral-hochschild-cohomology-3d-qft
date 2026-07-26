# F2 engine running log — BRST sl_5 (3,2), Kazhdan E_1-degeneration

## Stage 0 (done): reconnaissance
- Deliverable files ALREADY EXIST (commit bf71c4e "frontier compute scaffolding"):
  compute/lib/brst_sl5_subregular_engine.py (1929 lines, finite-dim pure-ghost toy, 47 tests pass),
  compute/tests/test_brst_sl5_subregular.py, test_brst_sl5_subregular_iv.py.
- Existing scaffold (Part A) facts: works with ad(h)-eigenvalue grading mislabeled "good even
  Dynkin" (true Dynkin grading ad(h)/2 is HALF-integral); n_+ := g_{>=1/2}, dim 10, LCS depth 3,
  10 brackets, chi inadmissible (3 violations) — that inadmissibility is the KNOWN reason for
  KRW neutral fermions / Lagrangian completion, NOT the Kazhdan E_1 verdict. Docstring sl_2-triple
  e = E12+2E23+E45 (1-indexed) FAILS [e,f]=h (block-3 gives diag(1,1,-2)).
- Manuscript thm:non-principal-ds-bar-cobar-obstruction-32 (thqg_bv_ht_extensions.tex:937-1006)
  mixes conventions: says dim n_+ = 10 / depth 3 AND "2-step with four nonzero commutators".
- Math worked out by hand (to be verified in code):
  h = diag(2,0,-2,1,-1), e = 2E01+2E12+E34, f = E10+E21+E43 (0-indexed): [e,f]=h checked.
  Dynkin (ad h)/2 grading dims: j = 0: 4; ±1/2: 4; ±1: 3; ±3/2: 2; ±2: 1. HALF-INTEGRAL.
  g^f: dim 8 (gl_5: sum of squares of conjugate partition (2,2,1) = 9, minus 1); sl_2-spins
  {2, 3/2, 3/2, 1, 1, 1/2, 1/2, 0} -> W-generator weights {1, 3/2, 3/2, 2, 2, 5/2, 5/2, 3}.
  Symplectic form on g_{1/2} = span{E03,E14,E31,E42}: omega(E03,E31)=1, omega(E14,E31)=-1,
  omega(E14,E42)=1, rest 0; nondegenerate. Lagrangian L = span{E03,E14}.
  m = L + g_{>=1} = {E03,E14,E01,E12,E34,E04,E32,E02}, dim 8: exactly FOUR nonzero commutators:
  [E03,E34]=E04, [E14,E01]=-E04, [E01,E12]=E02, [E03,E32]=E02; 2-step ([m,[m,m]]=0).
  chi = (f|.): chi(E01)=chi(E12)=chi(E34)=1, 0 else; chi([m,m])=0 (character) => Q^2=0.
  => FRONTIER F2 numbers (dim n_+=8, 2-step, 4 commutators, 17 sectors) are CORRECT for the
  KRW-Lagrangian m; scaffold's numbers describe g_{>=1/2} (different object).
  Expected H^0 graded dims (free generation, weights above): w=0:1, 1/2:0, 1:1, 3/2:2, 2:4,
  5/2:6, 3:11.
  Central charge identity: c_B(k) = 24k/(k+5) - 30k - 52 = KW formula; equals (A)-complex sum.
- Plan: Part B appended to engine: B1-B4 Lie layer; B5-B6 mode calculus + Q^2=0;
  B7-B8 cohomology per (weight, ghost) vs expected char; B9 Kazhdan SS E_1 both shears
  (gr d = d_0 = cJ+ccb for deg_K = Delta + charge; gr d = d_chi for sheared deg_K + gh),
  calibrated on principal sl_2; B10 report. k generic rational (avoid -5).

## Stage 1 (in progress): Lie layer B1-B4

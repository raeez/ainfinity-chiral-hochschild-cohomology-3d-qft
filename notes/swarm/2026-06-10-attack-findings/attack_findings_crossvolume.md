# CROSS-VOLUME ANCHOR VERIFICATION — attack journal
Date: 2026-06-10. Axis: cross-volume anchors (concordance.tex, dnp_identification_master.tex, Vol I, Vol III), N-to-weight Gritsenko table, fusion-conjecture scope, {5,5,5} independence, Vol I anchors, thm:G-eq-D-Yplus grade.

## F1. N-to-weight table (5,2,1,1,1/2) at N=(1,2,3,4,6) is a chimera — CONFIRMED
Severity: critical. Verdict: false (table as indexed); citation fabricated.
Anchors: compute/tests/test_kappa_tuple_iv.py:208,214 (c_table={1:10,2:4,3:2,4:2,6:1}, cites "Borcherds 1992 Frame-shape table", "Gritsenko 1999");
chapters/theory/introduction.tex:2344-2356 ("canonical Borcherds family Phi_N ... c_N(0)/2 in (5,2,1,1,1/2)" at "CHL index N in {1,2,3,4,6}");
chapters/connections/concordance.tex:164; chapters/connections/programme_climax_platonic.tex:1517-1531,1577-1585.
Evidence (first-principles sympy, /tmp/jacobi_phi0t2.py): canonical CY-type weak Jacobi forms:
phi_{0,3}=(theta(2z)/theta(z))^2 q^0 = r+2+1/r -> c_3(0,0)=2 -> Exp-lift weight 1 (= Delta_1 on paramodular Gamma_3);
phi_{0,4}=theta(3z)/theta(z) q^0 = r+1+1/r -> c_4(0,0)=1 -> weight 1/2 (= Delta_{1/2} on Gamma_4).
Gritsenko 1999 tower: t=(1,2,3,4) -> c=(10,4,2,1) -> weights (5,2,1,1/2). NO t=6 member; weight-1 form is t=3, weight-1/2 is t=4.
Repo asserts c_4(0)=2 (false: 1) and c_6(0)=1 (no such canonical form). Off-by-one stretch {1,2,3,4} -> {1,2,3,4,6} duplicating weight 1.
Also: actual K3xE/Z_N CHL Borcherds family (Jatkar-Sen, David-Jatkar-Sen, Govindarajan-Krishna): Phi_k(N) weights k=(10,6,4,3,2) at N=(1,2,3,4,6) via Frame shapes 1^24, 1^8 2^8, 1^6 3^6, 1^4 2^2 4^4, 1^2 2^2 3^2 6^2 (k = sum(a_d)/2 - 2); BKM square roots Delta_{k/2} weights (5,3,2,3/2,1). NEITHER family gives (5,2,1,1,1/2).
Citation fraud: "Borcherds Inventiones 109 (1992)" = Monstrous moonshine paper; contains NO paramodular c_N(0) Frame-shape table. Frame-shape weights at N={1,2,3,4,6} would give (5,3,2,3/2,1) not (5,2,1,1,1/2).
Internal inconsistency: conjecture scope says "N in {1,...,6}" (includes N=5) and "eight diagonal orbifolds N in {1,...,8}"; the table has five entries {1,2,3,4,6}.
Heal: replace table by the true Gritsenko tower indexed by paramodular polarization t in {1,2,3,4}: c_t(0)=(10,4,2,1), weights (5,2,1,1/2); OR if indexing K3xE/Z_N CHL orbifolds, use weights (5,3,2,3/2,1) at N=(1,2,3,4,6). Remove the fabricated "Borcherds 1992 Frame-shape table" citation. Note polarization t (paramodular Gamma_t) vs orbifold order N (Gamma_0(N)-level) are DIFFERENT parameters and cannot be identified.

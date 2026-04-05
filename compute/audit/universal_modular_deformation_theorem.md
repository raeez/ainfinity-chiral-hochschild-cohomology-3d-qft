# The Universal Modular Deformation Theorem

## The Problem

The slogan "every object is a projection of Theta^oc" is literally false: the curve X, Ran space, FM compactification, and Arnold relations are not projections of any MC element. The instinct is to downgrade to "every algebraic invariant of the genus expansion is a projection." But the stronger true theorem is available, and it is already proved in the manuscript -- it simply has not been stated in its correct categorical form.

## The Theorem

### Notation

- ChAlg_X^mk: the category of modular Koszul chiral algebras on a smooth projective curve X. Morphisms are maps compatible with augmentations, invariant forms, and genus towers (def:modular-koszul-chiral conditions (D1)-(D6)).
- dgLie_MC^pro: the category of pronilpotent completed dg Lie algebras equipped with a distinguished MC element. Morphisms (g, theta) -> (g', theta') are dg Lie maps f with f(theta) = theta'.
- G = {C_*(M-bar_{g,n})}: the modular cooperad, the chain-level geometric coefficient system encoding X, Ran(X), FM compactification, boundary stratification, and Arnold relations.

### Statement

**Theorem** (Universal modular deformation functor).

*(I) The functor.* The assignment

    Theta: ChAlg_X^mk --> dgLie_MC^pro,    A |--> (g_A^mod, Theta_A)

where g_A^mod = prod_{g,n} Hom_{Sigma_n}(C_*(M-bar_{g,n}), End_A(n)) is the modular convolution dg Lie algebra (def:modular-convolution-dg-lie) and Theta_A := D_A - d_0 is the bar-intrinsic MC element (thm:mc2-bar-intrinsic), is a functor. Concretely: for every morphism f: A -> B in ChAlg_X^mk, the induced map f_*: g_A^mod -> g_B^mod is a morphism of dg Lie algebras satisfying f_*(Theta_A) = Theta_B (prop:mc2-functoriality).

*(II) Arena separation.* The dg Lie algebra g_A^mod is the convolution internal Hom:

    g_A^mod = Hom_{Mod-op}(G, End_A)

between the geometric cooperad G (fixed) and the algebraic endomorphism operad End_A (variable). The geometric data -- the curve X, the Ran space, the FM compactification, the Arnold relations, the clutching maps, the planted-forest strata -- enters exclusively through G. The MC element Theta_A is a section of this Hom: it assigns to each geometric stratum of moduli space an algebraic operation on A. In particular:

- The five-component differential D = d_int + [tau,-] + d_sew + d_pf + hbar.Delta (eq:D-five-components) has each component determined by a structural map of G.
- D^2 = 0 (thm:convolution-d-squared-zero) is the codimension-2 face cancellation in M-bar_{g,n}.
- The bracket [-,-] is the convolution through separating clutching maps of G.

The geometric arena is not extracted from Theta_A; it is the base ring of the convolution algebra in which Theta_A lives.

*(III) Universality.* The functor Theta is universal in the following sense. The MC moduli MC_bullet(g_A^mod) is naturally isomorphic to the modular twisting morphism space:

    MC_bullet(g_A^mod) = Tw_alpha^mod

(cor:theta-twisting-morphism, thm:convolution-master-identification). The element Theta_A is the universal twisting morphism. Every natural derived invariant of the formal genus expansion is a natural transformation from Theta to a target functor on ChAlg_X^mk. The five main theorems are instances:

(A) Adjunction: the genus-0 restriction of Theta_A is the twisting morphism alpha_A inducing B |- Omega (thm:bar-cobar-adjunction).

(B) Inversion: bar-cobar inversion is the acyclicity of the genus-0 component of Theta_A, valid on the Koszul locus (thm:bar-cobar-inversion-qi). The Koszulness hypothesis is built into the domain category ChAlg_X^mk.

(C) Complementarity: the projection pi_{g,2} of the MC equation, composed with Verdier intertwining D(Theta_A) = Theta_{A^!}, yields the Lagrangian decomposition Q_g(A) + Q_g(A^!) = H*(M-bar_g, Z(A)) (thm:quantum-complementarity-main).

(D) Modular characteristic: kappa(A) = tr . pi_sc(Theta_A) is a natural transformation to the additive group (thm:modular-characteristic).

(E) Hochschild: ChirHoch*(A) is the Theta-twisted cohomology of the genus-0 convolution algebra (thm:hochschild-polynomial-growth).

And beyond the five:

(F) Shadow obstruction tower: pi_{r,bullet}(Theta_A) for r = 2, 3, 4, ... is the shadow Postnikov tower, each level a natural transformation (constr:shadow-extraction-explicit).

(G) Modular Koszul datum: the six-fold datum Pi_X(L) is canonically recovered from Theta_A (thm:platonic-recovery).

*(IV) Structural properties.* The functor Theta satisfies:

(a) Additivity: Theta_{A_1 otimes A_2} = Theta_{A_1} + Theta_{A_2} for independent tensor products (prop:independent-sum-factorization).

(b) Verdier intertwining: D_Ran(Theta_A) = Theta_{A^!} (thm:mc2-bar-intrinsic(iv)).

(c) Homotopy invariance: a quasi-isomorphism f: A -> A' induces a quasi-isomorphism f_*: g_A^mod -> g_{A'}^mod and a gauge equivalence Theta_A ~ Theta_{A'} in MC_bullet (thm:shadow-homotopy-invariance).

(d) Weight filtration: f_* preserves the weight filtration F^N, and the shadow algebra A^sh = H_bullet(Def_cyc^mod(A)) is a functor to graded commutative rings (prop:mc2-functoriality(ii)-(iii)).

### Proof

Every claim in (I)-(IV) is proved in the manuscript. The references are:

- (I): prop:mc2-functoriality (functoriality of A |-> Theta_A), thm:mc2-bar-intrinsic (construction of Theta_A), def:modular-convolution-dg-lie (construction of g_A^mod).
- (II): eq:modular-convolution (convolution Hom), eq:D-five-components (five-component differential), thm:convolution-d-squared-zero (D^2 = 0 from face cancellation), const:explicit-convolution-bracket (bracket from clutching).
- (III): thm:five-from-theta (five theorems as projections), thm:convolution-master-identification + cor:theta-twisting-morphism (MC = Tw representability), thm:platonic-recovery (modular Koszul datum recovery).
- (IV)(a): prop:independent-sum-factorization. (b): thm:mc2-bar-intrinsic(iv). (c): thm:shadow-homotopy-invariance. (d): prop:mc2-functoriality(ii)-(iii).

The categorical packaging into a single functor statement is new; all constituent results are proved.

## Why This is the Correct Statement

### Why stronger than "every algebraic invariant is a projection"

The downgraded version is a statement about a single MC element Theta_A for a fixed algebra A. It says: given Theta_A, you can extract invariants by projection. This is true but does not constrain the STRUCTURE of the extraction -- any ad hoc function of Theta_A would qualify.

The theorem above says: the assignment A |-> (g_A^mod, Theta_A) is a FUNCTOR with four structural properties (additivity, Verdier intertwining, homotopy invariance, weight filtration), and ALL derived invariants of the formal genus expansion are NATURAL TRANSFORMATIONS from this functor. Naturality is a falsifiable constraint: it says the extraction maps commute with morphisms of chiral algebras. Any proposed invariant that fails to be natural would falsify the theorem.

This is the difference between "pi appears in the area formula of every circle" and "area is a natural transformation from smooth manifolds to real-valued measures, and pi is the value at S^2." The second constrains what CAN be an invariant.

### Why stronger than the original "every object is a projection"

The original slogan is false because it conflates the arena with the invariant. The theorem cleanly separates them:

- The geometric ARENA (X, Ran, FM, Arnold) enters through the cooperad G, the SOURCE of the convolution Hom. It is the base ring, the coefficient system, the fixed infrastructure.
- The algebraic INVARIANT (kappa, shadows, partition function, R-matrix) is a section Theta_A of this Hom, a natural transformation evaluated at A.

The relationship is not "the arena is a projection of Theta" but "Theta is a section of the Hom over the arena." The arena DEFINES the space in which Theta lives; Theta is the universal element of that space.

This is Yoneda: the MC element Theta_A is the universal element of the representable functor MC_bullet(g_A^mod) = Tw_alpha^mod. Everything is extracted by evaluating natural transformations at the universal element. The geometric arena is the Yoneda embedding's target category; the MC element is the representing object's distinguished point.

### The scope boundary (honest)

The theorem controls: all derived invariants of the formal algebraic genus expansion that can be expressed as evaluations of the modular deformation complex on chains of M-bar_{g,n}.

The theorem does not control:
- Non-perturbative invariants (BTZ entropy, phase transitions, sum over geometries)
- Analytic invariants requiring the sewing envelope A^sew
- Invariants depending on real structure, unitarity, or Hilbert space data
- The geometric arena itself (which is input, not output)
- The open/closed extension (Theta^oc) at genus >= 2, where the clutching chain-level verification is incomplete beyond the annulus

## The Correct Slogan

**False:** "Every object is a projection of Theta^oc."

**Downgraded (true but weak):** "Every algebraic invariant of the genus expansion is a projection of Theta_A."

**Correct (true and strong):** "The modular deformation functor Theta: A |-> (g_A^mod, Theta_A) is the universal derived invariant of the formal genus expansion: every natural invariant factors through it. The geometric arena is the coefficient cooperad that defines the convolution algebra; the MC element is the universal section."

Or in one sentence for a treatise:

**"Chiral Koszul duality on a curve X is the study of a representable modular deformation functor whose universal element is the bar-intrinsic MC class, whose coefficient cooperad is the chain-level moduli of stable curves, and whose natural transformations are the derived invariants of the genus expansion."**

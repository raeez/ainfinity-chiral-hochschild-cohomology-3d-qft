# B1: local densities and integrated domains

The 003 review is preserved in full in `review-003.md`. It isolates the undefined integrated pairing and action on unrestricted smooth fields. The chapter remains a theory on that unrestricted smooth complex; the repair specifies where integration is defined.

For an ordinary homogeneous p-form field, define the local top-form density

`omega_loc(u,v)=(-1)^(p+1)(alpha wedge beta_prime + beta wedge alpha_prime)_[3]`.

The integrated pairing is its normalized integral if at least one input is compactly supported. The product density and all total-derivative primitives used in the Green identity are then compactly supported, because multiplication and differential operators preserve support in the required argument. Smooth coefficients are bounded on the compact set supporting each integral. Stokes therefore has no term at infinity. The boundary expression vanishes if one input has zero boundary pullback or both satisfy the stated chiral condition.

The action starts with the local density `L_k=b wedge dprime(a)+(k/2)a wedge partial(a)`. The full graded BV density is `S_loc=omega_loc(Phi,Q_k Phi)/2`. It is defined for unrestricted smooth fields; its integral is defined on compactly supported fields. For ghost-number-zero fields,

`S_loc = L_k - dprime(a wedge b)/2`.

For compact fields satisfying the boundary condition, the total derivative integrates to zero. The source explicitly defines `E_(k,L,c)`, keeps the unrestricted field complex for local equations and cohomology, and evaluates the variational Hamiltonian identity against compactly supported variations. The displayed first variation is finite even at an unrestricted smooth background because every term contains a compact variation or one of its derivatives. A compact gauge parameter preserves the compact field domain; the local gauge formula remains valid for general smooth gauge parameters satisfying the boundary condition. The master equation is an integrated functional identity on `E_(k,L,c)` and a local variational identity otherwise.

The nonzero field change and its inverse are pointwise linear maps with constant coefficients. They preserve compact support. The comparison now first identifies local pairing and Hamiltonian densities, then integrates on the precise domains. Thus no field complex, local cohomology statement, or compact observable theorem is silently replaced.

## Deciding configurations

Use the chapter's orientation, which is the negative of `dx wedge dy wedge dt`, and normalized integration `1/(2 pi i)`.

1. For nonzero real chi supported in `(1,2)`, set `u=(chi dt,0)` and `v=(0,chi dz wedge dbarz)`. Over the disk `|z|<=R`, their normalized pairing is `R^2 integral chi(t)^2 dt`. This diverges as R tends to infinity. Both fields vanish near the actual boundary. This refutes the missing implication from support away from the boundary to integrability.
2. For `a=t dbarz`, `b=dz`, the chiral boundary condition holds. The action density is `dz wedge dt wedge dbarz = 2i dx wedge dy wedge dt`. Its normalized integral over `|z|<=R`, `0<=t<=T`, is `-R^2 T`. The local field configuration and equations exist, but no unrestricted integrated action is assigned to it.

These are direct symbolic calculations, not numerical approximations. The source includes both configurations. No finite test is represented as a proof of compact-support analytic facts.

## Unchanged mathematics and residual targets

`unchanged-mathematics.json` records the exact unchanged suffix beginning with the zero-level Green operator. It includes the compact quantum collar comparison for every complex level, the full chiral extraction at zero, the nonzero scalar/mode normalization, and all auxiliary interval operations. Their existing finite checks remain in 003 and were not rerun as new evidence for this domain repair.

The full nonzero chiral extraction, an actual boundary Green kernel, intrinsic bulk collision-chain operations and interacting quantum higher operations remain outside the established claim, exactly as stated in003. This repair supplies no whole-book or integration acceptance.

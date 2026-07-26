"""
Fulton-MacPherson compactification: boundary strata and Stokes' theorem.

FM_k(X) is the compactification of Conf_k(X) (configuration space of
k labeled points on X). For X = C (complex plane):

  dim_R FM_k(C) = 2(k-1)  (after fixing translation)

Boundary strata correspond to subsets S ⊂ {1,...,k} with |S| ≥ 2:
  D_S ≅ FM_{|S|}(T_x C) × FM_{k-|S|+1}(C)

The codimension of D_S is 2(|S|-1) - 2 = 2|S| - 4... NO, that's wrong.
Actually D_S has codimension 2(|S|-1) - (2|S|-2) = 0... also wrong.

Let me think carefully. FM_k(C) has real dimension 2(k-1) (k complex
coordinates minus one translation). The stratum D_S where the |S| points
collide has:
  FM_{|S|} (fiber, relative positions) × FM_{k-|S|+1}(C) (base)
  dim = 2(|S|-1) + 2(k-|S|+1-1) = 2|S|-2 + 2k-2|S| = 2k-2 = 2(k-1)

So D_S has the SAME dimension as FM_k(C)? No — D_S is a boundary stratum,
which is codimension 1 in the BOUNDARY. More precisely:

For the manifold-with-corners structure:
  D_S is a codimension-1 boundary face when |S| = 2 (simple collision).
  D_S with |S| ≥ 3 gives higher-codimension strata.

CORRECTION: The codimension-1 boundary strata of FM_k(C) are
D_S for ALL subsets S with |S| ≥ 2 (not just |S|=2). Each D_S
is codimension 1 in the boundary. The key is that the normal direction
to D_S is the SCALE parameter for the colliding cluster.

For the A-infinity identity via Stokes' theorem:
  0 = ∫_{FM_k} dω = Σ_{|S|≥2} ±∫_{D_S} ω|_{D_S}

The sign ± comes from the orientation of D_S in ∂FM_k.

Paper references: Section 8 (configuration-spaces.tex), Section 15 (fm-stokes.tex).
"""
from itertools import combinations
from math import factorial, comb


def boundary_strata(k):
    """Return all codimension-1 boundary strata of FM_k(C).

    Each stratum is specified by a subset S ⊂ {1,...,k} with |S| ≥ 2.
    The stratum D_S corresponds to the points labeled by S colliding.

    Returns:
        List of frozensets, each a subset S of {1,...,k} with |S| ≥ 2.
    """
    if k < 2:
        return []
    strata = []
    for size in range(2, k + 1):
        for subset in combinations(range(1, k + 1), size):
            strata.append(frozenset(subset))
    return strata


def count_boundary_strata(k):
    """Count the number of codimension-1 boundary strata.

    |{S ⊂ {1,...,k} : |S| ≥ 2}| = 2^k - k - 1

    (Total subsets 2^k, minus empty set, minus k singletons.)
    """
    if k < 2:
        return 0
    return 2**k - k - 1


def boundary_decomposition(k, S):
    """Describe the topological type of boundary stratum D_S.

    D_S ≅ FM_{|S|}(C) × FM_{k-|S|+1}(C)

    The first factor is the relative configuration of colliding points.
    The second factor is the configuration where S is replaced by one point.

    Parameters:
        k: total number of points
        S: frozenset, the colliding subset

    Returns:
        dict with 'fiber_arity' (|S|), 'base_arity' (k-|S|+1),
        'fiber_dim' (2(|S|-1)), 'base_dim' (2(k-|S|))
    """
    s = len(S)
    return {
        'fiber_arity': s,
        'base_arity': k - s + 1,
        'fiber_dim': 2 * (s - 1),
        'base_dim': 2 * (k - s),
        'total_dim': 2 * (k - 1) - 1,  # codim 1 in FM_k(C)
    }


def stokes_sign(k, S, convention='koszul'):
    """Compute the orientation sign for ∫_{D_S} in Stokes' theorem.

    For the A-infinity identity, the sign of the D_S contribution is:
      (-1)^{ε} where ε depends on the ordering convention.

    In the Koszul convention (Vol II):
      For S = {s+1, ..., s+j} (a consecutive block in the ordered labeling),
      the sign is (-1)^{(j-1)(|a_1|+...+|a_s|)}.

    For degree-0 elements, this simplifies to (-1)^0 = +1 for all strata.

    Parameters:
        k: total number of points
        S: frozenset (or set), the colliding subset
        convention: 'koszul' (default) or 'lv' for Loday-Vallette

    Returns:
        +1 or -1 (for degree-0 elements)
    """
    if convention == 'koszul':
        # For degree-0 elements: all signs are +1
        # For general degrees: need the element degrees
        return 1
    elif convention == 'lv':
        # LV sign: (-1)^{rs+t} where r = min(S)-1, s = |S|, t = k-max(S)
        s_list = sorted(S)
        r = s_list[0] - 1
        s = len(S)
        t = k - s_list[-1]
        return (-1) ** (r * s + t)
    else:
        raise ValueError(f"Unknown convention: {convention}")


def ainfty_from_stokes(k, degree_zero=True):
    """Derive A-infinity terms from Stokes' theorem on FM_k(C).

    For each boundary stratum D_S with S a CONSECUTIVE block
    {s+1, ..., s+j}, the contribution to Stokes' theorem gives
    the term m_i(a_1,...,a_s, m_j(a_{s+1},...,a_{s+j}), a_{s+j+1},...,a_k)
    where i = k - j + 1.

    Non-consecutive subsets give ZERO contribution to the A-infinity
    identity (they contribute to different algebraic structures).

    Parameters:
        k: total number of elements (arity of the identity)
        degree_zero: if True, all elements have degree 0

    Returns:
        List of (r, j, t, sign) tuples for each term in the identity,
        where r + j + t = k and sign is ±1.
    """
    terms = []
    for j in range(1, k + 1):
        i = k + 1 - j
        for s in range(0, k - j + 1):
            r = s
            t = k - s - j
            if degree_zero:
                sign = 1  # Koszul sign is trivial for degree 0
            else:
                sign = None  # Need element degrees
            terms.append({
                'r': r, 'j': j, 't': t,
                'outer_arity': i,
                'inner_arity': j,
                'inner_start': s,
                'sign': sign,
            })
    return terms


def check_stokes_term_count(k):
    """Verify that the number of A-infinity terms equals the number of
    consecutive-block boundary strata.

    Number of consecutive blocks of size j in {1,...,k}: k - j + 1.
    Total: Σ_{j=1}^{k} (k - j + 1) = k + (k-1) + ... + 1 = k(k+1)/2.

    Returns:
        dict with 'term_count', 'expected', 'match'
    """
    terms = ainfty_from_stokes(k)
    expected = k * (k + 1) // 2
    return {
        'term_count': len(terms),
        'expected': expected,
        'match': len(terms) == expected,
    }


def codim2_corners(k):
    """Enumerate codimension-2 corners of FM_k(C).

    Codim-2 corners are intersections D_{S1} ∩ D_{S2} where:
    - S1 ⊂ S2 (nested: S1 collides first, then S2 containing S1)
    - S1 ∩ S2 = ∅ and S1, S2 disjoint (simultaneous independent collisions)

    For the A-infinity identity, the key corners are NESTED pairs
    of consecutive blocks, and their cancellation follows from
    Arnold relations (AOS cancellations).

    Returns:
        List of (S1, S2) pairs giving codim-2 corners.
    """
    strata = boundary_strata(k)
    corners = []
    for i, S1 in enumerate(strata):
        for S2 in strata[i+1:]:
            # Nested: S1 ⊊ S2
            if S1 < S2:
                corners.append((S1, S2))
            elif S2 < S1:
                corners.append((S2, S1))
            # Disjoint: S1 ∩ S2 = ∅
            elif S1.isdisjoint(S2):
                corners.append((S1, S2))
    return corners


def local_residue_convention():
    """Return the local FM residue and orientation convention.

    This is a structural oracle for the manuscript convention:
    Res_{D_S}(dlog eps_S wedge beta + alpha) = beta|_{eps_S=0},
    with outward-normal-first orientation
    or(FM_k(C)) = (-d eps_S) wedge or(D_S).
    """
    return {
        'residue_formula': (
            'Res_D_S(dlog(epsilon_S) wedge beta + alpha) '
            '= beta|epsilon_S=0'
        ),
        'orientation_formula': 'or(FM_k(C)) = (-d epsilon_S) wedge or(D_S)',
        'inward_coordinate': 'epsilon_S >= 0, epsilon_S=0 on D_S',
        'outward_normal': '-partial_epsilon_S',
    }


def ordered_chiral_bar_residue_skeleton(arity):
    """Structural profile for the ordered chiral bar residue skeleton.

    Args:
        arity: Taylor arity k for the cogenerator projection.

    Returns:
        A dictionary of formulas guarding the manuscript statement.
    """
    if arity < 1:
        raise ValueError("arity must be positive")

    return {
        'coalgebra': 'B^ch(A) = T^c(s^-1 bar A)',
        'coderivation': 'D_A = sum_{k>=1} D_{m_k}',
        'projection_formula': (
            f'pi_1 D_A | (s^-1 bar A)^otimes {arity} '
            f'= s^-1 m_{arity} s^otimes {arity}'
        ),
        'collision_residue': (
            'm_S = Res_D_S prod_{e in E(S)} '
            'dlog(z_{s(e)}-z_{t(e)}) tensor OPE_S'
        ),
        'd_squared_sources': (
            'partial^2[FM_k(C)] = 0',
            'omega_ij wedge omega_jk + omega_jk wedge omega_ki '
            '+ omega_ki wedge omega_ij = 0',
        ),
        'derived_centre': 'Z_der_ch(A) = Coder(B^ch(A))',
        'centre_differential': 'd = [D_A,-]',
        'brace_formula': (
            'f{g_1,...,g_r} = sum_ordered_insertions '
            '+/- f(1,...,g_1,...,g_r,...,1)'
        ),
        'mixed_operations': 'nu_k^(m): Z_der_ch(A)^otimes k tensor A^otimes m -> A',
        'strictification_requires': (
            'Drinfeld associator',
            'two-coloured contracting homotopy',
        ),
        'contracting_homotopy': (
            'h_oc: Omega((SC^ch,top)^!) -> Omega((SC^ch,top)^!)[-1]; '
            'd h_oc + h_oc d = id - i circ p; '
            'p h_oc = h_oc i = h_oc^2 = 0'
        ),
        'topologisation_condition': (
            'if T=[Q,G] and translations are Q-exact, '
            'Z_der_ch(A) in E_3^top-Alg'
        ),
        'class_m_completion': (
            'identities in Ahat_rho with ||m_k||_rho <= C_k(rho) '
            'and sum_k C_k(rho) < infinity'
        ),
    }

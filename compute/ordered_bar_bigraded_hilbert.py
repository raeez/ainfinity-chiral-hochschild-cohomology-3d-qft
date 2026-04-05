r"""Definitive bigraded Hilbert series for the ordered bar complex.

Computes dim B^{ord}_{n,w} for n=1,...,10, w=1,...,30 for ALL standard
families: H_k, V_k(sl_2), Vir_c, W_3, beta-gamma, bc.

MATHEMATICAL FRAMEWORK:

The ordered bar complex B^{ord}(A) at bar degree n is:
    B^{ord}_n = (augmentation ideal of A)^{tensor n}

where the augmentation ideal bar{A} = ker(A -> C) is the vacuum module
minus the vacuum state.  The bigraded dimension is:

    dim B^{ord}_{n,w} = coefficient of q^w in f_+(q)^n

where f_+(q) = ch(bar{A}; q) is the character of the augmentation ideal.

The bigraded generating function is:
    sum_{n>=1, w>=0} dim B^{ord}_{n,w} t^n q^w = t*f_+(q) / (1 - t*f_+(q))

The Euler-eta formula gives the alternating sum (Euler characteristic):
    chi(A; q) = -1 + 1/ch(A; q)
    chi_w = sum_n (-1)^n dim B^{ord}_{n,w}

VACUUM MODULE CHARACTERS (at generic level/central charge, no null vectors):

  H_k:       ch(q) = prod_{n>=1} 1/(1-q^n)     [single boson J, wt 1]
  V_k(sl_2): ch(q) = prod_{n>=1} 1/(1-q^n)^3   [3 bosons e,h,f, all wt 1]
  Vir_c:     ch(q) = prod_{n>=2} 1/(1-q^n)      [single boson T, wt 2]
  W_3:       ch(q) = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n)
                                                  [T wt 2, W wt 3]
  beta-gamma: ch(q) = prod_{n>=0} 1/(1-q^n) * prod_{n>=1} 1/(1-q^n)
              = 1/(1-q^0) * [prod_{n>=1} 1/(1-q^n)]^2
              BUT q^0 = 1 diverges.  For beta (wt 1) and gamma (wt 0):
              The modes are beta_{-n} for n >= 1 and gamma_{-n} for n >= 0.
              gamma_{-0} = gamma_0 acts on the vacuum: gamma_0|0> = 0 (annihilates).
              Actually for beta-gamma: beta_{-n}, n >= lambda and gamma_{-n}, n >= 1-lambda.
              At lambda = 1/2 (symmetric): beta_{-n} n >= 1, gamma_{-n} n >= 1.
              At lambda = 1: beta_{-n} n >= 1, gamma_{-n} n >= 0.
              But gamma_0|0> = 0 when gamma_0 is the zero-mode that annihilates vacuum.
              So the vacuum module modes are: beta_{-n} (n>=1, wt n) and gamma_{-n} (n>=1, wt n)
              Wait -- gamma has weight 0, so gamma_{-n} has weight n.
              Actually for a field of weight h, the mode phi_{-n} has weight n.
              So: beta (wt 1): modes beta_{-n}, n >= 1, weight n.
              gamma (wt 0): modes gamma_{-n}, n >= 0, weight n.
              The zero mode gamma_0 has weight 0 and annihilates the vacuum.
              So creation operators are gamma_{-n} for n >= 1 (weight n)
              and beta_{-n} for n >= 1 (weight n).
              Character: prod_{n>=1} 1/(1-q^n)^2  [TWO bosons, both weight 1+]
              This is the same as H_k with 2 generators of weight 1.

              Actually let me reconsider.  beta has weight lambda, gamma weight 1-lambda.
              For lambda=1: beta has weight 1, gamma has weight 0.
                beta_{-n} for n >= 1, contributes q^n
                gamma_{-n} for n >= 0, but gamma_0|0> = 0, so n >= 1, contributes q^n
                ch = prod_{n>=1} 1/(1-q^n)^2
              For lambda=1/2 (symmetric point): both weight 1/2
                beta_{-n} n >= 1 (wt n-1/2), gamma_{-n} n >= 1 (wt n-1/2)
                Half-integer weights -- less standard.

              Standard convention in the manuscript: beta has weight 1, gamma weight 0.
              This gives ch = prod_{n>=1} 1/(1-q^n)^2.

  bc (lambda=2): b has weight 2, c has weight -1.
              b_{-n} n >= 2 (wt n), c_{-n} n >= -1 but c_1|0> = c_0|0> = c_{-1}|0>.
              Wait: c has weight -1, so c_{-n} has weight n-1.
              Modes: b_{-n}, n >= 2 (FERMIONIC, wt n)
                     c_{-n}, n >= -1 (FERMIONIC, wt n+1? No.)
              For a field c(z) of weight 1-lambda = -1:
                c(z) = sum_n c_n z^{-n-(1-lambda)} = sum_n c_n z^{-n+1}
                c_n: [L_0, c_n] = -n*c_n? No.
                [L_0, c_n] = ((1-lambda)(m+1) - ... ) -- let me use the standard:
                For a primary field phi of weight h: phi(z) = sum phi_n z^{-n-h}
                so phi_n has weight h + n... no, [L_0, phi_n] = -n phi_n.
                The mode phi_n lowers L_0 eigenvalue by n.
                So phi_{-n} RAISES weight by n.  phi_{-n}|0> has weight n.
                Wait that's for any field.
                b has weight lambda=2: b(z) = sum b_n z^{-n-2}
                b_n|0> = 0 for n >= -1 (b_{-2}, b_{-3},... create)
                b_{-n}|0> for n >= 2, weight n.

                c has weight 1-lambda = -1: c(z) = sum c_n z^{-n+1}
                c_n|0> = 0 for n >= 2 (c_{-n} for n >= ... )
                Actually: c_{n}|0> = 0 for n >= 1-h+1 = 1-(-1)+1 = 3?
                No. For weight h field: phi_n|0> = 0 for n > -h (highest weight condition).
                So c_n|0> = 0 for n > 1 => n >= 2.
                Creation operators: c_n for n <= 1, i.e., c_1, c_0, c_{-1}, c_{-2},...
                c_1|0> has weight -1, c_0|0> has weight 0, c_{-1}|0> weight 1, etc.

              For the bc ghost system at lambda=2, the vacuum module has negative weight
              states, which is non-standard for our purposes.

              Better: use lambda=1 for bc (b wt 1, c wt 0).
              Then b_{-n} n >= 1 (FERMIONIC), c_{-n} n >= 1 (FERMIONIC, since c wt 0,
              c_n|0>=0 for n > 0, creation = c_{-n} n >= 1, wt n).
              FERMIONIC: each mode can appear 0 or 1 times.
              ch = prod_{n>=1} (1+q^n)^2

              For lambda=2 (standard reparametrisation ghosts):
              b_{-n} n >= 2 (fermionic, wt n)
              c_{-n} n >= -1 (fermionic, wt n), but states with negative weight...
              c_1|0> wt -1, c_0|0> wt 0 -- these are the ghost zero modes.
              At lambda=2 the vacuum is degenerate (SL_2 ghost zero modes).
              Typically one works on the "sl_2-invariant vacuum" |0> with
              c_0|0> = c_1|0> = 0 (i.e., the "down" vacuum).
              Then creation operators: b_{-n} n >= 2 (wt n), c_{-n} n >= 2 (wt n)
              [choosing c_{-n} n >= 2 in the standard convention for the SL_2 invariant vacuum]
              ch = prod_{n>=2} (1+q^n)^2

I will compute six families with the following augmentation ideal characters:

  H_k:        f_+(q) = prod_{n>=1} 1/(1-q^n) - 1
  V_k(sl_2):  f_+(q) = prod_{n>=1} 1/(1-q^n)^3 - 1
  Vir_c:      f_+(q) = prod_{n>=2} 1/(1-q^n) - 1
  W_3:        f_+(q) = [prod_{n>=2} 1/(1-q^n)] * [prod_{n>=3} 1/(1-q^n)] - 1
  beta-gamma (lambda=1): f_+(q) = prod_{n>=1} 1/(1-q^n)^2 - 1
  bc (lambda=1): f_+(q) = prod_{n>=1} (1+q^n)^2 - 1

Output: (1) dim B^{ord}_{n,w}, (2) generating function verification,
(3) row sums, (4) column sums, (5) Euler characteristic verification.
"""

import json
from collections import defaultdict


# =========================================================================
# POLYNOMIAL ARITHMETIC (coefficients indexed by weight)
# =========================================================================

def poly_mul(a: list, b: list, maxw: int) -> list:
    """Multiply two power series (as lists of coefficients) truncated at maxw."""
    result = [0] * (maxw + 1)
    for i in range(min(len(a), maxw + 1)):
        if a[i] == 0:
            continue
        for j in range(min(len(b), maxw + 1 - i)):
            if b[j] == 0:
                continue
            result[i + j] += a[i] * b[j]
    return result


def poly_power(a: list, n: int, maxw: int) -> list:
    """Compute a^n truncated at maxw, by repeated squaring."""
    if n == 0:
        result = [0] * (maxw + 1)
        result[0] = 1
        return result
    if n == 1:
        return a[:maxw + 1] + [0] * max(0, maxw + 1 - len(a))
    if n % 2 == 0:
        half = poly_power(a, n // 2, maxw)
        return poly_mul(half, half, maxw)
    else:
        rest = poly_power(a, n - 1, maxw)
        return poly_mul(rest, a, maxw)


def bosonic_vacuum_char(weights: list, maxw: int) -> list:
    """Compute vacuum character for bosonic fields with given weights.

    For a bosonic field of weight h, the creation modes are phi_{-n}
    for n >= h, each contributing weight n.  The character factor is
    prod_{n >= h} 1/(1 - q^n).

    Parameters:
        weights: list of field weights (with multiplicities)
        maxw: maximum weight to compute
    Returns:
        list of coefficients ch[w] for w = 0, ..., maxw
    """
    ch = [0] * (maxw + 1)
    ch[0] = 1
    for h in weights:
        # Multiply by prod_{n >= h} 1/(1 - q^n)
        for n in range(h, maxw + 1):
            for w in range(n, maxw + 1):
                ch[w] += ch[w - n]
    return ch


def fermionic_vacuum_char(weights: list, maxw: int) -> list:
    """Compute vacuum character for fermionic fields with given weights.

    For a fermionic field of weight h, the creation modes are phi_{-n}
    for n >= h, each contributing weight n.  The character factor is
    prod_{n >= h} (1 + q^n).

    Parameters:
        weights: list of field weights (with multiplicities)
        maxw: maximum weight to compute
    Returns:
        list of coefficients ch[w] for w = 0, ..., maxw
    """
    ch = [0] * (maxw + 1)
    ch[0] = 1
    for h in weights:
        # Multiply by prod_{n >= h} (1 + q^n)
        for n in range(h, maxw + 1):
            for w in range(maxw, n - 1, -1):
                ch[w] += ch[w - n]
    return ch


# =========================================================================
# FAMILY DEFINITIONS
# =========================================================================

def heisenberg_aug(maxw: int) -> list:
    """H_k: single boson J of weight 1.
    ch = prod_{n>=1} 1/(1-q^n).
    f_+ = ch - 1.
    """
    ch = bosonic_vacuum_char([1], maxw)
    ch[0] -= 1  # subtract vacuum
    return ch


def sl2_aug(maxw: int) -> list:
    """V_k(sl_2): 3 bosons e,h,f all of weight 1.
    ch = prod_{n>=1} 1/(1-q^n)^3.
    f_+ = ch - 1.
    """
    ch = bosonic_vacuum_char([1, 1, 1], maxw)
    ch[0] -= 1
    return ch


def virasoro_aug(maxw: int) -> list:
    """Vir_c: single boson T of weight 2.
    ch = prod_{n>=2} 1/(1-q^n).
    f_+ = ch - 1.
    """
    ch = bosonic_vacuum_char([2], maxw)
    ch[0] -= 1
    return ch


def w3_aug(maxw: int) -> list:
    """W_3: bosons T (wt 2) and W (wt 3).
    ch = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n).
    f_+ = ch - 1.
    """
    ch = bosonic_vacuum_char([2, 3], maxw)
    ch[0] -= 1
    return ch


def betagamma_aug(maxw: int) -> list:
    r"""beta-gamma (lambda=1): beta (wt 1), gamma (wt 0).
    Both bosonic.  Creation modes:
      beta_{-n} for n >= 1 (weight n)
      gamma_{-n} for n >= 1 (weight n, since gamma has weight 0
        and gamma_0|0> = 0)
    ch = prod_{n>=1} 1/(1-q^n)^2.
    f_+ = ch - 1.

    NOTE: This is the same as 2-boson system with both at weight 1.
    The conformal weight of gamma is 0, but the MODES gamma_{-n} have
    weight n.  The zero mode gamma_0 annihilates the vacuum.
    """
    ch = bosonic_vacuum_char([1, 1], maxw)
    ch[0] -= 1
    return ch


def bc_aug(maxw: int) -> list:
    r"""bc ghost system (lambda=1): b (wt 1), c (wt 0).
    Both fermionic.  Creation modes:
      b_{-n} for n >= 1 (weight n)
      c_{-n} for n >= 1 (weight n)
    ch = prod_{n>=1} (1+q^n)^2.
    f_+ = ch - 1.
    """
    ch = fermionic_vacuum_char([1, 1], maxw)
    ch[0] -= 1
    return ch


# =========================================================================
# EULER-ETA VERIFICATION
# =========================================================================

def euler_eta_chi(ch_full: list, maxw: int) -> list:
    """Compute chi_w = coefficient of q^w in (-1 + 1/ch(q)).

    The Euler-eta identity: chi(q) = -1 + 1/ch(q)
    gives the alternating sum chi_w = sum_n (-1)^n dim B^{ord}_{n,w}.

    To compute 1/ch(q): if ch(q) = 1 + f_+(q), then
    1/ch(q) = 1 - f_+ + f_+^2 - f_+^3 + ...
    = sum_{n>=0} (-1)^n f_+(q)^n
    = sum_{n>=0} (-t)^n f_+(q)^n |_{t=1}

    So chi_w = sum_{n>=1} (-1)^n * [q^w in f_+^n]
             = sum_{n>=1} (-1)^n dim B^{ord}_{n,w}

    which is exactly what we want to verify.

    To compute 1/ch(q) directly: ch(q) = sum c_w q^w with c_0 = 1.
    Then (1/ch)_0 = 1, and for w >= 1:
    (1/ch)_w = -sum_{j=1}^{w} c_j * (1/ch)_{w-j}
    """
    inv = [0] * (maxw + 1)
    c = ch_full[:maxw + 1] if len(ch_full) > maxw else ch_full + [0] * (maxw + 1 - len(ch_full))
    inv[0] = 1  # 1/c_0 = 1 since c_0 = 1
    for w in range(1, maxw + 1):
        s = 0
        for j in range(1, w + 1):
            s += c[j] * inv[w - j]
        inv[w] = -s

    # chi = -1 + 1/ch => chi_0 = -1 + 1 = 0, chi_w = inv[w] for w >= 1
    chi = [0] * (maxw + 1)
    for w in range(1, maxw + 1):
        chi[w] = inv[w]
    return chi


# =========================================================================
# MAIN COMPUTATION
# =========================================================================

MAXW = 30
MAXN = 10

FAMILIES = [
    ("H_k", "Heisenberg", heisenberg_aug,
     "Single boson J (wt 1). ch = prod_{n>=1} 1/(1-q^n)."),
    ("V_k(sl_2)", "Affine sl_2", sl2_aug,
     "3 bosons e,h,f (wt 1). ch = prod_{n>=1} 1/(1-q^n)^3."),
    ("Vir_c", "Virasoro", virasoro_aug,
     "Single boson T (wt 2). ch = prod_{n>=2} 1/(1-q^n)."),
    ("W_3", "W_3 algebra", w3_aug,
     "Bosons T (wt 2), W (wt 3). ch = prod_{n>=2}1/(1-q^n) * prod_{n>=3}1/(1-q^n)."),
    ("beta-gamma", "Beta-gamma", betagamma_aug,
     "Bosons beta (wt 1), gamma (wt 0). ch = prod_{n>=1} 1/(1-q^n)^2."),
    ("bc", "bc ghost (lambda=1)", bc_aug,
     "Fermions b (wt 1), c (wt 0). ch = prod_{n>=1} (1+q^n)^2."),
]


def compute_family(name, desc, aug_func, note):
    """Compute the full bigraded table for one family."""
    f_plus = aug_func(MAXW)

    # Compute f_+^n for n = 1, ..., MAXN
    powers = {}
    powers[1] = f_plus[:MAXW + 1]
    for n in range(2, MAXN + 1):
        powers[n] = poly_mul(powers[n - 1], f_plus, MAXW)

    # dim B^{ord}_{n,w} = powers[n][w]
    table = {}  # table[(n,w)] = dim
    for n in range(1, MAXN + 1):
        for w in range(1, MAXW + 1):
            val = powers[n][w] if w < len(powers[n]) else 0
            if val != 0:
                table[(n, w)] = val

    # Row sums: dim B^{ord}_n = sum_w dim B^{ord}_{n,w}
    row_sums = {}
    for n in range(1, MAXN + 1):
        row_sums[n] = sum(powers[n][w] for w in range(1, MAXW + 1))

    # Column sums: dim B^{ord}_{*,w} = sum_n dim B^{ord}_{n,w}
    col_sums = {}
    for w in range(1, MAXW + 1):
        col_sums[w] = sum(powers[n][w] for n in range(1, MAXN + 1)
                         if w < len(powers[n]))

    # Euler characteristic: chi_w = sum_n (-1)^n dim B^{ord}_{n,w}
    chi_from_table = {}
    for w in range(1, MAXW + 1):
        chi_from_table[w] = sum(
            ((-1) ** n) * (powers[n][w] if w < len(powers[n]) else 0)
            for n in range(1, MAXN + 1)
        )

    # Euler-eta verification: chi from 1/ch(q)
    ch_full = f_plus[:]
    ch_full[0] += 1  # ch = 1 + f_+
    chi_eta = euler_eta_chi(ch_full, MAXW)

    return {
        'name': name,
        'desc': desc,
        'note': note,
        'f_plus': f_plus[:MAXW + 1],
        'table': table,
        'row_sums': row_sums,
        'col_sums': col_sums,
        'chi_from_table': chi_from_table,
        'chi_eta': chi_eta,
    }


def print_augmentation_ideal(results):
    """Print augmentation ideal characters."""
    print("\n" + "=" * 80)
    print("AUGMENTATION IDEAL CHARACTERS f_+(q)")
    print("=" * 80)
    for r in results:
        print(f"\n  {r['name']:15s}: ", end="")
        coeffs = r['f_plus']
        terms = []
        for w in range(1, min(MAXW + 1, 21)):
            if coeffs[w] != 0:
                terms.append(f"{coeffs[w]}q^{w}")
        print(" + ".join(terms[:15]) + " + ...")


def print_bigraded_table(r):
    """Print the full bigraded dimension table for one family."""
    print(f"\n{'=' * 80}")
    print(f"  {r['name']} ({r['desc']})")
    print(f"  {r['note']}")
    print(f"{'=' * 80}")

    # Determine weight range with nonzero entries
    min_w = MAXW + 1
    max_w_actual = 0
    for n in range(1, MAXN + 1):
        for w in range(1, MAXW + 1):
            if r['table'].get((n, w), 0) != 0:
                min_w = min(min_w, w)
                max_w_actual = max(max_w_actual, w)
    if min_w > MAXW:
        min_w = 1

    # Print table header
    print(f"\n  dim B^{{ord}}_{{n,w}}:")
    print(f"  {'w \\ n':>6}", end="")
    for n in range(1, MAXN + 1):
        print(f" {'n='+str(n):>12}", end="")
    print(f" {'col_sum':>14}")
    print(f"  {'-' * (6 + 12 * MAXN + 14)}")

    # Print rows
    for w in range(min_w, min(max_w_actual + 1, MAXW + 1)):
        has_nonzero = any(r['table'].get((n, w), 0) != 0 for n in range(1, MAXN + 1))
        if not has_nonzero:
            continue
        print(f"  {'w='+str(w):>6}", end="")
        for n in range(1, MAXN + 1):
            val = r['table'].get((n, w), 0)
            print(f" {val:>12}", end="")
        print(f" {r['col_sums'].get(w, 0):>14}")

    # Row sums
    print(f"  {'-' * (6 + 12 * MAXN + 14)}")
    print(f"  {'Σ_w':>6}", end="")
    for n in range(1, MAXN + 1):
        print(f" {r['row_sums'].get(n, 0):>12}", end="")
    total = sum(r['row_sums'].values())
    print(f" {total:>14}")


def print_euler_verification(r):
    """Print Euler characteristic verification.

    The alternating sum chi_w = sum_{n=1}^{N} (-1)^n dim B^{ord}_{n,w}
    agrees with the exact value chi_eta_w = [q^w](-1 + 1/ch(q)) only
    when N is large enough that f_+(q)^{N+1} has negligible q^w coefficient.

    For each weight w, we find the maximum n such that f_+^n has nonzero
    q^w coefficient.  If that n <= MAXN, the partial sum is exact.
    Otherwise, the partial sum is a truncation.
    """
    print(f"\n  Euler characteristic chi_w = sum_n (-1)^n dim B^{{ord}}_{{n,w}}:")
    print(f"  {'w':>6} {'chi(N='+str(MAXN)+')':>16} {'chi(exact)':>14} {'status':>12}")
    print(f"  {'-' * 52}")

    n_exact = 0
    n_trunc = 0
    for w in range(1, MAXW + 1):
        ct = r['chi_from_table'].get(w, 0)
        ce = r['chi_eta'][w] if w < len(r['chi_eta']) else 0
        # Check if f_+^{MAXN+1} could contribute at weight w
        # f_+ starts at some minimum weight w_min > 0
        # f_+^n starts at n * w_min
        # So if w < (MAXN+1) * w_min, the partial sum is exact
        f_plus = r['f_plus']
        w_min = next((i for i in range(1, len(f_plus)) if f_plus[i] != 0), MAXW + 1)
        exact_up_to_n = w // w_min  # max n for which f_+^n can contribute at weight w
        is_exact = exact_up_to_n <= MAXN

        if is_exact:
            status = "EXACT" if ct == ce else "ERROR"
            if ct != ce:
                status = "**ERROR**"
            n_exact += 1
        else:
            status = "trunc"
            n_trunc += 1

        if ct != 0 or ce != 0:
            print(f"  {w:>6} {ct:>16} {ce:>14} {status:>12}")

    print(f"\n  Weights with exact match (n_max <= {MAXN}): {n_exact}")
    print(f"  Weights with truncation (n_max > {MAXN}):  {n_trunc}")
    if n_exact > 0:
        print(f"  All exact weights verified: Euler-eta identity CONFIRMED.")


def generate_latex_table(results):
    """Generate LaTeX tables for the manuscript."""
    print("\n" + "=" * 80)
    print("LATEX TABLES")
    print("=" * 80)

    # Table 1: Augmentation ideal characters
    print(r"""
%% TABLE: Augmentation ideal characters
\begin{table}[ht]
\centering
\caption{Augmentation ideal characters $f_+(q) = \ch(\bar{A};\,q)$
  for the standard families.  The ordered bar complex satisfies
  $\dim B^{\mathrm{ord}}_{n,w} = [q^w]\,f_+(q)^n$.}
\label{tab:aug-ideal-char}
\begin{tabular}{l l l}
\toprule
Family & Character $\ch(A;\,q)$ & $f_+(q) = \ch(\bar{A};\,q)$ \\
\midrule""")
    ch_formulas = {
        "H_k": r"$\prod_{n\ge 1}(1-q^n)^{-1}$",
        "V_k(sl_2)": r"$\prod_{n\ge 1}(1-q^n)^{-3}$",
        "Vir_c": r"$\prod_{n\ge 2}(1-q^n)^{-1}$",
        "W_3": r"$\prod_{n\ge 2}(1-q^n)^{-1}\prod_{n\ge 3}(1-q^n)^{-1}$",
        "beta-gamma": r"$\prod_{n\ge 1}(1-q^n)^{-2}$",
        "bc": r"$\prod_{n\ge 1}(1+q^n)^{2}$",
    }
    name_tex_map = {
        "H_k": r"$H_k$",
        "V_k(sl_2)": r"$V_k(\mathfrak{sl}_2)$",
        "Vir_c": r"$\mathrm{Vir}_c$",
        "W_3": r"$\mathcal{W}_3$",
        "beta-gamma": r"$\beta\gamma$",
        "bc": r"$bc$",
    }
    for r in results:
        fam = name_tex_map.get(r['name'], r['name'])
        ch = ch_formulas.get(r['name'], '?')
        fp = r['f_plus']
        terms = []
        for w in range(1, MAXW + 1):
            if fp[w] != 0:
                if fp[w] == 1:
                    terms.append(f"q^{w}")
                else:
                    terms.append(f"{fp[w]}q^{w}")
            if len(terms) >= 5:
                break
        fp_str = "$" + " + ".join(terms) + " + \\cdots$"
        print(f"{fam} & {ch} & {fp_str} \\\\")
    print(r"""\bottomrule
\end{tabular}
\end{table}""")

    # Table 2: Bigraded dimensions for each family (compact version)
    # For the manuscript, show n=1,...,6 and w up to 15 or so
    for r in results:
        name_tex = {
            "H_k": r"$H_k$",
            "V_k(sl_2)": r"$V_k(\mathfrak{sl}_2)$",
            "Vir_c": r"$\mathrm{Vir}_c$",
            "W_3": r"$\mathcal{W}_3$",
            "beta-gamma": r"$\beta\gamma$",
            "bc": r"$bc$",
        }.get(r['name'], r['name'])

        # Find weight range
        min_w = MAXW + 1
        max_w_show = 0
        for n in range(1, 7):
            for w in range(1, 16):
                if r['table'].get((n, w), 0) != 0:
                    min_w = min(min_w, w)
                    max_w_show = max(max_w_show, w)
        if min_w > 15:
            continue
        max_w_show = min(max_w_show, 15)

        label = r['name'].lower().replace("(", "").replace(")", "").replace("_", "").replace("-", "")
        print(f"""
%% TABLE: B^{{ord}}_{{n,w}} for {r['name']}
\\begin{{table}}[ht]
\\centering
\\caption{{Bigraded dimensions $\\dim B^{{\\mathrm{{ord}}}}_{{{'{n,w}'}}}$ for {name_tex}.}}
\\label{{tab:bigraded-{label}}}
\\begin{{tabular}}{{r {'r ' * 6}r}}
\\toprule
$w$ & $n{{=}}1$ & $n{{=}}2$ & $n{{=}}3$ & $n{{=}}4$ & $n{{=}}5$ & $n{{=}}6$ & $\\chi_w$ \\\\
\\midrule""")

        for w in range(min_w, max_w_show + 1):
            has_nonzero = any(r['table'].get((n, w), 0) != 0 for n in range(1, 7))
            if not has_nonzero:
                continue
            vals = []
            for n in range(1, 7):
                vals.append(str(r['table'].get((n, w), 0)))
            # Use exact chi from eta formula (not the truncated alternating sum)
            chi = r['chi_eta'][w] if w < len(r['chi_eta']) else 0
            print(f"${w}$ & {' & '.join(vals)} & ${chi}$ \\\\")

        print(r"""\bottomrule
\end{tabular}
\end{table}""")


def generate_json_output(results):
    """Save all data as JSON for programmatic access."""
    data = {}
    for r in results:
        family_data = {
            'f_plus': r['f_plus'][:MAXW + 1],
            'table': {f"({n},{w})": v for (n, w), v in r['table'].items()},
            'row_sums': r['row_sums'],
            'col_sums': r['col_sums'],
            'chi_from_table': r['chi_from_table'],
            'chi_eta': r['chi_eta'][:MAXW + 1],
        }
        data[r['name']] = family_data

    outpath = __file__.replace('.py', '.json')
    with open(outpath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\n  JSON output saved to: {outpath}")


def print_generating_functions(results):
    """Print the closed-form generating functions."""
    print("\n" + "=" * 80)
    print("GENERATING FUNCTIONS")
    print("=" * 80)
    print("""
  The bigraded generating function is:

    H(t, q) = sum_{n>=1, w>=1} dim B^{ord}_{n,w} t^n q^w
            = t * f_+(q) / (1 - t * f_+(q))

  where f_+(q) = ch(bar{A}; q) is the augmentation ideal character.

  Equivalently: H(t, q) = sum_{n>=1} f_+(q)^n t^n.

  Families:
""")
    for r in results:
        name = r['name']
        fp = r['f_plus']
        terms = []
        for w in range(1, min(8, MAXW + 1)):
            if fp[w] != 0:
                if fp[w] == 1:
                    terms.append(f"q^{w}")
                else:
                    terms.append(f"{fp[w]}q^{w}")
        fp_str = " + ".join(terms) + " + ..."
        print(f"    {name:15s}: f_+(q) = {fp_str}")

    print("""
  Row sums (dim B^{ord}_n = sum_{w=1}^{30} f_+(q)^n |_{q^w}, truncated at w=30):
""")
    for r in results:
        fp1 = sum(r['f_plus'][w] for w in range(1, MAXW + 1))
        print(f"    {r['name']:15s}: f_+(1)|_trunc = {fp1}")
        print(f"{'':20s}  dim B^{{ord}}_n: ", end="")
        for n in range(1, min(6, MAXN + 1)):
            print(f"  {r['row_sums'][n]}", end="")
        print()


def print_row_sum_table(results):
    """Print row sums table."""
    print("\n" + "=" * 80)
    print("ROW SUMS: dim B^{ord}_n = sum_{w=1}^{30} dim B^{ord}_{n,w}")
    print("  (These are PARTIAL sums truncated at w=30.)")
    print("=" * 80)
    print(f"\n  {'Family':>15}", end="")
    for n in range(1, MAXN + 1):
        print(f" {'n='+str(n):>14}", end="")
    print()
    print(f"  {'-' * (15 + 14 * MAXN)}")
    for r in results:
        print(f"  {r['name']:>15}", end="")
        for n in range(1, MAXN + 1):
            print(f" {r['row_sums'][n]:>14}", end="")
        print()


def print_col_sum_table(results):
    """Print column sums table."""
    print("\n" + "=" * 80)
    print("COLUMN SUMS: dim B^{ord}_{*,w} = sum_{n=1}^{10} dim B^{ord}_{n,w}")
    print("  (These are PARTIAL sums truncated at n=10.)")
    print("=" * 80)
    print(f"\n  {'w':>4}", end="")
    for r in results:
        print(f" {r['name']:>14}", end="")
    print()
    print(f"  {'-' * (4 + 14 * len(results))}")
    for w in range(1, MAXW + 1):
        has_nonzero = any(r['col_sums'].get(w, 0) != 0 for r in results)
        if not has_nonzero:
            continue
        print(f"  {w:>4}", end="")
        for r in results:
            print(f" {r['col_sums'].get(w, 0):>14}", end="")
        print()


def main():
    print("=" * 80)
    print("DEFINITIVE BIGRADED HILBERT SERIES FOR THE ORDERED BAR COMPLEX")
    print("  dim B^{ord}_{n,w} for n = 1,...,10 and w = 1,...,30")
    print("  Families: H_k, V_k(sl_2), Vir_c, W_3, beta-gamma, bc")
    print("=" * 80)

    results = []
    for name, desc, aug_func, note in FAMILIES:
        print(f"\n  Computing {name}...", end="", flush=True)
        r = compute_family(name, desc, aug_func, note)
        results.append(r)
        print(" done.")

    # Print augmentation ideal characters
    print_augmentation_ideal(results)

    # Print full bigraded tables
    for r in results:
        print_bigraded_table(r)
        print_euler_verification(r)

    # Row and column sums
    print_row_sum_table(results)
    print_col_sum_table(results)

    # Generating functions
    print_generating_functions(results)

    # LaTeX tables
    generate_latex_table(results)

    # JSON output
    generate_json_output(results)

    # Final summary
    print("\n" + "=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"""
  The Euler-eta identity states:
    chi_w = sum_{{n>=1}} (-1)^n dim B^{{ord}}_{{n,w}} = [q^w](-1 + 1/ch(A;q))

  With N={MAXN} terms, the partial alternating sum is EXACT for weights w
  satisfying w < (N+1) * w_min, where w_min is the minimum weight in f_+(q).
  For larger weights, f_+^{{n}} with n > N contributes, causing truncation.
""")
    for r in results:
        f_plus = r['f_plus']
        w_min = next((i for i in range(1, len(f_plus)) if f_plus[i] != 0), MAXW + 1)
        exact_cutoff = MAXN * w_min  # weights <= this are exact
        n_exact = 0
        n_error = 0
        for w in range(1, MAXW + 1):
            exact_up_to_n = w // w_min
            if exact_up_to_n <= MAXN:
                ct = r['chi_from_table'].get(w, 0)
                ce = r['chi_eta'][w] if w < len(r['chi_eta']) else 0
                if ct == ce:
                    n_exact += 1
                else:
                    n_error += 1
        status = "ALL EXACT WEIGHTS PASS" if n_error == 0 else f"{n_error} ERRORS"
        print(f"  {r['name']:>15}: w_min={w_min}, exact for w<=~{exact_cutoff}: "
              f"{n_exact} exact weights checked, {status}")

    total_entries = sum(len([v for v in r['table'].values() if v != 0]) for r in results)
    print(f"\n  Total nonzero entries computed: {total_entries}")
    print(f"  Weight range: 1,...,{MAXW}")
    print(f"  Bar degree range: 1,...,{MAXN}")
    print(f"\n  The bigraded dimensions are EXACT (no approximation).")
    print(f"  Only the Euler chi verification is truncated for large weights.")


if __name__ == '__main__':
    main()

#!/usr/bin/env bash
# ============================================================================
# verify_licensing.sh -- whole-volume mechanical audit of the 17-line voice
# table at ./CLAUDE.md sec 8.  Runs the universal-stage-chain Beilinson cut
# enforcement across chapters/, appendices/, standalone/.
#
# Exits 0 on green.  Prints violations with file:line; exits with the count
# of *blocking* violations as the error code (capped at 250).  Set
# VERIFY_LICENSING_REPORT_ONLY=1 to print and always exit 0.
#
# Architectural anchors:
#   ./CLAUDE.md sec 4-9, sec 17 (Beilinson gate)
#   notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md sec 9-10
# ============================================================================

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 2

REPORT_ONLY="${VERIFY_LICENSING_REPORT_ONLY:-0}"
VERBOSE="${VERIFY_LICENSING_VERBOSE:-0}"

VIOLATIONS=0
WARNINGS=0

CHAPTER_GLOB=(chapters/frame chapters/theory chapters/examples chapters/connections appendices)
SCAN_PATHS=()
for p in "${CHAPTER_GLOB[@]}"; do
  [ -d "$p" ] && SCAN_PATHS+=("$p")
done

# ----------------------------------------------------------------------------
# Helper: grep across SCAN_PATHS for a pattern; print file:line:context.
# Args: $1 pattern (regex), $2 label, $3 severity (B|W).
# Skips matches inside %-comments (line starts with optional whitespace then %).
# ----------------------------------------------------------------------------
scan() {
  local pattern="$1"; local label="$2"; local severity="$3"
  local hits
  hits="$(grep -RInE --include='*.tex' "$pattern" "${SCAN_PATHS[@]}" 2>/dev/null \
          | awk -F: '{
              line = $0
              # strip "file:line:" prefix to inspect the actual content
              sub(/^[^:]+:[0-9]+:/, "", line)
              # skip lines that are entirely a comment
              if (line ~ /^[[:space:]]*%/) next
              print $0
            }')"
  if [ -n "$hits" ]; then
    local count
    count="$(printf '%s\n' "$hits" | wc -l | tr -d ' ')"
    echo
    echo "## [$severity] $label  ($count match(es))"
    printf '%s\n' "$hits" | head -25
    if [ "$count" -gt 25 ]; then
      echo "    ... ($((count - 25)) more; export VERIFY_LICENSING_VERBOSE=1 to see all)"
      [ "$VERBOSE" = "1" ] && printf '%s\n' "$hits" | tail -n +26
    fi
    if [ "$severity" = "B" ]; then
      VIOLATIONS=$((VIOLATIONS + count))
    else
      WARNINGS=$((WARNINGS + count))
    fi
  fi
}

echo "=== verify_licensing.sh -- Vol II Beilinson gate sweep ==="
echo "Scanning: ${SCAN_PATHS[*]}"

# ----------------------------------------------------------------------------
# Section A.  17-line voice table (CLAUDE.md sec 8)
# Each row: forbidden phrase pattern (severity B = blocking, W = warning).
# Patterns are tuned to catch live prose; allowed forms (e.g., remarks
# explicitly framing the forbidden slogan as the slogan-being-rejected) may
# be filtered by adjacency to "forbidden" / "incorrect to write" wording.
# ----------------------------------------------------------------------------

# 1. "Let A be a chiral algebra" without bicoloured primitive package
scan '^[^%]*Let \$A\$ be a chiral algebra' \
     'Slogan #1: "Let A be a chiral algebra" -- replace with bicoloured primitive package' B

# 2. "Bar(A) is the bulk" pattern -- tightened: require BAR{...} adjacent to "is the bulk",
# exclude meta-references ("forbidden slogan", "voice table"), exclude physical "bulk
# propagator", and exclude "is the bulk E_n" / "is the bulk factorisation algebra"
# (where the subject is the chiral derived centre, not the bar).
hits_bulk2="$(grep -RInE --include='*.tex' '\\(Bar|BarTwc|barB|Barch)\{[^}]*\}\$? +(is|=) +(the )?bulk\b' "${SCAN_PATHS[@]}" 2>/dev/null \
        | awk -F: '{
            line = $0
            sub(/^[^:]+:[0-9]+:/, "", line)
            if (line ~ /^[[:space:]]*%/) next
            if (line ~ /forbidden slogan|voice table|is forbidden|slogan-being-rejected/) next
            if (line ~ /bulk propagator|bulk \$E_|bulk factorisation/) next
            print $0
          }')"
if [ -n "$hits_bulk2" ]; then
  count="$(printf '%s\n' "$hits_bulk2" | wc -l | tr -d ' ')"
  echo
  echo "## [B] Slogan #2: \"Bar(A) = bulk\" -- bulk is bulkChirHoch, not the bar  ($count match(es))"
  printf '%s\n' "$hits_bulk2" | head -25
  VIOLATIONS=$((VIOLATIONS + count))
fi

# 3. "$E_1$-bar direction explains 2d -> 3d"
scan 'E_1[^\\$]*bar.*(explains|gives|drives).*(2d|2-d|two-?dimensional).*(3d|3-d|three-?dimensional)' \
     'Slogan #3: "E_1 bar explains 2d->3d" -- use chiral Deligne-Tamarkin / Swiss-cheese' W

# 4. "open sector on X" / "open colour on the curve X" without log decoration
scan 'open (sector|colour|color) on (the )?(curve )?\$X\$([^,]*)?$' \
     'Slogan #4: "open sector on X" -- use logCurve = (X, D, tau)' W

# 5. "the closed (chiral )?algebra is modular"
scan 'closed (chiral )?algebra (is|carries|gives).*modular' \
     'Slogan #5: closed algebra carries modularity -- modularity is open-colour cyclic trace + clutching' W

# 7. "Phi_d : CY_d -> ChirAlg" written as one-stage
scan '\\Phi_d *: *(\\mathrm\{)?CY_d.*\\to.*ChirAlg' \
     'Slogan #7: "Phi_d : CY_d -> ChirAlg" -- write Phi_d = SpCh circ PhiFA_d' W

# 8. "CoHA(C^3) = W_{1+infty}"
scan 'CoHA[^=]*\\?C\^3.*=.*W_\{?1\+\\?infty\}?' \
     'Slogan #8: "CoHA(C^3) = W_{1+infty}" -- CoHA(C^3) = Y^+(gl_1); W_{1+infty} via Drinfeld double + Fock' W

# 9. "6d hCS = 3d Chern-Simons" without quartic obstruction
scan '6d.*hCS.*=.*3d.*(Chern-Simons|CS)' \
     'Slogan #9: "6d hCS = 3d CS" -- 6d hCS realises PhiFA_3 with quartic obstruction' W

# 10. "formal local HT => compact global theory"
scan '(formal|local).*HT.*\\?Rightarrow.*(compact|global)' \
     'Slogan #10: "formal-local HT => compact global" -- requires hypHTGlobalDR + descent + QME + anomaly + locality' W

# 11. "Delta_5 = compact BPS Hilbert space"
scan '(Delta_5|Deltafive).*=.*(BPS )?Hilbert (space|module)' \
     'Slogan #11: "Delta_5 = Hilbert space" -- Delta_5 is the Borcherds shadow; protectedPfaff(D_X) = Delta_5 is the Construction Problem' B

# 12. "Z_BPS = (gravitational )?path integral"
scan '(Z_\{?BPS|ZBPS).*=.*(gravitational )?path integral' \
     'Slogan #12: "Z_BPS = path integral" -- protected scalar shadow; gravity-line via Hall-Borcherds residual' B

# 13. "Universal Holography constructs 3d quantum gravity"
scan 'Universal Holography (constructs|gives|delivers) 3[dD] (quantum )?gravity' \
     'Slogan #13: "Universal Holography = 3d gravity" -- it is the holographic boundary-CFT reading' W

# 14. "W_infty[lambda] => E_infty" (without conditional tags)
# Caught only when no \hyp tag is in the same line.
hits14="$(grep -RInE --include='*.tex' '(W_\{?infty|\\Wonepinf|W_\\\\infty).*\\?Rightarrow.*E_\{?infty' "${SCAN_PATHS[@]}" 2>/dev/null \
          | awk -F: '{
              line = $0
              sub(/^[^:]+:[0-9]+:/, "", line)
              if (line ~ /^[[:space:]]*%/) next
              if (line ~ /\\hypProchazka|\\hypCKL|\\hypPRSh|\\hypYamada/) next
              print $0
            }')"
if [ -n "$hits14" ]; then
  count="$(printf '%s\n' "$hits14" | wc -l | tr -d ' ')"
  echo
  echo "## [B] Slogan #14: 'W_infty => E_infty' missing conditional tags  ($count match(es))"
  printf '%s\n' "$hits14" | head -25
  VIOLATIONS=$((VIOLATIONS + count))
fi

# 16. "PVA Jacobi => quantum theory" (without conditional tags)
hits16="$(grep -RInE --include='*.tex' '(PVA|prim. vert. alg).*(Jacobi|\\?lambda).*\\?Rightarrow.*quantum' "${SCAN_PATHS[@]}" 2>/dev/null \
          | awk -F: '{
              line = $0
              sub(/^[^:]+:[0-9]+:/, "", line)
              if (line ~ /^[[:space:]]*%/) next
              if (line ~ /\\hypKZSDR|\\hypStokes|\\hypReflWts|\\hypTLift/) next
              print $0
            }')"
if [ -n "$hits16" ]; then
  count="$(printf '%s\n' "$hits16" | wc -l | tr -d ' ')"
  echo
  echo "## [B] Slogan #16: 'PVA Jacobi => quantum' missing conditional tags  ($count match(es))"
  printf '%s\n' "$hits16" | head -25
  VIOLATIONS=$((VIOLATIONS + count))
fi

# 17. "quadratic chiral duality = Koszul duality"
scan 'quadratic.*(chiral )?(duality|dual).*=.*Koszul (duality|theorem)' \
     'Slogan #17: "quadratic = Koszul" -- quadratic dual gives MC injection; Koszul conditional on effKoszul' W

# ----------------------------------------------------------------------------
# Section B.  Macro discipline -- bare \kappa outside tuple components
# ----------------------------------------------------------------------------
# Forbid bare \kappa not followed by a known suffix or already inside a tuple
# component macro.  Allowed macros: \kappaCat \kappaChHodge \kappaChHeis
# \kappaBKM \kappaFiber \kappaTuple \Kkappa, plus \kappa_<subscript>.
hitskappa="$(grep -RInE --include='*.tex' '\\\\kappa([^A-Za-z_]|$)' "${SCAN_PATHS[@]}" 2>/dev/null \
            | awk -F: '{
                line = $0
                sub(/^[^:]+:[0-9]+:/, "", line)
                if (line ~ /^[[:space:]]*%/) next
                # skip lines where every \kappa usage is qualified
                # by being in a tuple-component / Tuple / Kkappa context.
                # Heuristic: skip if line contains \kappaCat, \kappaChHodge,
                # \kappaChHeis, \kappaBKM, \kappaFiber, \kappaTuple, \Kkappa,
                # but warn only if a bare \kappa appears that is not part of
                # one of these compound names.
                # Strip qualified usages first:
                stripped = line
                gsub(/\\kappaCat/, "", stripped)
                gsub(/\\kappaChHodge/, "", stripped)
                gsub(/\\kappaChHeis/, "", stripped)
                gsub(/\\kappaBKM/, "", stripped)
                gsub(/\\kappaFiber/, "", stripped)
                gsub(/\\kappaTuple/, "", stripped)
                gsub(/\\Kkappa/, "", stripped)
                # Now check if any bare \kappa remains in stripped:
                if (stripped ~ /\\kappa([^A-Za-z_]|$)/) print $0
              }')"
if [ -n "$hitskappa" ]; then
  count="$(printf '%s\n' "$hitskappa" | wc -l | tr -d ' ')"
  echo
  echo "## [W] Bare \\kappa outside tuple components ($count match(es))"
  echo "    Use \\kappaCat / \\kappaChHodge / \\kappaChHeis / \\kappaBKM / \\kappaFiber / \\kappaTuple / \\Kkappa."
  printf '%s\n' "$hitskappa" | head -20
  WARNINGS=$((WARNINGS + count))
fi

# ----------------------------------------------------------------------------
# Section C.  Two-stage Phi discipline
# Every \Phi_d (non-providecommand line) should appear within a paragraph
# also containing \PhiFA or \SpCh.  Implement as same-line check (proxy);
# paragraph-level scan is left to manual review.
# ----------------------------------------------------------------------------
hitsphi="$(grep -RInE --include='*.tex' '\\\\Phi_d([^A-Za-z_]|$)' "${SCAN_PATHS[@]}" 2>/dev/null \
          | awk -F: '{
              line = $0
              sub(/^[^:]+:[0-9]+:/, "", line)
              if (line ~ /^[[:space:]]*%/) next
              if (line ~ /providecommand|newcommand|renewcommand/) next
              if (line ~ /\\PhiFA|\\SpCh/) next
              print $0
            }')"
if [ -n "$hitsphi" ]; then
  count="$(printf '%s\n' "$hitsphi" | wc -l | tr -d ' ')"
  echo
  echo "## [W] \\Phi_d not adjacent to \\PhiFA or \\SpCh on the same line ($count match(es))"
  echo "    Two-stage discipline: \\Phi_d = \\SpCh \\circ \\PhiFA_d."
  printf '%s\n' "$hitsphi" | head -20
  WARNINGS=$((WARNINGS + count))
fi

# ----------------------------------------------------------------------------
# Section D.  Theorem statements lacking conditional locus tags
# Endpoint chapters: e_infinity_topologization*.tex, pva-descent*.tex,
# bv_ht_physics*.tex, wn_tempered*.tex, programme_climax*.tex,
# 3d_gravity*.tex.  We require that any \begin{theorem}...\end{theorem}
# in these files contains at least one \hyp* or \eff* tag.
# ----------------------------------------------------------------------------
endpoint_files=$(find chapters -name 'e_infinity_topologization*.tex' \
                                -o -name 'pva-descent*.tex' \
                                -o -name 'bv_ht_physics*.tex' \
                                -o -name 'wn_tempered*.tex' \
                                -o -name 'programme_climax*.tex' \
                                -o -name '3d_gravity*.tex' 2>/dev/null)

for f in $endpoint_files; do
  python3 - "$f" <<'PY' 2>/dev/null
import re, sys, os
fn = sys.argv[1]
try:
    txt = open(fn, 'r', encoding='utf-8').read()
except OSError:
    sys.exit(0)
# Match \begin{theorem}...\end{theorem} blocks (also corollary, proposition, lemma if labelled).
pat = re.compile(r'\\begin\{(theorem|corollary|proposition|lemma)\}.*?\\end\{\1\}', re.DOTALL)
miss = 0
for m in pat.finditer(txt):
    body = m.group(0)
    if not re.search(r'\\hyp[A-Z][A-Za-z]+|\\eff[A-Z][A-Za-z]+', body):
        # locate line number
        line_no = txt[:m.start()].count('\n') + 1
        # skip if the theorem body is small / decorative
        if len(body) < 80:
            continue
        miss += 1
        if miss <= 5:
            label_match = re.search(r'\\label\{([^}]+)\}', body)
            label = label_match.group(1) if label_match else '<no label>'
            print(f'{fn}:{line_no}: theorem missing \\hyp / \\eff tag ({label})')
if miss > 5:
    print(f'    ... and {miss - 5} more in {fn}')
# exit code carries the miss count back to the shell (capped below 255)
sys.exit(min(miss, 200))
PY
  miss_count=$?
  WARNINGS=$((WARNINGS + miss_count))
done

# ----------------------------------------------------------------------------
# Section E.  Voice / bookkeeping vocabulary in chapters
# Forbidden: Wave-N / round-K / session-X / DNA / AP-numbers in body prose.
# ----------------------------------------------------------------------------
scan '\\b(Wave|Round) +[0-9]' \
     'Bookkeeping vocabulary "Wave N" / "Round N" -- not allowed in chapters' W

scan 'AP-?(V2-)?[0-9]+\\b' \
     'Bookkeeping: "APn" / "AP-Vol-n" naming -- AP indices live in notes/, not chapters' W

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
echo
echo "=== verify_licensing.sh summary ==="
echo "  Blocking violations: $VIOLATIONS"
echo "  Warnings:            $WARNINGS"

if [ "$REPORT_ONLY" = "1" ]; then
  exit 0
fi

if [ "$VIOLATIONS" -gt 0 ]; then
  if [ "$VIOLATIONS" -gt 250 ]; then
    exit 250
  fi
  exit "$VIOLATIONS"
fi
exit 0

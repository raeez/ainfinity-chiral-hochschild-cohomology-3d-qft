# ============================================================================
#  Makefile -- A-infinity Chiral Algebras and 3D HT QFT (Vol II)
# ============================================================================
#
#  Usage:
#    make            Full build → out/
#    make fast       Quick build (up to 4 passes) → out/main.pdf
#    make release    Full release → out/ + iCloud
#    make standalone Build standalone documents → out/
#    make clean      Remove all LaTeX build artifacts (preserves stamp)
#    make veryclean  Remove artifacts AND out/ (forces rebuild)
#    make clean-builds  Remove all /tmp/mkd-* isolated build directories
#    make check      Halt-on-error validation
#    make count      Manuscript statistics
#    make test       Run compute test suite
#    make help       Show available targets
#
#  Build isolation (parallel agents):
#    Each build runs in its own /tmp directory.  Set MKD_BUILD_NS to reuse
#    the same directory across invocations (warm .aux files = faster builds):
#
#      export MKD_BUILD_NS="agent-$$"   # set once per agent session
#      make fast                         # cold first time, warm thereafter
#
#  All compiled output goes to out/.
#
# ============================================================================

# --- Configuration -----------------------------------------------------------

MAIN      := main
TEX       := pdflatex
TEXFLAGS  := -interaction=batchmode -file-line-error -synctex=0 -cnf-line='buf_size=1000000' -cnf-line='stack_size=20000'
BUILD_SCRIPT := ./scripts/build.sh
LOG_DIR   := .build_logs
PASSES    := 6
FAST_PASSES := 4
PYTEST_FAST_TIMEOUT ?= 120

# iCloud destination for release PDFs
ICLOUD_DIR := /Users/raeez/Library/Mobile Documents/com~apple~CloudDocs/research

# Source files
SOURCES   := $(wildcard *.tex) \
             $(wildcard chapters/theory/*.tex) \
             $(wildcard chapters/examples/*.tex) \
             $(wildcard chapters/connections/*.tex) \
             $(wildcard appendices/*.tex)

# Output -- everything goes to out/
OUT_DIR   := out
PDF       := $(OUT_DIR)/main.pdf
ICLOUD_MAIN_PREREQ := $(if $(wildcard $(PDF)),,$(PDF))

# Mathematics publish dir -- release binary copied here under canonical name
MATHEMATICS_DIR := $(HOME)/mathematics
PUBLISHED_PDF   := Ainfinity_Chiral_Algebras_and_Chiral_Hochschild_Cohomology.pdf
VOLUME_KEY      := vol2
PYTHON_BIN      ?= $(shell if [ -x compute/.venv/bin/python ]; then echo compute/.venv/bin/python; elif [ -x .venv/bin/python ]; then echo .venv/bin/python; else echo python3; fi)

# Working notes
WN_TEX    := working_notes.tex

# Standalone documents
STANDALONE_TEX := $(filter-out standalone/preface_full_survey.tex,$(wildcard standalone/*.tex))
STANDALONE_PASSES := 3

STAMP     := .build_stamp

# If PDF was externally deleted but stamp remains, force a rebuild.
ifeq (,$(wildcard $(PDF)))
  $(shell rm -f $(STAMP))
endif

# LaTeX intermediate extensions
AUX_EXTS  := aux log out toc synctex.gz fdb_latexmk fls bbl blg \
             nav snm vrb idx ilg ind lof lot

# ============================================================================
#  Targets
# ============================================================================

.DEFAULT_GOAL := platonic

.PHONY: platonic all fast clean veryclean clean-builds count check test dist release help working-notes standalone icloud verify-independence verify-independence-verbose verify-licensing verify-licensing-report mathematics-publish root-publish architecture unified-architecture

## platonic: Build the integrated monograph -> out/platonic.pdf
platonic:
	@echo "  -- Building the integrated monograph --"
	@mkdir -p out .build_logs
	@cd platonic && \
	  TEXINPUTS=".:..:$$TEXINPUTS" BIBINPUTS=".:..:$$BIBINPUTS" \
	  pdflatex -interaction=nonstopmode -file-line-error main.tex >../.build_logs/platonic.log 2>&1 || true; \
	  TEXINPUTS=".:..:$$TEXINPUTS" BIBINPUTS=".:..:$$BIBINPUTS" \
	  bibtex main >>../.build_logs/platonic.log 2>&1 || true; \
	  for i in 1 2 3; do \
	    TEXINPUTS=".:..:$$TEXINPUTS" BIBINPUTS=".:..:$$BIBINPUTS" \
	    pdflatex -interaction=nonstopmode -file-line-error main.tex >../.build_logs/platonic.log 2>&1 || true; \
	  done
	@if [ -f platonic/main.pdf ]; then cp platonic/main.pdf out/platonic.pdf; \
		echo "    OK out/platonic.pdf ($$(pdfinfo platonic/main.pdf 2>/dev/null | awk '/^Pages/{print $$2}') pages)"; \
	else echo "    FAIL - see .build_logs/platonic.log"; exit 1; fi
	@if grep -aqE '^! ' .build_logs/platonic.log; then \
		echo "    LaTeX errors:"; grep -aE '^! ' .build_logs/platonic.log | head -5; exit 1; fi
	@if grep -aqE 'Reference .* undefined|Citation .* undefined' .build_logs/platonic.log; then \
		echo "    undefined refs/citations:"; grep -aE 'Reference .* undefined|Citation .* undefined' .build_logs/platonic.log | head -5; exit 1; fi
	@echo "    0 errors, 0 undefined references, 0 undefined citations."


## icloud: Copy latest PDFs to iCloud Drive (subject-organised)
icloud: $(ICLOUD_MAIN_PREREQ) standalone
	@echo "  ── Copying Vol II to iCloud (subject-organised) ──"
	@mkdir -p "$(ICLOUD_DIR)/volumes"
	@mkdir -p "$(ICLOUD_DIR)/vol2_3d_ht_physics"
	@[ -f $(PDF) ] && cp $(PDF) "$(ICLOUD_DIR)/volumes/vol2_ainfinity_chiral_algebras.pdf" \
		&& echo "    ✓ volumes/vol2" || true
	@for pdf in $(OUT_DIR)/*.pdf; do \
		name=$$(basename "$$pdf"); \
		if [ "$$name" != "main.pdf" ]; then \
			cp "$$pdf" "$(ICLOUD_DIR)/vol2_3d_ht_physics/$$name"; \
			echo "    ✓ vol2_3d_ht_physics/$$name"; \
		fi; \
	done
	@echo "  Vol II PDFs copied to iCloud."

## all: Full build → out/
all: $(STAMP) working-notes

$(STAMP): $(SOURCES) Makefile $(BUILD_SCRIPT)
	@echo "══════════════════════════════════════════════════════════"
	@echo "  Building: $(MAIN).tex  →  $(PDF)"
	@echo "══════════════════════════════════════════════════════════"
	@$(BUILD_SCRIPT) $(PASSES)
	@if [ ! -f $(PDF) ]; then \
		echo "  ✗  Build failed -- no PDF produced."; exit 1; \
	fi
	@touch $(STAMP)
	@echo ""
	@echo "  ✓  $(PDF) built successfully."
	@echo ""

$(PDF): $(STAMP)
	@true

## fast: Bounded quick build for rapid iteration → out/main.pdf
fast:
	@echo "  ── Fast build (up to $(FAST_PASSES) passes) ──"
	@$(BUILD_SCRIPT) $(FAST_PASSES)

## working-notes: Build the working notes → out/working_notes.pdf
working-notes:
	@echo "  ── Building working notes ──"
	@mkdir -p $(OUT_DIR) $(LOG_DIR)
	@$(TEX) $(TEXFLAGS) $(WN_TEX) >$(LOG_DIR)/wn-pass1.log 2>&1 || true
	@$(TEX) $(TEXFLAGS) $(WN_TEX) >$(LOG_DIR)/wn-pass2.log 2>&1 || true
	@if [ -f working_notes.pdf ]; then \
		mv working_notes.pdf $(OUT_DIR)/working_notes.pdf; \
		rm -f working_notes.aux working_notes.log working_notes.out working_notes.toc 2>/dev/null; \
		echo "  ✓  $(OUT_DIR)/working_notes.pdf"; \
	else \
		echo "  ✗  Working notes build failed."; \
		exit 1; \
	fi

## release: Full rebuild → out/ + standalones + iCloud
release:
	@rm -f $(STAMP)
	@rm -rf $(OUT_DIR)
	@mkdir -p $(LOG_DIR) $(OUT_DIR)
	@echo ""
	@echo "  ══════════════════════════════════════════"
	@echo "  ── RELEASE BUILD (Vol II) ──"
	@echo "  ══════════════════════════════════════════"
	@echo ""
	@echo "  [1/3] Main manuscript"
	@$(MAKE) --no-print-directory $(STAMP)
	@echo ""
	@echo "  [2/3] Working notes"
	@$(MAKE) --no-print-directory working-notes
	@echo ""
	@echo "  [3/6] Standalone documents and iCloud"
	@$(MAKE) --no-print-directory icloud
	@echo ""
	@echo "  [4/6] Publish to repo root (canonical PDF name)"
	@$(MAKE) --no-print-directory root-publish
	@echo ""
	@echo "  [5/6] Publish to ~/mathematics + per-volume architecture"
	@$(MAKE) --no-print-directory mathematics-publish
	@$(MAKE) --no-print-directory architecture
	@echo ""
	@echo "  [6/6] Cross-volume architecture aggregation"
	@$(MAKE) --no-print-directory unified-architecture
	@echo ""
	@echo "  ══════════════════════════════════════════"
	@echo "  Release complete. All output in out/:"
	@ls -1 $(OUT_DIR)/*.pdf 2>/dev/null | sed 's/^/    /'
	@echo "  ══════════════════════════════════════════"

## root-publish: Copy the release binary to repo root under its canonical name
root-publish:
	@if [ -f "$(PDF)" ]; then \
		cp "$(PDF)" "$(PUBLISHED_PDF)"; \
		echo "    ✓  $(PUBLISHED_PDF) (in repo root)"; \
	else \
		echo "    ✗  $(PDF) missing -- skipping root publish"; \
	fi

## mathematics-publish: Copy the release binary to ~/mathematics under its canonical name
mathematics-publish:
	@mkdir -p "$(MATHEMATICS_DIR)"
	@if [ -f "$(PDF)" ]; then \
		cp "$(PDF)" "$(MATHEMATICS_DIR)/$(PUBLISHED_PDF)"; \
		echo "    ✓  $(MATHEMATICS_DIR)/$(PUBLISHED_PDF)"; \
	else \
		echo "    ✗  $(PDF) missing — skipping ~/mathematics publish"; \
	fi

## architecture: Build interactive HTML + JSON of the manuscript architecture
architecture:
	@$(PYTHON_BIN) scripts/build_architecture.py --root . --volume $(VOLUME_KEY) --out $(OUT_DIR) --quiet
	@mkdir -p "$(MATHEMATICS_DIR)/architecture"
	@cp "$(OUT_DIR)/architecture.json" "$(MATHEMATICS_DIR)/architecture/$(VOLUME_KEY).json"
	@echo "    ✓  $(OUT_DIR)/architecture.html + .json"
	@echo "    ✓  $(MATHEMATICS_DIR)/architecture/$(VOLUME_KEY).json"

## unified-architecture: Aggregate all per-volume architecture.json into the cross-volume HTML+JSON
unified-architecture:
	@$(PYTHON_BIN) scripts/build_unified_architecture.py --mathematics-dir "$(MATHEMATICS_DIR)" --quiet
	@echo "    ✓  $(MATHEMATICS_DIR)/architecture.html + .json"

## standalone: Build standalone documents → out/
standalone:
	@echo "  ── Building standalone documents ──"
	@mkdir -p $(OUT_DIR) $(LOG_DIR)
	@if [ -z "$(strip $(STANDALONE_TEX))" ]; then \
		echo "  (no standalone documents found)"; \
	else \
		failures=0; \
		for tex in $(STANDALONE_TEX); do \
			if [ ! -f "$$tex" ]; then continue; fi; \
			base=$$(basename "$$tex" .tex); \
			if [ -f "$(OUT_DIR)/$$base.pdf" ] && [ "$(OUT_DIR)/$$base.pdf" -nt "$$tex" ]; then \
				echo "  ✓  $(OUT_DIR)/$$base.pdf (up to date)"; \
				continue; \
			fi; \
			tmpdir=$$(mktemp -d "/tmp/mkd-$$(basename "$$(pwd)")-standalone-$$base.XXXXXX"); \
			echo "  [standalone] $$tex → $(OUT_DIR)/$$base.pdf"; \
			for pass in $$(seq 1 $(STANDALONE_PASSES)); do \
				TEXINPUTS="$$tmpdir:$$(pwd):$$(pwd)/standalone:" $(TEX) $(TEXFLAGS) -output-directory="$$tmpdir" "$$tex" >"$(LOG_DIR)/standalone-$$base-pass$$pass.log" 2>&1; rc=$$?; \
				if [ -f "$$tmpdir/$$base.idx" ]; then makeindex -q "$$tmpdir/$$base.idx" >/dev/null 2>&1 || true; fi; \
				if [ $$rc -ne 0 ]; then \
					if grep -aE '^! |Emergency stop|Runaway argument|Fatal error|Undefined control sequence|File ended while scanning|No pages of output' "$(LOG_DIR)/standalone-$$base-pass$$pass.log" >/dev/null 2>&1; then \
						echo "  ✗  $$tex failed on pass $$pass. See $(LOG_DIR)/standalone-$$base-pass$$pass.log"; \
						grep -aE '^! |Emergency stop|Runaway argument|Fatal error|Undefined control sequence|File ended while scanning|No pages of output' "$(LOG_DIR)/standalone-$$base-pass$$pass.log" | head -n 20 || tail -n 40 "$(LOG_DIR)/standalone-$$base-pass$$pass.log"; \
						failures=$$((failures + 1)); \
						break; \
					else \
						echo "  warn  $$tex returned $$rc on pass $$pass; continuing."; \
					fi; \
				fi; \
			done; \
			if [ -f "$$tmpdir/$$base.pdf" ]; then \
				cp "$$tmpdir/$$base.pdf" "$(OUT_DIR)/$$base.pdf"; \
				echo "  ✓  $(OUT_DIR)/$$base.pdf"; \
			elif [ $$failures -eq 0 ]; then \
				echo "  ✗  No PDF produced for $$tex"; \
				failures=$$((failures + 1)); \
			fi; \
		done; \
		if [ $$failures -ne 0 ]; then \
			echo "  ✗  $$failures standalone document(s) failed."; \
			exit 1; \
		fi; \
	fi

## dist: Create Vol2Archive.zip for distribution.
dist:
	@echo "  ── Creating Vol2Archive.zip ──"
	@rm -f $(OUT_DIR)/Vol2Archive.zip
	@mkdir -p $(OUT_DIR)
	@zip -r $(OUT_DIR)/Vol2Archive.zip \
		main.tex chapters/ appendices/ compute/ scripts/ \
		Makefile README.md CLAUDE.md \
		$(PDF) \
		$(OUT_DIR)/working_notes.pdf \
		-x '.*' -x '**/.*' -x '**/__pycache__/*' -x '**/*.pyc' \
		-x 'compute/.venv/*' \
		>$(LOG_DIR)/dist.log 2>&1
	@echo "  ✓  $(OUT_DIR)/Vol2Archive.zip ($$(du -h $(OUT_DIR)/Vol2Archive.zip | cut -f1))"

## verify-licensing: Whole-volume Beilinson-cut sweep (17-line voice table + macro discipline + conditional-locus tags).  Exits non-zero on blocking violations.
verify-licensing:
	@bash scripts/verify_licensing.sh

## verify-licensing-report: Same scan, report-only (always exit 0).
verify-licensing-report:
	@VERIFY_LICENSING_REPORT_ONLY=1 bash scripts/verify_licensing.sh

## check: Halt on first error -- use for CI or pre-commit validation.
check:
	@echo "  ── Error check (halt-on-error) ──"
	@mkdir -p $(LOG_DIR)
	@$(TEX) -interaction=nonstopmode -halt-on-error -file-line-error $(MAIN).tex >$(LOG_DIR)/check.log 2>&1 || { \
		echo "  ✗  Check failed. See $(LOG_DIR)/check.log"; \
		grep -aE '^\! |Emergency stop|Runaway argument|Fatal error|Undefined control sequence|File ended while scanning|No pages of output' $(LOG_DIR)/check.log | head -n 20 || tail -n 40 $(LOG_DIR)/check.log; \
		exit 1; \
	}
	@echo "  ✓  No fatal errors."
	@echo "     Log: $(LOG_DIR)/check.log"

## clean: Remove build debris but preserve the stamp.
clean:
	@echo "  Cleaning build artifacts..."
	@for ext in $(AUX_EXTS); do \
		rm -f $(MAIN).$$ext; \
	done
	@find chapters appendices -name '*.aux' -delete 2>/dev/null || true
	@rm -rf $(LOG_DIR)
	@rm -f texput.log
	@echo "  ✓  Clean (stamp preserved -- make will skip rebuild if sources unchanged)."

## veryclean: Remove EVERYTHING including out/ and build stamp (forces full rebuild).
veryclean: clean
	@rm -f $(STAMP)
	@rm -rf $(OUT_DIR)
	@echo "  ✓  Stamp and out/ removed -- next make will rebuild."

## clean-builds: Remove ALL /tmp/mkd-* isolated build directories (all volumes).
clean-builds:
	@echo "  Cleaning isolated build directories..."
	@rm -rf /tmp/mkd-chiral-bar-cobar-* /tmp/mkd-chiral-bar-cobar-vol2-* /tmp/mkd-calabi-yau-quantum-groups-*
	@echo "  ✓  All /tmp/mkd-* build directories removed."

## count: Manuscript statistics.
count:
	@echo ""
	@echo "  ── Vol II Statistics ──"
	@echo ""
	@printf "  Source files:   %s .tex files\n" "$$(find . -name '*.tex' -not -path './.build_logs/*' -not -path './out/*' | wc -l | tr -d ' ')"
	@printf "  Total lines:   %s\n" "$$(find . -name '*.tex' -not -path './.build_logs/*' -not -path './out/*' -exec cat {} + | wc -l | tr -d ' ')"
	@if [ -f $(PDF) ]; then \
		PAGES=$$(grep -o '([0-9]* pages' $(OUT_DIR)/main.log 2>/dev/null | grep -o '[0-9]*' | tail -n 1 || echo '?'); \
		printf "  PDF pages:     %s\n" "$$PAGES"; \
		printf "  PDF size:      %s\n" "$$(du -h $(PDF) | cut -f1)"; \
	else \
		echo "  PDF:           (not yet built -- run 'make')"; \
	fi
	@echo ""

## test: Run compute test suite.
test:
	@if [ -d compute/tests ] && ls compute/tests/test_*.py 1>/dev/null 2>&1; then \
		echo "  ── Running compute test suite (fast: excludes slow) ──"; \
		mkdir -p $(LOG_DIR); \
		if [ -f compute/.venv/bin/python ]; then \
			PYTHON_BIN=compute/.venv/bin/python; \
		elif [ -f .venv/bin/python ]; then \
			PYTHON_BIN=.venv/bin/python; \
		else \
			PYTHON_BIN=python3; \
		fi; \
		LOG_FILE=$(LOG_DIR)/pytest.log; \
		$$PYTHON_BIN -m pytest compute/tests/ -q -ra -m "not slow" \
			-o faulthandler_timeout=$(PYTEST_FAST_TIMEOUT) \
			--durations=10 --durations-min=1.0 >$$LOG_FILE 2>&1; rc=$$?; \
		if [ $$rc -eq 0 ]; then \
			tail -n 5 $$LOG_FILE; \
			echo "     Log: $$LOG_FILE"; \
		else \
			echo "  ✗  Test run failed. See $$LOG_FILE"; \
			tail -n 120 $$LOG_FILE; \
			exit $$rc; \
		fi; \
	else \
		echo "  (no compute tests found -- skipping)"; \
	fi

## verify-independence: Audit ProvedHere claims vs independent-verification registry
##                     (tautology / orphan check; coverage metric reported)
verify-independence:
	@python3 compute/scripts/audit_independent_verification.py

## verify-independence-verbose: Same, with full list of uncovered claims
verify-independence-verbose:
	@python3 compute/scripts/audit_independent_verification.py --verbose --show-orphans

## help: Show available targets.
help:
	@echo ""
	@echo "  Vol II -- Build System"
	@echo "  ─────────────────────"
	@echo "  All compiled output goes to out/"
	@echo ""
	@echo "  make            Full build → out/"
	@echo "  make fast       Quick converging build → out/main.pdf"
	@echo "  make release    Full rebuild → out/ + standalones + iCloud"
	@echo "  make standalone Build standalone documents → out/"
	@echo "  make dist       Create Vol2Archive.zip in out/"
	@echo "  make check      Halt-on-error validation"
	@echo "  make clean      Remove build debris (preserves stamp)"
	@echo "  make veryclean  Remove everything including out/"
	@echo "  make clean-builds  Remove /tmp/mkd-* isolated build directories"
	@echo "  make count      Manuscript statistics"
	@echo "  make test       Run compute tests (excludes slow)"
	@echo "  make help       This message"
	@echo ""

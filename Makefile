# ============================================================================
#  Makefile — A-infinity Chiral Algebras and 3D HT QFT (Vol II)
# ============================================================================
#
#  Usage:
#    make            Full build → out/
#    make fast       Quick build (up to 4 passes) for rapid iteration
#    make clean      Remove all LaTeX build artifacts (preserves stamp)
#    make veryclean  Remove artifacts AND compiled PDFs (forces rebuild)
#    make check      Halt-on-error validation
#    make count      Manuscript statistics
#    make test       Run compute test suite
#    make help       Show available targets
#
# ============================================================================

# --- Configuration -----------------------------------------------------------

MAIN      := main
TEX       := pdflatex
TEXFLAGS  := -interaction=batchmode -file-line-error -synctex=0 -cnf-line='buf_size=1000000' -cnf-line='stack_size=20000'
LOG_DIR   := .build_logs
PASSES    := 6
FAST_PASSES := 4
PYTEST_FAST_TIMEOUT ?= 120

# Source files
SOURCES   := $(wildcard *.tex) \
             $(wildcard chapters/theory/*.tex) \
             $(wildcard chapters/examples/*.tex) \
             $(wildcard chapters/connections/*.tex) \
             $(wildcard appendices/*.tex)

# Output
PDF       := $(MAIN).pdf
STAMP     := .build_stamp
OUT_DIR   := out
OUT_PDF   := $(OUT_DIR)/ainfinity_chiral_algebras.pdf

# Working notes
WN_DIR    := archive/source_tex
WN_TEX    := $(WN_DIR)/working_notes.tex
WN_PDF    := $(WN_DIR)/working_notes.pdf
OUT_WN    := $(OUT_DIR)/working_notes.pdf

# If PDF was externally deleted but stamp remains, force a rebuild.
ifeq (,$(wildcard $(PDF)))
  $(shell rm -f $(STAMP))
endif

# LaTeX intermediate extensions
AUX_EXTS  := aux log out toc synctex.gz fdb_latexmk fls bbl blg \
             nav snm vrb idx ilg ind lof lot

# ============================================================================
#  Helper functions
# ============================================================================

define count_matches
$(shell grep -aEc '$(1)' $(2) 2>/dev/null || echo 0)
endef

# ============================================================================
#  Targets
# ============================================================================

.PHONY: all fast clean veryclean count check test dist release help publish working-notes

## all: Full build → out/
##   Idempotent: no-op if no .tex files changed since last successful build.
all: $(STAMP) working-notes publish

$(STAMP): $(SOURCES)
	@echo "══════════════════════════════════════════════════════════"
	@echo "  Building: $(MAIN).tex  →  $(PDF)"
	@echo "  Engine:   $(TEX) (up to $(PASSES) passes)"
	@echo "══════════════════════════════════════════════════════════"
	@mkdir -p $(LOG_DIR)
	@converged=0; \
	for pass in $$(seq 1 $(PASSES)); do \
		echo "── Pass $$pass / $(PASSES) ──"; \
		$(TEX) $(TEXFLAGS) $(MAIN).tex >$(LOG_DIR)/pass$$pass.log 2>&1; \
		rc=$$?; \
		pages=$$(grep 'Output written' $(MAIN).log 2>/dev/null | sed 's/.*(\([0-9]*\) pages.*/\1/' | tail -n 1 || echo '?'); \
		undef_cit=$$(grep -aEc 'Citation.*undefined' $(MAIN).log 2>/dev/null || echo 0); \
		undef_ref=$$(grep -aEc 'Reference.*undefined' $(MAIN).log 2>/dev/null || echo 0); \
		overfull=$$(grep -aEc 'Overfull \\hbox' $(MAIN).log 2>/dev/null || echo 0); \
		underfull=$$(grep -aEc 'Underfull \\hbox|Underfull \\vbox' $(MAIN).log 2>/dev/null || echo 0); \
		rerun=$$(grep -aEc 'Label\(s\) may have changed|Package rerunfilecheck Warning' $(MAIN).log 2>/dev/null || echo 0); \
		echo "   $${pages}pp, $$undef_cit undef citations, $$undef_ref undef references, $$rerun rerun requests, $$overfull overfull, $$underfull underfull"; \
		if [ $$rc -ne 0 ] && { [ ! -f $(PDF) ] || ! grep -aq 'Output written' $(MAIN).log 2>/dev/null; }; then \
			echo "  ✗  Pass $$pass failed."; \
			grep -aE '^\! |Emergency stop|Runaway argument|Fatal error|Undefined control sequence|File ended while scanning|No pages of output' $(LOG_DIR)/pass$$pass.log | head -n 10; \
			exit $$rc; \
		fi; \
		if [ "$$pass" -ge 2 ] && [ "$$rerun" = "0" ] && [ "$$undef_cit" = "0" ] && [ "$$undef_ref" = "0" ]; then \
			echo "✓ Converged after $$pass passes."; \
			converged=1; \
			break; \
		fi; \
	done; \
	if [ "$$converged" = "0" ]; then \
		echo "⚠ Did not fully converge after $(PASSES) passes (Cit=$$undef_cit, Ref=$$undef_ref, Rerun=$$rerun)."; \
	fi
	@if [ ! -f $(PDF) ]; then \
		echo "  ✗  Build failed — no PDF produced."; exit 1; \
	fi
	@touch $(STAMP)
	@echo ""
	@echo "  ✓  $(PDF) built successfully."
	@echo "     Logs: $(LOG_DIR)/ and $(MAIN).log"
	@echo ""

## fast: Bounded quick build for rapid iteration.
fast:
	@echo "  ── Fast build (up to $(FAST_PASSES) passes) ──"
	@mkdir -p $(LOG_DIR)
	@converged=0; \
	for pass in $$(seq 1 $(FAST_PASSES)); do \
		echo "── Pass $$pass / $(FAST_PASSES) ──"; \
		$(TEX) $(TEXFLAGS) $(MAIN).tex >$(LOG_DIR)/fast-pass$$pass.log 2>&1; \
		rc=$$?; \
		pages=$$(grep 'Output written' $(MAIN).log 2>/dev/null | sed 's/.*(\([0-9]*\) pages.*/\1/' | tail -n 1 || echo '?'); \
		undef_cit=$$(grep -aEc 'Citation.*undefined' $(MAIN).log 2>/dev/null || echo 0); \
		undef_ref=$$(grep -aEc 'Reference.*undefined' $(MAIN).log 2>/dev/null || echo 0); \
		rerun=$$(grep -aEc 'Label\(s\) may have changed|Package rerunfilecheck Warning' $(MAIN).log 2>/dev/null || echo 0); \
		overfull=$$(grep -aEc 'Overfull \\hbox' $(MAIN).log 2>/dev/null || echo 0); \
		echo "   $${pages}pp, $$undef_cit undef cit, $$undef_ref undef ref, $$rerun rerun, $$overfull overfull"; \
		if [ $$rc -ne 0 ] && { [ ! -f $(PDF) ] || ! grep -aq 'Output written' $(MAIN).log 2>/dev/null; }; then \
			echo "  ✗  Pass $$pass failed."; \
			grep -aE '^\! |Emergency stop|Fatal error|Undefined control sequence' $(LOG_DIR)/fast-pass$$pass.log | head -n 10; \
			exit $$rc; \
		fi; \
		if [ "$$pass" -ge 2 ] && [ "$$rerun" = "0" ] && [ "$$undef_cit" = "0" ] && [ "$$undef_ref" = "0" ]; then \
			echo "✓ Converged after $$pass passes."; \
			converged=1; \
			break; \
		fi; \
	done; \
	if [ "$$converged" = "0" ] && [ $(FAST_PASSES) -gt 1 ]; then \
		echo "⚠ Did not converge in $(FAST_PASSES) passes (Cit=$$undef_cit, Ref=$$undef_ref, Rerun=$$rerun)."; \
	fi
	@echo "     Logs: $(LOG_DIR)/ and $(MAIN).log"

## publish: Copy final PDFs to out/ (does not trigger a rebuild).
publish:
	@mkdir -p $(OUT_DIR)
	@if [ -f $(PDF) ]; then cp $(PDF) $(OUT_PDF); echo "  ✓  $(OUT_PDF)"; \
	else echo "  ⚠  $(PDF) not found — run 'make fast' first."; fi
	@if [ -f $(WN_PDF) ]; then cp $(WN_PDF) $(OUT_WN); echo "  ✓  $(OUT_WN)"; fi

## working-notes: Build the working notes (standalone document).
working-notes: $(OUT_WN)

$(OUT_WN): $(WN_TEX)
	@echo "  ── Building working notes ──"
	@mkdir -p $(OUT_DIR) $(LOG_DIR)
	@cd $(WN_DIR) && \
		$(TEX) $(TEXFLAGS) working_notes.tex >../../$(LOG_DIR)/wn-pass1.log 2>&1 || true && \
		$(TEX) $(TEXFLAGS) working_notes.tex >../../$(LOG_DIR)/wn-pass2.log 2>&1 || true
	@if [ -f $(WN_PDF) ]; then \
		cp $(WN_PDF) $(OUT_WN); \
		echo "  ✓  $(OUT_WN)"; \
	else \
		echo "  ✗  Working notes build failed."; \
		exit 1; \
	fi

## release: Clean rebuild + named release PDF + working notes at root.
release:
	@rm -f $(STAMP) $(PDF) $(WN_PDF)
	@rm -rf $(OUT_DIR)
	@mkdir -p $(LOG_DIR) $(OUT_DIR)
	@$(MAKE) --no-print-directory $(STAMP) working-notes publish
	@cp $(PDF) Ainfinity_Chiral_Algebras_and_Chiral_Hochschild_Cohomology.pdf
	@echo "  ✓  Ainfinity_Chiral_Algebras_and_Chiral_Hochschild_Cohomology.pdf"

## dist: Create Vol2Archive.zip for distribution.
dist: publish
	@echo "  ── Creating Vol2Archive.zip ──"
	@rm -f $(OUT_DIR)/Vol2Archive.zip
	@mkdir -p $(OUT_DIR)
	@zip -r $(OUT_DIR)/Vol2Archive.zip \
		main.tex chapters/ appendices/ compute/ \
		Makefile README.md CLAUDE.md \
		$(OUT_DIR)/ainfinity_chiral_algebras.pdf \
		-x '.*' -x '**/.*' -x '**/__pycache__/*' -x '**/*.pyc' \
		-x 'compute/.venv/*' \
		>$(LOG_DIR)/dist.log 2>&1
	@echo "  ✓  $(OUT_DIR)/Vol2Archive.zip ($$(du -h $(OUT_DIR)/Vol2Archive.zip | cut -f1))"

## check: Halt on first error — use for CI or pre-commit validation.
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
	@echo "  ✓  Clean (stamp preserved — make will skip rebuild if sources unchanged)."

## veryclean: Remove EVERYTHING including PDF, out/, and build stamp (forces full rebuild).
veryclean: clean
	@rm -f $(PDF) $(STAMP)
	@rm -rf $(OUT_DIR)
	@echo "  ✓  Stamp, PDF, and out/ removed — next make will rebuild."

## count: Manuscript statistics.
count:
	@echo ""
	@echo "  ── Vol II Statistics ──"
	@echo ""
	@printf "  Source files:   %s .tex files\n" "$$(find . -name '*.tex' -not -path './.build_logs/*' | wc -l | tr -d ' ')"
	@printf "  Total lines:   %s\n" "$$(find . -name '*.tex' -not -path './.build_logs/*' -exec cat {} + | wc -l | tr -d ' ')"
	@if [ -f $(PDF) ]; then \
		PAGES=$$(grep 'Output written' $(MAIN).log 2>/dev/null | sed 's/.*(\([0-9]*\) pages.*/\1/' | tail -n 1 || echo '?'); \
		printf "  PDF pages:     %s\n" "$$PAGES"; \
		printf "  PDF size:      %s\n" "$$(du -h $(PDF) | cut -f1)"; \
	else \
		echo "  PDF:           (not yet built — run 'make')"; \
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
		echo "  (no compute tests found — skipping)"; \
	fi

## help: Show available targets.
help:
	@echo ""
	@echo "  Vol II — Build System"
	@echo "  ─────────────────────"
	@echo ""
	@echo "  make            Full build → out/ ($(PASSES) passes, stamp-based)"
	@echo "  make fast       Quick converging build (up to $(FAST_PASSES) passes)"
	@echo "  make release    Full rebuild + named release PDF at root"
	@echo "  make dist       Create Vol2Archive.zip in out/"
	@echo "  make check      Halt-on-error validation"
	@echo "  make clean      Remove build debris (preserves stamp)"
	@echo "  make veryclean  Remove everything (forces rebuild)"
	@echo "  make count      Manuscript statistics"
	@echo "  make test       Run compute tests (excludes slow)"
	@echo "  make help       This message"
	@echo ""

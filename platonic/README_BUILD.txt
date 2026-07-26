THE IGUSA SQUARE ROOT — PLATONIC INTEGRATED SOURCE

Entry point: main.tex

Required software:
  - TeX Live 2024 or newer
  - pdflatex
  - bibtex

Build:
  pdflatex -interaction=nonstopmode -halt-on-error main.tex
  bibtex main
  pdflatex -interaction=nonstopmode -halt-on-error main.tex
  pdflatex -interaction=nonstopmode -halt-on-error main.tex

Output: main.pdf

The source tree is self-contained.  It contains the complete integrated TeX,
shared style file, macros, bibliography, and all chapter sources used for the
released monograph.

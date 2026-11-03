.PHONY: install test check-reference figures paper

install:
	python -m pip install -e ".[dev]"

test:
	pytest

check-reference:
	npbc-water check-reference reference_results

figures:
	npbc-water rebuild-figures reference_results

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error SI_R15_R20_R30_analysis_methods.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error SI_R15_R20_R30_analysis_methods.tex

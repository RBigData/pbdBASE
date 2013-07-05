#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pbdBASE-guide.Rnw
bibtex pbdBASE-guide
pdflatex pbdBASE-guide.Rnw
pdflatex pbdBASE-guide.Rnw
pdflatex pbdBASE-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc *.dvi

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/

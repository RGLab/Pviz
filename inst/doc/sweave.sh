#!/bin/bash
R CMD BATCH '--args rnw="Pviz.Rnw"' sweave.R
pdflatex Pviz.tex
rm *.tex
rm *.toc
rm *.aux
rm *.log

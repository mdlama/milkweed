#!/bin/sh
pdflatex Figure1_Model.tex
pdfcrop Figure1_Model.pdf
rm Figure1_Model.pdf
mv Figure1_Model-crop.pdf Figure1_Model.pdf

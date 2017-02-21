#!/bin/sh
pdflatex Box1.tex
pdfcrop Box1.pdf
rm Box1.pdf
mv Box1-crop.pdf Box1.pdf
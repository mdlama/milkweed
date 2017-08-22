#!/bin/sh
pdflatex figure01_model.tex
pdfcrop figure01_model.pdf
rm figure01_model.pdf
mv figure01_model-crop.pdf figure01_model.pdf

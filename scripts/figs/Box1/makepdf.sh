#!/bin/sh
pdflatex $1
pdfcrop $1.pdf
rm $1.pdf
mv $1-crop.pdf $1.pdf

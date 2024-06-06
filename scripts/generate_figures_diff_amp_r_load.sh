#!/bin/sh

echo Generating Figures...

pdflatex ../doc/diff_amp_r_load.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 diff_amp_r_load.pdf diff_amp_r_load.svg
convert -flatten -density 600 diff_amp_r_load.pdf diff_amp_r_load.png

cp diff_amp_r_load.png ../doc/diff_amp_r_load.png
cp diff_amp_r_load.svg ../doc/diff_amp_r_load.svg

rm -rf diff_amp_r_load*

cd ../python/circuits/
python3 diff_amp_r_load_tb.py
cd -

#!/bin/sh

echo Generating Figures...

pdflatex ../doc/mosfet_bias_ex1.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 mosfet_bias_ex1.pdf mosfet_bias_ex1.svg
convert -flatten -density 600 mosfet_bias_ex1.pdf mosfet_bias_ex1.png

cp mosfet_bias_ex1.png ../doc/mosfet_bias_ex1.png
cp mosfet_bias_ex1.svg ../doc/mosfet_bias_ex1.svg

rm -rf mosfet_bias_ex1*

pdflatex ../doc/mosfet_bias_ex1_mirror.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 mosfet_bias_ex1_mirror.pdf mosfet_bias_ex1_mirror.svg
convert -flatten -density 600 mosfet_bias_ex1_mirror.pdf mosfet_bias_ex1_mirror.png

cp mosfet_bias_ex1_mirror.png ../doc/mosfet_bias_ex1_mirror.png
cp mosfet_bias_ex1_mirror.svg ../doc/mosfet_bias_ex1_mirror.svg

rm -rf mosfet_bias_ex1_mirror*

cd ../python/circuits/
python3 mosfet_bias_ex1.py
cd -

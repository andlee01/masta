#!/bin/sh

echo Generating Figures...

pdflatex ../doc/pmos_active_load.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 pmos_active_load.pdf pmos_active_load.svg
convert -flatten -density 600 pmos_active_load.pdf pmos_active_load.png

cp pmos_active_load.png ../doc/pmos_active_load.png
cp pmos_active_load.svg ../doc/pmos_active_load.svg

rm -rf pmos_active_load*

cd ../python/circuits/
python3 mosfet_pmos_active_load.py
cd -

#!/bin/sh

echo Generating Figures...

pdflatex ../doc/pmos_active_load_small_signal.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 pmos_active_load_small_signal.pdf pmos_active_load_small_signal.svg
convert -flatten -density 600 pmos_active_load_small_signal.pdf pmos_active_load_small_signal.png

cp pmos_active_load_small_signal.png ../doc/pmos_active_load_small_signal.png
cp pmos_active_load_small_signal.svg ../doc/pmos_active_load_small_signal.svg

rm -rf pmos_active_load_small_signal*

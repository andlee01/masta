#!/bin/sh

echo Generating Figures...

pdflatex ../doc/common_source_amp.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 common_source_amp.pdf common_source_amp.svg
convert -flatten -density 600 common_source_amp.pdf common_source_amp.png

cp common_source_amp.png ../doc/common_source_amp.png
cp common_source_amp.svg ../doc/common_source_amp.svg

rm -rf common_source_amp*

pdflatex ../doc/common_source_amp_small_signal.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 common_source_amp_small_signal.pdf common_source_amp_small_signal.svg
convert -flatten -density 600 common_source_amp_small_signal.pdf common_source_amp_small_signal.png

cp common_source_amp_small_signal.png ../doc/common_source_amp_small_signal.png
cp common_source_amp_small_signal.svg ../doc/common_source_amp_small_signal.svg

rm -rf common_source_amp_small_signal*

cd ../python/circuits/
python3 common_source_amp_tb.py
cd -

#!/bin/sh

echo Generating Figures...

pdflatex ../doc/beta_multiplier.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 beta_multiplier.pdf beta_multiplier.svg
convert -flatten -density 600 beta_multiplier.pdf beta_multiplier.png

cp beta_multiplier.png ../doc/beta_multiplier/beta_multiplier.png
cp beta_multiplier.svg ../doc/beta_multiplier/beta_multiplier.svg

pdflatex ../doc/beta_multiplier_small_signal.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 beta_multiplier_small_signal.pdf beta_multiplier_small_signal.svg
convert -flatten -density 600 beta_multiplier_small_signal.pdf beta_multiplier_small_signal.png

cp beta_multiplier_small_signal.png ../doc/beta_multiplier/beta_multiplier_small_signal.png
cp beta_multiplier_small_signal.svg ../doc/beta_multiplier/beta_multiplier_small_signal.svg

pdflatex ../doc/beta_multiplier_loop_gain_small_signal.tex > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 beta_multiplier_loop_gain_small_signal.pdf beta_multiplier_loop_gain_small_signal.svg
convert -flatten -density 600 beta_multiplier_loop_gain_small_signal.pdf beta_multiplier_loop_gain_small_signal.png

cp beta_multiplier_loop_gain_small_signal.png ../doc/beta_multiplier/beta_multiplier_loop_gain_small_signal.png
cp beta_multiplier_loop_gain_small_signal.svg ../doc/beta_multiplier/beta_multiplier_loop_gain_small_signal.svg

pdflatex "\def\EnableCstray{}\input{../doc/beta_multiplier_loop_gain_small_signal.tex}"  > /dev/null

echo Converting to SVG/PNG...
convert -flatten -density 600 beta_multiplier_loop_gain_small_signal.pdf beta_multiplier_loop_gain_small_signal_cstray.svg
convert -flatten -density 600 beta_multiplier_loop_gain_small_signal.pdf beta_multiplier_loop_gain_small_signal_cstray.png

cp beta_multiplier_loop_gain_small_signal_cstray.png ../doc/beta_multiplier/beta_multiplier_loop_gain_small_signal_cstray.png
cp beta_multiplier_loop_gain_small_signal_cstray.svg ../doc/beta_multiplier/beta_multiplier_loop_gain_small_signal_cstray.svg


rm -rf beta_multiplier*
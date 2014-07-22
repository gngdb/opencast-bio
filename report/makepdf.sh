#!/bin/bash
pdflatex -output-directory out/ main.tex 
evince out/main.pdf &

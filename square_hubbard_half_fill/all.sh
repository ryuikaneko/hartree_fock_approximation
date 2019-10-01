#!/bin/bash

gcc -O3 hf.c -llapack -lblas -lm -o hf.out -Wall -Wextra

./make_dat_n0.5_af.sh
./make_datene_n0.5_af.sh
./make_datgap_n0.5_af.sh
./make_datmag_n0.5_af.sh

gnuplot gnuplot_make_fig_1
gnuplot gnuplot_make_fig_2
gnuplot gnuplot_make_fig_3

rm hf.out

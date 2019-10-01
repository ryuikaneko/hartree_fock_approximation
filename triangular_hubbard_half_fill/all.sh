#!/bin/bash

gcc -O3 hf.c -llapack -lblas -lm -o hf.out -Wall -Wextra

./hf.out -U 10 -k 6 > dat_N6x6

./hf.out -U 10 -k 300 > dat_N300x300

rm ./hf.out

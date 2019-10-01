#!/bin/bash

n=0.5

for U in \
`seq -f %.1f 0.0 0.5 10.01`
do
  echo ${U}
  ./hf.out -U ${U} -k 240 -i 10000 -m 0.9 -f ${n} -c 1 -s 100 -b > dat_n${n}_af_U${U}
done

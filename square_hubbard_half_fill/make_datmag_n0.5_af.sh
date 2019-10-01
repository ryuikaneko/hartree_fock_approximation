#!/bin/bash

n=0.5

catdat=catdatmag_n${n}_af
echo -ne "# U SzA SzB NA NB\n" > ${catdat}

for U in \
`seq -f %.1f 0.0 0.5 10.01`
do
  echo ${U}
  dat=dat_n${n}_af_U${U}
  echo -ne "${U} " >> ${catdat}
  grep "# SzA SzB NA NB" -A1 ${dat} | tail -n 1 >> ${catdat}
done

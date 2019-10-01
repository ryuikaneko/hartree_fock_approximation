#!/bin/bash

n=0.5

catdat=catdatene_n${n}_af
echo -ne "# U iter Ene Ene_0 Ene_U ave(n) n_up n_dn ave(deln) deln_up deln_dn\n" > ${catdat}

for U in \
`seq -f %.1f 0.0 0.5 10.01`
do
  echo ${U}
  dat=dat_n${n}_af_U${U}
  echo -ne "${U} " >> ${catdat}
  grep "# converged" -A1 ${dat} | tail -n 1 >> ${catdat}
done

#!/bin/bash

n=0.5

catdat=catdatgap_n${n}_af
echo -ne "# U charge_gap\n" > ${catdat}

for U in \
`seq -f %.1f 0.0 0.5 10.01`
do
  echo ${U}
  dat=dat_n${n}_af_U${U}
  echo -ne "${U} " >> ${catdat}
  grep Egap ${dat} | sed 's/.*=//g' >> ${catdat}
done

#!/bin/bash -l

PW=pw.x
PH=ph.x
P2Y=p2y
YAMBO=yambo
LELPH=../../../../src/lelphc
NCORES=1
#run scf
cd scf
mpirun -np $NCORES $PW < scf.in | tee scf.out
cd ..
# run phonons
cd ph
cp -r ../scf/mos2.* .
mpirun -np $NCORES $PH < ph.in | tee ph.out
cd ..
# run nscf
cd nscf
cp -r ../scf/mos2.* .
mpirun -np $NCORES $PW < nscf.in | tee nscf.out
cd mos2.save
## create save dir
$P2Y
$YAMBO
cd ../..
######### lets start LetzElPhC 
## first create ph_save in ph directory
cd ph
$LELPH -pp --code=qe -F ph.in
cd ..
## now compute the elph matrix elements
cd elph
$LELPH -F elph.in
cd ..
cd bse
mpirun -np $NCORES $YAMBO -F bse.in -J BSE -C BSE -I ../nscf/hBN.save/SAVE
cd ..

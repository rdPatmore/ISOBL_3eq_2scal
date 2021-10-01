#!/bin/bash

JOBNO=516

# move to run directory
cd ../run

# make target directory
HOMEDIR=$HOMEROOT/ISOBLR_${JOBNO}/run
ssh -t $HOMEHOST "mkdir -p $HOMEDIR"

# rsync peripheral files
rsync -avzL hFac* $HOMEHOST:$HOMEDIR
rsync -avzL Depth* $HOMEHOST:$HOMEDIR
rsync -avzL DRF* $HOMEHOST:$HOMEDIR
rsync -avzL DRC* $HOMEHOST:$HOMEDIR
rsync -avzL DXC* $HOMEHOST:$HOMEDIR
rsync -avzL DYC* $HOMEHOST:$HOMEDIR
rsync -avzL DXG* $HOMEHOST:$HOMEDIR
rsync -avzL DYG* $HOMEHOST:$HOMEDIR
rsync -avzL RAC* $HOMEHOST:$HOMEDIR
rsync -avzL RC* $HOMEHOST:$HOMEDIR
rsync -avzL XC* $HOMEHOST:$HOMEDIR
rsync -avzL YC* $HOMEHOST:$HOMEDIR

# for each diagnostic file, make netcdf and rsync
VARS="
      vels 
      scalars 
      state2D 
      diffvisc
      uStar
      SHIDrag
     "

for VAR in $VARS
do
  rm -rf $VAR.nc
  echo 'seconds since 1971-01-01 00:00:00' > file_list
  ls $VAR.*.data >> file_list
  python $M_ROOT/utilities/mit2nc/mit2nc.py
  rsync -avzL $VAR.nc $HOMEHOST:$HOMEDIR
  rm -rf $VAR.nc
done

#! /bin/bash
infile=$HOME/Models/icon/build_lndonly/experiments/jsbalone_R2B4_sfa/jsbalone_R2B4_sfa_lnd_basic_ml_20000201T000000Z.nc
outfile=jsbalone_latlon.nc

cdo -r remapdis,t63grid ${infile} ${outfile}

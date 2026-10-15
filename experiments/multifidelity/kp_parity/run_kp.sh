#!/bin/bash
# usage: run_kp.sh TAG RUN_ID N DT_PY FIN [edit_dt_line(0/1) DT_M]
set -e
TAG=$1; RUN=$2; N=$3; DTPY=$4; FIN=$5; EDITDT=${6:-0}; DTM=${7:-}
KP=/oscar/scratch/eaguerov/kp_parity
REPO=/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
SRC=/oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/postprocessing/bio_mix_elvis.m
D=$KP/wd/D:/1_bioreactor/Data_save/bio_rec_new6_save/$TAG
mkdir -p $D/${TAG}_save
cp $SRC $KP/bio_mix_$TAG.m
sed -i "5s/.*/NAME    = '$TAG';/; 9s/.*/N       = $N;/; 18s/.*/Fin_t   = $FIN;/" $KP/bio_mix_$TAG.m
if [ "$EDITDT" = 1 ]; then sed -i "16s/.*/dt_File = $DTM;/" $KP/bio_mix_$TAG.m; fi
diff $SRC $KP/bio_mix_$TAG.m > $KP/edit_$TAG.diff || true
cd $REPO
uv run python scripts/kim_convert.py $REPO/runs/$RUN/fields $D/${TAG}_save $N 48.5863329 $DTPY $FIN 2>&1 | grep -v warning
cd $KP/wd
export PATH=/oscar/rt/sw/external/matlab/R2025a/bin:$PATH
cp $KP/bio_mix_$TAG.m ./bio_mix_$TAG.m
matlab -batch "bio_mix_$TAG" -singleCompThread > $KP/matlab_$TAG.log 2>&1
cp $D/out_mix_all.mat $KP/out_mix_all_$TAG.mat
rm -rf $D/${TAG}_save
echo done $TAG

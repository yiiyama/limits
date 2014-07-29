#!/bin/bash

cd $HOME/cmssw/Combine612
eval `scram runtime -sh`

cd $HOME/src/GammaL/limits

SUFFIX=_mt100

#for mglu in $(seq 400 50 1300); do
for mglu in $(seq 400 50 800); do
    [ $mglu -ge 800 -a $mglu -lt 1000 ] && continue
    for mchi in $(seq 25 50 $mglu); do
        echo $mglu $mchi
        python writeDataCard.py T5wg_${mglu}_${mchi} $SUFFIX || exit 1
        python computeLimits.py $HOME/work/datacards${SUFFIX}/T5wg_${mglu}_${mchi}.dat $HOME/work/logs/T5wg_${mglu}_${mchi}.log $HOME/work/limits${SUFFIX} || exit 1
        rm -r /tmp/yiiyama/T5wg_${mglu}_${mchi}
    done
done
       
hadd -f $HOME/output/GammaL/limits/T5wg${SUFFIX}.root $HOME/work/limits${SUFFIX}/T5wg*.root

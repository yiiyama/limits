#!/bin/bash

GRID=true

MODEL=$1
POINT=$2
SUFFIX=$3

DEST=$HOME/work/limits${SUFFIX}
PKLDIR=$HOME/output/GammaL/limits
RESPKL=result.pkl

source /afs/cern.ch/cms/cmsset_default.sh

cd $HOME/cmssw/Combine612
eval `scram runtime -sh`

cd $HOME/src/GammaL/limits

if $GRID; then
    python computeLimits.py -g $DEST/grid_${MODEL}_${POINT} $MODEL $POINT $RESPKL $PKLDIR $DEST
else
    python computeLimits.py $MODEL $POINT $RESPKL $PKLDIR $DEST
fi

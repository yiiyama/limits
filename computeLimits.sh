#!/bin/bash

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

python computeLimits.py $MODEL $POINT $RESPKL $PKLDIR $DEST

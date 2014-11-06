#!/bin/bash

MODEL=$1
POINT=$2
SUFFIX=$3

SOURCE=$HOME/work/datacards${SUFFIX}
DEST=$HOME/work/limits${SUFFIX}
PKL=$HOME/output/GammaL/limits/full.pkl

LOG=$HOME/work/logs/${MODEL}_${POINT}${SUFFIX}.log

source /afs/cern.ch/cms/cmsset_default.sh

cd $HOME/cmssw/Combine612
eval `scram runtime -sh`

cd $HOME/src/GammaL/limits

echo "" >> $LOG
date >> $LOG
echo "" >> $LOG

python computeLimits.py $MODEL $POINT $SOURCE $PKL $DEST >> $LOG 2>&1

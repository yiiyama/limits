#!/bin/bash

MODEL=$1
POINT=$2

source /afs/cern.ch/cms/cmsset_default.sh

cd $HOME/cmssw/SLC6Ntuplizer5314
eval `scram runtime -sh`

cd $HOME/src/GammaL/plotstack/GammaL

python countSignal.py -i rooth://ncmu40//store/countSignal/trees -p ${MODEL} -m ${POINT} > $HOME/work/logs/${MODEL}_${POINT}.log 2>&1

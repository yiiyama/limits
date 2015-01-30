#!/bin/bash

MODEL=$1
POINT=$2
INDEX=$3

DATADIR=$HOME/output/GammaL/limits
OUTPUTDIR=$HOME/work/limits_full
LOGDIR=$HOME/work/logs

CWD=$HOME/src/GammaL/limits

NSTEPS=100

if [ $INDEX ]; then
    source /afs/cern.ch/cms/cmsset_default.sh
    
    cd $HOME/cmssw/Combine612
    eval `scram runtime -sh`
    
    cd $CWD
    
    python clsgrid.py $MODEL $POINT $INDEX $DATADIR $OUTPUTDIR
else
    POINTS=""
    if [ $POINT ]; then
        POINTS=$POINT
    else
        if [ $MODEL = "TChiwg" ]; then
            for MCHI in $(seq 200 10 800); do
                POINTS="${POINTS}$MCHI "
            done
        elif [ $MODEL = "T5wg" ]; then
            for MGLU in $(seq 700 50 1500); do
                for MCHI in $(seq 25 50 $MGLU); do
                    POINTS="${POINTS}${MGLU}_${MCHI} "
                done
            done
        elif [ $MODEL = "Spectra_gW" ]; then
            for M3 in $(seq 715 50 1515); do
                for M2 in $(seq 205 50 $M3); do
                    POINTS="${POINTS}M3_${M3}_M2_${M2} "
                done
            done
        else
            exit 1
        fi
    fi

    for POINT in $POINTS; do
        POINTNAME=${MODEL}_${POINT}

        mkdir $LOGDIR/$POINTNAME 2> /dev/null

        for INDEX in $(seq 0 $(($NSTEPS-1))); do
            JOB=${POINTNAME}_${INDEX}

            [ $(ls $OUTPUTDIR/grid_$POINTNAME/higgsCombineGrid${INDEX}.* 2>/dev/null | wc -l) -ne 0 ] && continue

            echo bsub -q 8nh -J $JOB -o $LOGDIR/$POINTNAME/$INDEX.log "$CWD/clsgrid.sh $MODEL $POINT $INDEX"
            bsub -q 8nh -J $JOB -o $LOGDIR/$POINTNAME/$INDEX.log "$CWD/clsgrid.sh $MODEL $POINT $INDEX"
        done
    done
fi

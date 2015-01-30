import os
import pickle
import shutil
import array
import sys
import math
import subprocess
import ROOT

import datacard

SETENV = 'cd /afs/cern.ch/user/y/yiiyama/cmssw/Combine612; eval `scram runtime -sh`;'
NSTEPS = 100

def makeGrid(model, point, index, datadir, outputdir):

    pointName = model + '_' + point

    # get rvalue range from {model}.root

    print 'Calculating rvalue for index', index

    limitFile = ROOT.TFile.Open(datadir + '/' + model + '.root')
    limitTree = limitFile.Get('limitTree')

    vPointName = array.array('c', ' ' * 100)
    vLimM2s = array.array('d', [0.])
    vLimP2s = array.array('d', [0.])
            
    limitTree.SetBranchAddress('pointName', vPointName)
    limitTree.SetBranchAddress('m2s', vLimM2s)
    limitTree.SetBranchAddress('p2s', vLimP2s)

    iEntry = 0
    while limitTree.GetEntry(iEntry) > 0:
        iEntry += 1
        
        if vPointName.tostring().strip().strip('\0') == pointName:
            lnrlow = math.log(vLimM2s[0] * 0.7)
            lnrhigh = math.log(vLimP2s[0] * 1.5)
            break
    else:
        raise RuntimeError('Bounds for ' + pointName + ' not found')

    limitFile.Close()

    rval = math.exp(lnrlow + (lnrhigh - lnrlow) / NSTEPS * index)

    # prepare data card from result.pkl (observation & exp errors) and {model}.pkl (signal expectation)

    print 'Preparing data card'

    with open(datadir + '/result.pkl') as source:
        channels = pickle.load(source)

    with open(datadir + '/' + model + '.pkl') as source:
        signals = pickle.load(source)

    signalData = signals[pointName][0]
    for name, channel in channels.items():
        channel.processes['signal'] = signalData[name]

    workdir = os.environ['TMPDIR'] + '/' + pointName
    os.makedirs(workdir)

    cardPath = workdir + '/' + pointName + '.dat'

    datacard.writeDataCard(channels, cardPath)

    # run combine

    suffix = 'Grid' + str(index)
    command = 'combine ' + cardPath + ' -n ' + suffix + ' -M HybridNew -s -1 --freq --clsAcc 0 -T 500 -i 4 --fork 4 --saveToys --saveHybridResult --singlePoint ' + str(rval)

    print command

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; ' + command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    for line in proc.communicate()[0].split('\n'):
        print line
    
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    for fileName in os.listdir(workdir):
        if fileName.startswith('higgsCombine' + suffix + '.HybridNew') and fileName.endswith('.root'):

            sys.stderr.flush()
            sys.stderr.close()
            sys.stderr = sys.stdout
            os.close(2)

            err = os.open(workdir + '/errors.txt', os.O_WRONLY | os.O_CREAT)

            source = ROOT.TFile.Open(workdir + '/' + fileName)
            source.Get('toys')
            source.Get('ProcessID0')
            source.Get('limit')
            source.Close()

            os.close(err)

            with open(workdir + '/errors.txt') as errors:
                for line in errors:
                    if 'Error' in line:
                        raise RuntimeError('Error found in output')

            shutil.copy(workdir + '/' + fileName, outputdir + '/' + fileName)

if __name__ == '__main__':

    model = sys.argv[1]
    point = sys.argv[2]
    index = int(sys.argv[3])
    datadir = sys.argv[4]
    outputdir = sys.argv[5]

    makeGrid(model, point, index, datadir, outputdir + '/grid_' + model + '_' + point)

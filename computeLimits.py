import os
import sys
import shutil
import time
import math
import subprocess
import array
import re
import pickle
import ROOT

ROOT.gROOT.SetBatch(True)

import datacard

SETENV = 'cd /afs/cern.ch/user/y/yiiyama/cmssw/Combine612; eval `scram runtime -sh`;'
XSECDIR = '/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs'
FULLCLS = False
FORCEPROF = False

def getLimits(fileName, vLimits, calculate = False):
    """
    Extract median and +-1/2 sigma expected limits from a 6-entry tree output from combine. If calculate = True, tree is assumed to contain raw data (one toy per entry).
    """

    resFile = ROOT.TFile.Open(fileName)
    resTree = resFile.Get('limit')

    resTree.SetEstimate(resTree.GetEntries() + 1)
    nEntries = resTree.Draw('quantileExpected:limit:iToy', 'limit > 0.', 'goff') # condition to avoid nans
    quantiles = resTree.GetV1()
    rvalues = resTree.GetV2()
    toyIndices = resTree.GetV3()

    if calculate:
        rArr = [rvalues[i] for i in range(nEntries) if toyIndices[i] > 0]

        rArr.sort()

        # remove outliers
        binwidth = (rArr[int(0.5 * nEntries)] - rArr[0]) / 5.
        hist = ROOT.TH1D('hist', 'r value distribution', int((rArr[-1] - rArr[0] + 1.) / binwidth) + 1, rArr[0], rArr[-1] + 1.)
        for r in rArr: hist.Fill(r)

        for iX in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinContent(iX) == 0. and hist.GetBinContent(iX + 1) == 0.: break

        while rArr[-1] > hist.GetXaxis().GetBinLowEdge(iX): rArr.pop()

        n = len(rArr)

        vLimits['m2s'][0] = rArr[int(0.05 * n)]
        vLimits['m1s'][0] = rArr[int(0.32 * n)]
        vLimits['med'][0] = rArr[int(0.5 * n)]
        vLimits['p1s'][0] = rArr[int(0.68 * n)]
        vLimits['p2s'][0] = rArr[int(0.95 * n)]

        try:
            obsIndex = next(i for i in range(nEntries) if toyIndices[i] == 0)
            vLimits['obs'][0] = rvalues[obsIndex]
        except:
            return False

        return True

    else:
        limSet = set()
        
        for iEntry in range(nEntries):
            if quantiles[iEntry] < 0.: quant = 'obs'
            elif quantiles[iEntry] < 0.03: quant = 'm2s'
            elif quantiles[iEntry] < 0.2: quant = 'm1s'
            elif quantiles[iEntry] < 0.6: quant = 'med'
            elif quantiles[iEntry] < 0.9: quant = 'p1s'
            else: quant = 'p2s'

            limSet.add(quant)
            
            vLimits[quant][0] = rvalues[iEntry]

        return len(quant) == 6

        
def writeLog(header, content = ''):
    sys.stdout.write('-' * 20)
    sys.stdout.write(' ' + header + ' ')
    sys.stdout.write('-' * 20)
    sys.stdout.write('\n')
    sys.stdout.write(content + '\n')


def asymptotic(card, workdir, vLimits):
    """
    Run combine in asymptotic mode and return the resulting expected limits in the form of python dict (using getLimits)
    """

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -M Asymptotic -s -1 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = proc.communicate()[0]

    writeLog('Asymptotic', output)

    # read the last line of r value calculation
    rLine = ''
    for line in output.split('\n'):
        if re.match('^At r = ', line.strip()):
            rLine = line.strip()
        elif rLine:
            if 'q_mu = nan' in rLine:
                rLine = ''

    if not rLine:
        return {}

    for fileName in os.listdir(workdir):
        if 'higgsCombineTest.Asymptotic.' in fileName: break
    else:
        raise RuntimeError('No asymptotic result found')

    return getLimits(workdir + '/' + fileName, vLimits)


def profileLikelihood(card, workdir, vLimits):
    """
    Run combine in profile likelihood mode and return the resulting expected limits in the form of python dict (using getLimits)
    """

    writeLog('ProfileLikelihood Expected')

    procs = []
    seeds = [str(i) for i in range(1, 9)]
    while True:
        for seed, proc in procs:
            if proc.poll() is None: continue

            output = proc.communicate()[0]

            writeLog('Proc ' + seed, output)

            procs.remove((seed, proc))
        
        if len(procs) == 0 and len(seeds) == 0: break

        if len(procs) >= 4: continue

        try:
            seed = seeds.pop()
        except IndexError:
            continue

        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -M ProfileLikelihood -t 20 -s ' + seed + ' 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

        procs.append((seed, proc))

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Obs -M ProfileLikelihood -s -1 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = proc.communicate()[0]

    writeLog('ProfileLikelihood Observed', output)

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd ProfileLikelihood.root higgsCombine*.ProfileLikelihood.*.root 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = proc.communicate()[0]

    writeLog('Merge ProfileLikelihood', output)

    return getLimits(workdir + '/ProfileLikelihood.root', vLimits, calculate = True)


def makeGrid(bounds, card, workdir):
    """
    Generate toys for 100 r-values (signal strengths; mu-value) between the given bounds. To be used for frequentist expected limits. Resulting ROOT files are merged with hadd into [workdir]/HybridGrid.root
    """

    nSteps = 100
    rvalues = [bounds[0] + (bounds[1] - bounds[0]) / nSteps * i for i in range(nSteps)]
    procs = []
    while True:
        for rval, proc in procs:
            if proc.poll() is None: continue

            output = proc.communicate()[0]

            writeLog('Grid ' + rval, output)

            procs.remove((rval, proc))

        if len(procs) == 0 and len(rvalues) == 0: break

        if len(procs) >= 4: continue

        try:
            rval = '%.4f' % rvalues.pop()
        except IndexError:
            continue

        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Grid -M HybridNew -s -1 --freq --clsAcc 0 -T 1000 -i 1 --saveToys --saveHybridResult --singlePoint ' + rval, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

        procs.append((rval, proc))

    writeLog('Merge Grid')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridGrid.root higgsCombineGrid.HybridNew.*.root 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    proc.communicate()


def fullCLs(card, workdir):
    """
    Run combine in HybridNew mode and return the resulting expected limits in the form of python dict (using getLimits)    
    """

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n FullCLs -M HybridNew -s -1 --freq --grid HybridGrid.root 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = proc.communicate()[0]

    writeLog('Observed', output)

    for quant in ['0.025', '0.16', '0.5', '0.84', '0.975']:
        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n FullCLs -M HybridNew -s -1 --freq --grid HybridGrid.root --expectedFromGrid ' + quant + ' 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        output = proc.communicate()[0]
        
        writeLog('Expected ' + quant, output)

    writeLog('Merge FullCLs')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridFullCLs.root higgsCombineFullCLs.HybridNew.*.root 2>&1', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    proc.communicate()

    return getLimits(workdir + '/HybridFullCLs.root', vLimits)


def computeLimits(model, point, channels, outputdir):

    pointName = model + '_' + point

    try:
        workdir = os.environ['TMPDIR'] + '/' + pointName
    except KeyError:
        workdir = '/afs/cern.ch/user/y/yiiyama/work/tmp/' + pointName

    try:
        shutil.rmtree(workdir)
    except OSError:
        pass
    try:
        os.makedirs(workdir)
    except OSError:
        pass

    # Results are written into a tree (to be merged with calculations for all other signal points)

    limitPoints = ['obs', 'med', 'm2s', 'm1s', 'p1s', 'p2s']

    outputFile = ROOT.TFile.Open(workdir + '/' + pointName + '.root', 'recreate')
    limitTree = ROOT.TTree('limitTree', 'Limit Tree')

    vPointName = array.array('c', pointName + '\0')
    vMethod = array.array('c', '\0' * 64)
    vLimits = dict([(p, array.array('d', [0.])) for p in limitPoints])

    limitTree.Branch('pointName', vPointName, 'pointName/C')
    limitTree.Branch('method', vMethod, 'method/C')
    for p in limitPoints:
        limitTree.Branch(p, vLimits[p], p + '/D')

    cardPath = workdir + '/' + pointName + '.dat'

    datacard.writeDataCard(channels, cardPath)

    writeLog('Calculating asymptotic limits')
    method = 'asymptotic'

    converged = asymptotic(cardPath, workdir, vLimits)

    if converged:
        for iC in range(len(method)):
            vMethod[iC] = method[iC]
            vMethod[iC + 1] = '\0'
            
        limitTree.Fill()
    else:
        print 'Asymptotic method did not converge.'

    if FORCEPROF or not converged:
        writeLog('Using profile likelihood')
        method = 'profileLikelihood'

        if profileLikelihood(cardPath, workdir, vLimits):
            for iC in range(len(method)):
                vMethod[iC] = method[iC]
                vMethod[iC + 1] = '\0'
            
            limitTree.Fill()

    if FULLCLS:
        if limitTree.GetEntries() == 0:
            writeLog('No estimate of bounds for grid production.')
            sys.exit(1)
        
        writeLog('Calculating full CLs')
        method = 'fullCLs'

        makeGrid((vLimits['m2s'][0], vLimits['p2s'][0]), cardPath, workdir)

        if fullCLs(cardPath, workdir, vLimits)
            for iC in range(len(method)):
                vMethod[iC] = method[iC]
                vMethod[iC + 1] = '\0'
        
            limitTree.Fill()

    if limitTree.GetEntries() == 0:
        writeLog('Failed to calculate limits.')
        sys.exit(1)

    limitTree.Scan('*')

    limitTreeFile.Write()
    limitTreeFile.Close()

    shutil.copyfile(workdir + '/' + pointName + '.root', outputdir + '/' + pointName + '.root')


if __name__ == '__main__':

    model = sys.argv[1]
    point = sys.argv[2]
    result = sys.argv[3]
    pkldir = os.path.realpath(sys.argv[4])
    outputdir = sys.argv[5]

    with open(pkldir + '/' + result) as source:
        channels = pickle.load(source)

    with open(pkldir + '/' + model + '.pkl') as source:
        signals = pickle.load(source)

    for name, channel in channels.items():
        channel.processes['signal'] = signals[model + '_' + point][name]

    computeLimits(model, point, channels, outputdir)

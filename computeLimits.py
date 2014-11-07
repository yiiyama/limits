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

from datacard import *

SETENV = 'cd /afs/cern.ch/user/y/yiiyama/cmssw/Combine612; eval `scram runtime -sh`;'
FULLCLS = False

def getLimits(fileName, calculate = False):
    """
    Extract median and +-1/2 sigma expected limits from a 6-entry tree output from combine. If calculate = True, tree is assumed to contain raw data (one toy per entry).
    """

    resFile = ROOT.TFile.Open(fileName)
    resTree = resFile.Get('limit')

    resTree.SetEstimate(resTree.GetEntries() + 1)
    nEntries = resTree.Draw('quantileExpected:limit', 'limit > 0.', 'goff') # condition to avoid nans
    quantiles = resTree.GetV1()
    rvalues = resTree.GetV2()

    limits = {}

    if calculate:
        rArr = [rvalues[i] for i in range(nEntries)]

        rArr.sort()

        # remove outliers
        binwidth = (rArr[int(0.5 * nEntries)] - rArr[0]) / 5.
        hist = ROOT.TH1D('hist', 'r value distribution', int((rArr[-1] - rArr[0] + 1.) / binwidth) + 1, rArr[0], rArr[-1] + 1.)
        for r in rArr: hist.Fill(r)

        for iX in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinContent(iX) == 0. and hist.GetBinContent(iX + 1) == 0.: break

        while rArr[-1] > hist.GetXaxis().GetBinLowEdge(iX): rArr.pop()

        n = len(rArr)

        limits['m2s'] = rArr[int(0.05 * n)]
        limits['m1s'] = rArr[int(0.32 * n)]
        limits['med'] = rArr[int(0.5 * n)]
        limits['p1s'] = rArr[int(0.68 * n)]
        limits['p2s'] = rArr[int(0.95 * n)]
        limits['obs'] = 0.

    else:
        for iEntry in range(nEntries):
            if quantiles[iEntry] < 0.: quant = 'obs'
            elif quantiles[iEntry] < 0.03: quant = 'm2s'
            elif quantiles[iEntry] < 0.2: quant = 'm1s'
            elif quantiles[iEntry] < 0.6: quant = 'med'
            elif quantiles[iEntry] < 0.9: quant = 'p1s'
            else: quant = 'p2s'
    
            limits[quant] = rvalues[iEntry]

    return limits


def writeLog(header, content = ''):
    sys.stdout.write('-' * 20)
    sys.stdout.write(' ' + header + ' ')
    sys.stdout.write('-' * 20)
    sys.stdout.write('\n')
    sys.stdout.write(content)


def asymptotic(card, workdir):
    """
    Run combine in asymptotic mode and return the resulting expected limits in the form of python dict (using getLimits)
    """

    writeLog('Asymptotic')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -M Asymptotic -s -1 2>&1', shell = True, stdout = subprocess.PIPE)
    output = ''
    while proc.poll() is None:
        output += proc.communicate()[0]
        time.sleep(1)

    try:
        output += proc.communicate()[0]
    except:
        pass

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

    return getLimits(workdir + '/' + fileName)


def profileLikelihood(card, workdir):
    """
    Run combine in profile likelihood mode and return the resulting expected limits in the form of python dict (using getLimits)
    """

    writeLog('ProfileLikelihood Expected')

    outputs = {}
    procs = []
    seeds = [str(i) for i in range(1, 9)]
    while True:
        while len(procs):
            for seed, proc in procs:
                try:
                    outputs[seed] += proc.communicate()[0]
                except:
                    pass

            try:
                seed, proc = next(x for x in procs if x[1].poll() is not None)

                writeLog('Proc ' + seed, outputs[seed])

                procs.remove((seed, proc))

            except StopIteration:
                break

        else:
            if len(seeds) == 0: break
        
        if len(procs) >= 4 or len(seeds) == 0:
            time.sleep(10)
            continue

        seed = str(seeds.pop())

        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -M ProfileLikelihood -t 20 -s ' + seed + ' 2>&1', shell = True, stdout = subprocess.PIPE)

        procs.append((seed, proc))
        outputs[seed] = ''

    writeLog('Merge ProfileLikelihood')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd ProfileLikelihood.root higgsCombineTest.ProfileLikelihood.*.root 2>&1', shell = True)
    while proc.poll() is None:
        time.sleep(1)

    limits = getLimits(workdir + '/ProfileLikelihood.root', calculate = True)

    writeLog('ProfileLikelihood Observed')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Obs -M ProfileLikelihood -s -1 2>&1', shell = True, stdout = subprocess.PIPE)
    output = ''
    while proc.poll() is None:
        output += proc.communicate()[0]
        time.sleep(1)

    for fileName in os.listdir(workdir):
        if 'higgsCombineObs.ProfileLikelihood.' in fileName: break
    else:
        raise RuntimeError('No asymptotic result found')

    obsLim = getLimits(workdir + '/' + fileName)

    limits.update(obsLim)

    return limits


def makeGrid(bounds, card, workdir):
    """
    Generate toys for 100 r-values (signal strengths; mu-value) between the given bounds. To be used for frequentist expected limits. Resulting ROOT files are merged with hadd into [workdir]/HybridGrid.root
    """

    nSteps = 100
    rvalues = [bounds[0] + (bounds[1] - bounds[0]) / nSteps * i for i in range(nSteps)]
    procs = []
    outputs = {}
    while True:
        while len(procs):
            for rval, proc in procs:
                try:
                    outputs[rval] += proc.communicate()[0]
                except:
                    pass

            try:
                rval, proc = next(x for x in procs if x[1].poll() is not None)

                writeLog('Grid ' + rval, outputs[rval])

                procs.remove((rval, proc))

            except StopIteration:
                break

        else:
            if len(rvalues) == 0: break
        
        if len(procs) >= 4 or len(rvalues) == 0:
            time.sleep(10)
            continue

        rval = '%.4f' % rvalues.pop()

        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Grid -M HybridNew -s -1 --freq --clsAcc 0 -T 1000 -i 1 --saveToys --saveHybridResult --singlePoint ' + rval, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

        procs.append((rval, proc))
        outputs[rval] = ''

    writeLog('Merge Grid')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridGrid.root higgsCombineGrid.HybridNew.*.root 2>&1', shell = True)
    while proc.poll() is None:
        time.sleep(1)


def fullCLs(card, workdir):
    """
    Run combine in HybridNew mode and return the resulting expected limits in the form of python dict (using getLimits)    
    """

    writeLog('Observed')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n FullCLs -M HybridNew -s -1 --freq --grid HybridGrid.root 2>&1', shell = True)
    while proc.poll() is None:
        time.sleep(1)

    for quant in ['0.025', '0.16', '0.5', '0.84', '0.975']:
        writeLog('Expected ' + quant)

        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n FullCLs -M HybridNew -s -1 --freq --grid HybridGrid.root --expectedFromGrid ' + quant + ' 2>&1', shell = True)
        while proc.poll() is None:
            time.sleep(1)

    writeLog('Merge FullCLs')

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridFullCLs.root higgsCombineFullCLs.HybridNew.*.root 2>&1', shell = True)
    while proc.poll() is None:
        time.sleep(1)

    return getLimits(workdir + '/HybridFullCLs.root')


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

    cardPath = workdir + '/' + pointName + '.dat'

    writeDataCard(channels, cardPath)

    print 'Calculating asymptotic limits..'

    limits = asymptotic(cardPath, workdir)

    if len(limits) != 6:
        print 'Asymptotic method did not converge. Using profile likelihood..'

        limits = profileLikelihood(cardPath, workdir)

    if FULLCLS:
        print 'Calculating full CLs..'

        makeGrid((limits['m2s'], limits['p2s']), cardPath, workdir)
        limits = expected(cardPath, workdir)

    # Write result to a single-entry tree (to be merged with calculations for all other signal points)

    xsecs = {}
    nEvents = {}

    if model == 'T5wg+TChiwg':
        with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/T5wg.xsecs') as xsecsource:
            pname = 'T5wg_' + point
            for line in xsecsource:
                p, c, u, n = line.strip().split()
                if p != pname: continue

                xsecs[p] = (float(c), math.pow(float(u) * float(c), 2.))
                nEvents[p] = int(n)
                break
            else:
                raise RuntimeError('No point ' + pname + ' found')

        with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/TChiwg.xsecs') as xsecsource:
            pname = 'TChiwg_' + point[point.find('_') + 1:]
            for line in xsecsource:
                p, c, u, n = line.strip().split()
                if p != pname: continue

                xsecs[p] = (float(c) / 0.326, math.pow(float(u) * float(c) / 0.326, 2.)) # TChiwg_suppl generated with lepton filter
                nEvents[p] = int(n)
                break
            else:
                raise RuntimeError('No point ' + pname + ' found')
        
    elif model == 'Spectra_gW':
        with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/' + model + '.xsecs') as xsecsource:
            for line in xsecsource:
                p, c, u, n = line.strip().split()
                if pointName not in p: continue

                xsecs[p] = (float(c), math.pow(float(u) * float(c), 2.))
                nEvents[p] = int(n)

            if len(xsecs) == 0:
                raise RuntimeError('No point ' + pointName + ' found')

    else:
        with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/' + model + '.xsecs') as xsecsource:
            for line in xsecsource:
                p, c, u, n = line.strip().split()
                if p != pointName: continue

                xsecs[p] = (float(c), math.pow(float(u) * float(c), 2.))
                nEvents[p] = int(n)
                break
            else:
                raise RuntimeError('No point ' + pointName + ' found')
    
    outputFile = ROOT.TFile.Open(workdir + '/' + pointName + '.root', 'recreate')
    output = ROOT.TTree('limitTree', 'Limit Tree')

    vPointName = array.array('c', pointName + '\0')
    vXsec = array.array('d', [0.])
    vXsecErr = array.array('d', [0.])
    vLimObs = array.array('d', [0.])
    vLimMed = array.array('d', [0.])
    vLimM2s = array.array('d', [0.])
    vLimM1s = array.array('d', [0.])
    vLimP1s = array.array('d', [0.])
    vLimP2s = array.array('d', [0.])
    vNEvents = array.array('i', [0])
    vNEventsEff = array.array('d', [0.])
    vYield = array.array('i', [0] * len(channels))
    vYieldEff = array.array('d', [0.] * len(channels))

    channelNames = sorted(channels.keys())
                
    output.Branch('pointName', vPointName, 'pointName/C')
    output.Branch('xsec', vXsec, 'xsec/D')
    output.Branch('xsecErr', vXsecErr, 'xsecErr/D')
    output.Branch('limObs', vLimObs, 'limObs/D')
    output.Branch('limMed', vLimMed, 'limMed/D')
    output.Branch('limM2s', vLimM2s, 'limM2s/D')
    output.Branch('limM1s', vLimM1s, 'limM1s/D')
    output.Branch('limP1s', vLimP1s, 'limP1s/D')
    output.Branch('limP2s', vLimP2s, 'limP2s/D')
    output.Branch('nEvents', vNEvents, 'nEvents/I')
    output.Branch('nEventsEff', vNEventsEff, 'nEventsEff/D')
    output.Branch('yield', vYield, channelNames[0] + '/I:' + ':'.join(channelNames[1:]))
    output.Branch('yieldEff', vYieldEff, channelNames[0] + '/D:' + ':'.join(channelNames[1:]))

    vXsec[0] = sum(x[0] for x in xsecs.values())
    vXsecErr[0] = math.sqrt(sum(x[1] for x in xsecs.values())) / vXsec[0]
    vNEvents[0] = sum(n for n in nEvents.values())
    sumNSigma = sum(nEvents[p] * xsecs[p][0] for p in nEvents.keys())
    sumNSigma2 = sum(nEvents[p] * xsecs[p][0] * xsecs[p][0] for p in nEvents.keys())
    vNEventsEff[0] = math.pow(sumNSigma, 2.) / sumNSigma2

    for iCh in range(len(channelNames)):
        ch = channelNames[iCh]
        signal = channels[ch].processes['signal']
        vYield[iCh] = sum(r[1] for r in signal.rates.values())
        sumNSigma = 0.
        sumNSigma2 = 0.
        for p, (r, c) in signal.rates.items():
            sumNSigma += c * xsecs[p][0]
            sumNSigma2 += c * xsecs[p][0] * xsecs[p][0]
            
        if sumNSigma2 > 0.:
            vYieldEff[iCh] = math.pow(sumNSigma, 2.) / sumNSigma2
        else:
            vYieldEff[iCh] = 0.

    if len(limits) == 6:
        vLimObs[0] = limits['obs']
        vLimMed[0] = limits['med']
        vLimM2s[0] = limits['m2s']
        vLimM1s[0] = limits['m1s']
        vLimP1s[0] = limits['p1s']
        vLimP2s[0] = limits['p2s']
    
        output.Fill()
    else:
        print 'Failed to calculate limits.'

    outputFile.Write()
    outputFile.Close()

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

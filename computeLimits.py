import os
import shutil
import time
import subprocess
import array
import ROOT

SETENV = 'cd /afs/cern.ch/user/y/yiiyama/cmssw/Combine612; eval `scram runtime -sh`;'

def getExpected(fileName):
    """
    Extract median and +-1/2 sigma expected limits from a 5-entry tree output from combine
    """

    resFile = ROOT.TFile.Open(fileName)
    resTree = resFile.Get('limit')
    resTree.SetEstimate(resTree.GetEntries() + 1)
    nEntries = resTree.Draw('quantileExpected:limit', '', 'goff')
    quantiles = resTree.GetV1()
    rvalues = resTree.GetV2()

    limits = {}
    for i in range(nEntries):
        if quantiles[i] < 0.03: quant = 'm2s'
        elif quantiles[i] < 0.2: quant = 'm1s'
        elif quantiles[i] < 0.6: quant = 'med'
        elif quantiles[i] < 0.9: quant = 'p1s'
        else: quant = 'p2s'

        limits[quant] = rvalues[i]

    return limits


def writeLog(fileName, header, content):
    with open(fileName, 'a') as logFile:
        logFile.write('-' * 20)
        logFile.write(' ' + header + ' ')
        logFile.write('-' * 20)
        logFile.write('\n')
        logFile.write(content)


def asymptotic(card, workdir, log):
    """
    Run combine in asymptotic mode and return the resulting expected limits in the form of python dict (using getExpected)
    """

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -M Asymptotic -s -1 --run expected', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    while proc.poll() is None:
        time.sleep(1)

    writeLog(log, 'Asymptotic', proc.communicate()[0])

    for fileName in os.listdir(workdir):
        if 'higgsCombineTest.Asymptotic.' in fileName: break
    else:
        raise RuntimeError('No asymptotic result found')

    return getExpected(workdir + '/' + fileName)


def makeGrid(bounds, card, workdir, log):
    """
    Generate toys for 100 r-values (signal strengths; mu-value) between the given bounds. To be used for frequentist expected limits. Resulting ROOT files are merged with hadd into [workdir]/HybridGrid.root
    """

    nSteps = 100
    rvalues = [bounds[0] + (bounds[1] - bounds[0]) / nSteps * i for i in range(nSteps)]
    procs = []
    while True:
        while len(procs):
            try:
                rval, proc = next(x for x in procs if x[1].poll() is not None)

                writeLog(log, 'Grid ' + rval, proc.communicate()[0])

                procs.remove((rval, proc))

            except StopIteration:
                break

        else:
            if len(rvalues) == 0: break
        
        if len(procs) >= 4 or len(rvalues) == 0:
            time.sleep(10)
            continue

        rval = '%.4f' % rvalues.pop()
        procs.append((rval, subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Grid -M HybridNew -s -1 --freq --clsAcc 0 -T 1000 -i 1 --saveToys --saveHybridResult --singlePoint ' + rval, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)))

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridGrid.root higgsCombineGrid.HybridNew.*.root', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    while proc.poll() is None:
        time.sleep(1)

    writeLog(log, 'Merge Grid', proc.communicate()[0])


def expected(card, workdir, log):
    """
    Run combine in HybridNew mode and return the resulting expected limits in the form of python dict (using getExpected)    
    """

    for quant in ['0.025', '0.16', '0.5', '0.84', '0.975']:
        proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; combine ' + card + ' -n Expected -M HybridNew -s -1 --freq --grid HybridGrid.root --expectedFromGrid ' + quant, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        while proc.poll() is None:
            time.sleep(1)

        writeLog(log, 'Expected ' + quant, proc.communicate()[0])

    proc = subprocess.Popen(SETENV + ' cd ' + workdir + '; hadd HybridExpected.root higgsCombineExpected.HybridNew.*.root', shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    while proc.poll() is None:
        time.sleep(1)

    writeLog(log, 'Merge Expected', proc.communicate()[0])

    return getExpected(workdir + '/HybridExpected.root')
    

if __name__ == '__main__':

    import sys

    card = os.path.realpath(sys.argv[1])
    log = os.path.realpath(sys.argv[2])
    outputdir = sys.argv[3]

    cardName = os.path.basename(card)
    pointName = cardName[0:cardName.rfind('.')]

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

    with open(log, 'w') as logFile:
        pass

    # Calculate asymptotic expected limits

    limits = asymptotic(card, workdir, log)

    if len(limits) != 5:
        sys.exit(0)

    # Calculate frequentist expected limits

#    makeGrid((limits['m2s'], limits['p2s']), card, workdir, log)
    
#    limits = expected(card, workdir, log)

#    if len(limits) != 5:
#        sys.exit(0)

    # Write result to a single-entry tree (to be merged with calculations for all other signal points)

    model = pointName[0:pointName.find('_')]

    with open('/afs/cern.ch/user/y/yiiyama/src/GammaL/xsec/' + model + '.xsecs') as xsecsource:
        for line in xsecsource:
            p, c, u, n = line.strip().split()
            if p == pointName:
                xsec = float(c)
                uncert = float(u)
                nEvents = int(n)
                break
        else:
            raise RuntimeError('No point ' + pointName + ' found')
    
    outputFile = ROOT.TFile.Open(workdir + '/' + pointName + '.root', 'recreate')
    output = ROOT.TTree('limitTree', 'Limit Tree')

    vPointName = array.array('c', pointName + '\0')
    vXsec = array.array('d', [xsec])
    vXsecErr = array.array('d', [uncert])
    vLimit = array.array('d', [limits['med'] * xsec])
    vLimM2s = array.array('d', [limits['m2s'] * xsec])
    vLimM1s = array.array('d', [limits['m1s'] * xsec])
    vLimP1s = array.array('d', [limits['p1s'] * xsec])
    vLimP2s = array.array('d', [limits['p2s'] * xsec])
    vNEvents = array.array('i', [nEvents])
                
    output.Branch('pointName', vPointName, 'pointName/C')
    output.Branch('xsec', vXsec, 'xsec/D')
    output.Branch('xsecErr', vXsecErr, 'xsecErr/D')
    output.Branch('limit', vLimit, 'limit/D')
    output.Branch('limM2s', vLimM2s, 'limM2s/D')
    output.Branch('limM1s', vLimM1s, 'limM1s/D')
    output.Branch('limP1s', vLimP1s, 'limP1s/D')
    output.Branch('limP2s', vLimP2s, 'limP2s/D')
    output.Branch('nEvents', vNEvents, 'nEvents/I')

    yields = []
    with open(card) as cardFile:
        for line in cardFile:
            cols = line.strip().split()
            if cols[0].endswith('_signal'):
                channel = cols[0][0:cols[0].rfind('_')]
                y = array.array('i', [int(cols[2])])
                output.Branch(channel + 'Yield', y, channel + 'Yield/I')
                yields.append(y)

    output.Fill()
    outputFile.Write()
    outputFile.Close()

    shutil.copyfile(workdir + '/' + pointName + '.root', outputdir + '/' + pointName + '.root')

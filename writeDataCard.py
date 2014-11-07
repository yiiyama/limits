import sys
import os
import ROOT

from datacard import *

sys.path.append('/afs/cern.ch/user/y/yiiyama/src/GammaL/plotstack')
import rootconfig
import locations
from stack import Group
from dataset import Dataset
from GammaL.config import stackConfigs
        
def setupChannel(channel):

    groups = stackConfigs[channel.stackName].groups
    scaleSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + channel.stackName + '.root')

    for group in groups:
        if group.category == Group.OBSERVED:
            count = 0
            for sample in group.samples:
                count += getRateAndCount(sample, channel.cut)[1]
            channel.observed = count

        elif group.category == Group.BACKGROUND:
            rate = 0.
            for sample in group.samples:
                r, count = getRateAndCount(sample, channel.cut)
                if group.name == 'VGamma':
                    r *= scaleSource.Get('TemplateFitError/VGamma').GetY()[0]
                elif group.name == 'JLFake':
                    r *= scaleSource.Get('TemplateFitError/QCD').GetY()[0]

                rate += r

            if rate < 0.: rate = 0.
            if len(group.samples) != 1:
                count = 0

            channel.addProcess(group.name, rate, count)

    scaleSource.Close()

    channel.addProcess('signal', 0., 0, signal = True)


def setNuisances(nuisances, channel):

    nuisances['lumi'][channel.processes['EWK']] = 0.026
    nuisances['lumi'][channel.processes['signal']] = 0.026

    nuisances['effcorr'][channel.processes['EWK']] = 0.08
    nuisances['effcorr'][channel.processes['signal']] = 0.08

    nuisances['ewkxsec'][channel.processes['EWK']] = 0.5

    nuisances['proxyshape'][channel.processes['EGFake']] = 0.2
    if channel.lepton == 'Electron':
        nuisances['proxyshape'][channel.processes['JGFake']] = 0.004
        nuisances['proxyshape'][channel.processes['JLFake']] = 0.2
    elif channel.lepton == 'Muon':
        nuisances['proxyshape'][channel.processes['JGFake']] = 0.4
        nuisances['proxyshape'][channel.processes['JLFake']] = 0.2

    scaleSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + channel.stackName + '.root')
    scaleGr = scaleSource.Get('TemplateFitError/VGammaNoEff')
    nuisances['vgscale'][channel.processes['VGamma']] = boundVal(scaleGr.GetErrorY(0) / scaleGr.GetY()[0])

    groups = stackConfigs[channel.stackName].groups

    for group in groups:
        if group.name != 'EWK' and group.name != 'VGamma': continue

        process = channel.processes[group.name]
        if process.rate() == 0.: continue

        if group.name == 'VGamma':
            scale = scaleSource.Get('TemplateFitError/VGamma').GetY()[0]
        else:
            scale = 1.

        nuisances['jes'][process] = getJESUncert([(sample, scale) for sample in group.samples], channel.cut, process.rate())

    process = channel.processes['VGamma']
    group = next(group for group in groups if group.name == 'VGamma')
    if process.rate() != 0.:
        leptonPtSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/dimuonPt_binned.root')
        ptlWeight = leptonPtSource.Get('ratio')

        rate = 0.
        for sample in group.samples:
            if sample.name in treeStore:
                tree = treeStore[sample.name]
            else:
                sample.loadTree(locations.eventListDir)
                tree = sample.tree

            tree.SetEstimate(tree.GetEntries() + 1)
            if channel.lepton == 'Electron':
                ptVar = 'electron.pt[0]'
            elif channel.lepton == 'Muon':
                ptVar = 'muon.pt[0]'

            nEntries = tree.Draw(ptVar + ':eventSigma * ' + LUMI + ' * puWeight * effScale', channel.cut, 'goff')
            ptArr = tree.GetV1()
            wArr = tree.GetV2()
            for iE in range(nEntries):
                iBin = ptlWeight.FindFixBin(ptArr[iE])
                if iBin == 0: iBin = 1
                elif iBin == ptlWeight.GetNbinsX() + 1: iBin = ptlWeight.GetNbinsX()

                rate += ptlWeight.GetBinContent(iBin) * wArr[iE]

            if sample.name not in treeStore:
                tree = None
                sample.releaseTree()

        scaleGr = scaleSource.Get('TemplateFitError/VGamma')
        rate *= scaleGr.GetY()[0]

        diff = abs(rate - process.rate()) / process.rate()
        if diff > 5.e-4:
            nuisances['vgshape'][process] = boundVal(diff)

    scaleSource.Close()    


if __name__ == '__main__':

    import sys
    import pickle
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: writeDataCard.py [options] outputName')
    parser.add_option('-d', '--directory', dest = 'outputDir', default = '/afs/cern.ch/user/y/yiiyama/work/datacards', help = 'output directory')

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    outputName = args[0]
    if not outputName.endswith('.pkl'):
        print 'Output name must end with .pkl'
        sys.exit(1)

    if os.path.exists(os.environ['TMPDIR'] + '/writeDataCard_tmp.root'):
        print 'Using pre-made skim'
        tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/writeDataCard_tmp.root')
        for key in tmpFile.GetListOfKeys():
            name = key.GetName().replace('eventList_', '')
            treeStore[name] = key.ReadObj()

    else:
        tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/writeDataCard_tmp.root', 'recreate')
        samples = set([])
        for stackName in ['FloatingVGammaE', 'FloatingVGammaM']:
            for group in stackConfigs[stackName].groups:
                if group.category == Group.SIGNAL: continue
                for sample in group.samples:
                    samples.add(sample)
    
        for sample in samples:
            print 'Skimming', sample.name
            sample.loadTree(locations.eventListDir)
            tmpFile.cd()
            tree = sample.tree.CopyTree('(mt >= 100. && met >= 120.) || (mtJESUp >= 100. && metJESUp >= 120.) || (mtJESDown >= 100. && metJESDown >= 120.)')
            sample.releaseTree()
            tree.SetName('eventList_' + sample.name)
            treeStore[sample.name] = tree

        tmpFile.cd()
        tmpFile.Write()

    # setup background and observed

    channels = {}
    for lep, lepton, chCut, stack in [('el', 'Electron', '(mass2 < 81. || mass2 > 101.) && ', 'FloatingVGammaE'), ('mu', 'Muon', '', 'FloatingVGammaM')]:
        for ptCutName, ptCut in [('LowPt', 'photon.pt[0] < 80.'), ('HighPt', 'photon.pt[0] >= 80.')]:
            for htCutName, htCut in [('LowHt', 'ht < 100.'), ('MidHt', 'ht >= 100. && ht < 400.'), ('HighHt', 'ht >= 400.')]:
                for metLow, metHigh in [(120., 200.), (200., 300.), (300., 8000.)]:
                    channels[lep + ptCutName + htCutName + str(int(metLow))] = Channel(lepton, stack, chCut + 'mt >= 100. && ' + ptCut + ' && ' + htCut + ' && met >= ' + str(metLow) + ' && met < ' + str(metHigh))

    for name, ch in channels.items():
        print name
        setupChannel(ch)

# remove background process with 0 count in all channels - necessary?

    nuisances = {
        'lumi': {},
        'effcorr': {},
        'ewkxsec': {},
        'vgscale': {},
        'vgshape': {},
        'proxyshape': {},
        'jes': {},
        'isr': {}
    }

    for name, ch in channels.items():
        setNuisances(nuisances, ch)

    tmpFile.Close()

    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + outputName, 'w') as outputFile:
        pickle.dump((channels, nuisances), outputFile)

    writeDataCard(channels, nuisances, options.outputDir + '/' + outputName.replace('.pkl', '.dat'))

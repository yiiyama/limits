import sys
import os
import ROOT

import datacard

from stack import Group
import locations
from GammaL.config import stackConfigs
        
def setupChannel(channel):

    groups = stackConfigs[channel.stackName].groups

    scaleSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + channel.stackName + '.root')
    vgscale = scaleSource.Get('TemplateFitError/VGamma').GetY()[0]
    scaleGr = scaleSource.Get('TemplateFitError/VGammaNoEff')
    vgscaleRelErr = datacard.boundVal(scaleGr.GetErrorY(0) / scaleGr.GetY()[0])
    jlscale = scaleSource.Get('TemplateFitError/QCD').GetY()[0]
    scaleSource.Close()

    for group in groups:
        if group.category == Group.OBSERVED:
            count = 0
            for sample in group.samples:
                count += datacard.getRateAndCount(sample, channel.cut)[1]
            channel.observed = count

        elif group.category == Group.BACKGROUND:
            process = datacard.Process(group.name)

            for sample in group.samples:
                rate, count = datacard.getRateAndCount(sample, channel.cut)
                if group.name == 'VGamma':
                    rate *= vgscale
                elif group.name == 'JLFake':
                    rate *= jlscale

                if rate < 0.: rate = 0.

                process.addRate(sample.name, rate, count)

            channel.processes[group.name] = process

            if process.count() == 0: continue
            
            if group.name == 'EGFake':
                process.nuisances['proxyshape'] = 0.2

            elif group.name == 'JGFake':
                if channel.lepton == 'Electron':
                    process.nuisances['proxyshape'] = 0.004
                elif channel.lepton == 'Muon':
                    process.nuisances['proxyshape'] = 0.4

            elif group.name == 'JLFake':
                if channel.lepton == 'Electron':
                    process.nuisances['proxyshape'] = 0.2
                elif channel.lepton == 'Muon':
                    process.nuisances['proxyshape'] = 0.2

            else:
                if group.name == 'VGamma':
                    process.nuisances['vgscale'] = vgscaleRelErr
                    process.nuisances['vgshape'] = vgshape(group, channel, vgscale, process.rate())
                    scale = vgscale
    
                elif group.name == 'EWK':
                    process.nuisances['lumi'] = 0.026
                    process.nuisances['effcorr'] = 0.08
                    process.nuisances['ewkxsec'] = 0.5
                    scale = 1.

                process.nuisances['jes'] = datacard.getJESUncert([(sample, scale) for sample in group.samples], channel.cut, process.rate())
                process.nuisances['jer'] = datacard.getJERUncert([(sample, scale) for sample in group.samples], channel.cut, process.rate())


def vgshape(group, channel, vgscale, nominal):
    leptonPtSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/dimuonPt_binned.root')
    ptlWeight = leptonPtSource.Get('ratio')

    rate = 0.
    for sample in group.samples:
        if sample.name in datacard.treeStore:
            tree = datacard.treeStore[sample.name]
        else:
            sample.loadTree(locations.eventListDir)
            tree = sample.tree

        tree.SetEstimate(tree.GetEntries() + 1)
        if channel.lepton == 'Electron':
            ptVar = 'electron.pt[0]'
        elif channel.lepton == 'Muon':
            ptVar = 'muon.pt[0]'

        nEntries = tree.Draw(ptVar + ':eventSigma * ' + datacard.LUMI + ' * puWeight * effScale', channel.cut, 'goff')
        ptArr = tree.GetV1()
        wArr = tree.GetV2()
        for iE in range(nEntries):
            iBin = ptlWeight.FindFixBin(ptArr[iE])
            if iBin == 0: iBin = 1
            elif iBin == ptlWeight.GetNbinsX() + 1: iBin = ptlWeight.GetNbinsX()

            rate += ptlWeight.GetBinContent(iBin) * wArr[iE]

        if sample.name not in datacard.treeStore:
            tree = None
            sample.releaseTree()

    rate *= vgscale

    return datacard.boundVal(abs(rate - nominal) / nominal)


if __name__ == '__main__':

    import sys
    import pickle
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: prepareResultsData.py [options] outputName')

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    outputName = args[0]
    if not outputName.endswith('.pkl'):
        print 'Output name must end with .pkl'
        sys.exit(1)

    if os.path.exists(os.environ['TMPDIR'] + '/prepareResultsData_tmp.root'):
        print 'Using pre-made skim'
        tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/prepareResultsData_tmp.root')
        for key in tmpFile.GetListOfKeys():
            name = key.GetName().replace('eventList_', '')
            datacard.treeStore[name] = key.ReadObj()

    else:
        tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/prepareResultsData_tmp.root', 'recreate')
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
            suffices = ['', 'JESUp', 'JESDown', 'Smeared', 'SmearedUp', 'SmearedDown']
            selection = '||'.join(['(mt{suffix} >= 100. && met{suffix} >= 120.)'.format(suffix = s) for s in suffices])
            tree = sample.tree.CopyTree(selection)
            sample.releaseTree()
            tree.SetName('eventList_' + sample.name)
            datacard.treeStore[sample.name] = tree

        tmpFile.cd()
        tmpFile.Write()

    # setup background and observed

    channels = {}
    for lep, lepton, chCut, stack in [('el', 'Electron', '(mass2 < 81. || mass2 > 101.) && ', 'FloatingVGammaE'), ('mu', 'Muon', '', 'FloatingVGammaM')]:
        for ptCutName, ptCut in [('LowPt', 'photon.pt[0] < 80.'), ('HighPt', 'photon.pt[0] >= 80.')]:
            for htCutName, htCut in [('LowHt', 'ht < 100.'), ('MidHt', 'ht >= 100. && ht < 400.'), ('HighHt', 'ht >= 400.')]:
                for metLow, metHigh in [(120., 200.), (200., 300.), (300., 8000.)]:
                    cut = chCut + 'mt >= 100. && {ptCut} && {htCut} && met >= {metLow} && met < {metHigh}'.format(ptCut = ptCut, htCut = htCut, metLow = metLow, metHigh = metHigh)
                    channelName = lep + ptCutName + htCutName + str(int(metLow))
                    
                    print channelName
                    ch = datacard.Channel(channelName, lepton, stack, cut)
                    setupChannel(ch)
                    channels[channelName] = ch

    tmpFile.Close()

    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + outputName, 'wb') as outputFile:
        pickle.dump(channels, outputFile)

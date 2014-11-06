import sys
import os
import math
import array
from datacard import *

sys.path.append('/afs/cern.ch/user/y/yiiyama/src/GammaL/plotstack')
import rootconfig
from dataset import Dataset
from GammaL.countSignal import getDataset

tmpFile = None

class ISRDatabase(object):
    def __init__(self):
        self.source = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/isrWeights.root')
        self.tree = self.source.Get('isrWeights')

        self.vPointName = array.array('c', ['\0'] * 100)
        self.vNEvents = array.array('I', [0])
        self.vSumOfWeights = array.array('d', [0.])

        self.tree.SetBranchAddress('pointName', self.vPointName)
        self.tree.SetBranchAddress('nEvents', self.vNEvents)
        self.tree.SetBranchAddress('sumOfWeights', self.vSumOfWeights)

    def __del__(self):
        self.source.Close()

    def getScale(self, pointFullName):
        iEntry = 0
        while self.tree.GetEntry(iEntry) > 0:
            if pointFullName in self.vPointName.tostring().strip():
                return self.vNEvents[0] / self.vSumOfWeights[0]
            iEntry += 1

        raise RuntimeError('Point ' + pointFullName + ' not found')

isrDatabase = ISRDatabase()

def getISRUncert(model, pointName, samplesAndScales, cut, nominal):

    scale = isrDatabase.getScale(model + '_' + pointName)
    wstrs = []
    for low, high, w in [(0., 120., 1.), (120., 150., 0.95), (150., 250., 0.90), (250., 8000., 0.8)]:
        wstrs.append('(genBoost > %.0f && genBoost <= %.0f) * %f' % (low, high, w * scale))

    weight = ' + '.join(wstrs)

    rate = 0.
    for sample, s in samplesAndScales:
        rate += getRateAndCount(sample, cut, weight)[0] * s

    shift = (rate - nominal) / nominal

    if abs(shift) < 5.e-4: return 0.

    return boundVal(shift)


def setSignal(model, pointName, channels, nuisances, procName = 'signal'):

    dataset = getDataset(model, pointName)
    if not dataset: return

    for sample in dataset.samples.values():
        if sample.name not in treeStore:
            if tmpFile is not None:
                sample.loadTree('rooth://ncmu40//store/countSignal/' + model)
                tmpFile.cd()
                tree = sample.tree.CopyTree('(mt >= 100. && met >= 120.) || (mtJESUp >= 100. && metJESUp >= 120.) || (mtJESDown >= 100. && metJESDown >= 120.)')
                sample.releaseTree()
                tree.SetName('eventList_' + sample.name)
                treeStore[sample.name] = tree
            else:
                sample.loadTree('rooth://ncmu40//store/countSignal/' + model)
                treeStore[sample.name] = sample.tree

    for name, channel in channels.items():
        scaleSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + channel.stackName + '.root')
        jlScale = scaleSource.Get('TemplateFitError/QCD').GetY()[0]

        candSample = dataset.samples['PhotonAnd' + channel.lepton]
        egSample = dataset.samples['ElePhotonAnd' + channel.lepton]
        jgSample = dataset.samples['FakePhotonAnd' + channel.lepton]
        jlSample = dataset.samples['PhotonAndFake' + channel.lepton]

        rate, count = getRateAndCount(candSample, channel.cut)
        rate -= getRateAndCount(egSample, channel.cut)[0]
        rate -= getRateAndCount(jgSample, channel.cut)[0]
        rate -= getRateAndCount(jlSample, channel.cut)[0] * jlScale

        if procName not in channel.processes:
            channel.addProcess(procName, 0., 0, signal = True)
        
        process = channel.processes[procName]

        if rate < 0.: rate = 0.

        if rate != 0. and model != 'Spectra_gW':
            samplesAndScales = [(candSample, 1.), (egSample, -1.), (jgSample, -1.), (jlSample, -jlScale)]
            isrUncert = getISRUncert(model, pointName, samplesAndScales, channel.cut, rate)

            if process in nuisances['isr']:
                isrUncert *= rate
                isrUncert += nuisances['isr'][process] * process.rate
                isrUncert /= (rate + process.rate)

            nuisances['isr'][process] = isrUncert

        # adding count is technically incorrect; does not matter as long as acceptance is similar between the models
        channel.setProcess(procName, rate + process.rate, count + process.count, signal = True)

        scaleSource.Close()

    if tmpFile is None:
        for sample in dataset.samples.values():
            treeStore.pop(sample.name)
            sample.releaseTree()


def clearSignal(channels, nuisances):
    for channel in channels.values():
        process = channel.processes['signal']
    
        process.rate = 0.
        process.count = 0
        
        nuisances['isr'][process] = 0.


if __name__ == '__main__':
    import sys
    import pickle
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: writeSignalDataCards.py [options] inputName')
    parser.add_option('-d', '--directory', dest = 'outputDir', default = '/afs/cern.ch/user/y/yiiyama/work/datacards')
    parser.add_option('-p', '--point', dest = 'point', default = '')
    parser.add_option('-m', '--model', dest = 'model', default = '')

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    inputName = args[0]

    channels, nuisances = pickle.load(open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + inputName))

    tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/writeDataCard_signal_tmp.root', 'recreate')

    pointList = {}
    pointList['TChiwg'] = dict([('%d' % mchi, [('TChiwg', str(mchi))]) for mchi in range(100, 810, 10)])
    pointList['T5wg'] = dict([('%d_%d' % (mglu, mchi), [('T5wg', '%d_%d' % (mglu, mchi))]) for mglu in range(400, 1550, 50) for mchi in range(25, mglu, 50)])
    pointList['T5wg+TChiwg'] = {}
    for mglu in range(400, 1550, 50):
        for mchi in range(125, mglu, 50):
            pointList['T5wg+TChiwg']['%d_%d' % (mglu, mchi)] = [('T5wg', '%d_%d' % (mglu, mchi)), ('TChiwg', '%d' % mchi)]
    pointList['Spectra_gW'] = {}
    for m3 in range(715, 1565, 50):
        for m2 in range(205, m3, 50):
            pointList['Spectra_gW']['M3_%d_M2_%d' % (m3, m2)] = [('Spectra_gW', 'M3_%d_M2_%d_%s' % (m3, m2, proc)) for proc in ['gg', 'ncp', 'ncm']]

    if options.model:
        models = [options.model]
    else:
        if options.point:
            raise RuntimeError('Point given without model')

        models = pointList.keys()

    for model in models:
        for title in sorted(pointList[model].keys()):
            if options.point and title != options.point: continue

            fullTitle = model + '_' + title
            print fullTitle

            plist = pointList[model][title]

            clearSignal(channels, nuisances)
            for source, point in plist:
                setSignal(source, point, channels, nuisances)

            writeDataCard(channels, nuisances, options.outputDir + '/' + fullTitle + '.dat')

    tmpFile.Close()

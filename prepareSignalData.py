import sys
import os
import math
import array
import datacard

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
        rate += datacard.getRateAndCount(sample, cut, weight)[0] * s

    shift = (rate - nominal) / nominal

    if abs(shift) < 5.e-4: return 0.

    return datacard.boundVal(shift)


def setSignal(model, pointName, processes, channels, ratescale):

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

    stackNames = set([c.stackName for c in channels.values()])
    scaleSources = dict([(s, ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + s + '.root')) for s in stackNames])

    for channelName, channel in channels.items():
        jlScale = scaleSources[channel.stackName].Get('TemplateFitError/QCD').GetY()[0]

        candSample = dataset.samples['PhotonAnd' + channel.lepton]
        egSample = dataset.samples['ElePhotonAnd' + channel.lepton]
        jgSample = dataset.samples['FakePhotonAnd' + channel.lepton]
        jlSample = dataset.samples['PhotonAndFake' + channel.lepton]

        rate, count = getRateAndCount(candSample, channel.cut)
        rate -= getRateAndCount(egSample, channel.cut)[0]
        rate -= getRateAndCount(jgSample, channel.cut)[0]
        rate -= getRateAndCount(jlSample, channel.cut)[0] * jlScale

        if rate < 0.: rate = 0.

        rate *= ratescale

        if channelName not in processes:
            processes[channelName] = datacard.Process('signal', signal = True)

        process = processes[channelName]

        if rate != 0.:
            process.nuisances['lumi'] = 0.026
            process.nuisances['effcorr'] = 0.08

            if model == 'Spectra_gW':
                process.nuisances['isr'] = 0.05
            else:
                samplesAndScales = [(candSample, 1.), (egSample, -1.), (jgSample, -1.), (jlSample, -jlScale)]
                isrUncert = getISRUncert(model, pointName, samplesAndScales, channel.cut, rate)
    
                if 'isr' in process.nuisances:
                    isrUncert *= rate
                    isrUncert += process.nuisances['isr'] * process.rate()
                    isrUncert /= (rate + process.rate())
    
                process.nuisances['isr'] = isrUncert

            jesUncert = datacard.getJESUncert([(candSample, 1.)], channel.cut, rate)

            if 'jes' in process.nuisances:
                jesUncert *= rate
                jesUncert += process.nuisances['jes'] * process.rate()
                jesUncert /= (rate + process.rate())

            process.nuisances['jes'] = jesUncert

        process.addRate(model + '_' + pointName, rate, count, genInfo = dataset.GenInfo(dataset.sigma, dataset.sigmaRelErr, dataset.nEvents))

    for source in scaleSources.values():
        source.Close()

    if tmpFile is None:
        for sample in dataset.samples.values():
            treeStore.pop(sample.name)
            sample.releaseTree()


if __name__ == '__main__':
    import sys
    import pickle
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: writeSignalDataCards.py [options] inputName')
    parser.add_option('-p', '--point', dest = 'point', default = '')
    parser.add_option('-m', '--model', dest = 'model', default = '')

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    inputName = args[0]

    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + inputName) as source:
        channels = pickle.load(source)

    tmpFile = ROOT.TFile.Open(os.environ['TMPDIR'] + '/writeDataCard_signal_tmp.root', 'recreate')

    scale = {}
    pointList = {}
    
    pointList['TChiwg'] = {}
    for mchi in range(100, 810, 10):
        pointList['TChiwg']['%d' % mchi] = [('TChiwg', str(mchi))]

    pointList['T5wg'] = {}
    for mglu in range(700, 1550, 50):
        for mchi in range(25, mglu, 50):
            pointList['T5wg']['%d_%d' % (mglu, mchi)] = [('T5wg', '%d_%d' % (mglu, mchi))]

#    pointList['T5wg+TChiwg'] = {}
#    for mglu in range(400, 1550, 50):
#        for mchi in range(125, mglu, 50):
#            pointList['T5wg+TChiwg']['%d_%d' % (mglu, mchi)] = [('T5wg', '%d_%d' % (mglu, mchi)), ('TChiwg', '%d' % mchi)]

    pointList['Spectra_gW'] = {}
    for m3 in range(715, 1565, 50):
        for m2 in range(205, m3, 50):
            pointList['Spectra_gW']['M3_%d_M2_%d' % (m3, m2)] = [('Spectra_gW', 'M3_%d_M2_%d_%s' % (m3, m2, proc)) for proc in ['gg', 'ncp', 'ncm']]

    pointList['Spectra_gW_gg'] = {}
    for m3 in range(715, 1565, 50):
        for m2 in range(205, m3, 50):
            pointList['Spectra_gW_gg']['M3_%d_M2_%d' % (m3, m2)] = [('Spectra_gW', 'M3_%d_M2_%d_gg' % (m3, m2))]
    scale['Spectra_gW_gg'] = 1. / (4. / 9. * 0.23) * 0.5

    pointList['Spectra_gW_nc'] = {}
    for m3 in range(715, 1565, 50):
        for m2 in range(205, m3, 50):
            pointList['Spectra_gW_nc']['M3_%d_M2_%d' % (m3, m2)] = [('Spectra_gW', 'M3_%d_M2_%d_%s' % (m3, m2, proc)) for proc in ['ncp', 'ncm']]
    scale['Spectra_gW_nc'] = 1. / 0.23

    if options.model:
        models = [options.model]
    else:
        if options.point:
            raise RuntimeError('Point given without model')

        models = pointList.keys()

    for model in models:

        outputFileName = '/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + model
        if options.point:
            outputFileName += '_' + options.point
        outputFileName += '.pkl'

        if model in scale:
            ratescale = scale[model]
        else:
            ratescale = 1.

        outputFile = open(outputFileName, 'wb')
        data = {}

        for title in sorted(pointList[model].keys()):
            if options.point and title != options.point: continue

            fullTitle = model + '_' + title
            print fullTitle

            signals = pointList[model][title]
            processes = {}
            for source, point in signals:
                setSignal(source, point, processes, channels, ratescale)

            data[fullTitle] = processes

        pickle.dump(data, outputFile)

        outputFile.close()

    tmpFile.Close()

import sys
import math

import ROOT

class Channel(object):
    def __init__(self, name, sourceName, histoName, lepton, datasetName, x):
        self.name = name
        self.sourceName = sourceName
        self.histoName = histoName
        self.lepton = lepton
        self.datasetName = datasetName
        self.x = x

        self.observed = 0
        self.processes = {}

    def addProcess(self,  process):
        self.processes[process.name] = process


class Process(object):
    def __init__(self, name, channel, rate, count = 0, signal = False):
        self.name = name
        self.channel = channel
        self.rate = rate
        self.count = count
        self.signal = signal

 
def setupFromHistograms(pointName):

    def getBinContent(source, bin, plotName, rawPlotName = ''):
        hist = source.Get(plotName)
        rate = hist.GetBinContent(bin) * hist.GetXaxis().GetBinWidth(bin)
        if bin == hist.GetNbinsX() - 1:
            rate += hist.GetBinContent(bin + 1)

        if not rawPlotName: return rate, 0

        hist = source.Get(rawPlotName)
        count = hist.GetBinContent(bin)
        if bin == hist.GetNbinsX() - 1:
            count += hist.GetBinContent(bin + 1)
    
        return rate, int(count)


    elSourceName = '/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FixedScalesE.root'
    muSourceName = '/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FixedScalesM.root'

    lowPtPlot = 'MetMt100Pt4080'
    highPtPlot = 'MetMt100Pt80'

    channels = [
        Channel('el120LowPt', elSourceName, lowPtPlot, 'Electron', 'DataE', 120.),
        Channel('el200LowPt', elSourceName, lowPtPlot, 'Electron', 'DataE', 200.),
        Channel('el300LowPt', elSourceName, lowPtPlot, 'Electron', 'DataE', 300.),
        Channel('mu120LowPt', muSourceName, lowPtPlot, 'Muon', 'DataM', 120.),
        Channel('mu200LowPt', muSourceName, lowPtPlot, 'Muon', 'DataM', 200.),
        Channel('mu300LowPt', muSourceName, lowPtPlot, 'Muon', 'DataM', 300.),
        Channel('el120HighPt', elSourceName, highPtPlot, 'Electron', 'DataE', 120.),
        Channel('el200HighPt', elSourceName, highPtPlot, 'Electron', 'DataE', 200.),
        Channel('el300HighPt', elSourceName, highPtPlot, 'Electron', 'DataE', 300.),
        Channel('mu120HighPt', muSourceName, highPtPlot, 'Muon', 'DataM', 120.),
        Channel('mu200HighPt', muSourceName, highPtPlot, 'Muon', 'DataM', 200.),
        Channel('mu300HighPt', muSourceName, highPtPlot, 'Muon', 'DataM', 300.)
    ]

    bkgNames = ['EGFake', 'JGFake', 'JLFake', 'VGamma', 'EWK']

    signalSource = ROOT.TFile.Open('/tmp/yiiyama/countSignal/plots/' + pointName + '.root')

    for ch in channels:
        source = ROOT.TFile.Open(ch.sourceName)

        bin = source.Get('groups/' + ch.histoName + '/' + ch.histoName + '_Observed').GetXaxis().FindFixBin(ch.x)

        ch.observed = int(getBinContent(source, bin, 'groups/' + ch.histoName + '/' + ch.histoName + '_Observed')[0])
    
        for proc in bkgNames:
            rawPlotName = 'components/' + ch.histoName + '/' + ch.histoName + '_' + ch.datasetName
            if proc == 'EGFake':
                rawPlotName += '_ElePhotonAnd' + ch.lepton + '_Raw'
            elif proc == 'JGFake':
                rawPlotName += '_FakePhotonAnd' + ch.lepton + '_Raw'
            elif proc == 'JLFake':
                rawPlotName += '_PhotonAndFake' + ch.lepton + '_Raw'
            else:
                rawPlotName = ''

            rate, count = getBinContent(source, bin, 'groups/' + ch.histoName + '/' + ch.histoName + '_' + proc, rawPlotName)
            if rate < 0.: rate = 0.

            ch.addProcess(Process(proc, ch, rate, count = count))

        rate, count = getBinContent(signalSource, bin, ch.lepton + '/groups/' + ch.histoName + '/' + ch.histoName + '_' + pointName, ch.lepton + '/components/' + ch.histoName + '/' + ch.histoName + '_' + pointName + '_PhotonAnd' + ch.lepton + '_Raw')
        if rate < 0.: rate = 0.
        ch.addProcess(Process('signal', ch, rate, count = count, signal = True))

        source.Close()

    signalSource.Close()

    chList = dict([(ch.name, ch) for ch in channels])

    return chList


def getNuisances(channels):

    nuisances = {
        'lumi': {},
        'ewkxsec': {},
        'vgscale_Electron': {},
        'vgscale_Muon': {},
        'effscale': {}
    }

    for channel in channels.values():
        nuisances['lumi'][channel.processes['EWK']] = 0.026
        nuisances['lumi'][channel.processes['signal']] = 0.026

        nuisances['ewkxsec'][channel.processes['EWK']] = 0.5

        source = ROOT.TFile.Open(channel.sourceName)
        scaleGr = source.Get('TemplateFitError/VGamma')
        nuisances['vgscale_' + channel.lepton][channel.processes['VGamma']] = scaleGr.GetErrorY(0) / scaleGr.GetY()[0]
        source.Close()

        nuisances['effscale'][channel.processes['EWK']] = 0.1
        nuisances['effscale'][channel.processes['signal']] = 0.1

    return nuisances


def writeDataCard(channels, nuisances, cardName):

    HLINE = '----------------------------------------'

    def alignColumns(lines):
        maxWidth = 0
        for line in lines:
            if line == HLINE: continue
            if len(line[0]) > maxWidth: maxWidth = len(line[0])
        form = '%-' + str(maxWidth) + 's'

        for line in lines:
            if line == HLINE: continue
            line[0] = form % line[0]
    
        maxWidth = 0
        for line in lines:
            if line == HLINE: continue
            for word in line[1:]:
                if len(word) > maxWidth: maxWidth = len(word)
        form = '%' + str(maxWidth) + 's'

        for line in lines:
            if line == HLINE: continue
            for iW in range(1, len(line)):
                line[iW] = form % line[iW]

    bkgNames = []
    for channel in channels.values():
        names = set([name for name, process in channel.processes.items() if not process.signal])
        if len(bkgNames) == 0:
            bkgNames = sorted(list(names))
            continue
        elif len(set(bkgNames) - names) or len(names - set(bkgNames)):
            raise RuntimeError('Background names do not match')

    channelNames = sorted(channels.keys())

    nuisanceNames = sorted(nuisances.keys())

    lines = []

    lines.append(['imax', str(len(channels))])
    lines.append(['jmax', str(len(bkgNames))])
    lines.append(['kmax'])

    lines.append(HLINE)

    blockStart = len(lines)

    line = ['bin'] + channelNames
    lines.append(line)

    line = ['observation'] + [str(channels[ch].observed) for ch in channelNames]
    lines.append(line)

    alignColumns(lines[blockStart:])

    lines.append(HLINE)

    blockStart = len(lines)

    binLine = ['bin']
    procNameLine = ['process']
    procIDLine = ['process']
    rateLine = ['rate']
    for ch in channelNames:
        channel = channels[ch]
        binLine += ([ch] * (len(channel.processes)))
        procNameLine.append('signal')
        procIDLine.append('0')
        rateLine.append('%.3f' % channel.processes['signal'].rate)
        bkgID = 1
        for proc in bkgNames:
            procNameLine.append(proc)
            procIDLine.append(str(bkgID))
            rateLine.append('%.3f' % channel.processes[proc].rate)
            bkgID += 1

    lines.append(binLine)
    lines.append(procNameLine)
    lines.append(procIDLine)
    lines.append(rateLine)

    lines.append(HLINE)

    kmax = 0

    for name in nuisanceNames:
        nuisance = nuisances[name]
        line = [name + ' lnN']
        for ch in channelNames:
            channel = channels[ch]
            for proc in ['signal'] + bkgNames:
                process = channel.processes[proc]
                try:
                    val = 1. + nuisance[process]
                    line.append('%.3f' % val)
                except KeyError:
                    line.append('-')

        lines.append(line)
        kmax += 1

    for ch in channelNames:
        channel = channels[ch]
        for proc in ['signal'] + bkgNames:
            process = channel.processes[proc]
            line = [ch + '_' + proc + ' gmN ' + str(process.count)]
            for chproc in [(c, p) for c in channelNames for p in ['signal'] + bkgNames]:
                if chproc == (ch, proc):
                    if process.count == 0:
                        line.append('%.2e' % 0.)
                    else:
                        line.append('%.2e' % (process.rate / process.count))
                else:
                    line.append('-')
            
            lines.append(line)
            kmax += 1

    alignColumns(lines[blockStart:])

    kmaxLine = next(line for line in lines if line[0] == 'kmax')
    kmaxLine.append(str(kmax))

    with open(cardName, 'w') as datacard:
        for line in lines:
            if line == HLINE:
                datacard.write(line + '\n')
            else:
                datacard.write(' '.join(line) + '\n')


if __name__ == '__main__':

    pointName = sys.argv[1]
    if len(sys.argv) == 3:
        suffix = sys.argv[2]
    else:
        suffix = ''

    channels = setupFromHistograms(pointName)
    nuisances = getNuisances(channels)

    writeDataCard(channels, nuisances, '/afs/cern.ch/user/y/yiiyama/work/datacards' + suffix + '/' + pointName + '.dat')

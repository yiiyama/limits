import ROOT

class Process(object):
    def __init__(self, name, channel, rate, count = 0, signal = False):
        self.name = name
        self.channel = channel
        self.rate = rate
        self.count = count
        self.signal = signal


class Channel(object):
    def __init__(self, histoName, lepton, datasetName, xmin, xmax = 8000.):
        self.source = None
        self.histoName = histoName
        self.lepton = lepton
        self.datasetName = datasetName
        self.xmin = xmin
        self.xmax = xmax

        self.observed = 0
        self.processes = {}

    def addProcess(self, name, rate, count = 0, signal = False):
        self.processes[name] = Process(name, self, rate, count, signal)

    def setProcess(self, name, rate, count = 0, signal = False):
        process = self.processes[name]
        process.rate = rate
        process.count = count
        process.signal = signal


def getBinContent(hist, xmin, xmax, rawHist = None):

    iBin = hist.GetXaxis().FindFixBin(xmin)
    rate = 0.
    while hist.GetXaxis().GetBinLowEdge(iBin) < xmax:
        if iBin == hist.GetNbinsX():
            rate += hist.GetBinContent(iBin)
            break
        else:
            rate += hist.GetBinContent(iBin) * hist.GetXaxis().GetBinWidth(iBin)
            iBin += 1

    if rawHist is None: return rate

    iBin = rawHist.GetXaxis().FindFixBin(xmin)
    count = 0.
    while rawHist.GetXaxis().GetBinLowEdge(iBin) < xmax:
        count += rawHist.GetBinContent(iBin)
        if iBin == rawHist.GetNbinsX():
            break
        iBin += 1

    return rate, int(count)

        
def setupChannels(channels):

    bkgNames = ['EGFake', 'JGFake', 'JLFake', 'VGamma', 'EWK']
    bkgRates = dict([(proc, 0.) for proc in bkgNames])

    for ch in channels.values():
        histName = 'groups/' + ch.histoName + '/' + ch.histoName + '_Observed'
        ch.observed = int(getBinContent(ch.source.Get(histName), ch.xmin, ch.xmax))
    
        for proc in bkgNames:
            rawPlotName = 'components/' + ch.histoName + '/' + ch.histoName + '_' + ch.datasetName
            if proc == 'EGFake':
                rawPlotName += '_ElePhotonAnd' + ch.lepton + '_Raw'
            elif proc == 'JGFake':
                rawPlotName += '_FakePhotonAnd' + ch.lepton + '_Raw'
            elif proc == 'JLFake':
                rawPlotName += '_PhotonAndFake' + ch.lepton + '_Raw'
            else:
                rawPlotName += '_PhotonAnd' + ch.lepton + '_Raw'

            histName = 'groups/' + ch.histoName + '/' + ch.histoName + '_' + proc
            rate, count = getBinContent(ch.source.Get(histName), ch.xmin, ch.xmax, ch.source.Get(rawPlotName))
                
            if rate < 0.: rate = 0.

            ch.addProcess(proc, rate, count = count)

            bkgRates[proc] += rate

        ch.addProcess('signal', 0., 0, signal = True)

    for proc, rate in bkgRates.items():
        if rate != 0.: continue
        for ch in channels.values():
            ch.processes.pop(proc)


def boundVal(val, bound = 1.):
    if val > bound: return bound - 0.001
    if val < -bound: return -bound + 0.001
    return val


def getNuisances(channels):

    nuisances = {
        'lumi': {},
        'effcorr': {},
        'ewkxsec': {},
        'vgscale_Electron': {},
        'vgscale_Muon': {},
        'vgshape': {},
        'jes': {}
    }

    for ch in channels.values():
        nuisances['lumi'][ch.processes['EWK']] = 0.026
        nuisances['lumi'][ch.processes['signal']] = 0.026

        nuisances['effcorr'][ch.processes['EWK']] = 0.08
        nuisances['effcorr'][ch.processes['signal']] = 0.08

        nuisances['ewkxsec'][ch.processes['EWK']] = 0.5

    return nuisances


def writeDataCard(channels, nuisances, cardName, allowZeroSignal = False):

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
    for ch in channels.values():
        names = set([name for name, process in ch.processes.items() if not process.signal])
        if len(bkgNames) == 0:
            bkgNames = sorted(list(names))
            continue
        elif len(set(bkgNames) - names) or len(names - set(bkgNames)):
            raise RuntimeError('Background names do not match')

    if allowZeroSignal:
        channelNames = sorted(channels.keys())
    else:
        channelNames = sorted([name for name, channel in channels.items() if channel.processes['signal'].rate > 0.])

    nuisanceNames = sorted(nuisances.keys())

    lines = []

    lines.append(['imax', str(len(channelNames))])
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
        rateLine.append('%.3e' % channel.processes['signal'].rate)
        bkgID = 1
        for proc in bkgNames:
            procNameLine.append(proc)
            procIDLine.append(str(bkgID))
            rateLine.append('%.3e' % channel.processes[proc].rate)
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
                    line.append('%.3e' % val)
                except KeyError:
                    line.append('-')

        lines.append(line)
        kmax += 1

    for ch in channelNames:
        channel = channels[ch]
        for proc in ['signal'] + bkgNames:
            process = channel.processes[proc]
            if process.count == 0: continue

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

    # setup background and observed

    OUTPUTDIR = '/afs/cern.ch/user/y/yiiyama/output/GammaL/main'
    
    lowHtLowPtPlot = 'MetHighMtLowHtLowPhotonPt'
    highHtLowPtPlot = 'MetHighMtHighHtLowPhotonPt'
    noJetHighPtPlot = 'MetHighMtNoJetHighPhotonPt'
    midHtHighPtPlot = 'MetHighMtMidHtHighPhotonPt'
    highHtHighPtPlot = 'MetHighMtHighHtHighPhotonPt'
    
    channels = {
        'el120LowHtLowPt': Channel(lowHtLowPtPlot, 'Electron', 'DataE', 120., 200.),
        'el200LowHtLowPt': Channel(lowHtLowPtPlot, 'Electron', 'DataE', 200.),
        'mu120LowHtLowPt': Channel(lowHtLowPtPlot, 'Muon', 'DataM', 120., 200.),
        'mu200LowHtLowPt': Channel(lowHtLowPtPlot, 'Muon', 'DataM', 200.),
        'el120HighHtLowPt': Channel(highHtLowPtPlot, 'Electron', 'DataE', 120., 200.),
        'el200HighHtLowPt': Channel(highHtLowPtPlot, 'Electron', 'DataE', 200.),
        'mu120HighHtLowPt': Channel(highHtLowPtPlot, 'Muon', 'DataM', 120., 200.),
        'mu200HighHtLowPt': Channel(highHtLowPtPlot, 'Muon', 'DataM', 200.),
        'el120NoJetHighPt': Channel(noJetHighPtPlot, 'Electron', 'DataE', 120., 200.),
        'el200NoJetHighPt': Channel(noJetHighPtPlot, 'Electron', 'DataE', 200.),
        'mu120NoJetHighPt': Channel(noJetHighPtPlot, 'Muon', 'DataM', 120., 200.),
        'mu200NoJetHighPt': Channel(noJetHighPtPlot, 'Muon', 'DataM', 200.),
        'el120MidHtHighPt': Channel(midHtHighPtPlot, 'Electron', 'DataE', 120., 200.),
        'el200MidHtHighPt': Channel(midHtHighPtPlot, 'Electron', 'DataE', 200.),
        'mu120MidHtHighPt': Channel(midHtHighPtPlot, 'Muon', 'DataM', 120., 200.),
        'mu200MidHtHighPt': Channel(midHtHighPtPlot, 'Muon', 'DataM', 200.),
        'el120HighHtHighPt': Channel(highHtHighPtPlot, 'Electron', 'DataE', 120., 200.),
        'el200HighHtHighPt': Channel(highHtHighPtPlot, 'Electron', 'DataE', 200.),
        'mu120HighHtHighPt': Channel(highHtHighPtPlot, 'Muon', 'DataM', 120., 200.),
        'mu200HighHtHighPt': Channel(highHtHighPtPlot, 'Muon', 'DataM', 200.),
    }

#    noJetHighPtPlot = 'MetHighMtNoJetHighPhotonPt'
#    midHtHighPtPlot = 'MetHighMtMidHtHighPhotonPt'
#    highHtHighPtPlot = 'MetHighMtHighHtHighPhotonPt'
#    
#    channels = {
#        'el120NoJetHighPt': Channel(noJetHighPtPlot, 'Electron', 'DataE', 120., 200.),
#        'el200NoJetHighPt': Channel(noJetHighPtPlot, 'Electron', 'DataE', 200., 300.),
#        'el300NoJetHighPt': Channel(noJetHighPtPlot, 'Electron', 'DataE', 300.),
#        'mu120NoJetHighPt': Channel(noJetHighPtPlot, 'Muon', 'DataM', 120., 200.),
#        'mu200NoJetHighPt': Channel(noJetHighPtPlot, 'Muon', 'DataM', 200., 300.),
#        'mu300NoJetHighPt': Channel(noJetHighPtPlot, 'Muon', 'DataM', 300.),
#        'el120MidHtHighPt': Channel(midHtHighPtPlot, 'Electron', 'DataE', 120., 200.),
#        'el200MidHtHighPt': Channel(midHtHighPtPlot, 'Electron', 'DataE', 200., 300.),
#        'el300MidHtHighPt': Channel(midHtHighPtPlot, 'Electron', 'DataE', 300.),
#        'mu120MidHtHighPt': Channel(midHtHighPtPlot, 'Muon', 'DataM', 120., 200.),
#        'mu200MidHtHighPt': Channel(midHtHighPtPlot, 'Muon', 'DataM', 200., 300.),
#        'mu300MidHtHighPt': Channel(midHtHighPtPlot, 'Muon', 'DataM', 300.),
#        'el120HighHtHighPt': Channel(highHtHighPtPlot, 'Electron', 'DataE', 120., 200.),
#        'el200HighHtHighPt': Channel(highHtHighPtPlot, 'Electron', 'DataE', 200., 300.),
#        'el300HighHtHighPt': Channel(highHtHighPtPlot, 'Electron', 'DataE', 300.),
#        'mu120HighHtHighPt': Channel(highHtHighPtPlot, 'Muon', 'DataM', 120., 200.),
#        'mu200HighHtHighPt': Channel(highHtHighPtPlot, 'Muon', 'DataM', 200., 300.),
#        'mu300HighHtHighPt': Channel(highHtHighPtPlot, 'Muon', 'DataM', 300.),
#    }

    source = {
        'Electron': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaE.root'),
        'Muon': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaM.root')
    }

    for channel in channels.values():
        channel.source = source[channel.lepton]
    
    setupChannels(channels)
    nuisances = getNuisances(channels)

    sourceScale = source

    sourceJESUp = {
        'Electron': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaJESUpE.root'),
        'Muon': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaJESUpM.root')
    }
    sourceJESDown = {
        'Electron': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaJESDownE.root'),
        'Muon': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaJESDownM.root')
    }

    for channelName, channel in channels.items():
        scaleGr = sourceScale[channel.lepton].Get('TemplateFitError/VGammaNoEff')
        nuisances['vgscale_' + channel.lepton][channel.processes['VGamma']] = boundVal(scaleGr.GetErrorY(0) / scaleGr.GetY()[0])

        for procName, process in channel.processes.items():
            if procName == 'signal': continue
            if process.rate == 0.: continue
            histName = 'groups/' + channel.histoName + '/' + channel.histoName + '_' + procName
            up = (getBinContent(sourceJESUp[channel.lepton].Get(histName), channel.xmin, channel.xmax) - process.rate) / process.rate
            down = (getBinContent(sourceJESDown[channel.lepton].Get(histName), channel.xmin, channel.xmax) - process.rate) / process.rate

            if abs(up) < 5.e-4 and abs(down) < 5.e-4: continue

            # always use sign of upward JES shift
            if abs(up) > abs(down):
                nuisances['jes'][process] = boundVal(up)
                
            elif abs(down) > abs(up):
                if down * up < 0.:
                    nuisances['jes'][process] = boundVal(-down)
                else:
                    nuisances['jes'][process] = boundVal(down)

    for src in sourceJESUp.values():
        src.Close()
    for src in sourceJESDown.values():
        src.Close()

    sourcePtL = {
        'Electron': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaPtLE.root'),
        'Muon': ROOT.TFile.Open(OUTPUTDIR + '/FloatingVGammaPtLM.root')
    }

    for channelName, channel in channels.items():
        process = channel.processes['VGamma']

        if process.rate == 0.: continue

        histName = 'groups/' + channel.histoName + '/' + channel.histoName + '_VGamma'
        diff = abs(getBinContent(sourcePtL[channel.lepton].Get(histName), channel.xmin, channel.xmax) - process.rate) / process.rate

        if diff < 5.e-4: continue

        nuisances['vgshape'][process] = boundVal(diff)


    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + outputName, 'w') as outputFile:
        pickle.dump((channels, nuisances), outputFile)

    writeDataCard(channels, nuisances, options.outputDir + '/' + outputName.replace('.pkl', '.dat'), allowZeroSignal = True)

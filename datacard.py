import ROOT

LUMI = '19712.'

treeStore = {}

class Process(object):
    def __init__(self, name, channel, rate, count, signal = False):
        self.name = name
        self.channel = channel
        self.rate = rate
        self.count = count
        self.signal = signal


class Channel(object):
    def __init__(self, lepton, stackName, cut):
        self.lepton = lepton
        self.stackName = stackName
        self.cut = cut

        self.observed = 0
        self.processes = {}

    def addProcess(self, name, rate, count, signal = False):
        self.processes[name] = Process(name, self, rate, count, signal = signal)

    def setProcess(self, name, rate, count, signal = False):
        process = self.processes[name]
        process.rate = rate
        process.count = count
        process.signal = signal

    def bkgTotal(self):
        return reduce(lambda x, y: x + y, [p.rate for p in self.processes.values() if not p.signal])


def getRateAndCount(sample, cut, weight = ''):

    ROOT.gROOT.cd()
    counter = ROOT.TH1D('counter', 'counter', 1, 0., 1.)

    if sample.name in treeStore:
        tree = treeStore[sample.name]
    else:
        sample.loadTree(locations.eventListDir)
        tree = sample.tree

    ROOT.gROOT.cd()
    wstr = 'eventSigma * ' + LUMI + ' * puWeight * effScale'
    if weight:
        wstr += ' * (' + weight + ')'
    tree.Draw('0.5>>+counter', wstr + ' * (' + cut + ')', 'goff')

    if sample.name not in treeStore:        
        tree = None
        sample.releaseTree()

    rate = counter.GetBinContent(1)
    count = int(counter.GetEntries())
    counter.Delete()

    return rate, count


def boundVal(val, bound = 1.):

    if val > bound: return bound - 0.001
    if val < -bound: return -bound + 0.001
    return val


def getJESUncert(samplesAndScales, cut, nominal):

    up = 0.
    down = 0.
    for sample, scale in samplesAndScales:
        cutUp = cut.replace(' met ', ' metJESUp ').replace(' mt ', ' mtJESUp ').replace(' ht ', ' htJESUp ')
        up += getRateAndCount(sample, cutUp)[0] * scale

        cutDown = cut.replace(' met ', ' metJESDown ').replace(' mt ', ' mtJESDown ').replace(' ht ', ' htJESDown ')
        down += getRateAndCount(sample, cutDown)[0] * scale

    shiftUp = (nominal - up) / nominal
    shiftDown = (down - nominal) / nominal

    if abs(shiftUp) < 5.e-4 and abs(shiftDown) < 5.e-4:
        return 0.

    # always use sign of upward JES shift
    if abs(shiftUp) > abs(shiftDown):
        return boundVal(shiftUp)
        
    else:
        if shiftDown * shiftUp < 0.:
            return boundVal(-shiftDown)
        else:
            return boundVal(shiftDown)


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
    for ch in channels.values():
        names = set([name for name, process in ch.processes.items() if not process.signal])
        if len(bkgNames) == 0:
            bkgNames = sorted(list(names))
            continue
        elif len(set(bkgNames) - names) or len(names - set(bkgNames)):
            raise RuntimeError('Background names do not match')

    channelNames = sorted(channels.keys())
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
                if process in nuisance and nuisance[process] > 0.001:
                    line.append('%9.3f' % (1. + nuisance[process]))
                else:
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

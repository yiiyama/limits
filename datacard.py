import sys
import math
import collections
import ROOT

sys.path.append('/afs/cern.ch/user/y/yiiyama/src/GammaL/plotstack')
import locations

LUMI = '19712.'

treeStore = {}

GenInfo = collections.namedtuple('GenInfo', ['xsec', 'relErr', 'nEvents'])

class Process(object):
    def __init__(self, name, signal = False):
        self.name = name
        self.signal = signal
        self.nuisances = {}

        self.rates = {}

    def addRate(self, sample, rate, count):
        self.rates[sample] = (rate, count)

    def rate(self):
        return sum([r for r, c in self.rates.values()])

    def count(self):
        if len(self.rates) == 0: return 0

        # effective count = (sum{count * w})^2 / sum{count * w^2}
        # r = c * w
        R = self.rate()
        if R <= 0.: return 0

        R2 = math.pow(R, 2.)
        D = sum([r * r / c for r, c in self.rates.values() if c != 0])
        return int(round(R2 / D))


class Channel(object):
    def __init__(self, name, lepton, stackName, cut):
        self.name = name
        self.lepton = lepton
        self.stackName = stackName
        self.cut = cut

        self.observed = 0
        self.processes = {}

    def bkgTotal(self):
        return sum([p.rate() for p in self.processes.values() if not p.signal])


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
    
    if nominal <= 0.: return 0.

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


def getJERUncert(samplesAndScales, cut, nominal):
    
    if nominal <= 0.: return 0.

    center = 0.
    up = 0.
    down = 0.
    for sample, scale in samplesAndScales:
#        cut = cut.replace(' met ', ' metSmeared ').replace(' mt ', ' mtSmeared ').replace(' ht ', ' htSmeared ')
#        center += getRateAndCount(sample, cut)[0] * scale

        cutUp = cut.replace(' met ', ' metSmearedUp ').replace(' mt ', ' mtSmearedUp ').replace(' ht ', ' htSmearedUp ')
        up += getRateAndCount(sample, cutUp)[0] * scale

#        cutDown = cut.replace(' met ', ' metSmearedDown ').replace(' mt ', ' mtSmearedDown ').replace(' ht ', ' htSmearedDown ')
#        down += getRateAndCount(sample, cutDown)[0] * scale

    shiftUp = (nominal - up) / nominal
#    shiftDown = (down - nominal) / nominal
    shiftDown = 0.

    if abs(shiftUp) < 5.e-4 and abs(shiftDown) < 5.e-4:
        return 0.

    if abs(shiftUp) > abs(shiftDown):
        return boundVal(shiftUp)
        
    else:
        if shiftDown * shiftUp < 0.:
            return boundVal(-shiftDown)
        else:
            return boundVal(shiftDown)


def writeDataCard(channels, cardName):

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

    processNames = []
    for channel in channels.values():
        names = set([name for name, process in channel.processes.items()])
        if len(processNames) == 0:
            processNames = sorted(list(names))
            continue
        elif len(set(processNames) - names) or len(names - set(processNames)):
            raise RuntimeError('Background names do not match')

    nBkg = len(processNames)

    hasSignal = 'signal' in processNames
    if hasSignal:
        nBkg -= 1
        processNames.remove('signal')
        processNames.insert(0, 'signal')

    channelNames = sorted(channels.keys())
    nuisanceNamesAll = []
    for channel in channels.values():
        for process in channel.processes.values():
            nuisanceNamesAll += process.nuisances.keys()

    nuisanceNames = sorted(list(set(nuisanceNamesAll)))

    lines = []

    lines.append(['imax', str(len(channelNames))])
    lines.append(['jmax', str(nBkg)])
    lines.append(['kmax'])
    kmaxLine = lines[-1]

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
        binLine += ([ch] * len(channel.processes))
        procID = 0 if hasSignal else 1
        for proc in processNames:
            procNameLine.append(proc)
            procIDLine.append(str(procID))
            rateLine.append('%.3e' % channel.processes[proc].rate())
            procID += 1

    lines.append(binLine)
    lines.append(procNameLine)
    lines.append(procIDLine)
    lines.append(rateLine)

    lines.append(HLINE)

    kmax = 0

    for name in nuisanceNames:
        line = [name + ' lnN']
        for ch in channelNames:
            channel = channels[ch]
            for proc in processNames:
                process = channel.processes[proc]
                nuisances = process.nuisances
                if name in nuisances and abs(nuisances[name]) >= 0.001:
                    line.append('%9.3f' % (1. + nuisances[name]))
                else:
                    line.append('-')

        lines.append(line)
        kmax += 1

    for ch in channelNames:
        channel = channels[ch]
        for proc in processNames:
            process = channel.processes[proc]
            count = process.count()
            if count == 0: continue

            line = [ch + '_' + proc + ' gmN ' + str(count)]
            for chproc in [(c, p) for c in channelNames for p in processNames]:
                if chproc == (ch, proc):
                    line.append('%.2e' % (process.rate() / count))
                else:
                    line.append('-')
            
            lines.append(line)
            kmax += 1

    alignColumns(lines[blockStart:])

    kmaxLine.append(str(kmax))

    with open(cardName, 'w') as datacard:
        for line in lines:
            if line == HLINE:
                datacard.write(line + '\n')
            else:
                datacard.write(' '.join(line) + '\n')

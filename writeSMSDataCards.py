from writeDataCard import *

def getSignal(pointName, channels, plotsDir, proc = ''):

    source = ROOT.TFile.Open(plotsDir + '/' + pointName + '.root')

    for ch in channels.values():
        rateHistName = ch.lepton + '/groups/' + ch.histoName + '/' + ch.histoName + '_' + pointName
        rawHistName = ch.lepton + '/components/' + ch.histoName + '/' + ch.histoName + '_' + pointName + '_PhotonAnd' + ch.lepton + '_Raw'
        rate, count = getBinContent(source.Get(rateHistName), ch.xmin, ch.xmax, source.Get(rawHistName))
        if rate < 0.: rate = 0.
        if proc:
            ch.addProcess(proc, rate, count = count, signal = True)
        else:
            ch.setProcess('signal', rate, count = count, signal = True)

    source.Close()


if __name__ == '__main__':
    import sys
    import pickle
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: writeSMSDataCards.py [options] inputName')
    parser.add_option('-d', '--directory', dest = 'outputDir', default = '/afs/cern.ch/user/y/yiiyama/work/datacards')
    parser.add_option('-p', '--plots-directory', dest = 'plotsDir', default = 'rooth://ncmu40//store/countSignal/plots')

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    inputName = args[0]

    channels, nuisances = pickle.load(open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + inputName))

    for mchi in range(200, 810, 10):
        pointName = 'TChiwg_{mchi}'.format(mchi = mchi)
        print pointName

        getSignal(pointName, channels, options.plotsDir)
        writeDataCard(channels, nuisances, options.outputDir + '/' + pointName + '.dat')

    for mglu in range(700, 1350, 50):
        for mchi in range(25, mglu, 50):
            pointName = 'T5wg_{mglu}_{mchi}'.format(mglu = mglu, mchi = mchi)
            print pointName

            getSignal(pointName, channels, options.plotsDir)
            writeDataCard(channels, nuisances, options.outputDir + '/' + pointName + '.dat')


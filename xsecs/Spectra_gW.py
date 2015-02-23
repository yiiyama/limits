import ROOT
import math

ROOT.gROOT.Macro('/afs/cern.ch/user/y/yiiyama/work/rootlogon.C')

points = ['M3_{M3}_M2_{M2}'.format(M3 = M3, M2 = M2) for M3 in range(715, 1565, 50) for M2 in range(205, M3, 50)]

plist = [
    ('gg', ['gg'], 250000),
    ('nc', ['ncp', 'ncm'], 125000)
]

prospino = {}
with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/gen/Spectra_gW/prospino.dat') as source:
    for line in source:
        name, central, upS, downS = line.strip().split()
        up = float(upS)
        down = float(downS)
        prospino[name] = (float(central), math.sqrt(up * up + down * down) / float(central))

with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/Spectra_gW.xsecs', 'w') as outputFile:
    for point in points:
        for name, procs, events in plist:
            if name == 'nc':
                entries = {}

                for proc in procs:
                    source = ROOT.TFile.Open('/store/RA3Ntuples/SusyNtuples/cms538v1p2/PrivateMC/Spectra_gW/susyEvents_' + point + '_' + proc + '.root')
                    tree = source.Get('susyTree')
                    entries[proc] = float(tree.GetEntries())
                    source.Close()

                # Pythia generates chi^0 + chi^+-. Breakdown to positive and negative charginos unknown. Inferred from how many of each survived the gen-level filter.
                # Since BRs chi^0 -> photon and chi^+- -> W^+- -> leptons do not depend on the sign of the chargino, this should give a fairly accurate estimation of the
                # generated numbers.
                total = dict([(p, int(events * entries[p] / (entries['ncp'] + entries['ncm']))) for p in ['ncp', 'ncm']])

            else:
                total = {'gg': events}

            for proc in procs:
                xsec, relErr = prospino[point + '_' + proc]

                outputFile.write('Spectra_gW_%s_%s %.5e %.2e %d\n' % (point, proc, xsec, relErr, total[proc]))
            

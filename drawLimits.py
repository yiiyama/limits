import re
import math
import array
import pickle
import ROOT

ROOT.gROOT.SetBatch(True)
rootlogon = ROOT.gEnv.GetValue("Rint.Logon", "")
if rootlogon:
    ROOT.gROOT.Macro(rootlogon)

ROOT.gErrorIgnoreLevel = 2000

lumi = 19712.
imgform = '.pdf'

def suffix(sig):
    suf = ''
    if sig < 0:
        suf = '_m' + str(abs(sig))
    elif sig == 3:
        suf = '_obs'
    elif sig > 0:
        suf = '_p' + str(abs(sig))

    return suf


def setPalette(suppression):
    ROOT.gStyle.SetNumberContours(99)
    ROOT.TColor.InitializeColors()
    stops = [0.]
    red = [1. - (1. - 0.573) / suppression]
    green = [0.]
    blue = [1.]
    for iC in range(1, 11):
        stops.append(0.1 * iC)
        col = ROOT.gROOT.GetColor(50 + iC * 5)
        red.append(1. - (1. - col.GetRed()) / suppression)
        green.append(1. - (1. - col.GetGreen()) / suppression)
        blue.append(1. - (1. - col.GetBlue()) / suppression)
    
    ROOT.TColor.CreateGradientColorTable(11, array.array('d', stops), array.array('d', red), array.array('d', green), array.array('d', blue), 255)
    

def truncateContour(contour, base):
    x = ROOT.Double()
    y = ROOT.Double()

    iP = contour.GetN() - 1
    while iP >= 0:
        contour.GetPoint(iP, x, y)
        if base.GetBinContent(base.FindBin(x, y)) != 0.:
            break

        iP -= 1
        
    contour.Set(iP + 1)

    iP = 0
    while iP < contour.GetN():
        contour.GetPoint(iP, x, y)
        if base.GetBinContent(base.FindBin(x, y)) != 0.:
            break

        iP += 1

    shift = iP
    for iP in range(contour.GetN() - shift):
        contour.GetPoint(iP + shift, x, y)
        contour.SetPoint(iP, x, y)

    contour.Set(contour.GetN() - shift)


def closeContour(contour, base):
    # currently disabled

    x = ROOT.Double()
    y = ROOT.Double()

    contour.GetPoint(0, x, y)
    xBegin = float(x)
    yBegin = float(y)
    contour.GetPoint(contour.GetN() - 1, x, y)
    xEnd = float(x)
    yEnd = float(y)

    xmin = base.GetXaxis().GetXmin()
    ymin = base.GetYaxis().GetXmin()
    xmax = base.GetXaxis().GetXmax()
    ymax = base.GetYaxis().GetXmax()
    xw = base.GetXaxis().GetBinWidth(1)
    yw = base.GetYaxis().GetBinWidth(1)

    if abs(xBegin - xEnd) < xw and abs(yBegin - yEnd) < yw:
        return

    if xBegin - xmin < xw or yBegin - ymin < yw:
        contour.GetPoint(1, x, y)
        xNext = float(x)
        yNext = float(y)
        xExtr = max(xBegin + (xBegin - xNext) / (yBegin - yNext) * (ymin - yBegin), xmin)
        yExtr = max(yBegin + (yBegin - yNext) / (xBegin - xNext) * (xmin - xBegin), ymin)
        contour.Set(contour.GetN() + 1)
        for iP in range(contour.GetN() - 1, 0, -1):
            contour.GetPoint(iP - 1, x, y)
            contour.SetPoint(iP, x, y)
        x = ROOT.Double(xExtr)
        y = ROOT.Double(yExtr)
        contour.SetPoint(0, x, y)

    if xEnd - xmin < xw or yEnd - ymin < yw:
        contour.GetPoint(contour.GetN() - 2, x, y)
        xNext = float(x)
        yNext = float(y)
        x = ROOT.Double(max(xEnd + (xEnd - xNext) / (yEnd - yNext) * (ymin - yEnd), xmin))
        y = ROOT.Double(max(yEnd + (yEnd - yNext) / (xEnd - xNext) * (xmin - xEnd), ymin))
        contour.Set(contour.GetN() + 1)
        contour.SetPoint(contour.GetN() - 1, x, y)

    if xmax - xBegin < xw or ymax - yBegin < yw:
        contour.GetPoint(1, x, y)
        xNext = float(x)
        yNext = float(y)
        xExtr = min(xBegin + (xBegin - xNext) / (yBegin - yNext) * (ymax - yBegin), xmax)
        yExtr = min(yBegin + (yBegin - yNext) / (xBegin - xNext) * (xmax - xBegin), ymax)
        contour.Set(contour.GetN() + 1)
        for iP in range(contour.GetN() - 1, 0, -1):
            contour.GetPoint(iP - 1, x, y)
            contour.SetPoint(iP, x, y)
        x = ROOT.Double(xExtr)
        y = ROOT.Double(yExtr)
        contour.SetPoint(0, x, y)

    if xmax - xEnd < xw or ymax - yEnd < yw:
        contour.GetPoint(contour.GetN() - 2, x, y)
        xNext = float(x)
        yNext = float(y)
        x = ROOT.Double(max(xEnd + (xEnd - xNext) / (yEnd - yNext) * (ymax - yEnd), xmax))
        y = ROOT.Double(max(yEnd + (yEnd - yNext) / (xEnd - xNext) * (xmax - xEnd), ymax))
        contour.Set(contour.GetN() + 1)
        contour.SetPoint(contour.GetN() - 1, x, y)


def valueAxis(plot, ndim):
    if ndim == 1:
        return plot.GetYaxis()
    else:
        return plot.GetZaxis()


def drawLimits(model, sourceName, plotsDir, pointFormat, titles, order, axisRange, xsecScale = 1., outputName = ''):

    ndim = len(titles) - 1
    if ndim == 1:
        DRAWOPTION = 'CP'
    else:
        DRAWOPTION = 'COLZ'

    if not outputName:
        outputName = model + '_limits'

    chName = {'Electron': 'e#gamma', 'Muon': '#mu#gamma'}

    # Tree with information on one signal point per row. Columns are the name, calculation method, median and +-1/2 sigma signal strength upper bounds of the point.
    inputTree = ROOT.TChain('limitTree')
    inputTree.Add(sourceName)

    vPointName = array.array('c', ' ' * 100)
    vMethod = array.array('c', ' ' * 64)
    vLimObs = array.array('d', [0.])
    vLimMed = array.array('d', [0.])
    vLimM2s = array.array('d', [0.])
    vLimM1s = array.array('d', [0.])
    vLimP1s = array.array('d', [0.])
    vLimP2s = array.array('d', [0.])
            
    inputTree.SetBranchAddress('pointName', vPointName)
    inputTree.SetBranchAddress('method', vMethod)
    inputTree.SetBranchAddress('obs', vLimObs)
    inputTree.SetBranchAddress('med', vLimMed)
    inputTree.SetBranchAddress('m2s', vLimM2s)
    inputTree.SetBranchAddress('m1s', vLimM1s)
    inputTree.SetBranchAddress('p1s', vLimP1s)
    inputTree.SetBranchAddress('p2s', vLimP2s)

    methodPriority = {'asymptotic': 0, 'profileLikelihood': 1, 'fullCLs': 2}

    limits = dict([(i, {}) for i in range(-2, 4)])
   
    centers = [[], []]
    bestMethods = {}
    coords = {}

    iEntry = 0
    while inputTree.GetEntry(iEntry) > 0:
        iEntry += 1
        
        pointStr = vPointName.tostring()
        pointName = pointStr[0:pointStr.find('\0')]
        methodStr = vMethod.tostring()
        methodName = methodStr[0:methodStr.find('\0')]

        matches = re.match(pointFormat, pointName)
        if not matches:
            raise RuntimeError('Invalid point name ' + pointName)

        if pointName in bestMethods:
            if bestMethods[pointName] > methodPriority[methodName]:
                continue

        bestMethods[pointName] = methodPriority[methodName]

        coord = tuple(map(int, [matches.group(iD + 1) for iD in range(ndim)]))
        coords[pointName] = coord

        for iD in range(ndim):
            centers[iD].append(coord[iD])

        limits[3][coord] = vLimObs[0]
        limits[-2][coord] = vLimM2s[0]
        limits[-1][coord] = vLimM1s[0]
        limits[0][coord] = vLimMed[0]
        limits[1][coord] = vLimP1s[0]
        limits[2][coord] = vLimP2s[0]

    asymptotic = (min(bestMethods.values()) == 0)

    xsecs = {}
    xsecErrors = {}
    nEvents = {}
    yields = {'Electron': {}, 'Muon': {}}
       
    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + model + '.pkl', 'rb') as pickleSource:
        datacard = pickle.load(pickleSource)

    for pointName, coord in coords.items():
        processes, genInfo = datacard[pointName]
        
        xsecs[coord] = 0.

        xsecErr2 = 0.
        nEventsW = 0.
        nEventsW2 = 0.

        for info in genInfo.values():
            xsecs[coord] += info.xsec * xsecScale
            xsecErr2 += math.pow(info.relErr * info.xsec * xsecScale, 2.)
            nEventsW += info.nEvents * info.xsec
            nEventsW2 += info.nEvents * info.xsec * info.xsec

        yieldsW = {'Electron': 0., 'Muon': 0.}
        yieldsW2 = {'Electron': 0., 'Muon': 0.}

        for channelName, process in processes.items():
            if channelName.startswith('el'):
                lepton = 'Electron'
            elif channelName.startswith('mu'):
                lepton = 'Muon'
            
            for componentPoint, (rate, count) in process.rates.items():
                yieldsW[lepton] += count * genInfo[componentPoint].xsec
                yieldsW2[lepton] += count * genInfo[componentPoint].xsec * genInfo[componentPoint].xsec

        xsecErrors[coord] = math.sqrt(xsecErr2) / xsecs[coord]

        # effective
        #  . number of events = (sum_proc{N*sigma})^2 / sum_proc{N*sigma^2}
        #  . yield = (sum_proc_ch{y*sigma})^2 / sum_proc_ch{y*sigma^2}

        if nEventsW2 > 0.:
            nEvents[coord] = math.pow(nEventsW, 2.) / nEventsW2
        else:
            nEvents[coord] = 0.

        for lepton in ['Electron', 'Muon']:
            if yieldsW2[lepton] > 0.:
                yields[lepton][coord] = math.pow(yieldsW[lepton], 2.) / yieldsW2[lepton]
            else:
                yields[lepton][coord] = 0.

    # temporary - BR should be part of GenInfo in the future
    glEvents = {'Electron': {}, 'Muon': {}}
    if model == 'Spectra_gW':
        with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/brs/' + model + '.brs') as brSource:
            dump = dict([(lepton, dict([(coord, {}) for coord in coords.values()])) for lepton in ['Electron', 'Muon']])
            for line in brSource:
                componentPoint, nFiltered, nEl, nMu = line.strip().split()
                pointName = componentPoint[:componentPoint.rfind('_')]
                if pointName not in coords: continue
                coord = coords[pointName]

                dump['Electron'][coord][componentPoint] = int(nEl)
                dump['Muon'][coord][componentPoint] = int(nMu)

        for pointName, coord in coords.items():
            genInfo = datacard[pointName][1]
            for lepton in ['Electron', 'Muon']:
                sumW = 0.
                sumW2 = 0.
                for componentPoint, info in genInfo.items():
                    sumW += dump[lepton][coord][componentPoint] * info.xsec
                    sumW2 += dump[lepton][coord][componentPoint] * info.xsec * info.xsec

                glEvents[lepton][coord] = math.pow(sumW, 2.) / sumW2

    else:
        for lepton in ['Electron', 'Muon']:
            for coord in coords.values():
                glEvents[lepton][coord] = nEvents[coord] * 0.107

#    edges = [[], []]
#    for iD in range(ndim):
#        centers[iD] = sorted(list(set(centers[iD])))
#        for iX in range(len(centers[iD]) - 1):
#            edges[iD].append((centers[iD][iX + 1] + centers[iD][iX]) / 2.)
#
#        edges[iD].insert(0, centers[iD][0] - (centers[iD][1] - centers[iD][0]) / 2.)
#        edges[iD].append(centers[iD][-1] + (centers[iD][-1] - centers[iD][-2]) / 2.)
#
#    if ndim == 1:
#        template = ROOT.TH1D('template', ';' + titles[1], len(edges[0]) - 1, array.array('d', edges[0]))
#    else:
#        template = ROOT.TH2D('template', ';' + ';'.join(titles[1:]), len(edges[0]) - 1, array.array('d', edges[0]), len(edges[1]) - 1, array.array('d', edges[1]))

    width = [0.] * ndim
    nbins = [0] * ndim
    for iD in range(ndim):
        centers[iD] = sorted(list(set(centers[iD])))
        width[iD] = float(centers[iD][1] - centers[iD][0])
        nbins[iD] = int(float(centers[iD][-1] - centers[iD][0]) / width[iD]) + 1

# graphical elements

    for x in ['X', 'Y', 'Z']:
        ROOT.gStyle.SetLabelSize(0.04, x)
        ROOT.gStyle.SetTitleSize(0.05, x)
        ROOT.gStyle.SetNdivisions(210, x)

    if ndim == 1:
        ROOT.gStyle.SetTitleOffset(1.2, 'Y')
    else:
        ROOT.gStyle.SetTitleOffset(1.35, 'Y')
    ROOT.gStyle.SetTitleOffset(1.3, 'Z')

    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetHistLineColor(ROOT.kBlack)

    if ndim == 1:
        canvas = ROOT.TCanvas('limits', 'limits', 600, 600)
    else:
        canvas = ROOT.TCanvas('limits', 'limits', 650, 600)

    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetLeftMargin(0.14)
    if ndim == 2:
        canvas.SetRightMargin(0.17)
        
    canvas.SetTicks(1, 1)
        
    if ndim == 1:
        template = ROOT.TH1D('template', ';' + titles[1], nbins[0], centers[0][0] - width[0] / 2., centers[0][-1] + width[0] / 2.)
        template.SetMarkerSize(0)
        template.SetMarkerStyle(1)

    else:
        template = ROOT.TH2D('template', ';' + ';'.join(titles[1:]), nbins[0], centers[0][0] - width[0] / 2., centers[0][-1] + width[0] / 2., nbins[1], centers[1][0] - width[1] / 2., centers[1][-1] + width[1] / 2.)

    template.Sumw2()

    sqrtPave = ROOT.TPaveText()
    sqrtPave.SetBorderSize(0)
    sqrtPave.SetMargin(0.02)
    sqrtPave.SetTextAlign(32)
    sqrtPave.AddText('#sqrt{s} = 8 TeV')
    sqrtPave.SetX1NDC(0.5)
    sqrtPave.SetX2NDC(1. - canvas.GetRightMargin())
    sqrtPave.SetY1NDC(1. - canvas.GetTopMargin())
    sqrtPave.SetY2NDC(1.)

    lumiPave = ROOT.TPaveText()
    lumiPave.SetBorderSize(0)
    lumiPave.SetMargin(0.)
    lumiPave.SetTextAlign(32)
    lumiPave.AddText('#sqrt{s} = 8 TeV, L = %.1f fb^{-1}' % (lumi / 1000.))
    lumiPave.SetX1NDC(0.)
    lumiPave.SetX2NDC(0.95)
    lumiPave.SetY1NDC(1. - canvas.GetTopMargin())
    lumiPave.SetY2NDC(1.)

    cmsPave = ROOT.TPaveText()
    cmsPave.SetBorderSize(0)
    cmsPave.SetMargin(0.)
    cmsPave.SetTextAlign(12)
#    cmsPave.AddText('CMS #font[52]{Preliminary}')
    cmsPave.AddText('CMS')
    cmsPave.SetX1NDC(canvas.GetLeftMargin() + 0.03)
    cmsPave.SetX2NDC(1.)
    cmsPave.SetY1NDC(1. - canvas.GetTopMargin() - 0.08)
    cmsPave.SetY2NDC(1. - canvas.GetTopMargin() - 0.03)

    cmsuPave = ROOT.TPaveText()
    cmsuPave.SetBorderSize(0)
    cmsuPave.SetMargin(0.)
    cmsuPave.SetTextAlign(12)
    cmsuPave.AddText('CMS #font[52]{Unpublished}')
    cmsuPave.SetX1NDC(canvas.GetLeftMargin() + 0.03)
    cmsuPave.SetX2NDC(1.)
    cmsuPave.SetY1NDC(1. - canvas.GetTopMargin() - 0.08)
    cmsuPave.SetY2NDC(1. - canvas.GetTopMargin() - 0.03)

    output = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + outputName + '.root', 'recreate')

    if ndim == 1:
        canvas.SetLogy(True)
    else:
        canvas.SetLogz(True)
        setPalette(1.)

    # BR plot
    BR = []
    BRHisto = {}
    for lept, table in glEvents.items():
        if ndim == 1:
            numer = template.Clone('numer')
            denom = template.Clone('denom')
            for coord, y in table.items():
                bin = numer.FindBin(*coord)
                numer.SetBinContent(bin, y)
                denom.SetBinContent(bin, nEvents[coord])
                f = float(y) / nEvents[coord]
                if f > 0.:
                    BR.append(f)
            
            branching = ROOT.TGraphAsymmErrors(numer, denom)
            branching.SetName('branching_' + lept)
            branching.SetTitle('')
            branching.GetXaxis().SetTitle(template.GetXaxis().GetTitle())
            branching.SetLineColor(ROOT.kBlack)
            branching.SetMarkerColor(ROOT.kBlack)
            branching.SetMarkerStyle(8)
            branching.SetMarkerSize(0.5)

            BRHisto[lept] = numer.Clone('branching_' + lept)
            BRHisto[lept].Divide(denom)
            BRHisto[lept].Scale(100.)

            for iP in range(branching.GetN()):
                branching.GetY()[iP] *= 100.
                branching.SetPointEYhigh(iP, branching.GetErrorYhigh(iP) * 100.)
                branching.SetPointEYlow(iP, branching.GetErrorYlow(iP) * 100.)

            numer.Delete()
            denom.Delete()

            drawOption = 'AP'

        else:
            branching = template.Clone('branching_' + lept)

            for coord, y in table.items():
                try:
                    bin = branching.FindBin(*coord)
                    N = float(nEvents[coord])
                    f = y / N
                    branching.SetBinContent(bin, f * 100.)
                    branching.SetBinError(bin, f * math.sqrt(1. / y + 1. / N) * 100.)
                    if f > 0.:
                        BR.append(f)
                except:
                    pass

            BRHisto[lept] = branching

            drawOption = DRAWOPTION

        valueAxis(branching, ndim).SetTitle('BR(#gamma+#font[12]{l}+X) #times 100 (' + chName[lept] + ' channel)')

        branching.Write()

        if max(BR) / min(BR) < 10.:
            if ndim == 1:
                canvas.SetLogy(False)
            else:
                canvas.SetLogz(False)

        branching.Draw(drawOption)
        ROOT.gPad.Update()

        if ndim == 2:
            branching.GetZaxis().SetTickLength(0.02)
            palette = branching.GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)

        lumiPave.Draw()
        cmsuPave.Draw()
        canvas.Print(plotsDir + '/' + model + '_branching_' + lept + imgform)


    # acceptance plot
    Ae = []
    AeHisto = {}
    for lept, table in yields.items():
        if ndim == 1:
            numer = template.Clone('numer')
            denom = template.Clone('denom')
            for coord, y in table.items():
                bin = numer.FindBin(*coord)
                numer.SetBinContent(bin, y)
                denom.SetBinContent(bin, glEvents[lept][coord])
                f = float(y) / glEvents[lept][coord]
                if f > 0.:
                    Ae.append(f)
            
            acceptance = ROOT.TGraphAsymmErrors(numer, denom)
            acceptance.SetName('acceptance_' + lept)
            acceptance.SetTitle('')
            acceptance.GetXaxis().SetTitle(template.GetXaxis().GetTitle())
            acceptance.SetLineColor(ROOT.kBlack)
            acceptance.SetMarkerColor(ROOT.kBlack)
            acceptance.SetMarkerStyle(8)
            acceptance.SetMarkerSize(0.5)

            AeHisto[lept] = numer.Clone('acceptance_' + lept)
            AeHisto[lept].Divide(denom)

            numer.Delete()
            denom.Delete()

            drawOption = 'AP'

        else:
            acceptance = template.Clone('acceptance_' + lept)

            for coord, y in table.items():
                try:
                    bin = acceptance.FindBin(*coord)
                    N = float(glEvents[lept][coord])
                    f = y / N
                    acceptance.SetBinContent(bin, f)
                    acceptance.SetBinError(bin, f * math.sqrt(1. / y + 1. / N))
                    if f > 0.:
                        Ae.append(f)
                except:
                    pass

            AeHisto[lept] = acceptance

            drawOption = DRAWOPTION

        valueAxis(acceptance, ndim).SetTitle('A #times #varepsilon (' + chName[lept] + ' channel)')

        acceptance.Write()

        if max(Ae) / min(Ae) < 10.:
            if ndim == 1:
                canvas.SetLogy(False)
            else:
                canvas.SetLogz(False)

        acceptance.Draw(drawOption)
        ROOT.gPad.Update()

        if ndim == 2:
            acceptance.GetZaxis().SetTickLength(0.02)
            palette = acceptance.GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)

        lumiPave.Draw()
        cmsuPave.Draw()
        canvas.Print(plotsDir + '/' + model + '_acceptance_' + lept + imgform)


    if ndim == 1:
        canvas.SetLogy(True)
    else:
        canvas.SetLogz(True)


    # cross section plots
    crosssect = {}
    for sig in range(-2, 3):
        crosssect[sig] = template.Clone('crosssect_' + str(sig) + 'sigma')
        title = 'Cross section'
        if sig != 0:
            title += ' %+d #sigma_{%s}' % (sig, order)
        if xsecScale != 1.:
            title += ' (BR=%.2f)' % xsecScale

        valueAxis(crosssect[sig], ndim).SetTitle(title)

        minXS = 1000.
        for coord, xsec in xsecs.items():
            bin = crosssect[sig].FindBin(*coord)
            xs = xsec * (1. + sig * xsecErrors[coord])
            crosssect[sig].SetBinContent(bin, xs)
            if xs < minXS: minXS = xs

        crosssect[sig].SetMinimum(minXS * 0.5)

        crosssect[sig].Write()

        crosssect[sig].Draw(DRAWOPTION)
        ROOT.gPad.Update()

        if ndim == 2:
            crosssect[sig].GetZaxis().SetTickLength(0.02)
            palette = crosssect[sig].GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)

        sqrtPave.Draw()
        canvas.Print(plotsDir + '/' + model + '_crosssect' + suffix(sig) + imgform)

#    canvas.SetLogz(True)

    # expected events
    for lepton in ['Electron', 'Muon']:
        expN = template.Clone('expN_' + lepton)
        valueAxis(expN, ndim).SetTitle('Expected number of events (' + chName[lepton] + ' channel)')

        maxN = 0.
        minN = 1000.
        for coord, xsec in xsecs.items():
            bin = expN.FindBin(*coord)
            n = xsec * BRHisto[lepton].GetBinContent(bin) / 100. * AeHisto[lepton].GetBinContent(bin) * lumi
            expN.SetBinContent(bin, n)
            if n > maxN: maxN = n
            if n < minN: minN = n

        expN.SetMinimum(minN * 0.5)
        expN.SetMaximum(maxN * 1.5)
            
        expN.Write()
        expN.Draw(DRAWOPTION)
        ROOT.gPad.Update()

        if ndim == 2:
            expN.GetZaxis().SetTickLength(0.02)
            palette = expN.GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)

        lumiPave.Draw()
        cmsuPave.Draw()
        canvas.Print(plotsDir + '/' + model + '_expN_' + lepton + imgform)


    # upper limit plots
    upperlim = {}
    for sig in range(-2, 4):
        if sig == 3:
            upperlim[sig] = template.Clone('upperlim_obs')
        else:
            upperlim[sig] = template.Clone('upperlim_' + str(sig) + 'sigma')
        axisTitle = '#sigma_{95%}'
        if sig == 3:
            axisTitle += '^{observed}'
        else:
            axisTitle += '^{expected}'
            if sig != 0:
                quantile = ['2.5', '16', '50', '84', '97.5']
                axisTitle += quantile[sig + 2] + '% quantile'

        if asymptotic: axisTitle += ' asymptotic CL_{s}'
        else: axisTitle += ' CL_{s}'
        axisTitle += ' (pb)'
            
        valueAxis(upperlim[sig], ndim).SetTitle(axisTitle)
        
        for coord, limit in limits[sig].items():
            bin = upperlim[sig].FindBin(*coord)
            upperlim[sig].SetBinContent(bin, limit * xsecs[coord])

        upperlim[sig].Write()

        upperlim[sig].Draw(DRAWOPTION)
        ROOT.gPad.Update()

        if ndim == 2:
            upperlim[sig].GetZaxis().SetTickLength(0.02)
            palette = upperlim[sig].GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)

        lumiPave.Draw()
        cmsuPave.Draw()
        canvas.Print(plotsDir + '/' + model + '_upperlim' + suffix(sig) + imgform)


    if ndim == 2:
        canvas.SetLogz(False)

    # expected exclusion plots
    exclusion = {}
    contours = {}
    for sig in range(-1, 2) + [3]:
        for theo in range(-1, 2):
            if sig != 3 and theo != 0: continue

            key = (sig, theo)

            if sig == 3:
                exclusion[key] = template.Clone('exclusion_obs_' + str(theo) + 'sigmatheo')
            else:
                exclusion[key] = template.Clone('exclusion_expected_' + str(sig) + 'sigmaexp')

            theory = '#sigma_{%s}' % order
            limit = '#sigma_{95%}'

            if theo != 0:
                theory += '^{%+d#sigma}' % theo
            if sig == 3:
                limit += '^{observed}'
            elif sig == 0:
                limit += '^{expected}'
            else:
                quantile = ['2.5', '16', '50', '84', '97.5']
                limit += '^{expected %s%% quantile}' % quantile[sig + 2]

            axisTitle = theory + ' - ' + limit
            if xsecScale != 1.:
                axisTitle += ' (BR=%.2f)' % xsecScale

            if asymptotic: axisTitle += ' asymptotic CL_{s}'
            else: axisTitle += ' CL_{s}'
            axisTitle += ' (pb)'
    
            valueAxis(exclusion[key], ndim).SetTitle(axisTitle)

            exclusion[key].Add(crosssect[theo])
            exclusion[key].Add(upperlim[sig], -1.)
    
            if ndim == 2:
                contours[key] = []

                clevel = array.array('d', [0.])
                contSource = exclusion[key].Clone('contSource')
                contSource.SetContour(1, clevel)
                contSource.Draw('CONT LIST')
                canvas.Update()
                contList = ROOT.gROOT.GetListOfSpecials().FindObject('contours').At(0)
                for contour in contList:
                    cont = contour.Clone()
                    cont.Write(exclusion[key].GetName() + '_contour')
                    contours[key].append(cont)
    
                contSource.Delete()


#    expcol = ROOT.gROOT.GetListOfColors().GetSize()
#    col = ROOT.TColor(expcol, 196./255., 196./255., 196./255.)
    expcol = ROOT.kRed
    obscol = ROOT.kBlack

    if ndim == 1:
        ymax = 0.
        ymin = 100.

        theory = ROOT.TGraph(len(xsecs.items()))
        theoryUp = ROOT.TGraph(len(xsecs.items()))
        theoryDown = ROOT.TGraph(len(xsecs.items()))
        if model == 'TChiwg':
            theoryScaled = ROOT.TGraph(len(xsecs.items()))
            theoryScaledUp = ROOT.TGraph(len(xsecs.items()))
            theoryScaledDown = ROOT.TGraph(len(xsecs.items()))

        iP = 0
        for coord in sorted(xsecs.keys()):
            yUp = xsecs[coord] * (1. + xsecErrors[coord])
            yDown = xsecs[coord] * (1. - xsecErrors[coord])

            theory.SetPoint(iP, coord[0], xsecs[coord])
            theoryUp.SetPoint(iP, coord[0], yUp)
            theoryDown.SetPoint(iP, coord[0], yDown)
            if model == 'TChiwg':
                theoryScaled.SetPoint(iP, coord[0], xsecs[coord] * 0.22)
                theoryScaledUp.SetPoint(iP, coord[0], yUp * 0.22)
                theoryScaledDown.SetPoint(iP, coord[0], yDown * 0.22)

            iP += 1

            if yUp > ymax: ymax = yUp
            if yDown < ymin and yDown > 0.: ymin = yDown

        expected = ROOT.TGraph(upperlim[0].GetNbinsX())
        exp1s = ROOT.TGraphAsymmErrors(expected.GetN())
        exp2s = ROOT.TGraphAsymmErrors(expected.GetN())
        observed = ROOT.TGraph(upperlim[3].GetNbinsX())
        for iP in range(expected.GetN()):
            yExp = upperlim[0].GetBinContent(iP + 1)
            point = (upperlim[0].GetXaxis().GetBinCenter(iP + 1), yExp)
            expected.SetPoint(iP, *point)

            yUp = upperlim[2].GetBinContent(iP + 1)
            yDown = upperlim[-2].GetBinContent(iP + 1)

            exp1s.SetPoint(iP, *point)
            exp1s.SetPointEYhigh(iP, upperlim[1].GetBinContent(iP + 1) - yExp)
            exp1s.SetPointEYlow(iP, yExp - upperlim[-1].GetBinContent(iP + 1))

            exp2s.SetPoint(iP, *point)
            exp2s.SetPointEYhigh(iP, yUp - yExp)
            exp2s.SetPointEYlow(iP, yExp - yDown)

            yObs = upperlim[3].GetBinContent(iP + 1)
            observed.SetPoint(iP, upperlim[3].GetXaxis().GetBinCenter(iP + 1), yObs)

            if yUp > ymax: ymax = yUp
            if yDown < ymin and yDown > 0.: ymin = yDown

        theory.SetLineColor(ROOT.kBlue)
        theory.SetLineStyle(ROOT.kSolid)
        theory.SetLineWidth(2)
        theoryUp.SetLineColor(ROOT.kBlue)
        theoryUp.SetLineStyle(ROOT.kDashed)
        theoryUp.SetLineWidth(2)
        theoryDown.SetLineColor(ROOT.kBlue)
        theoryDown.SetLineStyle(ROOT.kDashed)
        theoryDown.SetLineWidth(2)
        if model == 'TChiwg':
            theoryScaled.SetLineColor(ROOT.kOrange + 10)
            theoryScaled.SetLineStyle(ROOT.kSolid)
            theoryScaled.SetLineWidth(2)
            theoryScaledUp.SetLineColor(ROOT.kOrange + 10)
            theoryScaledUp.SetLineStyle(ROOT.kDashed)
            theoryScaledUp.SetLineWidth(2)
            theoryScaledDown.SetLineColor(ROOT.kOrange + 10)
            theoryScaledDown.SetLineStyle(ROOT.kDashed)
            theoryScaledDown.SetLineWidth(2)

        expected.SetLineColor(ROOT.kBlack)
        expected.SetLineStyle(ROOT.kDashed)
        expected.SetLineWidth(2)
        exp1s.SetFillColor(ROOT.kGreen)
        exp1s.SetLineWidth(0)
        exp2s.SetFillColor(ROOT.kYellow)
        exp2s.SetLineWidth(0)
        observed.SetLineColor(obscol)
        observed.SetLineStyle(ROOT.kSolid)
        observed.SetLineWidth(2)
        observed.SetMarkerStyle(21)
        observed.SetMarkerSize(0.5)

        exp2s.Draw('AC3')
        exp1s.Draw('C3')
        theory.Draw('C')
        theoryUp.Draw('C')
        theoryDown.Draw('C')
        if model == 'TChiwg':
            theoryScaled.Draw('C')
            theoryScaledUp.Draw('C')
            theoryScaledDown.Draw('C')

        expected.Draw('C')
        observed.Draw('CP')

        if axisRange:
            rangeMin, rangeMax = axisRange
        else:
            rangeMin = ymin * 0.5
            rangeMax = math.exp(math.log(ymax / rangeMin) * 0.9 / 0.6 + math.log(rangeMin))

        exp2s.SetTitle('')
        exp2s.GetXaxis().SetTitle(upperlim[0].GetXaxis().GetTitle())
        exp2s.GetYaxis().SetTitle('')
        exp2s.GetYaxis().SetLabelSize(0)
        exp2s.GetYaxis().SetRangeUser(rangeMin, rangeMax)

        canvas.Update()

        bottom = canvas.GetBottomMargin()
        top = 1. - canvas.GetTopMargin()
        xmin = exp2s.GetXaxis().GetXmin()
        ymin = canvas.GetUymin()
        ymax = canvas.GetUymax()
        ymincoord = math.pow(10., ymin)
        newymaxcoord = math.pow(10., ymin + (ymax - ymin) / (top - bottom) * (0.7 - bottom))

        yaxis = ROOT.TGaxis(xmin, ymincoord, xmin, newymaxcoord, ymincoord, newymaxcoord, exp2s.GetYaxis().GetNdivisions(), 'GSB')
        yaxis.SetTickSize(0.)
        yaxis.SetLabelSize(ROOT.gStyle.GetLabelSize('Y'))
        yaxis.SetLabelFont(ROOT.gStyle.GetLabelFont('Y'))
        yaxis.SetTitleSize(ROOT.gStyle.GetTitleSize('Y'))
        yaxis.SetTitleFont(ROOT.gStyle.GetTitleFont('Y'))
        yaxis.SetTitleOffset(ROOT.gStyle.GetTitleOffset('Y'))
        yaxis.SetTitle('#sigma (pb)')
        yaxis.Draw()

    else:
        canvas.SetLogz(True)
        canvas.Clear()

        background = template.Clone('background')
        for iX in range(1, background.GetNbinsX() + 1):
            for iY in range(1, background.GetNbinsY() + 1):
                background.SetBinContent(iX, iY, upperlim[3].GetBinContent(iX, iY))

        setPalette(1.5)

        background.GetZaxis().SetTickLength(0.02)
        background.GetZaxis().SetTitle('95% CL cross section upper limit (pb)')
        background.Draw(DRAWOPTION)

        canvas.Update()

        maxY = 0.
        for contour in contours[(0, 0)]:
            truncateContour(contour, crosssect[0])
#            closeContour(contour, crosssect[0])
            contour.SetLineColor(expcol)
            contour.SetLineWidth(4)
            contour.Draw('C')
            
            for iP in range(contour.GetN()):
                y = contour.GetY()[iP]
                if y > maxY: maxY = y
    
        for sig in [-1, 1]:
            for contour in contours[(sig, 0)]:
                truncateContour(contour, crosssect[0])
#                closeContour(contour, crosssect[0])
                contour.SetLineColor(expcol)
                contour.SetLineWidth(2)
                contour.SetLineStyle(ROOT.kDashed)
                contour.Draw('C')
    
        for contour in contours[(3, 0)]:
            truncateContour(contour, crosssect[0])
#            closeContour(contour, crosssect[0])
            contour.SetLineColor(obscol)
            contour.SetLineWidth(4)
            contour.Draw('C')
            
            for iP in range(contour.GetN()):
                y = contour.GetY()[iP]
                if y > maxY: maxY = y
    
        for theo in [-1, 1]:
            for contour in contours[(3, theo)]:
                truncateContour(contour, crosssect[0])
#                closeContour(contour, crosssect[0])
                contour.SetLineColor(obscol)
                contour.SetLineWidth(2)
                contour.SetLineStyle(ROOT.kDashed)
                contour.Draw('C')
    
        yUp = (maxY - background.GetYaxis().GetXmin()) * 0.9 / 0.6 + background.GetYaxis().GetXmin()
        if yUp <= background.GetYaxis().GetXmax():
            background.GetYaxis().SetRangeUser(background.GetYaxis().GetXmin(), yUp)
        else:
            n = background.GetNbinsY()
            w = background.GetYaxis().GetBinWidth(1)
            n += int((yUp - background.GetYaxis().GetXmax()) / w)
            background.SetBins(background.GetNbinsX(), background.GetXaxis().GetXmin(), background.GetXaxis().GetXmax(), n, background.GetYaxis().GetXmin(), background.GetYaxis().GetXmin() + n * w)

        if axisRange:
            minZ, maxZ = axisRange
        else:
            minZ = 1000.
            maxZ = -1.
            for iX in range(1, background.GetNbinsX() + 1):
                for iY in range(1, background.GetNbinsY() + 1):
                    cont = background.GetBinContent(iX, iY)
                    if cont != 0. and cont < minZ:
                        minZ = cont
                    if cont > maxZ:
                        maxZ = cont

            minZ *= 0.8
            maxZ *= 1.2
    
        background.SetMinimum(minZ * 0.8)
        background.SetMaximum(maxZ * 0.8)

        palette = background.GetListOfFunctions().FindObject('palette')
        palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)
        palette.GetAxis().SetTitle(background.GetZaxis().GetTitle())

    header = ROOT.TPad()
    canvas.cd()
    header.Draw()
    header.cd()

    if type(titles[0]) is not tuple:
        titleLines = (titles[0],)
    else:
        titleLines = titles[0]

    nTitleLines = len(titleLines)
    if ndim == 1:
        nTotalLines = nTitleLines + 4.
    else:
        nTotalLines = nTitleLines + 3.
        
    x1 = canvas.GetLeftMargin()
    y2 = 1. - canvas.GetTopMargin()
    y1 = y2 - 0.045 * nTotalLines
    x2 = 1. - canvas.GetRightMargin()
    if ndim == 2:
        x1 += 0.01
        y2 -= 0.025
        y1 -= 0.025

    header.SetPad(x1, y1, x2, y2)
    if ndim == 1:
        header.SetFillStyle(1001)
        headerBack = ROOT.TBox(0., 0., 1., 1.)
        headerBack.SetLineWidth(1)
        headerBack.SetLineColor(ROOT.kBlack)
        headerBack.Draw()
    else:
        header.SetFillStyle(0)

    textSize = 0.72 / (nTotalLines - 1.)

    title = ROOT.TPaveText()
    title.SetMargin(0.)
    title.SetTextSize(textSize)
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextAlign(12)
    title.SetX1NDC(0.)
    title.SetX2NDC(1.)
    title.SetY1NDC(1. - nTitleLines / nTotalLines)
    title.SetY2NDC(1.)
    y = 0.93 - 0.495 / nTitleLines
    for line in titleLines:
        title.AddText(0.02, y, line)
        y -= 0.93 / nTitleLines

    title.Draw()

    legend = ROOT.TPaveText()
    legend.SetMargin(0.)
    legend.SetTextSize(textSize)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12)
    legend.SetX1NDC(0.)
    legend.SetX2NDC(1.)
    legend.SetY1NDC(0.)
    legend.SetY2NDC(1. - nTitleLines / nTotalLines)

    if ndim == 1:
        y = 0.82
        text = '95%'
        if asymptotic: text += ' asymptotic'
        text += ' CL upper limits'
        if xsecScale != 1.:
            text += ' (BR=%.2f)' % xsecScale
        legend.AddText(0.02, y, text)

        y = 0.61
        line = legend.AddLine(0.02, y, 0.08, y)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(obscol)
        box = legend.AddBox(0.046, y - 0.03, 0.054, y + 0.03)
        box.SetFillStyle(1001)
        box.SetFillColor(obscol)
        box.SetLineWidth(0)
        box.SetLineColor(obscol)
        legend.AddText(0.1, y, 'Observed')

        y = 0.36
        box = legend.AddBox(0.02, y - 0.06, 0.08, y + 0.06)
        box.SetFillStyle(1001)
        box.SetFillColor(ROOT.kGreen)
        box.SetLineWidth(0)
        box.SetLineColor(ROOT.kGreen)
        line = legend.AddLine(0.02, y, 0.08, y)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlack)
        legend.AddText(0.1, y, 'Expected #pm 1 #sigma')

        y = 0.11
        box = legend.AddBox(0.02, y - 0.06, 0.08, y + 0.06)
        box.SetFillStyle(1001)
        box.SetFillColor(ROOT.kYellow)
        box.SetLineWidth(0)
        box.SetLineColor(ROOT.kYellow)
        line = legend.AddLine(0.02, y, 0.08, y)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlack)
        legend.AddText(0.1, y, 'Expected #pm 2 #sigma')

        y = 0.82
        line = legend.AddLine(0.52, y + 0.06, 0.58, y + 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(ROOT.kBlue)
        line = legend.AddLine(0.52, y, 0.58, y)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(ROOT.kBlue)
        line = legend.AddLine(0.52, y - 0.06, 0.58, y - 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(ROOT.kBlue)
        legend.AddText(0.6, y, order + ' #pm 1 #sigma')

        if model == 'TChiwg':
            y = 0.58
            line = legend.AddLine(0.52, y + 0.06, 0.58, y + 0.06)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kDashed)
            line.SetLineColor(ROOT.kOrange + 10)
            line = legend.AddLine(0.52, y, 0.58, y)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kSolid)
            line.SetLineColor(ROOT.kOrange + 10)
            line = legend.AddLine(0.52, y - 0.06, 0.58, y - 0.06)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kDashed)
            line.SetLineColor(ROOT.kOrange + 10)
            legend.AddText(0.6, y, ' #times sin^{2}#theta_{W}')

    else:
        y = 0.83
        text = '%s exclusion' % order
        if xsecScale != 1.:
            text += ' (BR=%.2f)' % xsecScale
        legend.AddText(0.02, y, text)

        y = 0.51
        line = legend.AddLine(0.02, y + 0.06, 0.14, y + 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(obscol)
        line = legend.AddLine(0.02, y, 0.14, y)
        line.SetLineWidth(4)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(obscol)
        line = legend.AddLine(0.02, y - 0.06, 0.14, y - 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(obscol)
        legend.AddText(0.16, y - 0.02, 'Observed #pm 1 #sigma (theory)')

        y = 0.18
        line = legend.AddLine(0.02, y + 0.06, 0.14, y + 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(expcol)
        line = legend.AddLine(0.02, y, 0.14, y)
        line.SetLineWidth(4)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(expcol)
        line = legend.AddLine(0.02, y - 0.06, 0.14, y - 0.06)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(expcol)
        legend.AddText(0.16, y - 0.02, 'Expected #pm 1 #sigma (exp.)')

    legend.Draw()

    canvas.cd()

    cmsPave.SetX1NDC(canvas.GetLeftMargin())
    cmsPave.SetX2NDC(1.)
    cmsPave.SetY1NDC(1. - canvas.GetTopMargin())
    cmsPave.SetY2NDC(1.)
    cmsPave.SetTextAlign(12)

    lumiPave.Draw()
    cmsPave.Draw()
    canvas.Print(plotsDir + '/' + model + '_exclusion' + imgform)

    header.Delete()
        

if __name__ == '__main__':

    import sys
    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: writeDataCard.py [options] outputName')
    parser.add_option('-s', '--source', dest = 'sourceName', default = '', help = 'source root file')
    parser.add_option('-x', '--scale-xsec', dest = 'xsecScale', type = 'float', default = 1., help = 'cross section scale factor')
    parser.add_option('-o', '--output-name', dest = 'outputName', default = '', help = 'output ROOT file name')

    options, args = parser.parse_args()

    model = args[0]
    if options.sourceName:
        sourceName = options.sourceName
    else:
        sourceName = '/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + model + '.root'

    plotsDir = args[1]

    order = 'NLO+NLL'

    if model == 'T5wg':
        form = 'T5wg_([0-9]+)_([0-9]+)'
#        titles = (('pp #rightarrow #tilde{g}#kern[0.2]{#tilde{g}}, #tilde{g} #rightarrow q#kern[0.2]{#bar{q}}#kern[0.2]{#tilde{#chi}^{0/#pm}}', '#tilde{#chi}^{#pm} #rightarrow W^{#pm}#tilde{G}, #tilde{#chi}^{0} #rightarrow #gamma#tilde{G}'), 'm_{#tilde{g}} (GeV)', 'm_{#tilde{#chi}} (GeV)')
        titles = ('pp #rightarrow #tilde{g}#kern[0.2]{#tilde{g}}, #tilde{g} #rightarrow q#kern[0.2]{#bar{q}}#kern[0.2]{#tilde{#chi}^{0/#pm}}, #tilde{#chi}^{#pm} #rightarrow W^{#pm}#tilde{G}, #tilde{#chi}^{0} #rightarrow #gamma#tilde{G}', 'm_{#tilde{g}} (GeV)', 'm_{#tilde{#chi}} (GeV)')
        axisRange = (1.e-3, 5.e-2)
    elif model == 'TChiwg':
        form = 'TChiwg_([0-9]+)'
        titles = ('pp #rightarrow #tilde{#chi}^{#pm}#kern[0.2]{#tilde{#chi}^{0}}, #tilde{#chi}^{#pm} #rightarrow W^{#pm}#tilde{G}, #tilde{#chi}^{0} #rightarrow #gamma#tilde{G}', 'm_{#tilde{#chi}} (GeV)')
        axisRange = None
    elif model == 'T5wg+TChiwg':
        form = 'T5wg\\+TChiwg_([0-9]+)_([0-9]+)'
        titles = ('Simplified GMSB, #tilde{#chi}^{#pm} #rightarrow W^{#pm}, #tilde{#chi}^{0} #rightarrow #gamma', 'm_{#tilde{g}} (GeV)', 'm_{#tilde{#chi}} (GeV)')
        axisRange = None
    elif model == 'Spectra_gW':
        form = 'Spectra_gW_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP', 'gluino mass (GeV)', 'wino mass (GeV)')
        axisRange = (1.e-2, 2.e-1)
    elif model == 'Spectra_gW_nc':
        form = 'Spectra_gW_nc_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP EWK only', 'gluino mass (GeV)', 'wino mass (GeV)')
        axisRange = None
    elif model == 'Spectra_gW_gg':
        form = 'Spectra_gW_gg_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP strong only', 'gluino mass (GeV)', 'wino mass (GeV)')
        axisRange = None

    drawLimits(model, sourceName, plotsDir, form, titles, order, axisRange, xsecScale = options.xsecScale, outputName = options.outputName)

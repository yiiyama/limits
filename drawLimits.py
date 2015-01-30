import re
import math
import array
import pickle
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000
ROOT.gStyle.SetOptStat(0)

TITLEBASE = 'CMS Preliminary, 19.7 fb^{-1}, #sqrt{s} = 8 TeV'

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


def drawLimits(model, sourceName, plotsDir, pointFormat, titles, xsecScale = 1., outputName = ''):

    ndim = len(titles) - 1
    if ndim == 1:
        DRAWOPTION = 'CP'
    else:
        DRAWOPTION = 'COLZ'

    if not outputName:
        outputName = model + '_limits'

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

    if ndim == 1:
        template = ROOT.TH1D('template', ';' + titles[1], nbins[0], centers[0][0] - width[0] / 2., centers[0][-1] + width[0] / 2.)
        template.SetLineWidth(2)
        template.SetMarkerSize(0)
        template.SetMarkerStyle(1)
        template.SetLineColor(ROOT.kBlack)

    else:
        template = ROOT.TH2D('template', ';' + ';'.join(titles[1:]), nbins[0], centers[0][0] - width[0] / 2., centers[0][-1] + width[0] / 2., nbins[1], centers[1][0] - width[1] / 2., centers[1][-1] + width[1] / 2.)

    template.GetXaxis().SetTitleOffset(1.3)
    template.GetYaxis().SetTitleOffset(1.6)
    template.GetXaxis().SetLabelSize(0.03)
    template.GetYaxis().SetLabelSize(0.03)

    canvas = ROOT.TCanvas('limits', 'limits', 800, 800)
    canvas.SetLeftMargin(0.13)
    if ndim == 1:
        canvas.SetRightMargin(0.1)
    else:
        canvas.SetRightMargin(0.17)

    output = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + outputName + '.root', 'recreate')

    if ndim == 1:
        canvas.SetLogy(True)
    else:
        canvas.SetLogz(True)
        setPalette(1.)

    # acceptance plot
    Ae = []
    for lept, table in yields.items():
        if ndim == 1:
            numer = template.Clone('numer')
            denom = template.Clone('denom')
            for coord, y in table.items():
                bin = numer.FindBin(*coord)
                numer.SetBinContent(bin, y)
                denom.SetBinContent(bin, nEvents[coord])
                f = float(y) / nEvents[coord]
                if f > 0.:
                    Ae.append(f)
            
            acceptance = ROOT.TGraphAsymmErrors(numer, denom)
            acceptance.GetXaxis().SetTitle(template.GetXaxis().GetTitle())
            acceptance.SetLineColor(ROOT.kBlack)
            acceptance.SetMarkerColor(ROOT.kBlack)
            acceptance.SetMarkerStyle(8)
            acceptance.SetMarkerSize(0.5)

            numer.Delete()
            denom.Delete()

            drawOption = 'AP'

        else:
            acceptance = template.Clone('acceptance_' + lept)

            for coord, y in table.items():
                try:
                    bin = acceptance.FindBin(*coord)
                    N = float(nEvents[coord])
                    f = y / N
                    acceptance.SetBinContent(bin, f)
                    acceptance.SetBinError(bin, f * math.sqrt(1. / y + 1. / N))
                    if f > 0.:
                        Ae.append(f)
                except:
                    pass

            drawOption = DRAWOPTION

        acceptance.SetTitle('A #times #epsilon (' + lept + ' channel)')

        acceptance.Write()

        if max(Ae) / min(Ae) < 10.:
            if ndim == 1:
                canvas.SetLogy(False)
            else:
                canvas.SetLogz(False)

        acceptance.Draw(drawOption)
        ROOT.gPad.Update()

        if ndim == 2:
            zaxis = acceptance.GetZaxis()
            zaxis.SetTickLength(0.02)
            palette = acceptance.GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)
                palette.SetTitleOffset(1.5)
                palette.SetLabelSize(0.03)

        canvas.Print(plotsDir + '/' + model + '_acceptance_' + lept + '.pdf')


    if ndim == 1:
        canvas.SetLogy(True)
    else:
        canvas.SetLogz(True)


    # cross section plots
    crosssect = {}
    for sig in range(-2, 3):
        crosssect[sig] = template.Clone('crosssect_' + str(sig) + 'sigma')
        title = 'Cross Section'
        if sig != 0:
            title += ' %+d #sigma_{theory}' % sig
        if xsecScale != 1.:
            title += ' (BR=%.2f)' % xsecScale
        crosssect[sig].SetTitle(title)

        axisTitle = '#sigma (pb)'
        if ndim == 1:
            crosssect[sig].GetYaxis().SetTitle(axisTitle)
        else:
            crosssect[sig].GetZaxis().SetTitle(axisTitle)

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
            zaxis = crosssect[sig].GetZaxis()
            zaxis.SetTickLength(0.02)
            palette = crosssect[sig].GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)
                palette.SetTitleOffset(1.5)
                palette.SetLabelSize(0.03)

        canvas.Print(plotsDir + '/' + model + '_crosssect' + suffix(sig) + '.pdf')

#    canvas.SetLogz(True)

    # upper limit plots
    upperlim = {}
    for sig in range(-2, 4):
        upperlim[sig] = template.Clone('upperlim_' + str(sig) + 'sigma')
        title = TITLEBASE
        axisTitle = '#sigma_{95%}'
        if sig == 3:
            axisTitle += '^{observed}'
        elif sig == 0:
            title += ' Expected'
            axisTitle += '^{expected}'
        else:
            title += ' Expected %+d #sigma_{experiment}' % sig
            axisTitle += '^{expected}'

        upperlim[sig].SetTitle(title)

        if asymptotic: axisTitle += ' asymptotic CL_{s}'
        else: axisTitle += ' CL_{s}'
        axisTitle += ' (pb)'

        if ndim == 1:
            upperlim[sig].GetYaxis().SetTitle(axisTitle)
        else:
            upperlim[sig].GetZaxis().SetTitle(axisTitle)
        
        for coord, limit in limits[sig].items():
            bin = upperlim[sig].FindBin(*coord)
            upperlim[sig].SetBinContent(bin, limit * xsecs[coord])

        upperlim[sig].Write()

        upperlim[sig].Draw(DRAWOPTION)
        ROOT.gPad.Update()

        if ndim == 2:
            zaxis = upperlim[sig].GetZaxis()
            zaxis.SetTickLength(0.02)
            palette = upperlim[sig].GetListOfFunctions().FindObject('palette')
            if palette:
                palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)
                palette.SetTitleOffset(1.5)
                palette.SetLabelSize(0.03)

        canvas.Print(plotsDir + '/' + model + '_upperlim' + suffix(sig) + '.pdf')


    if ndim == 2:
        canvas.SetLogz(False)

    # expected exclusion plots
    exclusion = {}
    contours = {}
    for sig in range(-1, 2) + [3]:
        for theo in range(-1, 2):
            if sig != 3 and theo != 0: continue

            key = (sig, theo)

            exclusion[key] = template.Clone('exclusion_' + str(sig) + 'sigma')
            title = 'Exclusion'
            if sig == 3:
                if theo == 0:
                    title += ' (Observed)'
                else:
                    title += ' (Observed %+d #sigma_{theory})' % theo
            elif sig == 0:
                title += ' (Expected)'
            else:
                title += ' (Expected %+d #sigma_{experiment})' % sig
            if xsecScale != 1.:
                title += ' (BR=%.2f)' % xsecScale
    
            exclusion[key].SetTitle(title)

            axisTitle = '#sigma - #sigma_{95%CL_{s}}^{up} (pb)'
            if ndim == 1:
                exclusion[key].GetYaxis().SetTitle(axisTitle)
            else:
                exclusion[key].GetZaxis().SetTitle(axisTitle)

            exclusion[key].Add(crosssect[theo])
            exclusion[key].Add(upperlim[sig], -1.)
    
            exclusion[key].Write()
    
            if ndim == 2:
                contours[key] = []

                clevel = array.array('d', [0.])
                contSource = exclusion[key].Clone('contSource')
                contSource.SetContour(1, clevel)
                contSource.Draw('CONT LIST')
                canvas.Update()
                contList = ROOT.gROOT.GetListOfSpecials().FindObject('contours').At(0)
                for contour in contList:
                    contours[key].append(contour.Clone())
    
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

        rangeMin = ymin * 0.5
        rangeMax = math.exp(math.log(ymax / rangeMin) * 0.9 / 0.6 + math.log(rangeMin))

        exp2s.SetTitle(TITLEBASE)
        exp2s.GetXaxis().SetTitle(upperlim[0].GetXaxis().GetTitle())
        exp2s.GetYaxis().SetTitle('#sigma (pb)')
        exp2s.GetYaxis().SetRangeUser(rangeMin, rangeMax)
        exp2s.GetYaxis().SetTitleOffset(1.5)

    else:
        canvas.SetLogz(True)

        background = upperlim[3].Clone('background')
        setPalette(1.5)
    
        background.Draw(DRAWOPTION)

        zaxis = background.GetZaxis()
        zaxis.SetTickLength(0.02)
        palette = background.GetListOfFunctions().FindObject('palette')
        palette.SetX2NDC(palette.GetX1NDC() + (palette.GetX2NDC() - palette.GetX1NDC()) * 0.7)
        palette.SetTitleOffset(1.5)
        palette.SetLabelSize(0.03)

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
    
        minZ = 1000.
        for iX in range(1, background.GetNbinsX() + 1):
            for iY in range(1, background.GetNbinsY() + 1):
                cont = background.GetBinContent(iX, iY)
                if cont != 0. and cont < minZ:
                    minZ = cont
    
        background.SetMinimum(minZ * 0.8)

    header = ROOT.TPaveText()
    header.SetOption('')
    header.SetX1NDC(0.13)
    if ndim == 1:
        header.SetX2NDC(0.9)
    else:
        header.SetX2NDC(0.83)
    header.SetY2NDC(0.9)
    header.SetY1NDC(0.7)
    header.ConvertNDCtoPad()
    header.SetTextFont(62)
    header.SetTextSize(0.03)
    header.SetFillStyle(1001)
    header.SetFillColor(ROOT.kWhite)
    header.SetLineColor(ROOT.kBlack)
    header.SetLineWidth(1)
    header.SetBorderSize(1)

    if ndim == 1:
        header.AddText(0.02, 0.84, titles[0])
        text = '95%'
        if asymptotic: text += ' asymptotic'
        text += ' CL_{s} cross section upper limits'
        if xsecScale != 1.:
            text += ' (BR=%.2f)' % xsecScale
        header.AddText(0.02, 0.64, text)

        line = header.AddLine(0.02, 0.48, 0.08, 0.48)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(obscol)
        box = header.AddBox(0.046, 0.46, 0.054, 0.5)
        box.SetFillColor(obscol)
        box.SetLineWidth(0)
        header.AddText(0.1, 0.44, 'Observed')

        box = header.AddBox(0.02, 0.24, 0.08, 0.36)
        box.SetFillColor(ROOT.kGreen)
        box.SetLineWidth(0)
        line = header.AddLine(0.02, 0.3, 0.08, 0.3)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlack)
        header.AddText(0.1, 0.24, 'Expected (68%)')

        box = header.AddBox(0.02, 0.06, 0.08, 0.18)
        box.SetFillColor(ROOT.kYellow)
        box.SetLineWidth(0)
        line = header.AddLine(0.02, 0.12, 0.08, 0.12)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlack)
        header.AddText(0.1, 0.04, 'Expected (95%)')

        line = header.AddLine(0.52, 0.52, 0.58, 0.52)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(ROOT.kBlue)
        line = header.AddLine(0.52, 0.48, 0.58, 0.48)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(ROOT.kBlue)
        line = header.AddLine(0.52, 0.44, 0.58, 0.44)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(ROOT.kBlue)
        header.AddText(0.6, 0.44, 'Theory #pm 1 #sigma')

        if model == 'TChiwg':
            line = header.AddLine(0.52, 0.34, 0.58, 0.34)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kDashed)
            line.SetLineColor(ROOT.kOrange + 10)
            line = header.AddLine(0.52, 0.3, 0.58, 0.3)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kSolid)
            line.SetLineColor(ROOT.kOrange + 10)
            line = header.AddLine(0.52, 0.26, 0.58, 0.26)
            line.SetLineWidth(2)
            line.SetLineStyle(ROOT.kDashed)
            line.SetLineColor(ROOT.kOrange + 10)
            header.AddText(0.6, 0.26, ' #times sin^{2}#theta_{W}')

    else:
        header.AddText(0.02, 0.84, titles[0])
        text = '95%'
        if asymptotic: text += ' asymptotic'
        text += ' CL_{s} exclusion'
        if xsecScale != 1.:
            text += ' (BR=%.2f)' % xsecScale
        header.AddText(0.02, 0.58, text)

        line = header.AddLine(0.02, 0.43, 0.08, 0.43)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(obscol)
        line = header.AddLine(0.02, 0.39, 0.08, 0.39)
        line.SetLineWidth(4)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(obscol)
        line = header.AddLine(0.02, 0.35, 0.08, 0.35)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(obscol)
        header.AddText(0.1, 0.33, 'Observed #pm 1 #sigma_{theory}')

        line = header.AddLine(0.02, 0.18, 0.08, 0.18)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(expcol)
        line = header.AddLine(0.02, 0.14, 0.08, 0.14)
        line.SetLineWidth(4)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineColor(expcol)
        line = header.AddLine(0.02, 0.1, 0.08, 0.1)
        line.SetLineWidth(2)
        line.SetLineStyle(ROOT.kDashed)
        line.SetLineColor(expcol)
        header.AddText(0.1, 0.08, 'Expected #pm 1 #sigma_{experiment}')

    header.Draw()
    header.SetDrawOption(' ')
    
    canvas.Print(plotsDir + '/' + model + '_exclusion.pdf')
        

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

    if model == 'T5wg':
        form = 'T5wg_([0-9]+)_([0-9]+)'
        titles = ('pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{0/#pm}, #tilde{#chi}^{#pm} #rightarrow W^{#pm}, #tilde{#chi}^{0} #rightarrow #gamma NLO+NLL', 'M_{#tilde{g}} (GeV)', 'M_{#tilde{#chi}} (GeV)')
    elif model == 'TChiwg':
        form = 'TChiwg_([0-9]+)'
        titles = ('pp #rightarrow #tilde{#chi}^{#pm} #tilde{#chi}^{0}, #tilde{#chi}^{#pm} #rightarrow W^{#pm}, #tilde{#chi}^{0} #rightarrow #gamma NLO+NLL', 'M_{#tilde{#chi}} (GeV)')
    elif model == 'T5wg+TChiwg':
        form = 'T5wg\\+TChiwg_([0-9]+)_([0-9]+)'
        titles = ('Simplified GMSB, #tilde{#chi}^{#pm} #rightarrow W^{#pm}, #tilde{#chi}^{0} #rightarrow #gamma NLO+NLL', 'M_{#tilde{g}} (GeV)', 'M_{#tilde{#chi}} (GeV)')
    elif model == 'Spectra_gW':
        form = 'Spectra_gW_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP NLO+NLL', 'M_{3} (GeV)', 'M_{2} (GeV)')
    elif model == 'Spectra_gW_nc':
        form = 'Spectra_gW_nc_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP EWK only NLO+NLL', 'M_{3} (GeV)', 'M_{2} (GeV)')
    elif model == 'Spectra_gW_gg':
        form = 'Spectra_gW_gg_M3_([0-9]+)_M2_([0-9]+)'
        titles = ('GGM Wino-like NLSP strong only NLO+NLL', 'M_{3} (GeV)', 'M_{2} (GeV)')

    drawLimits(model, sourceName, plotsDir, form, titles, xsecScale = options.xsecScale, outputName = options.outputName)

import re
import math
import array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000
ROOT.gStyle.SetOptStat(0)

PLOTSDIR = '/afs/cern.ch/user/y/yiiyama/www/plots/limits'

def suffix(sig):
    suf = ''
    if sig < 0:
        suf = '_m' + str(abs(sig))
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
        x = xExtr
        y = yExtr
        contour.SetPoint(0, x, y)

    if xEnd - xmin < xw or yEnd - ymin < yw:
        contour.GetPoint(contour.GetN() - 2, x, y)
        xNext = float(x)
        yNext = float(y)
        x = max(xEnd + (xEnd - xNext) / (yEnd - yNext) * (ymin - yEnd), xmin)
        y = max(yEnd + (yEnd - yNext) / (xEnd - xNext) * (xmin - xEnd), ymin)
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
        x = xExtr
        y = yExtr
        contour.SetPoint(0, x, y)

    if xmax - xEnd < xw or ymax - yEnd < yw:
        contour.GetPoint(contour.GetN() - 2, x, y)
        xNext = float(x)
        yNext = float(y)
        x = max(xEnd + (xEnd - xNext) / (yEnd - yNext) * (ymax - yEnd), xmax)
        y = max(yEnd + (yEnd - yNext) / (xEnd - xNext) * (xmax - xEnd), ymax)
        contour.Set(contour.GetN() + 1)
        contour.SetPoint(contour.GetN() - 1, x, y)


def drawLimits(model, pointFormat, titles):

    ndim = len(titles) - 1
    if ndim == 1:
        drawOption = 'LP'
    else:
        drawOption = 'COLZ'

    # Tree with information on one signal point per row. Columns are the name, cross section, median and +-1/2 sigma xsec upper bounds of the point.
    source = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/' + model + '.root', )
    input = source.Get('limitTree')

    vPointName = array.array('c', ' ' * 100)
    vXsec = array.array('d', [0.])
    vXsecErr = array.array('d', [0.])
    vLimit = array.array('d', [0.])
    vLimM2s = array.array('d', [0.])
    vLimM1s = array.array('d', [0.])
    vLimP1s = array.array('d', [0.])
    vLimP2s = array.array('d', [0.])
    vNEvents = array.array('i', [0])
                
    input.SetBranchAddress('pointName', vPointName)
    input.SetBranchAddress('xsec', vXsec)
    input.SetBranchAddress('xsecErr', vXsecErr)
    input.SetBranchAddress('limit', vLimit)
    input.SetBranchAddress('limM2s', vLimM2s)
    input.SetBranchAddress('limM1s', vLimM1s)
    input.SetBranchAddress('limP1s', vLimP1s)
    input.SetBranchAddress('limP2s', vLimP2s)
    input.SetBranchAddress('nEvents', vNEvents)

    vYields = {}
    branches = input.GetListOfBranches()
    for branch in branches:
        if 'Yield' in branch.GetName():
            y = array.array('i', [0])
            branch.SetAddress(y)
            vYields[branch.GetName()] = y

    xsecs = {}
    xsecErrors = {}
    limits = dict([(i, {}) for i in range(-2, 3)])
    nEvents = {}
    
    yields = {'Electron': {}, 'Muon': {}}

    centers = [[], []]

    iEntry = 0
    while input.GetEntry(iEntry) > 0:
        iEntry += 1

        pointStr = vPointName.tostring()
        pointName = pointStr[0:pointStr.find('\0')]
        matches = re.match(pointFormat, pointName)
        if not matches:
            raise RuntimeError('Invalid point name ' + pointName)

        coord = tuple(map(int, [matches.group(iD + 1) for iD in range(ndim)]))
        for iD in range(ndim):
            centers[iD].append(coord[iD])

        xsecs[coord] = vXsec[0]
        xsecErrors[coord] = vXsecErr[0]
        limits[-2][coord] = vLimM2s[0]
        limits[-1][coord] = vLimM1s[0]
        limits[0][coord] = vLimit[0]
        limits[1][coord] = vLimP1s[0]
        limits[2][coord] = vLimP2s[0]
        nEvents[coord] = vNEvents[0]

        yields['Electron'][coord] = 0
        yields['Muon'][coord] = 0
        for name, v in vYields.items():
            if name.startswith('el'):
                yields['Electron'][coord] += v[0]
            elif name.startswith('mu'):
                yields['Muon'][coord] += v[0]

    source.Close()

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
    else:
        template = ROOT.TH2D('template', ';' + ';'.join(titles[1:]), nbins[0], centers[0][0] - width[0] / 2., centers[0][-1] + width[0] / 2., nbins[1], centers[1][0] - width[1] / 2., centers[1][-1] + width[1] / 2.)

    canvas = ROOT.TCanvas('limits', 'limits', 800, 800)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    output = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/' + model + '.root', 'recreate')

    canvas.SetLogz(True)
    setPalette(1.)

    # acceptance plot
    for lept, table in yields.items():
        acceptance = template.Clone('acceptance_' + lept)
        acceptance.SetTitle('A #times #epsilon (' + lept + ' channel)')
        acceptance.GetYaxis().SetTitleOffset(1.35)
        for coord, y in table.items():
            try:
                bin = acceptance.FindBin(*coord)
                N = float(nEvents[coord])
                acceptance.SetBinContent(bin, y / N)
                acceptance.SetBinError(bin, y / N * math.sqrt(1. / y + 1. / N))
            except:
                pass

        acceptance.Write()

        acceptance.Draw(drawOption)
        canvas.Print(PLOTSDIR + '/' + model + '_acceptance_' + lept + '.pdf')

#    canvas.SetLogz(False)

    # cross section plots
    crosssect = {}
    for sig in range(-2, 3):
        crosssect[sig] = template.Clone('crosssect_' + str(sig) + 'sigma')
        title = 'Cross Section'
        if sig != 0:
            title += ' (' + str(sig) + ' #sigma)'
        crosssect[sig].SetTitle(title)
        crosssect[sig].GetYaxis().SetTitleOffset(1.35)
        crosssect[sig].GetZaxis().SetTitle('#sigma (pb)')

        minZ = 1000.
        for coord, xsec in xsecs.items():
            bin = crosssect[sig].FindBin(*coord)
            z = xsec * (1. + sig * xsecErrors[coord])
            crosssect[sig].SetBinContent(bin, z)
            if z < minZ: minZ = z

        crosssect[sig].SetMinimum(minZ * 0.5)

        crosssect[sig].Write()

        crosssect[sig].Draw(drawOption)
        canvas.Print(PLOTSDIR + '/' + model + '_crosssect' + suffix(sig) + '.pdf')

#    canvas.SetLogz(True)

    # upper limit plots
    upperlim = {}
    for sig in range(-2, 3):
        upperlim[sig] = template.Clone('upperlim[sig]_' + str(sig) + 'sigma')
        title = 'CMS Preliminary, 19.7 fb^{-1}, #sqrt{s} = 8 TeV'
        if sig != 0:
            title += ' (' + str(sig) + ' #sigma)'
        upperlim[sig].SetTitle(title)
        upperlim[sig].GetYaxis().SetTitleOffset(1.35)
        upperlim[sig].GetZaxis().SetTitle('95% CL upper limit on cross section (pb)')
        
        for coord, limit in limits[sig].items():
            bin = upperlim[sig].FindBin(*coord)
            upperlim[sig].SetBinContent(bin, limit)

        upperlim[sig].Write()

        upperlim[sig].Draw(drawOption)
        canvas.Print(PLOTSDIR + '/' + model + '_upperlim' + suffix(sig) + '.pdf')

    canvas.SetLogz(False)

    # expected exclusion plots
    expected = {}
    contours = {}
    for sig in range(-1, 2):
        expected[sig] = template.Clone('expected_' + str(sig) + 'sigma')
        title = 'Expected Exclusion'
        if sig != 0:
            title += ' (' + str(sig) + ' #sigma)'
        expected[sig].SetTitle(title)
        expected[sig].GetZaxis().SetTitle('#sigma - #sigma_{95%CL}^{up} (pb)')

        expected[sig].Add(crosssect[sig])
        expected[sig].Add(upperlim[sig], -1.)

        expected[sig].Write()

        contours[sig] = []

        if ndim == 2:
            clevel = array.array('d', [0.])
            contSource = expected[sig].Clone('contSource')
            contSource.SetContour(1, clevel)
            contSource.Draw('CONT LIST')
            canvas.Update()
            contList = ROOT.gROOT.GetListOfSpecials().FindObject('contours').At(0)
            for contour in contList:
                contours[sig].append(contour.Clone())

            contSource.Delete()

    canvas.SetLogz(True)

#    min = expected[0].GetMinimum()
#    max = expected[0].GetMaximum()
#    for iY in range(expected[0].GetNbinsY() + 1):
#        for iX in range(expected[0].GetNbinsX() + 1):
#            if crosssect[0].GetBinContent(iX, iY) == 0.:
#                expected[0].SetBinContent(iX, iY, -1.e10)
#
#    expected[0].SetMinimum(min * 1.1)
#    expected[0].SetMaximum(max * 1.1)
#
#    expected[0].Draw(drawOption)

    background = upperlim[0].Clone('background')
    setPalette(1.5)

    background.Draw(drawOption)
    maxY = 0.
    for contour in contours[0]:
        truncateContour(contour, crosssect[0])
        closeContour(contour, crosssect[0])
        contour.SetLineColor(ROOT.kRed)
        contour.SetLineWidth(4)
        contour.Draw('L')
        
        for iP in range(contour.GetN()):
            y = contour.GetY()[iP]
            if y > maxY: maxY = y

    for sig in [-1, 1]:
        for contour in contours[sig]:
            truncateContour(contour, crosssect[0])
            closeContour(contour, crosssect[0])
            contour.SetLineColor(ROOT.kRed)
            contour.SetLineWidth(2)
            contour.SetLineStyle(ROOT.kDashed)
            contour.Draw('L')

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
    header.SetX1NDC(0.15)
    header.SetX2NDC(0.85)
    header.SetY1NDC(0.7)
    header.SetY2NDC(0.9)
    header.SetTextFont(62)
    header.SetTextSize(0.035)
    header.SetFillStyle(1001)
    header.SetFillColor(ROOT.kWhite)
    header.SetBorderSize(2)
    header.AddText(0.02, 0.7, titles[0])
    line = header.AddLine(0.02, 0.52, 0.08, 0.52)
    line.SetLineWidth(2)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineColor(ROOT.kRed)
    line = header.AddLine(0.02, 0.48, 0.08, 0.48)
    line.SetLineWidth(4)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineColor(ROOT.kRed)
    line = header.AddLine(0.02, 0.44, 0.08, 0.44)
    line.SetLineWidth(2)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineColor(ROOT.kRed)
    header.AddText(0.1, 0.4, 'Expected #pm 1 #sigma_{experiment}')

    header.Draw()
    
    canvas.Print(PLOTSDIR + '/' + model + '_expected.pdf')
        

if __name__ == '__main__':

    import sys

    model = sys.argv[1]

    if 'T5wg' in model:
        form = 'T5wg_([0-9]+)_([0-9]+)'
        titles = ('pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{0/#pm} NLO+NLL exclusion', 'M_{#tilde{g}} (GeV)', 'M_{#tilde{#chi}} (GeV)')

    drawLimits(model, form, titles)

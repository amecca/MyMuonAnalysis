#!/usr/bin/env python3

from os import makedirs
import numpy as np
from optparse import OptionParser
import ROOT

parser = OptionParser()
# parser.add_option("-b", "--batch-mode", action="store_true", dest="doBatchMode", default=False)
parser.add_option("-o", "--outputdir", action="store", default="last_trackSeedSA")

(options, args) = parser.parse_args()

# if(options.doBatchMode):
ROOT.gROOT.SetBatch(True)

if(options.outputdir.startswith(('/', '.'))):
    outdir = options.outputdir
else:
    outdir = '/eos/home-a/amecca/www/MuonReco/'+options.outputdir

makedirs(outdir, exist_ok=True)

theFile = ROOT.TFile("muonAnalysis.root")
if(not theFile.IsOpen()):
    exit(1)

dirOld = theFile.Get("muonAnalysisOLD")
dirNew = theFile.Get("muonAnalysisNEW")

if (not(dirOld and dirNew)):
    print('Could not get OLD and NEW directories from file "{:s}"'.format(theFile.GetName()))
    exit(1)

def needsLogPlot(plot, threshold=0.7):
    if(plot.Class().InheritsFrom("TH1")):
       x = [ plot.GetBinContent(i) for i in range(1, plot.GetNbinsX()) ]
    elif(plot.Class().InheritsFrom("TGraph")):
       x = plot.GetY()

    if(not x):
        print(">>>Warning: {:s} <{:s}>: could not get x axis!".format(plot.GetName(), plot.Class().GetName()) )
        return False
       
    x = [e for e in x if abs(e) > 1e-10]

    if(len(x) == 0):
        print(">>>Warning {:s}: len == 0!".format(plot.GetName()))
        return False
    
    sorted_x = np.sort(x)
    n = len(x)
    cumx = np.cumsum(sorted_x, dtype=float)
    gini = (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return gini > threshold


def getTGraphMax(g):
    if(g.GetN() == 0): return 1
    if  (g.Class().InheritsFrom("TGraphAsymmErrors")):
        def getPMax(g, i): return g.GetPointY(i) + g.GetErrorYhigh(i)
    elif(g.Class().InheritsFrom("TGraphErrors")):
        def getPMax(g, i): return g.GetPointY(i) + g.GetErrorY(i)
    else:
        def getPMax(g, i): return g.GetPointY(i)
    
    return max([ getPMax(g, i) for i in range(g.GetN()) ])

def getTGraphMin(g):
    if(g.GetN() == 0): return 0
    if  (g.Class().InheritsFrom("TGraphAsymmErrors")):
        def getPMin(g, i): return g.GetPointY(i) - g.GetErrorYlow(i)
    elif(g.Class().InheritsFrom("TGraphErrors")):
        def getPMin(g, i): return g.GetPointY(i) - g.GetErrorY(i)
    else:
        def getPMin(g, i): return g.GetPointY(i)
    return min([ getPMin(g, i) for i in range(g.GetN()) ])


c0 = ROOT.TCanvas("c0", "Canvas", 1600, 900)
c0.cd()
ROOT.gStyle.SetOptStat(0)

for key in dirOld.GetListOfKeys():
    name_ = key.GetName()
    old = key.ReadObj().Clone()
    class_ = old.Class()

    if(not (class_.InheritsFrom("TH1") or class_.InheritsFrom("TGraph")) ):
        print('Skipping "{:s}" of type {:s}.'.format(name_, class_.GetName()))
        continue

    new = dirNew.Get(name_).Clone()
    if(not new):
        print('Could not retrieve "{:s}" from NEW folder.')
        continue

    old.SetLineColor(ROOT.kRed)
    new.SetLineColor(ROOT.kBlue)

    legend = ROOT.TLegend(0.78, 0.75, 0.9, 0.9)
    
    if(class_.InheritsFrom("TH1")):
        print('TH1', name_)
        old.SetFillColor(ROOT.kRed)
        new.SetFillColor(ROOT.kBlue)
        old.SetFillStyle(3003)
        new.SetFillStyle(3005)
        yMax = max(old.GetMaximum(), new.GetMaximum())
        yMin = max(old.GetMinimum(), new.GetMinimum())
        old.SetMaximum( yMax + (yMax-yMin)*0.05 )
        old.Draw("hist")
        new.Draw("hist same")
        legend.AddEntry(old, "old", "lpf")
        legend.AddEntry(new, "new", "lpf")
    elif(class_.InheritsFrom("TGraph")):
        old.SetMarkerColor(ROOT.kRed)
        new.SetMarkerColor(ROOT.kBlue)
        old.SetMarkerSize(2)
        new.SetMarkerSize(2)
        old.SetMarkerStyle(ROOT.kFullTriangleUp)
        new.SetMarkerStyle(ROOT.kOpenCircle)
        yMax = max(getTGraphMax(old), getTGraphMax(new))
        yMin = min(getTGraphMin(old), getTGraphMin(new))
        old.SetMaximum(1.05)
        old.GetYaxis().SetRangeUser(yMin - abs(yMax-yMin)*0.05, 1.05)#yMax + abs(yMax-yMin)*0.05)
        old.Draw("APE")
        new.Draw("PE")
        legend.AddEntry(old, "old", "lpe")
        legend.AddEntry(new, "new", "lpe")
    else:
        print('Problem with "{:s}"'.format(name_))
    
    if(needsLogPlot(old)):
        ROOT.gPad.SetLogy(True)
    else:
        ROOT.gPad.SetLogy(False)

    legend.Draw("same")
    ROOT.gPad.Update()
    for extension in ['png']:#, 'pdf', 'root']:
        ROOT.gPad.SaveAs("{:s}/{:s}.{:s}".format(outdir, name_, extension))
    ROOT.gPad.Clear()

del c0

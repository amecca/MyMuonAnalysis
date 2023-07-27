#!/usr/bin/env python3
import sys

# Use python3 instead of the default python2.7.5 (which breaks in importing ROOT, for some reason...)
pyver = sys.version_info.major
if pyver < 3:
    print(' __________________________         __________________________')
    print('|  ________________________ WARNING _________________________ ')
    print('| |')
    print('| | Use python3 instead of default python version ('+pyver+')')
    print('| | For some reason, python'+pyver+' breaks when importing ROOT')
    print('| | Type instead:')
    print('| |   python3 '+' '.join([str(a) for a in sys.argv]))
    print('| |__________________________________________________________')
    print('|_____________________________________________________________')
    print('')
    exit(1)

import ROOT as r
import os, math
from subprocess import check_output
from argparse import ArgumentParser

#r.gROOT.LoadMacro('/afs/cern.ch/user/t/trocino/Utilities/danystyle.C')
#r.setTDRStyle()
r.gROOT.SetBatch(True)
#r.tdrStyle.SetLegendTextSize(0.034)

#folds = ['OLD', 'NEW']
#plots = ['n', 'pt', 'eta', 'dxy', 'dz']
#ctgrs = ['sa', 'saupd']

parser = ArgumentParser(
    description = 'Make pretty plots from the output of MuonGeneralAnalyzer'
)
parser.add_argument('-i', '--infile', default='muonAnalysis.root')
parser.add_argument('-0', '--hits0' , action='store_true', help='Use tracks with 0 valid hits')

args = parser.parse_args()
print(args)

infile = r.TFile.Open(args.infile,'read')
olddir = infile.Get('muonAnalysisOLD')
newdir = infile.Get('muonAnalysisNEW')

username = os.environ['USER']
site_basedir = '/afs/cern.ch/user/{initial}/{user}/www/'.format(initial=username[0], user=username)
#outputdir = os.path.join(site_basedir, 'MuonReco', check_output(['date', '+%Y-%m-%d']).decode().rstrip('\n')+'_trackSeedSA')
outputdir = os.path.join(site_basedir, 'MuonReco', 'last')

if outputdir[-1:] != '/': outputdir += '/'
if(not os.path.isdir(outputdir)):
    os.mkdir(outputdir)
if not os.path.exists(outputdir+'index.php'):
    os.system('cp {}/index.php {}'.format(site_basedir, outputdir))

labstr = 'sa_'
plotnames = [key.GetName()[len(labstr):] for key in olddir.GetListOfKeys() if labstr in key.GetName()]# and '_den_' not in key.GetName()]
# for p in plotnames:
#     print(p)
# exit(0)

r.gStyle.SetOptStat("")#('neiuo')

c = r.TCanvas('cc', 'cc', 1600, 1600)
c.SetLeftMargin  (0.14)
# c.SetRightMargin (0.1)
# c.SetTopMargin   (0.1)
# c.SetBottomMargin(0.1)
c.Draw()
c.SetGridx()
c.SetGridy()


for plotname in plotnames:
    # if(not 'err' in plotname):#'eff_pt_err' in plotname):
    #     continue    
    hits0 = '0hits' if args.hits0 else ''
    oldsa    = olddir.Get('{:s}{:s}_{:s}'.format('sa'   , hits0, plotname))
    oldsaupd = olddir.Get('{:s}{:s}_{:s}'.format('saupd', hits0, plotname))
    oldgl    = olddir.Get('{:s}{:s}_{:s}'.format('gl'   , hits0, plotname))
    newsa    = newdir.Get('{:s}{:s}_{:s}'.format('sa'   , hits0, plotname))
    newsaupd = newdir.Get('{:s}{:s}_{:s}'.format('saupd', hits0, plotname))
    newgl    = newdir.Get('{:s}{:s}_{:s}'.format('gl'   , hits0, plotname))
    
    firstplot = oldsaupd #oldsa # just a reference
    
    c.SetLogy(False)
    c.SetLogx(False)
        
    isth1 = firstplot.InheritsFrom('TH1')
    if isth1:
        # print('>>> {:32s}'.format(plotname), ' oldgl:', oldgl.GetEntries(), '- newgl:', newgl.GetEntries())
        if oldsaupd.GetEntries()+newsaupd.GetEntries()==0: continue  #oldsa.GetEntries()+newsa.GetEntries()+
    else:
        # print('>>> {:32s}'.format(plotname), ' oldgl:', oldgl.GetN()      , '- newgl:', newgl.GetN())
        if oldsaupd.GetN()+newsaupd.GetN()==0: continue  #oldsa.GetN()+newsa.GetN()+
    
    # OLD, SA
    # oldsa.SetLineColor(r.kBlack)
    # OLD, SA updated
    oldsaupd.SetLineColor(r.kGreen+2)
    # OLD, GLOBAL
    oldgl.SetLineColor(r.kBlue)
    # NEW, SA
    # newsa.SetLineColor(r.kYellow+1)
    # NEW, SA updated
    newsaupd.SetLineColor(r.kOrange-3)#(r.kMagenta)
    # NEW, GLOBAL
    newgl.SetLineColor(r.kRed) #r.kCyan+1

    for plot in [oldsa, oldsaupd, oldgl, newsa, newsaupd, newgl]:
        plot.SetLineWidth(2)
    if isth1:
        oldsa.SetLineStyle(2)
        oldsaupd.SetLineStyle(9)
        oldgl.SetLineStyle(1)
        newsa.SetLineStyle(2)
        newsaupd.SetLineStyle(9)
        newgl.SetLineStyle(1)
    else:
        oldsa.SetMarkerStyle(r.kFullCircle)
        oldsaupd.SetMarkerStyle(r.kOpenCircle)
        oldgl.SetMarkerStyle(r.kFullTriangleUp)
        newsa.SetMarkerStyle(r.kFullSquare)
        newsaupd.SetMarkerStyle(r.kOpenSquare)
        newgl.SetMarkerStyle(r.kFullTriangleDown)
        for graph in [oldsa, oldsaupd, oldgl, newsa, newsaupd, newgl]:
            graph.SetMarkerSize(1.)
            graph.SetMarkerColor(graph.GetLineColor())
            for b in range(0, graph.GetN()):
                graph.SetPointEXhigh(b, 0.)
                graph.SetPointEXlow (b, 0.)
    
    ##
    c.cd()
    if('_pt' in plotname):
        c.SetLogx(True)
        firstplot.GetXaxis().SetNoExponent()
        firstplot.GetXaxis().SetMoreLogLabels()
        firstplot.GetXaxis().SetRangeUser(2., 300.)
        
    drawopt = '' if isth1 else 'ALPE'
    firstplot.Draw(drawopt)
    
    firstplot.GetXaxis().SetTitleSize(0.04)
    firstplot.GetYaxis().SetTitleSize(0.04)
    # ymax = max(oldsa.GetMaximum(), oldsaupd.GetMaximum(), oldgl.GetMaximum(), newsa.GetMaximum(), newsaupd.GetMaximum(), newgl.GetMaximum(), 1.)
    # ymin = max(min(oldsa.GetMinimum(), oldsaupd.GetMinimum(), oldgl.GetMinimum(), newsa.GetMinimum(), newsaupd.GetMinimum(), newgl.GetMinimum()), 0.)
    if(isth1):
        ymax = max(oldsaupd.GetMaximum(), oldgl.GetMaximum(), newsaupd.GetMaximum(), newgl.GetMaximum())
        firstplot.SetMaximum(ymax*1.15)  # (firstplot.GetMaximum()*1.15)
        if(not 'res' in plotname):
            firstplot.SetMinimum(0.)
        else:
            c.SetLogy()
            firstplot.SetMaximum(firstplot.GetMaximum()*4)
            firstplot.GetYaxis().SetMoreLogLabels()
            firstplot.GetYaxis().SetNoExponent()
    else:
        ymin, ymax = 0.9, 1. #oldsa.GetYaxis().GetXmin(), oldsa.GetYaxis().GetXmax()
        if(args.hits0):
            ymin = 10**-4
            c.SetLogy()
            firstplot.GetYaxis().SetMoreLogLabels()
            firstplot.GetYaxis().SetNoExponent()
        if(not '_pt_' in plotname):
            firstplot.SetMinimum(ymin)
        firstplot.SetMaximum(ymax + min(ymax*0.1, (ymax-ymin)*0.2))
    
    drawopt = 'SAME' if isth1 else 'LPE'
    # oldsa.Draw(drawopt)
    # oldsaupd.Draw(drawopt)
    oldgl.Draw(drawopt)
    # newsa.Draw(drawopt)
    newsaupd.Draw(drawopt)
    newgl.Draw(drawopt)
        
    # Legend
    left_legend = 'pt_err' in plotname
    x0_leg = 0.60 if not left_legend else 0.20
    leg = r.TLegend(x0_leg, 0.78, x0_leg+0.30, 0.90)
    #leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetLineWidth(0)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    
    legopt = 'L' if isth1 else 'LEP'
    # leg.AddEntry(oldsa   , 'standard SA'           , legopt)
    leg.AddEntry(oldsaupd, 'standard SA'  , legopt) # + vertex
    leg.AddEntry(oldgl   , 'standard GLOBAL'       , legopt)
    # leg.AddEntry(newsa   , 'track-seed SA'         , legopt)
    leg.AddEntry(newsaupd, 'track-seed SA', legopt) # + vertex
    leg.AddEntry(newgl   , 'track-seed GLOBAL'     , legopt)
    leg.Draw()
    
    # Save canvas
    c.SaveAs(outputdir+plotname+('_'+hits0 if args.hits0 else '')+'.png')

infile.Close()

#!/usr/bin/env python3
import sys

# Use python3 instead of the default python2.7.5 (which breaks in importing ROOT on LXPLUS, since it was built with python3 bindings only)
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
import os, subprocess
from subprocess import check_output
from argparse import ArgumentParser
import logging
import itertools


plot_style = {
    'oldSda'   : {'title':'standard SdA'       , 'color':r.kBlack   , 'line':r.kDashed, 'marker':r.kOpenTriangleUp  },
    'oldSdaUpd': {'title':'standard SdA Upd'   , 'color':r.kGreen+2 , 'line':9        , 'marker':r.kOpenTriangleDown},
    'oldGl'    : {'title':'standard Global'    , 'color':r.kBlue    , 'line':r.kSolid , 'marker':r.kOpenSquare      },
    'outer'    : {'title':'outer'              , 'color':r.kGray    , 'line':r.kDotted, 'marker':r.kOpenCircle      },
    'global'   : {'title':'global'             , 'color':r.kGray+1  , 'line':r.kDotted, 'marker':r.kOpenDiamond     },
    'newSda'   : {'title':'trak-seeded SdA'    , 'color':r.kYellow+1, 'line':r.kDashed, 'marker':r.kFullTriangleUp  },
    'newSdaUpd': {'title':'trak-seeded SdA Upd', 'color':r.kOrange-3, 'line':9        , 'marker':r.kFullTriangleDown},
    'newGl'    : {'title':'trak-seeded Global' , 'color':r.kRed     , 'line':r.kSolid , 'marker':r.kFullSquare      }
}


def customize_graph(graph, numerator):
    graph.SetLineWidth(2)

    graph.SetLineColor(plot_style[numerator]['color'])
    graph.SetMarkerStyle(plot_style[numerator]['marker'])
    graph.SetMarkerColor(graph.GetLineColor())
    graph.SetMarkerSize(1.)

    for b in range(0, graph.GetN()):
        graph.SetPointEXhigh(b, 0.)
        graph.SetPointEXlow (b, 0.)


def customize_hist(hist, numerator):
    graph.SetLineWidth(2)
    graph.SetLineColor(plot_style[numerator]['color'])


def main():
    #r.gROOT.LoadMacro('/afs/cern.ch/user/t/trocino/Utilities/danystyle.C')
    #r.setTDRStyle()
    r.gROOT.SetBatch(True)
    #r.tdrStyle.SetLegendTextSize(0.034)

    parser = ArgumentParser(
        description = 'Make pretty plots from the output of MuonGeneralAnalyzer'
    )
    parser.add_argument('-i', '--infile', default='muonAnalysis.root')
    parser.add_argument('-0', '--hits0' , action='store_true', help='Use tracks with 0 valid hits')
    parser.add_argument('--log', dest='loglevel', type=str.upper, metavar='LEVEL', default='WARNING')
    parser.add_argument('-o', '--output', default=None)

    args = parser.parse_args()

    if(not hasattr(logging, args.loglevel)):
        raise ValueError('Invalid log level: %s' % args.loglevel)

    logging.basicConfig(level=getattr(logging, args.loglevel))
    logging.info('args = %s', args)

    infile = r.TFile.Open(args.infile,'read')
    analysisdir = infile.Get('muonAnalysis')
    assert analysisdir, 'The file "{:s}" does not contain a folder named "muonAnalysis"'.format(infile.GetName())

    username = os.environ['USER']
    site_basedir = '/afs/cern.ch/user/{initial}/{user}/www/'.format(initial=username[0], user=username)
    if(args.output is None):
        #outputdir = os.path.join(site_basedir, 'MuonReco', check_output(['date', '+%Y-%m-%d']).decode().rstrip('\n')+'_trackSeedSA')
        outputdir = os.path.join(site_basedir, 'MuonReco', 'last')
    else:
        outputdir = args.output
    logging.info('Output will be saved in %s', outputdir)

    if(not os.path.isdir(outputdir)):
        os.mkdir(outputdir)
    if not os.path.exists(os.path.join(outputdir, 'index.php')):
        default_index = os.path.join(site_basedir, 'index.php')
        if(os.path.exists(default_index)):
            subprocess.check_call(['cp', default_index, outputdir+'/'])


    r.gStyle.SetOptStat("")#('neiuo')

    c = r.TCanvas('cc', 'cc', 1600, 1600)
    c.SetLeftMargin  (0.14)
    # c.SetRightMargin (0.1)
    # c.SetTopMargin   (0.1)
    # c.SetBottomMargin(0.1)
    c.Draw()
    c.SetGridx()
    c.SetGridy()

    ### Efficiency ###

    plotnames = [k.GetName() for k in analysisdir.GetListOfKeys()]
    num_names = {k for k in plotnames if k.startswith('num_')}
    den_names = {k for k in plotnames if k.startswith('den_')}
    res_names = {k for k in plotnames if k.startswith('res_')}

    logging.info('efficiency num plots = %d', len(num_names))
    logging.info('efficiency den plots = %d', len(den_names))
    logging.debug('efficiency num/den  = %.1f', len(num_names)/len(den_names))
    logging.info('resolution plots     = %d', len(res_names))
    logging.info('total plots          = %d', len(plotnames))

    # Check that each numerator histogram can be associated to its denominator
    for num_name in num_names:
        num_split = num_name.split('_')
        den_name = 'den_{den:s}_{variable:s}'.format(den=num_split[1], variable=num_split[3])
        if(den_name not in den_names):
            print('ERROR: "{:s}" not found in file'.format(den_name))
            continue
        # else: print(num_name, '-->', den_name)

    denominators     = {n.split('_')[1] for n in den_names}
    variables        = {n.split('_')[2] for n in den_names}
    denominators_num = {n.split('_')[1] for n in num_names}
    numerators       = {n.split('_')[2] for n in num_names}
    variables_num    = {n.split('_')[3] for n in num_names}

    assert denominators == denominators_num, 'The two denominators sets are different'
    assert variables    == variables_num   , 'The two variables sets are different'

    logging.info('denominators(%d) = %s', len(denominators), denominators)
    logging.info('variables   (%d) = %s', len(variables   ), variables   )
    logging.info('numerators  (%d) = %s', len(numerators  ), numerators  )

    for denominator, variable in itertools.product(denominators, variables):
        # Retrieve denominator
        den_hist = analysisdir.Get(f'den_{denominator}_{variable}')
        den_entries = den_hist.GetEntries()
        logging.debug('den_hist = %s', den_hist)
        if(den_entries == 0):
            logging.warning('%s has 0 entries. Skipping', den_hist.GetName())
            continue

        # Legend
        left_legend = variable == 'pt'
        x0_leg = 0.60 if not left_legend else 0.20
        leg = r.TLegend(x0_leg, 0.78, x0_leg+0.30, 0.90)
        #leg.SetNColumns(2)
        leg.SetBorderSize(0)
        leg.SetLineWidth(0)
        leg.SetLineColor(0)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)

        # Retrieve numerators and create ratio
        err_graphs = []
        for numerator in numerators:
            num_hist =  analysisdir.Get(f'num_{denominator}_{numerator}_{variable}')
            assert num_hist, 'Could not retrieve %s' % (num_hist.GetName())
            num_entries = num_hist.GetEntries()
            if(num_entries == 0):
                logging.warning('%s has 0 entries. Skipping', num_hist.GetName())
                continue
            logging.debug( 'Num={:.0f}  Den={:.0f}  ratio={:.2f}'.format(num_entries, den_entries, num_entries/den_entries) )
            err_graph = r.TGraphAsymmErrors(num_hist, den_hist)
            customize_graph(err_graph, numerator=numerator)
            err_graphs.append(err_graph)
            leg.AddEntry(err_graph, 'standard GLOBAL', 'LEP')

        # Sanity check: if all the numerators were missing, skip this (denominator, variable) combination
        if(len(err_graphs) == 0):
            continue

        firstplot = err_graphs[0]
        c.cd()
        if(variable == 'pt'):
            c.SetLogx(True)
            firstplot.GetXaxis().SetNoExponent()
            firstplot.GetXaxis().SetMoreLogLabels()
            firstplot.GetXaxis().SetRangeUser(2., 300.)
        else:
            c.SetLogx(False)

        firstplot.Draw('ALPE')
        firstplot.GetXaxis().SetTitleSize(0.04)
        firstplot.GetYaxis().SetTitleSize(0.04)

        ymin, ymax = 0.9, 1. #oldsa.GetYaxis().GetXmin(), oldsa.GetYaxis().GetXmax()
        firstplot.SetMaximum(ymax + min(ymax*0.1, (ymax-ymin)*0.2))
        for g in err_graphs[1:]:
            g.Draw('LPE')

        leg.Draw()

        # Save canvas
        c.SaveAs(os.path.join(outputdir, f'eff_{denominator}_{numerator}_{variable}.png'))

    ### Resolution ###
    # Group plots based on reference object (the "denominator" if these were efficiencies) and variable
    references_res = {k.split('_')[1] for k in res_names}
    variables_res  = {k.split('_')[3] for k in res_names}

    for denominator, variable in itertools.product(references_res, variables_res):
        pass

    infile.Close()


def old_procedure():
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

if __name__ == '__main__':
    main()

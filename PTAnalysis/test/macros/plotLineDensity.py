#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT

import CMS_lumi, tdrstyle


#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "      Phase 2 Simulation"
#CMS_lumi.cmsTextSize = 0.35
CMS_lumi.lumi_sqrtS = "14 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPeriod = 0
iPos = 0

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
#ROOT.gStyle.SetLabelSize(0.04,'X')
#ROOT.gStyle.SetLabelSize(0.04,'Y')
#ROOT.gStyle.SetTitleSize(0.04,'X')
#ROOT.gStyle.SetTitleSize(0.04,'Y')
#ROOT.gStyle.SetTitleOffset(1.0,'X')
#ROOT.gStyle.SetTitleOffset(1.4,'Y')
#ROOT.gStyle.SetTextFont(42)
#ROOT.gStyle.SetLegendFont(42)


resol = [30, 40, 50, 60, 70]


f = {}

h_linedensity_noRelChIsoZCut = {}
h_linedensity_noRelChIsoZTCut = {}
h_linedensity_RelChIsoZCut = {}
h_linedensity_RelChIsoZTCut = {}

h_eff_RelChIsoZCut = {}
h_eff_RelChIsoZTCut = {}


nRe = 20

for i,res in enumerate(resol):
    f[res] = {}
    h_linedensity_noRelChIsoZCut[res] = {}
    h_linedensity_noRelChIsoZTCut[res] = {}
    h_linedensity_RelChIsoZCut[res] = {}
    h_linedensity_RelChIsoZTCut[res] = {}
    h_eff_RelChIsoZCut[res] = {}
    h_eff_RelChIsoZTCut[res] = {}
    for proc in 'sig', 'bkg':
        fname = '../output_muIso_DYToLL_%dps_prompt.root'%res
        if (proc == 'bkg'): fname = '../output_muIso_TTbar_%dps_fake.root'%res
        f[res][proc] = ROOT.TFile.Open(fname)

        #h_linedensity_noRelChIsoZCut[res][proc]  =  f[res][proc].Get('h_linedensity_noRelChIsoZCut')
        #h_linedensity_noRelChIsoZTCut[res][proc]  =  f[res][proc].Get('h_linedensity_noRelChIsoZTCut')
        #h_linedensity_RelChIsoZCut[res][proc]  =  f[res][proc].Get('h_linedensity_RelChIsoZCut_simVtx')
        #h_linedensity_RelChIsoZTCut[res][proc]  =  f[res][proc].Get('h_linedensity_RelChIsoZTCut_simVtx')

        h_linedensity_noRelChIsoZCut[res][proc]  =  f[res][proc].Get('h_linedensity_noRelChIsoZCut_endcap')
        h_linedensity_noRelChIsoZTCut[res][proc]  =  f[res][proc].Get('h_linedensity_noRelChIsoZTCut_endcap')
        h_linedensity_RelChIsoZCut[res][proc]  =  f[res][proc].Get('h_linedensity_RelChIsoZCut_simVtx_endcap')
        h_linedensity_RelChIsoZTCut[res][proc]  =  f[res][proc].Get('h_linedensity_RelChIsoZTCut_simVtx_endcap')

        h_linedensity_noRelChIsoZCut[res][proc].Rebin(nRe)
        h_linedensity_noRelChIsoZTCut[res][proc].Rebin(nRe)
        h_linedensity_RelChIsoZCut[res][proc].Rebin(nRe)
        h_linedensity_RelChIsoZTCut[res][proc].Rebin(nRe)
        
        h_linedensity_noRelChIsoZCut[res][proc].Sumw2()
        h_linedensity_noRelChIsoZTCut[res][proc].Sumw2()
        h_linedensity_RelChIsoZCut[res][proc].Sumw2()
        h_linedensity_RelChIsoZTCut[res][proc].Sumw2()
        
        h_eff_RelChIsoZCut[res][proc] =  h_linedensity_RelChIsoZCut[res][proc].Clone('h_eff_RelChIsoZCut_%s_%dps'%(proc,res))
        h_eff_RelChIsoZTCut[res][proc] =  h_linedensity_RelChIsoZTCut[res][proc].Clone('h_eff_RelChIsoZTCut_%s_%dps'%(proc,res))

        h_eff_RelChIsoZCut[res][proc].Divide(h_eff_RelChIsoZCut[res][proc],h_linedensity_noRelChIsoZCut[res][proc])
        h_eff_RelChIsoZTCut[res][proc].Divide(h_eff_RelChIsoZTCut[res][proc],h_linedensity_noRelChIsoZTCut[res][proc])

        h_eff_RelChIsoZCut[res][proc].SetLineColor(ROOT.kBlue)
        h_eff_RelChIsoZTCut[res][proc].SetLineColor(ROOT.kRed+i)

        h_eff_RelChIsoZCut[res][proc].SetMarkerColor(ROOT.kBlue)
        h_eff_RelChIsoZTCut[res][proc].SetMarkerColor(ROOT.kRed+i)

        h_eff_RelChIsoZCut[res][proc].SetMarkerStyle(20)
        h_eff_RelChIsoZTCut[res][proc].SetMarkerStyle(20)
        
        h_eff_RelChIsoZCut[res][proc].SetLineWidth(2)
        h_eff_RelChIsoZTCut[res][proc].SetLineWidth(2)

        
canvas = ROOT.TCanvas('eff_vs_linedensity','eff_vs_linedensity')
canvas.Divide(1,2)
canvas.cd(1)
canvas.cd(1).SetGridx()
canvas.cd(1).SetGridy()
h_eff_RelChIsoZCut[30]['sig'].GetYaxis().SetRangeUser(0.8,1.0)
h_eff_RelChIsoZCut[30]['sig'].GetYaxis().SetTitle('prompt efficiency')
h_eff_RelChIsoZCut[30]['sig'].GetXaxis().SetTitle('line density (mm^{-1})')
h_eff_RelChIsoZCut[30]['sig'].Draw('e')
h_eff_RelChIsoZTCut[30]['sig'].Draw('esame')
h_eff_RelChIsoZTCut[40]['sig'].Draw('esame')
h_eff_RelChIsoZTCut[60]['sig'].Draw('esame')

canvas.cd(2)
canvas.cd(2).SetGridx()
canvas.cd(2).SetGridy()
h_eff_RelChIsoZCut[30]['bkg'].GetYaxis().SetRangeUser(0.0,0.1)
h_eff_RelChIsoZCut[30]['bkg'].GetYaxis().SetTitle('non-prompt efficiency')
h_eff_RelChIsoZCut[30]['bkg'].GetXaxis().SetTitle('line density (mm^{-1})')
h_eff_RelChIsoZCut[30]['bkg'].Draw('e')
h_eff_RelChIsoZTCut[30]['bkg'].Draw('esame')
h_eff_RelChIsoZTCut[40]['bkg'].Draw('esame')
h_eff_RelChIsoZTCut[60]['bkg'].Draw('esame')

CMS_lumi.CMS_lumi(canvas.cd(1), iPeriod, iPos)

#canvas.SaveAs(canvas.GetName()+'_simVtx.pdf')
canvas.SaveAs(canvas.GetName()+'.pdf')


raw_input('ok?')


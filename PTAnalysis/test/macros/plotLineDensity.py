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




fSig = '../output_muIso_DYToLL_30ps_prompt.root '
fBkg = '../output_muIso_TTbar_30ps_fake.root'
#fBkg = '../output_muIso_QCD_30ps_fake.root'

f = {}
f['sig'] = ROOT.TFile.Open(fSig)
f['bkg'] = ROOT.TFile.Open(fBkg)

h_linedensity_noRelChIsoZCut = {}
h_linedensity_noRelChIsoZTCut = {}
h_linedensity_RelChIsoZCut = {}
h_linedensity_RelChIsoZTCut = {}

h_eff_RelChIsoZCut = {}
h_eff_RelChIsoZTCut = {}


nRe = 20

for proc in 'sig', 'bkg':
    h_linedensity_noRelChIsoZCut[proc]  =  f[proc].Get('h_linedensity_noRelChIsoZCut')
    h_linedensity_noRelChIsoZTCut[proc]  =  f[proc].Get('h_linedensity_noRelChIsoZTCut')
    h_linedensity_RelChIsoZCut[proc]  =  f[proc].Get('h_linedensity_RelChIsoZCut')
    h_linedensity_RelChIsoZTCut[proc]  =  f[proc].Get('h_linedensity_RelChIsoZTCut')
    #h_linedensity_RelChIsoZCut[proc]  =  f[proc].Get('h_linedensity_RelChIsoZCut_simVtx')
    #h_linedensity_RelChIsoZTCut[proc]  =  f[proc].Get('h_linedensity_RelChIsoZTCut_simVtx')

    h_linedensity_noRelChIsoZCut[proc].Rebin(nRe)
    h_linedensity_noRelChIsoZTCut[proc].Rebin(nRe)
    h_linedensity_RelChIsoZCut[proc].Rebin(nRe)
    h_linedensity_RelChIsoZTCut[proc].Rebin(nRe)

    h_linedensity_noRelChIsoZCut[proc].Sumw2()
    h_linedensity_noRelChIsoZTCut[proc].Sumw2()
    h_linedensity_RelChIsoZCut[proc].Sumw2()
    h_linedensity_RelChIsoZTCut[proc].Sumw2()
    
    h_eff_RelChIsoZCut[proc] =  h_linedensity_RelChIsoZCut[proc].Clone('h_eff_RelChIsoZCut_%s'%proc)
    h_eff_RelChIsoZTCut[proc] =  h_linedensity_RelChIsoZTCut[proc].Clone('h_eff_RelChIsoZTCut_%s'%proc)

    h_eff_RelChIsoZCut[proc].Divide(h_eff_RelChIsoZCut[proc],h_linedensity_noRelChIsoZCut[proc])
    h_eff_RelChIsoZTCut[proc].Divide(h_eff_RelChIsoZTCut[proc],h_linedensity_noRelChIsoZTCut[proc])

    h_eff_RelChIsoZCut[proc].SetLineColor(ROOT.kBlue)
    h_eff_RelChIsoZTCut[proc].SetLineColor(ROOT.kRed)

    h_eff_RelChIsoZCut[proc].SetMarkerColor(ROOT.kBlue)
    h_eff_RelChIsoZTCut[proc].SetMarkerColor(ROOT.kRed)

    h_eff_RelChIsoZCut[proc].SetMarkerStyle(20)
    h_eff_RelChIsoZTCut[proc].SetMarkerStyle(20)

    h_eff_RelChIsoZCut[proc].SetLineWidth(2)
    h_eff_RelChIsoZTCut[proc].SetLineWidth(2)

canvas = ROOT.TCanvas('eff_vs_linedensity','eff_vs_linedensity')
canvas.Divide(1,2)
canvas.cd(1)
canvas.cd(1).SetGridx()
canvas.cd(1).SetGridy()
h_eff_RelChIsoZCut['sig'].GetYaxis().SetRangeUser(0.8,1.0)
h_eff_RelChIsoZCut['sig'].GetYaxis().SetTitle('prompt efficiency')
h_eff_RelChIsoZCut['sig'].GetXaxis().SetTitle('line density (mm^{-1})')
h_eff_RelChIsoZCut['sig'].Draw('e')
h_eff_RelChIsoZTCut['sig'].Draw('esame')

canvas.cd(2)
canvas.cd(2).SetGridx()
canvas.cd(2).SetGridy()
h_eff_RelChIsoZCut['bkg'].GetYaxis().SetRangeUser(0.0,0.1)
h_eff_RelChIsoZCut['bkg'].GetYaxis().SetTitle('non-prompt efficiency')
h_eff_RelChIsoZCut['bkg'].GetXaxis().SetTitle('line density (mm^{-1})')
h_eff_RelChIsoZCut['bkg'].Draw('e')
h_eff_RelChIsoZTCut['bkg'].Draw('esame')

CMS_lumi.CMS_lumi(canvas.cd(1), iPeriod, iPos)

#canvas.SaveAs(canvas.GetName()+'_simVtx.pdf')
canvas.SaveAs(canvas.GetName()+'.pdf')


raw_input('ok?')


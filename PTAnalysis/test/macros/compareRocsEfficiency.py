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
CMS_lumi.cmsTextSize = 0.35
CMS_lumi.lumi_sqrtS = "14 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumiTextSize = 0.35

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


whichVtx = sys.argv[1]

resolutions = ['30']
efficiencies = [ 90, 100 ]

tl = ROOT.TLatex( 0.70, 0.84,'<PU> = 200')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, t#bar{t}')
#tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, QCD')
tl2.SetNDC()
tl2.SetTextSize(0.035)

f = {}
roc = {}
roc_dT =  {}

c = {}

regions = ['','barrel','endcap']


col = {'30':ROOT.kRed,
       '40':ROOT.kOrange,
       '50':ROOT.kGreen,
       '60':ROOT.kCyan,
       '70':ROOT.kBlue,
}

linestyle = { 100 : 1, 90 : 2}



leg = {}

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.13, '|#eta|<2.5')
tt['barrel'] = ROOT.TLatex( 0.15, 0.13, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.13, '1.5 < |#eta| < 2.5 ')
for reg in regions:
    tt[reg].SetNDC()
    tt[reg].SetTextSize(0.035)

    
gdummy = ROOT.TGraph()
gdummy_dT = ROOT.TGraph()

for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.65, 0.45, 0.89)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('muonIso_roc_comparison_eff_%s'%reg,'muonIso_roc_comparison_eff_%s'%reg,700,700)
    if (reg == ''):
        c[reg] = ROOT.TCanvas('muonIso_roc_comparison_eff','muonIso_roc_comparison_eff',700,700)
    c[reg].SetGridx()
    c[reg].SetGridy()
    
    for ires,res in enumerate(resolutions):
        f[res] = {}
        roc[reg][res] = {}
        roc_dT[reg][res] = {}
        
        for ieff, eff in enumerate(efficiencies):
            if (eff == 100):
                f[res][eff] = ROOT.TFile.Open('%sps_PU200/roc_%sps_PU200.root'%(res,res))
            else:
                f[res][eff] = ROOT.TFile.Open('%sps_PU200_eff%d/roc_%sps_PU200.root'%(res,eff,res))

            
                if (whichVtx == 'simVtx'):
                    gname  = 'g_roc_relChIso03_dZ1_simVtx_%s'%reg
                    gname_dT = 'g_roc_relChIso03_dZ1_simVtx_dT3s_%s'%reg
                else:
                    gname  = 'g_roc_relChIso03_dZ1_%s'%reg
                    gname_dT = 'g_roc_relChIso03_dZ1_dT3s_%s'%reg
                    #gname  = 'g_roc_relChIso03_dZ2_%s'%reg
                    #gname_dT = 'g_roc_relChIso03_dZ2_dT_%s'%reg

            
                if (reg == ''):
                    gname = gname[:-1]
                    gname_dT = gname_dT[:-1]


            roc[reg][res][eff] = f[res][eff].Get(gname) 
            roc_dT[reg][res][eff] = f[res][eff].Get(gname_dT)
            

            roc[reg][res][eff].SetLineColor(1)
            roc[reg][res][eff].SetLineStyle(1)
            roc[reg][res][eff].SetLineWidth(2)
            roc[reg][res][eff].SetFillStyle(1001)
            roc[reg][res][eff].SetFillColorAlpha(ROOT.kBlack,0.25)
        
            roc_dT[reg][res][eff].SetLineColor(col[res])
            roc_dT[reg][res][eff].SetLineStyle(linestyle[eff])
            roc_dT[reg][res][eff].SetFillStyle(0)
            roc_dT[reg][res][eff].SetFillColorAlpha(col[res],0.0)
            
            if (ires == 0 and ieff == 0):
                roc[reg][res][eff].GetXaxis().SetTitle('prompt efficiency')
                roc[reg][res][eff].GetYaxis().SetTitle('non-prompt efficiency')
                roc[reg][res][eff].GetYaxis().SetTitleOffset(1.4)
                roc[reg][res][eff].Draw('A C E3')
                leg[reg].AddEntry(roc[reg][res][eff],'no MTD','PL')
            
            roc_dT[reg][res][eff].Draw('C E3 same')
            leg[reg].AddEntry(roc_dT[reg][res][eff],'MTD, #sigma_{t} = %s ps, #epsilon = %d %%'%(res,eff),'PL')

    tt[reg].Draw()    
    leg[reg].Draw()
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    raw_input('ok?')


#save plots
dirname = 'RocsComparison/'
if (whichVtx == 'simVtx'):
    dirname = 'RocsComparison_simVtx/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for reg in regions:
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')

raw_input('ok?')

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

CMS_lumi.writeExtraText = True
iPeriod = 0
iPos = 0

#set the tdr style
tdrstyle.setTDRStyle()



pu = sys.argv[1]
inputDir = '93X/30ps_PU200_TTbar_minTkPtCut/'
if (pu == 'noPU'):
    inputDir = '93X/30ps_noPU_TTbar_minTkPtCut/'
    
dts = ['2s','3s','5s']
dtvals = {'2s': 2, '3s': 3, '5s': 5 }  

f = {}
roc = {}
roc_dT =  {}
c = {}

regions = ['barrel','endcap']


col = { '2s':ROOT.kRed,
        '3s':ROOT.kOrange,
        '5s':ROOT.kGreen}


tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
if (pu == 'noPU'):
    tl = ROOT.TLatex( 0.65, 0.88,'noPU')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
tl2.SetNDC()
tl2.SetTextSize(0.035)

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.16, '|#eta|<2.8')
tt['barrel'] = ROOT.TLatex( 0.15, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.16, '1.5 < |#eta| < 2.8')
for reg in regions:
    tt[reg].SetNDC()
    tt[reg].SetTextSize(0.035)
    
leg = {}

for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('roc_comparison_dt_%s_%s'%(reg,pu),'roc_comparison_dt_%s_%s'%(reg,pu))
    if (reg == ''):
        c[reg] = ROOT.TCanvas('roc_comparison_dt_%s'%pu,'roc_comparison_dt_%s'%pu)
    c[reg].SetGridx()
    c[reg].SetGridy()

    for idt,dt in enumerate(dts):
        if (pu == 'PU200'):
            f[dt] = ROOT.TFile.Open(inputDir+'/roc_30ps_PU200.root')
        if (pu == 'noPU'):
            f[dt] = ROOT.TFile.Open(inputDir+'/roc_30ps_noPU.root')

            
        if (reg == 'barrel'):
            gname = 'g_roc_relChIso03_dZ1_simVtx_%s'%(reg)
            gname_dT = 'g_roc_relChIso03_dZ1_simVtx_dT%s_%s'%(dt,reg)
        else:
            gname = 'g_roc_relChIso03_dZ2_simVtx_%s'%(reg)
            gname_dT = 'g_roc_relChIso03_dZ2_simVtx_dT%s_%s'%(dt,reg)
            
            
        roc[reg][dt] = f[dt].Get(gname) 
        roc_dT[reg][dt] = f[dt].Get(gname_dT)
                    
        roc[reg][dt].SetMarkerColor(col[dt])
        roc[reg][dt].SetLineColor(col[dt])
        roc[reg][dt].SetLineWidth(2)
        roc[reg][dt].SetFillStyle(0)
        roc[reg][dt].SetFillColorAlpha(col[dt],0.0)
    
        roc_dT[reg][dt].SetMarkerColor(col[dt])
        roc_dT[reg][dt].SetLineColor(col[dt])
        roc_dT[reg][dt].SetLineWidth(2)
        roc_dT[reg][dt].SetLineStyle(2)
        roc_dT[reg][dt].SetFillStyle(0)
        roc_dT[reg][dt].SetFillColorAlpha(col[dt],0.0)
    
        if (idt == 0):
            roc[reg][dt].GetXaxis().SetTitle('Prompt efficiency')
            roc[reg][dt].GetYaxis().SetTitle('Non-prompt efficiency')
            roc[reg][dt].GetXaxis().SetRangeUser(0.85,1.0)
            roc[reg][dt].GetYaxis().SetRangeUser(0.0,0.1)
            roc[reg][dt].Draw('A C E3')
            leg[reg].AddEntry(roc[reg][dt],'no MTD','L')
        #else:
        #    roc[reg][dt].Draw('C E3 SAME')

        roc_dT[reg][dt].Draw('C E3 same')
            

        #leg[reg].AddEntry(roc[reg][dt],'no MTD','PL')
        leg[reg].AddEntry(roc_dT[reg][dt],'with MTD - |#Deltat| < %d#sigma_{t}'%dtvals[dt],'L')

    tt[reg].Draw()    
    tl.Draw()
    tl2.Draw()
    leg[reg].Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)

#save plots
dirname = 'RocsComparisonDt/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for reg in regions:
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')

raw_input('ok?')

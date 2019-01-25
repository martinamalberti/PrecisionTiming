#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT

ROOT.gStyle.SetOptTitle(0)

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetStatX(0.6);
#ROOT.gStyle.SetStatY(0.6);
#ROOT.gStyle.SetStatW(0.3);
#ROOT.gStyle.SetStatH(0.3);

ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetFitFormat("4.4g")

ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.0,'Y')
#ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


pu = sys.argv[1]

dts = ['2s','3s','5s']
dtvals = {'2s': 2, '3s': 3, '5s': 5 }  

f = {}
roc = {}
roc_dT =  {}
c = {}

regions = ['','barrel','endcap']


col = { '2s':ROOT.kRed,
        '3s':ROOT.kOrange,
        '5s':ROOT.kGreen}

leg = {}

for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('roc_comparison_dt_%s_%s'%(reg,pu),'roc_comparison_dt_%s_%s'%(reg,pu),700,700)
    if (reg == ''):
        c[reg] = ROOT.TCanvas('roc_comparison_dt_%s'%pu,'roc_comparison_dt_%s'%pu,700,700)
    c[reg].SetGridx()
    c[reg].SetGridy()

    for idt,dt in enumerate(dts):
        f[dt] = ROOT.TFile.Open('30ps_PU200/roc_30ps_PU200.root')
        if (pu == 'noPU'):
            f[dt] = ROOT.TFile.Open('30ps_noPU/roc_30ps_noPU.root')

        #gname = 'g_roc_relChIso03_dZ1_%s'%(reg)
        #gname_dT = 'g_roc_relChIso03_dZ1_dT%s_%s'%(dt,reg)

        gname = 'g_roc_relChIso03_dZ1_simVtx_%s'%(reg)
        gname_dT = 'g_roc_relChIso03_dZ1_simVtx_dT%s_%s'%(dt,reg)

        if (reg == ''):
            gname = gname[:-1]
            gname_dT = gname_dT[:-1]

            
        roc[reg][dt] = f[dt].Get(gname) 
        roc_dT[reg][dt] = f[dt].Get(gname_dT)

                    
        roc[reg][dt].SetLineColor(col[dt])
        roc[reg][dt].SetLineWidth(2)
        roc[reg][dt].SetFillStyle(0)
        roc[reg][dt].SetFillColorAlpha(col[dt],0.0)
    
        roc_dT[reg][dt].SetLineColor(col[dt])
        roc_dT[reg][dt].SetLineWidth(2)
        roc_dT[reg][dt].SetLineStyle(2)
        roc_dT[reg][dt].SetFillStyle(0)
        roc_dT[reg][dt].SetFillColorAlpha(col[dt],0.0)
    
        if (idt == 0):
            roc[reg][dt].Draw('A C E3')
        else:
            roc[reg][dt].Draw('C E3 SAME')

        roc_dT[reg][dt].Draw('C E3 same')
            

        leg[reg].AddEntry(roc[reg][dt],'no MTD','PL')
        leg[reg].AddEntry(roc_dT[reg][dt],'with MTD - |#Delta t| < %d#sigma_{t}'%dtvals[dt],'PL')


    leg[reg].Draw()


#save plots
dirname = 'RocsComparisonDt/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for reg in regions:
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')

raw_input('ok?')

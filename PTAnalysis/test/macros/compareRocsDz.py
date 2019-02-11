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

inputDir = '93X/30ps_PU200_TTbar_minTkPtCut/'
if (pu == 'noPU'):
    inputDir = '93X/30ps_noPU_TTbar_minTkPtCut/'

print inputDir

dzs = ['05','1','2','3', '10']
dzvals = {'1': 1, '05': 0.5, '2':2, '3':3, '10':10 }

f = {}
roc = {}
roc_dT =  {}
c = {}

regions = ['','barrel','endcap']


col = { '05':ROOT.kOrange,
        '1' :ROOT.kRed,
        '2' :ROOT.kGreen,
        '3' :ROOT.kCyan,
        '5' :ROOT.kMagenta,
        '10' :ROOT.kBlue}

leg = {}

for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('roc_comparison_dz_%s_%s'%(reg,pu),'roc_comparison_dz_%s_%s'%(reg,pu),700,700)
    if (reg == ''):
        c[reg] = ROOT.TCanvas('roc_comparison_dz_%s'%pu,'roc_comparison_dz_%s'%pu,700,700)
    c[reg].SetGridx()
    c[reg].SetGridy()

    for idz,dz in enumerate(dzs):
        if (pu == 'PU200'):
            f[dz] = ROOT.TFile.Open(inputDir+'/roc_30ps_PU200.root')
        if (pu == 'noPU'):
            f[dz] = ROOT.TFile.Open(inputDir+'/roc_30ps_noPU.root')

        #gname = 'g_roc_relChIso03_dZ%s_%s'%(dz,reg)
        #gname_dT = 'g_roc_relChIso03_dZ%s_dT3s_%s'%(dz,reg)

        gname = 'g_roc_relChIso03_dZ%s_simVtx_%s'%(dz,reg)
        gname_dT = 'g_roc_relChIso03_dZ%s_simVtx_dT3s_%s'%(dz,reg)
        
        if (reg == ''):
            gname = gname[:-1]
            gname_dT = gname_dT[:-1]

                
        roc[reg][dz] = f[dz].Get(gname) 
        roc_dT[reg][dz] = f[dz].Get(gname_dT)

                    
        roc[reg][dz].SetLineColor(col[dz])
        roc[reg][dz].SetLineWidth(2)
        roc[reg][dz].SetFillStyle(0)
        roc[reg][dz].SetFillColorAlpha(col[dz],0.0)
    
        roc_dT[reg][dz].SetLineColor(col[dz])
        roc_dT[reg][dz].SetLineWidth(2)
        roc_dT[reg][dz].SetLineStyle(2)
        roc_dT[reg][dz].SetFillStyle(0)
        roc_dT[reg][dz].SetFillColorAlpha(col[dz],0.0)
    
        if (idz == 0):
            roc[reg][dz].Draw('A C E3')
        else:
            roc[reg][dz].Draw('C E3 SAME')

        roc_dT[reg][dz].Draw('C E3 same')
            

        leg[reg].AddEntry(roc[reg][dz],'no MTD - |d_{z}| < %.01f mm'%dzvals[dz],'PL')
        leg[reg].AddEntry(roc_dT[reg][dz],'with MTD - |d_{z}| < %.01f mm'%dzvals[dz],'PL')


    leg[reg].Draw()


#save plots
dirname = 'RocsComparisonDz/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for reg in regions:
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')

raw_input('ok?')

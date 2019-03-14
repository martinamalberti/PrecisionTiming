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
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)

pu = sys.argv[1]
release = '93X'

#inputDir = release+'/30ps_PU200_TTbar_minTkPtCut_noDxy/'
#if (pu == 'noPU'):
#    inputDir = release+'/30ps_noPU_TTbar_minTkPtCut_noDxy/'

inputDir = release+'/30ps_PU200_QCD_minTkPtCut_noDxy/'
if (pu == 'noPU'):
    inputDir = release+'/30ps_noPU_QCD_minTkPtCut_noDxy/'

#inputDir = '10_4_0_mtd5/35ps_PU200_TTbar/'
#if (pu == 'noPU'):
#    inputDir = '10_4_0_mtd5/35ps_noPU_TTbar/'

print inputDir

#dzs = ['05','1','2','3', '10', '50']
dzs = ['05','1','2','3','10']
dzvals = {'1': 1, '05': 0.5, '2':2, '3':3, '10':10, '50': 50}


minEff = 0.80

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
        '10' :ROOT.kBlue,
        '50' :ROOT.kBlue+2}

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
if (pu == 'noPU'):
    tl = ROOT.TLatex( 0.65, 0.88,'noPU')
tl.SetNDC()
tl.SetTextSize(0.035)

#tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, QCD')
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
    leg[reg] = ROOT.TLegend(0.15, 0.65, 0.50, 0.92)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('roc_comparison_dz_%s_%s'%(reg,pu),'roc_comparison_dz_%s_%s'%(reg,pu))
    if (reg == ''):
        c[reg] = ROOT.TCanvas('roc_comparison_dz_%s'%pu,'roc_comparison_dz_%s'%pu)
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


        roc[reg][dz].SetMarkerColor(col[dz])
        roc[reg][dz].SetLineColor(col[dz])
        roc[reg][dz].SetLineWidth(2)
        roc[reg][dz].SetFillStyle(0)
        roc[reg][dz].SetFillColorAlpha(col[dz],0.0)

        roc_dT[reg][dz].SetMarkerColor(col[dz])
        roc_dT[reg][dz].SetLineColor(col[dz])
        roc_dT[reg][dz].SetLineWidth(2)
        roc_dT[reg][dz].SetLineStyle(2)
        roc_dT[reg][dz].SetFillStyle(0)
        roc_dT[reg][dz].SetFillColorAlpha(col[dz],0.0)
    
        if (idz == 0):
            roc[reg][dz].GetXaxis().SetTitle('Prompt efficiency')
            roc[reg][dz].GetYaxis().SetTitle('Non-prompt efficiency')
            roc[reg][dz].GetYaxis().SetTitleOffset(1.1)
            roc[reg][dz].GetXaxis().SetRangeUser(minEff,1.0)
            roc[reg][dz].GetYaxis().SetRangeUser(0.0,0.1)
            roc[reg][dz].Draw('A C E3')
        else:
            roc[reg][dz].Draw('C E3 SAME')

        #roc_dT[reg][dz].Draw('C E3 same')
            

        leg[reg].AddEntry(roc[reg][dz],'no MTD - |#Deltaz| < %.01f mm'%dzvals[dz],'L')
        #leg[reg].AddEntry(roc_dT[reg][dz],'with MTD - |#Deltaz| < %.01f mm'%dzvals[dz],'L')

    leg[reg].Draw()
    tt[reg].Draw()    
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    

#save plots
dirname = release+'/RocsComparisonDz/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for reg in regions:
    if ('TTbar' in inputDir):
        c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
        c[reg].SaveAs(dirname+c[reg].GetName()+'.png')
        c[reg].SaveAs(dirname+c[reg].GetName()+'.C')
    if ('QCD' in inputDir):
        c[reg].SaveAs(dirname+c[reg].GetName()+'_QCD.pdf')
        c[reg].SaveAs(dirname+c[reg].GetName()+'_QCD.png')
        c[reg].SaveAs(dirname+c[reg].GetName()+'_QCD.C')

raw_input('ok?')

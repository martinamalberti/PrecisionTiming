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


release = '10_4_0_mtd5'

inputDir = release+'/tracksPuidMva_PU200_TTbar_minTrackPt/'
inputDirNoPU = '93X/30ps_noPU_TTbar_minTkPtCut/'

resolutions = ['30']

minEff = 0.85

dirname = release+'/RocsComparisonTrackMVA/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
tl2.SetNDC()
tl2.SetTextSize(0.035)

f = {}
roc = {}
roc_dT =  {}
roc_3dmva = {}
roc_4dmva =  {}


f_NoPU = {}
roc_NoPU = {}

c = {}

regions = ['barrel','endcap']


col = {'30':ROOT.kRed,
       '40':ROOT.kOrange+1,
       '50':ROOT.kGreen,
       '60':ROOT.kCyan,
       '70':ROOT.kBlue,
}

leg = {}

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.16, '|#eta|<2.8')
tt['barrel'] = ROOT.TLatex( 0.15, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.16, '1.5 < |#eta| < 2.8')
for reg in regions:
    tt[reg].SetNDC()
    tt[reg].SetTextSize(0.035)

   
for ireg,reg in enumerate(regions):
    roc_NoPU[reg] = {}
    roc[reg] = {}
    roc_dT[reg] = {}
    roc_3dmva[reg] = {}
    roc_4dmva[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.70, 0.45, 0.92)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('muonIso_roc_comparison_trackMVA_%s'%reg,'muonIso_roc_comparison_trackMVA_%s'%reg)
    c[reg].SetGridx()
    c[reg].SetGridy()

    for ires,res in enumerate(resolutions):
        f[res] = ROOT.TFile.Open( (inputDir+'/roc_%sps_PU200.root'%(res)).replace('30',res) )
        f_NoPU[res] = ROOT.TFile.Open(inputDirNoPU+'/roc_30ps_noPU.root')

        print f[res].GetName()
        
        gname  = 'g_roc_relChIso03_dZ1_%s'%reg
        gname_dT = 'g_roc_relChIso03_dZ1_dT3s_%s'%reg
        gname_3dmva = 'g_roc_relChIso03_mva3d_%s'%reg
        gname_4dmva = 'g_roc_relChIso03_mva4d_%s'%reg
        gnameNoPU  = 'g_roc_relChIso03_dZ10_simVtx_%s'%reg
            
        roc[reg][res] = f[res].Get(gname) 
        roc_dT[reg][res] = f[res].Get(gname_dT)
        roc_3dmva[reg][res] = f[res].Get(gname_3dmva) 
        roc_4dmva[reg][res] = f[res].Get(gname_4dmva)
        roc_NoPU[reg][res] = f_NoPU[res].Get(gnameNoPU) 

        roc[reg][res].SetLineColor(ROOT.kBlue)
        roc[reg][res].SetLineStyle(1)
        roc[reg][res].SetLineWidth(2)
        roc[reg][res].SetFillStyle(1001)
        roc[reg][res].SetFillColorAlpha(ROOT.kBlue,0.25)

        roc_NoPU[reg][res].SetLineColor(1)
        roc_NoPU[reg][res].SetLineStyle(2)
        roc_NoPU[reg][res].SetLineWidth(2)
        roc_NoPU[reg][res].SetFillStyle(1001)
        roc_NoPU[reg][res].SetFillColorAlpha(ROOT.kBlack,0.0)
        
        roc_dT[reg][res].SetLineColor(col[res])
        roc_dT[reg][res].SetFillStyle(1001)
        roc_dT[reg][res].SetFillColorAlpha(col[res],0.25)
        roc_dT[reg][res].SetLineWidth(2)

        roc_3dmva[reg][res].SetFillStyle(1001)
        roc_3dmva[reg][res].SetFillColorAlpha(ROOT.kGreen+2,0.25)
        roc_3dmva[reg][res].SetLineWidth(2)
        roc_3dmva[reg][res].SetLineStyle(1)

        roc_4dmva[reg][res].SetFillStyle(1001)
        roc_4dmva[reg][res].SetFillColorAlpha(ROOT.kMagenta,0.25)
        roc_4dmva[reg][res].SetLineWidth(2)
        roc_4dmva[reg][res].SetLineStyle(1)
        
        if (ires == 0):
            roc[reg][res].GetXaxis().SetTitle('Prompt efficiency')
            roc[reg][res].GetYaxis().SetTitle('Non-prompt efficiency')
            roc[reg][res].GetXaxis().SetRangeUser(minEff,1.0)
            roc[reg][res].GetYaxis().SetRangeUser(0.0,0.1)
            roc[reg][res].Draw('A C E3')
            leg[reg].AddEntry(roc[reg][res],'no MTD, dz < 1 mm','FL')
            roc_NoPU[reg][res].Draw('C E3 same')
            leg[reg].AddEntry(roc_NoPU[reg][res],'no MTD, noPU','FL')
        roc_dT[reg][res].Draw('C E3 same')
        roc_3dmva[reg][res].Draw('C E3 same')
        roc_4dmva[reg][res].Draw('C E3 same')
        leg[reg].AddEntry(roc_dT[reg][res],'MTD, dz < 1 mm, dt < 3#sigma_{t}','FL')
        leg[reg].AddEntry(roc_3dmva[reg][res],'no MTD, 3D MVA','FL')
        leg[reg].AddEntry(roc_4dmva[reg][res],'MTD, 4D MVA','FL')


    tt[reg].Draw()    
    leg[reg].Draw()
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    c[reg].Update()
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.C')
    raw_input('ok?')



#for reg in regions:
#    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
#    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')
#raw_input('ok?')

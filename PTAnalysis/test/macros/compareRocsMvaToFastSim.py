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

proc = sys.argv[1]

if (proc == 'TTbar'):
    inputDirFastSim = '93X/30ps_PU200_TTbar_minTkPtCut_noDxy/'
    inputDirFullSim = '10_4_0_mtd5/tracksPuidMva_PU200_TTbar_minTrackPt_noDxy_allVtxs/'
if (proc == 'QCD'):
    inputDirFastSim = '93X/30ps_PU200_QCD_minTkPtCut_noDxy/'
    inputDirFullSim = '10_4_0_mtd5/tracksPuidMva_PU200_QCD_minTrackPt_noDxy/'

resolutions = ['40']

minEffPrompt = 0.80
maxEffFake = 0.035

dirname = '10_4_0_mtd5/RocsComparisonTrackMVAFastSim/'
if ('QCD' in inputDirFullSim):
    dirname = '10_4_0_mtd5/RocsComparisonTrackMVAFastSim_QCD/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
#tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
if ('QCD' in inputDirFullSim): tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, QCD')
tl2.SetNDC()
#tl2.SetTextSize(0.035)

fFastSim = {}
fFullSim = {}

roc_fast = {}
roc_dT_fast =  {}
roc = {}
roc_dT = {}
roc_3dmva = {}
roc_4dmva =  {}
roc_3dmva_weighted = {}
roc_4dmva_weighted =  {}


c = {}


leg = {}

regions = ['barrel','endcap']

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.16, '|#eta|<2.8')
tt['barrel'] = ROOT.TLatex( 0.15, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.16, '1.5 < |#eta| < 2.8')
for reg in regions:
    tt[reg].SetNDC()
    #tt[reg].SetTextSize(0.035)

   
for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    roc_3dmva[reg] = {}
    roc_4dmva[reg] = {}
    roc_3dmva_weighted[reg] = {}
    roc_4dmva_weighted[reg] = {}
    roc_fast[reg] = {}
    roc_dT_fast[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.70, 0.45, 0.92)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('muonIso_roc_comparison_trackMVA_%s'%reg,'muonIso_roc_comparison_trackMVA_%s'%reg)
    c[reg].SetGridx()
    c[reg].SetGridy()

    for ires,res in enumerate(resolutions):
        fFastSim[res] = ROOT.TFile.Open( (inputDirFastSim+'/roc_%sps_PU200.root'%(res)).replace('30',res) )
        #fFastSim[res] = ROOT.TFile.Open(  inputDirFullSim+'/roc_35ps_PU200.root' )
        fFullSim[res] = ROOT.TFile.Open(  inputDirFullSim+'/roc_35ps_PU200.root' )
   
        if (reg == 'barrel'):
            gname_fast  = 'g_roc_relChIso03_dZ2_simVtx_%s'%reg
            gname_dT_fast = 'g_roc_relChIso03_dZ2_simVtx_dT3s_%s'%reg
            gname  = 'g_roc_relChIso03_dZ2_%s'%reg
            gname_dT = 'g_roc_relChIso03_dZ2_dT3s_%s'%reg
        else:
            gname_fast  = 'g_roc_relChIso03_dZ3_simVtx_%s'%reg
            gname_dT_fast = 'g_roc_relChIso03_dZ3_simVtx_dT3s_%s'%reg
            gname  = 'g_roc_relChIso03_dZ3_%s'%reg
            gname_dT = 'g_roc_relChIso03_dZ3_dT3s_%s'%reg
            
        gname_3dmva = 'g_roc_relChIso03_mva3d_%s'%reg
        gname_4dmva = 'g_roc_relChIso03_mva4d_%s'%reg

        gname_3dmva_weighted = 'g_roc_relChIso03_mva3d_weighted_%s'%reg
        gname_4dmva_weighted = 'g_roc_relChIso03_mva4d_weighted_%s'%reg
        
        roc_fast[reg][res] = fFastSim[res].Get(gname_fast) 
        roc_dT_fast[reg][res] = fFastSim[res].Get(gname_dT_fast)
        roc[reg][res] = fFullSim[res].Get(gname) 
        roc_dT[reg][res] = fFullSim[res].Get(gname_dT)
        roc_3dmva[reg][res] = fFullSim[res].Get(gname_3dmva) 
        roc_4dmva[reg][res] = fFullSim[res].Get(gname_4dmva)
        roc_3dmva_weighted[reg][res] = fFullSim[res].Get(gname_3dmva_weighted) 
        roc_4dmva_weighted[reg][res] = fFullSim[res].Get(gname_4dmva_weighted)

        roc[reg][res].SetLineColor(ROOT.kGreen+2)
        roc[reg][res].SetLineStyle(2)
        roc[reg][res].SetLineWidth(3)
        roc[reg][res].SetFillStyle(1001)
        roc[reg][res].SetFillColorAlpha(ROOT.kGreen+2,0.0)

        roc_dT[reg][res].SetLineColor(ROOT.kMagenta+2)
        roc_dT[reg][res].SetLineStyle(2)
        roc_dT[reg][res].SetLineWidth(3)
        roc_dT[reg][res].SetFillStyle(1001)
        roc_dT[reg][res].SetFillColorAlpha(ROOT.kMagenta+2,0.0)
        
        roc_fast[reg][res].SetLineColor(ROOT.kBlack)
        roc_fast[reg][res].SetLineStyle(1)
        roc_fast[reg][res].SetLineWidth(3)
        roc_fast[reg][res].SetFillStyle(1001)
        roc_fast[reg][res].SetFillColorAlpha(ROOT.kBlack,0.0)
        
        roc_dT_fast[reg][res].SetLineColor(ROOT.kOrange+7)
        roc_dT_fast[reg][res].SetFillStyle(1001)
        roc_dT_fast[reg][res].SetFillColorAlpha(ROOT.kOrange+7,0.0)
        roc_dT_fast[reg][res].SetLineWidth(3)

        roc_3dmva[reg][res].SetLineColor(ROOT.kGreen+2)
        roc_3dmva[reg][res].SetFillStyle(1001)
        roc_3dmva[reg][res].SetFillColorAlpha(ROOT.kGreen+2,0.0)
        roc_3dmva[reg][res].SetLineWidth(3)
        roc_3dmva[reg][res].SetLineStyle(1)

        roc_4dmva[reg][res].SetLineColor(ROOT.kMagenta+2)
        roc_4dmva[reg][res].SetFillStyle(1001)
        roc_4dmva[reg][res].SetFillColorAlpha(ROOT.kMagenta+2,0.0)
        roc_4dmva[reg][res].SetLineWidth(3)
        roc_4dmva[reg][res].SetLineStyle(1)

        roc_3dmva_weighted[reg][res].SetFillStyle(1001)
        roc_3dmva_weighted[reg][res].SetFillColorAlpha(ROOT.kGreen,0.0)
        roc_3dmva_weighted[reg][res].SetLineWidth(3)
        roc_3dmva_weighted[reg][res].SetLineStyle(2)

        roc_4dmva_weighted[reg][res].SetFillStyle(1001)
        roc_4dmva_weighted[reg][res].SetFillColorAlpha(ROOT.kMagenta,0.0)
        roc_4dmva_weighted[reg][res].SetLineWidth(3)
        roc_4dmva_weighted[reg][res].SetLineStyle(2)
        
        
        roc[reg][res].GetXaxis().SetTitle('Prompt efficiency')
        roc[reg][res].GetYaxis().SetTitle('Non-prompt efficiency')
        roc[reg][res].GetYaxis().SetTitleOffset(1.2)
        roc[reg][res].GetXaxis().SetRangeUser(minEffPrompt,1.0)
        roc[reg][res].GetYaxis().SetRangeUser(0.0,maxEffFake)


        # plot
        roc[reg][res].Draw('A C E3')
        roc_dT[reg][res].Draw('C E3 same')

        #roc_4dmva[reg][res].Draw('C E3 same')
        #roc_3dmva[reg][res].Draw('C E3 same')

        roc_fast[reg][res].Draw('C E3 same')
        roc_dT_fast[reg][res].Draw('C E3 same')
        #roc_3dmva_weighted[reg][res].Draw('C E3 same')
        #roc_4dmva_weighted[reg][res].Draw('C E3 same')

        leg[reg].AddEntry(roc_fast[reg][res],'no MTD, dz (fastSim, 9_3_2)','FL')
        leg[reg].AddEntry(roc_dT_fast[reg][res],'MTD, dz + dt (fastSim, 9_3_2)','FL')
        leg[reg].AddEntry(roc[reg][res],'no MTD, dz (fastSim, 10_4_0)','FL')
        leg[reg].AddEntry(roc_dT[reg][res],'MTD, dz+dt (fastSim, 10_4_0)','FL')
        #leg[reg].AddEntry(roc_3dmva[reg][res],'no MTD, 3D MVA (fullSim)','FL')
        #leg[reg].AddEntry(roc_4dmva[reg][res],'MTD, 4D MVA (fullSim)','FL')
        #leg[reg].AddEntry(roc_3dmva_weighted[reg][res],'no MTD, 3D MVA weighted','FL')
        #leg[reg].AddEntry(roc_4dmva_weighted[reg][res],'MTD, 4D MVA weighted','FL')


    tt[reg].Draw()    
    leg[reg].Draw()
    #tl.Draw()
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

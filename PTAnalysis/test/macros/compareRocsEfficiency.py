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

release = '93X'
inputDir = '93X/30ps_PU200_TTbar_minTkPtCut/'
#inputDirEff = '93X/30ps_PU200_TTbar_minTkPtCut/'
inputDirEff = '93X/30ps_PU200_TTbar_minTkPtCut_btlEff85_etlEff90/'

resolutions = ['40']
regions = ['barrel','endcap']
efficiencies = {}
efficiencies['barrel'] = [ 85, 100 ]
efficiencies['endcap'] = [ 90, 100 ]

dirname = release+'/RocsComparisonEfficiency/'
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

c = {}

col = {'30':ROOT.kRed,
       '40':ROOT.kOrange+1,
       '50':ROOT.kGreen,
       '60':ROOT.kCyan,
       '70':ROOT.kBlue,
}

linestyle = { 100 : 1, 90 : 2, 85 : 2}



leg = {}

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.16, '|#eta|<2.5')
tt['barrel'] = ROOT.TLatex( 0.15, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.16, '1.5 < |#eta| < 2.5 ')
for reg in regions:
    tt[reg].SetNDC()
    tt[reg].SetTextSize(0.035)


for ireg,reg in enumerate(regions):
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.65, 0.45, 0.92)
    leg[reg].SetBorderSize(0)

    c[reg] = ROOT.TCanvas('muonIso_roc_comparison_eff_%s'%reg,'muonIso_roc_comparison_eff_%s'%reg)
    c[reg].SetGridx()
    c[reg].SetGridy()
    
    for ires,res in enumerate(resolutions):
        f[res] = {}
        roc[reg][res] = {}
        roc_dT[reg][res] = {}
        
        for ieff, eff in enumerate(efficiencies[reg]):
            if (eff == 100):
                f[res][eff] = ROOT.TFile.Open( (inputDir+'/roc_%sps_PU200.root'%(res)).replace('30',res) )
            else:
                f[res][eff] = ROOT.TFile.Open( (inputDirEff+'/roc_%sps_PU200.root'%(res)).replace('30',res) )

                                              
            if (reg == 'barrel'):
                gname  = 'g_roc_relChIso03_dZ1_simVtx_%s'%reg
                gname_dT = 'g_roc_relChIso03_dZ1_simVtx_dT3s_%s'%reg
            else:
                gname  = 'g_roc_relChIso03_dZ2_simVtx_%s'%reg
                gname_dT = 'g_roc_relChIso03_dZ2_simVtx_dT3s_%s'%reg
                                              
                                              
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
                roc[reg][res][eff].GetXaxis().SetRangeUser(0.85, 1)
                roc[reg][res][eff].GetYaxis().SetRangeUser(0.0, 0.1)
                roc[reg][res][eff].GetXaxis().SetTitle('Prompt efficiency')
                roc[reg][res][eff].GetYaxis().SetTitle('Non-prompt efficiency')
                roc[reg][res][eff].Draw('A C E3')
                leg[reg].AddEntry(roc[reg][res][eff],'no MTD','FL')
            
            roc_dT[reg][res][eff].Draw('C E3 same')
            leg[reg].AddEntry(roc_dT[reg][res][eff],'MTD, #sigma_{t} = %s ps, #epsilon = %d %%'%(res,eff),'L')

    tt[reg].Draw()    
    leg[reg].Draw()
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    raw_input('ok?')


#save plots
for reg in regions:
    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')

raw_input('ok?')

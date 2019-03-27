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

whichVtx = sys.argv[1]

#inputDir = '93X/30ps_PU200_TTbar_minTkPtCut_muonPt10-20/'
#inputDirNoPU = '93X/30ps_noPU_TTbar_minTkPtCut_muonPt10-20/'
inputDir = '93X/30ps_PU200_TTbar_minTkPtCut/'
inputDirNoPU = '93X/30ps_noPU_TTbar_minTkPtCut/'

resolutions = ['30', '40', '50', '60', '70']
#resolutions = ['40']

canvasWithRatios = False
minEff = 0.85

dirname = release+'/RocsComparisonTimeResolution/'
if (whichVtx == 'simVtx'):
    dirname = release+'/RocsComparisonTimeResolution_simVtx/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
#tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, QCD')
tl2.SetNDC()
tl2.SetTextSize(0.035)

f = {}
roc = {}
roc_dT =  {}

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

    
gdummy = ROOT.TGraph()
gdummy_dT = ROOT.TGraph()

g_prompt_vs_fake = {}
g_prompt_vs_fake_dT = {}
g_ratio_prompt = {}
g_ratio_fake = {}

for ireg,reg in enumerate(regions):
    roc_NoPU[reg] = {}
    roc[reg] = {}
    roc_dT[reg] = {}
    leg[reg] = ROOT.TLegend(0.15, 0.70, 0.45, 0.92)
    leg[reg].SetBorderSize(0)
    
    c[reg] = ROOT.TCanvas('muonIso_roc_comparison_%s'%reg,'muonIso_roc_comparison_%s'%reg)
    if (reg == ''):
        c[reg] = ROOT.TCanvas('muonIso_roc_comparison','muonIso_roc_comparison')
    c[reg].SetGridx()
    c[reg].SetGridy()

    if (canvasWithRatios):
        ## -- inside this canvas, we create two pads
        pad1 = ROOT.TPad("pad1","pad1",0.05,0.25,0.75,0.95);
        pad2 = ROOT.TPad("pad2","pad2",0.71,0.25,0.95,0.95);
        pad3 = ROOT.TPad("pad3","pad3",0.05,0.10,0.75,0.25);
        pad1.Draw()
        pad2.Draw()
        pad3.Draw()

        for pad in pad1, pad2, pad3:
            pad.SetGridx()
            pad.SetGridy()
        
  
    g_prompt_vs_fake[reg] = {}
    g_prompt_vs_fake_dT[reg] = {}
    g_ratio_prompt[reg]  = {}
    g_ratio_fake[reg] = {}

    for ires,res in enumerate(resolutions):
        f[res] = ROOT.TFile.Open( (inputDir+'/roc_%sps_PU200.root'%(res)).replace('30',res) )
        print inputDirNoPU+'/roc_30ps_noPU.root'
        f_NoPU[res] = ROOT.TFile.Open(inputDirNoPU+'/roc_30ps_noPU.root')
        
        if (whichVtx == 'simVtx'):
            if (reg == 'barrel'):
                gname  = 'g_roc_relChIso03_dZ1_simVtx_%s'%reg
                gname_dT = 'g_roc_relChIso03_dZ1_simVtx_dT3s_%s'%reg
                gnameNoPU  = 'g_roc_relChIso03_dZ10_simVtx_%s'%reg
            else:
                gname  = 'g_roc_relChIso03_dZ2_simVtx_%s'%reg
                gname_dT = 'g_roc_relChIso03_dZ2_simVtx_dT3s_%s'%reg
                gnameNoPU  = 'g_roc_relChIso03_dZ10_simVtx_%s'%reg
        else:
            gname  = 'g_roc_relChIso03_dZ1_%s'%reg
            gname_dT = 'g_roc_relChIso03_dZ1_dT3s_%s'%reg

            
            
        roc[reg][res] = f[res].Get(gname) 
        roc_dT[reg][res] = f[res].Get(gname_dT)
        roc_NoPU[reg][res] = f_NoPU[res].Get(gnameNoPU) 

        g_prompt_vs_fake[reg][res] = ROOT.TGraph()
        g_prompt_vs_fake_dT[reg][res] = ROOT.TGraph()
        g_ratio_prompt[reg][res] = ROOT.TGraph()
        g_ratio_fake[reg][res] = ROOT.TGraph()
        
        eff_fake = ROOT.Double(0)
        eff_prompt = ROOT.Double(0)
        eff_fake_dT = ROOT.Double(0)
        eff_prompt_dT = ROOT.Double(0)
        n = 0
        for i in range(0,roc[reg][res].GetN() ):
            roc[reg][res].GetPoint(i,eff_prompt, eff_fake) 
            roc_dT[reg][res].GetPoint(i,eff_prompt_dT, eff_fake_dT)
            g_prompt_vs_fake[reg][res].SetPoint(i,eff_fake,eff_prompt )
            g_prompt_vs_fake_dT[reg][res].SetPoint(i,eff_fake_dT,eff_prompt_dT )
            if (eff_prompt > 0.65):
                g_ratio_prompt[reg][res].SetPoint(n,  g_prompt_vs_fake_dT[reg][res].Eval(eff_fake)/g_prompt_vs_fake[reg][res].Eval(eff_fake) , eff_fake)
                n = n+1
            g_ratio_fake[reg][res].SetPoint(i, eff_prompt,  roc_dT[reg][res].Eval(eff_prompt)/roc[reg][res].Eval(eff_prompt))

        roc[reg][res].SetLineColor(1)
        roc[reg][res].SetLineStyle(1)
        roc[reg][res].SetLineWidth(2)
        roc[reg][res].SetFillStyle(1001)
        roc[reg][res].SetFillColorAlpha(ROOT.kBlack,0.25)

        roc_NoPU[reg][res].SetLineColor(1)
        roc_NoPU[reg][res].SetLineStyle(2)
        roc_NoPU[reg][res].SetLineWidth(2)
        roc_NoPU[reg][res].SetFillStyle(1001)
        roc_NoPU[reg][res].SetFillColorAlpha(ROOT.kBlack,0.0)
        
        roc_dT[reg][res].SetLineColor(col[res])
        roc_dT[reg][res].SetFillStyle(0)
        roc_dT[reg][res].SetFillColorAlpha(col[res],0.0)
        if (res == '30'):
            roc_dT[reg][res].SetFillStyle(1001)
            roc_dT[reg][res].SetFillColorAlpha(col[res],0.0)
        roc_dT[reg][res].SetLineWidth(2)


        
        if (ires == 0):
            roc[reg][res].GetXaxis().SetTitle('Prompt efficiency')
            roc[reg][res].GetYaxis().SetTitle('Non-prompt efficiency')
            roc[reg][res].GetYaxis().SetTitleOffset(1.4)
            if (canvasWithRatios):
                pad1.cd()
            roc[reg][res].GetXaxis().SetRangeUser(minEff,1.0)
            roc[reg][res].GetYaxis().SetRangeUser(0.0,0.1)
            roc[reg][res].Draw('A C E3')
            leg[reg].AddEntry(roc[reg][res],'no MTD','FL')
            roc_NoPU[reg][res].Draw('C E3 same')
            leg[reg].AddEntry(roc_NoPU[reg][res],'no MTD, noPU','PL')
        if (canvasWithRatios):
            pad1.cd()  
        roc_dT[reg][res].Draw('C E3 same')
        leg[reg].AddEntry(roc_dT[reg][res],'MTD, #sigma_{t} = %s ps'%res,'PL')

        if (canvasWithRatios):
            pad2.cd()
            g_ratio_prompt[reg][res].SetLineColor(col[res])
        
            if (ires == 0):
                g_ratio_prompt[reg][res].Draw('ac')
                g_ratio_prompt[reg][res].GetXaxis().SetRangeUser(1.0,1.2)
                g_ratio_prompt[reg][res].GetYaxis().SetRangeUser(0.0,0.1)
                g_ratio_prompt[reg][res].GetXaxis().SetLabelFont(43)
                g_ratio_prompt[reg][res].GetXaxis().SetLabelSize(18)
                g_ratio_prompt[reg][res].GetXaxis().SetNdivisions(205)
            else: g_ratio_prompt[reg][res].Draw('csame')

            pad3.cd()
            g_ratio_fake[reg][res].SetLineColor(col[res])
            if (ires==0):
                g_ratio_fake[reg][res].Draw('ac')
                g_ratio_fake[reg][res].GetXaxis().SetRangeUser(minEff,1.01)
                g_ratio_fake[reg][res].GetYaxis().SetRangeUser(0.70,1.0)
                g_ratio_fake[reg][res].GetXaxis().SetLabelFont(43)
                g_ratio_fake[reg][res].GetXaxis().SetLabelSize(18)
                g_ratio_fake[reg][res].GetYaxis().SetLabelFont(43)
                g_ratio_fake[reg][res].GetYaxis().SetLabelSize(18)
                g_ratio_fake[reg][res].GetYaxis().SetNdivisions(205)
            else: g_ratio_fake[reg][res].Draw('csame')

    if (canvasWithRatios):
        pad1.cd()      

    tt[reg].Draw()    
    leg[reg].Draw()
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    c[reg].Update()
    #c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
    #c[reg].SaveAs(dirname+c[reg].GetName()+'.png')
    #c[reg].SaveAs(dirname+c[reg].GetName()+'.C')
    c[reg].SaveAs(dirname+c[reg].GetName()+'_timeResol.pdf')
    c[reg].SaveAs(dirname+c[reg].GetName()+'_timeResol.png')
    c[reg].SaveAs(dirname+c[reg].GetName()+'_timeResol.C')
#    raw_input('ok?')



#for reg in regions:
#    c[reg].SaveAs(dirname+c[reg].GetName()+'.pdf')
#    c[reg].SaveAs(dirname+c[reg].GetName()+'.png')
#raw_input('ok?')

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


inputDir = '93X/40ps_PU200_TTbar_minTkPtCut_noDxy/'
inputDirRef = '93X/30ps_PU200_TTbar_minTkPtCut_eff100_noDxy/'

f = ROOT.TFile.Open( inputDir+'/roc_40ps_PU200.root' )
f_ref  = ROOT.TFile.Open( inputDirRef+'/roc_30ps_PU200.root' )


canvasWithRatios = True
minEffPrompt = 0.80
maxEffFake = 0.05

dirname = release+'/RocsComparisonToIdealMTD/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
#tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, QCD')
tl2.SetNDC()
tl2.SetTextSize(0.035)

roc = {}
roc_dT =  {}
roc_dT_ref =  {}

c = {}

regions = ['barrel','endcap']


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
g_prompt_vs_fake_dT_ref = {}
g_ratio_prompt = {}
g_ratio_prompt_ref = {}
g_ratio_fake = {}
g_ratio_fake_ref = {}

g_doubleratio_prompt = {}
g_doubleratio_fake = {}
g_ratiofake_vs_derivative = {}


for ireg,reg in enumerate(regions):
    roc_dT[reg] = {}
    roc_dT_ref[reg] = {}
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
        
  
    if (reg == 'barrel'):
        gname = 'g_roc_relChIso03_dZ2_simVtx_%s'%reg
        gname_dT  = 'g_roc_relChIso03_dZ2_simVtx_dT3s_%s'%reg
        gname_dT_ref = 'g_roc_relChIso03_dZ2_simVtx_dT3s_%s'%reg
    else:
        gname  = 'g_roc_relChIso03_dZ3_simVtx_%s'%reg
        gname_dT  = 'g_roc_relChIso03_dZ3_simVtx_dT3s_%s'%reg
        gname_dT_ref = 'g_roc_relChIso03_dZ3_simVtx_dT3s_%s'%reg
                    
            
    roc[reg] = f.Get(gname)
    roc_dT[reg] = f.Get(gname_dT)
    roc_dT_ref[reg] = f_ref.Get(gname_dT_ref)
    
    g_prompt_vs_fake[reg] = ROOT.TGraph()
    g_prompt_vs_fake_dT[reg] = ROOT.TGraph()
    g_prompt_vs_fake_dT_ref[reg] = ROOT.TGraph()
    g_ratio_prompt_ref[reg] = ROOT.TGraph()
    g_ratio_prompt[reg] = ROOT.TGraph()
    g_ratio_fake[reg] = ROOT.TGraph()
    g_ratio_fake_ref[reg] = ROOT.TGraph()
    g_doubleratio_prompt[reg]  = ROOT.TGraph()
    g_doubleratio_prompt[reg].SetName('ratio_prompt_%s'%reg)
    g_doubleratio_fake[reg]  = ROOT.TGraph()
    g_doubleratio_fake[reg].SetName('ratio_fake_%s'%reg)
    g_ratiofake_vs_derivative[reg]= ROOT.TGraph()
    g_ratiofake_vs_derivative[reg].SetName('ratio_fake_vs_derivative_%s'%reg)
      
    eff_fake = ROOT.Double(0)
    eff_prompt = ROOT.Double(0)
    eff_fake_dT = ROOT.Double(0)
    eff_prompt_dT = ROOT.Double(0)
    eff_fake_dT_ref = ROOT.Double(0)
    eff_prompt_dT_ref = ROOT.Double(0)
    n = 0
    for i in range(1,roc_dT_ref[reg].GetN()-1 ):
        roc_dT_ref[reg].GetPoint(i,eff_prompt_dT_ref, eff_fake_dT_ref) 
        roc_dT[reg].GetPoint(i,eff_prompt_dT, eff_fake_dT)
        roc[reg].GetPoint(i,eff_prompt, eff_fake)
        g_prompt_vs_fake[reg].SetPoint(i,eff_fake,eff_prompt )
        g_prompt_vs_fake_dT[reg].SetPoint(i,eff_fake_dT,eff_prompt_dT )
        g_prompt_vs_fake_dT_ref[reg].SetPoint(i,eff_fake_dT_ref,eff_prompt_dT_ref )
        if (eff_prompt > minEffPrompt):
            g_ratio_prompt[reg].SetPoint(n,  g_prompt_vs_fake_dT[reg].Eval(eff_fake)/g_prompt_vs_fake[reg].Eval(eff_fake) , eff_fake)
            g_ratio_prompt_ref[reg].SetPoint(n,  g_prompt_vs_fake_dT_ref[reg].Eval(eff_fake)/g_prompt_vs_fake[reg].Eval(eff_fake) , eff_fake)
            g_doubleratio_prompt[reg].SetPoint(n, eff_fake, g_prompt_vs_fake_dT[reg].Eval(eff_fake)/g_prompt_vs_fake_dT_ref[reg].Eval(eff_fake) )
            n = n+1
        g_ratio_fake[reg].SetPoint(i, eff_prompt,  roc_dT[reg].Eval(eff_prompt)/roc[reg].Eval(eff_prompt))
        g_ratio_fake_ref[reg].SetPoint(i, eff_prompt,  roc_dT_ref[reg].Eval(eff_prompt)/roc[reg].Eval(eff_prompt))
        g_doubleratio_fake[reg].SetPoint(i, eff_prompt,  roc_dT[reg].Eval(eff_prompt)/roc_dT_ref[reg].Eval(eff_prompt))

        x1, x2, y1, y2 = [ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)]
        roc_dT_ref[reg].GetPoint(i-1, x1, y1)
        roc_dT_ref[reg].GetPoint(i+1, x2, y2)
        if (x1!=x2): derivative = (y2-y1)/(x2-x1)
        else: derivative = 9999.
        #print x1, y1, '   ', x2, y2, '  ', derivative
        g_ratiofake_vs_derivative[reg].SetPoint(i, derivative, roc_dT[reg].Eval(eff_prompt)/roc_dT_ref[reg].Eval(eff_prompt) )

        
        
    roc_dT_ref[reg].SetLineColor(1)
    roc_dT_ref[reg].SetLineStyle(1)
    roc_dT_ref[reg].SetLineWidth(2)
    roc_dT_ref[reg].SetFillStyle(1001)
    roc_dT_ref[reg].SetFillColorAlpha(ROOT.kBlack,0.0)
    
    roc_dT[reg].SetLineColor(2)
    roc_dT[reg].SetFillStyle(0)
    roc_dT[reg].SetFillColorAlpha(2,0.0)
    roc_dT[reg].SetLineWidth(2)
    
    roc_dT_ref[reg].GetXaxis().SetTitle('Prompt efficiency')
    roc_dT_ref[reg].GetYaxis().SetTitle('Non-prompt efficiency')
    roc_dT_ref[reg].GetYaxis().SetTitleOffset(1.12)

    if (canvasWithRatios):
        pad1.cd()
        roc[reg].GetXaxis().SetRangeUser(minEffPrompt,1.0)
        roc[reg].GetYaxis().SetRangeUser(0.0,maxEffFake)
        roc[reg].Draw('A C E3')
        leg[reg].AddEntry(roc[reg],'no MTD','FL')

    if (canvasWithRatios):
        pad1.cd()  
        roc_dT_ref[reg].Draw('C E3 same')
        roc_dT[reg].Draw('C E3 same')
        leg[reg].AddEntry(roc_dT_ref[reg],'MTD, #sigma_{t} = 30 ps, 100%','L')
        leg[reg].AddEntry(roc_dT[reg],'MTD, #sigma_{t} = 40 ps','L')
        

    if (canvasWithRatios):
        pad2.cd()
        g_ratio_prompt[reg].SetLineColor(1)
        
        g_ratio_prompt[reg].Draw('ac')
        g_ratio_prompt[reg].GetXaxis().SetRangeUser(0.5,1.1)
        g_ratio_prompt[reg].GetYaxis().SetRangeUser(0.0,maxEffFake)
        g_ratio_prompt[reg].GetXaxis().SetLabelFont(43)
        g_ratio_prompt[reg].GetXaxis().SetLabelSize(18)
        g_ratio_prompt[reg].GetXaxis().SetNdivisions(205)
        g_ratio_prompt_ref[reg].Draw('csame')

        pad3.cd()
        g_ratio_fake[reg].SetLineColor(1)
        g_ratio_fake[reg].Draw('ac')
        g_ratio_fake[reg].GetXaxis().SetRangeUser(minEffPrompt,1.01)
        g_ratio_fake[reg].GetYaxis().SetRangeUser(0.70,1.1)
        g_ratio_fake[reg].GetXaxis().SetLabelFont(43)
        g_ratio_fake[reg].GetXaxis().SetLabelSize(18)
        g_ratio_fake[reg].GetYaxis().SetLabelFont(43)
        g_ratio_fake[reg].GetYaxis().SetLabelSize(18)
        g_ratio_fake[reg].GetYaxis().SetNdivisions(205)
        g_ratio_fake_ref[reg].Draw('csame')
    
    if (canvasWithRatios):
        pad1.cd()      

    tt[reg].Draw()    
    leg[reg].Draw()
    tl.Draw()
    tl2.Draw()
    CMS_lumi.CMS_lumi(c[reg], iPeriod, iPos)
    c[reg].Update()
    
    raw_input('ok?')

outfile = ROOT.TFile('ratioToRescaleTaus.root','recreate')
for reg in 'barrel', 'endcap':
    g_doubleratio_prompt[reg].Write()
    g_doubleratio_fake[reg].Write()
    g_ratiofake_vs_derivative[reg].Write()
outfile.Close()



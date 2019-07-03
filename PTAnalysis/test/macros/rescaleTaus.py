
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


roc_dT =[
    [0.3,             0.006],
    [0.349714,        0.00943396],
    [0.375429,        0.0106132] , 
    [0.399429,        0.0117925] , 
    [0.425143,        0.0141509] , 
    [0.449143,        0.0159198] , 
    [0.473714,        0.0188679] , 
    [0.499429,        0.0235849] , 
    [0.525143,        0.0283019] , 
    [0.550286,        0.0347877] , 
    [0.575429,        0.0465802] , 
    [0.598857,        0.0642689] ,
    [0.605566,        0.0863636] ,
    #[0.608571,        0.0996462] ,
    [0.609082,        0.099697] ,
    [0.616113,        0.119394],
    [0.62,            0.14092], 
    [0.623429,        0.209906],
    [0.624,           0.25059],
    [0.624,           0.277019],
]  


roc = [
    #[0.3,             0.00884434],
    [0.3,             0.0084],
    [0.328,           0.0117925] ,    
    [0.349714,        0.0129717],
    [0.374857,        0.0153302],
    [0.4,             0.0176887],
    [0.426286,        0.021816],
    [0.449714,        0.0253538],
    [0.475429,        0.0294811],
    [0.498857,        0.035967] ,
    [0.525714,        0.0436321],
    [0.549714,        0.0542453],
    [0.574286,        0.0689858],
    [0.599429,        0.0925708],
    [0.605566,        0.104242], 
    [0.612571,        0.117925] ,
    [0.616,           0.149764],
    [0.621714,        0.200472],
    [0.624,           0.245283],
    [0.624,           0.277019],
]


g_roc = ROOT.TGraph()
g_roc_dT = ROOT.TGraph()
g_roc_dT_scaled = ROOT.TGraph()


f_scaling = ROOT.TFile.Open('ratioToRescaleTaus.root')
#f_scaling = ROOT.TFile.Open('ratioToRescaleTaus_50ps.root')
g_scaling = f_scaling.Get('ratio_fake_vs_derivative_endcap')

dirname = '93X/tauRescaling_40ps/'

print len(roc)
for i in range(0,len(roc)):
    eff_sig =  roc[i][0]
    eff_bkg =  roc[i][1]
    g_roc.SetPoint(i, eff_sig, eff_bkg)
    eff_sig_dT =  roc_dT[i][0]
    eff_bkg_dT =  roc_dT[i][1]
    g_roc_dT.SetPoint(i, eff_sig_dT, eff_bkg_dT)


# rescale
for i in range(0,len(roc_dT)):
    x1, x2, y1, y2 = [ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)]
    if (i == 0):
        g_roc_dT.GetPoint(i, x1, y1)
        g_roc_dT.GetPoint(i+1, x2, y2)
    elif (i == len(roc_dT)):
        g_roc_dT.GetPoint(i-1, x1, y1)
        g_roc_dT.GetPoint(i, x2, y2)
    else:
        g_roc_dT.GetPoint(i-1, x1, y1)
        g_roc_dT.GetPoint(i+1, x2, y2)
        
    if (x1!=x2): derivative = (y2-y1)/(x2-x1)
    else: derivative = 9999.

    eff_sig_dT, eff_bkg_dT = [ROOT.Double(0),ROOT.Double(0)]
    g_roc_dT.GetPoint(i, eff_sig_dT, eff_bkg_dT)
    eff_bkg_dT_scaled = eff_bkg_dT * g_scaling.Eval(derivative)
    g_roc_dT_scaled.SetPoint(i, eff_sig_dT, eff_bkg_dT_scaled)


for g in g_roc, g_roc_dT, g_roc_dT_scaled:
    g.SetLineWidth(3)
    
canvas = ROOT.TCanvas('roc_taus_rescaled','roc_taus_rescaled')
canvas.SetGridx()
canvas.SetGridy()
g_roc.SetLineColor(ROOT.kBlue)
g_roc_dT.SetLineColor(ROOT.kRed)
g_roc_dT_scaled.SetLineColor(ROOT.kOrange+7)

g_roc.GetXaxis().SetRangeUser(0.3,0.63)
g_roc.GetYaxis().SetRangeUser(0.0,0.10)
g_roc.GetXaxis().SetTitle('#tau_{h} efficiency')
g_roc.GetYaxis().SetTitle('Jet efficiency')
g_roc.Draw('al')
#g_roc_dT.Draw('lsame')
g_roc_dT_scaled.Draw('lsame')

leg = ROOT.TLegend(0.16, 0.70, 0.50, 0.92)
leg.SetBorderSize(0)
leg.AddEntry(g_roc,'no MTD','L')
#leg.AddEntry(g_roc_dT,'MTD, #sigma_{t} = 30 ps, MTD #epsilon = 100%','L')
#leg.AddEntry(g_roc_dT_scaled,'MTD, #sigma_{t} = 40 ps, #epsilon = 85(90)% BTL(ETL)','L')
leg.AddEntry(g_roc_dT_scaled,'MTD, #sigma_{t} = 40 ps','L')

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#tau#tau, QCD')
tl2.SetNDC()
#tl2.SetTextSize(0.035)

leg.Draw()
#tl.Draw()
tl2.Draw()

CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


canvas.SaveAs(dirname+canvas.GetName()+'.pdf')
canvas.SaveAs(dirname+canvas.GetName()+'.png')
canvas.SaveAs(dirname+canvas.GetName()+'.C')
raw_input('ok?')

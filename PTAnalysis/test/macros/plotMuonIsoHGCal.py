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

def makeRoc(hsig, hbkg, graph):
    nbins = hsig.GetNbinsX()
    for ibin in range(1, nbins+1):
        nPassSig = hsig.Integral(1,ibin)
        nSig     = hsig.Integral(1,nbins+1)
        nPassBkg = hbkg.Integral(1,ibin)
        nBkg     = hbkg.Integral(1,nbins+1)
        
        effSig = nPassSig/nSig
        effBkg = nPassBkg/nBkg

        effSigErrDown =  effSig - ROOT.TEfficiency.ClopperPearson(nSig, nPassSig, 0.68, False)
        effSigErrUp   =  ROOT.TEfficiency.ClopperPearson(nSig, nPassSig, 0.68, True) - effSig

        effBkgErrDown =  effBkg - ROOT.TEfficiency.ClopperPearson(nBkg, nPassBkg, 0.68, False)
        effBkgErrUp   =  ROOT.TEfficiency.ClopperPearson(nBkg, nPassBkg, 0.68, True) - effBkg

        graph.SetPoint(ibin-1, effSig, effBkg)    
        graph.SetPointError(ibin-1, effSigErrDown, effSigErrUp, effBkgErrDown, effBkgErrUp)    
        
    return


f = {}
fname = {}
#fname['sig'] = '../93X/testTracksHGCal_DYToLL_PU200.root'
#fname['bkg'] = '../93X/testTracksHGCal_TTbar_PU200.root'
fname['sig'] = '../93X/testTracksHGCal_DYToLL_PU200_noDxy.root'
fname['bkg'] = '../93X/testTracksHGCal_TTbar_PU200_noDxy.root'
#fname['sig'] = '../93X/testTracksHGCal_DYToLL_noPU.root'
#fname['bkg'] = '../93X/testTracksHGCal_TTbar_noPU.root'

configs = ['','HGCal','ETL2','ETL1','HGCal_ETL1','HGCal_ETL2'] 
col = {''          : ROOT.kBlack,
       'HGCal'     : ROOT.kBlue,
       'ETL2'      : ROOT.kRed,
       'ETL1'      : ROOT.kRed+2,
       'HGCal_ETL1': ROOT.kGreen+2,
       'HGCal_ETL2': ROOT.kGreen+1}

styl = {''          : 20,
        'HGCal'     : 21,
        'ETL2'      : 20,
        'ETL1'      : 20,
        'HGCal_ETL1': 21,
        'HGCal_ETL2': 22}

minEffPrompt = 0.75
maxEffFake = 0.05


h_relChIso = {}
p_efficiency_vs_linedensity = {}
p_efficiency_vs_muonpt = {}
g_roc ={}

for proc in 'sig', 'bkg':
    f[proc] = ROOT.TFile.Open(fname[proc])
    h_relChIso[proc] = {}
    p_efficiency_vs_linedensity[proc] = {}
    p_efficiency_vs_muonpt[proc] = {}
    for config in configs:
        if (config == ''):
            h_relChIso[proc][config] = f[proc].Get('h_muon_relChIso03_dZ2')
            p_efficiency_vs_linedensity[proc][config] = f[proc].Get('p_efficiency_relChIso03_dZ2_vs_linedensity')
            p_efficiency_vs_muonpt[proc][config] = f[proc].Get('p_efficiency_relChIso03_dZ2_vs_muonpt')
        else:
            h_relChIso[proc][config] = f[proc].Get('h_muon_relChIso03_dZ2_%s'%config)
            p_efficiency_vs_linedensity[proc][config] = f[proc].Get('p_efficiency_relChIso03_dZ2_%s_vs_linedensity'%config)
            p_efficiency_vs_muonpt[proc][config] = f[proc].Get('p_efficiency_relChIso03_dZ2_%s_vs_muonpt'%config)


for config in configs:
    g_roc[config] =  ROOT.TGraphAsymmErrors()
    g_roc[config].SetName('g_roc_%s'%config)
    makeRoc(h_relChIso['sig'][config], h_relChIso['bkg'][config], g_roc[config])
    
#plot
leg = ROOT.TLegend(0.15, 0.70, 0.45, 0.92)
leg.SetBorderSize(0)

tsig  = ROOT.TLatex( 0.60, 0.86,'<PU> = 200, Z#rightarrow#mu#mu')
tbkg  = ROOT.TLatex( 0.60, 0.86,'<PU> = 200, t#bar{t}')
tcut  = ROOT.TLatex( 0.60, 0.82,'rel chIso < 0.08')
for tt in tsig, tbkg, tcut:
    tt.SetNDC()
    tt.SetTextSize(0.035)


# vs line density
c_prompt = ROOT.TCanvas('eff_prompt_vs_linedensity_endcap_hgcal','eff_prompt_vs_linedensity_endcap_hgcal')
c_prompt.SetGridx()
c_prompt.SetGridy()
for config in configs:
    p_efficiency_vs_linedensity['sig'][config].SetLineColor(col[config])
    p_efficiency_vs_linedensity['sig'][config].SetMarkerColor(col[config])
    p_efficiency_vs_linedensity['sig'][config].SetMarkerStyle(20)
    p_efficiency_vs_linedensity['sig'][config].GetYaxis().SetRangeUser(minEffPrompt,1.1)
    p_efficiency_vs_linedensity['sig'][config].GetYaxis().SetTitle('Prompt efficiency')
    p_efficiency_vs_linedensity['sig'][config].GetXaxis().SetTitle('Line density (mm^{-1})')

    if (config==''):
        p_efficiency_vs_linedensity['sig'][config].Draw('e')
    else:
        p_efficiency_vs_linedensity['sig'][config].Draw('esame')

leg.AddEntry(p_efficiency_vs_linedensity['sig'][''],'no timing','PL')
leg.AddEntry(p_efficiency_vs_linedensity['sig']['HGCal'],'HGCal alone','PL')
leg.AddEntry(p_efficiency_vs_linedensity['sig']['ETL1'],'MTD 1-disk','PL')
leg.AddEntry(p_efficiency_vs_linedensity['sig']['ETL2'],'MTD 2-disks','PL')
leg.AddEntry(p_efficiency_vs_linedensity['sig']['HGCal_ETL1'],'HGCal + MTD 1-disk','PL')
leg.AddEntry(p_efficiency_vs_linedensity['sig']['HGCal_ETL2'],'HGCal + MTD 2-disks','PL')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
CMS_lumi.CMS_lumi(c_prompt, iPeriod, iPos)



c_fake = ROOT.TCanvas('eff_fake_vs_linedensity_endcap_hgcal','eff_fake_vs_linedensity_endcap_hgcal')
c_fake.SetGridx()
c_fake.SetGridy()
for config in configs:
    p_efficiency_vs_linedensity['bkg'][config].SetLineColor(col[config])
    p_efficiency_vs_linedensity['bkg'][config].SetMarkerColor(col[config])
    p_efficiency_vs_linedensity['bkg'][config].SetMarkerStyle(20)
    p_efficiency_vs_linedensity['bkg'][config].GetYaxis().SetRangeUser(0.0,maxEffFake)
    p_efficiency_vs_linedensity['bkg'][config].GetYaxis().SetTitle('Non-prompt efficiency')
    p_efficiency_vs_linedensity['bkg'][config].GetXaxis().SetTitle('Line density (mm^{-1})')

    if (config==''):
        p_efficiency_vs_linedensity['bkg'][config].Draw('e')
    else:
        p_efficiency_vs_linedensity['bkg'][config].Draw('esame')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
CMS_lumi.CMS_lumi(c_fake, iPeriod, iPos)


# vs pT
c_prompt_pt = ROOT.TCanvas('eff_prompt_vs_muonpt_endcap_hgcal','eff_prompt_vs_muonpt_endcap_hgcal')
c_prompt_pt.SetGridx()
c_prompt_pt.SetGridy()
for config in configs:
    p_efficiency_vs_muonpt['sig'][config].SetLineColor(col[config])
    p_efficiency_vs_muonpt['sig'][config].SetMarkerColor(col[config])
    p_efficiency_vs_muonpt['sig'][config].SetMarkerStyle(20)
    p_efficiency_vs_muonpt['sig'][config].GetXaxis().SetRangeUser(0.0,70.0)
    p_efficiency_vs_muonpt['sig'][config].GetYaxis().SetRangeUser(0.60,1.1)
    p_efficiency_vs_muonpt['sig'][config].GetYaxis().SetTitle('Prompt efficiency')
    p_efficiency_vs_muonpt['sig'][config].GetXaxis().SetTitle('Line density (mm^{-1})')
    if (config==''):
        p_efficiency_vs_muonpt['sig'][config].Draw('e')
    else:
        p_efficiency_vs_muonpt['sig'][config].Draw('esame')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
CMS_lumi.CMS_lumi(c_prompt_pt, iPeriod, iPos)


c_fake_pt = ROOT.TCanvas('eff_fake_vs_muonpt_endcap_hgcal','eff_fake_vs_muonpt_endcap_hgcal')
c_fake_pt.SetGridx()
c_fake_pt.SetGridy()
for config in configs:
    p_efficiency_vs_muonpt['bkg'][config].SetLineColor(col[config])
    p_efficiency_vs_muonpt['bkg'][config].SetMarkerColor(col[config])
    p_efficiency_vs_muonpt['bkg'][config].SetMarkerStyle(20)
    p_efficiency_vs_muonpt['bkg'][config].GetXaxis().SetRangeUser(0.0,70.)
    p_efficiency_vs_muonpt['bkg'][config].GetYaxis().SetRangeUser(0.0,0.1)
    p_efficiency_vs_muonpt['bkg'][config].GetYaxis().SetTitle('Non-prompt efficiency')
    p_efficiency_vs_muonpt['bkg'][config].GetXaxis().SetTitle('Line density (mm^{-1})')
    if (config==''):
        p_efficiency_vs_muonpt['bkg'][config].Draw('e')
    else:
        p_efficiency_vs_muonpt['bkg'][config].Draw('esame')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
CMS_lumi.CMS_lumi(c_fake_pt, iPeriod, iPos)



# plot ROCS
c_roc =  ROOT.TCanvas('muonIsolation_roc_endcap_hgcal','muonIsolation_roc_endcap_hgcal')
c_roc.SetGridx()
c_roc.SetGridy()
for config in configs:
    g_roc[config].SetLineWidth(2)
    g_roc[config].SetFillStyle(1001)
    g_roc[config].SetFillColorAlpha(ROOT.kBlack,0.0)
    g_roc[config].SetMarkerStyle(20)
    g_roc[config].SetMarkerSize(0.2)
    g_roc[config].SetMarkerColor(col[config])
    g_roc[config].SetLineColor(col[config])
    g_roc[config].GetXaxis().SetRangeUser(minEffPrompt, 1.0)
    g_roc[config].GetYaxis().SetRangeUser(0.0, maxEffFake)
    g_roc[config].GetXaxis().SetTitle('Prompt efficiency')
    g_roc[config].GetYaxis().SetTitle('Non-prompt efficiency')
    if (config==''):
        g_roc[config].Draw('A C E3')
    else:
        g_roc[config].Draw('C E3 same')
leg.Draw('same')


dirname = '93X/'

for canvas in c_prompt, c_fake, c_prompt_pt, c_fake_pt, c_roc:
    canvas.SaveAs(dirname+canvas.GetName()+'.png')
    canvas.SaveAs(dirname+canvas.GetName()+'.pdf')
    #canvas.SaveAs(dirname+canvas.GetName()+'.C')


raw_input('ok?')

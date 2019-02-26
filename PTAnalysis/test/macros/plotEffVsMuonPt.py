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



release = '93X'

#resol = [30, 40, 50, 60, 70]
resol = [30]




f = {}

h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel = {}
h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel = {}
h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel = {}

h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap= {}
h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap = {}
h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap = {}

nRe = 1

for i,res in enumerate(resol):
    f[res] = {}
    h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[res] = {}
    h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[res] = {}
    h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[res] = {}
    h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[res] = {}
    h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap[res] = {}
    h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[res] = {}

    for proc in 'sig', 'bkg':
        fname = '../93X/output_muIso_DYToLL_%dps_prompt_minTkPtCut_muonPt10-9999.root'%res
        if (proc == 'bkg'): fname = '../93X/output_muIso_TTbar_%dps_fake_minTkPtCut_muonPt10-9999.root'%res
        #fname = '../93X/output_muIso_DYToLL_%dps_prompt_minTkPtCut_muonPt10-9999_relIsoCut0.1.root'%res
        #if (proc == 'bkg'): fname = '../93X/output_muIso_TTbar_%dps_fake_minTkPtCut_muonPt10-9999_relIsoCut0.1.root'%res
        f[res][proc] = ROOT.TFile.Open(fname)

        h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ1_simVtx_barrel')
        h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ1_dT3s_simVtx_barrel')
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ2_dT3s_simVtx_barrel')
        h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ2_simVtx_endcap')
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ2_dT3s_simVtx_endcap')
        h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[res][proc] = f[res][proc].Get('p_relChIsoEfficiency_vs_muonpt_dZ3_dT3s_simVtx_endcap')

        
        h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[res][proc].SetLineColor(ROOT.kBlue)
        h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[res][proc].SetLineColor(ROOT.kRed+i)
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[res][proc].SetLineColor(ROOT.kRed+i)

        h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[res][proc].SetLineColor(ROOT.kBlue)
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap[res][proc].SetLineColor(ROOT.kRed+i)
        h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[res][proc].SetLineColor(ROOT.kRed+i)

        h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[res][proc].SetMarkerColor(ROOT.kBlue)
        h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[res][proc].SetMarkerColor(ROOT.kRed+i)
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[res][proc].SetMarkerColor(ROOT.kRed+i)

        h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[res][proc].SetMarkerColor(ROOT.kBlue)
        h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap[res][proc].SetMarkerColor(ROOT.kRed+i)
        h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[res][proc].SetMarkerColor(ROOT.kRed+i)


        for h in h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[res][proc], h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[res][proc], h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[res][proc], h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[res][proc],  h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_endcap[res][proc], h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[res][proc]:
            h.SetMarkerStyle(20)
            h.SetLineWidth(2)


#plot
leg = ROOT.TLegend(0.15, 0.70, 0.45, 0.92)
leg.SetBorderSize(0)

tsig  = ROOT.TLatex( 0.60, 0.86,'<PU> = 200, Z#rightarrow#mu#mu')
tbkg  = ROOT.TLatex( 0.60, 0.86,'<PU> = 200, t#bar{t}')
tcut  = ROOT.TLatex( 0.60, 0.82,'rel chIso < 0.05')
for tt in tsig, tbkg, tcut:
    tt.SetNDC()
    tt.SetTextSize(0.035)

tt = {}
tt['barrel'] = ROOT.TLatex( 0.65, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.65, 0.16, '1.5 < |#eta| < 2.8')
for reg in 'barrel','endcap':
    tt[reg].SetNDC()
    tt[reg].SetTextSize(0.035)
    
c_barrel_prompt = ROOT.TCanvas('eff_prompt_vs_muonpt_barrel','eff_prompt_vs_muonpt_barrel')
c_barrel_prompt.SetGridx()
c_barrel_prompt.SetGridy()
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'].GetXaxis().SetRangeUser(0.0,70.0)
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'].GetYaxis().SetRangeUser(0.60,1.1)
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'].GetYaxis().SetTitle('Prompt efficiency')
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'].GetXaxis().SetTitle('Muon p_{T} (GeV)')
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'].Draw('e')
h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[30]['sig'].Draw('esame')
leg.AddEntry(h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['sig'],'no MTD','PL')
leg.AddEntry(h_relChIsoEff_vs_muonpt_dZ1_dT3s_simVtx_barrel[30]['sig'],'MTD, #sigma_{t} = 30 ps','PL')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
tt['barrel'].Draw()
CMS_lumi.CMS_lumi(c_barrel_prompt, iPeriod, iPos)

c_barrel_fake = ROOT.TCanvas('eff_fake_vs_muonpt_barrel','eff_fake_vs_muonpt_barrel')
c_barrel_fake.SetGridx()
c_barrel_fake.SetGridy()
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['bkg'].GetXaxis().SetRangeUser(0.0,70.0)
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['bkg'].GetYaxis().SetRangeUser(0.0,0.1)
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['bkg'].GetYaxis().SetTitle('Non-prompt efficiency')
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['bkg'].GetXaxis().SetTitle('Muon p_{T} (GeV)')
h_relChIsoEff_vs_muonpt_dZ1_simVtx_barrel[30]['bkg'].Draw('e')
h_relChIsoEff_vs_muonpt_dZ2_dT3s_simVtx_barrel[30]['bkg'].Draw('esame')
leg.Draw('same')
tbkg.Draw()
tcut.Draw()
tt['barrel'].Draw()
CMS_lumi.CMS_lumi(c_barrel_fake, iPeriod, iPos)

c_endcap_prompt = ROOT.TCanvas('eff_prompt_vs_muonpt_endcap','eff_prompt_vs_muonpt_endcap')
c_endcap_prompt
c_endcap_prompt
c_endcap_prompt.SetGridx()
c_endcap_prompt.SetGridy()
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['sig'].GetXaxis().SetRangeUser(0.0,70.0)
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['sig'].GetYaxis().SetRangeUser(0.60,1.1)
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['sig'].GetYaxis().SetTitle('Prompt efficiency')
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['sig'].GetXaxis().SetTitle('Muon p_{T} (GeV)')
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['sig'].Draw('e')
h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[30]['sig'].Draw('esame')
leg.Draw('same')
tsig.Draw()
tcut.Draw()
tt['endcap'].Draw()
CMS_lumi.CMS_lumi(c_endcap_prompt, iPeriod, iPos)

c_endcap_fake = ROOT.TCanvas('eff_fake_vs_muonpt_endcap','eff_fake_vs_muonpt_endcap')
c_endcap_fake.SetGridx()
c_endcap_fake.SetGridy()
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['bkg'].GetXaxis().SetRangeUser(0.0,70.0)
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['bkg'].GetYaxis().SetRangeUser(0.0,0.1)
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['bkg'].GetYaxis().SetTitle('Non-prompt efficiency')
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['bkg'].GetXaxis().SetTitle('Muon p_{T} (GeV)')
h_relChIsoEff_vs_muonpt_dZ2_simVtx_endcap[30]['bkg'].Draw('e')
h_relChIsoEff_vs_muonpt_dZ3_dT3s_simVtx_endcap[30]['bkg'].Draw('esame')
leg.Draw('same')
tbkg.Draw()
tcut.Draw()
tt['endcap'].Draw()
CMS_lumi.CMS_lumi(c_endcap_fake, iPeriod, iPos)


dirname = release+'/'

for canvas in c_barrel_prompt, c_endcap_prompt, c_barrel_fake, c_endcap_fake:
    canvas.SaveAs(dirname+canvas.GetName()+'.png')
    canvas.SaveAs(dirname+canvas.GetName()+'.pdf')
    canvas.SaveAs(dirname+canvas.GetName()+'.C')


raw_input('ok?')


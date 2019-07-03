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
ROOT.gROOT.ForceStyle()

f = {}
f1 = {}
f2 = {}

f1['sig'] = ROOT.TFile.Open('../93X/testTracks_DYToLL_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy.root')
f1['bkg'] = ROOT.TFile.Open('../93X/testTracks_TTbar_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy.root')

f2['sig'] =  ROOT.TFile.Open('../10_4_0_mtd5/testTracksMva_DYToLL_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy_allVtxs.root')
f2['bkg'] =  ROOT.TFile.Open('../10_4_0_mtd5/testTracksMva_TTbar_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy_allVtxs.root')

dirname = 'ComparisonReleases/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

tl = ROOT.TLatex( 0.65, 0.88,'<PU> = 200')
tl.SetNDC()
#tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.65, 0.84,'Z#rightarrow#mu#mu, t#bar{t}')
tl2.SetNDC()
#tl2.SetTextSize(0.035)


hnames = ['h_muon_pt',
          'h_tracks_dz_vtx',
          'h_tracks_dt_vtx',
          'h_tracks_genmatched_dz_vtx',
          'h_tracks_genmatched_dt_vtx',
          'h_tracks_pt',
          'p_tracks_removed_pt_vs_linedensity',
          'p_tracks_removed_n_vs_linedensity',
          'p_tracks_removed_sumpt_vs_linedensity',
]


xtitle = [ 'muon p_{T} (GeV)',
           'dz(track, PV) (cm)',
           'dt(track, PV) (cm)',
           'dz(track, PV) (cm)',
           'dt(track, PV) (cm)',
           'track p_{T} (GeV)',
           'line density (mm^{-1})',
           'line density (mm^{-1})',
           'line density (mm^{-1})',          
]


c = {}
leg = ROOT.TLegend(0.16, 0.80, 0.45, 0.92)
leg.SetBorderSize(0)

regions = ['barrel','endcap']

tt = {}
tt[''] = ROOT.TLatex( 0.15, 0.16, '|#eta|<2.8')
tt['barrel'] = ROOT.TLatex( 0.15, 0.16, '|#eta|<1.5')
tt['endcap'] = ROOT.TLatex( 0.15, 0.16, '1.5 < |#eta| < 2.8')
for reg in regions:
    tt[reg].SetNDC()

h1 = {}
h2 = {}
c = {}

procs = ['sig','bkg']
#procs = ['bkg']

for proc in procs:
    h1[proc] = {}
    h2[proc] = {}
    c[proc] = {}
    for reg in ['barrel','endcap']:
        h1[proc][reg] = {}
        h2[proc][reg] = {}
        c[proc][reg] = {}
        for ih,hname in enumerate(hnames):
            print hname.replace('h_','c_')
            cname = (hname.replace('h_','c_').replace('p_','c_'))+'_%s'%reg+'_%s'%proc 
            c[proc][reg][hname] = ROOT.TCanvas(cname,cname)

            h1[proc][reg][hname] = f1[proc].Get(hname+'_%s'%reg)
            h2[proc][reg][hname] = f2[proc].Get(hname+'_%s'%reg)

            h1[proc][reg][hname].SetMarkerColor(1)
            h2[proc][reg][hname].SetMarkerColor(2)
            
            h1[proc][reg][hname].SetLineColor(1)
            h2[proc][reg][hname].SetLineColor(2)

            h1[proc][reg][hname].SetLineWidth(2)
            h2[proc][reg][hname].SetLineWidth(2)
            
            if ('dz' in hname):
                #h1[proc][reg][hname].Rebin(2)
                #h2[proc][reg][hname].Rebin(2)
                h2[proc][reg][hname].GetXaxis().SetRangeUser(-0.2, 0.2)
                if (reg == 'endcap'):
                    #h1[proc][reg][hname].Rebin(4)
                    #h2[proc][reg][hname].Rebin(4)
                    h2[proc][reg][hname].GetXaxis().SetRangeUser(-0.3, 0.3)
                
            h2[proc][reg][hname].GetYaxis().SetRangeUser(0.0001,h2[proc][reg][hname].GetMaximum()*1.5)
            h2[proc][reg][hname].GetXaxis().SetTitle(xtitle[ih])
            h2[proc][reg][hname].DrawNormalized('histo')
            h1[proc][reg][hname].DrawNormalized('histo same')
            if (hname.startswith('p_')):
                h2[proc][reg][hname].DrawNormalized('e')
                h1[proc][reg][hname].DrawNormalized('e same')
            
            if (proc == procs[0] and reg == 'barrel' and  hname == hnames[0]):
                leg.AddEntry(h1[proc][reg][hname],'9_3_2','L')
                leg.AddEntry(h2[proc][reg][hname],'10_4_0','L')
            
            tt[reg].Draw()    
            leg.Draw()
            CMS_lumi.CMS_lumi(c[proc][reg][hname], iPeriod, iPos)
            c[proc][reg][hname].Update()
            c[proc][reg][hname].SaveAs(dirname+c[proc][reg][hname].GetName()+'.pdf')
            c[proc][reg][hname].SaveAs(dirname+c[proc][reg][hname].GetName()+'.png')
            raw_input('ok?')



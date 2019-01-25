#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT

ROOT.gStyle.SetOptTitle(0)

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetStatX(0.6);
#ROOT.gStyle.SetStatY(0.6);
#ROOT.gStyle.SetStatW(0.3);
#ROOT.gStyle.SetStatH(0.3);

ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetFitFormat("4.4g")

ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.0,'Y')
#ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


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

def makeEff(h, heff):
    nbins = h.GetNbinsX()
    for ibin in range(1, nbins+1):
        x = h.GetBinCenter(ibin)
        nPassSig = h.Integral(1,ibin)
        nSig     = h.Integral(1,nbins+1)
        effSig = nPassSig/nSig
        heff.SetBinContent(ibin, effSig)
    return



#conesDR = ['02','03','04','05']
conesDR = ['03']


resolution = int(sys.argv[1])
pu = sys.argv[2]
bkgProc = sys.argv[3]
suffix = sys.argv[4]

if (pu == 'noPU'):
    fSig = '../output_muIso_DYToLL_noPU_%dps_prompt.root '%resolution
    fBkg = '../output_muIso_%s_noPU_%dps_fake.root' %(bkgProc,resolution)
    #fSig = '../10_4_0_mtd3/output_muIso_DYToLL_noPU_%dps_prompt_1040mtd3.root '%resolution
    #fBkg = '../10_4_0_mtd3/output_muIso_%s_noPU_%dps_fake_1040mtd3.root' %(bkgProc,resolution)
else:
    fSig = '../output_muIso_DYToLL_%dps_prompt.root '%resolution
    fBkg = '../output_muIso_%s_%dps_fake.root'  %(bkgProc,resolution)
    #fSig = '../output_muIso_DYToLL_%dps_prompt_eff90.root '%resolution
    #fBkg = '../output_muIso_%s_%dps_fake_eff90.root'  %(bkgProc,resolution)
    #fSig = '../output_muIso_DYToLL_%dps_prompt_ZTClosestVtx.root '%resolution
    #fBkg = '../output_muIso_%s_%dps_fake_ZTClosestVtx.root'  %(bkgProc,resolution)

    
    
print fSig
print fBkg

tl = ROOT.TLatex( 0.70, 0.84,'<PU> = %s'%pu.replace('PU',''))
if (pu == 'noPU'):
    tl = ROOT.TLatex( 0.70, 0.84,'<PU> = 0')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, t#bar{t}')
if ('QCD' in fBkg):
    tl2 = ROOT.TLatex( 0.70, 0.78,'Z#rightarrow#mu#mu, QCD')
tl2.SetNDC()
tl2.SetTextSize(0.035)


textProc = 't#bar{t}'
if ('QCD' in fBkg):
    textProc = 'QCD'


f = {}
f['sig'] = ROOT.TFile.Open(fSig)
f['bkg'] = ROOT.TFile.Open(fBkg)


# control plots
h_pt = {}
h_eta = {}
h_phi = {}
h_chIsoRatio = {}
h_chIsoRatio_barrel = {}
h_chIsoRatio_endcap = {}
h_vtx_dz3D = {}
h_vtx_dz4D = {}
h_vtx_dt4D = {}
h_vtx_dz3D_pull = {}
h_vtx_dz4D_pull = {}
h_vtx_dt4D_pull = {}


c = {}

for proc in 'sig', 'bkg':
    h_pt[proc]  =  f[proc].Get('h_muon_pt')
    h_eta[proc] =  f[proc].Get('h_muon_eta')
    h_phi[proc] =  f[proc].Get('h_muon_phi')
    h_chIsoRatio[proc] = f[proc].Get('h_muon_relChIso03_ratio')
    h_chIsoRatio_barrel[proc] = f[proc].Get('h_muon_relChIso03_ratio_barrel')
    h_chIsoRatio_endcap[proc] = f[proc].Get('h_muon_relChIso03_ratio_endcap')
    h_vtx_dz3D[proc] = f[proc].Get('h_vtx_dz3D')
    h_vtx_dz4D[proc] = f[proc].Get('h_vtx_dz4D')
    h_vtx_dt4D[proc] = f[proc].Get('h_vtx_dt4D')
    h_vtx_dz3D_pull[proc] = f[proc].Get('h_vtx_dz3D_pull')
    h_vtx_dz4D_pull[proc] = f[proc].Get('h_vtx_dz4D_pull')
    h_vtx_dt4D_pull[proc] = f[proc].Get('h_vtx_dt4D_pull')

    
leg1 = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
leg1.SetBorderSize(0)

for ih,h in  enumerate([h_pt, h_eta, h_phi, h_chIsoRatio,  h_chIsoRatio_barrel, h_chIsoRatio_endcap, h_vtx_dz3D, h_vtx_dz4D, h_vtx_dt4D]):
    cname = h['sig'].GetName().replace('h_','')
    c[cname] = ROOT.TCanvas(cname,cname,500,500)
    if ('relChIso' in cname or 'vtx' in cname):
        c[cname].SetLogy()
    h['sig'].SetLineColor(ROOT.kBlue)
    h['bkg'].SetLineColor(ROOT.kRed)
    h['sig'].GetYaxis().SetRangeUser(0,h['sig'].GetMaximum()*1.5)
    if ('relChIso' in cname or 'vtx' in cname):
        h['sig'].GetYaxis().SetRangeUser(1.0,h['sig'].GetMaximum()*10.)
    h['sig'].DrawNormalized()
    h['bkg'].DrawNormalized('histo same')
    if (ih == 0):
        leg1.AddEntry(h_pt['sig'],'prompt muons (Z#rightarrow#mu#mu)','L')
        leg1.AddEntry(h_pt['bkg'],'non-prompt muons (%s)'%textProc,'L')
        leg1.Draw('same')
    else:
        leg1.Draw('same')
    tl.Draw()
    #tl2.Draw()

    if ('relChIso' in cname):
        bin = h['sig'].FindBin(1)
        fracSig = h['sig'].Integral(0,bin-1) / h['sig'].Integral(0, h['sig'].GetNbinsX()+1 )
        fracBkg = h['bkg'].Integral(0,bin-1) / h['bkg'].Integral(0, h['bkg'].GetNbinsX()+1 )
        print 'Fraction of signal events with Ratio(chIso) < 1 : %.03f'%fracSig
        print 'Fraction of background events with Ratio(chIso) < 1 : %.03f'%fracBkg

        fracSig = h['sig'].Integral(bin+1, h['sig'].GetNbinsX()+1) / h['sig'].Integral(0, h['sig'].GetNbinsX()+1 )
        fracBkg = h['bkg'].Integral(bin+1, h['bkg'].GetNbinsX()+1) / h['bkg'].Integral(0, h['bkg'].GetNbinsX()+1 )
        print 'Fraction of signal events with Ratio(chIso) > 1 : %.03f'%fracSig
        print 'Fraction of background events with Ratio(chIso) > 1 : %.03f'%fracBkg

    #raw_input('ok?')


## --- vtx pulls
for ih,h in  enumerate([h_vtx_dz3D_pull, h_vtx_dz4D_pull, h_vtx_dt4D_pull]):
    cname = h['sig'].GetName().replace('h_','')
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptFit(1)
    c[cname] = ROOT.TCanvas(cname,cname,500,500)
    c[cname].SetLogy()
    h['sig'].GetYaxis().SetRangeUser(1.0,h['sig'].GetMaximum()*5.0)
    h['sig'].SetLineColor(ROOT.kBlue)
    h['bkg'].SetLineColor(ROOT.kRed)
    h['sig'].SetMarkerColor(ROOT.kBlue)
    h['bkg'].SetMarkerColor(ROOT.kRed)
    h['sig'].Draw('s')
    c[cname].Update()
    statsbox = h['sig'].GetListOfFunctions().FindObject('stats')
    statsbox.SetY1NDC(0.55)
    statsbox.SetY2NDC(0.74)
    statsbox.SetX1NDC(0.70)
    statsbox.SetX2NDC(0.89)
    statsbox.SetTextColor(ROOT.kBlue)
    h['bkg'].Draw('s')
    c[cname].Modified()
    c[cname].Update()
    statsbox2 = h['bkg'].GetListOfFunctions().FindObject('stats')
    statsbox2.SetY1NDC(0.35)
    statsbox2.SetY2NDC(0.54)
    statsbox2.SetX1NDC(0.70)
    statsbox2.SetX2NDC(0.89)
    statsbox2.SetTextColor(ROOT.kRed)
    h['sig'].DrawNormalized('s')
    h['bkg'].DrawNormalized('same s')
    statsbox.Draw()
    statsbox2.Draw()
    c[cname].Modified()
    c[cname].Update()
    leg1.Draw('same')
    tl.Draw()
    #tl2.Draw()

    

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
    
## chIso plots
chIsoNames =  {}
for ir,dr in enumerate(conesDR):
    chIsoNames[dr] = ['relChIso%s_dZ05_simVtx'%dr,
                      'relChIso%s_dZ1_simVtx'%dr,
                      'relChIso%s_dZ2_simVtx'%dr,
                      'relChIso%s_dZ05'%dr,
                      'relChIso%s_dZ1'%dr,
                      'relChIso%s_dZ2'%dr,
                      #'relChIso%s_reldZ'%dr,
                      'relChIso%s_dZmu05'%dr,
                      'relChIso%s_dZmu1'%dr,
                      'relChIso%s_dZmu2'%dr,
                      'relChIso%s_dZmu5'%dr,
                      'relChIso%s_dZmu10'%dr
                     ]


tcuts = ['2s','3s','5s']
#tcuts = ['3s', '5s']    
#tcuts = ['3s']

# load histograms
h_chIso = {}
h_chIso_barrel = {}
h_chIso_endcap = {}
h_chIso_dT = {}
h_chIso_dT_barrel = {}
h_chIso_dT_endcap = {}

for ir,dr in enumerate(conesDR):
    h_chIso[dr] = {}
    h_chIso_barrel[dr] = {}
    h_chIso_endcap[dr] = {}

    h_chIso_dT[dr] = {}
    h_chIso_dT_barrel[dr] = {}
    h_chIso_dT_endcap[dr] = {}
    
    for name in chIsoNames[dr]:
        h_chIso[dr][name] = {}
        h_chIso_barrel[dr][name] = {}
        h_chIso_endcap[dr][name] = {}
        h_chIso_dT[dr][name] = {}
        h_chIso_dT_barrel[dr][name]  = {}
        h_chIso_dT_endcap[dr][name] = {}

        for proc in 'sig','bkg':
            h_chIso[dr][name][proc] = f[proc].Get('h_muon_'+name)
            h_chIso_barrel[dr][name][proc] = f[proc].Get('h_muon_'+name+'_barrel')
            h_chIso_endcap[dr][name][proc] = f[proc].Get('h_muon_'+name+'_endcap')

            h_chIso_dT[dr][name][proc] = {}
            h_chIso_dT_barrel[dr][name][proc] = {} 
            h_chIso_dT_endcap[dr][name][proc] = {}
            
            for tcut in tcuts:  
                h_chIso_dT[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dT%s'%tcut)
                h_chIso_dT_barrel[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dT%s_barrel'%tcut)
                h_chIso_dT_endcap[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dT%s_endcap'%tcut)
                    
                if ('simVtx' in name):
                    h_chIso_dT[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name.replace('_simVtx','')+'_dT%s_simVtx'%tcut)
                    h_chIso_dT_barrel[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name.replace('_simVtx','')+'_dT%s_simVtx_barrel'%tcut)
                    h_chIso_dT_endcap[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name.replace('_simVtx','')+'_dT%s_simVtx_endcap'%tcut)

                if ( 'dZmu' in name ):
                    if ( tcut == '2s' or tcut == '5s' ): continue
                    h_chIso_dT[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dTmu')
                    h_chIso_dT_barrel[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dTmu_barrel')
                    h_chIso_dT_endcap[dr][name][proc][tcut] = f[proc].Get('h_muon_'+name+'_dTmu_endcap')
                    
                    
                    
                

# make efficiency plpots and rocs
h_eff_chIso = {}
h_eff_chIso_barrel = {}
h_eff_chIso_endcap = {}
h_eff_chIso_dT = {}
h_eff_chIso_dT_barrel = {}
h_eff_chIso_dT_endcap = {}

g_roc_chIso = {}
g_roc_chIso_barrel = {}
g_roc_chIso_endcap = {}
g_roc_chIso_dT = {}
g_roc_chIso_dT_barrel = {}
g_roc_chIso_dT_endcap = {}


c_roc = {}
c_roc_barrel = {}
c_roc_endcap = {}


leg2 = ROOT.TLegend(0.15,0.7, 0.5, 0.89)
leg2.SetBorderSize(0)


# efficiency plots
for ir,dr in enumerate(conesDR):
    h_eff_chIso[dr] = {}
    h_eff_chIso_barrel[dr] = {}
    h_eff_chIso_endcap[dr] = {}
    h_eff_chIso_dT[dr] = {}
    h_eff_chIso_dT_barrel[dr] = {}
    h_eff_chIso_dT_endcap[dr] = {}
    
    for name in chIsoNames[dr]:
        h_eff_chIso[dr][name] = {}
        h_eff_chIso_barrel[dr][name] = {}
        h_eff_chIso_endcap[dr][name] = {}
        h_eff_chIso_dT[dr][name] = {}
        h_eff_chIso_dT_barrel[dr][name] = {}
        h_eff_chIso_dT_endcap[dr][name] = {}
    
        for proc in 'sig','bkg':
            
            h_eff_chIso[dr][name][proc] = h_chIso[dr][name][proc].Clone('h_eff_%s_%s'%(name,proc))
            h_eff_chIso_barrel[dr][name][proc] = h_chIso_barrel[dr][name][proc].Clone('h_eff_%s_barrel_%s'%(name,proc))
            h_eff_chIso_endcap[dr][name][proc] = h_chIso_endcap[dr][name][proc].Clone('h_eff_%s_endcap_%s'%(name,proc))

            print 'Making efficiency plots for ', name
            makeEff(h_chIso[dr][name][proc],h_eff_chIso[dr][name][proc])
            makeEff(h_chIso_barrel[dr][name][proc],h_eff_chIso_barrel[dr][name][proc])
            makeEff(h_chIso_endcap[dr][name][proc],h_eff_chIso_endcap[dr][name][proc])
            
            h_eff_chIso_dT[dr][name][proc] = {}
            h_eff_chIso_dT_barrel[dr][name][proc] = {}
            h_eff_chIso_dT_endcap[dr][name][proc] = {}
            
            
            for tcut in tcuts:
                if ('dZmu' in name):
                    if (tcut == '2s' or tcut == '5s'): continue
                    h_eff_chIso_dT[dr][name][proc][tcut] = h_chIso_dT[dr][name][proc][tcut].Clone('h_eff_%s_dTmu_%s'%(name,proc))
                    h_eff_chIso_dT_barrel[dr][name][proc][tcut] = h_chIso_dT_barrel[dr][name][proc][tcut].Clone('h_eff_%s_dTmu_barrel_%s'%(name,proc))
                    h_eff_chIso_dT_endcap[dr][name][proc][tcut] = h_chIso_dT_endcap[dr][name][proc][tcut].Clone('h_eff_%s_dTmu_endcap_%s'%(name,proc))
                else:
                    h_eff_chIso_dT[dr][name][proc][tcut] = h_chIso_dT[dr][name][proc][tcut].Clone('h_eff_%s_dT%s_%s'%(name,tcut,proc))
                    h_eff_chIso_dT_barrel[dr][name][proc][tcut] = h_chIso_dT_barrel[dr][name][proc][tcut].Clone('h_eff_%s_dT%s_barrel_%s'%(name,tcut,proc))
                    h_eff_chIso_dT_endcap[dr][name][proc][tcut] = h_chIso_dT_endcap[dr][name][proc][tcut].Clone('h_eff_%s_dT%s_endcap_%s'%(name,tcut,proc))
                    
                print 'Making efficiency plots for ', name, tcut
                makeEff(h_chIso_dT[dr][name][proc][tcut],h_eff_chIso_dT[dr][name][proc][tcut])
                makeEff(h_chIso_dT_barrel[dr][name][proc][tcut],h_eff_chIso_dT_barrel[dr][name][proc][tcut])
                makeEff(h_chIso_dT_endcap[dr][name][proc][tcut],h_eff_chIso_dT_endcap[dr][name][proc][tcut])

            
# rocs
for ir,dr in enumerate(conesDR):
    g_roc_chIso[dr] = {}
    g_roc_chIso_barrel[dr] = {}
    g_roc_chIso_endcap[dr] = {}
    g_roc_chIso_dT[dr] = {}
    g_roc_chIso_dT_barrel[dr] = {}
    g_roc_chIso_dT_endcap[dr] = {}

    c_roc[dr] = {}
    c_roc_barrel[dr] = {}
    c_roc_endcap[dr] = {}

    for iname, name in enumerate(chIsoNames[dr]):
        g_roc_chIso[dr][name] = ROOT.TGraphAsymmErrors()
        g_roc_chIso[dr][name].SetName('g_roc_%s'%name)

        g_roc_chIso_barrel[dr][name] = ROOT.TGraphAsymmErrors()
        g_roc_chIso_barrel[dr][name].SetName('g_roc_%s_barrel'%name)
        
        g_roc_chIso_endcap[dr][name] = ROOT.TGraphAsymmErrors()
        g_roc_chIso_endcap[dr][name].SetName('g_roc_%s_endcap'%name)

        # make ROCs
        print 'Making ROCs for ', name
        makeRoc(h_chIso[dr][name]['sig'], h_chIso[dr][name]['bkg'], g_roc_chIso[dr][name])
        makeRoc(h_chIso_barrel[dr][name]['sig'], h_chIso_barrel[dr][name]['bkg'], g_roc_chIso_barrel[dr][name])
        makeRoc(h_chIso_endcap[dr][name]['sig'], h_chIso_endcap[dr][name]['bkg'], g_roc_chIso_endcap[dr][name])

        
        c_roc[dr][name] = {}
        c_roc_barrel[dr][name] = {}
        c_roc_endcap[dr][name] = {}

        g_roc_chIso_dT[dr][name] = {}
        g_roc_chIso_dT_barrel[dr][name] = {}
        g_roc_chIso_dT_endcap[dr][name] = {}
    
        for icut, tcut in enumerate(tcuts):
            g_roc_chIso_dT[dr][name][tcut] = ROOT.TGraphAsymmErrors()
            g_roc_chIso_dT[dr][name][tcut].SetName('g_roc_%s_dT%s'%(name,tcut))
        
            g_roc_chIso_dT_barrel[dr][name][tcut] = ROOT.TGraphAsymmErrors()
            g_roc_chIso_dT_barrel[dr][name][tcut].SetName('g_roc_%s_dT%s_barrel'%(name,tcut))
            
            g_roc_chIso_dT_endcap[dr][name][tcut] = ROOT.TGraphAsymmErrors()
            g_roc_chIso_dT_endcap[dr][name][tcut].SetName('g_roc_%s_dT%s_endcap'%(name,tcut))
            
            if ('dZmu' in name):
                if (tcut == '2s' or tcut == '5s'): continue
                g_roc_chIso_dT[dr][name][tcut].SetName('g_roc_%s_dTmu'%name)
                g_roc_chIso_dT_barrel[dr][name][tcut].SetName('g_roc_%s_dTmu_barrel'%name)
                g_roc_chIso_dT_endcap[dr][name][tcut].SetName('g_roc_%s_dTmu_endcap'%name)

            print 'Making ROCs for ', name, tcut
            makeRoc(h_chIso_dT[dr][name]['sig'][tcut], h_chIso_dT[dr][name]['bkg'][tcut], g_roc_chIso_dT[dr][name][tcut])
            makeRoc(h_chIso_dT_barrel[dr][name]['sig'][tcut], h_chIso_dT_barrel[dr][name]['bkg'][tcut], g_roc_chIso_dT_barrel[dr][name][tcut])
            makeRoc(h_chIso_dT_endcap[dr][name]['sig'][tcut], h_chIso_dT_endcap[dr][name]['bkg'][tcut], g_roc_chIso_dT_endcap[dr][name][tcut])

            
            #draw ROCs
            c_roc[dr][name][tcut] = ROOT.TCanvas('roc_%s_dT%s'%(name,tcut),'roc_%s_dT%s'%(name,tcut),500,500)
            c_roc[dr][name][tcut].SetGridx()
            c_roc[dr][name][tcut].SetGridy()
            c_roc[dr][name][tcut].SetTickx()
            c_roc[dr][name][tcut].SetTicky()
            g_roc_chIso[dr][name].GetXaxis().SetTitle('#epsilon_{prompt}')
            g_roc_chIso[dr][name].GetYaxis().SetTitle('#epsilon_{non-prompt}')
            g_roc_chIso[dr][name].GetXaxis().SetRangeUser(0.80,1.01)
            g_roc_chIso[dr][name].GetYaxis().SetRangeUser(0.0,0.1)
            g_roc_chIso[dr][name].SetMarkerColor(ROOT.kBlue)
            g_roc_chIso[dr][name].SetLineColor(ROOT.kBlue)
            g_roc_chIso[dr][name].SetLineWidth(2)
            g_roc_chIso[dr][name].SetFillColorAlpha(ROOT.kBlue, 0.35)
            g_roc_chIso[dr][name].SetFillStyle(1001)
            g_roc_chIso[dr][name].Draw('A L E3')
            g_roc_chIso_dT[dr][name][tcut].SetMarkerColor(ROOT.kRed)
            g_roc_chIso_dT[dr][name][tcut].SetLineColor(ROOT.kRed)
            g_roc_chIso_dT[dr][name][tcut].SetLineWidth(2)
            g_roc_chIso_dT[dr][name][tcut].SetFillColorAlpha(ROOT.kRed,0.35)
            g_roc_chIso_dT[dr][name][tcut].SetFillStyle(1001)
            g_roc_chIso_dT[dr][name][tcut].Draw('L E3 same')

            if (ir==0 and iname==0 and icut==0):
                leg2.AddEntry(g_roc_chIso[dr][name],'no MTD', 'L')
                leg2.AddEntry(g_roc_chIso_dT[dr][name][tcut],'MTD, #sigma_{t} = %d ps'%resolution, 'L')
            
            leg2.Draw('same')
            tl.Draw()
            tl2.Draw()
            #raw_input('ok?')



            c_roc_barrel[dr][name][tcut] = ROOT.TCanvas('roc_%s_dT%s_barrel'%(name,tcut),'roc_%s_dT%s_barrel'%(name,tcut),500,500)
            c_roc_barrel[dr][name][tcut].SetGridx()
            c_roc_barrel[dr][name][tcut].SetGridy()
            c_roc_barrel[dr][name][tcut].SetTickx()
            c_roc_barrel[dr][name][tcut].SetTicky()
            g_roc_chIso_barrel[dr][name].GetXaxis().SetTitle('#epsilon_{prompt}')
            g_roc_chIso_barrel[dr][name].GetYaxis().SetTitle('#epsilon_{non-prompt}')
            g_roc_chIso_barrel[dr][name].GetXaxis().SetRangeUser(0.80,1.01)
            g_roc_chIso_barrel[dr][name].GetYaxis().SetRangeUser(0.0,0.1) 
            g_roc_chIso_barrel[dr][name].SetMarkerColor(ROOT.kBlue)
            g_roc_chIso_barrel[dr][name].SetLineColor(ROOT.kBlue)
            g_roc_chIso_barrel[dr][name].SetLineWidth(2)
            g_roc_chIso_barrel[dr][name].SetFillColorAlpha(ROOT.kBlue,0.35)
            g_roc_chIso_barrel[dr][name].SetFillStyle(1001)
            g_roc_chIso_barrel[dr][name].Draw('A L E3')
            g_roc_chIso_dT_barrel[dr][name][tcut].SetMarkerColor(ROOT.kRed)
            g_roc_chIso_dT_barrel[dr][name][tcut].SetLineColor(ROOT.kRed)
            g_roc_chIso_dT_barrel[dr][name][tcut].SetLineWidth(2)
            g_roc_chIso_dT_barrel[dr][name][tcut].SetFillColorAlpha(ROOT.kRed,0.35)
            g_roc_chIso_dT_barrel[dr][name][tcut].SetFillStyle(1001)
            g_roc_chIso_dT_barrel[dr][name][tcut].Draw('L E3 same')
            leg2.Draw('same')
            tl.Draw()
            tl2.Draw()
            #raw_input('ok?')

            
            c_roc_endcap[dr][name][tcut] = ROOT.TCanvas('roc_%s_dT%s_endcap'%(name,tcut),'roc_%s_dT%s_endcap'%(name,tcut),500,500)
            c_roc_endcap[dr][name][tcut].SetGridx()
            c_roc_endcap[dr][name][tcut].SetGridy()
            c_roc_endcap[dr][name][tcut].SetTickx()
            c_roc_endcap[dr][name][tcut].SetTicky()
            g_roc_chIso_endcap[dr][name].GetXaxis().SetTitle('#epsilon_{prompt}')
            g_roc_chIso_endcap[dr][name].GetYaxis().SetTitle('#epsilon_{non-prompt}')
            g_roc_chIso_endcap[dr][name].GetXaxis().SetRangeUser(0.80,1.01)
            g_roc_chIso_endcap[dr][name].GetYaxis().SetRangeUser(0.0,0.1)
            g_roc_chIso_endcap[dr][name].SetMarkerColor(ROOT.kBlue)
            g_roc_chIso_endcap[dr][name].SetLineColor(ROOT.kBlue)
            g_roc_chIso_endcap[dr][name].SetLineWidth(2)
            g_roc_chIso_endcap[dr][name].SetFillColorAlpha(ROOT.kBlue,0.35)
            g_roc_chIso_endcap[dr][name].SetFillStyle(1001)
            g_roc_chIso_endcap[dr][name].Draw('A L E3')
            g_roc_chIso_dT_endcap[dr][name][tcut].SetMarkerColor(ROOT.kRed)
            g_roc_chIso_dT_endcap[dr][name][tcut].SetLineColor(ROOT.kRed)
            g_roc_chIso_dT_endcap[dr][name][tcut].SetLineWidth(2)
            g_roc_chIso_dT_endcap[dr][name][tcut].SetFillColorAlpha(ROOT.kRed,0.35)
            g_roc_chIso_dT_endcap[dr][name][tcut].SetFillStyle(1001)
            g_roc_chIso_dT_endcap[dr][name][tcut].Draw('L E3 same')
            leg2.Draw('same')
            tl.Draw()
            tl2.Draw()
            #raw_input('ok?')

            

c_chIso = {}
c_chIso_barrel = {}
c_chIso_endcap = {}

c_eff = {}
c_eff_barrel= {}
c_eff_endcap= {}

color = {'sig':ROOT.kBlue, 'bkg':ROOT.kRed}

leg3 = ROOT.TLegend(0.15, 0.7, 0.65, 0.89)
leg3.SetBorderSize(0)

nre = 1
tcut = '3s'
            
for ir,dr in enumerate(conesDR):
    c_chIso[dr] = {}
    c_chIso_barrel[dr] = {}
    c_chIso_endcap[dr] = {}
    
    c_eff[dr] = {}
    c_eff_barrel[dr] = {}
    c_eff_endcap[dr] = {}

    for iname,name in enumerate(chIsoNames[dr]):
        
        for proc in 'sig','bkg':
            
            h_chIso[dr][name][proc].Rebin(nre)
            h_chIso_barrel[dr][name][proc].Rebin(nre)
            h_chIso_endcap[dr][name][proc].Rebin(nre)

            h_chIso[dr][name][proc].SetLineColor(color[proc])
            h_chIso_barrel[dr][name][proc].SetLineColor(color[proc])
            h_chIso_endcap[dr][name][proc].SetLineColor(color[proc])

            h_eff_chIso[dr][name][proc].SetLineColor(color[proc])
            h_eff_chIso_barrel[dr][name][proc].SetLineColor(color[proc])
            h_eff_chIso_endcap[dr][name][proc].SetLineColor(color[proc])

            h_chIso_dT[dr][name][proc][tcut].Rebin(nre)
            h_chIso_dT_barrel[dr][name][proc][tcut].Rebin(nre)
            h_chIso_dT_endcap[dr][name][proc][tcut].Rebin(nre)
        
            h_chIso_dT[dr][name][proc][tcut].SetLineColor(color[proc])
            h_chIso_dT_barrel[dr][name][proc][tcut].SetLineColor(color[proc])
            h_chIso_dT_endcap[dr][name][proc][tcut].SetLineColor(color[proc])
            
            h_chIso_dT[dr][name][proc][tcut].SetLineStyle(2)
            h_chIso_dT_barrel[dr][name][proc][tcut].SetLineStyle(2)
            h_chIso_dT_endcap[dr][name][proc][tcut].SetLineStyle(2)
            
            h_eff_chIso_dT[dr][name][proc][tcut].SetLineColor(color[proc])
            h_eff_chIso_dT_barrel[dr][name][proc][tcut].SetLineColor(color[proc])
            h_eff_chIso_dT_endcap[dr][name][proc][tcut].SetLineColor(color[proc])
            
            h_eff_chIso_dT[dr][name][proc][tcut].SetLineStyle(2)
            h_eff_chIso_dT_barrel[dr][name][proc][tcut].SetLineStyle(2)
            h_eff_chIso_dT_endcap[dr][name][proc][tcut].SetLineStyle(2)


        #BTL+ETL        
        c_chIso[dr][name] = ROOT.TCanvas('muon_%s'%name,'muon_%s'%name,500,500)
        print  c_chIso[dr][name].GetName()
        c_chIso[dr][name].SetLogy()
        h_chIso[dr][name]['sig'].GetXaxis().SetRangeUser(0,1.2)    
        h_chIso[dr][name]['sig'].DrawNormalized()
        h_chIso_dT[dr][name]['sig'][tcut].DrawNormalized('same')
        h_chIso[dr][name]['bkg'].DrawNormalized('histo same')
        h_chIso_dT[dr][name]['bkg'][tcut].DrawNormalized('histo same')
        if (ir==0 and iname==0):
            leg3.AddEntry(h_chIso[dr][name]['sig'],'prompt muons (Z#rightarrow#mu#mu), no MTD', 'L')
            leg3.AddEntry(h_chIso_dT[dr][name]['sig'][tcut],'prompt muons (Z#rightarrow#mu#mu), with MTD (#sigma_{t} = %d ps)'%resolution, 'L')
            leg3.AddEntry(h_chIso[dr][name]['bkg'],'non-prompt muons (%s), no MTD'%textProc, 'L')
            leg3.AddEntry(h_chIso_dT[dr][name]['bkg'][tcut],'non-prompt muons (%s), with MTD (#sigma_{t} = %d ps)'%(textProc,resolution), 'L')
        leg3.Draw('same')
        tl.Draw()
        #tl2.Draw()

        #BTL only
        c_chIso_barrel[dr][name] = ROOT.TCanvas('muon_%s_barrel'%name,'muon_%s_barrel'%name,500,500)
        c_chIso_barrel[dr][name].SetLogy()
        h_chIso_barrel[dr][name]['sig'].GetXaxis().SetRangeUser(0,1.1)    
        h_chIso_barrel[dr][name]['sig'].DrawNormalized()
        h_chIso_dT_barrel[dr][name]['sig'][tcut].DrawNormalized('same')
        h_chIso_barrel[dr][name]['bkg'].DrawNormalized('histo same')
        h_chIso_dT_barrel[dr][name]['bkg'][tcut].DrawNormalized('histo same')
        leg3.Draw('same')
        tl.Draw()
        #tl2.Draw()
        
        #ETL only
        c_chIso_endcap[dr][name] = ROOT.TCanvas('muon_%s_endcap'%name,'muon_%s_endcap'%name,500,500)
        c_chIso_endcap[dr][name].SetLogy()
        h_chIso_endcap[dr][name]['sig'].GetXaxis().SetRangeUser(0,1.1)    
        h_chIso_endcap[dr][name]['sig'].DrawNormalized()
        h_chIso_dT_endcap[dr][name]['sig'][tcut].DrawNormalized('same')
        h_chIso_endcap[dr][name]['bkg'].DrawNormalized('histo same')
        h_chIso_dT_endcap[dr][name]['bkg'][tcut].DrawNormalized('histo same')
        leg3.Draw('same')
        tl.Draw()
        #tl2.Draw()

        c_eff[dr][name] = ROOT.TCanvas('efficiency_%s'%name,'eff_%s'%name,500,500)
        c_eff[dr][name].SetLogy()
        c_eff[dr][name].SetGridx()
        c_eff[dr][name].SetGridy()
        h_eff_chIso[dr][name]['sig'].GetXaxis().SetRangeUser(0,0.1)
        h_eff_chIso[dr][name]['sig'].GetYaxis().SetRangeUser(0.01,10)
        h_eff_chIso[dr][name]['sig'].GetXaxis().SetTitle('relative charged isolation')
        h_eff_chIso[dr][name]['sig'].GetYaxis().SetTitle('efficiency')
        h_eff_chIso[dr][name]['sig'].Draw()
        h_eff_chIso[dr][name]['bkg'].Draw('h same')
        h_eff_chIso_dT[dr][name]['sig'][tcut].Draw('h same')
        h_eff_chIso_dT[dr][name]['bkg'][tcut].Draw('h same')
        leg3.Draw('same')
        tl.Draw()
        #tl2.Draw()

        
    
#save plots
dirname = '%dps_%s_%s'%(resolution, pu, bkgProc)
if (suffix != ''):
    dirname = dirname+'_'+suffix
dirname = dirname+'/'
os.system('mkdir %s'%dirname)
shutil.copy('index.php', dirname)

for cname in c:
    c[cname].SaveAs(dirname+c[cname].GetName()+'.png')
    c[cname].SaveAs(dirname+c[cname].GetName()+'.pdf')


for dr in conesDR:
    for name in chIsoNames[dr]:
        for typ in '.png','.pdf':
            c_chIso[dr][name].SaveAs(dirname+c_chIso[dr][name].GetName()+typ)
            c_chIso_barrel[dr][name].SaveAs(dirname+c_chIso_barrel[dr][name].GetName()+typ)
            c_chIso_endcap[dr][name].SaveAs(dirname+c_chIso_endcap[dr][name].GetName()+typ)

            c_eff[dr][name].SaveAs(dirname+c_eff[dr][name].GetName()+typ)

            for tcut in tcuts:
                if ('dZmu' in name and (tcut != '3s') ): continue
                c_roc[dr][name][tcut].SaveAs(dirname+c_roc[dr][name][tcut].GetName()+typ)
                c_roc_barrel[dr][name][tcut].SaveAs(dirname+c_roc_barrel[dr][name][tcut].GetName()+typ)
                c_roc_endcap[dr][name][tcut].SaveAs(dirname+c_roc_endcap[dr][name][tcut].GetName()+typ)

            
fout = ROOT.TFile(dirname+'roc_%dps_%s.root'%(resolution, pu),'recreate')
for dr in conesDR:
    for name in chIsoNames[dr]:
        g_roc_chIso[dr][name].Write()
        g_roc_chIso_barrel[dr][name].Write()
        g_roc_chIso_endcap[dr][name].Write()
        for tcut in tcuts:
            if ('dZmu' in name and (tcut != '3s') ): continue
            g_roc_chIso_dT[dr][name][tcut].Write()
            g_roc_chIso_dT_barrel[dr][name][tcut].Write()
            g_roc_chIso_dT_endcap[dr][name][tcut].Write()

      

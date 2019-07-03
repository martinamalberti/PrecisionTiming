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
        #graph.SetPointError(ibin-1, effSigErrDown, effSigErrUp, effBkgErrDown, effBkgErrUp)
        graph.SetPointError(ibin-1, effSigErrDown, effSigErrUp)    
        
    return


# ==== MAIN  =====

release = '10_4_0_mtd5'
pu =  sys.argv[1]
bkgProc = sys.argv[2]
suffix = sys.argv[3]

fSig = '../'+release+'/testTracksMva_DYToLL_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy.root'
fBkg = '../'+release+'/testTracksMva_%s_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt_noDxy.root'%bkgProc
    
    
plotdir = release+'/optimizationTrackPuidMva_%s_%s'%(pu,bkgProc)
plotdir = plotdir+'_'+suffix
    
os.system('mkdir %s'%plotdir)
shutil.copy('index.php', plotdir)
    
tl = ROOT.TLatex( 0.70, 0.84,'<PU> = 200')
if ('noPU' in fSig):
    tl = ROOT.TLatex( 0.70, 0.84,'<PU> = 0')
tl.SetNDC()
tl.SetTextSize(0.035)

tl2 = {}
tl2['barrel'] = ROOT.TLatex( 0.69, 0.78,'barrel muons')
tl2['barrel'] .SetNDC()
tl2['barrel'] .SetTextSize(0.030)

tl2['endcap'] = ROOT.TLatex( 0.68, 0.78,'endcap muons')
tl2['endcap'] .SetNDC()
tl2['endcap'] .SetTextSize(0.030)

tl3 = {}
tl3['prompt'] = ROOT.TLatex( 0.65, 0.78,'Z#rightarrow#mu#mu')
tl3['prompt'] .SetNDC()
tl3['prompt'] .SetTextSize(0.030)
tl3['fake'] = ROOT.TLatex( 0.65, 0.78,'t#bar{t}')
tl3['fake'] .SetNDC()
tl3['fake'] .SetTextSize(0.030)


f = {}
f['prompt'] = ROOT.TFile.Open(fSig)
f['fake'] = ROOT.TFile.Open(fBkg)

h_relChIso03_mva3d = {}
h_relChIso03_mva4d = {}

g_roc_relChIso03_mva3d = {}
g_roc_relChIso03_mva4d = {}



cuts3d = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
cuts4d = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
                
col = []
for i in range(0, len(cuts3d)):
    if (i < 4 ): col.append(ROOT.kRed+i)
    if (i < 8 ): col.append(ROOT.kOrange+i)
    if (i < 12 ): col.append(ROOT.kGreen+i)
    if (i < 16 ): col.append(ROOT.kCyan+i)
    if (i < 20 ): col.append(ROOT.kBlue+i)
    if (i < 24 ): col.append(ROOT.kGray+i)
    

nRe = 4
    
for id, d in enumerate(['barrel','endcap']):
    h_relChIso03_mva3d[d] = {}
    h_relChIso03_mva4d[d] = {}
    g_roc_relChIso03_mva3d[d] = {}
    g_roc_relChIso03_mva4d[d] = {}

    # mva 3d
    for cut in cuts3d:
        h_relChIso03_mva3d[d][cut] = {}
        g_roc_relChIso03_mva3d[d][cut] = ROOT.TGraphErrors()
        g_roc_relChIso03_mva3d[d][cut].SetName('g_roc_relChIso03_mva3d_%f_%s'%(cut,d)) 
        for proc in ['prompt','fake']:
            h_relChIso03_mva3d[d][cut][proc] = f[proc].Get('h_muon_relChIso03_mva3D_%.03f_%s'%(cut,d))
            h_relChIso03_mva3d[d][cut][proc].Rebin(nRe)
        makeRoc(h_relChIso03_mva3d[d][cut]['prompt'], h_relChIso03_mva3d[d][cut]['fake'], g_roc_relChIso03_mva3d[d][cut] )

    # mva 4d
    for cut in cuts4d:
        h_relChIso03_mva4d[d][cut] = {}
        g_roc_relChIso03_mva4d[d][cut] = ROOT.TGraphErrors()
        g_roc_relChIso03_mva4d[d][cut].SetName('g_roc_relChIso03_mva4d_%f_%s'%(cut,d)) 
        for proc in ['prompt','fake']:
            h_relChIso03_mva4d[d][cut][proc] = f[proc].Get('h_muon_relChIso03_mva4D_%.03f_%s'%(cut,d))
            h_relChIso03_mva4d[d][cut][proc].Rebin(nRe)
        makeRoc(h_relChIso03_mva4d[d][cut]['prompt'], h_relChIso03_mva4d[d][cut]['fake'], g_roc_relChIso03_mva4d[d][cut] )



# now plot
c_roc_relChIso03_mva3d = {}
c_roc_relChIso03_mva4d = {}
leg3 = {}
leg4 = {}

for id, d in enumerate(['barrel','endcap']):
    c_roc_relChIso03_mva3d[d] =  ROOT.TCanvas('roc_relChIso03_mva3d_scan_%s'%d, 'roc_relChIso03_mva3d_scan_%s'%d)
    c_roc_relChIso03_mva3d[d].SetGridx()
    c_roc_relChIso03_mva3d[d].SetGridy()
    leg3[d] = ROOT.TLegend(0.15, 0.7, 0.55, 0.92)
    leg3[d].SetBorderSize(0)

    for icut,cut in enumerate(cuts3d):
        g_roc_relChIso03_mva3d[d][cut].GetXaxis().SetRangeUser(0.80,1.01)
        g_roc_relChIso03_mva3d[d][cut].GetYaxis().SetRangeUser(0.0,0.035)
        g_roc_relChIso03_mva3d[d][cut].GetXaxis().SetTitle("Prompt efficiency")
        g_roc_relChIso03_mva3d[d][cut].GetYaxis().SetTitle("Non-prompt efficiency")
        g_roc_relChIso03_mva3d[d][cut].SetLineColor(col[icut])
        g_roc_relChIso03_mva3d[d][cut].SetLineWidth(3)
        g_roc_relChIso03_mva3d[d][cut].SetFillColorAlpha(col[icut],0.0)
        g_roc_relChIso03_mva3d[d][cut].SetFillStyle(1001)
        if (icut == 0 ):
            g_roc_relChIso03_mva3d[d][cut].Draw('a l e3')
        else:
            g_roc_relChIso03_mva3d[d][cut].Draw('l e3 same')
            
        leg3[d].AddEntry( g_roc_relChIso03_mva3d[d][cut],'no MTD, 3DMVA > %f'%cut,'L')

    tl.Draw()
    tl2[d].Draw()
    leg3[d].Draw()
    CMS_lumi.CMS_lumi(c_roc_relChIso03_mva3d[d], iPeriod, iPos)



    c_roc_relChIso03_mva4d[d] =  ROOT.TCanvas('roc_relChIso03_mva4d_scan_%s'%d, 'roc_relChIso03_mva4d_scan_%s'%d)
    c_roc_relChIso03_mva4d[d].SetGridx()
    c_roc_relChIso03_mva4d[d].SetGridy()
    leg4[d] = ROOT.TLegend(0.15, 0.7, 0.55, 0.92)
    leg4[d].SetBorderSize(0)

    for icut,cut in enumerate(cuts4d):
        g_roc_relChIso03_mva4d[d][cut].GetXaxis().SetRangeUser(0.80,1.01)
        g_roc_relChIso03_mva4d[d][cut].GetYaxis().SetRangeUser(0.0,0.035)
        g_roc_relChIso03_mva4d[d][cut].GetXaxis().SetTitle("Prompt efficiency")
        g_roc_relChIso03_mva4d[d][cut].GetYaxis().SetTitle("Non-prompt efficiency")
        g_roc_relChIso03_mva4d[d][cut].SetLineColor(col[icut])
        g_roc_relChIso03_mva4d[d][cut].SetLineWidth(3)
        g_roc_relChIso03_mva4d[d][cut].SetFillColorAlpha(col[icut],0.0)
        g_roc_relChIso03_mva4d[d][cut].SetFillStyle(1001)
        if (icut == 0 ):
            g_roc_relChIso03_mva4d[d][cut].Draw('a l e3')
        else:
            g_roc_relChIso03_mva4d[d][cut].Draw('l e3 same')
            
        leg4[d].AddEntry( g_roc_relChIso03_mva4d[d][cut],'no MTD, 4DMVA > %.03f'%cut,'L')

    tl.Draw()
    tl2[d].Draw()
    leg4[d].Draw()
    CMS_lumi.CMS_lumi(c_roc_relChIso03_mva4d[d], iPeriod, iPos)


    
    raw_input('ok?')

# save canvases
for id, d in enumerate(['barrel','endcap']):
    c_roc_relChIso03_mva3d[d].SaveAs(plotdir+'/'+c_roc_relChIso03_mva3d[d].GetName()+'.png')
    c_roc_relChIso03_mva4d[d].SaveAs(plotdir+'/'+c_roc_relChIso03_mva4d[d].GetName()+'.png')
    c_roc_relChIso03_mva3d[d].SaveAs(plotdir+'/'+c_roc_relChIso03_mva3d[d].GetName()+'.pdf')
    c_roc_relChIso03_mva4d[d].SaveAs(plotdir+'/'+c_roc_relChIso03_mva4d[d].GetName()+'.pdf')

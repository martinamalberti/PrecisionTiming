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

pu = sys.argv[1]
bkgProc = sys.argv[2]


fSig = '../testTracks_DYToLL_PU200.root'
fBkg = '../testTracks_TTbar_PU200.root'
#fBkg = '../testTracks_QCD_PU200.root'

if (pu == 'noPU'):
    fSig = '../testTracks_DYToLL_noPU.root'
    fBkg = '../testTracks_TTbar_noPU.root'
    #fBkg = '../testTracks_QCD_noPU.root'

#fSig = '../testTracks_DYToLL_PU200_simVtxCut.root'
#fBkg = '../testTracks_TTbar_PU200_simVtxCut.root'
#fBkg = '../testTracks_QCD_PU200_simVtxCut.root'

plotdir = 'tracksInCone_%s_%s'%(pu,bkgProc)
if ('simVtx' in fBkg):
    plotdir = plotdir+'_SimVtxCut'
    
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


f = {}
f['prompt'] = ROOT.TFile.Open(fSig)
f['fake'] = ROOT.TFile.Open(fBkg)

h_tracks_dt_vtx = {}
h_tracks_dz_vtx = {}

h_tracks_dz_mu = {}
h_tracks_dt_mu = {}

h_tracks_pt_notiming = {}
h_tracks_pt_timing_PU = {}
h_tracks_pt_timing_PV = {}

h_tracks_eta_notiming = {}
h_tracks_eta_timing_PU = {}
h_tracks_eta_timing_PV = {}

h_tracks_dr_notiming = {}
h_tracks_dr_timing_PU = {}
h_tracks_dr_timing_PV = {}

h2_tracks_dzvtx_dtvtx = {}

for d in ['barrel','endcap']:
    h_tracks_dt_vtx[d] = {}
    h_tracks_dz_vtx[d] = {}
    h_tracks_dz_mu[d] = {}
    h_tracks_dt_mu[d] = {}
    h_tracks_pt_notiming[d] = {}
    h_tracks_pt_timing_PU[d] = {}
    h_tracks_pt_timing_PV[d] = {}
    
    h_tracks_eta_notiming[d] = {}
    h_tracks_eta_timing_PU[d] = {}
    h_tracks_eta_timing_PV[d] = {}
    
    h_tracks_dr_notiming[d] = {}
    h_tracks_dr_timing_PU[d] = {}
    h_tracks_dr_timing_PV[d] = {}

    h2_tracks_dzvtx_dtvtx[d] = {}
    
    for proc in 'prompt', 'fake':
        h_tracks_dz_vtx[d][proc]        =  f[proc].Get('h_tracks_dz_vtx_'+d)
        h_tracks_dt_vtx[d][proc]        =  f[proc].Get('h_tracks_dt_vtx_'+d)
        h_tracks_dz_mu[d][proc]         =  f[proc].Get('h_tracks_dz_mu_'+d)
        h_tracks_dt_mu[d][proc]         =  f[proc].Get('h_tracks_dt_mu_'+d)
        
        h_tracks_pt_notiming[d][proc]   = f[proc].Get('h_tracks_pt_notiming_'+d)
        h_tracks_pt_timing_PV[d][proc]  = f[proc].Get('h_tracks_pt_timing_PV_'+d)
        h_tracks_pt_timing_PU[d][proc]  = f[proc].Get('h_tracks_pt_timing_PU_'+d)
        
        h_tracks_eta_notiming[d][proc]  = f[proc].Get('h_tracks_eta_notiming_'+d)
        h_tracks_eta_timing_PV[d][proc] = f[proc].Get('h_tracks_eta_timing_PV_'+d)
        h_tracks_eta_timing_PU[d][proc] = f[proc].Get('h_tracks_eta_timing_PU_'+d)
    
        h_tracks_dr_notiming[d][proc] = f[proc].Get('h_tracks_dr_notiming_'+d)
        h_tracks_dr_timing_PV[d][proc] = f[proc].Get('h_tracks_dr_timing_PV_'+d)
        h_tracks_dr_timing_PU[d][proc] = f[proc].Get('h_tracks_dr_timing_PU_'+d)

        h2_tracks_dzvtx_dtvtx[d][proc] = f[proc].Get('h2_tracks_dzvtx_dtvtx_'+d)

leg1 = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
leg1.SetBorderSize(0)


c = {}

for id,d in enumerate(['barrel', 'endcap']):
    for ih,h in  enumerate([h_tracks_dz_vtx[d], h_tracks_dz_mu[d], h_tracks_pt_notiming[d], h_tracks_pt_timing_PV[d], h_tracks_pt_timing_PU[d], h_tracks_eta_notiming[d], h_tracks_eta_timing_PV[d], h_tracks_eta_timing_PU[d], h_tracks_dr_notiming[d], h_tracks_dr_timing_PV[d], h_tracks_dr_timing_PU[d] ]):
        cname = h['prompt'].GetName().replace('h_','')
        c[cname] = ROOT.TCanvas(cname,cname,500,500)
        h['prompt'].SetLineColor(ROOT.kBlue)
        h['fake'].SetLineColor(ROOT.kRed)
        h['prompt'].GetYaxis().SetRangeUser(0,h['prompt'].GetMaximum()*1.5)
        h['prompt'].DrawNormalized()
        h['fake'].DrawNormalized('histo same')
        if (id == 0 and ih == 0):
            leg1.AddEntry(h['prompt'],'prompt muons (Z#rightarrow#mu#mu)','L')
            leg1.AddEntry(h['fake'],'non-prompt muons (t#bar{t})','L')
            leg1.Draw('same')
        else:
            leg1.Draw('same')
        tl.Draw()
        tl2[d].Draw()
        #raw_input('ok?')


        
for id,d in enumerate(['barrel', 'endcap']):
    for proc in ['prompt','fake']:
        cname = h2_tracks_dzvtx_dtvtx[d][proc].GetName().replace('h2_','')+'_'+proc
        c[cname] = ROOT.TCanvas(cname,cname,500,500)
        h2_tracks_dzvtx_dtvtx[d][proc].GetXaxis().SetRangeUser(-0.1,0.1)
        h2_tracks_dzvtx_dtvtx[d][proc].GetYaxis().SetRangeUser(-0.2,0.2)
        h2_tracks_dzvtx_dtvtx[d][proc].GetXaxis().SetTitle('z_{track} - z_{vtx4D}')
        h2_tracks_dzvtx_dtvtx[d][proc].GetYaxis().SetTitle('t_{track} - t_{vtx4D}')
        h2_tracks_dzvtx_dtvtx[d][proc].Draw('colz')
        #raw_input('ok?')
 
cc  = {}
fitfun2 = {}
fpv = {}
fpu = {}

tpu = {}
tpv = {}


for d in 'barrel','endcap':
    cc[d] = {}
    fitfun2[d] = {}
    fpv[d] = {}
    fpu[d] = {}
    tpu[d] = {}
    tpv[d] = {}
    for proc in ['prompt','fake']:
        cname = h_tracks_dt_mu[d]['prompt'].GetName().replace('h_','') + '_' + proc
        cc[d][proc]= ROOT.TCanvas(cname,cname,500,500)
        h_tracks_dt_mu[d][proc].GetYaxis().SetRangeUser(0.0000001,h_tracks_dt_mu[d][proc].GetMaximum()*1.5)
        h_tracks_dt_mu[d][proc].GetXaxis().SetTitle('t_{track} - t_{#mu} (ns)')
        fitfun2[d][proc] = h_tracks_dt_mu[d][proc].GetFunction('fitfun2_%s'%d)
        fitfun2[d][proc].SetLineColor(ROOT.kGreen+1)
        h_tracks_dt_mu[d][proc].Fit( fitfun2[d][proc],'QRMS')
        h_tracks_dt_mu[d][proc].Draw()
        
        fpv[d][proc] = ROOT.TF1('fpv_%s_%s'%(proc,d),'gaus(0)',-10,10)
        fpv[d][proc].SetParameters(fitfun2[d][proc].GetParameter(0), fitfun2[d][proc].GetParameter(1), fitfun2[d][proc].GetParameter(2))
        
        fpu[d][proc] = ROOT.TF1('fpu_%s_%s'%(proc,d),'gaus(0)',-10,10)
        fpu[d][proc].SetParameters(fitfun2[d][proc].GetParameter(3), fitfun2[d][proc].GetParameter(4), fitfun2[d][proc].GetParameter(5))
        fpu[d][proc].SetLineColor(2)
        fpu[d][proc].Draw('same')

        sigmat_pv = fitfun2[d][proc].GetParameter(2)*1000.
        sigmat_pu = fitfun2[d][proc].GetParameter(5)*1000.

        sigmat_pv_err = fitfun2[d][proc].GetParError(2)*1000.
        sigmat_pu_err = fitfun2[d][proc].GetParError(5)*1000.
        
        ntracks_pv = fpv[d][proc].Integral(-10,10)/h_tracks_dt_mu[d][proc].GetBinWidth(1)
        ntracks_pu = fpu[d][proc].Integral(-10,10)/h_tracks_dt_mu[d][proc].GetBinWidth(1)

        ntracks_pv_err = fpv[d][proc].IntegralError(-10,10)/h_tracks_dt_mu[d][proc].GetBinWidth(1)
        ntracks_pu_err = fpu[d][proc].IntegralError(-10,10)/h_tracks_dt_mu[d][proc].GetBinWidth(1)

        print 'sigma_t PV tracks : %.01f +/- %.01f'%(sigmat_pv, sigmat_pv_err)
        print 'sigma_t PU tracks : %.01f +/- %.01f'%(sigmat_pu, sigmat_pu_err)
        
        print 'Number of PV tracks per muon : %.02f +/- %.02f'%(ntracks_pv, ntracks_pv_err )
        print 'Number of PU tracks per muon : %.02f +/- %.02f'%(ntracks_pu, ntracks_pu_err) 
    
        
        tpv[d][proc] = ROOT.TLatex( 0.12, 0.82, '#splitline{PV tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pv, ntracks_pv))
        tpv[d][proc].SetTextColor(ROOT.kGreen+1)
        tpv[d][proc].SetNDC()
        tpv[d][proc].SetTextSize(0.030)
        
        tpu[d][proc] = ROOT.TLatex( 0.12, 0.72, '#splitline{PU tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pu, ntracks_pu))
        tpu[d][proc].SetTextColor(2)
        tpu[d][proc].SetNDC()
        tpu[d][proc].SetTextSize(0.030)
        
        tpu[d][proc].Draw()
        tpv[d][proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        #raw_input('ok?')




# dt wrt vtx
ccc  = {}
fitfun = {}
ffpv = {}
ffpu = {}
ttpu = {}
ttpv = {}

for d in 'barrel','endcap':
    ccc[d] = {}
    fitfun[d] = {}
    ffpv[d] = {}
    ffpu[d] = {}
    ttpu[d] = {}
    ttpv[d] = {}
    for proc in ['prompt','fake']:
        cname = h_tracks_dt_vtx[d]['prompt'].GetName().replace('h_','') + '_' + proc
        ccc[d][proc]= ROOT.TCanvas(cname,cname,500,500)
        h_tracks_dt_vtx[d][proc].GetYaxis().SetRangeUser(0.0000001,h_tracks_dt_vtx[d][proc].GetMaximum()*1.5)
        h_tracks_dt_vtx[d][proc].GetXaxis().SetTitle('t_{track} - t_{vtx4D} (ns)')
        fitfun[d][proc] = h_tracks_dt_vtx[d][proc].GetFunction('fitfun_%s'%d)
        fitfun[d][proc].SetLineColor(ROOT.kGreen+1)
        fitfun[d][proc].SetRange(-1,1)
        #fitfun[d][proc].FixParameter(2, 0.030)
        #fitfun[d][proc].FixParameter(5, 0.270)
        fitfun[d][proc].SetParameter(2, 0.030)
        fitfun[d][proc].SetParameter(5, 0.270)
        #h_tracks_dt_vtx[d][proc].Rebin(4)
        h_tracks_dt_vtx[d][proc].Fit( fitfun[d][proc],'QRMS')
        h_tracks_dt_vtx[d][proc].Draw()
        
        ffpv[d][proc] = ROOT.TF1('ffpv_%s_%s'%(proc,d),'gaus(0)',-1,1)
        ffpv[d][proc].SetParameters(fitfun[d][proc].GetParameter(0), fitfun[d][proc].GetParameter(1), fitfun[d][proc].GetParameter(2))
        
        ffpu[d][proc] = ROOT.TF1('ffpu_%s_%s'%(proc,d),'gaus(0)',-1,1)
        ffpu[d][proc].SetParameters(fitfun[d][proc].GetParameter(3), fitfun[d][proc].GetParameter(4), fitfun[d][proc].GetParameter(5))
        ffpu[d][proc].SetLineColor(2)
        ffpu[d][proc].Draw('same')

        sigmat_pv = fitfun[d][proc].GetParameter(2)*1000.
        sigmat_pu = fitfun[d][proc].GetParameter(5)*1000.

        sigmat_pv_err = fitfun[d][proc].GetParError(2)*1000.
        sigmat_pu_err = fitfun[d][proc].GetParError(5)*1000.
        
        ntracks_pv = ffpv[d][proc].Integral(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)
        ntracks_pu = ffpu[d][proc].Integral(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)

        ntracks_pv_err = ffpv[d][proc].IntegralError(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)
        ntracks_pu_err = ffpu[d][proc].IntegralError(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)

        print 'sigma_t PV tracks : %.01f +/- %.01f'%(sigmat_pv, sigmat_pv_err)
        print 'sigma_t PU tracks : %.01f +/- %.01f'%(sigmat_pu, sigmat_pu_err)
        
        print 'Number of PV tracks per muon : %.02f +/- %.02f'%(ntracks_pv, ntracks_pv_err )
        print 'Number of PU tracks per muon : %.02f +/- %.02f'%(ntracks_pu, ntracks_pu_err) 
    
        
        ttpv[d][proc] = ROOT.TLatex( 0.12, 0.82, '#splitline{PV tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pv, ntracks_pv))
        ttpv[d][proc].SetTextColor(ROOT.kGreen+1)
        ttpv[d][proc].SetNDC()
        ttpv[d][proc].SetTextSize(0.030)
        
        ttpu[d][proc] = ROOT.TLatex( 0.12, 0.72, '#splitline{PU tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pu, ntracks_pu))
        #if (  ntracks_pu < 0.01 or sigmat_pu < 0 ) :
        #    ttpu[d][proc] = ROOT.TLatex( 0.12, 0.72, '#splitline{PU tracks}{#sigma_{t} = n.a. ps, <n/ev> = 0}')

        ttpu[d][proc].SetTextColor(2)
        ttpu[d][proc].SetNDC()
        ttpu[d][proc].SetTextSize(0.030)
        
        ttpu[d][proc].Draw()
        ttpv[d][proc].Draw()
        tl.Draw()
        tl2[d].Draw()

        #raw_input('ok?')

# save canvases
for cname in c:
    c[cname].SaveAs(plotdir+'/'+c[cname].GetName()+'.png')
    c[cname].SaveAs(plotdir+'/'+c[cname].GetName()+'.pdf')

for d in ['barrel', 'endcap']:
    for proc in ['prompt','fake']:
        cc[d][proc].SaveAs(plotdir+'/'+cc[d][proc].GetName()+'.png')
        cc[d][proc].SaveAs(plotdir+'/'+cc[d][proc].GetName()+'.pdf')
        ccc[d][proc].SaveAs(plotdir+'/'+ccc[d][proc].GetName()+'.png')
        ccc[d][proc].SaveAs(plotdir+'/'+ccc[d][proc].GetName()+'.pdf')
  
        

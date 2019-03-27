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

#release = '93X'
release = '10_4_0_mtd5'


pu = sys.argv[1]
bkgProc = sys.argv[2]

fSig = '../'+release+'/testTracksMva_DYToLL_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt.root'
fBkg = '../'+release+'/testTracksMva_%s_PU200_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt.root'%bkgProc


if (pu == 'noPU'):
    fSig = '../'+release+'/testTracksMva_DYToLL_noPU_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt.root'
    fBkg = '../'+release+'/testTracksMva_DYToLL_noPU_dT3sigma_minMuonPt20_maxMuonPt9999_minTrackPt.root'
    #fBkg = '../'+release+'/testTracks_%s_noPU_dT3sigma_minTrackPt.root'%bkgProc
    #fSig = '../'+release+'/testTracks_DYToLL_noPU_dT3sigma_minMuonPt10_maxMuonPt15_minTrackPt.root'
    #fBkg = '../'+release+'/testTracks_%s_noPU_dT3sigma_minMuonPt10_maxMuonPt15_minTrackPt.root'%bkgProc
    
    
plotdir = release+'/tracksPuidMva_%s_%s'%(pu,bkgProc)
if ('minTrackPt' in fSig):
    plotdir = plotdir+'_minTrackPt'
    
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

h_tracks_dt_vtx = {}
h_tracks_dz_vtx = {}
h2_tracks_dzvtx_dtvtx = {}
h_tracks_dz_mu = {}
h_tracks_dt_mu = {}

h_tracks_pt = {}
h_tracks_n = {}
h_tracks_sumpt = {}
p_tracks_pt = {}
p_tracks_n = {}
p_tracks_sumpt = {}

h_tracks_removed_pt = {}
h_tracks_removed_n = {}
h_tracks_removed_sumpt = {}
p_tracks_removed_pt = {}
p_tracks_removed_n = {}
p_tracks_removed_sumpt = {}

h_tracks_kept_pt = {}
h_tracks_kept_n = {}
h_tracks_kept_sumpt = {}
p_tracks_kept_pt = {}
p_tracks_kept_n = {}
p_tracks_kept_sumpt = {}



for d in ['barrel','endcap']:
    h_tracks_dt_vtx[d] = {}
    h_tracks_dz_vtx[d] = {}
    h_tracks_dz_mu[d] = {}
    h_tracks_dt_mu[d] = {}
    h2_tracks_dzvtx_dtvtx[d] = {}

    h_tracks_pt[d]  = {}
    h_tracks_n[d]  = {}
    h_tracks_sumpt[d] = {}
    p_tracks_pt[d]  = {}
    p_tracks_n[d]  = {}
    p_tracks_sumpt[d]  = {}

    h_tracks_removed_pt[d]  = {}
    h_tracks_removed_n[d]  = {}
    h_tracks_removed_sumpt[d]  = {}
    p_tracks_removed_pt[d]  = {}
    p_tracks_removed_n[d]  = {}
    p_tracks_removed_sumpt[d]  = {}

    h_tracks_kept_pt[d]  = {}
    h_tracks_kept_n[d]  = {}
    h_tracks_kept_sumpt[d]  = {}
    p_tracks_kept_pt[d]  = {}
    p_tracks_kept_n[d]  = {}
    p_tracks_kept_sumpt[d]  = {}

    for proc in ['prompt', 'fake']:
        h_tracks_dz_vtx[d][proc]        =  f[proc].Get('h_tracks_dz_vtx_'+d)
        h_tracks_dt_vtx[d][proc]        =  f[proc].Get('h_tracks_dt_vtx_'+d)
        h2_tracks_dzvtx_dtvtx[d][proc] = f[proc].Get('h2_tracks_dzvtx_dtvtx_'+d)

        h_tracks_pt[d][proc]    = f[proc].Get('h_tracks_pt_'+d)
        h_tracks_n[d][proc]     = f[proc].Get('h_tracks_n_'+d)
        h_tracks_sumpt[d][proc] = f[proc].Get('h_tracks_sumpt_'+d)
        p_tracks_pt[d][proc]    = f[proc].Get('p_tracks_pt_vs_linedensity_'+d)
        p_tracks_n[d][proc]     = f[proc].Get('p_tracks_n_vs_linedensity_'+d)
        p_tracks_sumpt[d][proc] = f[proc].Get('p_tracks_sumpt_vs_linedensity_'+d)

        h_tracks_removed_pt[d][proc]    = f[proc].Get('h_tracks_removed_pt_'+d)
        h_tracks_removed_n[d][proc]     = f[proc].Get('h_tracks_removed_n_'+d)
        h_tracks_removed_sumpt[d][proc] = f[proc].Get('h_tracks_removed_sumpt_'+d)
        p_tracks_removed_pt[d][proc]    = f[proc].Get('p_tracks_removed_pt_vs_linedensity_'+d)
        p_tracks_removed_n[d][proc]     = f[proc].Get('p_tracks_removed_n_vs_linedensity_'+d)
        p_tracks_removed_sumpt[d][proc] = f[proc].Get('p_tracks_removed_sumpt_vs_linedensity_'+d)

        h_tracks_kept_pt[d][proc]    = f[proc].Get('h_tracks_kept_pt_'+d)
        h_tracks_kept_n[d][proc]     = f[proc].Get('h_tracks_kept_n_'+d)
        h_tracks_kept_sumpt[d][proc] = f[proc].Get('h_tracks_kept_sumpt_'+d)
        p_tracks_kept_pt[d][proc]    = f[proc].Get('p_tracks_kept_pt_vs_linedensity_'+d)
        p_tracks_kept_n[d][proc]     = f[proc].Get('p_tracks_kept_n_vs_linedensity_'+d)
        p_tracks_kept_sumpt[d][proc] = f[proc].Get('p_tracks_kept_sumpt_vs_linedensity_'+d)

        

xtitle = {'h_tracks_dz_vtx' : 'z_{track} - z_{vtx} (cm)',
          'h_tracks_dz_mu' : 'z_{#mu} - z_{vtx} (cm)',
          'h_tracks_pt' : 'p^{T} (GeV)',
          'h_tracks_n' : 'number of tracks',
          'h_tracks_sumpt' : '#Sigma p^{T} in iso cone (GeV)',
          'p_tracks_pt_vs_linedensity' : 'line density (mm^{-1})',
          'p_tracks_n_vs_linedensity' : 'line density (mm^{-1})',
          'p_tracks_sumpt_vs_linedensity' : 'line density (mm^{-1})',
}

ytitle = {'h_tracks_dz_vtx' : 'a.u.',
          'h_tracks_dz_mu' : 'a.u.',
          'h_tracks_pt' : 'a.u.',
          'h_tracks_n' : 'a.u.',
          'h_tracks_sumpt' : 'a.u.',
          'p_tracks_pt_vs_linedensity' : '< p^{T} > (GeV)',
          'p_tracks_n_vs_linedensity' : ' < n_{tracks} >',
          'p_tracks_sumpt_vs_linedensity' : ' < #Sigma p^{T} > (GeV)',
          }


        
leg1 = ROOT.TLegend(0.15, 0.7, 0.45, 0.89)
leg1.SetBorderSize(0)



c = {}


# compare prompt and non-prompt
for id,d in enumerate(['barrel', 'endcap']):
    for ih,h in  enumerate([h_tracks_dz_vtx[d],
                            h_tracks_pt[d], h_tracks_n[d], h_tracks_sumpt[d],
                            h_tracks_removed_pt[d], h_tracks_removed_n[d], h_tracks_removed_sumpt[d],
                            h_tracks_kept_pt[d], h_tracks_kept_n[d], h_tracks_kept_sumpt[d],
                            p_tracks_pt[d], p_tracks_n[d], p_tracks_sumpt[d],
                            p_tracks_removed_pt[d], p_tracks_removed_n[d], p_tracks_removed_sumpt[d],
                            p_tracks_kept_pt[d], p_tracks_kept_n[d], p_tracks_kept_sumpt[d]   ]):
        print h['prompt'].GetName()
        cname = h['prompt'].GetName().replace('h_','')
        c[cname] = ROOT.TCanvas(cname,cname)
        h['prompt'].SetLineColor(ROOT.kBlue)
        h['fake'].SetLineColor(ROOT.kRed)
        h['prompt'].SetMarkerColor(ROOT.kBlue)
        h['fake'].SetMarkerColor(ROOT.kRed)
        h['prompt'].SetMarkerStyle(20)
        h['fake'].SetMarkerStyle(20)
        h['prompt'].GetYaxis().SetRangeUser(0.01,h['prompt'].GetMaximum()*1.3)
        h['prompt'].GetYaxis().SetRangeUser(0.01,h['prompt'].GetMaximum()*1.3)
        thisname = h['prompt'].GetName().replace('_barrel','').replace('_endcap','').replace('_removed','').replace('_kept','')
        h['prompt'].GetXaxis().SetTitle(xtitle[thisname])
        h['prompt'].GetYaxis().SetTitle(ytitle[thisname])
        h['fake'].GetXaxis().SetTitle(xtitle[thisname])
        h['fake'].GetYaxis().SetTitle(ytitle[thisname])
        #print thisname, xtitle[thisname], ytitle[thisname]
        if (thisname.startswith('h_tracks_sumpt') or thisname.startswith('h_tracks_removed_sumpt') or thisname.startswith('h_tracks_kept_sumpt') ):
            c[cname].SetLogy()
        if (h['prompt'].GetName().startswith('h_')):
            h['prompt'].DrawNormalized()
            h['fake'].DrawNormalized('histo same')
        elif (h['prompt'].GetName().startswith('p_')):
            h['fake'].GetYaxis().SetRangeUser(0.00001,h['fake'].GetMaximum()*1.3)
            h['fake'].GetYaxis().SetRangeUser(0.00001,h['fake'].GetMaximum()*1.3)
            h['fake'].Draw()
            h['prompt'].Draw('esame')
        if (id == 0 and ih == 0):
            leg1.AddEntry(h['prompt'],'prompt muons (Z#rightarrow#mu#mu)','PL')
            leg1.AddEntry(h['fake'],'non-prompt muons (t#bar{t})','PL')
            leg1.Draw('same')
        else:
            leg1.Draw('same')
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c[cname], iPeriod, iPos)
        c[cname].Update()
        #raw_input('ok?')



# compare all tracks in cone with removed tracks
c_tracks_pt = {} 
c_tracks_n = {}
c_tracks_sumpt = {}
c_tracks_pt_vs_linedensity = {}
c_tracks_n_vs_linedensity = {}
c_tracks_sumpt_vs_linedensity = {}

leg2 = {}
leg3 = {}

for id, d in enumerate(['barrel', 'endcap']):
    c_tracks_pt[d] = {}
    c_tracks_n[d] = {}
    c_tracks_sumpt[d] = {}
    c_tracks_pt_vs_linedensity[d] = {}
    c_tracks_n_vs_linedensity[d] = {}
    c_tracks_sumpt_vs_linedensity[d] = {}
    for ip,proc in enumerate(['prompt', 'fake']):
    
        c_tracks_pt[d][proc] = ROOT.TCanvas('tracks_pt_%s_%s'%(proc,d),'tracks_pt_%s_%s'%(proc,d))
        c_tracks_n[d][proc] = ROOT.TCanvas('tracks_n_%s_%s'%(proc,d), 'tracks_n_%s_%s'%(proc,d))
        c_tracks_sumpt[d][proc] = ROOT.TCanvas('tracks_sumpt_%s_%s'%(proc,d),'tracks_sumpt_%s_%s'%(proc,d),500,500 )

        c_tracks_pt_vs_linedensity[d][proc] = ROOT.TCanvas('tracks_pt_vs_linedensity_%s_%s'%(proc,d),'tracks_pt_vs_linedensity_%s_%s'%(proc,d))
        c_tracks_n_vs_linedensity[d][proc] = ROOT.TCanvas('tracks_n_vs_linedensity_%s_%s'%(proc,d), 'tracks_n_vs_linedensity_%s_%s'%(proc,d))
        c_tracks_sumpt_vs_linedensity[d][proc] = ROOT.TCanvas('tracks_sumpt_vs_linedensity_%s_%s'%(proc,d),'tracks_vs_linedensity_sumpt_%s_%s'%(proc,d),500,500 )
        
        c_tracks_pt[d][proc].cd()
        h_tracks_pt[d][proc].GetXaxis().SetRangeUser(0.0, 10.)
        h_tracks_pt[d][proc].GetYaxis().SetRangeUser(0.0001,h_tracks_pt[d][proc].GetMaximum()*1.2)
        h_tracks_pt[d][proc].Draw()
        h_tracks_removed_pt[d][proc].SetLineStyle(2)
        h_tracks_removed_pt[d][proc].Draw('same')
        if (id == 0):
            leg2[proc] = ROOT.TLegend(0.15, 0.7, 0.55, 0.89)
            leg2[proc].SetBorderSize(0)
            leg2[proc].AddEntry( h_tracks_pt[d][proc], 'all tracks in iso cone', 'PL')
            leg2[proc].AddEntry( h_tracks_removed_pt[d][proc], 'tracks removed from iso cone', 'PL')
        leg2[proc].Draw() 
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_pt[d][proc], iPeriod, iPos)
         
        c_tracks_n[d][proc].cd()
        h_tracks_n[d][proc].GetXaxis().SetRangeUser(0,10)
        h_tracks_n[d][proc].GetYaxis().SetRangeUser(0.0001,h_tracks_n[d][proc].GetMaximum()*1.2)
        h_tracks_n[d][proc].Draw()
        h_tracks_removed_n[d][proc].SetLineStyle(2)
        h_tracks_removed_n[d][proc].Draw('same')
        leg2[proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_n[d][proc], iPeriod, iPos)
         
        c_tracks_sumpt[d][proc].cd()
        c_tracks_sumpt[d][proc].SetLogy()
        h_tracks_sumpt[d][proc].GetXaxis().SetRangeUser(0.0, 50.)
        h_tracks_sumpt[d][proc].GetYaxis().SetRangeUser(0.01,h_tracks_sumpt[d][proc].GetMaximum()*1.2)
        h_tracks_sumpt[d][proc].Draw()
        h_tracks_removed_sumpt[d][proc].SetLineStyle(2)
        h_tracks_removed_sumpt[d][proc].Draw('same')
        leg2[proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_sumpt[d][proc], iPeriod, iPos)
         
        
        c_tracks_pt_vs_linedensity[d][proc].cd()
        p_tracks_pt[d][proc].GetYaxis().SetRangeUser(0.0001,p_tracks_pt[d][proc].GetMaximum()*1.2)
        p_tracks_pt[d][proc].Draw()
        p_tracks_removed_pt[d][proc].SetMarkerStyle(24)
        p_tracks_removed_pt[d][proc].Draw('same')
        if (id == 0 ):
            leg3[proc] = ROOT.TLegend(0.15, 0.7, 0.55, 0.89)
            leg3[proc].SetBorderSize(0)
            leg3[proc].AddEntry( p_tracks_pt[d][proc], 'all tracks in iso cone', 'PL')
            leg3[proc].AddEntry( p_tracks_removed_pt[d][proc], 'tracks removed from iso cone', 'PL')
        leg3[proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_pt_vs_linedensity[d][proc], iPeriod, iPos)
        
        c_tracks_n_vs_linedensity[d][proc].cd()
        p_tracks_n[d][proc].GetYaxis().SetRangeUser(0.0001,p_tracks_n[d][proc].GetMaximum()*1.2)
        p_tracks_n[d][proc].Draw()
        p_tracks_removed_n[d][proc].SetMarkerStyle(24)
        p_tracks_removed_n[d][proc].Draw('same')
        leg3[proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_n_vs_linedensity[d][proc], iPeriod, iPos)
        
        c_tracks_sumpt_vs_linedensity[d][proc].cd()
        p_tracks_sumpt[d][proc].GetYaxis().SetRangeUser(0.0001,p_tracks_sumpt[d][proc].GetMaximum()*1.2)
        p_tracks_sumpt[d][proc].Draw()
        p_tracks_removed_sumpt[d][proc].SetMarkerStyle(24)
        p_tracks_removed_sumpt[d][proc].Draw('same')
        leg3[proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        CMS_lumi.CMS_lumi(c_tracks_sumpt_vs_linedensity[d][proc], iPeriod, iPos)
        # raw_input('ok?')
        

# 2D plots
for id,d in enumerate(['barrel', 'endcap']):
    for proc in ['prompt','fake']:
        cname = h2_tracks_dzvtx_dtvtx[d][proc].GetName().replace('h2_','')+'_'+proc
        c[cname] = ROOT.TCanvas(cname,cname)
        h2_tracks_dzvtx_dtvtx[d][proc].GetXaxis().SetRangeUser(-0.1,0.1)
        h2_tracks_dzvtx_dtvtx[d][proc].GetYaxis().SetRangeUser(-0.2,0.2)
        h2_tracks_dzvtx_dtvtx[d][proc].GetXaxis().SetTitle('z_{track} - z_{vtx4D}')
        h2_tracks_dzvtx_dtvtx[d][proc].GetYaxis().SetTitle('t_{track} - t_{vtx4D}')
        h2_tracks_dzvtx_dtvtx[d][proc].Draw('colz')
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
        ccc[d][proc]= ROOT.TCanvas(cname,cname) 
        ccc[d][proc].SetLogy()
        ROOT.gStyle.SetOptFit(0)
        h_tracks_dt_vtx[d][proc].GetYaxis().SetRangeUser(0.00001,h_tracks_dt_vtx[d][proc].GetMaximum()*1.5)
        #h_tracks_dt_vtx[d][proc].GetXaxis().SetTitle('t_{track} - t_{vtx4D} (ns)')
        h_tracks_dt_vtx[d][proc].GetXaxis().SetTitle('t_{track} - t_{vtxGen} (ns)')
        fitfun[d][proc] = h_tracks_dt_vtx[d][proc].GetFunction('fitfun_%s'%d)
        #fitfun[d][proc] = ROOT.TF1('fitfun_%s'%d,'gaus(0)+gaus(3)+gaus(6)')
        fitfun[d][proc].SetLineColor(ROOT.kGreen+1)
        fitfun[d][proc].SetRange(-1,1)
        fitfun[d][proc].SetParameter(2, 0.030)
        fitfun[d][proc].SetParameter(5, 0.270)
        fitfun[d][proc].SetParameter(8, 0.050)
        fitfun[d][proc].FixParameter(1, 0.0)
        fitfun[d][proc].FixParameter(4, 0.0)
        fitfun[d][proc].FixParameter(7, 0.0)
        h_tracks_dt_vtx[d][proc].Fit( fitfun[d][proc],'QRMS')
        h_tracks_dt_vtx[d][proc].Draw()
        
        ffpv[d][proc] = ROOT.TF1('ffpv_%s_%s'%(proc,d),'gaus(0)',-1,1)
        ffpv[d][proc].SetParameters(fitfun[d][proc].GetParameter(0), fitfun[d][proc].GetParameter(1), fitfun[d][proc].GetParameter(2))
        
        ffpu[d][proc] = ROOT.TF1('ffpu_%s_%s'%(proc,d),'gaus(0)',-1,1)
        ffpu[d][proc].SetParameters(fitfun[d][proc].GetParameter(3), fitfun[d][proc].GetParameter(4), fitfun[d][proc].GetParameter(5))
        ffpu[d][proc].SetLineColor(2)
        ffpu[d][proc].Draw('same')

        # 3 gaus
        #ffpv[d][proc] = ROOT.TF1('ffpv_%s_%s'%(proc,d),'gaus(0)+gaus(3)',-1,1)
        #ffpv[d][proc].SetParameters(fitfun[d][proc].GetParameter(0), fitfun[d][proc].GetParameter(1), fitfun[d][proc].GetParameter(2),
        #                            fitfun[d][proc].GetParameter(6), fitfun[d][proc].GetParameter(7), fitfun[d][proc].GetParameter(8))
        
        #ffpu[d][proc] = ROOT.TF1('ffpu_%s_%s'%(proc,d),'gaus(0)',-1,1)
        #ffpu[d][proc].SetParameters(fitfun[d][proc].GetParameter(3), fitfun[d][proc].GetParameter(4), fitfun[d][proc].GetParameter(5))
        #ffpu[d][proc].SetLineColor(2)
        #ffpu[d][proc].Draw('same')

        sigmat_pv = fitfun[d][proc].GetParameter(2)*1000.
        sigmat_pu = fitfun[d][proc].GetParameter(5)*1000.

        sigmat_pv_err = fitfun[d][proc].GetParError(2)*1000.
        sigmat_pu_err = fitfun[d][proc].GetParError(5)*1000.
        
        ntracks_pv = ffpv[d][proc].Integral(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)
        ntracks_pu = ffpu[d][proc].Integral(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)
        #ntracks_pv = h_tracks_dt_vtx[d][proc].Integral() - ntracks_pu

        ntracks_pv_err = ffpv[d][proc].IntegralError(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)
        ntracks_pu_err = ffpu[d][proc].IntegralError(-10,10)/h_tracks_dt_vtx[d][proc].GetBinWidth(1)

        print 'sigma_t PV tracks : %.01f +/- %.01f'%(sigmat_pv, sigmat_pv_err)
        print 'sigma_t PU tracks : %.01f +/- %.01f'%(sigmat_pu, sigmat_pu_err)
        
        print 'Number of PV tracks per muon : %.02f +/- %.02f'%(ntracks_pv, ntracks_pv_err )
        print 'Number of PU tracks per muon : %.02f +/- %.02f'%(ntracks_pu, ntracks_pu_err) 
    
        
        ttpv[d][proc] = ROOT.TLatex( 0.15, 0.82, '#splitline{PV tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pv, ntracks_pv))
        ttpv[d][proc].SetTextColor(ROOT.kGreen+1)
        ttpv[d][proc].SetNDC()
        ttpv[d][proc].SetTextSize(0.030)
        
        ttpu[d][proc] = ROOT.TLatex( 0.15, 0.72, '#splitline{PU tracks}{#sigma_{t} = %.0f ps, <n/ev> = %.02f}'%(sigmat_pu, ntracks_pu))
        #if (  ntracks_pu < 0.01 or sigmat_pu < 0 ) :
        #    ttpu[d][proc] = ROOT.TLatex( 0.12, 0.72, '#splitline{PU tracks}{#sigma_{t} = n.a. ps, <n/ev> = 0}')

        ttpu[d][proc].SetTextColor(2)
        ttpu[d][proc].SetNDC()
        ttpu[d][proc].SetTextSize(0.030)
        
        ttpu[d][proc].Draw()
        ttpv[d][proc].Draw()
        tl.Draw()
        tl2[d].Draw()
        tl3[proc].Draw()
        CMS_lumi.CMS_lumi(ccc[d][proc], iPeriod, iPos)
        #raw_input('ok?')





# make rocs
h_relChIso03_dZ1 = {}
h_relChIso03_dZ1_dT3s = {}

g_roc_relChIso03_dZ1 = {}
g_roc_relChIso03_dZ1_dT3s = {}

c_roc_relChIso03_dZ1_dT3s = {}

h_relChIso03_mva3d= {}
h_relChIso03_mva4d = {}

g_roc_relChIso03_mva3d = {}
g_roc_relChIso03_mva4d = {}

c_roc_relChIso03_mva = {}

leg4 = {}
for id, d in enumerate(['barrel','endcap']):
    h_relChIso03_dZ1[d] = {}
    h_relChIso03_dZ1_dT3s[d] = {}
    h_relChIso03_mva3d[d] = {}
    h_relChIso03_mva4d[d] = {}
    for proc in ['prompt','fake']:
        h_relChIso03_dZ1[d][proc] = f[proc].Get('h_muon_relChIso03_dZ_%s'%d)
        h_relChIso03_dZ1_dT3s[d][proc] = f[proc].Get('h_muon_relChIso03_dZ_dT3s_%s'%d)
        h_relChIso03_mva3d[d][proc] = f[proc].Get('h_muon_relChIso03_mva3D_%s'%d)
        h_relChIso03_mva4d[d][proc] = f[proc].Get('h_muon_relChIso03_mva4D_%s'%d)
                 
    g_roc_relChIso03_dZ1[d] = ROOT.TGraphErrors()
    g_roc_relChIso03_dZ1[d].SetName('g_roc_relChIso03_dZ1_%s'%d)
    g_roc_relChIso03_dZ1_dT3s[d] = ROOT.TGraphErrors()
    g_roc_relChIso03_dZ1_dT3s[d].SetName('g_roc_relChIso03_dZ1_dT3s_%s'%d)

    g_roc_relChIso03_mva3d[d] = ROOT.TGraphErrors()
    g_roc_relChIso03_mva3d[d].SetName('g_roc_relChIso03_mva3d_%s'%d)
    g_roc_relChIso03_mva4d[d] = ROOT.TGraphErrors()
    g_roc_relChIso03_mva4d[d].SetName('g_roc_relChIso03_mva4d_%s'%d)
 
    makeRoc(h_relChIso03_dZ1[d]['prompt'], h_relChIso03_dZ1[d]['fake'], g_roc_relChIso03_dZ1[d] )
    makeRoc(h_relChIso03_dZ1_dT3s[d]['prompt'], h_relChIso03_dZ1_dT3s[d]['fake'],g_roc_relChIso03_dZ1_dT3s[d] )

    makeRoc(h_relChIso03_mva3d[d]['prompt'], h_relChIso03_mva3d[d]['fake'], g_roc_relChIso03_mva3d[d] )
    makeRoc(h_relChIso03_mva4d[d]['prompt'], h_relChIso03_mva4d[d]['fake'], g_roc_relChIso03_mva4d[d] )
        
 
    c_roc_relChIso03_dZ1_dT3s[d] = ROOT.TCanvas('roc_relChIso03_dZ1_dT3s_%s'%d, 'roc_relChIso03_dZ1_dT3s_%s'%d)
    c_roc_relChIso03_dZ1_dT3s[d].SetGridx()
    c_roc_relChIso03_dZ1_dT3s[d].SetGridy()
    g_roc_relChIso03_dZ1[d].GetXaxis().SetRangeUser(0.85,1.01)
    g_roc_relChIso03_dZ1[d].GetYaxis().SetRangeUser(0.0,0.1)
    g_roc_relChIso03_dZ1[d].GetXaxis().SetTitle("Prompt efficiency")
    g_roc_relChIso03_dZ1[d].GetYaxis().SetTitle("Non-prompt efficiency")
    g_roc_relChIso03_dZ1[d].SetLineColor(ROOT.kBlue)
    g_roc_relChIso03_dZ1[d].SetFillColorAlpha(ROOT.kBlue,0.35)
    g_roc_relChIso03_dZ1[d].SetFillStyle(1001)
    g_roc_relChIso03_dZ1_dT3s[d].SetLineColor(ROOT.kRed)
    g_roc_relChIso03_dZ1_dT3s[d].SetFillColorAlpha(ROOT.kRed,0.35)
    g_roc_relChIso03_dZ1_dT3s[d].SetFillStyle(1001)
    g_roc_relChIso03_dZ1[d].Draw('A L E3')
    g_roc_relChIso03_dZ1_dT3s[d].Draw('L E3 same')
    g_roc_relChIso03_mva3d[d].SetLineColor(ROOT.kGreen+2)
    g_roc_relChIso03_mva3d[d].SetFillColorAlpha(ROOT.kGreen+2,0.35)
    g_roc_relChIso03_mva3d[d].SetFillStyle(1001)
    g_roc_relChIso03_mva4d[d].SetLineColor(ROOT.kMagenta+2)
    g_roc_relChIso03_mva4d[d].SetFillColorAlpha(ROOT.kMagenta+2,0.35)
    g_roc_relChIso03_mva4d[d].SetFillStyle(1001)
    g_roc_relChIso03_mva3d[d].SetLineStyle(2)
    g_roc_relChIso03_mva4d[d].SetLineStyle(2)
    g_roc_relChIso03_mva3d[d].Draw('L E3 same')
    g_roc_relChIso03_mva4d[d].Draw('L E3 same')
    leg4[d] = ROOT.TLegend(0.15, 0.7, 0.55, 0.92)
    leg4[d].SetBorderSize(0)
    leg4[d].AddEntry( g_roc_relChIso03_mva3d[d],'no MTD, 3DMVA','L')
    leg4[d].AddEntry( g_roc_relChIso03_mva4d[d],'MTD, 4DMVA (#sigma_{t} = 35 ps)','L')
    leg4[d].AddEntry( g_roc_relChIso03_dZ1[d],'no MTD, dz < 1 mm','L')
    leg4[d].AddEntry( g_roc_relChIso03_dZ1_dT3s[d],'MTD, dz < 1 mm, dt < 3#sigma_{t} (#sigma_{t} = 35 ps)','L')
    tl.Draw()
    tl2[d].Draw()
    leg4[d].Draw()
    CMS_lumi.CMS_lumi(c_roc_relChIso03_dZ1_dT3s[d], iPeriod, iPos)
    
    #raw_input('ok?')


# save canvases
for cname in c:
    c[cname].SaveAs(plotdir+'/'+c[cname].GetName()+'.png')
    c[cname].SaveAs(plotdir+'/'+c[cname].GetName()+'.pdf')

for d in ['barrel', 'endcap']:
    c_roc_relChIso03_dZ1_dT3s[d].SaveAs(plotdir+'/'+c_roc_relChIso03_dZ1_dT3s[d].GetName()+'.png')
    c_roc_relChIso03_dZ1_dT3s[d].SaveAs(plotdir+'/'+c_roc_relChIso03_dZ1_dT3s[d].GetName()+'.pdf')
    
    for proc in ['prompt','fake']:
        ccc[d][proc].SaveAs(plotdir+'/'+ccc[d][proc].GetName()+'.png')
        ccc[d][proc].SaveAs(plotdir+'/'+ccc[d][proc].GetName()+'.pdf')

        c_tracks_pt[d][proc].SaveAs(plotdir+'/'+c_tracks_pt[d][proc].GetName()+'.png')
        c_tracks_pt[d][proc].SaveAs(plotdir+'/'+c_tracks_pt[d][proc].GetName()+'.pdf')

        c_tracks_n[d][proc].SaveAs(plotdir+'/'+c_tracks_n[d][proc].GetName()+'.png')
        c_tracks_n[d][proc].SaveAs(plotdir+'/'+c_tracks_n[d][proc].GetName()+'.pdf')

        c_tracks_sumpt[d][proc].SaveAs(plotdir+'/'+c_tracks_sumpt[d][proc].GetName()+'.png')
        c_tracks_sumpt[d][proc].SaveAs(plotdir+'/'+c_tracks_sumpt[d][proc].GetName()+'.pdf')

        c_tracks_pt_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_pt_vs_linedensity[d][proc].GetName()+'.png')
        c_tracks_pt_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_pt_vs_linedensity[d][proc].GetName()+'.pdf')

        c_tracks_n_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_n_vs_linedensity[d][proc].GetName()+'.png')
        c_tracks_n_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_n_vs_linedensity[d][proc].GetName()+'.pdf')

        c_tracks_sumpt_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_sumpt_vs_linedensity[d][proc].GetName()+'.png')
        c_tracks_sumpt_vs_linedensity[d][proc].SaveAs(plotdir+'/'+c_tracks_sumpt_vs_linedensity[d][proc].GetName()+'.pdf')



fout = ROOT.TFile(plotdir+'/roc_30ps_%s.root'%(pu),'recreate')
for d in ['barrel', 'endcap']:
    g_roc_relChIso03_dZ1[d].Write()
    g_roc_relChIso03_dZ1_dT3s[d].Write()
    g_roc_relChIso03_mva3d[d].Write()
    g_roc_relChIso03_mva4d[d].Write()
   

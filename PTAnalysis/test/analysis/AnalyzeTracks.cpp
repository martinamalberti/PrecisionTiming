// per compilare: g++ -Wall -o AnalyzeTracks `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit AnalyzeTracks.cpp                        

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
//#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>


using namespace std;

int main(int argc, char** argv){

  float timeResolution = 0.030;
  int nsigma = 3;

  float minMuonPt = 20.;
  float maxMuonPt = 9999.;
  float maxDz = 0.1;
  float maxDxy = 0.02;
  float maxDzMu  = 0.5;
  float btlMinTrackPt = 0.7;
  float etlMinTrackPt = 0.4;
  
  string process = argv[1];
  bool prompt = false;

  string pu = argv[2];
  
  // -- get TChain
  TChain* chain = new TChain("analysis/tree_30ps","tree");
  if (process.find("DYToLL") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_muIso/190111_164118/0000/muonIsolation_*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_noPU_muIso/190111_164205/0000/muonIsolation_*.root");
    prompt = true;
  }
  
  if (process.find("TTbar") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_muIso/190111_164134/0000/muonIsolation_*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_noPU_muIso/190111_164224/0000/muonIsolation_*.root");
    prompt = false;
  }
  
  if (process.find("QCD") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/test_QCD_muIso/190111_164150/0000/muonIsolation_*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/test_QCD_noPU_muIso/190111_164240/0000/muonIsolation_*.root");
    prompt = false; 
  }
  
  cout << "Using prompt muons = " << prompt <<endl;
  
  
  // -- tree vars
  int npu;
  float vtx3D_z;
  float vtx4D_z;
  float vtx4D_t;
  float vtx4D_tErr;
  int vtx4D_isFake;
  float vtxGen_t;
  float vtxGen_z;
  vector<int> *muon_isLoose;
  vector<int> *muon_isPrompt;
  vector<int> *muon_isMatchedToGenJet;
  vector<int>   *muon_isFromTauDecay;
  vector<float> *muon_pt;
  vector<float> *muon_eta;
  vector<float> *muon_phi;
  vector<float> *muon_t;
  vector<float> *muon_dz3D;
  vector<float> *muon_dz4D;
  vector<float> *muon_dxy3D;
  vector<float> *muon_dxy4D;
  vector<float> *muon_chIso03_dZ1_simVtx;
  vector<float> *track_pt;
  vector<float> *track_eta;
  vector<float> *track_phi;
  vector<float> *track_dz3D;
  vector<float> *track_dz4D;
  vector<float> *track_dxy3D;
  vector<float> *track_dxy4D;
  vector<float> *track_t;
  vector<int> *track_muIndex;
  
  muon_pt = 0;
  muon_eta = 0;
  muon_phi = 0;
  muon_dz3D = 0;
  muon_dz4D = 0;
  muon_dxy3D = 0;
  muon_dxy4D = 0;
  muon_t = 0;
  muon_isLoose= 0;
  muon_isMatchedToGenJet = 0;
  muon_isPrompt = 0;
  muon_isFromTauDecay = 0;
  muon_chIso03_dZ1_simVtx = 0;
  track_pt = 0;
  track_eta = 0;
  track_phi = 0;
  track_dz3D = 0;
  track_dz4D = 0;
  track_dxy3D = 0;
  track_dxy4D = 0;
  track_t = 0;
  track_muIndex = 0;

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("npu",1);                    chain->SetBranchAddress("npu",            &npu);
  chain->SetBranchStatus("vtx3D_z",1);                chain->SetBranchAddress("vtx3D_z",        &vtx3D_z);
  chain->SetBranchStatus("vtx4D_z",1);                chain->SetBranchAddress("vtx4D_z",        &vtx4D_z);
  chain->SetBranchStatus("vtx4D_t",1);                chain->SetBranchAddress("vtx4D_t",        &vtx4D_t);
  chain->SetBranchStatus("vtx4D_tErr",1);             chain->SetBranchAddress("vtx4D_tErr",     &vtx4D_tErr);
  chain->SetBranchStatus("vtx4D_isFake",1);           chain->SetBranchAddress("vtx4D_isFake",   &vtx4D_isFake);
  chain->SetBranchStatus("vtxGen_z",1);               chain->SetBranchAddress("vtxGen_z",       &vtxGen_z);
  chain->SetBranchStatus("vtxGen_t",1);               chain->SetBranchAddress("vtxGen_t",       &vtxGen_t);
  chain->SetBranchStatus("muon_pt",1);                chain->SetBranchAddress("muon_pt",        &muon_pt);
  chain->SetBranchStatus("muon_eta",1);               chain->SetBranchAddress("muon_eta",       &muon_eta);
  chain->SetBranchStatus("muon_phi",1);               chain->SetBranchAddress("muon_phi",       &muon_phi);
  chain->SetBranchStatus("muon_t",1);                 chain->SetBranchAddress("muon_t",         &muon_t);
  chain->SetBranchStatus("muon_dz3D",1);              chain->SetBranchAddress("muon_dz3D",      &muon_dz3D);
  chain->SetBranchStatus("muon_dz4D",1);              chain->SetBranchAddress("muon_dz4D",      &muon_dz4D);
  chain->SetBranchStatus("muon_dxy3D",1);             chain->SetBranchAddress("muon_dxy3D",     &muon_dxy3D);
  chain->SetBranchStatus("muon_dxy4D",1);             chain->SetBranchAddress("muon_dxy4D",     &muon_dxy4D);
  chain->SetBranchStatus("muon_isLoose",1);           chain->SetBranchAddress("muon_isLoose",   &muon_isLoose);
  chain->SetBranchStatus("muon_isPrompt",1);          chain->SetBranchAddress("muon_isPrompt",  &muon_isPrompt);
  chain->SetBranchStatus("muon_isMatchedToGenJet",1); chain->SetBranchAddress("muon_isMatchedToGenJet", &muon_isMatchedToGenJet);
  chain->SetBranchStatus("muon_isFromTauDecay",1); chain->SetBranchAddress("muon_isFromTauDecay", &muon_isFromTauDecay);
  chain->SetBranchStatus("muon_chIso03_dZ1_simVtx",1);     chain->SetBranchAddress("muon_chIso03_dZ1_simVtx",    &muon_chIso03_dZ1_simVtx);
  
  chain->SetBranchStatus("track_pt",1);               chain->SetBranchAddress("track_pt",       &track_pt);
  chain->SetBranchStatus("track_eta",1);              chain->SetBranchAddress("track_eta",      &track_eta);
  chain->SetBranchStatus("track_phi",1);              chain->SetBranchAddress("track_phi",      &track_phi);
  chain->SetBranchStatus("track_dz3D",1);             chain->SetBranchAddress("track_dz3D",     &track_dz3D);
  chain->SetBranchStatus("track_dz4D",1);             chain->SetBranchAddress("track_dz4D",     &track_dz4D);
  chain->SetBranchStatus("track_dxy3D",1);            chain->SetBranchAddress("track_dxy3D",    &track_dxy3D);
  chain->SetBranchStatus("track_dxy4D",1);            chain->SetBranchAddress("track_dxy4D",    &track_dxy4D);
  chain->SetBranchStatus("track_t",1);                chain->SetBranchAddress("track_t",        &track_t);
  chain->SetBranchStatus("track_muIndex",1);          chain->SetBranchAddress("track_muIndex",  &track_muIndex);

  TH1F *hsimvtx_z = new TH1F("hsimvtx_z","hsimvtx_z", 500,-25,25);
  TH1F *hsimvtx_t = new TH1F("hsimvtx_t","hsimvtx_t", 1000,-1,1);

  TH1F *h_muon_pt_barrel = new TH1F("h_muon_pt_barrel","h_muon_pt_barrel",100,0,100);
  TH1F *h_muon_pt_endcap = new TH1F("h_muon_pt_endcap","h_muon_pt_endcap",100,0,100);

  TH1F *h_muon_dz3D = new TH1F("h_muon_dz3D","h_muon_dz3D",1000,-1.,1.);
  TH1F *h_muon_dz4D = new TH1F("h_muon_dz4D","h_muon_dz4D",1000,-1.,1.);
  TH1F *h_muon_dt4D = new TH1F("h_muon_dt4D","h_muon_dt4D",1000,-1.,1.);
  
  TH1F *h_vtx_dz3D_sim = new TH1F("h_vtx_dz3D_sim","h_vtx_dz3D_sim",1000, -0.2, 0.2);
  TH1F *h_vtx_dz4D_sim = new TH1F("h_vtx_dz4D_sim","h_vtx_dz4D_sim",1000, -0.2, 0.2);
  TH1F *h_vtx_dt4D_sim = new TH1F("h_vtx_dt4D_sim","h_vtx_dt4D_sim",10000, -1.0, 1.0);
  
  TH1F *h_tracks_dt_vtx_barrel = new TH1F("h_tracks_dt_vtx_barrel","h_tracks_dt_vtx_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_dz_vtx_barrel = new TH1F("h_tracks_dz_vtx_barrel","h_tracks_dz_vtx_barrel",500, -1.0, 1.0);
  
  TH1F *h_tracks_dt_mu_barrel = new TH1F("h_tracks_dt_mu_barrel","h_tracks_dt_mu_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_dz_mu_barrel = new TH1F("h_tracks_dz_mu_barrel","h_tracks_dz_mu_barrel",500, -1.0, 1.0);

  TH1F *h_tracks_dt_vtx_endcap = new TH1F("h_tracks_dt_vtx_endcap","h_tracks_dt_vtx_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_dz_vtx_endcap = new TH1F("h_tracks_dz_vtx_endcap","h_tracks_dz_vtx_endcap",500, -1.0, 1.0);
  
  TH1F *h_tracks_dt_mu_endcap = new TH1F("h_tracks_dt_mu_endcap","h_tracks_dt_mu_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_dz_mu_endcap = new TH1F("h_tracks_dz_mu_endcap","h_tracks_dz_mu_endcap",500, -1.0, 1.0);
  
  TH2F *h2_tracks_dzvtx_dtvtx_barrel = new TH2F("h2_tracks_dzvtx_dtvtx_barrel","h2_tracks_dzvtx_dtvtx_barrel",1000,-1,1, 1000, -1, 1);
  TH2F *h2_tracks_dzvtx_dtvtx_endcap = new TH2F("h2_tracks_dzvtx_dtvtx_endcap","h2_tracks_dzvtx_dtvtx_endcap",1000,-1,1, 1000, -1, 1);

  TH2F *h2_tracks_dzmu_dtmu_barrel = new TH2F("h2_tracks_dzmu_dtmu_barrel","h2_tracks_dzmu_dtmu_barrel",1000,-1,1, 1000, -1, 1);
  TH2F *h2_tracks_dzmu_dtmu_endcap = new TH2F("h2_tracks_dzmu_dtmu_endcap","h2_tracks_dzmu_dtmu_endcap",1000,-1,1, 1000, -1, 1);


  // -- to debug tracks in cone 
  TH1F *h_tracks_pt_barrel = new TH1F("h_tracks_pt_barrel","h_tracks_pt_barrel",400,0.,20.);
  TH1F *h_tracks_removed_pt_barrel = new TH1F("h_tracks_removed_pt_barrel","h_tracks_removed_pt_barrel",400,0.,20.);
  TH1F *h_tracks_kept_pt_barrel = new TH1F("h_tracks_kept_pt_barrel","h_tracks_kept_pt_barrel",400,0.,20.);

  TH1F *h_tracks_n_barrel = new TH1F("h_tracks_n_barrel","h_ntracks_pt_barrel",50,0.,50.);
  TH1F *h_tracks_removed_n_barrel = new TH1F("h_tracks_removed_n_barrel","h_tracks_removed_n_barrel",50,0.,50.);
  TH1F *h_tracks_kept_n_barrel = new TH1F("h_tracks_kept_n_barrel","h_tracks_kept_n_barrel",50,0.,50.);

  TH1F *h_tracks_sumpt_barrel = new TH1F("h_tracks_sumpt_barrel","h_tracks_sumpt_barrel",1000,0.,100.);
  TH1F *h_tracks_removed_sumpt_barrel = new TH1F("h_tracks_removed_sumpt_barrel","h_tracks_removed_sumpt_barrel",1000,0.,100.);
  TH1F *h_tracks_kept_sumpt_barrel = new TH1F("h_tracks_kept_sumpt_barrel","h_tracks_kept_sumpt_barrel",1000,0.,100.);

  TProfile *p_tracks_pt_vs_linedensity_barrel = new TProfile("p_tracks_pt_vs_linedensity_barrel","p_tracks_pt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_removed_pt_vs_linedensity_barrel = new TProfile("p_tracks_removed_pt_vs_linedensity_barrel","p_tracks_removed_pt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_kept_pt_vs_linedensity_barrel = new TProfile("p_tracks_kept_pt_vs_linedensity_barrel","p_tracks_kept_pt_vs_linedensity_barrel",22,-0.05,2.05);

  TProfile *p_tracks_n_vs_linedensity_barrel = new TProfile("p_tracks_n_vs_linedensity_barrel","p_tracks_n_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_removed_n_vs_linedensity_barrel = new TProfile("p_tracks_removed_n_vs_linedensity_barrel","p_tracks_removed_n_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_kept_n_vs_linedensity_barrel = new TProfile("p_tracks_kept_n_vs_linedensity_barrel","p_tracks_kept_n_vs_linedensity_barrel",22,-0.05,2.05);
  
  TProfile *p_tracks_sumpt_vs_linedensity_barrel = new TProfile("p_tracks_sumpt_vs_linedensity_barrel","p_tracks_sumpt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_removed_sumpt_vs_linedensity_barrel = new TProfile("p_tracks_removed_sumpt_vs_linedensity_barrel","p_tracks_removed_sumpt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_tracks_kept_sumpt_vs_linedensity_barrel = new TProfile("p_tracks_kept_sumpt_vs_linedensity_barrel","p_tracks_kept_sumpt_vs_linedensity_barrel",22,-0.05,2.05);


  TH1F *h_tracks_pt_endcap = new TH1F("h_tracks_pt_endcap","h_tracks_pt_endcap",400,0.,20.);
  TH1F *h_tracks_removed_pt_endcap = new TH1F("h_tracks_removed_pt_endcap","h_tracks_removed_pt_endcap",400,0.,20.);
  TH1F *h_tracks_kept_pt_endcap = new TH1F("h_tracks_kept_pt_endcap","h_tracks_kept_pt_endcap",400,0.,20.);

  TH1F *h_tracks_n_endcap = new TH1F("h_tracks_n_endcap","h_ntracks_pt_endcap",50,0.,50.);
  TH1F *h_tracks_removed_n_endcap = new TH1F("h_tracks_removed_n_endcap","h_tracks_removed_n_endcap",50,0.,50.);
  TH1F *h_tracks_kept_n_endcap = new TH1F("h_tracks_kept_n_endcap","h_tracks_kept_n_endcap",50,0.,50.);

  TH1F *h_tracks_sumpt_endcap = new TH1F("h_tracks_sumpt_endcap","h_tracks_sumpt_endcap",1000,0.,100.);
  TH1F *h_tracks_removed_sumpt_endcap = new TH1F("h_tracks_removed_sumpt_endcap","h_tracks_removed_sumpt_endcap",1000,0.,100.);
  TH1F *h_tracks_kept_sumpt_endcap = new TH1F("h_tracks_kept_sumpt_endcap","h_tracks_kept_sumpt_endcap",1000,0.,100.);

  TProfile *p_tracks_pt_vs_linedensity_endcap = new TProfile("p_tracks_pt_vs_linedensity_endcap","p_tracks_pt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_removed_pt_vs_linedensity_endcap = new TProfile("p_tracks_removed_pt_vs_linedensity_endcap","p_tracks_removed_pt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_kept_pt_vs_linedensity_endcap = new TProfile("p_tracks_kept_pt_vs_linedensity_endcap","p_tracks_kept_pt_vs_linedensity_endcap",22,-0.05,2.05);

  TProfile *p_tracks_n_vs_linedensity_endcap = new TProfile("p_tracks_n_vs_linedensity_endcap","p_tracks_n_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_removed_n_vs_linedensity_endcap = new TProfile("p_tracks_removed_n_vs_linedensity_endcap","p_tracks_removed_n_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_kept_n_vs_linedensity_endcap = new TProfile("p_tracks_kept_n_vs_linedensity_endcap","p_tracks_kept_n_vs_linedensity_endcap",22,-0.05,2.05);
  
  TProfile *p_tracks_sumpt_vs_linedensity_endcap = new TProfile("p_tracks_sumpt_vs_linedensity_endcap","p_tracks_sumpt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_removed_sumpt_vs_linedensity_endcap = new TProfile("p_tracks_removed_sumpt_vs_linedensity_endcap","p_tracks_removed_sumpt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_tracks_kept_sumpt_vs_linedensity_endcap = new TProfile("p_tracks_kept_sumpt_vs_linedensity_endcap","p_tracks_kept_sumpt_vs_linedensity_endcap",22,-0.05,2.05);


  // chargedIso for rocs
  TH1F *h_muon_relChIso03_dZ1_barrel = new TH1F("h_muon_relChIso03_dZ1_barrel","h_muon_relChIso03_dZ1_barrel",5000, 0, 5); 
  TH1F *h_muon_relChIso03_dZ1_endcap = new TH1F("h_muon_relChIso03_dZ1_endcap","h_muon_relChIso03_dZ1_endcap",5000,0,5); 

  TH1F *h_muon_relChIso03_dZ1_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ1_dT3s_barrel","h_muon_relChIso03_dT3s_dZ1_barrel",5000, 0, 5); 
  TH1F *h_muon_relChIso03_dZ1_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ1_dT3s_endcap","h_muon_relChIso03_dT3s_dZ1_endcap",5000,0,5); 

  TH1F *h_muon_rawChIso03_dZ1_barrel = new TH1F("h_muon_rawChIso03_dZ1_barrel","h_muon_rawChIso03_dZ1_barrel",5000, 0, 500); 
  TH1F *h_muon_rawChIso03_dZ1_endcap = new TH1F("h_muon_rawChIso03_dZ1_endcap","h_muon_rawChIso03_dZ1_endcap",5000,0,500); 

  TH1F *h_muon_rawChIso03_dZ1_dT3s_barrel = new TH1F("h_muon_rawChIso03_dZ1_dT3s_barrel","h_muon_rawChIso03_dT3s_dZ1_barrel",5000, 0, 500); 
  TH1F *h_muon_rawChIso03_dZ1_dT3s_endcap = new TH1F("h_muon_rawChIso03_dZ1_dT3s_endcap","h_muon_rawChIso03_dT3s_dZ1_endcap",5000,0,500); 
  

  // efficiency vs pT for relChIsoCut at 0.1 
  TProfile *p_muonEff_relChIso03_dZ1_vs_pt_barrel = new TProfile("p_muonEff_relChIso03_dZ1_vs_pt_barrel","p_muonEff_relChIso03_dZ1_vs_pt_barrel",100, 0, 100);
  TProfile *p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel = new TProfile("p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel","p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel",100, 0, 100);

  TProfile *p_muonEff_relChIso03_dZ1_vs_pt_endcap = new TProfile("p_muonEff_relChIso03_dZ1_vs_pt_endcap","p_muonEff_relChIso03_dZ1_vs_pt_endcap",100, 0, 100);
  TProfile *p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap = new TProfile("p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap","p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap",100, 0, 100);
  


  cout << "Analyzing " << chain->GetEntries() << "  events" <<endl;
    
  int maxEntries = chain->GetEntries();
  //  maxEntries = 500000;
  for (int ientry = 0; ientry< maxEntries; ientry++) {
    
    chain -> GetEntry(ientry);
    
    if (ientry%1000==0) cout << "Analyzing event " << ientry << "\r" << flush;


    hsimvtx_z  -> Fill(vtxGen_z);    
    hsimvtx_t  -> Fill(vtxGen_t*1000000000.);    
    
    float linedensity = npu*TMath::Gaus(fabs(10.*vtxGen_z), 0, 42., 1);

    for (unsigned int imu = 0; imu < muon_pt->size(); imu++){
      if (muon_pt->at(imu) < minMuonPt ) continue;
      if (muon_pt->at(imu) > maxMuonPt ) continue;
      if (!muon_isLoose->at(imu) ) continue;
      if (fabs(muon_dz4D->at(imu)) > 0.5 ) continue;
      if (fabs(muon_dxy4D->at(imu)) > 0.2 ) continue;
      if (vtx4D_isFake) continue;

      bool isBarrel = fabs(muon_eta->at(imu)) < 1.5;

      bool pass = false;
      if ( prompt) pass = muon_isPrompt->at(imu);
      if (!prompt) pass = !muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu) && !muon_isFromTauDecay->at(imu);
      
      if (!pass) continue;

      h_vtx_dz3D_sim -> Fill( vtx3D_z - vtxGen_z);
      h_vtx_dz4D_sim -> Fill( vtx4D_z - vtxGen_z);
      h_vtx_dt4D_sim -> Fill( vtx4D_t - vtxGen_t*1000000000.);

      if ( isBarrel) h_muon_pt_barrel   -> Fill(muon_pt->at(imu));
      if ( !isBarrel) h_muon_pt_endcap   -> Fill(muon_pt->at(imu));
      h_muon_dz3D -> Fill(muon_dz3D->at(imu));
      h_muon_dz4D -> Fill(muon_dz4D->at(imu));
      h_muon_dt4D -> Fill(muon_t->at(imu) -  vtx4D_t);
       

      int ntracks = 0;
      int ntracks_removed = 0;
      int ntracks_kept = 0;
      float sumpt = 0.;
      float sumpt_removed = 0.;
      float sumpt_kept = 0.;

      // start loop over tracks in isolation cone
      for (unsigned int itk = 0; itk < track_pt->size(); itk++){
	
	// -- exclude tracks corresponding to the isolation cone of the other muon(s)
	if (track_muIndex ->at(itk) != int(imu)) continue;
	
	if (isnan(track_t->at(itk))) cout << " track time = " << track_t->at(itk) <<endl;

	float deta = fabs(muon_eta->at(imu) - track_eta->at(itk));
	float dphi = fabs(muon_phi->at(imu) - track_phi->at(itk));
	if (dphi > TMath::Pi()) dphi = 2*TMath::Pi()- dphi;
	if (dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
	float dr = sqrt(deta*deta+dphi*dphi);
	
	if ( dr > 0.3 ) continue;
	
	if ( fabs(track_eta->at(itk)) < 1.48 && track_pt->at(itk) < btlMinTrackPt ) continue;
	if ( fabs(track_eta->at(itk)) > 1.48 && track_pt->at(itk) < etlMinTrackPt ) continue;

	float dz = track_dz4D->at(itk);
	float dxy = track_dxy4D->at(itk);
	float dt = track_t->at(itk) - vtx4D_t;
	float dzsim = track_dz4D->at(itk) + vtx4D_z - vtxGen_z;
	float dtsim = track_t->at(itk) - vtxGen_t*1000000000.;
	
	dz = dzsim;
	dt = dtsim;

	// ---only tracks within dz < 0.1 cm from PV
	if (fabs(dz) < maxDz && fabs(dxy) < maxDxy) {

	  // -- barrel
          if (isBarrel){
	    h_tracks_dz_vtx_barrel -> Fill(dz);
	    if (track_t->at(itk) != -999){
	      h_tracks_dt_vtx_barrel->Fill(dt);
	      h2_tracks_dzvtx_dtvtx_barrel->Fill(dz, dt);
	    }
	    
	    // -- all tracks
	    h_tracks_pt_barrel -> Fill( track_pt->at(itk) );
	    p_tracks_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
	    ntracks++;
	    sumpt+=track_pt->at(itk);
	    // -- removed tracks
	    if ( track_t->at(itk) != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
	      h_tracks_removed_pt_barrel -> Fill( track_pt->at(itk) );
	      p_tracks_removed_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
	      ntracks_removed++;
	      sumpt_removed+=track_pt->at(itk);
	    }
	    else{
	      h_tracks_kept_pt_barrel -> Fill( track_pt->at(itk) );
	      p_tracks_kept_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
	      ntracks_kept++;
	      sumpt_kept+=track_pt->at(itk);
	    }
	  }// end if barrel
	  // -- endcap
          else{
	    h_tracks_dz_vtx_endcap -> Fill(dz);
            if (track_t->at(itk) != -999){
              h_tracks_dt_vtx_endcap->Fill(dt);
              h2_tracks_dzvtx_dtvtx_endcap->Fill(dz, dt);
            }
            // -- all tracks
            h_tracks_pt_endcap -> Fill( track_pt->at(itk) );
            p_tracks_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
            ntracks++;
            sumpt+=track_pt->at(itk);
            // -- removed tracks
            if ( track_t->at(itk) != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
              h_tracks_removed_pt_endcap -> Fill( track_pt->at(itk) );
              p_tracks_removed_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
              ntracks_removed++;
              sumpt_removed+=track_pt->at(itk);
            }
	    else{
	      h_tracks_kept_pt_endcap -> Fill( track_pt->at(itk) );
              p_tracks_kept_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
              ntracks_kept++;
              sumpt_kept+=track_pt->at(itk);
	    }
	  }// end endcap
	}


	// --- dz wrt to muon
	float dtmu = track_t->at(itk) - muon_t->at(imu);
	float dzmu = muon_dz4D->at(imu) - track_dz4D->at(itk) ;
	if (fabs(dzmu) < maxDzMu){
	  // -- barrel
	  if (isBarrel){
	    h_tracks_dz_mu_barrel -> Fill(dzmu);
	    // if timing info available
	    if (track_t->at(itk) != -999){	
	      h_tracks_dt_mu_barrel->Fill(dtmu);
	      h2_tracks_dzmu_dtmu_barrel->Fill(dzmu, dtmu);
	    }
	  }
	  // -- endcap
	  else {
	    h_tracks_dz_mu_endcap-> Fill(dzmu);
	    // if timing info available 
	    if (track_t->at(itk) != -999){
              h_tracks_dt_mu_endcap->Fill(dtmu);
              h2_tracks_dzmu_dtmu_endcap->Fill(dzmu, dtmu);
	    }
	  }
	}
	
      }// end loop over tracks

      float relChIso = sumpt/muon_pt->at(imu);
      float relChIso_dT = sumpt_kept/muon_pt->at(imu);
      
      if (isBarrel){
	// -- all tracks in isolation cone
	h_tracks_n_barrel -> Fill( ntracks );
	h_tracks_sumpt_barrel -> Fill( sumpt );
	p_tracks_n_vs_linedensity_barrel -> Fill( linedensity, ntracks );
	p_tracks_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt );

	// -- tracks removed from isolation cone
	h_tracks_removed_n_barrel -> Fill( ntracks_removed );
	h_tracks_removed_sumpt_barrel -> Fill( sumpt_removed );
	p_tracks_removed_n_vs_linedensity_barrel -> Fill( linedensity, ntracks_removed );
	p_tracks_removed_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_removed );

	// -- tracks kept from isolation cone
	h_tracks_kept_n_barrel -> Fill( ntracks_kept );
	h_tracks_kept_sumpt_barrel -> Fill( sumpt_kept );
	p_tracks_kept_n_vs_linedensity_barrel -> Fill( linedensity, ntracks_kept );
	p_tracks_kept_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_kept );

	// -- rel chIso
	h_muon_relChIso03_dZ1_barrel ->Fill(sumpt/muon_pt->at(imu));
	h_muon_relChIso03_dZ1_dT3s_barrel ->Fill(sumpt_kept/muon_pt->at(imu));

	// -- raw chIso
	h_muon_rawChIso03_dZ1_barrel ->Fill(sumpt);
	h_muon_rawChIso03_dZ1_dT3s_barrel ->Fill(sumpt_kept);

	// -- eff 
	if (relChIso < 0.05) p_muonEff_relChIso03_dZ1_vs_pt_barrel->Fill(muon_pt->at(imu), 1.);
	else p_muonEff_relChIso03_dZ1_vs_pt_barrel->Fill(muon_pt->at(imu), 0.);

	if (relChIso_dT < 0.05) p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel->Fill(muon_pt->at(imu), 1.);
	else p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel->Fill(muon_pt->at(imu), 0.);
      }
      else{
	h_tracks_n_endcap -> Fill( ntracks );
	h_tracks_sumpt_endcap -> Fill( sumpt );
	p_tracks_n_vs_linedensity_endcap -> Fill( linedensity, ntracks );
	p_tracks_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt );

	h_tracks_removed_n_endcap -> Fill( ntracks_removed );
	h_tracks_removed_sumpt_endcap -> Fill( sumpt_removed );
	p_tracks_removed_n_vs_linedensity_endcap -> Fill( linedensity, ntracks_removed );
	p_tracks_removed_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_removed );

	h_tracks_kept_n_endcap -> Fill( ntracks_kept );
	h_tracks_kept_sumpt_endcap -> Fill( sumpt_kept );
	p_tracks_kept_n_vs_linedensity_endcap -> Fill( linedensity, ntracks_kept );
	p_tracks_kept_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_kept );

	// -- rel chIso
	h_muon_relChIso03_dZ1_endcap ->Fill(sumpt/muon_pt->at(imu));
	h_muon_relChIso03_dZ1_dT3s_endcap ->Fill(sumpt_kept/muon_pt->at(imu));

	// -- raw chIso
	h_muon_rawChIso03_dZ1_endcap ->Fill(sumpt);
	h_muon_rawChIso03_dZ1_dT3s_endcap ->Fill(sumpt_kept);

	// -- eff 
	if (relChIso < 0.1) p_muonEff_relChIso03_dZ1_vs_pt_endcap->Fill(muon_pt->at(imu), 1.);
	else p_muonEff_relChIso03_dZ1_vs_pt_endcap->Fill(muon_pt->at(imu), 0.);

	if (relChIso_dT < 0.1) p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap->Fill(muon_pt->at(imu), 1.);
	else p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap->Fill(muon_pt->at(imu), 0.);
      }

    }// end loop over muons

  }//end loop over events
  
  // -- dt track in muon iso cone
  TF1 *fitfun_barrel = new TF1("fitfun_barrel","gaus(0)+gaus(3)", -10, 10);
  fitfun_barrel->SetParameter(1, 0.);
  fitfun_barrel->SetParameter(2, 0.030);
  fitfun_barrel->SetParameter(4, 0.);
  fitfun_barrel->SetParameter(5, 0.270);
  fitfun_barrel->SetNpx(1000);
  h_tracks_dt_vtx_barrel->Scale(1./h_muon_pt_barrel->GetSumOfWeights());
  h_tracks_dt_vtx_barrel->Fit("fitfun_barrel","QR");
  
  TF1 *fitfun_endcap = new TF1("fitfun_endcap","gaus(0)+gaus(3)", -10, 10);
  fitfun_endcap->SetParameter(1, 0.);
  fitfun_endcap->SetParameter(2, 0.030);
  fitfun_endcap->SetParameter(4, 0.);
  fitfun_endcap->SetParameter(5, 0.270);
  fitfun_endcap->SetNpx(1000);
  h_tracks_dt_vtx_endcap->Scale(1./h_muon_pt_endcap->GetSumOfWeights());
  h_tracks_dt_vtx_endcap->Fit("fitfun_endcap","QR");
  
  
  TF1 *fitfun2_barrel = new TF1("fitfun2_barrel","gaus(0)+gaus(3)", -10, 10);
  fitfun2_barrel->SetParameter(1, 0.);
  fitfun2_barrel->SetParameter(2, 0.042);
  fitfun2_barrel->SetParameter(4, 0.);
  fitfun2_barrel->SetParameter(5, 0.270);
  fitfun2_barrel->SetNpx(1000);
  h_tracks_dt_mu_barrel->Scale(1./h_muon_pt_barrel->GetSumOfWeights());
  h_tracks_dt_mu_barrel->Fit("fitfun2_barrel","QR");
  
  TF1 *fitfun2_endcap = new TF1("fitfun2_endcap","gaus(0)+gaus(3)", -10, 10);
  fitfun2_endcap->SetParameter(1, 0.);
  fitfun2_endcap->SetParameter(2, 0.042);
  fitfun2_endcap->SetParameter(4, 0.);
  fitfun2_endcap->SetParameter(5, 0.270);
  fitfun2_endcap->SetNpx(1000);
  h_tracks_dt_mu_endcap->Scale(1./h_muon_pt_endcap->GetSumOfWeights());
  h_tracks_dt_mu_endcap->Fit("fitfun2_endcap","QR");


  
  // -- save histograms in output file
  std::string foutName = "testTracks_" + process + "_" + pu + "_dT"+ std::to_string(nsigma)+"sigma"+
                         "_minMuonPt"+std::to_string(int(minMuonPt))+
                         "_maxMuonPt"+std::to_string(int(maxMuonPt))+
                         "_minTrackPt.root";

  TFile *fout = new TFile(foutName.c_str(),"recreate");

  hsimvtx_z->Write();
  hsimvtx_t->Write();

  h_vtx_dz3D_sim->Write();
  h_vtx_dz4D_sim->Write();
  h_vtx_dt4D_sim->Write();
  h_muon_pt_barrel->Write();
  h_muon_pt_endcap->Write();
  h_muon_dz3D->Write();
  h_muon_dz4D->Write();
  h_muon_dt4D->Write();
  
  h_tracks_dz_vtx_barrel->Write();
  h_tracks_dz_vtx_endcap->Write();

  h_tracks_dt_vtx_barrel->Write();
  h_tracks_dt_vtx_endcap->Write();
  
  h_tracks_dz_mu_barrel->Write();
  h_tracks_dt_mu_barrel->Write();
  
  h_tracks_dz_mu_endcap->Write();
  h_tracks_dt_mu_endcap->Write();
  
  h2_tracks_dzvtx_dtvtx_barrel-> Write();
  h2_tracks_dzvtx_dtvtx_endcap-> Write();

  h2_tracks_dzmu_dtmu_barrel->Write();  
  h2_tracks_dzmu_dtmu_endcap->Write();



  h_tracks_pt_barrel ->Write();
  h_tracks_removed_pt_barrel ->Write();
  h_tracks_kept_pt_barrel ->Write();

  h_tracks_n_barrel ->Write();
  h_tracks_removed_n_barrel ->Write();
  h_tracks_kept_n_barrel ->Write();

  h_tracks_sumpt_barrel ->Write();
  h_tracks_removed_sumpt_barrel ->Write();
  h_tracks_kept_sumpt_barrel ->Write();

  p_tracks_pt_vs_linedensity_barrel -> Write();
  p_tracks_removed_pt_vs_linedensity_barrel -> Write();
  p_tracks_kept_pt_vs_linedensity_barrel -> Write();

  p_tracks_n_vs_linedensity_barrel -> Write();
  p_tracks_removed_n_vs_linedensity_barrel -> Write();
  p_tracks_kept_n_vs_linedensity_barrel -> Write();

  p_tracks_sumpt_vs_linedensity_barrel -> Write();
  p_tracks_removed_sumpt_vs_linedensity_barrel -> Write();
  p_tracks_kept_sumpt_vs_linedensity_barrel -> Write();


  h_tracks_pt_endcap ->Write();
  h_tracks_removed_pt_endcap ->Write();
  h_tracks_kept_pt_endcap ->Write();

  h_tracks_n_endcap ->Write();
  h_tracks_removed_n_endcap ->Write();
  h_tracks_kept_n_endcap ->Write();

  h_tracks_sumpt_endcap ->Write();
  h_tracks_removed_sumpt_endcap ->Write();
  h_tracks_kept_sumpt_endcap ->Write();

  p_tracks_pt_vs_linedensity_endcap -> Write();
  p_tracks_removed_pt_vs_linedensity_endcap -> Write();
  p_tracks_kept_pt_vs_linedensity_endcap -> Write();

  p_tracks_n_vs_linedensity_endcap -> Write();
  p_tracks_removed_n_vs_linedensity_endcap -> Write();
  p_tracks_kept_n_vs_linedensity_endcap -> Write();

  p_tracks_sumpt_vs_linedensity_endcap -> Write();
  p_tracks_removed_sumpt_vs_linedensity_endcap -> Write();
  p_tracks_kept_sumpt_vs_linedensity_endcap -> Write();

  h_muon_relChIso03_dZ1_barrel -> Write();
  h_muon_relChIso03_dZ1_dT3s_barrel -> Write();

  h_muon_relChIso03_dZ1_endcap -> Write();
  h_muon_relChIso03_dZ1_dT3s_endcap -> Write();


  h_muon_rawChIso03_dZ1_barrel -> Write();
  h_muon_rawChIso03_dZ1_dT3s_barrel -> Write();

  h_muon_rawChIso03_dZ1_endcap -> Write();
  h_muon_rawChIso03_dZ1_dT3s_endcap -> Write();

  p_muonEff_relChIso03_dZ1_vs_pt_barrel->Write();
  p_muonEff_relChIso03_dZ1_dT3s_vs_pt_barrel->Write();

  p_muonEff_relChIso03_dZ1_vs_pt_endcap->Write();
  p_muonEff_relChIso03_dZ1_dT3s_vs_pt_endcap->Write();


  fout->Close();
  
}



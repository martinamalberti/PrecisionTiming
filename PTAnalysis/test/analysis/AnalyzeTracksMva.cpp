// per compilare: g++ -Wall -o AnalyzeTracksMva `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit AnalyzeTracksMva.cpp                        

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

  float timeResolution = 0.040;
  int nsigma = 3;

  float minMuonPt = 20.;
  float maxMuonPt = 9999.;
  float maxDzBtl = 0.2;
  float maxDzEtl = 0.3;
  float maxDxy = 9999.;
  float btlMinTrackPt = 0.7;
  float etlMinTrackPt = 0.4;
  float btlEfficiency = 0.85;
  float etlEfficiency = 0.90;
 
  string process = argv[1];
  bool prompt = false;

  string pu = argv[2];
  
  // -- get TChain
  TChain* chain = new TChain("MuonIsolationAnalyzer/tree_35ps","tree");
  if (process.find("DYToLL") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/DYToLL_M-50_14TeV_TuneCP5_pythia8/test_DYToLL_muIso/190320_223157/0000/muonIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/DYToLL_M-50_14TeV_TuneCP5_pythia8/test_DYToLL_noPU_muIso/190321_100924/0000/muonIsolation*.root");
    prompt = true;
  }
  
  if (process.find("TTbar") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/TTbar_14TeV_TuneCP5_Pythia8/test_TTbar_muIso/190320_223713/0000/muonIsolation*.root");
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/TTbar_14TeV_TuneCP5_Pythia8/test_TTbar_muIso/190320_223713/0001/muonIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/RelValTTbar_Tauola_14TeV/test_TTbar_noPU_muIso/190321_080855/0000/muonIsolation*.root");
    prompt = false;
  }
  
  if (process.find("QCD") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/test_QCD_muIso/190320_224211/0000/muonIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/1040mtd5/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/test_QCD_muIso/190320_224211/0000/muonIsolation*.root");
    prompt = false; 
    maxDzBtl = 0.1;
    maxDzEtl = 0.1;
  }
  
  cout << "Using prompt muons = " << prompt <<endl;
    
  // -- tree vars
  int npu;
  float vtx3D_z;
  float vtx4D_z;
  float vtx4D_t;
  float vtx4D_tErr;
  int vtx3D_isFake;
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
  vector<float> *track_pt;
  vector<float> *track_eta;
  vector<float> *track_phi;
  vector<float> *track_dz3D;
  vector<float> *track_dz4D;
  vector<float> *track_dxy3D;
  vector<float> *track_dxy4D;
  vector<float> *track_t;
  vector<float> *track_tFastSim;
  vector<int> *track_muIndex;
  vector<int> *track_isMatchedToGenParticle;
  vector<float> *track_puid3dmva;
  vector<float> *track_puid4dmva;
  
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
  track_pt = 0;
  track_eta = 0;
  track_phi = 0;
  track_dz3D = 0;
  track_dz4D = 0;
  track_dxy3D = 0;
  track_dxy4D = 0;
  track_t = 0;
  track_tFastSim = 0;
  track_muIndex = 0;
  track_isMatchedToGenParticle = 0;
  track_puid3dmva = 0;
  track_puid4dmva = 0;

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("npu",1);                    chain->SetBranchAddress("npu",            &npu);
  chain->SetBranchStatus("vtx3D_z",1);                chain->SetBranchAddress("vtx3D_z",        &vtx3D_z);
  chain->SetBranchStatus("vtx4D_z",1);                chain->SetBranchAddress("vtx4D_z",        &vtx4D_z);
  chain->SetBranchStatus("vtx4D_t",1);                chain->SetBranchAddress("vtx4D_t",        &vtx4D_t);
  chain->SetBranchStatus("vtx4D_tErr",1);             chain->SetBranchAddress("vtx4D_tErr",     &vtx4D_tErr);
  chain->SetBranchStatus("vtx4D_isFake",1);           chain->SetBranchAddress("vtx4D_isFake",   &vtx4D_isFake);
  chain->SetBranchStatus("vtx3D_isFake",1);           chain->SetBranchAddress("vtx3D_isFake",   &vtx3D_isFake);
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
  
  chain->SetBranchStatus("track_pt",1);               chain->SetBranchAddress("track_pt",       &track_pt);
  chain->SetBranchStatus("track_eta",1);              chain->SetBranchAddress("track_eta",      &track_eta);
  chain->SetBranchStatus("track_phi",1);              chain->SetBranchAddress("track_phi",      &track_phi);
  chain->SetBranchStatus("track_dz3D",1);             chain->SetBranchAddress("track_dz3D",     &track_dz3D);
  chain->SetBranchStatus("track_dz4D",1);             chain->SetBranchAddress("track_dz4D",     &track_dz4D);
  chain->SetBranchStatus("track_dxy3D",1);            chain->SetBranchAddress("track_dxy3D",    &track_dxy3D);
  chain->SetBranchStatus("track_dxy4D",1);            chain->SetBranchAddress("track_dxy4D",    &track_dxy4D);
  chain->SetBranchStatus("track_t",1);                chain->SetBranchAddress("track_t",        &track_t);
  chain->SetBranchStatus("track_tFastSim",1);         chain->SetBranchAddress("track_tFastSim", &track_tFastSim);
  chain->SetBranchStatus("track_muIndex",1);          chain->SetBranchAddress("track_muIndex",  &track_muIndex);
  chain->SetBranchStatus("track_isMatchedToGenParticle",1);  chain->SetBranchAddress("track_isMatchedToGenParticle",  &track_isMatchedToGenParticle);

  chain->SetBranchStatus("track_puid3dmva",1);               chain->SetBranchAddress("track_puid3dmva",       &track_puid3dmva);
  chain->SetBranchStatus("track_puid4dmva",1);               chain->SetBranchAddress("track_puid4dmva",       &track_puid4dmva);


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
  TH1F *h_tracks_dz_vtx_barrel = new TH1F("h_tracks_dz_vtx_barrel","h_tracks_dz_vtx_barrel",2000, -1.0, 1.0);

  TH1F *h_tracks_dt_vtx_endcap = new TH1F("h_tracks_dt_vtx_endcap","h_tracks_dt_vtx_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_dz_vtx_endcap = new TH1F("h_tracks_dz_vtx_endcap","h_tracks_dz_vtx_endcap",2000, -1.0, 1.0);
  
  TH1F *h_tracks_genmatched_dt_vtx_barrel = new TH1F("h_tracks_genmatched_dt_vtx_barrel","h_tracks_genmatched_dt_vtx_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_genmatched_dz_vtx_barrel = new TH1F("h_tracks_genmatched_dz_vtx_barrel","h_tracks_genmatched_dz_vtx_barrel",2000, -1.0, 1.0);

  TH1F *h_tracks_genmatched_dt_vtx_endcap = new TH1F("h_tracks_genmatched_dt_vtx_endcap","h_tracks_genmatched_dt_vtx_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_genmatched_dz_vtx_endcap = new TH1F("h_tracks_genmatched_dz_vtx_endcap","h_tracks_genmatched_dz_vtx_endcap",2000, -1.0, 1.0);


  TH2F *h2_tracks_dzvtx_dtvtx_barrel = new TH2F("h2_tracks_dzvtx_dtvtx_barrel","h2_tracks_dzvtx_dtvtx_barrel",1000,-1,1, 1000, -1, 1);
  TH2F *h2_tracks_dzvtx_dtvtx_endcap = new TH2F("h2_tracks_dzvtx_dtvtx_endcap","h2_tracks_dzvtx_dtvtx_endcap",1000,-1,1, 1000, -1, 1);

  TH1F *h_tracks_puid3dmva_barrel = new TH1F("h_tracks_puid3dmva_barrel","h_tracks_puid3dmva_barrel", 100, -1, 1);
  TH1F *h_tracks_puid4dmva_barrel = new TH1F("h_tracks_puid4dmva_barrel","h_tracks_puid4dmva_barrel", 100, -1, 1);

  TH1F *h_tracks_puid3dmva_endcap = new TH1F("h_tracks_puid3dmva_endcap","h_tracks_puid3dmva_endcap", 100, -1, 1);
  TH1F *h_tracks_puid4dmva_endcap = new TH1F("h_tracks_puid4dmva_endcap","h_tracks_puid4dmva_endcap", 100, -1, 1);

  TH1F *h_tracks_genmatched_puid3dmva_barrel = new TH1F("h_tracks_genmatched_puid3dmva_barrel","h_tracks_genmatched_puid3dmva_barrel", 100, -1, 1);
  TH1F *h_tracks_genmatched_puid4dmva_barrel = new TH1F("h_tracks_genmatched_puid4dmva_barrel","h_tracks_genmatched_puid4dmva_barrel", 100, -1, 1);

  TH1F *h_tracks_genmatched_puid3dmva_endcap = new TH1F("h_tracks_genmatched_puid3dmva_endcap","h_tracks_genmatched_puid3dmva_endcap", 100, -1, 1);
  TH1F *h_tracks_genmatched_puid4dmva_endcap = new TH1F("h_tracks_genmatched_puid4dmva_endcap","h_tracks_genmatched_puid4dmva_endcap", 100, -1, 1);


  TH1F *h_tracks_pass3dmva_dt_vtx_fullSim_barrel = new TH1F("h_tracks_pass3dmva_dt_vtx_fullSim_barrel","h_tracks_pass3dmva_dt_vtx_fullSim_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_pass3dmva_dz_vtx_fullSim_barrel = new TH1F("h_tracks_pass3dmva_dz_vtx_fullSim_barrel","h_tracks_pass3dmva_dz_vtx_fullSim_barrel",500, -1.0, 1.0);

  TH1F *h_tracks_pass3dmva_dt_vtx_fullSim_endcap = new TH1F("h_tracks_pass3dmva_dt_vtx_fullSim_endcap","h_tracks_pass3dmva_dt_vtx_fullSim_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_pass3dmva_dz_vtx_fullSim_endcap = new TH1F("h_tracks_pass3dmva_dz_vtx_fullSim_endcap","h_tracks_pass3dmva_dz_vtx_fullSim_endcap",500, -1.0, 1.0);

  TH1F *h_tracks_pass4dmva_dt_vtx_fullSim_barrel = new TH1F("h_tracks_pass4dmva_dt_vtx_fullSim_barrel","h_tracks_pass4dmva_dt_vtx_fullSim_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_pass4dmva_dz_vtx_fullSim_barrel = new TH1F("h_tracks_pass4dmva_dz_vtx_fullSim_barrel","h_tracks_pass4dmva_dz_vtx_fullSim_barrel",500, -1.0, 1.0);

  TH1F *h_tracks_pass4dmva_dt_vtx_fullSim_endcap = new TH1F("h_tracks_pass4dmva_dt_vtx_fullSim_endcap","h_tracks_pass4dmva_dt_vtx_fullSim_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_pass4dmva_dz_vtx_fullSim_endcap = new TH1F("h_tracks_pass4dmva_dz_vtx_fullSim_endcap","h_tracks_pass4dmva_dz_vtx_fullSim_endcap",500, -1.0, 1.0);


  TH1F *h_tracks_genmatched_dt_vtx_fullSim_barrel = new TH1F("h_tracks_genmatched_dt_vtx_fullSim_barrel","h_tracks_genmatched_dt_vtx_fullSim_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_genmatched_dz_vtx_fullSim_barrel = new TH1F("h_tracks_genmatched_dz_vtx_fullSim_barrel","h_tracks_genmatched_dz_vtx_fullSim_barrel",500, -1.0, 1.0);

  TH1F *h_tracks_genmatched_dt_vtx_fullSim_endcap = new TH1F("h_tracks_genmatched_dt_vtx_fullSim_endcap","h_tracks_genmatched_dt_vtx_fullSim_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_genmatched_dz_vtx_fullSim_endcap = new TH1F("h_tracks_genmatched_dz_vtx_fullSim_endcap","h_tracks_genmatched_dz_vtx_fullSim_endcap",500, -1.0, 1.0);


  TH1F *h_tracks_pass3dmva_pt_barrel = new TH1F("h_tracks_pass3dmva_pt_barrel","h_tracks_pass3dmva_pt_barrel",400, 0, 20);
  TH1F *h_tracks_pass3dmva_pt_endcap = new TH1F("h_tracks_pass3dmva_pt_endcap","h_tracks_pass3dmva_pt_endcap",400, 0, 20);

  TH1F *h_tracks_pass4dmva_pt_barrel = new TH1F("h_tracks_pass4dmva_pt_barrel","h_tracks_pass4dmva_pt_barrel",400, 0, 20);
  TH1F *h_tracks_pass4dmva_pt_endcap = new TH1F("h_tracks_pass4dmva_pt_endcap","h_tracks_pass4dmva_pt_endcap",400, 0, 20);

  TH1F *h_tracks_pass3dmva_eta_barrel = new TH1F("h_tracks_pass3dmva_eta_barrel","h_tracks_pass3dmva_eta_barrel",100, -4, 4);
  TH1F *h_tracks_pass3dmva_eta_endcap = new TH1F("h_tracks_pass3dmva_eta_endcap","h_tracks_pass3dmva_eta_endcap",100, -4, 4);
  TH1F *h_tracks_pass4dmva_eta_barrel = new TH1F("h_tracks_pass4dmva_eta_barrel","h_tracks_pass4dmva_eta_barrel",100, -4, 4);
  TH1F *h_tracks_pass4dmva_eta_endcap = new TH1F("h_tracks_pass4dmva_eta_endcap","h_tracks_pass4dmva_eta_endcap",100, -4, 4);

  


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



  TH1F *h_tracks_pt_withtiming_barrel = new TH1F("h_tracks_pt_withtiming_barrel","h_tracks_withtiming_pt_barrel",400,0.,20.);
  TH1F *h_tracks_pt_withtiming_endcap = new TH1F("h_tracks_pt_withtiming_endcap","h_tracks_withtiming_pt_endcap",400,0.,20.);


  // chargedIso for rocs
  const int NCUTSDZ = 4;
  float dzCut[NCUTSDZ] = {0.1, 0.2, 0.3, 1.0};
  TH1F *h_muon_relChIso03_dZ_barrel[NCUTSDZ];
  TH1F *h_muon_relChIso03_dZ_endcap[NCUTSDZ];
  TH1F *h_muon_relChIso03_dZ_dT3s_barrel[NCUTSDZ];
  TH1F *h_muon_relChIso03_dZ_dT3s_endcap[NCUTSDZ];
  for (int i = 0; i < NCUTSDZ; i++){
    h_muon_relChIso03_dZ_barrel[i] = new TH1F(Form("h_muon_relChIso03_dZ%.00f_barrel", dzCut[i]*10) , Form("h_muon_relChIso03_dZ%.00f_barrel", dzCut[i]*10) ,5000, 0, 5);
    h_muon_relChIso03_dZ_endcap[i] = new TH1F(Form("h_muon_relChIso03_dZ%.00f_endcap", dzCut[i]*10) , Form("h_muon_relChIso03_dZ%.00f_endcap", dzCut[i]*10) ,5000, 0, 5);
    h_muon_relChIso03_dZ_dT3s_barrel[i] = new TH1F(Form("h_muon_relChIso03_dZ%.00f_dT3s_barrel", dzCut[i]*10) , Form("h_muon_relChIso03_dZ%.00f_dT3s_barrel", dzCut[i]*10) ,5000, 0, 5);
    h_muon_relChIso03_dZ_dT3s_endcap[i] = new TH1F(Form("h_muon_relChIso03_dZ%.00f_dT3s_endcap", dzCut[i]*10) , Form("h_muon_relChIso03_dZ%.00f_dT3s_endcap", dzCut[i]*10) ,5000, 0, 5);
  }

  const int NCUTSMVA = 13;
  TH1F *h_muon_relChIso03_mva3D_barrel[NCUTSMVA];
  TH1F *h_muon_relChIso03_mva3D_endcap[NCUTSMVA];
  TH1F *h_muon_relChIso03_mva4D_barrel[NCUTSMVA];
  TH1F *h_muon_relChIso03_mva4D_endcap[NCUTSMVA];
  float mvaCut3D[NCUTSMVA] = {0.7228171, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.80, 0.85, 0.90};
  float mvaCut4D[NCUTSMVA] = {0.7534892, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.80, 0.85, 0.90};
  for (int i = 0; i < NCUTSMVA; i++){
    h_muon_relChIso03_mva3D_barrel[i] = new TH1F(Form("h_muon_relChIso03_mva3D_%.03f_barrel", mvaCut3D[i]) , Form("h_muon_relChIso03_mva3D_%.03f_barrel", mvaCut3D[i]) ,5000, 0, 5);
    h_muon_relChIso03_mva3D_endcap[i] = new TH1F(Form("h_muon_relChIso03_mva3D_%.03f_endcap", mvaCut3D[i]) , Form("h_muon_relChIso03_mva3D_%.03f_endcap", mvaCut3D[i]) ,5000, 0, 5);
    h_muon_relChIso03_mva4D_barrel[i] = new TH1F(Form("h_muon_relChIso03_mva4D_%.03f_barrel", mvaCut4D[i]) , Form("h_muon_relChIso03_mva4D_%.03f_barrel", mvaCut4D[i]) ,5000, 0, 5);
    h_muon_relChIso03_mva4D_endcap[i] = new TH1F(Form("h_muon_relChIso03_mva4D_%.03f_endcap", mvaCut4D[i]) , Form("h_muon_relChIso03_mva4D_%.03f_endcap", mvaCut4D[i]) ,5000, 0, 5);    
  }

  TH1F *h_muon_relChIso03_mva3D_weighted_barrel = new TH1F("h_muon_relChIso03_mva3D_weighted_barrel","h_muon_relChIso03_mva3D_weighted_barrel",5000, 0, 5); 
  TH1F *h_muon_relChIso03_mva3D_weighted_endcap = new TH1F("h_muon_relChIso03_mva3D_weighted_endcap","h_muon_relChIso03_mva3D_weighted_endcap",5000,0,5); 

  TH1F *h_muon_relChIso03_mva4D_weighted_barrel = new TH1F("h_muon_relChIso03_mva4D_weighted_barrel","h_muon_relChIso03_mva4D_weighted_barrel",5000, 0, 5); 
  TH1F *h_muon_relChIso03_mva4D_weighted_endcap = new TH1F("h_muon_relChIso03_mva4D_weighted_endcap","h_muon_relChIso03_mva4D_weighted_endcap",5000,0,5); 




  cout << "Analyzing " << chain->GetEntries() << "  events" <<endl;
    


  TRandom *gRandom  = new TRandom();
  TRandom *gRandom2  = new TRandom();

  int maxEntries = chain->GetEntries();
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
      if (vtx4D_isFake || vtx3D_isFake) continue;
      
      bool pass = false;
      if ( prompt) pass = muon_isPrompt->at(imu);
      if (!prompt) pass = !muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu) && !muon_isFromTauDecay->at(imu);
      
      if (!pass) continue;

      h_vtx_dz3D_sim -> Fill( vtx3D_z - vtxGen_z);
      h_vtx_dz4D_sim -> Fill( vtx4D_z - vtxGen_z);
      h_vtx_dt4D_sim -> Fill( vtx4D_t - vtxGen_t*1000000000.);


      // use only evenst with vtx0 matched to gen vtx
      if (fabs(vtx4D_z - vtxGen_z) > 0.01) continue;
      if (fabs(vtx3D_z - vtxGen_z) > 0.01) continue;


      bool isBarrel = fabs(muon_eta->at(imu)) < 1.5;
      if ( isBarrel) h_muon_pt_barrel   -> Fill(muon_pt->at(imu));
      if ( !isBarrel) h_muon_pt_endcap   -> Fill(muon_pt->at(imu));
      h_muon_dz3D -> Fill(muon_dz3D->at(imu));
      h_muon_dz4D -> Fill(muon_dz4D->at(imu));
      h_muon_dt4D -> Fill(muon_t->at(imu) -  vtx4D_t);
       

      // dz < X mm
      int ntracks_dz[NCUTSDZ];
      int ntracks_dz_removed[NCUTSDZ];
      int ntracks_dz_kept[NCUTSDZ];
      float sumpt_dz[NCUTSDZ];
      float sumpt_dz_removed[NCUTSDZ];
      float sumpt_dz_kept[NCUTSDZ];

      for (int icut = 0; icut < NCUTSDZ; icut++){
	ntracks_dz[icut] = 0;
	ntracks_dz_removed[icut] = 0;
	ntracks_dz_kept[icut] = 0;
	sumpt_dz[icut] = 0;
	sumpt_dz_removed[icut] = 0;
	sumpt_dz_kept[icut] = 0;
      }

      // mva
      int ntracks_mva3d_kept[NCUTSMVA] ;
      int ntracks_mva4d_kept[NCUTSMVA] ;
      float sumpt_mva3d_kept[NCUTSMVA] ;
      float sumpt_mva4d_kept[NCUTSMVA] ;
      float sumpt_mva3d_weighted = 0.;
      float sumpt_mva4d_weighted = 0.;

      for (int icut = 0; icut < NCUTSMVA; icut++){
	ntracks_mva3d_kept[icut] = 0;
	ntracks_mva4d_kept[icut] = 0;
	sumpt_mva3d_kept[icut] = 0.;
	sumpt_mva4d_kept[icut] = 0.;
      }


      // start loop over tracks in isolation cone
      for (unsigned int itk = 0; itk < track_pt->size(); itk++){
	
	// -- exclude tracks corresponding to the isolation cone of the other muon(s)
	if (track_muIndex ->at(itk) != int(imu)) continue;
	
	float deta = fabs(muon_eta->at(imu) - track_eta->at(itk));
	float dphi = fabs(muon_phi->at(imu) - track_phi->at(itk));
	if (dphi > TMath::Pi()) dphi = 2*TMath::Pi()- dphi;
	if (dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
	float dr = sqrt(deta*deta+dphi*dphi);
	
	if ( dr > 0.3 ) continue;
	
	if ( fabs(track_eta->at(itk)) < 1.5 && track_pt->at(itk) < btlMinTrackPt ) continue;
	if ( fabs(track_eta->at(itk)) > 1.5 && track_pt->at(itk) < etlMinTrackPt ) continue;

	// -- emulate BTL and ETL efficiency
	bool useTrackTime = true;
	double rndEff = gRandom2->Uniform(0.,1.); 
	if ( fabs(track_eta->at(itk)) < 1.5 && rndEff > btlEfficiency ) useTrackTime = false; 
	if ( fabs(track_eta->at(itk)) > 1.5 && rndEff > etlEfficiency ) useTrackTime = false; 


	float dxy   = track_dxy4D->at(itk);
	float dz    = track_dz4D->at(itk);
	float dzsim = track_dz4D->at(itk) + vtx4D_z - vtxGen_z;
	float trackTime = track_tFastSim->at(itk);
	// -- extra-smearing not needed beacuase it's already added at ntuple prod.
	if (trackTime!=-999)  {
	  float extra_smearing = sqrt( timeResolution*timeResolution - 0.035*0.035); // smear up to match 40 ps time resolution
	  trackTime= trackTime + gRandom->Gaus(0, extra_smearing );
	}

	float dt    = trackTime - vtx4D_t;
	float dtsim = trackTime - vtxGen_t*1000000000.;
	float dtFullSim = track_t->at(itk) - vtxGen_t*1000000000.;

	dt = dtsim;
	dz = dzsim;
	
	// --- dz, dt bix cuts with timing FastSim
	if (fabs(dxy) < maxDxy) {

	  for (int icut = 0; icut < NCUTSDZ; icut++){
	    
	    // -- barrel
	    if (isBarrel && fabs(dz) < dzCut[icut]){
	      
	      if (icut == 1){ // dz = 2 mm 
		h_tracks_dz_vtx_barrel -> Fill(dz);
		if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_dz_vtx_barrel->Fill(dz);	
		if ( trackTime!= -999){
		  h_tracks_dt_vtx_barrel->Fill(dt);
		  h_tracks_pt_withtiming_barrel -> Fill(track_pt->at(itk));
		  if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_dt_vtx_barrel->Fill(dt);
		  h2_tracks_dzvtx_dtvtx_barrel->Fill(dz, dt);
		}
		h_tracks_pt_barrel -> Fill( track_pt->at(itk) );
		p_tracks_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
	      }
	      ntracks_dz[icut]++;
	      sumpt_dz[icut]+=track_pt->at(itk);

	      // -- removed tracks
	      if ( useTrackTime && trackTime != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
		if (icut == 1){
		  h_tracks_removed_pt_barrel -> Fill( track_pt->at(itk) );
		  p_tracks_removed_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
		}
		ntracks_dz_removed[icut]++;
		sumpt_dz_removed[icut]+=track_pt->at(itk);
	      }
	      // -- kept tracks
	      else{
		if (icut == 1){
		  h_tracks_kept_pt_barrel -> Fill( track_pt->at(itk) );
		  p_tracks_kept_pt_vs_linedensity_barrel -> Fill( linedensity, track_pt->at(itk) );
		}
		ntracks_dz_kept[icut]++;
		sumpt_dz_kept[icut]+=track_pt->at(itk);
	      }
	    }// end if barrel
	    // -- endcap
	    if (!isBarrel && fabs(dz) < dzCut[icut]){
	      if (icut == 2){ 
		h_tracks_dz_vtx_endcap -> Fill(dz);
		if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_dz_vtx_endcap->Fill(dz);
		if ( trackTime!= -999){
		  h_tracks_dt_vtx_endcap->Fill(dt);
		  h_tracks_pt_withtiming_endcap-> Fill(track_pt->at(itk));
		  if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_dt_vtx_endcap->Fill(dt);
		  h2_tracks_dzvtx_dtvtx_endcap->Fill(dz, dt);
		}
		h_tracks_pt_endcap -> Fill( track_pt->at(itk) );
		p_tracks_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
	      }
	      ntracks_dz[icut]++;
	      sumpt_dz[icut]+=track_pt->at(itk);
	      
	      // -- removed tracks
	      if ( useTrackTime && trackTime != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
		if (icut == 2){
		  h_tracks_removed_pt_endcap -> Fill( track_pt->at(itk) );
		  p_tracks_removed_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
		}
		ntracks_dz_removed[icut]++;
		sumpt_dz_removed[icut]+=track_pt->at(itk);
	      }
	      // -- kept tracks
	      else{
		if (icut == 2){
		  h_tracks_kept_pt_endcap -> Fill( track_pt->at(itk) );
		  p_tracks_kept_pt_vs_linedensity_endcap -> Fill( linedensity, track_pt->at(itk) );
		}
		ntracks_dz_kept[icut]++;
		sumpt_dz_kept[icut]+=track_pt->at(itk);
	      }
	    }// end endcap
	  }// end loop over dz cuts
	}

	
	// BDT 3D and 4D
	float dz3D = fabs(track_dz3D->at(itk));
	float dz4D = fabs(track_dz4D->at(itk));
	if (isBarrel) { 
	  if ( dz3D < 1.0 ) {
	    h_tracks_puid3dmva_barrel ->Fill(track_puid3dmva->at(itk));
	    if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_puid3dmva_barrel ->Fill(track_puid3dmva->at(itk));
	  }
	  if ( dz4D < 1.0 ) {
	    h_tracks_puid4dmva_barrel ->Fill(track_puid4dmva->at(itk));
	    if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_puid4dmva_barrel ->Fill(track_puid4dmva->at(itk));
	  }
	}
	else{
	  if ( dz3D < 1.0 ) {
	    h_tracks_puid3dmva_endcap ->Fill(track_puid3dmva->at(itk));
	    if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_puid3dmva_endcap ->Fill(track_puid3dmva->at(itk));
	  }
	  if ( dz4D < 1.0 ) {
	    h_tracks_puid4dmva_endcap ->Fill(track_puid4dmva->at(itk));
	    if (track_isMatchedToGenParticle->at(itk)) h_tracks_genmatched_puid4dmva_endcap ->Fill(track_puid4dmva->at(itk));
	  }
	}
	
	if (dz3D < 1.0){
	  sumpt_mva3d_weighted+=track_pt->at(itk)*track_puid3dmva->at(itk);
	} 

	if (dz4D < 1.0){
	  sumpt_mva4d_weighted+=track_pt->at(itk)*track_puid4dmva->at(itk);
	} 

	if (dz3D < 1.0 && dz4D < 1.0){
	  if (track_isMatchedToGenParticle->at(itk)){
	    if (isBarrel){
	      h_tracks_genmatched_dz_vtx_fullSim_barrel ->Fill(dz);
	      h_tracks_genmatched_dt_vtx_fullSim_barrel ->Fill(dtFullSim);
	    } 
	    else{
	      h_tracks_genmatched_dz_vtx_fullSim_endcap ->Fill(dz);
	      h_tracks_genmatched_dt_vtx_fullSim_endcap ->Fill(dtFullSim);
	    }
	  }
	}

	// 3D MVA
	for (int icut = 0; icut < NCUTSMVA; icut++){
	  //if ( dz3D < 1.0  && track_puid3dmva->at(itk) > 0.7228171 ){
	  if ( dz3D < 1.0  && (track_puid3dmva->at(itk) == -1 || track_puid3dmva->at(itk) > mvaCut3D[icut]) ){
	    // -- all tracks kept by mva cut
	    ntracks_mva3d_kept[icut]++;
	    sumpt_mva3d_kept[icut]+=track_pt->at(itk);
	    if (isBarrel){
	      if (icut == 2){
		h_tracks_pass3dmva_dz_vtx_fullSim_barrel ->Fill(dz);
		h_tracks_pass3dmva_dt_vtx_fullSim_barrel ->Fill(dtFullSim);
		h_tracks_pass3dmva_pt_barrel->Fill(track_pt->at(itk));
		h_tracks_pass3dmva_eta_barrel->Fill(track_eta->at(itk));
	      }
	    }
	    else{
	      if (icut == 2){
		h_tracks_pass3dmva_dz_vtx_fullSim_endcap ->Fill(dz);
		h_tracks_pass3dmva_dt_vtx_fullSim_endcap ->Fill(dtFullSim);
		h_tracks_pass3dmva_pt_endcap->Fill(track_pt->at(itk));
		h_tracks_pass3dmva_eta_endcap->Fill(track_eta->at(itk));
	      }
	    }
	  }
	}


	// 4D MVA
        for (int icut = 0; icut < NCUTSMVA; icut++){
	  //if ( dz4D < 1.0 && track_puid4dmva->at(itk) > 0.7534892 ){
	  if ( dz4D < 1.0  && (track_puid4dmva->at(itk) == -1 ||  track_puid4dmva->at(itk) > mvaCut4D[icut] ) ){	  
	    // -- all tracks kept by mva cut 
	    ntracks_mva4d_kept[icut]++;
	    sumpt_mva4d_kept[icut]+=track_pt->at(itk);
	    if (isBarrel){
	      if (icut == 2){
		h_tracks_pass4dmva_dz_vtx_fullSim_barrel ->Fill(dz);
		h_tracks_pass4dmva_dt_vtx_fullSim_barrel ->Fill(dtFullSim);
		h_tracks_pass4dmva_pt_barrel->Fill(track_pt->at(itk));
		h_tracks_pass4dmva_eta_barrel->Fill(track_eta->at(itk));
	      }
	    }
	    else{
	      if (icut == 2){
		h_tracks_pass4dmva_dz_vtx_fullSim_endcap ->Fill(dz);
		h_tracks_pass4dmva_dt_vtx_fullSim_endcap ->Fill(dtFullSim);
		h_tracks_pass4dmva_pt_endcap->Fill(track_pt->at(itk));
		h_tracks_pass4dmva_eta_endcap->Fill(track_eta->at(itk));
	      }
	    }
	  }
	}
		  
      }// end loop over tracks
      

      // --- Fill histograms
 
      if (isBarrel){
	// -- all tracks in isolation cone
	h_tracks_n_barrel -> Fill( ntracks_dz[1] );
	h_tracks_sumpt_barrel -> Fill( sumpt_dz[1] );
	p_tracks_n_vs_linedensity_barrel -> Fill( linedensity, ntracks_dz[1] );
	p_tracks_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_dz[1] );

	// -- tracks removed from isolation cone
	h_tracks_removed_n_barrel -> Fill( ntracks_dz_removed[1] );
	h_tracks_removed_sumpt_barrel -> Fill( sumpt_dz_removed[1] );
	p_tracks_removed_n_vs_linedensity_barrel -> Fill( linedensity, ntracks_dz_removed[1] );
	p_tracks_removed_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_dz_removed[1] );

	// -- tracks kept from isolation cone
	h_tracks_kept_n_barrel -> Fill( ntracks_dz_kept[1] );
	h_tracks_kept_sumpt_barrel -> Fill( sumpt_dz_kept[1] );
	p_tracks_kept_n_vs_linedensity_barrel -> Fill( linedensity, ntracks_dz_kept[1] );
	p_tracks_kept_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_dz_kept[1] );

	// -- rel chIso
	for (int icut = 0; icut < NCUTSDZ; icut++){
	  h_muon_relChIso03_dZ_barrel[icut] ->Fill(sumpt_dz[icut]/muon_pt->at(imu));
	  h_muon_relChIso03_dZ_dT3s_barrel[icut] ->Fill(sumpt_dz_kept[icut]/muon_pt->at(imu));
	}

	h_muon_relChIso03_mva3D_weighted_barrel ->Fill(sumpt_mva3d_weighted/muon_pt->at(imu));
	h_muon_relChIso03_mva4D_weighted_barrel ->Fill(sumpt_mva4d_weighted/muon_pt->at(imu));
      
	for (int icut = 0; icut < NCUTSMVA; icut++){
	  h_muon_relChIso03_mva3D_barrel[icut] ->Fill(sumpt_mva3d_kept[icut]/muon_pt->at(imu));
	  h_muon_relChIso03_mva4D_barrel[icut] ->Fill(sumpt_mva4d_kept[icut]/muon_pt->at(imu));
	}

      }
      else{
	h_tracks_n_endcap -> Fill( ntracks_dz[2] );
	h_tracks_sumpt_endcap -> Fill( sumpt_dz[2] );
	p_tracks_n_vs_linedensity_endcap -> Fill( linedensity, ntracks_dz[2] );
	p_tracks_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_dz[2] );

	h_tracks_removed_n_endcap -> Fill( ntracks_dz_removed[2] );
	h_tracks_removed_sumpt_endcap -> Fill( sumpt_dz_removed[2]);
	p_tracks_removed_n_vs_linedensity_endcap -> Fill( linedensity, ntracks_dz_removed[2] );
	p_tracks_removed_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_dz_removed[2] );

	h_tracks_kept_n_endcap -> Fill( ntracks_dz_kept[2] );
	h_tracks_kept_sumpt_endcap -> Fill( sumpt_dz_kept[2] );
	p_tracks_kept_n_vs_linedensity_endcap -> Fill( linedensity, ntracks_dz_kept[2] );
	p_tracks_kept_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_dz_kept[2] );

	// -- rel chIso
	for (int icut = 0; icut < NCUTSDZ; icut++){
	  h_muon_relChIso03_dZ_endcap[icut] ->Fill(sumpt_dz[icut]/muon_pt->at(imu));
	  h_muon_relChIso03_dZ_dT3s_endcap[icut] ->Fill(sumpt_dz_kept[icut]/muon_pt->at(imu));
	}	

	h_muon_relChIso03_mva3D_weighted_endcap ->Fill(sumpt_mva3d_weighted/muon_pt->at(imu));
	h_muon_relChIso03_mva4D_weighted_endcap ->Fill(sumpt_mva4d_weighted/muon_pt->at(imu));

	for (int icut = 0; icut < NCUTSMVA; icut++){
	  h_muon_relChIso03_mva3D_endcap[icut] ->Fill(sumpt_mva3d_kept[icut]/muon_pt->at(imu));
	  h_muon_relChIso03_mva4D_endcap[icut] ->Fill(sumpt_mva4d_kept[icut]/muon_pt->at(imu));
	}

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
  
  
  // -- save histograms in output file
  std::string foutName = "testTracksMva_" + process + "_" + pu + "_dT"+ std::to_string(nsigma)+"sigma"+
                         "_minMuonPt"+std::to_string(int(minMuonPt))+
                         "_maxMuonPt"+std::to_string(int(maxMuonPt))+
                         "_minTrackPt_noDxy.root";

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
  
  h2_tracks_dzvtx_dtvtx_barrel-> Write();
  h2_tracks_dzvtx_dtvtx_endcap-> Write();

  h_tracks_genmatched_dz_vtx_barrel->Write();
  h_tracks_genmatched_dz_vtx_endcap->Write();

  h_tracks_genmatched_dt_vtx_barrel->Write();
  h_tracks_genmatched_dt_vtx_endcap->Write();

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


  h_tracks_pt_withtiming_barrel -> Write();
  h_tracks_pt_withtiming_endcap -> Write();


  for (int icut = 0; icut < NCUTSDZ; icut++){
    h_muon_relChIso03_dZ_barrel[icut] -> Write();
    h_muon_relChIso03_dZ_dT3s_barrel[icut] -> Write();

    h_muon_relChIso03_dZ_endcap[icut] -> Write();
    h_muon_relChIso03_dZ_dT3s_endcap[icut] -> Write();
  }


  h_tracks_puid3dmva_barrel ->Write();
  h_tracks_puid4dmva_barrel ->Write();

  h_tracks_puid3dmva_endcap ->Write();
  h_tracks_puid4dmva_endcap ->Write();


  h_tracks_genmatched_puid3dmva_barrel ->Write();
  h_tracks_genmatched_puid4dmva_barrel ->Write();

  h_tracks_genmatched_puid3dmva_endcap ->Write();
  h_tracks_genmatched_puid4dmva_endcap ->Write();

  h_muon_relChIso03_mva3D_weighted_barrel -> Write();
  h_muon_relChIso03_mva4D_weighted_barrel -> Write();

  h_muon_relChIso03_mva3D_weighted_endcap -> Write();
  h_muon_relChIso03_mva4D_weighted_endcap -> Write();

  for (int icut = 0; icut < NCUTSMVA; icut++){
    h_muon_relChIso03_mva3D_barrel[icut]->Write();
    h_muon_relChIso03_mva4D_barrel[icut]->Write();
    h_muon_relChIso03_mva3D_endcap[icut]->Write();
    h_muon_relChIso03_mva4D_endcap[icut]->Write();
  }


  h_tracks_pass3dmva_dz_vtx_fullSim_barrel->Write();
  h_tracks_pass3dmva_dz_vtx_fullSim_endcap->Write();

  h_tracks_pass3dmva_dt_vtx_fullSim_barrel->Write();
  h_tracks_pass3dmva_dt_vtx_fullSim_endcap->Write();

  h_tracks_pass4dmva_dz_vtx_fullSim_barrel->Write();
  h_tracks_pass4dmva_dz_vtx_fullSim_endcap->Write();

  h_tracks_pass4dmva_dt_vtx_fullSim_barrel->Write();
  h_tracks_pass4dmva_dt_vtx_fullSim_endcap->Write();


  h_tracks_genmatched_dz_vtx_fullSim_barrel->Write();
  h_tracks_genmatched_dz_vtx_fullSim_endcap->Write();

  h_tracks_genmatched_dt_vtx_fullSim_barrel->Write();
  h_tracks_genmatched_dt_vtx_fullSim_endcap->Write();

  h_tracks_genmatched_dz_vtx_fullSim_barrel->Write();
  h_tracks_genmatched_dz_vtx_fullSim_endcap->Write();

  h_tracks_genmatched_dt_vtx_fullSim_barrel->Write();
  h_tracks_genmatched_dt_vtx_fullSim_endcap->Write();



  h_tracks_pass3dmva_pt_barrel->Write();
  h_tracks_pass3dmva_pt_endcap->Write();
  h_tracks_pass4dmva_pt_barrel->Write();
  h_tracks_pass4dmva_pt_endcap->Write();
  h_tracks_pass3dmva_eta_barrel->Write();
  h_tracks_pass3dmva_eta_endcap->Write();
  h_tracks_pass4dmva_eta_barrel->Write();
  h_tracks_pass4dmva_eta_endcap->Write();

  fout->Close();
  
}



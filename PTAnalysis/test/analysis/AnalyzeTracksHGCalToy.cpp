// per compilare: g++ -Wall -o AnalyzeTracksHGCalToy `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit AnalyzeTracksHGCalToy.cpp                        

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

// -- parametrization of HGCal efficiency
float HGCalEfficiency(float pt) {
  float eff = 1.-exp(-(pt-0.2)/1.4) ;
  //float eff = 1.;
  //if (pt < 2.)  eff = 0.;
  return eff;
}

// -- parametrization of HGCal time resolution
float HGCalTimeResolution(float pt) {
  float resol = sqrt(3600./pt+600.) / 1000. ; ///in ns
  //float resol = 0.030;
  return resol;
}

// weighted average
void combineTimes(float t1, float t2, float t1res, float t2res, float &tcomb, float &tres){
  float w1 = 1./t1res/t1res;
  float w2 = 1./t2res/t2res;
  float wsum = w1+w2;
  w1/=wsum;
  w2/=wsum;
  tcomb = w1*t1+w2*t2;
  tres = sqrt(1./wsum);
}




int main(int argc, char** argv){

  float defaultTimeResolution = 0.030;
  // ---- ETL resolution and efficiency
  // -- 2 layers: 35 ps with 90%, 1 layer 50 ps with 90%x85% efficiency
  float timeResolutionETL2 = 0.035;
  float effETL2 = 0.90;
  float timeResolutionETL1 = 0.050;
  float effETL1 = 0.90 * 0.85;

  float timeResolutionETL1opt = 0.035;
  float minTkEtaETLsmall = 2.3; // this correspons to an ETL radius of 60 cm instead of 127 cm

  int nsigma = 3;

  float minMuonPt = 20.;
  float maxMuonPt = 9999.;
  float maxDz = 0.3; // best dz cut for etl tracks is 2 mm;
  float maxDxy = 9999.; // was 0.02
  float btlMinTrackPt = 0.7;
  float etlMinTrackPt = 0.4;
  float relChIsoCut = 0.08;

  string process = argv[1];
  bool prompt = false;

  string pu = argv[2];
  
  // -- get TChain
  TChain* chain = new TChain("analysis/tree_30ps","tree");
  if (process.find("DYToLL") != std::string::npos) {
    //if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_muIso/190111_164118/0000/muonIsolation*.root");
    //if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_noPU_muIso/190111_164205/0000/muonIsolation*.root");
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_muIso_minPt/190208_131558/0000/muonIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_noPU_muIso_minPt/190208_131709/0000/muonIsolation*.root");
        prompt = true;
  }
  
  if (process.find("TTbar") != std::string::npos) {
    //if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_muIso/190111_164134/0000/muonIsolation*.root");
    //if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_noPU_muIso/190111_164224/0000/muonIsolation*.root");
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_muIso_minPt/190208_131612/0000/muonIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_noPU_muIso_minPt/190208_131723/0000/muonIsolation*.root");
    prompt = false;
  }
  
  
  cout << "Using prompt muons = " << prompt <<endl;

  cout << "Analyzing " << chain->GetEntries() << "  events" <<endl;
    
  
  
  //if PU=0, don't cut on dz
  //if ( pu.find("noPU")  != std::string::npos ) {
  //  maxDz = 9999.;
  //  maxDxy = 9999.;
  //}
  
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
  vector<int> *track_muIndex;
  vector<int> *track_isMatchedToGenParticle;
  
  muon_pt = 0;
  muon_eta = 0;
  muon_phi = 0;
  muon_dz3D = 0;
  muon_dz4D = 0;
  muon_dxy3D = 0;
  muon_dxy4D = 0;
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
  track_muIndex = 0;
  track_isMatchedToGenParticle = 0;

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
  chain->SetBranchStatus("track_muIndex",1);          chain->SetBranchAddress("track_muIndex",  &track_muIndex);
  chain->SetBranchStatus("track_isMatchedToGenParticle",1);          chain->SetBranchAddress("track_isMatchedToGenParticle",  &track_isMatchedToGenParticle);

  
  // control plot for HGCal model
  TProfile *p_hgcal_efficiency_vs_pt = new TProfile("p_hgcal_efficiency_vs_pt","p_hgcal_efficiency_vs_pt", 100, 0., 20.);
  TProfile *p_hgcal_timeResolution_vs_pt = new TProfile("p_hgcal_timeResolution_vs_pt","p_hgcal_timeResolution_vs_pt", 100, 0., 20.);
  



  // chargedIso for rocs
  // -- no timing
  TH1F *h_muon_relChIso03_dZ2 = new TH1F("h_muon_relChIso03_dZ2","h_muon_relChIso03_dZ2",5000,0,5); 
  // -- HGCal alone
  TH1F *h_muon_relChIso03_dZ2_HGCal = new TH1F("h_muon_relChIso03_dZ2_HGCal","h_muon_relChIso03_dZ2_HGCal",5000, 0, 5); 
  // -- 2 layer MTD (40 ps 90% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_ETL2 = new TH1F("h_muon_relChIso03_dZ2_ETL2","h_muon_relChIso03_dZ2_ETL2",5000, 0, 5);
  // -- 1 layer MTD (50 ps 85*90% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_ETL1 = new TH1F("h_muon_relChIso03_dZ2_ETL1","h_muon_relChIso03_dZ2_ETL1",5000, 0, 5);
  // -- 2 layer MTD, smalle radius (40 ps 90% efficiency 60 cm radius --> eta 2.3-3.0) 
  TH1F *h_muon_relChIso03_dZ2_ETL2small = new TH1F("h_muon_relChIso03_dZ2_ETL2small","h_muon_relChIso03_dZ2_ETL2small",5000, 0, 5);
  // -- 1 layer MTD (35 ps 85*90% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_ETL1opt = new TH1F("h_muon_relChIso03_dZ2_ETL1opt","h_muon_relChIso03_dZ2_ETL1opt",5000, 0, 5);
  // -- HGCal (+) 1 layer MTD (50 ps 90*85% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_HGCal_ETL1 = new TH1F("h_muon_relChIso03_dZ2_HGCal_ETL1","h_muon_relChIso03_dZ2_HGCal_ETL1",5000, 0, 5);
  // -- HGCal (+) 2 layers MTD (40 ps 90% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_HGCal_ETL2 = new TH1F("h_muon_relChIso03_dZ2_HGCal_ETL2","h_muon_relChIso03_dZ2_HGCal_ETL2",5000, 0, 5);
  // -- HGCal (+) 1 layer MTD (35 ps 90*85% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_HGCal_ETL1opt = new TH1F("h_muon_relChIso03_dZ2_HGCal_ETL1opt","h_muon_relChIso03_dZ2_HGCal_ETL1opt",5000, 0, 5);
  // -- HGCal (+) 2 layers MTD small radius(40 ps 90% efficiency) 
  TH1F *h_muon_relChIso03_dZ2_HGCal_ETL2small = new TH1F("h_muon_relChIso03_dZ2_HGCal_ETL2small","h_muon_relChIso03_dZ2_HGCal_ETL2small",5000, 0, 5);


  // eff vs line density
  TProfile *p_efficiency_relChIso03_dZ2_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_vs_linedensity","p_efficiency_relChIso03_dZ2_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity","p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity","p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity","p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity","p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity","p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity","p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity","p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity","p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity", 21, -0.05, 2.05);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity","p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity", 21, -0.05, 2.05);

  // eff vs pT
  TProfile *p_efficiency_relChIso03_dZ2_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_vs_muonpt","p_efficiency_relChIso03_dZ2_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt","p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt","p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt","p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt","p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt","p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt","p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt","p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt","p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt", 100, 0., 100.);
  TProfile *p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt = new TProfile("p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt","p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt", 100, 0., 100.);


  // control plot
  TH1F *h_dt_HGCal = new TH1F("h_dt_HGCal","h_dt_HGCal", 500, -0.5, 0.5);
  TH1F *h_dt_ETL1  = new TH1F("h_dt_ETL1","h_dt_ETL1", 500, -0.5, 0.5);
  TH1F *h_dt_ETL2  = new TH1F("h_dt_ETL2","h_dt_ETL2", 500, -0.5, 0.5);
  TH1F *h_dt_HGCal_ETL2  = new TH1F("h_dt_HGCal_ETL2","h_dt_HGCal_ETL2", 500, -0.5, 0.5);



  TRandom *gRandom = new TRandom();
  TRandom *gRandomHGCal = new TRandom();
  TRandom *gRandomETL1 = new TRandom();
  TRandom *gRandomETL2 = new TRandom();


  int maxEntries = chain->GetEntries();
  for (int ientry = 0; ientry< maxEntries; ientry++) {
    
    chain -> GetEntry(ientry);
    
    if (ientry%1000==0) cout << "Analyzing event " << ientry << "\r" << flush;
    
    float linedensity = npu*TMath::Gaus(fabs(10.*vtxGen_z), 0, 42., 1);

    for (unsigned int imu = 0; imu < muon_pt->size(); imu++){
      
      // -- prompt or fake muons
      bool pass = false;
      if ( prompt) pass = muon_isPrompt->at(imu);
      if (!prompt) pass = !muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu) && !muon_isFromTauDecay->at(imu);
      if (!pass) continue;


      // -- reco muon selections
      if (muon_pt->at(imu) < minMuonPt ) continue;
      if (muon_pt->at(imu) > maxMuonPt ) continue;
      if (!muon_isLoose->at(imu)) continue;
      if (fabs(muon_eta->at(imu)) < 1.5 ) continue; // only endcap muons
    
      
      bool pass3D = fabs(muon_dz3D->at(imu)) < 0.5 && fabs(muon_dxy3D->at(imu)) < 0.2 && !vtx3D_isFake;
      bool pass4D = fabs(muon_dz4D->at(imu)) < 0.5 && fabs(muon_dxy4D->at(imu)) < 0.2 && !vtx4D_isFake;
      
      float sumpt = 0.;
      float sumpt_hgcal = 0.;
      float sumpt_hgcal_etl1 = 0.;
      float sumpt_hgcal_etl2 = 0.;
      float sumpt_hgcal_etl1opt = 0.;
      float sumpt_hgcal_etl2small = 0.;
      float sumpt_etl2 = 0.;
      float sumpt_etl1 = 0.;
      float sumpt_etl2small = 0.;
      float sumpt_etl1opt = 0.;
      
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
	
	if ( fabs(track_eta->at(itk)) < 1.5 && track_pt->at(itk) < btlMinTrackPt ) continue;
	if ( fabs(track_eta->at(itk)) > 1.5 && track_pt->at(itk) < etlMinTrackPt ) continue;
	//if ( fabs(track_eta->at(itk)) > 3.0 ) continue;

	float dxy3D = track_dxy3D->at(itk);
	float dxy4D = track_dxy4D->at(itk);
	//float dz3D  = track_dz3D->at(itk); 
	//float dz4D  = track_dz4D->at(itk);
	float dzsim = track_dz4D->at(itk) + vtx4D_z - vtxGen_z;
	float t0 =  vtxGen_t*1000000000.;
	
	// -- emulate HGCal efficiency
	float effHGCal  = HGCalEfficiency(track_pt->at(itk));
	bool trackHasHGCalTime = true;
	double rndEffHGCal = gRandomHGCal->Uniform(0.,1.);
	if (rndEffHGCal > effHGCal) trackHasHGCalTime = false;
	p_hgcal_efficiency_vs_pt ->Fill(track_pt->at(itk),effHGCal);	

	// -- emulate ETL 1-layer efficiency
	bool trackHasETL1Time = true;
        double rndEffETL1 = gRandomETL1->Uniform(0.,1.);
	if (rndEffETL1 > effETL1) trackHasETL1Time = false;
	
	// -- emulate ETL 2-layers efficiency
	bool trackHasETL2Time = true;
        double rndEffETL2 = gRandomETL2->Uniform(0.,1.);
	if (rndEffETL2 > effETL2) trackHasETL2Time = false;
	
	
	if ( fabs(track_eta->at(itk)) < 1.59 || fabs(track_eta->at(itk)) > 3.0 ) {
	  trackHasHGCalTime = false;
	  trackHasETL1Time  = false;
	  trackHasETL2Time  = false;
	}

	if ( track_t->at(itk) == -999){
	  trackHasHGCalTime = false;
          trackHasETL1Time  = false;
          trackHasETL2Time  = false;
	}
	
	// -- emulate HGCal time resolution
	float timeResolutionHGCal = HGCalTimeResolution(track_pt->at(itk));
	float extra_resol = 0.;
	if ( timeResolutionHGCal > defaultTimeResolution  ) extra_resol = sqrt(timeResolutionHGCal*timeResolutionHGCal - defaultTimeResolution *defaultTimeResolution );
	else extra_resol = 0.;
	double rnd0 = gRandom->Gaus(0., extra_resol);
	float tHGCal = track_t->at(itk) + rnd0;
	p_hgcal_timeResolution_vs_pt ->Fill(track_pt->at(itk),timeResolutionHGCal);	
	
	// -- emulate ETL 1-layer time resolution with extra smearing to emulate different time resolution
	if (timeResolutionETL1 > defaultTimeResolution ) extra_resol = sqrt(timeResolutionETL1*timeResolutionETL1 - defaultTimeResolution *defaultTimeResolution );
	else extra_resol =0.;
	double rnd1 = gRandom->Gaus(0., extra_resol);
	float tETL1 = track_t->at(itk) + rnd1;
	
	// -- emulate ETL 2-layers time resolution with extra smearing to emulate different time resolution
	if (timeResolutionETL2 > defaultTimeResolution ) extra_resol = sqrt(timeResolutionETL2*timeResolutionETL2 - defaultTimeResolution *defaultTimeResolution );
	else extra_resol =0.;
	double rnd2 = gRandom->Gaus(0., extra_resol);
	float tETL2 = track_t->at(itk) + rnd2;
		
	// -- emulate ETL 1-layer time resolution with extra smearing to emulate 35 ps time resolution
	if (timeResolutionETL1opt > defaultTimeResolution ) extra_resol = sqrt(timeResolutionETL1opt*timeResolutionETL1opt - defaultTimeResolution *defaultTimeResolution );
	else extra_resol =0.;
	double rnd1opt = gRandom->Gaus(0., extra_resol);
	float tETL1opt = track_t->at(itk) + rnd1opt;
	
	// -- no timing 
	if ( fabs(dzsim) < maxDz && fabs(dxy3D) < maxDxy ){
	  sumpt+=track_pt->at(itk);
	}
	
	// -- with timing
	if ( fabs(dzsim) < maxDz && fabs(dxy4D) < maxDxy){
	  
	  // -- HGCal alone
	  if ( !trackHasHGCalTime || (trackHasHGCalTime && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal) ) {
	      sumpt_hgcal+=track_pt->at(itk);
	  }
	  
	  // -- ETL 1 layer alone
	  if ( !trackHasETL1Time || (trackHasETL1Time && fabs( tETL1 - t0 ) < nsigma * timeResolutionETL1)) {
	    sumpt_etl1+=track_pt->at(itk);
	  }

	  // -- ETL 2 layers alone
	  if ( !trackHasETL2Time || (trackHasETL2Time && fabs( tETL2 - t0 ) < nsigma * timeResolutionETL2)) {
	    sumpt_etl2+=track_pt->at(itk);
	  }
	
	  // -- ETL 1 layer alone with optimized (35 ps) time resolution
          if ( !trackHasETL1Time || (trackHasETL1Time && fabs( tETL1opt - t0 ) < nsigma * timeResolutionETL1opt)) {
            sumpt_etl1opt+=track_pt->at(itk);
          }

	  // -- ETL 2 layers with smaller radius (eta 2.3-3.0)
	  if ( !trackHasETL2Time || (trackHasETL2Time && fabs(track_eta->at(itk)) < minTkEtaETLsmall) || (trackHasETL2Time && fabs(track_eta->at(itk)) > minTkEtaETLsmall && fabs( tETL2 - t0 ) < nsigma * timeResolutionETL2)) {
	    sumpt_etl2small+=track_pt->at(itk);
	  }

  	  
	  // -- HGCal  + 1 MTD layer
	  float tComb, timeResolutionComb;
	  combineTimes(tHGCal, tETL1, timeResolutionHGCal, timeResolutionETL1, tComb, timeResolutionComb);
	  if (timeResolutionComb > defaultTimeResolution ) {
	    extra_resol = sqrt(timeResolutionComb*timeResolutionComb - defaultTimeResolution *defaultTimeResolution );
	  }
	  else {
	    extra_resol =0.;
	    timeResolutionComb = 0.030;	  
	  }
	  tComb = track_t->at(itk) + gRandom->Gaus(0., extra_resol);
	  if ( !trackHasHGCalTime && !trackHasETL1Time) sumpt_hgcal_etl1+=track_pt->at(itk);
	  if ( !trackHasETL1Time  && trackHasHGCalTime && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal) sumpt_hgcal_etl1+=track_pt->at(itk); 
	  if ( !trackHasHGCalTime && trackHasETL1Time  && fabs( tETL1  - t0 ) < nsigma * timeResolutionETL1 ) sumpt_hgcal_etl1+=track_pt->at(itk);
	  if (  trackHasHGCalTime && trackHasETL1Time  && fabs( tComb  - t0 ) < nsigma * timeResolutionComb ) sumpt_hgcal_etl1+=track_pt->at(itk);
	  		  
	  // -- HGCal  + 2 MTD layers
	  combineTimes(tHGCal, tETL2, timeResolutionHGCal, timeResolutionETL2, tComb, timeResolutionComb );
	  if (timeResolutionComb > defaultTimeResolution ) {
	    extra_resol = sqrt(timeResolutionComb*timeResolutionComb - defaultTimeResolution *defaultTimeResolution );
	  }
	  else {
	    extra_resol = 0.;
	    timeResolutionComb = 0.030;
	  }
	  tComb = track_t->at(itk) + gRandom->Gaus(0., extra_resol);
	  if ( !trackHasHGCalTime && !trackHasETL2Time) sumpt_hgcal_etl2+=track_pt->at(itk);
	  if ( !trackHasETL2Time  && trackHasHGCalTime && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal) sumpt_hgcal_etl2+=track_pt->at(itk); 
	  if ( !trackHasHGCalTime && trackHasETL2Time  && fabs( tETL2  - t0 ) < nsigma * timeResolutionETL2 ) sumpt_hgcal_etl2+=track_pt->at(itk);
	  if (  trackHasHGCalTime && trackHasETL2Time  && fabs( tComb  - t0 ) < nsigma * timeResolutionComb ) sumpt_hgcal_etl2+=track_pt->at(itk);
	  


	  // -- HGCal  + 1 MTD layer optimized to 35 ps
          combineTimes(tHGCal, tETL1opt, timeResolutionHGCal, timeResolutionETL1opt, tComb, timeResolutionComb);
          if (timeResolutionComb > defaultTimeResolution ) {
            extra_resol = sqrt(timeResolutionComb*timeResolutionComb - defaultTimeResolution *defaultTimeResolution );
          }
          else {
            extra_resol =0.;
            timeResolutionComb = 0.030;
          }
          tComb = track_t->at(itk) + gRandom->Gaus(0., extra_resol);
          if ( !trackHasHGCalTime && !trackHasETL1Time) sumpt_hgcal_etl1opt+=track_pt->at(itk);
          if ( !trackHasETL1Time  && trackHasHGCalTime && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal) sumpt_hgcal_etl1opt+=track_pt->at(itk);
          if ( !trackHasHGCalTime && trackHasETL1Time  && fabs( tETL1opt  - t0 ) < nsigma * timeResolutionETL1opt ) sumpt_hgcal_etl1opt+=track_pt->at(itk);
          if (  trackHasHGCalTime && trackHasETL1Time  && fabs( tComb  - t0 ) < nsigma * timeResolutionComb ) sumpt_hgcal_etl1opt+=track_pt->at(itk);


	  // -- HGCal  + 2 MTD layers at smaller radius
          combineTimes(tHGCal, tETL2, timeResolutionHGCal, timeResolutionETL2, tComb, timeResolutionComb );
          if (timeResolutionComb > defaultTimeResolution ) {
            extra_resol = sqrt(timeResolutionComb*timeResolutionComb - defaultTimeResolution *defaultTimeResolution );
          }
          else {
            extra_resol = 0.;
            timeResolutionComb = 0.030;
          }
          tComb = track_t->at(itk) + gRandom->Gaus(0., extra_resol);
          if ( !trackHasHGCalTime && !trackHasETL2Time) sumpt_hgcal_etl2small+=track_pt->at(itk);
          if ( !trackHasETL2Time && trackHasHGCalTime && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal) sumpt_hgcal_etl2small+=track_pt->at(itk);
          if ( !trackHasHGCalTime && trackHasETL2Time && fabs(track_eta->at(itk))> minTkEtaETLsmall && fabs( tETL2  - t0 ) < nsigma * timeResolutionETL2 ) sumpt_hgcal_etl2small+=track_pt->at(itk);
	  if ( !trackHasHGCalTime && trackHasETL2Time && fabs(track_eta->at(itk))< minTkEtaETLsmall) sumpt_hgcal_etl2small+=track_pt->at(itk);
          if (  trackHasHGCalTime && trackHasETL2Time && fabs(track_eta->at(itk))> minTkEtaETLsmall  && fabs( tComb  - t0 ) < nsigma * timeResolutionComb ) sumpt_hgcal_etl2small+=track_pt->at(itk);
          if (  trackHasHGCalTime && trackHasETL2Time && fabs(track_eta->at(itk))< minTkEtaETLsmall  && fabs( tHGCal - t0 ) < nsigma * timeResolutionHGCal ) sumpt_hgcal_etl2small+=track_pt->at(itk);


	  // control plots
	  if (track_isMatchedToGenParticle->at(itk)){
	    if (trackHasHGCalTime) h_dt_HGCal ->Fill(tHGCal - t0 );
	    if (trackHasETL1Time) h_dt_ETL1 ->Fill(tETL1 - t0 );
	    if (trackHasETL2Time) h_dt_ETL2 ->Fill(tETL2 - t0 );
	    if (trackHasETL2Time && trackHasHGCalTime) h_dt_HGCal_ETL2 ->Fill(tComb - t0 );
	  }
	}
      
      }// end loop over tracks
      
      
      // --  fill relChIso histograms and eff vs line density
      
      float muonpt = muon_pt->at(imu);
      
      if (pass3D) {
	h_muon_relChIso03_dZ2 -> Fill(sumpt/muonpt);
      }

      if (pass4D){
	h_muon_relChIso03_dZ2_HGCal -> Fill(sumpt_hgcal/muonpt);
	h_muon_relChIso03_dZ2_ETL2 -> Fill(sumpt_etl2/muonpt);
	h_muon_relChIso03_dZ2_ETL1 -> Fill(sumpt_etl1/muonpt);
	h_muon_relChIso03_dZ2_ETL2small -> Fill(sumpt_etl2small/muonpt);
	h_muon_relChIso03_dZ2_ETL1opt -> Fill(sumpt_etl1opt/muonpt);
	h_muon_relChIso03_dZ2_HGCal_ETL1 -> Fill(sumpt_hgcal_etl1/muonpt);
	h_muon_relChIso03_dZ2_HGCal_ETL2 -> Fill(sumpt_hgcal_etl2/muonpt);
	h_muon_relChIso03_dZ2_HGCal_ETL1opt -> Fill(sumpt_hgcal_etl1opt/muonpt);
	h_muon_relChIso03_dZ2_HGCal_ETL2small -> Fill(sumpt_hgcal_etl2small/muonpt);
      }

      
      if (pass3D) {
	if ( sumpt/muonpt < relChIsoCut) {
	  p_efficiency_relChIso03_dZ2_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_vs_muonpt->Fill(muonpt, 0.);
	}
      }


      if (pass4D){

	if ( sumpt_hgcal/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt->Fill(muonpt, 0.);
	}
	

	if ( sumpt_etl1/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt->Fill(muonpt, 0.);
	}


	if ( sumpt_etl2/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt->Fill(muonpt, 0.);
	}



	if ( sumpt_etl1opt/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt->Fill(muonpt, 0.);
	}


	if ( sumpt_etl2small/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt->Fill(muonpt, 0.);
	}


	
	if ( sumpt_hgcal_etl1/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt->Fill(muonpt, 0.);
	}

      
	if ( sumpt_hgcal_etl2/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt->Fill(muonpt, 0.);
	}


	if ( sumpt_hgcal_etl1opt/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt->Fill(muonpt, 0.);
	}

      
	if ( sumpt_hgcal_etl2small/muonpt < relChIsoCut){
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity->Fill(linedensity, 1.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt->Fill(muonpt, 1.);
	}
	else{
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity->Fill(linedensity, 0.);
	  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt->Fill(muonpt, 0.);
	}

      }

    }// end loop over muons
    
  }// end loop over events


  // -- save histograms in output file
  std::string foutName = "testTracksHGCal_" + process + "_" + pu + "_noDxy_new.root";
  
  TFile *fout = new TFile(foutName.c_str(),"recreate");
  
  p_hgcal_efficiency_vs_pt->Write();
  p_hgcal_timeResolution_vs_pt->Write();

  h_muon_relChIso03_dZ2 -> Write();
  h_muon_relChIso03_dZ2_HGCal -> Write();
  h_muon_relChIso03_dZ2_ETL1 -> Write();
  h_muon_relChIso03_dZ2_ETL2 -> Write();
  h_muon_relChIso03_dZ2_ETL1opt -> Write();
  h_muon_relChIso03_dZ2_ETL2small -> Write();
  h_muon_relChIso03_dZ2_HGCal_ETL1 -> Write();
  h_muon_relChIso03_dZ2_HGCal_ETL2 -> Write();
  h_muon_relChIso03_dZ2_HGCal_ETL1opt -> Write();
  h_muon_relChIso03_dZ2_HGCal_ETL2small -> Write();

  p_efficiency_relChIso03_dZ2_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_ETL1_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_ETL2_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_ETL1opt_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_ETL2small_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_linedensity -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_linedensity -> Write();

  p_efficiency_relChIso03_dZ2_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_ETL1_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_ETL2_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_ETL1opt_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_ETL2small_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL1_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL2_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL1opt_vs_muonpt -> Write();
  p_efficiency_relChIso03_dZ2_HGCal_ETL2small_vs_muonpt -> Write();


  h_dt_HGCal -> Write();
  h_dt_ETL1 -> Write();
  h_dt_ETL2 -> Write();
  h_dt_HGCal_ETL2 -> Write();


  fout->Close();
  
}



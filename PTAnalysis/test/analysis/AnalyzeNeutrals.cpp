// per compilare: g++ -Wall -o AnalyzeNeutrals `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit AnalyzeNeutrals.cpp                        

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

  float timeResolution = 0.035;
  int nsigma = 3;

  float minMuonPt = 20.;
  float maxMuonPt = 9999.;
  
  string process = argv[1];
  bool prompt = false;

  string pu = argv[2];
  
  // -- get TChain
  TChain* chain = new TChain("analysis/tree_35ps","tree");
  if (process.find("DYToLL") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/NeutralIso/DYToLL_M-50_14TeV_TuneCP5_pythia8/test_DYToLL_neutralMuIso/190212_133305/0000/muonNeutrIsolation*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/NeutralIso/DYToLL_M-50_14TeV_TuneCP5_pythia8/test_DYToLL_noPU_neutralMuIso/190212_133029/0000/muonNeutrIsolation*.root");
    prompt = true;
  }
  
  if (process.find("TTbar") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/MTD/NeutralIso/TTbar_14TeV_TuneCP5_Pythia8/test_TTbar_neutralMuIso/190212_145344/0000/muonNeutrIsolation_*.root"); 
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_noPU_muIso/190111_164224/0000/muonNeutrIsolation_*.root");
    prompt = false;
  }
  
  if (process.find("QCD") != std::string::npos) {
    if (pu.find("PU200") != std::string::npos) chain->Add("/eos/cms/store/user/malberti/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/test_QCD_muIso/190206_095346/0000/muonNeutrIsolation_*.root");
    if (pu.find("noPU")  != std::string::npos) chain->Add("/eos/cms/store/user/malberti/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/test_QCD_noPU_muIso/190111_164240/0000/muonNeutrIsolation_*.root");
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
  vector<float> *neutrPfCand_particleId;
  vector<float> *neutrPfCand_pt;
  vector<float> *neutrPfCand_eta;
  vector<float> *neutrPfCand_phi;
  vector<float> *neutrPfCand_tCluster;
  vector<float> *neutrPfCand_dRcluster;
  vector<float> *neutrPfCand_dRmu;
  vector<int> *neutrPfCand_muIndex;
  
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
  neutrPfCand_pt = 0;
  neutrPfCand_eta = 0;
  neutrPfCand_phi = 0;
  neutrPfCand_particleId= 0;
  neutrPfCand_dRcluster = 0;
  neutrPfCand_tCluster = 0;
  neutrPfCand_dRmu = 0;
  neutrPfCand_muIndex = 0;

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
  chain->SetBranchStatus("muon_isFromTauDecay",1);    chain->SetBranchAddress("muon_isFromTauDecay", &muon_isFromTauDecay);
  
  chain->SetBranchStatus("neutrPfCand_pt",1);               chain->SetBranchAddress("neutrPfCand_pt",        &neutrPfCand_pt);
  chain->SetBranchStatus("neutrPfCand_eta",1);              chain->SetBranchAddress("neutrPfCand_eta",       &neutrPfCand_eta);
  chain->SetBranchStatus("neutrPfCand_phi",1);              chain->SetBranchAddress("neutrPfCand_phi",       &neutrPfCand_phi);
  chain->SetBranchStatus("neutrPfCand_particleId",1);       chain->SetBranchAddress("neutrPfCand_particleId",&neutrPfCand_particleId);
  chain->SetBranchStatus("neutrPfCand_tCluster",1);         chain->SetBranchAddress("neutrPfCand_tCluster",  &neutrPfCand_tCluster);
  chain->SetBranchStatus("neutrPfCand_dRcluster",1);        chain->SetBranchAddress("neutrPfCand_dRcluster", &neutrPfCand_dRcluster);
  chain->SetBranchStatus("neutrPfCand_dRmu",1);             chain->SetBranchAddress("neutrPfCand_dRmu",      &neutrPfCand_dRmu);
  chain->SetBranchStatus("neutrPfCand_muIndex",1);          chain->SetBranchAddress("neutrPfCand_muIndex",   &neutrPfCand_muIndex);

  TH1F *hsimvtx_z = new TH1F("hsimvtx_z","hsimvtx_z", 500,-25,25);
  TH1F *hsimvtx_t = new TH1F("hsimvtx_t","hsimvtx_t", 1000,-1,1);

  TH1F *h_muon_pt_barrel = new TH1F("h_muon_pt_barrel","h_muon_pt_barrel",100,0,100);
  TH1F *h_muon_pt_endcap = new TH1F("h_muon_pt_endcap","h_muon_pt_endcap",100,0,100);

  TH1F *h_neutrals_dt_vtx_barrel = new TH1F("h_neutrals_dt_vtx_barrel","h_neutrals_dt_vtx_barrel",5000, -5.0, 20.0);
  TH1F *h_neutrals_dt_vtx_endcap = new TH1F("h_neutrals_dt_vtx_endcap","h_neutrals_dt_vtx_endcap",5000, -5.0, 20.0);

  TH1F *h_neutrals_matchingToMtd_barrel = new TH1F("h_neutrals_matchingToMtd_barrel","h_neutrals_matchingToMtd_barrel",2, -0.5, 1.5);
  TH1F *h_neutrals_matchingToMtd_endcap = new TH1F("h_neutrals_matchingToMtd_endcap","h_neutrals_matchingToMtd_endcap",2, -0.5, 1.5);

  TH1F *h_neutrals_dRcluster_barrel = new TH1F("h_neutrals_dRcluster_barrel","h_neutrals_dRcluster_barrel",100, 0., 0.2);
  TH1F *h_neutrals_dRcluster_endcap = new TH1F("h_neutrals_dRcluster_endcap","h_neutrals_dRcluster_endcap",100, 0., 0.2);


  // -- to debug tracks in cone 
  TH1F *h_neutrals_pt_barrel = new TH1F("h_neutrals_pt_barrel","h_neutrals_pt_barrel",400,0.,20.);
  TH1F *h_neutrals_removed_pt_barrel = new TH1F("h_neutrals_removed_pt_barrel","h_neutrals_removed_pt_barrel",400,0.,20.);
  TH1F *h_neutrals_kept_pt_barrel = new TH1F("h_neutrals_kept_pt_barrel","h_neutrals_kept_pt_barrel",400,0.,20.);

  TH1F *h_neutrals_n_barrel = new TH1F("h_neutrals_n_barrel","h_npfcands_pt_barrel",50,0.,50.);
  TH1F *h_neutrals_removed_n_barrel = new TH1F("h_neutrals_removed_n_barrel","h_neutrals_removed_n_barrel",50,0.,50.);
  TH1F *h_neutrals_kept_n_barrel = new TH1F("h_neutrals_kept_n_barrel","h_neutrals_kept_n_barrel",50,0.,50.);

  TH1F *h_neutrals_sumpt_barrel = new TH1F("h_neutrals_sumpt_barrel","h_neutrals_sumpt_barrel",1000,0.,100.);
  TH1F *h_neutrals_removed_sumpt_barrel = new TH1F("h_neutrals_removed_sumpt_barrel","h_neutrals_removed_sumpt_barrel",1000,0.,100.);
  TH1F *h_neutrals_kept_sumpt_barrel = new TH1F("h_neutrals_kept_sumpt_barrel","h_neutrals_kept_sumpt_barrel",1000,0.,100.);

  TProfile *p_neutrals_pt_vs_linedensity_barrel = new TProfile("p_neutrals_pt_vs_linedensity_barrel","p_neutrals_pt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_removed_pt_vs_linedensity_barrel = new TProfile("p_neutrals_removed_pt_vs_linedensity_barrel","p_neutrals_removed_pt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_kept_pt_vs_linedensity_barrel = new TProfile("p_neutrals_kept_pt_vs_linedensity_barrel","p_neutrals_kept_pt_vs_linedensity_barrel",22,-0.05,2.05);

  TProfile *p_neutrals_n_vs_linedensity_barrel = new TProfile("p_neutrals_n_vs_linedensity_barrel","p_neutrals_n_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_removed_n_vs_linedensity_barrel = new TProfile("p_neutrals_removed_n_vs_linedensity_barrel","p_neutrals_removed_n_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_kept_n_vs_linedensity_barrel = new TProfile("p_neutrals_kept_n_vs_linedensity_barrel","p_neutrals_kept_n_vs_linedensity_barrel",22,-0.05,2.05);
  
  TProfile *p_neutrals_sumpt_vs_linedensity_barrel = new TProfile("p_neutrals_sumpt_vs_linedensity_barrel","p_neutrals_sumpt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_removed_sumpt_vs_linedensity_barrel = new TProfile("p_neutrals_removed_sumpt_vs_linedensity_barrel","p_neutrals_removed_sumpt_vs_linedensity_barrel",22,-0.05,2.05);
  TProfile *p_neutrals_kept_sumpt_vs_linedensity_barrel = new TProfile("p_neutrals_kept_sumpt_vs_linedensity_barrel","p_neutrals_kept_sumpt_vs_linedensity_barrel",22,-0.05,2.05);

  TH1F *h_neutrals_pt_endcap = new TH1F("h_neutrals_pt_endcap","h_neutrals_pt_endcap",400,0.,20.);
  TH1F *h_neutrals_removed_pt_endcap = new TH1F("h_neutrals_removed_pt_endcap","h_neutrals_removed_pt_endcap",400,0.,20.);
  TH1F *h_neutrals_kept_pt_endcap = new TH1F("h_neutrals_kept_pt_endcap","h_neutrals_kept_pt_endcap",400,0.,20.);

  TH1F *h_neutrals_n_endcap = new TH1F("h_neutrals_n_endcap","h_npfcands_pt_endcap",50,0.,50.);
  TH1F *h_neutrals_removed_n_endcap = new TH1F("h_neutrals_removed_n_endcap","h_neutrals_removed_n_endcap",50,0.,50.);
  TH1F *h_neutrals_kept_n_endcap = new TH1F("h_neutrals_kept_n_endcap","h_neutrals_kept_n_endcap",50,0.,50.);

  TH1F *h_neutrals_sumpt_endcap = new TH1F("h_neutrals_sumpt_endcap","h_neutrals_sumpt_endcap",1000,0.,100.);
  TH1F *h_neutrals_removed_sumpt_endcap = new TH1F("h_neutrals_removed_sumpt_endcap","h_neutrals_removed_sumpt_endcap",1000,0.,100.);
  TH1F *h_neutrals_kept_sumpt_endcap = new TH1F("h_neutrals_kept_sumpt_endcap","h_neutrals_kept_sumpt_endcap",1000,0.,100.);

  TProfile *p_neutrals_pt_vs_linedensity_endcap = new TProfile("p_neutrals_pt_vs_linedensity_endcap","p_neutrals_pt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_removed_pt_vs_linedensity_endcap = new TProfile("p_neutrals_removed_pt_vs_linedensity_endcap","p_neutrals_removed_pt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_kept_pt_vs_linedensity_endcap = new TProfile("p_neutrals_kept_pt_vs_linedensity_endcap","p_neutrals_kept_pt_vs_linedensity_endcap",22,-0.05,2.05);

  TProfile *p_neutrals_n_vs_linedensity_endcap = new TProfile("p_neutrals_n_vs_linedensity_endcap","p_neutrals_n_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_removed_n_vs_linedensity_endcap = new TProfile("p_neutrals_removed_n_vs_linedensity_endcap","p_neutrals_removed_n_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_kept_n_vs_linedensity_endcap = new TProfile("p_neutrals_kept_n_vs_linedensity_endcap","p_neutrals_kept_n_vs_linedensity_endcap",22,-0.05,2.05);
  
  TProfile *p_neutrals_sumpt_vs_linedensity_endcap = new TProfile("p_neutrals_sumpt_vs_linedensity_endcap","p_neutrals_sumpt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_removed_sumpt_vs_linedensity_endcap = new TProfile("p_neutrals_removed_sumpt_vs_linedensity_endcap","p_neutrals_removed_sumpt_vs_linedensity_endcap",22,-0.05,2.05);
  TProfile *p_neutrals_kept_sumpt_vs_linedensity_endcap = new TProfile("p_neutrals_kept_sumpt_vs_linedensity_endcap","p_neutrals_kept_sumpt_vs_linedensity_endcap",22,-0.05,2.05);


  TProfile *p_neutrals_n_vs_npu_barrel = new TProfile("p_neutrals_n_vs_npu_barrel","p_neutrals_n_vs_npu_barrel",100, 150, 250);
  TProfile *p_neutrals_n_vs_npu_endcap = new TProfile("p_neutrals_n_vs_npu_endcap","p_neutrals_n_vs_npu_endcap",100, 150, 250);


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

      if ( isBarrel) h_muon_pt_barrel   -> Fill(muon_pt->at(imu));
      if ( !isBarrel) h_muon_pt_endcap   -> Fill(muon_pt->at(imu));

      int npfcands = 0;
      int npfcands_removed = 0;
      int npfcands_kept = 0;
      float sumpt = 0.;
      float sumpt_removed = 0.;
      float sumpt_kept = 0.;

      // start loop over pfcands in isolation cone
      for (unsigned int icand = 0; icand < neutrPfCand_pt->size(); icand++){
	
	// -- exclude pfcands corresponding to the isolation cone of the other muon(s)
	if (neutrPfCand_muIndex ->at(icand) != int(imu)) continue;
	
	float drclus = neutrPfCand_dRcluster->at(icand);
	float dt = neutrPfCand_tCluster->at(icand) - vtxGen_t*1000000000.;
	float pt = neutrPfCand_pt->at(icand);

	// -- barrel
	if (isBarrel){
	  h_neutrals_dRcluster_barrel->Fill( drclus );
	  if (neutrPfCand_tCluster->at(icand) != -999){
	    h_neutrals_matchingToMtd_barrel->Fill(1.);
	    h_neutrals_dt_vtx_barrel->Fill(dt);
	  }
	  else{
	    h_neutrals_matchingToMtd_barrel->Fill(0.);
	  }
	  // -- all neutrals
	  h_neutrals_pt_barrel -> Fill( pt );
	  p_neutrals_pt_vs_linedensity_barrel -> Fill( linedensity, pt );
	  npfcands++;
	  sumpt+=pt;
	  // -- removed pfcands
	  if ( neutrPfCand_tCluster->at(icand) != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
	    h_neutrals_removed_pt_barrel -> Fill( pt );
	    p_neutrals_removed_pt_vs_linedensity_barrel -> Fill( linedensity, pt );
	    npfcands_removed++;
	    sumpt_removed+=pt;
	  }
	  else{
	    h_neutrals_kept_pt_barrel -> Fill( pt );
	    p_neutrals_kept_pt_vs_linedensity_barrel -> Fill( linedensity, pt );
	    npfcands_kept++;
	    sumpt_kept+=pt;
	  }
	}// end if barrel

	// endcap
	else{
	  h_neutrals_dRcluster_endcap->Fill( drclus );
	  if (neutrPfCand_tCluster->at(icand) != -999){
	    h_neutrals_matchingToMtd_endcap->Fill(1.);
	    h_neutrals_dt_vtx_endcap->Fill(dt);
	  }
	  else{
	    h_neutrals_matchingToMtd_endcap->Fill(0.);
	  }
	  // -- all neutrals
	  h_neutrals_pt_endcap -> Fill( pt );
	  p_neutrals_pt_vs_linedensity_endcap -> Fill( linedensity, pt );
	  npfcands++;
	  sumpt+=pt;
	  // -- removed pfcands
	  if ( neutrPfCand_tCluster->at(icand) != -999 && fabs(dt) > float(nsigma) * timeResolution ) {
	    h_neutrals_removed_pt_endcap -> Fill( pt );
	    p_neutrals_removed_pt_vs_linedensity_endcap -> Fill( linedensity, pt );
	    npfcands_removed++;
	    sumpt_removed+=pt;
	  }
	  else{
	    h_neutrals_kept_pt_endcap -> Fill( pt );
	    p_neutrals_kept_pt_vs_linedensity_endcap -> Fill( linedensity, pt );
	    npfcands_kept++;
	    sumpt_kept+=pt;
	  }
	}	

      } // end loop over neutral pf cands

      

      //float relChIso = sumpt/muon_pt->at(imu);
      //float relChIso_dT = sumpt_kept/muon_pt->at(imu);
      
      if (isBarrel){
	// -- all pfcands in isolation cone
	h_neutrals_n_barrel -> Fill( npfcands );
	h_neutrals_sumpt_barrel -> Fill( sumpt );
	p_neutrals_n_vs_linedensity_barrel -> Fill( linedensity, npfcands );
	p_neutrals_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt );
	p_neutrals_n_vs_npu_barrel  -> Fill( npu, npfcands );
	// -- pfcands removed from isolation cone
	h_neutrals_removed_n_barrel -> Fill( npfcands_removed );
	h_neutrals_removed_sumpt_barrel -> Fill( sumpt_removed );
	p_neutrals_removed_n_vs_linedensity_barrel -> Fill( linedensity, npfcands_removed );
	p_neutrals_removed_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_removed );

	// -- pfcands kept from isolation cone
	h_neutrals_kept_n_barrel -> Fill( npfcands_kept );
	h_neutrals_kept_sumpt_barrel -> Fill( sumpt_kept );
	p_neutrals_kept_n_vs_linedensity_barrel -> Fill( linedensity, npfcands_kept );
	p_neutrals_kept_sumpt_vs_linedensity_barrel -> Fill( linedensity, sumpt_kept );
      }
      else{
	h_neutrals_n_endcap -> Fill( npfcands );
	h_neutrals_sumpt_endcap -> Fill( sumpt );
	p_neutrals_n_vs_linedensity_endcap -> Fill( linedensity, npfcands );
	p_neutrals_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt );
	p_neutrals_n_vs_npu_endcap  -> Fill( npu, npfcands );

	h_neutrals_removed_n_endcap -> Fill( npfcands_removed );
	h_neutrals_removed_sumpt_endcap -> Fill( sumpt_removed );
	p_neutrals_removed_n_vs_linedensity_endcap -> Fill( linedensity, npfcands_removed );
	p_neutrals_removed_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_removed );

	h_neutrals_kept_n_endcap -> Fill( npfcands_kept );
	h_neutrals_kept_sumpt_endcap -> Fill( sumpt_kept );
	p_neutrals_kept_n_vs_linedensity_endcap -> Fill( linedensity, npfcands_kept );
	p_neutrals_kept_sumpt_vs_linedensity_endcap -> Fill( linedensity, sumpt_kept );
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
  h_neutrals_dt_vtx_barrel->Scale(1./h_muon_pt_barrel->GetSumOfWeights());
  h_neutrals_dt_vtx_barrel->Fit("fitfun_barrel","QR");
  
  TF1 *fitfun_endcap = new TF1("fitfun_endcap","gaus(0)+gaus(3)", -10, 10);
  fitfun_endcap->SetParameter(1, 0.);
  fitfun_endcap->SetParameter(2, 0.030);
  fitfun_endcap->SetParameter(4, 0.);
  fitfun_endcap->SetParameter(5, 0.270);
  fitfun_endcap->SetNpx(1000);
  h_neutrals_dt_vtx_endcap->Scale(1./h_muon_pt_endcap->GetSumOfWeights());
  h_neutrals_dt_vtx_endcap->Fit("fitfun_endcap","QR");
  
  
  
  // -- save histograms in output file
  std::string foutName = "testNeutrals_" + process + "_" + pu +
                         "_minMuonPt"+std::to_string(int(minMuonPt))+
                         "_maxMuonPt"+std::to_string(int(maxMuonPt))+
                         ".root";

  TFile *fout = new TFile(foutName.c_str(),"recreate");

  hsimvtx_z->Write();
  hsimvtx_t->Write();

  h_muon_pt_barrel->Write();
  h_muon_pt_endcap->Write();

  h_neutrals_matchingToMtd_barrel->Write();
  h_neutrals_matchingToMtd_endcap->Write();

  h_neutrals_dRcluster_barrel->Write();
  h_neutrals_dRcluster_endcap->Write();

  h_neutrals_dt_vtx_barrel->Write();
  h_neutrals_dt_vtx_endcap->Write();
  
  h_neutrals_pt_barrel ->Write();
  h_neutrals_removed_pt_barrel ->Write();
  h_neutrals_kept_pt_barrel ->Write();

  h_neutrals_n_barrel ->Write();
  h_neutrals_removed_n_barrel ->Write();
  h_neutrals_kept_n_barrel ->Write();

  h_neutrals_sumpt_barrel ->Write();
  h_neutrals_removed_sumpt_barrel ->Write();
  h_neutrals_kept_sumpt_barrel ->Write();

  p_neutrals_pt_vs_linedensity_barrel -> Write();
  p_neutrals_removed_pt_vs_linedensity_barrel -> Write();
  p_neutrals_kept_pt_vs_linedensity_barrel -> Write();

  p_neutrals_n_vs_linedensity_barrel -> Write();
  p_neutrals_removed_n_vs_linedensity_barrel -> Write();
  p_neutrals_kept_n_vs_linedensity_barrel -> Write();

  p_neutrals_sumpt_vs_linedensity_barrel -> Write();
  p_neutrals_removed_sumpt_vs_linedensity_barrel -> Write();
  p_neutrals_kept_sumpt_vs_linedensity_barrel -> Write();

  h_neutrals_pt_endcap ->Write();
  h_neutrals_removed_pt_endcap ->Write();
  h_neutrals_kept_pt_endcap ->Write();

  h_neutrals_n_endcap ->Write();
  h_neutrals_removed_n_endcap ->Write();
  h_neutrals_kept_n_endcap ->Write();

  h_neutrals_sumpt_endcap ->Write();
  h_neutrals_removed_sumpt_endcap ->Write();
  h_neutrals_kept_sumpt_endcap ->Write();

  p_neutrals_pt_vs_linedensity_endcap -> Write();
  p_neutrals_removed_pt_vs_linedensity_endcap -> Write();
  p_neutrals_kept_pt_vs_linedensity_endcap -> Write();

  p_neutrals_n_vs_linedensity_endcap -> Write();
  p_neutrals_removed_n_vs_linedensity_endcap -> Write();
  p_neutrals_kept_n_vs_linedensity_endcap -> Write();

  p_neutrals_sumpt_vs_linedensity_endcap -> Write();
  p_neutrals_removed_sumpt_vs_linedensity_endcap -> Write();
  p_neutrals_kept_sumpt_vs_linedensity_endcap -> Write();

  p_neutrals_n_vs_npu_barrel-> Write();
  p_neutrals_n_vs_npu_endcap -> Write();

  fout->Close();
  
}



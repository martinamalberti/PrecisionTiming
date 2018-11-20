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

  float maxDzVtx = 0.1;
  float maxDzMu  = 0.1;
  
  string process = argv[1];
  bool prompt = false;
  
  // -- get TChain
  TChain* chain = new TChain("analysis/tree_30ps","tree");
  if (process.find("DYToLL") != std::string::npos) {
    chain->Add("/eos/cms/store/user/malberti/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/test_DYToLL_muIso/181117_083340/0000/muonIsolation_*.root");
    prompt = true;
  }
  
  if (process.find("TTbar") != std::string::npos) {
    chain->Add("/eos/cms/store/user/malberti/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/test_TTbar_muIso/181117_083704/0000/muonIsolation_*.root");
    prompt = false;
  }
  
  if (process.find("QCD") != std::string::npos) {
    chain->Add("/eos/cms/store/user/malberti/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/test_QCD_muIso/181119_120017/0000/muonIsolation_*.root");
    prompt = false; 
  }
  
  cout << "Using prompt muons = " << prompt <<endl;
  
  
  // -- tree vars
  float vtx3D_z;
  float vtx4D_z;
  float vtx4D_t;
  int vtx4D_isFake;
  float vtxGen_t;
  float vtxGen_z;
  vector<int> *muon_isLoose;
  vector<int> *muon_isPrompt;
  vector<int> *muon_isMatchedToGenJet;
  vector<float> *muon_pt;
  vector<float> *muon_eta;
  vector<float> *muon_phi;
  vector<float> *muon_t;
  vector<float> *muon_dz3D;
  vector<float> *muon_dz4D;
  vector<float> *track_pt;
  vector<float> *track_eta;
  vector<float> *track_phi;
  vector<float> *track_dz3D;
  vector<float> *track_dz4D;
  vector<float> *track_t;
  vector<int> *track_muIndex;
  
  muon_pt = 0;
  muon_eta = 0;
  muon_phi = 0;
  muon_dz3D = 0;
  muon_dz4D = 0;
  muon_t = 0;
  muon_isLoose= 0;
  muon_isMatchedToGenJet = 0;
  muon_isPrompt = 0;
  track_pt = 0;
  track_eta = 0;
  track_phi = 0;
  track_dz3D = 0;
  track_dz4D = 0;
  track_t = 0;
  track_muIndex = 0;
  
  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("vtx3D_z",1);                chain->SetBranchAddress("vtx3D_z",        &vtx3D_z);
  chain->SetBranchStatus("vtx4D_z",1);                chain->SetBranchAddress("vtx4D_z",        &vtx4D_z);
  chain->SetBranchStatus("vtx4D_t",1);                chain->SetBranchAddress("vtx4D_t",        &vtx4D_t);
  chain->SetBranchStatus("vtx4D_isFake",1);           chain->SetBranchAddress("vtx4D_isFake",   &vtx4D_isFake);
  chain->SetBranchStatus("vtxGen_z",1);               chain->SetBranchAddress("vtxGen_z",       &vtxGen_z);
  chain->SetBranchStatus("vtxGen_t",1);               chain->SetBranchAddress("vtxGen_t",       &vtxGen_t);
  chain->SetBranchStatus("muon_pt",1);                chain->SetBranchAddress("muon_pt",        &muon_pt);
  chain->SetBranchStatus("muon_eta",1);               chain->SetBranchAddress("muon_eta",       &muon_eta);
  chain->SetBranchStatus("muon_phi",1);               chain->SetBranchAddress("muon_phi",       &muon_phi);
  chain->SetBranchStatus("muon_t",1);                 chain->SetBranchAddress("muon_t",         &muon_t);
  chain->SetBranchStatus("muon_dz3D",1);              chain->SetBranchAddress("muon_dz3D",      &muon_dz3D);
  chain->SetBranchStatus("muon_dz4D",1);              chain->SetBranchAddress("muon_dz4D",      &muon_dz4D);
  chain->SetBranchStatus("muon_isLoose",1);           chain->SetBranchAddress("muon_isLoose",   &muon_isLoose);
  chain->SetBranchStatus("muon_isPrompt",1);          chain->SetBranchAddress("muon_isPrompt",  &muon_isPrompt);
  chain->SetBranchStatus("muon_isMatchedToGenJet",1); chain->SetBranchAddress("muon_isMatchedToGenJet", &muon_isMatchedToGenJet);
  
  chain->SetBranchStatus("track_pt",1);               chain->SetBranchAddress("track_pt",       &track_pt);
  chain->SetBranchStatus("track_eta",1);              chain->SetBranchAddress("track_eta",      &track_eta);
  chain->SetBranchStatus("track_phi",1);              chain->SetBranchAddress("track_phi",      &track_phi);
  chain->SetBranchStatus("track_dz3D",1);             chain->SetBranchAddress("track_dz3D",     &track_dz3D);
  chain->SetBranchStatus("track_dz4D",1);             chain->SetBranchAddress("track_dz4D",     &track_dz4D);
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
  TH1F *h_vtx_dt4D_sim = new TH1F("h_vtx_dt4D_sim","h_vtx_dt4D_sim",1000, -0.5, 0.5);
  
  TH1F *h_tracks_dt_vtx_barrel = new TH1F("h_tracks_dt_vtx_barrel","h_tracks_dt_vtx_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_dz_vtx_barrel = new TH1F("h_tracks_dz_vtx_barrel","h_tracks_dz_vtx_barrel",500, -1.0, 1.0);
  
  TH1F *h_tracks_dt_mu_barrel = new TH1F("h_tracks_dt_mu_barrel","h_tracks_dt_mu_barrel",500, -1.0, 1.0);
  TH1F *h_tracks_dz_mu_barrel = new TH1F("h_tracks_dz_mu_barrel","h_tracks_dz_mu_barrel",500, -1.0, 1.0);
  

  TH1F *h_tracks_pt_notiming_barrel = new TH1F("h_tracks_pt_notiming_barrel","h_tracks_pt_notiming_barrel",1000, 0, 20);
  TH1F *h_tracks_pt_timing_PV_barrel = new TH1F("h_tracks_pt_timing_PV_barrel","h_tracks_pt_timing_PV_barrel",1000, 0, 20);
  TH1F *h_tracks_pt_timing_PU_barrel = new TH1F("h_tracks_pt_timing_PU_barrel","h_tracks_pt_timing_PU_barrel",1000, 0, 20);
  
  TH1F *h_tracks_eta_notiming_barrel = new TH1F("h_tracks_eta_notiming_barrel","h_tracks_eta_notiming_barrel",300, -3, 3);
  TH1F *h_tracks_eta_timing_PV_barrel = new TH1F("h_tracks_eta_timing_PV_barrel","h_tracks_eta_timing_PV_barrel",300, -3, 3);
  TH1F *h_tracks_eta_timing_PU_barrel = new TH1F("h_tracks_eta_timing_PU_barrel","h_tracks_eta_timing_PU_barrel",300, -3, 3);
  
  TH1F *h_tracks_dr_notiming_barrel = new TH1F("h_tracks_dr_notiming_barrel","h_tracks_dr_notiming_barrel",1000, 0, 1);
  TH1F *h_tracks_dr_timing_PV_barrel = new TH1F("h_tracks_dr_timing_PV_barrel","h_tracks_dr_timing_PV_barrel",1000, 0, 1);
  TH1F *h_tracks_dr_timing_PU_barrel = new TH1F("h_tracks_dr_timing_PU_barrel","h_tracks_dr_timing_PU_barrel",1000, 0, 1);


  TH1F *h_tracks_dt_vtx_endcap = new TH1F("h_tracks_dt_vtx_endcap","h_tracks_dt_vtx_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_dz_vtx_endcap = new TH1F("h_tracks_dz_vtx_endcap","h_tracks_dz_vtx_endcap",500, -1.0, 1.0);
  
  TH1F *h_tracks_dt_mu_endcap = new TH1F("h_tracks_dt_mu_endcap","h_tracks_dt_mu_endcap",500, -1.0, 1.0);
  TH1F *h_tracks_dz_mu_endcap = new TH1F("h_tracks_dz_mu_endcap","h_tracks_dz_mu_endcap",500, -1.0, 1.0);
  
  TH1F *h_tracks_pt_notiming_endcap = new TH1F("h_tracks_pt_notiming_endcap","h_tracks_pt_notiming_endcap",1000, 0, 20);
  TH1F *h_tracks_pt_timing_PV_endcap = new TH1F("h_tracks_pt_timing_PV_endcap","h_tracks_pt_timing_PV_endcap",1000, 0, 20);
  TH1F *h_tracks_pt_timing_PU_endcap = new TH1F("h_tracks_pt_timing_PU_endcap","h_tracks_pt_timing_PU_endcap",1000, 0, 20);
  
  TH1F *h_tracks_eta_notiming_endcap = new TH1F("h_tracks_eta_notiming_endcap","h_tracks_eta_notiming_endcap",300, -3, 3);
  TH1F *h_tracks_eta_timing_PV_endcap = new TH1F("h_tracks_eta_timing_PV_endcap","h_tracks_eta_timing_PV_endcap",300, -3, 3);
  TH1F *h_tracks_eta_timing_PU_endcap = new TH1F("h_tracks_eta_timing_PU_endcap","h_tracks_eta_timing_PU_endcap",300, -3, 3);
  
  TH1F *h_tracks_dr_notiming_endcap = new TH1F("h_tracks_dr_notiming_endcap","h_tracks_dr_notiming_endcap",1000, 0, 1);
  TH1F *h_tracks_dr_timing_PV_endcap = new TH1F("h_tracks_dr_timing_PV_endcap","h_tracks_dr_timing_PV_endcap",1000, 0, 1);
  TH1F *h_tracks_dr_timing_PU_endcap = new TH1F("h_tracks_dr_timing_PU_endcap","h_tracks_dr_timing_PU_endcap",1000, 0, 1);
  
  TH2F *h2_tracks_dzvtx_dtvtx_barrel = new TH2F("h2_tracks_dzvtx_dtvtx_barrel","h2_tracks_dzvtx_dtvtx_barrel",1000,-1,1, 1000, -1, 1);
  TH2F *h2_tracks_dzvtx_dtvtx_endcap = new TH2F("h2_tracks_dzvtx_dtvtx_endcap","h2_tracks_dzvtx_dtvtx_endcap",1000,-1,1, 1000, -1, 1);

  TH2F *h2_tracks_dzmu_dtmu_barrel = new TH2F("h2_tracks_dzmu_dtmu_barrel","h2_tracks_dzmu_dtmu_barrel",1000,-1,1, 1000, -1, 1);
  TH2F *h2_tracks_dzmu_dtmu_endcap = new TH2F("h2_tracks_dzmu_dtmu_endcap","h2_tracks_dzmu_dtmu_endcap",1000,-1,1, 1000, -1, 1);



  cout << "Analyzing " << chain->GetEntries() << "  events" <<endl;
  
  
  //int maxEntries = 500000;
  int maxEntries = chain->GetEntries();
  for (int ientry = 0; ientry< maxEntries; ientry++) {
    
    chain -> GetEntry(ientry);
    
    if (ientry%1000==0) cout << "Analyzing event " << ientry << "\r" << flush;
    
    for (unsigned int imu = 0; imu < muon_pt->size(); imu++){
      if (muon_pt->at(imu)< 20 ) continue;
      if (!muon_isLoose->at(imu) ) continue;
     
      bool isBarrel = fabs(muon_eta->at(imu)) < 1.5;

      bool pass = false;
      if ( prompt) pass = muon_isPrompt->at(imu);
      if (!prompt) pass = !muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu);
      
      if (!pass) continue;

      hsimvtx_z  -> Fill(vtxGen_z);    
      hsimvtx_t  -> Fill(vtxGen_t*1000000000.);    

      h_vtx_dz3D_sim -> Fill( vtx3D_z - vtxGen_z);
      h_vtx_dz4D_sim -> Fill( vtx4D_z - vtxGen_z);
      h_vtx_dt4D_sim -> Fill( vtx4D_t - vtxGen_t*1000000000.);




      //if (fabs(vtx4D_z - vtxGen_z) > 0.005 ) continue;
      //if (fabs(vtx4D_t - vtxGen_t*1000000000) > 0.015 ) continue;


      // test
      if (fabs(vtx4D_z - vtxGen_z) < 0.01 && fabs(vtx4D_t - vtxGen_t*1000000000) < 0.02) continue;

      
      if ( isBarrel) h_muon_pt_barrel   -> Fill(muon_pt->at(imu));
      if ( !isBarrel) h_muon_pt_endcap   -> Fill(muon_pt->at(imu));
      h_muon_dz3D -> Fill(muon_dz3D->at(imu));
      h_muon_dz4D -> Fill(muon_dz4D->at(imu));
      h_muon_dt4D -> Fill(muon_t->at(imu) -  vtx4D_t);
      
      for (unsigned int itk = 0; itk < track_eta->size(); itk++){
	
	// -- exclude tracks corresponding to the isolation cone of the other muon(s)
	if (track_muIndex ->at(itk) != int(imu)) continue;
	
	float deta = fabs(muon_eta->at(imu) - track_eta->at(itk));
	float dphi = fabs(muon_phi->at(imu) - track_phi->at(itk));
	if (dphi > TMath::Pi()) dphi = 2*TMath::Pi()- dphi;
	if (dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
	float dr = sqrt(deta*deta+dphi*dphi);
	
	if ( dr > 0.3 ) continue;
	
	float dz = track_dz4D->at(itk);
	float dt = track_t->at(itk) - vtx4D_t;
	if ( isBarrel) h_tracks_dz_vtx_barrel -> Fill(dz);
        if (!isBarrel) h_tracks_dz_vtx_endcap-> Fill(dz);

	if (fabs(dz) < maxDzVtx) {
	  // -- barrel
          if (isBarrel){
	    // if timing info available
	    if (track_t->at(itk) > -99){
	      h_tracks_dt_vtx_barrel->Fill(dt);
	      h2_tracks_dzvtx_dtvtx_barrel->Fill(dz, dt);
	      if ( fabs(dt) > 0.030 * 3) {
		h_tracks_pt_timing_PU_barrel->Fill(track_pt->at(itk));
                h_tracks_eta_timing_PU_barrel->Fill(track_eta->at(itk));
                h_tracks_dr_timing_PU_barrel->Fill(dr);
	      }
	      else {
	        h_tracks_pt_timing_PV_barrel->Fill(track_pt->at(itk));
                h_tracks_eta_timing_PV_barrel->Fill(track_eta->at(itk));
                h_tracks_dr_timing_PV_barrel->Fill(dr);
	      }
	    }
	    // if timing info not available
            else{
	      h_tracks_pt_notiming_barrel->Fill(track_pt->at(itk));
              h_tracks_eta_notiming_barrel->Fill(track_eta->at(itk));
              h_tracks_dr_notiming_barrel->Fill(dr); 
	    }
	  }// end if barrel
	  
          // -- endcap
          else {
            // if timing info available
            if (track_t->at(itk) > -99){
              h_tracks_dt_vtx_endcap->Fill(dt);
              h2_tracks_dzvtx_dtvtx_endcap->Fill(dz, dt);
              if ( fabs(dt) > 0.030 * 3) {
                h_tracks_pt_timing_PU_endcap->Fill(track_pt->at(itk));
                h_tracks_eta_timing_PU_endcap->Fill(track_eta->at(itk));
                h_tracks_dr_timing_PU_endcap->Fill(dr);
              }
              else {
                h_tracks_pt_timing_PV_endcap->Fill(track_pt->at(itk));
                h_tracks_eta_timing_PV_endcap->Fill(track_eta->at(itk));
                h_tracks_dr_timing_PV_endcap->Fill(dr);
              }
            }
            // if timing info not available
            else{
              h_tracks_pt_notiming_endcap->Fill(track_pt->at(itk));
              h_tracks_eta_notiming_endcap->Fill(track_eta->at(itk));
              h_tracks_dr_notiming_endcap->Fill(dr);
            }
          }
	}
      

	// wrt to muon
	float dtmu = track_t->at(itk) - muon_t->at(imu);
	float dzmu = muon_dz4D->at(imu) - track_dz4D->at(itk) ;
	if ( isBarrel) h_tracks_dz_mu_barrel -> Fill(dzmu);
	if (!isBarrel) h_tracks_dz_mu_endcap-> Fill(dzmu);
	
	if (fabs(dzmu) < maxDzMu){
	  // -- barrel
	  if (isBarrel){
	    // if timing info available
	    if (track_t->at(itk) > -99){	
	      h_tracks_dt_mu_barrel->Fill(dtmu);
	      h2_tracks_dzmu_dtmu_barrel->Fill(dzmu, dtmu);
	      /*
	      if ( fabs(dtmu) > 0.030 * sqrt(2) * 3) { 
		h_tracks_pt_timing_PU_barrel->Fill(track_pt->at(itk));
		h_tracks_eta_timing_PU_barrel->Fill(track_eta->at(itk));
		h_tracks_dr_timing_PU_barrel->Fill(dr);
	      }
	      else { 
		h_tracks_pt_timing_PV_barrel->Fill(track_pt->at(itk));
		h_tracks_eta_timing_PV_barrel->Fill(track_eta->at(itk));
		h_tracks_dr_timing_PV_barrel->Fill(dr);
	      }
	      */
	    }
	    // if timing info not available    
	    /*
	    else{
	      h_tracks_pt_notiming_barrel->Fill(track_pt->at(itk));
	      h_tracks_eta_notiming_barrel->Fill(track_eta->at(itk));
	      h_tracks_dr_notiming_barrel->Fill(dr);
	    }
	    */
	  }

	  // -- endcap
	  else {
	    // if timing info available 
	    if (track_t->at(itk) > -99){
              h_tracks_dt_mu_endcap->Fill(dtmu);
              h2_tracks_dzmu_dtmu_endcap->Fill(dzmu, dtmu);
	      /*
	      if ( fabs(dtmu) > 0.030 * sqrt(2) * 3) {
		h_tracks_pt_timing_PU_endcap->Fill(track_pt->at(itk));
		h_tracks_eta_timing_PU_endcap->Fill(track_eta->at(itk));
		h_tracks_dr_timing_PU_endcap->Fill(dr);
	      }
	      else {
		h_tracks_pt_timing_PV_endcap->Fill(track_pt->at(itk));
		h_tracks_eta_timing_PV_endcap->Fill(track_eta->at(itk));
		h_tracks_dr_timing_PV_endcap->Fill(dr);
	      }
	      */
            }
            // if timing info not available
            /*
	    else{
	      h_tracks_pt_notiming_endcap->Fill(track_pt->at(itk));
	      h_tracks_eta_notiming_endcap->Fill(track_eta->at(itk));
              h_tracks_dr_notiming_endcap->Fill(dr);
	    }
	    */
	  }
	}

      }// end loop over tracks
    }// end loop over muons
  }//end loop over events
  
  

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
  
  
  std::string foutName = "testTracks_" + process + ".root";
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
  
  h_tracks_pt_notiming_barrel->Write();
  h_tracks_pt_timing_PU_barrel->Write();
  h_tracks_pt_timing_PV_barrel->Write();
  
  h_tracks_eta_notiming_barrel->Write();
  h_tracks_eta_timing_PU_barrel->Write();
  h_tracks_eta_timing_PV_barrel->Write();
  
  h_tracks_dr_notiming_barrel->Write();
  h_tracks_dr_timing_PU_barrel->Write();
  h_tracks_dr_timing_PV_barrel->Write();
  

  
  h_tracks_dt_vtx_endcap->Write();
  h_tracks_dt_vtx_endcap->Write();
  
  h_tracks_dz_mu_endcap->Write();
  h_tracks_dt_mu_endcap->Write();
  
  h_tracks_pt_notiming_endcap->Write();
  h_tracks_pt_timing_PU_endcap->Write();
  h_tracks_pt_timing_PV_endcap->Write();
  
  h_tracks_eta_notiming_endcap->Write();
  h_tracks_eta_timing_PU_endcap->Write();
  h_tracks_eta_timing_PV_endcap->Write();
  
  h_tracks_dr_notiming_endcap->Write();
  h_tracks_dr_timing_PU_endcap->Write();
  h_tracks_dr_timing_PV_endcap->Write();

  h2_tracks_dzvtx_dtvtx_barrel->Write();  
  h2_tracks_dzvtx_dtvtx_endcap->Write();

  h2_tracks_dzmu_dtmu_barrel->Write();  
  h2_tracks_dzmu_dtmu_endcap->Write();
  
  fout->Close();
  
}



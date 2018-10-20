// per compilare: g++ -Wall -o MuonIsolationAnalysis `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit MuonIsolationAnalysis.cpp

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

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>


using namespace std;


int main(int argc, char** argv)
{
  // Check the number of parameters
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " [folder]  [timeResol]  [usePromptMuons] [ptReweighting] [outputFile]" << std::endl;
    return 1;
  }
  
  string folderName(argv[1]);
  string timeResolution(argv[2]);
  int usePromptMuons = atoi(argv[3]); 
  int applyPtReweighting = atoi(argv[4]); 

  TH2F *hratio2D;
  if (applyPtReweighting){
    TFile *fWeights;
    if (folderName.find("noPU")!= std::string::npos){
      fWeights= TFile::Open("ptRatio_mu_noPU_TTbar.root");
      if ( folderName.find("QCD")!= std::string::npos )
	fWeights= TFile::Open("ptRatio_mu_noPU_QCD.root");
    }
    else{
      fWeights= TFile::Open("ptRatio_mu_TTbar.root");
      if ( folderName.find("QCD")!= std::string::npos )
	fWeights= TFile::Open("ptRatio_mu_QCD.root");
    }
    std::cout << "Applying pT/eta weights from " << fWeights->GetName()<<std::endl;
    hratio2D = (TH2F*)fWeights->Get("hratio2D");
  }
  

  cout <<endl;
  cout << "Start analyzing " << folderName <<endl;

  // -- get TChain
  TChain* chain = new TChain(("analysis/tree_"+timeResolution+"ps").c_str(),"tree");
  chain->Add((folderName+"/0000/muonIsolation_*.root").c_str());

  // -- tree vars
  int npu;
  vector<float> *muon_pt;
  vector<float> *muon_eta;
  vector<float> *muon_phi;
  vector<float> *muon_dz3D;
  vector<float> *muon_dxy3D;
  vector<float> *muon_dz4D;
  vector<float> *muon_dxy4D;
  vector<int  > *muon_isLoose;
  vector<int>   *muon_isPrompt;
  vector<int>   *muon_isTight3D;
  vector<int>   *muon_isTight4D;
  vector<int>   *muon_isMatchedToGenJet;
  vector<float> *muon_chIso02;
  vector<float> *muon_chIso02_dT;
  vector<float> *muon_chIso03;
  vector<float> *muon_chIso03_dT;
  vector<float> *muon_chIso04;
  vector<float> *muon_chIso04_dT;
  vector<float> *muon_chIso05;
  vector<float> *muon_chIso05_dT;
  
  muon_pt = 0;
  muon_eta = 0;
  muon_phi = 0;
  muon_dz3D = 0;
  muon_dxy3D = 0;
  muon_dz4D = 0;
  muon_dxy4D = 0;
  muon_isLoose = 0;
  muon_isTight3D = 0;
  muon_isTight4D = 0;
  muon_isPrompt = 0;
  muon_isMatchedToGenJet = 0;
  muon_chIso02 = 0;
  muon_chIso02_dT = 0;
  muon_chIso03 = 0;
  muon_chIso03_dT = 0;
  muon_chIso04 = 0;
  muon_chIso04_dT = 0;
  muon_chIso05 = 0;
  muon_chIso05_dT = 0;


  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("npu",1);                     chain->SetBranchAddress("npu",              &npu);
  chain->SetBranchStatus("muon_pt",1);                 chain->SetBranchAddress("muon_pt",        &muon_pt);
  chain->SetBranchStatus("muon_eta",1);                chain->SetBranchAddress("muon_eta",       &muon_eta);
  chain->SetBranchStatus("muon_phi",1);                chain->SetBranchAddress("muon_phi",       &muon_phi);
  chain->SetBranchStatus("muon_dz3D",1);                chain->SetBranchAddress("muon_dz3D",       &muon_dz3D);
  chain->SetBranchStatus("muon_dxy3D",1);                chain->SetBranchAddress("muon_dxy3D",       &muon_dxy3D);
  chain->SetBranchStatus("muon_dz4D",1);                chain->SetBranchAddress("muon_dz4D",       &muon_dz4D);
  chain->SetBranchStatus("muon_dxy4D",1);                chain->SetBranchAddress("muon_dxy4D",       &muon_dxy4D);
  chain->SetBranchStatus("muon_isLoose",1);            chain->SetBranchAddress("muon_isLoose",       &muon_isLoose);
  chain->SetBranchStatus("muon_isTight3D",1);            chain->SetBranchAddress("muon_isTight3D",       &muon_isTight3D);
  chain->SetBranchStatus("muon_isTight4D",1);            chain->SetBranchAddress("muon_isTight4D",       &muon_isTight4D);
  chain->SetBranchStatus("muon_isPrompt",1);           chain->SetBranchAddress("muon_isPrompt",       &muon_isPrompt);
  chain->SetBranchStatus("muon_isMatchedToGenJet",1);  chain->SetBranchAddress("muon_isMatchedToGenJet",  &muon_isMatchedToGenJet);
  chain->SetBranchStatus("muon_chIso02",1);            chain->SetBranchAddress("muon_chIso02",   &muon_chIso02);
  chain->SetBranchStatus("muon_chIso02_dT",1);         chain->SetBranchAddress("muon_chIso02_dT",&muon_chIso02_dT);
  chain->SetBranchStatus("muon_chIso03",1);            chain->SetBranchAddress("muon_chIso03",   &muon_chIso03);
  chain->SetBranchStatus("muon_chIso03_dT",1);         chain->SetBranchAddress("muon_chIso03_dT",&muon_chIso03_dT);
  chain->SetBranchStatus("muon_chIso04",1);            chain->SetBranchAddress("muon_chIso04",   &muon_chIso04);
  chain->SetBranchStatus("muon_chIso04_dT",1);         chain->SetBranchAddress("muon_chIso04_dT",&muon_chIso04_dT);
  chain->SetBranchStatus("muon_chIso05",1);            chain->SetBranchAddress("muon_chIso05",   &muon_chIso05);
  chain->SetBranchStatus("muon_chIso05_dT",1);         chain->SetBranchAddress("muon_chIso05_dT",&muon_chIso05_dT);


  // -- book histograms
  TH1F *h_npu   = new TH1F("h_npu","h_npu",200,50,250);
  h_npu->GetXaxis()->SetTitle("number of pileup vertices");
  
  TH1F *h_muon_pt = new TH1F("h_muon_pt","h_muon_pt",200,0,200);
  h_muon_pt->GetXaxis()->SetTitle("muon p_{T} (GeV)");
  
  TH1F *h_muon_eta = new TH1F("h_muon_eta","h_muon_eta",100,-3,3);
  h_muon_eta->GetXaxis()->SetTitle("muon eta");

  TH1F *h_muon_phi = new TH1F("h_muon_phi","h_muon_phi",100,-4,4);
  h_muon_phi->GetXaxis()->SetTitle("muon phi");

  TH2F *h2_muon_pt_vs_eta = new TH2F("h2_muon_pt_vs_eta","h2_muon_pt_vs_eta",150, -3, 3, 500, 0, 1000);
  
  // -- relative isolation
  TH1F *h_muon_relChIso02 = new TH1F("h_muon_relChIso02","h_muon_relChIso02",5000,0,5);
  h_muon_relChIso02->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso02_dT = new TH1F("h_muon_relChIso02_dT","h_muon_relChIso02_dT",5000,0,5);
  h_muon_relChIso02_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03 = new TH1F("h_muon_relChIso03","h_muon_relChIso03",5000,0,5);
  h_muon_relChIso03->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT = new TH1F("h_muon_relChIso03_dT","h_muon_relChIso03_dT",5000,0,5);
  h_muon_relChIso03_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04 = new TH1F("h_muon_relChIso04","h_muon_relChIso04",5000,0,5);
  h_muon_relChIso04->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04_dT = new TH1F("h_muon_relChIso04_dT","h_muon_relChIso04_dT",5000,0,5);
  h_muon_relChIso04_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05 = new TH1F("h_muon_relChIso05","h_muon_relChIso05",5000,0,5);
  h_muon_relChIso05->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05_dT = new TH1F("h_muon_relChIso05_dT","h_muon_relChIso05_dT",5000,0,5);
  h_muon_relChIso05_dT->GetXaxis()->SetTitle("relative charged isolation");


  // -- relative isolation barrel
  TH1F *h_muon_relChIso02_barrel = new TH1F("h_muon_relChIso02_barrel","h_muon_relChIso02_barrel",5000,0,5);
  h_muon_relChIso02_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso02_dT_barrel = new TH1F("h_muon_relChIso02_dT_barrel","h_muon_relChIso02_dT_barrel",5000,0,5);
  h_muon_relChIso02_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_barrel = new TH1F("h_muon_relChIso03_barrel","h_muon_relChIso03_barrel",5000,0,5);
  h_muon_relChIso03_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_barrel = new TH1F("h_muon_relChIso03_dT_barrel","h_muon_relChIso03_dT_barrel",5000,0,5);
  h_muon_relChIso03_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04_barrel = new TH1F("h_muon_relChIso04_barrel","h_muon_relChIso04_barrel",5000,0,5);
  h_muon_relChIso04_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04_dT_barrel = new TH1F("h_muon_relChIso04_dT_barrel","h_muon_relChIso04_dT_barrel",5000,0,5);
  h_muon_relChIso04_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05_barrel = new TH1F("h_muon_relChIso05_barrel","h_muon_relChIso05_barrel",5000,0,5);
  h_muon_relChIso05_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05_dT_barrel = new TH1F("h_muon_relChIso05_dT_barrel","h_muon_relChIso05_dT_barrel",5000,0,5);
  h_muon_relChIso05_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  // -- relative isolation endcap
  TH1F *h_muon_relChIso02_endcap = new TH1F("h_muon_relChIso02_endcap","h_muon_relChIso02_endcap",5000,0,5);
  h_muon_relChIso02_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso02_dT_endcap = new TH1F("h_muon_relChIso02_dT_endcap","h_muon_relChIso02_dT_endcap",5000,0,5);
  h_muon_relChIso02_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_endcap = new TH1F("h_muon_relChIso03_endcap","h_muon_relChIso03_endcap",5000,0,5);
  h_muon_relChIso03_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_endcap = new TH1F("h_muon_relChIso03_dT_endcap","h_muon_relChIso03_dT_endcap",5000,0,5);
  h_muon_relChIso03_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04_endcap = new TH1F("h_muon_relChIso04_endcap","h_muon_relChIso04_endcap",5000,0,5);
  h_muon_relChIso04_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso04_dT_endcap = new TH1F("h_muon_relChIso04_dT_endcap","h_muon_relChIso04_dT_endcap",5000,0,5);
  h_muon_relChIso04_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05_endcap = new TH1F("h_muon_relChIso05_endcap","h_muon_relChIso05_endcap",5000,0,5);
  h_muon_relChIso05_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso05_dT_endcap = new TH1F("h_muon_relChIso05_dT_endcap","h_muon_relChIso05_dT_endcap",5000,0,5);
  h_muon_relChIso05_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");


  //
  TH1F *h_muon_relChIso03_diff = new TH1F("h_muon_relChIso03_diff","h_muon_relChIso03_diff",1000,0,1);
  h_muon_relChIso03_diff->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_diff_barrel = new TH1F("h_muon_relChIso03_diff_barrel","h_muon_relChIso03_diff_barrel",1000,0,1);
  h_muon_relChIso03_diff_barrel->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_diff_endcap = new TH1F("h_muon_relChIso03_diff_endcap","h_muon_relChIso03_diff_endcap",1000,0,1);
  h_muon_relChIso03_diff_endcap ->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");


  // -- loop over events
  //int maxEntries = std::min(int(chain ->GetEntries()),1000000);  
  int maxEntries = chain ->GetEntries();
  cout << "Number of events to be analyzed : " << maxEntries <<endl;
  float w = 1;
  float pt = 0;
  float maxdz = 0.5;
  float maxdxy = 0.2;

  for (int i = 0; i < maxEntries; i++) {
    
    chain -> GetEntry(i);
    
    if (i%1000==0) cout << "Analyzing event " << i << "\r" << flush;
    
    h_npu->Fill(npu);

    for (unsigned int imu = 0; imu < muon_pt->size(); imu++){
      
      pt = muon_pt->at(imu);
      
      if ( pt < 20. )  continue;
      if ( muon_isLoose->at(imu) == 0 ) continue;


      // -- pt/eta reweighting
      if (applyPtReweighting){
	int bin = hratio2D->FindBin(muon_eta->at(imu), pt ); 
	w = hratio2D -> GetBinContent(bin); 
      }
      else {
	w = 1;
      }

      
      // -- prompt muons or non-prompt muons
      bool pass = false;
      if (usePromptMuons && muon_isPrompt->at(imu)) pass = true; 
      if (!usePromptMuons && (!muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu)) ) pass = true;

      if (!pass) continue;

      
      h_muon_pt->Fill(pt, w);
      h_muon_eta->Fill(muon_eta->at(imu), w);
      h_muon_phi->Fill(muon_phi->at(imu), w);
      h2_muon_pt_vs_eta->Fill(muon_eta->at(imu), pt, w);

      // control plots
      if ( fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy && fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy){
	float chIsoDiff = ( muon_chIso03->at(imu) - muon_chIso03_dT->at(imu) ) / muon_chIso03->at(imu);
	if ( muon_chIso03->at(imu) == 0 && chIsoDiff ) chIsoDiff = 0.; 
	h_muon_relChIso03_diff -> Fill( chIsoDiff, w );
	if (fabs(muon_eta->at(imu))<1.5){ 
	  h_muon_relChIso03_diff_barrel -> Fill( chIsoDiff, w );
	}
        else{ 
	  h_muon_relChIso03_diff_endcap -> Fill( chIsoDiff, w );
	}
      }

      // all
      if ( fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy){
	h_muon_relChIso02    -> Fill(muon_chIso02->at(imu)/pt, w);
	h_muon_relChIso03    -> Fill(muon_chIso03->at(imu)/pt, w);
	h_muon_relChIso04    -> Fill(muon_chIso04->at(imu)/pt, w);
	h_muon_relChIso05    -> Fill(muon_chIso05->at(imu)/pt, w);
      }
      if ( fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy){
	h_muon_relChIso02_dT -> Fill(muon_chIso02_dT->at(imu)/pt, w);
	h_muon_relChIso03_dT -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	h_muon_relChIso04_dT -> Fill(muon_chIso04_dT->at(imu)/pt, w);
	h_muon_relChIso05_dT -> Fill(muon_chIso05_dT->at(imu)/pt, w);
      }

      // barrel 
      if (fabs(muon_eta->at(imu))<1.5){
	if ( fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy){
	  h_muon_relChIso02_barrel    -> Fill(muon_chIso02->at(imu)/pt, w);
	  h_muon_relChIso03_barrel    -> Fill(muon_chIso03->at(imu)/pt, w);
	  h_muon_relChIso04_barrel    -> Fill(muon_chIso04->at(imu)/pt, w);
	  h_muon_relChIso05_barrel    -> Fill(muon_chIso05->at(imu)/pt, w);
	}
	if ( fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy ){
	  h_muon_relChIso02_dT_barrel -> Fill(muon_chIso02_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dT_barrel -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	  h_muon_relChIso04_dT_barrel -> Fill(muon_chIso04_dT->at(imu)/pt, w);
	  h_muon_relChIso05_dT_barrel -> Fill(muon_chIso05_dT->at(imu)/pt, w);
	}
      }
      // endcap
      else{
	if ( fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy){
	  h_muon_relChIso02_endcap    -> Fill(muon_chIso02->at(imu)/pt, w);
	  h_muon_relChIso03_endcap    -> Fill(muon_chIso03->at(imu)/pt, w);
	  h_muon_relChIso04_endcap    -> Fill(muon_chIso04->at(imu)/pt, w);
	  h_muon_relChIso05_endcap    -> Fill(muon_chIso05->at(imu)/pt, w);
	}
	if ( fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy ){
	  h_muon_relChIso02_dT_endcap -> Fill(muon_chIso02_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dT_endcap -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	  h_muon_relChIso04_dT_endcap -> Fill(muon_chIso04_dT->at(imu)/pt, w);
	  h_muon_relChIso05_dT_endcap -> Fill(muon_chIso05_dT->at(imu)/pt, w);
	}
      }
        
    }// end loop over muons
  }// end loop over events

  cout << endl;

  cout << "Saving histograms in file " << argv[5] << endl;

  // -- save output file
  TFile *fout = new TFile(argv[5],"recreate");
  h_npu->Write();
  h_muon_pt->Write();
  h_muon_eta->Write();
  h_muon_phi->Write();
  h2_muon_pt_vs_eta->Write();

  h_muon_relChIso03_diff ->Write();
  h_muon_relChIso03_diff_barrel ->Write();
  h_muon_relChIso03_diff_endcap ->Write();

  h_muon_relChIso02->Write();
  h_muon_relChIso02_barrel->Write();
  h_muon_relChIso02_endcap->Write();
  h_muon_relChIso02_dT->Write();
  h_muon_relChIso02_dT_barrel->Write();
  h_muon_relChIso02_dT_endcap->Write();

  h_muon_relChIso03->Write();
  h_muon_relChIso03_barrel->Write();
  h_muon_relChIso03_endcap->Write();
  h_muon_relChIso03_dT->Write();
  h_muon_relChIso03_dT_barrel->Write();
  h_muon_relChIso03_dT_endcap->Write();

  h_muon_relChIso04->Write();
  h_muon_relChIso04_barrel->Write();
  h_muon_relChIso04_endcap->Write();
  h_muon_relChIso04_dT->Write();
  h_muon_relChIso04_dT_barrel->Write();
  h_muon_relChIso04_dT_endcap->Write();

  h_muon_relChIso05->Write();
  h_muon_relChIso05_barrel->Write();
  h_muon_relChIso05_endcap->Write();
  h_muon_relChIso05_dT->Write();
  h_muon_relChIso05_dT_barrel->Write();
  h_muon_relChIso05_dT_endcap->Write();

  fout->Close();
  
}


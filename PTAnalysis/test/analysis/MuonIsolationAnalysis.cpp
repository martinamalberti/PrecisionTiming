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
  float vtxGen_t;
  float vtxGen_z;
  float vtx3D_z;
  float vtx4D_z;
  float vtx4D_t;
  float vtx3D_zErr;
  float vtx4D_zErr;
  float vtx4D_tErr;
  int vtx3D_isFake;
  int vtx4D_isFake;
  vector<float> *muon_pt;
  vector<float> *muon_eta;
  vector<float> *muon_phi;
  vector<float> *muon_dz3D;
  vector<float> *muon_dxy3D;
  vector<float> *muon_dz4D;
  vector<float> *muon_dxy4D;
  vector<int  > *muon_isLoose;
  vector<int>   *muon_isTight3D;
  vector<int>   *muon_isTight4D;
  vector<int>   *muon_isPrompt;
  vector<int>   *muon_isMatchedToGenJet;
  vector<int>   *muon_isFromTauDecay;

  vector<float> *muon_chIso03;
  vector<float> *muon_chIso03_dT;

  vector<float> *muon_chIso03_dZmu;
  vector<float> *muon_chIso03_dZmu_dTmu;

  vector<float> *muon_chIso03_simVtx;
  vector<float> *muon_chIso03_dT_simVtx;

  vector<float> *muon_chIso03_dZ05;
  vector<float> *muon_chIso03_dZ2;
  vector<float> *muon_chIso03_dZ05_dT;
  vector<float> *muon_chIso03_dZ2_dT;
  vector<float> *muon_chIso03_reldZ;
  vector<float> *muon_chIso03_reldZ_dT;  

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
  muon_isFromTauDecay = 0;
  muon_chIso03 = 0;
  muon_chIso03_dT = 0;
  muon_chIso03_dZmu = 0;
  muon_chIso03_dZmu_dTmu = 0;
  muon_chIso03_simVtx = 0;
  muon_chIso03_dT_simVtx = 0;
  muon_chIso03_dZ05 = 0;
  muon_chIso03_dZ2 = 0;
  muon_chIso03_dZ05_dT = 0;
  muon_chIso03_dZ2_dT = 0;
  muon_chIso03_reldZ = 0;
  muon_chIso03_reldZ_dT = 0;

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("npu",1);                     chain->SetBranchAddress("npu",            &npu);

  chain->SetBranchStatus("vtxGen_z",1);                chain->SetBranchAddress("vtxGen_z",       &vtxGen_z);
  chain->SetBranchStatus("vtxGen_t",1);                chain->SetBranchAddress("vtxGen_t",       &vtxGen_t);
  chain->SetBranchStatus("vtx3D_isFake",1);            chain->SetBranchAddress("vtx3D_isFake",   &vtx3D_isFake);
  chain->SetBranchStatus("vtx4D_isFake",1);            chain->SetBranchAddress("vtx4D_isFake",   &vtx4D_isFake);
  chain->SetBranchStatus("vtx3D_z",1);                 chain->SetBranchAddress("vtx3D_z",        &vtx3D_z);
  chain->SetBranchStatus("vtx4D_z",1);                 chain->SetBranchAddress("vtx4D_z",        &vtx4D_z);
  chain->SetBranchStatus("vtx4D_t",1);                 chain->SetBranchAddress("vtx4D_t",        &vtx4D_t);
  chain->SetBranchStatus("vtx3D_zErr",1);              chain->SetBranchAddress("vtx3D_zErr",     &vtx3D_zErr);
  chain->SetBranchStatus("vtx4D_zErr",1);              chain->SetBranchAddress("vtx4D_zErr",     &vtx4D_zErr);
  chain->SetBranchStatus("vtx4D_tErr",1);              chain->SetBranchAddress("vtx4D_tErr",     &vtx4D_tErr);

  chain->SetBranchStatus("muon_pt",1);                 chain->SetBranchAddress("muon_pt",        &muon_pt);
  chain->SetBranchStatus("muon_eta",1);                chain->SetBranchAddress("muon_eta",       &muon_eta);
  chain->SetBranchStatus("muon_phi",1);                chain->SetBranchAddress("muon_phi",       &muon_phi);
  chain->SetBranchStatus("muon_dz3D",1);               chain->SetBranchAddress("muon_dz3D",      &muon_dz3D);
  chain->SetBranchStatus("muon_dxy3D",1);              chain->SetBranchAddress("muon_dxy3D",     &muon_dxy3D);
  chain->SetBranchStatus("muon_dz4D",1);               chain->SetBranchAddress("muon_dz4D",      &muon_dz4D);
  chain->SetBranchStatus("muon_dxy4D",1);              chain->SetBranchAddress("muon_dxy4D",     &muon_dxy4D);
  chain->SetBranchStatus("muon_isLoose",1);            chain->SetBranchAddress("muon_isLoose",   &muon_isLoose);
  chain->SetBranchStatus("muon_isTight3D",1);          chain->SetBranchAddress("muon_isTight3D", &muon_isTight3D);
  chain->SetBranchStatus("muon_isTight4D",1);          chain->SetBranchAddress("muon_isTight4D", &muon_isTight4D);
  chain->SetBranchStatus("muon_isPrompt",1);           chain->SetBranchAddress("muon_isPrompt",  &muon_isPrompt);
  chain->SetBranchStatus("muon_isMatchedToGenJet",1);  chain->SetBranchAddress("muon_isMatchedToGenJet",  &muon_isMatchedToGenJet);
  chain->SetBranchStatus("muon_isFromTauDecay",1);     chain->SetBranchAddress("muon_isFromTauDecay",  &muon_isFromTauDecay);

  chain->SetBranchStatus("muon_chIso03",1);            chain->SetBranchAddress("muon_chIso03",   &muon_chIso03);
  chain->SetBranchStatus("muon_chIso03_dT",1);         chain->SetBranchAddress("muon_chIso03_dT",&muon_chIso03_dT);

  chain->SetBranchStatus("muon_chIso03_dZmu",1);       chain->SetBranchAddress("muon_chIso03_dZmu",&muon_chIso03_dZmu);
  chain->SetBranchStatus("muon_chIso03_dZmu_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu_dTmu",&muon_chIso03_dZmu_dTmu);

  chain->SetBranchStatus("muon_chIso03_simVtx",1);     chain->SetBranchAddress("muon_chIso03_simVtx",   &muon_chIso03_simVtx);
  chain->SetBranchStatus("muon_chIso03_dT_simVtx",1);  chain->SetBranchAddress("muon_chIso03_dT_simVtx",&muon_chIso03_dT_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ05",1);       chain->SetBranchAddress("muon_chIso03_dZ05",   &muon_chIso03_dZ05);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT",1);    chain->SetBranchAddress("muon_chIso03_dZ05_dT",&muon_chIso03_dZ05_dT);

  chain->SetBranchStatus("muon_chIso03_dZ2",1);       chain->SetBranchAddress("muon_chIso03_dZ2",   &muon_chIso03_dZ2);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT",1);    chain->SetBranchAddress("muon_chIso03_dZ2_dT",&muon_chIso03_dZ2_dT);

  chain->SetBranchStatus("muon_chIso03_reldZ",1);       chain->SetBranchAddress("muon_chIso03_reldZ",   &muon_chIso03_reldZ);
  chain->SetBranchStatus("muon_chIso03_reldZ_dT",1);    chain->SetBranchAddress("muon_chIso03_reldZ_dT",&muon_chIso03_reldZ_dT);


  // -- book histograms
  TH1F *h_npu   = new TH1F("h_npu","h_npu",200,50,250);
  h_npu->GetXaxis()->SetTitle("number of pileup vertices");

  TH1F *h_vtx_dz3D = new TH1F("h_vtx_dz3D","h_vtx_dz3D",1000, -0.2, 0.2);
  h_vtx_dz3D -> GetXaxis()->SetTitle("z_{3Dvtx} - z_{gen} (mm)");

  TH1F *h_vtx_dz4D = new TH1F("h_vtx_dz4D","h_vtx_dz4D",1000, -0.2, 0.2);
  h_vtx_dz4D -> GetXaxis()->SetTitle("z_{4Dvtx} - z_{gen} (mm)");

  TH1F *h_vtx_dt4D = new TH1F("h_vtx_dt4D","h_vtx_dt4D",1000, -0.5, 0.5);
  h_vtx_dt4D -> GetXaxis()->SetTitle("t_{4Dvtx} - t_{gen} (ns)");

  TH1F *h_vtx_dz3D_pull = new TH1F("h_vtx_dz3D_pull","h_vtx_dz3D_pull",200, -10, 10);
  h_vtx_dz3D_pull -> GetXaxis()->SetTitle("(z_{3Dvtx} - z_{gen})/#sigma_{z}");

  TH1F *h_vtx_dz4D_pull = new TH1F("h_vtx_dz4D_pull","h_vtx_dz4D_pull",200, -10, 10);
  h_vtx_dz4D_pull -> GetXaxis()->SetTitle("(z_{4Dvtx} - z_{gen})/#sigma_{z}");

  TH1F *h_vtx_dt4D_pull = new TH1F("h_vtx_dt4D_pull","h_vtx_dt4D_pull",200, -10, 10);
  h_vtx_dt4D_pull -> GetXaxis()->SetTitle("(t_{4Dvtx} - t_{gen})/#sigma_{t}");


  TH1F *h_muon_pt = new TH1F("h_muon_pt","h_muon_pt",200,0,200);
  h_muon_pt->GetXaxis()->SetTitle("muon p_{T} (GeV)");
  
  TH1F *h_muon_eta = new TH1F("h_muon_eta","h_muon_eta",100,-3,3);
  h_muon_eta->GetXaxis()->SetTitle("muon eta");

  TH1F *h_muon_phi = new TH1F("h_muon_phi","h_muon_phi",100,-4,4);
  h_muon_phi->GetXaxis()->SetTitle("muon phi");

  TH2F *h2_muon_pt_vs_eta = new TH2F("h2_muon_pt_vs_eta","h2_muon_pt_vs_eta",150, -3, 3, 500, 0, 1000);
  
  // -- relative isolation -- dZ = 1 mm
  TH1F *h_muon_relChIso03 = new TH1F("h_muon_relChIso03","h_muon_relChIso03",5000,0,5);
  h_muon_relChIso03->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT = new TH1F("h_muon_relChIso03_dT","h_muon_relChIso03_dT",5000,0,5);
  h_muon_relChIso03_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_barrel = new TH1F("h_muon_relChIso03_barrel","h_muon_relChIso03_barrel",5000,0,5);
  h_muon_relChIso03_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_barrel = new TH1F("h_muon_relChIso03_dT_barrel","h_muon_relChIso03_dT_barrel",5000,0,5);
  h_muon_relChIso03_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_endcap = new TH1F("h_muon_relChIso03_endcap","h_muon_relChIso03_endcap",5000,0,5);
  h_muon_relChIso03_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_endcap = new TH1F("h_muon_relChIso03_dT_endcap","h_muon_relChIso03_dT_endcap",5000,0,5);
  h_muon_relChIso03_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // chIso wrt to sim vertex
  TH1F *h_muon_relChIso03_simVtx = new TH1F("h_muon_relChIso03_simVtx","h_muon_relChIso03_simVtx",5000,0,5);
  h_muon_relChIso03_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_simVtx = new TH1F("h_muon_relChIso03_dT_simVtx","h_muon_relChIso03_dT_simVtx",5000,0,5);
  h_muon_relChIso03_dT_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_simVtx_barrel = new TH1F("h_muon_relChIso03_simVtx_barrel","h_muon_relChIso03_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_simVtx_barrel = new TH1F("h_muon_relChIso03_dT_simVtx_barrel","h_muon_relChIso03_dT_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dT_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_simVtx_endcap = new TH1F("h_muon_relChIso03_simVtx_endcap","h_muon_relChIso03_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dT_simVtx_endcap = new TH1F("h_muon_relChIso03_dT_simVtx_endcap","h_muon_relChIso03_dT_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dT_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // dZ = 0.5 mm
  TH1F *h_muon_relChIso03_dZ05 = new TH1F("h_muon_relChIso03_dZ05","h_muon_relChIso03_dZ05",5000,0,5);
  h_muon_relChIso03_dZ05->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT = new TH1F("h_muon_relChIso03_dZ05_dT","h_muon_relChIso03_dZ05_dT",5000,0,5);
  h_muon_relChIso03_dZ05_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_barrel = new TH1F("h_muon_relChIso03_dZ05_barrel","h_muon_relChIso03_dZ05_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT_barrel = new TH1F("h_muon_relChIso03_dZ05_dT_barrel","h_muon_relChIso03_dZ05_dT_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_endcap = new TH1F("h_muon_relChIso03_dZ05_endcap","h_muon_relChIso03_dZ05_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT_endcap = new TH1F("h_muon_relChIso03_dZ05_dT_endcap","h_muon_relChIso03_dZ05_dT_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // dZ = 2.0 mm
  TH1F *h_muon_relChIso03_dZ2 = new TH1F("h_muon_relChIso03_dZ2","h_muon_relChIso03_dZ2",5000,0,5);
  h_muon_relChIso03_dZ2->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT = new TH1F("h_muon_relChIso03_dZ2_dT","h_muon_relChIso03_dZ2_dT",5000,0,5);
  h_muon_relChIso03_dZ2_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_barrel = new TH1F("h_muon_relChIso03_dZ2_barrel","h_muon_relChIso03_dZ2_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT_barrel = new TH1F("h_muon_relChIso03_dZ2_dT_barrel","h_muon_relChIso03_dZ2_dT_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_endcap = new TH1F("h_muon_relChIso03_dZ2_endcap","h_muon_relChIso03_dZ2_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT_endcap = new TH1F("h_muon_relChIso03_dZ2_dT_endcap","h_muon_relChIso03_dZ2_dT_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

  // rel dZ ( 3*sigma_Z)
  TH1F *h_muon_relChIso03_reldZ = new TH1F("h_muon_relChIso03_reldZ","h_muon_relChIso03_reldZ",5000,0,5);
  h_muon_relChIso03_reldZ->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_reldZ_dT = new TH1F("h_muon_relChIso03_reldZ_dT","h_muon_relChIso03_reldZ_dT",5000,0,5);
  h_muon_relChIso03_reldZ_dT->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_reldZ_barrel = new TH1F("h_muon_relChIso03_reldZ_barrel","h_muon_relChIso03_reldZ_barrel",5000,0,5);
  h_muon_relChIso03_reldZ_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_reldZ_dT_barrel = new TH1F("h_muon_relChIso03_reldZ_dT_barrel","h_muon_relChIso03_reldZ_dT_barrel",5000,0,5);
  h_muon_relChIso03_reldZ_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_reldZ_endcap = new TH1F("h_muon_relChIso03_reldZ_endcap","h_muon_relChIso03_reldZ_endcap",5000,0,5);
  h_muon_relChIso03_reldZ_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_reldZ_dT_endcap = new TH1F("h_muon_relChIso03_reldZ_dT_endcap","h_muon_relChIso03_reldZ_dT_endcap",5000,0,5);
  h_muon_relChIso03_reldZ_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // dz, dt wrt muon
  TH1F *h_muon_relChIso03_dZmu = new TH1F("h_muon_relChIso03_dZmu","h_muon_relChIso03_dZmu",5000,0,5);
  h_muon_relChIso03_dZmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu_dTmu = new TH1F("h_muon_relChIso03_dZmu_dTmu","h_muon_relChIso03_dZmu_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu_barrel = new TH1F("h_muon_relChIso03_dZmu_barrel","h_muon_relChIso03_dZmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu_dTmu_barrel","h_muon_relChIso03_dZmu_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu_endcap = new TH1F("h_muon_relChIso03_dZmu_endcap","h_muon_relChIso03_dZmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu_dTmu_endcap","h_muon_relChIso03_dZmu_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // ratio iso
  TH1F *h_muon_relChIso03_ratio = new TH1F("h_muon_relChIso03_ratio","h_muon_relChIso03_ratio",1000,0,2);
  h_muon_relChIso03_ratio->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_ratio_barrel = new TH1F("h_muon_relChIso03_ratio_barrel","h_muon_relChIso03_ratio_barrel",1000,0,2);
  h_muon_relChIso03_ratio_barrel->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_ratio_endcap = new TH1F("h_muon_relChIso03_ratio_endcap","h_muon_relChIso03_ratio_endcap",1000,0,2);
  h_muon_relChIso03_ratio_endcap->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");





  // -- loop over events
  //int maxEntries = std::min(int(chain ->GetEntries()),1000000);  
  int maxEntries = chain ->GetEntries();
  cout << "Number of events to be analyzed : " << maxEntries <<endl;
  float w = 1;
  float pt = 0;
  float maxdz = 0.5;
  float maxdxy = 0.2;

  int nMuonsInEvent = 0;

  for (int i = 0; i < maxEntries; i++) {
    
    chain -> GetEntry(i);
    
    if (i%1000==0) cout << "Analyzing event " << i << "\r" << flush;
    
    h_npu->Fill(npu);

    nMuonsInEvent = 0;

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
      if (!usePromptMuons && (!muon_isPrompt->at(imu) && muon_isMatchedToGenJet->at(imu) && !muon_isFromTauDecay->at(imu)) ) pass = true;

      if (!pass) continue;

      nMuonsInEvent++;
      
      bool pass3D = fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy && !vtx3D_isFake;
      bool pass4D = fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy && !vtx4D_isFake;

      if (nMuonsInEvent == 1){

      	if (pass3D) h_vtx_dz3D->Fill( vtx3D_z - vtxGen_z);
      	if (pass4D) h_vtx_dz4D->Fill( vtx4D_z - vtxGen_z);
      	if (pass4D) h_vtx_dt4D->Fill( vtx4D_t - vtxGen_t*1000000000.);

      	if (pass3D) h_vtx_dz3D_pull->Fill( (vtx3D_z - vtxGen_z) / vtx3D_zErr );
	if (pass4D) h_vtx_dz4D_pull->Fill( (vtx4D_z - vtxGen_z) / vtx4D_zErr );
      	if (pass4D) h_vtx_dt4D_pull->Fill( (vtx4D_t - vtxGen_t*1000000000.) / vtx4D_tErr);

      }

      h_muon_pt->Fill(pt, w);
      h_muon_eta->Fill(muon_eta->at(imu), w);
      h_muon_phi->Fill(muon_phi->at(imu), w);
      h2_muon_pt_vs_eta->Fill(muon_eta->at(imu), pt, w);


      // control plots
      if (pass3D && pass4D){	
	float chIsoRatio = muon_chIso03_dT->at(imu)/muon_chIso03->at(imu);
	if ( muon_chIso03->at(imu) == 0 ) chIsoRatio = 1;
	h_muon_relChIso03_ratio -> Fill( chIsoRatio, w );
	if (fabs(muon_eta->at(imu))<1.5){ 
	  h_muon_relChIso03_ratio_barrel -> Fill( chIsoRatio, w );
	}
        else{ 
	  h_muon_relChIso03_ratio_endcap -> Fill( chIsoRatio, w );
	}
      }
      
      // chiIso plots - all
      if ( pass3D ){
	h_muon_relChIso03 -> Fill(muon_chIso03->at(imu)/pt, w);
	h_muon_relChIso03_simVtx -> Fill(muon_chIso03_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ05 -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
	h_muon_relChIso03_dZ2 -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
	h_muon_relChIso03_reldZ -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
      }
      if ( pass4D ){
	h_muon_relChIso03_dT -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	h_muon_relChIso03_dT_simVtx -> Fill(muon_chIso03_dT_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ05_dT -> Fill(muon_chIso03_dZ05_dT->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_dT -> Fill(muon_chIso03_dZ2_dT->at(imu)/pt, w);
	h_muon_relChIso03_reldZ_dT -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
      }
      // barrel 
      if (fabs(muon_eta->at(imu))<1.5){
	if ( pass3D ) {
	  h_muon_relChIso03_barrel    -> Fill(muon_chIso03->at(imu)/pt, w);
	  h_muon_relChIso03_simVtx_barrel -> Fill(muon_chIso03_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_barrel -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_barrel -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
	  h_muon_relChIso03_reldZ_barrel -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
	}
	if ( pass4D ){
	  h_muon_relChIso03_dT_barrel -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dT_simVtx_barrel -> Fill(muon_chIso03_dT_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT_barrel -> Fill(muon_chIso03_dZ05_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT_barrel -> Fill(muon_chIso03_dZ2_dT->at(imu)/pt, w);
	  h_muon_relChIso03_reldZ_dT_barrel -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
	}
      }
      // endcap
      else{
	if ( pass3D ){
	  h_muon_relChIso03_endcap    -> Fill(muon_chIso03->at(imu)/pt, w);
	  h_muon_relChIso03_simVtx_endcap -> Fill(muon_chIso03_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_endcap -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_endcap -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
	  h_muon_relChIso03_reldZ_endcap -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
	}
	if ( pass4D ) {
	  h_muon_relChIso03_dT_endcap -> Fill(muon_chIso03_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dT_simVtx_endcap -> Fill(muon_chIso03_dT_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT_endcap -> Fill(muon_chIso03_dZ05_dT->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT_endcap -> Fill(muon_chIso03_dZ2_dT->at(imu)/pt, w);
	  h_muon_relChIso03_reldZ_dT_endcap -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
	}
      }


      // dz, dt cuts wrt muon
      h_muon_relChIso03_dZmu -> Fill(muon_chIso03_dZmu->at(imu)/pt, w);
      h_muon_relChIso03_dZmu_dTmu -> Fill(muon_chIso03_dZmu_dTmu->at(imu)/pt, w);
      if (fabs(muon_eta->at(imu))<1.5){ 
	h_muon_relChIso03_dZmu_barrel -> Fill(muon_chIso03_dZmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu_dTmu_barrel -> Fill(muon_chIso03_dZmu_dTmu->at(imu)/pt, w);
      }
      else{
	h_muon_relChIso03_dZmu_endcap -> Fill(muon_chIso03_dZmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu_dTmu_endcap -> Fill(muon_chIso03_dZmu_dTmu->at(imu)/pt, w);
      }
        
    }// end loop over muons
  }// end loop over events

  cout << endl;

  TF1 *f1 = new TF1("f1","gaus",-10., 10.);
  f1->SetParameter(2, h_vtx_dz3D_pull->GetRMS() );
  f1->SetRange(h_vtx_dz3D_pull->GetMean()-1.5*h_vtx_dz3D_pull->GetRMS(), h_vtx_dz3D_pull->GetMean()+1.5*h_vtx_dz3D_pull->GetRMS() );
  h_vtx_dz3D_pull->Scale(1./h_vtx_dz3D_pull->GetSumOfWeights());
  h_vtx_dz3D_pull->Fit("f1", "QR");

  TF1 *f2 = new TF1("f2","gaus",-10., 10.);
  f2->SetParameter(2, h_vtx_dz4D_pull->GetRMS() );
  f2->SetRange(h_vtx_dz4D_pull->GetMean()-1.5*h_vtx_dz4D_pull->GetRMS(), h_vtx_dz4D_pull->GetMean()+1.5*h_vtx_dz4D_pull->GetRMS() );
  h_vtx_dz4D_pull->Scale(1./h_vtx_dz4D_pull->GetSumOfWeights());
  h_vtx_dz4D_pull->Fit("f1", "QR");

  TF1 *f3 = new TF1("f3","gaus",-10., 10.);
  f3->SetParameter(2, h_vtx_dt4D_pull->GetRMS() );
  f3->SetRange(h_vtx_dt4D_pull->GetMean()-1.5*h_vtx_dt4D_pull->GetRMS(), h_vtx_dt4D_pull->GetMean()+1.5*h_vtx_dt4D_pull->GetRMS() );
  h_vtx_dt4D_pull->Scale(1./h_vtx_dt4D_pull->GetSumOfWeights());
  h_vtx_dt4D_pull->Fit("f1", "QR");


  cout << "Saving histograms in file " << argv[5] << endl;

  // -- save output file
  TFile *fout = new TFile(argv[5],"recreate");
  h_npu->Write();

  h_vtx_dz3D-> Write();
  h_vtx_dz4D-> Write();
  h_vtx_dt4D-> Write();

  h_vtx_dz3D_pull-> Write();
  h_vtx_dz4D_pull-> Write();
  h_vtx_dt4D_pull-> Write();

  h_muon_pt->Write();
  h_muon_eta->Write();
  h_muon_phi->Write();
  h2_muon_pt_vs_eta->Write();

  h_muon_relChIso03_ratio ->Write();
  h_muon_relChIso03_ratio_barrel ->Write();
  h_muon_relChIso03_ratio_endcap ->Write();

  h_muon_relChIso03->Write();
  h_muon_relChIso03_barrel->Write();
  h_muon_relChIso03_endcap->Write();
  h_muon_relChIso03_dT->Write();
  h_muon_relChIso03_dT_barrel->Write();
  h_muon_relChIso03_dT_endcap->Write();

  h_muon_relChIso03_simVtx->Write();
  h_muon_relChIso03_simVtx_barrel->Write();
  h_muon_relChIso03_simVtx_endcap->Write();
  h_muon_relChIso03_dT_simVtx->Write();
  h_muon_relChIso03_dT_simVtx_barrel->Write();
  h_muon_relChIso03_dT_simVtx_endcap->Write();

  h_muon_relChIso03_dZmu->Write();
  h_muon_relChIso03_dZmu_barrel->Write();
  h_muon_relChIso03_dZmu_endcap->Write();

  h_muon_relChIso03_dZmu_dTmu->Write();
  h_muon_relChIso03_dZmu_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu_dTmu_endcap->Write();


  h_muon_relChIso03_dZ05->Write();
  h_muon_relChIso03_dZ05_barrel->Write();
  h_muon_relChIso03_dZ05_endcap->Write();
  h_muon_relChIso03_dZ05_dT->Write();
  h_muon_relChIso03_dZ05_dT_barrel->Write();
  h_muon_relChIso03_dZ05_dT_endcap->Write();

  h_muon_relChIso03_dZ2->Write();
  h_muon_relChIso03_dZ2_barrel->Write();
  h_muon_relChIso03_dZ2_endcap->Write();
  h_muon_relChIso03_dZ2_dT->Write();
  h_muon_relChIso03_dZ2_dT_barrel->Write();
  h_muon_relChIso03_dZ2_dT_endcap->Write();


  h_muon_relChIso03_reldZ->Write();
  h_muon_relChIso03_reldZ_barrel->Write();
  h_muon_relChIso03_reldZ_endcap->Write();
  h_muon_relChIso03_reldZ_dT->Write();
  h_muon_relChIso03_reldZ_dT_barrel->Write();
  h_muon_relChIso03_reldZ_dT_endcap->Write();

  fout->Close();
  
}


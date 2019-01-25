// per compilare: g++ -Wall -o PhotonIsolationAnalysis `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit PhotonIsolationAnalysis.cpp

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
    std::cerr << "Usage: " << argv[0] << " [folder]  [timeResol]  [usePromptPhotons] [ptReweighting] [outputFile]" << std::endl;
    return 1;
  }
  
  string folderName(argv[1]);
  string timeResolution(argv[2]);
  int usePromptPhotons = atoi(argv[3]); 
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
  chain->Add((folderName+"/0000/photonIsolation_*.root").c_str());

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
  vector<float> *photon_t;
  vector<float> *photon_pt;
  vector<float> *photon_eta;
  vector<float> *photon_phi;
  vector<float> *photon_r9;
  vector<float> *photon_sigmaIetaIeta;
  vector<float> *photon_hasPixelSeed;
  vector<int>   *photon_isPrompt;
  vector<int>   *photon_isMatchedToGenJet;

  vector<float> *photon_chIso03_dZ05_simVtx;
  vector<float> *photon_chIso03_dZ05_dT2s_simVtx;
  vector<float> *photon_chIso03_dZ05_dT3s_simVtx;
  vector<float> *photon_chIso03_dZ05_dT5s_simVtx;
  vector<float> *photon_chIso03_dZ1_simVtx;
  vector<float> *photon_chIso03_dZ1_dT2s_simVtx;
  vector<float> *photon_chIso03_dZ1_dT3s_simVtx;
  vector<float> *photon_chIso03_dZ1_dT5s_simVtx;
  vector<float> *photon_chIso03_dZ2_simVtx;
  vector<float> *photon_chIso03_dZ2_dT2s_simVtx;
  vector<float> *photon_chIso03_dZ2_dT3s_simVtx;
  vector<float> *photon_chIso03_dZ2_dT5s_simVtx;

  vector<float> *photon_chIso03_dZ05;
  vector<float> *photon_chIso03_dZ05_dT2s;
  vector<float> *photon_chIso03_dZ05_dT3s;
  vector<float> *photon_chIso03_dZ05_dT5s;
  vector<float> *photon_chIso03_dZ1;
  vector<float> *photon_chIso03_dZ1_dT2s;
  vector<float> *photon_chIso03_dZ1_dT3s;
  vector<float> *photon_chIso03_dZ1_dT5s;
  vector<float> *photon_chIso03_dZ2;
  vector<float> *photon_chIso03_dZ2_dT2s;
  vector<float> *photon_chIso03_dZ2_dT3s;
  vector<float> *photon_chIso03_dZ2_dT5s;

  photon_t = 0;
  photon_pt = 0;
  photon_eta = 0;
  photon_phi = 0;
  photon_r9 = 0;
  photon_sigmaIetaIeta = 0;
  photon_hasPixelSeed = 0;
  photon_isPrompt = 0;
  photon_isMatchedToGenJet = 0;

  photon_chIso03_dZ05_simVtx = 0;
  photon_chIso03_dZ05_dT2s_simVtx = 0;
  photon_chIso03_dZ05_dT3s_simVtx = 0;
  photon_chIso03_dZ05_dT5s_simVtx = 0;
  photon_chIso03_dZ1_simVtx = 0;
  photon_chIso03_dZ1_dT2s_simVtx = 0;
  photon_chIso03_dZ1_dT3s_simVtx = 0;
  photon_chIso03_dZ1_dT5s_simVtx = 0;
  photon_chIso03_dZ2_simVtx = 0;
  photon_chIso03_dZ2_dT2s_simVtx = 0;
  photon_chIso03_dZ2_dT3s_simVtx = 0;
  photon_chIso03_dZ2_dT5s_simVtx = 0;

  photon_chIso03_dZ05 = 0;
  photon_chIso03_dZ05_dT2s = 0;
  photon_chIso03_dZ05_dT3s = 0;
  photon_chIso03_dZ05_dT5s = 0;
  photon_chIso03_dZ1 = 0;
  photon_chIso03_dZ1_dT2s = 0;
  photon_chIso03_dZ1_dT3s = 0;
  photon_chIso03_dZ1_dT5s = 0;
  photon_chIso03_dZ2 = 0;
  photon_chIso03_dZ2_dT2s = 0;
  photon_chIso03_dZ2_dT3s = 0;
  photon_chIso03_dZ2_dT5s = 0;

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

  chain->SetBranchStatus("photon_t",1);                  chain->SetBranchAddress("photon_t",        &photon_t);
  chain->SetBranchStatus("photon_pt",1);                 chain->SetBranchAddress("photon_pt",        &photon_pt);
  chain->SetBranchStatus("photon_eta",1);                chain->SetBranchAddress("photon_eta",       &photon_eta);
  chain->SetBranchStatus("photon_phi",1);                chain->SetBranchAddress("photon_phi",       &photon_phi);
  chain->SetBranchStatus("photon_r9",1);                 chain->SetBranchAddress("photon_r9",        &photon_r9);
  chain->SetBranchStatus("photon_sigmaIetaIeta",1);      chain->SetBranchAddress("photon_sigmaIetaIeta", &photon_sigmaIetaIeta);
  chain->SetBranchStatus("photon_hasPixelSeed",1);       chain->SetBranchAddress("photon_hasPixelSeed", &photon_hasPixelSeed);
  chain->SetBranchStatus("photon_isPrompt",1);           chain->SetBranchAddress("photon_isPrompt",  &photon_isPrompt);
  chain->SetBranchStatus("photon_isMatchedToGenJet",1);  chain->SetBranchAddress("photon_isMatchedToGenJet",  &photon_isMatchedToGenJet);

  chain->SetBranchStatus("photon_chIso03_dZ05_simVtx",1);    chain->SetBranchAddress("photon_chIso03_dZ05_simVtx",   &photon_chIso03_dZ05_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT2s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ05_dT2s_simVtx",&photon_chIso03_dZ05_dT2s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT3s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ05_dT3s_simVtx",&photon_chIso03_dZ05_dT3s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT5s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ05_dT5s_simVtx",&photon_chIso03_dZ05_dT5s_simVtx);

  chain->SetBranchStatus("photon_chIso03_dZ1_simVtx",1);     chain->SetBranchAddress("photon_chIso03_dZ1_simVtx",    &photon_chIso03_dZ1_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT2s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ1_dT2s_simVtx",&photon_chIso03_dZ1_dT2s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT3s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ1_dT3s_simVtx",&photon_chIso03_dZ1_dT3s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT5s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ1_dT5s_simVtx",&photon_chIso03_dZ1_dT5s_simVtx);

  chain->SetBranchStatus("photon_chIso03_dZ2_simVtx",1);     chain->SetBranchAddress("photon_chIso03_dZ2_simVtx",    &photon_chIso03_dZ2_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT2s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ2_dT2s_simVtx",&photon_chIso03_dZ2_dT2s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT3s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ2_dT3s_simVtx",&photon_chIso03_dZ2_dT3s_simVtx);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT5s_simVtx",1); chain->SetBranchAddress("photon_chIso03_dZ2_dT5s_simVtx",&photon_chIso03_dZ2_dT5s_simVtx);

  chain->SetBranchStatus("photon_chIso03_dZ05",1);      chain->SetBranchAddress("photon_chIso03_dZ05",   &photon_chIso03_dZ05);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT2s",1);   chain->SetBranchAddress("photon_chIso03_dZ05_dT2s",&photon_chIso03_dZ05_dT2s);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT3s",1);   chain->SetBranchAddress("photon_chIso03_dZ05_dT3s",&photon_chIso03_dZ05_dT3s);
  chain->SetBranchStatus("photon_chIso03_dZ05_dT5s",1);   chain->SetBranchAddress("photon_chIso03_dZ05_dT5s",&photon_chIso03_dZ05_dT5s);

  chain->SetBranchStatus("photon_chIso03_dZ1",1);       chain->SetBranchAddress("photon_chIso03_dZ1",    &photon_chIso03_dZ1);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT2s",1);   chain->SetBranchAddress("photon_chIso03_dZ1_dT2s",&photon_chIso03_dZ1_dT2s);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT3s",1);   chain->SetBranchAddress("photon_chIso03_dZ1_dT3s",&photon_chIso03_dZ1_dT3s);
  chain->SetBranchStatus("photon_chIso03_dZ1_dT5s",1);   chain->SetBranchAddress("photon_chIso03_dZ1_dT5s",&photon_chIso03_dZ1_dT5s);

  chain->SetBranchStatus("photon_chIso03_dZ2",1);       chain->SetBranchAddress("photon_chIso03_dZ2",    &photon_chIso03_dZ2);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT2s",1);   chain->SetBranchAddress("photon_chIso03_dZ2_dT2s",&photon_chIso03_dZ2_dT2s);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT3s",1);   chain->SetBranchAddress("photon_chIso03_dZ2_dT3s",&photon_chIso03_dZ2_dT3s);
  chain->SetBranchStatus("photon_chIso03_dZ2_dT5s",1);   chain->SetBranchAddress("photon_chIso03_dZ2_dT5s",&photon_chIso03_dZ2_dT5s);


  // -- book histograms
  TH1F *h_npu   = new TH1F("h_npu","h_npu",200,50,250);
  h_npu->GetXaxis()->SetTitle("number of pileup vertices");

  TH1F *h_vtx_dz3D = new TH1F("h_vtx_dz3D","h_vtx_dz3D",1000, -0.2, 0.2);
  h_vtx_dz3D -> GetXaxis()->SetTitle("z_{3Dvtx} - z_{gen} (cm)");

  TH1F *h_vtx_dz4D = new TH1F("h_vtx_dz4D","h_vtx_dz4D",1000, -0.2, 0.2);
  h_vtx_dz4D -> GetXaxis()->SetTitle("z_{4Dvtx} - z_{gen} (cm)");

  TH1F *h_vtx_dt4D = new TH1F("h_vtx_dt4D","h_vtx_dt4D",1000, -0.5, 0.5);
  h_vtx_dt4D -> GetXaxis()->SetTitle("t_{4Dvtx} - t_{gen} (ns)");

  TH1F *h_vtx_dz3D_pull = new TH1F("h_vtx_dz3D_pull","h_vtx_dz3D_pull",200, -10, 10);
  h_vtx_dz3D_pull -> GetXaxis()->SetTitle("(z_{3Dvtx} - z_{gen})/#sigma_{z}");

  TH1F *h_vtx_dz4D_pull = new TH1F("h_vtx_dz4D_pull","h_vtx_dz4D_pull",200, -10, 10);
  h_vtx_dz4D_pull -> GetXaxis()->SetTitle("(z_{4Dvtx} - z_{gen})/#sigma_{z}");

  TH1F *h_vtx_dt4D_pull = new TH1F("h_vtx_dt4D_pull","h_vtx_dt4D_pull",200, -10, 10);
  h_vtx_dt4D_pull -> GetXaxis()->SetTitle("(t_{4Dvtx} - t_{gen})/#sigma_{t}");


  TH1F *h_photon_time = new TH1F("h_photon_time","h_photon_time",500,-1,1);
  h_photon_time->GetXaxis()->SetTitle("photon time (ns)");

  TH1F *h_photon_pt = new TH1F("h_photon_pt","h_photon_pt",200,0,200);
  h_photon_pt->GetXaxis()->SetTitle("photon p_{T} (GeV)");
  
  TH1F *h_photon_eta = new TH1F("h_photon_eta","h_photon_eta",100,-3,3);
  h_photon_eta->GetXaxis()->SetTitle("photon eta");

  TH1F *h_photon_phi = new TH1F("h_photon_phi","h_photon_phi",100,-4,4);
  h_photon_phi->GetXaxis()->SetTitle("photon phi");

  TH2F *h2_photon_pt_vs_eta = new TH2F("h2_photon_pt_vs_eta","h2_photon_pt_vs_eta",150, -3, 3, 500, 0, 1000);
  
  // --- chIso wrt to sim vertex
 
  // -- dz = 0.05 
  TH1F *h_photon_relChIso03_dZ05_simVtx = new TH1F("h_photon_relChIso03_dZ05_simVtx","h_photon_relChIso03_dZ05_simVtx",5000,0,5);
  h_photon_relChIso03_dZ05_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s_simVtx = new TH1F("h_photon_relChIso03_dZ05_dT2s_simVtx","h_photon_relChIso03_dZ05_dT2s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s_simVtx = new TH1F("h_photon_relChIso03_dZ05_dT3s_simVtx","h_photon_relChIso03_dZ05_dT3s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s_simVtx = new TH1F("h_photon_relChIso03_dZ05_dT5s_simVtx","h_photon_relChIso03_dZ05_dT5s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ05_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ05_simVtx_barrel","h_photon_relChIso03_dZ05_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ05_dT2s_simVtx_barrel","h_photon_relChIso03_dZ05_dT2s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ05_dT3s_simVtx_barrel","h_photon_relChIso03_dZ05_dT3s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ05_dT5s_simVtx_barrel","h_photon_relChIso03_dZ05_dT5s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ05_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ05_simVtx_endcap","h_photon_relChIso03_dZ05_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ05_dT2s_simVtx_endcap","h_photon_relChIso03_dZ05_dT2s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ05_dT3s_simVtx_endcap","h_photon_relChIso03_dZ05_dT3s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ05_dT5s_simVtx_endcap","h_photon_relChIso03_dZ05_dT5s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  // -- dz = 1. mm
  TH1F *h_photon_relChIso03_dZ1_simVtx = new TH1F("h_photon_relChIso03_dZ1_simVtx","h_photon_relChIso03_dZ1_simVtx",5000,0,5);
  h_photon_relChIso03_dZ1_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s_simVtx = new TH1F("h_photon_relChIso03_dZ1_dT2s_simVtx","h_photon_relChIso03_dZ1_dT2s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s_simVtx = new TH1F("h_photon_relChIso03_dZ1_dT3s_simVtx","h_photon_relChIso03_dZ1_dT3s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s_simVtx = new TH1F("h_photon_relChIso03_dZ1_dT5s_simVtx","h_photon_relChIso03_dZ1_dT5s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ1_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ1_simVtx_barrel","h_photon_relChIso03_dZ1_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ1_dT2s_simVtx_barrel","h_photon_relChIso03_dZ1_dT2s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ1_dT3s_simVtx_barrel","h_photon_relChIso03_dZ1_dT3s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ1_dT5s_simVtx_barrel","h_photon_relChIso03_dZ1_dT5s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ1_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ1_simVtx_endcap","h_photon_relChIso03_dZ1_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ1_dT2s_simVtx_endcap","h_photon_relChIso03_dZ1_dT2s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ1_dT3s_simVtx_endcap","h_photon_relChIso03_dZ1_dT3s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ1_dT5s_simVtx_endcap","h_photon_relChIso03_dZ1_dT5s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 2. mm
  TH1F *h_photon_relChIso03_dZ2_simVtx = new TH1F("h_photon_relChIso03_dZ2_simVtx","h_photon_relChIso03_dZ2_simVtx",5000,0,5);
  h_photon_relChIso03_dZ2_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s_simVtx = new TH1F("h_photon_relChIso03_dZ2_dT2s_simVtx","h_photon_relChIso03_dZ2_dT2s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s_simVtx = new TH1F("h_photon_relChIso03_dZ2_dT3s_simVtx","h_photon_relChIso03_dZ2_dT3s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s_simVtx = new TH1F("h_photon_relChIso03_dZ2_dT5s_simVtx","h_photon_relChIso03_dZ2_dT5s_simVtx",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ2_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ2_simVtx_barrel","h_photon_relChIso03_dZ2_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ2_dT2s_simVtx_barrel","h_photon_relChIso03_dZ2_dT2s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ2_dT3s_simVtx_barrel","h_photon_relChIso03_dZ2_dT3s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s_simVtx_barrel = new TH1F("h_photon_relChIso03_dZ2_dT5s_simVtx_barrel","h_photon_relChIso03_dZ2_dT5s_simVtx_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ2_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ2_simVtx_endcap","h_photon_relChIso03_dZ2_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ2_dT2s_simVtx_endcap","h_photon_relChIso03_dZ2_dT2s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ2_dT3s_simVtx_endcap","h_photon_relChIso03_dZ2_dT3s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s_simVtx_endcap = new TH1F("h_photon_relChIso03_dZ2_dT5s_simVtx_endcap","h_photon_relChIso03_dZ2_dT5s_simVtx_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");
  

  // --- reco vtx closest to sim vertex
 
  // -- dz = 0.05 
  TH1F *h_photon_relChIso03_dZ05 = new TH1F("h_photon_relChIso03_dZ05","h_photon_relChIso03_dZ05",5000,0,5);
  h_photon_relChIso03_dZ05->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s = new TH1F("h_photon_relChIso03_dZ05_dT2s","h_photon_relChIso03_dZ05_dT2s",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s = new TH1F("h_photon_relChIso03_dZ05_dT3s","h_photon_relChIso03_dZ05_dT3s",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s = new TH1F("h_photon_relChIso03_dZ05_dT5s","h_photon_relChIso03_dZ05_dT5s",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ05_barrel = new TH1F("h_photon_relChIso03_dZ05_barrel","h_photon_relChIso03_dZ05_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s_barrel = new TH1F("h_photon_relChIso03_dZ05_dT2s_barrel","h_photon_relChIso03_dZ05_dT2s_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s_barrel = new TH1F("h_photon_relChIso03_dZ05_dT3s_barrel","h_photon_relChIso03_dZ05_dT3s_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s_barrel = new TH1F("h_photon_relChIso03_dZ05_dT5s_barrel","h_photon_relChIso03_dZ05_dT5s_barrel",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ05_endcap = new TH1F("h_photon_relChIso03_dZ05_endcap","h_photon_relChIso03_dZ05_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT2s_endcap = new TH1F("h_photon_relChIso03_dZ05_dT2s_endcap","h_photon_relChIso03_dZ05_dT2s_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT3s_endcap = new TH1F("h_photon_relChIso03_dZ05_dT3s_endcap","h_photon_relChIso03_dZ05_dT3s_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ05_dT5s_endcap = new TH1F("h_photon_relChIso03_dZ05_dT5s_endcap","h_photon_relChIso03_dZ05_dT5s_endcap",5000,0,5);
  h_photon_relChIso03_dZ05_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  // -- dz = 1. mm
  TH1F *h_photon_relChIso03_dZ1 = new TH1F("h_photon_relChIso03_dZ1","h_photon_relChIso03_dZ1",5000,0,5);
  h_photon_relChIso03_dZ1->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s = new TH1F("h_photon_relChIso03_dZ1_dT2s","h_photon_relChIso03_dZ1_dT2s",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s = new TH1F("h_photon_relChIso03_dZ1_dT3s","h_photon_relChIso03_dZ1_dT3s",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s = new TH1F("h_photon_relChIso03_dZ1_dT5s","h_photon_relChIso03_dZ1_dT5s",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ1_barrel = new TH1F("h_photon_relChIso03_dZ1_barrel","h_photon_relChIso03_dZ1_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s_barrel = new TH1F("h_photon_relChIso03_dZ1_dT2s_barrel","h_photon_relChIso03_dZ1_dT2s_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s_barrel = new TH1F("h_photon_relChIso03_dZ1_dT3s_barrel","h_photon_relChIso03_dZ1_dT3s_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s_barrel = new TH1F("h_photon_relChIso03_dZ1_dT5s_barrel","h_photon_relChIso03_dZ1_dT5s_barrel",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ1_endcap = new TH1F("h_photon_relChIso03_dZ1_endcap","h_photon_relChIso03_dZ1_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT2s_endcap = new TH1F("h_photon_relChIso03_dZ1_dT2s_endcap","h_photon_relChIso03_dZ1_dT2s_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT3s_endcap = new TH1F("h_photon_relChIso03_dZ1_dT3s_endcap","h_photon_relChIso03_dZ1_dT3s_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ1_dT5s_endcap = new TH1F("h_photon_relChIso03_dZ1_dT5s_endcap","h_photon_relChIso03_dZ1_dT5s_endcap",5000,0,5);
  h_photon_relChIso03_dZ1_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 2. mm
  TH1F *h_photon_relChIso03_dZ2 = new TH1F("h_photon_relChIso03_dZ2","h_photon_relChIso03_dZ2",5000,0,5);
  h_photon_relChIso03_dZ2->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s = new TH1F("h_photon_relChIso03_dZ2_dT2s","h_photon_relChIso03_dZ2_dT2s",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s = new TH1F("h_photon_relChIso03_dZ2_dT3s","h_photon_relChIso03_dZ2_dT3s",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s = new TH1F("h_photon_relChIso03_dZ2_dT5s","h_photon_relChIso03_dZ2_dT5s",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ2_barrel = new TH1F("h_photon_relChIso03_dZ2_barrel","h_photon_relChIso03_dZ2_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s_barrel = new TH1F("h_photon_relChIso03_dZ2_dT2s_barrel","h_photon_relChIso03_dZ2_dT2s_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s_barrel = new TH1F("h_photon_relChIso03_dZ2_dT3s_barrel","h_photon_relChIso03_dZ2_dT3s_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s_barrel = new TH1F("h_photon_relChIso03_dZ2_dT5s_barrel","h_photon_relChIso03_dZ2_dT5s_barrel",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_photon_relChIso03_dZ2_endcap = new TH1F("h_photon_relChIso03_dZ2_endcap","h_photon_relChIso03_dZ2_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT2s_endcap = new TH1F("h_photon_relChIso03_dZ2_dT2s_endcap","h_photon_relChIso03_dZ2_dT2s_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT3s_endcap = new TH1F("h_photon_relChIso03_dZ2_dT3s_endcap","h_photon_relChIso03_dZ2_dT3s_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_photon_relChIso03_dZ2_dT5s_endcap = new TH1F("h_photon_relChIso03_dZ2_dT5s_endcap","h_photon_relChIso03_dZ2_dT5s_endcap",5000,0,5);
  h_photon_relChIso03_dZ2_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // ratio iso
  TH1F *h_photon_relChIso03_ratio = new TH1F("h_photon_relChIso03_ratio","h_photon_relChIso03_ratio",1000,0,2);
  h_photon_relChIso03_ratio->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_photon_relChIso03_ratio_barrel = new TH1F("h_photon_relChIso03_ratio_barrel","h_photon_relChIso03_ratio_barrel",1000,0,2);
  h_photon_relChIso03_ratio_barrel->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_photon_relChIso03_ratio_endcap = new TH1F("h_photon_relChIso03_ratio_endcap","h_photon_relChIso03_ratio_endcap",1000,0,2);
  h_photon_relChIso03_ratio_endcap->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");


  // line density
  TH1F *h_linedensity_noRelChIsoZCut = new TH1F("h_linedensity_noRelChIsoZCut","h_linedensity_noRelChIsoZCut", 100, 0, 2);
  TH1F *h_linedensity_noRelChIsoZTCut = new TH1F("h_linedensity_noRelChIsoZTCut","h_linedensity_noRelChIsoZTCut", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut = new TH1F("h_linedensity_RelChIsoZCut ","h_linedensity_RelChIsoZCut ", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut  = new TH1F("h_linedensity_RelChIsoZTCut","h_linedensity_RelChIsoZTCut", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_simVtx = new TH1F("h_linedensity_RelChIsoZCut_simVtx","h_linedensity_RelChIsoZCut_simVtx", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_simVtx  = new TH1F("h_linedensity_RelChIsoZTCut_simVtx","h_linedensity_RelChIsoZTCut_simVtx", 100, 0, 2);
  


  // -- loop over events
  //int maxEntries = std::min(int(chain ->GetEntries()),1000000);  
  int maxEntries = chain ->GetEntries();
  cout << "Number of events to be analyzed : " << maxEntries <<endl;
  float w = 1;
  float pt = 0;

  int nPhotonsInEvent = 0;

  for (int i = 0; i < maxEntries; i++) {
    
    chain -> GetEntry(i);
    
    if (i%1000==0) cout << "Analyzing event " << i << "\r" << flush;
    
    h_npu->Fill(npu);

    nPhotonsInEvent = 0;

    for (unsigned int ipho = 0; ipho < photon_pt->size(); ipho++){
      

      pt = photon_pt->at(ipho);
      
      if ( pt < 15. )  continue;
      // - pixel seed veto
      //if ( photon_hasPixelSeed->at(ipho) ) continue;


      // -- pt/eta reweighting
      if (applyPtReweighting){
	int bin = hratio2D->FindBin(photon_eta->at(ipho), pt ); 
	w = hratio2D -> GetBinContent(bin); 
      }
      else {
	w = 1;
      }
      
      // -- prompt photons or non-prompt photons
      bool pass = false;
      if (usePromptPhotons && photon_isPrompt->at(ipho)) pass = true; 
      if (!usePromptPhotons && (!photon_isPrompt->at(ipho) && photon_isMatchedToGenJet->at(ipho) ) ) pass = true;
      
      if (!pass) continue;

      nPhotonsInEvent++;
      
      bool pass3D = !vtx3D_isFake;
      bool pass4D = !vtx4D_isFake;
      
      
      if (nPhotonsInEvent == 1){

      	if (pass3D) h_vtx_dz3D->Fill( vtx3D_z - vtxGen_z);
      	if (pass4D) h_vtx_dz4D->Fill( vtx4D_z - vtxGen_z);
      	if (pass4D) h_vtx_dt4D->Fill( vtx4D_t - vtxGen_t*1000000000.);

      	if (pass3D) h_vtx_dz3D_pull->Fill( (vtx3D_z - vtxGen_z) / vtx3D_zErr );
	if (pass4D) h_vtx_dz4D_pull->Fill( (vtx4D_z - vtxGen_z) / vtx4D_zErr );
      	if (pass4D) h_vtx_dt4D_pull->Fill( (vtx4D_t - vtxGen_t*1000000000.) / vtx4D_tErr);

      }

      h_photon_time->Fill(photon_t->at(ipho), w);
      h_photon_pt->Fill(pt, w);
      h_photon_eta->Fill(photon_eta->at(ipho), w);
      h_photon_phi->Fill(photon_phi->at(ipho), w);
      h2_photon_pt_vs_eta->Fill(photon_eta->at(ipho), pt, w);


      // control plots
      if (pass3D && pass4D){	
	float chIsoRatio = photon_chIso03_dZ1_dT3s->at(ipho)/photon_chIso03_dZ1->at(ipho);
	if ( photon_chIso03_dZ1->at(ipho) == 0 ) chIsoRatio = 1;
	h_photon_relChIso03_ratio -> Fill( chIsoRatio, w );
	if (fabs(photon_eta->at(ipho))<1.48){ 
	  h_photon_relChIso03_ratio_barrel -> Fill( chIsoRatio, w );
	}
        else{ 
	  h_photon_relChIso03_ratio_endcap -> Fill( chIsoRatio, w );
	}
      }
      
      // chiIso plots - all
      if ( pass3D ){
	h_photon_relChIso03_dZ05_simVtx -> Fill(photon_chIso03_dZ05_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1_simVtx -> Fill(photon_chIso03_dZ1_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2_simVtx -> Fill(photon_chIso03_dZ2_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ05 -> Fill(photon_chIso03_dZ05->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1  -> Fill(photon_chIso03_dZ1->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2  -> Fill(photon_chIso03_dZ2->at(ipho)/pt, w);
      }
      if ( pass4D ){
	h_photon_relChIso03_dZ05_dT2s_simVtx -> Fill(photon_chIso03_dZ05_dT2s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ05_dT3s_simVtx -> Fill(photon_chIso03_dZ05_dT3s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ05_dT5s_simVtx -> Fill(photon_chIso03_dZ05_dT5s_simVtx->at(ipho)/pt, w);

	h_photon_relChIso03_dZ1_dT2s_simVtx -> Fill(photon_chIso03_dZ1_dT2s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1_dT3s_simVtx -> Fill(photon_chIso03_dZ1_dT3s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1_dT5s_simVtx -> Fill(photon_chIso03_dZ1_dT5s_simVtx->at(ipho)/pt, w);

	h_photon_relChIso03_dZ2_dT2s_simVtx -> Fill(photon_chIso03_dZ2_dT2s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2_dT3s_simVtx -> Fill(photon_chIso03_dZ2_dT3s_simVtx->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2_dT5s_simVtx -> Fill(photon_chIso03_dZ2_dT5s_simVtx->at(ipho)/pt, w);

	h_photon_relChIso03_dZ05_dT2s -> Fill(photon_chIso03_dZ05_dT2s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ05_dT3s -> Fill(photon_chIso03_dZ05_dT3s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ05_dT5s -> Fill(photon_chIso03_dZ05_dT5s->at(ipho)/pt, w);

	h_photon_relChIso03_dZ1_dT2s -> Fill(photon_chIso03_dZ1_dT2s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1_dT3s -> Fill(photon_chIso03_dZ1_dT3s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ1_dT5s -> Fill(photon_chIso03_dZ1_dT5s->at(ipho)/pt, w);

	h_photon_relChIso03_dZ2_dT2s -> Fill(photon_chIso03_dZ2_dT2s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2_dT3s -> Fill(photon_chIso03_dZ2_dT3s->at(ipho)/pt, w);
	h_photon_relChIso03_dZ2_dT5s -> Fill(photon_chIso03_dZ2_dT5s->at(ipho)/pt, w);
      }

      // barrel 
      if (fabs(photon_eta->at(ipho))<1.48){
	if ( pass3D ) {
	  h_photon_relChIso03_dZ05_simVtx_barrel -> Fill(photon_chIso03_dZ05_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_simVtx_barrel -> Fill(photon_chIso03_dZ1_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_simVtx_barrel -> Fill(photon_chIso03_dZ2_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_barrel -> Fill(photon_chIso03_dZ05->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_barrel  -> Fill(photon_chIso03_dZ1->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_barrel  -> Fill(photon_chIso03_dZ2->at(ipho)/pt, w);
	}
	if ( pass4D ){
	  h_photon_relChIso03_dZ05_dT2s_simVtx_barrel -> Fill(photon_chIso03_dZ05_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT3s_simVtx_barrel -> Fill(photon_chIso03_dZ05_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT5s_simVtx_barrel -> Fill(photon_chIso03_dZ05_dT5s_simVtx->at(ipho)/pt, w);
	  
	  h_photon_relChIso03_dZ1_dT2s_simVtx_barrel -> Fill(photon_chIso03_dZ1_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT3s_simVtx_barrel -> Fill(photon_chIso03_dZ1_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT5s_simVtx_barrel -> Fill(photon_chIso03_dZ1_dT5s_simVtx->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ2_dT2s_simVtx_barrel -> Fill(photon_chIso03_dZ2_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT3s_simVtx_barrel -> Fill(photon_chIso03_dZ2_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT5s_simVtx_barrel -> Fill(photon_chIso03_dZ2_dT5s_simVtx->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ05_dT2s_barrel -> Fill(photon_chIso03_dZ05_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT3s_barrel -> Fill(photon_chIso03_dZ05_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT5s_barrel -> Fill(photon_chIso03_dZ05_dT5s->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ1_dT2s_barrel -> Fill(photon_chIso03_dZ1_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT3s_barrel -> Fill(photon_chIso03_dZ1_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT5s_barrel -> Fill(photon_chIso03_dZ1_dT5s->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ2_dT2s_barrel -> Fill(photon_chIso03_dZ2_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT3s_barrel -> Fill(photon_chIso03_dZ2_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT5s_barrel -> Fill(photon_chIso03_dZ2_dT5s->at(ipho)/pt, w);

	}
      }
      // endcap
      else{
	if ( pass3D ){
	  h_photon_relChIso03_dZ05_simVtx_endcap -> Fill(photon_chIso03_dZ05_simVtx->at(ipho)/pt, w);
          h_photon_relChIso03_dZ1_simVtx_endcap -> Fill(photon_chIso03_dZ1_simVtx->at(ipho)/pt, w);
          h_photon_relChIso03_dZ2_simVtx_endcap -> Fill(photon_chIso03_dZ2_simVtx->at(ipho)/pt, w);
          h_photon_relChIso03_dZ05_endcap -> Fill(photon_chIso03_dZ05->at(ipho)/pt, w);
          h_photon_relChIso03_dZ1_endcap  -> Fill(photon_chIso03_dZ1->at(ipho)/pt, w);
          h_photon_relChIso03_dZ2_endcap  -> Fill(photon_chIso03_dZ2->at(ipho)/pt, w);
	}
	if ( pass4D ) {
	  h_photon_relChIso03_dZ05_dT2s_simVtx_endcap -> Fill(photon_chIso03_dZ05_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT3s_simVtx_endcap -> Fill(photon_chIso03_dZ05_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT5s_simVtx_endcap -> Fill(photon_chIso03_dZ05_dT5s_simVtx->at(ipho)/pt, w);
	  
	  h_photon_relChIso03_dZ1_dT2s_simVtx_endcap -> Fill(photon_chIso03_dZ1_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT3s_simVtx_endcap -> Fill(photon_chIso03_dZ1_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT5s_simVtx_endcap -> Fill(photon_chIso03_dZ1_dT5s_simVtx->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ2_dT2s_simVtx_endcap -> Fill(photon_chIso03_dZ2_dT2s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT3s_simVtx_endcap -> Fill(photon_chIso03_dZ2_dT3s_simVtx->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT5s_simVtx_endcap -> Fill(photon_chIso03_dZ2_dT5s_simVtx->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ05_dT2s_endcap -> Fill(photon_chIso03_dZ05_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT3s_endcap -> Fill(photon_chIso03_dZ05_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ05_dT5s_endcap -> Fill(photon_chIso03_dZ05_dT5s->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ1_dT2s_endcap -> Fill(photon_chIso03_dZ1_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT3s_endcap -> Fill(photon_chIso03_dZ1_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ1_dT5s_endcap -> Fill(photon_chIso03_dZ1_dT5s->at(ipho)/pt, w);

	  h_photon_relChIso03_dZ2_dT2s_endcap -> Fill(photon_chIso03_dZ2_dT2s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT3s_endcap -> Fill(photon_chIso03_dZ2_dT3s->at(ipho)/pt, w);
	  h_photon_relChIso03_dZ2_dT5s_endcap -> Fill(photon_chIso03_dZ2_dT5s->at(ipho)/pt, w);
	}
      }
      
      // -- line density
      float linedensity = 200.*TMath::Gaus(fabs(10*vtxGen_z), 0, 42., 1);
      if (pass3D){
	h_linedensity_noRelChIsoZCut->Fill(linedensity,w); 
	if ( (photon_chIso03_dZ1->at(ipho)/photon_pt->at(ipho)) < 0.05 ) {
	  h_linedensity_RelChIsoZCut->Fill(linedensity,w); 
	}
	if ( (photon_chIso03_dZ1_simVtx->at(ipho)/photon_pt->at(ipho)) < 0.05 ) {
	  h_linedensity_RelChIsoZCut_simVtx->Fill(linedensity,w); 
	}
      }
      if (pass4D){
	h_linedensity_noRelChIsoZTCut->Fill(linedensity,w);
	if ( (photon_chIso03_dZ1_dT3s->at(ipho)/photon_pt->at(ipho)) < 0.05 ) {
	  h_linedensity_RelChIsoZTCut->Fill(linedensity,w); 
	}
	if ( (photon_chIso03_dZ1_dT3s_simVtx->at(ipho)/photon_pt->at(ipho)) < 0.05 ) {
	  h_linedensity_RelChIsoZTCut_simVtx->Fill(linedensity,w); 
	}

      }
    

    }// end loop over photons
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

  h_photon_time->Write();
  h_photon_pt->Write();
  h_photon_eta->Write();
  h_photon_phi->Write();
  h2_photon_pt_vs_eta->Write();

  h_photon_relChIso03_ratio ->Write();
  h_photon_relChIso03_ratio_barrel ->Write();
  h_photon_relChIso03_ratio_endcap ->Write();


  h_photon_relChIso03_dZ05_simVtx->Write();
  h_photon_relChIso03_dZ05_simVtx_barrel->Write();
  h_photon_relChIso03_dZ05_simVtx_endcap->Write();

  h_photon_relChIso03_dZ05_dT2s_simVtx->Write();
  h_photon_relChIso03_dZ05_dT2s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ05_dT2s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ05_dT3s_simVtx->Write();
  h_photon_relChIso03_dZ05_dT3s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ05_dT3s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ05_dT5s_simVtx->Write();
  h_photon_relChIso03_dZ05_dT5s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ05_dT5s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ1_simVtx->Write();
  h_photon_relChIso03_dZ1_simVtx_barrel->Write();
  h_photon_relChIso03_dZ1_simVtx_endcap->Write();

  h_photon_relChIso03_dZ1_dT2s_simVtx->Write();
  h_photon_relChIso03_dZ1_dT2s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ1_dT2s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ1_dT3s_simVtx->Write();
  h_photon_relChIso03_dZ1_dT3s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ1_dT3s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ1_dT5s_simVtx->Write();
  h_photon_relChIso03_dZ1_dT5s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ1_dT5s_simVtx_endcap->Write();
  
  h_photon_relChIso03_dZ2_simVtx->Write();
  h_photon_relChIso03_dZ2_simVtx_barrel->Write();
  h_photon_relChIso03_dZ2_simVtx_endcap->Write();

  h_photon_relChIso03_dZ2_dT2s_simVtx->Write();
  h_photon_relChIso03_dZ2_dT2s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ2_dT2s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ2_dT3s_simVtx->Write();
  h_photon_relChIso03_dZ2_dT3s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ2_dT3s_simVtx_endcap->Write();

  h_photon_relChIso03_dZ2_dT5s_simVtx->Write();
  h_photon_relChIso03_dZ2_dT5s_simVtx_barrel->Write();
  h_photon_relChIso03_dZ2_dT5s_simVtx_endcap->Write();



  h_photon_relChIso03_dZ05->Write();
  h_photon_relChIso03_dZ05_barrel->Write();
  h_photon_relChIso03_dZ05_endcap->Write();

  h_photon_relChIso03_dZ05_dT2s->Write();
  h_photon_relChIso03_dZ05_dT2s_barrel->Write();
  h_photon_relChIso03_dZ05_dT2s_endcap->Write();

  h_photon_relChIso03_dZ05_dT3s->Write();
  h_photon_relChIso03_dZ05_dT3s_barrel->Write();
  h_photon_relChIso03_dZ05_dT3s_endcap->Write();

  h_photon_relChIso03_dZ05_dT5s->Write();
  h_photon_relChIso03_dZ05_dT5s_barrel->Write();
  h_photon_relChIso03_dZ05_dT5s_endcap->Write();

  h_photon_relChIso03_dZ1->Write();
  h_photon_relChIso03_dZ1_barrel->Write();
  h_photon_relChIso03_dZ1_endcap->Write();

  h_photon_relChIso03_dZ1_dT2s->Write();
  h_photon_relChIso03_dZ1_dT2s_barrel->Write();
  h_photon_relChIso03_dZ1_dT2s_endcap->Write();

  h_photon_relChIso03_dZ1_dT3s->Write();
  h_photon_relChIso03_dZ1_dT3s_barrel->Write();
  h_photon_relChIso03_dZ1_dT3s_endcap->Write();

  h_photon_relChIso03_dZ1_dT5s->Write();
  h_photon_relChIso03_dZ1_dT5s_barrel->Write();
  h_photon_relChIso03_dZ1_dT5s_endcap->Write();
  
  h_photon_relChIso03_dZ2->Write();
  h_photon_relChIso03_dZ2_barrel->Write();
  h_photon_relChIso03_dZ2_endcap->Write();

  h_photon_relChIso03_dZ2_dT2s->Write();
  h_photon_relChIso03_dZ2_dT2s_barrel->Write();
  h_photon_relChIso03_dZ2_dT2s_endcap->Write();

  h_photon_relChIso03_dZ2_dT3s->Write();
  h_photon_relChIso03_dZ2_dT3s_barrel->Write();
  h_photon_relChIso03_dZ2_dT3s_endcap->Write();

  h_photon_relChIso03_dZ2_dT5s->Write();
  h_photon_relChIso03_dZ2_dT5s_barrel->Write();
  h_photon_relChIso03_dZ2_dT5s_endcap->Write();


  h_linedensity_noRelChIsoZCut->Write();
  h_linedensity_noRelChIsoZTCut->Write();

  h_linedensity_RelChIsoZCut->Write();
  h_linedensity_RelChIsoZTCut->Write();

  h_linedensity_RelChIsoZCut_simVtx->Write();
  h_linedensity_RelChIsoZTCut_simVtx->Write();

  fout->Close();
  
}


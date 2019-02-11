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
  chain->Add((folderName+"/0000/muonIsolation*.root").c_str());

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
  vector<float> *muon_t;
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

  vector<float> *muon_chIso03_dZ05_simVtx;
  vector<float> *muon_chIso03_dZ05_dT2s_simVtx;
  vector<float> *muon_chIso03_dZ05_dT3s_simVtx;
  vector<float> *muon_chIso03_dZ05_dT5s_simVtx;
  vector<float> *muon_chIso03_dZ1_simVtx;
  vector<float> *muon_chIso03_dZ1_dT2s_simVtx;
  vector<float> *muon_chIso03_dZ1_dT3s_simVtx;
  vector<float> *muon_chIso03_dZ1_dT5s_simVtx;
  vector<float> *muon_chIso03_dZ2_simVtx;
  vector<float> *muon_chIso03_dZ2_dT2s_simVtx;
  vector<float> *muon_chIso03_dZ2_dT3s_simVtx;
  vector<float> *muon_chIso03_dZ2_dT5s_simVtx;
  vector<float> *muon_chIso03_dZ3_simVtx;
  vector<float> *muon_chIso03_dZ3_dT2s_simVtx;
  vector<float> *muon_chIso03_dZ3_dT3s_simVtx;
  vector<float> *muon_chIso03_dZ3_dT5s_simVtx;
  vector<float> *muon_chIso03_dZ10_simVtx;
  vector<float> *muon_chIso03_dZ10_dT2s_simVtx;
  vector<float> *muon_chIso03_dZ10_dT3s_simVtx;
  vector<float> *muon_chIso03_dZ10_dT5s_simVtx;

  vector<float> *muon_chIso03_dZ05;
  vector<float> *muon_chIso03_dZ05_dT2s;
  vector<float> *muon_chIso03_dZ05_dT3s;
  vector<float> *muon_chIso03_dZ05_dT5s;
  vector<float> *muon_chIso03_dZ1;
  vector<float> *muon_chIso03_dZ1_dT2s;
  vector<float> *muon_chIso03_dZ1_dT3s;
  vector<float> *muon_chIso03_dZ1_dT5s;
  vector<float> *muon_chIso03_dZ2;
  vector<float> *muon_chIso03_dZ2_dT2s;
  vector<float> *muon_chIso03_dZ2_dT3s;
  vector<float> *muon_chIso03_dZ2_dT5s;
  vector<float> *muon_chIso03_dZ3;
  vector<float> *muon_chIso03_dZ3_dT2s;
  vector<float> *muon_chIso03_dZ3_dT3s;
  vector<float> *muon_chIso03_dZ3_dT5s;
  vector<float> *muon_chIso03_dZ10;
  vector<float> *muon_chIso03_dZ10_dT2s;
  vector<float> *muon_chIso03_dZ10_dT3s;
  vector<float> *muon_chIso03_dZ10_dT5s;

  vector<float> *muon_chIso03_reldZ;
  vector<float> *muon_chIso03_reldZ_dT;  

  vector<float> *muon_chIso03_dZmu05;
  vector<float> *muon_chIso03_dZmu05_dTmu;

  vector<float> *muon_chIso03_dZmu1;
  vector<float> *muon_chIso03_dZmu1_dTmu;

  vector<float> *muon_chIso03_dZmu2;
  vector<float> *muon_chIso03_dZmu2_dTmu;

  vector<float> *muon_chIso03_dZmu5;
  vector<float> *muon_chIso03_dZmu5_dTmu;

  vector<float> *muon_chIso03_dZmu10;
  vector<float> *muon_chIso03_dZmu10_dTmu;

  muon_t = 0;
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

  muon_chIso03_dZ05_simVtx = 0;
  muon_chIso03_dZ05_dT2s_simVtx = 0;
  muon_chIso03_dZ05_dT3s_simVtx = 0;
  muon_chIso03_dZ05_dT5s_simVtx = 0;
  muon_chIso03_dZ1_simVtx = 0;
  muon_chIso03_dZ1_dT2s_simVtx = 0;
  muon_chIso03_dZ1_dT3s_simVtx = 0;
  muon_chIso03_dZ1_dT5s_simVtx = 0;
  muon_chIso03_dZ2_simVtx = 0;
  muon_chIso03_dZ2_dT2s_simVtx = 0;
  muon_chIso03_dZ2_dT3s_simVtx = 0;
  muon_chIso03_dZ2_dT5s_simVtx = 0;
  muon_chIso03_dZ3_simVtx = 0;
  muon_chIso03_dZ3_dT2s_simVtx = 0;
  muon_chIso03_dZ3_dT3s_simVtx = 0;
  muon_chIso03_dZ3_dT5s_simVtx = 0;
  muon_chIso03_dZ10_simVtx = 0;
  muon_chIso03_dZ10_dT2s_simVtx = 0;
  muon_chIso03_dZ10_dT3s_simVtx = 0;
  muon_chIso03_dZ10_dT5s_simVtx = 0;

  muon_chIso03_dZ05 = 0;
  muon_chIso03_dZ05_dT2s = 0;
  muon_chIso03_dZ05_dT3s = 0;
  muon_chIso03_dZ05_dT5s = 0;
  muon_chIso03_dZ1 = 0;
  muon_chIso03_dZ1_dT2s = 0;
  muon_chIso03_dZ1_dT3s = 0;
  muon_chIso03_dZ1_dT5s = 0;
  muon_chIso03_dZ2 = 0;
  muon_chIso03_dZ2_dT2s = 0;
  muon_chIso03_dZ2_dT3s = 0;
  muon_chIso03_dZ2_dT5s = 0;
  muon_chIso03_dZ3 = 0;
  muon_chIso03_dZ3_dT2s = 0;
  muon_chIso03_dZ3_dT3s = 0;
  muon_chIso03_dZ3_dT5s = 0;
  muon_chIso03_dZ10 = 0;
  muon_chIso03_dZ10_dT2s = 0;
  muon_chIso03_dZ10_dT3s = 0;
  muon_chIso03_dZ10_dT5s = 0;

  muon_chIso03_dZmu05 = 0;
  muon_chIso03_dZmu05_dTmu = 0;

  muon_chIso03_dZmu1 = 0;
  muon_chIso03_dZmu1_dTmu = 0;

  muon_chIso03_dZmu2 = 0;
  muon_chIso03_dZmu2_dTmu = 0;

  muon_chIso03_dZmu5 = 0;
  muon_chIso03_dZmu5_dTmu = 0;

  muon_chIso03_dZmu10 = 0;
  muon_chIso03_dZmu10_dTmu = 0;

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

  chain->SetBranchStatus("muon_t",1);                  chain->SetBranchAddress("muon_t",        &muon_t);
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

  chain->SetBranchStatus("muon_chIso03_dZ05_simVtx",1);    chain->SetBranchAddress("muon_chIso03_dZ05_simVtx",   &muon_chIso03_dZ05_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT2s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ05_dT2s_simVtx",&muon_chIso03_dZ05_dT2s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT3s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ05_dT3s_simVtx",&muon_chIso03_dZ05_dT3s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT5s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ05_dT5s_simVtx",&muon_chIso03_dZ05_dT5s_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ1_simVtx",1);     chain->SetBranchAddress("muon_chIso03_dZ1_simVtx",    &muon_chIso03_dZ1_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT2s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ1_dT2s_simVtx",&muon_chIso03_dZ1_dT2s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT3s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ1_dT3s_simVtx",&muon_chIso03_dZ1_dT3s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT5s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ1_dT5s_simVtx",&muon_chIso03_dZ1_dT5s_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ2_simVtx",1);     chain->SetBranchAddress("muon_chIso03_dZ2_simVtx",    &muon_chIso03_dZ2_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT2s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ2_dT2s_simVtx",&muon_chIso03_dZ2_dT2s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT3s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ2_dT3s_simVtx",&muon_chIso03_dZ2_dT3s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT5s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ2_dT5s_simVtx",&muon_chIso03_dZ2_dT5s_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ3_simVtx",1);     chain->SetBranchAddress("muon_chIso03_dZ3_simVtx",    &muon_chIso03_dZ3_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT2s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ3_dT2s_simVtx",&muon_chIso03_dZ3_dT2s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT3s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ3_dT3s_simVtx",&muon_chIso03_dZ3_dT3s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT5s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ3_dT5s_simVtx",&muon_chIso03_dZ3_dT5s_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ10_simVtx",1);     chain->SetBranchAddress("muon_chIso03_dZ10_simVtx",    &muon_chIso03_dZ10_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT2s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ10_dT2s_simVtx",&muon_chIso03_dZ10_dT2s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT3s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ10_dT3s_simVtx",&muon_chIso03_dZ10_dT3s_simVtx);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT5s_simVtx",1); chain->SetBranchAddress("muon_chIso03_dZ10_dT5s_simVtx",&muon_chIso03_dZ10_dT5s_simVtx);

  chain->SetBranchStatus("muon_chIso03_dZ05",1);      chain->SetBranchAddress("muon_chIso03_dZ05",   &muon_chIso03_dZ05);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT2s",1);   chain->SetBranchAddress("muon_chIso03_dZ05_dT2s",&muon_chIso03_dZ05_dT2s);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT3s",1);   chain->SetBranchAddress("muon_chIso03_dZ05_dT3s",&muon_chIso03_dZ05_dT3s);
  chain->SetBranchStatus("muon_chIso03_dZ05_dT5s",1);   chain->SetBranchAddress("muon_chIso03_dZ05_dT5s",&muon_chIso03_dZ05_dT5s);

  chain->SetBranchStatus("muon_chIso03_dZ1",1);       chain->SetBranchAddress("muon_chIso03_dZ1",    &muon_chIso03_dZ1);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT2s",1);   chain->SetBranchAddress("muon_chIso03_dZ1_dT2s",&muon_chIso03_dZ1_dT2s);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT3s",1);   chain->SetBranchAddress("muon_chIso03_dZ1_dT3s",&muon_chIso03_dZ1_dT3s);
  chain->SetBranchStatus("muon_chIso03_dZ1_dT5s",1);   chain->SetBranchAddress("muon_chIso03_dZ1_dT5s",&muon_chIso03_dZ1_dT5s);

  chain->SetBranchStatus("muon_chIso03_dZ2",1);       chain->SetBranchAddress("muon_chIso03_dZ2",    &muon_chIso03_dZ2);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT2s",1);   chain->SetBranchAddress("muon_chIso03_dZ2_dT2s",&muon_chIso03_dZ2_dT2s);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT3s",1);   chain->SetBranchAddress("muon_chIso03_dZ2_dT3s",&muon_chIso03_dZ2_dT3s);
  chain->SetBranchStatus("muon_chIso03_dZ2_dT5s",1);   chain->SetBranchAddress("muon_chIso03_dZ2_dT5s",&muon_chIso03_dZ2_dT5s);

  chain->SetBranchStatus("muon_chIso03_dZ3",1);       chain->SetBranchAddress("muon_chIso03_dZ3",    &muon_chIso03_dZ3);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT2s",1);   chain->SetBranchAddress("muon_chIso03_dZ3_dT2s",&muon_chIso03_dZ3_dT2s);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT3s",1);   chain->SetBranchAddress("muon_chIso03_dZ3_dT3s",&muon_chIso03_dZ3_dT3s);
  chain->SetBranchStatus("muon_chIso03_dZ3_dT5s",1);   chain->SetBranchAddress("muon_chIso03_dZ3_dT5s",&muon_chIso03_dZ3_dT5s);

  chain->SetBranchStatus("muon_chIso03_dZ10",1);       chain->SetBranchAddress("muon_chIso03_dZ10",    &muon_chIso03_dZ10);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT2s",1);   chain->SetBranchAddress("muon_chIso03_dZ10_dT2s",&muon_chIso03_dZ10_dT2s);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT3s",1);   chain->SetBranchAddress("muon_chIso03_dZ10_dT3s",&muon_chIso03_dZ10_dT3s);
  chain->SetBranchStatus("muon_chIso03_dZ10_dT5s",1);   chain->SetBranchAddress("muon_chIso03_dZ10_dT5s",&muon_chIso03_dZ10_dT5s);


  chain->SetBranchStatus("muon_chIso03_reldZ",1);       chain->SetBranchAddress("muon_chIso03_reldZ",   &muon_chIso03_reldZ);
  chain->SetBranchStatus("muon_chIso03_reldZ_dT",1);    chain->SetBranchAddress("muon_chIso03_reldZ_dT",&muon_chIso03_reldZ_dT);

  chain->SetBranchStatus("muon_chIso03_dZmu05",1);       chain->SetBranchAddress("muon_chIso03_dZmu05",&muon_chIso03_dZmu05);
  chain->SetBranchStatus("muon_chIso03_dZmu05_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu05_dTmu",&muon_chIso03_dZmu05_dTmu);

  chain->SetBranchStatus("muon_chIso03_dZmu1",1);       chain->SetBranchAddress("muon_chIso03_dZmu1",&muon_chIso03_dZmu1);
  chain->SetBranchStatus("muon_chIso03_dZmu1_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu1_dTmu",&muon_chIso03_dZmu1_dTmu);

  chain->SetBranchStatus("muon_chIso03_dZmu2",1);       chain->SetBranchAddress("muon_chIso03_dZmu2",&muon_chIso03_dZmu2);
  chain->SetBranchStatus("muon_chIso03_dZmu2_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu2_dTmu",&muon_chIso03_dZmu2_dTmu);

  chain->SetBranchStatus("muon_chIso03_dZmu5",1);       chain->SetBranchAddress("muon_chIso03_dZmu5",&muon_chIso03_dZmu5);
  chain->SetBranchStatus("muon_chIso03_dZmu5_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu5_dTmu",&muon_chIso03_dZmu5_dTmu);

  chain->SetBranchStatus("muon_chIso03_dZmu10",1);       chain->SetBranchAddress("muon_chIso03_dZmu10",&muon_chIso03_dZmu10);
  chain->SetBranchStatus("muon_chIso03_dZmu10_dTmu",1);  chain->SetBranchAddress("muon_chIso03_dZmu10_dTmu",&muon_chIso03_dZmu10_dTmu);


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


  TH1F *h_muon_time = new TH1F("h_muon_time","h_muon_time",500,-1,1);
  h_muon_time->GetXaxis()->SetTitle("muon time (ns)");

  TH1F *h_muon_pt = new TH1F("h_muon_pt","h_muon_pt",200,0,200);
  h_muon_pt->GetXaxis()->SetTitle("muon p_{T} (GeV)");
  
  TH1F *h_muon_eta = new TH1F("h_muon_eta","h_muon_eta",100,-3,3);
  h_muon_eta->GetXaxis()->SetTitle("muon eta");

  TH1F *h_muon_phi = new TH1F("h_muon_phi","h_muon_phi",100,-4,4);
  h_muon_phi->GetXaxis()->SetTitle("muon phi");

  TH2F *h2_muon_pt_vs_eta = new TH2F("h2_muon_pt_vs_eta","h2_muon_pt_vs_eta",150, -3, 3, 500, 0, 1000);
  
  // --- chIso wrt to sim vertex
 
  // -- dz = 0.05 
  TH1F *h_muon_relChIso03_dZ05_simVtx = new TH1F("h_muon_relChIso03_dZ05_simVtx","h_muon_relChIso03_dZ05_simVtx",5000,0,5);
  h_muon_relChIso03_dZ05_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s_simVtx = new TH1F("h_muon_relChIso03_dZ05_dT2s_simVtx","h_muon_relChIso03_dZ05_dT2s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s_simVtx = new TH1F("h_muon_relChIso03_dZ05_dT3s_simVtx","h_muon_relChIso03_dZ05_dT3s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s_simVtx = new TH1F("h_muon_relChIso03_dZ05_dT5s_simVtx","h_muon_relChIso03_dZ05_dT5s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ05_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ05_simVtx_barrel","h_muon_relChIso03_dZ05_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ05_dT2s_simVtx_barrel","h_muon_relChIso03_dZ05_dT2s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ05_dT3s_simVtx_barrel","h_muon_relChIso03_dZ05_dT3s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ05_dT5s_simVtx_barrel","h_muon_relChIso03_dZ05_dT5s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ05_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ05_simVtx_endcap","h_muon_relChIso03_dZ05_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ05_dT2s_simVtx_endcap","h_muon_relChIso03_dZ05_dT2s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ05_dT3s_simVtx_endcap","h_muon_relChIso03_dZ05_dT3s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ05_dT5s_simVtx_endcap","h_muon_relChIso03_dZ05_dT5s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  // -- dz = 1. mm
  TH1F *h_muon_relChIso03_dZ1_simVtx = new TH1F("h_muon_relChIso03_dZ1_simVtx","h_muon_relChIso03_dZ1_simVtx",5000,0,5);
  h_muon_relChIso03_dZ1_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s_simVtx = new TH1F("h_muon_relChIso03_dZ1_dT2s_simVtx","h_muon_relChIso03_dZ1_dT2s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s_simVtx = new TH1F("h_muon_relChIso03_dZ1_dT3s_simVtx","h_muon_relChIso03_dZ1_dT3s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s_simVtx = new TH1F("h_muon_relChIso03_dZ1_dT5s_simVtx","h_muon_relChIso03_dZ1_dT5s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ1_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ1_simVtx_barrel","h_muon_relChIso03_dZ1_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ1_dT2s_simVtx_barrel","h_muon_relChIso03_dZ1_dT2s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ1_dT3s_simVtx_barrel","h_muon_relChIso03_dZ1_dT3s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ1_dT5s_simVtx_barrel","h_muon_relChIso03_dZ1_dT5s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ1_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ1_simVtx_endcap","h_muon_relChIso03_dZ1_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ1_dT2s_simVtx_endcap","h_muon_relChIso03_dZ1_dT2s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ1_dT3s_simVtx_endcap","h_muon_relChIso03_dZ1_dT3s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ1_dT5s_simVtx_endcap","h_muon_relChIso03_dZ1_dT5s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 2. mm
  TH1F *h_muon_relChIso03_dZ2_simVtx = new TH1F("h_muon_relChIso03_dZ2_simVtx","h_muon_relChIso03_dZ2_simVtx",5000,0,5);
  h_muon_relChIso03_dZ2_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s_simVtx = new TH1F("h_muon_relChIso03_dZ2_dT2s_simVtx","h_muon_relChIso03_dZ2_dT2s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s_simVtx = new TH1F("h_muon_relChIso03_dZ2_dT3s_simVtx","h_muon_relChIso03_dZ2_dT3s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s_simVtx = new TH1F("h_muon_relChIso03_dZ2_dT5s_simVtx","h_muon_relChIso03_dZ2_dT5s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ2_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ2_simVtx_barrel","h_muon_relChIso03_dZ2_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ2_dT2s_simVtx_barrel","h_muon_relChIso03_dZ2_dT2s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ2_dT3s_simVtx_barrel","h_muon_relChIso03_dZ2_dT3s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ2_dT5s_simVtx_barrel","h_muon_relChIso03_dZ2_dT5s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ2_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ2_simVtx_endcap","h_muon_relChIso03_dZ2_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ2_dT2s_simVtx_endcap","h_muon_relChIso03_dZ2_dT2s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ2_dT3s_simVtx_endcap","h_muon_relChIso03_dZ2_dT3s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ2_dT5s_simVtx_endcap","h_muon_relChIso03_dZ2_dT5s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 3. mm
  TH1F *h_muon_relChIso03_dZ3_simVtx = new TH1F("h_muon_relChIso03_dZ3_simVtx","h_muon_relChIso03_dZ3_simVtx",5000,0,5);
  h_muon_relChIso03_dZ3_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s_simVtx = new TH1F("h_muon_relChIso03_dZ3_dT2s_simVtx","h_muon_relChIso03_dZ3_dT2s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s_simVtx = new TH1F("h_muon_relChIso03_dZ3_dT3s_simVtx","h_muon_relChIso03_dZ3_dT3s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s_simVtx = new TH1F("h_muon_relChIso03_dZ3_dT5s_simVtx","h_muon_relChIso03_dZ3_dT5s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ3_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ3_simVtx_barrel","h_muon_relChIso03_dZ3_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ3_dT2s_simVtx_barrel","h_muon_relChIso03_dZ3_dT2s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ3_dT3s_simVtx_barrel","h_muon_relChIso03_dZ3_dT3s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ3_dT5s_simVtx_barrel","h_muon_relChIso03_dZ3_dT5s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ3_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ3_simVtx_endcap","h_muon_relChIso03_dZ3_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ3_dT2s_simVtx_endcap","h_muon_relChIso03_dZ3_dT2s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ3_dT3s_simVtx_endcap","h_muon_relChIso03_dZ3_dT3s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ3_dT5s_simVtx_endcap","h_muon_relChIso03_dZ3_dT5s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 10. mm
  TH1F *h_muon_relChIso03_dZ10_simVtx = new TH1F("h_muon_relChIso03_dZ10_simVtx","h_muon_relChIso03_dZ10_simVtx",5000,0,5);
  h_muon_relChIso03_dZ10_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s_simVtx = new TH1F("h_muon_relChIso03_dZ10_dT2s_simVtx","h_muon_relChIso03_dZ10_dT2s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s_simVtx = new TH1F("h_muon_relChIso03_dZ10_dT3s_simVtx","h_muon_relChIso03_dZ10_dT3s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s_simVtx->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s_simVtx = new TH1F("h_muon_relChIso03_dZ10_dT5s_simVtx","h_muon_relChIso03_dZ10_dT5s_simVtx",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s_simVtx->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ10_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ10_simVtx_barrel","h_muon_relChIso03_dZ10_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ10_dT2s_simVtx_barrel","h_muon_relChIso03_dZ10_dT2s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ10_dT3s_simVtx_barrel","h_muon_relChIso03_dZ10_dT3s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s_simVtx_barrel = new TH1F("h_muon_relChIso03_dZ10_dT5s_simVtx_barrel","h_muon_relChIso03_dZ10_dT5s_simVtx_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ10_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ10_simVtx_endcap","h_muon_relChIso03_dZ10_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ10_dT2s_simVtx_endcap","h_muon_relChIso03_dZ10_dT2s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ10_dT3s_simVtx_endcap","h_muon_relChIso03_dZ10_dT3s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s_simVtx_endcap = new TH1F("h_muon_relChIso03_dZ10_dT5s_simVtx_endcap","h_muon_relChIso03_dZ10_dT5s_simVtx_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");
  

  // --- reco vtx closest to sim vertex
 
  // -- dz = 0.05 
  TH1F *h_muon_relChIso03_dZ05 = new TH1F("h_muon_relChIso03_dZ05","h_muon_relChIso03_dZ05",5000,0,5);
  h_muon_relChIso03_dZ05->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s = new TH1F("h_muon_relChIso03_dZ05_dT2s","h_muon_relChIso03_dZ05_dT2s",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s = new TH1F("h_muon_relChIso03_dZ05_dT3s","h_muon_relChIso03_dZ05_dT3s",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s = new TH1F("h_muon_relChIso03_dZ05_dT5s","h_muon_relChIso03_dZ05_dT5s",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ05_barrel = new TH1F("h_muon_relChIso03_dZ05_barrel","h_muon_relChIso03_dZ05_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s_barrel = new TH1F("h_muon_relChIso03_dZ05_dT2s_barrel","h_muon_relChIso03_dZ05_dT2s_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ05_dT3s_barrel","h_muon_relChIso03_dZ05_dT3s_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s_barrel = new TH1F("h_muon_relChIso03_dZ05_dT5s_barrel","h_muon_relChIso03_dZ05_dT5s_barrel",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ05_endcap = new TH1F("h_muon_relChIso03_dZ05_endcap","h_muon_relChIso03_dZ05_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT2s_endcap = new TH1F("h_muon_relChIso03_dZ05_dT2s_endcap","h_muon_relChIso03_dZ05_dT2s_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ05_dT3s_endcap","h_muon_relChIso03_dZ05_dT3s_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ05_dT5s_endcap = new TH1F("h_muon_relChIso03_dZ05_dT5s_endcap","h_muon_relChIso03_dZ05_dT5s_endcap",5000,0,5);
  h_muon_relChIso03_dZ05_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  // -- dz = 1. mm
  TH1F *h_muon_relChIso03_dZ1 = new TH1F("h_muon_relChIso03_dZ1","h_muon_relChIso03_dZ1",5000,0,5);
  h_muon_relChIso03_dZ1->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s = new TH1F("h_muon_relChIso03_dZ1_dT2s","h_muon_relChIso03_dZ1_dT2s",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s = new TH1F("h_muon_relChIso03_dZ1_dT3s","h_muon_relChIso03_dZ1_dT3s",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s = new TH1F("h_muon_relChIso03_dZ1_dT5s","h_muon_relChIso03_dZ1_dT5s",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ1_barrel = new TH1F("h_muon_relChIso03_dZ1_barrel","h_muon_relChIso03_dZ1_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s_barrel = new TH1F("h_muon_relChIso03_dZ1_dT2s_barrel","h_muon_relChIso03_dZ1_dT2s_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ1_dT3s_barrel","h_muon_relChIso03_dZ1_dT3s_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s_barrel = new TH1F("h_muon_relChIso03_dZ1_dT5s_barrel","h_muon_relChIso03_dZ1_dT5s_barrel",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ1_endcap = new TH1F("h_muon_relChIso03_dZ1_endcap","h_muon_relChIso03_dZ1_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT2s_endcap = new TH1F("h_muon_relChIso03_dZ1_dT2s_endcap","h_muon_relChIso03_dZ1_dT2s_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ1_dT3s_endcap","h_muon_relChIso03_dZ1_dT3s_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ1_dT5s_endcap = new TH1F("h_muon_relChIso03_dZ1_dT5s_endcap","h_muon_relChIso03_dZ1_dT5s_endcap",5000,0,5);
  h_muon_relChIso03_dZ1_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 2. mm
  TH1F *h_muon_relChIso03_dZ2 = new TH1F("h_muon_relChIso03_dZ2","h_muon_relChIso03_dZ2",5000,0,5);
  h_muon_relChIso03_dZ2->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s = new TH1F("h_muon_relChIso03_dZ2_dT2s","h_muon_relChIso03_dZ2_dT2s",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s = new TH1F("h_muon_relChIso03_dZ2_dT3s","h_muon_relChIso03_dZ2_dT3s",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s = new TH1F("h_muon_relChIso03_dZ2_dT5s","h_muon_relChIso03_dZ2_dT5s",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ2_barrel = new TH1F("h_muon_relChIso03_dZ2_barrel","h_muon_relChIso03_dZ2_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s_barrel = new TH1F("h_muon_relChIso03_dZ2_dT2s_barrel","h_muon_relChIso03_dZ2_dT2s_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ2_dT3s_barrel","h_muon_relChIso03_dZ2_dT3s_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s_barrel = new TH1F("h_muon_relChIso03_dZ2_dT5s_barrel","h_muon_relChIso03_dZ2_dT5s_barrel",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ2_endcap = new TH1F("h_muon_relChIso03_dZ2_endcap","h_muon_relChIso03_dZ2_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT2s_endcap = new TH1F("h_muon_relChIso03_dZ2_dT2s_endcap","h_muon_relChIso03_dZ2_dT2s_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ2_dT3s_endcap","h_muon_relChIso03_dZ2_dT3s_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ2_dT5s_endcap = new TH1F("h_muon_relChIso03_dZ2_dT5s_endcap","h_muon_relChIso03_dZ2_dT5s_endcap",5000,0,5);
  h_muon_relChIso03_dZ2_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");



  // -- dz = 3. mm
  TH1F *h_muon_relChIso03_dZ3 = new TH1F("h_muon_relChIso03_dZ3","h_muon_relChIso03_dZ3",5000,0,5);
  h_muon_relChIso03_dZ3->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s = new TH1F("h_muon_relChIso03_dZ3_dT2s","h_muon_relChIso03_dZ3_dT2s",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s = new TH1F("h_muon_relChIso03_dZ3_dT3s","h_muon_relChIso03_dZ3_dT3s",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s = new TH1F("h_muon_relChIso03_dZ3_dT5s","h_muon_relChIso03_dZ3_dT5s",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ3_barrel = new TH1F("h_muon_relChIso03_dZ3_barrel","h_muon_relChIso03_dZ3_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s_barrel = new TH1F("h_muon_relChIso03_dZ3_dT2s_barrel","h_muon_relChIso03_dZ3_dT2s_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ3_dT3s_barrel","h_muon_relChIso03_dZ3_dT3s_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s_barrel = new TH1F("h_muon_relChIso03_dZ3_dT5s_barrel","h_muon_relChIso03_dZ3_dT5s_barrel",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ3_endcap = new TH1F("h_muon_relChIso03_dZ3_endcap","h_muon_relChIso03_dZ3_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT2s_endcap = new TH1F("h_muon_relChIso03_dZ3_dT2s_endcap","h_muon_relChIso03_dZ3_dT2s_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ3_dT3s_endcap","h_muon_relChIso03_dZ3_dT3s_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ3_dT5s_endcap = new TH1F("h_muon_relChIso03_dZ3_dT5s_endcap","h_muon_relChIso03_dZ3_dT5s_endcap",5000,0,5);
  h_muon_relChIso03_dZ3_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // -- dz = 10. mm
  TH1F *h_muon_relChIso03_dZ10 = new TH1F("h_muon_relChIso03_dZ10","h_muon_relChIso03_dZ10",5000,0,5);
  h_muon_relChIso03_dZ10->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s = new TH1F("h_muon_relChIso03_dZ10_dT2s","h_muon_relChIso03_dZ10_dT2s",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s = new TH1F("h_muon_relChIso03_dZ10_dT3s","h_muon_relChIso03_dZ10_dT3s",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s = new TH1F("h_muon_relChIso03_dZ10_dT5s","h_muon_relChIso03_dZ10_dT5s",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ10_barrel = new TH1F("h_muon_relChIso03_dZ10_barrel","h_muon_relChIso03_dZ10_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s_barrel = new TH1F("h_muon_relChIso03_dZ10_dT2s_barrel","h_muon_relChIso03_dZ10_dT2s_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s_barrel = new TH1F("h_muon_relChIso03_dZ10_dT3s_barrel","h_muon_relChIso03_dZ10_dT3s_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s_barrel = new TH1F("h_muon_relChIso03_dZ10_dT5s_barrel","h_muon_relChIso03_dZ10_dT5s_barrel",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s_barrel->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZ10_endcap = new TH1F("h_muon_relChIso03_dZ10_endcap","h_muon_relChIso03_dZ10_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT2s_endcap = new TH1F("h_muon_relChIso03_dZ10_dT2s_endcap","h_muon_relChIso03_dZ10_dT2s_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT2s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT3s_endcap = new TH1F("h_muon_relChIso03_dZ10_dT3s_endcap","h_muon_relChIso03_dZ10_dT3s_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT3s_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZ10_dT5s_endcap = new TH1F("h_muon_relChIso03_dZ10_dT5s_endcap","h_muon_relChIso03_dZ10_dT5s_endcap",5000,0,5);
  h_muon_relChIso03_dZ10_dT5s_endcap->GetXaxis()->SetTitle("relative charged isolation");


  // --- rel dZ ( 3*sigma_Z)
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



  // --- dz, dt wrt muon
  TH1F *h_muon_relChIso03_dZmu05 = new TH1F("h_muon_relChIso03_dZmu05","h_muon_relChIso03_dZmu05",5000,0,5);
  h_muon_relChIso03_dZmu05->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu05_dTmu = new TH1F("h_muon_relChIso03_dZmu05_dTmu","h_muon_relChIso03_dZmu05_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu05_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu05_barrel = new TH1F("h_muon_relChIso03_dZmu05_barrel","h_muon_relChIso03_dZmu05_barrel",5000,0,5);
  h_muon_relChIso03_dZmu05_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu05_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu05_dTmu_barrel","h_muon_relChIso03_dZmu05_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu05_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu05_endcap = new TH1F("h_muon_relChIso03_dZmu05_endcap","h_muon_relChIso03_dZmu05_endcap",5000,0,5);
  h_muon_relChIso03_dZmu05_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu05_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu05_dTmu_endcap","h_muon_relChIso03_dZmu05_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu05_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZmu1 = new TH1F("h_muon_relChIso03_dZmu1","h_muon_relChIso03_dZmu1",5000,0,5);
  h_muon_relChIso03_dZmu1->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu1_dTmu = new TH1F("h_muon_relChIso03_dZmu1_dTmu","h_muon_relChIso03_dZmu1_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu1_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu1_barrel = new TH1F("h_muon_relChIso03_dZmu1_barrel","h_muon_relChIso03_dZmu1_barrel",5000,0,5);
  h_muon_relChIso03_dZmu1_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu1_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu1_dTmu_barrel","h_muon_relChIso03_dZmu1_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu1_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu1_endcap = new TH1F("h_muon_relChIso03_dZmu1_endcap","h_muon_relChIso03_dZmu1_endcap",5000,0,5);
  h_muon_relChIso03_dZmu1_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu1_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu1_dTmu_endcap","h_muon_relChIso03_dZmu1_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu1_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");



  TH1F *h_muon_relChIso03_dZmu2 = new TH1F("h_muon_relChIso03_dZmu2","h_muon_relChIso03_dZmu2",5000,0,5);
  h_muon_relChIso03_dZmu2->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu2_dTmu = new TH1F("h_muon_relChIso03_dZmu2_dTmu","h_muon_relChIso03_dZmu2_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu2_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu2_barrel = new TH1F("h_muon_relChIso03_dZmu2_barrel","h_muon_relChIso03_dZmu2_barrel",5000,0,5);
  h_muon_relChIso03_dZmu2_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu2_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu2_dTmu_barrel","h_muon_relChIso03_dZmu2_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu2_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu2_endcap = new TH1F("h_muon_relChIso03_dZmu2_endcap","h_muon_relChIso03_dZmu2_endcap",5000,0,5);
  h_muon_relChIso03_dZmu2_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu2_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu2_dTmu_endcap","h_muon_relChIso03_dZmu2_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu2_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");



  TH1F *h_muon_relChIso03_dZmu5 = new TH1F("h_muon_relChIso03_dZmu5","h_muon_relChIso03_dZmu5",5000,0,5);
  h_muon_relChIso03_dZmu5->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu5_dTmu = new TH1F("h_muon_relChIso03_dZmu5_dTmu","h_muon_relChIso03_dZmu5_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu5_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu5_barrel = new TH1F("h_muon_relChIso03_dZmu5_barrel","h_muon_relChIso03_dZmu5_barrel",5000,0,5);
  h_muon_relChIso03_dZmu5_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu5_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu5_dTmu_barrel","h_muon_relChIso03_dZmu5_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu5_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu5_endcap = new TH1F("h_muon_relChIso03_dZmu5_endcap","h_muon_relChIso03_dZmu5_endcap",5000,0,5);
  h_muon_relChIso03_dZmu5_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu5_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu5_dTmu_endcap","h_muon_relChIso03_dZmu5_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu5_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");


  TH1F *h_muon_relChIso03_dZmu10 = new TH1F("h_muon_relChIso03_dZmu10","h_muon_relChIso03_dZmu10",5000,0,5);
  h_muon_relChIso03_dZmu10->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu10_dTmu = new TH1F("h_muon_relChIso03_dZmu10_dTmu","h_muon_relChIso03_dZmu10_dTmu",5000,0,5);
  h_muon_relChIso03_dZmu10_dTmu->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu10_barrel = new TH1F("h_muon_relChIso03_dZmu10_barrel","h_muon_relChIso03_dZmu10_barrel",5000,0,5);
  h_muon_relChIso03_dZmu10_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu10_dTmu_barrel = new TH1F("h_muon_relChIso03_dZmu10_dTmu_barrel","h_muon_relChIso03_dZmu10_dTmu_barrel",5000,0,5);
  h_muon_relChIso03_dZmu10_dTmu_barrel->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu10_endcap = new TH1F("h_muon_relChIso03_dZmu10_endcap","h_muon_relChIso03_dZmu10_endcap",5000,0,5);
  h_muon_relChIso03_dZmu10_endcap->GetXaxis()->SetTitle("relative charged isolation");

  TH1F *h_muon_relChIso03_dZmu10_dTmu_endcap = new TH1F("h_muon_relChIso03_dZmu10_dTmu_endcap","h_muon_relChIso03_dZmu10_dTmu_endcap",5000,0,5);
  h_muon_relChIso03_dZmu10_dTmu_endcap->GetXaxis()->SetTitle("relative charged isolation");




  // ratio iso
  TH1F *h_muon_relChIso03_ratio = new TH1F("h_muon_relChIso03_ratio","h_muon_relChIso03_ratio",1000,0,2);
  h_muon_relChIso03_ratio->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_ratio_barrel = new TH1F("h_muon_relChIso03_ratio_barrel","h_muon_relChIso03_ratio_barrel",1000,0,2);
  h_muon_relChIso03_ratio_barrel->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

  TH1F *h_muon_relChIso03_ratio_endcap = new TH1F("h_muon_relChIso03_ratio_endcap","h_muon_relChIso03_ratio_endcap",1000,0,2);
  h_muon_relChIso03_ratio_endcap->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");


  // line density
  TH1F *h_linedensity_noRelChIsoZCut = new TH1F("h_linedensity_noRelChIsoZCut","h_linedensity_noRelChIsoZCut", 100, 0, 2);
  TH1F *h_linedensity_noRelChIsoZTCut = new TH1F("h_linedensity_noRelChIsoZTCut","h_linedensity_noRelChIsoZTCut", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut = new TH1F("h_linedensity_RelChIsoZCut","h_linedensity_RelChIsoZCut ", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut  = new TH1F("h_linedensity_RelChIsoZTCut","h_linedensity_RelChIsoZTCut", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_simVtx = new TH1F("h_linedensity_RelChIsoZCut_simVtx","h_linedensity_RelChIsoZCut_simVtx", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_simVtx  = new TH1F("h_linedensity_RelChIsoZTCut_simVtx","h_linedensity_RelChIsoZTCut_simVtx", 100, 0, 2);


  TH1F *h_linedensity_noRelChIsoZCut_barrel = new TH1F("h_linedensity_noRelChIsoZCut_barrel","h_linedensity_noRelChIsoZCut_barrel", 100, 0, 2);
  TH1F *h_linedensity_noRelChIsoZTCut_barrel = new TH1F("h_linedensity_noRelChIsoZTCut_barrel","h_linedensity_noRelChIsoZTCut_barrel", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_barrel = new TH1F("h_linedensity_RelChIsoZCut_barrel","h_linedensity_RelChIsoZCut_barrel", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_barrel  = new TH1F("h_linedensity_RelChIsoZTCut_barrel","h_linedensity_RelChIsoZTCut_barrel", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_simVtx_barrel = new TH1F("h_linedensity_RelChIsoZCut_simVtx_barrel","h_linedensity_RelChIsoZCut_simVtx_barrel", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_simVtx_barrel  = new TH1F("h_linedensity_RelChIsoZTCut_simVtx_barrel","h_linedensity_RelChIsoZTCut_simVtx_barrel", 100, 0, 2);


  TH1F *h_linedensity_noRelChIsoZCut_endcap = new TH1F("h_linedensity_noRelChIsoZCut_endcap","h_linedensity_noRelChIsoZCut_endcap", 100, 0, 2);
  TH1F *h_linedensity_noRelChIsoZTCut_endcap = new TH1F("h_linedensity_noRelChIsoZTCut_endcap","h_linedensity_noRelChIsoZTCut_endcap", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_endcap = new TH1F("h_linedensity_RelChIsoZCut_endcap","h_linedensity_RelChIsoZCut_endcap", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_endcap  = new TH1F("h_linedensity_RelChIsoZTCut_endcap","h_linedensity_RelChIsoZTCut_endcap", 100, 0, 2);

  TH1F *h_linedensity_RelChIsoZCut_simVtx_endcap = new TH1F("h_linedensity_RelChIsoZCut_simVtx_endcap","h_linedensity_RelChIsoZCut_simVtx_endcap", 100, 0, 2);
  TH1F *h_linedensity_RelChIsoZTCut_simVtx_endcap  = new TH1F("h_linedensity_RelChIsoZTCut_simVtx_endcap","h_linedensity_RelChIsoZTCut_simVtx_endcap", 100, 0, 2);


  


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
      
      // -- test !!!!!!!!!!!!!!!
      //bool pass3D = fabs(muon_dz3D->at(imu)) < maxdz && fabs(muon_dxy3D->at(imu)) < maxdxy && !vtx3D_isFake  && fabs(vtx3D_z) < 2.0;
      //bool pass4D = fabs(muon_dz4D->at(imu)) < maxdz && fabs(muon_dxy4D->at(imu)) < maxdxy && !vtx4D_isFake  && fabs(vtx4D_z) < 2.0;

      if (nMuonsInEvent == 1){

      	if (pass3D) h_vtx_dz3D->Fill( vtx3D_z - vtxGen_z);
      	if (pass4D) h_vtx_dz4D->Fill( vtx4D_z - vtxGen_z);
      	if (pass4D) h_vtx_dt4D->Fill( vtx4D_t - vtxGen_t*1000000000.);

      	if (pass3D) h_vtx_dz3D_pull->Fill( (vtx3D_z - vtxGen_z) / vtx3D_zErr );
	if (pass4D) h_vtx_dz4D_pull->Fill( (vtx4D_z - vtxGen_z) / vtx4D_zErr );
      	if (pass4D) h_vtx_dt4D_pull->Fill( (vtx4D_t - vtxGen_t*1000000000.) / vtx4D_tErr);

      }


      h_muon_time->Fill(muon_t->at(imu), w);
      h_muon_pt->Fill(pt, w);
      h_muon_eta->Fill(muon_eta->at(imu), w);
      h_muon_phi->Fill(muon_phi->at(imu), w);
      h2_muon_pt_vs_eta->Fill(muon_eta->at(imu), pt, w);


      // control plots
      if (pass3D && pass4D){	
	float chIsoRatio = muon_chIso03_dZ1_dT3s->at(imu)/muon_chIso03_dZ1->at(imu);
	if ( muon_chIso03_dZ1->at(imu) == 0 ) chIsoRatio = 1;
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
	h_muon_relChIso03_dZ05_simVtx -> Fill(muon_chIso03_dZ05_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ1_simVtx -> Fill(muon_chIso03_dZ1_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_simVtx -> Fill(muon_chIso03_dZ2_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ3_simVtx -> Fill(muon_chIso03_dZ3_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ10_simVtx -> Fill(muon_chIso03_dZ10_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ05 -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
	h_muon_relChIso03_dZ1  -> Fill(muon_chIso03_dZ1->at(imu)/pt, w);
	h_muon_relChIso03_dZ2  -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
	h_muon_relChIso03_dZ3  -> Fill(muon_chIso03_dZ3->at(imu)/pt, w);
	h_muon_relChIso03_dZ10 -> Fill(muon_chIso03_dZ10->at(imu)/pt, w);
	h_muon_relChIso03_reldZ -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
      }
      if ( pass4D ){
	h_muon_relChIso03_dZ05_dT2s_simVtx -> Fill(muon_chIso03_dZ05_dT2s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ05_dT3s_simVtx -> Fill(muon_chIso03_dZ05_dT3s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ05_dT5s_simVtx -> Fill(muon_chIso03_dZ05_dT5s_simVtx->at(imu)/pt, w);

	h_muon_relChIso03_dZ1_dT2s_simVtx -> Fill(muon_chIso03_dZ1_dT2s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ1_dT3s_simVtx -> Fill(muon_chIso03_dZ1_dT3s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ1_dT5s_simVtx -> Fill(muon_chIso03_dZ1_dT5s_simVtx->at(imu)/pt, w);

	h_muon_relChIso03_dZ2_dT2s_simVtx -> Fill(muon_chIso03_dZ2_dT2s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_dT3s_simVtx -> Fill(muon_chIso03_dZ2_dT3s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_dT5s_simVtx -> Fill(muon_chIso03_dZ2_dT5s_simVtx->at(imu)/pt, w);

	h_muon_relChIso03_dZ3_dT2s_simVtx -> Fill(muon_chIso03_dZ3_dT2s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ3_dT3s_simVtx -> Fill(muon_chIso03_dZ3_dT3s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ3_dT5s_simVtx -> Fill(muon_chIso03_dZ3_dT5s_simVtx->at(imu)/pt, w);

	h_muon_relChIso03_dZ10_dT2s_simVtx -> Fill(muon_chIso03_dZ10_dT2s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ10_dT3s_simVtx -> Fill(muon_chIso03_dZ10_dT3s_simVtx->at(imu)/pt, w);
	h_muon_relChIso03_dZ10_dT5s_simVtx -> Fill(muon_chIso03_dZ10_dT5s_simVtx->at(imu)/pt, w);

	h_muon_relChIso03_dZ05_dT2s -> Fill(muon_chIso03_dZ05_dT2s->at(imu)/pt, w);
	h_muon_relChIso03_dZ05_dT3s -> Fill(muon_chIso03_dZ05_dT3s->at(imu)/pt, w);
	h_muon_relChIso03_dZ05_dT5s -> Fill(muon_chIso03_dZ05_dT5s->at(imu)/pt, w);

	h_muon_relChIso03_dZ1_dT2s -> Fill(muon_chIso03_dZ1_dT2s->at(imu)/pt, w);
	h_muon_relChIso03_dZ1_dT3s -> Fill(muon_chIso03_dZ1_dT3s->at(imu)/pt, w);
	h_muon_relChIso03_dZ1_dT5s -> Fill(muon_chIso03_dZ1_dT5s->at(imu)/pt, w);

	h_muon_relChIso03_dZ2_dT2s -> Fill(muon_chIso03_dZ2_dT2s->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_dT3s -> Fill(muon_chIso03_dZ2_dT3s->at(imu)/pt, w);
	h_muon_relChIso03_dZ2_dT5s -> Fill(muon_chIso03_dZ2_dT5s->at(imu)/pt, w);

	h_muon_relChIso03_dZ3_dT2s -> Fill(muon_chIso03_dZ3_dT2s->at(imu)/pt, w);
	h_muon_relChIso03_dZ3_dT3s -> Fill(muon_chIso03_dZ3_dT3s->at(imu)/pt, w);
	h_muon_relChIso03_dZ3_dT5s -> Fill(muon_chIso03_dZ3_dT5s->at(imu)/pt, w);

	h_muon_relChIso03_dZ10_dT2s -> Fill(muon_chIso03_dZ10_dT2s->at(imu)/pt, w);
	h_muon_relChIso03_dZ10_dT3s -> Fill(muon_chIso03_dZ10_dT3s->at(imu)/pt, w);
	h_muon_relChIso03_dZ10_dT5s -> Fill(muon_chIso03_dZ10_dT5s->at(imu)/pt, w);

	h_muon_relChIso03_reldZ_dT -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
      }

      // barrel 
      if (fabs(muon_eta->at(imu))<1.5){
	if ( pass3D ) {
	  h_muon_relChIso03_dZ05_simVtx_barrel -> Fill(muon_chIso03_dZ05_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_simVtx_barrel -> Fill(muon_chIso03_dZ1_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_simVtx_barrel -> Fill(muon_chIso03_dZ2_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_simVtx_barrel -> Fill(muon_chIso03_dZ3_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_simVtx_barrel -> Fill(muon_chIso03_dZ10_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_barrel -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_barrel  -> Fill(muon_chIso03_dZ1->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_barrel  -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_barrel  -> Fill(muon_chIso03_dZ3->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_barrel -> Fill(muon_chIso03_dZ10->at(imu)/pt, w);
	  h_muon_relChIso03_reldZ_barrel -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
	}
	if ( pass4D ){
	  h_muon_relChIso03_dZ05_dT2s_simVtx_barrel -> Fill(muon_chIso03_dZ05_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT3s_simVtx_barrel -> Fill(muon_chIso03_dZ05_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT5s_simVtx_barrel -> Fill(muon_chIso03_dZ05_dT5s_simVtx->at(imu)/pt, w);
	  
	  h_muon_relChIso03_dZ1_dT2s_simVtx_barrel -> Fill(muon_chIso03_dZ1_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT3s_simVtx_barrel -> Fill(muon_chIso03_dZ1_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT5s_simVtx_barrel -> Fill(muon_chIso03_dZ1_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ2_dT2s_simVtx_barrel -> Fill(muon_chIso03_dZ2_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT3s_simVtx_barrel -> Fill(muon_chIso03_dZ2_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT5s_simVtx_barrel -> Fill(muon_chIso03_dZ2_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ3_dT2s_simVtx_barrel -> Fill(muon_chIso03_dZ3_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT3s_simVtx_barrel -> Fill(muon_chIso03_dZ3_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT5s_simVtx_barrel -> Fill(muon_chIso03_dZ3_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ10_dT2s_simVtx_barrel -> Fill(muon_chIso03_dZ10_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT3s_simVtx_barrel -> Fill(muon_chIso03_dZ10_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT5s_simVtx_barrel -> Fill(muon_chIso03_dZ10_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ05_dT2s_barrel -> Fill(muon_chIso03_dZ05_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT3s_barrel -> Fill(muon_chIso03_dZ05_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT5s_barrel -> Fill(muon_chIso03_dZ05_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ1_dT2s_barrel -> Fill(muon_chIso03_dZ1_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT3s_barrel -> Fill(muon_chIso03_dZ1_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT5s_barrel -> Fill(muon_chIso03_dZ1_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ2_dT2s_barrel -> Fill(muon_chIso03_dZ2_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT3s_barrel -> Fill(muon_chIso03_dZ2_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT5s_barrel -> Fill(muon_chIso03_dZ2_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ3_dT2s_barrel -> Fill(muon_chIso03_dZ3_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT3s_barrel -> Fill(muon_chIso03_dZ3_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT5s_barrel -> Fill(muon_chIso03_dZ3_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ10_dT2s_barrel -> Fill(muon_chIso03_dZ10_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT3s_barrel -> Fill(muon_chIso03_dZ10_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT5s_barrel -> Fill(muon_chIso03_dZ10_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_reldZ_dT_barrel -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
	}
      }
      // endcap
      else{
	if ( pass3D ){
	  h_muon_relChIso03_dZ05_simVtx_endcap -> Fill(muon_chIso03_dZ05_simVtx->at(imu)/pt, w);
          h_muon_relChIso03_dZ1_simVtx_endcap -> Fill(muon_chIso03_dZ1_simVtx->at(imu)/pt, w);
          h_muon_relChIso03_dZ2_simVtx_endcap -> Fill(muon_chIso03_dZ2_simVtx->at(imu)/pt, w);
          h_muon_relChIso03_dZ3_simVtx_endcap -> Fill(muon_chIso03_dZ3_simVtx->at(imu)/pt, w);
          h_muon_relChIso03_dZ10_simVtx_endcap -> Fill(muon_chIso03_dZ10_simVtx->at(imu)/pt, w);
          h_muon_relChIso03_dZ05_endcap -> Fill(muon_chIso03_dZ05->at(imu)/pt, w);
          h_muon_relChIso03_dZ1_endcap  -> Fill(muon_chIso03_dZ1->at(imu)/pt, w);
          h_muon_relChIso03_dZ2_endcap  -> Fill(muon_chIso03_dZ2->at(imu)/pt, w);
          h_muon_relChIso03_dZ3_endcap  -> Fill(muon_chIso03_dZ3->at(imu)/pt, w);
          h_muon_relChIso03_dZ10_endcap  -> Fill(muon_chIso03_dZ10->at(imu)/pt, w);
          h_muon_relChIso03_reldZ_endcap -> Fill(muon_chIso03_reldZ->at(imu)/pt, w);
	}
	if ( pass4D ) {
	  h_muon_relChIso03_dZ05_dT2s_simVtx_endcap -> Fill(muon_chIso03_dZ05_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT3s_simVtx_endcap -> Fill(muon_chIso03_dZ05_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT5s_simVtx_endcap -> Fill(muon_chIso03_dZ05_dT5s_simVtx->at(imu)/pt, w);
	  
	  h_muon_relChIso03_dZ1_dT2s_simVtx_endcap -> Fill(muon_chIso03_dZ1_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT3s_simVtx_endcap -> Fill(muon_chIso03_dZ1_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT5s_simVtx_endcap -> Fill(muon_chIso03_dZ1_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ2_dT2s_simVtx_endcap -> Fill(muon_chIso03_dZ2_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT3s_simVtx_endcap -> Fill(muon_chIso03_dZ2_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT5s_simVtx_endcap -> Fill(muon_chIso03_dZ2_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ3_dT2s_simVtx_endcap -> Fill(muon_chIso03_dZ3_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT3s_simVtx_endcap -> Fill(muon_chIso03_dZ3_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT5s_simVtx_endcap -> Fill(muon_chIso03_dZ3_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ10_dT2s_simVtx_endcap -> Fill(muon_chIso03_dZ10_dT2s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT3s_simVtx_endcap -> Fill(muon_chIso03_dZ10_dT3s_simVtx->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT5s_simVtx_endcap -> Fill(muon_chIso03_dZ10_dT5s_simVtx->at(imu)/pt, w);

	  h_muon_relChIso03_dZ05_dT2s_endcap -> Fill(muon_chIso03_dZ05_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT3s_endcap -> Fill(muon_chIso03_dZ05_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ05_dT5s_endcap -> Fill(muon_chIso03_dZ05_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ1_dT2s_endcap -> Fill(muon_chIso03_dZ1_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT3s_endcap -> Fill(muon_chIso03_dZ1_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ1_dT5s_endcap -> Fill(muon_chIso03_dZ1_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ2_dT2s_endcap -> Fill(muon_chIso03_dZ2_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT3s_endcap -> Fill(muon_chIso03_dZ2_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ2_dT5s_endcap -> Fill(muon_chIso03_dZ2_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ3_dT2s_endcap -> Fill(muon_chIso03_dZ3_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT3s_endcap -> Fill(muon_chIso03_dZ3_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ3_dT5s_endcap -> Fill(muon_chIso03_dZ3_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_dZ10_dT2s_endcap -> Fill(muon_chIso03_dZ10_dT2s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT3s_endcap -> Fill(muon_chIso03_dZ10_dT3s->at(imu)/pt, w);
	  h_muon_relChIso03_dZ10_dT5s_endcap -> Fill(muon_chIso03_dZ10_dT5s->at(imu)/pt, w);

	  h_muon_relChIso03_reldZ_dT_endcap -> Fill(muon_chIso03_reldZ_dT->at(imu)/pt, w);
	}
      }
      

      // dz, dt cuts wrt muon
      //if (muon_t->at(imu) == 0.) { continue; }
      h_muon_relChIso03_dZmu05 -> Fill(muon_chIso03_dZmu05->at(imu)/pt, w);
      h_muon_relChIso03_dZmu05_dTmu -> Fill(muon_chIso03_dZmu05_dTmu->at(imu)/pt, w);
      h_muon_relChIso03_dZmu1 -> Fill(muon_chIso03_dZmu1->at(imu)/pt, w);
      h_muon_relChIso03_dZmu1_dTmu -> Fill(muon_chIso03_dZmu1_dTmu->at(imu)/pt, w);
      h_muon_relChIso03_dZmu2 -> Fill(muon_chIso03_dZmu2->at(imu)/pt, w);
      h_muon_relChIso03_dZmu2_dTmu -> Fill(muon_chIso03_dZmu2_dTmu->at(imu)/pt, w);
      h_muon_relChIso03_dZmu5 -> Fill(muon_chIso03_dZmu5->at(imu)/pt, w);
      h_muon_relChIso03_dZmu5_dTmu -> Fill(muon_chIso03_dZmu5_dTmu->at(imu)/pt, w);
      h_muon_relChIso03_dZmu10 -> Fill(muon_chIso03_dZmu10->at(imu)/pt, w);
      h_muon_relChIso03_dZmu10_dTmu -> Fill(muon_chIso03_dZmu10_dTmu->at(imu)/pt, w);
      if (fabs(muon_eta->at(imu))<1.5){ 
	h_muon_relChIso03_dZmu05_barrel -> Fill(muon_chIso03_dZmu05->at(imu)/pt, w);
	h_muon_relChIso03_dZmu05_dTmu_barrel -> Fill(muon_chIso03_dZmu05_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu1_barrel -> Fill(muon_chIso03_dZmu1->at(imu)/pt, w);
	h_muon_relChIso03_dZmu1_dTmu_barrel -> Fill(muon_chIso03_dZmu1_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu2_barrel -> Fill(muon_chIso03_dZmu2->at(imu)/pt, w);
	h_muon_relChIso03_dZmu2_dTmu_barrel -> Fill(muon_chIso03_dZmu2_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu5_barrel -> Fill(muon_chIso03_dZmu5->at(imu)/pt, w);
	h_muon_relChIso03_dZmu5_dTmu_barrel -> Fill(muon_chIso03_dZmu5_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu10_barrel -> Fill(muon_chIso03_dZmu10->at(imu)/pt, w);
	h_muon_relChIso03_dZmu10_dTmu_barrel -> Fill(muon_chIso03_dZmu10_dTmu->at(imu)/pt, w);
      }
      else{
	h_muon_relChIso03_dZmu05_endcap -> Fill(muon_chIso03_dZmu05->at(imu)/pt, w);
	h_muon_relChIso03_dZmu05_dTmu_endcap -> Fill(muon_chIso03_dZmu05_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu1_endcap -> Fill(muon_chIso03_dZmu1->at(imu)/pt, w);
	h_muon_relChIso03_dZmu1_dTmu_endcap -> Fill(muon_chIso03_dZmu1_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu2_endcap -> Fill(muon_chIso03_dZmu2->at(imu)/pt, w);
	h_muon_relChIso03_dZmu2_dTmu_endcap -> Fill(muon_chIso03_dZmu2_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu5_endcap -> Fill(muon_chIso03_dZmu5->at(imu)/pt, w);
	h_muon_relChIso03_dZmu5_dTmu_endcap -> Fill(muon_chIso03_dZmu5_dTmu->at(imu)/pt, w);
	h_muon_relChIso03_dZmu10_endcap -> Fill(muon_chIso03_dZmu10->at(imu)/pt, w);
	h_muon_relChIso03_dZmu10_dTmu_endcap -> Fill(muon_chIso03_dZmu10_dTmu->at(imu)/pt, w);
      }
              

      // -- line density
      float linedensity = 200.*TMath::Gaus(fabs(10*vtxGen_z), 0, 42., 1);
      if (pass3D){
	h_linedensity_noRelChIsoZCut->Fill(linedensity,w); 
	if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_noRelChIsoZCut_barrel->Fill(linedensity,w);}
	else{ h_linedensity_noRelChIsoZCut_endcap->Fill(linedensity,w); }
	if ( (muon_chIso03_dZ1->at(imu)/muon_pt->at(imu)) < 0.05 ) {
	  h_linedensity_RelChIsoZCut->Fill(linedensity,w); 
	  if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_RelChIsoZCut_barrel->Fill(linedensity,w);}
	  else {h_linedensity_RelChIsoZCut_endcap->Fill(linedensity,w);}
	}
	if ( (muon_chIso03_dZ1_simVtx->at(imu)/muon_pt->at(imu)) < 0.05 ) {
	  h_linedensity_RelChIsoZCut_simVtx->Fill(linedensity,w); 
	  if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_RelChIsoZCut_simVtx_barrel->Fill(linedensity,w);}
	  else {h_linedensity_RelChIsoZCut_simVtx_endcap->Fill(linedensity,w);}
	}
      }
      if (pass4D){
	h_linedensity_noRelChIsoZTCut->Fill(linedensity,w);
	if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_noRelChIsoZTCut_barrel->Fill(linedensity,w);}
        else{ h_linedensity_noRelChIsoZTCut_endcap->Fill(linedensity,w); }
	if ( (muon_chIso03_dZ1_dT3s->at(imu)/muon_pt->at(imu)) < 0.05 ) {
	  h_linedensity_RelChIsoZTCut->Fill(linedensity,w); 
	  if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_RelChIsoZTCut_barrel->Fill(linedensity,w); }
	  else {h_linedensity_RelChIsoZTCut_endcap->Fill(linedensity,w);}
	}
	if ( (muon_chIso03_dZ1_dT3s_simVtx->at(imu)/muon_pt->at(imu)) < 0.05 ) {
	  h_linedensity_RelChIsoZTCut_simVtx->Fill(linedensity,w); 
	  if (fabs(muon_eta->at(imu)) <1.5){ h_linedensity_RelChIsoZTCut_simVtx_barrel->Fill(linedensity,w); }
          else {h_linedensity_RelChIsoZTCut_simVtx_endcap->Fill(linedensity,w);}
	}

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

  h_muon_time->Write();
  h_muon_pt->Write();
  h_muon_eta->Write();
  h_muon_phi->Write();
  h2_muon_pt_vs_eta->Write();

  h_muon_relChIso03_ratio ->Write();
  h_muon_relChIso03_ratio_barrel ->Write();
  h_muon_relChIso03_ratio_endcap ->Write();


  h_muon_relChIso03_dZ05_simVtx->Write();
  h_muon_relChIso03_dZ05_simVtx_barrel->Write();
  h_muon_relChIso03_dZ05_simVtx_endcap->Write();

  h_muon_relChIso03_dZ05_dT2s_simVtx->Write();
  h_muon_relChIso03_dZ05_dT2s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ05_dT2s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ05_dT3s_simVtx->Write();
  h_muon_relChIso03_dZ05_dT3s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ05_dT3s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ05_dT5s_simVtx->Write();
  h_muon_relChIso03_dZ05_dT5s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ05_dT5s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ1_simVtx->Write();
  h_muon_relChIso03_dZ1_simVtx_barrel->Write();
  h_muon_relChIso03_dZ1_simVtx_endcap->Write();

  h_muon_relChIso03_dZ1_dT2s_simVtx->Write();
  h_muon_relChIso03_dZ1_dT2s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ1_dT2s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ1_dT3s_simVtx->Write();
  h_muon_relChIso03_dZ1_dT3s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ1_dT3s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ1_dT5s_simVtx->Write();
  h_muon_relChIso03_dZ1_dT5s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ1_dT5s_simVtx_endcap->Write();
  
  h_muon_relChIso03_dZ2_simVtx->Write();
  h_muon_relChIso03_dZ2_simVtx_barrel->Write();
  h_muon_relChIso03_dZ2_simVtx_endcap->Write();

  h_muon_relChIso03_dZ2_dT2s_simVtx->Write();
  h_muon_relChIso03_dZ2_dT2s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ2_dT2s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ2_dT3s_simVtx->Write();
  h_muon_relChIso03_dZ2_dT3s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ2_dT3s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ2_dT5s_simVtx->Write();
  h_muon_relChIso03_dZ2_dT5s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ2_dT5s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ3_simVtx->Write();
  h_muon_relChIso03_dZ3_simVtx_barrel->Write();
  h_muon_relChIso03_dZ3_simVtx_endcap->Write();

  h_muon_relChIso03_dZ3_dT2s_simVtx->Write();
  h_muon_relChIso03_dZ3_dT2s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ3_dT2s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ3_dT3s_simVtx->Write();
  h_muon_relChIso03_dZ3_dT3s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ3_dT3s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ3_dT5s_simVtx->Write();
  h_muon_relChIso03_dZ3_dT5s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ3_dT5s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ10_simVtx->Write();
  h_muon_relChIso03_dZ10_simVtx_barrel->Write();
  h_muon_relChIso03_dZ10_simVtx_endcap->Write();

  h_muon_relChIso03_dZ10_dT2s_simVtx->Write();
  h_muon_relChIso03_dZ10_dT2s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ10_dT2s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ10_dT3s_simVtx->Write();
  h_muon_relChIso03_dZ10_dT3s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ10_dT3s_simVtx_endcap->Write();

  h_muon_relChIso03_dZ10_dT5s_simVtx->Write();
  h_muon_relChIso03_dZ10_dT5s_simVtx_barrel->Write();
  h_muon_relChIso03_dZ10_dT5s_simVtx_endcap->Write();



  h_muon_relChIso03_dZ05->Write();
  h_muon_relChIso03_dZ05_barrel->Write();
  h_muon_relChIso03_dZ05_endcap->Write();

  h_muon_relChIso03_dZ05_dT2s->Write();
  h_muon_relChIso03_dZ05_dT2s_barrel->Write();
  h_muon_relChIso03_dZ05_dT2s_endcap->Write();

  h_muon_relChIso03_dZ05_dT3s->Write();
  h_muon_relChIso03_dZ05_dT3s_barrel->Write();
  h_muon_relChIso03_dZ05_dT3s_endcap->Write();

  h_muon_relChIso03_dZ05_dT5s->Write();
  h_muon_relChIso03_dZ05_dT5s_barrel->Write();
  h_muon_relChIso03_dZ05_dT5s_endcap->Write();

  h_muon_relChIso03_dZ1->Write();
  h_muon_relChIso03_dZ1_barrel->Write();
  h_muon_relChIso03_dZ1_endcap->Write();

  h_muon_relChIso03_dZ1_dT2s->Write();
  h_muon_relChIso03_dZ1_dT2s_barrel->Write();
  h_muon_relChIso03_dZ1_dT2s_endcap->Write();

  h_muon_relChIso03_dZ1_dT3s->Write();
  h_muon_relChIso03_dZ1_dT3s_barrel->Write();
  h_muon_relChIso03_dZ1_dT3s_endcap->Write();

  h_muon_relChIso03_dZ1_dT5s->Write();
  h_muon_relChIso03_dZ1_dT5s_barrel->Write();
  h_muon_relChIso03_dZ1_dT5s_endcap->Write();
  
  h_muon_relChIso03_dZ2->Write();
  h_muon_relChIso03_dZ2_barrel->Write();
  h_muon_relChIso03_dZ2_endcap->Write();

  h_muon_relChIso03_dZ2_dT2s->Write();
  h_muon_relChIso03_dZ2_dT2s_barrel->Write();
  h_muon_relChIso03_dZ2_dT2s_endcap->Write();

  h_muon_relChIso03_dZ2_dT3s->Write();
  h_muon_relChIso03_dZ2_dT3s_barrel->Write();
  h_muon_relChIso03_dZ2_dT3s_endcap->Write();

  h_muon_relChIso03_dZ2_dT5s->Write();
  h_muon_relChIso03_dZ2_dT5s_barrel->Write();
  h_muon_relChIso03_dZ2_dT5s_endcap->Write();

  h_muon_relChIso03_dZ3->Write();
  h_muon_relChIso03_dZ3_barrel->Write();
  h_muon_relChIso03_dZ3_endcap->Write();

  h_muon_relChIso03_dZ3_dT2s->Write();
  h_muon_relChIso03_dZ3_dT2s_barrel->Write();
  h_muon_relChIso03_dZ3_dT2s_endcap->Write();

  h_muon_relChIso03_dZ3_dT3s->Write();
  h_muon_relChIso03_dZ3_dT3s_barrel->Write();
  h_muon_relChIso03_dZ3_dT3s_endcap->Write();

  h_muon_relChIso03_dZ3_dT5s->Write();
  h_muon_relChIso03_dZ3_dT5s_barrel->Write();
  h_muon_relChIso03_dZ3_dT5s_endcap->Write();

  h_muon_relChIso03_dZ10->Write();
  h_muon_relChIso03_dZ10_barrel->Write();
  h_muon_relChIso03_dZ10_endcap->Write();

  h_muon_relChIso03_dZ10_dT2s->Write();
  h_muon_relChIso03_dZ10_dT2s_barrel->Write();
  h_muon_relChIso03_dZ10_dT2s_endcap->Write();

  h_muon_relChIso03_dZ10_dT3s->Write();
  h_muon_relChIso03_dZ10_dT3s_barrel->Write();
  h_muon_relChIso03_dZ10_dT3s_endcap->Write();

  h_muon_relChIso03_dZ10_dT5s->Write();
  h_muon_relChIso03_dZ10_dT5s_barrel->Write();
  h_muon_relChIso03_dZ10_dT5s_endcap->Write();



  h_muon_relChIso03_reldZ->Write();
  h_muon_relChIso03_reldZ_barrel->Write();
  h_muon_relChIso03_reldZ_endcap->Write();

  h_muon_relChIso03_reldZ_dT->Write();
  h_muon_relChIso03_reldZ_dT_barrel->Write();
  h_muon_relChIso03_reldZ_dT_endcap->Write();




  h_muon_relChIso03_dZmu05->Write();
  h_muon_relChIso03_dZmu05_barrel->Write();
  h_muon_relChIso03_dZmu05_endcap->Write();

  h_muon_relChIso03_dZmu05_dTmu->Write();
  h_muon_relChIso03_dZmu05_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu05_dTmu_endcap->Write();

  h_muon_relChIso03_dZmu1->Write();
  h_muon_relChIso03_dZmu1_barrel->Write();
  h_muon_relChIso03_dZmu1_endcap->Write();

  h_muon_relChIso03_dZmu1_dTmu->Write();
  h_muon_relChIso03_dZmu1_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu1_dTmu_endcap->Write();

  h_muon_relChIso03_dZmu2->Write();
  h_muon_relChIso03_dZmu2_barrel->Write();
  h_muon_relChIso03_dZmu2_endcap->Write();

  h_muon_relChIso03_dZmu2_dTmu->Write();
  h_muon_relChIso03_dZmu2_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu2_dTmu_endcap->Write();

  h_muon_relChIso03_dZmu5->Write();
  h_muon_relChIso03_dZmu5_barrel->Write();
  h_muon_relChIso03_dZmu5_endcap->Write();

  h_muon_relChIso03_dZmu5_dTmu->Write();
  h_muon_relChIso03_dZmu5_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu5_dTmu_endcap->Write();

  h_muon_relChIso03_dZmu10->Write();
  h_muon_relChIso03_dZmu10_barrel->Write();
  h_muon_relChIso03_dZmu10_endcap->Write();

  h_muon_relChIso03_dZmu10_dTmu->Write();
  h_muon_relChIso03_dZmu10_dTmu_barrel->Write();
  h_muon_relChIso03_dZmu10_dTmu_endcap->Write();



  h_linedensity_noRelChIsoZCut->Write();
  h_linedensity_noRelChIsoZTCut->Write();

  h_linedensity_RelChIsoZCut->Write();
  h_linedensity_RelChIsoZTCut->Write();

  h_linedensity_RelChIsoZCut_simVtx->Write();
  h_linedensity_RelChIsoZTCut_simVtx->Write();


  h_linedensity_noRelChIsoZCut_barrel->Write();
  h_linedensity_noRelChIsoZTCut_barrel->Write();

  h_linedensity_RelChIsoZCut_barrel->Write();
  h_linedensity_RelChIsoZTCut_barrel->Write();

  h_linedensity_RelChIsoZCut_simVtx_barrel->Write();
  h_linedensity_RelChIsoZTCut_simVtx_barrel->Write();


  h_linedensity_noRelChIsoZCut_endcap->Write();
  h_linedensity_noRelChIsoZTCut_endcap->Write();

  h_linedensity_RelChIsoZCut_endcap->Write();
  h_linedensity_RelChIsoZTCut_endcap->Write();

  h_linedensity_RelChIsoZCut_simVtx_endcap->Write();
  h_linedensity_RelChIsoZTCut_simVtx_endcap->Write();

  fout->Close();
  
}


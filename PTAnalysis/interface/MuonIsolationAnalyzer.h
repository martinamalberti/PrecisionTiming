#ifndef _MUON_ISOLATION_H
#define _MUON_ISOLATION_H
// -*- C++ -*-
//
// Package:    PrecisionTiming/MuonIsolationAnalyzer
// Class:      MuonIsolationAnalyzer
// 
/**\class MuonIsolationAnalyzer MuonIsolationAnalyzer.cc PrecisionTiming/MuonIsolationAnalyzer/plugins/MuonIsolationAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Martina Malberti
//         Created:  Mon, 10 Oct 2016 14:06:02 GMT
//
//


// system include files
#include <memory>
#include <cstdlib>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/Common/interface/View.h"

#include "PrecisionTiming/PTAnalysis/interface/Utils.h"

#include <vector>
#include <map>
#include "TTree.h"
#include <TRandom.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


using namespace std;
using namespace edm;
using namespace reco;
using namespace math;


struct eventInfo
{
  int npu;
  vector<float> track_pt;
  vector<float> track_eta;
  vector<float> track_phi;
  vector<float> track_dz3D;
  vector<float> track_dz4D;
  vector<float> track_dxy3D;
  vector<float> track_dxy4D;
  vector<float> track_t;
  vector<float> track_dRmu;
  vector<int> track_muIndex;
  vector<int> track_isMatchedToGenParticle;
  vector<int> track_isUnmatchedToGenParticle;

  float vtxGen_z;
  float vtxGen_t;
  float vtx4D_z;
  float vtx3D_z;
  float vtx3D_zErr;
  float vtx4D_t;
  float vtx4D_tErr;
  float vtx4D_zErr;
  int vtx3D_isFake;
  int vtx4D_isFake;
  vector<float> muon_pt;
  vector<float> muon_eta;
  vector<float> muon_phi;
  vector<int> muon_isLoose;
  vector<int> muon_isMedium;
  vector<int> muon_isTight3D;
  vector<int> muon_isTight4D;
  vector<float> muon_dz3D;
  vector<float> muon_dxy3D;
  vector<float> muon_dz4D;
  vector<float> muon_dxy4D;
  vector<float> muon_t;
  vector<int> muon_isPrompt;
  vector<int> muon_isMatchedToGenJet;
  vector<int> muon_isFromTauDecay;

  vector<float> muon_chIso_dZ05_simVtx;
  vector<float> muon_chIso_dZ05_dT2s_simVtx;
  vector<float> muon_chIso_dZ05_dT3s_simVtx;
  vector<float> muon_chIso_dZ05_dT5s_simVtx;

  vector<float> muon_chIso_dZ1_simVtx;
  vector<float> muon_chIso_dZ1_dT2s_simVtx;
  vector<float> muon_chIso_dZ1_dT3s_simVtx;
  vector<float> muon_chIso_dZ1_dT5s_simVtx;

  vector<float> muon_chIso_dZ2_simVtx;
  vector<float> muon_chIso_dZ2_dT2s_simVtx;
  vector<float> muon_chIso_dZ2_dT3s_simVtx;
  vector<float> muon_chIso_dZ2_dT5s_simVtx;

  vector<float> muon_chIso_dZ3_simVtx;
  vector<float> muon_chIso_dZ3_dT2s_simVtx;
  vector<float> muon_chIso_dZ3_dT3s_simVtx;
  vector<float> muon_chIso_dZ3_dT5s_simVtx;

  vector<float> muon_chIso_dZ10_simVtx;
  vector<float> muon_chIso_dZ10_dT2s_simVtx;
  vector<float> muon_chIso_dZ10_dT3s_simVtx;
  vector<float> muon_chIso_dZ10_dT5s_simVtx;

  vector<float> muon_chIso_dZ05;
  vector<float> muon_chIso_dZ05_dT2s;
  vector<float> muon_chIso_dZ05_dT3s;
  vector<float> muon_chIso_dZ05_dT5s;

  vector<float> muon_chIso_dZ1;
  vector<float> muon_chIso_dZ1_dT2s;
  vector<float> muon_chIso_dZ1_dT3s;
  vector<float> muon_chIso_dZ1_dT5s;

  vector<float> muon_chIso_dZ2;
  vector<float> muon_chIso_dZ2_dT2s;
  vector<float> muon_chIso_dZ2_dT3s;
  vector<float> muon_chIso_dZ2_dT5s;

  vector<float> muon_chIso_dZ3;
  vector<float> muon_chIso_dZ3_dT2s;
  vector<float> muon_chIso_dZ3_dT3s;
  vector<float> muon_chIso_dZ3_dT5s;

  vector<float> muon_chIso_dZ10;
  vector<float> muon_chIso_dZ10_dT2s;
  vector<float> muon_chIso_dZ10_dT3s;
  vector<float> muon_chIso_dZ10_dT5s;

  vector<float> muon_chIso_reldZ;
  vector<float> muon_chIso_reldZ_dT;

  vector<float> muon_chIso_dZmu05;
  vector<float> muon_chIso_dZmu05_dTmu;

  vector<float> muon_chIso_dZmu1;
  vector<float> muon_chIso_dZmu1_dTmu;

  vector<float> muon_chIso_dZmu2;
  vector<float> muon_chIso_dZmu2_dTmu;

  vector<float> muon_chIso_dZmu5;
  vector<float> muon_chIso_dZmu5_dTmu;

  vector<float> muon_chIso_dZmu10;
  vector<float> muon_chIso_dZmu10_dTmu;

};


//class MuonIsolationAnalyzer : public edm::EDAnalyzer  
class MuonIsolationAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit MuonIsolationAnalyzer(const edm::ParameterSet&);
  ~MuonIsolationAnalyzer();

  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
    
  void initEventStructure();


  //---inputs  
  EDGetTokenT<vector<PileupSummaryInfo> > PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken3D_;
  EDGetTokenT<View<reco::Vertex> > vertexToken4D_;
  EDGetTokenT<edm::View<reco::PFCandidate> >      pfcandToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  EDGetTokenT<vector<SimVertex> >  genVertexToken_;
  EDGetTokenT<genXYZ> genXYZToken_;
  EDGetTokenT<float>  genT0Token_;
  EDGetTokenT<View<reco::GenJet> > genJetsToken_;
  EDGetTokenT<View<reco::Muon> > muonsToken_; 
  
  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree[10];
  eventInfo *evInfo[10];
  
  //--- options
  vector<double> timeResolutions_;
  bool mtd5sample_;
  double isoConeDR_;
  bool saveTracks_;
  double maxDz_;
  double minDr_;
  double btlMinTrackPt_;
  double etlMinTrackPt_;
  bool useVertexClosestToGenZ_;
  bool useVertexClosestToGenZT_;
  double btlEfficiency_;
  double etlEfficiency_;

  // -- 
  TRandom *gRandom;
  TRandom *gRandom2;

};


/*
bool isPromptMuon(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles);
bool isMatchedToGenJet(const reco::Muon &muon, const edm::View<reco::GenJet>& genJet);
bool isFromTau(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles);
bool isMatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles);
bool isUnmatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles);
*/
#endif

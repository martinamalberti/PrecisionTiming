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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/Common/interface/View.h"

#include "PrecisionTiming/PTAnalysis/interface/Utils.h"

#include <vector>
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
  vector<int> track_phoIndex;
  vector<int> track_isMatchedToGenParticle;

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
  vector<float> photon_pt;
  vector<float> photon_eta;
  vector<float> photon_phi;
  vector<float> photon_t;
  vector<float> photon_hasConversionTracks;
  vector<float> photon_hasPixelSeed;
  vector<float> photon_sigmaIetaIeta;
  vector<float> photon_r9;
  vector<int> photon_isPrompt;
  vector<int> photon_isMatchedToGenJet;

  vector<float> photon_chIso_dZ05_simVtx;
  vector<float> photon_chIso_dZ05_dT2s_simVtx;
  vector<float> photon_chIso_dZ05_dT3s_simVtx;
  vector<float> photon_chIso_dZ05_dT5s_simVtx;

  vector<float> photon_chIso_dZ1_simVtx;
  vector<float> photon_chIso_dZ1_dT2s_simVtx;
  vector<float> photon_chIso_dZ1_dT3s_simVtx;
  vector<float> photon_chIso_dZ1_dT5s_simVtx;

  vector<float> photon_chIso_dZ2_simVtx;
  vector<float> photon_chIso_dZ2_dT2s_simVtx;
  vector<float> photon_chIso_dZ2_dT3s_simVtx;
  vector<float> photon_chIso_dZ2_dT5s_simVtx;

  vector<float> photon_chIso_dZ05;
  vector<float> photon_chIso_dZ05_dT2s;
  vector<float> photon_chIso_dZ05_dT3s;
  vector<float> photon_chIso_dZ05_dT5s;

  vector<float> photon_chIso_dZ1;
  vector<float> photon_chIso_dZ1_dT2s;
  vector<float> photon_chIso_dZ1_dT3s;
  vector<float> photon_chIso_dZ1_dT5s;

  vector<float> photon_chIso_dZ2;
  vector<float> photon_chIso_dZ2_dT2s;
  vector<float> photon_chIso_dZ2_dT3s;
  vector<float> photon_chIso_dZ2_dT5s;


};


//class PhotonIsolationAnalyzer : public edm::EDAnalyzer  
class PhotonIsolationAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit PhotonIsolationAnalyzer(const edm::ParameterSet&);
  ~PhotonIsolationAnalyzer();

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
  EDGetTokenT<View<reco::Photon> > barrelPhotonsToken_;
  EDGetTokenT<View<reco::Photon> > endcapPhotonsToken_;
  EDGetTokenT<ValueMap<float> > trackTimeToken_;
  EDGetTokenT<ValueMap<float> > trackTimeErrToken_;
  
  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree[6];
  eventInfo *evInfo[6];
  
  //--- options
  vector<double> timeResolutions_;
  double isoConeDR_;
  bool saveTracks_;
  double maxDz_;
  double maxDxy_;
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

//bool isPromptPhoton(const reco::Photon &photon, const edm::View<reco::GenParticle>& genParticles);
//bool isMatchedToGenJet(const reco::Photon &photon, const edm::View<reco::GenJet>& genJet);

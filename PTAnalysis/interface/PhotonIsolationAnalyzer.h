// -*- C++ -*-
//
// Package:    PrecisionTiming/PhotonIsolationAnalyzer
// Class:      PhotonIsolationAnalyzer
// 
/**\class PhotonIsolationAnalyzer PhotonIsolationAnalyzer.cc PrecisionTiming/PhotonIsolationAnalyzer/plugins/PhotonIsolationAnalyzer.cc

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

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
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>
#include "TTree.h"

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
  vector<float> track_dz;
  vector<float> track_t;

  float vtxGen_z;
  float vtxGen_t;
  float vtx_z;
  float vtx_t;
  float vtx3D_z;
  vector<float> photon_pt;
  vector<float> photon_eta;
  vector<float> photon_phi;
  vector<float> photon_isPrompt;
  vector<float> photon_isMatchedToGenJet;
  vector<float> photon_hasConversionTracks;
  vector<float> photon_hasPixelSeed;
  vector<float> photon_sigmaIetaIeta;
  vector<float> photon_r9;
  vector<float> photon_chIso[10];
  vector<float> photon_chIso_dT[10][10];
};


class PhotonIsolationAnalyzer : public edm::EDAnalyzer  
{
public:
  explicit PhotonIsolationAnalyzer(const edm::ParameterSet&);
  ~PhotonIsolationAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  //virtual void beginJob() override;
  //virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  //virtual void endJob() override;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
    
  void initEventStructure();


  //---inputs  
  EDGetTokenT<vector<PileupSummaryInfo> > PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken3D_;
  EDGetTokenT<View<reco::Vertex> > vertexToken4D_;
  EDGetTokenT<edm::View<reco::PFCandidate> >      pfcandToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  EDGetTokenT<vector<SimVertex> >  genVertexToken_;
  EDGetTokenT<View<reco::GenJet> > genJetsToken_;
  EDGetTokenT<View<reco::Photon> > barrelPhotonsToken_; 
  EDGetTokenT<View<reco::Photon> > endcapPhotonsToken_; 
  
  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree[10];
  eventInfo evInfo[10];

  //--- options
  vector<double> timeResolutions_;
  vector<double> isoConeDR_;
  bool saveTracks_;
  float maxDz_;
  float minDr_;
};

bool isPromptPhoton(const reco::Photon &photon, const edm::View<reco::GenParticle>& genParticles);
bool isMatchedToGenJet(const reco::Photon &photon, const edm::View<reco::GenJet>& genJet);
math::XYZTLorentzVector correctP4(const reco::Photon &photon, const reco::Vertex& vtx);

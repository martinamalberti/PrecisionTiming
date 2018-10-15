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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/Common/interface/View.h"

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
  vector<float> track_dz3D;
  vector<float> track_dz;
  vector<float> track_t;

  float vtxGen_z;
  float vtxGen_t;
  float vtx_z;
  float vtx_t;
  float vtx3D_z;
  vector<float> muon_pt;
  vector<float> muon_eta;
  vector<float> muon_phi;
  vector<int> muon_isLoose;
  vector<int> muon_isMedium;
  vector<int> muon_isTight3D;
  vector<int> muon_isTight4D;
  vector<int> muon_isPrompt;
  vector<int> muon_isMatchedToGenJet;
  vector<float> muon_chIso[10];
  vector<float> muon_chIso_dT[10][10];
};


class MuonIsolationAnalyzer : public edm::EDAnalyzer  
{
public:
  explicit MuonIsolationAnalyzer(const edm::ParameterSet&);
  ~MuonIsolationAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
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
  EDGetTokenT<View<reco::Muon> > muonsToken_; 
  
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

bool isPromptMuon(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles);
bool isMatchedToGenJet(const reco::Muon &muon, const edm::View<reco::GenJet>& genJet);
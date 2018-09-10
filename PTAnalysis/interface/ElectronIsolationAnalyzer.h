// -*- C++ -*-
//
// Package:    PrecisionTiming/ElectronIsolationAnalyzer
// Class:      ElectronIsolationAnalyzer
// 
/**\class ElectronIsolationAnalyzer ElectronIsolationAnalyzer.cc PrecisionTiming/ElectronIsolationAnalyzer/plugins/ElectronIsolationAnalyzer.cc

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
#include "DataFormats/PatCandidates/interface/Electron.h"

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
  float vtx3D_z;
  float vtx_z;
  float vtx_t;
  vector<float> electron_pt;
  vector<float> electron_eta;
  vector<float> electron_phi;
  vector<float> electron_isMatchedToGen;
  vector<float> electron_r9;
  vector<float> electron_chIso[10];
  vector<float> electron_chIso_dT[10][10];
};


class ElectronIsolationAnalyzer : public edm::EDAnalyzer  
{
public:
  explicit ElectronIsolationAnalyzer(const edm::ParameterSet&);
  ~ElectronIsolationAnalyzer();

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
  EDGetTokenT<View<reco::Track> > tracksToken_;
  EDGetTokenT<ValueMap<float> > trackTimeToken_;
  EDGetTokenT<edm::View<reco::PFCandidate> >      pfcandToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  EDGetTokenT<vector<SimVertex> >  genVertexToken_;
  EDGetTokenT<View<reco::GsfElectron> > electronsToken_; 
  
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

bool isMatchedToGen(const reco::GsfElectron &electron, const edm::View<reco::GenParticle>& genParticles);

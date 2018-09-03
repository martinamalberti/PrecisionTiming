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


struct eventInfo
{
  int npu;
  vector<float> tkEta;
  vector<float> tkPhi;
  vector<float> tkPt;
  vector<float> tkTime;
  /*vector<float> tkOuterR;
  vector<float> tkOuterX;
  vector<float> tkOuterY;
  vector<float> tkOuterZ;
  */
  
  //int index_closestToGen;
  float vtx_t;
  float vtx_z;
  float vtx_nTracks;
  float vtx_nTracks_dT;
  vector<float> photon_pt;
  vector<float> photon_eta;
  vector<float> photon_phi;
  vector<float> photon_chIso[4];
  vector<float> photon_chIso_dT[4][3];
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
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  EDGetTokenT<View<reco::Track> > tracksToken_;
  EDGetTokenT<ValueMap<float> > trackTimeToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  EDGetTokenT<vector<SimVertex> >  genVertexToken_;
  EDGetTokenT<View<reco::Photon> > photonsToken_; 
  
  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree[3];
  eventInfo evInfo[4];

  //--- options
  vector<double> timeResolutions_;
  vector<double> isoConeDR_;
  bool saveTracks_;
  double BTLEfficiency_;
  double ETLEfficiency_;
};


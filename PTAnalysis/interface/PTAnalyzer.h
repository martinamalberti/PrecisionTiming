// -*- C++ -*-
//
// Package:    PrecisionTiming/PTAnalyzer
// Class:      PTAnalyzer
// 
/**\class PTAnalyzer PTAnalyzer.cc PrecisionTiming/PTAnalyzer/plugins/PTAnalyzer.cc

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

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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
  vector<float> tkOuterR;
  vector<float> tkOuterX;
  vector<float> tkOuterY;
  vector<float> tkOuterZ;

  vector<float> vtx_z;
  vector<int> vtx_nTks;
  vector<float> vtx1D_z;
  vector<int> vtx1D_nTks;
  vector<float> vtx4D_z;
  vector<int> vtx4D_nTks;
  vector<float> vtx4D_t;


  vector<float> pfEta;
  vector<float> pfPhi;
  vector<float> pfPt;
  vector<float> pfType;

};


class PTAnalyzer : public edm::EDAnalyzer  
{
public:
  explicit PTAnalyzer(const edm::ParameterSet&);
  ~PTAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
    
  void initEventStructure();

      
  //---inputs  
  EDGetTokenT<vector<PileupSummaryInfo> > PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  EDGetTokenT<View<reco::Vertex> > vertex1DToken_;
  EDGetTokenT<View<reco::Vertex> > vertex4DToken_;
  EDGetTokenT<View<reco::Track> > tracksToken_;
  EDGetTokenT<ValueMap<float> > trackTimeToken_;
  EDGetTokenT<View<reco::PFCandidate> > pfCandToken_;

  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree;
  eventInfo evInfo;

};

#ifndef _MUON_NEUTR_ISOLATION_H
#define _MUON_NEUTR_ISOLATION_H

// -*- C++ -*-
//
// Package:    PrecisionTiming/MuonNeutralIsolationAnalyzer
// Class:      MuonNeutralIsolationAnalyzer
// 
/**\class MuonNeutralIsolationAnalyzer MuonNeutralIsolationAnalyzer.cc PrecisionTiming/MuonNeutralIsolationAnalyzer/plugins/MuonNeutralIsolationAnalyzer.cc

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
//#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
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

#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
//#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
//#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
//#include "DataFormats/FTLRecHit/interface/FTLUncalibratedRecHit.h"
//#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
//#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
//#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
//#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
//#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
//#include "DataFormats/GeometrySurface/interface/Plane.h"
//#include "DataFormats/GeometrySurface/interface/Cylinder.h"

//#include "TrackingTools/PatternTools/interface/Trajectory.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
//#include "RecoMTD/TransientTrackingRecHit/interface/MTDTransientTrackingRecHitBuilder.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
//#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
//#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
//#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
//#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
/*
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
*/

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
//using namespace math;



struct eventInfo
{
  int npu;
  vector<float> neutrPfCand_pt;
  vector<float> neutrPfCand_eta;
  vector<float> neutrPfCand_phi;
  vector<float> neutrPfCand_tCluster;
  vector<float> neutrPfCand_tClusterSeed;
  vector<float> neutrPfCand_dRmu;
  vector<int> neutrPfCand_muIndex;
  vector<int> neutrPfCand_isMatchedToGenParticle;
  vector<int> neutrPfCand_isUnmatchedToGenParticle;

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

  vector<float> muon_neutrIso;
  vector<float> muon_neutrIso_dT2s;
  vector<float> muon_neutrIso_dT3s;
  vector<float> muon_neutrIso_dT5s;

};


//class MuonNeutralIsolationAnalyzer : public edm::EDAnalyzer  
class MuonNeutralIsolationAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit MuonNeutralIsolationAnalyzer(const edm::ParameterSet&);
  ~MuonNeutralIsolationAnalyzer();

  //  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
    
  void initEventStructure();


  //---inputs  
  const MTDGeometry* mtdGeometry_;

  EDGetTokenT<vector<PileupSummaryInfo> > PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken3D_;
  EDGetTokenT<View<reco::Vertex> > vertexToken4D_;
  EDGetTokenT<edm::View<reco::PFCandidate> >      pfcandToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  EDGetTokenT<vector<SimVertex> >  genVertexToken_;
  //EDGetTokenT<genXYZ> genXYZToken_;
  //EDGetTokenT<float>  genT0Token_;
  EDGetTokenT<View<reco::GenJet> > genJetsToken_;
  EDGetTokenT<View<reco::Muon> > muonsToken_; 
  EDGetTokenT<FTLClusterCollection> clustersBTLToken_;    
  EDGetTokenT<FTLClusterCollection> clustersETLToken_;    
  
  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree[10];
  eventInfo *evInfo[10];
  
  //--- options
  vector<double> timeResolutions_;
  bool mtd5sample_;
  double isoConeDR_;
  bool savePfCands_;
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
bool isPromptMuon(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles) ;
bool isMatchedToGenJet(const reco::Muon &muon, const edm::View<reco::GenJet>& genJet) ;
bool isFromTau(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles) ;
bool isMatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles) ;
bool isUnmatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles) ;
bool isMatchedToMTDCluster(const reco::PFCandidate &pfcand, FTLClusterCollection &clustersBTL, const MTDGeometry* mtdGeometry_, double &clusterTime) ;
*/


#endif

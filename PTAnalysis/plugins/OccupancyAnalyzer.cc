// -*- C++ -*-
//
// Package:    PrecisionTiming/PTAnalysis
// Class:      OccupancyAnalyzer
// 
/**\class OccupancyAnalyzer OccupancyAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/OccupancyAnalyzer.cc

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
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "PrecisionTiming/PTAnalysis/interface/OccupancyAnalyzer.h"

#include <TMath.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//
// constructors and destructor
//
OccupancyAnalyzer::OccupancyAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  //vertex1DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex1DTag" ) ) ),
  //vertex4DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex4DTag" ) ) ),
  tracksToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag>( "TracksTag" ) ) ),
  //trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  pfCandToken_( consumes<View<reco::PFCandidate> >( iConfig.getParameter<InputTag>( "PFCandidateTag" ) ) )
{
   //Now do what ever initialization is needed
  eventTree = fs_->make<TTree>( "event", "event" );
}


OccupancyAnalyzer::~OccupancyAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
OccupancyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  // -- get the 3D vertex collection
  Handle<View<reco::Vertex> > VertexCollectionH;
  iEvent.getByToken( vertexToken_, VertexCollectionH );
  const edm::View<reco::Vertex>& vertices = *VertexCollectionH;

  //// -- get the 3D vertex collection with same pT cut (0.7 GeV) as 4D
  //Handle<View<reco::Vertex> > Vertex1DCollectionH;
  //iEvent.getByToken( vertex1DToken_, Vertex1DCollectionH );
  //const edm::View<reco::Vertex>& vertices1D = *Vertex1DCollectionH;

  // -- get the 4D vertex collection
  //Handle<View<reco::Vertex> > Vertex4DCollectionH;
  //iEvent.getByToken( vertex4DToken_, Vertex4DCollectionH );
  //const edm::View<reco::Vertex>& vertices4D = *Vertex4DCollectionH;

  // -- get the PU 
  Handle<vector<PileupSummaryInfo> > PileupInfos;
  if( !iEvent.isRealData() ){
    iEvent.getByToken( PileUpToken_, PileupInfos );
   } else return;

  // -- get the track collection
  Handle<View<reco::Track> > TrackCollectionH;
  iEvent.getByToken(tracksToken_, TrackCollectionH);
  const edm::View<reco::Track>& tracks = *TrackCollectionH;

  //// -- get the trackTimeValueMap
  //Handle<ValueMap<float> > trackTimeValueMap;
  //iEvent.getByToken( trackTimeToken_, trackTimeValueMap );

  // -- get the PFCandidate collection
  Handle<View<reco::PFCandidate> > PFCandidateCollectionH;
  iEvent.getByToken(pfCandToken_, PFCandidateCollectionH);
  const edm::View<reco::PFCandidate>& pfCandidates = *PFCandidateCollectionH;

  

   // -- initialize output tree
  initEventStructure();
  
  
  // -- number of pileup events
  if( ! iEvent.isRealData() ) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI){
      Int_t pu_bunchcrossing = PVI->getBunchCrossing();
      if( pu_bunchcrossing == 0 ) {
	evInfo.npu = PVI->getPU_NumInteractions();
      }
    }
  }

  // -- number of reco vertices
  for(unsigned int ivtx=0; ivtx < vertices.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices[ivtx];
    evInfo.vtx_z.push_back(vtx.z());
    evInfo.vtx_nTks.push_back(vtx.tracksSize());
  }

  /*  for(unsigned int ivtx=0; ivtx < vertices1D.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices1D[ivtx];
    evInfo.vtx1D_z.push_back(vtx.z());
    evInfo.vtx1D_nTks.push_back(vtx.tracksSize());
  }

  for(unsigned int ivtx=0; ivtx < vertices4D.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices4D[ivtx];
    evInfo.vtx4D_z.push_back(vtx.z());
    evInfo.vtx4D_t.push_back(vtx.t());
    evInfo.vtx4D_nTks.push_back(vtx.tracksSize());
  }
  */

  // -- tracks
  for(unsigned int itk =0; itk < tracks.size(); itk++ ){
    const reco::Track& tk = tracks[itk];
    const auto tkRef = TrackCollectionH->refAt(itk);
    //cout << "Track pt, eta, phi, time : " << tk.pt() << " "<< tk.eta()<< " " << tk.phi() << " " << (*trackTimeValueMap)[tkRef]<<endl;  
  
    evInfo.tkPt.push_back(tk.pt());
    evInfo.tkEta.push_back(tk.eta());
    evInfo.tkPhi.push_back(tk.phi());
    evInfo.tkOuterR.push_back(tk.outerRadius());
    evInfo.tkOuterX.push_back(tk.outerX());
    evInfo.tkOuterY.push_back(tk.outerY());
    evInfo.tkOuterZ.push_back(tk.outerZ());
    //    evInfo.tkTime.push_back( (*trackTimeValueMap)[tkRef] );

  }


  // -- PF Candidates
  for(unsigned int ipf =0; ipf < pfCandidates.size(); ipf++ ){
    const reco::PFCandidate& cand = pfCandidates[ipf];
    //reco::PFCandidateRef pfRef = pfCandidates.refAt(ipf).castTo<reco::PFCandidateRef>(); 
    //reco::PFCandidate cand = *pfRef;
    //math::XYZPoint v(0., 0., 25.);
    //cand.setVertex(v);

    evInfo.pfPt.push_back(cand.pt());
    evInfo.pfEta.push_back(cand.eta());
    evInfo.pfPhi.push_back(cand.phi());
    evInfo.pfType.push_back(cand.particleId());

		   //    if (cand.pt()>0.7)    cout << cand.particleId()<< "  " << cand.pt()<< "   "<<cand.eta()<<"  " << cand.positionAtECALEntrance().x()<< "  " << cand.positionAtECALEntrance().z() <<endl;
  }
  
  // --- fill the tree
  eventTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
OccupancyAnalyzer::beginJob()
{

  eventTree->Branch( "npu",        &evInfo.npu);
  eventTree->Branch( "tkPt",       &evInfo.tkPt);
  eventTree->Branch( "tkEta",      &evInfo.tkEta);
  eventTree->Branch( "tkPhi",      &evInfo.tkPhi);
  //eventTree->Branch( "tkTime",     &evInfo.tkTime);
  eventTree->Branch( "tkOuterR",   &evInfo.tkOuterR);
  eventTree->Branch( "tkOuterX",   &evInfo.tkOuterX);
  eventTree->Branch( "tkOuterY",   &evInfo.tkOuterY);
  eventTree->Branch( "tkOuterZ",   &evInfo.tkOuterZ);
  eventTree->Branch( "vtx_z",      &evInfo.vtx_z);
  eventTree->Branch( "vtx_nTks",   &evInfo.vtx_nTks);
  //eventTree->Branch( "vtx1D_z",    &evInfo.vtx1D_z);
  //eventTree->Branch( "vtx1D_nTks", &evInfo.vtx1D_nTks);
  //eventTree->Branch( "vtx4D_z",    &evInfo.vtx4D_z);
  //eventTree->Branch( "vtx4D_nTks", &evInfo.vtx4D_nTks);
  //eventTree->Branch( "vtx4D_t",    &evInfo.vtx4D_t);
  eventTree->Branch( "pfPt",       &evInfo.pfPt);
  eventTree->Branch( "pfEta",      &evInfo.pfEta);
  eventTree->Branch( "pfPhi",      &evInfo.pfPhi);
  eventTree->Branch( "pfType",       &evInfo.pfType);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
OccupancyAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
OccupancyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
OccupancyAnalyzer::initEventStructure()
{
  // per-event tree:
  evInfo.npu = -1;
  evInfo.tkPt.clear();
  evInfo.tkEta.clear();
  evInfo.tkPhi.clear();
  // evInfo.tkTime.clear();
  evInfo.tkOuterR.clear();
  evInfo.tkOuterX.clear();
  evInfo.tkOuterY.clear();
  evInfo.tkOuterZ.clear();
  evInfo.vtx_z.clear();
  evInfo.vtx_nTks.clear();
  //evInfo.vtx1D_z.clear();
  //evInfo.vtx1D_nTks.clear();
  //evInfo.vtx4D_z.clear();
  //evInfo.vtx4D_nTks.clear();
  //evInfo.vtx4D_t.clear();
  evInfo.pfPt.clear();
  evInfo.pfEta.clear();
  evInfo.pfPhi.clear();
  evInfo.pfType.clear();


}

//define this as a plug-in
DEFINE_FWK_MODULE(OccupancyAnalyzer);

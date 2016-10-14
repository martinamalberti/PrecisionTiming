// -*- C++ -*-
//
// Package:    PTAnalyzer/PTAnalyzer
// Class:      PTAnalyzer
// 
/**\class PTAnalyzer PTAnalyzer.cc PTAnalysis/PTAnalysis/plugins/PTAnalyzer.cc

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

#include "PTAnalysis/PTAnalysis/interface/PTAnalyzer.h"


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
PTAnalyzer::PTAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  vertex4DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex4DTag" ) ) ),
  tracksToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag>( "TracksTag" ) ) ),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) )
{
   //Now do what ever initialization is needed
  eventTree = fs_->make<TTree>( "event", "event" );
}


PTAnalyzer::~PTAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  // -- get the 3D vertex collection
  Handle<View<reco::Vertex> > VertexCollectionH;
  iEvent.getByToken( vertexToken_, VertexCollectionH );
  const edm::View<reco::Vertex>& vertices = *VertexCollectionH;

  // -- get the 4D vertex collection
  Handle<View<reco::Vertex> > Vertex4DCollectionH;
  iEvent.getByToken( vertex4DToken_, Vertex4DCollectionH );
  const edm::View<reco::Vertex>& vertices4D = *Vertex4DCollectionH;

  // -- get the PU 
  Handle<vector<PileupSummaryInfo> > PileupInfos;
  if( !iEvent.isRealData() ){
    iEvent.getByToken( PileUpToken_, PileupInfos );
   } else return;

  // -- get the track collection
  //get track collections
  Handle<View<reco::Track> > TrackCollectionH;
  iEvent.getByToken(tracksToken_, TrackCollectionH);
  const edm::View<reco::Track>& tracks = *TrackCollectionH;


  // -- get the trackTimeValueMap
  Handle<ValueMap<float> > trackTimeValueMap;
  iEvent.getByToken( trackTimeToken_, trackTimeValueMap );

  /*
  // -- get the SimHits collection
  Handle<PSimHitContainer> ftHits;
  iEvent.getByToken(SimHitsToken_,ftHits);
  if (!ftHits.isValid()) {
    edm::LogWarning("PTAnalyzer  ") << "Error! can't get FastTimerHits product!";
    return ;
  }
  const edm::PSimHitContainer * ftHit = ftHits.product () ;
  */
  

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
  evInfo.nvtx   = vertices.size() ;
  for(unsigned int ivtx=0; ivtx < vertices.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices[ivtx];
    evInfo.vtx_z.push_back(vtx.z());
    evInfo.vtx_nTks.push_back(vtx.tracksSize());
  }

  evInfo.nvtx4D   = vertices4D.size() ;
  for(unsigned int ivtx=0; ivtx < vertices4D.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices4D[ivtx];
    evInfo.vtx4D_z.push_back(vtx.z());
    evInfo.vtx4D_t.push_back(vtx.t());
    evInfo.vtx4D_nTks.push_back(vtx.tracksSize());
  }

  // -- tracks
  evInfo.nTracks = tracks.size();
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
    evInfo.tkTime.push_back( (*trackTimeValueMap)[tkRef] );

  }


  // -- loop over simhits
  /* evInfo.nftHits = ftHit->size();
  
  for (vector<PSimHit>::const_iterator ftItr = ftHit->begin(); ftItr != ftHit->end(); ++ftItr) {
    FastTimeDetId id(ftItr->detUnitId());
    cout << " FT hit momentum    = " << ftItr->pabs() << endl;
    cout << " FT hit TOF         = " << ftItr->timeOfFlight() << endl;
    cout << " FT hit energy loss = " << ftItr->energyLoss() << endl;
    cout << " FT hit detUnitId() = " << ftItr->detUnitId()  << endl;
    cout << " FT hit id.ix, iy   = " << id.ix() << "," << id.iy()  << endl;
    cout << " FT isFastTime      = " << id.isFastTime()<<endl;
    cout << endl;
    
    evInfo.ftHitEnergyLoss.push_back(ftItr->energyLoss());
    evInfo.ftHitIx.push_back(id.ix());
    evInfo.ftHitIy.push_back(id.iy());
  }
  */
  
  // --- fill the tree
  eventTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
PTAnalyzer::beginJob()
{

  eventTree->Branch( "npu",     &evInfo.npu);
  eventTree->Branch( "nvtx",    &evInfo.nvtx);
  eventTree->Branch( "nvtx4D",  &evInfo.nvtx4D);
  eventTree->Branch( "nTracks", &evInfo.nTracks);
  eventTree->Branch( "tkPt",    &evInfo.tkPt);
  eventTree->Branch( "tkEta",   &evInfo.tkEta);
  eventTree->Branch( "tkPhi",   &evInfo.tkPhi);
  eventTree->Branch( "tkTime",  &evInfo.tkTime);
  eventTree->Branch( "tkOuterR",&evInfo.tkOuterR);
  eventTree->Branch( "tkOuterX",&evInfo.tkOuterX);
  eventTree->Branch( "tkOuterY",&evInfo.tkOuterY);
  eventTree->Branch( "tkOuterZ",&evInfo.tkOuterZ);
  eventTree->Branch( "vtx_z",   &evInfo.vtx_z);
  eventTree->Branch( "vtx_nTks",   &evInfo.vtx_nTks);
  eventTree->Branch( "vtx4D_z", &evInfo.vtx4D_z);
  eventTree->Branch( "vtx4D_nTks",   &evInfo.vtx4D_nTks);
  eventTree->Branch( "vtx4D_t", &evInfo.vtx4D_t);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PTAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
PTAnalyzer::initEventStructure()
{
  // per-event tree:
  evInfo.npu = -1;
  evInfo.nvtx = -1;
  evInfo.nvtx4D = -1;
  evInfo.nTracks = -1;
  evInfo.tkPt.clear();
  evInfo.tkEta.clear();
  evInfo.tkPhi.clear();
  evInfo.tkTime.clear();
  evInfo.tkOuterR.clear();
  evInfo.tkOuterX.clear();
  evInfo.tkOuterY.clear();
  evInfo.tkOuterZ.clear();
  evInfo.vtx_z.clear();
  evInfo.vtx_nTks.clear();
  evInfo.vtx4D_z.clear();
  evInfo.vtx4D_t.clear();
  evInfo.vtx4D_nTks.clear();


}

//define this as a plug-in
DEFINE_FWK_MODULE(PTAnalyzer);

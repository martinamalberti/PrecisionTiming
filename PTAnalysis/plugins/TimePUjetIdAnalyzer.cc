// -*- C++ -*-
//
// Package:    PrecisionTiming/PTAnalysis 
// Class:      TimePUJetIdAnalyzer
// 
/**\class TimePUJetIdAnalyzer TimePUJetIdAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/TimePUJetIdAnalyzer.cc

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

#include "DataFormats/Math/interface/deltaR.h"

#include "PrecisionTiming/PTAnalysis/interface/TimePUJetIdAnalyzer.h"


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
TimePUJetIdAnalyzer::TimePUJetIdAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  vertex4DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex4DTag" ) ) ),
  chsjetsToken_( consumes<View<reco::PFJet> >( iConfig.getParameter<InputTag>( "ChsJetsTag" ) ) ),
  genjetsToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag>( "GenJetsTag" ) ) ),
  muonsToken_( consumes<View<reco::Muon> >( iConfig.getParameter<InputTag>( "MuonsTag" ) ) ),
  genparticlesToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "GenParticlesTag" ) ) ),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) )
{
   //Now do what ever initialization is needed
  eventTree = fs_->make<TTree>( "event", "event" );
}


TimePUJetIdAnalyzer::~TimePUJetIdAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TimePUJetIdAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  // -- get the trackTimeValueMap
  Handle<ValueMap<float> > trackTimeValueMap;
  iEvent.getByToken( trackTimeToken_, trackTimeValueMap );

  // -- get the CHS jets collection
  Handle<View<reco::PFJet> > ChsJetsCollectionH;
  iEvent.getByToken(chsjetsToken_, ChsJetsCollectionH);
  const edm::View<reco::PFJet>& chsjets = *ChsJetsCollectionH;

  // -- get the GEN jets collection
  Handle<View<reco::GenJet> > GenJetsCollectionH;
  iEvent.getByToken(genjetsToken_, GenJetsCollectionH);
  const edm::View<reco::GenJet>& genjets = *GenJetsCollectionH;

  // -- get the genparticles
  Handle<View<reco::GenParticle> > GenParticlesCollectionH;
  iEvent.getByToken(genparticlesToken_, GenParticlesCollectionH);
  const edm::View<reco::GenParticle>& genparticles = *GenParticlesCollectionH;
  
  //   // -- initialize output tree
  ///initEventStructure();  
  
  // -- number of pileup events
  float pu = 0;
  if( ! iEvent.isRealData() ) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI){
      Int_t pu_bunchcrossing = PVI->getBunchCrossing();
      if( pu_bunchcrossing == 0 ) {
	pu = PVI->getPU_NumInteractions();
      }
    }
  }
  // -- vertex
  const reco::Vertex& vtx4D = vertices4D[0];
  

  const reco::Vertex& vtx = vertices[0];

  // -- gen vertex
  float genVtxZ = -999;
  for(unsigned int ip =0; ip < genparticles.size(); ip++ ){
    const reco::GenParticle& p = genparticles[ip];
    if (p.pdgId()==23) {
      genVtxZ = p.vertex().z();
      break;
    }
  }
  

  // -- CHS jets
  for(unsigned int ijet =0; ijet < chsjets.size(); ijet++ ){
    const reco::PFJet& jet = chsjets[ijet];
    
    // -- initialize output tree
    initEventStructure();

    evInfo.npu = pu;
    evInfo.vtx4D_z = vtx4D.z();
    evInfo.vtx4D_t = vtx4D.t();
    evInfo.vtx_z = vtx.z();
    evInfo.genvtx_z = genVtxZ;

    // --- check gen matching
    float mindr = 999999;
    bool genMatched = false;
    for (unsigned int ig =0; ig < genjets.size(); ig++ ){
      const reco::GenJet& genjet = genjets[ig];
      float dr = deltaR(jet, genjet);
      if (dr<0.4 && dr<mindr) {
        mindr = dr;
        genMatched = true;
	break;
      }
    }

    // --- remove overlap with muons 
    // ... to be added ...

    // --- loop over charged jet constituents
    reco::TrackRefVector tkRefs = jet.getTrackRefs() ;
    for (unsigned int ii = 0; ii < jet.getTrackRefs().size() ; ii ++) {
      reco::TrackRef tkRef = jet.getTrackRefs().at(ii);
      if ((*tkRef).pt() > 0.7 ){ 
	evInfo.chTime.push_back( (*trackTimeValueMap)[tkRef] );
	evInfo.chPt.push_back( tkRef->pt() );
	evInfo.chEta.push_back( tkRef->eta() );	
	evInfo.chPhi.push_back( tkRef->phi() );
      }
    }

    // -- tree variables
    evInfo.jetPt = jet.pt();
    evInfo.jetEta = jet.eta();
    evInfo.jetPhi = jet.phi();
    evInfo.jetIsMatchedToGen = genMatched;
    
    // --- fill the tree
    eventTree->Fill();
    
    
  }// end loop over jets
      
}


// ------------ method called once each job just before starting event loop  ------------
void 
TimePUJetIdAnalyzer::beginJob()
{

  eventTree->Branch( "npu",     &evInfo.npu);
  eventTree->Branch( "jetPt",   &evInfo.jetPt);
  eventTree->Branch( "jetEta",  &evInfo.jetEta);
  eventTree->Branch( "jetPhi",  &evInfo.jetPhi);
  eventTree->Branch( "jetIsMatchedToGen",   &evInfo.jetIsMatchedToGen);
  eventTree->Branch( "chTime",  &evInfo.chTime);
  eventTree->Branch( "chPt",    &evInfo.chPt);
  eventTree->Branch( "chEta",   &evInfo.chEta);
  eventTree->Branch( "chPhi",   &evInfo.chPhi);
  eventTree->Branch( "vtx4D_z", &evInfo.vtx4D_z);
  eventTree->Branch( "vtx4D_t", &evInfo.vtx4D_t);
  eventTree->Branch( "vtx_z",   &evInfo.vtx_z);
  eventTree->Branch( "genvtx_z",&evInfo.genvtx_z);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TimePUJetIdAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TimePUJetIdAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
TimePUJetIdAnalyzer::initEventStructure()
{
  // per-event tree:
  evInfo.npu = -1;
  evInfo.vtx4D_t   = -999;
  evInfo.vtx4D_z   = -999;
  evInfo.vtx_z   = -999;
  evInfo.genvtx_z   = -999;
  evInfo.jetPt = -999.;
  evInfo.jetEta = -999.;
  evInfo.jetPhi = -999.;
  evInfo.jetIsMatchedToGen = false;
  evInfo.chTime.clear();
  evInfo.chPt.clear();
  evInfo.chEta.clear();
  evInfo.chPhi.clear();




}

//define this as a plug-in
DEFINE_FWK_MODULE(TimePUJetIdAnalyzer);

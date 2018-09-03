// -*- C++ -*-
//
// Package:    PrecisionTiming/PhotonIsolationAnalyzer
// Class:      PhotonIsolationAnalyzer
// 
/**\class PhotonIsolationAnalyzer PhotonIsolationAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/PhotonIsolationAnalyzer.cc

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
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PrecisionTiming/PTAnalysis/interface/PhotonIsolationAnalyzer.h"

#include <TMath.h>
#include <TRandom.h>
#include "DataFormats/Math/interface/deltaR.h"

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
PhotonIsolationAnalyzer::PhotonIsolationAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  tracksToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag>( "TracksTag" ) ) ),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<InputTag>("genVtxTag"))),
  photonsToken_(consumes<View<reco::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonsTag")))
{
  timeResolutions_ = iConfig.getUntrackedParameter<vector<double> >("timeResolutions");
  isoConeDR_       = iConfig.getUntrackedParameter<vector<double> >("isoConeDR");
  BTLEfficiency_   = iConfig.getUntrackedParameter<double>("BTLEfficiency");
  ETLEfficiency_   = iConfig.getUntrackedParameter<double>("ETLEfficiency");
  saveTracks_      = iConfig.getUntrackedParameter<bool>("saveTracks");

   //Now do what ever initialization is needed
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    eventTree[iRes] = fs_->make<TTree>( Form("tree_%dps",int(timeResolutions_[iRes]*100) ), Form("tree_%dps", int(timeResolutions_[iRes]*100)) );
  }

}


PhotonIsolationAnalyzer::~PhotonIsolationAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  // -- get the 3D vertex collection
  Handle<View<reco::Vertex> > VertexCollectionH;
  iEvent.getByToken( vertexToken_, VertexCollectionH );
  const edm::View<reco::Vertex>& vertices = *VertexCollectionH;

  // -- get the PU 
  Handle<vector<PileupSummaryInfo> > PileupInfos;
  if( !iEvent.isRealData() ){
    iEvent.getByToken( PileUpToken_, PileupInfos );
   } else return;

  // -- get the photons
  Handle<View<reco::Photon> > PhotonCollectionH;
  iEvent.getByToken(photonsToken_, PhotonCollectionH);
  const edm::View<reco::Photon>& photons = *PhotonCollectionH;

  // -- get the track collection
  Handle<View<reco::Track> > TrackCollectionH;
  iEvent.getByToken(tracksToken_, TrackCollectionH);
  const edm::View<reco::Track>& tracks = *TrackCollectionH;

  // -- get the trackTimeValueMap
  Handle<ValueMap<float> > trackTimeValueMap;
  iEvent.getByToken( trackTimeToken_, trackTimeValueMap );

  // -- get the gen particles collection
  Handle<View<reco::GenParticle> > GenParticleCollectionH;
  iEvent.getByToken(genPartToken_, GenParticleCollectionH);
  const edm::View<reco::GenParticle>& genParticles = *GenParticleCollectionH;

  // -- get the gen vertex collection
  Handle<vector<SimVertex> > GenVertexCollectionH;
  iEvent.getByToken( genVertexToken_, GenVertexCollectionH );
  const vector<SimVertex>& genVertices = *GenVertexCollectionH;

  
  // -- initialize output tree
  initEventStructure();
  
  
  // -- number of pileup events
  int nPU = 0;
  if( ! iEvent.isRealData() ) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI){
      Int_t pu_bunchcrossing = PVI->getBunchCrossing();
      if( pu_bunchcrossing == 0 ) {
	nPU = PVI->getPU_NumInteractions();
      }
    }
  }

  //---get truth PV
  SimVertex genPV = genVertices.at(0);
  double mindz = 999999.;
  int pv_index = 0;

  TRandom *gRandom = new TRandom();


  // -- find the reco vertex closest to the gen vertex
  for(unsigned int ivtx=0; ivtx < vertices.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices[ivtx];
    const float dz = std::abs(vtx.z() - genPV.position().z());
    if( dz < mindz )
      {
	mindz = dz;
	pv_index = ivtx;
      }
  }

  // -- get isolation around a candidate photon 
  // --- using only vtx closest to gen vtx
  const reco::Vertex& vtx = vertices[pv_index];

  for(unsigned int ipho=0; ipho < photons.size(); ipho++ ){

    const reco::Photon& photon = photons[ipho];

    // -- minimal checks
    if(photon.pt() < 15.) continue;
    
    // -- real photons or fake photons
    


    // -- compute charged isolations
    float chIso[4] = {};
    float chIso_dT[4][3] = {{}};
    
    for(reco::Vertex::trackRef_iterator ti = vtx.tracks_begin(); ti!=vtx.tracks_end(); ++ti){

      // AAAAAAAAAAAAAAAA !!!!!! need to recompute photon eta, phi wrt to the PV
      float dr = deltaR(photon.eta(), photon.phi(), (*ti) -> eta(), (*ti) -> phi());

      // get track ref to generalTracks and track time
      reco::TrackRef vtxtkRef = ti->castTo<reco::TrackRef>();
      float trkTime = (*trackTimeValueMap)[vtxtkRef];
      //cout << "Track pt, eta, phi, time : " << (*ti)->pt() << " "<< (*ti)->eta()<< " " << (*ti)->phi() << "  "<<(*trackTimeValueMap)[vtxtkRef]<<endl; 


      // -- loop over cone sizes
      for (unsigned int iCone = 0 ; iCone < isoConeDR_.size(); iCone++){
      
	if (dr < isoConeDR_[iCone]){

	  chIso[iCone]+= (*ti) -> pt();

	  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
	    float time_resol = timeResolutions_[iRes]; 
	    float extra_resol = sqrt(time_resol*time_resol - 0.030*0.030);
	    trkTime+=gRandom->Gaus(0., extra_resol);
	    float dt = std::abs(trkTime - vtx.t());
	    if ( dt < 3. * time_resol){
	      chIso_dT[iCone][iRes]+= (*ti) -> pt();
	    }
	  }// end loop over time resolutions

	}
 
      }// end loop over cone sizes

      if (saveTracks_){
	for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
	  evInfo[iRes].tkTime.push_back(trkTime);
	  evInfo[iRes].tkPt.push_back((*ti)->pt());
	  evInfo[iRes].tkEta.push_back((*ti)->eta());
	  evInfo[iRes].tkPhi.push_back((*ti)->phi());
	}
      }

    }// end loop over PV tracks
    
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes].npu = nPU;
      evInfo[iRes].vtx_z = vtx.z();
      evInfo[iRes].vtx_t = vtx.t();
      evInfo[iRes].vtx_nTracks = vtx.nTracks();
      evInfo[iRes].photon_pt.push_back(photon.pt());  
      evInfo[iRes].photon_eta.push_back(photon.eta());  
      evInfo[iRes].photon_phi.push_back(photon.phi());  
      for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
	evInfo[iRes].photon_chIso[iCone].push_back(chIso[iCone]);  
	evInfo[iRes].photon_chIso_dT[iCone][iRes].push_back(chIso_dT[iCone][iRes]);
      }  
    
    }

  }// end loop over photons

  // --- fill the tree
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++)
    eventTree[iRes]->Fill();
  
}// end analyze event
  

// ------------ method called once each job just before starting event loop  ------------
void 
PhotonIsolationAnalyzer::beginJob()
{
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
  
    eventTree[iRes]->Branch( "npu",               &evInfo[iRes].npu);
    eventTree[iRes]->Branch( "vtx_z",             &evInfo[iRes].vtx_z);
    eventTree[iRes]->Branch( "vtx_nTracks",       &evInfo[iRes].vtx_nTracks);
    eventTree[iRes]->Branch( "vtx_nTracks_dT",    &evInfo[iRes].vtx_nTracks_dT);
    eventTree[iRes]->Branch( "photon_pt",         &evInfo[iRes].photon_pt);  
    eventTree[iRes]->Branch( "photon_eta",        &evInfo[iRes].photon_eta);  
    eventTree[iRes]->Branch( "photon_phi",        &evInfo[iRes].photon_phi);  

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      eventTree[iRes]->Branch( Form("photon_chIso%.2d",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].photon_chIso[iCone]);  
      eventTree[iRes]->Branch( Form("photon_chIso%.2d_dT",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].photon_chIso[iCone][iRes]);  
    }

    if (saveTracks_){
      eventTree[iRes]->Branch( "tkTime",   &evInfo[iRes].tkTime);
      eventTree[iRes]->Branch( "tkPt",     &evInfo[iRes].tkPt);
      eventTree[iRes]->Branch( "tkEta",    &evInfo[iRes].tkEta);
      eventTree[iRes]->Branch( "tkPhi",    &evInfo[iRes].tkPhi);
    }
  
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonIsolationAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
PhotonIsolationAnalyzer::initEventStructure()
{
  // per-event trees:
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    evInfo[iRes].npu = -1;
    evInfo[iRes].vtx_z = -999;
    evInfo[iRes].vtx_nTracks = -999.;
    evInfo[iRes].vtx_nTracks_dT = -999.;

    evInfo[iRes].photon_pt.clear();
    evInfo[iRes].photon_eta.clear();
    evInfo[iRes].photon_phi.clear();

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      evInfo[iRes].photon_chIso[iCone].clear();
      evInfo[iRes].photon_chIso_dT[iCone][iRes].clear();
    }

    if (saveTracks_){
      evInfo[iRes].tkTime.clear();
      evInfo[iRes].tkPt.clear();
      evInfo[iRes].tkEta.clear();
      evInfo[iRes].tkPhi.clear();
    }
  }
}



//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIsolationAnalyzer);

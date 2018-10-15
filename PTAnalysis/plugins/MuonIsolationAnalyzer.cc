// -*- C++ -*-%
//
// Package:    PrecisionTiming/MuonIsolationAnalyzer
// Class:      MuonIsolationAnalyzer
// 
/**\class MuonIsolationAnalyzer MuonIsolationAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/MuonIsolationAnalyzer.cc

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

#include "PrecisionTiming/PTAnalysis/interface/MuonIsolationAnalyzer.h"

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
MuonIsolationAnalyzer::MuonIsolationAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken3D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag3D" ) ) ),
  vertexToken4D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag4D" ) ) ),
  pfcandToken_( consumes<View<reco::PFCandidate> >( iConfig.getParameter<InputTag>( "PFCandidateTag" ) ) ),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<InputTag>("genVtxTag"))),
  genJetsToken_(consumes<View<reco::GenJet> >(iConfig.getUntrackedParameter<InputTag>("genJetsTag"))),
  muonsToken_(consumes<View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonsTag")))
{
  timeResolutions_ = iConfig.getUntrackedParameter<vector<double> >("timeResolutions");
  isoConeDR_       = iConfig.getUntrackedParameter<vector<double> >("isoConeDR");
  saveTracks_      = iConfig.getUntrackedParameter<bool>("saveTracks");
  maxDz_           = iConfig.getUntrackedParameter<double>("maxDz");
  minDr_           = iConfig.getUntrackedParameter<double>("minDr");

  //Now do what ever initialization is needed
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    eventTree[iRes] = fs_->make<TTree>( Form("tree_%dps",int(timeResolutions_[iRes]*1000) ), Form("tree_%dps", int(timeResolutions_[iRes]*1000)) );
    cout << iRes << "  " << timeResolutions_[iRes] << "  " << eventTree[iRes]->GetName() <<endl;
  }

}


MuonIsolationAnalyzer::~MuonIsolationAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
  // -- get the vertex 3D collection
  Handle<View<reco::Vertex> > Vertex3DCollectionH;
  iEvent.getByToken( vertexToken3D_, Vertex3DCollectionH );
  const edm::View<reco::Vertex>& vertices3D = *Vertex3DCollectionH;

  // -- get the vertex 4D collection
  Handle<View<reco::Vertex> > Vertex4DCollectionH;
  iEvent.getByToken( vertexToken4D_, Vertex4DCollectionH );
  const edm::View<reco::Vertex>& vertices4D = *Vertex4DCollectionH;

  // -- get the PU 
  Handle<vector<PileupSummaryInfo> > PileupInfos;
  if( !iEvent.isRealData() ){
    iEvent.getByToken( PileUpToken_, PileupInfos );
   } else return;

  // -- get the muons
  Handle<View<reco::Muon> > MuonCollectionH;
  iEvent.getByToken(muonsToken_, MuonCollectionH);
  const edm::View<reco::Muon>& muons = *MuonCollectionH;

  // -- get the PFCandidate collection
  Handle<View<reco::PFCandidate> > PFCandidateCollectionH;
  iEvent.getByToken(pfcandToken_, PFCandidateCollectionH);
  const edm::View<reco::PFCandidate>& pfcands = *PFCandidateCollectionH;

  // -- get the gen particles collection
  Handle<View<reco::GenParticle> > GenParticleCollectionH;
  iEvent.getByToken(genPartToken_, GenParticleCollectionH);
  const edm::View<reco::GenParticle>& genParticles = *GenParticleCollectionH;

  // -- get the gen vertex collection
  Handle<vector<SimVertex> > GenVertexCollectionH;
  iEvent.getByToken( genVertexToken_, GenVertexCollectionH );
  const vector<SimVertex>& genVertices = *GenVertexCollectionH;

  // -- get the gen jets collection
  Handle<View<reco::GenJet> > GenJetCollectionH;
  iEvent.getByToken(genJetsToken_, GenJetCollectionH);
  const edm::View<reco::GenJet>& genJets = *GenJetCollectionH;
  
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
  int pv_index_3D = 0;
  int pv_index_4D = 0;

  TRandom *gRandom = new TRandom();

  // -- find the reco vertex closest to the gen vertex (3D)
  for(unsigned int ivtx=0; ivtx < vertices3D.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices3D[ivtx];
    const float dz = std::abs(vtx.z() - genPV.position().z());
    if( dz < mindz )
      {
	mindz = dz;
	pv_index_3D = ivtx;
      }
  }
  
  // -- find the reco vertex closest to the gen vertex (4D)
  mindz = 999999.;
  for(unsigned int ivtx=0; ivtx < vertices4D.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices4D[ivtx];
    const float dz = std::abs(vtx.z() - genPV.position().z());
    if( dz < mindz )
      {
	mindz = dz;
	pv_index_4D = ivtx;
      }
  }



  // -- get isolation around a candidate photon 
  // --- using only vtx closest to gen vtx
  const reco::Vertex& vtx   = vertices4D[pv_index_4D];
  const reco::Vertex& vtx3D = vertices3D[pv_index_3D];


  // --- start loop over muons
  for(unsigned int imu=0; imu < muons.size(); imu++ ){

    const reco::Muon& muon = muons[imu];

    // -- minimal checks
    if (muon.track().isNull()) continue;
    if (muon.pt() < 10.) continue;
  
    // -- check if prompt or fake muon
    bool isPromptMu = isPromptMuon(muon,genParticles);
    bool isMatchedGenJet  = isMatchedToGenJet(muon, genJets);


    // -- compute charged isolations
    const int nCones = isoConeDR_.size();
    const int nResol = timeResolutions_.size();
    float chIso[nCones] = {0.};
    float chIso_dT[nCones][nResol] = {{0.}} ;
    float time[nCones][nResol] = {{0.}}; 
    
    // -- loop over charged pf candidates
    for(unsigned icand = 0; icand < pfcands.size(); ++icand) {
      const reco::PFCandidate& pfcand = pfcands[icand];
      if (pfcand.charge() == 0 ) continue;
      
      // -- skip the track if it is the muon track
      auto pfcandRef = pfcands.refAt(icand);
      reco::TrackRef trackRef = pfcandRef->trackRef();                   
      if ( trackRef.isNull() ) continue;
      if ( trackRef == muon.track() ) continue;
      

      
      // -- skip tracks if dz(track, vertex) > dzCut 
      float dz4D = std::abs( trackRef->dz(vtx.position()) );
      float dz3D = std::abs( trackRef->dz(vtx3D.position()) );
      //std::cout << "*** 3D ****  dz = " << dz3D << "  dz = " <<  std::abs(trackRef->dz( vtx3D.position()  )) << "   diff = " << dz3D-std::abs(trackRef->dz( vtx3D.position()  ))<<std::endl;;
      //std::cout << "*** 4D ****  dz = " << dz4D << "  dz = " <<  std::abs(trackRef->dz( vtx.position()  )) << "   diff = " << dz4D-std::abs(trackRef->dz( vtx.position()  ))<<std::endl;;

      float dxy4D = std::abs( trackRef->dxy(vtx.position()) );
      float dxy3D = std::abs( trackRef->dxy(vtx3D.position()) );
     
      float dr  = deltaR(muon.eta(), muon.phi(), pfcand.eta(), pfcand.phi());
      
      //std::cout << pfcand.eta() << "  " << trackRef->eta() << "  " << pfcand.phi() << "  " << trackRef->phi() <<std::endl;

      // --- no timing 
      if (dz3D < maxDz_  && dxy3D < 0.02){
	for (unsigned int iCone = 0 ; iCone < isoConeDR_.size(); iCone++){
	  if (dr > minDr_ && dr < isoConeDR_[iCone]){
	    chIso[iCone]+= pfcand.pt();
	  }
	}
      }


      // --- with timing
      if ( dz4D < maxDz_  && dxy4D < 0.02 ){
      
	double dt = 0;
	for (unsigned int iCone = 0 ; iCone < isoConeDR_.size(); iCone++){
	  if (dr > minDr_ && dr < isoConeDR_[iCone]){ 	  
	    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){                                                                                                    
	      double time_resol = timeResolutions_[iRes];
	      double extra_resol = sqrt(time_resol*time_resol - 0.03*0.03);                                                                                                      
	      if ( pfcand.isTimeValid() ) {
		time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
		dt = std::abs(time[iCone][iRes] - vtx.t());
	      }
	      else{
		dt = 0;
	      }
	      if ( dt < 3*time_resol){                                                                                                                                         
		chIso_dT[iCone][iRes]+= pfcand.pt();                                                                                                                            
	      }                                                                                                                                                                   
	    }// end loop over time resolutions                                                                                                                                    
	
	    // save info for tracks in the isolation cone
	    if ( isoConeDR_[iCone] == 0.3 && saveTracks_){
	      for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){                                                                                                      
		evInfo[iRes].track_t.push_back(time[iCone][iRes]);                                                                                                                               
		evInfo[iRes].track_dz.push_back(trackRef->dz( vtx.position() )); 
		evInfo[iRes].track_dz3D.push_back(trackRef->dz( vtx3D.position() )); 
		evInfo[iRes].track_pt.push_back(pfcand.pt());
		evInfo[iRes].track_eta.push_back(pfcand.eta());                                                                                                                           
		evInfo[iRes].track_phi.push_back(pfcand.phi());                                                                                                                           
	      }
	    }
	  }                                                                                                                                                                      
	}// end loop over cone sizes                                                                                                                                              
      }
    }// end loop over tracks
        
    // -- fill muon info
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes].muon_pt.push_back(muon.pt());  
      evInfo[iRes].muon_eta.push_back(muon.eta());  
      evInfo[iRes].muon_phi.push_back(muon.phi()); 
      evInfo[iRes].muon_isLoose.push_back(muon::isLooseMuon(muon)); 
      evInfo[iRes].muon_isMedium.push_back(muon::isMediumMuon(muon));
      int isTight3D = muon::isTightMuon(muon, vtx3D);
      int isTight4D = muon::isTightMuon(muon, vtx);
      evInfo[iRes].muon_isTight3D.push_back(isTight3D); 
      evInfo[iRes].muon_isTight4D.push_back(isTight4D); 
      evInfo[iRes].muon_isPrompt.push_back(isPromptMu);
      evInfo[iRes].muon_isMatchedToGenJet.push_back(isMatchedGenJet);
        
      for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
	(evInfo[iRes].muon_chIso[iCone]).push_back(chIso[iCone]);  
	(evInfo[iRes].muon_chIso_dT[iCone][iRes]).push_back(chIso_dT[iCone][iRes]);
      }
    }
  }// end loop over muon



  // -- fill event info
  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
    evInfo[iRes].npu = nPU;
    evInfo[iRes].vtx_t = vtx.t();
    evInfo[iRes].vtx_z = vtx.z();
    evInfo[iRes].vtx3D_z = vtx3D.z();
    evInfo[iRes].vtxGen_z = genPV.position().z();
    evInfo[iRes].vtxGen_t = genPV.position().t();
  }
  
  
  // --- fill the tree
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
    eventTree[iRes]->Fill();
  }

  
}// end analyze event
  

// ------------ method called once each job just before starting event loop  ------------
void 
MuonIsolationAnalyzer::beginJob()
{
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
  
    eventTree[iRes]->Branch( "npu",             &evInfo[iRes].npu);
    eventTree[iRes]->Branch( "vtxGen_z",        &evInfo[iRes].vtxGen_z);
    eventTree[iRes]->Branch( "vtxGen_t",        &evInfo[iRes].vtxGen_t);
    eventTree[iRes]->Branch( "vtx3D_z",         &evInfo[iRes].vtx3D_z);
    eventTree[iRes]->Branch( "vtx_z",           &evInfo[iRes].vtx_z);
    eventTree[iRes]->Branch( "vtx_t",           &evInfo[iRes].vtx_t);
    eventTree[iRes]->Branch( "muon_pt",         &evInfo[iRes].muon_pt);  
    eventTree[iRes]->Branch( "muon_eta",        &evInfo[iRes].muon_eta);  
    eventTree[iRes]->Branch( "muon_phi",        &evInfo[iRes].muon_phi);  
    eventTree[iRes]->Branch( "muon_isLoose",    &evInfo[iRes].muon_isLoose);  
    eventTree[iRes]->Branch( "muon_isMedium",   &evInfo[iRes].muon_isMedium);  
    eventTree[iRes]->Branch( "muon_isTight3D",  &evInfo[iRes].muon_isTight3D);  
    eventTree[iRes]->Branch( "muon_isTight4D",  &evInfo[iRes].muon_isTight4D);  
    eventTree[iRes]->Branch( "muon_isPrompt",   &evInfo[iRes].muon_isPrompt);  
    eventTree[iRes]->Branch( "muon_isMatchedToGenJet",   &evInfo[iRes].muon_isMatchedToGenJet);  

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      eventTree[iRes]->Branch( Form("muon_chIso%.2d",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].muon_chIso[iCone]);  
      eventTree[iRes]->Branch( Form("muon_chIso%.2d_dT",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].muon_chIso_dT[iCone][iRes]);  
    }

    if (saveTracks_){
      eventTree[iRes]->Branch( "track_t",      &evInfo[iRes].track_t);
      eventTree[iRes]->Branch( "track_dz",     &evInfo[iRes].track_dz);
      eventTree[iRes]->Branch( "track_dz3D",   &evInfo[iRes].track_dz3D);
      eventTree[iRes]->Branch( "track_pt",     &evInfo[iRes].track_pt);
      eventTree[iRes]->Branch( "track_eta",    &evInfo[iRes].track_eta);
      eventTree[iRes]->Branch( "track_phi",    &evInfo[iRes].track_phi);
    }
  
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonIsolationAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
MuonIsolationAnalyzer::initEventStructure()
{
  // per-event trees:
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    evInfo[iRes].npu = -1;
    evInfo[iRes].vtxGen_t = -999;
    evInfo[iRes].vtxGen_z = -999;
    evInfo[iRes].vtx3D_z = -999;
    evInfo[iRes].vtx_z = -999;
    evInfo[iRes].vtx_t = -999;

    evInfo[iRes].muon_pt.clear();
    evInfo[iRes].muon_eta.clear();
    evInfo[iRes].muon_phi.clear();
    evInfo[iRes].muon_isLoose.clear();
    evInfo[iRes].muon_isMedium.clear();
    evInfo[iRes].muon_isTight3D.clear();    
    evInfo[iRes].muon_isTight4D.clear();
    evInfo[iRes].muon_isPrompt.clear();
    evInfo[iRes].muon_isMatchedToGenJet.clear();

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      evInfo[iRes].muon_chIso[iCone].clear();
      evInfo[iRes].muon_chIso_dT[iCone][iRes].clear();
    }

    if (saveTracks_){
      evInfo[iRes].track_t.clear();
      evInfo[iRes].track_dz.clear();
      evInfo[iRes].track_dz3D.clear();
      evInfo[iRes].track_pt.clear();
      evInfo[iRes].track_eta.clear();
      evInfo[iRes].track_phi.clear();
    }
  }
}



// --- matching to gen muon
bool isPromptMuon(const reco::Muon& muon, const edm::View<reco::GenParticle>& genParticles)
{
  bool isPrompt = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 13) continue;
    if ( !genp.isPromptFinalState() ) continue;
    if ( genp.pt() < 10.0 ) continue;
    double dr = deltaR(muon,genp);
    if (dr > 0.1){ 
      continue;
    }
    else{ 
      isPrompt=true;
      break;
    }
  }
  
  return isPrompt;
}


// --- matching to gen jet
bool isMatchedToGenJet(const reco::Muon& muon, const edm::View<reco::GenJet>& genJets)
{
  bool isMatched = false;

  for(unsigned int ip=0; ip < genJets.size(); ip++ ){
    const reco::GenJet& genj = genJets[ip];
    if ( genj.pt() < 15.0  || genj.hadEnergy()/genj.energy() < 0.3) continue;
    double dr = deltaR(muon,genj);
    if (dr > 0.3){ 
      continue;
    }
    else{ 
      isMatched=true;
      break;
    }
  }
  
  return isMatched;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationAnalyzer);

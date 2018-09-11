// -*- C++ -*-%
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
  vertexToken3D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag3D" ) ) ),
  vertexToken4D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag4D" ) ) ),
  pfcandToken_( consumes<View<reco::PFCandidate> >( iConfig.getParameter<InputTag>( "PFCandidateTag" ) ) ),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<InputTag>("genVtxTag"))),
  barrelPhotonsToken_(consumes<View<reco::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("barrelPhotonsTag"))),
  endcapPhotonsToken_(consumes<View<reco::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapPhotonsTag")))
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

  // -- get the barrel photons
  Handle<View<reco::Photon> > BarrelPhotonCollectionH;
  iEvent.getByToken(barrelPhotonsToken_, BarrelPhotonCollectionH);
  const edm::View<reco::Photon>& barrelPhotons = *BarrelPhotonCollectionH;

  // -- get the endcap photons
  Handle<View<reco::Photon> > EndcapPhotonCollectionH;
  iEvent.getByToken(endcapPhotonsToken_, EndcapPhotonCollectionH);
  const edm::View<reco::Photon>& endcapPhotons = *EndcapPhotonCollectionH;

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


  // --- start loop over barrel 
  for(unsigned int ipho=0; ipho < barrelPhotons.size(); ipho++ ){

    const reco::Photon& photon = barrelPhotons[ipho];

    // -- minimal checks
    if(photon.pt() < 15.) continue;
  
    // -- check if prompt or fake photon
    bool isPromptPho = isPromptPhoton(photon,genParticles);


    // -- recompute p4 wrt to the chosen vertex (taken from flashgg)
    math::XYZTLorentzVector pho_p4 = correctP4(photon, vtx);
    math::XYZTLorentzVector pho_p4_3Dvtx = correctP4(photon, vtx3D);

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

      // -- skip tracks if dz(track, vertex) > dzCut 
      float dz = std::abs(pfcand.vz() - vtx.z());
      if (dz  > maxDz_) continue; // 1 mm

      float dxy = sqrt( pow(pfcand.vx() - vtx.x(),2) + pow(pfcand.vy() - vtx.y(), 2)); 
      if (dxy > 0.02) continue;
     

      // --- no timing 
      float dr  = deltaR(pho_p4_3Dvtx.eta(), pho_p4_3Dvtx.phi(), pfcand.eta(), pfcand.phi());
      for (unsigned int iCone = 0 ; iCone < isoConeDR_.size(); iCone++){
        if (dr > minDr_ && dr < isoConeDR_[iCone]){
	  chIso[iCone]+= pfcand.pt();
	}
      }


      // --- with timing
      dr  = deltaR(pho_p4.eta(), pho_p4.phi(), pfcand.eta(), pfcand.phi());     
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
	      evInfo[iRes].track_dz.push_back(dz);                                                                                                                               
	      evInfo[iRes].track_pt.push_back(pfcand.pt());
	      evInfo[iRes].track_eta.push_back(pfcand.eta());                                                                                                                           
	      evInfo[iRes].track_phi.push_back(pfcand.phi());                                                                                                                           
	    }                                                                                                                                                                       
	  }
	}                                                                                                                                                                      
      }// end loop over cone sizes                                                                                                                                              
    }// end loop over tracks

    
    // -- fill photon info
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes].photon_pt.push_back(pho_p4.pt());  
      evInfo[iRes].photon_eta.push_back(pho_p4.eta());  
      evInfo[iRes].photon_phi.push_back(pho_p4.phi()); 
      evInfo[iRes].photon_isPrompt.push_back(isPromptPho);
      evInfo[iRes].photon_hasConversionTracks.push_back(photon.hasConversionTracks());
      evInfo[iRes].photon_sigmaIetaIeta.push_back(photon.sigmaIetaIeta());
      evInfo[iRes].photon_r9.push_back(photon.r9());
      
      for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
	(evInfo[iRes].photon_chIso[iCone]).push_back(chIso[iCone]);  
	(evInfo[iRes].photon_chIso_dT[iCone][iRes]).push_back(chIso_dT[iCone][iRes]);
      }
    }
  }// end loop over barrel photons



  // --- start loop over endcap
  for(unsigned int ipho=0; ipho < endcapPhotons.size(); ipho++ ){
    
    const reco::Photon& photon = endcapPhotons[ipho];

    // -- minimal checks
    if (photon.pt()  < 15.) continue;
    if (fabs(photon.eta()) < 1.5) continue; // only endcap photons

    // -- check if prompt or fake photon
    bool isPromptPho = isPromptPhoton(photon,genParticles);

    // -- recompute p4 wrt to the chosen vertex (taken from flashgg)
    math::XYZTLorentzVector pho_p4 = correctP4(photon, vtx);
    math::XYZTLorentzVector pho_p4_3Dvtx = correctP4(photon, vtx3D);

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

      // -- skip tracks if dz(track, vertex) > dzCut
      float dz = std::abs(pfcand.vz() - vtx.z());
      if (dz  > maxDz_) continue; // 1 mm

      float dxy = sqrt( pow(pfcand.vx() - vtx.x(),2) + pow(pfcand.vy() - vtx.y(), 2));
      if (dxy > 0.02) continue;



      // --- no timing
      float dr  = deltaR(pho_p4_3Dvtx.eta(), pho_p4_3Dvtx.phi(), pfcand.eta(), pfcand.phi());
      for (unsigned int iCone = 0 ; iCone < isoConeDR_.size(); iCone++){
        if (dr > minDr_ && dr < isoConeDR_[iCone]){
          chIso[iCone]+= pfcand.pt();
        }
      }



      // --- with timing
      dr  = deltaR(pho_p4.eta(), pho_p4.phi(), pfcand.eta(), pfcand.phi());
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
              evInfo[iRes].track_dz.push_back(dz);
              evInfo[iRes].track_pt.push_back(pfcand.pt());
              evInfo[iRes].track_eta.push_back(pfcand.eta());
              evInfo[iRes].track_phi.push_back(pfcand.phi());
            }
          }
        }
      }// end loop over cone sizes
    }// end loop over tracks                                                                                                                                                                                                                
   
    // -- fill photon info
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes].photon_pt.push_back(pho_p4.pt());
      evInfo[iRes].photon_eta.push_back(pho_p4.eta());
      evInfo[iRes].photon_phi.push_back(pho_p4.phi());
      evInfo[iRes].photon_isPrompt.push_back(isPromptPho);
      evInfo[iRes].photon_hasConversionTracks.push_back(photon.hasConversionTracks());
      evInfo[iRes].photon_sigmaIetaIeta.push_back(photon.sigmaIetaIeta());
      evInfo[iRes].photon_r9.push_back(photon.r9());
      
      for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
        (evInfo[iRes].photon_chIso[iCone]).push_back(chIso[iCone]);
        (evInfo[iRes].photon_chIso_dT[iCone][iRes]).push_back(chIso_dT[iCone][iRes]);
      }
    }
  }// end loop over endcap photons


  // -- fill event info
  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
    evInfo[iRes].npu = nPU;
    evInfo[iRes].vtx_t = vtx.t();
    evInfo[iRes].vtx_z = vtx.z();
    evInfo[iRes].vtx3D_z = vtx3D.z();
    evInfo[iRes].vtxGen_z = genPV.position().z();
  }
  
  
  // --- fill the tree
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
    eventTree[iRes]->Fill();
  }

  
}// end analyze event
  

// ------------ method called once each job just before starting event loop  ------------
void 
PhotonIsolationAnalyzer::beginJob()
{
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
  
    eventTree[iRes]->Branch( "npu",               &evInfo[iRes].npu);
    eventTree[iRes]->Branch( "vtxGen_z",          &evInfo[iRes].vtxGen_z);
    eventTree[iRes]->Branch( "vtx3D_z",           &evInfo[iRes].vtx3D_z);
    eventTree[iRes]->Branch( "vtx_z",             &evInfo[iRes].vtx_z);
    eventTree[iRes]->Branch( "vtx_t",             &evInfo[iRes].vtx_t);
    eventTree[iRes]->Branch( "photon_pt",         &evInfo[iRes].photon_pt);  
    eventTree[iRes]->Branch( "photon_eta",        &evInfo[iRes].photon_eta);  
    eventTree[iRes]->Branch( "photon_phi",        &evInfo[iRes].photon_phi);  
    eventTree[iRes]->Branch( "photon_isPrompt",   &evInfo[iRes].photon_isPrompt);  
    eventTree[iRes]->Branch( "photon_hasConversionTracks",   &evInfo[iRes].photon_hasConversionTracks);  
    eventTree[iRes]->Branch( "photon_sigmaIetaIeta",   &evInfo[iRes].photon_sigmaIetaIeta);  
    eventTree[iRes]->Branch( "photon_r9",   &evInfo[iRes].photon_r9);  

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      eventTree[iRes]->Branch( Form("photon_chIso%.2d",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].photon_chIso[iCone]);  
      eventTree[iRes]->Branch( Form("photon_chIso%.2d_dT",int(isoConeDR_[iCone]*10) ), &evInfo[iRes].photon_chIso_dT[iCone][iRes]);  
    }

    if (saveTracks_){
      eventTree[iRes]->Branch( "track_t",   &evInfo[iRes].track_t);
      eventTree[iRes]->Branch( "track_dz",   &evInfo[iRes].track_dz);
      eventTree[iRes]->Branch( "track_pt",     &evInfo[iRes].track_pt);
      eventTree[iRes]->Branch( "track_eta",    &evInfo[iRes].track_eta);
      eventTree[iRes]->Branch( "track_phi",    &evInfo[iRes].track_phi);
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
    evInfo[iRes].vtxGen_z = -999;
    evInfo[iRes].vtx3D_z = -999;
    evInfo[iRes].vtx_z = -999;
    evInfo[iRes].vtx_t = -999;

    evInfo[iRes].photon_pt.clear();
    evInfo[iRes].photon_eta.clear();
    evInfo[iRes].photon_phi.clear();
    evInfo[iRes].photon_isPrompt.clear();
    evInfo[iRes].photon_hasConversionTracks.clear();
    evInfo[iRes].photon_sigmaIetaIeta.clear();
    evInfo[iRes].photon_r9.clear();

    for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++){
      evInfo[iRes].photon_chIso[iCone].clear();
      evInfo[iRes].photon_chIso_dT[iCone][iRes].clear();
    }

    if (saveTracks_){
      evInfo[iRes].track_t.clear();
      evInfo[iRes].track_dz.clear();
      evInfo[iRes].track_pt.clear();
      evInfo[iRes].track_eta.clear();
      evInfo[iRes].track_phi.clear();
    }
  }
}




bool isPromptPhoton(const reco::Photon& photon, const edm::View<reco::GenParticle>& genParticles)
{
  bool isPrompt = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 22) continue;
    if ( !genp.isPromptFinalState() ) continue;
    if ( genp.pt() < 5.0 ) continue;
    double dr = deltaR(photon,genp);
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





XYZTLorentzVector correctP4(const reco::Photon &photon, const reco::Vertex& vtx){
  // -- recompute p4 wrt to the chosen vertex (taken from flashgg)                                                                                                                                      
  math::XYZVector vtx_Pos( vtx.x(), vtx.y(), vtx.z() );
  math::XYZVector sc_Pos( photon.superCluster()->x(), photon.superCluster()->y(), photon.superCluster()->z());
  math::XYZVector direction = sc_Pos - vtx_Pos;
  math::XYZVector p = ( direction.Unit() ) * ( photon.energy() );
  math::XYZTLorentzVector corr_p4( p.x(), p.y(), p.z(), photon.energy() );
  
  return corr_p4;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIsolationAnalyzer);

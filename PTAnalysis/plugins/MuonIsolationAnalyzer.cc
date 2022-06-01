// -*- C++ -*-%//
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
  pfcandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidateTag"))),
  genPartToken_(consumes<edm::View<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genPartTag"))),
  genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
  genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))),      
  genJetsToken_(consumes<View<reco::GenJet> >(iConfig.getUntrackedParameter<InputTag>("genJetsTag"))),
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonsTag")))
{
  timeResolutions_ = iConfig.getUntrackedParameter<vector<double> >("timeResolutions");
  isoConeDR_       = iConfig.getUntrackedParameter<double>("isoConeDR");
  saveTracks_      = iConfig.getUntrackedParameter<bool>("saveTracks");
  maxDz_           = iConfig.getUntrackedParameter<double>("maxDz");
  minDr_           = iConfig.getUntrackedParameter<double>("minDr");
  minTrackPt_      = iConfig.getUntrackedParameter<double>("minTrackPt");
  useVertexClosestToGen_ = iConfig.getUntrackedParameter<bool>("useVertexClosestToGen");
  btlEfficiency_ = iConfig.getUntrackedParameter<double>("btlEfficiency");
  etlEfficiency_ = iConfig.getUntrackedParameter<double>("etlEfficiency");

  //Now do what ever initialization is needed
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    if (timeResolutions_[iRes] == -1){
      cout << iRes << " --> Using default time resolution " << endl;
      eventTree[iRes] = fs_->make<TTree>("tree", "tree");
    }
    else {
      eventTree[iRes] = fs_->make<TTree>( Form("tree_%dps",int(timeResolutions_[iRes]*1000) ), Form("tree_%dps", int(timeResolutions_[iRes]*1000)) );
    }
  }


  gRandom  = new TRandom();
  gRandom2 = new TRandom();

}


MuonIsolationAnalyzer::~MuonIsolationAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  cout << "CIAO CIAO" << endl; 
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
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_, muons);

  // -- get the PFCandidate collection
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandToken_, pfcands);

  // -- get the gen particles collection
  edm::Handle<edm::View<pat::PackedGenParticle>> hPackedGenParticles;
  iEvent.getByToken(genPartToken_, hPackedGenParticles);
  auto genParticles = hPackedGenParticles.isValid() ? hPackedGenParticles.product() : nullptr;


  // -- get the genXYZ
  Handle<genXYZ> genXYZH;
  iEvent.getByToken(genXYZToken_, genXYZH);
  
  // -- get the genT0
  Handle<float> genT0H;
  iEvent.getByToken(genT0Token_, genT0H);
    
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
  SimVertex genPV; 
  auto xyz = genXYZH.product();
  auto t = *genT0H.product();
  auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
  genPV = SimVertex(v, t);
    
  float mindz = 999999.;
  int pv_index_3D = -1;
  int pv_index_4D = -1;

  if (useVertexClosestToGen_){
    // -- find the reco vertex closest to the gen vertex (3D)
    for(unsigned int ivtx=0; ivtx < vertices3D.size(); ivtx++ ){
      const reco::Vertex& vtx = vertices3D[ivtx];
      float dz = std::abs(vtx.z() - genPV.position().z());
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
      float dz = std::abs(vtx.z() - genPV.position().z());
      if( dz < mindz )
	{
	  mindz = dz;
	  pv_index_4D = ivtx;
	}
    }
  }

  // -- if PV index == -1, use highest ranked vertex 
  if (pv_index_3D==-1) pv_index_3D = 0;
  if (pv_index_4D==-1) pv_index_4D = 0;

  //cout << "pv_index_3D = " << pv_index_3D <<endl;
  //cout << "pv_index_4D = " << pv_index_4D <<endl;


  // -- get isolation around a candidate photon 
  // --- using only vtx closest to gen vtx
  const reco::Vertex& vtx4D = vertices4D[pv_index_4D];
  const reco::Vertex& vtx3D = vertices3D[pv_index_3D];

  int muonIndex =  -1;

  const int nResol = timeResolutions_.size();

  // --- start loop over muons
  for (const pat::Muon &muon : *muons){ 
    
    // -- minimal checks
    if (muon.track().isNull()) continue;
    if (muon.pt() < 5.) continue;
    
    muonIndex++;

    // -- check if prompt or fake muon
    bool isPromptMu = isPromptMuon(muon, *genParticles);
    bool isMatchedGenJet  = isMatchedToGenJet(muon, genJets);
    bool isFromTauDecay = isFromTau(muon, *genParticles);

    // -- compute charged isolations
    float chIso_dZ05_simVtx = 0.;
    float chIso_dZ1_simVtx = 0.;
    float chIso_dZ2_simVtx = 0.;
    float chIso_dZ05 = 0.;
    float chIso_dZ1 = 0.;
    float chIso_dZ2 = 0.;
    float chIso_reldZ = 0.;
    float chIso_dZmu05 = 0.;
    float chIso_dZmu1 = 0.;
    float chIso_dZmu2 = 0.;

    float chIso_dZ05_dT_simVtx[nResol];
    float chIso_dZ1_dT_simVtx[nResol];
    float chIso_dZ2_dT_simVtx[nResol];
    float chIso_dZ05_dT[nResol];
    float chIso_dZ1_dT[nResol];
    float chIso_dZ2_dT[nResol];
    float chIso_reldZ_dT[nResol];
    float chIso_dZmu05_dTmu[nResol];
    float chIso_dZmu1_dTmu[nResol];
    float chIso_dZmu2_dTmu[nResol];
    float time[nResol]; 
    
    // -- initialize 
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      chIso_dZ05_dT_simVtx[iRes] = 0.;
      chIso_dZ1_dT_simVtx[iRes]  = 0.;
      chIso_dZ2_dT_simVtx[iRes]  = 0.;
      chIso_dZ05_dT[iRes] = 0.;
      chIso_dZ1_dT[iRes] = 0.;
      chIso_dZ2_dT[iRes] = 0.;
      chIso_reldZ_dT[iRes] = 0.;
      chIso_dZmu05_dTmu[iRes] = 0.;
      chIso_dZmu1_dTmu[iRes] = 0.;
      chIso_dZmu2_dTmu[iRes] = 0.;
      time[iRes] = 0.;
    }

    
    float muonTime =  0;
      
    // -- loop over charged pf candidates
    for (const pat::PackedCandidate &pfcand : *pfcands){
      //for(unsigned icand = 0; icand < pfcands->size(); ++icand) {
      //const reco::PFCandidate& pfcand = *pfcands[icand];
      if (pfcand.charge() == 0 ) continue;
      
      // -- skip the track if it is the muon track
      //reco::TrackRef trackRef = pfcand.bestTrack();
      //if ( trackRef.isNull() ) continue;
      //if ( !(trackRef->quality(reco::TrackBase::highPurity)) ) continue;
      //if ( trackRef == muon.track() ) muonTime = pfcand.time();
      //if ( trackRef == muon.track() ) continue;

      if (pfcand.hasTrackDetails() == false) continue;
      if (pfcand.trackHighPurity() == false) continue;
      const reco::Track *track = pfcand.bestTrack();
      if ( track == muon.bestTrack() ) muonTime = pfcand.time();
      //if ( &track == muon.track() ) continue;
      

      //-- use tracks with pT above thtreshold
      if ( pfcand.pt() < minTrackPt_) continue;

      // -- compute dz, dxy 
      float dz4D  = std::abs( pfcand.dz(vtx4D.position()) );
      float dz3D  = std::abs( pfcand.dz(vtx3D.position()) );
      float dxy4D = std::abs( pfcand.dxy(vtx4D.position()) );
      float dxy3D = std::abs( pfcand.dxy(vtx3D.position()) );
   
      float dz4Drel = std::abs(dz4D/std::sqrt(pfcand.dzError()*pfcand.dzError() + vtx4D.zError()*vtx4D.zError()));
      float dz3Drel = std::abs(dz3D/std::sqrt(pfcand.dzError()*pfcand.dzError() + vtx3D.zError()*vtx3D.zError()));

      //float dzsim   = std::abs( pfcand.dz(Point(v)) ); 
      //float dxysim  = std::abs( pfcand.dxy(Point(v)) ); 

      float dzsim   = 0; 
      float dxysim  = 0;

      //float dzmu  = std::abs( pfcand.dz(vtx4D.position()) - muon.track()->dz(vtx4D.position()) );
      //float dxymu = std::abs( pfcand.dxy(vtx4D.position()) - muon.track()->dxy(vtx4D.position()) );
      float dzmu  =  0;
      float dxymu  =  0;


      float dr  = deltaR(muon.eta(), muon.phi(), pfcand.eta(), pfcand.phi());

      // --- no timing 
      if (dr > minDr_ && dr < isoConeDR_){
	// -- sim vertex 
	if (dzsim < 0.05 && dxysim < 0.02) { chIso_dZ05_simVtx+= pfcand.pt(); }
	if (dzsim < 0.1  && dxysim < 0.02) { chIso_dZ1_simVtx+= pfcand.pt(); }
	if (dzsim < 0.2  && dxysim < 0.02) { chIso_dZ2_simVtx+= pfcand.pt(); }
	
	// -- reco vtx closest to the sim one
	if (dz3D < 0.05  &&  dxy3D < 0.02) { chIso_dZ05+= pfcand.pt(); }
	if (dz3D < 0.1   &&  dxy3D < 0.02) { chIso_dZ1+= pfcand.pt(); }
	if (dz3D < 0.2   &&  dxy3D < 0.02) { chIso_dZ2+= pfcand.pt(); }
	
	// -- using reco vtx closest to the sim one and cut on relative dz  
	if (dz3Drel < 3.0  && dxy3D < 0.02) { chIso_reldZ+= pfcand.pt(); } 
	
	// -- dz wrt to muon
	if (dzmu < 0.05) { chIso_dZmu05+= pfcand.pt(); }
	if (dzmu < 0.1 ) { chIso_dZmu1+= pfcand.pt(); }
	if (dzmu < 0.2 ) { chIso_dZmu2+= pfcand.pt(); }
      }
      

      // --- with timing
      if (dr > minDr_ && dr < isoConeDR_){ 	  

	// -- emulate BTL and ETL efficiency
	double rndEff = gRandom2->Uniform(0.,1.); 
	if ( std::abs(pfcand.eta()) < 1.5 && rndEff > btlEfficiency_ ) continue; 	
	if ( std::abs(pfcand.eta()) > 1.5 && rndEff > etlEfficiency_ ) continue; 	

	for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){                                                                                                    
	  double defaultTimeResol  = 0.;
	  double extra_resol = 0.;

	  double targetTimeResol = timeResolutions_[iRes];
	  if ( pfcand.timeError()>-1 ) defaultTimeResol  = double(pfcand.timeError()); 

	  if (targetTimeResol == -1) extra_resol = 0;
	  else {
	    if ( targetTimeResol > defaultTimeResol) extra_resol = sqrt(targetTimeResol*targetTimeResol - defaultTimeResol*defaultTimeResol);
	  }
    
	  double dtsim = 0.;
	  double dt = 0.;
	  double dtmu = 0.;
	  time[iRes] = -999.;
	  if ( pfcand.timeError()>-1 ) {
	    double rnd = gRandom->Gaus(0., extra_resol);
	    //cout << "target time resol = "<< targetTimeResol << "  extra_resol = "<< extra_resol << "  rnd = " << rnd <<endl;
	    time[iRes] = pfcand.time() + rnd;
	    dtsim = std::abs(time[iRes] - genPV.position().t()*1000000000.);
	    dt    = std::abs(time[iRes] - vtx4D.t());
	    dtmu  = std::abs(time[iRes] - muonTime);
	  }
	  else{
	    dtsim = 0.;
	    dt    = 0.;
	    dtmu  = 0.;
	  }
	    
	  // -- sim vertex
	  if (dtsim < 3*targetTimeResol && dzsim < 0.05 && dxysim < 0.02 ) { chIso_dZ05_dT_simVtx[iRes]+= pfcand.pt();}
	  if (dtsim < 3*targetTimeResol && dzsim < 0.1  && dxysim < 0.02 ) { chIso_dZ1_dT_simVtx[iRes]+= pfcand.pt();}
	  if (dtsim < 3*targetTimeResol && dzsim < 0.2  && dxysim < 0.02 ) { chIso_dZ2_dT_simVtx[iRes]+= pfcand.pt();}
	  
	  // -- reco vtx closest to the sim one
	  if (dt < 3*targetTimeResol && dz4D < 0.05 &&  dxy4D < 0.02) { chIso_dZ05_dT[iRes]+= pfcand.pt(); }
	  if (dt < 3*targetTimeResol && dz4D < 0.1  &&  dxy4D < 0.02) { chIso_dZ1_dT[iRes]+= pfcand.pt(); }
	  if (dt < 3*targetTimeResol && dz4D < 0.2  &&  dxy4D < 0.02) { chIso_dZ2_dT[iRes]+= pfcand.pt(); }
	  
	  // -- using reco vtx closest to the sim one   and cut on relative dz
	  if (dt < 3*targetTimeResol && dz4Drel < 3.0  && dxy4D < 0.02) { chIso_reldZ_dT[iRes]+= pfcand.pt(); }
	  
	  // -- dT wrt to muon 
	  if ( dtmu < sqrt(2)*3*targetTimeResol && dzmu < 0.05) { chIso_dZmu05_dTmu[iRes]+= pfcand.pt(); }
	  if ( dtmu < sqrt(2)*3*targetTimeResol && dzmu < 0.1 ) { chIso_dZmu1_dTmu[iRes]+= pfcand.pt(); }
	  if ( dtmu < sqrt(2)*3*targetTimeResol && dzmu < 0.2 ) { chIso_dZmu2_dTmu[iRes]+= pfcand.pt(); }
	  
	  
	  // -- save info for tracks in the isolation cone (only for DR = 0.3)
	  if (saveTracks_ && (dz4D < 1.0 || dz3D < 1. || dzmu < 1. ) ) { // save a subset of tracks with loose dz selection
	    evInfo[iRes].track_t.push_back(time[iRes]); 
	    evInfo[iRes].track_dz4D.push_back(pfcand.dz( vtx4D.position() )); 
	    evInfo[iRes].track_dz3D.push_back(pfcand.dz( vtx3D.position() )); 
	    evInfo[iRes].track_dxy4D.push_back(pfcand.dxy( vtx4D.position() )); 
	    evInfo[iRes].track_dxy3D.push_back(pfcand.dxy( vtx3D.position() )); 
	    evInfo[iRes].track_pt.push_back(pfcand.pt());
	    evInfo[iRes].track_eta.push_back(pfcand.eta());
	    evInfo[iRes].track_phi.push_back(pfcand.phi());
	    evInfo[iRes].track_muIndex.push_back(muonIndex);
	  }
	  
	}// end loop over time resolutions                                                                                                                                    
      }
    
    }// end loop over tracks
        

    // -- fill muon info for each resolution scenario
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes].muon_pt.push_back(muon.pt());  
      evInfo[iRes].muon_eta.push_back(muon.eta());  
      evInfo[iRes].muon_phi.push_back(muon.phi()); 
      evInfo[iRes].muon_isLoose.push_back(muon::isLooseMuon(muon)); 
      evInfo[iRes].muon_isMedium.push_back(muon::isMediumMuon(muon));
      int isTight3D = muon::isTightMuon(muon, vtx3D);
      int isTight4D = muon::isTightMuon(muon, vtx4D);
      evInfo[iRes].muon_isTight3D.push_back(isTight3D); 
      evInfo[iRes].muon_isTight4D.push_back(isTight4D); 
      evInfo[iRes].muon_dz4D.push_back( muon.track()->dz(vtx4D.position()) );
      evInfo[iRes].muon_dxy4D.push_back( muon.track()->dxy(vtx4D.position()) );
      evInfo[iRes].muon_dz3D.push_back( muon.track()->dz(vtx3D.position()) );
      evInfo[iRes].muon_dxy3D.push_back( muon.track()->dxy(vtx3D.position()) );
      evInfo[iRes].muon_t.push_back(muonTime); 
      evInfo[iRes].muon_isPrompt.push_back(isPromptMu);
      evInfo[iRes].muon_isMatchedToGenJet.push_back(isMatchedGenJet);
      evInfo[iRes].muon_isFromTauDecay.push_back(isFromTauDecay);
        
      evInfo[iRes].muon_chIso_dZ05_simVtx.push_back(chIso_dZ05_simVtx);
      evInfo[iRes].muon_chIso_dZ05_dT_simVtx.push_back(chIso_dZ05_dT_simVtx[iRes]);
      
      evInfo[iRes].muon_chIso_dZ1_simVtx.push_back(chIso_dZ1_simVtx);
      evInfo[iRes].muon_chIso_dZ1_dT_simVtx.push_back(chIso_dZ1_dT_simVtx[iRes]);
      
      evInfo[iRes].muon_chIso_dZ2_simVtx.push_back(chIso_dZ2_simVtx);
      evInfo[iRes].muon_chIso_dZ2_dT_simVtx.push_back(chIso_dZ2_dT_simVtx[iRes]);
      
      evInfo[iRes].muon_chIso_dZ05.push_back(chIso_dZ05);  
      evInfo[iRes].muon_chIso_dZ05_dT.push_back(chIso_dZ05_dT[iRes]);
      
      evInfo[iRes].muon_chIso_dZ1.push_back(chIso_dZ1);  
      evInfo[iRes].muon_chIso_dZ1_dT.push_back(chIso_dZ1_dT[iRes]);
      
      evInfo[iRes].muon_chIso_dZ2.push_back(chIso_dZ2);  
      evInfo[iRes].muon_chIso_dZ2_dT.push_back(chIso_dZ2_dT[iRes]);
      
      evInfo[iRes].muon_chIso_reldZ.push_back(chIso_reldZ);  
      evInfo[iRes].muon_chIso_reldZ_dT.push_back(chIso_reldZ_dT[iRes]);
      
      evInfo[iRes].muon_chIso_dZmu05.push_back(chIso_dZmu05);
      evInfo[iRes].muon_chIso_dZmu05_dTmu.push_back(chIso_dZmu05_dTmu[iRes]);
      
      evInfo[iRes].muon_chIso_dZmu1.push_back(chIso_dZmu1);
      evInfo[iRes].muon_chIso_dZmu1_dTmu.push_back(chIso_dZmu1_dTmu[iRes]);
      
      evInfo[iRes].muon_chIso_dZmu2.push_back(chIso_dZmu2);
      evInfo[iRes].muon_chIso_dZmu2_dTmu.push_back(chIso_dZmu2_dTmu[iRes]);
      
    } // end loop over time resolutions

  }// end loop over muons


  // -- fill event info
  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
    evInfo[iRes].npu = nPU;
    evInfo[iRes].vtx4D_t = vtx4D.t();
    evInfo[iRes].vtx4D_tErr = vtx4D.tError();
    evInfo[iRes].vtx4D_z = vtx4D.z();
    evInfo[iRes].vtx4D_zErr = vtx4D.zError();
    evInfo[iRes].vtx3D_z = vtx3D.z();
    evInfo[iRes].vtx3D_zErr = vtx3D.zError();
    evInfo[iRes].vtxGen_z = genPV.position().z();
    evInfo[iRes].vtxGen_t = genPV.position().t();
    evInfo[iRes].vtx3D_isFake = vtx3D.isFake();
    evInfo[iRes].vtx4D_isFake = vtx4D.isFake();
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
    eventTree[iRes]->Branch( "vtx3D_zErr",      &evInfo[iRes].vtx3D_zErr);
    eventTree[iRes]->Branch( "vtx4D_z",         &evInfo[iRes].vtx4D_z);
    eventTree[iRes]->Branch( "vtx4D_zErr",      &evInfo[iRes].vtx4D_zErr);
    eventTree[iRes]->Branch( "vtx4D_t",         &evInfo[iRes].vtx4D_t);
    eventTree[iRes]->Branch( "vtx4D_tErr",      &evInfo[iRes].vtx4D_tErr);
    eventTree[iRes]->Branch( "vtx4D_isFake",    &evInfo[iRes].vtx4D_isFake);
    eventTree[iRes]->Branch( "vtx3D_isFake",    &evInfo[iRes].vtx3D_isFake);
    eventTree[iRes]->Branch( "muon_pt",         &evInfo[iRes].muon_pt);  
    eventTree[iRes]->Branch( "muon_eta",        &evInfo[iRes].muon_eta);  
    eventTree[iRes]->Branch( "muon_phi",        &evInfo[iRes].muon_phi);  
    eventTree[iRes]->Branch( "muon_isLoose",    &evInfo[iRes].muon_isLoose);  
    eventTree[iRes]->Branch( "muon_isMedium",   &evInfo[iRes].muon_isMedium);  
    eventTree[iRes]->Branch( "muon_isTight3D",  &evInfo[iRes].muon_isTight3D);  
    eventTree[iRes]->Branch( "muon_isTight4D",  &evInfo[iRes].muon_isTight4D);  
    eventTree[iRes]->Branch( "muon_dz3D",       &evInfo[iRes].muon_dz3D);  
    eventTree[iRes]->Branch( "muon_dz4D",       &evInfo[iRes].muon_dz4D);  
    eventTree[iRes]->Branch( "muon_dxy3D",      &evInfo[iRes].muon_dxy3D);  
    eventTree[iRes]->Branch( "muon_dxy4D",      &evInfo[iRes].muon_dxy4D);  
    eventTree[iRes]->Branch( "muon_t",          &evInfo[iRes].muon_t);  
    eventTree[iRes]->Branch( "muon_isPrompt",   &evInfo[iRes].muon_isPrompt);  
    eventTree[iRes]->Branch( "muon_isMatchedToGenJet",&evInfo[iRes].muon_isMatchedToGenJet);  
    eventTree[iRes]->Branch( "muon_isFromTauDecay",   &evInfo[iRes].muon_isFromTauDecay);  

    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ05_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ05_simVtx);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ05_dT_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ05_dT_simVtx);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ1_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ1_simVtx);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ1_dT_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ1_dT_simVtx);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ2_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ2_simVtx);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ2_dT_simVtx",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ2_dT_simVtx);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ05",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ05);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ05_dT",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ05_dT);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ1",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ1);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ1_dT",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ1_dT);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ2",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ2);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZ2_dT",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZ2_dT);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_reldZ",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_reldZ);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_reldZ_dT",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_reldZ_dT);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu05",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu05);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu05_dTmu",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu05_dTmu);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu1",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu1);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu1_dTmu",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu1_dTmu);  
    
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu2",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu2);  
    eventTree[iRes]->Branch( Form("muon_chIso%.2d_dZmu2_dTmu",int(isoConeDR_*10) ), &evInfo[iRes].muon_chIso_dZmu2_dTmu);  


    if (saveTracks_){
      eventTree[iRes]->Branch( "track_t",      &evInfo[iRes].track_t);
      eventTree[iRes]->Branch( "track_dz4D",   &evInfo[iRes].track_dz4D);
      eventTree[iRes]->Branch( "track_dz3D",   &evInfo[iRes].track_dz3D);
      eventTree[iRes]->Branch( "track_dxy4D",  &evInfo[iRes].track_dxy4D);
      eventTree[iRes]->Branch( "track_dxy3D",  &evInfo[iRes].track_dxy3D);
      eventTree[iRes]->Branch( "track_pt",     &evInfo[iRes].track_pt);
      eventTree[iRes]->Branch( "track_eta",    &evInfo[iRes].track_eta);
      eventTree[iRes]->Branch( "track_phi",    &evInfo[iRes].track_phi);
      eventTree[iRes]->Branch( "track_muIndex",&evInfo[iRes].track_muIndex);
    }
  
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonIsolationAnalyzer::endJob() 
{

  cout << "Analysis finished " <<endl;
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
    evInfo[iRes].vtx3D_zErr = -999;
    evInfo[iRes].vtx4D_z = -999;
    evInfo[iRes].vtx4D_zErr = -999;
    evInfo[iRes].vtx4D_t = -999;
    evInfo[iRes].vtx4D_tErr = -999;
    evInfo[iRes].vtx3D_isFake = -999;
    evInfo[iRes].vtx4D_isFake = -999;

    evInfo[iRes].muon_pt.clear();
    evInfo[iRes].muon_eta.clear();
    evInfo[iRes].muon_phi.clear();
    evInfo[iRes].muon_isLoose.clear();
    evInfo[iRes].muon_isMedium.clear();
    evInfo[iRes].muon_isTight3D.clear();    
    evInfo[iRes].muon_isTight4D.clear();
    evInfo[iRes].muon_dz3D.clear();
    evInfo[iRes].muon_dxy3D.clear();
    evInfo[iRes].muon_dz4D.clear();
    evInfo[iRes].muon_dxy4D.clear();
    evInfo[iRes].muon_t.clear();
    evInfo[iRes].muon_isPrompt.clear();
    evInfo[iRes].muon_isMatchedToGenJet.clear();
    evInfo[iRes].muon_isFromTauDecay.clear();

    evInfo[iRes].muon_chIso_dZ05_simVtx.clear();
    evInfo[iRes].muon_chIso_dZ05_dT_simVtx.clear();
    evInfo[iRes].muon_chIso_dZ1_simVtx.clear();
    evInfo[iRes].muon_chIso_dZ1_dT_simVtx.clear();
    evInfo[iRes].muon_chIso_dZ2_simVtx.clear();
    evInfo[iRes].muon_chIso_dZ2_dT_simVtx.clear();

    evInfo[iRes].muon_chIso_dZ05.clear();
    evInfo[iRes].muon_chIso_dZ05_dT.clear();
    evInfo[iRes].muon_chIso_dZ1.clear();
    evInfo[iRes].muon_chIso_dZ1_dT.clear();
    evInfo[iRes].muon_chIso_dZ2.clear();
    evInfo[iRes].muon_chIso_dZ2_dT.clear();

    evInfo[iRes].muon_chIso_reldZ.clear();
    evInfo[iRes].muon_chIso_reldZ_dT.clear();

    evInfo[iRes].muon_chIso_dZmu05.clear();
    evInfo[iRes].muon_chIso_dZmu05_dTmu.clear();
    evInfo[iRes].muon_chIso_dZmu1.clear();
    evInfo[iRes].muon_chIso_dZmu1_dTmu.clear();
    evInfo[iRes].muon_chIso_dZmu2.clear();
    evInfo[iRes].muon_chIso_dZmu2_dTmu.clear();


    if (saveTracks_){
      evInfo[iRes].track_t.clear();
      evInfo[iRes].track_dz4D.clear();
      evInfo[iRes].track_dz3D.clear();
      evInfo[iRes].track_dxy4D.clear();
      evInfo[iRes].track_dxy3D.clear();
      evInfo[iRes].track_pt.clear();
      evInfo[iRes].track_eta.clear();
      evInfo[iRes].track_phi.clear();
      evInfo[iRes].track_muIndex.clear();
    }
  }
}



// --- matching to gen muon
bool isPromptMuon(const pat::Muon& muon, const edm::View<pat::PackedGenParticle>& genParticles)
//bool isPromptMuon(const pat::Muon& muon, const std::vector<pat::PackedGenParticle>* genParticles)
{
  bool isPrompt = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const pat::PackedGenParticle& genp = genParticles.at(ip);
    //    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 13) continue;
    if (genp.status() != 1 ) continue; 
    if ( !genp.isPromptFinalState() ) continue;
    double dr = deltaR(muon,genp);
    if (dr > 0.2){ 
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
bool isMatchedToGenJet(const pat::Muon& muon, const edm::View<reco::GenJet>& genJets)
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



// --- matching to muons from tau decays
bool isFromTau(const pat::Muon& muon, const edm::View<pat::PackedGenParticle>& genParticles)
{
  
  bool fromTau = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const pat::PackedGenParticle& genp = genParticles.at(ip);
    //const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 13) continue;
    if ( !genp.isDirectPromptTauDecayProductFinalState() ) continue;
    double dr = deltaR(muon,genp);
    if (dr > 0.2){
      continue;
    }
    else{
      fromTau=true;
      break;
    }
  }

  return fromTau;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationAnalyzer);

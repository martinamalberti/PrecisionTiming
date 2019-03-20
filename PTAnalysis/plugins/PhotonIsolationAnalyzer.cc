// -*- C++ -*-%//
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
#include "PrecisionTiming/PTAnalysis/interface/PhotonIsolationAnalyzer.h"

#include <TMath.h>
#include <TRandom.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <limits>

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
  genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
  genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))),      
  genJetsToken_(consumes<View<reco::GenJet> >(iConfig.getUntrackedParameter<InputTag>("genJetsTag"))),
  barrelPhotonsToken_(consumes<View<reco::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("barrelPhotonsTag"))),
  endcapPhotonsToken_(consumes<View<reco::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapPhotonsTag"))),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  trackTimeErrToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeErrValueMapTag" ) ) )
{
  timeResolutions_ = iConfig.getUntrackedParameter<vector<double> >("timeResolutions");
  isoConeDR_       = iConfig.getUntrackedParameter<double>("isoConeDR");
  saveTracks_      = iConfig.getUntrackedParameter<bool>("saveTracks");
  maxDz_           = iConfig.getUntrackedParameter<double>("maxDz");
  maxDxy_          = iConfig.getUntrackedParameter<double>("maxDxy");
  minDr_           = iConfig.getUntrackedParameter<double>("minDr");
  btlMinTrackPt_   = iConfig.getUntrackedParameter<double>("btlMinTrackPt");
  etlMinTrackPt_   = iConfig.getUntrackedParameter<double>("etlMinTrackPt");
  useVertexClosestToGenZ_ = iConfig.getUntrackedParameter<bool>("useVertexClosestToGenZ");
  useVertexClosestToGenZT_ = iConfig.getUntrackedParameter<bool>("useVertexClosestToGenZT");
  btlEfficiency_ = iConfig.getUntrackedParameter<double>("btlEfficiency");
  etlEfficiency_ = iConfig.getUntrackedParameter<double>("etlEfficiency");

  //Now do what ever initialization is needed
  usesResource("TFileService");

  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
    evInfo[iRes] = new eventInfo;
    eventTree[iRes] = fs_->make<TTree>( Form("tree_%dps",int(timeResolutions_[iRes]*1000) ), Form("tree_%dps", int(timeResolutions_[iRes]*1000)) );
    cout << iRes << "  " << timeResolutions_[iRes] << "  " << eventTree[iRes]->GetName() <<endl;
  }


  gRandom  = new TRandom();
  gRandom2 = new TRandom();

}


PhotonIsolationAnalyzer::~PhotonIsolationAnalyzer()
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
  
  // -- get the trackTimeValueMap
  Handle<ValueMap<float> > trackTimeValueMap;
  iEvent.getByToken( trackTimeToken_, trackTimeValueMap );

  // -- get the trackTimeErrValueMap
  Handle<ValueMap<float> > trackTimeErrValueMap;
  iEvent.getByToken( trackTimeErrToken_, trackTimeErrValueMap );

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
  if ( GenVertexCollectionH.isValid()) {
    const vector<SimVertex>& genVertices = *GenVertexCollectionH;   
    genPV = genVertices.at(0);
  }
  else{
    auto xyz = genXYZH.product();
    auto t = *genT0H.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    genPV = SimVertex(v, t);
  }
  
  float mindz = std::numeric_limits<float>::max();
  int pv_index_3D = -1;
  int pv_index_4D = -1;

  if (useVertexClosestToGenZ_){
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
    mindz = std::numeric_limits<float>::max();
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


  if (useVertexClosestToGenZT_){
    // -- find the reco vertex closest to the gen vertex (3D)
    mindz = std::numeric_limits<float>::max();
    for(unsigned int ivtx=0; ivtx < vertices3D.size(); ivtx++ ){
      const reco::Vertex& vtx = vertices3D[ivtx];
      float dz = std::abs(vtx.z() - genPV.position().z())/vtx.zError() ;
      if( dz < mindz )
        {
          mindz = dz;
          pv_index_3D = ivtx;
        }
    }
    // -- find the reco vertex closest to the gen vertex (4D)
    float mindist = std::numeric_limits<float>::max();
    for(unsigned int ivtx=0; ivtx < vertices4D.size(); ivtx++ ){
      const reco::Vertex& vtx = vertices4D[ivtx];
      float dz = std::abs(vtx.z() - genPV.position().z()) / vtx.zError();
      if ( vtx.tError() <=0 ) continue;
      float dt = std::abs(vtx.t() - genPV.position().t()*1000000000.) / vtx.tError();
      float dist = sqrt(dz*dz+dt*dt);
      //cout << ivtx << "dist = " << dist << "   dz = "<< dz <<  "   dt = "<< dt << "vtx time = " << vtx.t() << "  tErr = " << vtx.tError() <<endl;
      if( dist < mindist )
        {
          mindist = dist;
          pv_index_4D = ivtx;
        }
    }
  }


  // -- if PV index == -1, use highest ranked vertex 
  if (pv_index_3D==-1) {
    cout << "pv_index_3D = " << pv_index_3D <<endl;
    pv_index_3D = 0;
  }
  if (pv_index_4D==-1) {
    cout << "pv_index_4D = " << pv_index_4D <<endl;
    pv_index_4D = 0;
  }


  // -- get isolation around a candidate photon 
  // --- using only vtx closest to gen vtx
  const reco::Vertex& vtx4D = vertices4D[pv_index_4D];
  const reco::Vertex& vtx3D = vertices3D[pv_index_3D];

  int photonIndex =  -1;

  const int nResol = timeResolutions_.size();

  // --- start loop over photon candidates
  edm::View<reco::Photon> photons ; 

  for (int iPhotonColl = 0; iPhotonColl < 2; iPhotonColl++){
    
    if (iPhotonColl==0) photons = barrelPhotons;
    if (iPhotonColl==1) photons = endcapPhotons;

    for(unsigned int ipho=0; ipho < photons.size(); ipho++ ){

      const reco::Photon& photon = photons[ipho];
      
      // -- minimal checks
      if(photon.pt() < 15.) continue;
      
      photonIndex++;
      
      // -- check if prompt or fake photon
      bool isPromptPho = isPromptPhoton(photon,genParticles);
      bool isMatched   = isPhotonMatchedToGenJet(photon, genJets);
      
      // -- compute charged isolations
      float chIso_dZ05_simVtx = 0.;
      float chIso_dZ1_simVtx = 0.;
      float chIso_dZ2_simVtx = 0.;
      float chIso_dZ05 = 0.;
      float chIso_dZ1 = 0.;
      float chIso_dZ2 = 0.;
      
      float chIso_dZ05_dT2s_simVtx[nResol];
      float chIso_dZ1_dT2s_simVtx[nResol];
      float chIso_dZ2_dT2s_simVtx[nResol];
      float chIso_dZ05_dT2s[nResol];
      float chIso_dZ1_dT2s[nResol];
      float chIso_dZ2_dT2s[nResol];
      
      float chIso_dZ05_dT3s_simVtx[nResol];
      float chIso_dZ1_dT3s_simVtx[nResol];
      float chIso_dZ2_dT3s_simVtx[nResol];
      float chIso_dZ05_dT3s[nResol];
      float chIso_dZ1_dT3s[nResol];
      float chIso_dZ2_dT3s[nResol];
      
      float chIso_dZ05_dT5s_simVtx[nResol];
      float chIso_dZ1_dT5s_simVtx[nResol];
      float chIso_dZ2_dT5s_simVtx[nResol];
      float chIso_dZ05_dT5s[nResol];
      float chIso_dZ1_dT5s[nResol];
      float chIso_dZ2_dT5s[nResol];
      
      float time[nResol]; 
      
      // -- initialize 
      for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
	chIso_dZ05_dT2s_simVtx[iRes] = 0.;
	chIso_dZ1_dT2s_simVtx[iRes]  = 0.;
	chIso_dZ2_dT2s_simVtx[iRes]  = 0.;
	chIso_dZ05_dT2s[iRes] = 0.;
	chIso_dZ1_dT2s[iRes] = 0.;
	chIso_dZ2_dT2s[iRes] = 0.;
	
	chIso_dZ05_dT3s_simVtx[iRes] = 0.;
	chIso_dZ1_dT3s_simVtx[iRes]  = 0.;
	chIso_dZ2_dT3s_simVtx[iRes]  = 0.;
	chIso_dZ05_dT3s[iRes] = 0.;
	chIso_dZ1_dT3s[iRes] = 0.;
	chIso_dZ2_dT3s[iRes] = 0.;
	
	chIso_dZ05_dT5s_simVtx[iRes] = 0.;
	chIso_dZ1_dT5s_simVtx[iRes]  = 0.;
	chIso_dZ2_dT5s_simVtx[iRes]  = 0.;
	chIso_dZ05_dT5s[iRes] = 0.;
	chIso_dZ1_dT5s[iRes] = 0.;
	chIso_dZ2_dT5s[iRes] = 0.;
	
	time[iRes] = 0.;
      }
      
      
      float photonTime =  -999.;
      
      // -- loop over charged pf candidates
      for(unsigned icand = 0; icand < pfcands.size(); ++icand) {
	const reco::PFCandidate& pfcand = pfcands[icand];
	if (pfcand.charge() == 0 ) continue;
	
	// -- skip the track if it is the muon track
	reco::TrackRef trackRef = pfcand.trackRef();
	if ( trackRef.isNull() ) continue;
	if ( !(trackRef->quality(reco::TrackBase::highPurity)) ) continue;
	
	//-- use tracks with pT above thtreshold
	if ( std::abs(pfcand.eta()) < 1.48 && pfcand.pt() < btlMinTrackPt_) continue;
	if ( std::abs(pfcand.eta()) > 1.48 && pfcand.pt() < etlMinTrackPt_) continue;

	// -- compute dz, dxy 
	float dz4D  = std::abs( trackRef->dz(vtx4D.position()) );
	float dz3D  = std::abs( trackRef->dz(vtx3D.position()) );
	float dxy4D = std::abs( trackRef->dxy(vtx4D.position()) );
	float dxy3D = std::abs( trackRef->dxy(vtx3D.position()) );
	
	float dzsim   = std::abs( trackRef->vz() - genPV.position().z() ); 
	float dxysim  = sqrt ( pow(trackRef->vx() - genPV.position().x(),2) + pow(trackRef->vy() - genPV.position().y(),2) ); 
	
	float dr  = deltaR(photon.eta(), photon.phi(), pfcand.eta(), pfcand.phi());
	

	// -- time 
	float pfcandtime    = (*trackTimeValueMap)[trackRef] ;
        float pfcandtimeErr = (*trackTimeErrValueMap)[trackRef] ;



	// --- no timing 
	if (dr > minDr_ && dr < isoConeDR_){
	  // -- sim vertex 
	  if (dzsim < 0.05 && dxysim < maxDxy_) { chIso_dZ05_simVtx+= pfcand.pt(); }
	  if (dzsim < 0.1  && dxysim < maxDxy_) { chIso_dZ1_simVtx+= pfcand.pt(); }
	  if (dzsim < 0.2  && dxysim < maxDxy_) { chIso_dZ2_simVtx+= pfcand.pt(); }
	  
	  // -- reco vtx closest to the sim one
	  if (dz3D < 0.05  &&  dxy3D < maxDxy_) { chIso_dZ05+= pfcand.pt(); }
	  if (dz3D < 0.1   &&  dxy3D < maxDxy_) { chIso_dZ1+= pfcand.pt(); }
	  if (dz3D < 0.2   &&  dxy3D < maxDxy_) { chIso_dZ2+= pfcand.pt(); }
	  
	}
	
	// --- with timing
	if (dr > minDr_ && dr < isoConeDR_){ 	  
	  
	  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){                                                                                                    
	    double targetTimeResol = timeResolutions_[iRes];
	    double defaultTimeResol  = 0.;
	    //if ( pfcand.isTimeValid() ) defaultTimeResol  = double(pfcand.timeError()); 
	    if ( pfcandtimeErr ) defaultTimeResol  = double(pfcandtimeErr); 
	    double extra_resol = sqrt(targetTimeResol*targetTimeResol - defaultTimeResol*defaultTimeResol);                 
	    double dtsim = 0.;
	    double dt = 0.;
	    time[iRes] = -999.;
	    //if ( pfcand.isTimeValid() && !isnan( pfcand.time() )) {
	    if ( pfcandtimeErr != -1 && !isnan( pfcandtime )) {
	      // -- emulate BTL and ETL efficiency
	      bool keepTrack = true;
	      double rndEff = gRandom2->Uniform(0.,1.); 
	      if ( std::abs(pfcand.eta()) < 1.5 && rndEff > btlEfficiency_ ) keepTrack = false; 	
	      if ( std::abs(pfcand.eta()) > 1.5 && rndEff > etlEfficiency_ ) keepTrack = false; 	
	      if (keepTrack) {
		// -- extra smearing to emulate different time resolution
		double rnd   = gRandom->Gaus(0., extra_resol);
		time[iRes] = pfcandtime + rnd;
		dtsim = std::abs(time[iRes] - genPV.position().t()*1000000000.);
		dt    = std::abs(time[iRes] - vtx4D.t());
	      }
	      else{
		dtsim = 0.;
		dt    = 0.;
	      }
	    }
	    else{
	      dtsim = 0.;
	      dt    = 0.;
	    }
	    
	    // -- sim vertex
	    if (dtsim < 2.*targetTimeResol && dzsim < 0.05 && dxysim < maxDxy_) { chIso_dZ05_dT2s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 2.*targetTimeResol && dzsim < 0.1  && dxysim < maxDxy_) { chIso_dZ1_dT2s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 2.*targetTimeResol && dzsim < 0.2  && dxysim < maxDxy_) { chIso_dZ2_dT2s_simVtx[iRes]+= pfcand.pt();}
	    
	    if (dtsim < 3.*targetTimeResol && dzsim < 0.05 && dxysim < maxDxy_) { chIso_dZ05_dT3s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 3.*targetTimeResol && dzsim < 0.1  && dxysim < maxDxy_) { chIso_dZ1_dT3s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 3.*targetTimeResol && dzsim < 0.2  && dxysim < maxDxy_) { chIso_dZ2_dT3s_simVtx[iRes]+= pfcand.pt();}
	    
	    if (dtsim < 5.*targetTimeResol && dzsim < 0.05 && dxysim < maxDxy_) { chIso_dZ05_dT5s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 5.*targetTimeResol && dzsim < 0.1  && dxysim < maxDxy_) { chIso_dZ1_dT5s_simVtx[iRes]+= pfcand.pt();}
	    if (dtsim < 5.*targetTimeResol && dzsim < 0.2  && dxysim < maxDxy_) { chIso_dZ2_dT5s_simVtx[iRes]+= pfcand.pt();}
	    
	    // -- reco vtx closest to the sim one
	    if (dt < 2.*targetTimeResol && dz4D < 0.05 &&  dxy4D < maxDxy_) { chIso_dZ05_dT2s[iRes]+= pfcand.pt(); }
	    if (dt < 2.*targetTimeResol && dz4D < 0.1  &&  dxy4D < maxDxy_) { chIso_dZ1_dT2s[iRes]+= pfcand.pt(); }
	    if (dt < 2.*targetTimeResol && dz4D < 0.2  &&  dxy4D < maxDxy_) { chIso_dZ2_dT2s[iRes]+= pfcand.pt(); }
	    
	    if (dt < 3.*targetTimeResol && dz4D < 0.05 &&  dxy4D < maxDxy_) { chIso_dZ05_dT3s[iRes]+= pfcand.pt(); }
	    if (dt < 3.*targetTimeResol && dz4D < 0.1  &&  dxy4D < maxDxy_) { chIso_dZ1_dT3s[iRes]+= pfcand.pt(); }
	    if (dt < 3.*targetTimeResol && dz4D < 0.2  &&  dxy4D < maxDxy_) { chIso_dZ2_dT3s[iRes]+= pfcand.pt(); }
	    
	    if (dt < 5.*targetTimeResol && dz4D < 0.05 &&  dxy4D < maxDxy_) { chIso_dZ05_dT5s[iRes]+= pfcand.pt(); }
	    if (dt < 5.*targetTimeResol && dz4D < 0.1  &&  dxy4D < maxDxy_) { chIso_dZ1_dT5s[iRes]+= pfcand.pt(); }
	    if (dt < 5.*targetTimeResol && dz4D < 0.2  &&  dxy4D < maxDxy_) { chIso_dZ2_dT5s[iRes]+= pfcand.pt(); }
	    
	    // -- save info for tracks in the isolation cone (only for DR = 0.3)
	    if (saveTracks_ && (dz4D < 1.0 || dz3D < 1. || dzsim < 1.) ) { // save a subset of tracks with loose dz selection
	      bool genMatching  = isMatchedToGenParticle(pfcand, genParticles);
	      evInfo[iRes]->track_t.push_back(time[iRes]); 
	      evInfo[iRes]->track_dz4D.push_back(trackRef->dz( vtx4D.position() )); 
	      evInfo[iRes]->track_dz3D.push_back(trackRef->dz( vtx3D.position() )); 
	      evInfo[iRes]->track_dxy4D.push_back(trackRef->dxy( vtx4D.position() )); 
	      evInfo[iRes]->track_dxy3D.push_back(trackRef->dxy( vtx3D.position() )); 
	      evInfo[iRes]->track_pt.push_back(pfcand.pt());
	      evInfo[iRes]->track_eta.push_back(pfcand.eta());
	      evInfo[iRes]->track_phi.push_back(pfcand.phi());
	      evInfo[iRes]->track_phoIndex.push_back(photonIndex);
	      evInfo[iRes]->track_isMatchedToGenParticle.push_back(genMatching);
	    }
	    
	  }// end loop over time resolutions                                                                                                                                    
	}
	
      }// end loop over tracks
        

      // -- fill photon info for each resolution scenario
      for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
	evInfo[iRes]->photon_pt.push_back(photon.pt());  
	evInfo[iRes]->photon_eta.push_back(photon.eta());  
	evInfo[iRes]->photon_phi.push_back(photon.phi()); 
	evInfo[iRes]->photon_t.push_back(photonTime); 
	evInfo[iRes]->photon_isPrompt.push_back(isPromptPho);
	evInfo[iRes]->photon_isMatchedToGenJet.push_back(isMatched);
	evInfo[iRes]->photon_hasConversionTracks.push_back(photon.hasConversionTracks());
	evInfo[iRes]->photon_hasPixelSeed.push_back(photon.hasPixelSeed());
	evInfo[iRes]->photon_sigmaIetaIeta.push_back(photon.sigmaIetaIeta());
	evInfo[iRes]->photon_r9.push_back(photon.r9());
	
	evInfo[iRes]->photon_chIso_dZ05_simVtx.push_back(chIso_dZ05_simVtx);
	evInfo[iRes]->photon_chIso_dZ05_dT2s_simVtx.push_back(chIso_dZ05_dT2s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ05_dT3s_simVtx.push_back(chIso_dZ05_dT3s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ05_dT5s_simVtx.push_back(chIso_dZ05_dT5s_simVtx[iRes]);
	
	evInfo[iRes]->photon_chIso_dZ1_simVtx.push_back(chIso_dZ1_simVtx);
	evInfo[iRes]->photon_chIso_dZ1_dT2s_simVtx.push_back(chIso_dZ1_dT2s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ1_dT3s_simVtx.push_back(chIso_dZ1_dT3s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ1_dT5s_simVtx.push_back(chIso_dZ1_dT5s_simVtx[iRes]);
	
	evInfo[iRes]->photon_chIso_dZ2_simVtx.push_back(chIso_dZ2_simVtx);
	evInfo[iRes]->photon_chIso_dZ2_dT2s_simVtx.push_back(chIso_dZ2_dT2s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ2_dT3s_simVtx.push_back(chIso_dZ2_dT3s_simVtx[iRes]);
	evInfo[iRes]->photon_chIso_dZ2_dT5s_simVtx.push_back(chIso_dZ2_dT5s_simVtx[iRes]);
	
	evInfo[iRes]->photon_chIso_dZ05.push_back(chIso_dZ05);  
	evInfo[iRes]->photon_chIso_dZ05_dT2s.push_back(chIso_dZ05_dT2s[iRes]);
	evInfo[iRes]->photon_chIso_dZ05_dT3s.push_back(chIso_dZ05_dT3s[iRes]);
	evInfo[iRes]->photon_chIso_dZ05_dT5s.push_back(chIso_dZ05_dT5s[iRes]);
	
	evInfo[iRes]->photon_chIso_dZ1.push_back(chIso_dZ1);  
	evInfo[iRes]->photon_chIso_dZ1_dT2s.push_back(chIso_dZ1_dT2s[iRes]);
	evInfo[iRes]->photon_chIso_dZ1_dT3s.push_back(chIso_dZ1_dT3s[iRes]);
	evInfo[iRes]->photon_chIso_dZ1_dT5s.push_back(chIso_dZ1_dT5s[iRes]);
	
	evInfo[iRes]->photon_chIso_dZ2.push_back(chIso_dZ2);  
	evInfo[iRes]->photon_chIso_dZ2_dT2s.push_back(chIso_dZ2_dT2s[iRes]);
	evInfo[iRes]->photon_chIso_dZ2_dT3s.push_back(chIso_dZ2_dT3s[iRes]);
	evInfo[iRes]->photon_chIso_dZ2_dT5s.push_back(chIso_dZ2_dT5s[iRes]);
	
      } // end loop over time resolutions
      
    }// end loop over photons

  }// end loop over photon collections (barrel/endcap)

  // -- fill event info
  for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
    evInfo[iRes]->npu = nPU;
    evInfo[iRes]->vtx4D_t = vtx4D.t();
    evInfo[iRes]->vtx4D_tErr = vtx4D.tError();
    evInfo[iRes]->vtx4D_z = vtx4D.z();
    evInfo[iRes]->vtx4D_zErr = vtx4D.zError();
    evInfo[iRes]->vtx3D_z = vtx3D.z();
    evInfo[iRes]->vtx3D_zErr = vtx3D.zError();
    evInfo[iRes]->vtxGen_z = genPV.position().z();
    evInfo[iRes]->vtxGen_t = genPV.position().t() * 1.E09;
    evInfo[iRes]->vtx3D_isFake = vtx3D.isFake();
    evInfo[iRes]->vtx4D_isFake = vtx4D.isFake();
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
  //  fout = new TFile("test.root","recreate");
  
  for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
  
    eventTree[iRes]->Branch( "npu",             &evInfo[iRes]->npu);
    eventTree[iRes]->Branch( "vtxGen_z",        &evInfo[iRes]->vtxGen_z);
    eventTree[iRes]->Branch( "vtxGen_t",        &evInfo[iRes]->vtxGen_t);
    eventTree[iRes]->Branch( "vtx3D_z",         &evInfo[iRes]->vtx3D_z);
    eventTree[iRes]->Branch( "vtx3D_zErr",      &evInfo[iRes]->vtx3D_zErr);
    eventTree[iRes]->Branch( "vtx4D_z",         &evInfo[iRes]->vtx4D_z);
    eventTree[iRes]->Branch( "vtx4D_zErr",      &evInfo[iRes]->vtx4D_zErr);
    eventTree[iRes]->Branch( "vtx4D_t",         &evInfo[iRes]->vtx4D_t);
    eventTree[iRes]->Branch( "vtx4D_tErr",      &evInfo[iRes]->vtx4D_tErr);
    eventTree[iRes]->Branch( "vtx4D_isFake",    &evInfo[iRes]->vtx4D_isFake);
    eventTree[iRes]->Branch( "vtx3D_isFake",    &evInfo[iRes]->vtx3D_isFake);
    eventTree[iRes]->Branch( "photon_pt",         &evInfo[iRes]->photon_pt);  
    eventTree[iRes]->Branch( "photon_eta",        &evInfo[iRes]->photon_eta);  
    eventTree[iRes]->Branch( "photon_phi",        &evInfo[iRes]->photon_phi);  
    eventTree[iRes]->Branch( "photon_t",          &evInfo[iRes]->photon_t);  
    eventTree[iRes]->Branch( "photon_isPrompt",   &evInfo[iRes]->photon_isPrompt);  
    eventTree[iRes]->Branch( "photon_isMatchedToGenJet",&evInfo[iRes]->photon_isMatchedToGenJet);  
    eventTree[iRes]->Branch( "photon_hasConversionTracks",   &evInfo[iRes]->photon_hasConversionTracks);
    eventTree[iRes]->Branch( "photon_hasPixelSeed",   &evInfo[iRes]->photon_hasPixelSeed);
    eventTree[iRes]->Branch( "photon_sigmaIetaIeta",   &evInfo[iRes]->photon_sigmaIetaIeta);
    eventTree[iRes]->Branch( "photon_r9",   &evInfo[iRes]->photon_r9);

    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT2s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT2s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT3s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT3s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT5s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT5s_simVtx);  
    
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT2s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT2s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT3s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT3s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT5s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT5s_simVtx);  
    
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT2s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT2s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT3s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT3s_simVtx);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT5s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT5s_simVtx);  
    
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT2s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT2s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT3s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT3s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ05_dT5s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ05_dT5s);  
    
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT2s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT2s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT3s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT3s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ1_dT5s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ1_dT5s);  
    
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT2s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT2s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT3s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT3s);  
    eventTree[iRes]->Branch( Form("photon_chIso%.2d_dZ2_dT5s",int(isoConeDR_*10) ), &evInfo[iRes]->photon_chIso_dZ2_dT5s);  
    
    if (saveTracks_){
      eventTree[iRes]->Branch( "track_t",      &evInfo[iRes]->track_t);
      eventTree[iRes]->Branch( "track_dz4D",   &evInfo[iRes]->track_dz4D);
      eventTree[iRes]->Branch( "track_dz3D",   &evInfo[iRes]->track_dz3D);
      eventTree[iRes]->Branch( "track_dxy4D",  &evInfo[iRes]->track_dxy4D);
      eventTree[iRes]->Branch( "track_dxy3D",  &evInfo[iRes]->track_dxy3D);
      eventTree[iRes]->Branch( "track_pt",     &evInfo[iRes]->track_pt);
      eventTree[iRes]->Branch( "track_eta",    &evInfo[iRes]->track_eta);
      eventTree[iRes]->Branch( "track_phi",    &evInfo[iRes]->track_phi);
      eventTree[iRes]->Branch( "track_phoIndex",&evInfo[iRes]->track_phoIndex);
      eventTree[iRes]->Branch( "track_isMatchedToGenParticle",&evInfo[iRes]->track_isMatchedToGenParticle);
    }
  
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonIsolationAnalyzer::endJob() 
{
  cout << "Analysis finished " <<endl;
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
    evInfo[iRes]->npu = -1;
    evInfo[iRes]->vtxGen_t = -999;
    evInfo[iRes]->vtxGen_z = -999;
    evInfo[iRes]->vtx3D_z = -999;
    evInfo[iRes]->vtx3D_zErr = -999;
    evInfo[iRes]->vtx4D_z = -999;
    evInfo[iRes]->vtx4D_zErr = -999;
    evInfo[iRes]->vtx4D_t = -999;
    evInfo[iRes]->vtx4D_tErr = -999;
    evInfo[iRes]->vtx3D_isFake = -999;
    evInfo[iRes]->vtx4D_isFake = -999;

    evInfo[iRes]->photon_pt.clear();
    evInfo[iRes]->photon_eta.clear();
    evInfo[iRes]->photon_phi.clear();
    evInfo[iRes]->photon_t.clear();
    evInfo[iRes]->photon_isPrompt.clear();
    evInfo[iRes]->photon_isMatchedToGenJet.clear();
    evInfo[iRes]->photon_hasConversionTracks.clear();
    evInfo[iRes]->photon_hasPixelSeed.clear();
    evInfo[iRes]->photon_sigmaIetaIeta.clear();
    evInfo[iRes]->photon_r9.clear();


    evInfo[iRes]->photon_chIso_dZ05_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT2s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT3s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT5s_simVtx.clear();

    evInfo[iRes]->photon_chIso_dZ1_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT2s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT3s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT5s_simVtx.clear();
 
    evInfo[iRes]->photon_chIso_dZ2_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT2s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT3s_simVtx.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT5s_simVtx.clear();

    evInfo[iRes]->photon_chIso_dZ05.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT2s.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT3s.clear();
    evInfo[iRes]->photon_chIso_dZ05_dT5s.clear();

    evInfo[iRes]->photon_chIso_dZ1.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT2s.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT3s.clear();
    evInfo[iRes]->photon_chIso_dZ1_dT5s.clear();

    evInfo[iRes]->photon_chIso_dZ2.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT2s.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT3s.clear();
    evInfo[iRes]->photon_chIso_dZ2_dT5s.clear();

    if (saveTracks_){
      evInfo[iRes]->track_t.clear();
      evInfo[iRes]->track_dz4D.clear();
      evInfo[iRes]->track_dz3D.clear();
      evInfo[iRes]->track_dxy4D.clear();
      evInfo[iRes]->track_dxy3D.clear();
      evInfo[iRes]->track_pt.clear();
      evInfo[iRes]->track_eta.clear();
      evInfo[iRes]->track_phi.clear();
      evInfo[iRes]->track_phoIndex.clear();
      evInfo[iRes]->track_isMatchedToGenParticle.clear();
    }
  }
}


/*
// --- matching to gen photon
bool isPromptPhoton(const reco::Photon& photon, const edm::View<reco::GenParticle>& genParticles)
{
  bool isPrompt = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 22) continue;
    if (genp.status() != 1 || !genp.isLastCopy() ) continue; // -- from Simone
    if ( !genp.isPromptFinalState() ) continue;
    double dr = deltaR(photon,genp);
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
bool isMatchedToGenJet(const reco::Photon& photon, const edm::View<reco::GenJet>& genJets)
{
  bool isMatched = false;

  for(unsigned int ip=0; ip < genJets.size(); ip++ ){
    const reco::GenJet& genj = genJets[ip];
    if ( genj.pt() < 15.0  || genj.hadEnergy()/genj.energy() < 0.3) continue;
    double dr = deltaR(photon,genj);
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


*/

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIsolationAnalyzer);

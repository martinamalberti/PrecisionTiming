// -*- C++ -*-
//
// Package:    PrecisionTiming/PTAnalyzer
// Class:      PTAnalyzer
// 
/**\class PTAnalyzer PTAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/PTAnalyzer.cc

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

#include "PrecisionTiming/PTAnalysis/interface/PTAnalyzer.h"

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
PTAnalyzer::PTAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  vertex1DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex1DTag" ) ) ),
  vertex4DToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "Vertex4DTag" ) ) ),
  tracksToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag>( "TracksTag" ) ) ),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  pfCandToken_( consumes<View<reco::PFCandidate> >( iConfig.getParameter<InputTag>( "PFCandidateTag" ) ) ),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<InputTag>("genVtxTag")))
{
   //Now do what ever initialization is needed
  eventTree = fs_->make<TTree>( "event", "event" );
  keepMuons_ = iConfig.getUntrackedParameter<bool>("keepMuons");
  BTLEfficiency_ = iConfig.getUntrackedParameter<double>("BTLEfficiency");
  ETLEfficiency_ = iConfig.getUntrackedParameter<double>("ETLEfficiency");
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

  // -- get the 3D vertex collection with same pT cut (0.7 GeV) as 4D
  Handle<View<reco::Vertex> > Vertex1DCollectionH;
  iEvent.getByToken( vertex1DToken_, Vertex1DCollectionH );
  const edm::View<reco::Vertex>& vertices1D = *Vertex1DCollectionH;

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
  Handle<View<reco::Track> > TrackCollectionH;
  iEvent.getByToken(tracksToken_, TrackCollectionH);
  const edm::View<reco::Track>& tracks = *TrackCollectionH;

  // -- get the trackTimeValueMap
  Handle<ValueMap<float> > trackTimeValueMap;
  iEvent.getByToken( trackTimeToken_, trackTimeValueMap );


  // -- get the muons collection
  //  iEvent.getByToken(muonsToken_, muonsHandle_);

  // -- get the gen particles collection
  Handle<View<reco::GenParticle> > GenParticleCollectionH;
  iEvent.getByToken(genPartToken_, GenParticleCollectionH);
  const edm::View<reco::GenParticle>& genParticles = *GenParticleCollectionH;

  // -- get the gen vertex collection
  Handle<vector<SimVertex> > GenVertexCollectionH;
  iEvent.getByToken( genVertexToken_, GenVertexCollectionH );
  const vector<SimVertex>& genVertices = *GenVertexCollectionH;


  // -- get the PFCandidate collection
  //  Handle<View<reco::PFCandidate> > PFCandidateCollectionH;
  //iEvent.getByToken(pfCandToken_, PFCandidateCollectionH);
  //const edm::View<reco::PFCandidate>& pfCandidates = *PFCandidateCollectionH;

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

  //---get truth PV
  SimVertex genPV = genVertices.at(0);
  double mindz = 999999.;
  int vtx_index = 0;

  TRandom *gRandom = new TRandom();

  // -- reco vertices
  for(unsigned int ivtx=0; ivtx < vertices.size(); ivtx++ ){
    const reco::Vertex& vtx = vertices[ivtx];
    evInfo.vtx_z.push_back(vtx.z());

    // -- find teh vertex closest to gen vtx
    const float dz = std::abs(vtx.z() - genPV.position().z());
    if( dz < mindz )
      {
	mindz = dz;
	vtx_index = ivtx;
      }


    // -- count tracks multiplicity
    int nTks =0;
    int nTks_09_B = 0;
    int nTks_09_E = 0;
    int nTks_2_E = 0;

    bool isBTL;
    bool isETL;
    bool isHGCal;

    for(reco::Vertex::trackRef_iterator ti = vtx.tracks_begin(); ti!=vtx.tracks_end(); ++ti){

      // drop muon tracks
      const reco::Track & track=*(ti->get()); 
      bool isMuon = isMuonTrack(track, genParticles);
      if ( !keepMuons_ && isMuon )  { continue;}

      isBTL = ( (*ti) -> pt() > 0.9 && fabs( (*ti) -> eta() ) < 1.5 );
      isETL = ( (*ti) -> p()  > 0.9 && fabs( (*ti) -> eta() ) > 1.5 && fabs( (*ti) -> eta() ) < 3.0 );
      isHGCal = ( (*ti) -> pt() > 2.0 && fabs( (*ti) -> eta() ) > 1.5 && fabs( (*ti) -> eta() ) < 3.0 );

      // emulate BTL and ETL efficiency
      double rndB = gRandom->Uniform(0.,1.); 
      double rndE = gRandom->Uniform(0.,1.); 
      
      nTks++;
      if ( isBTL && rndB <= BTLEfficiency_) nTks_09_B++; 
      if ( isETL && rndE <= ETLEfficiency_) nTks_09_E++;
      if ( isHGCal ) nTks_2_E++;
      
    }

    evInfo.vtx_nTks.push_back(nTks);
    evInfo.vtx_nTks_pt09_B.push_back(nTks_09_B);
    evInfo.vtx_nTks_p09_E.push_back(nTks_09_E);
    evInfo.vtx_nTks_pt2_E.push_back(nTks_2_E);

  }

  evInfo.index_closestToGen = vtx_index;


  for(unsigned int ivtx=0; ivtx < vertices1D.size(); ivtx++ ){
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
    evInfo.tkTime.push_back( (*trackTimeValueMap)[tkRef] );

  }

  /*
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
  */

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

  eventTree->Branch( "npu",        &evInfo.npu);
  eventTree->Branch( "tkPt",       &evInfo.tkPt);
  eventTree->Branch( "tkEta",      &evInfo.tkEta);
  eventTree->Branch( "tkPhi",      &evInfo.tkPhi);
  eventTree->Branch( "tkTime",     &evInfo.tkTime);
  eventTree->Branch( "tkOuterR",   &evInfo.tkOuterR);
  eventTree->Branch( "tkOuterX",   &evInfo.tkOuterX);
  eventTree->Branch( "tkOuterY",   &evInfo.tkOuterY);
  eventTree->Branch( "tkOuterZ",   &evInfo.tkOuterZ);
  eventTree->Branch( "index_closestToGen",      &evInfo.index_closestToGen);
  eventTree->Branch( "vtx_z",      &evInfo.vtx_z);
  eventTree->Branch( "vtx_nTks",   &evInfo.vtx_nTks);
  eventTree->Branch( "vtx_nTks_pt09_B",  &evInfo.vtx_nTks_pt09_B);
  eventTree->Branch( "vtx_nTks_p09_E",   &evInfo.vtx_nTks_p09_E);
  eventTree->Branch( "vtx_nTks_pt2_E",   &evInfo.vtx_nTks_pt2_E);

  eventTree->Branch( "vtx1D_z",    &evInfo.vtx1D_z);
  eventTree->Branch( "vtx1D_nTks", &evInfo.vtx1D_nTks);
  eventTree->Branch( "vtx4D_z",    &evInfo.vtx4D_z);
  eventTree->Branch( "vtx4D_nTks", &evInfo.vtx4D_nTks);
  eventTree->Branch( "vtx4D_t",    &evInfo.vtx4D_t);
  eventTree->Branch( "pfPt",       &evInfo.pfPt);
  eventTree->Branch( "pfEta",      &evInfo.pfEta);
  eventTree->Branch( "pfPhi",      &evInfo.pfPhi);
  eventTree->Branch( "pfType",       &evInfo.pfType);

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
  evInfo.tkPt.clear();
  evInfo.tkEta.clear();
  evInfo.tkPhi.clear();
  evInfo.tkTime.clear();
  evInfo.tkOuterR.clear();
  evInfo.tkOuterX.clear();
  evInfo.tkOuterY.clear();
  evInfo.tkOuterZ.clear();
  evInfo.index_closestToGen = 0;
  evInfo.vtx_z.clear();
  evInfo.vtx_nTks.clear();
  evInfo.vtx_nTks_pt09_B.clear();
  evInfo.vtx_nTks_p09_E.clear();
  evInfo.vtx_nTks_pt2_E.clear();

  evInfo.vtx1D_z.clear();
  evInfo.vtx1D_nTks.clear();
  evInfo.vtx4D_z.clear();
  evInfo.vtx4D_nTks.clear();
  evInfo.vtx4D_t.clear();
  evInfo.pfPt.clear();
  evInfo.pfEta.clear();
  evInfo.pfPhi.clear();
  evInfo.pfType.clear();


}



bool isMuonTrack(const reco::Track & track, const edm::View<reco::GenParticle>& genParticles)
{
  bool isMuon = false;
  double mindr = 9999999.;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& p = genParticles[ip];
    if (p.status() != 1 || !p.isLastCopy() ) continue;
    if (std::abs(p.pdgId()) == 13)
      {      
	double dr = deltaR(track,p);
	//cout<<"dr="<<dr <<endl;
	if (dr < 0.02  && dr<mindr){
	  isMuon=true;
	}
      }
  }
  
  return isMuon;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PTAnalyzer);

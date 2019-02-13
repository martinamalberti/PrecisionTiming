// -*- C++ -*-%//
// Package:    PrecisionTiming/MuonNeutralIsolationAnalyzer
// Class:      MuonNeutralIsolationAnalyzer
// 
/**\class MuonNeutralIsolationAnalyzer MuonNeutralIsolationAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/MuonNeutralIsolationAnalyzer.cc

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
#include "PrecisionTiming/PTAnalysis/interface/MuonNeutralIsolationAnalyzer.h"

//#include <TMath.h>
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
MuonNeutralIsolationAnalyzer::MuonNeutralIsolationAnalyzer(const edm::ParameterSet& iConfig):
  PileUpToken_( consumes<vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken3D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag3D" ) ) ),
  vertexToken4D_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag4D" ) ) ),
  pfcandToken_( consumes<View<reco::PFCandidate> >( iConfig.getParameter<InputTag>( "PFCandidateTag" ) ) ),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<InputTag>("genVtxTag"))),
  //genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
  //genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))),      
  genJetsToken_(consumes<View<reco::GenJet> >(iConfig.getUntrackedParameter<InputTag>("genJetsTag"))),
  muonsToken_(consumes<View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonsTag"))),
  clustersBTLToken_(consumes<FTLClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("clustersBTLTag"))),
  clustersETLToken_(consumes<FTLClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("clustersETLTag")))
{
  timeResolutions_ = iConfig.getUntrackedParameter<vector<double> >("timeResolutions");
  mtd5sample_      = iConfig.getUntrackedParameter<bool>("mtd5sample");
  isoConeDR_       = iConfig.getUntrackedParameter<double>("isoConeDR");
  savePfCands_      = iConfig.getUntrackedParameter<bool>("savePfCands");
  maxDz_           = iConfig.getUntrackedParameter<double>("maxDz");
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


MuonNeutralIsolationAnalyzer::~MuonNeutralIsolationAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  //cout << "MuonNeutralIsolationAnalyzer destructor" <<endl;
  //for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++){
  //cout << iRes << "  " << evInfo[iRes]->npu <<endl;
  //cout << iRes << "  " << eventTree[iRes]->GetEntries()<<endl;
  ////  delete eventTree[iRes];
  ////delete evInfo[iRes];
  //}
  cout << "CIAO CIAO" << endl; 
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonNeutralIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //---get the MTD geometry
  edm::ESHandle<MTDGeometry> geoHandle;
  iSetup.get<MTDDigiGeometryRecord>().get(geoHandle);
  mtdGeometry_ = geoHandle.product();
  
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
  //const vector<SimVertex>& genVertices = *GenVertexCollectionH;

  // -- get the genXYZ
  //Handle<genXYZ> genXYZH;
  //iEvent.getByToken(genXYZToken_, genXYZH);
  
  // -- get the genT0
  //Handle<float> genT0H;
  //iEvent.getByToken(genT0Token_, genT0H);
    
  // -- get the gen jets collection
  Handle<View<reco::GenJet> > GenJetCollectionH;
  iEvent.getByToken(genJetsToken_, GenJetCollectionH);
  const edm::View<reco::GenJet>& genJets = *GenJetCollectionH;
  

  // -- get the MTD clusters collections
  Handle<FTLClusterCollection> clustersBTLHandle_;
  iEvent.getByToken(clustersBTLToken_, clustersBTLHandle_);
  auto clustersBTL = *clustersBTLHandle_.product();

  Handle<FTLClusterCollection> clustersETLHandle_;
  iEvent.getByToken(clustersETLToken_, clustersETLHandle_);
  auto clustersETL = *clustersETLHandle_.product();

  

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
  // else{
  //  auto xyz = genXYZH.product();
  //  auto t = *genT0H.product();
  //  auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
  //  genPV = SimVertex(v, t);
  // }
  
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

  int muonIndex =  -1;

  const int nResol = timeResolutions_.size();

  // --- start loop over muons
  for(unsigned int imu=0; imu < muons.size(); imu++ ){

    const reco::Muon& muon = muons[imu];

    // -- minimal checks
    if (muon.track().isNull()) continue;
    if (muon.pt() < 5.) continue;
    
    muonIndex++;

    // -- check if prompt or fake muon
    bool isPromptMu = isPromptMuon(muon,genParticles);
    bool isMatchedGenJet  = isMatchedToGenJet(muon, genJets);
    bool isFromTauDecay = isFromTau(muon,genParticles);

    // -- compute neutral isolations
    float neutrIso = 0.;
    float neutrIso_dT2s[nResol];
    float neutrIso_dT3s[nResol];
    float neutrIso_dT5s[nResol];
    float neutrIso_dT2s_simVtx[nResol];
    float neutrIso_dT3s_simVtx[nResol];
    float neutrIso_dT5s_simVtx[nResol];
    float time[nResol]; 
    
    // -- initialize 
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      neutrIso_dT2s[iRes] = 0 ;
      neutrIso_dT3s[iRes] = 0 ;
      neutrIso_dT5s[iRes] = 0 ;
      neutrIso_dT2s_simVtx[iRes] = 0 ;
      neutrIso_dT3s_simVtx[iRes] = 0 ;
      neutrIso_dT5s_simVtx[iRes] = 0 ;
      time[iRes] = 0.;
    }

    float muonTime =  -999.;
      
    // -- loop over neutral pf candidates
    for(unsigned icand = 0; icand < pfcands.size(); ++icand) {
      const reco::PFCandidate& pfcand = pfcands[icand];

      //if (pfcand.charge() != 0 ){
      //	auto pfcandRef = pfcands.refAt(icand);
      //	reco::TrackRef trackRef = pfcandRef->trackRef();
      //	if (trackRef.isNull()) continue;
      //	cout << icand << " pfcand.time = " << pfcand.time() << " trackRef.t0 = " << trackRef -> t0() <<endl; 
      //	cout << icand << " pfcand.timeError = " << pfcand.timeError() << " trackRef.t0Error = " << trackRef -> t0Error() <<endl; 
      //}

      if (pfcand.charge() != 0 ) continue;
      
      float dr  = deltaR(muon.eta(), muon.phi(), pfcand.eta(), pfcand.phi());

      // --- no timing 
      if (dr > minDr_ && dr < isoConeDR_){
	neutrIso+=pfcand.pt();
      }
      
      // --- with timing
      if (dr > minDr_ && dr < isoConeDR_){ 	  

	for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){                                                                                                    
	  float targetTimeResol = timeResolutions_[iRes];
	  float defaultTimeResol = 0.030;
	  if (mtd5sample_) defaultTimeResol = 0.035;
	  float extra_resol = 0.;
	  if ( targetTimeResol > defaultTimeResol) { extra_resol = sqrt(targetTimeResol*targetTimeResol - defaultTimeResol*defaultTimeResol); }

	  // check is it's matched to MTD cluster
	  bool matchedToMTD = false;
	  float clusterTime = -999.;
	  float clusterSeedTime = -999.;
	  float clusterX = -999.;
	  float clusterY = -999.;
	  float clusterZ = -999.;
	  float clusterR = -999.;
	  float dRcluster = -999.;
	  time[iRes] = -999.;
	  float dt = 0;
	  float dtsim = 0;
	  if (std::abs(pfcand.eta() < 1.5)) matchedToMTD = isMatchedToMTDCluster(pfcand, clustersBTL, mtdGeometry_, clusterTime, clusterSeedTime, dRcluster, clusterX, clusterY, clusterZ, clusterR);
	  if (std::abs(pfcand.eta() > 1.5)) matchedToMTD = isMatchedToMTDCluster(pfcand, clustersETL, mtdGeometry_, clusterTime, clusterSeedTime, dRcluster, clusterX, clusterY, clusterZ, clusterR);
	  if (matchedToMTD){
	    // -- extra smearing to emulate different time resolution
	    float rnd   = gRandom->Gaus(0., extra_resol);
	    time[iRes] = clusterTime + rnd;
	    float tof    = sqrt( pow(clusterZ - vtx4D.z(), 2)  + clusterR*clusterR ) / (TMath::C()*1.E-07) ; // m/s --> cm/ns
	    float tofsim = sqrt( pow(clusterZ - genPV.position().z(), 2) + clusterR*clusterR ) / (TMath::C()*1.E-07) ; // m/s --> cm/ns
	    dt = std::abs(time[iRes] - tof - vtx4D.t());
	    dtsim = std::abs(time[iRes] - tofsim - genPV.position().t()*1.E09);
	  }
	  else {
	    dt = 0. ;
	    dtsim = 0. ;
	  }
	  	   
	  // -- reco vertex
	  if (dt < 2.*targetTimeResol) { neutrIso_dT2s[iRes]+= pfcand.pt();}
	  if (dt < 3.*targetTimeResol) { neutrIso_dT3s[iRes]+= pfcand.pt();}
	  if (dt < 5.*targetTimeResol) { neutrIso_dT5s[iRes]+= pfcand.pt();}

	  // -- sim vertex
	  if (dtsim < 2.*targetTimeResol) { neutrIso_dT2s_simVtx[iRes]+= pfcand.pt();}
	  if (dtsim < 3.*targetTimeResol) { neutrIso_dT3s_simVtx[iRes]+= pfcand.pt();}
	  if (dtsim < 5.*targetTimeResol) { neutrIso_dT5s_simVtx[iRes]+= pfcand.pt();}
	  	  
	  // -- save info for neutral pf cands in the isolation cone 
	  if ( savePfCands_ )  {

	    bool genMatching  = isMatchedToGenParticle(pfcand, genParticles);
	    bool genUnmatching  = isUnmatchedToGenParticle(pfcand, genParticles);

	    evInfo[iRes]->neutrPfCand_particleId.push_back(pfcand.translateTypeToPdgId(pfcand.particleId())); 
	    evInfo[iRes]->neutrPfCand_cluster_t.push_back(clusterTime); 
	    evInfo[iRes]->neutrPfCand_clusterSeed_t.push_back(clusterSeedTime); 
	    evInfo[iRes]->neutrPfCand_cluster_x.push_back(clusterX); 
	    evInfo[iRes]->neutrPfCand_cluster_y.push_back(clusterY); 
	    evInfo[iRes]->neutrPfCand_cluster_z.push_back(clusterZ); 
	    evInfo[iRes]->neutrPfCand_cluster_R.push_back(clusterR); 
	    evInfo[iRes]->neutrPfCand_dRcluster.push_back(dRcluster); 
	    evInfo[iRes]->neutrPfCand_pt.push_back(pfcand.pt());
	    evInfo[iRes]->neutrPfCand_eta.push_back(pfcand.eta());
	    evInfo[iRes]->neutrPfCand_phi.push_back(pfcand.phi());
	    evInfo[iRes]->neutrPfCand_dRmu.push_back(dr);
	    evInfo[iRes]->neutrPfCand_muIndex.push_back(muonIndex);
	    evInfo[iRes]->neutrPfCand_isMatchedToGenParticle.push_back(genMatching);
	    evInfo[iRes]->neutrPfCand_isUnmatchedToGenParticle.push_back(genUnmatching);
	  }
	  
	}// end loop over time resolutions                                                                                                                                    
      }
    
    }// end loop over tracks
        

    // -- fill muon info for each resolution scenario
    for (unsigned int iRes = 0; iRes<timeResolutions_.size(); iRes++){
      evInfo[iRes]->muon_pt.push_back(muon.pt());  
      evInfo[iRes]->muon_eta.push_back(muon.eta());  
      evInfo[iRes]->muon_phi.push_back(muon.phi()); 
      evInfo[iRes]->muon_isLoose.push_back(muon::isLooseMuon(muon)); 
      evInfo[iRes]->muon_isMedium.push_back(muon::isMediumMuon(muon));
      int isTight3D = muon::isTightMuon(muon, vtx3D);
      int isTight4D = muon::isTightMuon(muon, vtx4D);
      evInfo[iRes]->muon_isTight3D.push_back(isTight3D); 
      evInfo[iRes]->muon_isTight4D.push_back(isTight4D); 
      evInfo[iRes]->muon_dz4D.push_back( muon.track()->dz(vtx4D.position()) );
      evInfo[iRes]->muon_dxy4D.push_back( muon.track()->dxy(vtx4D.position()) );
      evInfo[iRes]->muon_dz3D.push_back( muon.track()->dz(vtx3D.position()) );
      evInfo[iRes]->muon_dxy3D.push_back( muon.track()->dxy(vtx3D.position()) );
      evInfo[iRes]->muon_t.push_back(muonTime); 
      evInfo[iRes]->muon_isPrompt.push_back(isPromptMu);
      evInfo[iRes]->muon_isMatchedToGenJet.push_back(isMatchedGenJet);
      evInfo[iRes]->muon_isFromTauDecay.push_back(isFromTauDecay);
        
      evInfo[iRes]->muon_neutrIso.push_back(neutrIso);
      evInfo[iRes]->muon_neutrIso_dT2s.push_back(neutrIso_dT2s[iRes]);
      evInfo[iRes]->muon_neutrIso_dT3s.push_back(neutrIso_dT3s[iRes]);
      evInfo[iRes]->muon_neutrIso_dT5s.push_back(neutrIso_dT5s[iRes]);
      evInfo[iRes]->muon_neutrIso_dT2s_simVtx.push_back(neutrIso_dT2s_simVtx[iRes]);
      evInfo[iRes]->muon_neutrIso_dT3s_simVtx.push_back(neutrIso_dT3s_simVtx[iRes]);
      evInfo[iRes]->muon_neutrIso_dT5s_simVtx.push_back(neutrIso_dT5s_simVtx[iRes]);
    } // end loop over time resolutions

  }// end loop over muons


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
    evInfo[iRes]->vtxGen_t = genPV.position().t()*1.E09; // ns
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
MuonNeutralIsolationAnalyzer::beginJob()
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
    eventTree[iRes]->Branch( "muon_pt",         &evInfo[iRes]->muon_pt);  
    eventTree[iRes]->Branch( "muon_eta",        &evInfo[iRes]->muon_eta);  
    eventTree[iRes]->Branch( "muon_phi",        &evInfo[iRes]->muon_phi);  
    eventTree[iRes]->Branch( "muon_isLoose",    &evInfo[iRes]->muon_isLoose);  
    eventTree[iRes]->Branch( "muon_isMedium",   &evInfo[iRes]->muon_isMedium);  
    eventTree[iRes]->Branch( "muon_isTight3D",  &evInfo[iRes]->muon_isTight3D);  
    eventTree[iRes]->Branch( "muon_isTight4D",  &evInfo[iRes]->muon_isTight4D);  
    eventTree[iRes]->Branch( "muon_dz3D",       &evInfo[iRes]->muon_dz3D);  
    eventTree[iRes]->Branch( "muon_dz4D",       &evInfo[iRes]->muon_dz4D);  
    eventTree[iRes]->Branch( "muon_dxy3D",      &evInfo[iRes]->muon_dxy3D);  
    eventTree[iRes]->Branch( "muon_dxy4D",      &evInfo[iRes]->muon_dxy4D);  
    eventTree[iRes]->Branch( "muon_t",          &evInfo[iRes]->muon_t);  
    eventTree[iRes]->Branch( "muon_isPrompt",   &evInfo[iRes]->muon_isPrompt);  
    eventTree[iRes]->Branch( "muon_isMatchedToGenJet",&evInfo[iRes]->muon_isMatchedToGenJet);  
    eventTree[iRes]->Branch( "muon_isFromTauDecay",   &evInfo[iRes]->muon_isFromTauDecay);  

    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT2s",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT2s);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT3s",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT3s);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT5s",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT5s);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT2s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT2s_simVtx);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT3s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT3s_simVtx);  
    eventTree[iRes]->Branch( Form("muon_neutrIso%.2d_dT5s_simVtx",int(isoConeDR_*10) ), &evInfo[iRes]->muon_neutrIso_dT5s_simVtx);  
    
    if (savePfCands_){
      eventTree[iRes]->Branch( "neutrPfCand_particleId",    &evInfo[iRes]->neutrPfCand_particleId);
      eventTree[iRes]->Branch( "neutrPfCand_cluster_t",     &evInfo[iRes]->neutrPfCand_cluster_t);
      eventTree[iRes]->Branch( "neutrPfCand_clusterSeed_t", &evInfo[iRes]->neutrPfCand_clusterSeed_t);
      eventTree[iRes]->Branch( "neutrPfCand_cluster_x",     &evInfo[iRes]->neutrPfCand_cluster_x);
      eventTree[iRes]->Branch( "neutrPfCand_cluster_y",     &evInfo[iRes]->neutrPfCand_cluster_y);
      eventTree[iRes]->Branch( "neutrPfCand_cluster_z",     &evInfo[iRes]->neutrPfCand_cluster_z);
      eventTree[iRes]->Branch( "neutrPfCand_cluster_R",     &evInfo[iRes]->neutrPfCand_cluster_R);
      eventTree[iRes]->Branch( "neutrPfCand_dRcluster",     &evInfo[iRes]->neutrPfCand_dRcluster);
      eventTree[iRes]->Branch( "neutrPfCand_pt",     &evInfo[iRes]->neutrPfCand_pt);
      eventTree[iRes]->Branch( "neutrPfCand_eta",    &evInfo[iRes]->neutrPfCand_eta);
      eventTree[iRes]->Branch( "neutrPfCand_phi",    &evInfo[iRes]->neutrPfCand_phi);
      eventTree[iRes]->Branch( "neutrPfCand_dRmu",   &evInfo[iRes]->neutrPfCand_dRmu);
      eventTree[iRes]->Branch( "neutrPfCand_muIndex",&evInfo[iRes]->neutrPfCand_muIndex);
      eventTree[iRes]->Branch( "neutrPfCand_isMatchedToGenParticle",&evInfo[iRes]->neutrPfCand_isMatchedToGenParticle);
      eventTree[iRes]->Branch( "neutrPfCand_isUnmatchedToGenParticle",&evInfo[iRes]->neutrPfCand_isUnmatchedToGenParticle);
    }
  
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonNeutralIsolationAnalyzer::endJob() 
{
  cout << "Analysis finished " <<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonNeutralIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ method initialize tree structure ------------
void
MuonNeutralIsolationAnalyzer::initEventStructure()
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

    evInfo[iRes]->muon_pt.clear();
    evInfo[iRes]->muon_eta.clear();
    evInfo[iRes]->muon_phi.clear();
    evInfo[iRes]->muon_isLoose.clear();
    evInfo[iRes]->muon_isMedium.clear();
    evInfo[iRes]->muon_isTight3D.clear();    
    evInfo[iRes]->muon_isTight4D.clear();
    evInfo[iRes]->muon_dz3D.clear();
    evInfo[iRes]->muon_dxy3D.clear();
    evInfo[iRes]->muon_dz4D.clear();
    evInfo[iRes]->muon_dxy4D.clear();
    evInfo[iRes]->muon_t.clear();
    evInfo[iRes]->muon_isPrompt.clear();
    evInfo[iRes]->muon_isMatchedToGenJet.clear();
    evInfo[iRes]->muon_isFromTauDecay.clear();

    evInfo[iRes]->muon_neutrIso.clear();
    evInfo[iRes]->muon_neutrIso_dT2s.clear();
    evInfo[iRes]->muon_neutrIso_dT3s.clear();
    evInfo[iRes]->muon_neutrIso_dT5s.clear();
    evInfo[iRes]->muon_neutrIso_dT2s_simVtx.clear();
    evInfo[iRes]->muon_neutrIso_dT3s_simVtx.clear();
    evInfo[iRes]->muon_neutrIso_dT5s_simVtx.clear();


    if (savePfCands_){
      evInfo[iRes]->neutrPfCand_particleId.clear();
      evInfo[iRes]->neutrPfCand_cluster_t.clear();
      evInfo[iRes]->neutrPfCand_clusterSeed_t.clear();
      evInfo[iRes]->neutrPfCand_cluster_x.clear();
      evInfo[iRes]->neutrPfCand_cluster_y.clear();
      evInfo[iRes]->neutrPfCand_cluster_z.clear();
      evInfo[iRes]->neutrPfCand_cluster_R.clear();
      evInfo[iRes]->neutrPfCand_dRcluster.clear();
      evInfo[iRes]->neutrPfCand_pt.clear();
      evInfo[iRes]->neutrPfCand_eta.clear();
      evInfo[iRes]->neutrPfCand_phi.clear();
      evInfo[iRes]->neutrPfCand_dRmu.clear();
      evInfo[iRes]->neutrPfCand_muIndex.clear();
      evInfo[iRes]->neutrPfCand_isMatchedToGenParticle.clear();
      evInfo[iRes]->neutrPfCand_isUnmatchedToGenParticle.clear();
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonNeutralIsolationAnalyzer);

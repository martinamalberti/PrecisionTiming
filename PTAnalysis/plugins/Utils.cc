// user include files                                                                                                                                                                                                                  
#include "PrecisionTiming/PTAnalysis/interface/Utils.h"

// --- matching to gen muon ----------------------------------------------------------------
bool isPromptMuon(const reco::Muon& muon, const edm::View<reco::GenParticle>& genParticles)
{
  bool isPrompt = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 13) continue;
    if (genp.status() != 1 || !genp.isLastCopy() ) continue; // -- from Simone                                                                                                                                                          
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
//-------------------------------------------------------------------------------------------



// --- matching to gen jet ------------------------------------------------------------------
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
//-------------------------------------------------------------------------------------------



// --- matching to muons from tau decays ----------------------------------------------------
bool isFromTau(const reco::Muon& muon, const edm::View<reco::GenParticle>& genParticles)
{

  bool fromTau = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
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
//-------------------------------------------------------------------------------------------




//--- matching pfcands to gen level particles (from PV) ------------------------------------------------------- 
bool isMatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles){

  bool isMatched = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    // -- stable particle
    if (genp.status() !=1 ) continue;
    // -- charged or neutral particle
    if (genp.charge()!=pfcand.charge()) continue;
    // -- pt matching 
    float dpt = std::abs((pfcand.pt()-genp.pt())/genp.pt());
    if ( dpt > 0.1) continue;
    // -- deltaR matching
    double dr = deltaR(pfcand,genp);
    if (dr > 0.05){
      continue;
    }
    else{
      isMatched=true;
      break;
    }
  }

  return isMatched;

}
//-------------------------------------------------------------------------------------------



// -- not matching to gen particles (looser deltaR) -------------------------------------------------------------
bool isUnmatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles){

  bool isUnmatched = true;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    // -- stable particle
    if (genp.status() !=1 ) continue;
    // -- charged or neutral particle
    if (genp.charge()!=pfcand.charge()) continue;
    // -- not a muon
    if (genp.pdgId() == 13 ) continue;
    // -- pt matching
    //float dpt = std::abs((pfcand.pt()-genp.pt())/genp.pt());
    //if ( dpt > 0.2) continue; 
    // -- deltaR matching
    double dr = deltaR(pfcand,genp);
    if (dr > 0.15){
      continue;
    }
    else{
      isUnmatched=false;
      break;
    }
  }

  return isUnmatched;

}
//-------------------------------------------------------------------------------------------




// --- matching pfcand to MTD cluster -------------------------------------------------------
bool isMatchedToMTDCluster(const reco::PFCandidate &pfcand, FTLClusterCollection &clustersBTL, const MTDGeometry* mtdGeometry_, double &clusterTime,  double &clusterSeedTime){
  bool ismatched = false;

  for(auto clusIt : clustersBTL)
    {
      DetId id = clusIt.detId();
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      for(auto cluster : clusIt)
        {
          MeasurementPoint mp(cluster.x(),cluster.y());
          LocalPoint lp = topo.localPosition(mp);
          GlobalPoint gp = det->toGlobal(lp);

          float cluster_eta = gp.eta();
          float cluster_phi = gp.phi();

          float deta  = pfcand.eta() - cluster_eta;
          float dphi  = deltaPhi(pfcand.phi(), cluster_phi);
          float DR    = sqrt(deta*deta+dphi*dphi);

          if (DR<0.05){
            ismatched = true;
            clusterTime = cluster.time();
            clusterSeedTime = cluster.seed().time();
            break;
          }
        }
    }// end loop over clusters

  return ismatched;

}

//-------------------------------------------------------------------------------------------

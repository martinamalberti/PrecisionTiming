// user include files                                                                                                                                                                                                                  
#include "PrecisionTiming/PTAnalysis/interface/Utils.h"

// --- matching to gen muon ----------------------------------------------------------------
bool isPromptMuon(const reco::Muon& muon, const edm::View<reco::GenParticle>& genParticles, float & genpt)
{
  bool isPrompt = false;
  float drmin = 9999.;
  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genp.pdgId()) != 13) continue;
    if (genp.status() != 1 || !genp.isLastCopy() ) continue; // -- from Simone                                                                                                                                                          
    if ( !genp.isPromptFinalState() ) continue;
    double dr = deltaR(muon,genp);
    if (dr > 0.2){
      continue;
    }
    //else{
    //  isPrompt=true;
    //  break;
    //}
    if (dr < drmin && dr < 0.2){
      drmin = dr;
      genpt = genp.pt();
      isPrompt=true;   
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





#ifndef flashgg_HHWWggCandidate
#define flashgg_HHWWggCandidate

// https://root.cern.ch/doc/v608/TLorentzVector_8h_source.html
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

namespace flashgg {

  // HHWWggCandidate is a sub class or derived class of WeightedObject 
  class HHWWggCandidate : public WeightedObject
  {
  // access specifier 
  public:
    //---ctors---
    // when constructor overloading, each must have different number or specific types of input variables 
    HHWWggCandidate() ;
    // before adding jets
    //HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector);
    
    // After adding jets 
    //diphoVector_, phoVector, vertex_zero, genVertex, goodElectrons_, goodMuons_, theMET_, genParticlesVector, tagJets_, SLW_tag
    HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector, std::vector<flashgg::Jet> JetVector, bool SLW_tag, bool Pass_PS); 
    //HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector, std::vector<edm::Ptr<Jet>> tagJets); 
    // Testing with new constructor to make plot of variable from new item such as diphoton vector size 

    //---dtor---
    ~HHWWggCandidate();

    //---utils---
    const std::vector<flashgg::DiPhotonCandidate> diphoVector() const { return diphoVector_; };
    const std::vector<flashgg::Photon> phoVector() const { return phoVector_; };
    const std::vector<flashgg::Electron> electronVector() const {return electronVector_;} 
    const std::vector<flashgg::Muon> muonVector() const {return muonVector_;}
    const std::vector<flashgg::Met> METVector() const {return METVector_;}
    const std::vector<reco::GenParticle> GenParticlesVector() const {return GenParticlesVector_;}
    const std::vector<flashgg::Jet> JetVector() const {return JetVector_;}
    //const std::vector<edm::Ptr<Jet>> tagJets() const {return tagJets_;}
    //const reco::Candidate::LorentzVector& jet1() const { return jet1_; };
    const edm::Ptr<reco::Vertex> & vertex() const { return vertex_;  };
    const reco::GenParticle::Point & genVertex() const { return genVertex_;  };
    const std::vector<flashgg::Photon> phoP4Corrected() const { return phoP4Corrected_; };
    const float pho1_MVA() const { return pho1_MVA_; };
    const float pho2_MVA() const { return pho2_MVA_; };
    const float pho3_MVA() const { return pho3_MVA_; };
    const float pho4_MVA() const { return pho4_MVA_; };
    const reco::Candidate::LorentzVector& MET_fourvec() const { return MET_fourvec_; };
    const reco::Candidate::LorentzVector& leading_dpho() const { return leading_dpho_; };
    const reco::Candidate::LorentzVector& leading_pho() const { return leading_pho_; };
    const reco::Candidate::LorentzVector& leading_elec() const { return leading_elec_; };
    const reco::Candidate::LorentzVector& subleading_elec() const { return subleading_elec_; };
    const reco::Candidate::LorentzVector& leading_muon() const { return leading_muon_; };
    const reco::Candidate::LorentzVector& subleading_muon() const { return subleading_muon_; };
    const reco::Candidate::LorentzVector& muon1() const { return muon1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho1() const { return dp1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho2() const { return dp2_; };
    const reco::Candidate::LorentzVector& HHWWggFourVect() const { return tp_; };     
    const reco::Candidate::LorentzVector& MatchingDiJet() const { return mdij_; }; 
    const reco::Candidate::LorentzVector& NonMatchingDiJet() const { return nmdij_; }; 
    const reco::Candidate::LorentzVector& MatchingDiQuark() const { return mdiq_; }; 
    const reco::Candidate::LorentzVector& NonMatchingDiQuark() const { return nmdiq_; }; 
    //float getCosThetaStar_CS(float ebeam) const;
    //std::vector<float> CosThetaAngles() const;
    //float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
    const float theMETcorpt() const { return theMETcorpt_; };
    const float W1_TM() const { return W1_TM_; };
    const reco::Candidate::LorentzVector& gen_leading_elec() const { return gen_leading_elec_; };
    const reco::Candidate::LorentzVector& gen_subleading_elec() const { return gen_subleading_elec_; };
    const reco::Candidate::LorentzVector& gen_leading_muon() const { return gen_leading_muon_; };
    const reco::Candidate::LorentzVector& gen_subleading_muon() const { return gen_subleading_muon_; };
    bool test() const { return test_; };
    bool SLW_tag() const { return SLW_tag_; }; 
    bool Pass_PS() const { return Pass_PS_; }; 
    const reco::Candidate::LorentzVector& lsl_dij() const { return lsl_dij_; };
  private:

    std::vector<flashgg::DiPhotonCandidate> diphoVector_;
    std::vector<flashgg::Photon> phoVector_;
    edm::Ptr<reco::Vertex>               vertex_;
    reco::GenParticle::Point genVertex_;
    std::vector<flashgg::Electron> electronVector_;
    std::vector<flashgg::Muon> muonVector_;
    std::vector<flashgg::Met> METVector_;
    std::vector<reco::GenParticle> GenParticlesVector_;
    //std::vector<edm::Ptr<Jet>> tagJets_;
    std::vector<flashgg::Jet> JetVector_;
    std::vector<flashgg::Photon> phoP4Corrected_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
    reco::Candidate::LorentzVector MET_fourvec_;
    reco::Candidate::LorentzVector leading_dpho_;
    reco::Candidate::LorentzVector leading_pho_;
    reco::Candidate::LorentzVector leading_elec_;
    reco::Candidate::LorentzVector subleading_elec_;
    reco::Candidate::LorentzVector leading_muon_;
    reco::Candidate::LorentzVector subleading_muon_;
    reco::Candidate::LorentzVector muon1_;
    reco::Candidate::LorentzVector dp1_;
    reco::Candidate::LorentzVector dp2_;
    reco::Candidate::LorentzVector tp_;
    reco::Candidate::LorentzVector mdij_;
    reco::Candidate::LorentzVector nmdij_;
    reco::Candidate::LorentzVector mdiq_;
    reco::Candidate::LorentzVector nmdiq_;
    float theMETcorpt_;
    float W1_TM_;
    reco::Candidate::LorentzVector gen_leading_elec_;
    reco::Candidate::LorentzVector gen_subleading_elec_;
    reco::Candidate::LorentzVector gen_leading_muon_;
    reco::Candidate::LorentzVector gen_subleading_muon_;
    bool test_;
    bool SLW_tag_;
    bool Pass_PS_;
    reco::Candidate::LorentzVector lsl_dij_;

  };
  typedef std::vector<HHWWggCandidate> HHWWggCandidateCollection; // define new type: vector of HHWWggCandidates 

}

#endif

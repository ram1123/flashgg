#ifndef flashgg_HHWWggCandidate
#define flashgg_HHWWggCandidate

// https://root.cern.ch/doc/v608/TLorentzVector_8h_source.html
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
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
    //HHWWggCandidate( std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex);
    //HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET);
    //HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET);
    HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector);

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
    const edm::Ptr<reco::Vertex> & vertex() const { return vertex_;  };
    const reco::GenParticle::Point & genVertex() const { return genVertex_;  };
    const std::vector<flashgg::Photon> phoP4Corrected() const { return phoP4Corrected_; };
    const float pho1_MVA() const { return pho1_MVA_; };
    const float pho2_MVA() const { return pho2_MVA_; };
    const float pho3_MVA() const { return pho3_MVA_; };
    const float pho4_MVA() const { return pho4_MVA_; };
    const reco::Candidate::LorentzVector& MET_fourvec() const { return MET_fourvec_; };
    const reco::Candidate::LorentzVector& Abe_HHWWggDiPho() const { return abe_dp_; };
    const reco::Candidate::LorentzVector& elec1() const { return elec1_; };
    const reco::Candidate::LorentzVector& elec2() const { return elec2_; };
    const reco::Candidate::LorentzVector& muon1() const { return muon1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho1() const { return dp1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho2() const { return dp2_; };
    const reco::Candidate::LorentzVector& HHWWggFourVect() const { return tp_; };
    //float getCosThetaStar_CS(float ebeam) const;
    //std::vector<float> CosThetaAngles() const;
    //float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
    const float theMETcorpt() const { return theMETcorpt_; };
    const float W1_TM() const { return W1_TM_; };
    const float gen_lepton_pt() const { return gen_lepton_pt_; };
    const float gen_neutrino_pt() const {return gen_neutrino_pt_;};
    //const reco::Candidate::LorentzVector& gen() const {return gen_;};

  private:

    std::vector<flashgg::DiPhotonCandidate> diphoVector_;
    std::vector<flashgg::Photon> phoVector_;
    edm::Ptr<reco::Vertex>               vertex_;
    reco::GenParticle::Point genVertex_;
    std::vector<flashgg::Electron> electronVector_;
    std::vector<flashgg::Muon> muonVector_;
    std::vector<flashgg::Met> METVector_;
    std::vector<reco::GenParticle> GenParticlesVector_;
    std::vector<flashgg::Photon> phoP4Corrected_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
    reco::Candidate::LorentzVector MET_fourvec_;
    reco::Candidate::LorentzVector abe_dp_;
    reco::Candidate::LorentzVector elec1_;
    reco::Candidate::LorentzVector elec2_;
    reco::Candidate::LorentzVector muon1_;
    reco::Candidate::LorentzVector dp1_;
    reco::Candidate::LorentzVector dp2_;
    reco::Candidate::LorentzVector tp_;
    float theMETcorpt_;
    float W1_TM_;
    float gen_lepton_pt_;
    float gen_neutrino_pt_;
    //reco::Candidate::LorentzVector gen_;

  };
  typedef std::vector<HHWWggCandidate> HHWWggCandidateCollection; // define new type: vector of HHWWggCandidates 


}

#endif

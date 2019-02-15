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
#include "flashgg/DataFormats/interface/Met.h"

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
    HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Met> METVector);

    // Testing with new constructor to make plot of variable from new item such as diphoton vector size 

    //---dtor---
    ~HHWWggCandidate();

    //---utils---
    const std::vector<flashgg::DiPhotonCandidate> diphoVector() const { return diphoVector_; };
    const std::vector<flashgg::Photon> phoVector() const { return phoVector_; };
    const std::vector<flashgg::Electron> electronVector() const {return electronVector_;}
    const std::vector<flashgg::Met> METVector() const {return METVector_;}
    //const edm::Ptr<flashgg::Met> theMET() const {return theMET_;}
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
    const reco::Candidate::LorentzVector& HHWWggDiPho1() const { return dp1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho2() const { return dp2_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho1_Pho1() const { return dp1_pho1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho1_Pho2() const { return dp1_pho2_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho2_Pho1() const { return dp2_pho1_; };
    const reco::Candidate::LorentzVector& HHWWggDiPho2_Pho2() const { return dp2_pho2_; };
    const int& HHWWggDiPho1_iPho1() const { return dp1_ipho1_; };
    const int& HHWWggDiPho1_iPho2() const { return dp1_ipho2_; };
    const int& HHWWggDiPho2_iPho1() const { return dp2_ipho1_; };
    const int& HHWWggDiPho2_iPho2() const { return dp2_ipho2_; };
    const reco::Candidate::LorentzVector& HHWWggPho12() const { return pho12_; };
    const reco::Candidate::LorentzVector& HHWWggPho13() const { return pho13_; };
    const reco::Candidate::LorentzVector& HHWWggPho14() const { return pho14_; };
    const reco::Candidate::LorentzVector& HHWWggPho23() const { return pho23_; };
    const reco::Candidate::LorentzVector& HHWWggPho24() const { return pho24_; };
    const reco::Candidate::LorentzVector& HHWWggPho34() const { return pho34_; };
    const reco::Candidate::LorentzVector& HHWWggFourVect() const { return tp_; };
    float getCosThetaStar_CS(float ebeam) const;
    std::vector<float> CosThetaAngles() const;
    float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
    const float theMETcorpt() const { return theMETcorpt_; };
    const float W1_TM() const { return W1_TM_; };

  private:

    //std::vector<edm::Ptr<Electron> > Electrons_;
    //std::vector<flashgg::Electron> Electrons_;
    std::vector<flashgg::DiPhotonCandidate> diphoVector_;
    std::vector<flashgg::Photon> phoVector_;
    edm::Ptr<reco::Vertex>               vertex_;
    reco::GenParticle::Point genVertex_;
    std::vector<flashgg::Electron> electronVector_;
    std::vector<flashgg::Met> METVector_;
    //edm::Ptr<flashgg::Met> theMET_;
    //edm::Ptr<flashgg::Met> theMET_;
    std::vector<flashgg::Photon> phoP4Corrected_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
    reco::Candidate::LorentzVector MET_fourvec_;
    reco::Candidate::LorentzVector abe_dp_;
    reco::Candidate::LorentzVector elec1_;
    reco::Candidate::LorentzVector elec2_;
    reco::Candidate::LorentzVector dp1_;
    reco::Candidate::LorentzVector dp2_;
    reco::Candidate::LorentzVector dp1_pho1_;
    reco::Candidate::LorentzVector dp1_pho2_;
    reco::Candidate::LorentzVector dp2_pho1_;
    reco::Candidate::LorentzVector dp2_pho2_;
    int dp1_ipho1_;
    int dp1_ipho2_;
    int dp2_ipho1_;
    int dp2_ipho2_;
    reco::Candidate::LorentzVector pho12_;
    reco::Candidate::LorentzVector pho13_;
    reco::Candidate::LorentzVector pho14_;
    reco::Candidate::LorentzVector pho23_;
    reco::Candidate::LorentzVector pho24_;
    reco::Candidate::LorentzVector pho34_;
    reco::Candidate::LorentzVector tp_;
    float theMETcorpt_;
    float W1_TM_;

  };
  typedef std::vector<HHWWggCandidate> HHWWggCandidateCollection; // define new type: vector of HHWWggCandidates 


}

#endif

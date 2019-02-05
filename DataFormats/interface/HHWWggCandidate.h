#ifndef flashgg_HHWWggCandidate
#define flashgg_HHWWggCandidate

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

namespace flashgg {

  class HHWWggCandidate : public WeightedObject
  {
  public:
    //---ctors---
    HHWWggCandidate() ;
    HHWWggCandidate( std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex);

    //---dtor---
    ~HHWWggCandidate();

    //---utils---
    const std::vector<flashgg::Photon> phoVector() const { return phoVector_; };
    const edm::Ptr<reco::Vertex> & vertex() const { return vertex_;  };
    const reco::GenParticle::Point & genVertex() const { return genVertex_;  };
    const std::vector<flashgg::Photon> phoP4Corrected() const { return phoP4Corrected_; };
    const float pho1_MVA() const { return pho1_MVA_; };
    const float pho2_MVA() const { return pho2_MVA_; };
    const float pho3_MVA() const { return pho3_MVA_; };
    const float pho4_MVA() const { return pho4_MVA_; };
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

  private:

    std::vector<flashgg::Photon> phoVector_;
    edm::Ptr<reco::Vertex>               vertex_;
    reco::GenParticle::Point genVertex_;
    std::vector<flashgg::Photon> phoP4Corrected_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
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

  };
  typedef std::vector<HHWWggCandidate> HHWWggCandidateCollection;


}

#endif

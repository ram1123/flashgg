#ifndef flashgg_HHWWggTag
#define flashgg_HHWWggTag

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
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

namespace flashgg {

  // HHWWggTag is a sub class or derived class of WeightedObject 
  class HHWWggTag: public DiPhotonTagBase, public reco::LeafCandidate // inherit from tag producer  
  {
  // access specifier 
  public:
    //---ctors--- 
    // when constructor overloading, each must have different number or specific types of input variables 
    
    HHWWggTag();

    //---dtor---
    ~HHWWggTag();

    //-- SemiLeptonic final state constructors

    //- Without cut flow analysis and b tag variables 
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>); // HHWWggTag_0
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Muon>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>); // HHWWggTag_1

    //- With Cutflow analysis
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, std::vector<edm::Ptr<flashgg::Electron>>, std::vector<edm::Ptr<flashgg::Electron>>,
              std::vector<edm::Ptr<flashgg::Muon>>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, 
              std::vector<edm::Ptr<flashgg::Jet>>, std::vector<edm::Ptr<flashgg::Jet>>,
              std::vector<double>); 

    HHWWggTag(edm::Ptr<DiPhotonCandidate>, std::vector<edm::Ptr<flashgg::Electron>>, edm::Ptr<flashgg::Muon>, std::vector<edm::Ptr<flashgg::Muon>>, std::vector<edm::Ptr<flashgg::Muon>>,
              edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, 
              std::vector<edm::Ptr<flashgg::Jet>>, std::vector<edm::Ptr<flashgg::Jet>>,
              std::vector<double>); 

    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<double>); 
    // dipho, tag_electron, electrons, goodElectrons, theMET, jet1, jet2, Jets_, tagJets_
    // dipho, tag_muon, muons, goodMuons, theMET, jet1, jet2, Jets_, tagJets_

    //- With cut flow analysis and b tag variables 
    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<flashgg::Jet>, std::vector<double>); 
    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, std::vector<flashgg::Jet>, std::vector<double>); 
    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Muon>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<flashgg::Jet>, std::vector<double>); 
    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Muon>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, std::vector<flashgg::Jet>, std::vector<double>); 
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<double>); 
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, std::vector<double>); 
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Muon>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<double>); 
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Muon>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, std::vector<double>);     

    // Untagged, with cutflow analysis
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, 
              std::vector<edm::Ptr<flashgg::Electron>>, std::vector<edm::Ptr<flashgg::Electron>>,
              std::vector<edm::Ptr<flashgg::Muon>>, std::vector<edm::Ptr<flashgg::Muon>>,
              edm::Ptr<flashgg::Met>,
              std::vector<edm::Ptr<flashgg::Jet>>, std::vector<edm::Ptr<flashgg::Jet>>,
              std::vector<double>
              );
    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Electron>, edm::Ptr<flashgg::Met>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, std::vector<double>); 
    // dipho, electrons, goodElectrons, muons, goodMuons, theMET, Jets_, tagJets_

    // HHWWggTag(edm::Ptr<DiPhotonCandidate>, std::vector<flashgg::Jet>, edm::Ptr<flashgg::Met>, std::vector<double>);
    HHWWggTag(edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Met>, std::vector<double>);

    // no jets 

    virtual HHWWggTag *clone() const override; // You need this because HHWWggTag is derived from the DiPhotonTagBase

    // void GetPhotons(edm::Ptr<DiPhotonCandidate> dipho);
    void GetObjects(edm::Ptr<DiPhotonCandidate> dipho);
    void GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Met>);
    void GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET, edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2);
    void GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET, edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2);

    // Get Leptons and Jets 
    
    std::vector<flashgg::Electron> GetElectrons(std::vector<edm::Ptr<flashgg::Electron>>);
    std::vector<flashgg::Muon> GetMuons(std::vector<edm::Ptr<flashgg::Muon>>);
    std::vector<flashgg::Jet> GetJets(std::vector<edm::Ptr<flashgg::Jet>>);

    void setMVA(double x) { mva_ = x; }
    double MVA() const { return mva_; }
    void setEventNumber(double x) { eventNumber_ = x; }

    const flashgg::Photon* Leading_Photon() const { return Leading_Photon_; };
    const flashgg::Photon* Subleading_Photon() const { return Subleading_Photon_; };
    const float lp_Hgg_MVA() const {return lp_Hgg_MVA_;};
    const float slp_Hgg_MVA() const {return slp_Hgg_MVA_;};
    const flashgg::Electron Electron() const { return Electron_; };
    const flashgg::Muon Muon() const { return Muon_; };
    const flashgg::Jet Leading_Jet() const { return Leading_Jet_; };
    const flashgg::Jet Subleading_Jet() const { return Subleading_Jet_; };
    const flashgg::Met MET() const {return MET_; }; 
    //---utils---
    // const std::vector<flashgg::Jet> JetVector() const {return JetVector_;}
    const std::vector<double> Cut_Variables() const { return Cut_Variables_; };
    const std::vector<flashgg::Electron> allElectrons() const {return allElectrons_;}
    const std::vector<flashgg::Electron> goodElectrons() const {return goodElectrons_;}
    const std::vector<flashgg::Muon> allMuons() const {return allMuons_;}
    const std::vector<flashgg::Muon> goodMuons() const {return goodMuons_;}
    const std::vector<flashgg::Jet> allJets() const {return allJets_;}
    const std::vector<flashgg::Jet> goodJets() const {return goodJets_;}

    // const std::vector<flashgg::DiPhotonCandidate> diphoVector() const { return diphoVector_; };
    // const std::vector<flashgg::Electron> electronVector() const {return electronVector_;} 
    // const std::vector<flashgg::Muon> muonVector() const {return muonVector_;}
    // const std::vector<flashgg::Met> METVector() const {return METVector_;}
    // const std::vector<reco::GenParticle> GenParticlesVector() const {return GenParticlesVector_;}
    // const std::vector<reco::Candidate::LorentzVector> GenElectrons() const {return GenElectrons_;}
    // const std::vector<reco::Candidate::LorentzVector> GenMuons() const {return GenMuons_;}
    // const std::vector<reco::Candidate::LorentzVector> GenPhotons() const {return GenPhotons_;}
    // const std::vector<reco::Candidate::LorentzVector> GenQuarks() const {return GenQuarks_;}
    // const std::vector<reco::Candidate::LorentzVector> GenNeutrinos() const {return GenNeutrinos_;}
    // const std::vector<flashgg::Jet> JetVector() const {return JetVector_;}
    // const reco::Candidate::LorentzVector& MET_fourvec() const { return MET_fourvec_; };
    // const reco::Candidate::LorentzVector& leading_dpho() const { return leading_dpho_; };
    // const reco::Candidate::LorentzVector& leading_elec() const { return leading_elec_; };
    // const reco::Candidate::LorentzVector& subleading_elec() const { return subleading_elec_; };
    // const reco::Candidate::LorentzVector& leading_muon() const { return leading_muon_; };
    // const reco::Candidate::LorentzVector& subleading_muon() const { return subleading_muon_; };
    //float getCosThetaStar_CS(float ebeam) const;
    //std::vector<float> CosThetaAngles() const;
    //float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
    // const float theMETcorpt() const { return theMETcorpt_; };
    // const std::vector<double> Vertex_Variables() const { return Vertex_Variables_; };
    // const float dipho_MVA() const {return dipho_MVA_;};
    // const float CMS_hgg_mass() const {return CMS_hgg_mass_;};
    // const float dZ() const {return dZ_;};

  private:
    double mva_;
    long eventNumber_;
    const flashgg::Photon* Leading_Photon_;
    const flashgg::Photon* Subleading_Photon_;
    float lp_Hgg_MVA_;
    float slp_Hgg_MVA_;
    // std::vector<flashgg::Jet> JetVector_;
    std::vector<double> Cut_Variables_;
    flashgg::Electron Electron_;
    flashgg::Muon Muon_; 
    flashgg::Jet Leading_Jet_;
    flashgg::Jet Subleading_Jet_;
    flashgg::Met MET_;
    
    // std::vector<flashgg::DiPhotonCandidate> diphoVector_;
    std::vector<flashgg::Electron> allElectrons_;
    std::vector<flashgg::Electron> goodElectrons_;
    std::vector<flashgg::Muon> allMuons_;
    std::vector<flashgg::Muon> goodMuons_;
    std::vector<flashgg::Jet> allJets_;
    std::vector<flashgg::Jet> goodJets_;

    // std::vector<flashgg::Met> METVector_;
    // std::vector<reco::GenParticle> GenParticlesVector_;
    // std::vector<reco::Candidate::LorentzVector> GenElectrons_;
    // std::vector<reco::Candidate::LorentzVector> GenMuons_;
    // std::vector<reco::Candidate::LorentzVector> GenPhotons_;
    // std::vector<reco::Candidate::LorentzVector> GenQuarks_;
    // std::vector<reco::Candidate::LorentzVector> GenNeutrinos_;
    // reco::Candidate::LorentzVector MET_fourvec_;
    // reco::Candidate::LorentzVector leading_dpho_;
    // reco::Candidate::LorentzVector leading_elec_;
    // reco::Candidate::LorentzVector subleading_elec_;
    // reco::Candidate::LorentzVector leading_muon_;
    // reco::Candidate::LorentzVector subleading_muon_;
    // float theMETcorpt_;
    // std::vector<double> Vertex_Variables_;
    // float dipho_MVA_;
    // float CMS_hgg_mass_;
    // float dZ_;
    // float lp_Hgg_MVA_;
    // float slp_Hgg_MVA_;
    
  };

}

#endif

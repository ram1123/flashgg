// include path /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/src/
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/HHWWggTag.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace flashgg; // makes flashgg sub members visible 
HHWWggTag::HHWWggTag() : DiPhotonTagBase::DiPhotonTagBase(), mva_(-2.), JetVector_ (), Cut_Variables_ ()
{

}

HHWWggTag::~HHWWggTag() {}

void HHWWggTag::GetPhotons(edm::Ptr<DiPhotonCandidate> dipho)
{
  // dipho->makePhotonsPersistent();
  // auto diphop4 = dipho_.p4();

  // Get photons 
  // flashgg::Photon leading_photon = dipho->getLeadingPhoton();
  // flashgg::Photon subleading_photon = dipho->getSubLeadingPhoton();

  // Save as dumper objects 
  Leading_Photon_ = dipho->leadingPhoton();
  Subleading_Photon_ = dipho->subLeadingPhoton();

  lp_Hgg_MVA_ = Leading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() ); 
  slp_Hgg_MVA_ = Subleading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() ); 
}

// 2 jets, electron 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET,
                      edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  // flashgg::DiPhotonCandidate dipho_ = diphoVector_[0];
  // // if (!dipho_zero_vtx){
  // //   cout << "Diphoton vertex is not zero. Need to recompute" << endl; 
  // // }

  // leading_dpho_ = dipho;

  // // Get photons 

  dipho_ = dipho;
  GetPhotons(dipho);
}

// 1 jet, electron
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET,
                      edm::Ptr<flashgg::Jet> jet, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;  
  GetPhotons(dipho);
}

// two jets, muon 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,
                      edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;
  GetPhotons(dipho);
}

// 1 jet, muon 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,
                      edm::Ptr<flashgg::Jet> jet, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;
  GetPhotons(dipho);
}

// Untagged
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;
  GetPhotons(dipho);
}

//FullLep
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron1, edm::Ptr<flashgg::Electron> electron2, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables):Cut_Variables_(Cut_Variables)
{
   dipho_ = dipho;
   GetPhotons(dipho);
}

HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon1, edm::Ptr<flashgg::Muon> muon2, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables):Cut_Variables_(Cut_Variables)
{
   dipho_ = dipho;
   GetPhotons(dipho);
}

HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables):Cut_Variables_(Cut_Variables)
{
   dipho_ = dipho;
   GetPhotons(dipho);
}



HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho,std::vector<double> Cut_Variables):Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;                                    
  GetPhotons(dipho);
}
// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// You need this because HHWWggTag is derived from another class 
HHWWggTag *HHWWggTag::clone() const
{
    HHWWggTag *result = new HHWWggTag( *this );
    return result;
}


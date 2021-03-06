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
//---ctors---
// HHWWggTag::HHWWggTag() : DiPhotonTagBase::DiPhotonTagBase(), mva_(-2.), Cut_Variables_ (), MuonVars_ (), JetVars_ (), genMhh_(0.), genCosThetaStar_CS_(0.)
HHWWggTag::HHWWggTag() : DiPhotonTagBase::DiPhotonTagBase(), mva_(-2.), Cut_Variables_(), MuonVars_(), JetVars_(), genMhh_(0.), genCosThetaStar_CS_(0.), GenVtx_z_(-99), HggVtx_z_(-99), ZeroVtx_z_(-99)
{

}

//---dtor---
HHWWggTag::~HHWWggTag() {}

void HHWWggTag::GetPhoAtt(edm::Ptr<DiPhotonCandidate> dipho)
{
  Leading_Photon_ = dipho->leadingPhoton();
  Subleading_Photon_ = dipho->subLeadingPhoton();

  lp_Hgg_MVA_ = Leading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() );
  slp_Hgg_MVA_ = Subleading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() );
  lp_pt_ = Leading_Photon_->pt();
  slp_pt_ = Subleading_Photon_->pt();
}

void HHWWggTag::GetObjects(edm::Ptr<DiPhotonCandidate> dipho)
{
  GetPhoAtt(dipho);
}

void HHWWggTag::GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Met> MET)
{
  GetPhoAtt(dipho);
  MET_ = *MET;
}

void HHWWggTag::GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET, edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2)
{
  GetPhoAtt(dipho);
  Electron_ = *electron;
  MET_ = *MET;
  Leading_Jet_ = *jet1;
  Subleading_Jet_ = *jet2;

  dijet_ = (*jet1).p4() + (*jet2).p4();

  HWW_ = Leading_Jet_.p4() + Subleading_Jet_.p4() + Electron_.p4() + MET_.p4();
  HGG_ = dipho->p4();
  HH_ = Leading_Jet_.p4() + Subleading_Jet_.p4() + Electron_.p4() + MET_.p4() + dipho->p4();

}

void HHWWggTag::GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET, edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2)
{
  GetPhoAtt(dipho);
  Muon_ = *muon;
  MET_ = *MET;
  Leading_Jet_ = *jet1;
  Subleading_Jet_ = *jet2;

  dijet_ = (*jet1).p4() + (*jet2).p4();
  HWW_ = Leading_Jet_.p4() + Subleading_Jet_.p4() + Muon_.p4() + MET_.p4();
  HGG_ = dipho->p4();
  HH_ = Leading_Jet_.p4() + Subleading_Jet_.p4() + Muon_.p4() + MET_.p4() + dipho->p4();
}

void HHWWggTag::GetObjects(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Met> MET, edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2, edm::Ptr<flashgg::Jet> jet3, edm::Ptr<flashgg::Jet> jet4)
{
  GetPhoAtt(dipho);
  MET_ = *MET;
  Leading_Jet_ = *jet1;
  Subleading_Jet_ = *jet2;
  Sub2leading_Jet_ = *jet3;
  Sub3leading_Jet_ = *jet4;

  dijet_ = (*jet1).p4() + (*jet2).p4();
  dijet2_ = (*jet3).p4() + (*jet4).p4();

  HWW_ = dijet_ + dijet2_;
  HGG_ = dipho->p4();
  HH_ = dijet_ + dijet2_ + dipho->p4();
}

// Get Leptons and Jets 
std::vector<flashgg::Electron> HHWWggTag::GetElectrons(std::vector<edm::Ptr<flashgg::Electron>> electrons)
{
  std::vector<flashgg::Electron> savedElectrons;
  for (unsigned int i = 0; i < electrons.size(); i++){
    auto electron = electrons[i];
    flashgg::Electron * thisElectronPointer = const_cast<flashgg::Electron *>(electron.get());
    savedElectrons.push_back(*thisElectronPointer);
  }

  return savedElectrons;
}

void HHWWggTag::SetAllElectrons(std::vector<edm::Ptr<flashgg::Electron>> allElectrons){
  allElectrons_ = GetElectrons(allElectrons);
}

void HHWWggTag::SetGoodElectrons(std::vector<edm::Ptr<flashgg::Electron>> goodElectrons){
  goodElectrons_ = GetElectrons(goodElectrons);
}

void HHWWggTag::SetAllMuons(std::vector<edm::Ptr<flashgg::Muon>> allMuons){
  allMuons_ = GetMuons(allMuons); 
}

void HHWWggTag::SetGoodMuons(std::vector<edm::Ptr<flashgg::Muon>> goodMuons){
  goodMuons_ = GetMuons(goodMuons);
}

std::vector<flashgg::Muon> HHWWggTag::GetMuons(std::vector<edm::Ptr<flashgg::Muon>> muons)
{
  std::vector<flashgg::Muon> savedMuons;
  for (unsigned int i = 0; i < muons.size(); i++){
    auto muon = muons[i];
    flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(muon.get());
    savedMuons.push_back(*thisMuonPointer);
  }

  return savedMuons;
}

std::vector<flashgg::Jet> HHWWggTag::GetJets(std::vector<edm::Ptr<flashgg::Jet>> jets)
{
  std::vector<flashgg::Jet> savedJets;
  for (unsigned int i = 0; i < jets.size(); i++){
    auto jet = jets[i];
    flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(jet.get());
    savedJets.push_back(*thisJetPointer);
  }

  return savedJets;
}

// To call in tag producer without constructor 

// void HHWWggTag::setZeroVertex(edm::Ptr<reco::Vertex> ZeroVertex){
  // ZeroVertex_ = ZeroVertex;
// }

void HHWWggTag::SetGoodJets(std::vector<edm::Ptr<flashgg::Jet>> jets){
  goodJets_ = GetJets(jets);
}

void HHWWggTag::SetAllJets(std::vector<edm::Ptr<flashgg::Jet>> jets){
  allJets_ = GetJets(jets);
}

void HHWWggTag::SetDiPhoMVA(double dipho_MVA){
  dipho_MVA_ = dipho_MVA;
}

void HHWWggTag::SetDiPhoPt(double dipho_pt){
  dipho_pt_ = dipho_pt;
}

void HHWWggTag::SetDiphoCentralWeight(double DiphoCentralWeight){
  DiphoCentralWeight_ = DiphoCentralWeight;
}

void HHWWggTag::SetJetIDs(vector<vector<double>> passesIDs){
  goodJets_passJetPUID_ = passesIDs; 
}

void HHWWggTag::SetGenObjs(std::vector<edm::Ptr<reco::GenParticle> > genHiggsBosons, std::vector<edm::Ptr<reco::GenParticle> > genWBosons, std::vector<edm::Ptr<reco::GenParticle> > genPhotons,
                           std::vector<edm::Ptr<reco::GenParticle> > genQuarks, std::vector<edm::Ptr<reco::GenParticle> >  genLeptons,  std::vector<edm::Ptr<reco::GenParticle> > genNeutrinos){

  int pdgId; 

  for(unsigned int i = 0; i < genHiggsBosons.size(); i++){
    reco::Candidate::LorentzVector genFourVec = genHiggsBosons.at(i)->p4();
    genHiggsBosons_.push_back(genFourVec);
  }

  for(unsigned int i = 0; i < genWBosons.size(); i++){
    reco::Candidate::LorentzVector genFourVec = genWBosons.at(i)->p4();
    genWBosons_.push_back(genFourVec);    
  }

  for(unsigned int i = 0; i < genPhotons.size(); i++){
    reco::Candidate::LorentzVector genFourVec = genPhotons.at(i)->p4();
    genPhotons_.push_back(genFourVec);    
  }    

  for(unsigned int i = 0; i < genQuarks.size(); i++){
    pdgId = genQuarks.at(i)->pdgId();
    reco::Candidate::LorentzVector genFourVec = genQuarks.at(i)->p4();
    genQuarks_.push_back(genFourVec);    
    genQuarksPdgIds_.push_back(pdgId);    
  } 

  for(unsigned int i = 0; i < genLeptons.size(); i++){
    pdgId = genLeptons.at(i)->pdgId();
    reco::Candidate::LorentzVector genFourVec = genLeptons.at(i)->p4();
    genLeptons_.push_back(genFourVec);    
    genLeptonsPdgIds_.push_back(pdgId);    
  } 

  for(unsigned int i = 0; i < genNeutrinos.size(); i++){
    pdgId = genNeutrinos.at(i)->pdgId();
    reco::Candidate::LorentzVector genFourVec = genNeutrinos.at(i)->p4();
    genNeutrinos_.push_back(genFourVec);    
    genNeutrinosPdgIds_.push_back(pdgId);    
  }       

}

// void SavePDFInfo(GenEventInfoProduct::PDF const * pdf){

// }

//-- Fully Leptonic Leptons 

void HHWWggTag::GetFLElectrons(edm::Ptr<flashgg::Electron> Ele1,edm::Ptr<flashgg::Electron> Ele2)
{
  Leading_lepton_ = Ele1->p4();
  Subleading_lepton_ = Ele2->p4();
}

void HHWWggTag::GetFLMuons(edm::Ptr<flashgg::Muon> muon1,edm::Ptr<flashgg::Muon> muon2)
{
  Leading_lepton_ = muon1->p4();
  Subleading_lepton_ = muon2->p4();
}

void HHWWggTag::GetDiffLeptons(edm::Ptr<flashgg::Electron> Ele,edm::Ptr<flashgg::Muon> muon)
{
  if (Ele->p4().pt()>muon->p4().pt()){
    Leading_lepton_ = Ele->p4();
    Subleading_lepton_ = muon->p4();
  }
  else {
    Leading_lepton_ = muon->p4();
    Subleading_lepton_ = Ele->p4();
  }
}

void HHWWggTag::GetMET(edm::Ptr<flashgg::Met> MET) 
{
  MET_ = *MET; 
}

//-- Categories - No cutFlowAnalysis. Minimal variables only 

HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho)
{
  dipho_ = dipho; // need to save DiPhotonTagBase dipho_ object for CMS_hgg_mass 
  GetObjects(dipho);
}

//-- With cutFlowAnalysis. Save extra variables for final state objects 

// HHWWggTag_0 - Semileptonic Electron with cutFlowAnalysis
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron,
                     std::vector<edm::Ptr<flashgg::Electron>> allElectrons,
                     std::vector<edm::Ptr<flashgg::Electron>> goodElectrons,
                     std::vector<edm::Ptr<flashgg::Muon>> allMuons,
                     edm::Ptr<flashgg::Met> MET,
                     edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2,
                     std::vector<edm::Ptr<flashgg::Jet>> allJets, std::vector<edm::Ptr<flashgg::Jet>> goodJets,
                     std::vector<double> Cut_Variables, std::vector<double> MuonVars, std::vector<double> JetVars): Cut_Variables_(Cut_Variables), MuonVars_(MuonVars), JetVars_(JetVars)
{
  dipho_ = dipho;
  GetObjects(dipho, electron, MET, jet1, jet2);
  // GetObjects(dipho, electron, allElectrons, goodElectrons, MET, jet1, jet2);
  allElectrons_ = GetElectrons(allElectrons);
  goodElectrons_ = GetElectrons(goodElectrons);
  allMuons_ = GetMuons(allMuons);
  // goodMuons_ = GetMuons(goodMuons); // no good muons in this category
  allJets_ = GetJets(allJets);
  goodJets_ = GetJets(goodJets);
}

// HHWWggTag_1 - Semileptonic Muon with cutFlowAnalysis
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, std::vector<edm::Ptr<flashgg::Electron>> allElectrons,
                     edm::Ptr<flashgg::Muon> muon, std::vector<edm::Ptr<flashgg::Muon>> allMuons,
                     std::vector<edm::Ptr<flashgg::Muon>> goodMuons,
                     edm::Ptr<flashgg::Met> MET,
                     edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2,
                     std::vector<edm::Ptr<flashgg::Jet>> allJets, std::vector<edm::Ptr<flashgg::Jet>> goodJets,
                     std::vector<double> Cut_Variables, std::vector<double> MuonVars, std::vector<double> JetVars) : Cut_Variables_(Cut_Variables), MuonVars_(MuonVars), JetVars_(JetVars)
{
  dipho_ = dipho;
  GetObjects(dipho, muon, MET, jet1, jet2);
  allElectrons_ = GetElectrons(allElectrons);
  // goodElectrons_ = GetElectrons(goodElectrons); // no good electrons in this category
  allMuons_ = GetMuons(allMuons);
  goodMuons_ = GetMuons(goodMuons);
  allJets_ = GetJets(allJets);
  goodJets_ = GetJets(goodJets);
}

// HHWWggTag_1 - Fully Hadronic with cutFlowAnalysis
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho,
              edm::Ptr<flashgg::Met> MET,
              edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2,
              edm::Ptr<flashgg::Jet> jet3, edm::Ptr<flashgg::Jet> jet4,
              std::vector<edm::Ptr<flashgg::Jet>> allJets, std::vector<edm::Ptr<flashgg::Jet>> goodJets, 
              std::vector<double> Cut_Variables, std::vector<double> JetVars) : Cut_Variables_(Cut_Variables), JetVars_(JetVars)
{
  //std::cout << "DEBUG: Inside HHWWggtags class constructor for 4 jets only and other vars" << std::endl;
  dipho_ = dipho;
  GetObjects(dipho, MET, jet1, jet2, jet3, jet4);
  allJets_ = GetJets(allJets);
  goodJets_ = GetJets(goodJets);
}

// HHWWggTag_2 - Fully Leptonic, ee
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron1, edm::Ptr<flashgg::Electron> electron2, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables, double dipho_MVA):Cut_Variables_(Cut_Variables),dipho_MVA_(dipho_MVA)
{
   dipho_ = dipho;
   GetPhoAtt(dipho);
   GetFLElectrons(electron1,electron2);
   GetMET(MET);
   //std::cout<<Cut_Variables[0]<<std::endl;
}

// HHWWggTag_2 - Fully Leptonic, mumu
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon1, edm::Ptr<flashgg::Muon> muon2, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables, double dipho_MVA ):Cut_Variables_(Cut_Variables),dipho_MVA_(dipho_MVA)
{
   dipho_ = dipho;
   GetPhoAtt(dipho);
   GetFLMuons(muon1,muon2);
   GetMET(MET);
}

// HHWWggTag_2 - Fully Leptonic, emu 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,std::vector<double> Cut_Variables, double dipho_MVA ):Cut_Variables_(Cut_Variables), dipho_MVA_(dipho_MVA)
{
   dipho_ = dipho;
   GetPhoAtt(dipho);
   GetDiffLeptons(electron,muon);
   GetMET(MET);
}

//-- HHWWggTag_3 - Untagged with cutFlowAnalysis
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho,
                     std::vector<edm::Ptr<flashgg::Electron>> allElectrons, std::vector<edm::Ptr<flashgg::Electron>> goodElectrons,
                     std::vector<edm::Ptr<flashgg::Muon>> allMuons, std::vector<edm::Ptr<flashgg::Muon>> goodMuons,
                     edm::Ptr<flashgg::Met> MET,
                     std::vector<edm::Ptr<flashgg::Jet>> allJets, std::vector<edm::Ptr<flashgg::Jet>> goodJets,
                     std::vector<double> Cut_Variables, std::vector<double> MuonVars, std::vector<double> JetVars) : Cut_Variables_(Cut_Variables), MuonVars_(MuonVars), JetVars_(JetVars)
{
  dipho_ = dipho;
  GetObjects(dipho, MET);
  allElectrons_ = GetElectrons(allElectrons);
  goodElectrons_ = GetElectrons(goodElectrons);
  allMuons_ = GetMuons(allMuons);
  goodMuons_ = GetMuons(goodMuons);
  allJets_ = GetJets(allJets);
  goodJets_ = GetJets(goodJets);
}

//-- SL systematic tree. Good objects only 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho,
                     std::vector<edm::Ptr<flashgg::Electron>> goodElectrons,
                     std::vector<edm::Ptr<flashgg::Muon>> goodMuons,
                     edm::Ptr<flashgg::Met> MET,
                     std::vector<edm::Ptr<flashgg::Jet>> goodJets,
                     std::vector<double> Cut_Variables) : Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;
  GetObjects(dipho, MET);
  goodElectrons_ = GetElectrons(goodElectrons);
  goodMuons_ = GetMuons(goodMuons);
  goodJets_ = GetJets(goodJets);
}

//-- FH systematic tree 
HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho,
              edm::Ptr<flashgg::Met> MET,
              edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2,
              edm::Ptr<flashgg::Jet> jet3, edm::Ptr<flashgg::Jet> jet4,
              std::vector<edm::Ptr<flashgg::Jet>> goodJets, 
              std::vector<double> Cut_Variables) : Cut_Variables_(Cut_Variables)
{
  dipho_ = dipho;
  GetObjects(dipho, MET, jet1, jet2, jet3, jet4);
  goodJets_ = GetJets(goodJets); // Note that goodJets are pT ordered. jet1-4 above are ordered based on FH algorithm which is not necessarily pT ordered 
}

// Require this because HHWWggTag is derived from another class
HHWWggTag *HHWWggTag::clone() const
{
    HHWWggTag *result = new HHWWggTag( *this );
    return result;
}

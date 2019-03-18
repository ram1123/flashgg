#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/HHWWggCandidate.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"

using namespace flashgg; // makes flashgg sub members visible 
HHWWggCandidate::HHWWggCandidate():
diphoVector_ (),
phoVector_ (),
vertex_ (),
electronVector_ (),
muonVector_ (),
METVector_ (),
GenParticlesVector_ (),
//tagJets_ (),
JetVector_ (),
phoP4Corrected_ (),
pho1_MVA_ (),
pho2_MVA_ (),
pho3_MVA_ (),
pho4_MVA_ (),
MET_fourvec_ (),
abe_dp_(),
leading_elec_(),
subleading_elec_(),
leading_muon_(),
subleading_muon_(),
muon1_(),
dp1_ (),
dp2_ (),
tp_ (),
theMETcorpt_(),
W1_TM_ (),
leading_gen_elec_pt_ (),
subleading_gen_elec_pt_ (),
leading_gen_muon_pt_ (),
subleading_gen_muon_pt_ ()

 // Need absence of comma on last variable 

{}

  HHWWggCandidate::~HHWWggCandidate() {}

  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex):
  //phoVector_(phoVector), vertex_(vertex)

  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET):
  //diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), theMET_(theMET)

  // HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET):
  // diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), theMET_(theMET)

  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector):
  //diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), muonVector_(muonVector), METVector_(METVector), GenParticlesVector_(GenParticlesVector)

  HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector, std::vector<flashgg::Jet> JetVector):
  diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), muonVector_(muonVector), METVector_(METVector), GenParticlesVector_(GenParticlesVector), JetVector_(JetVector)
 
  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector, std::vector<edm::Ptr<Jet>> tagJets):
  //diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), muonVector_(muonVector), METVector_(METVector), GenParticlesVector_(GenParticlesVector), tagJets_(tagJets)

  {

    // void test_function(){

    //   cout << "Inside test function" << endl;
    // }

    // Call function finding leading and subleading momentum for each 

    // If there is only one diphoton candidate, plot its invariant mass
    // going to come from diphoVector[]
    // if size diphoVector == 1, save diphoton object as TLorentzvector to plot invariant mass 
    // Save Diphoton object from flashgg::DiPhotonCandidate 

    //bool one_PS_dpho = false; // at least one diphoton that passed preselection 
    //if (diphoVector_.size() > 0) one_PS_dpho = true;

    if (diphoVector_.size() == 1)
    {
      flashgg::DiPhotonCandidate dipho_ = diphoVector_[0];
      auto dipho = dipho_.p4();
      abe_dp_ = dipho;

    } 

    // MET 
    if (METVector_.size() == 1)
    {
      flashgg::Met met__ = METVector_[0];
      auto met_ = met__.p4();
      MET_fourvec_ = met_;

      //auto W = met_ + elec1;
      //W1_TM_ = W.Mt();

    } 

    // ---
    // Get leading and subleading leptons
    // --- 

    // Get leading electron 
    float l_elec_pt = -99;
    for (unsigned int i = 0; i < electronVector_.size(); i++ )
    {
      flashgg::Electron current_elec = electronVector_[i];
      auto current_elec_4vec = current_elec.p4();
      float current_elec_pt = current_elec_4vec.pt();

      // If current electron pt is greater than maximum so far, make it the leading pt electron 
      if (current_elec_pt > l_elec_pt){ 
        l_elec_pt = current_elec_pt; 
        auto leading_elec = current_elec_4vec;
        leading_elec_ = leading_elec;

      }
    } 

    // Get Subleading electron 
    float sl_elec_pt = -99;
    for (unsigned int i = 0; i < electronVector_.size(); i++ )
    {
      flashgg::Electron current_elec = electronVector_[i];
      auto current_elec_4vec = current_elec.p4();
      float current_elec_pt = current_elec_4vec.pt();

      // If current electron pt is greater than max subleading pt, and isn't the leading pt electron, make it the subleading electron 
      if ( (current_elec_pt > sl_elec_pt) && (current_elec_pt != l_elec_pt) ){ 
        sl_elec_pt = current_elec_pt; 
        auto subleading_elec = current_elec_4vec;
        subleading_elec_ = subleading_elec;

      }
    } 

    // Get leading muon
    float l_muon_pt = -99;

    //if (muonVector_.size() == 0 ) cout << "muonVector_.size() == 0" << endl;

    for (unsigned int i = 0; i < muonVector_.size(); i++ )
    {
      flashgg::Muon current_muon = muonVector_[i];
      auto current_muon_4vec = current_muon.p4();
      float current_muon_pt = current_muon_4vec.pt();

      // If current muon pt is greater than maximum so far, make it the leading pt muon
      if (current_muon_pt > l_muon_pt){ 
        l_muon_pt = current_muon_pt; 
        auto leading_muon = current_muon_4vec;
        leading_muon_ = leading_muon;

      }
    } 

    //if (l_muon_pt == 0.) l_muon_pt = -999;

    // Get Subleading muon
    float sl_muon_pt = -99;
    for (unsigned int i = 0; i < muonVector_.size(); i++ )
    {
      flashgg::Muon current_muon = muonVector_[i];
      auto current_muon_4vec = current_muon.p4();
      float current_muon_pt = current_muon_4vec.pt();

      // If current muon pt is greater than max subleading pt, and isn't the leading pt muon, make it the subleading muon 
      if ( (current_muon_pt > sl_muon_pt) && (current_muon_pt != l_muon_pt) ){ 
        sl_muon_pt = current_muon_pt; 
        auto subleading_muon = current_muon_4vec;
        subleading_muon_ = subleading_muon;

      }
    } 

    // if (electronVector_.size() == 2)
    // {
    //   flashgg::Electron firstelec = electronVector_[0];
    //   flashgg::Electron secondelec = electronVector_[1];
    //   auto leading_elec = firstelec.p4();
    //   auto subleading_elec = secondelec.p4();
      
    //   // if the first elec pt is higher, call it elec1 
    //   if (firstelec.pt() > secondelec.pt())
    //   {
    //     leading_elec = firstelec.p4();
    //     elec2 = secondelec.p4();
    //   }
    //   // if the secobnd elec pt is higher, call it elec1 
    //   else if (secondelec.pt() > firstelec.pt())
    //   {
    //     leading_elec = secondelec.p4();
    //     elec2 = firstelec.p4();
    //   }

    //   // Same pt
    //   else
    //   {
    //     leading_elec = firstelec.p4();
    //     elec2 = secondelec.p4();  
    //   }

    //   leading_elec_ = leading_elec; // leading pt electron
    //   subleading_elec_ = subleading_elec; // subleading pt electron 

    //   // MET 
    //   if (METVector_.size() == 1)
    //   {
    //     flashgg::Met met__ = METVector_[0];
    //     auto met_ = met__.p4();
    //     MET_fourvec_ = met_;

    //     auto W = met_ + elec1;
    //     W1_TM_ = W.Mt();

    //   } 

    // } 

    // // semileptonic signal 
    // else if (electronVector_.size() == 1)
    // {
    //   flashgg::Electron firstelec = electronVector_[0];
    //   auto elec1 = firstelec.p4();
    //   elec1_ = elec1;

    //   // MET 
    //   if (METVector_.size() == 1)
    //   {
    //     flashgg::Met met__ = METVector_[0];
    //     auto met_ = met__.p4();
    //     MET_fourvec_ = met_;

    //     auto W = met_ + elec1;
    //     W1_TM_ = W.Mt();

    //   } 


    // }

    // if (muonVector_.size() == 1)
    // {
    //   flashgg::Muon firstmuon = muonVector_[0];
    //   auto muon1 = firstmuon.p4();
    //   muon1_ = muon1;

    //   // MET 
    //   if (METVector_.size() == 1)
    //   {
    //     flashgg::Met met__ = METVector_[0];
    //     auto met_ = met__.p4();
    //     MET_fourvec_ = met_;

    //     auto W = met_ + muon1;
    //     W1_TM_ = W.Mt();

    //   } 


    // }

    // MET 

    //theMETcorpt_ = theMET_.getCorPt();

    //edm::Ptr<flashgg::Met> met = theMET;
    // auto dipho = dipho_.p4();
    // theMET_ = 

    // Save Higgs object (coming from fully leptonic, semi-leptonic, or fully hadronic)
    // First iteration: H->qqenu (where q's are specifially charm strange)
    
    // Jets 

    // Create all possible jet pairs 
    // Match to gen quarks 
    // Split into Wqq pair or not Wqq pairs 
    // Plot variables 

    // if (JetVector_.size() > 0)
    // {
    //   flashgg::Jet firstjet = JetVector_[0];
    //   auto jet1 = firstjet.p4();
    //   jet1_ = jet1;

    // }

    // GEN 

    //std::vector<int> quark_pdgids = {1,2,3,4,5}; // Include B quarks
    std::vector<int> quark_pdgids = {1,2,3,4};
    vector<reco::GenParticle> quarkVector;
    vector<reco::GenParticle> gen_electronVector;
    vector<reco::GenParticle> gen_muonVector;
    //vector<> 

    //cout << "GenParticlesVector_.size() = " << GenParticlesVector_.size() << endl; 
    //     float l_elec_pt = 0.;
    for (unsigned int i = 0; i < GenParticlesVector_.size(); i++){
    
      reco::GenParticle gen_ = GenParticlesVector_[i];  
      //cout << "i = " << i << endl;
      //cout << "gen_.pdgId() = " << gen_.pdgId() << endl;

      // Get gen_leading_electron_pt 

      if (abs(gen_.pdgId()) == 11)
      {
        gen_electronVector.push_back(gen_);

      } 

      if (abs(gen_.pdgId()) == 13)
      {
        gen_muonVector.push_back(gen_);
      } 

      // if (abs(gen_.pdgId()) == 11 || abs(gen_.pdgId()) == 13 ){
      //   //cout << "Found electron or muon in HHWWggCandidate.cc" << endl;
      //   //cout << "electronVector.size() = " << electronVector.size() << endl;
      //   //cout << "gen_.pdgId() = " << gen_.pdgId() << endl;
      //   //cout << "gen_.eta() = " << gen_.eta() << endl;
      //   //cout << "gen_.eta() = " << gen_.phi() << endl;
      //   //cout << "gen_.pt() = " << gen_.pt() << endl;
      //   //auto gen = gen_.p4();
      //   //cout << "typeid(gen_.pt()).name() = " << typeid(gen_.pt()).name() << endl;
      //   auto gen_lepton_pt_val = gen_.pt();
      //   auto gen_lepton_eta_val = gen_.eta();
      //   //cout << "gen_pt_val = " << gen_pt_val << endl;
      //   gen_lepton_pt_ = gen_lepton_pt_val;
      //   //gen_ = gen;
      // }

      if (abs(gen_.pdgId()) == 12 || abs(gen_.pdgId()) == 14 ){
        //cout << "Found electron neutrino or muon neutrino" << endl;
        //cout << "electronVector.size() = " << electronVector.size() << endl;
        //cout << "gen_.pdgId() = " << gen_.pdgId() << endl;
        //cout << "gen_.eta() = " << gen_.eta() << endl;
        //cout << "gen_.phi() = " << gen_.phi() << endl;
        //cout << "gen_.pt() = " << gen_.pt() << endl;
        //auto gen = gen_.p4();
        //cout << "typeid(gen_.pt()).name() = " << typeid(gen_.pt()).name() << endl;
        //auto gen_neutrino_pt = gen_.pt();
        //cout << "gen_neutrino_pt = " << gen_neutrino_pt << endl;
        //gen_neutrino_pt_ = gen_neutrino_pt;
        //gen_ = gen;
      }


      // Create quark vector (gen) 
      for (unsigned int i = 0; i < quark_pdgids.size(); i++){
        int val = quark_pdgids[i];

        // all gen_ quark objects in genparticles vector came from W bosons, as dictated by HHWWggCandidateProducer.cc 
        // Should break into W+ quarks and W- quarks so you know if two quarks came from the same W and can accurately call a pair a match 
        if (abs(gen_.pdgId()) == val ){
            cout << "quark mother = " << gen_.mother(0) << endl;
            quarkVector.push_back(gen_);


            //will want a vector of quarks 
            //if 2 quarks: order by leading, subleading  
            //Want to order by leading, subleading then save pt values 
            //std::cout << "found a quark from W boson" << endl;
            //reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen_.get());
            //genParticlesVector.push_back(*thisGENPointer);
        } 

      }

    }

    // Now that we have gen quark vector, use to match to jets 
    unsigned int JVSize = JetVector_.size(), QVSize = quarkVector.size();
    bool match1 = false, match2 = false;
    double dr = 0.;
    int mii = 0; 

    //vector<reco::GenParticle> qVecCopy = quarkVector;
    for (unsigned int i = 0; i < JVSize; i++){
      flashgg::Jet ijet = JetVector_[i];

      auto ij = ijet.p4(); 

      // Loop with all GEN jets to see if it has exactly one match 
      // if a quark matches a jet, remove it from the temporary vector 
      double mdr = 100.; // minimum dr 
      // Loop 
      //double dr = 0.,; 
      vector<reco::GenParticle> qVecCopy = quarkVector; // make new copy with all quarks in
      for (unsigned int ii = 0; ii < QVSize; ii++){
        dr = 0.;
        dr = sqrt( (qVecCopy[ii].eta()-ij.eta())*(qVecCopy[ii].eta()-ij.eta()) + (qVecCopy[ii].phi()-ij.phi())*(qVecCopy[ii].phi()-ij.phi()) );
        cout << "dr = " << dr << endl;
        if (dr < mdr) {
          mdr = dr;
          mii = ii; // index of minimum dr quark 
        }
      // Save whatever you need before deleting the matching quark 
      qVecCopy.erase(qVecCopy.begin() + mii);

      //if (match1) break; // if match found, stop looping quarks 
      }
      // remove matched quark from cloned vector 

      //quarkVector

      for (unsigned int j = i; j < JVSize; j++){
        flashgg::Jet jjet = JetVector_[j];

        auto jj = jjet.p4(); // get eta and phi from this 

        mdr = 100.;
        mii = 0;
        for (unsigned int jjj = 0; jjj < QVSize - 1; jjj++){
          dr = 0.;
          dr = sqrt( (qVecCopy[jjj].eta()-jj.eta())*(qVecCopy[jjj].eta()-jj.eta()) + (qVecCopy[jjj].phi()-jj.phi())*(qVecCopy[jjj].phi()-jj.phi()) );
          cout << "dr = " << dr << endl;
          if (dr < mdr) {
            mdr = dr;
            mii = jjj; // index of minimum dr quark 
          }
        // Save whatever you need before deleting the matching quark 
        qVecCopy.erase(qVecCopy.begin() + jjj);

        //if (match1) break; // if match found, stop looping quarks 
        }

      }

    }



    // Get leading and subleading GEN leptons 

    // Get gen electrons

    // Get leading pt gen_electron 
    float l_gen_elec_pt = -99;
    for (unsigned int i = 0; i < gen_electronVector.size(); i++ )
    {
      reco::GenParticle current_gen_elec = gen_electronVector[i];
      auto current_gen_elec_pt = current_gen_elec.pt();

      // If current gen electron pt is greater than maximum so far, make it the leading pt gen electron 
      if (current_gen_elec_pt > l_gen_elec_pt){ 
        l_gen_elec_pt = current_gen_elec_pt; 
        leading_gen_elec_pt_ = l_gen_elec_pt;

      }
    }     

    if (l_gen_elec_pt == -99) leading_gen_elec_pt_ = -99;

    // Get subleading pt gen_electron 
    float sl_gen_elec_pt = -99;
    for (unsigned int i = 0; i < gen_electronVector.size(); i++ )
    {
      reco::GenParticle current_gen_elec = gen_electronVector[i];
      auto current_gen_elec_pt = current_gen_elec.pt();

      // If current gen electron pt is greater than max subleading gen pt, and isn't the leading pt gen electron, make it the subleading gen electron 
      if ( (current_gen_elec_pt > sl_gen_elec_pt) && (current_gen_elec_pt != l_gen_elec_pt) ){ 
        sl_gen_elec_pt = current_gen_elec_pt; 
        subleading_gen_elec_pt_ = sl_gen_elec_pt;

      }
    } 

    if (sl_gen_elec_pt == -99) subleading_gen_elec_pt_ = -99;

    // Get gen muons 

    //cout << "gen_muonVector.size() = " << gen_muonVector.size() << endl;
    // Get leading pt gen_muon
    float l_gen_muon_pt = -99;
    for (unsigned int i = 0; i < gen_muonVector.size(); i++ )
    {
      reco::GenParticle current_gen_muon = gen_muonVector[i];
      auto current_gen_muon_pt = current_gen_muon.pt();

      // If current gen muon pt is greater than maximum so far, make it the leading pt gen muon
      if (current_gen_muon_pt > l_gen_muon_pt){ 
        l_gen_muon_pt = current_gen_muon_pt; 
        leading_gen_muon_pt_ = l_gen_muon_pt;
        //cout << "leading_gen_muon_pt_ = " << leading_gen_muon_pt_ << endl;

      }
    }     

    if (l_gen_muon_pt == -99) leading_gen_muon_pt_ = -99;

    // Get subleading pt gen_muon
    float sl_gen_muon_pt = -99;
    for (unsigned int i = 0; i < gen_muonVector.size(); i++ )
    {
      reco::GenParticle current_gen_muon = gen_muonVector[i];
      auto current_gen_muon_pt = current_gen_muon.pt();

      // If current gen muon pt is greater than max subleading gen pt, and isn't the leading pt gen muon, make it the subleading gen muon
      if ( (current_gen_muon_pt > sl_gen_muon_pt) && (current_gen_muon_pt != l_gen_muon_pt) ){ 
        sl_gen_muon_pt = current_gen_muon_pt; 
        subleading_gen_muon_pt_ = sl_gen_muon_pt;

      }
    } 

    //cout << "sl_gen_elec_pt = " << sl_gen_elec_pt << endl;
    if (sl_gen_elec_pt == -99) subleading_gen_muon_pt_ = -99;

    //cout << "quarkVector.size() = " << quarkVector.size() << endl;

    // if (quarkVector.size() == 2){
    //   //cout << ""
    //   // order particles by pt 

    // }

    // W transverse mass

    // W1_TM = 
    // W2_TM = 

    // (MET + elec1).Mt()

    // Eventually, Create Higgs object (from lep, neutrinos and MET) as TLorentzVector and plot invariant mass. 
    // Then Create di-Higgs object from first higgs and diphoton (Or should higgs object be made from diphoton?) and plot invariant mass. Radion/Graviton. 




  //   float vtx_X = vertex_->x();
  //   float vtx_Y = vertex_->y();
  //   float vtx_Z = vertex_->z();
  //   math::XYZVector vtx_Pos( vtx_X, vtx_Y, vtx_Z );
  //   if (phoVector_.size() > 0)
  //   {
  //   for( int p = 0; p < (int) phoVector_.size(); p++ )
  //   {
  //     float sc_X = phoVector_[p].superCluster()->x();
  //     float sc_Y = phoVector_[p].superCluster()->y();
  //     float sc_Z = phoVector_[p].superCluster()->z();
  //     math::XYZVector sc_Pos( sc_X, sc_Y, sc_Z );
  //     math::XYZVector direction = sc_Pos - vtx_Pos;
  //     math::XYZVector pho = ( direction.Unit() ) * ( phoVector_[p].energy() );
  //     math::XYZTLorentzVector corrected_p4( pho.x(), pho.y(), pho.z(), phoVector_[p].energy() );
  //     phoVector_[p].setP4(corrected_p4);
  //     phoP4Corrected_.push_back(phoVector_[p]);
  //   }
  // }


  //std::cout << "phoP4Corrected_.size() = " << phoP4Corrected_.size() << std::endl;
  //std::cout << "a" << std::endl;
  // pho1_MVA_ = phoP4Corrected_.size() > 0 ? phoP4Corrected_[0].phoIdMvaDWrtVtx(vertex_) : -999;
  // pho2_MVA_ = phoP4Corrected_.size() > 0 ? phoP4Corrected_[1].phoIdMvaDWrtVtx(vertex_) : -999;
  // pho3_MVA_ = phoP4Corrected_.size() > 2 ? phoP4Corrected_[2].phoIdMvaDWrtVtx(vertex_) : -999;
  // pho4_MVA_ = phoP4Corrected_.size() > 3 ? phoP4Corrected_[3].phoIdMvaDWrtVtx(vertex_) : -999;
  //std::cout << "b" << std::endl;
  
    //float minDM = 1000000;



    // if (phoP4Corrected_.size() > 3)
    // {
    //   for (int i1=0; i1 < (int) phoP4Corrected_.size(); i1++)
    //   {
    //     flashgg::Photon pho1 = phoP4Corrected_[i1];
    //     for (int i2=0; i2 < (int) phoP4Corrected_.size(); i2++)
    //     {
    //       if (i2 <= i1 ){continue;}
    //       flashgg::Photon pho2 = phoP4Corrected_[i2];
    //       for (int i3=0; i3 < (int) phoP4Corrected_.size(); i3++)
    //       {
    //         if (i3 == i2 || i3 == i1){continue;}
    //         flashgg::Photon pho3 = phoP4Corrected_[i3];
    //         for (int i4=0; i4 < (int) phoP4Corrected_.size(); i4++)
    //         {
    //           if (i4 <= i3){continue;}
    //           if (i4 == i1 || i4 == i2){continue;}
    //           flashgg::Photon pho4 = phoP4Corrected_[i4];
    //           auto dipho1 = pho1.p4() + pho2.p4();
    //           auto dipho2 = pho3.p4() + pho4.p4();
    //           float deltaM = fabs( dipho1.mass() - dipho2.mass());
    //           if (deltaM < minDM){
    //             minDM = deltaM;
    //             dp1_pho1_ = pho1.p4();
    //             dp1_ipho1_ = i1;
    //             dp1_pho2_ = pho2.p4();
    //             dp1_ipho2_ = i2;
    //             dp2_pho1_ = pho3.p4();
    //             dp2_ipho1_ = i3;
    //             dp2_pho2_ = pho4.p4();
    //             dp2_ipho2_ = i4;
    //             if ((pho1.pt() + pho2.pt()) > (pho3.pt() + pho4.pt()) )
    //             {
    //               dp1_ = dipho1;
    //               dp2_ = dipho2;
    //             }
    //             else if ((pho1.pt() + pho2.pt()) < (pho3.pt() + pho4.pt()) )
    //             {
    //               dp1_ = dipho2;
    //               dp2_ = dipho1;
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // }

    // if (phoP4Corrected_.size() == 2)
    // {
    //   tp_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4();
    //   pho12_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4(); // invariant mass of these two photons 
    // }
    // else if (phoP4Corrected_.size() == 3)
    // {
    //   tp_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4() + phoP4Corrected_[2].p4();
    //   pho12_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4();
    //   pho13_ = phoP4Corrected_[0].p4() + phoP4Corrected_[2].p4();
    //   pho23_ = phoP4Corrected_[1].p4() + phoP4Corrected_[2].p4();
    // }
    // else if (phoP4Corrected_.size() > 3 )
    // {
    //   tp_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4() + phoP4Corrected_[2].p4() + phoP4Corrected_[3].p4();
    //   pho12_ = phoP4Corrected_[0].p4() + phoP4Corrected_[1].p4();
    //   pho13_ = phoP4Corrected_[0].p4() + phoP4Corrected_[2].p4();
    //   pho14_ = phoP4Corrected_[0].p4() + phoP4Corrected_[3].p4();
    //   pho23_ = phoP4Corrected_[1].p4() + phoP4Corrected_[2].p4();
    //   pho24_ = phoP4Corrected_[1].p4() + phoP4Corrected_[3].p4();
    //   pho34_ = phoP4Corrected_[2].p4() + phoP4Corrected_[3].p4();
    // }

  }

  // float HHWWggCandidate::getCosThetaStar_CS(float ebeam) const {
  //   TLorentzVector p1, p2;
  //   p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
  //   p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

  //   reco::Candidate::LorentzVector h_lor = dp1_ + dp2_;
  //   TLorentzVector h;
  //   h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

  //   TVector3 boost = - h.BoostVector();
  //   p1.Boost(boost);
  //   p2.Boost(boost);
  //   reco::Candidate::LorentzVector a1_lor = dp1_;
  //   TLorentzVector a_1;
  //   a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;
  //   a_1.Boost(boost);

  //   TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
  //   CSaxis.Unit();

  //   return 0; 
  //   return cos(   CSaxis.Angle( a_1.Vect().Unit() )    );
  // // }
  // std::vector<float> HHWWggCandidate::CosThetaAngles() const {
  //   std::vector<float> helicityThetas;

  //   TLorentzVector Boosted_a1(0,0,0,0);
  //   Boosted_a1.SetPxPyPzE(dp1_.px(),dp1_.py(),dp1_.pz(),dp1_.energy()) ;
  //   TLorentzVector BoostedLeadingPhoton_a1(0,0,0,0);
  //   BoostedLeadingPhoton_a1.SetPxPyPzE(dp1_pho1_.px(),dp1_pho1_.py(),dp1_pho1_.pz(),dp1_pho1_.energy()) ;

  //   helicityThetas.push_back( HelicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1));

  //   TLorentzVector Boosted_a2(0,0,0,0);
  //   Boosted_a2.SetPxPyPzE(dp2_.px(),dp2_.py(),dp2_.pz(),dp2_.energy()) ;
  //   TLorentzVector BoostedLeadingPhoton_a2(0,0,0,0);
  //   BoostedLeadingPhoton_a2.SetPxPyPzE(dp2_pho1_.px(),dp2_pho1_.py(),dp2_pho1_.pz(),dp2_pho1_.energy()) ;

  //   helicityThetas.push_back( HelicityCosTheta(Boosted_a2, BoostedLeadingPhoton_a2));

  //   return helicityThetas;

  // }

  // float HHWWggCandidate::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
  // {
  //   TVector3 BoostVector = Booster.BoostVector();
  //   Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
  //   return Boosted.CosTheta();
  // }

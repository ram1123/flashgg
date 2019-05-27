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
leading_dpho_ (),
leading_pho_ (),
leading_elec_(),
subleading_elec_(),
leading_muon_(),
subleading_muon_(),
muon1_(),
dp1_ (),
dp2_ (),
tp_ (),
mdij_ (),
nmdij_ (),
mdiq_ (),
nmdiq_ (),
theMETcorpt_(),
W1_TM_ (),
gen_leading_elec_ (),
gen_subleading_elec_ (),
gen_leading_muon_ (),
gen_subleading_muon_ (), 
test_ (),
SLW_tag_ (),
Pass_PS_ (),
lsl_dij_ () // Need absence of comma on last variable 

{}

  HHWWggCandidate::~HHWWggCandidate() {}

  HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Muon> muonVector, std::vector<flashgg::Met> METVector, std::vector<reco::GenParticle> GenParticlesVector, std::vector<flashgg::Jet> JetVector, bool SLW_tag, bool Pass_PS):
  diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), muonVector_(muonVector), METVector_(METVector), GenParticlesVector_(GenParticlesVector), JetVector_(JetVector), SLW_tag_(SLW_tag), Pass_PS_(Pass_PS)
 
  {

    
    // if (SLW_tag_){
    // cout << "sizes:" << endl;
    // cout << "" << endl;
    // cout << "diphoVector_.size() = " << diphoVector_.size() << endl;
    // cout << "phoVector_.size() = " << phoVector_.size() << endl;
    // cout << "electronVector_.size() = " << electronVector_.size() << endl;
    // cout << "muonVector_.size() = " << muonVector_.size() << endl;
    // cout << "METVector_.size() = " << METVector_.size() << endl;
    // cout << "GenParticlesVector_.size() = " << GenParticlesVector_.size() << endl;
    // cout << "JetVector_.size() = " << JetVector_.size() << endl;
    // }
    
    // if (diphoVector_.size() >= 1) test = 1;
    // else test = 0;

    // Might eventually sort by object, not GEN/RECO 

    //-- GEN and Jets 

    //std::vector<int> quark_pdgids = {1,2,3,4,5}; // Include b quarks
    std::vector<int> quark_pdgids = {1,2,3,4}; // Don't look for b quarks 
    vector<reco::GenParticle> quarkVector; // quarks
    vector<reco::GenParticle> gen_electronVector; // gen electrons
    vector<reco::GenParticle> gen_muonVector; // gen muons 

    //cout << "GenParticlesVector_.size() = " << GenParticlesVector_.size() << endl; 

    // Look through gen particles vector 
    for (unsigned int i = 0; i < GenParticlesVector_.size(); i++){

      // cout << "Hello, in genparticlesvector lop" << endl;

    
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
            //cout << "quark mother = " << gen_->mother(0)->pdgId() << endl;
            //cout << "quark mother = " << gen_.mother(0).pdgId() << endl;
            auto m = gen_.mother(0);
            //int mid = m->pdgId();
            //cout << "motherID = " << m->pdgId() << endl;
            //quarkVector.push_back( make_pair( gen_,mid ) );
            quarkVector.push_back( gen_ );
        } 

      }

    }

    //-- Jet/Quark Matching

    // The options for number of gen quarks from H->WW is -> qqqq (4), ->qqlnu (2), and ->lnulnu (0)
    // 4,2 or 0 Quarks. 
    // 2, 1 or 0 diquark matches (from H->WW)
    // 
    // Number of jets can be many. 
    // Max number of dijet matched from H->WW: 2 (b/c 2, 1 or 0)

    // Now that we have gen quark vector, we use it to match to jets 
    unsigned int JVSize = JetVector_.size(), QVSize = quarkVector.size();
    //bool match1 = false, match2 = false;
    bool is_pair = 0;
    // if (diphoVector_.size() >= 1) {
    //   test = 1;
    //   }
    // else {
    //   test = 0;
    //   }
    double dr = 0., dpt = 0.;
    int mii = 0; 
    int qim = 0, qjm = 0; // i and j quark mothers 
    //float inv_mass = 0;

    std::vector<flashgg::Jet> jVecCopy = JetVector_; // make new copy with all jets in it so you can remove a jet after it's matched to a quark 

    //cout << "QVSize = " << QVSize << endl;
    //cout << "JVSize = " << JVSize << endl;

    // Before using gen info, can already form dijet objects based on pT
    // Remember, JetVector_ is ordered by pT
    // Save object of leading/subleading jets

    if (JVSize >= 2){
      auto j0 = JetVector_[0].p4();
      auto j1 = JetVector_[1].p4();
      auto lsl_dij = j0 + j1;
      lsl_dij_ = lsl_dij;
    }
    // If there are more quarks than jets, there will be quarks without matching jets 
    // If there are extra jets, we don't care. 

    // for some reason there are 4 quarks coming from the 2 quark signal file. For now want to skip these events. 
    if (QVSize == 2){

    if ( (QVSize < JVSize) && (JVSize > 0) && (QVSize > 0) ){

      //cout << "Hello, in quark/jet matching loop" << endl;

      //vector<reco::GenParticle> qVecCopy = quarkVector;
      // Loop 'i' quarks 
      for (unsigned int i = 0; i < QVSize - 1; i++){ // last quark cannot be ith quark because then there will be no j quark 

        cout << "Looking at quark " << i << endl;

        //flashgg::Jet ijet = JetVector_[i];
        reco::GenParticle quarki = quarkVector[i]; //quark i 
        auto qim0 = quarki.mother(0);
        qim = qim0->pdgId();
        //cout << "qim = " << qim << endl;

        auto qi = quarki.p4(); // quark i four vector 

        // Loop all jets to find smallest dr 
        // if a jet matches the quark, remove the jet from the temporary vector 
        double mdr = 100.; // minimum dr 
        // Loop 
        //double dr = 0.,; 
        // vector<reco::GenParticle> qVecCopy = quarkVector; // make new copy with all quarks in
        for (unsigned int ii = 0; ii < jVecCopy.size(); ii++){
          dr = 0.;
          //cout << "jVecCopy[ii].eta() = " << jVecCopy[ii].eta() << endl;
          //cout << "jVecCopy[ii].p4().eta() = " << jVecCopy[ii].p4().eta() << endl;

          dr = sqrt( (jVecCopy[ii].eta()-qi.eta())*(jVecCopy[ii].eta()-qi.eta()) + (jVecCopy[ii].phi()-qi.phi())*(jVecCopy[ii].phi()-qi.phi()) );
          //cout << "dr = " << dr << endl;
          if (dr < mdr) {
            mdr = dr;
            mii = ii; // index of minimum dr jet 
          }

        }

        // First jet matched 
        // Save whatever you need before deleting the matching quark from qVecCopy
        //cout << "Matched quark mother = " << qVecCopy[mii].mother(0)->pdgId() << endl;
        //cout << "Matched quark/jet dr = " << mdr << endl;
        //cout << "Matched jet index = " << mii << endl;
        auto jet1 = jVecCopy[mii].p4();
          // for (unsigned int iii = 0; iii < jVecCopy.size(); iii++){
          //   cout << "jVecCopy[" << iii << "] = " << jVecCopy[iii] << endl;
          // }
        jVecCopy.erase(jVecCopy.begin() + mii); // remove matched jet from cloned vector 

          // for (unsigned int iii = 0; iii < jVecCopy.size(); iii++){
          //   cout << "jVecCopy[" << iii << "] = " << jVecCopy[iii] << endl;
          // }
        // use same cloned vector to match remaining jets to second quark

        // loop 'j' quarks. The ones remaining after 'i'
        //for (unsigned int j = i; j < JVSize; j++){
        for (unsigned int j = i + 1; j < QVSize; j++){
          reco::GenParticle quarkj = quarkVector[j]; //quark j
          auto qjm0 = quarkj.mother(0);
          qjm = qjm0->pdgId();
          //cout << "qjm = " << qjm << endl;
          auto qj = quarkj.p4(); // quark j four vector 
          mdr = 100.;
          mii = 0;
          for (unsigned int jjj = 0; jjj < jVecCopy.size(); jjj++){ //jVecCopy size should be different 
            dr = 0.;
            dr = sqrt( (jVecCopy[jjj].eta()-qj.eta())*(jVecCopy[jjj].eta()-qj.eta()) + (jVecCopy[jjj].phi()-qj.phi())*(jVecCopy[jjj].phi()-qj.phi()) );
            //cout << "dr = " << dr << endl;
            if (dr < mdr) {
              mdr = dr;
              mii = jjj; // index of minimum dr quark 
            }

          }
          // Save whatever you need before deleting the matching quark 
          //cout << "Matched quark/jet dr = " << mdr << endl;
          //cout << "Matched jet index = " << mii << endl;
          auto jet2 = jVecCopy[mii].p4();

          // for (unsigned int iii = 0; iii < jVecCopy.size(); iii++){
          //   cout << "jVecCopy[" << iii << "] = " << jVecCopy[iii] << endl;
          // }

          jVecCopy.erase(jVecCopy.begin() + mii); // remove matched jet from cloned vector 

          // for (unsigned int iii = 0; iii < jVecCopy.size(); iii++){
          //   cout << "jVecCopy[" << iii << "] = " << jVecCopy[iii] << endl;
          // }

          // After each second (or 'j') quark is matched, you have a dijet pair (i jet and j jet)
          // which are matched to quarks based on dr.
          // The mothers of the quarks tell you if the quarks come from the same W or not.

          // If quarks have same mother particle type *** This may need to be adjusted when there is pileup or backgrounds
          // because in that case, one quark may come from a signal W and the other from a background/pileup W 
          if (qim == qjm) is_pair = 1; // only have this info when checking quarks 
          else is_pair = 0;
        
          // Might need to add more matching and non matching jet objects when there are 4 quarks 
          // This is because only one of each object is saved per event 

          if (is_pair){
          // matching dijet object 
          auto mdij = jet1 + jet2; // does 'auto' use more memory than giving the proper type? Might want to change this to speed things up in the future. 
          mdij_ = mdij;

          // matching diquark object 
          auto mdiq = qi + qj;
          mdiq_ = mdiq;
          }

          if (!is_pair){
          // non-matching diquark object 
          auto nmdiq = qi + qj;
          nmdiq_ = nmdiq;
          
          // non-matching dijet object 
          auto nmdij = jet1 + jet2;
          nmdij_ = nmdij;
          }

        }

      }

    } // Loop performed if there are more jets than quarks 

    }

    else cout << "There are more quarks than jets, and/or no jets at all. Skipping." << endl;

    //cout << "Hello, outside of matching loop" << endl;

    // Get leading and subleading GEN leptons 

    // Get gen electrons

    // Get leading pt gen_electron 
    float l_gen_elec_pt = -99;
    for (unsigned int i = 0; i < gen_electronVector.size(); i++ )
    {
      reco::GenParticle current_gen_elec = gen_electronVector[i];
      auto current_gen_elec_4vec = current_gen_elec.p4();
      auto current_gen_elec_pt = current_gen_elec.pt();

      // If current gen electron pt is greater than maximum so far, make it the leading pt gen electron 
      if (current_gen_elec_pt > l_gen_elec_pt){ 
        l_gen_elec_pt = current_gen_elec_pt; 
        gen_leading_elec_ = current_gen_elec_4vec;

      }
    }     

    //if (l_gen_elec_pt == -99) leading_gen_elec_pt_ = -99;

    // Get subleading pt gen_electron 
    float sl_gen_elec_pt = -99;
    for (unsigned int i = 0; i < gen_electronVector.size(); i++ )
    {
      reco::GenParticle current_gen_elec = gen_electronVector[i];
      auto current_gen_elec_4vec = current_gen_elec.p4();
      auto current_gen_elec_pt = current_gen_elec.pt();

      // If current gen electron pt is greater than max subleading gen pt, and isn't the leading pt gen electron, make it the subleading gen electron 
      if ( (current_gen_elec_pt > sl_gen_elec_pt) && (current_gen_elec_pt != l_gen_elec_pt) ){ 
        sl_gen_elec_pt = current_gen_elec_pt; 
        gen_subleading_elec_ = current_gen_elec_4vec;

      }
    } 

    //if (sl_gen_elec_pt == -99) subleading_gen_elec_pt_ = -99;

    // Get gen muons 

    //cout << "gen_muonVector.size() = " << gen_muonVector.size() << endl;
    // Get leading pt gen_muon
    float l_gen_muon_pt = -99;
    for (unsigned int i = 0; i < gen_muonVector.size(); i++ )
    {
      reco::GenParticle current_gen_muon = gen_muonVector[i];
      auto current_gen_muon_4vec = current_gen_muon.p4();
      auto current_gen_muon_pt = current_gen_muon.pt();

      // If current gen muon pt is greater than maximum so far, make it the leading pt gen muon
      if (current_gen_muon_pt > l_gen_muon_pt){ 
        //l_gen_muon_pt = current_gen_muon_pt; 
        //leading_gen_muon_pt_ = l_gen_muon_pt;
        l_gen_muon_pt = current_gen_muon_pt; 
        gen_leading_muon_ = current_gen_muon_4vec;

      }
    }     

    //if (l_gen_muon_pt == -99) leading_gen_muon_pt_ = -99;

    // Get subleading pt gen_muon
    float sl_gen_muon_pt = -99;
    for (unsigned int i = 0; i < gen_muonVector.size(); i++ )
    {
      reco::GenParticle current_gen_muon = gen_muonVector[i];
      auto current_gen_muon_pt = current_gen_muon.pt();
      auto current_gen_muon_4vec = current_gen_muon.p4();

      // If current gen muon pt is greater than max subleading gen pt, and isn't the leading pt gen muon, make it the subleading gen muon
      if ( (current_gen_muon_pt > sl_gen_muon_pt) && (current_gen_muon_pt != l_gen_muon_pt) ){ 
        //sl_gen_muon_pt = current_gen_muon_pt; 
        //subleading_gen_muon_pt_ = sl_gen_muon_pt;
        sl_gen_muon_pt = current_gen_muon_pt; 
        gen_subleading_muon_ = current_gen_muon_4vec;

      }
    } 

    //-- RECO (besides jets) 

    //-- Diphotons 
    unsigned int ndpho = diphoVector_.size(); // number of diphotons 
    double tmp_dp_pt = 0, max_dp_pt = -99; // temporary diphoton pt 
    //bool test = 0;

    // If only one diphoton, save its four vector 
    if (ndpho == 1)
    {
      test_ = 1; 
      flashgg::DiPhotonCandidate dipho_ = diphoVector_[0];
      auto dipho = dipho_.p4();
      leading_dpho_ = dipho;

    } 

    // If more than one diphoton, take the highest pt diphoton 
    else if (ndpho > 1)
    {
      test_ = 1;
      
      for (unsigned int i = 0; i < ndpho; i ++)
      {
        flashgg::DiPhotonCandidate tmp_dipho_ = diphoVector_[i];
        auto dipho_ = tmp_dipho_.p4();
        tmp_dp_pt = dipho_.pt();
        if (tmp_dp_pt > max_dp_pt) 
        {
          max_dp_pt = tmp_dp_pt;
          leading_dpho_ = dipho_;
        }
      }

    }

    //-- Photons

      unsigned int npho = phoVector_.size(); // number of photons 
      double tmp_p_pt = 0, max_p_pt = -99; // temporary photon pt 

      
      for (unsigned int i = 0; i < npho; i ++)
      {
        flashgg::Photon tmp_pho_ = phoVector_[i];
        auto pho_ = tmp_pho_.p4();
        tmp_p_pt = pho_.pt();
        if (tmp_p_pt > max_p_pt) 
        {
          max_p_pt = tmp_p_pt;
          leading_pho_ = pho_;
        }
      }

    //cout << "Right before MET" << endl;

    // MET 
    if (METVector_.size() == 1)
    {
      flashgg::Met met__ = METVector_[0];
      auto met_ = met__.p4();
      MET_fourvec_ = met_;

      //auto W = met_ + elec1;
      //W1_TM_ = W.Mt();

    } 

    //cout << "Right after MET" << endl;

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

  }

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/HHWWggCandidate.h"
#include "flashgg/DataFormats/interface/Met.h"

using namespace flashgg; // makes flashgg sub members visible 
HHWWggCandidate::HHWWggCandidate():
diphoVector_ (),
phoVector_ (),
vertex_ (),
electronVector_ (),
METVector_ (),
//theMET_ (),
phoP4Corrected_ (),
pho1_MVA_ (),
pho2_MVA_ (),
pho3_MVA_ (),
pho4_MVA_ (),
MET_fourvec_ (),
abe_dp_(),
elec1_(),
elec2_(),
dp1_ (),
dp2_ (),
dp1_pho1_ (),
dp1_pho2_ (),
dp2_pho1_ (),
dp2_pho2_ (),
dp1_ipho1_ (),
dp1_ipho2_ (),
dp2_ipho1_ (),
dp2_ipho2_ (),
pho12_ (),
pho13_ (),
pho14_ (),
pho23_ (),
pho24_ (),
pho34_ (),
tp_ (),
theMETcorpt_(),
W1_TM_ ()
{}

  HHWWggCandidate::~HHWWggCandidate() {}

  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex):
  //phoVector_(phoVector), vertex_(vertex)

  //HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET):
  //diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), theMET_(theMET)

  // HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, edm::Ptr<flashgg::Met> theMET):
  // diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), theMET_(theMET)

  HHWWggCandidate::HHWWggCandidate( std::vector<flashgg::DiPhotonCandidate> diphoVector, std::vector<flashgg::Photon> phoVector, edm::Ptr<reco::Vertex> vertex, reco::GenParticle::Point genVertex, std::vector<flashgg::Electron> electronVector, std::vector<flashgg::Met> METVector):
  diphoVector_(diphoVector), phoVector_(phoVector), vertex_(vertex), electronVector_(electronVector), METVector_(METVector)

  {

    // If there is only one diphoton candidate, plot its invariant mass

    // going to come from diphoVector[]

    // if size diphoVector == 1, save diphoton object as TLorentzvector to plot invariant mass 

    // Save Diphoton object from flashgg::DiPhotonCandidate 


    //std::cout << "Got here" << std::endl;

    if (diphoVector_.size() == 1)
    {
      flashgg::DiPhotonCandidate dipho_ = diphoVector_[0];
      auto dipho = dipho_.p4();
      abe_dp_ = dipho;

    } 

    // Save Electrons 1 and 2 four momentum **IF** there are exactly two electrons. If more than 2, take leading pT? 

    //if (electronVector_.size() == 2)
    // Check for muons as well 

    if (electronVector_.size() == 2)
    {
      flashgg::Electron firstelec = electronVector_[0];
      flashgg::Electron secondelec = electronVector_[1];
      auto elec1 = firstelec.p4();
      auto elec2 = secondelec.p4();
      
      // if the first elec pt is higher, call it elec1 
      if (firstelec.pt() > secondelec.pt())
      {
        elec1 = firstelec.p4();
        elec2 = secondelec.p4();
      }
      // if the secobnd elec pt is higher, call it elec1 
      else if (secondelec.pt() > firstelec.pt())
      {
        elec1 = secondelec.p4();
        elec2 = firstelec.p4();
      }

      // Same pt
      else
      {
        elec1 = firstelec.p4();
        elec2 = secondelec.p4();  
      }

      elec1_ = elec1; // leading pt electron
      elec2_ = elec2; // subleading pt electron 

      // MET 
      if (METVector_.size() == 1)
      {
        flashgg::Met met__ = METVector_[0];
        auto met_ = met__.p4();
        MET_fourvec_ = met_;

        auto W = met_ + elec1;
        W1_TM_ = W.Mt();

      } 

    } 

    // semileptonic signal 
    else if (electronVector_.size() == 1)
    {
      flashgg::Electron firstelec = electronVector_[0];
      auto elec1 = firstelec.p4();
      elec1_ = elec1;

      // MET 
      if (METVector_.size() == 1)
      {
        flashgg::Met met__ = METVector_[0];
        auto met_ = met__.p4();
        MET_fourvec_ = met_;

        auto W = met_ + elec1;
        W1_TM_ = W.Mt();

      } 


    }

    // MET 

    //theMETcorpt_ = theMET_.getCorPt();

    //edm::Ptr<flashgg::Met> met = theMET;
    // auto dipho = dipho_.p4();
    // theMET_ = 

    // Save Higgs object (coming from fully leptonic, semi-leptonic, or fully hadronic)
    // First iteration: H->qqenu (where q's are specifially charm strange)
    



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

  float HHWWggCandidate::getCosThetaStar_CS(float ebeam) const {
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    reco::Candidate::LorentzVector h_lor = dp1_ + dp2_;
    TLorentzVector h;
    h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

    TVector3 boost = - h.BoostVector();
    p1.Boost(boost);
    p2.Boost(boost);
    reco::Candidate::LorentzVector a1_lor = dp1_;
    TLorentzVector a_1;
    a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;
    a_1.Boost(boost);

    TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
    CSaxis.Unit();

    return cos(   CSaxis.Angle( a_1.Vect().Unit() )    );
  }
  std::vector<float> HHWWggCandidate::CosThetaAngles() const {
    std::vector<float> helicityThetas;

    TLorentzVector Boosted_a1(0,0,0,0);
    Boosted_a1.SetPxPyPzE(dp1_.px(),dp1_.py(),dp1_.pz(),dp1_.energy()) ;
    TLorentzVector BoostedLeadingPhoton_a1(0,0,0,0);
    BoostedLeadingPhoton_a1.SetPxPyPzE(dp1_pho1_.px(),dp1_pho1_.py(),dp1_pho1_.pz(),dp1_pho1_.energy()) ;

    helicityThetas.push_back( HelicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1));

    TLorentzVector Boosted_a2(0,0,0,0);
    Boosted_a2.SetPxPyPzE(dp2_.px(),dp2_.py(),dp2_.pz(),dp2_.energy()) ;
    TLorentzVector BoostedLeadingPhoton_a2(0,0,0,0);
    BoostedLeadingPhoton_a2.SetPxPyPzE(dp2_pho1_.px(),dp2_pho1_.py(),dp2_pho1_.pz(),dp2_pho1_.energy()) ;

    helicityThetas.push_back( HelicityCosTheta(Boosted_a2, BoostedLeadingPhoton_a2));

    return helicityThetas;

  }

  float HHWWggCandidate::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
  {
    TVector3 BoostVector = Booster.BoostVector();
    Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
    return Boosted.CosTheta();
  }

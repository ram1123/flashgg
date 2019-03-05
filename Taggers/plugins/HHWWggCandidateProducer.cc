// Abe Tishelman-Charny
// 5 Feb 2019
//
// Starting with H4G, converting to HHWWgg 
// Currently have: photons, diphotons, vertex, gen particle
// Want to add:
//
// Electron
// Jets
// Met

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/HHWWggCandidate.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"


using namespace std;
using namespace edm;

//
// For ggqqlnu, I need a diphoton, lepton and neutrino. 
// Do I need a vector of photons? 
// Check VH and ZH for W info. Will need two W's: One on shell one off shell 
// Fully leptonic, semileptonic, fully hadronic 

namespace flashgg {
  // HHWWggCandidateProducer is a sub class or derived class of EDProducer 
  class HHWWggCandidateProducer : public EDProducer
  {
  public:
    //---typedef
    typedef math::XYZTLorentzVector LorentzVector;

    //---ctors
    HHWWggCandidateProducer();
    HHWWggCandidateProducer( const ParameterSet & );
  private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Photon> > photonToken_;
    Handle<View<flashgg::Photon> > photons;

    EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;

    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    Handle<View<reco::Vertex> > vertex;

    EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    Handle<View<reco::GenParticle> > genParticle;

    EDGetTokenT<View<Electron> > electronToken_;
    Handle<View<flashgg::Electron> > electrons;

    EDGetTokenT<View<Muon> > muonToken_;
    Handle<View<flashgg::Muon> > muons;

    EDGetTokenT<View<flashgg::Met> > METToken_;
    Handle<View<flashgg::Met> > METs;

    //---ID selector
    ConsumesCollector cc_;
    CutBasedDiPhotonObjectSelector idSelector_;

    //----output collection
    auto_ptr<vector<HHWWggCandidate> > HHWWggColl_;

  };

  //---constructors
  HHWWggCandidateProducer::HHWWggCandidateProducer( ):
  photonToken_(),
  diphotonToken_(),
  genParticleToken_(),
  electronToken_(),
  muonToken_(),
  METToken_(),
  cc_( consumesCollector() ),
  idSelector_( ParameterSet(), cc_ )

  {}

    //---standard
    HHWWggCandidateProducer::HHWWggCandidateProducer( const ParameterSet & pSet):
    photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
    diphotonToken_( consumes<View<DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    electronToken_( consumes<View<Electron> >( pSet.getParameter<InputTag> ( "ElectronTag" ) ) ), // prev iConfig.get..
    muonToken_( consumes<View<Muon> >( pSet.getParameter<InputTag> ( "MuonTag" ) ) ),
    //METToken_( consumes<View<flashgg::Met> >( pSet.getParameter<InputTag> ( "METTag" ) ) ),
    METToken_( consumes<View<Met> >( pSet.getParameter<InputTag> ( "METTag" ) ) ),

    cc_( consumesCollector() ),
    idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

    {
      produces<vector<HHWWggCandidate> > ();
    }

    // I'm not sure where HHWWggCandidateProducer::produce is called, but this is called event by event 
    void HHWWggCandidateProducer::produce( Event &event, const EventSetup & )
    {
      event.getByToken( photonToken_, photons );
      event.getByToken( diphotonToken_, diphotons );
      event.getByToken( vertexToken_, vertex );
      event.getByToken( genParticleToken_, genParticle );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );

      //---output collection
      std::unique_ptr<vector<HHWWggCandidate> > HHWWggColl_( new vector<HHWWggCandidate> );

      edm::Ptr<reco::Vertex> vertex_zero = vertex->ptrAt(0);
      //std::cout << "vertex_zero = " << vertex_zero << std::endl;
      reco::GenParticle::Point genVertex;
      

      //std::cout << "diphotons->size() = " << diphotons->size() << std::endl;

      std::vector<const flashgg::Photon*> phosTemp;
      // for( unsigned int dpIndex = 0; dpIndex < diphotons->size(); dpIndex++ )
      // {
      //   edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diphotons->ptrAt( dpIndex );
      //   flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
      //   atLeastOneDiphoPass |= idSelector_(*thisDPPointer, event);
      // }
      int n_electrons = electrons->size();
      int n_muons = muons->size();
      int n_photons = photons->size();
      //std::cout << "n_photons = " << n_photons << std::endl;
      //if (n_photons == 0) std::cout << "***************************n_photons = " << n_photons << std::endl;
      int n_diphotons = diphotons->size();
      int n_METs = METs->size(); // Should be 1..But using this as a way to obtain met four vector 
      std::vector<flashgg::Photon> phoVector;
      std::vector<flashgg::DiPhotonCandidate> diphoVector;
      std::vector<flashgg::Electron> electronVector;
      std::vector<flashgg::Muon> muonVector;
      std::vector<flashgg::Met> METVector;

      // Save gen particles 
      // Save all gen particles, sort out later 
      //std::vector
      //std::vector<> GENelectronVector;

      //flashgg::Met
      // if (atLeastOneDiphoPass)
      // {
        // Append electronVector
        for( int electronIndex = 0; electronIndex < n_electrons; electronIndex++ )
        {
          edm::Ptr<flashgg::Electron> elec = electrons->ptrAt( electronIndex );
          flashgg::Electron * thisElecPointer = const_cast<flashgg::Electron *>(elec.get());
          electronVector.push_back(*thisElecPointer);
        }

        // Append muonVector
        for( int muonIndex = 0; muonIndex < n_muons; muonIndex++ )
        {
          edm::Ptr<flashgg::Muon> mlep = muons->ptrAt( muonIndex );
          flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(mlep.get());
          muonVector.push_back(*thisMuonPointer);
        }


        // Append phoVector
        for( int phoIndex = 0; phoIndex < n_photons; phoIndex++ )
        {
          edm::Ptr<flashgg::Photon> pho = photons->ptrAt( phoIndex );
          flashgg::Photon * thisPPointer = const_cast<flashgg::Photon *>(pho.get());
          phoVector.push_back(*thisPPointer);
        }
        // Append diphoVector if it passes preselection 
        bool PassPS = false;
        for( int diphoIndex = 0; diphoIndex < n_diphotons; diphoIndex++ )
        {
          edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho.get());
          //---at least one diphoton should pass the low mass hgg pre-selection
          PassPS = false;
          PassPS |= idSelector_(*thisDPPointer, event);
          if (PassPS) diphoVector.push_back(*thisDPPointer);
          
        }

        // Trying MET vector
        if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
        for( int METIndex = 0; METIndex < n_METs; METIndex++ )
        {
          edm::Ptr<flashgg::Met> m_entry = METs->ptrAt( METIndex );
          flashgg::Met * thisMETPointer = const_cast<flashgg::Met *>(m_entry.get());
          METVector.push_back(*thisMETPointer);
        }

        // Want to save GEN particles to compare variables to RECO 
        // Only want to save GEN particles coming from mother particle of interest

        // Mother Daughter pdgid pairs 
        // if a particle of pdgid [0] came from pdgid [1], save it 
        std::vector<std::vector<int>> md_pairs = {};
        md_pairs.push_back({24,25}); // W from H 
        md_pairs.push_back({22,25}); // Photon from H 
        md_pairs.push_back({11,24}); // Electron from W
        md_pairs.push_back({12,24}); // Electron neutrino from W
        md_pairs.push_back({13,24}); // Muon from W
        md_pairs.push_back({14,24}); // Muon neutrino from W 
        
        std::vector<int> quark_pdgids = {1,2,3,4,5};

        for (unsigned int i = 0; i < quark_pdgids.size(); i++){
          int qid = quark_pdgids[i];
          md_pairs.push_back({qid,24}); // Quark from W
        }

        // vector to store genparticles in 
        vector<reco::GenParticle> genParticlesVector;

        // If MC event 
        if (! event.isRealData() ){
          // For each gen particle in event 
          for(size_t g=0; g < genParticle->size(); g++){
            auto gen = genParticle->ptrAt(g);

            // If the particle has a mother 
            if (gen->mother(0) != 0){
              int pid = gen->pdgId();
              int pmotid = gen->mother(0)->pdgId();              

              for(unsigned int i = 0; i < md_pairs.size(); i++){
                  int vec_id = md_pairs[i][0];
                  int vec_mid = md_pairs[i][1]; 

                  // if event gen particle and mother are on list of desired particles, add to genParticlesVector 
                  if ( (abs(pid) == abs(vec_id)) && (abs(pmotid) == abs(vec_mid))){ 
                    cout << "Found " << abs(pid) << " from " << abs(pmotid) << endl;
                    reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
                    genParticlesVector.push_back(*thisGENPointer);

                  }


              }
            }
          }
        }


      //   // If MC event
      //   if (! event.isRealData() ) 
      //   {

      //     for(size_t g=0; g < genParticle->size(); g++) {

      //         auto gen = genParticle->ptrAt(g);

      //         // 2212: proton 
      //         // If the particle has a mother 
      //         if (gen->mother(0) != 0){
      //           int pid = gen->pdgId();
      //           int pmotid = gen->mother(0)->pdgId();

      //           if (abs(pid) == 24 && pmotid == 25 ){
      //             cout << "Found W boson from Higgs boson" << endl;
      //             reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
      //             genParticlesVector.push_back(*thisGENPointer);

      //           }

      //           if (abs(pid) == 22 && pmotid == 25 ){
      //             cout << "Found photon from Higgs boson" << endl;
      //             reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
      //             genParticlesVector.push_back(*thisGENPointer);

      //           }

      //           if (abs(pid) == 11 && abs(pmotid) == 24){
      //             cout << "Found electron from W boson" << endl;
      //             reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
      //             genParticlesVector.push_back(*thisGENPointer);
      //           }

      //           if (abs(pid) == 12 && abs(pmotid) == 24 ){
      //             cout << "Found neutrino from W boson" << endl;
      //             reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
      //             genParticlesVector.push_back(*thisGENPointer);

      //           }

      //           for (unsigned int i = 0; i < quark_pdgids.size(); i++){
      //               int val = quark_pdgids[i];
      //               if (abs(pid) == val && abs(pmotid) == 24 ){
      //                   std::cout << "found a quark from W boson" << endl;
      //                   reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
      //                   genParticlesVector.push_back(*thisGENPointer);
      //               }

      //           }

                

      //         }

      //     }


      //   //   // for each gen particle in event 
      //   //   for( auto &part : *genParticle ) { 

      //   //     if( part.pdgId() != 2212 || part.vertex().z() != 0. ) 
      //   //     {
      //   //       genVertex = part.vertex();
      //   //     }
            
      //   // }

      // }


        //HHWWggCandidate HHWWgg(phoVector, vertex_zero, genVertex); // I think one for each event 
        //HHWWggCandidate HHWWgg(diphoVector, phoVector, vertex_zero, genVertex); 
        //HHWWggCandidate HHWWgg(diphoVector, phoVector, vertex_zero, genVertex, electronVector, theMET); 
        //HHWWggCandidate HHWWgg(diphoVector, phoVector, vertex_zero, genVertex, electronVector, METVector);  genParticlesVector
        HHWWggCandidate HHWWgg(diphoVector, phoVector, vertex_zero, genVertex, electronVector, muonVector, METVector, genParticlesVector);
        HHWWggColl_->push_back(HHWWgg);

      //}
      event.put( std::move(HHWWggColl_) );
    }
  }
  typedef flashgg::HHWWggCandidateProducer FlashggHHWWggCandidateProducer;
  DEFINE_FWK_MODULE( FlashggHHWWggCandidateProducer );

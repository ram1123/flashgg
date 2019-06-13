// Abe Tishelman-Charny
// 5 Feb 2019
//
// Started with H4G, converting to HHWWgg 

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/HHWWggCandidate.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

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

    //---Outtree 
    edm::Service<TFileService> fs;  
    TH1F* cutFlow;

  private:
    double genTotalWeight;
    void produce( Event &, const EventSetup & ) override;

    // Adding Jets 
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;

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

    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;

    Handle<View<reco::Vertex> > vertices;

    EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<edm::TriggerResults> triggerRECO_;
    edm::EDGetTokenT<edm::TriggerResults> triggerPAT_;
    edm::EDGetTokenT<edm::TriggerResults> triggerFLASHggMicroAOD_;
    edm::Handle<double>  rho;

    std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > JetToken_;

    //---ID selector
    ConsumesCollector cc_;
    CutBasedDiPhotonObjectSelector idSelector_;

    //----output collection
    auto_ptr<vector<HHWWggCandidate> > HHWWggColl_;

    // variables from WHLeptonicTagProducer

    double leptonPtThreshold_;
    double muonEtaThreshold_;
    double leadPhoOverMassThreshold_;
    double subleadPhoOverMassThreshold_;
    double MVAThreshold_;
    double deltaRMuonPhoThreshold_;
    double jetsNumberThreshold_;
    double jetPtThreshold_;
    double jetEtaThreshold_;
    double muPFIsoSumRelThreshold_;
    double PhoMVAThreshold_;
    double METThreshold_;
    bool useVertex0only_;
    double deltaRJetMuonThreshold_;
    double deltaRPhoLeadJet_;
    double deltaRPhoSubLeadJet_;

    double DeltaRTrkElec_;
    double TransverseImpactParam_;
    double LongitudinalImpactParam_;

    double deltaRPhoElectronThreshold_;
    double deltaMassElectronZThreshold_;

    bool hasGoodElec = false;
    bool hasGoodMuons = false;

    vector<double> nonTrigMVAThresholds_;
    vector<double> nonTrigMVAEtaCuts_;

    double electronIsoThreshold_;
    double electronNumOfHitsThreshold_;
    vector<double> electronEtaThresholds_;
    bool useElectronMVARecipe_;
    bool useElectronLooseID_;

    edm::InputTag genInfo_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
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
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    electronToken_( consumes<View<Electron> >( pSet.getParameter<InputTag> ( "ElectronTag" ) ) ), 
    muonToken_( consumes<View<Muon> >( pSet.getParameter<InputTag> ( "MuonTag" ) ) ),
    METToken_( consumes<View<Met> >( pSet.getParameter<InputTag> ( "METTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( pSet.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    rhoTag_( consumes<double>( pSet.getParameter<InputTag>( "rhoTag" ) ) ),
    triggerRECO_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("RECOfilters") ) ),
    triggerPAT_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("PATfilters") ) ),
    triggerFLASHggMicroAOD_( consumes<edm::TriggerResults>( pSet.getParameter<InputTag>("FLASHfilters") ) ),
    cc_( consumesCollector() ),
    idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

    {
      genInfo_ = pSet.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );
      cutFlow = fs->make<TH1F> ("cutFlow","Cut Flow",10,0,10);

      leptonPtThreshold_ = pSet.getParameter<double>( "leptonPtThreshold");
      muonEtaThreshold_ = pSet.getParameter<double>( "muonEtaThreshold");
      leadPhoOverMassThreshold_ = pSet.getParameter<double>( "leadPhoOverMassThreshold");
      subleadPhoOverMassThreshold_ = pSet.getParameter<double>( "subleadPhoOverMassThreshold");
      MVAThreshold_ = pSet.getParameter<double>( "MVAThreshold");
      deltaRMuonPhoThreshold_ = pSet.getParameter<double>( "deltaRMuonPhoThreshold");
      jetsNumberThreshold_ = pSet.getParameter<double>( "jetsNumberThreshold");
      jetPtThreshold_ = pSet.getParameter<double>( "jetPtThreshold");
      jetEtaThreshold_ = pSet.getParameter<double>( "jetEtaThreshold");
      muPFIsoSumRelThreshold_ = pSet.getParameter<double>( "muPFIsoSumRelThreshold");
      PhoMVAThreshold_ = pSet.getParameter<double>( "PhoMVAThreshold");
      METThreshold_ = pSet.getParameter<double>( "METThreshold");
      useVertex0only_              = pSet.getParameter<bool>("useVertex0only");
      deltaRJetMuonThreshold_ = pSet.getParameter<double>( "deltaRJetMuonThreshold");
      deltaRPhoLeadJet_ = pSet.getParameter<double>( "deltaRPhoLeadJet");
      deltaRPhoSubLeadJet_ = pSet.getParameter<double>( "deltaRPhoSubLeadJet");

      DeltaRTrkElec_ = pSet.getParameter<double>( "DeltaRTrkElec");
      TransverseImpactParam_ = pSet.getParameter<double>( "TransverseImpactParam");
      LongitudinalImpactParam_ = pSet.getParameter<double>( "LongitudinalImpactParam");

      deltaRPhoElectronThreshold_ = pSet.getParameter<double>( "deltaRPhoElectronThreshold");
      deltaMassElectronZThreshold_ = pSet.getParameter<double>( "deltaMassElectronZThreshold");

      nonTrigMVAThresholds_ =  pSet.getParameter<vector<double > >( "nonTrigMVAThresholds");
      nonTrigMVAEtaCuts_ =  pSet.getParameter<vector<double > >( "nonTrigMVAEtaCuts");
      electronIsoThreshold_ = pSet.getParameter<double>( "electronIsoThreshold");
      electronNumOfHitsThreshold_ = pSet.getParameter<double>( "electronNumOfHitsThreshold");
      electronEtaThresholds_ = pSet.getParameter<vector<double > >( "electronEtaThresholds");
      useElectronMVARecipe_=pSet.getParameter<bool>("useElectronMVARecipe");
      useElectronLooseID_=pSet.getParameter<bool>("useElectronLooseID");

      auto jetTags = pSet.getParameter<std::vector<edm::InputTag> > ( "JetTags" ); 
      for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }

      produces<vector<HHWWggCandidate> > ();
    }

    void HHWWggCandidateProducer::produce( Event &event, const EventSetup & )
    {
      // Get particle objects
      event.getByToken( photonToken_, photons );
      event.getByToken( diphotonToken_, diphotons );
      event.getByToken( vertexToken_, vertex );
      event.getByToken( genParticleToken_, genParticle );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );
      event.getByToken( mvaResultToken_, mvaResults );
      event.getByToken( vertexToken_, vertices );
      event.getByToken( rhoTag_, rho);

      double rho_    = *rho;

      // Set cut booleans
      std::vector<double> Cut_Results = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Cut_Results[i] = 1: Event Passed Cut i 
      std::vector<double> Cut_Variables = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Cut_Results[i] = 1: Event Passed Cut i 
      bool SLW_Tag = 0; // one SL W tagged  
      bool PS_dipho_tag = 0; // preselected diphoton 
      bool Pass_PS = 0; // Preselected Diphoton (may also pass SLW_tag)

      //---output collection
      std::unique_ptr<vector<HHWWggCandidate> > HHWWggColl_( new vector<HHWWggCandidate> );
      //edm::Ptr<reco::Vertex> vertex_zero = vertex->ptrAt(0);
      //reco::GenParticle::Point genVertex;
      int n_electrons = electrons->size();
      int n_muons = muons->size();
      int n_photons = photons->size();
      int n_diphotons = diphotons->size();
      int n_METs = METs->size(); // Should be 1, but using as a way to obtain met four vector 
      // int n_jets = 0;
      int n_good_electrons = 0;
      int n_good_muons = 0;
      int n_good_photons = 0;
      int n_good_diphotons = 0;
      int num_dipho = 0;
      // int n_METs = METs->size(); // Should be 1, but using as a way to obtain met four vector 
      int n_good_jets = 0;
      std::vector<flashgg::Photon> phoVector;
      std::vector<flashgg::DiPhotonCandidate> diphoVector;
      std::vector<flashgg::Electron> electronVector;
      std::vector<flashgg::Muon> muonVector;
      std::vector<flashgg::Met> METVector;
      std::vector<flashgg::Jet> JetVector;

      // Saved Objects after selections
      std::vector<flashgg::Jet> tagJets_;
      std::vector<flashgg::Muon> goodMuons_;
      std::vector<flashgg::Electron> goodElectrons_; 
      std::vector<flashgg::Met> theMET_;
      std::vector<flashgg::DiPhotonCandidate> diphoVector_; // diphoton from W tagged events  
      // std::vector<edm::Ptr<Jet> > tagJets_;
      // std::vector<edm::Ptr<flashgg::Muon> > goodMuons_;
      // std::vector<edm::Ptr<Electron> > goodElectrons_; 
      // Ptr<flashgg::Met> theMET_;

//-----------------------------------------------------------------------------------------------------------

      // For BSM Search, take 0th diphoton because expect everything to come from same vertex 

      bool passMETfilters=1;
      //Get trigger results relevant to MET filters                                                                                                                                              

      edm::Handle<edm::TriggerResults> triggerBits;
      if(! event.isRealData() )
          event.getByToken( triggerPAT_, triggerBits );
      else
          event.getByToken( triggerRECO_, triggerBits );

      edm::Handle<edm::TriggerResults> triggerFLASHggMicroAOD;
      event.getByToken( triggerFLASHggMicroAOD_, triggerFLASHggMicroAOD );
      const edm::TriggerNames &triggerNames = event.triggerNames( *triggerBits );

      //check if passMETfilters 
      std::vector<std::string> flagList {"Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_goodVertices","Flag_eeBadScFilter"};
      for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ )
          {
              if(!triggerBits->accept(i))
                  for(size_t j=0;j<flagList.size();j++)
                      {
                          if(flagList[j]==triggerNames.triggerName(i))
                              {
                                  passMETfilters=0;
                                  break;
                              }
                      }
          }

      std::vector<std::string> flashggFlagList {"flag_BadChargedCandidateFilter","flag_BadPFMuonFilter","flag_globalTightHalo2016Filter"};
      const edm::TriggerNames &flashggtriggerNames = event.triggerNames( *triggerFLASHggMicroAOD );
      for( unsigned int i = 0; i < flashggtriggerNames.triggerNames().size(); i++ )
          {
              if(!triggerFLASHggMicroAOD->accept(i))
                  for(size_t j=0;j<flagList.size();j++)
                      {
                          if(flagList[j]==flashggtriggerNames.triggerName(i))
                              {
                                  passMETfilters=0;
                                  break;
                              }
                      }
          }

      bool photonSelection = false;
      double idmva1 = 0.;
      double idmva2 = 0.;

      // Get first (zeroeth) diphoton
      // if (n_diphotons != 0){

      if (diphotons->size() > 0){
        edm::Ptr<flashgg::DiPhotonCandidate> first_dipho = diphotons->ptrAt( 0 ); 
        flashgg::DiPhotonCandidate * firstDPPointer = const_cast<flashgg::DiPhotonCandidate *>(first_dipho.get());
        Pass_PS |= idSelector_(*firstDPPointer, event);

        if ( Pass_PS ) {
          PS_dipho_tag = 1;
        } 
      }

      //else{ num_dipho = 0;}

      for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) {
        if (diphoIndex != 0) continue; // only look at 0th diphoton 
        hasGoodElec = false;
        hasGoodMuons = false;

        // Check MET Filters
        //if(!passMETfilters) {continue;}

        // Check if useVertex0only_
        if(useVertex0only_) // If only 0th vertex of diphotons is used, continue on all other vertex indices 
            if(diphotons->ptrAt(diphoIndex)->vertexIndex()!=0)
                {continue;}

        // Diphoton 

        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex ); 
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

        // Diphoton selections 

        // Preselection
        // flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho.get());
        // Pass_PS |= idSelector_(*thisDPPointer, event);
        // if ( Pass_PS ) {
        //   PS_dipho_tag = 1;
        // }    

        // leading/subleading photon pt 
        if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) { continue; }
        if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { continue; }

        // leading/subleading photon MVA
        idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) { continue; }

        // Diphoton MVA 
        if( mvares->result < MVAThreshold_ ) { continue; } // diphotonmva 
        photonSelection = true;

        // Electrons 
        std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho,vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
                                                                          useElectronMVARecipe_,useElectronLooseID_,
                                                                          deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_,
                                                                          rho_, event.isRealData() );

        
        // Muons                                                                   
        std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_,
        muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );

        n_good_electrons = goodElectrons.size();
        n_good_muons = goodMuons.size();

        hasGoodElec = ( goodElectrons.size() > 0 );
        hasGoodMuons = ( goodMuons.size() > 0 );
                    
        if( !hasGoodElec && !hasGoodMuons ) { continue; }

        // Jets 
        unsigned int jetCollectionIndex = diphotons->ptrAt( diphoIndex )->jetCollectionIndex();
        edm::Handle<edm::View<flashgg::Jet> > Jets_;
        event.getByToken( jetTokens_[jetCollectionIndex], Jets_);
        std::vector<edm::Ptr<Jet> > tagJets;

        // Jet Selection 
        for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ ) 
            {
                bool keepJet=true;
                edm::Ptr<flashgg::Jet> thejet = Jets_->ptrAt( candIndex_outer );
                // cout << " jet number = " << candIndex_outer << endl;
                // cout << " jet pt = " << thejet->pt() << endl;
                if(!thejet->passesJetID  ( flashgg::Tight2017 ) ) { continue; }
                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { keepJet=false; }

                if( thejet->pt() < jetPtThreshold_ ) { keepJet=false; }
                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                dipho->subLeadingPhoton()->superCluster()->phi() );
                
                if( dRPhoLeadJet < deltaRPhoLeadJet_ || dRPhoSubLeadJet < deltaRPhoSubLeadJet_ ) { keepJet=false; }
                if( hasGoodElec ) 
                    for( unsigned int electronIndex = 0; electronIndex < goodElectrons.size(); electronIndex++ ) 
                        {
                            Ptr<flashgg::Electron> electron = goodElectrons[electronIndex];
                            float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), electron->eta(), electron->phi() ) ;
                            if( dRJetElectron < deltaRJetMuonThreshold_ ) { keepJet=false; }
                        }
                if( hasGoodMuons ) 
                    for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) 
                        {
                            Ptr<flashgg::Muon> muon = goodMuons[muonIndex];
                            float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;
                            if( dRJetMuon < deltaRJetMuonThreshold_ ) { keepJet=false; }
                        }
                if(keepJet)
                    tagJets.push_back( thejet );
            }
        n_good_jets = tagJets.size();

        // MET 
        if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
        Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );


        //n_good_jets = tagJets.size();
        // Semileptonic W Tag 
        if( (tagJets.size() < jetsNumberThreshold_) && photonSelection && theMET->getCorPt()>METThreshold_) { 
          SLW_Tag = 1;
          // n_good_jets = tagJets.size();
          // n_good_electrons = goodElectrons.size();
          // n_good_muons = goodMuons.size();

          // how to go from pointer 'jet' to object '*thisJetPointer'
          //flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(jet.get());
          //JetVector.push_back(*thisJetPointer);

          for (unsigned int i = 0; i < tagJets.size(); i++){
            auto tag_jet = tagJets[i];
            flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(tag_jet.get());
            tagJets_.push_back(*thisJetPointer);
          }
          for (unsigned int i = 0; i < goodElectrons.size(); i++){
            auto good_elec = goodElectrons[i];
            flashgg::Electron * thisElectronPointer = const_cast<flashgg::Electron *>(good_elec.get());
            goodElectrons_.push_back(*thisElectronPointer);
          }
          for (unsigned int i = 0; i < goodMuons.size(); i++){
            auto good_muon = goodMuons[i];
            flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(good_muon.get());
            goodMuons_.push_back(*thisMuonPointer);
          }
          // Store MET as one element in MET vector 
          if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          for( int METIndex = 0; METIndex < n_METs; METIndex++ )
          {
            edm::Ptr<flashgg::Met> m_entry = METs->ptrAt( METIndex );
            flashgg::Met * thisMETPointer = const_cast<flashgg::Met *>(m_entry.get());
            theMET_.push_back(*thisMETPointer);
          }

          // dipho
          edm::Ptr<flashgg::DiPhotonCandidate> dipho_ = dipho;
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho_.get());
          // diphoton preselection 
          //Pass_PS = false;
          // Pass_PS |= idSelector_(*thisDPPointer, event);
          // if ( (Pass_PS) && (!one_ps_dp) ) {
          //   one_ps_dp = 1;
          //   Pass_PS = 1;
          // }
          // if (PassPS) {
          diphoVector_.push_back(*thisDPPointer);

        }

        // size_t vtx = (size_t)dipho->jetCollectionIndex();
        // edm::Handle<edm::View<flashgg::Jet> > jets;
        // event.getByToken( jetTokens_[vtx], jets);
        // n_jets = jets->size();
        // for( size_t ijet=0; ijet < jets->size(); ++ijet ) { //jets are ordered in pt
        //     auto jet = jets->ptrAt(ijet);

        //     flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(jet.get());
        //     JetVector.push_back(*thisJetPointer);

        // }

      }


      // -- Save good objects 
      // only care about events with a preselected diphoton because this is what we'll likely find in a signal event 
        // bool photonSelection = false;
        // double idmva1 = 0.;
        // double idmva2 = 0.;

        // for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) {
            // cout << "diphoindex = " << diphoIndex << endl;
            // cout << "diphotons->ptrAt(diphoIndex)->vertexIndex() = " << diphotons->ptrAt(diphoIndex)->vertexIndex() << endl;
            // hasGoodElec = false;
            // hasGoodMuons = false;
            // if(!passMETfilters) {continue;}
            // if(useVertex0only_) // If only 0th vertex of diphotons is used, continue on all other vertex indices 
            //     if(diphotons->ptrAt(diphoIndex)->vertexIndex()!=0)
            //         {continue;}
            // unsigned int jetCollectionIndex = diphotons->ptrAt( diphoIndex )->jetCollectionIndex();
            // edm::Handle<edm::View<flashgg::Jet> > Jets_;
            // event.getByToken( jetTokens_[jetCollectionIndex], Jets_);
            // cout << "-----------------------------------\n";
            // cout << "Before selections:" << endl;
            // cout << "Jets_->size() = " << Jets_->size() << endl;    
            // cout << "(jet number, pt)" << endl;
            // for (unsigned int i = 0; i < Jets_->size(); i ++ ){
            //         edm::Ptr<flashgg::Jet> thejet__ = Jets_->ptrAt( i );
            //         cout << "(" << i << ", " << thejet__->pt() << ")" << endl;
            // }
            // std::vector<edm::Ptr<Jet> > tagJets;
            // edm::Ptr<flashgg::Met>  tagMETs;

            // edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
            // edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );
            
            // if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) { continue; }
            // if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { continue; }
          
            // idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            // idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        
            // if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) { continue; }
            // if( mvares->result < MVAThreshold_ ) { continue; } // diphotonmva 
            
            // photonSelection = true;
            // std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_,
            //         muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );
            
            // std::vector<edm::Ptr<Electron> >goodElectrons = selectStdElectrons( electrons->ptrs(), dipho,vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
            //                                                                     useElectronMVARecipe_,useElectronLooseID_,
            //                                                                     deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_,
            //                                                                     rho_, event.isRealData() );
            
            // hasGoodElec = ( goodElectrons.size() > 0 );
            // hasGoodMuons = ( goodMuons.size() > 0 );
                        
            // if( !hasGoodElec && !hasGoodMuons ) { continue; }
            //including SFs for leading muon or electron 
            // if(goodMuons.size()>0){
            //     //whleptonictags_obj.includeWeightsByLabel( *goodMuons.at(0), "MuonWeight");
            // } else if (goodElectrons.size() > 0){
            //     //whleptonictags_obj.includeWeights( *goodElectrons.at(0));
            // }
            // cout << "-----------------------------------\n";
            // cout << "After selections:" << endl;
            // for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ ) 
            //     {
            //         bool keepJet=true;
            //         edm::Ptr<flashgg::Jet> thejet = Jets_->ptrAt( candIndex_outer );
            //         cout << " jet number = " << candIndex_outer << endl;
            //         cout << " jet pt = " << thejet->pt() << endl;
            //         if(!thejet->passesJetID  ( flashgg::Tight2017 ) ) { continue; }
            //         if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { keepJet=false; }

            //         if( thejet->pt() < jetPtThreshold_ ) { keepJet=false; }
            //         float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
            //         float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
            //                                         dipho->subLeadingPhoton()->superCluster()->phi() );
                    
            //         if( dRPhoLeadJet < deltaRPhoLeadJet_ || dRPhoSubLeadJet < deltaRPhoSubLeadJet_ ) { keepJet=false; }
            //         if( hasGoodElec ) 
            //             for( unsigned int electronIndex = 0; electronIndex < goodElectrons.size(); electronIndex++ ) 
            //                 {
            //                     Ptr<flashgg::Electron> electron = goodElectrons[electronIndex];
            //                     float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), electron->eta(), electron->phi() ) ;
            //                     if( dRJetElectron < deltaRJetMuonThreshold_ ) { keepJet=false; }
            //                 }
            //         if( hasGoodMuons ) 
            //             for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) 
            //                 {
            //                     Ptr<flashgg::Muon> muon = goodMuons[muonIndex];
            //                     float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;
            //                     if( dRJetMuon < deltaRJetMuonThreshold_ ) { keepJet=false; }
            //                 }
            //         if(keepJet)
            //             tagJets.push_back( thejet );
            //     }
            //------>MET info
            // if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
            // Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

            //has Nj<thresh; photonSelection; muons>=1||electrons>=1; metPt>thresh
            // if( (tagJets.size() < jetsNumberThreshold_) && photonSelection && theMET->getCorPt()>METThreshold_) {
            //     //cout << "Tagging semi-leptonic W decay " << endl;
            //     if (!SLW_Tag){
            //         SLW_Tag = 1;
            //         //SLW_tag = 1; 
            //       }
            //     // how to go from pointer 'jet' to object '*thisJetPointer'
            //     //flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(jet.get());
            //     //JetVector.push_back(*thisJetPointer);

            //     for (unsigned int i = 0; i < tagJets.size(); i++){
            //       auto tag_jet = tagJets[i];
            //       flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(tag_jet.get());
            //       tagJets_.push_back(*thisJetPointer);
            //     }
            //     for (unsigned int i = 0; i < goodElectrons.size(); i++){
            //       auto good_elec = goodElectrons[i];
            //       flashgg::Electron * thisElectronPointer = const_cast<flashgg::Electron *>(good_elec.get());
            //       goodElectrons_.push_back(*thisElectronPointer);
            //     }
            //     for (unsigned int i = 0; i < goodMuons.size(); i++){
            //       auto good_muon = goodMuons[i];
            //       flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(good_muon.get());
            //       goodMuons_.push_back(*thisMuonPointer);
            //     }
            //     // Store MET as one element in MET vector 
            //     if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
            //     for( int METIndex = 0; METIndex < n_METs; METIndex++ )
            //     {
            //       edm::Ptr<flashgg::Met> m_entry = METs->ptrAt( METIndex );
            //       flashgg::Met * thisMETPointer = const_cast<flashgg::Met *>(m_entry.get());
            //       theMET_.push_back(*thisMETPointer);
            //     }

            //     // dipho
            //     edm::Ptr<flashgg::DiPhotonCandidate> dipho_ = dipho;
            //     flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho_.get());
            //     // diphoton preselection 
            //     //Pass_PS = false;
            //     // Pass_PS |= idSelector_(*thisDPPointer, event);
            //     // if ( (Pass_PS) && (!one_ps_dp) ) {
            //     //   one_ps_dp = 1;
            //     //   Pass_PS = 1;
            //     // }
            //     // if (PassPS) {
            //     diphoVector_.push_back(*thisDPPointer);

            // }
        // }

        // // Check that at least one diphoton passes all preselections 

        // for( int diphoIndex = 0; diphoIndex < n_diphotons; diphoIndex++ )
        // {
        //   edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
        //   flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho.get());
        //   // Pass_PS |= idSelector_(*thisDPPointer, event);
        //   // if ( (Pass_PS) && (!one_ps_dp) ) {
        //   //   one_ps_dp = 1;
        //   //   Pass_PS = 1;
        //   // }          
        // }

        // If at least one diphoton passed preselection, save all objects:
        // if(one_ps_dp){
          // Diphotons
          for( int diphoIndex = 0; diphoIndex < n_diphotons; diphoIndex++ )
          {
            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
            flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho.get());
            diphoVector.push_back(*thisDPPointer);          
          }

          // Photons   
          for( int phoIndex = 0; phoIndex < n_photons; phoIndex++ )
          {
            edm::Ptr<flashgg::Photon> pho = photons->ptrAt( phoIndex );
            flashgg::Photon * thisPPointer = const_cast<flashgg::Photon *>(pho.get());
            phoVector.push_back(*thisPPointer);
          }

          // Electrons
          for( int electronIndex = 0; electronIndex < n_electrons; electronIndex++ )
          {
            edm::Ptr<flashgg::Electron> elec = electrons->ptrAt( electronIndex );
            flashgg::Electron * thisElecPointer = const_cast<flashgg::Electron *>(elec.get());
            //if ( elec->pt() >= leptonPtThreshold_) lep_thres_cut = 1; // check lepton pt cut 
            electronVector.push_back(*thisElecPointer);
          }

          // Muons 
          for( int muonIndex = 0; muonIndex < n_muons; muonIndex++ )
          {
            edm::Ptr<flashgg::Muon> mlep = muons->ptrAt( muonIndex );
            flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(mlep.get());
            //if ( mlep->pt() >= leptonPtThreshold_) lep_thres_cut = 1; // check lepton pt cut 
            muonVector.push_back(*thisMuonPointer); 
          }

          // MET
          // Store MET as one element in MET vector 
          if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          for( int METIndex = 0; METIndex < n_METs; METIndex++ )
          {
            edm::Ptr<flashgg::Met> m_entry = METs->ptrAt( METIndex );
            flashgg::Met * thisMETPointer = const_cast<flashgg::Met *>(m_entry.get());
            METVector.push_back(*thisMETPointer);
          }

        // } // Done saving all objects 

//////////////////////////////////////////////////////////////////

        //-- Get GEN particles 

        // Want to save GEN particles to compare variables to RECO 
        // Only want to save GEN particles coming from mother particles of interest

        // Mother Daughter pdgid pairs 
        // if a particle of abs(pdgid) = [0] came from abs(pdgid) = [1], store it 
        std::vector<std::vector<int>> md_pairs = {}; // mother daughter pairs 
        md_pairs.push_back({24,25}); // W from H 
        md_pairs.push_back({22,25}); // Photon from H 
        md_pairs.push_back({11,24}); // Electron from W
        md_pairs.push_back({12,24}); // Electron neutrino from W
        md_pairs.push_back({13,24}); // Muon from W
        md_pairs.push_back({14,24}); // Muon neutrino from W 
        
        //std::vector<int> quark_pdgids = {1,2,3,4,5}; // if you want to look for b quarks coming from W's 
        std::vector<int> quark_pdgids = {1,2,3,4}; // down, up, strange, charm 

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

              if ( ( abs(pid) == 5) && (abs(pmotid) == 24) ) cout << "***B quark Found!***" << endl;           

              for(unsigned int i = 0; i < md_pairs.size(); i++){
                int vec_id = md_pairs[i][0];
                int vec_mid = md_pairs[i][1]; 

                // if event gen particle and mother are on list of desired particles, add to genParticlesVector 
                if ( (abs(pid) == abs(vec_id)) && (abs(pmotid) == abs(vec_mid))){ 
                  //cout << "Found " << abs(pid) << " from " << abs(pmotid) << endl;
                  reco::GenParticle * thisGENPointer = const_cast<reco::GenParticle *>(gen.get());
                  genParticlesVector.push_back(*thisGENPointer);

                }
              }
            }
          }
        }

        Cut_Variables[0] = n_good_electrons;
        Cut_Variables[1] = n_good_muons;
        Cut_Variables[2] = n_good_jets;
        Cut_Variables[3] = passMETfilters; 

        // CutFlow for Semi-Leptonic Signal
        if (n_photons >= 2) Cut_Results[1] = 1.0;
        if (n_good_electrons >= 1) Cut_Results[2] = 1.0;
        if (n_good_muons >= 1) Cut_Results[3] = 1.0;
        if (n_good_jets >= 1) Cut_Results[4] = 1.0;
        if (PS_dipho_tag) Cut_Results[5] = 1.0;
        if (SLW_Tag) Cut_Results[6] = 1.0;
        if ( (PS_dipho_tag) and (SLW_Tag) ) Cut_Results[7] = 1.0;
        bool pass_cut = 0;
        double gen_value = 0;
        // Data 
        if (event.isRealData()){
          for (unsigned int i = 0; i < Cut_Results.size(); i++){
            pass_cut = Cut_Results[i]; // 1 = pass, 0 = fail 
            cutFlow->Fill(i,pass_cut);
          }
        }

        // MC 
        if (!event.isRealData()){
          Handle<GenEventInfoProduct> genInfo;
          event.getByToken(genInfoToken_, genInfo);
          genTotalWeight = genInfo->weight();
          for (unsigned int i = 0; i < Cut_Results.size(); i++){
            pass_cut = Cut_Results[i]; // 1 = pass, 0 = fail 
            gen_value = pass_cut*genTotalWeight;
            cutFlow->Fill(i,gen_value);
          }
        }

        // Create HHWWggCandidate Object 
        HHWWggCandidate HHWWgg(diphoVector, phoVector, electronVector, muonVector, theMET_, genParticlesVector, JetVector, Cut_Results, Cut_Variables);
        //HHWWggCandidate HHWWgg(diphoVector, phoVector, vertex_zero, genVertex, electronVector, muonVector, theMET_, genParticlesVector, JetVector, Cut_Results);
        HHWWggColl_->push_back(HHWWgg);

      event.put( std::move(HHWWggColl_) );
    }
  }
  typedef flashgg::HHWWggCandidateProducer FlashggHHWWggCandidateProducer;
  DEFINE_FWK_MODULE( FlashggHHWWggCandidateProducer );

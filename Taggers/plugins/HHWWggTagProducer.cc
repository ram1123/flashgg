// Abe Tishelman-Charny
// November 2019
// Derived from HH->WWgg event dumper and HH->bbgg tagger 

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
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/HHWWggTag.h"
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
  class HHWWggTagProducer : public EDProducer
  {
  public:
    //---typedef
    typedef math::XYZTLorentzVector LorentzVector;

    //---ctors
    // HHWWggTagProducer();
    HHWWggTagProducer( const ParameterSet & );

    //---Outtree 
    edm::Service<TFileService> fs;
    
    TH1F* numDiphoCand; 
    TH1F* diphoton_idx_h;
    TH1F* diPhotons_size_h;

  private:
    double genTotalWeight;
    void produce( Event &, const EventSetup & ) override;
    std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;
    std::string inputDiPhotonName_;

    std::string inputJetsName_;
    std::vector<std::string> inputJetsSuffixes_;
    unsigned int inputJetsCollSize_;

    // Adding Jets 
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;

    EDGetTokenT<View<Photon> > photonToken_;
    Handle<View<flashgg::Photon> > photons;

    // EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    // Handle<View<flashgg::DiPhotonCandidate> > diphotons;

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
    string systLabel_;
    edm::Handle<double>  rho;

    // std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > JetToken_;


    std::vector< std::string > systematicsLabels;
    std::vector<std::string> inputDiPhotonSuffixes_;

    //---ID selector
    ConsumesCollector cc_;
    GlobalVariablesComputer globalVariablesComputer_;
    // CutBasedDiPhotonObjectSelector idSelector_;

    //----output collection
    // auto_ptr<vector<HHWWggCandidate> > HHWWggColl_;

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
  // HHWWggTagProducer::HHWWggTagProducer( ):
  // photonToken_(),
  // diphotonToken_(),
  // genParticleToken_(),
  // electronToken_(),
  // muonToken_(),
  // METToken_(),
  // cc_( consumesCollector() )
  // // idSelector_( ParameterSet(), cc_ )

  // {}

    //---standard
    HHWWggTagProducer::HHWWggTagProducer( const ParameterSet & pSet):
    photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
    // diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
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
    systLabel_( pSet.getParameter<string> ( "SystLabel" ) ),
    cc_( consumesCollector() ), // need absence of comma on last entry 
    globalVariablesComputer_(pSet.getParameter<edm::ParameterSet>("globalVariables"), cc_)
    // idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

    {

      inputDiPhotonName_= pSet.getParameter<std::string > ( "DiPhotonName" );
      inputDiPhotonSuffixes_= pSet.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
      std::vector<edm::InputTag>  diPhotonTags;
      for (auto & suffix : inputDiPhotonSuffixes_){ 
          systematicsLabels.push_back(suffix);
          std::string inputName = inputDiPhotonName_;
          inputName.append(suffix);
          if (!suffix.empty()) diPhotonTags.push_back(edm::InputTag(inputName));
          else  diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_));
      }
      for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

      inputJetsName_= pSet.getParameter<std::string> ( "JetsName" );
      inputJetsCollSize_= pSet.getParameter<unsigned int> ( "JetsCollSize" );
      inputJetsSuffixes_= pSet.getParameter<std::vector<std::string> > ( "JetsSuffixes" );
      // cout << "inputJetsCollSize_ = " << inputJetsCollSize_ << endl;
      // std::vector<edm::InputTag>  jetTags; // With bregression on 
      for (auto & suffix : inputJetsSuffixes_) {
          if (!suffix.empty()) systematicsLabels.push_back(suffix);  //nominal is already put in the diphoton loop
          for (unsigned int i = 0; i < inputJetsCollSize_ ; i++) {
                std::string bregtag = suffix;
                bregtag.append(std::to_string(i));
                // jetTags.push_back(edm::InputTag(inputJetsName_,bregtag)); // With bregression on 
          }         
      }
      // for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); } // With bregression on 

      // Jets without bregression 
      auto jetTags = pSet.getParameter<std::vector<edm::InputTag> > ( "JetTags" ); 
      for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }

      genInfo_ = pSet.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );
      numDiphoCand = fs->make<TH1F> ("numDiphoCand","numDiphoCand",10,0,10); 
      diphoton_idx_h = fs->make<TH1F> ("diphoton_idx_h","diphoton_idx_h",20,0,20); 
      diPhotons_size_h = fs->make<TH1F> ("diPhotons_size_h","diPhotons_size_h",20,0,20); 
      // indexes = fs->make<TH1F> ("indexes","indexes",5,0,5);
      // numEvents = fs->make<TH1F> ("numEvents","numEvents",1,0,10);

      // gen_weights = fs->make<TH1F> ("gen_weights","gen_weights",1000,-2,2);
      // vars = fs->make<TH1F> ("vars","vars",10,0,10);
      // cutFlow = fs->make<TH1F> ("cutFlow","Cut Flow",10,0,10);
      // WTags = fs->make<TH1F> ("WTags","W Tags",3,0,3);

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

      // produces<vector<HHWWggTag>>();
      for (auto & systname : systematicsLabels) {
          produces<vector<HHWWggTag>>(systname);
      }
      produces<vector<TagTruthBase>>();
      // cout << "**************************** in HHWWggTagProducer.cc **********************************************" << endl;
    }

    void HHWWggTagProducer::produce( Event &event, const EventSetup & )
    {

      // update global variables
      globalVariablesComputer_.update(event);

      // Get particle objects
      event.getByToken( photonToken_, photons );
      // event.getByToken( diphotonToken_, diphotons );
      // event.getByToken( genParticleToken_, genParticle );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );
      event.getByToken( mvaResultToken_, mvaResults );
      event.getByToken( vertexToken_, vertices );
      event.getByToken( rhoTag_, rho);

      double rho_    = *rho;

      // Set cut booleans
      // std::vector<double> Cut_Results = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Cut_Results[i] = 1: Event Passed Cut i 
      std::vector<double> Cut_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i 
      std::vector<double> Vertex_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i 

      // Cut Variables 
      // double has_PS_Dipho = 0, pass_METfilters = 0, dipho_vertex_is_zero = 0, pass_leadPhoOverMassThreshold = 0, pass_subleadPhoOverMassThreshold = 0,
      //   pass_LeadPhoton_MVA = 0, pass_SubLeadPhoton_MVA = 0, pass_dipho_MVA = 0, number_passed_jetid = 0;
      // double dipho_vertex_is_zero = -999;
      // double SLW_Tag = 0.; // Semi-Leptonic W Tag  
      // double FLW_Tag = 0.; // Fully-Leptonic W Tag
      // double FHW_Tag = 0.; // Fully-Hadronic W Tag 
      // bool PS_dipho_tag = 0; // preselected diphoton 

      //---output collection
      // std::unique_ptr<vector<HHWWggCandidate> > HHWWggColl_( new vector<HHWWggCandidate> );
      // std::unique_ptr<vector<HHWWggTag> > tags( new vector<HHWWggTag> );
      int n_METs = METs->size(); // Should be 1, but using as a way to obtain met four vector 
      int n_good_electrons = 0;
      int n_good_muons = 0;
      int n_good_leptons = 0;
      int n_good_jets = 0;
      double dipho_MVA = -99;
      double lead_pho_Hgg_MVA = -99, sublead_pho_Hgg_MVA = -99;
      double CMS_hgg_mass = -99;

      // Saved Objects after selections
      std::vector<flashgg::Jet> tagJets_;
      std::vector<flashgg::Muon> goodMuons_;
      std::vector<flashgg::Electron> goodElectrons_; 
      std::vector<flashgg::Met> theMET_;
      std::vector<flashgg::DiPhotonCandidate> diphoVector_;
      reco::GenParticle::Point genVertex;

      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();

//-----------------------------------------------------------------------------------------------------------

      // Vertex variables
      // double gen_vertex_z = -999;
      // double hgg_vertex_z = -999;
      double zero_vertex_z = -999;
      // double vertex_diff_zeroeth = -999;
      // double vertex_diff_hgg = -999;
      // double num_vertices = -999;

      int diphoton_vertex_index;
      edm::Ptr<reco::Vertex> diphoton_vertex;
      edm::Ptr<reco::Vertex> zero_vertex;
      // num_vertices = (double)vertices->size();
      // cout << "vertices->size() = " << vertices->size() << endl;
      if (vertices->size() > 0){
        zero_vertex = vertices->ptrAt( 0 );
      }

      // MC truth
      TagTruthBase truth_obj;
      // double genMhh=0.;
      if( ! event.isRealData() ) {
          Handle<View<reco::GenParticle> > genParticles;
          std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
          event.getByToken( genParticleToken_, genParticles );
          reco::GenParticle::Point higgsVtx(0.,0.,0.);
          for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
              int pdgid = genParticles->ptrAt( genLoop )->pdgId(); 
              // if( pdgid == 25 || pdgid == 22 ) { // not so sure if this is correct for HHWWgg because of potential photons from hadronization 
              if( pdgid == 25 ) { 
                  higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                  // gen_vertex_z = higgsVtx.z();
                  break;
              }
          }
          for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
              edm::Ptr<reco::GenParticle> genPar = genParticles->ptrAt(genLoop);
              if (selHiggses.size()>1) break;
            if (genPar->pdgId()==25 && genPar->isHardProcess()){
                selHiggses.push_back(genPar);
            }   
          }
          // if (selHiggses.size()==2){
          //     TLorentzVector H1,H2;
          //     H1.SetPtEtaPhiE(selHiggses[0]->p4().pt(),selHiggses[0]->p4().eta(),selHiggses[0]->p4().phi(),selHiggses[0]->p4().energy());
          //     H2.SetPtEtaPhiE(selHiggses[1]->p4().pt(),selHiggses[1]->p4().eta(),selHiggses[1]->p4().phi(),selHiggses[1]->p4().energy());
          //     genMhh  = (H1+H2).M();
          // }
          truth_obj.setGenPV( higgsVtx );
          truths->push_back( truth_obj );
      }

      // // Get Gen vertex 
      // bool got_gen_vertex = 0;
      // if (! event.isRealData()){
      //   for( auto &part : *genParticle ) {
      //     if( ( part.pdgId() != 2212 || part.vertex().z() != 0.) && (!got_gen_vertex) ){
      //       genVertex = part.vertex();
      //       gen_vertex_z = genVertex.z();
      //       got_gen_vertex = 1;
      //       // cout << "Gen vertex z: " << genVertex.z() << endl;
      //     }
      //   }
      //   if (!got_gen_vertex){
      //     cout << "**********WARNING: Did not obtain non-zero GEN vertex from GEN particles" << endl;
      //   }
      // }

      // METfilters 
      // bool passMETfilters=1;
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
                                  // passMETfilters=0;  
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
                                  // passMETfilters=0;
                                  break;
                              }
                      }
          }

      // bool photonSelection = false;
      // double idmva1 = 0.;
      // double idmva2 = 0.;
      // bool checked_first = false; 
      // Pass_PS = false;
      // bool one_FH_dr = false;
      // bool one_FL_dr = false;
      double num_FL_dr = 0;
      double num_FH_dr = 0;
      float dr_ll = 0;
      // int n_ps_dpho = diphotons->size(); // number of preselected diphotons in event 

      // for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) { // look at all diphotons 

      // If there is one diphoton, and its vertex is not the zeroeth, recompute photon quantities relative to zeroeth vertex 
      // Then create recomputed diphoton object 

      // flashgg::Photon* lpho;
      // flashgg::Photon* slpho;
      // edm::Ptr<flashgg::Photon>
      // edm::Ptr<flashgg::Photon>
      
      // if ( (n_ps_dpho == 1)){
      //   edm::Ptr<flashgg::DiPhotonCandidate> unCorr_diphoton = diphotons->ptrAt( 0 );
      //   // cout << "photons->size() = " << photons->size() << endl;

      //   cout << "++++++++++++++++++++++++++++++++++++" << endl;
      //   for (unsigned int i = 0; i < photons->size(); i++){
      //     edm::Ptr<flashgg::Photon> pton = photons->ptrAt(i);
      //     // cout << "photons->ptrAt(" << i << ") = " << pton << endl;
      //   }
      //   cout << "++++++++++++++++++++++++++++++++++++" << endl;

      //   // edm::Ptr<flashgg::Photon> a = photons->ptrAt(0);
        



      // for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) { // only look at highest pt dipho
      //   edm::Ptr<flashgg::DiPhotonCandidate> dipho_ = diphotons->ptrAt( diphoIndex );
      //   diphoton_vertex_index = dipho_->vertexIndex();
      //   // cout << "vertex index = " << diphoton_vertex_index << endl;
      //   indexes->Fill(diphoton_vertex_index);
      //   // if (diphoton_vertex_index != 0){
      //   //   cout << "********************************************" << endl;
      //   //   cout << "********************************************" << endl;
      //   //   cout << "********************************************" << endl;
      //   //   cout << "diphoton vertex index not 0." << endl;
      //   //   cout << "Index: " << diphoton_vertex_index << endl;
      //   //   cout << "********************************************" << endl;
      //   //   cout << "********************************************" << endl;
      //   //   cout << "********************************************" << endl;
      //   // }
      // }

      // cout << "in HHWWggTagProducer.cc: Right before diphoton loop" << endl;

      // cout << "diphotons->size() = " << diphotons->size() << endl;

      // read diphotons
      for (unsigned int diphoton_idx = 0; diphoton_idx < diPhotonTokens_.size(); diphoton_idx++) { //looping over all diphoton systematics
        // cout << "diphoton_idx = " << diphoton_idx << endl;
        diphoton_idx_h->Fill(diphoton_idx);
        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        event.getByToken( diPhotonTokens_[diphoton_idx], diPhotons ); // for each diphoton systematic 

        // event.getByToken( diphotonToken_, diphotons );
        unsigned int loopOverJets = 1;
        if (inputDiPhotonSuffixes_[diphoton_idx].empty()) loopOverJets = inputJetsSuffixes_.size();
        for (unsigned int jet_col_idx = 0; jet_col_idx < loopOverJets; jet_col_idx++) {//looping over all jet systematics, only for nominal diphotons
          // cout << "jet_col_idx = " << jet_col_idx << endl;
          std::unique_ptr<vector<HHWWggTag> > tags( new vector<HHWWggTag> );

        diPhotons_size_h->Fill(diPhotons->size());
        if (diPhotons->size() > 0){ // for each systematic 
        // if (diphotons->size() > 0){
        for( unsigned int diphoIndex = 0; diphoIndex < 1; diphoIndex++ ) { // only look at highest pt dipho
          // std::unique_ptr<vector<HHWWggTag> > tags( new vector<HHWWggTag> );
          // edm::Ptr<flashgg::DiPhotonCandidate> dipho = Corrdiphoton; 
          edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex ); // systematic loop 
          // edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex ); 
          edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );    
          // cout << "dipho energy1 = " << dipho->genP4().E() << endl; 



          // MVA selections 
          // kinematic cuts on diphotons
          // auto leadPho = dipho->leadingPhoton();
          // auto subleadPho = dipho->subLeadingPhoton();



          // if (diphotons2->size() > 0){
          // edm::Ptr<flashgg::DiPhotonCandidate> dipho2 = diphotons2->ptrAt( diphoIndex ); 
          // // cout << "dipho energy2 = " << dipho2->genP4().E() << endl; 
          //   if ( dipho->genP4().mass() != dipho2->genP4().mass()){

          //     cout << "****************************************************************************************************************************************" << endl;
          //     cout << "different invariant masses" << endl;
          //     cout << "****************************************************************************************************************************************" << endl;
          //     cout << "****************************************************************************************************************************************" << endl;
          //     cout << "dipho mass1 = " << dipho->genP4().mass() << endl; 
          //     cout << "dipho mass2 = " << dipho2->genP4().mass() << endl; 
          //   }
          // }

          diphoton_vertex = dipho->vtx();        
          // diphoton_vertex_index = dipho->vertexIndex();
          // if (diphoton_vertex_index != 0){
          //   cout << "********************************************" << endl;
          //   cout << "diphoton vertex index not 0." << endl;
          //   cout << "Index: " << diphoton_vertex_index << endl;
          //   cout << "********************************************" << endl;
          // }
          // l_eta = dipho->leadingPhoton().superCluster()->eta();
          // l_r9 = dipho->leadingPhoton()->full5x5_r9();
          // sl_eta = dipho->leadingPhoton()->superCluster()->eta();
          // sl_r9 = dipho->leadingPhoton()->full5x5_r9();
          // Check that diphoton is preselected 
          // flashgg::DiPhotonCandidate * DPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho.get());

          // Pass_PS |= idSelector_(*initialDPPointer, event);
          // if (Pass_PS == 0){
          //   cout << "******************************************************************************************" << endl;
          //   cout << "Supposedly passed preselection but not idselector" << endl;
          //   cout << "******************************************************************************************" << endl;
          // }

          //if (!Pass_PS) continue;

          // else if ( Pass_PS ){ 
          //   // abs(dipho->leadingPhoton()->superCluster()->eta()) > 1.4;
          //   // (dipho->leadingPhoton()->full5x5_r9()) < 0.9;
          //   // abs(dipho->subLeadingPhoton()->superCluster()->eta());
          //   // (dipho->subLeadingPhoton()->full5x5_r9())< 0.9;
          //   //if ((abs(dipho->leadingPhoton()->superCluster()->eta()) > 1.4) && ((dipho->leadingPhoton()->full5x5_r9()) < 0.9) && (abs(dipho->subLeadingPhoton()->superCluster()->eta())) && ((dipho->subLeadingPhoton()->full5x5_r9())< 0.9))
          //   //{
          //   PS_dipho_tag = true;
          //   checked_first = true;
          //   //}
          // }

          hasGoodElec = false;
          hasGoodMuons = false;

          // Check MET Filters
          // if(passMETfilters){
          //   pass_METfilters = 1;
          // }

          // if(diphotons->ptrAt(diphoIndex)->vertexIndex()==0){
          //   dipho_vertex_is_zero = 1;
          // }

          // Check if useVertex0only_
          // if(useVertex0only_) // If only 0th vertex of diphotons is used, continue on all other vertex indices 
          //     if(diphotons->ptrAt(diphoIndex)->vertexIndex()==0){
          //       dipho_vertex_is_zero = 1;
          //     }
                      

          // leading/subleading photon pt 
          // if( dipho->leadingPhoton()->pt() > ( dipho->mass() )*leadPhoOverMassThreshold_ ){ 
          //     pass_leadPhoOverMassThreshold = 1;
          //   }
          // if( dipho->subLeadingPhoton()->pt() > ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { 
          //     pass_subleadPhoOverMassThreshold = 1;
          //   }

          // leading/subleading photon MVA
          // idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ); // can choose the Hgg MVA score with respect to a vertex 
          // idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ); // 
          // lead_pho_Hgg_MVA = idmva1;
          // sublead_pho_Hgg_MVA = idmva2;
          // if (idmva1 > PhoMVAThreshold_) pass_LeadPhoton_MVA = 1;
          // if (idmva2 > PhoMVAThreshold_) pass_SubLeadPhoton_MVA = 1;

          // if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) {
          //     pass_LeadPhoton_MVA = 1;
          //     // continue; 
          //    } // isn't this already applied in preselection? 

          // Diphoton MVA 
          // if ( mvares->result >= MVAThreshold_ ){
          //   pass_dipho_MVA = 1;
          // }
          // if( mvares->result < MVAThreshold_ ) { continue; }


          // dipho_MVA = mvares->result; 


          // cout << "diphomva = " << mvares->result << endl;
          // photonSelection = true;

          // Electrons 
          std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
                                                                            useElectronMVARecipe_,useElectronLooseID_,
                                                                            deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_,
                                                                            rho_, event.isRealData() );

          // Muons                                                                   
          std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, leptonPtThreshold_,
          muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );

          n_good_electrons = goodElectrons.size();
          n_good_muons = goodMuons.size();
          n_good_leptons = n_good_electrons + n_good_muons;
          hasGoodElec = ( goodElectrons.size() > 0 );
          hasGoodMuons = ( goodMuons.size() > 0 );
                
          // FL: Require at dr(l,l) > 0.4 
          // For which pairs of the >=2 good leptons should dr be greater than 4?

          if (hasGoodElec && hasGoodMuons){
            for (unsigned int ei = 0; ei < goodElectrons.size(); ei++){
              Ptr<flashgg::Electron> electron = goodElectrons[ei];
              for (unsigned int mi = 0; mi < goodMuons.size(); mi++){
                Ptr<flashgg::Muon> muon = goodMuons[mi];
                dr_ll = deltaR(electron->eta(), electron->phi(), muon->eta(), muon->phi()); 
                if (dr_ll > 0.4){ 
                  // one_FL_dr = true;
                  num_FL_dr += 1.0;
                  // break;
                }
              }
            }
          }

          else if (hasGoodElec && !hasGoodMuons){
            for (unsigned int ei = 0; ei < goodElectrons.size() - 1; ei++){ // the last electron cannot be the first one in the dR calculation 
              Ptr<flashgg::Electron> electroni = goodElectrons[ei];
              for (unsigned int ej = ei + 1; ej < goodElectrons.size(); ej++){
                Ptr<flashgg::Electron> electronj = goodElectrons[ej];
                dr_ll = deltaR(electroni->eta(), electroni->phi(), electronj->eta(), electronj->phi()); 
                if (dr_ll > 0.4){ 
                  // one_FL_dr = true;
                  num_FL_dr += 1.0;
                  // break;
                }
              }
            }
          }

          else if (!hasGoodElec && hasGoodMuons){
            for (unsigned int mi = 0; mi < goodMuons.size() - 1; mi++){
              Ptr<flashgg::Muon> muoni = goodMuons[mi];
              for (unsigned int mj = mi + 1; mj < goodMuons.size(); mj++){
                Ptr<flashgg::Muon> muonj = goodMuons[mj];
                dr_ll = deltaR(muoni->eta(), muoni->phi(), muonj->eta(), muonj->phi()); 
                if (dr_ll > 0.4){ 
                  // one_FL_dr = true;
                  num_FL_dr += 1.0;
                  // break;
                }
              }
            }
          }

          // Jets 
          unsigned int jetCollectionIndex = diPhotons->at( diphoIndex ).jetCollectionIndex(); // diphoton collection for each systematic 
          // unsigned int jetCollectionIndex = diphotons->at( diphoIndex ).jetCollectionIndex();
          size_t vtx = (size_t)dipho->jetCollectionIndex();
          edm::Handle<edm::View<flashgg::Jet> > Jets_;

          unsigned int jet_token_index = jet_col_idx*inputJetsCollSize_+vtx;
          // cout << "jet_token_index = " << jet_token_index << endl;
          // cout << "jetCollectionIndex = " << jetCollectionIndex << endl;

          // event.getByToken( jetTokens_[jetCollectionIndex], Jets_); // testing 
          event.getByToken( jetTokens_[jet_col_idx*inputJetsCollSize_+vtx], Jets_);  //take the corresponding vertex of current systematic //WORKS 
          // cout << "right after getbytoken jettokens jets" << endl;

          std::vector<edm::Ptr<Jet> > tagJets;

          // Jet Selections
          for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ ) 
              {
                  bool keepJet=true;
                  edm::Ptr<flashgg::Jet> thejet = Jets_->ptrAt( candIndex_outer );
                  // if(!thejet->passesJetID  ( flashgg::Tight2017 ) ) { continue; }
                  if(thejet->passesJetID  ( flashgg::Tight2017 ) ) {
                      // number_passed_jetid += 1;
                    }
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
          // if (n_good_jets == 1){
          //   edm::Ptr<flashgg::Jet> theOnejet = tagJets[0];
          //   j_mass = theOnejet->mass(); 
          //   // cout << "pt = " << theOnejet->pt() << endl;      
          // }

          // FH: Require at least one delta r pair with dr > 0.4 
          float dr_jj = 0;
          // bool one_FH_dr = false;
          if (tagJets.size() >= 2){
            for (unsigned int ji = 0; ji < tagJets.size() - 1; ji++){
              Ptr<flashgg::Jet> jeti = tagJets[ji];
              for (unsigned int jj = ji + 1; jj < tagJets.size(); jj++){
                Ptr<flashgg::Jet> jetj = tagJets[jj];
                dr_jj = deltaR(jeti->eta(), jeti->phi(), jetj->eta(), jetj->phi()); 
                if (dr_jj > 0.4){ 
                  // one_FH_dr = true;
                  num_FH_dr += 1.0;
                  // break; 
                } 
              }
            }
          }

          // MET 
          if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

          // Check W tags 
          // if ( (n_good_leptons == 1) && (theMET->getCorPt() >= 45) && (n_good_jets >= 1) ) SLW_Tag = 1.0; // Passed SL selections 
          // if ( (n_good_leptons >= 2) && (theMET->getCorPt() >= 70) && (n_good_jets < 2) && (one_FL_dr) ) FLW_Tag = 1.0; // Passed FL selections 
          // if ( (n_good_leptons == 0) && (theMET->getCorPt() < 45) && (n_good_jets >= 3) && (one_FH_dr) ) FHW_Tag = 1.0; // Passed FH selections 

          // how to go from pointer 'jet' to object '*thisJetPointer'
          //flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(jet.get());
          //JetVector.push_back(*thisJetPointer);

          // for (unsigned int i = 0; i < tagJets.size(); i++){
          //   auto tag_jet = tagJets[i]; 
          //   flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(tag_jet.get());
          //   tagJets_.push_back(*thisJetPointer);
          // }
          // for (unsigned int i = 0; i < goodElectrons.size(); i++){
          //   auto good_elec = goodElectrons[i];
          //   flashgg::Electron * thisElectronPointer = const_cast<flashgg::Electron *>(good_elec.get());
          //   goodElectrons_.push_back(*thisElectronPointer);
          // }
          // for (unsigned int i = 0; i < goodMuons.size(); i++){
          //   auto good_muon = goodMuons[i];
          //   flashgg::Muon * thisMuonPointer = const_cast<flashgg::Muon *>(good_muon.get());
          //   goodMuons_.push_back(*thisMuonPointer);
          // }
          // Store MET as one element in MET vector 
          // if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          // for( int METIndex = 0; METIndex < n_METs; METIndex++ )
          // {
          //   edm::Ptr<flashgg::Met> m_entry = METs->ptrAt( METIndex );
          //   flashgg::Met * thisMETPointer = const_cast<flashgg::Met *>(m_entry.get());
          //   theMET_.push_back(*thisMETPointer);
          // }

          // dipho
          edm::Ptr<flashgg::DiPhotonCandidate> dipho_ = dipho;
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(dipho_.get());
          diphoVector_.push_back(*thisDPPointer);

          //-- Tag object 
          if ( (n_good_leptons == 1) && ( (n_good_jets == 1) ) ){
          // if ( (n_good_leptons == 1) && ( (n_good_jets == 1) || (n_good_jets == 2) ) ){
          // if ( (n_good_leptons == 1) ){
            int catnum = 0;
            Ptr<flashgg::Jet> jet1 = tagJets[0];
            Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

            if (n_good_electrons == 1 && n_good_jets == 1){
              Ptr<flashgg::Electron> tag_electron = goodElectrons[0];
              HHWWggTag tag_obj(dipho, tag_electron, theMET, jet1);
              if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);
              tag_obj.setDiPhotonIndex( diphoIndex ); 
              tag_obj.setMVA( -0.9 );
              tag_obj.setCategoryNumber( catnum );
              // tag_obj.setEventNumber(event.id().event() );
              tags->push_back( tag_obj );
              if( ! event.isRealData() ) {
                tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              }
              // if (n_good_jets == 1){
              //   HHWWggTag tag_obj(dipho, tag_electron, theMET, jet1);
                
              //   if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              //   else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

              //   tag_obj.setDiPhotonIndex( diphoIndex ); 
              //   tag_obj.setMVA( -0.9 );
              //   tag_obj.setCategoryNumber( catnum );
              //   // tag_obj.setEventNumber(event.id().event() );
              //   tags->push_back( tag_obj );

              //   if( ! event.isRealData() ) {
              //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              //   }
              // }
              // if (n_good_jets == 2){

              //   Ptr<flashgg::Jet> jet2 = tagJets[1];
              //   HHWWggTag tag_obj(dipho, tag_electron, theMET, jet1, jet2); // electron, MET, jet1, jet2 
              //   if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              //   else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

              //   tag_obj.setDiPhotonIndex( diphoIndex );           
              //   tag_obj.setMVA( -0.9 );
              //   tag_obj.setCategoryNumber( catnum );
              //   // tag_obj.setEventNumber(event.id().event() );
              //   tags->push_back( tag_obj ); 

              //   if( ! event.isRealData() ) {
              //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              //   }  
              // }
            }
            
            if (n_good_muons == 1 && n_good_jets == 1){
              Ptr<flashgg::Muon> tag_muon = goodMuons[0];
              HHWWggTag tag_obj(dipho, tag_muon, theMET, jet1);
              tag_obj.setDiPhotonIndex( diphoIndex ); 
              tag_obj.setMVA( -0.9 );
              tag_obj.setCategoryNumber( catnum );
              // tag_obj.setEventNumber(event.id().event() );
              tags->push_back( tag_obj );   
              if( ! event.isRealData() ) {
                tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              }           
              // if (n_good_jets == 1){
              //   HHWWggTag tag_obj(dipho, tag_muon, theMET, jet1);
                
              //   if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              //   else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

              //   tag_obj.setDiPhotonIndex( diphoIndex ); 
              //   tag_obj.setMVA( -0.9 );
              //   tag_obj.setCategoryNumber( catnum );
              //   // tag_obj.setEventNumber(event.id().event() );
              //   tags->push_back( tag_obj );

              //   if( ! event.isRealData() ) {
              //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              //   }
              // }
              // else if (n_good_jets == 2){

              //   Ptr<flashgg::Jet> jet2 = tagJets[1];
              //   HHWWggTag tag_obj(dipho, tag_muon, theMET, jet1, jet2); // electron, MET, jet1, jet2 
              //   if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              //   else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

              //   tag_obj.setDiPhotonIndex( diphoIndex );           
              //   tag_obj.setMVA( -0.9 );
              //   tag_obj.setCategoryNumber( catnum );
              //   // tag_obj.setEventNumber(event.id().event() );
              //   tags->push_back( tag_obj ); 

              //   if( ! event.isRealData() ) {
              //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
              //   }  
              // }
            }

            // if (n_good_jets == 1){
            // // cout << "in HHWWggTagProducer.cc: saving tag object 1 jet" << endl;

            //   HHWWggTag tag_obj(dipho, tag_electron, theMET, jet1); // electron, MET, jet1, jet2 
            //   if (loopOverJets == 1) 
            //       tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
            //   else  
            //       tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);
            //   // cout << "HHWWgg diphoton candIndex = " << diphoIndex << endl;
            //   tag_obj.setDiPhotonIndex( diphoIndex ); 
            //   tag_obj.setMVA( -0.9 );
            //   tag_obj.setCategoryNumber( catnum );
            //   // tag_obj.setEventNumber(event.id().event() );
            //   tags->push_back( tag_obj );

            //   if( ! event.isRealData() ) {
            //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
            //   }
            // }
            // else if (n_good_jets == 2){
            // // cout << "in HHWWggTagProducer.cc: saving tag object 2 jets" << endl;

            //   Ptr<flashgg::Jet> jet2 = tagJets[1];
            //   HHWWggTag tag_obj(dipho, tag_electron, theMET, jet1, jet2); // electron, MET, jet1, jet2 
            //   if (loopOverJets == 1) 
            //       tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
            //   else  
            //       tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);
            //   // cout << "HHWWgg diphoton candIndex = " << diphoIndex << endl;
            //   tag_obj.setDiPhotonIndex( diphoIndex );           
            //   tag_obj.setMVA( -0.9 );
            //   tag_obj.setCategoryNumber( catnum );
            //   // tag_obj.setEventNumber(event.id().event() );
            //   tags->push_back( tag_obj ); 

            //   if( ! event.isRealData() ) {
            //     tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
            //   }        
            // }
            
          } // if ( (n_good_leptons == 1) && (n_good_jets >= 1) )

          // else{
          //   cout << "event does not have at least 1 good electron and at least 1 good jet " << endl;
          // }


        } // Diphoton loop 

        // cout << "Just left diphoton loop" << endl;

        // Compare diphoton vertex to gen vertex
        
        // if (! event.isRealData()){
        //   // hgg_vertex_z = diphoton_vertex->z();
        //   zero_vertex_z = zero_vertex->z();
        //   // cout << "hello" << endl;
        //   // cout << "genVertex.z() = " << genVertex.z() << endl;
        //   // cout << "diphoton_vertex->z() = " << diphoton_vertex->z() << endl;
        //   // vertex_diff_hgg = fabs(gen_vertex_z - hgg_vertex_z);
        //   vertex_diff_zeroeth = fabs(gen_vertex_z - zero_vertex_z);
        //   // cout << "vertex difference Hgg = " << vertex_diff_hgg << endl;
        //   // cout << "vertex difference Zero = " << vertex_diff_zeroeth << endl;
        // } 

        // }

          // tag_obj.setDiPhotonIndex( candIndex );
          // if (loopOverJets == 1) 
          //     tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
          // else  
          //     tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

          // if (tag_obj.dijet().mass()<mjjBoundaries_[0] || tag_obj.dijet().mass()>mjjBoundaries_[1]) continue;

          // // compute extra variables here
          // tag_obj.setMX( tag_obj.p4().mass() - tag_obj.dijet().mass() - tag_obj.diPhoton()->mass() + 250. );
          // tag_obj.setGenMhh( genMhh );
          // if (doReweight_>0) tag_obj.setBenchmarkReweight( reweight_values );
          
          // if(doSigmaMDecorr_){
          //     tag_obj.setSigmaMDecorrTransf(transfEBEB_,transfNotEBEB_);
          // }


          // // eval MVA discriminant
          // std::vector<float> mva_vector = mvaComputer_(tag_obj);
          // double mva = mva_vector[multiclassSignalIdx_];
          // if(doMVAFlattening_){
          //     double mvaScaled = mva/(mva*(1.-MVAscaling_)+MVAscaling_);
          //     mva = MVAFlatteningCumulative_->Eval(mvaScaled);
          // }

          // tag_obj.setEventNumber(evt.id().event() );
          // tag_obj.setMVA( mva );




      } // if at least 1 PS diphoton 
    // cout << "Just left if at least 1 PS diphoton" << endl;
      if (loopOverJets == 1){
        event.put( std::move( tags ),inputDiPhotonSuffixes_[diphoton_idx] );
        // cout << "inputDiPhotonSuffixes_[" << diphoton_idx << "] = " << inputDiPhotonSuffixes_[diphoton_idx] << endl; // Trying to debug HHWWgg 
      }
      else{  
        event.put( std::move( tags ),inputJetsSuffixes_[jet_col_idx] );
        // cout << "inputJetsSuffixes_[" << jet_col_idx << "] = " << inputJetsSuffixes_[jet_col_idx] << endl; // Trying to debug HHWWgg 
      }
    } // looping over jet systematics 
    // cout << "Just left looping over jet systematics" << endl;

    }  //looping over all diphoton systematics
    // cout << "Just left looping over all diphoton systematics" << endl;

      event.put( std::move( truths ) );

    } // HHWWggTagProducer::produce

  } // namespace flashgg

  typedef flashgg::HHWWggTagProducer FlashggHHWWggTagProducer;
  DEFINE_FWK_MODULE( FlashggHHWWggTagProducer );
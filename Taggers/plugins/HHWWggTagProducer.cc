// The purpose of this plugin is to tag HH->WWgg events

#include "stdlib.h"
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
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "flashgg/Taggers/interface/HHWWggTagProducer.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

namespace flashgg {

    void HHWWggTagProducer::produce( Event &event, const EventSetup & )
    {

      if (doHHWWggDebug_) cout << "[HHWWggTagProducer.cc - HHWWggTagProducer::produce] - systLabel: " << systLabel_ << endl;  

      // Set cut booleans
      std::vector<double> Cut_Variables(30,0.0); // Cut_Variables[i] = 1.0: Event Passed Cut i
      std::vector<double> MuonVars; // For saving Muon ID's and isolation
      std::vector<double> JetVars;

      // Cut Variables
      double pass_leadPhoOverMassThreshold = 0, pass_subleadPhoOverMassThreshold = 0;

      // Tag Variable
      bool FilledTag = 0; 
      string Tag_; 

      // Get particle objects
      event.getByToken( diphotonToken_, diphotons );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );
      event.getByToken( mvaResultToken_, mvaResults );
      event.getByToken( vertexToken_, vertices );

      //---output collection
      int n_good_electrons = 0;
      int n_good_muons = 0;
      int n_good_leptons = 0;
      int n_good_jets = 0;
      bool hasHighbTag = 0;
      float btagVal = 0;

      // Vertex 
      double GenVtx_z = -99; 
      double HggVtx_z = -99;
      double ZeroVtx_z = -99; 

      // Gen information 
      int pdgId; 
      std::vector<edm::Ptr<reco::GenParticle> > genHiggsBosons;
      std::vector<edm::Ptr<reco::GenParticle> > genWBosons;
      std::vector<edm::Ptr<reco::GenParticle> > genPhotons;    
      std::vector<edm::Ptr<reco::GenParticle> > genQuarks;    
      std::vector<edm::Ptr<reco::GenParticle> > genLeptons;    
      std::vector<edm::Ptr<reco::GenParticle> > genNeutrinos; 
      // edm::Ptr<reco::GenParticle& a;
      // edm::Ptr<reco::GenParticle& b;         

      // Misc
      int num_FL_dr = 0;
      int catnum = 3; // category number. default to untagged
      double DiphoCentralWeight = 1;
      // double prefireWeight = 1; // If you are setting prefire weight by hand 
      // Saved Objects after selections
      std::vector<edm::Ptr<flashgg::Jet> > allJets;
      std::vector<edm::Ptr<flashgg::Jet> > FHJets; // four jets for fully hadronic tag
      std::vector<edm::Ptr<flashgg::Electron> > allElectrons;
      std::vector<edm::Ptr<flashgg::Muon> > allMuons;
      reco::GenParticle::Point genVertex;

      // Four jets for the fully hadronic analysis. First two variables also used for SL final state jets
      Ptr<flashgg::Jet> jet1;
      Ptr<flashgg::Jet> jet2;
      Ptr<flashgg::Jet> jet3;
      Ptr<flashgg::Jet> jet4;
      bool keepJet = 1;
      bool passPUJetID = 1;

      // Saved Objects after selections
      std::vector<edm::Ptr<flashgg::Electron> > goodElectrons;
      std::vector<edm::Ptr<flashgg::Muon> > goodMuons; 

      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();

//-----------------------------------------------------------------------------------------------------------

      // Vertex variables
      edm::Ptr<reco::Vertex> diphoton_vertex;
      edm::Ptr<reco::Vertex> zero_vertex;

      if (vertices->size() > 0){
        zero_vertex = vertices->ptrAt( 0 );
        ZeroVtx_z = zero_vertex->z();
      }

      //-- MC truth
      TagTruthBase truth_obj;
      double genMhh = 0.;
      double genCosThetaStar_CS = 0.;   
      bool isData = event.isRealData();   
      if( ! isData ) {
          Handle<View<reco::GenParticle> > genParticles;
          std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
          event.getByToken( genParticleToken_, genParticles );

          edm::Handle<GenEventInfoProduct> genEvtInfo; 
          event.getByToken(genInfoToken_, genEvtInfo);

          // To obtain Higgs vertex information 
          reco::GenParticle::Point higgsVtx(0.,0.,0.);

          for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
               edm::Ptr<reco::GenParticle> genPar = genParticles->ptrAt(genLoop);
               if (selHiggses.size()>1) break;
              if (genPar->pdgId()==25 && genPar->isHardProcess()){
                  selHiggses.push_back(genPar); // to obtain diHiggs 
                  higgsVtx = genParticles->ptrAt( genLoop )->vertex(); // for dZ calculation
                  GenVtx_z = higgsVtx.z();
              } 
          }

          if (selHiggses.size()==2){
              TLorentzVector H1,H2;
              H1.SetPtEtaPhiE(selHiggses[0]->p4().pt(),selHiggses[0]->p4().eta(),selHiggses[0]->p4().phi(),selHiggses[0]->p4().energy());
              H2.SetPtEtaPhiE(selHiggses[1]->p4().pt(),selHiggses[1]->p4().eta(),selHiggses[1]->p4().phi(),selHiggses[1]->p4().energy());
              genMhh  = (H1+H2).M();
              genCosThetaStar_CS = getGenCosThetaStar_CS(H1,H2);   
          }

          // Save gen particles in nominal tree
          if(systLabel_ == ""){
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                edm::Ptr<reco::GenParticle> genPar = genParticles->ptrAt(genLoop);
                if(genPar->isHardProcess()){
                  pdgId = genPar->pdgId();
                  if(doHHWWggDebug_) std::cout << "==> Hard Process Gen Particle pdgID: " << pdgId << std::endl; 
                  if(pdgId==25){
                      genHiggsBosons.push_back(genPar);
                  }
                  else if(abs(pdgId)==24){
                      genWBosons.push_back(genPar);
                  }

                  else if(pdgId == 22){
                      genPhotons.push_back(genPar);
                  }

                  else if(abs(pdgId) >= 1 && abs(pdgId) <= 5){ // can check number of b quarks here 
                      genQuarks.push_back(genPar);
                  }

                  else if(abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15){ // can check number of taus here 
                      genLeptons.push_back(genPar);
                  }     

                  else if(abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16){ // can check number of tau neutrinos here 
                      genNeutrinos.push_back(genPar);
                  }                                    

                }
            }

            // PDF info 
            // GenEventInfoProduct::PDF const *pdf = genEvtInfo->pdf();
            // alternative PS event weights
            // std::vector<double> ps_weights = genEvtInfo->weights();
            // for (unsigned int i = 0; i < ps_weights.size(); i++){
              // std::cout << "ps_weights[" << i << "] = " << ps_weights[i] << std::endl;
            // }

            if(genHiggsBosons.size() >= 2) std::sort(genHiggsBosons.begin(), genHiggsBosons.end(), compPt);
            if(genWBosons.size() >= 2) std::sort(genWBosons.begin(), genWBosons.end(), compPt);
            if(genPhotons.size() >= 2) std::sort(genPhotons.begin(), genPhotons.end(), compPt);
            if(genQuarks.size() >= 2)    std::sort(genQuarks.begin(), genQuarks.end(), compPt);  
            if(genLeptons.size() >= 2)   std::sort(genLeptons.begin(), genLeptons.end(), compPt);   
            if(genNeutrinos.size() >= 2) std::sort(genNeutrinos.begin(), genNeutrinos.end(), compPt);    
                  
            if(doHHWWggDebug_){
              std::cout << "======================" << std::endl;
              std::cout << "==> N Gen Higgs: " << genHiggsBosons.size() << std::endl;
              std::cout << "==> N Gen W Bosons: " << genWBosons.size() << std::endl;
              std::cout << "==> N Gen Photons: " << genPhotons.size() << std::endl;
              std::cout << "==> N Gen Quarks: " << genQuarks.size() << std::endl;
              std::cout << "==> N Gen Leptons: " << genLeptons.size() << std::endl;
              std::cout << "==> N Gen Neutrinos: " << genNeutrinos.size() << std::endl;
            }  

          }          


          truth_obj.setGenPV( higgsVtx );
          truths->push_back( truth_obj );
      }

        std::unique_ptr<vector<HHWWggTag> > HHWWggtags( new vector<HHWWggTag> );

        //-- Loop through preselected diphotons. If more than one, it's possible the 2nd preselected diphoton will be tagged but not the 1st.
        for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ )
        { 

          //-- Get Diphoton and Diphoton MVA 
          edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex ); 
          edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

          //-- Get Photons
          // const flashgg::Photon* leadPho = dipho->leadingPhoton();
          // const flashgg::Photon* subleadPho = dipho->subLeadingPhoton();

           //-- Vertex Info 
          diphoton_vertex = dipho->vtx();
          if(HHWWgguseZeroVtx_) HggVtx_z = -99; 
          else HggVtx_z = diphoton_vertex->z();
          
          // leading/subleading photon pt
          if( dipho->leadingPhoton()->pt() > ( dipho->mass() )*leadPhoOverMassThreshold_ ){
              pass_leadPhoOverMassThreshold = 1;
            }
          if( dipho->subLeadingPhoton()->pt() > ( dipho->mass() )*subleadPhoOverMassThreshold_ ) {
              pass_subleadPhoOverMassThreshold = 1;
            }

          // // Doing cut flow analysis: Don't skip event, but check if photon selections are passed
          if(doHHWWggTagCutFlowAnalysis_){
            if( !pass_leadPhoOverMassThreshold || !pass_subleadPhoOverMassThreshold ) Cut_Variables[1] = 0.0; // failed photon selections
            else Cut_Variables[1] = 1.0; // passed photon selections
          }

          // // Not doing cut flow analysis: Skip event if photon selections are not passed
          else{
            if( !pass_leadPhoOverMassThreshold || !pass_subleadPhoOverMassThreshold ) continue; // Do not save event if leading and subleading photons don't pass pt/mgg cuts
          }

          hasGoodElec = false;
          hasGoodMuons = false;
          dipho_MVA = mvares->result;
          diPho_pT=dipho->pt();

          // Save information for all electrons and muons without selections 
          if(doHHWWggTagCutFlowAnalysis_ && systLabel_ == ""){
            for( unsigned int ei = 0; ei <  electrons->size() ; ei++ ){
              edm::Ptr<flashgg::Electron> theElectron = electrons->ptrAt( ei );
              allElectrons.push_back(theElectron);
            }
            for( unsigned int mi = 0; mi <  muons->size() ; mi++ ){
              edm::Ptr<flashgg::Muon> theMuon = muons->ptrAt( mi );
              allMuons.push_back(theMuon);
            }
          }

          // Save Good electron and muon vectors 

          // Electrons
          goodElectrons = LeptonSelection2018::selectElectrons( electrons->ptrs(), dipho, leptonPtThreshold_, electronEtaThresholds_, deltaRPhoElectronThreshold_,
                                                                deltaMassElectronZThreshold_, DeltaRTrkElec_, doHHWWggDebug_, ElectronID_ 
          );


          // Muons 
          goodMuons = LeptonSelection2018::selectMuons( muons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, muonEtaThreshold_, muPFIsoSumRelThreshold_, 
                                                        deltaRMuonPhoThreshold_, doHHWWggDebug_, MuonID_

          );

          // If doing cut flow analysis, save Muon IDs
          if(doHHWWggTagCutFlowAnalysis_ && systLabel_ == "") MuonVars = GetMuonVars(muons->ptrs(), vertices->ptrs());

          n_good_electrons = goodElectrons.size();
          n_good_muons = goodMuons.size();
          n_good_leptons = n_good_electrons + n_good_muons;
          hasGoodElec = ( goodElectrons.size() > 0 );
          hasGoodMuons = ( goodMuons.size() > 0 );

          // Jets
          unsigned int jetCollectionIndex = diphotons->at( diphoIndex ).jetCollectionIndex(); 
          edm::Handle<edm::View<flashgg::Jet> > Jets_;
          event.getByToken( jetTokens_[jetCollectionIndex], Jets_); 
          std::vector<edm::Ptr<Jet> > tagJets;

          // If doing cut flow analysis, save Jet IDs
          if(doHHWWggTagCutFlowAnalysis_ && systLabel_ == "") JetVars = GetJetVars(Jets_->ptrs(), dipho);

          // Jet Selections
          if(doHHWWggDebug_){
            cout << "applyPUJetID: " << applyPUJetID_ << endl; 
          }

          for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ )
              {
                  keepJet = true;
                  edm::Ptr<flashgg::Jet> thejet = Jets_->ptrAt( candIndex_outer );
                  if(systLabel_ == "") allJets.push_back(thejet);

                  // JetID selection
                  // Do not save jet in goodJets if it doesn't pass the recommended ID 
                  if( JetIDLevel_ == "Loose" && !thejet->passesJetID  ( flashgg::Loose ) ) continue;
                  if( JetIDLevel_ == "Tight" && !thejet->passesJetID  ( flashgg::Tight ) ) continue;
                  if( JetIDLevel_ == "Tight2017" && !thejet->passesJetID  ( flashgg::Tight2017 ) ) continue;
                  if( JetIDLevel_ == "Tight2018" && !thejet->passesJetID  ( flashgg::Tight2018 ) ) continue;                   

                  // Eta selection                 
                  if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }

                  // pT Selection
                  if( thejet->pt() < jetPtThreshold_ ) { continue; }

                  // dR Selections
                  float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                  float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                  dipho->subLeadingPhoton()->superCluster()->phi() );

                  if( dRPhoLeadJet < deltaRPhoLeadJet_ || dRPhoSubLeadJet < deltaRPhoSubLeadJet_ ) { continue; }

                  if(applyPUJetID_){
                    // PUJet ID Selection 
                    std::vector<std::pair<double,double> > eta_cuts_(4);
                    eta_cuts_[0] = std::make_pair (0    , 2.50 );
                    eta_cuts_[1] = std::make_pair (2.50 , 2.75 );
                    eta_cuts_[2] = std::make_pair (2.75 , 3.00 );
                    eta_cuts_[3] = std::make_pair (3.00 , 10);
                    
                    if ( (!pujid_wp_pt_bin_1_.empty())  &&
                        (!pujid_wp_pt_bin_2_.empty())  ){
                        passPUJetID = 0;
                        for (UInt_t eta_bin = 0; eta_bin < pujid_wp_pt_bin_1_.size(); eta_bin++ ){
                            if ( fabs( thejet->eta() ) >  eta_cuts_[eta_bin].first &&
                                fabs( thejet->eta() ) <= eta_cuts_[eta_bin].second){
                                if ( thejet->pt() >  20 &&
                                    thejet->pt() <= 30 && thejet->puJetIdMVA() > pujid_wp_pt_bin_1_[eta_bin] )
                                    passPUJetID = 1;
                                if ( thejet->pt() >  30 &&
                                    thejet->pt() <= 50 && thejet->puJetIdMVA() > pujid_wp_pt_bin_2_[eta_bin] )
                                    passPUJetID = 1;
                                if (thejet->pt() > 50) passPUJetID = 1;
                            }
                        }
                        if (!passPUJetID){
                          keepJet = false;
                        }  
                    }
                  }

                  if( hasGoodElec )
                      for( unsigned int electronIndex = 0; electronIndex < goodElectrons.size(); electronIndex++ )
                          {
                              Ptr<flashgg::Electron> electron = goodElectrons[electronIndex];
                              float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), electron->eta(), electron->phi() ) ;
                              if( dRJetElectron < deltaRJetMuonThreshold_ ){ keepJet = false; }
                          }

                  if( hasGoodMuons )
                      for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ )
                          {
                              Ptr<flashgg::Muon> muon = goodMuons[muonIndex];
                              float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;
                              if( dRJetMuon < deltaRJetMuonThreshold_ ) { keepJet = false; }
                          }

                  if(keepJet)
                      tagJets.push_back( thejet );

              }

          n_good_jets = tagJets.size(); 

          // MET
          if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

          //-- Categorize Events
          // In the semi-leptonic final state, because the categorization is done with MVAs, will just save all events orthogonal to the other two final states 
          if (doHHWWggDebug_) std::cout << "[HHWWggTagProducer] n_good_leptons = " << n_good_leptons << ",\t n_good_jets = " << n_good_jets << std::endl;

          // Three final state tags should be orthogonal because don't want to use same event in multiple background models because then it's incorrect to combine them 

          //-- Semi-Leptonic Final state tags
          if( (HHWWggAnalysisChannel_ == "SL" || HHWWggAnalysisChannel_ == "all") && FilledTag == 0)
          {
            Tag_ = "SL";
            if ( n_good_leptons == 1 ) {
              catnum = 0;
              if(doHHWWggDebug_) cout << "Filling Semi-Leptonic category..." << endl;
              HHWWggTag tag_obj;
              if(doHHWWggTagCutFlowAnalysis_){
                if(systLabel_ == ""){
                  HHWWggTag tag_obj_(dipho, allElectrons, goodElectrons, allMuons, goodMuons, theMET, allJets, tagJets, Cut_Variables, MuonVars, JetVars);                 
                  tag_obj = tag_obj_;
                } 
                else{
                  HHWWggTag tag_obj_(dipho, goodElectrons, goodMuons, theMET, tagJets, Cut_Variables); // If systematics tree, only save "good" objects (current analysis)
                  tag_obj = tag_obj_;
                }                 

                tag_obj.setGenVtx_z(GenVtx_z);
                tag_obj.setHggVtx_z(HggVtx_z);
                tag_obj.setZeroVtx_z(ZeroVtx_z);                                  

              }

              // If not doing cutflow analysis, save the minimum which is just the dipho information for the final fit, to keep process and output lightweight 
              else{
                HHWWggTag tag_obj_(dipho); // diphoton
                tag_obj = tag_obj_;
              } 

              // Set tag object attributes 
              tag_obj.setSystLabel( systLabel_);
              // tag_obj.setZeroVertex( zero_vertex );
              tag_obj.setDiPhotonIndex( diphoIndex );
              tag_obj.setCategoryNumber( catnum );
              tag_obj.setGenMhh( genMhh );
              tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS ); 
              DiphoCentralWeight = dipho->centralWeight();
              // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
              // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; // just used here to set diphocentral weight. Done separately for tag_obj central weight in SetCentralUpDownWeights()
              tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
              tag_obj.SetDiPhoMVA(dipho_MVA);
              tag_obj.SetDiPhoPt(diPho_pT);
              tag_obj = ComputePUJetIDs(tagJets, tag_obj);
              
              //-- Include Scale Factors 
              // Set CentralWeight values for each SF to access in trees 
              if(!isData){
                tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_);
                tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                // tag_obj.SavePDFInfo(pdf);
              } 

              // if(doHHWWggDebug_){
              //   cout << "*********************************************************" << endl; 
              //   cout << "HHWWgg Object Scale Factors:" << endl; 
              //   PrintScaleFactorsObj(tag_obj);
              // }

              // Push back tag object 
              HHWWggtags->push_back( tag_obj );
              FilledTag = 1; 
              if( ! isData ) {
                HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }

            } // if (n_good_leptons == 1)
          } // Semileptonic cateogries

          //-- Fully Hadronic Final State Tags
          if( (HHWWggAnalysisChannel_ == "FH" || HHWWggAnalysisChannel_ == "all") && (FilledTag == 0))
          {
            Tag_ = "FH";
            // double diPho_pT=dipho->pt();
            // if (n_good_leptons==0 && n_good_jets>=4 && diPho_pT > FH_Dipho_pT_Thre_)
            if (n_good_leptons==0 && n_good_jets>=4)
            {
              catnum = 1;
              if(doHHWWggDebug_) cout << "\n\nFilling Fully-Hadronic category..." << endl;
              // Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

              // Define four jets with WH min method, or just take four leading pT
              if(doHHWWggFHptOrdered_){
                if (doHHWWggDebug_) std::cout << "=============> doHHWWggFHptOrdered_ ==================" << std::endl;
                FHJets = GetFHPtOrderedJets(doHHWWggDebug_, tagJets);
                jet1 = FHJets[0];
                jet2 = FHJets[1];
                jet3 = FHJets[2];
                jet4 = FHJets[3];
              }
              else if (doHHWWggFHminWHJets_)
              {
                if (doHHWWggDebug_) std::cout << "=============> doHHWWggFHminWHJets_ ==================" << std::endl;
                FHJets = GetFHminWHJets(doHHWWggDebug_, tagJets);
                jet1 = FHJets[0];
                jet2 = FHJets[1];
                jet3 = FHJets[2];
                jet4 = FHJets[3];
              }
              else if (doHHWWggFHminWHLead2Jet_)
              {
                if (doHHWWggDebug_) std::cout << "=============> doHHWWggFHminWHLead2Jet_ ==================" << std::endl;
                FHJets = GetFHminWHLead2Jet(doHHWWggDebug_, tagJets);
                jet1 = FHJets[0];
                jet2 = FHJets[1];
                jet3 = FHJets[2];
                jet4 = FHJets[3];
              }
              else if (doHHWWggFHminHiggsMassOnly_) // kept doHHWWggFHminHiggsMassOnly_ to true in Taggers/python/flashggHHWWggTag_cfi.py; this ensures that this won't eneter in the last else
              {
                if (doHHWWggDebug_) std::cout << "=============> doHHWWggFHminHiggsMassOnly_ ==================" << std::endl;
                FHJets = GetFHminHiggsMassOnly(doHHWWggDebug_, tagJets);
                jet1 = FHJets[0];
                jet2 = FHJets[1];
                jet3 = FHJets[2];
                jet4 = FHJets[3];
              }
              else {
                std::cout << "You need to switch on any one hadronic selections..." << std::endl;
                std::cout << "==> HELP: doHHWWggFHptOrdered=1  or "                 << std::endl;
                std::cout << "          doHHWWggFHminWHJets=1  or "                 << std::endl;
                std::cout << "          doHHWWggFHminWHLead2Jet=1  or "             << std::endl;
                std::cout << "          doHHWWggFHminHiggsMassOnly=1 "              << std::endl;
                std::cout << "==> Please you any one options while running "        << std::endl;
                exit(0);
              }
              HHWWggTag tag_obj;

              if (doHHWWggTagCutFlowAnalysis_){
                // Save cut flow variables
                if ( (jet1->p4() + jet2->p4()).M() > 40 && (jet1->p4() + jet2->p4()).M() < 160 ) Cut_Variables[7] = 1.0;
                if ( (jet1->p4() + jet2->p4()).M() > 65 && (jet1->p4() + jet2->p4()).M() < 105 ) Cut_Variables[8] = 1.0;
                if ( (jet3->p4() + jet4->p4()).M() > 0 && (jet3->p4() + jet4->p4()).M() < 160 ) Cut_Variables[9] = 1.0;
                if ( (jet1->p4() + jet2->p4() + jet3->p4() + jet4->p4()).M() > 105 && (jet1->p4() + jet2->p4() + jet3->p4() + jet4->p4()).M() < 160 ) Cut_Variables[10] = 1.0;
                if ( (jet1->p4() + jet2->p4() + jet3->p4() + jet4->p4()).M() > 40 && (jet1->p4() + jet2->p4() + jet3->p4() + jet4->p4()).M() < 210 ) Cut_Variables[11] = 1.0;

                // save tag
                if(systLabel_ == ""){
                  HHWWggTag tag_obj_(dipho, theMET, jet1, jet2, jet3, jet4, allJets, tagJets, Cut_Variables, JetVars); 
                  tag_obj = tag_obj_;   
                } 
                else{
                  HHWWggTag tag_obj_(dipho, theMET, jet1, jet2, jet3, jet4, tagJets, Cut_Variables); // If systematics tree, only save "good" objects (current analysis)
                  tag_obj = tag_obj_;   
                } 
                tag_obj.setGenVtx_z(GenVtx_z);
                tag_obj.setHggVtx_z(HggVtx_z);
                tag_obj.setZeroVtx_z(ZeroVtx_z);                 
                            
              }
              else{
                HHWWggTag tag_obj_(dipho); // diphoton, electron, MET, jet1, jet2
                tag_obj = tag_obj_;
              }
              tag_obj.setSystLabel( systLabel_);
              // tag_obj.setZeroVertex( zero_vertex );
              tag_obj.setDiPhotonIndex( diphoIndex );
              tag_obj.setCategoryNumber( catnum );
              tag_obj.setGenMhh( genMhh );
              tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS );      

              //-- Include Scale Factors 
              DiphoCentralWeight = dipho->centralWeight();
              // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
              // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; 
              tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
              tag_obj.SetDiPhoMVA(dipho_MVA);
              tag_obj.SetDiPhoPt(diPho_pT);

              tag_obj = ComputePUJetIDs(tagJets, tag_obj);
                 
              //-- Include Scale Factors 
              // Set CentralWeight values for each SF to access in trees   
              if(!isData){
                tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_);  
                tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                // tag_obj.SavePDFInfo(pdf);
              } 

              HHWWggtags->push_back( tag_obj );
              FilledTag = 1;
              if( ! isData ) {
                HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }
            }  // if (n_good_leptons==0 && n_good_jets>=4)
          } // Fully Hadronic Categories

          //-- Fully-Leptonic Final State 

          if( (HHWWggAnalysisChannel_ == "FL" || HHWWggAnalysisChannel_ == "all") && FilledTag == 0)
          {  

            Tag_ = "FL";

            // For the Fully-Leptonic Final State, apply the Medium bVeto working point b/c retain high sig eff, remove vast ttHToGG
            // If jet collection has a jet suspected to be a bjet, don't save the event
            hasHighbTag = 0; // collection has jet with high b score 
            for (unsigned int j = 0; j < tagJets.size(); j++){
              Ptr<flashgg::Jet> jet_ = tagJets[j];
              btagVal = 0; 
              for(unsigned int BTagTypes_i = 0; BTagTypes_i < BTagTypes_.size(); BTagTypes_i ++){
                btagVal += jet_->bDiscriminator(BTagTypes_[BTagTypes_i]);
              }
              if (btagVal > btagThresh_){
                hasHighbTag = 1;
                if(doHHWWggDebug_){
                  cout << "Jet with score: " << btagVal << " is considered HIGH" << endl; 
                }
              } 
            }

            if( n_good_leptons >= 2 && !hasHighbTag ){

            num_FL_dr = GetNumFLDR(goodElectrons, goodMuons, deltaRLeps_);

            // if ( (n_good_leptons >=2 ) && (theMET->getCorPt() >= FL_METThreshold_) && num_FL_dr>=1 && (leadPho->p4().pt()+subleadPho->p4().pt())>0){
            if ( (n_good_leptons >=2 ) && (theMET->getCorPt() >= FL_METThreshold_) && num_FL_dr>=1){
              if(doHHWWggDebug_) cout << "Filling Fully-Leptonic category..." << endl;
              catnum = 2; 
              double Save=0.;
              double MassT;
              double MassT_l2;
              // double lep3pt = 0.;
              Ptr<flashgg::Electron> tag_electron1;
              Ptr<flashgg::Electron> tag_electron3;
              Ptr<flashgg::Electron> tag_electron2;
              Ptr<flashgg::Muon> tag_muon1;
              Ptr<flashgg::Muon> tag_muon3;
              Ptr<flashgg::Muon> tag_muon2;
              sort(goodElectrons.begin(),goodElectrons.end(),compEle);
              sort(goodMuons.begin(),goodMuons.end(),compMu);
              if (n_good_electrons >= 2){
                //std::cout<<"ele1:"<<goodElectrons[0]->pt()<<"ele2:"<<goodElectrons[1]->pt()<<"ele3:"<<goodElectrons[2]->pt()<<std::endl;
                if(n_good_electrons == 2 && n_good_muons < 1){
                tag_electron1 = goodElectrons[0];
                tag_electron2 = goodElectrons[1];
                Save = 1.;
                }
                else if(n_good_electrons > 2 && n_good_muons < 1){
                tag_electron1 = goodElectrons[0];
                tag_electron2 = goodElectrons[1];
                tag_electron3 = goodElectrons[2];
                Save = 1.;
                }
                else if(n_good_electrons ==2 && n_good_muons >=1){
                if(goodMuons[0]->pt()<lep3ptThre_){
                tag_electron1 = goodElectrons[0];
                tag_electron2 = goodElectrons[1];
                Save=1.;
                }
                else {
                  tag_electron1 = goodElectrons[0];
                  tag_electron2 = goodElectrons[1];
                  Save=0.;}
                }
                else if(n_good_electrons >2 && n_good_muons >=1){
                if(goodMuons[0]->pt()<lep3ptThre_){
                tag_electron1 = goodElectrons[0];
                tag_electron2 = goodElectrons[1];
                tag_electron3 = goodElectrons[2];
                Save = 1.;
                }
                else{
                tag_electron1 = goodElectrons[0];
                tag_electron2 = goodElectrons[1];
                tag_electron3 = goodElectrons[2];
                Save = 0.;
                }
                }
                reco::Candidate::LorentzVector Ele1LV = tag_electron1->p4();
                reco::Candidate::LorentzVector Ele2LV = tag_electron2->p4();
                reco::Candidate::LorentzVector DiEleLV = Ele1LV +  Ele2LV;
                MassT= sqrt(2*DiEleLV.Pt()*theMET->pt()*(1-cos(abs(DiEleLV.Phi()-theMET->phi()))));
                MassT_l2 = sqrt(2*Ele2LV.Pt()*theMET->pt()*(1-cos(abs(Ele2LV.Phi()-theMET->phi()))));
                  if((tag_electron1->charge()*tag_electron2->charge()==-1) && tag_electron1->pt()>=lep1ptThre_ && tag_electron2->pt()>=lep1ptThre_  && DiEleLV.Pt()>DiLepPtThre_ && (DiEleLV.M() > DiLepMassHigThre_ || DiEleLV.M() < DiLepMassLowThre_ )  && MassT>MassTThre_ && MassT_l2 >MassT_l2Thre_ && Save == 1.){
                    HHWWggTag tag_obj;

                    if(doHHWWggTagCutFlowAnalysis_){
                      Cut_Variables[18]=0.;
                      HHWWggTag tag_obj_(dipho, tag_electron1, tag_electron2, theMET, Cut_Variables, dipho_MVA); 
                      tag_obj = tag_obj_;    
                      tag_obj.SetGoodJets(tagJets); 
                      tag_obj.SetGoodElectrons(goodElectrons); 
                      tag_obj.SetGoodMuons(goodMuons); 
                      if(systLabel_ == ""){
                        tag_obj.SetAllJets(allJets);
                        tag_obj.SetAllElectrons(allElectrons);
                        tag_obj.SetAllMuons(allMuons);   
                      }                                          
                      tag_obj.setGenVtx_z(GenVtx_z);
                      tag_obj.setHggVtx_z(HggVtx_z);
                      tag_obj.setZeroVtx_z(ZeroVtx_z);          
                    }

                    else{
                      HHWWggTag tag_obj_(dipho);
                      tag_obj = tag_obj_;
                    }

                  tag_obj.setSystLabel(systLabel_);
                  // tag_obj.setZeroVertex( zero_vertex );
                  tag_obj.setDiPhotonIndex( diphoIndex );
                  tag_obj.setMVA( -0.9 );
                  tag_obj.setCategoryNumber( catnum );
                  DiphoCentralWeight = dipho->centralWeight();
                  // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
                  // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; 
                  tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
                  tag_obj.SetDiPhoMVA(dipho_MVA);
                  tag_obj.SetDiPhoPt(diPho_pT);
                  //-- Include Scale Factors 
                  // Set CentralWeight values for each SF to access in trees   
                  if(!isData){
                    tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_);
                    tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                    // tag_obj.SavePDFInfo(pdf);
                  } 

                  tag_obj.setGenMhh( genMhh );
                  tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS );  
                  tag_obj = ComputePUJetIDs(tagJets, tag_obj);                
                  HHWWggtags->push_back( tag_obj );
                  FilledTag = 1;
                  if( ! isData ) {
                    HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
                  }

                }
              }

              if (n_good_muons >= 2  ){
                if (n_good_muons == 2 && n_good_electrons <1){
                  tag_muon1 = goodMuons[0];
                  tag_muon2 = goodMuons[1];
                  Save =1.;
                }
                else if(n_good_muons >2 && n_good_electrons <1 ){
                tag_muon1 = goodMuons[0];
                tag_muon2 = goodMuons[1];
                tag_muon3 = goodMuons[2];
                Save = 1.;
                }
                else if(n_good_muons ==2 && n_good_electrons >=1){
                  if(goodElectrons[0]->pt()<lep3ptThre_){
                  tag_muon1 = goodMuons[0];
                  tag_muon2 = goodMuons[1];
                  Save = 1.;
                }
                  else{
                    tag_muon1 = goodMuons[0];
                    tag_muon2 = goodMuons[1];
                  Save = 0.;}
                }

                else if(n_good_muons >2 && n_good_electrons >=1){
                if(goodElectrons[0]->pt()<lep3ptThre_){
                tag_muon1 = goodMuons[0];
                tag_muon2 = goodMuons[1];
                tag_muon3 = goodMuons[2];
                Save = 1.;
                }
                else{
                  tag_muon1 = goodMuons[0];
                  tag_muon2 = goodMuons[1];
                  tag_muon3 = goodMuons[2];
                  Save = 0.;
                  }
                }
                reco::Candidate::LorentzVector Mu1LV = tag_muon1->p4();
                reco::Candidate::LorentzVector Mu2LV = tag_muon2->p4();
                reco::Candidate::LorentzVector DiMuLV = Mu1LV +  Mu2LV;
                MassT= sqrt(2*DiMuLV.Pt()*theMET->pt()*(1-cos(abs(DiMuLV.Phi()-theMET->phi()))));
                MassT_l2 = sqrt(2*Mu2LV.Pt()*theMET->pt()*(1-cos(abs(Mu2LV.Phi()-theMET->phi()))));
                if((tag_muon1->charge()*tag_muon2->charge() == -1) && tag_muon1->pt()>=lep1ptThre_ && tag_muon2->pt()>=lep2ptThre_  && DiMuLV.Pt()>DiLepPtThre_ && (DiMuLV.M() > DiLepMassHigThre_ || DiMuLV.M() < DiLepMassLowThre_ ) && MassT>MassTThre_ && MassT_l2>MassT_l2Thre_ && Save == 1.){
                  HHWWggTag tag_obj;

                  if(doHHWWggTagCutFlowAnalysis_){
                    Cut_Variables[18]=1.;
                    HHWWggTag tag_obj_(dipho, tag_muon1, tag_muon2, theMET, Cut_Variables, dipho_MVA);
                    tag_obj = tag_obj_;
                    tag_obj.SetGoodJets(tagJets); 
                    tag_obj.SetGoodElectrons(goodElectrons); 
                    tag_obj.SetGoodMuons(goodMuons); 
                    if(systLabel_ == ""){
                      tag_obj.SetAllJets(allJets);
                      tag_obj.SetAllElectrons(allElectrons);
                      tag_obj.SetAllMuons(allMuons);   
                    }                                        
                    tag_obj.setGenVtx_z(GenVtx_z);
                    tag_obj.setHggVtx_z(HggVtx_z);
                    tag_obj.setZeroVtx_z(ZeroVtx_z);        
                  }
                  else{
                    HHWWggTag tag_obj_(dipho);
                    tag_obj = tag_obj_;
                  }
                  tag_obj.setSystLabel(systLabel_);
                  // tag_obj.setZeroVertex( zero_vertex );
                  tag_obj.setDiPhotonIndex( diphoIndex );
                  tag_obj.setMVA( -0.9 );
                  tag_obj.setCategoryNumber( catnum );
                  //-- Include Scale Factors 
                  DiphoCentralWeight = dipho->centralWeight();
                  // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
                  // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; 
                  tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
                  tag_obj.SetDiPhoMVA(dipho_MVA);
                  tag_obj.SetDiPhoPt(diPho_pT);
                  //-- Include Scale Factors 
                  // Set CentralWeight values for each SF to access in trees 
                  if(!isData){
                    tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_);
                    tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                    // tag_obj.SavePDFInfo(pdf);
                  } 

                  tag_obj.setGenMhh( genMhh );
                  tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS );    
                  tag_obj = ComputePUJetIDs(tagJets, tag_obj);              
                  HHWWggtags->push_back( tag_obj );
                  FilledTag = 1;
                  if( ! isData ) {
                    HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
                  }
              }
            }

              if (n_good_muons >=1 && n_good_electrons >= 1){
              if (n_good_electrons == 1 && n_good_muons == 1){
                tag_electron1 = goodElectrons[0];
                tag_muon1 = goodMuons[0];
                Save=1.;
                }
              else if (n_good_electrons>1 && n_good_muons ==1){
              if(goodElectrons[1]->pt()<lep3ptThre_){
                tag_electron1 = goodElectrons[0];
                tag_muon1 = goodMuons[0];
                Save=1.;
              }
              else{
              tag_electron1 = goodElectrons[0];
              tag_muon1 = goodMuons[0];
              Save=0.;
              }
              }
              else if(n_good_electrons == 1 && n_good_muons > 1){
              if(goodMuons[1]->pt()<lep3ptThre_){
              tag_electron1 = goodElectrons[0];
              tag_muon1 = goodMuons[0];
              Save=1.;
              }
              else{
              tag_electron1 = goodElectrons[0];
              tag_muon1 = goodMuons[0];
              Save=0.;}
              }
              else if(n_good_electrons > 1 && n_good_muons > 1){
              if(goodMuons[1]->pt()<lep3ptThre_ && goodElectrons[1]->pt()<lep3ptThre_){
              tag_electron1 = goodElectrons[0];
              tag_muon1 = goodMuons[0];
              Save=1.;
              }
              else{tag_electron1 = goodElectrons[0];
                tag_muon1 = goodMuons[0];
                Save=0.;}
              }
              reco::Candidate::LorentzVector MuLV = tag_muon1->p4();
                reco::Candidate::LorentzVector EleLV = tag_electron1->p4();
                reco::Candidate::LorentzVector DiLepLV = EleLV +  MuLV;
                MassT= sqrt(2*DiLepLV.Pt()*theMET->pt()*(1-cos(abs(DiLepLV.Phi()-theMET->phi()))));
                MassT_l2= sqrt(2*MuLV.Pt()*theMET->pt()*(1-cos(abs(MuLV.Phi()-theMET->phi()))));
                if((((tag_electron1->pt()>=lep1ptThre_) && (tag_muon1->pt()>=lep2ptThre_))||((tag_muon1->pt()>=lep1ptThre_) && (tag_electron1->pt()>=lep2ptThre_))) && (tag_muon1->charge()*tag_electron1->charge()==-1) && (DiLepLV.M() > DiLepMassHigThre_ || DiLepLV.M() < DiLepMassLowThre_ )  && (DiLepLV.Pt()>DiLepPtThre_) && (MassT_l2>MassTThre_) && (MassT>MassT_l2Thre_) && (Save==1.)){
                  HHWWggTag tag_obj;
                  if(doHHWWggTagCutFlowAnalysis_){
                    if(tag_electron1->pt()>tag_muon1->pt()) Cut_Variables[18]=2.;//e mu
                    else Cut_Variables[18]=3.; //mu e
                    HHWWggTag tag_obj_(dipho, tag_electron1, tag_muon1, theMET, Cut_Variables, dipho_MVA);
                    tag_obj = tag_obj_; 
                    tag_obj.SetGoodJets(tagJets); 
                    tag_obj.SetGoodElectrons(goodElectrons); 
                    tag_obj.SetGoodMuons(goodMuons); 
                    if(systLabel_ == ""){
                      tag_obj.SetAllJets(allJets);
                      tag_obj.SetAllElectrons(allElectrons);
                      tag_obj.SetAllMuons(allMuons);   
                    }                                                            
                    tag_obj.setGenVtx_z(GenVtx_z);
                    tag_obj.setHggVtx_z(HggVtx_z);
                    tag_obj.setZeroVtx_z(ZeroVtx_z);             
                  }
                  else{
                    HHWWggTag tag_obj_(dipho);
                    tag_obj = tag_obj_;
                  }
                  tag_obj.setSystLabel(systLabel_);
                  // tag_obj.setZeroVertex( zero_vertex );
                  tag_obj.setDiPhotonIndex( diphoIndex );
                  tag_obj.setMVA( -0.9 );
                  tag_obj.setCategoryNumber( catnum );
                  //-- Include Scale Factors 
                  DiphoCentralWeight = dipho->centralWeight();
                  // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
                  // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; 
                  tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
                  tag_obj.SetDiPhoMVA(dipho_MVA);
                  tag_obj.SetDiPhoPt(diPho_pT);                    
                  //-- Include Scale Factors 
                  // Set CentralWeight values for each SF to access in trees 
                  if(!isData){
                    tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_); 
                    tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                    // tag_obj.SavePDFInfo(pdf);
                  }   

                  tag_obj.setGenMhh( genMhh );
                  tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS );   
                  tag_obj = ComputePUJetIDs(tagJets, tag_obj);               
                  HHWWggtags->push_back( tag_obj );
                  FilledTag = 1;
                  if( ! isData ) {
                    HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
                  }
              }
              }
              } // if n_good_leptons >= 2 and other FL cat selections met
            } // if n_good_leptons >= 2
          } // FL selection

          // Untagged category
          // Only fill this tag for the nominal tree. Do not fill for the systematics trees. 
          if(FilledTag == 0 && FillUntagged_ && systLabel_ == ""){
            Tag_ = "Untagged";
            if(doHHWWggTagCutFlowAnalysis_){
              if(doHHWWggDebug_) cout << "Filling untagged category..." << endl;
              catnum = 3;
              HHWWggTag tag_obj(dipho, allElectrons, goodElectrons, allMuons, goodMuons, theMET, allJets, tagJets, Cut_Variables, MuonVars, JetVars);
              tag_obj.setSystLabel(systLabel_);
              // tag_obj.setZeroVertex( zero_vertex );
              tag_obj.setDiPhotonIndex( diphoIndex );
              tag_obj.setCategoryNumber( catnum ); // Untagged category. Does not meet any selection criteria but want to save event
              tag_obj.setGenMhh( genMhh );
              tag_obj.setGenCosThetaStar_CS( genCosThetaStar_CS ); 
              tag_obj.setGenVtx_z(GenVtx_z);
              tag_obj.setHggVtx_z(HggVtx_z);
              tag_obj.setZeroVtx_z(ZeroVtx_z);                

              //-- Include Scale Factors 
              DiphoCentralWeight = dipho->centralWeight();
              // prefireWeight = dipho->weight("prefireWeightCentral"); // if setting pre fire weight by hand 
              // DiphoCentralWeight = DiphoCentralWeight * prefireWeight; 
              tag_obj.SetDiphoCentralWeight(DiphoCentralWeight);
              tag_obj.SetDiPhoMVA(dipho_MVA);
              tag_obj.SetDiPhoPt(diPho_pT);      

              //-- Include Scale Factors 
              // Set CentralWeight values for each SF to access in trees 
              if(!isData){
                tag_obj = SetCentralUpDownWeights(tag_obj, goodElectrons, goodMuons, tagJets, dipho, doHHWWggDebug_, MuonID_, muPFIsoSumRelThreshold_, Tag_);  
                tag_obj.SetGenObjs(genHiggsBosons, genWBosons, genPhotons, genQuarks, genLeptons, genNeutrinos);
                // tag_obj.SavePDFInfo(pdf);
              } 
              tag_obj = ComputePUJetIDs(tagJets, tag_obj);

              // if(doHHWWggDebug_){
              //   cout << "*********************************************************" << endl; 
              //   cout << "HHWWgg Object Scale Factors:" << endl; 
              //   PrintScaleFactorsObj(tag_obj);
              // }

              // Push back tag object 
              HHWWggtags->push_back( tag_obj );
              FilledTag = 1;  // actually need to configure this so that untagged isn't filled unless other preselected diphotons checked 
              if( ! isData ) {
                HHWWggtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }
            } // doHHWWggTagCutFlowAnalysis_
          } // Untagged category: if FilledTag == 0
        } // Preselected Diphoton loop 
      event.put( std::move( HHWWggtags ) );
      event.put( std::move( truths ) );
    } // HHWWggTagProducer::produce

  } // namespace flashgg

  typedef flashgg::HHWWggTagProducer FlashggHHWWggTagProducer;
  DEFINE_FWK_MODULE( FlashggHHWWggTagProducer );

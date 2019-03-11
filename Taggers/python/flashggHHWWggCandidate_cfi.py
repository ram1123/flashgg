import FWCore.ParameterSet.Config as cms
from HHWWggCandidateDumper_cfi import HHWWggCandidateDumper
from flashggHHWWggTag_cfi import flashggHHWWggTag

# cfi = configuration fragment include
# Clone these params into _cfg 
FlashggHHWWggCandidate = cms.EDProducer("FlashggHHWWggCandidateProducer", 
                                     PhotonTag              = cms.InputTag('flashggRandomizedPhotons'),
                                     DiPhotonTag            = cms.InputTag('flashggDiPhotons'),
                                     VertexTag              = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     GenParticleTag         = cms.InputTag('flashggPrunedGenParticles'),
                                     ElectronTag            = cms.InputTag('flashggSelectedElectrons'),
                                     #METTag                 = cms.InputTag('slimmedMETs'),
                                     MuonTag                = cms.InputTag('flashggSelectedMuons'),
                                     METTag                 = cms.InputTag('flashggMets'),
                                     # METTag                 = cms.InputTag('flashggMetsCorr'),
                                     JetTag                 = cms.InputTag('JetCollectionVInputTag')     
                                     # flashggSelectedJets
                                     # slimmedJets
                                     # JetCollectionVInputTag
                                     # UnpackedJetCollectionVInputTag

                                     )
flashggHHWWggTagSequence = cms.Sequence( flashggHHWWggTag )

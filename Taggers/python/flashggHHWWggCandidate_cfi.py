import FWCore.ParameterSet.Config as cms
from HHWWggCandidateDumper_cfi import HHWWggCandidateDumper
from flashggHHWWggTag_cfi import flashggHHWWggTag


FlashggHHWWggCandidate = cms.EDProducer("FlashggHHWWggCandidateProducer",
                                     PhotonTag              = cms.InputTag('flashggRandomizedPhotons'),
                                     DiPhotonTag            = cms.InputTag('flashggDiPhotons'),
                                     VertexTag              = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     GenParticleTag         = cms.InputTag('flashggPrunedGenParticles')
                                     )
flashggHHWWggTagSequence = cms.Sequence( flashggHHWWggTag )

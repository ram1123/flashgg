import FWCore.ParameterSet.Config as cms
#from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, flashggDeepCSV, maxJetCollections
from HHWWggCandidateDumper_cfi import HHWWggCandidateDumper
from flashggHHWWggTag_cfi import flashggHHWWggTag

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

# flashggUnpackedJetsaa = cms.EDProducer("FlashggVectorVectorJetUnpacker",
#                                     JetsTag = cms.InputTag("flashggFinalJets"),
#                                     #inputTagJets = cms.InputTag("flashggFinalJets"),
#                                     NCollections = cms.uint32(maxJetCollections)
#                                     )

# HTXSInputTags = cms.PSet(stage0cat = cms.InputTag("rivetProducerHTXS","stage0cat"), #2016
#                          stage1cat = cms.InputTag("rivetProducerHTXS","stage1cat"), #2016
#                          njets     = cms.InputTag("rivetProducerHTXS","njets"), #2016
#                          pTH       = cms.InputTag("rivetProducerHTXS","pTH"), #2016
#                          pTV       = cms.InputTag("rivetProducerHTXS","pTV"), #2016
#                          ClassificationObj = cms.InputTag("rivetProducerHTXS","HiggsClassification") # 2017
#                          )

# UnpackedJetCollectionVInputTag = cms.VInputTag()
# for i in range(0,maxJetCollections):
#         print 'i = ',i
#         UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

#inputTagJets= UnpackedJetCollectionVInputTag

# cfi = configuration fragment include
# Clone these params into _cfg 
FlashggHHWWggCandidate = cms.EDProducer("FlashggHHWWggCandidateProducer", 
                                     PhotonTag              = cms.InputTag('flashggRandomizedPhotons'),
                                     DiPhotonTag            = cms.InputTag('flashggDiPhotons'),
                                     #DiPhotonTag            = cms.InputTag('flashggPreselectedDiPhotons'),
                                     VertexTag              = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     GenParticleTag         = cms.InputTag('flashggPrunedGenParticles'),
                                     ElectronTag            = cms.InputTag('flashggSelectedElectrons'),
                                     #METTag                 = cms.InputTag('slimmedMETs'),
                                     MuonTag                = cms.InputTag('flashggSelectedMuons'),
                                     METTag                 = cms.InputTag('flashggMets'),
                                     #JetTag                 = cms.InputTag('flashggJets'),
                                     JetTags                = UnpackedJetCollectionVInputTag, # one jet per vertex
                                     #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                     #HTXSTags               = HTXSInputTags

                                     )
flashggHHWWggTagSequence = cms.Sequence( flashggHHWWggTag )

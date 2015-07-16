import FWCore.ParameterSet.Config as cms

selectedFlashggElectrons = cms.EDFilter("FLASHggElectronSelector",
    src = cms.InputTag("flashggElectrons"),
    cut = cms.string("pt > 9.")
)

selectedFlashggMuons = cms.EDFilter("FLASHggMuonSelector",
    src = cms.InputTag("flashggMuons"),
    cut = cms.string("pt > 9.")
)

selectedFlashggPhotons = cms.EDFilter("FLASHggPhotonSelector",
    src = cms.InputTag("flashggPhotons"),
    cut = cms.string("pt > 9.")
)

selectedFlashggJets = cms.EDFilter("FLASHggJetSelector",
    src = cms.InputTag("flashggJets"),
    cut = cms.string("pt > 15.")
)

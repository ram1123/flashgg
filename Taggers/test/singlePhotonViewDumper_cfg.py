#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

# process.source = cms.Source("PoolSource",
#                             fileNames=cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-2_7_7/2_7_7/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17-2_7_7-2_7_7-v0-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180117_155902/0000/myMicroAODOutputFile_1.root"
# ))
# process.source = cms.Source("PoolSource",
                        #     fileNames=cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2017_RR-31Mar2018_v2/legacyRun2FullV1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/Era2017_RR-31Mar2018_v2-legacyRun2FullV1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/190606_095939/0000/myMicroAODOutputFile_141.root") # GJet MicroAOD
                        #     secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/94X_mc2017_realistic_v10-v1/70000/FEFF944B-BAD7-E711-82EF-0025904A8ECE.root"
# ))

process.source = cms.Source ("PoolSource",
                        # fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190605_224652/0000/myMicroAODOutputFile_2.root"),
                        # secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv3/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/110000/CA97CD76-CD37-E911-ADBD-0090FAA57E64.root")
        # fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2017_RR-31Mar2018_v2/legacyRun2FullV1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/Era2017_RR-31Mar2018_v2-legacyRun2FullV1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/190606_095939/0000/myMicroAODOutputFile_141.root"),
        fileNames = cms.untracked.vstring("file:abeMicroAODTest.root"),  
        # fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/2EBBB73F-2ABA-E811-B295-0025905B85F6.root"),  
        # secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/FE600B1C-94B9-E811-A662-0CC47A4C8E46.root")
        #fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/DiPhotonJetsBox_M40_80-Sherpa/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/190715_222009/0000/myMicroAODOutputFile_62.root"),
        #secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/DiPhotonJetsBox_M40_80-Sherpa/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/60000/36573481-50CD-E811-95FB-0242AC130004.root")
        #fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/DoubleEG/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-Run2016B-17Jul2018_ver2-v1/190605_220256/0000/myMicroAODOutputFile_932.root"),
        #secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/20000/D03AED69-308D-E811-AFFC-008CFA197CD0.root","root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/20000/FCFDB07D-378D-E811-89A0-008CFAE45144.root")
            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
# process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)

process.load("flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi") 
process.flashggUpdatedIdMVADiPhotons.reRunRegression = cms.bool(False)
process.flashggUpdatedIdMVADiPhotons.doNon5x5transformation = cms.bool(False)
process.flashggUpdatedIdMVADiPhotons.do5x5correction = cms.bool(False)
process.flashggUpdatedIdMVADiPhotons.doIsoCorrection = cms.bool(False)

from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
process.kinPreselDiPhotons = flashggPreselectedDiPhotons.clone(
src = cms.InputTag("flashggUpdatedIdMVADiPhotons"),
cut=cms.string(
        "mass > 95"
        " && leadingPhoton.pt > 18 && subLeadingPhoton.pt > 18"
        " && abs(leadingPhoton.superCluster.eta)<2.5 && abs(subLeadingPhoton.superCluster.eta)<2.5 "
        " && ( abs(leadingPhoton.superCluster.eta)<1.4442 || abs(leadingPhoton.superCluster.eta)>1.566)"
        " && ( abs(subLeadingPhoton.superCluster.eta)<1.4442 || abs(subLeadingPhoton.superCluster.eta)>1.566)"
        " && (leadingPhoton.pt > 14 && leadingPhoton.hadTowOverEm < 0.15 && (leadingPhoton.full5x5_r9>0.8 || leadingPhoton.chargedHadronIso<20 || leadingPhoton.chargedHadronIso<(0.3*leadingPhoton.pt)))"
        " && (subLeadingPhoton.pt > 14 && subLeadingPhoton.hadTowOverEm < 0.15 && (subLeadingPhoton.full5x5_r9>0.8 || subLeadingPhoton.chargedHadronIso<20 || subLeadingPhoton.chargedHadronIso<(0.3*subLeadingPhoton.pt)))"
        )
)


process.flashggSinglePhotonViews = cms.EDProducer("FlashggSinglePhotonViewProducer",
                                                  DiPhotonTag=cms.InputTag('kinPreselDiPhotons'),                                         
                                                  maxCandidates = cms.int32(1),
                                                  EBreducedEcalRecHits                 = cms.InputTag("reducedEcalRecHitsEB"),
                                                  EEreducedEcalRecHits                 = cms.InputTag("reducedEcalRecHitsEE")
                                                  )

process.load("flashgg.Taggers.photonViewDumper_cfi") ##  import diphotonDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.photonViewDumper.src = "flashggSinglePhotonViews"
process.photonViewDumper.dumpTrees = True
process.photonViewDumper.dumpWorkspace = False
process.photonViewDumper.quietRooFit = True


## list of variables to be dumped in trees/datasets. Same variables for all categories
variables=[
           "pt                     := photon.pt",
        #    "energy                 := photon.energy",
           "eta                    := photon.eta", # What's the difference between photon eta and photon supercluster eta?
           "phi                    := photon.phi",
           "scEta                  := photon.superCluster.eta", 
           "scPhi                  := photon.superCluster.phi", 
        #    "SCRawE                 := photon.superCluster.rawEnergy",
        #    "etaWidth               := photon.superCluster.etaWidth",
        #    "phiWidth               := photon.superCluster.phiWidth",
        #    "covIphiIphi            := photon.sipip",
        #    "chgIsoWrtWorstVtx      := photon.pfChgIsoWrtWorstVtx03",
        #    "phoIso03               := photon.pfPhoIso03",
        #    "phoIsoCorr             := photon.pfPhoIso03Corr",
        #    "chgIsoWrtChosenVtx     := pfChIso03WrtChosenVtx",
        #    "hcalTowerSumEtConeDR03 := photon.hcalTowerSumEtConeDR03",
        #    "trkSumPtHollowConeDR03 := photon.trkSumPtHollowConeDR03",
        #    "hadTowOverEm           := photon.hadTowOverEm",
        #    "idMVA                  := phoIdMvaWrtChosenVtx",
        #    #"genIso                 := photon.userFloat('genIso')", 
        #    "eTrue                  := ? photon.hasMatchedGenPhoton ? photon.matchedGenPhoton.energy : 0",
        #    "sigmaIetaIeta          := photon.full5x5_sigmaIetaIeta",
        #    "r9                     := photon.full5x5_r9", 
        #    "esEffSigmaRR           := photon.esEffSigmaRR",
        #    "s4                     := photon.s4",
        #    "covIEtaIPhi            := photon.sieip",
        #    "esEnergy               := photon.superCluster.preshowerEnergy",
        #    "esEnergyOverRawE       := photon.superCluster.preshowerEnergy/photon.superCluster.rawEnergy",
        #    "ieta_0                 := ietas[0]",
        #    "iphi_0                 := iphis[0]",
        #    "recHit_0                := recHits[0]",
        #    "ietas                  := map(ietas::170,-85,85::ietas[0],ietas[1])"
        #    "recH"



           #"rho                    := global.rho",
           #"esEnergyPlane1         := photon.esEnergyPlane1",
           #"esEnergyPlane2         := photon.esEnergyPlane2",
           #"e1x3                   := photon.e1x3",
           #"e2x5max                := photon.e2x5max",
           #"e5x5                   := photon.e5x5"
           ]

# Doing this because don't know how to save a vector in variables

num_rec_hits = 100

# "jet0_pt                         := ? JetVector.size() >= 1 ? JetVector[0].pt() : -99 "
# "ietas_0                         := ? ietas.size()     >= 1 ? ietas[0]          : -9999 "
# variables.append("ietas_size := ietas.size()")
# variables.append("iphis_size := iphis.size()")
# variables.append("recHits_size := recHits.size()")
# variables.append("ieta_0 := ? ietas.size() >= 1 ? ietas[0] : -9999")
for i in range(num_rec_hits):
        variables.append("DOF1s_" + str(i) + " := DOF1s[" + str(i) + "] ")
        variables.append("DOF2s_" + str(i) + " := DOF2s[" + str(i) + "] ")
        variables.append("DOF3s_" + str(i) + " := DOF3s[" + str(i) + "] ")
        variables.append("recHit_" + str(i) + " := recHits[" + str(i) + "] ")
                
        # print"var = ","ietas_" + str(i) + " := ? ietas.size() >= " + str(i+1) + " ? ietas[" + str(i) + "] : -9999 "
        # variables.append("ietas_" + str(i) + " := ? ietas.size() >= " + str(i+1) + " ? ietas[" + str(i) + "] : -9999 ")
        # variables.append("iphis_" + str(i) + " := ? iphis.size() >= " + str(i+1) + " ? iphis[" + str(i) + "] : -9999 ")
        # variables.append("recHit_" + str(i) + " := ? recHits.size() >= " + str(i+1) + " ? recHits[" + str(i) + "] : -9999 ")

# variables.append("ietas_0 := ietas[0]")
print'variables = ',variables 

## list of histograms to be plotted
histograms=[
            "r9>>r9(110,0,1.1)",
            "scEta>>scEta(100,-2.5,2.5)",
        #     "ietas>>ietas(170,-85,85)",
        #     "iphis>>iphis(720,-360,360)",
        #     "recHits>>recHits(1000,0,100)",
            ]

## define categories and associated objects to dump
cfgTools.addCategory(process.photonViewDumper,
                     "Reject",
                     "abs(photon.superCluster.eta)>=1.4442&&abs(photon.superCluster.eta)<=1.566||abs(photon.superCluster.eta)>=2.5",
                     -1 ## if nSubcat is -1 do not store anythings
                     )

# interestng categories 
cfgTools.addCategories(process.photonViewDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [("promptPhotons","photon.genMatchType == 1",0), 
                        ("fakePhotons",  "photon.genMatchType != 1",0),
                        ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=variables,
                       ## histograms to be plotted. 
                       ## the variables need to be defined first
                       histograms=histograms,
                       )

process.p1 = cms.Path(process.flashggUpdatedIdMVADiPhotons*
                      process.kinPreselDiPhotons*
                      process.flashggSinglePhotonViews*
                      process.photonViewDumper
                      )


from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",10000)
customize(process)


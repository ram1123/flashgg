import importlib
import FWCore.ParameterSet.Config as cms

flashggDiPhotonSystematics = cms.EDProducer('FlashggDiPhotonSystematicProducer',
                src = cms.InputTag("flashggDifferentialPhoIdInputsCorrection"),
                SystMethods2D = cms.VPSet(),
                SystMethods = cms.VPSet()
)

def setupDiPhotonSystematics( process, options ):
   print'[setupDiPhotonSystematics] - Choosing diphoton systematics'
   process.load("flashgg.Systematics."+options.metaConditions["flashggDiPhotonSystematics"])
   print'[DEBUG1: setupDiPhotonSystematics] - Choosing diphoton systematics'
   sysmodule = importlib.import_module("flashgg.Systematics."+options.metaConditions["flashggDiPhotonSystematics"])
   print'[DEBUG2: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleHighR9EB)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleLowR9EB)
   print'[DEBUG3: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleHighR9EE)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleLowR9EE)
   print'[DEBUG4: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleGain6EB_EGM)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCScaleGain1EB_EGM)
   print'[DEBUG5: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MaterialCentralBarrel)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MaterialOuterBarrel)
   print'[DEBUG6: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MaterialForward)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.ShowerShapeHighR9EB)
   print'[DEBUG7: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.ShowerShapeHighR9EE)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.ShowerShapeLowR9EB)
   print'[DEBUG8: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.ShowerShapeLowR9EE)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.FNUFEB)
   print'[DEBUG9: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.FNUFEE)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCSmearHighR9EE)
   print'[DEBUG10: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCSmearLowR9EE)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCSmearHighR9EB)
   print'[DEBUG11: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MCSmearLowR9EB)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.MvaShift)
   print'[DEBUG12: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.PreselSF)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.electronVetoSF)
   print'[DEBUG13: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.TriggerWeight)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.LooseMvaSF)
   print'[DEBUG14: setupDiPhotonSystematics] - Choosing diphoton systematics'
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.SigmaEOverEShift)
   flashggDiPhotonSystematics.SystMethods.append(sysmodule.SigmaEOverESmearing)
   print'[DEBUG15: setupDiPhotonSystematics] - Choosing diphoton systematics'

   if (options.processId.count('ggh') or options.processId.count('GluGluHToGG')) and not options.processId.count('gghh'):
      flashggDiPhotonSystematics.SystMethods.append(sysmodule.FracRVWeight)
   #flashggDiPhotonSystematics.SystMethods.append(sysmodule.FracRVNvtxWeight)
   print'[DEBUG16: setupDiPhotonSystematics] - Choosing diphoton systematics'

   if options.ignoreNegR9:
      for syst_method in flashggDiPhotonSystematics.SystMethods:
         if hasattr(syst_method, "PhotonMethodName"):
            syst_method.OverallRange = syst_method.OverallRange._value + " && full5x5_r9>0."
   print'[DEBUG17: setupDiPhotonSystematics] - Choosing diphoton systematics'

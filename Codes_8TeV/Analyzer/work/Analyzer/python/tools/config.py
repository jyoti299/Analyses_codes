import FWCore.ParameterSet.Config as cms

process = cms.Process("NTuples")

process.load("work.Analyzer.analyzer_cfi")

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

process.GlobalTag.globaltag = cms.string('FT_53_V21_AN3::All')

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'],
                outputModules = [])


###########Required for using rechitsflagstobeExcluded defined in these for calculation of sigmaIetaIphi, sigmaIphiIphi......
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *                        ############ In barrel
process.My_Analyzer.RecHitFlagToBeExcludedEB = cleanedHybridSuperClusters.RecHitFlagToBeExcluded ######## Hybrid algorithm is used for barrel
process.My_Analyzer.RecHitSeverityToBeExcludedEB = cleanedHybridSuperClusters.RecHitSeverityToBeExcluded
from RecoEcal.EgammaClusterProducers.multi5x5BasicClusters_cfi import *                      ############ In endcap 
process.My_Analyzer.RecHitFlagToBeExcludedEE = multi5x5BasicClustersCleaned.RecHitFlagToBeExcluded ###### Multi5x5 alogrithm is used for endcap
process.My_Analyzer.RecHitSeverityToBeExcludedEE = cleanedHybridSuperClusters.RecHitSeverityToBeExcluded  ##### using severityflagstobeexcluded
                                                                                                          ##### from hybrid as this is not 
                                                                                                          ##### defined in multi5x5





process.My_Analyzer.debug = 0

#Uncomment it while running over MC
#process.My_Analyzer.isItAOD = True

if process.My_Analyzer.isItAOD == True:
     process.My_Analyzer.OutFile == "NTuples_MC.root"
     process.My_Analyzer.runMCPileUp = True
else:
     process.My_Analyzer.OutFile = "NTuples_Data.root"
     process.My_Analyzer.runPhotons = True    
     process.My_Analyzer.saveCrystalInfo = True
     process.My_Analyzer.savePFIsolation = True
     process.My_Analyzer.saveRecHitsInfo = True

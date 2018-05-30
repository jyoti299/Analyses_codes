import FWCore.ParameterSet.Config as cms

process = cms.Process("NTuples")

process.load("work.Analyzer.analyzer_cfi")

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

#Importing the Reco Jet Correction Services required by JetCorrectionRecord
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.GlobalTag.globaltag = cms.string('FT_53_V21_AN3::All') #for 2012 ABCD 22Jan2013 ReReco
 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://cmssrv32.fnal.gov//store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/82E6D4C4-0069-E211-86C1-002618943983.root',
        'root://xrootd.unl.edu//store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/82E6D4C4-0069-E211-86C1-002618943983.root',
       # '/store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/8079344F-FE68-E211-A139-003048FF9AA6.root',
#       'root://xrootd.unl.edu//store/mc/Summer12/G_Pt-120to170_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v1/0000/766E8556-B691-E111-9857-002481E0DC64.root',
       # '/store/mc/Summer12/G_Pt-120to170_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v1/0000/BA228BB5-7791-E111-9FC1-003048D43942.root',
    )
)

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'],
                 outputModules = [])


#Changing the level of corrections and applying pt and eta cut in the default pat::ak5CaloJet collection
process.patJetCorrFactors.levels = cms.vstring(['L1Offset','L2Relative','L3Absolute','L2L3Residual'])
process.selectedPatJets.cut = cms.string("pt > 10 & abs(eta) < 3.0")

#Adding AK5 PFJet collection in the pat::Jet collections (as default patSequence produce only the AK5 CaloJets)
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                 jetCollection    = cms.InputTag('ak5PFJets'),
                 algoLabel        = 'AK5',
                 typeLabel        = 'PF',
                 doJTA            = True,    
                 doBTagging       = True,
                 jetCorrLabel     = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])),
                 doType1MET       = True,
                 doL1Cleaning     = True,
                 doL1Counters     = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID          = True,
                 jetIdLabel       = "ak5"
                 )
process.selectedPatJetsAK5PF.cut = cms.string("pt > 10")


###########Required for using rechitsflagstobeExcluded defined in these for calculation of sigmaIetaIphi, sigmaIphiIphi......
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *                        ############ In barrel
process.My_Analyzer.RecHitFlagToBeExcludedEB = cleanedHybridSuperClusters.RecHitFlagToBeExcluded ######## Hybrid algorithm is used for barrel
process.My_Analyzer.RecHitSeverityToBeExcludedEB = cleanedHybridSuperClusters.RecHitSeverityToBeExcluded
from RecoEcal.EgammaClusterProducers.multi5x5BasicClusters_cfi import *                      ############ In endcap 
process.My_Analyzer.RecHitFlagToBeExcludedEE = multi5x5BasicClustersCleaned.RecHitFlagToBeExcluded ###### Multi5x5 alogrithm is used for endcap
process.My_Analyzer.RecHitSeverityToBeExcludedEE = cleanedHybridSuperClusters.RecHitSeverityToBeExcluded  ##### using severityflagstobeexclude
                                                                                                          ##### from hybrid as this is not 
                                                                                                          ##### defined in multi5x5

################################################################
### Adding and Rerunning B-tagging algorithms for reco::Jets ###
################################################################

###For AK5CaloJets (simply adding into the path as by default ak5CaloJets have been taken in the code)

from RecoJets.JetAssociationProducers.ak5JTA_cff import *

from RecoBTag.Configuration.RecoBTag_cff import *

process.ak5CaloJetsBTaggingSeq = cms.Sequence(ak5JetTracksAssociatorAtVertex *
                                              btagging *
                                              softElectronBJetTags  ###Explicitly running this algo as it is commented out in btagging seqence
                                              )

###For AK5PFJets

#----------------------------------------------------------------------------------------------------------------------------------------------------
#If b tagging for only AK5PF is required then do not need to do all this, simply import the above cff files and change the name of jets in          |
#ak5JetTracksAssociatorAtVertex by simply doing process.ak5JetTracksAssociatorAtVertex.jets = "ak5PFJets" and run the same seq as is done for calo. |
#This all needs to be done when b tagging for more than one kind of jet is required                                                                 |
#----------------------------------------------------------------------------------------------------------------------------------------------------

#Defining the PFJet in the JetTrack Associator used in impactParameterTagInfos required for b tagging (defined in ak5JTA_cff)
process.ak5PFJetTracksAssociatorAtVertex = process.ak5JetTracksAssociatorAtVertex.clone(
     jets = cms.InputTag("ak5PFJets")
     )

#Adding the newly defined JetTrackAssociator in ImpactParameterTagInfos
process.ak5PFImpactParameterTagInfos = process.impactParameterTagInfos.clone(
     jetTracks = cms.InputTag("ak5PFJetTracksAssociatorAtVertex")
     )

#Adding the new ImpactParameterTagInfos to SecondaryVertexTagInfos
process.ak5PFSecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone(
     trackIPTagInfos = cms.InputTag("ak5PFImpactParameterTagInfos")
     )

#Adding PFJets to SoftElectronTagInfos
process.ak5PFSoftElectronTagInfos = process.softElectronTagInfos.clone(
     jets = cms.InputTag("ak5PFJets")
     )

#Adding PFJets to SoftMuonTagInfos
process.ak5PFSoftMuonTagInfos = process.softMuonTagInfos.clone(
     jets = cms.InputTag("ak5PFJets")
     )

####Adding the newly defined algorithms (above) in all the b tagging algorithms wherever required

#TrackCountingHighEffBJetTags
process.ak5PFTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"))
     )

#TrackCountingHighPurBJetTags
process.ak5PFTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"))
     )

#SimpleSecondaryVertexHighEffBJetTags
process.ak5PFSimpleSecondaryVertexHighEffBJetTags = process.simpleSecondaryVertexHighEffBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFSecondaryVertexTagInfos"))
     )

#SimpleSecondaryVertexHighPurBJetTags
process.ak5PFSimpleSecondaryVertexHighPurBJetTags = process.simpleSecondaryVertexHighPurBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFSecondaryVertexTagInfos"))
     )

#CombinedSecondaryVertexBJetTags
process.ak5PFCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"),
                              cms.InputTag("ak5PFSecondaryVertexTagInfos"))
     )

#CombinedSecondaryVertexMVABJetTags
process.ak5PFCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"),
                              cms.InputTag("ak5PFSecondaryVertexTagInfos"))
     )
     
#JetProbabilityBJetTags
process.ak5PFJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"))
     )

#JetBProbabilityBJetTags
process.ak5PFJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFImpactParameterTagInfos"))
     )

#SoftElectronBJetTags
process.ak5PFSoftElectronBJetTags = process.softElectronBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFSoftElectronTagInfos"))
     )

#SoftMuonBJetTags
process.ak5PFSoftMuonBJetTags = process.softMuonBJetTags.clone(
     tagInfos = cms.VInputTag(cms.InputTag("ak5PFSoftMuonTagInfos"))
     )

#Adding all algorithms to the sequence
process.ak5PFJetsBTaggingSeq = cms.Sequence(process.ak5PFJetTracksAssociatorAtVertex *
                                            process.ak5PFImpactParameterTagInfos *
                                            process.ak5PFSecondaryVertexTagInfos *
                                            process.ak5PFSoftElectronTagInfos *
                                            process.ak5PFSoftMuonTagInfos *
                                            process.ak5PFTrackCountingHighEffBJetTags *
                                            process.ak5PFTrackCountingHighPurBJetTags *
                                            process.ak5PFSimpleSecondaryVertexHighEffBJetTags *
                                            process.ak5PFSimpleSecondaryVertexHighPurBJetTags *
                                            process.ak5PFCombinedSecondaryVertexBJetTags *
                                            process.ak5PFCombinedSecondaryVertexMVABJetTags *
                                            process.ak5PFJetProbabilityBJetTags *
                                            process.ak5PFJetBProbabilityBJetTags *
                                            process.ak5PFSoftElectronBJetTags *
                                            process.ak5PFSoftMuonBJetTags
                                            )


##Load the PU Jet Id Sequence for PF Pat Jets
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets=cms.InputTag("selectedPatJetsAK5PF")
process.puJetMva.jets=cms.InputTag("selectedPatJetsAK5PF")


process.My_Analyzer.debug = 0

#Uncomment it while running over MC
#process.My_Analyzer.isItAOD = True

if process.My_Analyzer.isItAOD == True:
     process.My_Analyzer.OutFile == "NTuples_MC.root"
     process.My_Analyzer.saveLuminosityInfo = True
     process.My_Analyzer.runMCPileUp = True
else:

     process.My_Analyzer.OutFile = "NTuples_Data.root"
     process.My_Analyzer.runPhotons = True    
     process.My_Analyzer.savePhotonCrystalInfo = True
     process.My_Analyzer.savePhotonPFIsolation = True
     process.My_Analyzer.savePhotonRecHitsInfo = True
     process.My_Analyzer.runJets = True
     process.My_Analyzer.runCaloPatJets = True
     process.My_Analyzer.runPFPatJets = True
     process.My_Analyzer.runCaloRecoJets = True
     process.My_Analyzer.runPFRecoJets = True
     process.My_Analyzer.savePUJetIdInfo_PFPat = True
     process.My_Analyzer.saveBTaggingInfo_CaloPat = True
     process.My_Analyzer.saveBTaggingInfo_PFPat = True
     process.My_Analyzer.saveBTaggingInfo_CaloReco = True
     process.My_Analyzer.saveBTaggingInfo_PFReco = True
     process.My_Analyzer.runVertex = True
     process.My_Analyzer.runTracks = True
     process.My_Analyzer.runScraping = True
     process.My_Analyzer.runHLT = True
process.My_Analyzer.CaloRecoJetCorrectionService = "ak5CaloL1L2L3Residual"
process.My_Analyzer.PFRecoJetCorrectionService = "ak5PFL1FastL2L3Residual"
      
process.p = cms.Path(process.patDefaultSequence *
                     process.puJetIdSqeuence *
                     process.ak5CaloJetsBTaggingSeq *
                     process.ak5PFJetsBTaggingSeq *
                     process.My_Analyzer)


# to write out events
process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string('outfile.root')
                               )
# write out file
process.output = cms.EndPath(
     process.out
     )

#print process.dumpPython()

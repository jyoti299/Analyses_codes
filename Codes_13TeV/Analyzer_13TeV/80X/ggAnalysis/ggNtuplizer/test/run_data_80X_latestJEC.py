import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v4')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
 connect = cms.string('sqlite:Spring16_23Sep2016AllV2_DATA.db'),
 toGet = cms.VPSet(
 cms.PSet(
  record = cms.string('JetCorrectionsRecord'),
  tag = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4PFchs'),
  label = cms.untracked.string('AK4PFchs')
 ),
 cms.PSet(
  record = cms.string('JetCorrectionsRecord'),
  tag = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK8PFchs'),
  label = cms.untracked.string('AK8PFchs')
  )))
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/data/Run2016B/SinglePhoton/MINIAOD/23Sep2016-v3/110000/06B8D6E3-4298-E611-BC45-0CC47A78A3EC.root'
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/04DC0281-C89F-E611-81C6-02163E0141E6.root'
        #'/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v2/80000/5A4402F5-638C-E611-A471-0025905A60AA.root'
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Winter_2016_reReco_v1_ele'),
    electrons = cms.InputTag("slimmedElectrons"),
    gbrForestName = cms.string('gedelectron_p4combination_25ns'),
    isMC = cms.bool(False),
    isSynchronization = cms.bool(False)
)

process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Winter_2016_reReco_v1_ele'),
    photons = cms.InputTag("slimmedPhotons"),
    isMC = cms.bool(False),
    isSynchronization = cms.bool(False)
)

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

# this loads all available b-taggers
#process.load("RecoBTag.Configuration.RecoBTag_cff")
#process.load("RecoBTag.SecondaryVertex.pfBoostedDoubleSecondaryVertexAK8BJetTags_cfi")
#process.pfImpactParameterTagInfosAK8.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.pfImpactParameterTagInfosAK8.candidates = cms.InputTag("packedPFCandidates")
#process.pfImpactParameterTagInfosAK8.jets = cms.InputTag("slimmedJetsAK8")
#process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfosAK8_cfi")
#process.pfInclusiveSecondaryVertexFinderTagInfosAK8.extSVCollection = cms.InputTag("slimmedSecondaryVertices")

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

jecLevels = [
  'Spring16_25nsV10p2_DATA_L2Relative_AK8PFchs.txt',
  'Spring16_25nsV10p2_DATA_L3Absolute_AK8PFchs.txt',
  'Spring16_25nsV10p2_DATA_L2L3Residual_AK8PFchs.txt'
]

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')
    )
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsAK8'),
    labelName = 'UpdatedJECAK8',
    jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')
    )

# MET correction and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=True                           
                           )

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.load("ggAnalysis.ggNtuplizer.ggPhotonIso_CITK_PUPPI_cff")
process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.jecAK8PayloadNames=cms.vstring(jecLevels)
process.ggNtuplizer.runHFElectrons=cms.bool(True)
process.ggNtuplizer.isAOD=cms.bool(False)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.doSyst=cms.bool(False)
process.ggNtuplizer.EleScale_Smearing_file="EgammaAnalysis/ElectronTools/data/ScalesSmearings/Winter_2016_reReco_v1_ele"

#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
                                                                
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
    
    process.p = cms.Path(
        ###process.reapplyJEC*
        ###process.pfImpactParameterTagInfosAK8 *
        ###process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
        ###process.pfBoostedDoubleSecondaryVertexAK8BJetTags *
        process.ggMETFiltersSequence* 
        process.calibratedPatElectrons*
        process.calibratedPatPhotons*
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggNtuplizer
        )
    
#print process.dumpPython()

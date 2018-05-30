import FWCore.ParameterSet.Config as cms

My_Analyzer = cms.EDAnalyzer("Analyzer",
                             debug                               = cms.untracked.int32( 0 ),                          
                             OutFile                             = cms.untracked.string(""),
                           
                             #All bools initialized to False
                             isItAOD                             = cms.untracked.bool(False),
                             saveLuminosityInfo                  = cms.untracked.bool(False),
                             runMCPileUp                         = cms.untracked.bool(False),
                             runPhotons                          = cms.untracked.bool(False),
                             savePhotonCrystalInfo               = cms.untracked.bool(False),
                             savePhotonPFIsolation               = cms.untracked.bool(False),
                             savePhotonRecHitsInfo               = cms.untracked.bool(False),
                             runJets                             = cms.untracked.bool(False),
                             runCaloPatJets                      = cms.untracked.bool(False),
                             runPFPatJets                        = cms.untracked.bool(False),
                             runCaloRecoJets                     = cms.untracked.bool(False),
                             runPFRecoJets                       = cms.untracked.bool(False),
                             savePUJetIdInfo_PFPat               = cms.untracked.bool(False),
                             saveBTaggingInfo_CaloPat            = cms.untracked.bool(False),
                             saveBTaggingInfo_PFPat              = cms.untracked.bool(False),
                             saveBTaggingInfo_CaloReco           = cms.untracked.bool(False),
                             saveBTaggingInfo_PFReco             = cms.untracked.bool(False),
                             runVertex                           = cms.untracked.bool(False),
                             runTracks                           = cms.untracked.bool(False),
                             runScraping                         = cms.untracked.bool(False),
                             runHLT                              = cms.untracked.bool(False),

 
                             #InputTags
                             MCpileupTag                         = cms.untracked.InputTag("addPileupInfo"),
                             patPhoTag                           = cms.untracked.InputTag("selectedPatPhotons"),
                             recoPhoTag                          = cms.untracked.InputTag("photons"),
                             pfCandidateTag                      = cms.untracked.InputTag("particleFlow"),
                             vertexTag                           = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BSTag                               = cms.untracked.InputTag("offlineBeamSpot"), 
                             ConvPhoTag                          = cms.untracked.InputTag("conversions"),
                             GsfEleTag                           = cms.untracked.InputTag("gsfElectrons"),
                             BarrelrechitTag                     = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                             EndcaprechitTag                     = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                             CaloPatJetTag                       = cms.untracked.InputTag("selectedPatJets"),
                             PFPatJetTag                         = cms.untracked.InputTag("selectedPatJetsAK5PF"),
                             CaloRecoJetTag                      = cms.untracked.InputTag("ak5CaloJets"),
                             PFRecoJetTag                        = cms.untracked.InputTag("ak5PFJets"),
                             #b tagging Input tags
                             CaloReco_TrkCountingHighEffBJetTag  = cms.untracked.InputTag("trackCountingHighEffBJetTags","","NTuples"),
                             CaloReco_TrkCountingHighPurBJetTag  = cms.untracked.InputTag("trackCountingHighPurBJetTags","","NTuples"),
                             CaloReco_SimpleSecVtxHighEffBJetTag = cms.untracked.InputTag("simpleSecondaryVertexHighEffBJetTags","","NTuples"),
                             CaloReco_SimpleSecVtxHighPurBJetTag = cms.untracked.InputTag("simpleSecondaryVertexHighPurBJetTags","","NTuples"),
                             CaloReco_CombinedSecVtxBJetTag      = cms.untracked.InputTag("combinedSecondaryVertexBJetTags","","NTuples"),
                             CaloReco_CombinedSecVtxMVABJetTag   = cms.untracked.InputTag("combinedSecondaryVertexMVABJetTags","","NTuples"),
                             CaloReco_JetProbBJetTag             = cms.untracked.InputTag("jetProbabilityBJetTags","","NTuples"),
                             CaloReco_JetBProbBJetTag            = cms.untracked.InputTag("jetBProbabilityBJetTags","","NTuples"),
                             CaloReco_SoftElectronBJetTag        = cms.untracked.InputTag("softElectronBJetTags","","NTuples"), 
                             CaloReco_SoftMuonBJetTag            = cms.untracked.InputTag("softMuonBJetTags","","NTuples"),
                             PFReco_TrkCountingHighEffBJetTag    = cms.untracked.InputTag("ak5PFTrackCountingHighEffBJetTags","","NTuples"),
                             PFReco_TrkCountingHighPurBJetTag    = cms.untracked.InputTag("ak5PFTrackCountingHighPurBJetTags","","NTuples"),
                             PFReco_SimpleSecVtxHighEffBJetTag   = cms.untracked.InputTag("ak5PFSimpleSecondaryVertexHighEffBJetTags","","NTuples"),
                             PFReco_SimpleSecVtxHighPurBJetTag   = cms.untracked.InputTag("ak5PFSimpleSecondaryVertexHighPurBJetTags","","NTuples"),
                             PFReco_CombinedSecVtxBJetTag        = cms.untracked.InputTag("ak5PFCombinedSecondaryVertexBJetTags","","NTuples"),
                             PFReco_CombinedSecVtxMVABJetTag     = cms.untracked.InputTag("ak5PFCombinedSecondaryVertexMVABJetTags","","NTuples"),
                             PFReco_JetProbBJetTag               = cms.untracked.InputTag("ak5PFJetProbabilityBJetTags","","NTuples"),
                             PFReco_JetBProbBJetTag              = cms.untracked.InputTag("ak5PFJetBProbabilityBJetTags","","NTuples"),
                             PFReco_SoftElectronBJetTag          = cms.untracked.InputTag("ak5PFSoftElectronBJetTags","","NTuples"),
                             PFReco_SoftMuonBJetTag              = cms.untracked.InputTag("ak5PFSoftMuonBJetTags","","NTuples"),
                             #Tracks Tag
                             TracksTag                           = cms.untracked.InputTag("generalTracks"),
                             #TriggerEvent Tag
                             triggerEventTag                     = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
                             #TriggerResults Tag
                             triggerResultsTag                   = cms.untracked.InputTag("TriggerResults","","HLT"),



                             #severity flags for Shower Shape
                             RecHitFlagToBeExcludedEB            = cms.untracked.vstring(""),
                             RecHitFlagToBeExcludedEE            = cms.untracked.vstring(""),
                             RecHitSeverityToBeExcludedEB        = cms.untracked.vstring(""),
                             RecHitSeverityToBeExcludedEE        = cms.untracked.vstring(""),
 
                             #Reco Jet Correction Service
                             CaloRecoJetCorrectionService        = cms.untracked.string(""),
                             PFRecoJetCorrectionService          = cms.untracked.string(""),


)


import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.coreTools import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("trig_analysis.Trig_analyzer.trig_analyzer_cfi")
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' ## same as in setup_cff.py

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
          'root://cmssrv32.fnal.gov//store/data/Run2012D/DoubleElectron/RAW-RECO/ZElectron-22Jan2013-v1/10000/4E810F67-698F-E211-B9CA-0026189438D3.root',
#         'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/0012FD53-AD92-E211-83BB-0025905964BC.root',
#         'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/0076C476-BE92-E211-9437-002354EF3BE3.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/00D660B1-C392-E211-981B-002618943953.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/028E9AC7-0E93-E211-B783-002618943923.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/02D65770-E392-E211-9DB3-003048678E2A.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/02D726B7-0393-E211-8E65-0025905964B6.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/02FFD6F1-E192-E211-B384-0025905964C4.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/040075CE-0E93-E211-8BD8-003048FF9AA6.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/041DB199-F992-E211-9995-0025905938B4.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/044F4260-CA92-E211-935B-003048679236.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/04655A49-F492-E211-AE4B-002618943974.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/04F25E7B-FE92-E211-9C0F-003048679076.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/0479472F-E592-E211-83D3-00304867904E.root',
#          'root://cmssrv32.fnal.gov//store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/002F22E7-AA92-E211-82E2-003048FF9AA6.root'
        )
)

#---------Fast Rho calculation-------------------------
#Rho for eta= 2.5
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets25.Ghost_EtaMax = cms.double(2.5)
process.kt6PFJets25.doAreaFastjet = True
process.kt6PFJets25.voronoiRfact = 0.9


#For calculating effective area give EAtarget any of the values (EleEANoCorr, EleEAData2011, EleEASummer11MC,EleEAFall11MC, EleEAData2012)
#corresponding to different data and MC. This is defined in this way in ElectronEffectiveArea.h
process.Trig_analyzer.EAtarget = "EleEAData2012"


process.Trig_analyzer.tnp_trigger = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50"
#process.Trig_analyzer.tnp_leadinglegfilter = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"
#process.Trig_analyzer.tnp_trailinglegfilter = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter" 
process.Trig_analyzer.Ele27_WP80 = "HLT_Ele27_WP80"


process.Trig_analyzer.pTOfflineProbeCut = 0.0

process.p = cms.Path(process.kt6PFJets25 * process.Trig_analyzer)

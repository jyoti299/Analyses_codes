# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Projects/GamJet/python/sherpa_13TeV_gam_2j_Pt-20-100_80X_MASTER_cff.py -s GEN -n 100 --no_exec --conditions auto:mc --eventcontent RAWSIM
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Projects/GamJet/python/sherpa_13TeV_gam_2j_Pt-20-100_80X_MASTER_cff.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('sherpa_13TeV_gam_2j_Pt-20-100_80X_MASTER_cff_py_GEN.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("SherpaGeneratorFilter",
    FetchSherpack = cms.bool(False),
    SherpaDefaultWeight = cms.double(1.0),
    SherpaParameters = cms.PSet(
        MPI_Cross_Sections = cms.vstring(' MPIs in Sherpa, Model = Amisic:', 
            ' semihard xsec = 43.2111 mb,', 
            ' non-diffractive xsec = 17.0318 mb with nd factor = 0.3142.'),
        Run = cms.vstring(' (run){', 
            ' EVENTS = 100;', 
            ' EVENT_MODE = HepMC;', 
            ' # avoid comix re-init after runcard modification', 
            ' WRITE_MAPPING_FILE 3;', 
            '}(run)', 
            ' (beam){', 
            ' BEAM_1 = 2212; BEAM_ENERGY_1 = 6500.;', 
            ' BEAM_2 = 2212; BEAM_ENERGY_2 = 6500.;', 
            '}(beam)', 
            ' (processes){', 
            ' Process 93 93 -> 22 93 93{1};', 
            ' Order (*,1);', 
            ' Enhance_Factor 2 {3};', 
            ' CKKW sqr(20./E_CMS);', 
            ' Integration_Error 0.005 {3}', 
            ' End process;', 
            '}(processes)', 
            ' (selector){', 
            ' PT 22 20. 100.;', 
            '}(selector)', 
            ' (shower){', 
            ' CSS_EW_MODE = 1', 
            '}(shower)', 
            ' (integration){', 
            ' FINISH_OPTIMIZATION = Off', 
            '}(integration)', 
            ' (isr){', 
            ' PDF_LIBRARY     = LHAPDFSherpa', 
            ' PDF_SET         = CT10', 
            ' PDF_SET_VERSION = 0', 
            ' PDF_GRID_PATH   = PDFsets', 
            '}(isr)', 
            ' (me){', 
            ' ME_SIGNAL_GENERATOR = Internal Comix', 
            ' EVENT_GENERATION_MODE = Unweighted;', 
            '}(me)', 
            ' (mi){', 
            ' MI_HANDLER = Amisic  # None or Amisic', 
            '}(mi)'),
        parameterSets = cms.vstring('MPI_Cross_Sections', 
            'Run')
    ),
    SherpaPath = cms.string('./'),
    SherpaPathPiece = cms.string('./'),
    SherpaProcess = cms.string('13TeV_gam_2j_Pt-20-100_80X'),
    SherpaResultDir = cms.string('Result'),
    SherpackChecksum = cms.string('a54966e19c6b611f9cf1f3bdd6152ee1'),
    SherpackLocation = cms.string('./'),
    crossSection = cms.untracked.double(-1),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.int32(0)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 



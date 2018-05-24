import FWCore.ParameterSet.Config as cms
import os

source = cms.Source("EmptySource")

generator = cms.EDFilter("SherpaGeneratorFilter",
  maxEventsToPrint = cms.int32(0),
  filterEfficiency = cms.untracked.double(1.0),
  crossSection = cms.untracked.double(-1),
  SherpaProcess = cms.string('8TeV_gamgam_2j_loop'),
  SherpackLocation = cms.string('./'),
  SherpackChecksum = cms.string('8c38277a1b5875b987f180129ace4a95'),
  FetchSherpack = cms.bool(False),
  SherpaPath = cms.string('./'),
  SherpaPathPiece = cms.string('./'),
  SherpaResultDir = cms.string('Result'),
  SherpaDefaultWeight = cms.double(1.0),
  SherpaParameters = cms.PSet(parameterSets = cms.vstring(
                             "MPI_Cross_Sections",
                             "Run"),
                              MPI_Cross_Sections = cms.vstring(
				" MPIs in Sherpa, Model = Amisic:",
				" semihard xsec = 38.3171 mb,",
				" non-diffractive xsec = 15.6401 mb with nd factor = 0.3142."
                                                  ),
                              Run = cms.vstring(
				" (run){",
				" EVENTS = 100;",
				" EVENT_MODE = HepMC;",
				" # avoid comix re-init after runcard modification",
				" WRITE_MAPPING_FILE 3;",
				"}(run)",
				" (beam){",
				" BEAM_1 = 2212; BEAM_ENERGY_1 = 4000.;",
				" BEAM_2 = 2212; BEAM_ENERGY_2 = 4000.;",
				"}(beam)",
				" (processes){",
				" Process 21 21 -> 22 22",
				" Scales VAR{Abs2(p[2]+p[3])}",
				" Loop_Generator gg_yy",
				" End process;",
				" Process 93 93 -> 22 22 93{2};",
				" Order (*,2);",
				" Enhance_Factor 2 {3};",
				" Enhance_Factor 5 {4};",
				" CKKW sqr(20./E_CMS);",
				" Integration_Error 0.005 {3};",
				" Integration_Error 0.03 {4};",
				" End process;",
				"}(processes)",
				" (selector){",
				" Mass  22 22 60. E_CMS;",
				" PT 22 10. E_CMS;",
				"}(selector)",
				" (shower){",
				" CSS_EW_MODE = 1",
				"}(shower)",
				" (integration){",
				" FINISH_OPTIMIZATION = Off",
				"}(integration)",
				" (isr){",
				" PDF_LIBRARY     = LHAPDFSherpa",
				" PDF_SET         = CT10",
				" PDF_SET_VERSION = 0",
				" PDF_GRID_PATH   = PDFsets",
				"}(isr)",
				" (me){",
				" ME_SIGNAL_GENERATOR = Internal Comix",
				" EVENT_GENERATION_MODE = Unweighted;",
				"}(me)",
				" (mi){",
				" MI_HANDLER = Amisic  # None or Amisic",
				"}(mi)"
                                                  ),
                             )
)

ProductionFilterSequence = cms.Sequence(generator)


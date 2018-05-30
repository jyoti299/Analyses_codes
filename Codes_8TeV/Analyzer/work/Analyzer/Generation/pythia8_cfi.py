import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    crossSection = cms.untracked.double(4.18), 
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring(

            'Main:timesAllowErrors = 10000',
            'ExcitedFermion:dg2dStar = on',  
            'ExcitedFermion:ug2uStar= on',
            '4000001:m0 = 1000.',
            '4000001:onMode = off',   
            '4000001:onIfMatch = 22 1',
            '4000002:m0 = 1000.',
            '4000002:onMode = off',
            '4000002:onIfMatch = 22 2',
            'ExcitedFermion:Lambda = 1000.',
            'ExcitedFermion:coupFprime = 1.',
            'ExcitedFermion:coupF = 1.',
            'ExcitedFermion:coupFcol = 1.',
            'Tune:pp 5'

            ),
        parameterSets = cms.vstring('processParameters')
    )
)

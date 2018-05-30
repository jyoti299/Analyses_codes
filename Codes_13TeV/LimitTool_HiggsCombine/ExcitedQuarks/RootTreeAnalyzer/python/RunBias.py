from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys
from RunCombine import massIterable

        
def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/dijet_bias.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="CaloDijet2016",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="qq",type="string",
                  help="signal model name")
    parser.add_option('--mass',dest="mass", default='750',type="string",
                  help="mass of resonance")
    parser.add_option('-l','--lumi',dest="lumi", default="12.910",type="string",
                  help="lumi in fb^-1, e.g.: 12.910")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    parser.add_option('--rMax',dest="rMax",default=5,type="float",
                  help="maximum r value (for better precision)")
    parser.add_option('-r',dest="r",default=1,type="float",
                  help="expect signal r value")
    parser.add_option('--rMin',dest="rMin",default=-5,type="float",
                  help="minimum r value (for better precision)")
    parser.add_option('--xsec',dest="xsec",default=10,type="float",
                  help="xsec for signal in pb (r = 1)")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='inputs/DijetFitResults.root',type="string",
                  help="input fit file")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store everything")
    parser.add_option('-t','--toys',dest="toys",default=1000,type="int",
                  help="number of toys")    
    parser.add_option('--gen-pdf',dest="genPdf", default="dijet", choices=['dijet','f1','f2','f3','f4','f5','f6'],
                  help="pdf for generating")
    parser.add_option('--fit-pdf',dest="fitPdf", default="dijet", choices=['dijet','f1','f2','f3','f4','f5','f6'],
                  help="pdf for fitting")
    parser.add_option('--asymptotic-file',dest="asymptoticFile",default=None,type="string",
                  help="load asymptotic cross section results file")
    parser.add_option('--computeBias',dest="computeBias",default=False,action='store_true',
                  help="Adding bias term in background component")
    
    (options,args) = parser.parse_args()

    pdfIndexMap = {'dijet': 0,
                   'f1': 1,
                   'f2': 2,
                   'f3': 3,
                   'f4': 4,
                   'f5': 5,
                   'f6': 6
                   }

    box = options.box
    lumi = float(options.lumi)
    model = options.model
    
    backgroundDsName = {'ExcitedQuarks2016':'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/Bias_study/TotalMC_Run2016_ReReco-BCDEFG_PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_Nodeta_Nodphi_CSVM_Mass695_QstarInvtMass.root',
                        'Excited1btagQuarks2016':'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/Bias_study/PseudoData_FromMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_1tagBstarInvtMass.root'
#'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/TotalMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_QstarInvtMass_No1p3scale.root'

#'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/TotalMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_QstarInvtMass.root'
#
                        }
            
    if box=='CaloDijet2015' or box=='CaloDijet20152016':
        signalDsName = 'inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15.root'%model
    elif box=='ExcitedQuarks2016':
        signalDsName = '/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/Bias_study/ResonanceShapes_%s_f-1p0_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16.root'%model
    elif box=='Excited1btagQuarks2016':
        signalDsName = '/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/Bias_study/ResonanceShapes_Bstar_1bTag_f1p0_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass700_80X_Summer16.root'
#'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/ResonanceShapes_%s_f-1p0_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass695_80X_Spring16.root'%model 
#'/uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/inputs/ResonanceShapes_%s_f-1p0_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16.root'%model 


#    xsecTree = None
    rDict = {}
    if options.asymptoticFile is not None:
        print "INFO: Input ref xsec file!"
        asymptoticRootFile = rt.TFile.Open(options.asymptoticFile,"READ")
        xsecTree = asymptoticRootFile.Get("xsecTree")        
        xsecTree.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = -1
        while True:
            entry = elist.Next()
            if entry == -1: break
            xsecTree.GetEntry(entry)
            rDict[int(eval('xsecTree.mass'))] = eval('xsecTree.xsecULExp_%s'%box)/options.xsec
    else:        
        for massPoint in massIterable(options.mass):   
            rDict[int(massPoint)] = options.r
    print rDict
        
    xsecString = '--xsec %f'%options.xsec
    rRangeString =  '--setPhysicsModelParameterRanges r=%.3f,%.3f'%(options.rMin,options.rMax)

    fixStringGen = '--setPhysicsModelParameters pdf_index=%i,r=%.5f'%(pdfIndexMap[options.genPdf],rDict[int(massPoint)])
    freezeStringGen = '--freezeNuisances pdf_index'
    if options.genPdf != 'dijet':
        freezeStringGen += ',p1_ExcitedQuarks2016,p2_ExcitedQuarks2016,p3_ExcitedQuarks2016'
#        freezeStringGen += ',p1_Excited1btagQuarks2016,p2_Excited1btagQuarks2016,p3_Excited1btagQuarks2016'
    if options.genPdf != 'f1':
        freezeStringGen += ',pF11_ExcitedQuarks2016,pF12_ExcitedQuarks2016,pF13_ExcitedQuarks2016'
#        freezeStringGen += ',pF11_Excited1btagQuarks2016,pF12_Excited1btagQuarks2016,pF13_Excited1btagQuarks2016'
    if options.genPdf != 'f2':
        freezeStringGen += ',pF21_ExcitedQuarks2016,pF22_ExcitedQuarks2016,pF23_ExcitedQuarks2016'
#        freezeStringGen += ',pF21_Excited1btagQuarks2016,pF22_Excited1btagQuarks2016,pF23_Excited1btagQuarks2016'
    if options.genPdf != 'f3':
        freezeStringGen += ',pF31_ExcitedQuarks2016,pF32_ExcitedQuarks2016'
#        freezeStringGen += ',pF31_Excited1btagQuarks2016,pF32_Excited1btagQuarks2016'
    if options.genPdf != 'f4':
        freezeStringGen += ',pF41_ExcitedQuarks2016,pF42_ExcitedQuarks2016'
#        freezeStringGen += ',pF41_Excited1btagQuarks2016,pF42_Excited1btagQuarks2016'
    if options.genPdf != 'f5':
        freezeStringGen += ',pF51_ExcitedQuarks2016,pF52_ExcitedQuarks2016,pF53_ExcitedQuarks2016'
#        freezeStringGen += ',pF51_Excited1btagQuarks2016,pF52_Excited1btagQuarks2016,pF53_Excited1btagQuarks2016'
    if options.genPdf != 'f6':
        freezeStringGen += ',pF61_ExcitedQuarks2016,pF62_ExcitedQuarks2016,pF63_ExcitedQuarks2016'
#        freezeStringGen += ',pF61_Excited1btagQuarks2016,pF62_Excited1btagQuarks2016,pF63_Excited1btagQuarks2016'
        
    fixStringFit = '--setPhysicsModelParameters pdf_index=%i,r=%.5f'%(pdfIndexMap[options.fitPdf],rDict[int(massPoint)])
    freezeStringFit = '--freezeNuisances pdf_index'
    if options.fitPdf != 'dijet':
        freezeStringFit += ',p1_ExcitedQuarks2016,p2_ExcitedQuarks2016,p3_ExcitedQuarks2016'
#        freezeStringFit += ',p1_Excited1btagQuarks2016,p2_Excited1btagQuarks2016,p3_Excited1btagQuarks2016'
    if options.fitPdf != 'f1':
        freezeStringFit += ',pF11_ExcitedQuarks2016,pF12_ExcitedQuarks2016,pF13_ExcitedQuarks2016'
#        freezeStringFit += ',pF11_Excited1btagQuarks2016,pF12_Excited1btagQuarks2016,pF13_Excited1btagQuarks2016'
    if options.fitPdf != 'f2':
        freezeStringFit += ',pF21_ExcitedQuarks2016,pF22_ExcitedQuarks2016,pF23_ExcitedQuarks2016'
#        freezeStringFit += ',pF21_Excited1btagQuarks2016,pF22_Excited1btagQuarks2016,pF23_Excited1btagQuarks2016'
    if options.fitPdf != 'f3':
        freezeStringFit += ',pF31_ExcitedQuarks2016,pF32_ExcitedQuarks2016'
#        freezeStringFit += ',pF31_Excited1btagQuarks2016,pF32_Excited1btagQuarks2016'
    if options.fitPdf != 'f4':
        freezeStringFit += ',pF41_ExcitedQuarks2016,pF42_ExcitedQuarks2016'
#        freezeStringFit += ',pF41_Excited1btagQuarks2016,pF42_Excited1btagQuarks2016'
    if options.fitPdf != 'f5':
        freezeStringFit += ',pF51_ExcitedQuarks2016,pF52_ExcitedQuarks2016,pF53_ExcitedQuarks2016'
#        freezeStringFit += ',pF51_Excited1btagQuarks2016,pF52_Excited1btagQuarks2016,pF53_Excited1btagQuarks2016'
    if options.fitPdf != 'f6':
        freezeStringFit += ',pF61_ExcitedQuarks2016,pF62_ExcitedQuarks2016,pF63_ExcitedQuarks2016'
#        freezeStringFit += ',pF61_Excited1btagQuarks2016,pF62_Excited1btagQuarks2016,pF63_Excited1btagQuarks2016'

    BiasString = ''
    if options.computeBias:
        BiasString = '--computeBias'
            
    for massPoint in massIterable(options.mass):        
        exec_me('python /uscms_data/d3/rocky86/13TeV/HiggsCombineTool/HiggsCombine_fromdijet/CMSSW_7_4_14/src/ExcitedQuarks/RootTreeAnalyzer/python/WriteDataCard.py -m %s --mass %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s --multi'%(model, massPoint, options.inputFitFile,1000*lumi,options.config,box,options.outDir,signalDsName,backgroundDsName[box],xsecString,BiasString),options.dryRun)
        exec_me('combine -M GenerateOnly %s/ExcitedQuarks_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_r-%.5f_%s_%s_%s %s %s %s --toysFrequentist --saveToys --expectSignal %.5f -t %i'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,rRangeString,fixStringGen,freezeStringGen,rDict[int(massPoint)],options.toys),options.dryRun)
        exec_me('combine -M MaxLikelihoodFit  %s/ExcitedQuarks_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_r-%.5f_%s_%s_%s --toysFile higgsCombine%s_%s_lumi-%.3f_r-%.5f_%s_%s_%s.GenerateOnly.mH120.123456.root -t %i %s %s %s --minimizerTolerance 0.01 --minimizerStrategy 2 --minos poi --saveWorkspace --robustFit 1 --minimizerAlgoForMinos Minuit2,Migrad'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,options.toys,rRangeString,fixStringFit,freezeStringFit),options.dryRun)
        exec_me('mv higgsCombine%s_%s_lumi-%.3f_r-%.5f_%s_%s_%s.GenerateOnly.mH120.123456.root %s/'%(model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,options.outDir),options.dryRun)
        exec_me('mv higgsCombine%s_%s_lumi-%.3f_r-%.5f_%s_%s_%s.MaxLikelihoodFit.mH120.123456.root %s/'%(model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,options.outDir),options.dryRun)
        exec_me('mv mlfit%s_%s_lumi-%.3f_r-%.5f_%s_%s_%s.root %s/'%(model,massPoint,lumi,rDict[int(massPoint)],box,options.genPdf,options.fitPdf,options.outDir),options.dryRun)

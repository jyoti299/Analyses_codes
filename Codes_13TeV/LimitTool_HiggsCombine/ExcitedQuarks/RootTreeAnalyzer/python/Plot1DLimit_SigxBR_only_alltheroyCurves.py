##### Script is done for sigma X BR X acc X eff, (the q* xs is getting divided by eff)
##### script not done for sigma X BR, the obs and exp limits are yet to be divided by eff, not done yet. 
#! /usr/bin/env python
import ROOT as rt
import os.path
import sys, glob, re
from array import *
from optparse import OptionParser

def getThyXsecDict():    
    thyXsecDict = {}
#    xsecFiles = ['data/Qstar_f1p0_xsec.txt', 'data/Qstar_f0p5_xsec.txt', 'data/Qstar_f0p1_xsec.txt', 'data/Bstar_f1p0_xsec.txt', 'data/Bstar_f0p5_xsec.txt', 'data/Bstar_f0p1_xsec.txt']
    xsecFiles = ['data/Qstar_f1p0_xsec.txt', 'data/Qstar_f0p5_xsec.txt', 'data/Qstar_f0p1_xsec.txt']
#    xsecFiles = ['data/Bstar_f0p1_xsec.txt']
#    xsecFiles = ['data/Bstar_1Tag_f1p0_xsec.txt']
#    xsecFiles = ['data/Optimization_80X/CSVDisc_Opti/Bstar_f1p0_xsec_CSVM.txt']
    print xsecFiles
    for xsecFile in xsecFiles:
        moreThyModels = []
        f = open(xsecFile)
        for i,line in enumerate(f.readlines()):
            if line[0]=='#': continue
            line = line.replace('\n','')
            line = line.replace('\t','')
            line = line.replace('\r','')
            lineList = [l for l in line.split(" ") if l!='']

            if lineList[0]=='Mass':
                for l in lineList:
                    if l=='Mass': continue
                    if l.find('Eff')!=-1: continue
                    thyXsecDict[l] = {}
                    moreThyModels.append(l)
            else:
                for j, thyModel in enumerate(moreThyModels):
                    thyXsecDict[thyModel][int(float(lineList[0]))] = float(lineList[j+1]) 
                    #print thyXsecDict[thyModel][int(float(lineList[0]))]
        f.close()

        #thyXsecDict['AxigluonkNLO'] = {}
    #for (mass,thyXsec) in thyXsecDict['Axigluon'].iteritems():
        #thyXsecDict['AxigluonkNLO'][mass] = 1.08 * thyXsec
    #for (mass,thyXsec) in thyXsecDict['DM1GeV'].iteritems():
        #thyXsecDict['DM1GeV'][mass] = (5./6.) * thyXsec
    return thyXsecDict


def getThyEffDict():    
    thyEffDict = {}
    effFiles = ['data/Qstar_f1p0_eff.txt']
#    effFiles = ['data/Qstar_f1p0_eff_EnvlCal_2sigma.txt']
#    effFiles = ['data/Bstar_f0p1_eff.txt']
    print effFiles
    for effFile in effFiles:
        f = open(effFile)
        for i,line in enumerate(f.readlines()):
            if line[0]=='#': continue
            line = line.replace('\n','')
            line = line.replace('\t','')
            line = line.replace('\r','')
            lineList = [l for l in line.split(" ") if l!='']

            if lineList[0]=='Mass': continue
            else:
                thyEffDict[int(float(lineList[0]))] = float(lineList[1])           
                #print thyEffDict[int(float(lineList[0]))]
        f.close()
        #print thyEffDict

    return thyEffDict

#def file_key(filename):
    #massPoint = re.findall("[0-9]+.000000",filename)
    #gluinoMass    = massPoint[0]
    #LSPMass  = massPoint[1]
    #return float(gluinoMass)
    
def getHybridCLsArrays(directory, model, Box, bayes):
    if bayes:
        tfile = rt.TFile.Open("%s/xsecUL_MarkovChainMC_%s_%s.root"%(directory,model,Box))
    else:
        tfile = rt.TFile.Open("%s/xsecUL_Asymptotic_%s_%s.root"%(directory,model,Box))
    xsecTree = tfile.Get("xsecTree")
    
    gluinoMassArray = array('d')
    gluinoMassArray_er = array('d')
    observedLimit = array('d')
    observedLimit_er = array('d')
    expectedLimit = array('d')
    expectedLimit_minus1sigma = array('d')
    expectedLimit_plus1sigma = array('d')
    expectedLimit_minus2sigma = array('d')
    expectedLimit_plus2sigma = array('d')

    
    xsecTree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        xsecTree.GetEntry(entry)

        gluinoMassArray.append(xsecTree.mass/1000)
        gluinoMassArray_er.append(0.0)
        
        exec 'xsecULObs = xsecTree.xsecULObs_%s'%Box
        exec 'xsecULExp = xsecTree.xsecULExp_%s'%Box
        exec 'xsecULExpPlus = xsecTree.xsecULExpPlus_%s'%Box
        exec 'xsecULExpMinus = xsecTree.xsecULExpMinus_%s'%Box
        exec 'xsecULExpPlus2 = xsecTree.xsecULExpPlus2_%s'%Box
        exec 'xsecULExpMinus2 = xsecTree.xsecULExpMinus2_%s'%Box
                       
        xsecULObs = xsecULObs
        xsecULExp = xsecULExp
        observedLimit.append(xsecULObs)#*crossSections[i])
        observedLimit_er.append(0.0)#*crossSections[i])

        expectedLimit.append(xsecULExp)#*crossSections[i])
                    
        xsecULExpPlus = max(xsecULExpPlus,xsecULExp)
        xsecULExpMinus = min(xsecULExpMinus,xsecULExp)
        xsecULExpPlus2 = max(xsecULExpPlus2,xsecULExpPlus)
        xsecULExpMinus2 = min(xsecULExpMinus2,xsecULExpMinus)

        expectedLimit_minus1sigma.append(xsecULExp - xsecULExpMinus)#*crossSections[i])
        expectedLimit_plus1sigma.append(xsecULExpPlus - xsecULExp)#*crossSections[i])
        expectedLimit_minus2sigma.append(xsecULExp - xsecULExpMinus2)#*crossSections[i])
        expectedLimit_plus2sigma.append(xsecULExpPlus2 - xsecULExp)#*crossSections[i])

    return gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma


def getSignificanceArrays(directory, model, Box):
    tfile = rt.TFile.Open("%s/xsecUL_ProfileLikelihood_%s_%s.root"%(directory,model,Box))
    xsecTree = tfile.Get("xsecTree")
    
    gluinoMassArray = array('d')
    gluinoMassArray_er = array('d')
    observedLimit = array('d')
    observedLimit_er = array('d')
    expectedLimit = array('d')
    expectedLimit_minus1sigma = array('d')
    expectedLimit_plus1sigma = array('d')
    expectedLimit_minus2sigma = array('d')
    expectedLimit_plus2sigma = array('d')

    
    xsecTree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        xsecTree.GetEntry(entry)

        gluinoMassArray.append(xsecTree.mass)
        gluinoMassArray_er.append(0.0)
        
        exec 'xsecULObs = xsecTree.xsecULObs_%s'%Box
        exec 'xsecULExp = xsecTree.xsecULExp_%s'%Box
        exec 'xsecULExpPlus = xsecTree.xsecULExpPlus_%s'%Box
        exec 'xsecULExpMinus = xsecTree.xsecULExpMinus_%s'%Box
        exec 'xsecULExpPlus2 = xsecTree.xsecULExpPlus2_%s'%Box
        exec 'xsecULExpMinus2 = xsecTree.xsecULExpMinus2_%s'%Box

            
            
        xsecULObs = xsecULObs
        xsecULExp = xsecULExp
        observedLimit.append(xsecULObs)#*crossSections[i])
        observedLimit_er.append(0.0)#*crossSections[i])

        expectedLimit.append(xsecULExp)#*crossSections[i])
        
        xsecULExpPlus = max(xsecULExpPlus,xsecULExp)
        xsecULExpMinus = min(xsecULExpMinus,xsecULExp)
        xsecULExpPlus2 = max(xsecULExpPlus2,xsecULExpPlus)
        xsecULExpMinus2 = min(xsecULExpMinus2,xsecULExpMinus)

        expectedLimit_minus1sigma.append(xsecULExp - xsecULExpMinus)#*crossSections[i])
        expectedLimit_plus1sigma.append(xsecULExpPlus - xsecULExp)#*crossSections[i])
        expectedLimit_minus2sigma.append(xsecULExp - xsecULExpMinus2)#*crossSections[i])
        expectedLimit_plus2sigma.append(xsecULExpPlus2 - xsecULExp)#*crossSections[i])
    

    return gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma
    
def setstyle():
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(600) #Width of canvas
    rt.gStyle.SetCanvasDefX(0)   #POsition on screen
    rt.gStyle.SetCanvasDefY(0)
    
    # For the Pad:
    rt.gStyle.SetPadBorderMode(0)
    # rt.gStyle.SetPadBorderSize(Width_t size = 1)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
    
    # For the frame:
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    # set the paper & margin sizes
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.09)
    rt.gStyle.SetPadRightMargin(0.065)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.06,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.06,"xyz")   # set the pad title size
    rt.gStyle.SetTitleSize(0.052,"xy")   # set the pad title size
    rt.gStyle.SetTitleOffset(1.2,"xy")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.03,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    rt.gStyle.SetMarkerStyle(8)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(11111111)
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    
    ncontours = 999
    
    stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
    green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
    blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
        
    npoints = len(s)
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)
   
    rt.gStyle.cd()
        

if __name__ == '__main__':
    
    rt.gROOT.SetBatch()
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="ExcitedQuarks",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")    
    parser.add_option('-l','--lumi',dest="lumi", default=1.,type="float",
                  help="integrated luminosity in fb^-1")
    parser.add_option('--massMin',dest="massMin", default=500.,type="float",
                  help="minimum mass")
    parser.add_option('--massMax',dest="massMax", default=8000.,type="float",
                  help="maximum mass")
    parser.add_option('--xsecMin',dest="xsecMin", default=1e-4,type="float",
                  help="minimum mass")
    parser.add_option('--xsecMax',dest="xsecMax", default=1e4,type="float",
                  help="maximum mass")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--bayes',dest="bayes",default=False,action='store_true',
                  help="for bayesian limits")
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="for no systematics limits")
    parser.add_option('-f','--coup',dest="SigCoup",default="f1p0",type="string",
                  help="SM coupling of signal")
    
    (options,args) = parser.parse_args()
    Boxes = options.box.split('_')
    models = options.model.split('_')
    model = models[0]
    directory      = options.outDir
    Box = Boxes[0]
    box = Box.lower()
    coupling = options.SigCoup
        
    thyXsecDict = getThyXsecDict() 
    thyEffDict  = getThyEffDict()
    thyModels = thyXsecDict.keys()

    thyModelsToDraw = thyXsecDict.keys()
#    thyModelsToDraw = []
#    if options.model=='Qstar':
#        if options.SigCoup=='f1p0':
#            thyModelsToDraw = ['Qstar_f1p0']        
#        elif options.SigCoup=='f0p5':
#            thyModelsToDraw = ['Qstar_f0p5']
#        elif options.SigCoup=='f0p1':
#            thyModelsToDraw = ['Qstar_f0p1']
#    if options.model=='Bstar':
#        if options.SigCoup=='f1p0':
#            thyModelsToDraw = ['Bstar_f1p0']        
#        elif options.SigCoup=='f0p5':
#            thyModelsToDraw = ['Bstar_f0p5']
#        elif options.SigCoup=='f0p1':
#            thyModelsToDraw = ['Bstar_f0p1']

    lineStyle = {'Qstar_f1p0':4,
                 'Qstar_f0p5':5,
                 'Qstar_f0p1':8,
                 'Bstar_f1p0':2,
                 'Bstar_f0p5':3,
                 'Bstar_f0p1':4,
                 'None':1               
                 }
        
    lineColor = {'Qstar_f1p0':rt.kRed+1,
                 'Qstar_f0p5':rt.kMagenta+2,
                 'Qstar_f0p1':rt.kCyan+2,
                 'Bstar_f1p0':rt.kRed+1,
                 'Bstar_f0p5':rt.kOrange+2,
                 'Bstar_f0p1':rt.kGreen+1,
                 'None':1
                 }
        
    markerStyle = {'Qstar_f1p0':20,
                   'Qstar_f0p5':29,
                   'Qstar_f0p1':21,
                   'Bstar_f1p0':24,
                   'Bstar_f0p5':30,
                   'Bstar_f0p1':25
                 }
        
    legendLabel = {'Qstar_f1p0':'Excited quark (#font[12]{f}=1.0)  ',
                   'Qstar_f0p5':'Excited quark (#font[12]{f}=0.5)  ',
                   'Qstar_f0p1':'Excited quark (#font[12]{f}=0.1)  ',
                   'Bstar_f1p0':'Excited b-quark (#font[12]{f}=1.0)',
                   'Bstar_f0p5':'Excited b-quark (#font[12]{f}=0.5)',
                   'Bstar_f0p1':'Excited b-quark (#font[12]{f}=0.1)'                 
                   }
    
    mass_xsec_gev = {}
    div = {}
    mass_xsec = {}
    sig_xsec = {}
    N_g_xsec = {}
    xsec_gr_nom = {}
    for thyModel in thyModelsToDraw:        
        mass_xsec_gev[thyModel] = array('d')
        div[thyModel] = array('d')
        mass_xsec[thyModel] = array('d')
        sig_xsec[thyModel] = array('d')
        for mg in sorted(thyXsecDict[thyModel].keys()):
            mass_xsec_gev[thyModel].append(mg)
            div[thyModel].append(1000)
            sig_xsec[thyModel].append(thyXsecDict[thyModel][mg])

        mass_xsec[thyModel] = array('d', [float(b) / float(m) for b,m in zip(mass_xsec_gev[thyModel], div[thyModel])])
        #print sig_xsec[thyModel] 
        N_g_xsec[thyModel] = len(mass_xsec[thyModel])
        xsec_gr_nom[thyModel] = rt.TGraph(N_g_xsec[thyModel], mass_xsec[thyModel], sig_xsec[thyModel])
        xsec_gr_nom[thyModel].SetMarkerSize(0)
        xsec_gr_nom[thyModel].SetLineWidth(2)
        xsec_gr_nom[thyModel].SetLineStyle(lineStyle[thyModel])
        xsec_gr_nom[thyModel].SetLineColor(lineColor[thyModel])

    N_g_eff = {}
    mass_eff = array('d')
    sig_eff = array('d')
    for mgg in sorted(thyEffDict.keys()):
        mass_eff.append(mgg)
        sig_eff.append(thyEffDict[mgg])
    N_g_eff = len(mass_eff)

    #print sig_eff
            
    setstyle()
    rt.gStyle.SetOptStat(0)
    c = rt.TCanvas("c","c",800,800)
    if options.doSignificance:
        c.SetLogy(0)
    else:        
        c.SetLogy()

    h_limit = rt.TMultiGraph()
    gr_observedLimit = {}
    gr_expectedLimit = {}
    gr_expectedLimit2sigma = {}
    gr_expectedLimit1sigma = {}
    gluinoMassArray = {}
    gluinoMassArray_er = {}
    observedLimit = {}
    observedLimit_er = {}
    expectedLimit = {}
    expectedLimit_minus1sigma = {}
    expectedLimit_plus1sigma = {}
    expectedLimit_minus2sigma = {}
    expectedLimit_plus2sigma = {}
    
    if options.doSignificance:
        h_limit.SetTitle(" ;Resonance mass [TeV];Local significance n#sigma")
    else:
        h_limit.SetTitle(" ;Resonance mass [TeV]; #sigma #times #bf{#it{#Beta}} [pb]")

    for Box in Boxes:
        for model in models:
            if len(models)>1:
                #directory =  options.outDir+'/%s_IntermediateRange'%model
                directory =  options.outDir+'/%s'%model
            if options.doSignificance:
                gluinoMassArray[(Box,model)], gluinoMassArray_er[(Box,model)], observedLimit[(Box,model)], observedLimit_er[(Box,model)], expectedLimit[(Box,model)], expectedLimit_minus1sigma[(Box,model)], expectedLimit_plus1sigma[(Box,model)], expectedLimit_minus2sigma[(Box,model)], expectedLimit_plus2sigma[(Box,model)] = getSignificanceArrays(directory, model, Box)
            else:        
                gluinoMassArray[(Box,model)], gluinoMassArray_er[(Box,model)], observedLimit[(Box,model)], observedLimit_er[(Box,model)], expectedLimit[(Box,model)], expectedLimit_minus1sigma[(Box,model)], expectedLimit_plus1sigma[(Box,model)], expectedLimit_minus2sigma[(Box,model)], expectedLimit_plus2sigma[(Box,model)] = getHybridCLsArrays(directory, model, Box, options.bayes)

            print "+++++++++++++++++++++++"
            print gluinoMassArray[(Box,model)]
            print "+++++++++++++++++++++++"
            print "+++++++++++++++++++++++"
            print observedLimit[(Box,model)]
            print "+++++++++++++++++++++++"

            nPoints = len(observedLimit[(Box,model)])

            if nPoints == N_g_eff:
                observedLimit[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(observedLimit[(Box,model)], sig_eff)])
                expectedLimit[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(expectedLimit[(Box,model)], sig_eff)])
                expectedLimit_minus1sigma[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(expectedLimit_minus1sigma[(Box,model)], sig_eff)])
                expectedLimit_plus1sigma[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(expectedLimit_plus1sigma[(Box,model)], sig_eff)])
                expectedLimit_minus2sigma[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(expectedLimit_minus2sigma[(Box,model)], sig_eff)])
                expectedLimit_plus2sigma[(Box,model)] = array('d', [float(b) / float(m) for b,m in zip(expectedLimit_plus2sigma[(Box,model)], sig_eff)])
            else:
                 print "Error:::: limit and efficiency arrays are of different lengths"

            #print "+++++++++++++++++++++++"
            #print expectedLimit[(Box,model)]
            #print "+++++++++++++++++++++++"


            gr_observedLimit[(Box,model)] = rt.TGraph(nPoints, gluinoMassArray[(Box,model)], observedLimit[(Box,model)])
            gr_observedLimit[(Box,model)].SetMarkerColor(1)
            gr_observedLimit[(Box,model)].SetMarkerStyle(22)
            gr_observedLimit[(Box,model)].SetMarkerSize(0.8)
            gr_observedLimit[(Box,model)].SetLineWidth(3)
            gr_observedLimit[(Box,model)].SetLineColor(rt.kBlack)
            gr_observedLimit[(Box,model)].SetMarkerStyle(20)
            if len(models)>1:
                gr_observedLimit[(Box,model)].SetLineColor(lineColor[model])
                gr_observedLimit[(Box,model)].SetMarkerStyle(markerStyle[model])
                gr_observedLimit[(Box,model)].SetMarkerColor(lineColor[model])


            gr_expectedLimit[(Box,model)] = rt.TGraph(nPoints, gluinoMassArray[(Box,model)], expectedLimit[(Box,model)])
            gr_expectedLimit[(Box,model)].SetLineWidth(3)
            gr_expectedLimit[(Box,model)].SetLineStyle(2)
            gr_expectedLimit[(Box,model)].SetLineColor(rt.kBlue)
            if len(models)>1:
                gr_expectedLimit[(Box,model)].SetLineColor(lineColor[model])

            gr_expectedLimit2sigma[(Box,model)] = rt.TGraphAsymmErrors(nPoints, gluinoMassArray[(Box,model)], expectedLimit[(Box,model)], gluinoMassArray_er[(Box,model)], gluinoMassArray_er[(Box,model)], expectedLimit_minus2sigma[(Box,model)], expectedLimit_plus2sigma[(Box,model)])
            #gr_expectedLimit2sigma[(Box,model)].SetLineColor(5)
            #gr_expectedLimit2sigma[(Box,model)].SetFillColor(5)
            gr_expectedLimit2sigma[(Box,model)].SetLineColor(rt.kOrange)
            gr_expectedLimit2sigma[(Box,model)].SetFillColor(rt.kOrange)
            #gr_expectedLimit2sigma[(Box,model)].SetLineColor(rt.kYellow)
            #gr_expectedLimit2sigma[(Box,model)].SetFillColor(rt.kYellow)
            gr_expectedLimit2sigma[(Box,model)].SetFillStyle(1001)

            gr_expectedLimit1sigma[(Box,model)] = rt.TGraphAsymmErrors(nPoints, gluinoMassArray[(Box,model)], expectedLimit[(Box,model)], gluinoMassArray_er[(Box,model)], gluinoMassArray_er[(Box,model)], expectedLimit_minus1sigma[(Box,model)], expectedLimit_plus1sigma[(Box,model)])

            #gr_expectedLimit1sigma[(Box,model)].SetLineColor(rt.kGreen-7)
            #gr_expectedLimit1sigma[(Box,model)].SetFillColor(rt.kGreen-7)
            gr_expectedLimit1sigma[(Box,model)].SetLineColor(rt.kGreen+1)
            gr_expectedLimit1sigma[(Box,model)].SetFillColor(rt.kGreen+1)

            if len(models)==1:
                h_limit.Add(gr_expectedLimit2sigma[(Box,model)])
                h_limit.Add(gr_expectedLimit1sigma[(Box,model)])
                h_limit.Add(gr_observedLimit[(Box,model)])

        
    for thyModel in thyModelsToDraw:
        h_limit.Add(xsec_gr_nom[thyModel])
        
    h_limit.Draw("a3")
    h_limit.GetXaxis().SetLimits(options.massMin,options.massMax)
    h_limit.GetXaxis().CenterTitle()
    h_limit.GetYaxis().CenterTitle()
    h_limit.GetXaxis().SetLabelFont(42)
    h_limit.GetYaxis().SetLabelFont(42)
    h_limit.GetXaxis().SetLabelSize(0.04)
    h_limit.GetYaxis().SetLabelSize(0.04)
    h_limit.GetXaxis().SetLabelOffset(0.008)
    h_limit.GetYaxis().SetLabelOffset(0.008)

    h_limit.GetXaxis().SetTitleFont(42)
    h_limit.GetYaxis().SetTitleFont(42)
    h_limit.GetXaxis().SetTitleSize(0.045)
    h_limit.GetYaxis().SetTitleSize(0.045)
    h_limit.GetXaxis().SetTitleOffset(1.0)
    h_limit.GetYaxis().SetTitleOffset(1.4)

    if options.doSignificance:
        h_limit.SetMaximum(4)
        h_limit.SetMinimum(0)
    else:
        h_limit.SetMaximum(options.xsecMax)
        h_limit.SetMinimum(options.xsecMin)
            
    h_limit.Draw("a3")
    if options.doSignificance:
        h_limit.GetYaxis().SetNdivisions(405,True)

    for Box in Boxes:
        for model in models:    
            if options.doSignificance:
                gr_observedLimit[(Box,model)].SetMarkerStyle(21)
                gr_observedLimit[(Box,model)].SetMarkerSize(1)
                gr_observedLimit[(Box,model)].SetLineColor(rt.kRed)
                gr_observedLimit[(Box,model)].SetMarkerColor(rt.kBlue)
                gr_observedLimit[(Box,model)].Draw("lp SAME")
            else:
                if len(models)==1:
                    gr_expectedLimit[(Box,model)].Draw("c same")
                for thyModel in thyModelsToDraw:
                    xsec_gr_nom[thyModel].Draw("c same")
            gr_observedLimit[(Box,model)].Draw("lp SAME")

            gr_expectedLimit1sigma[(Box,model)].SetLineStyle(2)
            gr_expectedLimit1sigma[(Box,model)].SetLineWidth(3)
            gr_expectedLimit1sigma[(Box,model)].SetLineColor(rt.kBlack)
            gr_expectedLimit2sigma[(Box,model)].SetLineStyle(2)
            gr_expectedLimit2sigma[(Box,model)].SetLineWidth(3)
            gr_expectedLimit2sigma[(Box,model)].SetLineColor(rt.kBlack)
    
    l = rt.TLatex()
    l.SetTextAlign(33)
    l.SetTextSize(0.045)
    l.SetNDC()
    l.SetTextFont(61)
    #l.DrawLatex(0.17,0.92,"CMS")    
    if len(Boxes)>1 and len(models)>1:
        l.DrawLatex(0.3,0.77,"CMS")
    elif len(Boxes)>1:
        l.DrawLatex(0.41,0.835,"CMS")
    else:
        l.DrawLatex(0.905,0.88,"CMS")
        
    l.SetTextFont(52)
    #l.DrawLatex(0.905,0.835,"Preliminary")
    l.SetTextFont(42)
    #l.DrawLatex(0.65,0.92,"%.0f pb^{-1} (13 TeV)"%(options.lumi*1000))
    l.DrawLatex(0.93,0.96,"%.1f fb^{-1} (13 TeV)"%(options.lumi))
    
    if options.model=="Qstar":
        if len(Boxes)>1:
            l.DrawLatex(0.4,0.74,"q* #rightarrow q#gamma")
        else:
            l.DrawLatex(0.89,0.75,"q* #rightarrow q#gamma")
    elif options.model=="Bstar":        
        if len(Boxes)>1:
            l.DrawLatex(0.4,0.74,"b* #rightarrow b#gamma")
        else:
            l.DrawLatex(0.4,0.74,"b* #rightarrow b#gamma")

    #if options.bayes:
    #    if options.noSys:        
    #        l.DrawLatex(0.2,0.85,"Bayesian, no syst.")
    #    else:
    #        l.DrawLatex(0.2,0.85,"Bayesian, with syst.")
    #else:        
    #    if options.noSys:        
    #        l.DrawLatex(0.2,0.85,"Frequentist, no syst.")
    #    else:
    #        l.DrawLatex(0.2,0.85,"Frequentist, with syst.")

    if options.doSignificance:
        c.SetGridy()
        leg = rt.TLegend(0.55,0.79,0.92,0.87)      
    else:        
        leg = rt.TLegend(0.21,0.164,0.58,0.48)
    
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)
    leg.SetFillColorAlpha(0,0)
    leg.SetLineColor(0)
    #leg.SetLineWidth(0.01)
    if not options.doSignificance:
        leg.SetHeader("95% CL upper limits")
    if len(models)==1:
        if options.doSignificance:
            leg.AddEntry(gr_observedLimit[(Box,model)], "Observed","lp")
        else:
            #leg.AddEntry(None,"95% CL limits","")
            #leg.AddEntry(None,"90% CL limits","")
            leg.AddEntry(gr_observedLimit[(Box,model)], "Observed","lp")
            leg.AddEntry(gr_expectedLimit[(Box,model)], "Median expected","l")
        if not options.doSignificance:
            leg.AddEntry(gr_expectedLimit1sigma[(Box,model)], "68% expected","lf")    
        if not options.doSignificance:
            leg.AddEntry(gr_expectedLimit2sigma[(Box,model)], "95% expected","lf")
        for thyModel in thyModelsToDraw:
            leg.AddEntry(xsec_gr_nom[thyModel], legendLabel[thyModel],'l')
    else:
        #leg.AddEntry(None,"95% CL limits","")
        #leg.AddEntry(None,"90% CL limits","")
        for model in models:
            leg.AddEntry(gr_observedLimit[(Box,model)], legendLabel[model],"lp")
            
    leg.Draw("SAME")
        
    if len(thyModelsToDraw)>0 and not options.doSignificance:        
        #legThyModel = rt.TLegend(0.2,0.17,0.55,0.45)
        legThyModel = rt.TLegend(0.45,0.7,0.9,0.92)
        legThyModel2 = rt.TLegend(0.55,0.54,0.9,0.7)            
        legThyModel2.SetTextFont(42)
        legThyModel2.SetFillColor(rt.kWhite)
        legThyModel2.SetLineColor(rt.kWhite)
        legThyModel2.SetFillColorAlpha(0,0)
        legThyModel2.SetLineColorAlpha(0,0)
        legThyModel.SetTextFont(42)
        legThyModel.SetFillColor(rt.kWhite)
        legThyModel.SetLineColor(rt.kWhite)
        legThyModel.SetFillColorAlpha(0,0)
        legThyModel.SetLineColorAlpha(0,0)
            
        #for i, thyModel in enumerate(thyModelsToDraw):
            #if i>4:
                #try:
                    #legThyModel2.AddEntry(xsec_gr_nom[thyModel], legendLabel[thyModel],'l')
                #except:
                    #pass
           # else:
                #legThyModel.AddEntry(xsec_gr_nom[thyModel], legendLabel[thyModel],'l')
        #legThyModel.Draw("same")
        #try:
            #legThyModel2.Draw("same")
        #except:
            #pass
            

    for Box in Boxes:
        for model in models:    
            if options.doSignificance:
                gr_observedLimit[(Box,model)].Draw("lp SAME")
            else:
                if len(models)==1:
                    gr_expectedLimit[(Box,model)].Draw("c same")
                for thyModel in thyModelsToDraw:
                    xsec_gr_nom[thyModel].Draw("c same")
                gr_observedLimit[(Box,model)].Draw("lp SAME")


    #if 'PF' in Box or options.massMax>1600:
        #h_limit.GetXaxis().SetTitle('Resonance mass [TeV]')
        #h_limit.GetXaxis().SetLabelOffset(1000)
        #h_fit_residual_vs_mass.GetXaxis().SetNoExponent()
        #h_fit_residual_vs_mass.GetXaxis().SetMoreLogLabels()    
        #xLab = rt.TLatex()
        #xLab.SetTextAlign(22)
        #xLab.SetTextSize(0.05)
        #xLab.SetTextFont(42)
        #xLab.SetTextSize(0.05)
        #if options.doSignificance:
            #yOffset = -0.138
        #else:
            #yOffset = 6.5e-5 # for 1e-4 min
            #yOffset = 5.25e-6 # for 1e-5 min
        #for i in range(1,8):
            #if i*1000>=options.massMin:
                #xLab.DrawLatex(i*1000, yOffset, "%g"%i)

    #else:
        #h_limit.GetXaxis().SetNdivisions(408,True)

    #if options.box=="CaloDijet2016_PFDijet2016":
        #line1 = rt.TLine(1600,1e4,1600,options.xsecMax)
        #line1 = rt.TLine(1600,1e-1,1600,options.xsecMax)
        #line1.SetLineStyle(2)
        #line1.SetLineWidth(2)
        #line1.SetLineColor(rt.kGray+1)
        #line1.Draw()
        #line2 = rt.TLine(1600,1e-1,1600,2)
        #line2.SetLineStyle(2)
        #line2.SetLineWidth(2)
        #line2.SetLineColor(rt.kGray+1)
        #line2.Draw()
        
          
        #lab = rt.TLatex()
        #lab.SetTextSize(0.035)
        #lab.SetTextFont(42)
        #lab.SetTextColor(rt.kGray+1)
        #lab.SetTextAlign(33)
        #lab.DrawLatex(1600-10,6e4,"#leftarrow")
        #lab.SetTextAlign(13)
        #lab.DrawLatex(1600+10,6e4,"#rightarrow") 
        #lab.SetTextAlign(23)
        #lab.DrawLatex(1600-400,3.5e4,"Low")
        #lab.DrawLatex(1600-400,1.2e4,"mass")
        #lab.DrawLatex(1600+400,3.5e4,"High")
        #lab.DrawLatex(1600+400,1.2e4,"mass")
        
    c.Update()
    #c.SetLogx()    
    c.RedrawAxis() # request from David

    if options.doSignificance:
        c.SaveAs(options.outDir+"/signif_"+options.model+"_"+options.box.lower()+".pdf")
        c.SaveAs(options.outDir+"/signif_"+options.model+"_"+options.box.lower()+".C")
    else:
        if options.bayes:
            if options.noSys:
                c.SaveAs(options.outDir+"/limits_bayes_nosys_"+options.model+"_"+options.box.lower()+".pdf")
                c.SaveAs(options.outDir+"/limits_bayes_nosys_"+options.model+"_"+options.box.lower()+".C")
            else:
                c.SaveAs(options.outDir+"/limits_bayes_"+options.model+"_"+options.box.lower()+".pdf")
                c.SaveAs(options.outDir+"/limits_bayes_"+options.model+"_"+options.box.lower()+".C")
        else:
            if options.noSys:
                c.SaveAs(options.outDir+"/limits_freq_nosys_"+options.model+"_"+options.box.lower()+".pdf")
                c.SaveAs(options.outDir+"/limits_freq_nosys_"+options.model+"_"+options.box.lower()+".C")
            else:
                c.SaveAs(options.outDir+"/limits_freq_"+options.model+"_"+options.SigCoup+"_"+options.box.lower()+".pdf")
                c.SaveAs(options.outDir+"/limits_freq_"+options.model+"_"+options.SigCoup+"_"+options.box.lower()+".png")
                c.SaveAs(options.outDir+"/limits_freq_"+options.model+"_"+options.SigCoup+"_"+options.box.lower()+".root")


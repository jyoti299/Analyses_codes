#! /usr/bin/env python
import ROOT as rt
import os.path
import sys, glob, re
from array import *
from optparse import OptionParser

def getThyXsecDict():    
    thyXsecDict = {}
#    xsecFiles = ['data/Qstar_f1p0_xsec.txt', 'data/Qstar_f0p5_xsec.txt', 'data/Qstar_f0p1_xsec.txt', 'data/Bstar_f1p0_xsec.txt', 'data/Bstar_f0p5_xsec.txt', 'data/Bstar_f0p1_xsec.txt']
####DEta files
    xsecFiles = ['data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_1.0.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_1.2.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_1.5.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_1.8.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_2.0.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_2.2.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_DEta_2.5.txt', 'data/Optimization_80X_ReReco/DEta_Opti/Qstar_f1p0_xsec_NoDEta.txt']
###DPhi Files
    #xsecFiles = ['data/Optimization_80X/DPhi_Opti/Qstar_f1p0_xsec_DPhi_1.5.txt', 'data/Optimization_80X/DPhi_Opti/Qstar_f1p0_xsec_DPhi_2.0.txt', 'data/Optimization_80X/DPhi_Opti/Qstar_f1p0_xsec_DPhi_2.5.txt', 'data/Optimization_80X/DPhi_Opti/Qstar_f1p0_xsec_NoDPhi.txt']
###Phid files
    #xsecFiles = ['data/Optimization_80X/PhId_Opti/Qstar_f1p0_xsec_PhId_Loose.txt', 'data/Optimization_80X/PhId_Opti/Qstar_f1p0_xsec_PhId_Medium.txt', 'data/Optimization_80X/PhId_Opti/Qstar_f1p0_xsec_PhId_Tight.txt' ]
###Jet Id files
    #xsecFiles = ['data/Optimization_80X/JetId_Opti/Qstar_f1p0_xsec_JetId_Tight.txt', 'data/Optimization_80X/JetId_Opti/Qstar_f1p0_xsec_JetId_TightLepVeto.txt']
##CSV files
    #xsecFiles = ['data/Optimization_80X/CSVDisc_Opti/Bstar_f1p0_xsec_CSVL.txt', 'data/Optimization_80X/CSVDisc_Opti/Bstar_f1p0_xsec_CSVM.txt', 'data/Optimization_80X/CSVDisc_Opti/Bstar_f1p0_xsec_CSVT.txt']

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
                    if options.DetectorIndependentLimit:
                        thyXsecDict[thyModel][int(float(lineList[0]))] = float(lineList[j+1]) 
                        #print thyXsecDict[thyModel][int(float(lineList[0]))]
                    else:
                        thyXsecDict[thyModel][int(float(lineList[0]))] = float(lineList[j+1]) * float(lineList[j+2])
                        #print thyXsecDict[thyModel][int(float(lineList[0]))]
        f.close()

        #thyXsecDict['AxigluonkNLO'] = {}
    #for (mass,thyXsec) in thyXsecDict['Axigluon'].iteritems():
        #thyXsecDict['AxigluonkNLO'][mass] = 1.08 * thyXsec
    #for (mass,thyXsec) in thyXsecDict['DM1GeV'].iteritems():
        #thyXsecDict['DM1GeV'][mass] = (5./6.) * thyXsec
    return thyXsecDict

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
    expectedLimit = array('d')
    
    xsecTree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        xsecTree.GetEntry(entry)

        gluinoMassArray.append(xsecTree.mass)
        gluinoMassArray_er.append(0.0)
        
        exec 'xsecULExp = xsecTree.xsecULExp_%s'%Box
                       
        xsecULExp = xsecULExp
        expectedLimit.append(xsecULExp)#*crossSections[i])
                    
    return gluinoMassArray, gluinoMassArray_er, expectedLimit


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
   # parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
    #              help="Input/Output directory to store output")    
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
    parser.add_option('--Independent',dest="DetectorIndependentLimit",default=False,action='store_true',
                  help="Get (sigma X BR) limit instead of (sigma X BR X Acc X Eff)")
    
    (options,args) = parser.parse_args()
    Boxes = options.box.split('_')
    models = options.model.split('_')
    model = models[0]
    #directory      = options.outDir
    Box = Boxes[0]
    box = Box.lower()
    coupling = options.SigCoup
        
    #thyXsecDict = getThyXsecDict() 
    #thyModels = thyXsecDict.keys()

    ### If we want to draw q* xsec also along with expected and obs limits, then the models in thyXsecDict and thyModelsToDraw should be same.
    thyModelsToDraw = []
    thyModelsToDraw = ['Exp. Limit with Bias', 'Exp. Limit without Bias']

    lineStyle = {'Exp. Limit with Bias':2,
                 'Exp. Limit without Bias':2,
                 'None':1               
                 }
        
    lineColor = {'Exp. Limit with Bias':rt.kRed+1,
                 'Exp. Limit without Bias':rt.kBlue+1,
                 'None':1
                 }

    markerColor = {'Exp. Limit with Bias':rt.kRed+1,
                 'Exp. Limit without Bias':rt.kBlue+1,
                 'None':1
                 }
        
    markerStyle = {'Exp. Limit with Bias':1,
                 'Exp. Limit without Bias':1
                 }
        
    legendLabel = {'Exp. Limit with Bias':'Exp. Limit with Bias',
                 'Exp. Limit without Bias':'Exp. Limit without Bias'
                 }

    Directory = {'Exp. Limit with Bias':'Cards_Qstar_f-1p0_ExcitedQuarks2016_ExpLimitFromData_withBias/',
                 'Exp. Limit without Bias':'Cards_Qstar_f-1p0_ExcitedQuarks2016_ExpLimitFromData/'
                 }
    
    setstyle()
    rt.gStyle.SetOptStat(0)
    c = rt.TCanvas("c","c",800,800)
    if options.doSignificance:
        c.SetLogy(0)
    else:        
        c.SetLogy()

    h_limit = rt.TMultiGraph()
    gr_expectedLimit = {}
    gluinoMassArray = {}
    gluinoMassArray_er = {}
    expectedLimit = {}
    
    if options.doSignificance:
        h_limit.SetTitle(" ;Resonance mass [GeV];Local significance n#sigma")
    else:
        if options.DetectorIndependentLimit:
            h_limit.SetTitle(" ;Resonance mass [GeV]; #sigma #times #bf{#it{#Beta}} [pb]")
        else:
            h_limit.SetTitle(" ;Resonance mass [GeV]; #sigma #times #bf{#it{#Beta}} #times #bf{#it{#Alpha}} #times #bf{#it{#epsilon}} [pb]")

    for thyModel in thyModelsToDraw:
        gluinoMassArray[thyModel], gluinoMassArray_er[thyModel], expectedLimit[thyModel] = getHybridCLsArrays(Directory[thyModel], model, Box, options.bayes)
        print expectedLimit[thyModel]        
        nPoints = len(expectedLimit[thyModel])

        gr_expectedLimit[thyModel] = rt.TGraph(nPoints, gluinoMassArray[thyModel], expectedLimit[thyModel])
        gr_expectedLimit[thyModel].SetMarkerColor(markerColor[thyModel])
        gr_expectedLimit[thyModel].SetMarkerStyle(markerStyle[thyModel])
        gr_expectedLimit[thyModel].SetLineColor(lineColor[thyModel])
        gr_expectedLimit[thyModel].SetLineStyle(lineStyle[thyModel])
        gr_expectedLimit[thyModel].SetMarkerSize(0.8)
        gr_expectedLimit[thyModel].SetLineWidth(2)

        h_limit.Add(gr_expectedLimit[thyModel])
                
    h_limit.Draw("a3")
    h_limit.GetXaxis().SetLimits(options.massMin,options.massMax)

    if options.doSignificance:
        h_limit.SetMaximum(4)
        h_limit.SetMinimum(0)
    else:
        h_limit.SetMaximum(options.xsecMax)
        h_limit.SetMinimum(options.xsecMin)
            
    h_limit.Draw("a3")
    if options.doSignificance:
        h_limit.GetYaxis().SetNdivisions(405,True)

    for thyModel in thyModelsToDraw:
        gr_expectedLimit[thyModel].Draw("lp SAME")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.045)
    l.SetNDC()
    l.SetTextFont(62)
    #l.DrawLatex(0.17,0.92,"CMS")    
    if len(Boxes)>1 and len(models)>1:
        l.DrawLatex(0.3,0.77,"CMS")
    elif len(Boxes)>1:
        l.DrawLatex(0.41,0.835,"CMS")
    else:
        l.DrawLatex(0.18,0.918,"CMS")
        
    l.SetTextFont(52)
    #l.DrawLatex(0.28,0.92,"Preliminary")
    l.SetTextFont(42)
    #l.DrawLatex(0.65,0.92,"%.0f pb^{-1} (13 TeV)"%(options.lumi*1000))
#    l.DrawLatex(0.63,0.92,"%.1f fb^{-1} (13 TeV)"%(options.lumi))
    
    if options.model=="Qstar":
        if len(Boxes)>1:
            l.DrawLatex(0.28,0.77,"q* #rightarrow q#gamma")
        else:
            l.DrawLatex(0.28,0.77,"q* #rightarrow q#gamma")
    elif options.model=="Bstar":        
        if len(Boxes)>1:
            l.DrawLatex(0.28,0.77,"b* #rightarrow b#gamma")
        else:
            l.DrawLatex(0.28,0.77,"b* #rightarrow b#gamma")

    if options.doSignificance:
        c.SetGridy()
        leg = rt.TLegend(0.51,0.79,0.90,0.90)      
    else:        
        leg = rt.TLegend(0.55,0.68,0.92,0.87)#0.51,0.79,0.90,0.90)

    leg.SetTextFont(42)
    leg.SetFillColorAlpha(0,0)
    leg.SetLineColor(0)
    leg.SetHeader("95% CL expected limits")
    for thyModel in thyModelsToDraw:
        leg.AddEntry(gr_expectedLimit[thyModel], legendLabel[thyModel],'lp')
            
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
            
    for thyModel in thyModelsToDraw:
        gr_expectedLimit[thyModel].Draw("lp same")        

    c.Update()
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
                c.SaveAs("Cards_Qstar_f-1p0_ExcitedQuarks2016_ExpLimitFromData_withBias/Explimits_comp_Bias_freq_"+options.model+"_"+options.box.lower()+".pdf")
                c.SaveAs("Cards_Qstar_f-1p0_ExcitedQuarks2016_ExpLimitFromData_withBias/Explimits_comp_Bias_freq_"+options.model+"_"+options.box.lower()+".png")
                c.SaveAs("Cards_Qstar_f-1p0_ExcitedQuarks2016_ExpLimitFromData_withBias/Explimits_comp_Bias_freq_"+options.model+"_"+options.box.lower()+".root")


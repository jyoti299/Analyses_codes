#!/bin/tcsh

################################################################
## Postanalyzer for q* -> gamma + q study at 13 TeV
## author : Varun Sharma
## Email : varun.sharma@cern.ch
################################################################

setenv pwd $PWD

set isM = 1 ## set 1 for MC
set isD = 0 ## set 1 for Data
#-----------------------------------------------------------

## 0=Qstarf1p0,  1=Qstarf0p5,  2=Qstarf0p1,  3=Bkg GJ,   4=Bkg Dijet,   5=Bkg GJ Pythia8, 6=EWK, 7 = TESTING
foreach runCase ( 0 1 2 3 4 6)
#foreach runCase ( 7)

set sampleIndex = 0

###----------------------------------------------------------------------------------------------------
if( ${runCase} == 0 ) then
echo "******** Running for Signal -- Qstar, f=1.0 ***********"
setenv tmp QstarToGJ_
foreach i (${tmp}M500_f1p0 ${tmp}M1000_f1p0 ${tmp}M2000_f1p0 ${tmp}M3000_f1p0 ${tmp}M4000_f1p0 ${tmp}M5000_f1p0 ${tmp}M6000_f1p0 ${tmp}M7000_f1p0 ${tmp}M8000_f1p0 ${tmp}M9000_f1p0)
set XS = ( 3.033E2         1.632E1          5.213E-1         4.272E-2         4.8E-3           5.835E-4         7.076E-5         8.66E-6          1.283E-6         2.985E-7 )
set totEvnts = (298953     299072           100000           74879            74271            74802            74822            49544            49628            50000 )
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/Signal/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 1 ) then
echo "******** Running for Signal -- Qstar, f=0.5 ***********"
setenv tmp QstarToGJ_
foreach i (${tmp}M500_f0p5 ${tmp}M1000_f0p5 ${tmp}M2000_f0p5 ${tmp}M3000_f0p5 ${tmp}M4000_f0p5 ${tmp}M5000_f0p5 ${tmp}M6000_f0p5 ${tmp}M7000_f0p5 ${tmp}M8000_f0p5 ${tmp}M9000_f0p5)
set XS = ( 7.378E1         4.129E0          1.328E-1         1.095E-2         1.212E-3          1.437E-4        1.62E-5          1.672E-6         1.647E-7         2.329E-8)
set totEvnts = (297624     97970            98740            74232            74857             74896           74328            49888            48200            49403)
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/Signal/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 2 ) then
echo "******** Running for Signal -- Qstar, f=0.1 ***********"
setenv tmp QstarToGJ_
foreach i (${tmp}M500_f0p1 ${tmp}M1000_f0p1 ${tmp}M2000_f0p1 ${tmp}M3000_f0p1 ${tmp}M4000_f0p1 ${tmp}M5000_f0p1 ${tmp}M6000_f0p1 ${tmp}M7000_f0p1 ${tmp}M8000_f0p1 ${tmp}M9000_f0p1)
set XS = ( 2.955E0         1.655E-1         5.315E-3         4.356E-4         4.861E-5         5.715E-6         6.241E-7         5.973E-8         4.515E-9         2.655E-10   )
set totEvnts = (97208      98290            73128            74286 	      74836            75000            49918            49679            49316            49782 )
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/Signal/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 3 ) then
echo "******** Running for Background -- photon + jet, madgraphMLM-pythia8, HT binned ***********"
foreach i (GJets_HT40To100  GJets_HT100To200  GJets_HT200To400  GJets_HT400To600  GJets_HT600ToInf)
set XS = ( 23080.0          9110.0            2281.0            273.0             94.5)
set totEvnts = (4816232     5027618           10424189          2476770           2550765)
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/GJets/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 4 ) then
echo "******** Running for Background -- dijet, pythia8, PT binned ***********"
foreach i (QCD_Pt_120to170 QCD_Pt_170to300 QCD_Pt_300to470 QCD_Pt_470to600 QCD_Pt_600to800 QCD_Pt_800to1000 QCD_Pt_1000to1400 QCD_Pt_1400to1800 QCD_Pt_1800to2400 QCD_Pt_2400to3200 QCD_Pt_3200toInf)
set XS = ( 471100.0        117276.0        7823.0          648.2           186.9               32.293       9.4183            0.84265           0.114943   	  0.00682981        0.000165445 )
set totEvnts = (3458385    3364368         2935633         1937537         1964128         1937216          1487136           197959            193608            194456            192944 )
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/DiJet/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 5 ) then
echo "******** Running for Background -- photon + jet, Pt-15-6000 Flat Pythia8 ***********"
foreach i ( GJet_Flat15To6000_Pythia )
set XS = ( 365896.0)
set totEvnts = (9891424)
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 6 ) then
echo "******** Running for EWK -- WGToLNuG, WJetsToLNu, ZGTo2LG ***********"
foreach i (     WGToLNuG   WJetsToLNu  ZGTo2LG )
set XS = (      377.5      50360       123.7   )
set totEvnts = (6099599    72220670    4451319 )
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 5000 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 7 ) then
echo "******** Running for Testing -- Qstar, f=1.0 ***********"
setenv tmp QstarToGJ_
foreach i ( ${tmp}M1000_f1p0 )
set XS = (      1.632E1      )
set totEvnts = (299072       )
setenv sourceDir /store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/Signal/${i}/
setenv InputFilesPath root://cmseos.fnal.gov/${sourceDir}
set filesPerJob = 20 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------

setenv FileNameTag ${i}

@ sampleIndex = ${sampleIndex} + 1

##_______________________________________________________________________________
## Read each file in the directory and submit seperate job for each file OR
## merge "r" files to submit one job
##-----------------------------------------------------------------
       
#Fill the range for input to code e.g 1-10 then r = 10
set r = ${filesPerJob}
set sf = 1          ## start file
set ef = ${r}       ## end file
 
eosls ${sourceDir} > ${FileNameTag}_dataset.txt
 
setenv datafile  ${FileNameTag}_dataset.txt
set file_length=`wc -l $datafile | cut -c1-4`     # Count Lines i.e. no. of files in dir
set p = 0
 
set Tot = ${file_length}
#set Tot = 1
 
#run till the last line of the input files
while (${sf} <= ${Tot})
 
@ p = ${p} + 1
set DataName=`tail -n +$p ${datafile} | head -1`  # Get name of dataset
 
#####More changes in constructor, in the outputfile and at the end of file--------------------------------


cat>PhoJet_Analyzer_MC.C<<EOF
#define PhoJet_Analyzer_MC_cxx
#include "PhoJet_Analyzer_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;
using namespace ROOT;

int main()
{
    gROOT->ProcessLine("#include <vector>");
    gROOT->ProcessLine("#include <map>");
    PhoJet_Analyzer_MC a;
    a.Loop();
}

void PhoJet_Analyzer_MC::Loop()
{
    //----------------------------------------------------------------
    // Give values for all cuts at the beginning
    //----------------------------------------------------------------
    // Luminosity
    Lumi = 2670.555; // pb^-1
    isMC = ${isM};
    isDATA = ${isD};

    // Fiducial Cuts
    cut_photonPt  = 190.0;  // (GeV)
    cut_photonEta = 1.4442; 

    cut_jetPt     = 190.0;  // (GeV)
    cut_jetEta    = 2.4;

    cut_gjDPhi = 2.5;
    cut_gjDEta = 2.0;

    cut_gjMass = 560.0;

    //Output File :
    TString OutputFile = "${FileNameTag}_${p}";
    f1 = new TFile(OutputFile+".root","recreate");
    f1->cd();
    //-------------------------------------------------------------------  


    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    std::cout<< " Will Analyze = " << nentries << " events" <<std::endl;
    std::cout<< " Total Events = $totEvnts[${sampleIndex}].0 events" <<std::endl;

    BookHistos();
    if(isMC) LumiReWeighting();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;

	Double_t pu_weight = 1.0;
	Int_t    indexBX0  = -1;

	if(isMC){
	    //onlyEvtWeight = Lumi*($XS[${sampleIndex}]/$totEvnts[${sampleIndex}].0);
	    onlyEvtWeight = Lumi*($XS[${sampleIndex}]/(double)nentries);

	    for(int iPU = 0; iPU < nPUInfo; ++iPU)
		if( puBX->at(iPU) == 0) indexBX0 = iPU;

	    pu_weight = puweight(puTrue->at(indexBX0));
	}
	if(isDATA) onlyEvtWeight = 1;
	EvtWeight = pu_weight*onlyEvtWeight ;

	PC = -1;
	JC = -1;
	goodVertex = 0;

	Bool_t passHLT = false;
	Bool_t primaryVtx = false;
	primaryVtx  = PrimaryVertex(goodVertex);
	
	h_Vertices[0]->Fill(goodVertex);
	h_Vertices[1]->Fill(goodVertex, pu_weight);
	h_trueInteractions[0]->Fill(puTrue->at(indexBX0));
	h_trueInteractions[1]->Fill(puTrue->at(indexBX0), pu_weight);

	//++++++++++++++++++++++++++++++++ PHOTON SELECTION ++++++++++++++++++++++++++++++++++++++++
	// Selecting photon candidate using loose/medium/tight photon ID.
	foundIsoPhoton.clear();
	for(int ipho=0; ipho<nPho; ++ipho){
	    if(CutBasedPFPhotonID(ipho, "loose") && NoSpike(ipho) &&(*phoEt)[ipho] > 30.0)
		foundIsoPhoton.push_back(ipho);
	}

	foundPhoton.clear();
	for( int isoPho = 0; isoPho < foundIsoPhoton.size(); ++isoPho){
	    if( (*phoEt)[foundIsoPhoton[isoPho]] > cut_photonPt  && fabs((*phoSCEta)[foundIsoPhoton[isoPho]]) < cut_photonEta)
		foundPhoton.push_back(foundIsoPhoton[isoPho]);
	}
	if(foundPhoton.size() != 0)PC = foundPhoton[0];

	h_nIsoPhoton[0]->Fill(foundIsoPhoton.size(), EvtWeight);  // Tot # of isolated photons
	h_nPhoton[0]->Fill(foundPhoton.size(), EvtWeight);        // Tot # of isolated photons with pt > cut and eta < cut
	for(Int_t i_p = 0 ; i_p < foundIsoPhoton.size() ; ++i_p)
	    h_nIsoPhotonPt[0]->Fill((*phoEt)[foundIsoPhoton[i_p]], i_p+1, EvtWeight);
	for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p)
	    h_nPhotonPt[0]->Fill((*phoEt)[foundPhoton[i_p]], i_p+1, EvtWeight);
	
	// ++++++++++++++++++++++++++++++++++++= JET SELECTION ++++++++++++++++++++++++++++++++
	// Selecting Jet candidate using loose/tight/tightLepVeto Jet ID    [tightLepVeto not there in current version]
	foundJet.clear();
	if(foundPhoton.size() != 0){
	    for(int ijet=0; ijet<nJet; ++ijet){
		if( getDR((*phoSCEta)[PC], (*jetEta)[ijet], (*phoSCPhi)[PC], (*jetPhi)[ijet] ) > 0.5 && jetID(ijet, "tight") && (*jetPt)[ijet] > 30.0 )
		    foundJet.push_back(ijet);
	    }
	}
	if(foundJet.size() != 0)JC = foundJet[0];

	h_nJet[0]->Fill(foundJet.size(), EvtWeight);
	for(Int_t i_j = 0 ; i_j < foundJet.size() ; ++i_j)
	    h_nJetPt[0]->Fill((*jetPt)[foundJet[i_j]], i_j+1, EvtWeight);
	
	Bool_t c_jetPtEta = false;
	Bool_t c_gjDPhi = false;
	Bool_t c_gjDEta = false;
	Bool_t c_InvMass = false;
	Bool_t c_InvMass_1p2to2p2 = false;

	if(PC > -1 && JC > -1){
	    c_jetPtEta = ( (*jetPt)[JC] > cut_jetPt  && fabs((*jetEta)[JC]) < cut_jetEta );
	    c_gjDPhi   = ( getDPhi( (*phoSCPhi)[PC], (*jetPhi)[JC] ) > cut_gjDPhi );
	    c_gjDEta   = ( getDEta( (*phoSCEta)[PC], (*jetEta)[JC] ) < cut_gjDEta );
	    c_InvMass  = ( getMass(PC,JC) > cut_gjMass );
	    
	    c_InvMass_1p2to2p2  = ( getMass(PC,JC) > 1807.0 && getMass(PC,JC) < 2102.0 );
	}


	h_CutFlow->Fill(0.5);
	h_CutExpFlow->Fill(0.5,EvtWeight);

	if(passHLT || isMC ){
	    h_CutFlow->Fill(1.5);
	    h_CutExpFlow->Fill(1.5,EvtWeight);

	    if(primaryVtx){
		h_CutFlow->Fill(2.5);
		h_CutExpFlow->Fill(2.5,EvtWeight);

		if(foundIsoPhoton.size() > 0){
		    h_CutFlow->Fill(3.5);
		    h_CutExpFlow->Fill(3.5,EvtWeight);

		    if(PC > -1){
			h_CutFlow->Fill(4.5);
			h_CutExpFlow->Fill(4.5,EvtWeight);

			if(JC > -1){
			    h_CutFlow->Fill(5.5);
			    h_CutExpFlow->Fill(5.5,EvtWeight);

			    if(c_jetPtEta){
				h_CutFlow->Fill(6.5);
				h_CutExpFlow->Fill(6.5,EvtWeight);

				if(c_gjDPhi){
				    h_CutFlow->Fill(7.5);
				    h_CutExpFlow->Fill(7.5,EvtWeight);

				    if(c_gjDEta){
					h_CutFlow->Fill(8.5);
					h_CutExpFlow->Fill(8.5,EvtWeight);

					h_PC             ->Fill( PC, EvtWeight);
					h_JC             ->Fill( JC, EvtWeight);

					h_nIsoPhoton[1]   ->Fill(foundIsoPhoton.size(), EvtWeight);
					h_nPhoton[1]      ->Fill(foundPhoton.size(), EvtWeight);
					for(Int_t i_p = 0 ; i_p < foundIsoPhoton.size() ; ++i_p)
					    h_nIsoPhotonPt[1]->Fill((*phoEt)[foundIsoPhoton[i_p]], i_p+1, EvtWeight);
					for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p)
					    h_nPhotonPt[1]->Fill((*phoEt)[foundPhoton[i_p]], i_p+1, EvtWeight);
	
					h_nJet[1]         ->Fill(foundJet.size(), EvtWeight);
					for(Int_t i_j = 1 ; i_j < foundJet.size() ; ++i_j)
					    h_nJetPt[1]     ->Fill((*jetPt)[foundJet[i_j]], i_j+1, EvtWeight);

					h_ptPhoton_noPU[0]    ->Fill((*phoEt)[PC] ,onlyEvtWeight);
					h_ptJet_noPU[0]       ->Fill((*jetPt)[JC] ,onlyEvtWeight);
					h_mass_VarBin_noPU[0] ->Fill(getMass(PC,JC), onlyEvtWeight);

					h_ptPhoton[0]    ->Fill((*phoEt)[PC] ,EvtWeight);
					h_ptJet[0]       ->Fill((*jetPt)[JC] ,EvtWeight);
					h_etaPhoton[0]   ->Fill((*phoSCEta)[PC] ,EvtWeight);
					h_etaJet[0]      ->Fill((*jetEta)[JC] ,EvtWeight);
					h_phiPhoton[0]   ->Fill((*phoSCPhi)[PC] ,EvtWeight);
					h_phiJet[0]      ->Fill((*jetPhi)[JC] ,EvtWeight);
					h_mass_VarBin[0] ->Fill(getMass(PC,JC), EvtWeight);
					h_mass_bin1[0]   ->Fill(getMass(PC,JC), EvtWeight);

					h_PFMet[0]       ->Fill(pfMET,EvtWeight);
					h_SumPFMet[0]    ->Fill(pfMETsumEt,EvtWeight);
					h_MetBySumMET[0] ->Fill(pfMET/pfMETsumEt, EvtWeight);

					h_pfMetgjMass[0] ->Fill(getMass(PC,JC),pfMET, EvtWeight);
					h_PtPhotJet[0]   ->Fill((*phoEt)[PC], (*jetPt)[JC], EvtWeight);
					h_etaPhotJet[0]  ->Fill((*phoSCEta)[PC], (*jetEta)[JC], EvtWeight);
					h_DR_PhotonJet[0]->Fill( getDR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC] ), EvtWeight);
										
					h_cosThetaStar[0]->Fill( getCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC] ), EvtWeight);
					h_dEta[0]        ->Fill( getDEta( (*phoSCEta)[PC], (*jetEta)[JC] ) , EvtWeight);
					h_dphi[0]        ->Fill( getDPhi( (*phoSCPhi)[PC], (*jetPhi)[JC] ) , EvtWeight);

					h_Photon_SigmaIetaIeta[0]->Fill( (*phoSigmaIEtaIEtaFull5x5)[PC], EvtWeight);
					h_HoE[0]                 ->Fill( (*phoHoverE)[PC], EvtWeight);
					h_CorrPFiso_Charged[0]   ->Fill( TMath::Max( ( (*phoPFChIso)[PC]  - rho*EAcharged((*phoSCEta)[PC]) ), 0.0), EvtWeight);
					h_CorrPFiso_Neutral[0]   ->Fill( TMath::Max( ( (*phoPFNeuIso)[PC] - rho*EAneutral((*phoSCEta)[PC]) ), 0.0), EvtWeight);
					h_CorrPFiso_Photon[0]    ->Fill( TMath::Max( ( (*phoPFPhoIso)[PC] - rho*EAphoton((*phoSCEta)[PC])  ), 0.0), EvtWeight);
					h_PFiso_Electronveto[0]  ->Fill( (*phoEleVeto)[PC], EvtWeight);

					h_jet_NEF[0]             ->Fill( (*jetNEF)[JC], EvtWeight);
					h_jet_NHF[0]             ->Fill( (*jetNHF)[JC], EvtWeight);
					h_jet_CEF[0]             ->Fill( (*jetCEF)[JC], EvtWeight);
					h_jet_CHF[0]             ->Fill( (*jetCHF)[JC], EvtWeight);
					h_jet_NConstituents[0]   ->Fill( (*jetNConstituents)[JC], EvtWeight);
					h_jet_ChargeMultiplicity[0]->Fill( (*jetNCH)[JC], EvtWeight);


					if(c_InvMass){
					    h_CutFlow->Fill(9.5);
					    h_CutExpFlow->Fill(9.5,EvtWeight);

					    h_Vertices[2]     ->Fill(goodVertex);
					    h_Vertices[3]     ->Fill(goodVertex, pu_weight);
					    h_nIsoPhoton[2]   ->Fill(foundIsoPhoton.size(), EvtWeight);
					    h_nPhoton[2]      ->Fill(foundPhoton.size(), EvtWeight);
					    for(Int_t i_p = 0 ; i_p < foundIsoPhoton.size() ; ++i_p)
						h_nIsoPhotonPt[2]->Fill((*phoEt)[foundIsoPhoton[i_p]], i_p+1, EvtWeight);
					    for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p)
						h_nPhotonPt[2]->Fill((*phoEt)[foundPhoton[i_p]], i_p+1, EvtWeight);
					    
					    h_nJet[2]         ->Fill(foundJet.size(), EvtWeight);
					    for(Int_t i_j = 1 ; i_j < foundJet.size() ; ++i_j)
						h_nJetPt[2]     ->Fill((*jetPt)[foundJet[i_j]], i_j+1, EvtWeight);

					    h_ptPhoton_noPU[1]    ->Fill((*phoEt)[PC] ,onlyEvtWeight);
					    h_ptJet_noPU[1]       ->Fill((*jetPt)[JC] ,onlyEvtWeight);
					    h_mass_VarBin_noPU[1] ->Fill(getMass(PC,JC), onlyEvtWeight);

					    h_ptPhoton[1]    ->Fill((*phoEt)[PC] ,EvtWeight);
					    h_ptJet[1]       ->Fill((*jetPt)[JC] ,EvtWeight);
					    h_etaPhoton[1]   ->Fill((*phoSCEta)[PC] ,EvtWeight);
					    h_etaJet[1]      ->Fill((*jetEta)[JC] ,EvtWeight);
					    h_phiPhoton[1]   ->Fill((*phoSCPhi)[PC] ,EvtWeight);
					    h_phiJet[1]      ->Fill((*jetPhi)[JC] ,EvtWeight);
					    h_mass_VarBin[1] ->Fill(getMass(PC,JC), EvtWeight);
					    h_mass_bin1[1]   ->Fill(getMass(PC,JC), EvtWeight);

					    h_PFMet[1]       ->Fill(pfMET,EvtWeight);
					    h_SumPFMet[1]    ->Fill(pfMETsumEt,EvtWeight);
					    h_MetBySumMET[1] ->Fill(pfMET/pfMETsumEt, EvtWeight);

					    h_pfMetgjMass[1] ->Fill(getMass(PC,JC),pfMET, EvtWeight);
					    h_PtPhotJet[1]   ->Fill((*phoEt)[PC], (*jetPt)[JC], EvtWeight);
					    h_etaPhotJet[1]  ->Fill((*phoSCEta)[PC], (*jetEta)[JC], EvtWeight);
					    h_DR_PhotonJet[1]->Fill( getDR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC] ), EvtWeight);

					    h_cosThetaStar[1]->Fill( getCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC] ), EvtWeight);
					    h_dEta[1]        ->Fill(getDEta( (*phoSCEta)[PC], (*jetEta)[JC] ) , EvtWeight);
					    h_dphi[1]        ->Fill(getDPhi( (*phoSCPhi)[PC], (*jetPhi)[JC] ) , EvtWeight);

					    h_Photon_SigmaIetaIeta[1]->Fill( (*phoSigmaIEtaIEtaFull5x5)[PC], EvtWeight);
					    h_HoE[1]                 ->Fill( (*phoHoverE)[PC], EvtWeight);
					    h_CorrPFiso_Charged[1]   ->Fill( TMath::Max( ( (*phoPFChIso)[PC]  - rho*EAcharged((*phoSCEta)[PC]) ), 0.0), EvtWeight);
					    h_CorrPFiso_Neutral[1]   ->Fill( TMath::Max( ( (*phoPFNeuIso)[PC] - rho*EAneutral((*phoSCEta)[PC]) ), 0.0), EvtWeight);
					    h_CorrPFiso_Photon[1]    ->Fill( TMath::Max( ( (*phoPFPhoIso)[PC] - rho*EAphoton((*phoSCEta)[PC])  ), 0.0), EvtWeight);
					    h_PFiso_Electronveto[1]  ->Fill( (*phoEleVeto)[PC], EvtWeight);

					    h_jet_NEF[1]             ->Fill( (*jetNEF)[JC], EvtWeight);
					    h_jet_NHF[1]             ->Fill( (*jetNHF)[JC], EvtWeight);
					    h_jet_CEF[1]             ->Fill( (*jetCEF)[JC], EvtWeight);
					    h_jet_CHF[1]             ->Fill( (*jetCHF)[JC], EvtWeight);
					    h_jet_NConstituents[1]   ->Fill( (*jetNConstituents)[JC], EvtWeight);
					    h_jet_ChargeMultiplicity[1]->Fill( (*jetNCH)[JC], EvtWeight);

					    if(c_InvMass_1p2to2p2){
						h_CutFlow->Fill(10.5);
						h_CutExpFlow->Fill(10.5,EvtWeight);

						h_nIsoPhoton[3]   ->Fill(foundIsoPhoton.size(), EvtWeight);
						h_nPhoton[3]      ->Fill(foundPhoton.size(), EvtWeight);
						for(Int_t i_p = 0 ; i_p < foundIsoPhoton.size() ; ++i_p)
						    h_nIsoPhotonPt[3]->Fill((*phoEt)[foundIsoPhoton[i_p]], i_p+1, EvtWeight);
						for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p)
						    h_nPhotonPt[3]->Fill((*phoEt)[foundPhoton[i_p]], i_p+1, EvtWeight);

						h_nJet[3]         ->Fill(foundJet.size(), EvtWeight);
						for(Int_t i_j = 1 ; i_j < foundJet.size() ; ++i_j)
						    h_nJetPt[3]     ->Fill((*jetPt)[foundJet[i_j]], i_j+1, EvtWeight);

						h_ptPhoton[2]    ->Fill((*phoEt)[PC] ,EvtWeight);
						h_ptJet[2]       ->Fill((*jetPt)[JC] ,EvtWeight);
						h_etaPhoton[2]   ->Fill((*phoSCEta)[PC] ,EvtWeight);
						h_etaJet[2]      ->Fill((*jetEta)[JC] ,EvtWeight);
						h_phiPhoton[2]   ->Fill((*phoSCPhi)[PC] ,EvtWeight);
						h_phiJet[2]      ->Fill((*jetPhi)[JC] ,EvtWeight);
						h_mass_VarBin[2] ->Fill(getMass(PC,JC), EvtWeight);
						h_mass_bin1[2]   ->Fill(getMass(PC,JC), EvtWeight);

						h_PFMet[2]       ->Fill(pfMET,EvtWeight);
						h_SumPFMet[2]    ->Fill(pfMETsumEt,EvtWeight);
						h_MetBySumMET[2] ->Fill(pfMET/pfMETsumEt, EvtWeight);

						h_pfMetgjMass[2] ->Fill(getMass(PC,JC),pfMET, EvtWeight);
						h_PtPhotJet[2]   ->Fill((*phoEt)[PC], (*jetPt)[JC], EvtWeight);
						h_etaPhotJet[2]  ->Fill((*phoSCEta)[PC], (*jetEta)[JC], EvtWeight);
						h_DR_PhotonJet[2]->Fill( getDR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC] ), EvtWeight);

						h_cosThetaStar[2]->Fill( getCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC] ), EvtWeight);
						h_dEta[2]        ->Fill(getDEta( (*phoSCEta)[PC], (*jetEta)[JC] ) , EvtWeight);
						h_dphi[2]        ->Fill(getDPhi( (*phoSCPhi)[PC], (*jetPhi)[JC] ) , EvtWeight);

						h_Photon_SigmaIetaIeta[2]->Fill( (*phoSigmaIEtaIEtaFull5x5)[PC], EvtWeight);
						h_HoE[2]                 ->Fill( (*phoHoverE)[PC], EvtWeight);
						h_CorrPFiso_Charged[2]   ->Fill( TMath::Max( ( (*phoPFChIso)[PC]  - rho*EAcharged((*phoSCEta)[PC]) ), 0.0), EvtWeight);
						h_CorrPFiso_Neutral[2]   ->Fill( TMath::Max( ( (*phoPFNeuIso)[PC] - rho*EAneutral((*phoSCEta)[PC]) ), 0.0), EvtWeight);
						h_CorrPFiso_Photon[2]    ->Fill( TMath::Max( ( (*phoPFPhoIso)[PC] - rho*EAphoton((*phoSCEta)[PC])  ), 0.0), EvtWeight);
						h_PFiso_Electronveto[2]  ->Fill( (*phoEleVeto)[PC], EvtWeight);

						h_jet_NEF[2]             ->Fill( (*jetNEF)[JC], EvtWeight);
						h_jet_NHF[2]             ->Fill( (*jetNHF)[JC], EvtWeight);
						h_jet_CEF[2]             ->Fill( (*jetCEF)[JC], EvtWeight);
						h_jet_CHF[2]             ->Fill( (*jetCHF)[JC], EvtWeight);
						h_jet_NConstituents[2]   ->Fill( (*jetNConstituents)[JC], EvtWeight);
						h_jet_ChargeMultiplicity[2]->Fill( (*jetNCH)[JC], EvtWeight);

					    }//--c_InvMass_1p2to2p2 
					}//--c_InvMass
				    }//--c_gjDEta
				}//--c_gjDPhi
			    }//--c_jetPtEta
			}//--JC
		    }//--PC
		}//--foundIsoPhoton
	    }//--primaryVtx
	}//--passHLT

    }//-- jentry for loop

    //WriteHistos(); 
}// -- Loop()
EOF
##################################################
###---------- .C File ends --------------------###

cat>PhoJet_Analyzer_MC.h<<EOF
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 28 18:07:04 2016 by ROOT version 6.02/13
// from TTree EventTree/Event data
// found on file: /eos/uscms/store/user/leptonjets/varun/13TeV/Ntuples/MC_76X/Signal/QstarToGJ_M1000_f1p0/ntuple_QstarToGJ_M-1000_f-1p0_76X_49.root
//////////////////////////////////////////////////////////

#ifndef PhoJet_Analyzer_MC_h
#define PhoJet_Analyzer_MC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h> 
#include <TGraphAsymmErrors.h>
#include <map>
#include "TRFIOFile.h"
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <Riostream.h>
#include <set>
#include "TKDE.h"
#include <TLorentzVector.h>

#ifdef __MAKECINT__ 
#pragma link C++ class vector<bool>+;
#pragma link C++ class map<TString,TH1D*>+;
#pragma link C++ class vector<unsigned long long>+;
#pragma link C++ class vector<vector<unsigned long long> >+;
#pragma link C++ class vector<ULong64_t>+;
#pragma link C++ class vector<vector<ULong64_t> >+;
#pragma link C++ class vector<long long>+;
#pragma link C++ class vector<vector<long long> >+;
#pragma link C++ class vector<Long64_t>+;
#pragma link C++ class vector<vector<Long64_t> >+;
#endif

using namespace std;        
using namespace ROOT;


class PhoJet_Analyzer_MC {
    public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Variables defined by me :
	TFile *f1;

	Bool_t isMC, isDATA;
	Bool_t _debug;

	Double_t EvtWeight;
	Double_t onlyEvtWeight;
	Double_t Lumi;

	Double_t cut_photonPt;
	Double_t cut_photonEta;
	Double_t cut_jetPt;
	Double_t cut_jetEta;

	Double_t cut_gjDPhi;
	Double_t cut_gjDEta;

	Double_t cut_gjMass;

	std::vector<int> foundIsoPhoton;
	std::vector<int> foundPhoton;
	std::vector<int> foundJet;

	Int_t goodVertex;
	Int_t PC;
	Int_t JC;

	// Declare user defined histograms below :
	TH1F * h_ptPhoton_noPU[2];
	TH1F * h_ptJet_noPU[2];
	TH1F * h_mass_VarBin_noPU[2];

	TH1F * h_ptPhoton[3];
	TH1F * h_ptJet[3];
	TH1F * h_mass_VarBin[3];
	TH1F * h_mass_bin1[3];
	TH1F * h_etaPhoton[3];
	TH1F * h_phiPhoton[3];
	TH1F * h_etaJet[3];
	TH1F * h_phiJet[3];

	TH1F * h_PFMet[3];
	TH1F * h_SumPFMet[3];
	TH1F * h_MetBySumMET[3];

	TH2F * h_pfMetgjMass[3];
	TH2F * h_PtPhotJet[3];
	TH2F * h_etaPhotJet[3];
	TH1F * h_DR_PhotonJet[3];
	TH1F * h_cosThetaStar[3];
	TH1F * h_dEta[3];
	TH1F * h_dphi[3];

	TH1F * h_Photon_SigmaIetaIeta[3];
	TH1F * h_HoE[3];
	TH1F * h_CorrPFiso_Charged[3];
	TH1F * h_CorrPFiso_Neutral[3];
	TH1F * h_CorrPFiso_Photon[3];
	TH1F * h_PFiso_Electronveto[3];

	TH1F * h_jet_NEF[3];
	TH1F * h_jet_NHF[3];
	TH1F * h_jet_CEF[3];
	TH1F * h_jet_CHF[3];
	TH1F * h_jet_NConstituents[3];
	TH1F * h_jet_ChargeMultiplicity[3];

	TH1F * h_nIsoPhoton[4];
	TH1F * h_nPhoton[4];
	TH2F * h_nIsoPhotonPt[4];
	TH2F * h_nPhotonPt[4];
	TH1F * h_nJet[4];
	TH2F * h_nJetPt[4];

	TH1F * h_Vertices[4];
	TH1F * h_trueInteractions[2];


	TH1F * h_PC;
	TH1F * h_JC;

	TH1F * h_CutFlow;
	TH1F * h_CutExpFlow;

	TH1F * h_dataPU;
	TH1F * h_mcPU;
	TH1F * h_weights; 

	// Fixed size dimensions of array or collections stored in the TTree if any.

	// Declaration of leaf types
	Int_t           run;
	Long64_t        event;
	Int_t           lumis;
	Bool_t          isData;
	Int_t           nVtx;
	Int_t           nTrksPV;
	Bool_t          isPVGood;
	Bool_t          hasGoodVtx;
	Float_t         vtx;
	Float_t         vty;
	Float_t         vtz;
	Float_t         rho;
	Float_t         rhoCentral;
	ULong64_t       HLTEleMuX;
	ULong64_t       HLTPho;
	ULong64_t       HLTJet;
	ULong64_t       HLTEleMuXIsPrescaled;
	ULong64_t       HLTPhoIsPrescaled;
	ULong64_t       HLTJetIsPrescaled;
	vector<float>   *pdf;
	Float_t         pthat;
	Float_t         processID;
	Float_t         genWeight;
	Float_t         genHT;
	Int_t           nPUInfo;
	vector<int>     *nPU;
	vector<int>     *puBX;
	vector<float>   *puTrue;
	Int_t           nMC;
	vector<int>     *mcPID;
	vector<float>   *mcVtx;
	vector<float>   *mcVty;
	vector<float>   *mcVtz;
	vector<float>   *mcPt;
	vector<float>   *mcMass;
	vector<float>   *mcEta;
	vector<float>   *mcPhi;
	vector<float>   *mcE;
	vector<float>   *mcEt;
	vector<int>     *mcGMomPID;
	vector<int>     *mcMomPID;
	vector<float>   *mcMomPt;
	vector<float>   *mcMomMass;
	vector<float>   *mcMomEta;
	vector<float>   *mcMomPhi;
	vector<int>     *mcIndex;
	vector<unsigned short> *mcStatusFlag;
	vector<int>     *mcParentage;
	vector<int>     *mcStatus;
	vector<float>   *mcCalIsoDR03;
	vector<float>   *mcTrkIsoDR03;
	vector<float>   *mcCalIsoDR04;
	vector<float>   *mcTrkIsoDR04;
	Int_t           metFilters;
	Float_t         genMET;
	Float_t         genMETPhi;
	Float_t         pfMET;
	Float_t         pfMETPhi;
	Float_t         pfMETsumEt;
	Float_t         pfMETmEtSig;
	Float_t         pfMETSig;
	Float_t         pfMET_T1JERUp;
	Float_t         pfMET_T1JERDo;
	Float_t         pfMET_T1JESUp;
	Float_t         pfMET_T1JESDo;
	Float_t         pfMET_T1MESUp;
	Float_t         pfMET_T1MESDo;
	Float_t         pfMET_T1EESUp;
	Float_t         pfMET_T1EESDo;
	Float_t         pfMET_T1PESUp;
	Float_t         pfMET_T1PESDo;
	Float_t         pfMET_T1TESUp;
	Float_t         pfMET_T1TESDo;
	Float_t         pfMET_T1UESUp;
	Float_t         pfMET_T1UESDo;
	Int_t           nPho;
	vector<float>   *phoE;
	vector<float>   *phoEt;
	vector<float>   *phoEta;
	vector<float>   *phoPhi;
	vector<float>   *phoSCE;
	vector<float>   *phoSCRawE;
	vector<float>   *phoESEn;
	vector<float>   *phoESEnP1;
	vector<float>   *phoESEnP2;
	vector<float>   *phoSCEta;
	vector<float>   *phoSCPhi;
	vector<float>   *phoSCEtaWidth;
	vector<float>   *phoSCPhiWidth;
	vector<float>   *phoSCBrem;
	vector<int>     *phohasPixelSeed;
	vector<int>     *phoEleVeto;
	vector<float>   *phoR9;
	vector<float>   *phoHoverE;
	vector<float>   *phoSigmaIEtaIEta;
	vector<float>   *phoSigmaIEtaIPhi;
	vector<float>   *phoSigmaIPhiIPhi;
	vector<float>   *phoE1x3;
	vector<float>   *phoE1x5;
	vector<float>   *phoE2x2;
	vector<float>   *phoE2x5Max;
	vector<float>   *phoE5x5;
	vector<float>   *phoESEffSigmaRR;
	vector<float>   *phoSigmaIEtaIEtaFull5x5;
	vector<float>   *phoSigmaIEtaIPhiFull5x5;
	vector<float>   *phoSigmaIPhiIPhiFull5x5;
	vector<float>   *phoE1x3Full5x5;
	vector<float>   *phoE1x5Full5x5;
	vector<float>   *phoE2x2Full5x5;
	vector<float>   *phoE2x5MaxFull5x5;
	vector<float>   *phoE5x5Full5x5;
	vector<float>   *phoR9Full5x5;
	vector<float>   *phoSeedBCE;
	vector<float>   *phoSeedBCEta;
	vector<float>   *phoPFChIso;
	vector<float>   *phoPFPhoIso;
	vector<float>   *phoPFNeuIso;
	vector<float>   *phoPFChWorstIso;
	vector<float>   *phoPFChIsoFrix1;
	vector<float>   *phoPFChIsoFrix2;
	vector<float>   *phoPFChIsoFrix3;
	vector<float>   *phoPFChIsoFrix4;
	vector<float>   *phoPFChIsoFrix5;
	vector<float>   *phoPFChIsoFrix6;
	vector<float>   *phoPFChIsoFrix7;
	vector<float>   *phoPFChIsoFrix8;
	vector<float>   *phoPFPhoIsoFrix1;
	vector<float>   *phoPFPhoIsoFrix2;
	vector<float>   *phoPFPhoIsoFrix3;
	vector<float>   *phoPFPhoIsoFrix4;
	vector<float>   *phoPFPhoIsoFrix5;
	vector<float>   *phoPFPhoIsoFrix6;
	vector<float>   *phoPFPhoIsoFrix7;
	vector<float>   *phoPFPhoIsoFrix8;
	vector<float>   *phoPFNeuIsoFrix1;
	vector<float>   *phoPFNeuIsoFrix2;
	vector<float>   *phoPFNeuIsoFrix3;
	vector<float>   *phoPFNeuIsoFrix4;
	vector<float>   *phoPFNeuIsoFrix5;
	vector<float>   *phoPFNeuIsoFrix6;
	vector<float>   *phoPFNeuIsoFrix7;
	vector<float>   *phoPFNeuIsoFrix8;
	vector<float>   *phoEcalRecHitSumEtConeDR03;
	vector<float>   *phohcalDepth1TowerSumEtConeDR03;
	vector<float>   *phohcalDepth2TowerSumEtConeDR03;
	vector<float>   *phohcalTowerSumEtConeDR03;
	vector<float>   *photrkSumPtHollowConeDR03;
	vector<float>   *photrkSumPtSolidConeDR03;
	vector<float>   *phoIDMVA;
	vector<int>     *phoFiredSingleTrgs;
	vector<int>     *phoFiredDoubleTrgs;
	vector<int>     *phoIEta;
	vector<int>     *phoIPhi;
	vector<float>   *phomaxXtalenergyFull5x5;
	vector<float>   *phoseedTimeFull5x5;
	vector<float>   *phomaxXtalenergy;
	vector<float>   *phoseedTime;
	vector<unsigned short> *phoxtalBits;
	vector<unsigned short> *phoIDbit;
	Int_t           nEle;
	vector<int>     *eleCharge;
	vector<int>     *eleChargeConsistent;
	vector<float>   *eleEn;
	vector<float>   *eleSCEn;
	vector<float>   *eleESEn;
	vector<float>   *eleESEnP1;
	vector<float>   *eleESEnP2;
	vector<float>   *eleD0;
	vector<float>   *eleDz;
	vector<float>   *elePt;
	vector<float>   *eleEta;
	vector<float>   *elePhi;
	vector<float>   *eleR9;
	vector<float>   *eleSCEta;
	vector<float>   *eleSCPhi;
	vector<float>   *eleSCRawEn;
	vector<float>   *eleSCEtaWidth;
	vector<float>   *eleSCPhiWidth;
	vector<float>   *eleHoverE;
	vector<float>   *eleEoverP;
	vector<float>   *eleEoverPout;
	vector<float>   *eleEoverPInv;
	vector<float>   *eleBrem;
	vector<float>   *eledEtaAtVtx;
	vector<float>   *eledPhiAtVtx;
	vector<float>   *eledEtaAtCalo;
	vector<float>   *eleSigmaIEtaIEta;
	vector<float>   *eleSigmaIEtaIPhi;
	vector<float>   *eleSigmaIPhiIPhi;
	vector<float>   *eleSigmaIEtaIEtaFull5x5;
	vector<float>   *eleSigmaIPhiIPhiFull5x5;
	vector<int>     *eleConvVeto;
	vector<int>     *eleMissHits;
	vector<float>   *eleESEffSigmaRR;
	vector<float>   *elePFChIso;
	vector<float>   *elePFPhoIso;
	vector<float>   *elePFNeuIso;
	vector<float>   *elePFPUIso;
	vector<float>   *elePFClusEcalIso;
	vector<float>   *elePFClusHcalIso;
	vector<float>   *elePFMiniIso;
	vector<float>   *eleIDMVANonTrg;
	vector<float>   *eleIDMVATrg;
	vector<float>   *eledEtaseedAtVtx;
	vector<float>   *eleE1x5;
	vector<float>   *eleE2x5;
	vector<float>   *eleE5x5;
	vector<float>   *eleE1x5Full5x5;
	vector<float>   *eleE2x5Full5x5;
	vector<float>   *eleE5x5Full5x5;
	vector<float>   *eleR9Full5x5;
	vector<int>     *eleEcalDrivenSeed;
	vector<float>   *eleDr03EcalRecHitSumEt;
	vector<float>   *eleDr03HcalDepth1TowerSumEt;
	vector<float>   *eleDr03HcalDepth2TowerSumEt;
	vector<float>   *eleDr03HcalTowerSumEt;
	vector<float>   *eleDr03TkSumPt;
	vector<float>   *elecaloEnergy;
	vector<float>   *eleTrkdxy;
	vector<float>   *eleKFHits;
	vector<float>   *eleKFChi2;
	vector<vector<float> > *eleGSFPt;
	vector<vector<float> > *eleGSFEta;
	vector<vector<float> > *eleGSFPhi;
	vector<vector<float> > *eleGSFCharge;
	vector<vector<int> > *eleGSFHits;
	vector<vector<int> > *eleGSFMissHits;
	vector<vector<int> > *eleGSFNHitsMax;
	vector<vector<float> > *eleGSFVtxProb;
	vector<vector<float> > *eleGSFlxyPV;
	vector<vector<float> > *eleGSFlxyBS;
	vector<vector<float> > *eleBCEn;
	vector<vector<float> > *eleBCEta;
	vector<vector<float> > *eleBCPhi;
	vector<vector<float> > *eleBCS25;
	vector<vector<float> > *eleBCS15;
	vector<vector<float> > *eleBCSieie;
	vector<vector<float> > *eleBCSieip;
	vector<vector<float> > *eleBCSipip;
	vector<int>     *eleFiredTrgs;
	vector<unsigned short> *eleIDbit;
	vector<float>   *eleESEnP1Raw;
	vector<float>   *eleESEnP2Raw;
	vector<vector<float> > *eleESEnEta;
	vector<vector<float> > *eleESEnPhi;
	vector<vector<int> > *eleESEnZ;
	vector<vector<int> > *eleESEnP;
	vector<vector<int> > *eleESEnX;
	vector<vector<int> > *eleESEnY;
	vector<vector<int> > *eleESEnS;
	vector<vector<float> > *eleESEnE;
	Int_t           nGSFTrk;
	vector<float>   *gsfPt;
	vector<float>   *gsfEta;
	vector<float>   *gsfPhi;
	Int_t           npfHF;
	vector<float>   *pfHFEn;
	vector<float>   *pfHFECALEn;
	vector<float>   *pfHFHCALEn;
	vector<float>   *pfHFPt;
	vector<float>   *pfHFEta;
	vector<float>   *pfHFPhi;
	vector<float>   *pfHFIso;
	Int_t           nMu;
	vector<float>   *muPt;
	vector<float>   *muEn;
	vector<float>   *muEta;
	vector<float>   *muPhi;
	vector<int>     *muCharge;
	vector<int>     *muType;
	vector<bool>    *muIsLooseID;
	vector<bool>    *muIsMediumID;
	vector<bool>    *muIsTightID;
	vector<bool>    *muIsSoftID;
	vector<bool>    *muIsHighPtID;
	vector<float>   *muD0;
	vector<float>   *muDz;
	vector<float>   *muChi2NDF;
	vector<float>   *muInnerD0;
	vector<float>   *muInnerDz;
	vector<int>     *muTrkLayers;
	vector<int>     *muPixelLayers;
	vector<int>     *muPixelHits;
	vector<int>     *muMuonHits;
	vector<int>     *muStations;
	vector<int>     *muMatches;
	vector<int>     *muTrkQuality;
	vector<float>   *muIsoTrk;
	vector<float>   *muPFChIso;
	vector<float>   *muPFPhoIso;
	vector<float>   *muPFNeuIso;
	vector<float>   *muPFPUIso;
	vector<float>   *muPFMiniIso;
	vector<int>     *muFiredTrgs;
	vector<float>   *muInnervalidFraction;
	vector<float>   *musegmentCompatibility;
	vector<float>   *muchi2LocalPosition;
	vector<float>   *mutrkKink;
	vector<float>   *muBestTrkPtError;
	vector<float>   *muBestTrkPt;
	vector<bool>    *muIsPFMuon;
	vector<bool>    *muIsGlobalMuon;
	vector<bool>    *muIsTrackerMuon;
	Int_t           nJet;
	vector<float>   *jetPt;
	vector<float>   *jetEn;
	vector<float>   *jetEta;
	vector<float>   *jetPhi;
	vector<float>   *jetRawPt;
	vector<float>   *jetRawEn;
	vector<float>   *jetMt;
	vector<float>   *jetArea;
	vector<float>   *jetLeadTrackPt;
	vector<float>   *jetLeadTrackEta;
	vector<float>   *jetLeadTrackPhi;
	vector<int>     *jetLepTrackPID;
	vector<float>   *jetLepTrackPt;
	vector<float>   *jetLepTrackEta;
	vector<float>   *jetLepTrackPhi;
	vector<float>   *jetpfCombinedInclusiveSecondaryVertexV2BJetTags;
	vector<float>   *jetJetProbabilityBJetTags;
	vector<float>   *jetpfCombinedMVABJetTags;
	vector<int>     *jetPartonID;
	vector<int>     *jetGenJetIndex;
	vector<float>   *jetGenJetEn;
	vector<float>   *jetGenJetPt;
	vector<float>   *jetGenJetEta;
	vector<float>   *jetGenJetPhi;
	vector<int>     *jetGenPartonID;
	vector<float>   *jetGenEn;
	vector<float>   *jetGenPt;
	vector<float>   *jetGenEta;
	vector<float>   *jetGenPhi;
	vector<int>     *jetGenPartonMomID;
	vector<bool>    *jetPFLooseId;
	vector<float>   *jetPUidFullDiscriminant;
	vector<float>   *jetJECUnc;
	vector<int>     *jetFiredTrgs;
	vector<float>   *jetCHF;
	vector<float>   *jetNHF;
	vector<float>   *jetCEF;
	vector<float>   *jetNEF;
	vector<int>     *jetNCH;
	vector<float>   *jetMUF;
	vector<float>   *jetMass;
	vector<float>   *jetPx;
	vector<float>   *jetPy;
	vector<float>   *jetPz;
	vector<float>   *jetVtxPt;
	vector<float>   *jetVtxMass;
	vector<float>   *jetVtxNtrks;
	vector<float>   *jetVtx3DVal;
	vector<float>   *jetVtx3DSig;
	vector<float>   *jetHFHAE;
	vector<float>   *jetHFEME;
	vector<int>     *jetNConstituents;
	Int_t           nAK8Jet;
	vector<float>   *AK8JetPt;
	vector<float>   *AK8JetEn;
	vector<float>   *AK8JetRawPt;
	vector<float>   *AK8JetRawEn;
	vector<float>   *AK8JetEta;
	vector<float>   *AK8JetPhi;
	vector<float>   *AK8JetMass;
	vector<float>   *AK8Jet_tau1;
	vector<float>   *AK8Jet_tau2;
	vector<float>   *AK8Jet_tau3;
	vector<float>   *AK8JetCHF;
	vector<float>   *AK8JetNHF;
	vector<float>   *AK8JetCEF;
	vector<float>   *AK8JetNEF;
	vector<int>     *AK8JetNCH;
	vector<int>     *AK8Jetnconstituents;
	vector<float>   *AK8JetMUF;
	vector<float>   *AK8JetPx;
	vector<float>   *AK8JetPy;
	vector<float>   *AK8JetPz;
	vector<bool>    *AK8JetPFLooseId;
	vector<bool>    *AK8JetPFTightLepVetoId;
	vector<float>   *AK8CHSSoftDropJetMass;
	vector<float>   *AK8CHSSoftDropJetMassCorr;
	vector<float>   *AK8CHSPrunedJetMass;
	vector<float>   *AK8CHSPrunedJetMassCorr;
	vector<float>   *AK8JetpfBoostedDSVBTag;
	vector<float>   *AK8JetJECUnc;
	vector<float>   *AK8JetL2L3corr;
	vector<int>     *AK8JetPartonID;
	vector<int>     *AK8JetGenJetIndex;
	vector<float>   *AK8JetGenJetEn;
	vector<float>   *AK8JetGenJetPt;
	vector<float>   *AK8JetGenJetEta;
	vector<float>   *AK8JetGenJetPhi;
	vector<int>     *AK8JetGenPartonID;
	vector<float>   *AK8JetGenEn;
	vector<float>   *AK8JetGenPt;
	vector<float>   *AK8JetGenEta;
	vector<float>   *AK8JetGenPhi;
	vector<int>     *AK8JetGenPartonMomID;
	vector<int>     *nAK8softdropSubjet;
	vector<vector<float> > *AK8softdropSubjetPt;
	vector<vector<float> > *AK8softdropSubjetEta;
	vector<vector<float> > *AK8softdropSubjetPhi;
	vector<vector<float> > *AK8softdropSubjetMass;
	vector<vector<float> > *AK8softdropSubjetE;
	vector<vector<int> > *AK8softdropSubjetCharge;
	vector<vector<int> > *AK8softdropSubjetFlavour;
	vector<vector<float> > *AK8softdropSubjetCSV;

	// List of branches
	TBranch        *b_run;   //!
	TBranch        *b_event;   //!
	TBranch        *b_lumis;   //!
	TBranch        *b_isData;   //!
	TBranch        *b_nVtx;   //!
	TBranch        *b_nTrksPV;   //!
	TBranch        *b_isPVGood;   //!
	TBranch        *b_hasGoodVtx;   //!
	TBranch        *b_vtx;   //!
	TBranch        *b_vty;   //!
	TBranch        *b_vtz;   //!
	TBranch        *b_rho;   //!
	TBranch        *b_rhoCentral;   //!
	TBranch        *b_HLTEleMuX;   //!
	TBranch        *b_HLTPho;   //!
	TBranch        *b_HLTJet;   //!
	TBranch        *b_HLTEleMuXIsPrescaled;   //!
	TBranch        *b_HLTPhoIsPrescaled;   //!
	TBranch        *b_HLTJetIsPrescaled;   //!
	TBranch        *b_pdf;   //!
	TBranch        *b_pthat;   //!
	TBranch        *b_processID;   //!
	TBranch        *b_genWeight;   //!
	TBranch        *b_genHT;   //!
	TBranch        *b_nPUInfo;   //!
	TBranch        *b_nPU;   //!
	TBranch        *b_puBX;   //!
	TBranch        *b_puTrue;   //!
	TBranch        *b_nMC;   //!
	TBranch        *b_mcPID;   //!
	TBranch        *b_mcVtx;   //!
	TBranch        *b_mcVty;   //!
	TBranch        *b_mcVtz;   //!
	TBranch        *b_mcPt;   //!
	TBranch        *b_mcMass;   //!
	TBranch        *b_mcEta;   //!
	TBranch        *b_mcPhi;   //!
	TBranch        *b_mcE;   //!
	TBranch        *b_mcEt;   //!
	TBranch        *b_mcGMomPID;   //!
	TBranch        *b_mcMomPID;   //!
	TBranch        *b_mcMomPt;   //!
	TBranch        *b_mcMomMass;   //!
	TBranch        *b_mcMomEta;   //!
	TBranch        *b_mcMomPhi;   //!
	TBranch        *b_mcIndex;   //!
	TBranch        *b_mcStatusFlag;   //!
	TBranch        *b_mcParentage;   //!
	TBranch        *b_mcStatus;   //!
	TBranch        *b_mcCalIsoDR03;   //!
	TBranch        *b_mcTrkIsoDR03;   //!
	TBranch        *b_mcCalIsoDR04;   //!
	TBranch        *b_mcTrkIsoDR04;   //!
	TBranch        *b_metFilters;   //!
	TBranch        *b_genMET;   //!
	TBranch        *b_genMETPhi;   //!
	TBranch        *b_pfMET;   //!
	TBranch        *b_pfMETPhi;   //!
	TBranch        *b_pfMETsumEt;   //!
	TBranch        *b_pfMETmEtSig;   //!
	TBranch        *b_pfMETSig;   //!
	TBranch        *b_pfMET_T1JERUp;   //!
	TBranch        *b_pfMET_T1JERDo;   //!
	TBranch        *b_pfMET_T1JESUp;   //!
	TBranch        *b_pfMET_T1JESDo;   //!
	TBranch        *b_pfMET_T1MESUp;   //!
	TBranch        *b_pfMET_T1MESDo;   //!
	TBranch        *b_pfMET_T1EESUp;   //!
	TBranch        *b_pfMET_T1EESDo;   //!
	TBranch        *b_pfMET_T1PESUp;   //!
	TBranch        *b_pfMET_T1PESDo;   //!
	TBranch        *b_pfMET_T1TESUp;   //!
	TBranch        *b_pfMET_T1TESDo;   //!
	TBranch        *b_pfMET_T1UESUp;   //!
	TBranch        *b_pfMET_T1UESDo;   //!
	TBranch        *b_nPho;   //!
	TBranch        *b_phoE;   //!
	TBranch        *b_phoEt;   //!
	TBranch        *b_phoEta;   //!
	TBranch        *b_phoPhi;   //!
	TBranch        *b_phoSCE;   //!
	TBranch        *b_phoSCRawE;   //!
	TBranch        *b_phoESEn;   //!
	TBranch        *b_phoESEnP1;   //!
	TBranch        *b_phoESEnP2;   //!
	TBranch        *b_phoSCEta;   //!
	TBranch        *b_phoSCPhi;   //!
	TBranch        *b_phoSCEtaWidth;   //!
	TBranch        *b_phoSCPhiWidth;   //!
	TBranch        *b_phoSCBrem;   //!
	TBranch        *b_phohasPixelSeed;   //!
	TBranch        *b_phoEleVeto;   //!
	TBranch        *b_phoR9;   //!
	TBranch        *b_phoHoverE;   //!
	TBranch        *b_phoSigmaIEtaIEta;   //!
	TBranch        *b_phoSigmaIEtaIPhi;   //!
	TBranch        *b_phoSigmaIPhiIPhi;   //!
	TBranch        *b_phoE1x3;   //!
	TBranch        *b_phoE1x5;   //!
	TBranch        *b_phoE2x2;   //!
	TBranch        *b_phoE2x5Max;   //!
	TBranch        *b_phoE5x5;   //!
	TBranch        *b_phoESEffSigmaRR;   //!
	TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
	TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
	TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
	TBranch        *b_phoE1x3Full5x5;   //!
	TBranch        *b_phoE1x5Full5x5;   //!
	TBranch        *b_phoE2x2Full5x5;   //!
	TBranch        *b_phoE2x5MaxFull5x5;   //!
	TBranch        *b_phoE5x5Full5x5;   //!
	TBranch        *b_phoR9Full5x5;   //!
	TBranch        *b_phoSeedBCE;   //!
	TBranch        *b_phoSeedBCEta;   //!
	TBranch        *b_phoPFChIso;   //!
	TBranch        *b_phoPFPhoIso;   //!
	TBranch        *b_phoPFNeuIso;   //!
	TBranch        *b_phoPFChWorstIso;   //!
	TBranch        *b_phoPFChIsoFrix1;   //!
	TBranch        *b_phoPFChIsoFrix2;   //!
	TBranch        *b_phoPFChIsoFrix3;   //!
	TBranch        *b_phoPFChIsoFrix4;   //!
	TBranch        *b_phoPFChIsoFrix5;   //!
	TBranch        *b_phoPFChIsoFrix6;   //!
	TBranch        *b_phoPFChIsoFrix7;   //!
	TBranch        *b_phoPFChIsoFrix8;   //!
	TBranch        *b_phoPFPhoIsoFrix1;   //!
	TBranch        *b_phoPFPhoIsoFrix2;   //!
	TBranch        *b_phoPFPhoIsoFrix3;   //!
	TBranch        *b_phoPFPhoIsoFrix4;   //!
	TBranch        *b_phoPFPhoIsoFrix5;   //!
	TBranch        *b_phoPFPhoIsoFrix6;   //!
	TBranch        *b_phoPFPhoIsoFrix7;   //!
	TBranch        *b_phoPFPhoIsoFrix8;   //!
	TBranch        *b_phoPFNeuIsoFrix1;   //!
	TBranch        *b_phoPFNeuIsoFrix2;   //!
	TBranch        *b_phoPFNeuIsoFrix3;   //!
	TBranch        *b_phoPFNeuIsoFrix4;   //!
	TBranch        *b_phoPFNeuIsoFrix5;   //!
	TBranch        *b_phoPFNeuIsoFrix6;   //!
	TBranch        *b_phoPFNeuIsoFrix7;   //!
	TBranch        *b_phoPFNeuIsoFrix8;   //!
	TBranch        *b_phoEcalRecHitSumEtConeDR03;   //!
	TBranch        *b_phohcalDepth1TowerSumEtConeDR03;   //!
	TBranch        *b_phohcalDepth2TowerSumEtConeDR03;   //!
	TBranch        *b_phohcalTowerSumEtConeDR03;   //!
	TBranch        *b_photrkSumPtHollowConeDR03;   //!
	TBranch        *b_photrkSumPtSolidConeDR03;   //!
	TBranch        *b_phoIDMVA;   //!
	TBranch        *b_phoFiredSingleTrgs;   //!
	TBranch        *b_phoFiredDoubleTrgs;   //!
	TBranch        *b_phoIEta;   //!
	TBranch        *b_phoIPhi;   //!
	TBranch        *b_phomaxXtalenergyFull5x5;   //!
	TBranch        *b_phoseedTimeFull5x5;   //!
	TBranch        *b_phomaxXtalenergy;   //!
	TBranch        *b_phoseedTime;   //!
	TBranch        *b_phoxtalBits;   //!
	TBranch        *b_phoIDbit;   //!
	TBranch        *b_nEle;   //!
	TBranch        *b_eleCharge;   //!
	TBranch        *b_eleChargeConsistent;   //!
	TBranch        *b_eleEn;   //!
	TBranch        *b_eleSCEn;   //!
	TBranch        *b_eleESEn;   //!
	TBranch        *b_eleESEnP1;   //!
	TBranch        *b_eleESEnP2;   //!
	TBranch        *b_eleD0;   //!
	TBranch        *b_eleDz;   //!
	TBranch        *b_elePt;   //!
	TBranch        *b_eleEta;   //!
	TBranch        *b_elePhi;   //!
	TBranch        *b_eleR9;   //!
	TBranch        *b_eleSCEta;   //!
	TBranch        *b_eleSCPhi;   //!
	TBranch        *b_eleSCRawEn;   //!
	TBranch        *b_eleSCEtaWidth;   //!
	TBranch        *b_eleSCPhiWidth;   //!
	TBranch        *b_eleHoverE;   //!
	TBranch        *b_eleEoverP;   //!
	TBranch        *b_eleEoverPout;   //!
	TBranch        *b_eleEoverPInv;   //!
	TBranch        *b_eleBrem;   //!
	TBranch        *b_eledEtaAtVtx;   //!
	TBranch        *b_eledPhiAtVtx;   //!
	TBranch        *b_eledEtaAtCalo;   //!
	TBranch        *b_eleSigmaIEtaIEta;   //!
	TBranch        *b_eleSigmaIEtaIPhi;   //!
	TBranch        *b_eleSigmaIPhiIPhi;   //!
	TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
	TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
	TBranch        *b_eleConvVeto;   //!
	TBranch        *b_eleMissHits;   //!
	TBranch        *b_eleESEffSigmaRR;   //!
	TBranch        *b_elePFChIso;   //!
	TBranch        *b_elePFPhoIso;   //!
	TBranch        *b_elePFNeuIso;   //!
	TBranch        *b_elePFPUIso;   //!
	TBranch        *b_elePFClusEcalIso;   //!
	TBranch        *b_elePFClusHcalIso;   //!
	TBranch        *b_elePFMiniIso;   //!
	TBranch        *b_eleIDMVANonTrg;   //!
	TBranch        *b_eleIDMVATrg;   //!
	TBranch        *b_eledEtaseedAtVtx;   //!
	TBranch        *b_eleE1x5;   //!
	TBranch        *b_eleE2x5;   //!
	TBranch        *b_eleE5x5;   //!
	TBranch        *b_eleE1x5Full5x5;   //!
	TBranch        *b_eleE2x5Full5x5;   //!
	TBranch        *b_eleE5x5Full5x5;   //!
	TBranch        *b_eleR9Full5x5;   //!
	TBranch        *b_eleEcalDrivenSeed;   //!
	TBranch        *b_eleDr03EcalRecHitSumEt;   //!
	TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
	TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
	TBranch        *b_eleDr03HcalTowerSumEt;   //!
	TBranch        *b_eleDr03TkSumPt;   //!
	TBranch        *b_elecaloEnergy;   //!
	TBranch        *b_eleTrkdxy;   //!
	TBranch        *b_eleKFHits;   //!
	TBranch        *b_eleKFChi2;   //!
	TBranch        *b_eleGSFPt;   //!
	TBranch        *b_eleGSFEta;   //!
	TBranch        *b_eleGSFPhi;   //!
	TBranch        *b_eleGSFCharge;   //!
	TBranch        *b_eleGSFHits;   //!
	TBranch        *b_eleGSFMissHits;   //!
	TBranch        *b_eleGSFNHitsMax;   //!
	TBranch        *b_eleGSFVtxProb;   //!
	TBranch        *b_eleGSFlxyPV;   //!
	TBranch        *b_eleGSFlxyBS;   //!
	TBranch        *b_eleBCEn;   //!
	TBranch        *b_eleBCEta;   //!
	TBranch        *b_eleBCPhi;   //!
	TBranch        *b_eleBCS25;   //!
	TBranch        *b_eleBCS15;   //!
	TBranch        *b_eleBCSieie;   //!
	TBranch        *b_eleBCSieip;   //!
	TBranch        *b_eleBCSipip;   //!
	TBranch        *b_eleFiredTrgs;   //!
	TBranch        *b_eleIDbit;   //!
	TBranch        *b_eleESEnP1Raw;   //!
	TBranch        *b_eleESEnP2Raw;   //!
	TBranch        *b_eleESEnEta;   //!
	TBranch        *b_eleESEnPhi;   //!
	TBranch        *b_eleESEnZ;   //!
	TBranch        *b_eleESEnP;   //!
	TBranch        *b_eleESEnX;   //!
	TBranch        *b_eleESEnY;   //!
	TBranch        *b_eleESEnS;   //!
	TBranch        *b_eleESEnE;   //!
	TBranch        *b_nGSFTrk;   //!
	TBranch        *b_gsfPt;   //!
	TBranch        *b_gsfEta;   //!
	TBranch        *b_gsfPhi;   //!
	TBranch        *b_npfHF;   //!
	TBranch        *b_pfHFEn;   //!
	TBranch        *b_pfHFECALEn;   //!
	TBranch        *b_pfHFHCALEn;   //!
	TBranch        *b_pfHFPt;   //!
	TBranch        *b_pfHFEta;   //!
	TBranch        *b_pfHFPhi;   //!
	TBranch        *b_pfHFIso;   //!
	TBranch        *b_nMu;   //!
	TBranch        *b_muPt;   //!
	TBranch        *b_muEn;   //!
	TBranch        *b_muEta;   //!
	TBranch        *b_muPhi;   //!
	TBranch        *b_muCharge;   //!
	TBranch        *b_muType;   //!
	TBranch        *b_muIsLooseID;   //!
	TBranch        *b_muIsMediumID;   //!
	TBranch        *b_muIsTightID;   //!
	TBranch        *b_muIsSoftID;   //!
	TBranch        *b_muIsHighPtID;   //!
	TBranch        *b_muD0;   //!
	TBranch        *b_muDz;   //!
	TBranch        *b_muChi2NDF;   //!
	TBranch        *b_muInnerD0;   //!
	TBranch        *b_muInnerDz;   //!
	TBranch        *b_muTrkLayers;   //!
	TBranch        *b_muPixelLayers;   //!
	TBranch        *b_muPixelHits;   //!
	TBranch        *b_muMuonHits;   //!
	TBranch        *b_muStations;   //!
	TBranch        *b_muMatches;   //!
	TBranch        *b_muTrkQuality;   //!
	TBranch        *b_muIsoTrk;   //!
	TBranch        *b_muPFChIso;   //!
	TBranch        *b_muPFPhoIso;   //!
	TBranch        *b_muPFNeuIso;   //!
	TBranch        *b_muPFPUIso;   //!
	TBranch        *b_muPFMiniIso;   //!
	TBranch        *b_muFiredTrgs;   //!
	TBranch        *b_muInnervalidFraction;   //!
	TBranch        *b_musegmentCompatibility;   //!
	TBranch        *b_muchi2LocalPosition;   //!
	TBranch        *b_mutrkKink;   //!
	TBranch        *b_muBestTrkPtError;   //!
	TBranch        *b_muBestTrkPt;   //!
	TBranch        *b_muIsPFMuon;   //!
	TBranch        *b_muIsGlobalMuon;   //!
	TBranch        *b_muIsTrackerMuon;   //!
	TBranch        *b_nJet;   //!
	TBranch        *b_jetPt;   //!
	TBranch        *b_jetEn;   //!
	TBranch        *b_jetEta;   //!
	TBranch        *b_jetPhi;   //!
	TBranch        *b_jetRawPt;   //!
	TBranch        *b_jetRawEn;   //!
	TBranch        *b_jetMt;   //!
	TBranch        *b_jetArea;   //!
	TBranch        *b_jetLeadTrackPt;   //!
	TBranch        *b_jetLeadTrackEta;   //!
	TBranch        *b_jetLeadTrackPhi;   //!
	TBranch        *b_jetLepTrackPID;   //!
	TBranch        *b_jetLepTrackPt;   //!
	TBranch        *b_jetLepTrackEta;   //!
	TBranch        *b_jetLepTrackPhi;   //!
	TBranch        *b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
	TBranch        *b_jetJetProbabilityBJetTags;   //!
	TBranch        *b_jetpfCombinedMVABJetTags;   //!
	TBranch        *b_jetPartonID;   //!
	TBranch        *b_jetGenJetIndex;   //!
	TBranch        *b_jetGenJetEn;   //!
	TBranch        *b_jetGenJetPt;   //!
	TBranch        *b_jetGenJetEta;   //!
	TBranch        *b_jetGenJetPhi;   //!
	TBranch        *b_jetGenPartonID;   //!
	TBranch        *b_jetGenEn;   //!
	TBranch        *b_jetGenPt;   //!
	TBranch        *b_jetGenEta;   //!
	TBranch        *b_jetGenPhi;   //!
	TBranch        *b_jetGenPartonMomID;   //!
	TBranch        *b_jetPFLooseId;   //!
	TBranch        *b_jetPUidFullDiscriminant;   //!
	TBranch        *b_jetJECUnc;   //!
	TBranch        *b_jetFiredTrgs;   //!
	TBranch        *b_jetCHF;   //!
	TBranch        *b_jetNHF;   //!
	TBranch        *b_jetCEF;   //!
	TBranch        *b_jetNEF;   //!
	TBranch        *b_jetNCH;   //!
	TBranch        *b_jetMUF;   //!
	TBranch        *b_jetMass;   //!
	TBranch        *b_jetPx;   //!
	TBranch        *b_jetPy;   //!
	TBranch        *b_jetPz;   //!
	TBranch        *b_jetVtxPt;   //!
	TBranch        *b_jetVtxMass;   //!
	TBranch        *b_jetVtxNtrks;   //!
	TBranch        *b_jetVtx3DVal;   //!
	TBranch        *b_jetVtx3DSig;   //!
	TBranch        *b_jetHFHAE;   //!
	TBranch        *b_jetHFEME;   //!
	TBranch        *b_jetNConstituents;   //!
	TBranch        *b_nAK8Jet;   //!
	TBranch        *b_AK8JetPt;   //!
	TBranch        *b_AK8JetEn;   //!
	TBranch        *b_AK8JetRawPt;   //!
	TBranch        *b_AK8JetRawEn;   //!
	TBranch        *b_AK8JetEta;   //!
	TBranch        *b_AK8JetPhi;   //!
	TBranch        *b_AK8JetMass;   //!
	TBranch        *b_AK8Jet_tau1;   //!
	TBranch        *b_AK8Jet_tau2;   //!
	TBranch        *b_AK8Jet_tau3;   //!
	TBranch        *b_AK8JetCHF;   //!
	TBranch        *b_AK8JetNHF;   //!
	TBranch        *b_AK8JetCEF;   //!
	TBranch        *b_AK8JetNEF;   //!
	TBranch        *b_AK8JetNCH;   //!
	TBranch        *b_AK8Jetnconstituents;   //!
	TBranch        *b_AK8JetMUF;   //!
	TBranch        *b_AK8JetPx;   //!
	TBranch        *b_AK8JetPy;   //!
	TBranch        *b_AK8JetPz;   //!
	TBranch        *b_AK8JetPFLooseId;   //!
	TBranch        *b_AK8JetPFTightLepVetoId;   //!
	TBranch        *b_AK8CHSSoftDropJetMass;   //!
	TBranch        *b_AK8CHSSoftDropJetMassCorr;   //!
	TBranch        *b_AK8CHSPrunedJetMass;   //!
	TBranch        *b_AK8CHSPrunedJetMassCorr;   //!
	TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
	TBranch        *b_AK8JetJECUnc;   //!
	TBranch        *b_AK8JetL2L3corr;   //!
	TBranch        *b_AK8JetPartonID;   //!
	TBranch        *b_AK8JetGenJetIndex;   //!
	TBranch        *b_AK8JetGenJetEn;   //!
	TBranch        *b_AK8JetGenJetPt;   //!
	TBranch        *b_AK8JetGenJetEta;   //!
	TBranch        *b_AK8JetGenJetPhi;   //!
	TBranch        *b_AK8JetGenPartonID;   //!
	TBranch        *b_AK8JetGenEn;   //!
	TBranch        *b_AK8JetGenPt;   //!
	TBranch        *b_AK8JetGenEta;   //!
	TBranch        *b_AK8JetGenPhi;   //!
	TBranch        *b_AK8JetGenPartonMomID;   //!
	TBranch        *b_nAK8softdropSubjet;   //!
	TBranch        *b_AK8softdropSubjetPt;   //!
	TBranch        *b_AK8softdropSubjetEta;   //!
	TBranch        *b_AK8softdropSubjetPhi;   //!
	TBranch        *b_AK8softdropSubjetMass;   //!
	TBranch        *b_AK8softdropSubjetE;   //!
	TBranch        *b_AK8softdropSubjetCharge;   //!
	TBranch        *b_AK8softdropSubjetFlavour;   //!
	TBranch        *b_AK8softdropSubjetCSV;   //!

	//PhoJet_Analyzer_MC(TTree *tree=0);
	PhoJet_Analyzer_MC();
	virtual ~PhoJet_Analyzer_MC();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	// User Defined functions :
	void     BookHistos();
	void     WriteHistos();
	Bool_t   NoSpike(Int_t ipho);
	Bool_t   PrimaryVertex(Int_t &goodVertex);
	Double_t getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
	Double_t getDEta(Double_t eta1, Double_t eta2);
	Double_t getCosThetaStar(Double_t eta1, Double_t eta2);
	Double_t getDPhi(Double_t phi1, Double_t phi2);
	Double_t getMass(Int_t p, Int_t j);
	Int_t    getPhotonCand(TString phoWP);
	Bool_t   CutBasedPFPhotonID(Int_t ipho, TString phoWP);
	Double_t EAcharged(Double_t eta);
	Double_t EAneutral(Double_t eta);
	Double_t EAphoton(Double_t eta);
	Bool_t   jetID(Int_t ijet, TString jetWP);
	Bool_t   AK8jetID(Int_t ijet, TString jetWP);
	void     LumiReWeighting();
	Double_t puweight(Float_t npv);

};

#endif

#ifdef PhoJet_Analyzer_MC_cxx
PhoJet_Analyzer_MC::PhoJet_Analyzer_MC()
{

    TChain *chain = new TChain("ggNtuplizer/EventTree");

    //add input files here
    ifstream datafile;
    datafile.open("${datafile}", ifstream::in);
    char datafilename[500];

    for(Int_t ii = 1; ii <= ${ef} && ii <= ${Tot}; ++ii){
	datafile >> datafilename;
	string fname(datafilename);
	string main_path       = "${InputFilesPath}";

	if(ii >= ${sf} && ii <= ${Tot}){
	    chain->Add((main_path+fname).c_str());
	    cout<<((main_path+fname).c_str())<<endl;
	    cout<<chain->GetEntries()<<endl;
	}
    }

    //    chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/Signal_Qstar/QstarToGJ_M1000_f1p0/AOD_QstarToGJ_M1000_f1p0_1.root");
    Init(chain);

}

PhoJet_Analyzer_MC::~PhoJet_Analyzer_MC()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();

    f1->cd();
    f1->Write();
    h_dataPU->Write();
    h_mcPU->Write();
   
    f1->Close();
}

Int_t PhoJet_Analyzer_MC::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t PhoJet_Analyzer_MC::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
	fCurrent = fChain->GetTreeNumber();
	Notify();
    }
    return centry;
}

void PhoJet_Analyzer_MC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdf = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoSigmaIEtaIEta = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoE1x3 = 0;
   phoE1x5 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE5x5 = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE1x3Full5x5 = 0;
   phoE1x5Full5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE2x5MaxFull5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoSeedBCE = 0;
   phoSeedBCEta = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoPFChIsoFrix1 = 0;
   phoPFChIsoFrix2 = 0;
   phoPFChIsoFrix3 = 0;
   phoPFChIsoFrix4 = 0;
   phoPFChIsoFrix5 = 0;
   phoPFChIsoFrix6 = 0;
   phoPFChIsoFrix7 = 0;
   phoPFChIsoFrix8 = 0;
   phoPFPhoIsoFrix1 = 0;
   phoPFPhoIsoFrix2 = 0;
   phoPFPhoIsoFrix3 = 0;
   phoPFPhoIsoFrix4 = 0;
   phoPFPhoIsoFrix5 = 0;
   phoPFPhoIsoFrix6 = 0;
   phoPFPhoIsoFrix7 = 0;
   phoPFPhoIsoFrix8 = 0;
   phoPFNeuIsoFrix1 = 0;
   phoPFNeuIsoFrix2 = 0;
   phoPFNeuIsoFrix3 = 0;
   phoPFNeuIsoFrix4 = 0;
   phoPFNeuIsoFrix5 = 0;
   phoPFNeuIsoFrix6 = 0;
   phoPFNeuIsoFrix7 = 0;
   phoPFNeuIsoFrix8 = 0;
   phoEcalRecHitSumEtConeDR03 = 0;
   phohcalDepth1TowerSumEtConeDR03 = 0;
   phohcalDepth2TowerSumEtConeDR03 = 0;
   phohcalTowerSumEtConeDR03 = 0;
   photrkSumPtHollowConeDR03 = 0;
   photrkSumPtSolidConeDR03 = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoIEta = 0;
   phoIPhi = 0;
   phomaxXtalenergyFull5x5 = 0;
   phoseedTimeFull5x5 = 0;
   phomaxXtalenergy = 0;
   phoseedTime = 0;
   phoxtalBits = 0;
   phoIDbit = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIPhi = 0;
   eleSigmaIPhiIPhi = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   elePFMiniIso = 0;
   eleIDMVANonTrg = 0;
   eleIDMVATrg = 0;
   eledEtaseedAtVtx = 0;
   eleE1x5 = 0;
   eleE2x5 = 0;
   eleE5x5 = 0;
   eleE1x5Full5x5 = 0;
   eleE2x5Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   eleDr03EcalRecHitSumEt = 0;
   eleDr03HcalDepth1TowerSumEt = 0;
   eleDr03HcalDepth2TowerSumEt = 0;
   eleDr03HcalTowerSumEt = 0;
   eleDr03TkSumPt = 0;
   elecaloEnergy = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFPt = 0;
   eleGSFEta = 0;
   eleGSFPhi = 0;
   eleGSFCharge = 0;
   eleGSFHits = 0;
   eleGSFMissHits = 0;
   eleGSFNHitsMax = 0;
   eleGSFVtxProb = 0;
   eleGSFlxyPV = 0;
   eleGSFlxyBS = 0;
   eleBCEn = 0;
   eleBCEta = 0;
   eleBCPhi = 0;
   eleBCS25 = 0;
   eleBCS15 = 0;
   eleBCSieie = 0;
   eleBCSieip = 0;
   eleBCSipip = 0;
   eleFiredTrgs = 0;
   eleIDbit = 0;
   eleESEnP1Raw = 0;
   eleESEnP2Raw = 0;
   eleESEnEta = 0;
   eleESEnPhi = 0;
   eleESEnZ = 0;
   eleESEnP = 0;
   eleESEnX = 0;
   eleESEnY = 0;
   eleESEnS = 0;
   eleESEnE = 0;
   gsfPt = 0;
   gsfEta = 0;
   gsfPhi = 0;
   pfHFEn = 0;
   pfHFECALEn = 0;
   pfHFHCALEn = 0;
   pfHFPt = 0;
   pfHFEta = 0;
   pfHFPhi = 0;
   pfHFIso = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIsLooseID = 0;
   muIsMediumID = 0;
   muIsTightID = 0;
   muIsSoftID = 0;
   muIsHighPtID = 0;
   muD0 = 0;
   muDz = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muPFMiniIso = 0;
   muFiredTrgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   muIsPFMuon = 0;
   muIsGlobalMuon = 0;
   muIsTrackerMuon = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetpfCombinedInclusiveSecondaryVertexV2BJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetpfCombinedMVABJetTags = 0;
   jetPartonID = 0;
   jetGenJetIndex = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetPFLooseId = 0;
   jetPUidFullDiscriminant = 0;
   jetJECUnc = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetMUF = 0;
   jetMass = 0;
   jetPx = 0;
   jetPy = 0;
   jetPz = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtxNtrks = 0;
   jetVtx3DVal = 0;
   jetVtx3DSig = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   AK8JetPt = 0;
   AK8JetEn = 0;
   AK8JetRawPt = 0;
   AK8JetRawEn = 0;
   AK8JetEta = 0;
   AK8JetPhi = 0;
   AK8JetMass = 0;
   AK8Jet_tau1 = 0;
   AK8Jet_tau2 = 0;
   AK8Jet_tau3 = 0;
   AK8JetCHF = 0;
   AK8JetNHF = 0;
   AK8JetCEF = 0;
   AK8JetNEF = 0;
   AK8JetNCH = 0;
   AK8Jetnconstituents = 0;
   AK8JetMUF = 0;
   AK8JetPx = 0;
   AK8JetPy = 0;
   AK8JetPz = 0;
   AK8JetPFLooseId = 0;
   AK8JetPFTightLepVetoId = 0;
   AK8CHSSoftDropJetMass = 0;
   AK8CHSSoftDropJetMassCorr = 0;
   AK8CHSPrunedJetMass = 0;
   AK8CHSPrunedJetMassCorr = 0;
   AK8JetpfBoostedDSVBTag = 0;
   AK8JetJECUnc = 0;
   AK8JetL2L3corr = 0;
   AK8JetPartonID = 0;
   AK8JetGenJetIndex = 0;
   AK8JetGenJetEn = 0;
   AK8JetGenJetPt = 0;
   AK8JetGenJetEta = 0;
   AK8JetGenJetPhi = 0;
   AK8JetGenPartonID = 0;
   AK8JetGenEn = 0;
   AK8JetGenPt = 0;
   AK8JetGenEta = 0;
   AK8JetGenPhi = 0;
   AK8JetGenPartonMomID = 0;
   nAK8softdropSubjet = 0;
   AK8softdropSubjetPt = 0;
   AK8softdropSubjetEta = 0;
   AK8softdropSubjetPhi = 0;
   AK8softdropSubjetMass = 0;
   AK8softdropSubjetE = 0;
   AK8softdropSubjetCharge = 0;
   AK8softdropSubjetFlavour = 0;
   AK8softdropSubjetCSV = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("hasGoodVtx", &hasGoodVtx, &b_hasGoodVtx);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1MESUp", &pfMET_T1MESUp, &b_pfMET_T1MESUp);
   fChain->SetBranchAddress("pfMET_T1MESDo", &pfMET_T1MESDo, &b_pfMET_T1MESDo);
   fChain->SetBranchAddress("pfMET_T1EESUp", &pfMET_T1EESUp, &b_pfMET_T1EESUp);
   fChain->SetBranchAddress("pfMET_T1EESDo", &pfMET_T1EESDo, &b_pfMET_T1EESDo);
   fChain->SetBranchAddress("pfMET_T1PESUp", &pfMET_T1PESUp, &b_pfMET_T1PESUp);
   fChain->SetBranchAddress("pfMET_T1PESDo", &pfMET_T1PESDo, &b_pfMET_T1PESDo);
   fChain->SetBranchAddress("pfMET_T1TESUp", &pfMET_T1TESUp, &b_pfMET_T1TESUp);
   fChain->SetBranchAddress("pfMET_T1TESDo", &pfMET_T1TESDo, &b_pfMET_T1TESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   fChain->SetBranchAddress("phoE1x5Full5x5", &phoE1x5Full5x5, &b_phoE1x5Full5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE2x5MaxFull5x5", &phoE2x5MaxFull5x5, &b_phoE2x5MaxFull5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoSeedBCE", &phoSeedBCE, &b_phoSeedBCE);
   fChain->SetBranchAddress("phoSeedBCEta", &phoSeedBCEta, &b_phoSeedBCEta);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFChIsoFrix1", &phoPFChIsoFrix1, &b_phoPFChIsoFrix1);
   fChain->SetBranchAddress("phoPFChIsoFrix2", &phoPFChIsoFrix2, &b_phoPFChIsoFrix2);
   fChain->SetBranchAddress("phoPFChIsoFrix3", &phoPFChIsoFrix3, &b_phoPFChIsoFrix3);
   fChain->SetBranchAddress("phoPFChIsoFrix4", &phoPFChIsoFrix4, &b_phoPFChIsoFrix4);
   fChain->SetBranchAddress("phoPFChIsoFrix5", &phoPFChIsoFrix5, &b_phoPFChIsoFrix5);
   fChain->SetBranchAddress("phoPFChIsoFrix6", &phoPFChIsoFrix6, &b_phoPFChIsoFrix6);
   fChain->SetBranchAddress("phoPFChIsoFrix7", &phoPFChIsoFrix7, &b_phoPFChIsoFrix7);
   fChain->SetBranchAddress("phoPFChIsoFrix8", &phoPFChIsoFrix8, &b_phoPFChIsoFrix8);
   fChain->SetBranchAddress("phoPFPhoIsoFrix1", &phoPFPhoIsoFrix1, &b_phoPFPhoIsoFrix1);
   fChain->SetBranchAddress("phoPFPhoIsoFrix2", &phoPFPhoIsoFrix2, &b_phoPFPhoIsoFrix2);
   fChain->SetBranchAddress("phoPFPhoIsoFrix3", &phoPFPhoIsoFrix3, &b_phoPFPhoIsoFrix3);
   fChain->SetBranchAddress("phoPFPhoIsoFrix4", &phoPFPhoIsoFrix4, &b_phoPFPhoIsoFrix4);
   fChain->SetBranchAddress("phoPFPhoIsoFrix5", &phoPFPhoIsoFrix5, &b_phoPFPhoIsoFrix5);
   fChain->SetBranchAddress("phoPFPhoIsoFrix6", &phoPFPhoIsoFrix6, &b_phoPFPhoIsoFrix6);
   fChain->SetBranchAddress("phoPFPhoIsoFrix7", &phoPFPhoIsoFrix7, &b_phoPFPhoIsoFrix7);
   fChain->SetBranchAddress("phoPFPhoIsoFrix8", &phoPFPhoIsoFrix8, &b_phoPFPhoIsoFrix8);
   fChain->SetBranchAddress("phoPFNeuIsoFrix1", &phoPFNeuIsoFrix1, &b_phoPFNeuIsoFrix1);
   fChain->SetBranchAddress("phoPFNeuIsoFrix2", &phoPFNeuIsoFrix2, &b_phoPFNeuIsoFrix2);
   fChain->SetBranchAddress("phoPFNeuIsoFrix3", &phoPFNeuIsoFrix3, &b_phoPFNeuIsoFrix3);
   fChain->SetBranchAddress("phoPFNeuIsoFrix4", &phoPFNeuIsoFrix4, &b_phoPFNeuIsoFrix4);
   fChain->SetBranchAddress("phoPFNeuIsoFrix5", &phoPFNeuIsoFrix5, &b_phoPFNeuIsoFrix5);
   fChain->SetBranchAddress("phoPFNeuIsoFrix6", &phoPFNeuIsoFrix6, &b_phoPFNeuIsoFrix6);
   fChain->SetBranchAddress("phoPFNeuIsoFrix7", &phoPFNeuIsoFrix7, &b_phoPFNeuIsoFrix7);
   fChain->SetBranchAddress("phoPFNeuIsoFrix8", &phoPFNeuIsoFrix8, &b_phoPFNeuIsoFrix8);
   fChain->SetBranchAddress("phoEcalRecHitSumEtConeDR03", &phoEcalRecHitSumEtConeDR03, &b_phoEcalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03, &b_phohcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03, &b_phohcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("phohcalTowerSumEtConeDR03", &phohcalTowerSumEtConeDR03, &b_phohcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("photrkSumPtHollowConeDR03", &photrkSumPtHollowConeDR03, &b_photrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("photrkSumPtSolidConeDR03", &photrkSumPtSolidConeDR03, &b_photrkSumPtSolidConeDR03);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoIEta", &phoIEta, &b_phoIEta);
   fChain->SetBranchAddress("phoIPhi", &phoIPhi, &b_phoIPhi);
   fChain->SetBranchAddress("phomaxXtalenergyFull5x5", &phomaxXtalenergyFull5x5, &b_phomaxXtalenergyFull5x5);
   fChain->SetBranchAddress("phoseedTimeFull5x5", &phoseedTimeFull5x5, &b_phoseedTimeFull5x5);
   fChain->SetBranchAddress("phomaxXtalenergy", &phomaxXtalenergy, &b_phomaxXtalenergy);
   fChain->SetBranchAddress("phoseedTime", &phoseedTime, &b_phoseedTime);
   fChain->SetBranchAddress("phoxtalBits", &phoxtalBits, &b_phoxtalBits);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("elePFMiniIso", &elePFMiniIso, &b_elePFMiniIso);
   fChain->SetBranchAddress("eleIDMVANonTrg", &eleIDMVANonTrg, &b_eleIDMVANonTrg);
   fChain->SetBranchAddress("eleIDMVATrg", &eleIDMVATrg, &b_eleIDMVATrg);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   fChain->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   fChain->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   fChain->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   fChain->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   fChain->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   fChain->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   fChain->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   fChain->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   fChain->SetBranchAddress("eleBCEn", &eleBCEn, &b_eleBCEn);
   fChain->SetBranchAddress("eleBCEta", &eleBCEta, &b_eleBCEta);
   fChain->SetBranchAddress("eleBCPhi", &eleBCPhi, &b_eleBCPhi);
   fChain->SetBranchAddress("eleBCS25", &eleBCS25, &b_eleBCS25);
   fChain->SetBranchAddress("eleBCS15", &eleBCS15, &b_eleBCS15);
   fChain->SetBranchAddress("eleBCSieie", &eleBCSieie, &b_eleBCSieie);
   fChain->SetBranchAddress("eleBCSieip", &eleBCSieip, &b_eleBCSieip);
   fChain->SetBranchAddress("eleBCSipip", &eleBCSipip, &b_eleBCSipip);
   fChain->SetBranchAddress("eleFiredTrgs", &eleFiredTrgs, &b_eleFiredTrgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleESEnP1Raw", &eleESEnP1Raw, &b_eleESEnP1Raw);
   fChain->SetBranchAddress("eleESEnP2Raw", &eleESEnP2Raw, &b_eleESEnP2Raw);
   fChain->SetBranchAddress("eleESEnEta", &eleESEnEta, &b_eleESEnEta);
   fChain->SetBranchAddress("eleESEnPhi", &eleESEnPhi, &b_eleESEnPhi);
   fChain->SetBranchAddress("eleESEnZ", &eleESEnZ, &b_eleESEnZ);
   fChain->SetBranchAddress("eleESEnP", &eleESEnP, &b_eleESEnP);
   fChain->SetBranchAddress("eleESEnX", &eleESEnX, &b_eleESEnX);
   fChain->SetBranchAddress("eleESEnY", &eleESEnY, &b_eleESEnY);
   fChain->SetBranchAddress("eleESEnS", &eleESEnS, &b_eleESEnS);
   fChain->SetBranchAddress("eleESEnE", &eleESEnE, &b_eleESEnE);
   fChain->SetBranchAddress("nGSFTrk", &nGSFTrk, &b_nGSFTrk);
   fChain->SetBranchAddress("gsfPt", &gsfPt, &b_gsfPt);
   fChain->SetBranchAddress("gsfEta", &gsfEta, &b_gsfEta);
   fChain->SetBranchAddress("gsfPhi", &gsfPhi, &b_gsfPhi);
   fChain->SetBranchAddress("npfHF", &npfHF, &b_npfHF);
   fChain->SetBranchAddress("pfHFEn", &pfHFEn, &b_pfHFEn);
   fChain->SetBranchAddress("pfHFECALEn", &pfHFECALEn, &b_pfHFECALEn);
   fChain->SetBranchAddress("pfHFHCALEn", &pfHFHCALEn, &b_pfHFHCALEn);
   fChain->SetBranchAddress("pfHFPt", &pfHFPt, &b_pfHFPt);
   fChain->SetBranchAddress("pfHFEta", &pfHFEta, &b_pfHFEta);
   fChain->SetBranchAddress("pfHFPhi", &pfHFPhi, &b_pfHFPhi);
   fChain->SetBranchAddress("pfHFIso", &pfHFIso, &b_pfHFIso);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIsLooseID", &muIsLooseID, &b_muIsLooseID);
   fChain->SetBranchAddress("muIsMediumID", &muIsMediumID, &b_muIsMediumID);
   fChain->SetBranchAddress("muIsTightID", &muIsTightID, &b_muIsTightID);
   fChain->SetBranchAddress("muIsSoftID", &muIsSoftID, &b_muIsSoftID);
   fChain->SetBranchAddress("muIsHighPtID", &muIsHighPtID, &b_muIsHighPtID);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muPFMiniIso", &muPFMiniIso, &b_muPFMiniIso);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("muIsPFMuon", &muIsPFMuon, &b_muIsPFMuon);
   fChain->SetBranchAddress("muIsGlobalMuon", &muIsGlobalMuon, &b_muIsGlobalMuon);
   fChain->SetBranchAddress("muIsTrackerMuon", &muIsTrackerMuon, &b_muIsTrackerMuon);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetpfCombinedInclusiveSecondaryVertexV2BJetTags", &jetpfCombinedInclusiveSecondaryVertexV2BJetTags, &b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetpfCombinedMVABJetTags", &jetpfCombinedMVABJetTags, &b_jetpfCombinedMVABJetTags);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetPUidFullDiscriminant", &jetPUidFullDiscriminant, &b_jetPUidFullDiscriminant);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetPx", &jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", &jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", &jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
   fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
   fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("nAK8Jet", &nAK8Jet, &b_nAK8Jet);
   fChain->SetBranchAddress("AK8JetPt", &AK8JetPt, &b_AK8JetPt);
   fChain->SetBranchAddress("AK8JetEn", &AK8JetEn, &b_AK8JetEn);
   fChain->SetBranchAddress("AK8JetRawPt", &AK8JetRawPt, &b_AK8JetRawPt);
   fChain->SetBranchAddress("AK8JetRawEn", &AK8JetRawEn, &b_AK8JetRawEn);
   fChain->SetBranchAddress("AK8JetEta", &AK8JetEta, &b_AK8JetEta);
   fChain->SetBranchAddress("AK8JetPhi", &AK8JetPhi, &b_AK8JetPhi);
   fChain->SetBranchAddress("AK8JetMass", &AK8JetMass, &b_AK8JetMass);
   fChain->SetBranchAddress("AK8Jet_tau1", &AK8Jet_tau1, &b_AK8Jet_tau1);
   fChain->SetBranchAddress("AK8Jet_tau2", &AK8Jet_tau2, &b_AK8Jet_tau2);
   fChain->SetBranchAddress("AK8Jet_tau3", &AK8Jet_tau3, &b_AK8Jet_tau3);
   fChain->SetBranchAddress("AK8JetCHF", &AK8JetCHF, &b_AK8JetCHF);
   fChain->SetBranchAddress("AK8JetNHF", &AK8JetNHF, &b_AK8JetNHF);
   fChain->SetBranchAddress("AK8JetCEF", &AK8JetCEF, &b_AK8JetCEF);
   fChain->SetBranchAddress("AK8JetNEF", &AK8JetNEF, &b_AK8JetNEF);
   fChain->SetBranchAddress("AK8JetNCH", &AK8JetNCH, &b_AK8JetNCH);
   fChain->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
   fChain->SetBranchAddress("AK8JetMUF", &AK8JetMUF, &b_AK8JetMUF);
   fChain->SetBranchAddress("AK8JetPx", &AK8JetPx, &b_AK8JetPx);
   fChain->SetBranchAddress("AK8JetPy", &AK8JetPy, &b_AK8JetPy);
   fChain->SetBranchAddress("AK8JetPz", &AK8JetPz, &b_AK8JetPz);
   fChain->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
   fChain->SetBranchAddress("AK8JetPFTightLepVetoId", &AK8JetPFTightLepVetoId, &b_AK8JetPFTightLepVetoId);
   fChain->SetBranchAddress("AK8CHSSoftDropJetMass", &AK8CHSSoftDropJetMass, &b_AK8CHSSoftDropJetMass);
   fChain->SetBranchAddress("AK8CHSSoftDropJetMassCorr", &AK8CHSSoftDropJetMassCorr, &b_AK8CHSSoftDropJetMassCorr);
   fChain->SetBranchAddress("AK8CHSPrunedJetMass", &AK8CHSPrunedJetMass, &b_AK8CHSPrunedJetMass);
   fChain->SetBranchAddress("AK8CHSPrunedJetMassCorr", &AK8CHSPrunedJetMassCorr, &b_AK8CHSPrunedJetMassCorr);
   fChain->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
   fChain->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
   fChain->SetBranchAddress("AK8JetL2L3corr", &AK8JetL2L3corr, &b_AK8JetL2L3corr);
   fChain->SetBranchAddress("AK8JetPartonID", &AK8JetPartonID, &b_AK8JetPartonID);
   fChain->SetBranchAddress("AK8JetGenJetIndex", &AK8JetGenJetIndex, &b_AK8JetGenJetIndex);
   fChain->SetBranchAddress("AK8JetGenJetEn", &AK8JetGenJetEn, &b_AK8JetGenJetEn);
   fChain->SetBranchAddress("AK8JetGenJetPt", &AK8JetGenJetPt, &b_AK8JetGenJetPt);
   fChain->SetBranchAddress("AK8JetGenJetEta", &AK8JetGenJetEta, &b_AK8JetGenJetEta);
   fChain->SetBranchAddress("AK8JetGenJetPhi", &AK8JetGenJetPhi, &b_AK8JetGenJetPhi);
   fChain->SetBranchAddress("AK8JetGenPartonID", &AK8JetGenPartonID, &b_AK8JetGenPartonID);
   fChain->SetBranchAddress("AK8JetGenEn", &AK8JetGenEn, &b_AK8JetGenEn);
   fChain->SetBranchAddress("AK8JetGenPt", &AK8JetGenPt, &b_AK8JetGenPt);
   fChain->SetBranchAddress("AK8JetGenEta", &AK8JetGenEta, &b_AK8JetGenEta);
   fChain->SetBranchAddress("AK8JetGenPhi", &AK8JetGenPhi, &b_AK8JetGenPhi);
   fChain->SetBranchAddress("AK8JetGenPartonMomID", &AK8JetGenPartonMomID, &b_AK8JetGenPartonMomID);
   fChain->SetBranchAddress("nAK8softdropSubjet", &nAK8softdropSubjet, &b_nAK8softdropSubjet);
   fChain->SetBranchAddress("AK8softdropSubjetPt", &AK8softdropSubjetPt, &b_AK8softdropSubjetPt);
   fChain->SetBranchAddress("AK8softdropSubjetEta", &AK8softdropSubjetEta, &b_AK8softdropSubjetEta);
   fChain->SetBranchAddress("AK8softdropSubjetPhi", &AK8softdropSubjetPhi, &b_AK8softdropSubjetPhi);
   fChain->SetBranchAddress("AK8softdropSubjetMass", &AK8softdropSubjetMass, &b_AK8softdropSubjetMass);
   fChain->SetBranchAddress("AK8softdropSubjetE", &AK8softdropSubjetE, &b_AK8softdropSubjetE);
   fChain->SetBranchAddress("AK8softdropSubjetCharge", &AK8softdropSubjetCharge, &b_AK8softdropSubjetCharge);
   fChain->SetBranchAddress("AK8softdropSubjetFlavour", &AK8softdropSubjetFlavour, &b_AK8softdropSubjetFlavour);
   fChain->SetBranchAddress("AK8softdropSubjetCSV", &AK8softdropSubjetCSV, &b_AK8softdropSubjetCSV);
   Notify();
}

Bool_t PhoJet_Analyzer_MC::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void PhoJet_Analyzer_MC::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t PhoJet_Analyzer_MC::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}


void PhoJet_Analyzer_MC::BookHistos(){

    // Define Histograms here +++++++++++++++++++++++++++++
    f1->cd(); 

    char name[100];
    const Int_t nMassBins = 119;
    const Double_t MassBin[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

    std::string cut0[2] = {"noMassCut", "MassCut"};
    for( Int_t hist = 0; hist < 2; ++hist){
	sprintf(name, "h_ptPhoton_noPU_%s",cut0[hist].c_str());
	h_ptPhoton_noPU[hist]  = new TH1F(name,"pt of photon",120,0.0,4800.0);  //
	h_ptPhoton_noPU[hist]->GetYaxis()->SetTitle("Events/40 GeV");           h_ptPhoton_noPU[hist]->GetYaxis()->CenterTitle();
	h_ptPhoton_noPU[hist]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton_noPU[hist]->GetXaxis()->CenterTitle();
	h_ptPhoton_noPU[hist]->Sumw2();

	sprintf(name, "h_ptJet_noPU_%s",cut0[hist].c_str());
	h_ptJet_noPU[hist]     =new TH1F(name,"pt of jet",120,0.0,4800.0);  // 
	h_ptJet_noPU[hist]->GetYaxis()->SetTitle("Events/40 GeV");       h_ptJet_noPU[hist]->GetYaxis()->CenterTitle();
	h_ptJet_noPU[hist]->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet_noPU[hist]->GetXaxis()->CenterTitle();
	h_ptJet_noPU[hist]->Sumw2();

	sprintf(name, "h_mass_VarBin_noPU_%s",cut0[hist].c_str());
	h_mass_VarBin_noPU[hist]      =new TH1F(name,"mass of photon+jet",nMassBins,MassBin);
	h_mass_VarBin_noPU[hist]->GetYaxis()->SetTitle("Events");                      h_mass_VarBin_noPU[hist]->GetYaxis()->CenterTitle();
	h_mass_VarBin_noPU[hist]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_VarBin_noPU[hist]->GetXaxis()->CenterTitle();
	h_mass_VarBin_noPU[hist]->Sumw2(); 
    }

    std::string cut1[3] = {"noMassCut", "MassCut", "M_1p2to2p2"};
    for( Int_t hist = 0; hist < 3; ++hist){
	sprintf(name, "h_ptPhoton_%s",cut1[hist].c_str());
	h_ptPhoton[hist]  = new TH1F(name,"pt of photon",120,0.0,4800.0);  //
	h_ptPhoton[hist]->GetYaxis()->SetTitle("Events/40 GeV");           h_ptPhoton[hist]->GetYaxis()->CenterTitle();
	h_ptPhoton[hist]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton[hist]->GetXaxis()->CenterTitle();
	h_ptPhoton[hist]->Sumw2();

	sprintf(name, "h_ptJet_%s",cut1[hist].c_str());
	h_ptJet[hist]     =new TH1F(name,"pt of jet",120,0.0,4800.0);  // 
	h_ptJet[hist]->GetYaxis()->SetTitle("Events/40 GeV");       h_ptJet[hist]->GetYaxis()->CenterTitle();
	h_ptJet[hist]->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet[hist]->GetXaxis()->CenterTitle();
	h_ptJet[hist]->Sumw2();

	sprintf(name, "h_mass_VarBin_%s",cut1[hist].c_str());
	h_mass_VarBin[hist]      =new TH1F(name,"mass of photon+jet",nMassBins,MassBin);
	h_mass_VarBin[hist]->GetYaxis()->SetTitle("Events");                      h_mass_VarBin[hist]->GetYaxis()->CenterTitle();
	h_mass_VarBin[hist]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_VarBin[hist]->GetXaxis()->CenterTitle();
	h_mass_VarBin[hist]->Sumw2(); 

	sprintf(name, "h_mass_bin1_%s",cut1[hist].c_str());
	h_mass_bin1[hist]     = new TH1F(name,"mass of photon+jet",14000,0.0,14000.0);
	h_mass_bin1[hist]->GetYaxis()->SetTitle("Events/1 GeV");                h_mass_bin1[hist]->GetYaxis()->CenterTitle();
	h_mass_bin1[hist]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_bin1[hist]->GetXaxis()->CenterTitle();
	h_mass_bin1[hist]->Sumw2();

	sprintf(name, "h_etaPhoton_%s",cut1[hist].c_str());
	h_etaPhoton[hist] =new TH1F(name,"eta of photon",50,-2.5,2.5);
	h_etaPhoton[hist]->GetYaxis()->SetTitle("Events");        h_etaPhoton[hist]->GetYaxis()->CenterTitle();
	h_etaPhoton[hist]->GetXaxis()->SetTitle("#eta^{#gamma}"); h_etaPhoton[hist]->GetXaxis()->CenterTitle();
	h_etaPhoton[hist]->Sumw2();

	sprintf(name, "h_etaJet_%s",cut1[hist].c_str());
	h_etaJet[hist]    =new TH1F(name,"eta of jet      ",60,-3.0,3.0);
	h_etaJet[hist]->GetYaxis()->SetTitle("Events");      h_etaJet[hist]->GetYaxis()->CenterTitle();
	h_etaJet[hist]->GetXaxis()->SetTitle("#eta^{jet}");  h_etaJet[hist]->GetXaxis()->CenterTitle();
	h_etaJet[hist]->Sumw2();

	sprintf(name, "h_phiPhoton_%s",cut1[hist].c_str());
	h_phiPhoton[hist] =new TH1F(name,"phi of photon",50,-3.2,3.2);
	h_phiPhoton[hist]->GetYaxis()->SetTitle("Events");        h_phiPhoton[hist]->GetYaxis()->CenterTitle();
	h_phiPhoton[hist]->GetXaxis()->SetTitle("#phi^{#gamma}"); h_phiPhoton[hist]->GetXaxis()->CenterTitle();
	h_phiPhoton[hist]->Sumw2();

	sprintf(name, "h_phiJet_%s",cut1[hist].c_str());
	h_phiJet[hist]    =new TH1F(name,"phi of jet",50,-3.2,3.2);
	h_phiJet[hist]->GetYaxis()->SetTitle("Events");      h_phiJet[hist]->GetYaxis()->CenterTitle();
	h_phiJet[hist]->GetXaxis()->SetTitle("#phi^{jet}");  h_phiJet[hist]->GetXaxis()->CenterTitle();
	h_phiJet[hist]->Sumw2();

	sprintf(name, "h_PFMet_%s",cut1[hist].c_str());
	h_PFMet[hist]  = new TH1F(name,"pt of PFMet",200,0.0,1000.0);
	h_PFMet[hist]->GetYaxis()->SetTitle("Events/5 GeV");           h_PFMet[hist]->GetYaxis()->CenterTitle();
	h_PFMet[hist]->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");    h_PFMet[hist]->GetXaxis()->CenterTitle();
	h_PFMet[hist]->Sumw2();

	sprintf(name, "h_SumPFMet_%s",cut1[hist].c_str());
	h_SumPFMet[hist]  = new TH1F(name,"SumET PF Met",80,0.0,4000.0);
	h_SumPFMet[hist]->GetYaxis()->SetTitle("Events/5 GeV");               h_SumPFMet[hist]->GetYaxis()->CenterTitle();
	h_SumPFMet[hist]->GetXaxis()->SetTitle("#sum#slash{E}_{T} (GeV)");    h_SumPFMet[hist]->GetXaxis()->CenterTitle();
	h_SumPFMet[hist]->Sumw2();

	sprintf(name, "h_MetBySumMET_%s",cut1[hist].c_str());
	h_MetBySumMET[hist]  = new TH1F(name,"MET / SumET PF Met",50,0.0,1.0);
	h_MetBySumMET[hist]->GetYaxis()->SetTitle("Events");                             h_MetBySumMET[hist]->GetYaxis()->CenterTitle();
	h_MetBySumMET[hist]->GetXaxis()->SetTitle("#slash{E}_{T}/#sum#slash{E}_{T}");    h_MetBySumMET[hist]->GetXaxis()->CenterTitle();
	h_MetBySumMET[hist]->Sumw2();

	sprintf(name, "h_pfMetgjMass_%s",cut1[hist].c_str());
	h_pfMetgjMass[hist]    =new TH2F(name,"MET vs GJ Mass ",120,0.0,4800.0,100,0.0,1000.0);
	h_pfMetgjMass[hist]->GetYaxis()->SetTitle("PFMET");                         h_pfMetgjMass[hist]->GetYaxis()->CenterTitle();
	h_pfMetgjMass[hist]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");     h_pfMetgjMass[hist]->GetXaxis()->CenterTitle();
	h_pfMetgjMass[hist]->Sumw2();

	sprintf(name, "h_PtPhotJet_%s",cut1[hist].c_str());
	h_PtPhotJet[hist]    =new TH2F(name,"Pt of Photon vs Jet ",120,0.0,4800.0,120,0.0,4800.0);
	h_PtPhotJet[hist]->GetYaxis()->SetTitle("P_{T}^{Jet}");     h_PtPhotJet[hist]->GetYaxis()->CenterTitle();
	h_PtPhotJet[hist]->GetXaxis()->SetTitle("P_{T}^{#gamma}");  h_PtPhotJet[hist]->GetXaxis()->CenterTitle();
	h_PtPhotJet[hist]->Sumw2();

	sprintf(name, "h_etaPhotJet_%s",cut1[hist].c_str());
	h_etaPhotJet[hist]    =new TH2F(name,"eta of Photon vs Jet ",100,-2.5,2.5,100,-2.5,2.5);
	h_etaPhotJet[hist]->GetYaxis()->SetTitle("#eta^{Jet}");     h_etaPhotJet[hist]->GetYaxis()->CenterTitle();
	h_etaPhotJet[hist]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_etaPhotJet[hist]->GetXaxis()->CenterTitle();
	h_etaPhotJet[hist]->Sumw2();

	sprintf(name, "h_DR_PhotonJet_%s",cut1[hist].c_str());
	h_DR_PhotonJet[hist] = new TH1F(name,"DeltaR between photon n jet",100,0.0,10.0);
	h_DR_PhotonJet[hist]->GetYaxis()->SetTitle("Events");        h_DR_PhotonJet[hist]->GetYaxis()->CenterTitle();
	h_DR_PhotonJet[hist]->GetXaxis()->SetTitle("#Delta R");      h_DR_PhotonJet[hist]->GetXaxis()->CenterTitle();
	h_DR_PhotonJet[hist]->Sumw2();

	sprintf(name, "h_cosThetaStar_%s",cut1[hist].c_str());
	h_cosThetaStar[hist]      =new TH1F(name,"Cos theta star of photon+jet",50,0,1);
	h_cosThetaStar[hist]->GetYaxis()->SetTitle("Events");        h_cosThetaStar[hist]->GetYaxis()->CenterTitle();
	h_cosThetaStar[hist]->GetXaxis()->SetTitle("cos#theta*");    h_cosThetaStar[hist]->GetXaxis()->CenterTitle();
	h_cosThetaStar[hist]->Sumw2();

	sprintf(name, "h_dEta_%s",cut1[hist].c_str());
	h_dEta[hist]      =new TH1F(name,"dEta of photon+jet",120,0,6);
	h_dEta[hist]->GetYaxis()->SetTitle("Events");        h_dEta[hist]->GetYaxis()->CenterTitle();
	h_dEta[hist]->GetXaxis()->SetTitle("#Delta #eta");   h_dEta[hist]->GetXaxis()->CenterTitle();
	h_dEta[hist]->Sumw2();

	sprintf(name, "h_dphi_%s",cut1[hist].c_str());
	h_dphi[hist]      =new TH1F(name,"dphi of photon+jet",64,0,3.2);
	h_dphi[hist]->GetYaxis()->SetTitle("Events");        h_dphi[hist]->GetYaxis()->CenterTitle();
	h_dphi[hist]->GetXaxis()->SetTitle("#Delta #phi");   h_dphi[hist]->GetXaxis()->CenterTitle();
	h_dphi[hist]->Sumw2();

	sprintf(name, "h_Photon_SigmaIetaIeta_%s",cut1[hist].c_str());
	h_Photon_SigmaIetaIeta[hist]  = new TH1F(name,"Photon SigmaIetaIeta",100,0.0,0.05);
	h_Photon_SigmaIetaIeta[hist]->GetYaxis()->SetTitle("Events");                 h_Photon_SigmaIetaIeta[hist]->GetYaxis()->CenterTitle();
	h_Photon_SigmaIetaIeta[hist]->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");   h_Photon_SigmaIetaIeta[hist]->GetXaxis()->CenterTitle(); 
	h_Photon_SigmaIetaIeta[hist]->Sumw2();

	sprintf(name, "h_HoE_%s",cut1[hist].c_str());
	h_HoE[hist]  = new TH1F(name,"PFIso HoE Barrel ",50,0,0.1);
	h_HoE[hist]->GetYaxis()->SetTitle("Events");              h_HoE[hist]->GetYaxis()->CenterTitle();
	h_HoE[hist]->GetXaxis()->SetTitle("PFIso HoE");           h_HoE[hist]->GetXaxis()->CenterTitle();
	h_HoE[hist]->Sumw2();

	sprintf(name, "h_CorrPFiso_Charged_%s",cut1[hist].c_str());
	h_CorrPFiso_Charged[hist] = new TH1F(name,"Rho Corrected PFIso Charged  ",300,0,24);
	h_CorrPFiso_Charged[hist]->GetYaxis()->SetTitle("Events");                        h_CorrPFiso_Charged[hist]->GetYaxis()->CenterTitle();
	h_CorrPFiso_Charged[hist]->GetXaxis()->SetTitle("Corrected PFIso Charged");       h_CorrPFiso_Charged[hist]->GetXaxis()->CenterTitle();
	h_CorrPFiso_Charged[hist]->Sumw2();

	sprintf(name, "h_CorrPFiso_Neutral_%s",cut1[hist].c_str());
	h_CorrPFiso_Neutral[hist] = new TH1F(name,"Rho Corrected PFIso Neutral barrel",200,0,300);
	h_CorrPFiso_Neutral[hist]->GetYaxis()->SetTitle("Events");                        h_CorrPFiso_Neutral[hist]->GetYaxis()->CenterTitle();
	h_CorrPFiso_Neutral[hist]->GetXaxis()->SetTitle("Corrected PFIso Neutral");       h_CorrPFiso_Neutral[hist]->GetXaxis()->CenterTitle();
	h_CorrPFiso_Neutral[hist]->Sumw2();

	sprintf(name, "h_CorrPFiso_Photon_%s",cut1[hist].c_str());
	h_CorrPFiso_Photon[hist] = new TH1F(name,"Rho Corrected PFIso Photon barrel",125,0,25);
	h_CorrPFiso_Photon[hist]->GetYaxis()->SetTitle("Events");                         h_CorrPFiso_Photon[hist]->GetYaxis()->CenterTitle();
	h_CorrPFiso_Photon[hist]->GetXaxis()->SetTitle("Corrected PFIso Photon");         h_CorrPFiso_Photon[hist]->GetXaxis()->CenterTitle();
	h_CorrPFiso_Photon[hist]->Sumw2();

	sprintf(name, "h_PFiso_Electronveto_%s",cut1[hist].c_str());
	h_PFiso_Electronveto[hist]  = new TH1F(name,"PFIso Electronveto",3,0,3);
	h_PFiso_Electronveto[hist]->GetYaxis()->SetTitle("Events");              h_PFiso_Electronveto[hist]->GetYaxis()->CenterTitle();
	h_PFiso_Electronveto[hist]->GetXaxis()->SetTitle("PFIso Electron veto"); h_PFiso_Electronveto[hist]->GetXaxis()->CenterTitle();                                   
	h_PFiso_Electronveto[hist]->Sumw2();


	sprintf(name, "h_jet_NEF_%s",cut1[hist].c_str());
	h_jet_NEF[hist] = new TH1F(name,"Neutral EM Fraction",25,0,1);
	h_jet_NEF[hist]->GetYaxis()->SetTitle("");              h_jet_NEF[hist]->GetYaxis()->CenterTitle();
	h_jet_NEF[hist]->GetXaxis()->SetTitle("jet_NEF");       h_jet_NEF[hist]->GetXaxis()->CenterTitle();
	h_jet_NEF[hist]->Sumw2();

	sprintf(name, "h_jet_NHF_%s",cut1[hist].c_str());
	h_jet_NHF[hist] = new TH1F(name,"Neutral hadron Fraction",25,0,1);
	h_jet_NHF[hist]->GetYaxis()->SetTitle("");              h_jet_NHF[hist]->GetYaxis()->CenterTitle();
	h_jet_NHF[hist]->GetXaxis()->SetTitle("jet_NHF");       h_jet_NHF[hist]->GetXaxis()->CenterTitle();
	h_jet_NHF[hist]->Sumw2();

	sprintf(name, "h_jet_CEF_%s",cut1[hist].c_str());
	h_jet_CEF[hist] = new TH1F(name,"Charged EM Fraction",25,0,1);
	h_jet_CEF[hist]->GetYaxis()->SetTitle("");              h_jet_CEF[hist]->GetYaxis()->CenterTitle();
	h_jet_CEF[hist]->GetXaxis()->SetTitle("jet_CEF");       h_jet_CEF[hist]->GetXaxis()->CenterTitle();
	h_jet_CEF[hist]->Sumw2();

	sprintf(name, "h_jet_CHF_%s",cut1[hist].c_str());
	h_jet_CHF[hist] = new TH1F(name,"Charged hadron Fraction",25,0,1);
	h_jet_CHF[hist]->GetYaxis()->SetTitle("");              h_jet_CHF[hist]->GetYaxis()->CenterTitle();
	h_jet_CHF[hist]->GetXaxis()->SetTitle("jet_CHF");       h_jet_CHF[hist]->GetXaxis()->CenterTitle();
	h_jet_CHF[hist]->Sumw2();

	sprintf(name, "h_jet_NConstituents_%s",cut1[hist].c_str());
	h_jet_NConstituents[hist] = new TH1F(name,"NConstituents",50,0,50);
	h_jet_NConstituents[hist]->GetYaxis()->SetTitle("");                   h_jet_NConstituents[hist]->GetYaxis()->CenterTitle();
	h_jet_NConstituents[hist]->GetXaxis()->SetTitle("jet_NConstituents");  h_jet_NConstituents[hist]->GetXaxis()->CenterTitle();
	h_jet_NConstituents[hist]->Sumw2();

	sprintf(name, "h_jet_ChargeMultiplicity_%s",cut1[hist].c_str());
	h_jet_ChargeMultiplicity[hist] = new TH1F(name,"ChargeMultiplicity",50,0,50);
	h_jet_ChargeMultiplicity[hist]->GetYaxis()->SetTitle("");                           h_jet_ChargeMultiplicity[hist]->GetYaxis()->CenterTitle();
	h_jet_ChargeMultiplicity[hist]->GetXaxis()->SetTitle("jet_ChargeMultiplicity");     h_jet_ChargeMultiplicity[hist]->GetXaxis()->CenterTitle();
	h_jet_ChargeMultiplicity[hist]->Sumw2(); 

    }//--cut1

    // Histogram fornumber for photons and jets
    std::string cut2[4] ={ "noCut", "noMassCut", "MassCut", "M_1p2to2p2"};
    std::string cut3[4] ={ "noPU", "PU", "MCut_noPU", "MCut_PU"}; // for h_Vertices
    for(Int_t hist = 0 ; hist < 4 ; ++hist){
	sprintf(name, "h_nIsoPhoton_%s",cut2[hist].c_str());
	h_nIsoPhoton[hist] = new TH1F(name,"no of Iso Photons",10,0,10);
	h_nIsoPhoton[hist]->GetXaxis()->SetTitle("No. of Iso photons");

	sprintf(name, "h_nPhoton_%s",cut2[hist].c_str());
	h_nPhoton[hist] = new TH1F(name,"no of Photons",10,0,10);
	h_nPhoton[hist]->GetXaxis()->SetTitle("No. of photons");

	sprintf(name, "h_nIsoPhotonPt_%s",cut2[hist].c_str());
	h_nIsoPhotonPt[hist] = new TH2F(name,"no of Iso Photons vs Pt",120,0.0,4800.0,10,0.0,10.0);

	sprintf(name, "h_nPhotonPt_%s",cut2[hist].c_str());
	h_nPhotonPt[hist] = new TH2F(name,"no of Photons vs Pt",120,0.0,4800.0,10,0.0,10.0);

	sprintf(name, "h_nJet_%s",cut2[hist].c_str());
	h_nJet[hist] = new TH1F(name,"no of Jets",20,0,20);
	h_nJet[hist]->GetXaxis()->SetTitle("No. of jets");

	sprintf(name, "h_nJetPt_%s",cut2[hist].c_str());
	h_nJetPt[hist] = new TH2F(name,"no of selected jets vs Pt",120,0.0,4800.0,20,0.0,20.0);

	sprintf(name, "h_Vertices_%s",cut3[hist].c_str());
	h_Vertices[hist] = new TH1F(name,"Vertices", 100, 0, 100);
    }

    std::string cut4[2] ={ "noPU", "PU"};
    for(Int_t hist = 0 ; hist < 2 ; ++hist){
	sprintf(name, "h_trueInteractions_%s",cut4[hist].c_str());
	h_trueInteractions[hist] = new TH1F(name,"true interactions", 50, 0, 50);
    }

    h_PC = new TH1F ("h_PC","Photon Candidate", 10, 0, 10);
    h_JC = new TH1F ("h_JC","Jet Candidate", 20, 0, 20);


    h_CutFlow = new TH1F("h_CutFlow", "cut flow of selection", 11, 0, 11);
    h_CutExpFlow = new TH1F("h_CutExpFlow", "exp cut flow of selection", 11, 0, 11);

    TString cutFlowLabel[11] = {"Total", "HLT", "PrimaryVtx", "PhotonID", "PhotonPtEta", "JetID", "JetPtEta", "Dphi", "DEta", "MassCut", "M_1p2to2p2"};
    for( Int_t bin = 1; bin <= h_CutFlow->GetNbinsX(); ++bin){
	h_CutFlow->GetXaxis()->SetBinLabel(bin, cutFlowLabel[bin-1]);
	h_CutExpFlow->GetXaxis()->SetBinLabel(bin, cutFlowLabel[bin-1]);
    }

} //-- BookHistos

void PhoJet_Analyzer_MC::WriteHistos(){  // WriteHistos

    for( Int_t hist = 0; hist < 3; ++hist){

	h_ptPhoton[hist]                 ->Write();  
	h_ptJet[hist]                    ->Write();     
	h_mass_VarBin[hist]              ->Write();     
	h_mass_bin1[hist]                ->Write();     
	h_etaPhoton[hist]                ->Write(); 
	h_etaJet[hist]                   ->Write();    
	h_phiPhoton[hist]                ->Write(); 
	h_phiJet[hist]                   ->Write();    
	h_PFMet[hist]                    ->Write();  
	h_SumPFMet[hist]                 ->Write();  
	h_MetBySumMET[hist]              ->Write();  
	h_pfMetgjMass[hist]              ->Write();    
	h_PtPhotJet[hist]                ->Write();    
	h_etaPhotJet[hist]               ->Write();    
	h_DR_PhotonJet[hist]             ->Write(); 
	h_cosThetaStar[hist]             ->Write();      
	h_dEta[hist]                     ->Write();      
	h_dphi[hist]                     ->Write();      
	h_Photon_SigmaIetaIeta[hist]     ->Write();  
	h_HoE[hist]                      ->Write();  
	h_CorrPFiso_Charged[hist]        ->Write(); 
	h_CorrPFiso_Neutral[hist]        ->Write(); 
	h_CorrPFiso_Photon[hist]         ->Write(); 
	h_PFiso_Electronveto[hist]       ->Write();  
	h_jet_NEF[hist]                  ->Write(); 
	h_jet_NHF[hist]                  ->Write(); 
	h_jet_CEF[hist]                  ->Write(); 
	h_jet_CHF[hist]                  ->Write(); 
	h_jet_NConstituents[hist]        ->Write(); 
	h_jet_ChargeMultiplicity[hist]   ->Write(); 

    }//--cut1

    // Histogram fornumber for photons and jets
    for(Int_t hist = 0 ; hist < 4 ; ++hist){
	h_nIsoPhoton[hist]        ->Write();
	h_nPhoton[hist]           ->Write();
	h_nIsoPhotonPt[hist]      ->Write();
	h_nPhotonPt[hist]         ->Write();
	h_nJet[hist]              ->Write();
	h_nJetPt[hist]            ->Write();
	h_Vertices[hist]          ->Write();
    }

    h_PC          ->Write();
    h_JC          ->Write();
    h_CutFlow     ->Write();
    h_CutExpFlow  ->Write();

}

//----------------------
// Spike Cut
Bool_t PhoJet_Analyzer_MC::NoSpike(Int_t ipho){
  Bool_t passSpike = false;
  if( (*phoSigmaIEtaIEtaFull5x5)[ipho]        > 0.001  && 
      (*phoSigmaIPhiIPhiFull5x5)[ipho]        > 0.001  &&
      // fabs(getLICTD(ipho))              < 5.0    &&	  
      fabs( (*phoseedTimeFull5x5)[ipho] )     < 3.0    &&
      (*phoR9Full5x5)[ipho]                   < 1.0){
      passSpike = true ; 
  } 

  return passSpike;
}

//-------------------
//Primary Vertex
Bool_t PhoJet_Analyzer_MC::PrimaryVertex(Int_t &goodVertex){

    Bool_t passVertex = false;
    goodVertex = 0;

    for(Int_t i=0; i<nVtx; ++i){
       if(hasGoodVtx)
	  goodVertex++;
    }
    if(goodVertex > 0) passVertex = true;

    return passVertex;
}

//----------------------
// Compute deltaR
Double_t PhoJet_Analyzer_MC::getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){
    Double_t DR = 0.0;
    DR = pow( ( pow((eta1 - eta2), 2.0) + pow((phi1 - phi2), 2.0) ), 0.5);
    return DR;
}

//----------------------
// Compute deltaEta
Double_t PhoJet_Analyzer_MC::getDEta(Double_t eta1, Double_t eta2){
    Double_t dEta = fabs(eta1 - eta2);
    return dEta;
}

//----------------------
// Compute cosThetaStar
Double_t PhoJet_Analyzer_MC::getCosThetaStar(Double_t eta1, Double_t eta2){
    Double_t theta = tanh( getDEta(eta1,eta2)/2.0 );
    return theta;
}

//----------------------
// Compute deltaPhi
Double_t PhoJet_Analyzer_MC::getDPhi(Double_t phi1, Double_t phi2){
    Double_t dPhi  = fabs(phi1 - phi2);
    Double_t twopi = 2.0*(TMath::Pi());

    if(dPhi < 0) dPhi = - dPhi;
    if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;

    return dPhi;
}

//------------------------
//Compute Invariant mass
Double_t PhoJet_Analyzer_MC::getMass(Int_t p, Int_t j){

    Double_t mass = 0.0;

    TLorentzVector pho;
    pho.SetPtEtaPhiE( (*phoEt)[p], (*phoSCEta)[p], (*phoSCPhi)[p], (*phoE)[p] );

    TLorentzVector jet;
    jet.SetPtEtaPhiE( (*jetPt)[j], (*jetEta)[j], (*jetPhi)[j], (*jetEn)[j] );

    mass = (pho+jet).M();

    return mass;
}

//get Photon Candidate
Int_t PhoJet_Analyzer_MC::getPhotonCand(TString phoWP){
    Int_t pc = -1;
    Bool_t ID = false;
    for(Int_t ipho = 0; ipho < nPho; ++ipho){
	ID = CutBasedPFPhotonID(ipho, phoWP);
	if(ID){
	    pc = ipho;
	    break;
	}
    }
    return pc;
}

//-------------------------------------------
// Cut based Photon ID SPRING15 selection 25ns 
//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns
Bool_t PhoJet_Analyzer_MC::CutBasedPFPhotonID(Int_t ipho, TString phoWP){

    Bool_t photonId = false;

    if(phoWP == "tight"){  // Tight
	if( fabs((*phoSCEta)[ipho]) < 1.4442){ // EB
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0100 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               <  0.76   ) &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (0.97 + (0.014 * (*phoEt)[ipho]) + (0.000019 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.08 + (0.0053 * (*phoEt)[ipho])) ) ); 
	}
	if( fabs((*phoSCEta)[ipho]) > 1.4442){ // EE
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0268 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               < 0.56    )  &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (2.09 + (0.139 * (*phoEt)[ipho]) + (0.000025 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.16 + (0.0034 * (*phoEt)[ipho])) ) ); 
	}
    }

    if(phoWP == "medium"){ // Medium
	if( fabs((*phoSCEta)[ipho]) < 1.4442){ // EB
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0102 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               <  1.37   )  &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (1.06 + (0.014 * (*phoEt)[ipho]) + (0.000019 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.28 + (0.0053 * (*phoEt)[ipho])) ) ); 
	}
	if( fabs((*phoSCEta)[ipho]) > 1.4442){ // EE
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0268 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               <  1.10   ) &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (2.69 + (0.0139 * (*phoEt)[ipho]) + (0.000025 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.39 + (0.0034 * (*phoEt)[ipho])) ) ); 
	}
    }

    if(phoWP == "loose"){ // Loose
	if( fabs((*phoSCEta)[ipho]) < 1.4442){ // EB
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0102 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               <  3.32   )  &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (1.92 + (0.014 * (*phoEt)[ipho]) + (0.000019 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.81 + (0.0053 * (*phoEt)[ipho])) ) ); 
	}
	if( fabs((*phoSCEta)[ipho]) > 1.4442){ // EE
	    photonId = ( 
		    ((*phoHoverE)[ipho]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[ipho]  <  0.0274 ) &&
		    ((*phoEleVeto)[ipho]              ==  1      ) && 
		    ((*phoPFChIso)[ipho]               <  1.97   )  &&
		    (TMath::Max(((*phoPFNeuIso)[ipho] - rho*EAneutral((*phoSCEta)[ipho])), 0.0) < (11.86 + (0.0139 * (*phoEt)[ipho]) + (0.000025 * pow((*phoEt)[ipho], 2.0))) )  &&
		    (TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAphoton( (*phoSCEta)[ipho])), 0.0) < (0.83 + (0.0034 * (*phoEt)[ipho])) ) ); 
	}
    }

    return photonId;
}
// Effective area to be needed in PF Iso for photon ID -- SPRING15 bx25ns
Double_t PhoJet_Analyzer_MC::EAcharged(Double_t eta){
    Float_t EffectiveArea = 0.0;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0;

    return EffectiveArea;
}

Double_t PhoJet_Analyzer_MC::EAneutral(Double_t eta){
    Float_t EffectiveArea = 0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0599;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0819;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0696;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0360;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0360;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0462;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0656;

    return EffectiveArea;
}

Double_t PhoJet_Analyzer_MC::EAphoton(Double_t eta){
    Float_t EffectiveArea = 0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1271;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1101;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0756;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.1175;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.1498;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.1857;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.2183;

    return EffectiveArea;
}


//-------------------------------------------
// Jet ID
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
Bool_t PhoJet_Analyzer_MC::jetID(Int_t ijet, TString jetWP){

    Bool_t jetID = false;

    if(jetWP == "loose"){ // loose
	jetID = ( ( (*jetNHF)[ijet] < 0.99 && (*jetNEF)[ijet] < 0.99 &&  (*jetNConstituents)[ijet] > 1 ) &&
		  ( (fabs((*jetEta)[ijet]) <= 2.4 &&  (*jetCHF)[ijet] > 0.0 &&  (*jetNCH)[ijet] > 0 && (*jetCEF)[ijet] < 0.99 ) || fabs((*jetEta)[ijet]) > 2.4) && 
		  fabs((*jetEta)[ijet]) <= 3.0);
    }

    if(jetWP == "tight"){ // tight
    	jetID = ( ( (*jetNHF)[ijet] < 0.90 && (*jetNEF)[ijet] < 0.90 &&  (*jetNConstituents)[ijet] > 1 ) &&
		  ( (fabs((*jetEta)[ijet]) <= 2.4 &&  (*jetCHF)[ijet] > 0.0 &&  (*jetNCH)[ijet] > 0 && (*jetCEF)[ijet] < 0.99 ) || fabs((*jetEta)[ijet]) > 2.4) && 
		  fabs((*jetEta)[ijet]) <= 3.0);
    }

    if(jetWP == "tightLepVeto"){ //tightLepVeto  
       jetID = ( ( (*jetNHF)[ijet] < 0.90 && (*jetNEF)[ijet] < 0.90 &&  (*jetNConstituents)[ijet] > 1 && (*jetMUF)[ijet] < 0.8) &&
	     ( (fabs((*jetEta)[ijet]) <= 2.4 &&  (*jetCHF)[ijet] > 0.0 &&  (*jetNCH)[ijet] > 0 && (*jetCEF)[ijet] < 0.99 ) || fabs((*jetEta)[ijet]) > 2.4) && 
	     fabs((*jetEta)[ijet]) <= 3.0);
    }

    return jetID;
}

//-------------------------------------------
// AK8 Jet ID
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
Bool_t PhoJet_Analyzer_MC::AK8jetID(Int_t ijet, TString jetWP){

    Bool_t jetID = false;

    if(jetWP == "loose"){ // loose
	jetID = ( ( (*AK8JetNHF)[ijet] < 0.99 && (*AK8JetNEF)[ijet] < 0.99 &&  (*AK8Jetnconstituents)[ijet] > 1 ) &&
		  ( (fabs((*AK8JetEta)[ijet]) <= 2.4 &&  (*AK8JetCHF)[ijet] > 0.0 &&  (*AK8JetNCH)[ijet] > 0 && (*AK8JetCEF)[ijet] < 0.99 ) || fabs((*AK8JetEta)[ijet]) > 2.4) && 
		  fabs((*AK8JetEta)[ijet]) <= 3.0);
    }

    if(jetWP == "tight"){ // tight
    	jetID = ( ( (*AK8JetNHF)[ijet] < 0.90 && (*AK8JetNEF)[ijet] < 0.90 &&  (*AK8Jetnconstituents)[ijet] > 1 ) &&
		  ( (fabs((*AK8JetEta)[ijet]) <= 2.4 &&  (*AK8JetCHF)[ijet] > 0.0 &&  (*AK8JetNCH)[ijet] > 0 && (*AK8JetCEF)[ijet] < 0.99 ) || fabs((*AK8JetEta)[ijet]) > 2.4) && 
		  fabs((*AK8JetEta)[ijet]) <= 3.0);
    }

    if(jetWP == "tightLepVeto"){ //tightLepVeto
       jetID = ( ( (*AK8JetNHF)[ijet] < 0.90 && (*AK8JetNEF)[ijet] < 0.90 &&  (*AK8Jetnconstituents)[ijet] > 1 && (*AK8JetMUF)[ijet] < 0.8) &&
	     ( (fabs((*AK8JetEta)[ijet]) <= 2.4 &&  (*AK8JetCHF)[ijet] > 0.0 &&  (*AK8JetNCH)[ijet] > 0 && (*AK8JetCEF)[ijet] < 0.99 ) || fabs((*AK8JetEta)[ijet]) > 2.4) && 
	     fabs((*AK8JetEta)[ijet]) <= 3.0);
    }

    return jetID;
}

void PhoJet_Analyzer_MC::LumiReWeighting(){


    TFile *fData = TFile::Open("/uscms_data/d3/varun/13TeV/Qstar/PileUpReweighting/PromptReco/hist_DataPU.root");
    h_dataPU = (TH1F*)fData->Get("pileup");
    h_dataPU->Scale(1.0/h_dataPU->Integral());

    TFile *fMC = TFile::Open("/uscms_data/d3/varun/13TeV/Qstar/PileUpReweighting/MC/SigBkg_MC_puHist.root");
    h_mcPU = (TH1F*)fMC->Get("h_mcPileUp50"); 
    h_mcPU->Scale(1.0/h_mcPU->Integral());

    if( h_dataPU->GetNbinsX() != h_mcPU->GetNbinsX() ){
	std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
	return;
    }

    h_weights = (TH1F*) h_dataPU->Clone("h_weights");
    
    h_weights->Divide(h_mcPU);

}

Double_t PhoJet_Analyzer_MC::puweight(Float_t npv){
    Int_t bin = h_weights->GetXaxis()->FindBin( npv );
    return h_weights->GetBinContent( bin );
}





#endif // #ifdef PhoJet_Analyzer_MC_cxx
EOF


###Now compilation
     
g++ -Wno-deprecated PhoJet_Analyzer_MC.C -o ${FileNameTag}_${p}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`
     
echo "--------------------Submitting Job for ${FileNameTag}_${p} Files"   
echo "------------------------Submitting Job #-> ${p}  for  ${sf} to ${ef} Files  -----Total = ${Tot}"

./runCondorFiles.csh ${FileNameTag}_${p} ${FileNameTag}_dataset.txt

##change for next file
@ sf = ${sf} + ${r}
@ ef = ${ef} + ${r}


rm PhoJet_Analyzer_MC.C
rm PhoJet_Analyzer_MC.h


end ## while loop
     echo " ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo "
end  ## foreach i loop
     echo " <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> "
     echo "  "
end ## runCase foreach

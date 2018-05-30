#!/bin/tcsh

##Blocks of code for splitting are enclosed in ##********

setenv pwd $PWD
setenv eosDir /eos/uscms/store/user/rocky86
setenv lpcqstarDir /eos/uscms/store/user/lpcqstar
setenv leptonjetsDir /eos/uscms/store/user/leptonjets

#Run2016BCDEFG_PromptReco for JetHT dataset
setenv DatanTuplesDir root://cmseos.fnal.gov//eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/reMiniAOD03Feb2017/

foreach i (singlePhotonRun2016Bv2  singlePhotonRun2016Cv1  singlePhotonRun2016Dv1  singlePhotonRun2016Ev1  singlePhotonRun2016Fv1  singlePhotonRun2016Gv1  singlePhotonRun2016Hv2  singlePhotonRun2016Hv3)

setenv sourceDir ${DatanTuplesDir}/${i}/

##*************************************************
if ( -f Dataset_${i}.txt ) then
echo "++++++++++++++ Deleting Dataset_${i}.txt ++++++++++++++"
rm Dataset_${i}.txt
endif

eosls ${sourceDir} > Dataset_${i}.txt
setenv dataset Dataset_${i}.txt

##Total no. of files per job
set r = 50
echo "Max Files per Job = ${r}"
set sf = 1
set maxf = ${r}

if( ${i} == "singlePhotonRun2016Hv3" ) then
set Total_files = `wc -l ${dataset} | cut -c1-3`
else
set Total_files = `wc -l ${dataset} | cut -c1-4`
endif

echo "Total number of files in Dataset ${i} are ${Total_files}"

#set Total_jobs = `((((${Total_files}/${r}))+1)`
#echo "Total Jobs to be submitted for ${i} = ${Total_jobs}"

while ( ${sf} <= ${Total_files} )

if ( ${maxf} <= ${Total_files} ) then
setenv filenameTag ${i}_${sf}_${maxf}
else
setenv filenameTag ${i}_${sf}_${Total_files}
endif
##**************************************************

setenv destinationDir ${sourceDir}

if ( -f PostAnalyzer_Data.C ) then
echo "++++++++++++++ Deleting PostAnalyzer_Data.C ++++++++++++++"
rm PostAnalyzer_Data.C
endif

if ( -f PostAnalyzer_Data.h ) then
echo "++++++++++++++ Deleting PostAnalyzer_Data.h ++++++++++++++"
rm PostAnalyzer_Data.h
endif

echo "Filename = ${filenameTag}"
echo "Source Dir = ${sourceDir}"

########## Making PostAnalyzer_Data.C ############
cat > PostAnalyzer_Data.C <<EOF
#define PostAnalyzer_Data_cxx
#include "PostAnalyzer_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PostAnalyzer_Data::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PostAnalyzer_Data.C
//      root> PostAnalyzer_Data t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Initializing various parameters here
   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 190.0; // GeV
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 190.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.8;

   Cut_GJInvtMass = 695.0;

   Cut_PhId = "loose";
   Cut_JetId = "tight";

  //List of single photon triggers in data (Listed here just for information, not used anywhere in the study)
   SinglePhotonTriggers.push_back("HLT_Photon75_v3"); //prescaled
   SinglePhotonTriggers.push_back("HLT_Photon90_v3"); //prescaled
   SinglePhotonTriggers.push_back("HLT_Photon120_v3"); //prescaled
   SinglePhotonTriggers.push_back("HLT_Photon165_HE10_v3"); //unprescaled
   SinglePhotonTriggers.push_back("HLT_Photon175_v3"); //unprescaled
   SinglePhotonTriggers.push_back("HLT_Photon250_NoHE_v2"); //unprescaled
   SinglePhotonTriggers.push_back("HLT_Photon300_NoHE_v2"); //unprescaled
   SinglePhotonTriggers.push_back("HLT_Photon500_v1"); //unprescaled
   SinglePhotonTriggers.push_back("HLT_Photon600_v1"); //unprescaled

   //triggers used in the study (efficiency of these three triggers is same as using all triggers)
   vector<ULong64_t> HLT;
   HLT.clear();
   HLT.push_back(HLT_Photon165_HE10_v);

   //triggers for the denominator of trigger turn on
   vector<ULong64_t> HLT_deno;
   HLT_deno.clear();
   HLT_deno.push_back(HLT_Photon75_v);
   HLT_deno.push_back(HLT_Photon90_v);
   HLT_deno.push_back(HLT_Photon120_v);

   //Uncomment this in script
   //Define Output file here
   //TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   //file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   file = new TFile(OutputFile+".root", "RECREATE");

   //Define Histograms here
   BookHistograms();

   //Defining CSVv2 bTag Operating Point (LOOSE, MEDIUM, TIGHT OR RESHAPING)  
   std::string CSV_WP = "M"; // required for Tagger (L, M or T)              

   //*********************************************************************************************************//
   //Get Event No., Lumi No. and Run no. for high invariant mass events
   Double_t htmass    = 4000.0;
   Int_t    RunNo     = 0;
   Long64_t EvtNo     = 0;
   Int_t    LumiNo    = 0;
   Double_t htPhoPt   = 0.0;
   Double_t htbJetPt  = 0.0;
   Double_t htPhoEta  = 0.0;
   Double_t htbJetEta = 0.0;
   Double_t htPhoPhi  = 0.0;
   Double_t htbJetPhi = 0.0;
   //********************************************************************************************************//

   //Event For loop starts from here
   Long64_t nentries = fChain->GetEntries();
   cout << "<Total entries" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "Analyzing entry:" << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      PC = -1;
      JC = -1;
      GoodVertex = 0;

      //Initially making all bools to false
      Pass_HLT = false;
      HasPrimaryVtx = false;
      Pass_PhoPt = false;
      Pass_PhoEtaEB = false;
      Pass_JetPt = false;
      Pass_JetEta = false;
      //No dphi and deta cut
      Pass_GJdPhi = true;
      Pass_GJdEta = true;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      GoodVertex = nGoodVtx;
      if(GoodVertex > 0) HasPrimaryVtx = true;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons(Cut_PhId); //All photons passing loose id, residual spikes and pt > 30

      GoodIsoBarrelPhotons.clear();
      if(GoodIsoPhotons.size() != 0){
	for(Int_t ip = 0; ip <  GoodIsoPhotons.size(); ip++){
	  if( (*phoEt)[GoodIsoPhotons[ip]] > Cut_Photon_pt && fabs((*phoSCEta)[GoodIsoPhotons[ip]]) < Cut_Photon_eta){
	  GoodIsoBarrelPhotons.push_back(GoodIsoPhotons[ip]);
	  }
	}
      }

      if(GoodIsoBarrelPhotons.size() != 0) PC = GoodIsoBarrelPhotons[0];

      GoodIsoJets.clear();
      GoodIsoJets = GoodJets(Cut_JetId); // All jets passing dR > 0.5, tight jet id and pt > 30.0 
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) < Cut_Jet_eta);
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      //      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);
      
      h_CutFlow_qstar->Fill(0.5);
      h_CutFlow_bstar->Fill(0.5);

      //------------------------------------------------------------
      //Photon distributions noCut
      if(PC > -1 && JC > -1){
	h_PhotonPt[0]               ->Fill((*phoEt)[0]);
	h_PhotonCalibPt[0]          ->Fill((*phoCalibEt)[0]);
	h_PhotonEta[0]              ->Fill((*phoSCEta)[0]);
	h_PhotonPhi[0]              ->Fill((*phoSCPhi)[0]);
	h_Photon_SigmaIEtaIEta[0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[0]);
	h_Photon_R9[0]              ->Fill((*phoR9)[0]);
	h_Photon_HoverE[0]          ->Fill((*phoHoverE)[0]);
	h_Photon_EleVeto[0]         ->Fill((*phoEleVeto)[0]);
	h_Photon_CorrPFChIso[0]     ->Fill(TMath::Max(((*phoPFChIso)[0] - rho*EAChargedHadrons((*phoSCEta)[0])), 0.0));
	h_Photon_CorrPFNeuIso[0]    ->Fill(TMath::Max(((*phoPFNeuIso)[0] - rho*EANeutralHadrons((*phoSCEta)[0])), 0.0));
	h_Photon_CorrPFPhoIso[0]    ->Fill(TMath::Max(((*phoPFPhoIso)[0] - rho*EAPhotons((*phoSCEta)[0])), 0.0));
		      
      //Jet distributions noCut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
	h_bJetPt[0]                 ->Fill((*jetPt)[0]);
	h_bJetEta[0]                ->Fill((*jetEta)[0]);
	h_bJetPhi[0]                ->Fill((*jetPhi)[0]);
	h_bJet_Mt[0]                ->Fill((*jetMt)[0]);
	h_bJet_area[0]              ->Fill((*jetArea)[0]);
	h_bJet_Mass[0]              ->Fill((*jetMass)[0]);
	h_bJet_NHEF[0]              ->Fill((*jetNHF)[0]);
	h_bJet_NEEF[0]              ->Fill((*jetNEF)[0]);
        h_bJet_NConst[0]            ->Fill((*jetNConstituents)[0]);
	h_bJet_CHEF[0]              ->Fill((*jetCHF)[0]);
	h_bJet_ChMult[0]            ->Fill((*jetNCH)[0]);
	h_bJet_CEEF[0]              ->Fill((*jetCEF)[0]);
	h_bJet_MUF[0]               ->Fill((*jetMUF)[0]);
	h_bJet_NNP[0]               ->Fill((*jetNNP)[0]);
		      
      //Photon+Jet distributions noCut
	h_GbJetInvtMass_VarBin[0]   ->Fill(GetInvtMass(0, 0));
	h_GbJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(0, 0));
	h_GbJet_dEta[0]             ->Fill(GetdEta((*phoSCEta)[0], (*jetEta)[0]));
	h_GbJet_dPhi[0]             ->Fill(GetdPhi((*phoSCPhi)[0], (*jetPhi)[0]));
	h_GbJet_dR[0]               ->Fill(GetdR((*phoSCEta)[0], (*jetEta)[0], (*phoSCPhi)[0], (*jetPhi)[0]));
	h_cosThetaStar[0]           ->Fill(GetCosThetaStar((*phoSCEta)[0], (*jetEta)[0]));
		      
      //PFMet distributions for noCut
	h_PFMet[0]                  ->Fill(pfMET);
	h_SumPFMet[0]               ->Fill(pfMETsumEt);
	h_MetBySumMET[0]            ->Fill(pfMET/pfMETsumEt);
	h_PFMetVsGJmass[0]          ->Fill(GetInvtMass(0, 0), pfMET);
	h_PFMetOverSumEtVsGJmass[0] ->Fill(GetInvtMass(0, 0), pfMET/pfMETsumEt);
		      
      //Photon vs Jet dist for noCut
	h_PhPt_vs_bJetPt[0]         ->Fill((*phoEt)[0], (*jetPt)[0]);
	h_PhEta_vs_bJetEta[0]       ->Fill((*phoSCEta)[0], (*jetEta)[0]);
		      
      //CSVv2 discriminator distributions for noCut
	h_CSVv2Dist[0]              ->Fill((*jetCSV2BJetTags)[0]);
	h_CSVv2_vs_bJetPt[0]        ->Fill((*jetPt)[0], (*jetCSV2BJetTags)[0]);
	h_CSVv2_vs_bJetEta[0]       ->Fill((*jetEta)[0], (*jetCSV2BJetTags)[0]);

      //Primary vertex and number of photon and jets for noCut
	h_goodPV[0]                      ->Fill(GoodVertex);
	h_nIsoPhotons[0]                 ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
        h_nGoodPhotons[0]                ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of isolated photons with pt > cut and eta < cut 
        for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
	  h_IsoPhotonIdxVsPt[0]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
	}
	for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
	  h_GoodPhotonIdxVsPt[0]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
	}				     
	h_nJets[0]                       ->Fill(GoodIsoJets.size());
	for(int ij = 0; ij < GoodIsoJets.size(); ij++){
	  h_JetIdxVsPt[0]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
	}
      }
      //------------------------------------------------------------

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	h_CutFlow_qstar->Fill(1.5);
         
	if(HasPrimaryVtx){ // && metFilters == 0){
	  h_CutFlow_qstar->Fill(2.5);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_qstar->Fill(3.5);
	  
	    if(PC > -1){
	      h_CutFlow_qstar->Fill(4.5);
	      
	      if(JC > -1){
		h_CutFlow_qstar->Fill(5.5);
		  
		if(Pass_JetPt){
		  h_CutFlow_qstar->Fill(6.5);

		  if(Pass_JetEta){
		    h_CutFlow_qstar->Fill(7.5);
		  
 		    if(Pass_GJdPhi){
		      h_CutFlow_qstar->Fill(8.5);

		      if(Pass_GJdEta){
			h_CutFlow_qstar->Fill(9.5);

			//----------------------------------------------------------
			//[1]
			//Photon Distributions noBTag_noMasscut
			h_PhotonPt[1]               ->Fill((*phoEt)[PC]);
                        h_PhotonCalibPt[1]          ->Fill((*phoCalibEt)[PC]);
			h_PhotonEta[1]              ->Fill((*phoSCEta)[PC]);
			h_PhotonPhi[1]              ->Fill((*phoSCPhi)[PC]);
			h_Photon_SigmaIEtaIEta[1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			h_Photon_R9[1]              ->Fill((*phoR9)[PC]);
			h_Photon_HoverE[1]          ->Fill((*phoHoverE)[PC]);
			h_Photon_EleVeto[1]         ->Fill((*phoEleVeto)[PC]);
			h_Photon_CorrPFChIso[1]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			h_Photon_CorrPFNeuIso[1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			h_Photon_CorrPFPhoIso[1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			//Jet Distributions noBTag_noMasscut 
			h_bJetPt[1]                 ->Fill((*jetPt)[JC]);
			h_bJetEta[1]                ->Fill((*jetEta)[JC]);
			h_bJetPhi[1]                ->Fill((*jetPhi)[JC]);
			h_bJet_Mt[1]                ->Fill((*jetMt)[JC]);
			h_bJet_area[1]              ->Fill((*jetArea)[JC]);
			h_bJet_Mass[1]              ->Fill((*jetMass)[JC]);
			h_bJet_NHEF[1]              ->Fill((*jetNHF)[JC]);
			h_bJet_NEEF[1]              ->Fill((*jetNEF)[JC]);
			h_bJet_NConst[1]            ->Fill((*jetNConstituents)[JC]);
			h_bJet_CHEF[1]              ->Fill((*jetCHF)[JC]);
			h_bJet_ChMult[1]            ->Fill((*jetNCH)[JC]);
			h_bJet_CEEF[1]              ->Fill((*jetCEF)[JC]);
			h_bJet_MUF[1]               ->Fill((*jetMUF)[JC]);
			h_bJet_NNP[1]               ->Fill((*jetNNP)[JC]);

			//Photon+Jet Distributions noBTag_noMasscut
			h_GbJetInvtMass_VarBin[1]   ->Fill(GetInvtMass(PC, JC));
			h_GbJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC));
			h_GbJet_dEta[1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			h_GbJet_dPhi[1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			h_GbJet_dR[1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			h_cosThetaStar[1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			//PFMet Distributions for noBTag_noMasscut
			h_PFMet[1]                  ->Fill(pfMET);
			h_SumPFMet[1]               ->Fill(pfMETsumEt);
                        h_MetBySumMET[1]            ->Fill(pfMET/pfMETsumEt);
			h_PFMetVsGJmass[1]          ->Fill(GetInvtMass(PC, JC), pfMET);
                        h_PFMetOverSumEtVsGJmass[1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			//Photon vs Jet Dists for noBTag_noMasscut
			h_PhPt_vs_bJetPt[1]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			h_PhEta_vs_bJetEta[1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			//CSVv2 discriminator Distributions for noBTag_noMasscut
			h_CSVv2Dist[1]              ->Fill((*jetCSV2BJetTags)[JC]);
			h_CSVv2_vs_bJetPt[1]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			h_CSVv2_vs_bJetEta[1]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			//Primary vertex and number of photon and jets for noBTag_noMasscut
			h_goodPV[1]                    ->Fill(GoodVertex);
			h_nIsoPhotons[1]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			h_nGoodPhotons[1]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			  h_IsoPhotonIdxVsPt[1]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			}
			for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			  h_GoodPhotonIdxVsPt[1]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			}
			h_nJets[1]                     ->Fill(GoodIsoJets.size());
			for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			  h_JetIdxVsPt[1]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			}
			//----------------------------------------------------------

			//Photon and Jet index for noBTag_noMasscut only
			h_PC                          ->Fill(PC);
			h_JC                          ->Fill(JC);
		      
			if(Pass_GJInvtMass){
			  h_CutFlow_qstar->Fill(10.5);
			  
			  //----------------------------------------------------------
			  //[2]
			  //Photon Distributions noBTag_Masscut
			  h_PhotonPt[2]               ->Fill((*phoEt)[PC]);
			  h_PhotonCalibPt[2]          ->Fill((*phoCalibEt)[PC]);
			  h_PhotonEta[2]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[2]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[2]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[2]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[2]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[2]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[2]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFNeuIso[2]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[2]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			  //Jet Distributions noBTag_Masscut 
			  h_bJetPt[2]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[2]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[2]                ->Fill((*jetPhi)[JC]);
			  h_bJet_Mt[2]                ->Fill((*jetMt)[JC]);
			  h_bJet_area[2]              ->Fill((*jetArea)[JC]);
			  h_bJet_Mass[2]              ->Fill((*jetMass)[JC]);
			  h_bJet_NHEF[2]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[2]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[2]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[2]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[2]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[2]              ->Fill((*jetCEF)[JC]);
			  h_bJet_MUF[2]               ->Fill((*jetMUF)[JC]);
			  h_bJet_NNP[2]               ->Fill((*jetNNP)[JC]);

			  //Photon+Jet Distributions noBTag_Masscut
			  h_GbJetInvtMass_VarBin[2]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[2]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[2]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[2]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[2]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[2]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			  //PFMet Distributions for noBTag_Masscut
			  h_PFMet[2]                  ->Fill(pfMET);
			  h_SumPFMet[2]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[2]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[2]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[2] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			  //Photon vs Jet Dists for noBTag_Masscut
			  h_PhPt_vs_bJetPt[2]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[2]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			  //CSVv2 discriminator Distributions for noBTag_Masscut
			  h_CSVv2Dist[2]              ->Fill((*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[2]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[2]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			  //Primary vertex and number of photon and jets for noBTag_Masscut
			  h_goodPV[2]                    ->Fill(GoodVertex);
			  h_nIsoPhotons[2]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			  h_nGoodPhotons[2]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[2]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[2]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			  }
			  h_nJets[2]                     ->Fill(GoodIsoJets.size());
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[2]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			  }
			  //----------------------------------------------------------
			
			}//if(Pass_GJInvtMass)			
		      }//if(Pass_GJdEta)
		    }//if(Pass_GJdPhi)
		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//PC > -1
	  }//PhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     BSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	h_CutFlow_bstar->Fill(1.5);
         
	if(HasPrimaryVtx && metFilters == 0){
	  h_CutFlow_bstar->Fill(2.5);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_bstar->Fill(3.5);
	  
	    if(PC > -1){
	      h_CutFlow_bstar->Fill(4.5);
	      
	      if(JC > -1){
		h_CutFlow_bstar->Fill(5.5);
		  
		if(Pass_JetPt){
		  h_CutFlow_bstar->Fill(6.5);

		  if(Pass_JetEta){
		    h_CutFlow_bstar->Fill(7.5);
		   
		    if(Pass_CSVv2bTag){
		      h_CutFlow_bstar->Fill(8.5);

		      if(Pass_GJdPhi){
			h_CutFlow_bstar->Fill(9.5);

			if(Pass_GJdEta){
			  h_CutFlow_bstar->Fill(10.5);
			      
			  //----------------------------------------------------------
			  //[3]
			  //Photon Distributions 1BTag_noMasscut
			  h_PhotonPt[3]               ->Fill((*phoEt)[PC]);
			  h_PhotonCalibPt[3]          ->Fill((*phoCalibEt)[PC]);
			  h_PhotonEta[3]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[3]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[3]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[3]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[3]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[3]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[3]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFNeuIso[3]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[3]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		     
			  //Jet Distributions 1BTag_noMasscut 
			  h_bJetPt[3]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[3]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[3]                ->Fill((*jetPhi)[JC]);
			  h_bJet_Mt[3]                ->Fill((*jetMt)[JC]);
			  h_bJet_area[3]              ->Fill((*jetArea)[JC]);
			  h_bJet_Mass[3]              ->Fill((*jetMass)[JC]);
			  h_bJet_NHEF[3]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[3]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[3]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[3]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[3]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[3]              ->Fill((*jetCEF)[JC]);
			  h_bJet_MUF[3]               ->Fill((*jetMUF)[JC]);
			  h_bJet_NNP[3]               ->Fill((*jetNNP)[JC]);

			  //Photon+Jet Distributions 1BTag_noMasscut
			  h_GbJetInvtMass_VarBin[3]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[3]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[3]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[3]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[3]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[3]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			  //PFMet Distributions for 1BTag_noMasscut
			  h_PFMet[3]                  ->Fill(pfMET);
			  h_SumPFMet[3]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[3]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[3]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[3] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			  //Photon vs Jet Dists for 1BTag_noMasscut
			  h_PhPt_vs_bJetPt[3]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[3]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			  //CSVv2 discriminator Distributions for 1BTag_noMasscut
			  h_CSVv2Dist[3]              ->Fill((*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[3]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[3]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			  //Primary vertex and number of photon and jets for 1BTag_noMasscut
			  h_goodPV[3]                    ->Fill(GoodVertex);
			  h_nIsoPhotons[3]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			  h_nGoodPhotons[3]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[3]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[3]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			  }
			  h_nJets[3]                     ->Fill(GoodIsoJets.size());
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[3]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			  }
			  //----------------------------------------------------------

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(11.5);
				
			    //----------------------------------------------------------
			    //[4]
			    //Photon Distributions 1BTag_Masscut
			    h_PhotonPt[4]               ->Fill((*phoEt)[PC]);
			    h_PhotonCalibPt[4]          ->Fill((*phoCalibEt)[PC]);
			    h_PhotonEta[4]              ->Fill((*phoSCEta)[PC]);
			    h_PhotonPhi[4]              ->Fill((*phoSCPhi)[PC]);
			    h_Photon_SigmaIEtaIEta[4]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			    h_Photon_R9[4]              ->Fill((*phoR9)[PC]);
			    h_Photon_HoverE[4]          ->Fill((*phoHoverE)[PC]);
			    h_Photon_EleVeto[4]         ->Fill((*phoEleVeto)[PC]);
			    h_Photon_CorrPFChIso[4]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFNeuIso[4]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFPhoIso[4]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			    //Jet Distributions 1BTag_Masscut 
			    h_bJetPt[4]                 ->Fill((*jetPt)[JC]);
			    h_bJetEta[4]                ->Fill((*jetEta)[JC]);
			    h_bJetPhi[4]                ->Fill((*jetPhi)[JC]);
			    h_bJet_Mt[4]                ->Fill((*jetMt)[JC]);
			    h_bJet_area[4]              ->Fill((*jetArea)[JC]);
			    h_bJet_Mass[4]              ->Fill((*jetMass)[JC]);
			    h_bJet_NHEF[4]              ->Fill((*jetNHF)[JC]);
			    h_bJet_NEEF[4]              ->Fill((*jetNEF)[JC]);
			    h_bJet_NConst[4]            ->Fill((*jetNConstituents)[JC]);
			    h_bJet_CHEF[4]              ->Fill((*jetCHF)[JC]);
			    h_bJet_ChMult[4]            ->Fill((*jetNCH)[JC]);
			    h_bJet_CEEF[4]              ->Fill((*jetCEF)[JC]);
			    h_bJet_MUF[4]               ->Fill((*jetMUF)[JC]);
			    h_bJet_NNP[4]               ->Fill((*jetNNP)[JC]);

			    //Photon+Jet Distributions 1BTag_Masscut
			    h_GbJetInvtMass_VarBin[4]   ->Fill(GetInvtMass(PC, JC));
			    h_GbJetInvtMass_UnitBin[4]  ->Fill(GetInvtMass(PC, JC));
			    h_GbJet_dEta[4]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			    h_GbJet_dPhi[4]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_GbJet_dR[4]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_cosThetaStar[4]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			    //PFMet Distributions for 1BTag_Masscut
			    h_PFMet[4]                  ->Fill(pfMET);
			    h_SumPFMet[4]               ->Fill(pfMETsumEt);
			    h_MetBySumMET[4]            ->Fill(pfMET/pfMETsumEt);
			    h_PFMetVsGJmass[4]          ->Fill(GetInvtMass(PC, JC), pfMET);
			    h_PFMetOverSumEtVsGJmass[4] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			    //Photon vs Jet Dists for 1BTag_Masscut
			    h_PhPt_vs_bJetPt[4]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			    h_PhEta_vs_bJetEta[4]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			    //CSVv2 discriminator Distributions for 1BTag_Masscut
			    h_CSVv2Dist[4]              ->Fill((*jetCSV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetPt[4]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetEta[4]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			    //Primary vertex and number of photon and jets for 1BTag_Masscut
			    h_goodPV[4]                    ->Fill(GoodVertex);
			    h_nIsoPhotons[4]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			    h_nGoodPhotons[4]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[4]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[4]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			    }
			    h_nJets[4]                     ->Fill(GoodIsoJets.size());
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[4]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			    }
			    //----------------------------------------------------------

			  }//if(Pass_GJInvtMass) inside if(Pass_CSVv2bTag)
			}//if(Pass_GJdEta) inside if(Pass_CSVv2bTag)
		      }//if(Pass_GJdPhi) inside if(Pass_CSVv2bTag)		       
		      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    }//if(Pass_CSVv2bTag)		     
		    else{
		      h_CutFlow_bstar->Fill(12.5);

		      if(Pass_GJdPhi){
			h_CutFlow_bstar->Fill(13.5);

			if(Pass_GJdEta){
			  h_CutFlow_bstar->Fill(14.5);
			  
			  //----------------------------------------------------------
			  //[5]
			  //Photon Distributions 0BTag_noMasscut
			  h_PhotonPt[5]               ->Fill((*phoEt)[PC]);
			  h_PhotonCalibPt[5]          ->Fill((*phoCalibEt)[PC]);
			  h_PhotonEta[5]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[5]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[5]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[5]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[5]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[5]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[5]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFNeuIso[5]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[5]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			  //Jet Distributions 0BTag_noMasscut 
			  h_bJetPt[5]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[5]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[5]                ->Fill((*jetPhi)[JC]);
			  h_bJet_Mt[5]                ->Fill((*jetMt)[JC]);
			  h_bJet_area[5]              ->Fill((*jetArea)[JC]);
			  h_bJet_Mass[5]              ->Fill((*jetMass)[JC]);
			  h_bJet_NHEF[5]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[5]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[5]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[5]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[5]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[5]              ->Fill((*jetCEF)[JC]);
			  h_bJet_MUF[5]               ->Fill((*jetMUF)[JC]);
			  h_bJet_NNP[5]               ->Fill((*jetNNP)[JC]);

			  //Photon+Jet Distributions 0BTag_noMasscut
			  h_GbJetInvtMass_VarBin[5]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[5]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[5]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[5]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[5]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[5]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			  //PFMet Distributions for 0BTag_noMasscut
			  h_PFMet[5]                  ->Fill(pfMET);
			  h_SumPFMet[5]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[5]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[5]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[5] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			  //Photon vs Jet Dists for 0BTag_noMasscut
			  h_PhPt_vs_bJetPt[5]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[5]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			  //CSVv2 discriminator Distributions for 0BTag_noMasscut
			  h_CSVv2Dist[5]              ->Fill((*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[5]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[5]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			  //Primary vertex and number of photon and jets for 0BTag_noMasscut
			  h_goodPV[5]                    ->Fill(GoodVertex);
			  h_nIsoPhotons[5]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			  h_nGoodPhotons[5]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[5]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[5]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			  }
			  h_nJets[5]                     ->Fill(GoodIsoJets.size());
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[5]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			  }
			  //----------------------------------------------------------

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(15.5);
			    
			    //----------------------------------------------------------
			    //[6]
			    //Photon Distributions 0BTag_Masscut
			    h_PhotonPt[6]               ->Fill((*phoEt)[PC]);
			    h_PhotonCalibPt[6]          ->Fill((*phoCalibEt)[PC]);
			    h_PhotonEta[6]              ->Fill((*phoSCEta)[PC]);
			    h_PhotonPhi[6]              ->Fill((*phoSCPhi)[PC]);
			    h_Photon_SigmaIEtaIEta[6]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			    h_Photon_R9[6]              ->Fill((*phoR9)[PC]);
			    h_Photon_HoverE[6]          ->Fill((*phoHoverE)[PC]);
			    h_Photon_EleVeto[6]         ->Fill((*phoEleVeto)[PC]);
			    h_Photon_CorrPFChIso[6]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFNeuIso[6]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFPhoIso[6]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			    //Jet Distributions 0BTag_Masscut 
			    h_bJetPt[6]                 ->Fill((*jetPt)[JC]);
			    h_bJetEta[6]                ->Fill((*jetEta)[JC]);
			    h_bJetPhi[6]                ->Fill((*jetPhi)[JC]);
			    h_bJet_Mt[6]                ->Fill((*jetMt)[JC]);
			    h_bJet_area[6]              ->Fill((*jetArea)[JC]);
			    h_bJet_Mass[6]              ->Fill((*jetMass)[JC]);
			    h_bJet_NHEF[6]              ->Fill((*jetNHF)[JC]);
			    h_bJet_NEEF[6]              ->Fill((*jetNEF)[JC]);
			    h_bJet_NConst[6]            ->Fill((*jetNConstituents)[JC]);
			    h_bJet_CHEF[6]              ->Fill((*jetCHF)[JC]);
			    h_bJet_ChMult[6]            ->Fill((*jetNCH)[JC]);
			    h_bJet_CEEF[6]              ->Fill((*jetCEF)[JC]);
			    h_bJet_MUF[6]               ->Fill((*jetMUF)[JC]);
			    h_bJet_NNP[6]               ->Fill((*jetNNP)[JC]);

			    //Photon+Jet Distributions 0BTag_Masscut
			    h_GbJetInvtMass_VarBin[6]   ->Fill(GetInvtMass(PC, JC));
			    h_GbJetInvtMass_UnitBin[6]  ->Fill(GetInvtMass(PC, JC));
			    h_GbJet_dEta[6]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			    h_GbJet_dPhi[6]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_GbJet_dR[6]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_cosThetaStar[6]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			    //PFMet Distributions for 0BTag_Masscut
			    h_PFMet[6]                  ->Fill(pfMET);
			    h_SumPFMet[6]               ->Fill(pfMETsumEt);
			    h_MetBySumMET[6]            ->Fill(pfMET/pfMETsumEt);
			    h_PFMetVsGJmass[6]          ->Fill(GetInvtMass(PC, JC), pfMET);
			    h_PFMetOverSumEtVsGJmass[6] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			    //Photon vs Jet Dists for 0BTag_Masscut
			    h_PhPt_vs_bJetPt[6]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			    h_PhEta_vs_bJetEta[6]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			    //CSVv2 discriminator Distributions for 0BTag_Masscut
			    h_CSVv2Dist[6]              ->Fill((*jetCSV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetPt[6]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetEta[6]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC]);

			    //Primary vertex and number of photon and jets for 0BTag_Masscut
			    h_goodPV[6]                    ->Fill(GoodVertex);
			    h_nIsoPhotons[6]               ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
			    h_nGoodPhotons[6]              ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of iso photons with pt>cut and eta<cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[6]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[6]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
			    }
			    h_nJets[6]                     ->Fill(GoodIsoJets.size());
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[6]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1);
			    }
			    //----------------------------------------------------------

			  }//if(Pass_GJInvtMass) inside else(Pass_CSVv2bTag)
			}//if(Pass_GJdEta) inside else(Pass_CSVv2bTag)
		      }//if(Pass_GJdPhi) inside else(Pass_CSVv2bTag)
		    }//else(Pass_CSVv2bTag)
		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//PC > -1
	  }//PhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      //---------------------------------------------------------------------------------   
      //Trigger Turn-on
      Bool_t PassHLTNum = false;
      Bool_t PassHLTDeno = false;
      Bool_t PassHLT_Pre = false;
      //      cout << PassHLT_Pre << endl;
      ULong64_t HLT_pre = HLT[0];

      PassHLTNum = PassHLT(HLT);
      PassHLTDeno = PassHLT(HLT_deno);
      PassHLT_Pre = PassHLT_Prescale(HLT_pre);

      std::vector<Int_t> TrigPhotons;
      std::vector<Int_t> TrigJets;

      TrigPhotons.clear();
      for(Int_t i = 0; i < nPho; i++){
	if(CutBasedPhotonID(i, "loose") && fabs((*phoSCEta)[i]) < Cut_Photon_eta) TrigPhotons.push_back(i);
      }

      TrigJets.clear();
      for(Int_t j = 0; j < nJet; j++){
	if(TrigPhotons.size() > 0 && JetId(j, "tight") && (*jetEta)[j] < Cut_Jet_eta)
	  if(GetdR((*phoSCEta)[TrigPhotons[0]], (*jetEta)[j], (*phoSCPhi)[TrigPhotons[0]], (*jetPhi)[j]) > 0.5)
	    TrigJets.push_back(j);
      }

      if(TrigPhotons.size() > 0){
        if(PassHLTDeno){
	  h_TrigPhotonPt[0]->Fill((*phoEt)[TrigPhotons[0]]);	
	  if(PassHLTNum && !(PassHLT_Pre)){
	    h_TrigPhotonPt[1]->Fill((*phoEt)[TrigPhotons[0]]);
	  }
	}
      }

      if(TrigPhotons.size() > 0 && TrigJets.size() > 0){
	if((*phoEt)[TrigPhotons[0]] > 190.0 && (*jetPt)[TrigJets[0]] > 60.0){
	  h_TrigGJmass_Jetpt60->Fill(GetInvtMass(TrigPhotons[0], TrigJets[0]));
	  if((*phoEt)[TrigPhotons[0]] > 190.0 && (*jetPt)[TrigJets[0]] > 100){
	    h_TrigGJmass_Jetpt100->Fill(GetInvtMass(TrigPhotons[0], TrigJets[0]));
	  }
	  if((*phoEt)[TrigPhotons[0]] > 190.0 && (*jetPt)[TrigJets[0]] > 120){
	    h_TrigGJmass_Jetpt120->Fill(GetInvtMass(TrigPhotons[0], TrigJets[0]));
	  }
	  if((*phoEt)[TrigPhotons[0]] > 190.0 && (*jetPt)[TrigJets[0]] > 150){
	    h_TrigGJmass_Jetpt150->Fill(GetInvtMass(TrigPhotons[0], TrigJets[0]));
	  }
	  if((*phoEt)[TrigPhotons[0]] > 190.0 && (*jetPt)[TrigJets[0]] > 190){
	    h_TrigGJmass_Jetpt190->Fill(GetInvtMass(TrigPhotons[0], TrigJets[0]));
	  }
	}
      }

      //--------------------------------------------------------------------------------------------------
      //Getting High mass events
      if(Pass_HLT){
	if(HasPrimaryVtx && metFilters == 0){
	  if(GoodIsoPhotons.size() > 0)
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
		    Double_t InvtMass = GetInvtMass(PC, JC);
		    if(InvtMass > htmass ){ 
		      if(Pass_CSVv2bTag){

			cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			cout <<"                                        1 BTAG EVENT                                           "<< endl;
			cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;

			RunNo = run;
			EvtNo = event;
			LumiNo = lumis;
			htPhoPt = (*phoEt)[PC];
			htbJetPt = (*jetPt)[JC];
			htPhoEta = (*phoSCEta)[PC];
			htbJetEta = (*jetEta)[JC];
			htPhoPhi = (*phoSCPhi)[PC];
			htbJetPhi = (*jetPhi)[JC];

		      }else{

			cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			cout <<"                                        0 BTAG EVENT                                           "<< endl;
			cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;

			RunNo = run;
			EvtNo = event;
			LumiNo = lumis;
			htPhoPt = (*phoEt)[PC];
			htbJetPt = (*jetPt)[JC];
			htPhoEta = (*phoSCEta)[PC];
			htbJetEta = (*jetEta)[JC];
			htPhoPhi = (*phoSCPhi)[PC];
			htbJetPhi = (*jetPhi)[JC];

		      }

		      cout <<"                                                                                               "<< endl;
		      cout <<"+++++++++++++++++++++++++++++++++++ EVENT INFO ==> ENTRY = " << jentry << "+++++++++++++++++++++++++++++++"<< endl;
		      cout <<"| InvtMass (GeV/c2) | Ph_Pt (GeV/c) | Jet_Pt (GeV/c) |  Ph_Eta  |  Jet_Eta  |  Ph_Phi  |  Jet_Phi  |  Run_No  |  Event_No  |  Lumi_No  |" << endl;
		      cout << "|   " << InvtMass << "   |   " << htPhoPt << "   |   " << htbJetPt << "   |   " << htPhoEta << "   |   " << htbJetEta << "   |   " << htPhoPhi << "   |   " << htbJetPhi << "   |   " << RunNo << "   |   " << EvtNo << "   |   " << LumiNo << "   |   " << endl;
		      cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
		      cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
		      cout <<"                                                                                               "<< endl;


		    }//if(InvtMass > htmass )
		  }
		}
	      }
	    }
	}
      }
      //--------------------------------------------------------------------------------------------------
      
   }//for jentry
}//Loop() 



EOF



########## Making PostAnalyzer_Data.h ############
cat > PostAnalyzer_Data.h <<EOF
/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 00:35:43 2015 by ROOT version 6.02/05
// from TChain ggNtuplizer/EventTree/
//////////////////////////////////////////////////////////

#ifndef PostAnalyzer_Data_h
#define PostAnalyzer_Data_h

//ROOT include files
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include "TRFIOFile.h"
#include "TKDE.h" 

// Header file for the classes stored in the TTree if any.
#include "vector"

//c++ include files
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <map>
#include <set>

using namespace std;
using namespace ROOT;

class PostAnalyzer_Data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   //-----------------------
   //Variables defined by me
   //-----------------------
   //Outut root file
   TFile *file;

   //Bools
   Bool_t Pass_HLT;
   Bool_t HasPrimaryVtx;
   Bool_t Pass_PhoPt;
   Bool_t Pass_PhoEtaEB;
   Bool_t Pass_JetPt;
   Bool_t Pass_JetEta;
   Bool_t Pass_GJdPhi;
   Bool_t Pass_GJdEta;
   Bool_t Pass_GJInvtMass;
   Bool_t Pass_CSVv2bTag;

   //Ints, Doubles etc.
   Int_t GoodVertex;
   Int_t PC, JC;

   Double_t Cut_Vtx_z; //(cm)
   Double_t Cut_Vtx_ndof;
   Double_t Cut_Vtx_rho; //(cm)

   Double_t Cut_Photon_pt; //(GeV)
   Double_t Cut_Photon_eta;

   Double_t Cut_Jet_pt; //(GeV)
   Double_t Cut_Jet_eta;

   Double_t Cut_GJdPhi;
   Double_t Cut_GJdEta;
   Double_t Cut_GJInvtMass;
   TString Cut_PhId;
   TString Cut_JetId;

   //Photon and Jet Containers
   vector<Int_t> GoodIsoPhotons;
   vector<Int_t> GoodIsoBarrelPhotons;
   vector<Int_t> GoodIsoJets;

   //Trigger Info (Required for data only)
   std::vector<std::string> SinglePhotonTriggers;
   //Unprescaled Triggers (Bits assigned according to ggAnalysis/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc)
   const ULong64_t HLT_Photon175_v = 7;
   const ULong64_t HLT_Photon250_NoHE_v = 8;
   const ULong64_t HLT_Photon300_NoHE_v = 9;
   const ULong64_t HLT_Photon500_v = 10;
   const ULong64_t HLT_Photon600_v = 11;
   const ULong64_t HLT_Photon165_HE10_v = 12;
   //Prescaled Triggers (Bits assigned according to ggAnalysis/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc)
   const ULong64_t HLT_Photon75_v = 4;
   const ULong64_t HLT_Photon90_v = 5;
   const ULong64_t HLT_Photon120_v = 6;

   //*****************************************************************
   //Histograms
   TH1F *h_PhotonPt[7];
   TH1F *h_PhotonCalibPt[7];
   TH1F *h_PhotonEta[7];
   TH1F *h_PhotonPhi[7];
   TH1F *h_Photon_SigmaIEtaIEta[7];
   TH1F *h_Photon_R9[7];
   TH1F *h_Photon_HoverE[7];
   TH1F *h_Photon_EleVeto[7];
   TH1F *h_Photon_CorrPFChIso[7];
   TH1F *h_Photon_CorrPFNeuIso[7];
   TH1F *h_Photon_CorrPFPhoIso[7];

   TH1F *h_bJetPt[7];
   TH1F *h_bJetEta[7];
   TH1F *h_bJetPhi[7];
   TH1F *h_bJet_Mt[7];
   TH1F *h_bJet_area[7];
   TH1F *h_bJet_Mass[7];
   TH1F *h_bJet_NHEF[7]; //Neutral Hadron energy fraction
   TH1F *h_bJet_NEEF[7]; //Neutral EM energy fraction
   TH1F *h_bJet_NConst[7]; //Number of constituents
   TH1F *h_bJet_CHEF[7];  //Charged Hadron energy fraction
   TH1F *h_bJet_ChMult[7]; //Charged Multiplicity
   TH1F *h_bJet_CEEF[7]; //Charged EM energy fraction
   TH1F *h_bJet_MUF[7]; //Muon energy fraction
   TH1F *h_bJet_NNP[7]; //Number of neutral particles

   TH1F *h_GbJetInvtMass_VarBin[7];
   TH1F *h_GbJetInvtMass_UnitBin[7];
   TH1F *h_GbJet_dEta[7];
   TH1F *h_GbJet_dPhi[7];
   TH1F *h_GbJet_dR[7];
   TH1F *h_cosThetaStar[7];

   TH1F *h_PFMet[7];
   TH1F *h_SumPFMet[7];
   TH1F *h_MetBySumMET[7];
   TH2F *h_PFMetVsGJmass[7];
   TH2F *h_PFMetOverSumEtVsGJmass[7];

   TH2F *h_PhPt_vs_bJetPt[7];
   TH2F *h_PhEta_vs_bJetEta[7];

   TH1F *h_CSVv2Dist[7];
   TH2F *h_CSVv2_vs_bJetPt[7];
   TH2F *h_CSVv2_vs_bJetEta[7];

   //Pileup
   TH1F *h_goodPV[7];

   //Number of photons and jets
   TH1F *h_nIsoPhotons[7];
   TH1F *h_nGoodPhotons[7];
   TH2F *h_IsoPhotonIdxVsPt[7];
   TH2F *h_GoodPhotonIdxVsPt[7];
   TH1F *h_nJets[7];
   TH2F *h_JetIdxVsPt[7];

   //Trigger Turn-on
   TH1F *h_TrigPhotonPt[2];
   TH1F *h_TrigGJmass_Jetpt60;
   TH1F *h_TrigGJmass_Jetpt100;
   TH1F *h_TrigGJmass_Jetpt120;
   TH1F *h_TrigGJmass_Jetpt150;
   TH1F *h_TrigGJmass_Jetpt190;

   TH1F *h_PC;
   TH1F *h_JC;

   //Cut Flow
   TH1F *h_CutFlow_qstar;
   TH1F *h_CutFlow_bstar;

   //***********************************************************

   //Variables of TTree::EventTree
   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Int_t           nTrksPV;
   Bool_t          isPVGood;
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
   vector<int>     *phoPrescale;
   Int_t           metFilters;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoPx;
   vector<float>   *phoPy;
   vector<float>   *phoPz;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
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
   vector<float>   *phoE1x3;
   vector<float>   *phoE1x5;
   vector<float>   *phoE2x2;
   vector<float>   *phoE2x5Max;
   vector<float>   *phoE5x5;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phomaxXtalenergyFull5x5;
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
   vector<float>   *phoCITKChIso;
   vector<float>   *phoCITKPhoIso;
   vector<float>   *phoCITKNeuIso;
   vector<float>   *phoIDMVA;
   vector<unsigned int> *phoFiredSingleTrgs;
   vector<unsigned int> *phoFiredDoubleTrgs;
   vector<unsigned int> *phoFiredL1Trgs;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<unsigned short> *phoxtalBits;
   vector<unsigned short> *phoIDbit;
   Int_t           npfPho;
   vector<float>   *pfphoEt;
   vector<float>   *pfphoEta;
   vector<float>   *pfphoPhi;
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
   vector<float>   *eleSIP;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
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
   vector<float>   *eleIDMVA;
   vector<float>   *eleIDMVAHZZ;
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
   vector<float>   *eleGSFChi2;
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
   vector<unsigned int> *eleFiredSingleTrgs;
   vector<unsigned int> *eleFiredDoubleTrgs;
   vector<unsigned int> *eleFiredL1Trgs;
   vector<unsigned short> *eleIDbit;
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
   vector<unsigned short> *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
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
   vector<unsigned int> *muFiredTrgs;
   vector<unsigned int> *muFiredL1Trgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetMass;
   vector<float>   *jetPx;
   vector<float>   *jetPy;
   vector<float>   *jetPz;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetLeadTrackEta;
   vector<float>   *jetLeadTrackPhi;
   vector<int>     *jetLepTrackPID;
   vector<float>   *jetLepTrackPt;
   vector<float>   *jetLepTrackEta;
   vector<float>   *jetLepTrackPhi;
   vector<float>   *jetCSV2BJetTags;
   vector<float>   *jetJetProbabilityBJetTags;
   vector<float>   *jetpfCombinedMVAV2BJetTags;
   vector<bool>    *jetPFLooseId;
   vector<int>     *jetID;
   vector<float>   *jetPUID;
   vector<float>   *jetJECUnc;
   vector<float>   *jetJERSmearing;
   vector<float>   *jetJERSmearingUp;
   vector<float>   *jetJERSmearingDown;
   vector<unsigned int> *jetFiredTrgs;
   vector<float>   *jetCHF; //chargedHadronEnergyFraction
   vector<float>   *jetNHF; //neutralHadronEnergyFraction
   vector<float>   *jetCEF; //chargedEmEnergyFraction
   vector<float>   *jetNEF; //neutralEmEnergyFraction
   vector<int>     *jetNCH; //chargedMultiplicity
   vector<int>     *jetNNP; //NumberofNeutralParticles or neutal Multiplicity
   vector<float>   *jetMUF; //Muon energy fraction
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtxNtrks;
   vector<float>   *jetVtx3DVal;
   vector<float>   *jetVtx3DSig;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConstituents; //Number of constituents
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
   vector<int>     *AK8JetNNP;
   vector<float>   *AK8JetMUF;
   vector<int>     *AK8Jetnconstituents;
   vector<bool>    *AK8JetPFLooseId;
   vector<bool>    *AK8JetPFTightLepVetoId;
   vector<float>   *AK8JetSoftDropMass;
   vector<float>   *AK8JetSoftDropMassCorr;
   vector<float>   *AK8JetPrunedMass;
   vector<float>   *AK8JetPrunedMassCorr;
   vector<float>   *AK8JetpfBoostedDSVBTag;
   vector<float>   *AK8JetDSVnewV4;
   vector<float>   *AK8JetCSV;
   vector<float>   *AK8JetJECUnc;
   vector<float>   *AK8JetL2L3corr;
   vector<float>   *AK8puppiPt;
   vector<float>   *AK8puppiMass;
   vector<float>   *AK8puppiEta;
   vector<float>   *AK8puppiPhi;
   vector<float>   *AK8puppiTau1;
   vector<float>   *AK8puppiTau2;
   vector<float>   *AK8puppiTau3;
   vector<float>   *AK8puppiSDL2L3corr;
   vector<float>   *AK8puppiSDMass;
   vector<float>   *AK8puppiSDMassL2L3Corr;
   vector<int>     *nAK8SDSJ;
   vector<vector<float> > *AK8SDSJPt;
   vector<vector<float> > *AK8SDSJEta;
   vector<vector<float> > *AK8SDSJPhi;
   vector<vector<float> > *AK8SDSJMass;
   vector<vector<float> > *AK8SDSJE;
   vector<vector<int> > *AK8SDSJCharge;
   vector<vector<int> > *AK8SDSJFlavour;
   vector<vector<float> > *AK8SDSJCSV;
   vector<int>     *nAK8puppiSDSJ;
   vector<vector<float> > *AK8puppiSDSJPt;
   vector<vector<float> > *AK8puppiSDSJEta;
   vector<vector<float> > *AK8puppiSDSJPhi;
   vector<vector<float> > *AK8puppiSDSJMass;
   vector<vector<float> > *AK8puppiSDSJE;
   vector<vector<int> > *AK8puppiSDSJCharge;
   vector<vector<int> > *AK8puppiSDSJFlavour;
   vector<vector<float> > *AK8puppiSDSJCSV;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_isPVGood;   //!
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
   TBranch        *b_phoPrescale;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
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
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phomaxXtalenergyFull5x5;   //!
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
   TBranch        *b_phoCITKChIso;   //!
   TBranch        *b_phoCITKPhoIso;   //!
   TBranch        *b_phoCITKNeuIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoxtalBits;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_npfPho;   //!
   TBranch        *b_pfphoEt;   //!
   TBranch        *b_pfphoEta;   //!
   TBranch        *b_pfphoPhi;   //!
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
   TBranch        *b_eleSIP;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
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
   TBranch        *b_eleIDMVA;   //!
   TBranch        *b_eleIDMVAHZZ;   //!
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
   TBranch        *b_eleGSFChi2;   //!
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
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleIDbit;   //!
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
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
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
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetPx;   //!
   TBranch        *b_jetPy;   //!
   TBranch        *b_jetPz;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetLeadTrackEta;   //!
   TBranch        *b_jetLeadTrackPhi;   //!
   TBranch        *b_jetLepTrackPID;   //!
   TBranch        *b_jetLepTrackPt;   //!
   TBranch        *b_jetLepTrackEta;   //!
   TBranch        *b_jetLepTrackPhi;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetpfCombinedMVAV2BJetTags;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetJERSmearing;   //!
   TBranch        *b_jetJERSmearingUp;   //!
   TBranch        *b_jetJERSmearingDown;   //!
   TBranch        *b_jetFiredTrgs;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetNNP;   //!
   TBranch        *b_jetMUF;   //!
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
   TBranch        *b_AK8JetNNP;   //!
   TBranch        *b_AK8JetMUF;   //!
   TBranch        *b_AK8Jetnconstituents;   //!
   TBranch        *b_AK8JetPFLooseId;   //!
   TBranch        *b_AK8JetPFTightLepVetoId;   //!
   TBranch        *b_AK8JetSoftDropMass;   //!
   TBranch        *b_AK8JetSoftDropMassCorr;   //!
   TBranch        *b_AK8JetPrunedMass;   //!
   TBranch        *b_AK8JetPrunedMassCorr;   //!
   TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
   TBranch        *b_AK8JetDSVnewV4;   //!
   TBranch        *b_AK8JetCSV;   //!
   TBranch        *b_AK8JetJECUnc;   //!
   TBranch        *b_AK8JetL2L3corr;   //!
   TBranch        *b_AK8puppiPt;   //!
   TBranch        *b_AK8puppiMass;   //!
   TBranch        *b_AK8puppiEta;   //!
   TBranch        *b_AK8puppiPhi;   //!
   TBranch        *b_AK8puppiTau1;   //!
   TBranch        *b_AK8puppiTau2;   //!
   TBranch        *b_AK8puppiTau3;   //!
   TBranch        *b_AK8puppiSDL2L3corr;   //!
   TBranch        *b_AK8puppiSDMass;   //!
   TBranch        *b_AK8puppiSDMassL2L3Corr;   //!
   TBranch        *b_nAK8SDSJ;   //!
   TBranch        *b_AK8SDSJPt;   //!
   TBranch        *b_AK8SDSJEta;   //!
   TBranch        *b_AK8SDSJPhi;   //!
   TBranch        *b_AK8SDSJMass;   //!
   TBranch        *b_AK8SDSJE;   //!
   TBranch        *b_AK8SDSJCharge;   //!
   TBranch        *b_AK8SDSJFlavour;   //!
   TBranch        *b_AK8SDSJCSV;   //!
   TBranch        *b_nAK8puppiSDSJ;   //!
   TBranch        *b_AK8puppiSDSJPt;   //!
   TBranch        *b_AK8puppiSDSJEta;   //!
   TBranch        *b_AK8puppiSDSJPhi;   //!
   TBranch        *b_AK8puppiSDSJMass;   //!
   TBranch        *b_AK8puppiSDSJE;   //!
   TBranch        *b_AK8puppiSDSJCharge;   //!
   TBranch        *b_AK8puppiSDSJFlavour;   //!
   TBranch        *b_AK8puppiSDSJCSV;   //!

   PostAnalyzer_Data(TTree *tree=0);
   virtual ~PostAnalyzer_Data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //User-Defined functions
   virtual Bool_t   PassHLT(vector<ULong64_t> trigBits);
   virtual Bool_t   PassHLT_Prescale(ULong64_t trigBit);
   //virtual Bool_t   GoodPrimaryVtx(Int_t &GoodVertex);

   virtual Bool_t   CutBasedPhotonID(Int_t ipho, TString phoWP);
   virtual Double_t EAChargedHadrons(Double_t eta);
   virtual Double_t EANeutralHadrons(Double_t eta);
   virtual Double_t EAPhotons(Double_t eta);
   virtual Int_t    FirstGoodPhoton(TString phoWP); 
   virtual vector<Int_t> GoodPhotons(TString phoWP); 

   virtual Bool_t   ResSpikes(Int_t);
   virtual Bool_t   JetId(Int_t iJet, TString jetWP);
   virtual Int_t    FirstGoodJet(TString jetWP);
   virtual vector<Int_t> GoodJets(TString jetWP);

   //For 80X (WP cuts need to change for 76X)
   virtual Bool_t   CSVv2bTag(Int_t ijet, TString WP);

   virtual Double_t GetdEta(Double_t eta1, Double_t eta2);
   virtual Double_t GetdPhi(Double_t phi1, Double_t phi2);
   virtual Double_t GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
   virtual Double_t GetCosThetaStar(Double_t eta1, Double_t eta2);
   virtual Double_t GetInvtMass(Int_t ph, Int_t jet);

   virtual void     BookHistograms();

};

#endif

#ifdef PostAnalyzer_Data_cxx
PostAnalyzer_Data::PostAnalyzer_Data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ggNtuplizer/EventTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ggNtuplizer/EventTree","");
      //      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/ReReco-BCDEFG_PromptReco-H/Run2016B-23Sep2016_ReReco/Data_Run2016B-23Sep2016_ReReco_4151.root/ggNtuplizer/EventTree");
      //chain->Add("/eos/uscms/store/user/rocky86/13TeV/Ntuples/80X/Data/ReReco-BCDEFG_PromptReco-H/Run2016G-23Sep2016_ReReco-v1/Data_Run2016G-23Sep2016_ReReco-v1_1.root/ggNtuplizer/EventTree");
      //      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/PromptReco/Run2016B_v2/Data_Run2016B_v2_1515.root/ggNtuplizer/EventTree");
      //      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/PromptReco/Run2016B_v2/Data_Run2016B_v2_1517.root/ggNtuplizer/EventTree");
      //      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/PromptReco/Run2016B_v2/Data_Run2016B_v2_1521.root/ggNtuplizer/EventTree");

/*
      //Uncomment this part in script

      ///--------------------Use this part while submitting job in one go for a dataset--------------///
      TString main_path = "${sourceDir}";

      TSystemDirectory sourceDir("sysDir",main_path);
      TList* fileList = sourceDir.GetListOfFiles();
      TIter next(fileList);
      TSystemFile* fileName;

      int fileNumber = 1;
      int maxFiles = -1;

      while ((fileName = (TSystemFile*)next()) && fileNumber > maxFiles){
        if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."){continue;}

        TString FullPathInputFile = (main_path+fileName->GetName());

        //cout << FullPathInputFile << endl;

        chain->Add(FullPathInputFile+"/ggNtuplizer/EventTree");

        fileNumber++;

      }

      cout << "Total files in this set = " << fileNumber - 1 << endl;
      ///-------------------------------------------------------------------------------------------------///
*/
      ///-----------Only change Tstring part with this part for job submission in parts for a dataset-----///
      ifstream Dataset;
      Dataset.open("${pwd}//${dataset}", ifstream::in);
      char datafilename[500];

      if(Dataset) cout << "${dataset} Opened" << endl;

      int a = 0;
      for(Int_t i = 1; i <= ${maxf} && i <= ${Total_files}; i++){      
	Dataset >> datafilename;   ////dataset >> datafilename always starts from the 1st line if the file is just opened and will start from the 
                                   //// next to last line if already opened.
	string fname(datafilename);
	string main_path = "${sourceDir}";
	string myevt = "/ggNtuplizer/EventTree";

        if(i >= ${sf} && i <= ${Total_files}){
          chain->Add((main_path+fname+myevt).c_str());
          a++;
        }

      //cout << (FullPathInputFile+myevt).c_str() << endl;

      }
      cout << "Total Files in this job are = " << a << endl;
      ///------------------------------------------------------------------------------------------------///

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PostAnalyzer_Data::~PostAnalyzer_Data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   file->cd();
   file->Write();
   file->Close();
}

Int_t PostAnalyzer_Data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PostAnalyzer_Data::LoadTree(Long64_t entry)
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

void PostAnalyzer_Data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   phoPrescale = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
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
   phoE1x3 = 0;
   phoE1x5 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE5x5 = 0;
   phoESEffSigmaRR = 0;
   phomaxXtalenergyFull5x5 = 0;
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
   phoCITKChIso = 0;
   phoCITKPhoIso = 0;
   phoCITKNeuIso = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoxtalBits = 0;
   phoIDbit = 0;
   pfphoEt = 0;
   pfphoEta = 0;
   pfphoPhi = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
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
   eleIDMVA = 0;
   eleIDMVAHZZ = 0;
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
   eleGSFChi2 = 0;
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
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleIDbit = 0;
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
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
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
   muFiredL1Trgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetMass = 0;
   jetPx = 0;
   jetPy = 0;
   jetPz = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetCSV2BJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetpfCombinedMVAV2BJetTags = 0;
   jetPFLooseId = 0;
   jetID = 0;
   jetPUID = 0;
   jetJECUnc = 0;
   jetJERSmearing = 0;
   jetJERSmearingUp = 0;
   jetJERSmearingDown = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetNNP = 0;
   jetMUF = 0;
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
   AK8JetNNP = 0;
   AK8JetMUF = 0;
   AK8Jetnconstituents = 0;
   AK8JetPFLooseId = 0;
   AK8JetPFTightLepVetoId = 0;
   AK8JetSoftDropMass = 0;
   AK8JetSoftDropMassCorr = 0;
   AK8JetPrunedMass = 0;
   AK8JetPrunedMassCorr = 0;
   AK8JetpfBoostedDSVBTag = 0;
   AK8JetDSVnewV4 = 0;
   AK8JetCSV = 0;
   AK8JetJECUnc = 0;
   AK8JetL2L3corr = 0;
   AK8puppiPt = 0;
   AK8puppiMass = 0;
   AK8puppiEta = 0;
   AK8puppiPhi = 0;
   AK8puppiTau1 = 0;
   AK8puppiTau2 = 0;
   AK8puppiTau3 = 0;
   AK8puppiSDL2L3corr = 0;
   AK8puppiSDMass = 0;
   AK8puppiSDMassL2L3Corr = 0;
   nAK8SDSJ = 0;
   AK8SDSJPt = 0;
   AK8SDSJEta = 0;
   AK8SDSJPhi = 0;
   AK8SDSJMass = 0;
   AK8SDSJE = 0;
   AK8SDSJCharge = 0;
   AK8SDSJFlavour = 0;
   AK8SDSJCSV = 0;
   nAK8puppiSDSJ = 0;
   AK8puppiSDSJPt = 0;
   AK8puppiSDSJEta = 0;
   AK8puppiSDSJPhi = 0;
   AK8puppiSDSJMass = 0;
   AK8puppiSDSJE = 0;
   AK8puppiSDSJCharge = 0;
   AK8puppiSDSJFlavour = 0;
   AK8puppiSDSJCSV = 0;
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
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
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
   fChain->SetBranchAddress("phoPrescale", &phoPrescale, &b_phoPrescale);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
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
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phomaxXtalenergyFull5x5", &phomaxXtalenergyFull5x5, &b_phomaxXtalenergyFull5x5);
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
   fChain->SetBranchAddress("phoCITKChIso", &phoCITKChIso, &b_phoCITKChIso);
   fChain->SetBranchAddress("phoCITKPhoIso", &phoCITKPhoIso, &b_phoCITKPhoIso);
   fChain->SetBranchAddress("phoCITKNeuIso", &phoCITKNeuIso, &b_phoCITKNeuIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoxtalBits", &phoxtalBits, &b_phoxtalBits);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("npfPho", &npfPho, &b_npfPho);
   fChain->SetBranchAddress("pfphoEt", &pfphoEt, &b_pfphoEt);
   fChain->SetBranchAddress("pfphoEta", &pfphoEta, &b_pfphoEta);
   fChain->SetBranchAddress("pfphoPhi", &pfphoPhi, &b_pfphoPhi);
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
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
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
   fChain->SetBranchAddress("eleIDMVA", &eleIDMVA, &b_eleIDMVA);
   fChain->SetBranchAddress("eleIDMVAHZZ", &eleIDMVAHZZ, &b_eleIDMVAHZZ);
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
   fChain->SetBranchAddress("eleGSFChi2", &eleGSFChi2, &b_eleGSFChi2);
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
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
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
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
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
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetPx", &jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", &jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", &jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags, &b_jetpfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetJERSmearing", &jetJERSmearing, &b_jetJERSmearing);
   fChain->SetBranchAddress("jetJERSmearingUp", &jetJERSmearingUp, &b_jetJERSmearingUp);
   fChain->SetBranchAddress("jetJERSmearingDown", &jetJERSmearingDown, &b_jetJERSmearingDown);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetNNP", &jetNNP, &b_jetNNP);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
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
   fChain->SetBranchAddress("AK8JetNNP", &AK8JetNNP, &b_AK8JetNNP);
   fChain->SetBranchAddress("AK8JetMUF", &AK8JetMUF, &b_AK8JetMUF);
   fChain->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
   fChain->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
   fChain->SetBranchAddress("AK8JetPFTightLepVetoId", &AK8JetPFTightLepVetoId, &b_AK8JetPFTightLepVetoId);
   fChain->SetBranchAddress("AK8JetSoftDropMass", &AK8JetSoftDropMass, &b_AK8JetSoftDropMass);
   fChain->SetBranchAddress("AK8JetSoftDropMassCorr", &AK8JetSoftDropMassCorr, &b_AK8JetSoftDropMassCorr);
   fChain->SetBranchAddress("AK8JetPrunedMass", &AK8JetPrunedMass, &b_AK8JetPrunedMass);
   fChain->SetBranchAddress("AK8JetPrunedMassCorr", &AK8JetPrunedMassCorr, &b_AK8JetPrunedMassCorr);
   fChain->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
   fChain->SetBranchAddress("AK8JetDSVnewV4", &AK8JetDSVnewV4, &b_AK8JetDSVnewV4);
   fChain->SetBranchAddress("AK8JetCSV", &AK8JetCSV, &b_AK8JetCSV);
   fChain->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
   fChain->SetBranchAddress("AK8JetL2L3corr", &AK8JetL2L3corr, &b_AK8JetL2L3corr);
   fChain->SetBranchAddress("AK8puppiPt", &AK8puppiPt, &b_AK8puppiPt);
   fChain->SetBranchAddress("AK8puppiMass", &AK8puppiMass, &b_AK8puppiMass);
   fChain->SetBranchAddress("AK8puppiEta", &AK8puppiEta, &b_AK8puppiEta);
   fChain->SetBranchAddress("AK8puppiPhi", &AK8puppiPhi, &b_AK8puppiPhi);
   fChain->SetBranchAddress("AK8puppiTau1", &AK8puppiTau1, &b_AK8puppiTau1);
   fChain->SetBranchAddress("AK8puppiTau2", &AK8puppiTau2, &b_AK8puppiTau2);
   fChain->SetBranchAddress("AK8puppiTau3", &AK8puppiTau3, &b_AK8puppiTau3);
   fChain->SetBranchAddress("AK8puppiSDL2L3corr", &AK8puppiSDL2L3corr, &b_AK8puppiSDL2L3corr);
   fChain->SetBranchAddress("AK8puppiSDMass", &AK8puppiSDMass, &b_AK8puppiSDMass);
   fChain->SetBranchAddress("AK8puppiSDMassL2L3Corr", &AK8puppiSDMassL2L3Corr, &b_AK8puppiSDMassL2L3Corr);
   fChain->SetBranchAddress("nAK8SDSJ", &nAK8SDSJ, &b_nAK8SDSJ);
   fChain->SetBranchAddress("AK8SDSJPt", &AK8SDSJPt, &b_AK8SDSJPt);
   fChain->SetBranchAddress("AK8SDSJEta", &AK8SDSJEta, &b_AK8SDSJEta);
   fChain->SetBranchAddress("AK8SDSJPhi", &AK8SDSJPhi, &b_AK8SDSJPhi);
   fChain->SetBranchAddress("AK8SDSJMass", &AK8SDSJMass, &b_AK8SDSJMass);
   fChain->SetBranchAddress("AK8SDSJE", &AK8SDSJE, &b_AK8SDSJE);
   fChain->SetBranchAddress("AK8SDSJCharge", &AK8SDSJCharge, &b_AK8SDSJCharge);
   fChain->SetBranchAddress("AK8SDSJFlavour", &AK8SDSJFlavour, &b_AK8SDSJFlavour);
   fChain->SetBranchAddress("AK8SDSJCSV", &AK8SDSJCSV, &b_AK8SDSJCSV);
   fChain->SetBranchAddress("nAK8puppiSDSJ", &nAK8puppiSDSJ, &b_nAK8puppiSDSJ);
   fChain->SetBranchAddress("AK8puppiSDSJPt", &AK8puppiSDSJPt, &b_AK8puppiSDSJPt);
   fChain->SetBranchAddress("AK8puppiSDSJEta", &AK8puppiSDSJEta, &b_AK8puppiSDSJEta);
   fChain->SetBranchAddress("AK8puppiSDSJPhi", &AK8puppiSDSJPhi, &b_AK8puppiSDSJPhi);
   fChain->SetBranchAddress("AK8puppiSDSJMass", &AK8puppiSDSJMass, &b_AK8puppiSDSJMass);
   fChain->SetBranchAddress("AK8puppiSDSJE", &AK8puppiSDSJE, &b_AK8puppiSDSJE);
   fChain->SetBranchAddress("AK8puppiSDSJCharge", &AK8puppiSDSJCharge, &b_AK8puppiSDSJCharge);
   fChain->SetBranchAddress("AK8puppiSDSJFlavour", &AK8puppiSDSJFlavour, &b_AK8puppiSDSJFlavour);
   fChain->SetBranchAddress("AK8puppiSDSJCSV", &AK8puppiSDSJCSV, &b_AK8puppiSDSJCSV);
   Notify();
}

Bool_t PostAnalyzer_Data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PostAnalyzer_Data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PostAnalyzer_Data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t PostAnalyzer_Data::PassHLT(vector<ULong64_t> trigBits){

  int count = 0;
  bool trigFired = false;

  for(Int_t i = 0; i < trigBits.size(); i++){
    if(HLTPho>>(trigBits[i]) & 1){            //HLTPho is a 64 bit integer having only those bits as one whose corresponding triggers have been fired
      count++;                                //by this particular event. see the correspondence between triggers and bits in ggAnalysis/ggNtuplizer/
    }                                         //plugins/ggNtuplizer_globalEvent.cc. trigBit in input represent the bit for trigger we want to check.
  }                                           //So if a trigger has been fired by evt, then its bit in HLTPho is 1 and right shift operator (>>) will
  if(count > 0) trigFired = true;                 //shift this bit to the first place and its bitwise AND (&) with 1 will be 1, otherwise not.
  return trigFired;
}

Bool_t PostAnalyzer_Data::PassHLT_Prescale(ULong64_t trigBit){ //HLTPhoIsPrescaled is 1 for trig bit whose prescale value > 1 i.e trig is prescaled
                                                               //So PassHLT_Prescale() will retrun true for the trigger with prescale > 1.
  bool trigPre = false;                                        //As for HLT_Photon165_HE10, prescale = 1, so this fun will return 0 for that.

  if((HLTPhoIsPrescaled >> trigBit) & 1){
    trigPre = true;
  }                          
                                         
  return trigPre;
}

/*
Bool_t PostAnalyzer_Data::GoodPrimaryVtx(Int_t &GoodVertex){

  Bool_t passVtx = false;
  GoodVertex = 0;

  for(Int_t i=0; i < nVtx; ++i){
    if( (fabs((*vtz)[i])) <= Cut_Vtx_z &&
        (*vndof)[i] >= Cut_Vtx_ndof    &&
        !((*isFake)[i])                &&
        (fabs((*vrho)[i])) <= Cut_Vtx_rho )
      GoodVertex++;
  }

  if(GoodVertex > 0) passVtx = true;

  return passVtx;

}
*/

//Cut Based Ph ID for 13 TeV 2016data(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2)
//Latest ID for full 2016 Rereco data (Spring16 selection)
Bool_t PostAnalyzer_Data::CutBasedPhotonID(Int_t ipho, TString phoWP){

  Bool_t PhID = false;

  if(phoWP == "loose"){ //loose 
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.0597)                   &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.01031)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 1.295)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 10.910 + 0.0148*(*phoEt)[ipho] + 0.000017*(*phoEt)[ipho]*(*phoEt)[ipho])                                      &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 3.630+0.0047*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0481)                     &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.03013)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 1.011)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 5.931 + 0.0163*(*phoEt)[ipho] + 0.000014*(*phoEt)[ipho]*(*phoEt)[ipho])                                         &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 6.641 + 0.0034*(*phoEt)[ipho]);

    }
  }
  if(phoWP == "medium"){ //medium
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.0396)                   &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.01022)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.441)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 2.725 + 0.0148*(*phoEt)[ipho] + 0.000017*(*phoEt)[ipho]*(*phoEt)[ipho])                                       &&
	((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.571 + 0.0047*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0219)                     &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.03001)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.442)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 1.715 + 0.0163*(*phoEt)[ipho] + 0.000014*(*phoEt)[ipho]*(*phoEt)[ipho])                                         &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 3.863 + 0.0034*(*phoEt)[ipho]);

    }
  }        
  if(phoWP == "tight"){ //tight
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.0269)                   &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.00994)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.202)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 0.264 + 0.0148*(*phoEt)[ipho] + 0.000017*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.362 + 0.0047*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0213)                     &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.03000)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.034)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 0.586 + 0.0163*(*phoEt)[ipho] + 0.000014*(*phoEt)[ipho]*(*phoEt)[ipho])                                           &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.617 + 0.0034*(*phoEt)[ipho]);

    }
  }
  return PhID;
}

//(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Selection implementation details for SPRING16)
Double_t PostAnalyzer_Data::EAChargedHadrons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.0360;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.0377;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0306;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0283;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.0254;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.0217;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.0167;

  return EffArea;

}

Double_t PostAnalyzer_Data::EANeutralHadrons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.0597;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.0807;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0629;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0197;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.0184;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.0284;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.0591;

  return EffArea;

}

Double_t PostAnalyzer_Data::EAPhotons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.1210;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1107;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0699;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.1056;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1457;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.1719;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.1998;

  return EffArea;

}

Int_t PostAnalyzer_Data::FirstGoodPhoton(TString phoWP){

  Int_t pc = -1;
  Bool_t ID = false;
  for(Int_t i = 0; i < nPho; i++){
    ID = CutBasedPhotonID(i, phoWP);
    if(ID){
      pc = i;
      break;
    }
  }
  return pc;
}


vector<Int_t> PostAnalyzer_Data::GoodPhotons(TString phoWP){

  vector<Int_t> goodphs;
  goodphs.clear();

  for(Int_t i = 0; i < nPho; i++){
    if(CutBasedPhotonID(i, phoWP) && ResSpikes(i) && (*phoEt)[i] > 30.0){
      goodphs.push_back(i);
    }
  }
  return goodphs;
}

Bool_t PostAnalyzer_Data::ResSpikes(Int_t i){
  Bool_t spikes = false;
  if( fabs((*phoSeedTime)[i]) < 3.0    &&  //time of arrival of ith photon at seed crystal
      (*phoSigmaIEtaIEtaFull5x5)[i] > 0.001   &&
      (*phoSigmaIPhiIPhiFull5x5)[i] > 0.001   &&
      //fabs(GetLICTD(i)) < 5.0               &&   //LICTD is the largest time difference between the seed crystal and the any other crystal 
      (*phoR9Full5x5)[i] < 1.0){
    spikes = true;
  }
  return spikes;
}

//Recommended JetID for 13 TeV 2016 data(https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016)
Bool_t PostAnalyzer_Data::JetId(Int_t iJet, TString jetWP){

  Bool_t JetID = false;

  if(fabs((*jetEta)[iJet]) <= 2.7){
    if(jetWP == "loose"){

      JetID = ((*jetNHF)[iJet] < 0.99 && (*jetNEF)[iJet] < 0.99 && (*jetNConstituents)[iJet] > 1 ) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0 && (*jetCEF)[iJet] < 0.99) || fabs((*jetEta)[iJet]) > 2.4);
	 
    }
    if(jetWP == "tight"){

      JetID = ((*jetNHF)[iJet] < 0.90 && (*jetNEF)[iJet] < 0.90 && (*jetNConstituents)[iJet] > 1 ) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0 && (*jetCEF)[iJet] < 0.99) || fabs((*jetEta)[iJet]) > 2.4);

    }  
    if(jetWP == "tightLepVeto"){

      JetID = ((*jetNHF)[iJet] < 0.90 && (*jetNEF)[iJet] < 0.90 && (*jetNConstituents)[iJet] > 1 && (*jetMUF)[iJet] < 0.8) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0 && (*jetCEF)[iJet] < 0.90) || fabs((*jetEta)[iJet]) > 2.4);

    }
  }

  if(fabs((*jetEta)[iJet]) > 2.7 && fabs((*jetEta)[iJet]) <= 3.0){
    if(jetWP == "loose"){
      JetID = (*jetNEF)[iJet] > 0.01 && (*jetNHF)[iJet] < 0.98 && (*jetNNP)[iJet] > 2;
    }
    if(jetWP == "tight"){
      JetID = (*jetNEF)[iJet] > 0.01 && (*jetNHF)[iJet] < 0.98 && (*jetNNP)[iJet] > 2;
    }
  }
  
  if(fabs((*jetEta)[iJet]) > 3.0){
    if(jetWP == "loose"){
      JetID = (*jetNEF)[iJet] < 0.90 && (*jetNNP)[iJet] > 10;
    }
    if(jetWP == "tight"){
      JetID = (*jetNEF)[iJet] < 0.90 && (*jetNNP)[iJet]> 10;
    }
  }

  return JetID;
} 


Int_t PostAnalyzer_Data::FirstGoodJet(TString jetWP){

  Int_t jc = -1;
  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    ID = JetId(i, jetWP);
    if(ID && (*jetPt)[i] > 30.0){
      double minDR = 99.0;
      for(Int_t ph = 0; ph < nPho; ph++){
        if(CutBasedPhotonID(ph, "loose") && ResSpikes(ph) && (*phoEt)[ph] > 30.0){
          double temp_dR = GetdR((*phoSCEta)[ph], (*jetEta)[i], (*phoSCPhi)[ph], (*jetPhi)[i]);
          if(temp_dR < minDR) minDR = temp_dR;
        }
      }
      if(minDR > 0.5 && minDR < 99.0){
	jc = i;
	break;
      }
    }
  }
  return jc;
}

vector<Int_t> PostAnalyzer_Data::GoodJets(TString jetWP){

  vector<Int_t> goodjets;
  goodjets.clear();

  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    ID = JetId(i, jetWP);
    if(ID && (*jetPt)[i] > 30.0){
      double minDR = 99.0;
      for(Int_t ph = 0; ph < nPho; ph++){
        if(CutBasedPhotonID(ph, "loose") && ResSpikes(ph) && (*phoEt)[ph] > 30.0){
          double temp_dR = GetdR((*phoSCEta)[ph], (*jetEta)[i], (*phoSCPhi)[ph], (*jetPhi)[i]);
          if(temp_dR < minDR) minDR = temp_dR;
        }
      }
      if(minDR > 0.5 && minDR < 99.0) goodjets.push_back(i);
    }
  }
  return goodjets;
}


//80XReReco recommendations (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco)
Bool_t PostAnalyzer_Data::CSVv2bTag(Int_t ijet, TString WP){

  Bool_t passTag = false;

  if(WP == "L"){ //loose
    if((*jetCSV2BJetTags)[ijet] > 0.5426) passTag = true;
  }
  if(WP == "M"){ //medium
    if((*jetCSV2BJetTags)[ijet] > 0.8484) passTag = true;
  }
  if(WP == "T"){ //tight
    if((*jetCSV2BJetTags)[ijet] > 0.9535) passTag = true;
  }

  return passTag;
}


Double_t PostAnalyzer_Data::GetdEta(Double_t eta1, Double_t eta2){

  Double_t dEta = fabs(eta1 - eta2);
  return dEta;
}

Double_t PostAnalyzer_Data::GetdPhi(Double_t phi1, Double_t phi2){

  Double_t dphi = (phi1 - phi2);
  Double_t twoPi = 2.0*(TMath::Pi());

  if(dphi < 0) dphi = - dphi;
  if(dphi >= (twoPi - dphi)) dphi = twoPi - dphi;

  return dphi;
}

Double_t PostAnalyzer_Data::GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

  Double_t dEta = GetdEta(eta1, eta2);
  Double_t dPhi = GetdPhi(phi1, phi2);

  Double_t dR = 0.0;
  dR = sqrt(dEta*dEta + dPhi*dPhi);

  return dR;
}

//----------------------                                                          
// Compute cosThetaStar                                                                                                                            
//---------------------
Double_t PostAnalyzer_Data::GetCosThetaStar(Double_t eta1, Double_t eta2){ 
  Double_t theta = tanh( GetdEta(eta1,eta2)/2.0 );                                                                                                  
  return theta;                                                                                                                                     
}      


Double_t PostAnalyzer_Data::GetInvtMass(Int_t pho, Int_t jet){

  Double_t mass = 0.0;

  TLorentzVector Pho;
  Pho.SetPtEtaPhiE((*phoEt)[pho], (*phoSCEta)[pho], (*phoSCPhi)[pho], (*phoE)[pho]);

  TLorentzVector Jet;
  Jet.SetPtEtaPhiE((*jetPt)[jet], (*jetEta)[jet], (*jetPhi)[jet], (*jetEn)[jet] );

  mass = (Pho+Jet).M();

  return mass;
}

void PostAnalyzer_Data::BookHistograms(){
  file->cd();

  char name[100];
  /*
  const Int_t nMassBins_qstar = 119;
  const Double_t MassBin_qstar[nMassBins_qstar+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};
  */
  const Int_t nMassBins = 119;
  const Double_t MassBin[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

  std::string cut0[7] = {"noCut", "noMassCut", "MassCut", "1BTag_noMassCut", "1BTag_MassCut", "0BTag_noMassCut", "0BTag_MassCut"};
  for(Int_t hist0 = 0; hist0 < 7; ++hist0){

    sprintf(name, "h_PhotonPt_%s",cut0[hist0].c_str());
    h_PhotonPt[hist0] = new TH1F(name,"Pt distribution of photons",150,0.0,6000.0);
    h_PhotonPt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_PhotonPt[hist0]->GetYaxis()->CenterTitle();
    h_PhotonPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_PhotonPt[hist0]->GetXaxis()->CenterTitle();
    h_PhotonPt[hist0]->Sumw2();
  
    sprintf(name, "h_PhotonCalibPt_%s",cut0[hist0].c_str());
    h_PhotonCalibPt[hist0] = new TH1F(name,"Pt distribution of Calibrated photons",150,0.0,6000.0);
    h_PhotonCalibPt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_PhotonCalibPt[hist0]->GetYaxis()->CenterTitle();
    h_PhotonCalibPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_PhotonCalibPt[hist0]->GetXaxis()->CenterTitle();
    h_PhotonCalibPt[hist0]->Sumw2();

    sprintf(name, "h_PhotonEta_%s",cut0[hist0].c_str());
    h_PhotonEta[hist0] = new TH1F(name,"Eta distribution of photons",100,-2.5,2.5);
    h_PhotonEta[hist0]->GetYaxis()->SetTitle("Events");  h_PhotonEta[hist0]->GetYaxis()->CenterTitle();
    h_PhotonEta[hist0]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_PhotonEta[hist0]->GetXaxis()->CenterTitle();
    h_PhotonEta[hist0]->Sumw2();

    sprintf(name, "h_PhotonPhi_%s",cut0[hist0].c_str());
    h_PhotonPhi[hist0] = new TH1F(name,"Phi distribution of photons",100,-4.0,4.0);
    h_PhotonPhi[hist0]->GetYaxis()->SetTitle("Events");  h_PhotonPhi[hist0]->GetYaxis()->CenterTitle();
    h_PhotonPhi[hist0]->GetXaxis()->SetTitle("#phi^{#gamma}");   h_PhotonPhi[hist0]->GetXaxis()->CenterTitle();
    h_PhotonPhi[hist0]->Sumw2();

    sprintf(name, "h_Photon_SigmaIEtaIEta_%s",cut0[hist0].c_str());
    h_Photon_SigmaIEtaIEta[hist0] = new TH1F(name,"Photon SigmaIetaIeta Distribution",100,0.0,0.05);
    h_Photon_SigmaIEtaIEta[hist0]->GetYaxis()->SetTitle("Events");  h_Photon_SigmaIEtaIEta[hist0]->GetYaxis()->CenterTitle();
    h_Photon_SigmaIEtaIEta[hist0]->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");  h_Photon_SigmaIEtaIEta[hist0]->GetXaxis()->CenterTitle();
    h_Photon_SigmaIEtaIEta[hist0]->Sumw2();

    sprintf(name, "h_Photon_R9_%s",cut0[hist0].c_str());
    h_Photon_R9[hist0] = new TH1F(name,"Photon R9 Distribution",100,0.0,10.0);
    h_Photon_R9[hist0]->GetYaxis()->SetTitle("Events");       h_Photon_R9[hist0]->GetYaxis()->CenterTitle();
    h_Photon_R9[hist0]->GetXaxis()->SetTitle("Photon_r9");    h_Photon_R9[hist0]->GetXaxis()->CenterTitle();
    h_Photon_R9[hist0]->Sumw2();

    sprintf(name, "h_Photon_HoverE_%s",cut0[hist0].c_str());
    h_Photon_HoverE[hist0] = new TH1F(name,"Photon HoverE Distribution",50,0.0,0.1);
    h_Photon_HoverE[hist0]->GetYaxis()->SetTitle("Events");   h_Photon_HoverE[hist0]->GetYaxis()->CenterTitle();
    h_Photon_HoverE[hist0]->GetXaxis()->SetTitle("H/E");      h_Photon_HoverE[hist0]->GetXaxis()->CenterTitle();
    h_Photon_HoverE[hist0]->Sumw2();

    sprintf(name, "h_Photon_EleVeto_%s",cut0[hist0].c_str());
    h_Photon_EleVeto[hist0] = new TH1F(name,"Photon ElectronVeto",3,0,3);
    h_Photon_EleVeto[hist0]->GetYaxis()->SetTitle("Events");     h_Photon_EleVeto[hist0]->GetYaxis()->CenterTitle();
    h_Photon_EleVeto[hist0]->GetXaxis()->SetTitle("EleVeto");    h_Photon_EleVeto[hist0]->GetXaxis()->CenterTitle();
    h_Photon_EleVeto[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFChIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFChIso[hist0] = new TH1F(name,"Rho Corrected PF Charged Isolation",625,0,25);
    h_Photon_CorrPFChIso[hist0]->GetYaxis()->SetTitle("Events");        h_Photon_CorrPFChIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFChIso[hist0]->GetXaxis()->SetTitle("CorrPFChIso");   h_Photon_CorrPFChIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFChIso[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFNeuIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFNeuIso[hist0] = new TH1F(name,"Rho Corrected PF Neutral Isolation",100,0,300);
    h_Photon_CorrPFNeuIso[hist0]->GetYaxis()->SetTitle("Events");    h_Photon_CorrPFNeuIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFNeuIso[hist0]->GetXaxis()->SetTitle("CorrPFNeuIso");    h_Photon_CorrPFNeuIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFNeuIso[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFPhoIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFPhoIso[hist0] = new TH1F(name,"Rho Corrected PF Photon Isolation",125,0,25);
    h_Photon_CorrPFPhoIso[hist0]->GetYaxis()->SetTitle("Events");    h_Photon_CorrPFPhoIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFPhoIso[hist0]->GetXaxis()->SetTitle("CorrPFPhoIso"); h_Photon_CorrPFPhoIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFPhoIso[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetPt_%s",cut0[hist0].c_str());
    h_bJetPt[hist0] = new TH1F(name,"Pt distribution of jets",150,0.0,6000.0);
    h_bJetPt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{Jet} (GeV)");  h_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetEta_%s",cut0[hist0].c_str());
    h_bJetEta[hist0] = new TH1F(name,"Eta distribution of jets",120,-3.0,3.0);
    h_bJetEta[hist0]->GetYaxis()->SetTitle("Events");      h_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{Jet}");  h_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_bJetEta[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetPhi_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetPhi_%s",cut0[hist0].c_str());
    h_bJetPhi[hist0] = new TH1F(name,"Phi distribution of jets",100,-4.0,4.0);
    h_bJetPhi[hist0]->GetYaxis()->SetTitle("Events");      h_bJetPhi[hist0]->GetYaxis()->CenterTitle();
    h_bJetPhi[hist0]->GetXaxis()->SetTitle("#phi^{Jet}");  h_bJetPhi[hist0]->GetXaxis()->CenterTitle();
    h_bJetPhi[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_Mt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_Mt_%s",cut0[hist0].c_str());
    h_bJet_Mt[hist0] = new TH1F(name,"Transverse mass distribution of jets",150,0.0,6000.0);
    h_bJet_Mt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_Mt[hist0]->GetYaxis()->CenterTitle();
    h_bJet_Mt[hist0]->GetXaxis()->SetTitle("M_{T}^{Jet}");  h_bJet_Mt[hist0]->GetXaxis()->CenterTitle();
    h_bJet_Mt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_area_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_area_%s",cut0[hist0].c_str());
    h_bJet_area[hist0] = new TH1F(name,"Jets area",25,0.0,1000.0);
    h_bJet_area[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_area[hist0]->GetYaxis()->CenterTitle();
    h_bJet_area[hist0]->GetXaxis()->SetTitle("JetArea");  h_bJet_area[hist0]->GetXaxis()->CenterTitle();
    h_bJet_area[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_Mass_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_Mass_%s",cut0[hist0].c_str());
    h_bJet_Mass[hist0] = new TH1F(name,"Mass distribution of jets",150,0.0,6000.0);
    h_bJet_Mass[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_Mass[hist0]->GetYaxis()->CenterTitle();
    h_bJet_Mass[hist0]->GetXaxis()->SetTitle("M^{Jet}");  h_bJet_Mass[hist0]->GetXaxis()->CenterTitle();
    h_bJet_Mass[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NHEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NHEF_%s",cut0[hist0].c_str());
    h_bJet_NHEF[hist0] = new TH1F(name,"Neutral Hadron Energy Fraction",25,0,1);
    h_bJet_NHEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NHEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NHEF[hist0]->GetXaxis()->SetTitle("Jet_NHEF");    h_bJet_NHEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NHEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NEEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NEEF_%s",cut0[hist0].c_str());
    h_bJet_NEEF[hist0] = new TH1F(name,"Neutral Em Energy Fraction",25,0,1);
    h_bJet_NEEF[hist0]->GetYaxis()->SetTitle("Events");       h_bJet_NEEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NEEF[hist0]->GetXaxis()->SetTitle("Jet_NEEF");     h_bJet_NEEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NEEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NConst_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NConst_%s",cut0[hist0].c_str());
    h_bJet_NConst[hist0] = new TH1F(name,"No. of Constituents",100,0,100);
    h_bJet_NConst[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NConst[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NConst[hist0]->GetXaxis()->SetTitle("Jet_NConst");  h_bJet_NConst[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NConst[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_CHEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_CHEF_%s",cut0[hist0].c_str());
    h_bJet_CHEF[hist0] = new TH1F(name,"Charged Hadron Energy Fraction",25,0,1);
    h_bJet_CHEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_CHEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_CHEF[hist0]->GetXaxis()->SetTitle("Jet_CHEF");    h_bJet_CHEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_CHEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_ChMult_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_ChMult_%s",cut0[hist0].c_str());
    h_bJet_ChMult[hist0] = new TH1F(name,"Charged Multiplicity",100,0,100);
    h_bJet_ChMult[hist0]->GetYaxis()->SetTitle("Events");       h_bJet_ChMult[hist0]->GetYaxis()->CenterTitle();
    h_bJet_ChMult[hist0]->GetXaxis()->SetTitle("Jet_ChMult");   h_bJet_ChMult[hist0]->GetXaxis()->CenterTitle();
    h_bJet_ChMult[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_CEEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_CEEF_%s",cut0[hist0].c_str());
    h_bJet_CEEF[hist0] = new TH1F(name,"Charged Em Energy Fraction",25,0,1);
    h_bJet_CEEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_CEEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_CEEF[hist0]->GetXaxis()->SetTitle("Jet_CEEF");    h_bJet_CEEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_CEEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_MUF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_MUF_%s",cut0[hist0].c_str());
    h_bJet_MUF[hist0] = new TH1F(name,"Muon Fraction",25,0,1);
    h_bJet_MUF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_MUF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_MUF[hist0]->GetXaxis()->SetTitle("Jet_MUF");    h_bJet_MUF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_MUF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NNP_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NNP_%s",cut0[hist0].c_str());
    h_bJet_NNP[hist0] = new TH1F(name,"Number of neutral particles",100,0,100);
    h_bJet_NNP[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NNP[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NNP[hist0]->GetXaxis()->SetTitle("Jet_NNP");    h_bJet_NNP[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NNP[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJetInvtMass_VarBin_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJetInvtMass_VarBin_%s",cut0[hist0].c_str());
    h_GbJetInvtMass_VarBin[hist0] = new TH1F(name, "Invt mass distribution", nMassBins, MassBin);
    h_GbJetInvtMass_VarBin[hist0]->GetYaxis()->SetTitle("Events/VarBin"); h_GbJetInvtMass_VarBin[hist0]->GetYaxis()->CenterTitle();
    h_GbJetInvtMass_VarBin[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");    h_GbJetInvtMass_VarBin[hist0]->GetXaxis()->CenterTitle();
    h_GbJetInvtMass_VarBin[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJetInvtMass_UnitBin_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJetInvtMass_UnitBin_%s",cut0[hist0].c_str());
    h_GbJetInvtMass_UnitBin[hist0] = new TH1F(name,"Invt mass distribution", 14000, 0.0, 14000.0);
    h_GbJetInvtMass_UnitBin[hist0]->GetYaxis()->SetTitle("Events/UnitBin");   h_GbJetInvtMass_UnitBin[hist0]->GetYaxis()->CenterTitle();
    h_GbJetInvtMass_UnitBin[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)"); h_GbJetInvtMass_UnitBin[hist0]->GetXaxis()->CenterTitle();
    h_GbJetInvtMass_UnitBin[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dEta_%s",cut0[hist0].c_str());
    h_GbJet_dEta[hist0] = new TH1F(name, "delta Eta dist", 120, 0, 6);
    h_GbJet_dEta[hist0]->GetYaxis()->SetTitle("Events");         h_GbJet_dEta[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dEta[hist0]->GetXaxis()->SetTitle("#Delta #eta");    h_GbJet_dEta[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dEta[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dPhi_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dPhi_%s",cut0[hist0].c_str());
    h_GbJet_dPhi[hist0] = new TH1F(name, "delta Phi dist", 64, 0, 3.2);
    h_GbJet_dPhi[hist0]->GetYaxis()->SetTitle("Events");       h_GbJet_dPhi[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dPhi[hist0]->GetXaxis()->SetTitle("#Delta #phi");  h_GbJet_dPhi[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dPhi[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dR_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dR_%s",cut0[hist0].c_str());
    h_GbJet_dR[hist0] = new TH1F(name, "delta R dist", 100, 0.0, 10.0);
    h_GbJet_dR[hist0]->GetYaxis()->SetTitle("Events");          h_GbJet_dR[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dR[hist0]->GetXaxis()->SetTitle("#Delta R");        h_GbJet_dR[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dR[hist0]->Sumw2();

    sprintf(name, "h_cosThetaStar_%s",cut0[hist0].c_str());
    h_cosThetaStar[hist0] = new TH1F(name,"Cos theta star of photon+jet",50,0,1);
    h_cosThetaStar[hist0]->GetYaxis()->SetTitle("Events");        h_cosThetaStar[hist0]->GetYaxis()->CenterTitle();
    h_cosThetaStar[hist0]->GetXaxis()->SetTitle("cos#theta*");    h_cosThetaStar[hist0]->GetXaxis()->CenterTitle();
    h_cosThetaStar[hist0]->Sumw2();

    sprintf(name, "h_PFMet_%s",cut0[hist0].c_str());
    h_PFMet[hist0] = new TH1F(name, "PFMet distribution", 200,0.0,1000.0);
    h_PFMet[hist0]->GetYaxis()->SetTitle("Events/5 GeV");           h_PFMet[hist0]->GetYaxis()->CenterTitle();
    h_PFMet[hist0]->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");    h_PFMet[hist0]->GetXaxis()->CenterTitle();
    h_PFMet[hist0]->Sumw2();

    sprintf(name, "h_SumPFMet_%s",cut0[hist0].c_str());
    h_SumPFMet[hist0]  = new TH1F(name,"SumET PF Met distribution",80,0.0,4000.0);
    h_SumPFMet[hist0]->GetYaxis()->SetTitle("Events/50 GeV");            h_SumPFMet[hist0]->GetYaxis()->CenterTitle();
    h_SumPFMet[hist0]->GetXaxis()->SetTitle("#sum#slash{E}_{T} (GeV)");  h_SumPFMet[hist0]->GetXaxis()->CenterTitle();
    h_SumPFMet[hist0]->Sumw2();

    sprintf(name, "h_MetBySumMET_%s",cut0[hist0].c_str());
    h_MetBySumMET[hist0]  = new TH1F(name,"MET / SumET PF Met",50,0.0,1.0);
    h_MetBySumMET[hist0]->GetYaxis()->SetTitle("Events");        h_MetBySumMET[hist0]->GetYaxis()->CenterTitle();
    h_MetBySumMET[hist0]->GetXaxis()->SetTitle("#slash{E}_{T}/#sum#slash{E}_{T}"); h_MetBySumMET[hist0]->GetXaxis()->CenterTitle();
    h_MetBySumMET[hist0]->Sumw2();

    sprintf(name, "h_PFMetVsGJmass_%s",cut0[hist0].c_str());
    h_PFMetVsGJmass[hist0] = new TH2F(name, "PFMet Vs InvtMass", 1400, 0.0, 14000.0, 150, 0.0, 6000.0);
    h_PFMetVsGJmass[hist0]->GetYaxis()->SetTitle("PFMET (GeV)");    h_PFMetVsGJmass[hist0]->GetYaxis()->CenterTitle();
    h_PFMetVsGJmass[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");  h_PFMetVsGJmass[hist0]->GetXaxis()->CenterTitle();
    h_PFMetVsGJmass[hist0]->Sumw2();

    sprintf(name, "h_PFMetOverSumEtVsGJmass_%s",cut0[hist0].c_str());
    h_PFMetOverSumEtVsGJmass[hist0] = new TH2F(name, "PFMet/SumEt Vs InvtMass", 1400, 0.0, 14000.0, 1000, 0.0, 10.0);
    h_PFMetOverSumEtVsGJmass[hist0]->GetYaxis()->SetTitle("PFMET/SumEt (GeV)");  h_PFMetOverSumEtVsGJmass[hist0]->GetYaxis()->CenterTitle();
    h_PFMetOverSumEtVsGJmass[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");  h_PFMetOverSumEtVsGJmass[hist0]->GetXaxis()->CenterTitle();
    h_PFMetOverSumEtVsGJmass[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_PhPt_vs_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_PhPt_vs_bJetPt_%s",cut0[hist0].c_str());
    h_PhPt_vs_bJetPt[hist0] = new TH2F(name,"Pt of Photon vs Jet ",120,0.0,4800.0,120,0.0,4800.0);
    h_PhPt_vs_bJetPt[hist0]->GetYaxis()->SetTitle("P_{T}^{Jet}");   h_PhPt_vs_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_PhPt_vs_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma}");  h_PhPt_vs_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_PhPt_vs_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_PhEta_vs_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_PhEta_vs_bJetEta_%s",cut0[hist0].c_str());
    h_PhEta_vs_bJetEta[hist0] = new TH2F(name,"eta of Photon vs Jet ",100,-2.5,2.5,100,-2.5,2.5);
    h_PhEta_vs_bJetEta[hist0]->GetYaxis()->SetTitle("#eta^{Jet}");     h_PhEta_vs_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_PhEta_vs_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_PhEta_vs_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_PhEta_vs_bJetEta[hist0]->Sumw2();

    sprintf(name, "h_CSVv2Dist_%s",cut0[hist0].c_str());
    h_CSVv2Dist[hist0] = new TH1F(name, "CSVv2 bTagger Distribution", 500, 0.0, 5.0);
    h_CSVv2Dist[hist0]->GetYaxis()->SetTitle("Events");      h_CSVv2Dist[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2Dist[hist0]->GetXaxis()->SetTitle("CSVv2");       h_CSVv2Dist[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2Dist[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_CSVv2_vs_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_CSVv2_vs_bJetPt_%s",cut0[hist0].c_str());
    h_CSVv2_vs_bJetPt[hist0] = new TH2F(name,"CSVv2 vs Jet Pt", 120, 0.0, 4800.0, 500, 0.0, 5.0);
    h_CSVv2_vs_bJetPt[hist0]->GetYaxis()->SetTitle("CSVv2");         h_CSVv2_vs_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2_vs_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{Jet}");   h_CSVv2_vs_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2_vs_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_CSVv2_vs_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_CSVv2_vs_bJetEta_%s",cut0[hist0].c_str());
    h_CSVv2_vs_bJetEta[hist0] = new TH2F(name,"CSVv2 vs Jet Eta", 100, -2.5, 2.5, 500, 0.0, 5.0);
    h_CSVv2_vs_bJetEta[hist0]->GetYaxis()->SetTitle("CSVv2");         h_CSVv2_vs_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2_vs_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{Jet}");    h_CSVv2_vs_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2_vs_bJetEta[hist0]->Sumw2();

    sprintf(name, "h_goodPV_%s",cut0[hist0].c_str());
    h_goodPV[hist0] = new TH1F(name,"No. of Good Primary Vertices", 50, 0, 50);
    h_goodPV[hist0]->GetYaxis()->SetTitle("Events");   h_goodPV[hist0]->GetYaxis()->CenterTitle();
    h_goodPV[hist0]->GetXaxis()->SetTitle("nPV");      h_goodPV[hist0]->GetXaxis()->CenterTitle();
    h_goodPV[hist0]->Sumw2();

    sprintf(name, "h_nIsoPhotons_%s",cut0[hist0].c_str());
    h_nIsoPhotons[hist0] = new TH1F(name,"No. of Isolated Photons", 400, 0, 200);
    h_nIsoPhotons[hist0]->GetYaxis()->SetTitle("Events");         h_nIsoPhotons[hist0]->GetYaxis()->CenterTitle();
    h_nIsoPhotons[hist0]->GetXaxis()->SetTitle("nIsoPhotons");    h_nIsoPhotons[hist0]->GetXaxis()->CenterTitle();
    h_nIsoPhotons[hist0]->Sumw2();

    sprintf(name, "h_nGoodPhotons_%s",cut0[hist0].c_str());
    h_nGoodPhotons[hist0] = new TH1F(name,"No. of Good Isolated Photons", 400, 0, 200);
    h_nGoodPhotons[hist0]->GetYaxis()->SetTitle("Events");          h_nGoodPhotons[hist0]->GetYaxis()->CenterTitle();
    h_nGoodPhotons[hist0]->GetXaxis()->SetTitle("nGoodPhotons");    h_nGoodPhotons[hist0]->GetXaxis()->CenterTitle();
    h_nGoodPhotons[hist0]->Sumw2();

    sprintf(name, "h_IsoPhotonIdxVsPt_%s",cut0[hist0].c_str());
    h_IsoPhotonIdxVsPt[hist0] = new TH2F(name,"Isolated Photons Index Vs Pt", 120, 0.0, 4800.0, 10, 0.0, 10.0);
    h_IsoPhotonIdxVsPt[hist0]->GetYaxis()->SetTitle("IsoPhotonIdx");      h_IsoPhotonIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_IsoPhotonIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{#gamma}");    h_IsoPhotonIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_IsoPhotonIdxVsPt[hist0]->Sumw2();

    sprintf(name, "h_GoodPhotonIdxVsPt_%s",cut0[hist0].c_str());
    h_GoodPhotonIdxVsPt[hist0] = new TH2F(name,"Good Isolated Photons Idx Vs Pt", 120, 0.0, 4800.0, 10, 0.0, 10.0);
    h_GoodPhotonIdxVsPt[hist0]->GetYaxis()->SetTitle("GoodPhotonIdx");         h_GoodPhotonIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_GoodPhotonIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{#gamma}");       h_GoodPhotonIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_GoodPhotonIdxVsPt[hist0]->Sumw2();

    sprintf(name, "h_nJets_%s",cut0[hist0].c_str());
    h_nJets[hist0] = new TH1F(name,"No. of Jets", 400, 0, 200);
    h_nJets[hist0]->GetYaxis()->SetTitle("Events");      h_nJets[hist0]->GetYaxis()->CenterTitle();
    h_nJets[hist0]->GetXaxis()->SetTitle("nJets");       h_nJets[hist0]->GetXaxis()->CenterTitle();
    h_nJets[hist0]->Sumw2();

    sprintf(name, "h_JetIdxVsPt_%s",cut0[hist0].c_str());
    h_JetIdxVsPt[hist0] = new TH2F(name,"Jet Idx vs Pt", 120, 0.0, 4800.0, 20, 0.0, 20.0);
    h_JetIdxVsPt[hist0]->GetYaxis()->SetTitle("Jet Idx");           h_JetIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_JetIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{Jet}");       h_JetIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_JetIdxVsPt[hist0]->Sumw2();

  }

  //Position of photon and jet
  h_PC = new TH1F("h_PC", "Photon Candidate", 10, 0, 10);
  h_PC->GetYaxis()->SetTitle("Events");                       h_PC->GetYaxis()->CenterTitle();
  h_PC->GetXaxis()->SetTitle("Position of Photon");           h_PC->GetXaxis()->CenterTitle();
  h_PC->Sumw2();

  h_JC = new TH1F("h_JC", "Jet Candidate", 20, 0, 20);
  h_JC->GetYaxis()->SetTitle("Events");                       h_JC->GetYaxis()->CenterTitle();
  h_JC->GetXaxis()->SetTitle("Position of Jet");              h_JC->GetXaxis()->CenterTitle();
  h_JC->Sumw2();
  
  //Trigger Turn-on
  const Int_t nTrigPtBins = 31; 
  const Double_t TrigPtBins[nTrigPtBins+1] = { 0, 50, 100, 120, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 250, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000};       

  std::string cut1[2] = {"deno", "num"};
  for(Int_t hist1 =0; hist1 < 2; ++hist1){
  
    sprintf(name, "h_TrigPhotonPt_%s",cut1[hist1].c_str());
    h_TrigPhotonPt[hist1] = new TH1F(name, "pt of trigger photon", nTrigPtBins, TrigPtBins);
    h_TrigPhotonPt[hist1]->GetYaxis()->SetTitle("");      h_TrigPhotonPt[hist1]->GetYaxis()->CenterTitle();                      
    h_TrigPhotonPt[hist1]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");   h_TrigPhotonPt[hist1]->GetXaxis()->CenterTitle();
    h_TrigPhotonPt[hist1]->Sumw2();
  }

  h_TrigGJmass_Jetpt60 = new TH1F("h_TrigGJmass_Jetpt60", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt60->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt60->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt60->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt60->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt60->Sumw2();

  h_TrigGJmass_Jetpt100 = new TH1F("h_TrigGJmass_Jetpt100", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt100->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt100->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt100->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt100->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt100->Sumw2();

  h_TrigGJmass_Jetpt120 = new TH1F("h_TrigGJmass_Jetpt120", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt120->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt120->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt120->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt120->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt120->Sumw2();

  h_TrigGJmass_Jetpt150 = new TH1F("h_TrigGJmass_Jetpt150", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt150->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt150->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt150->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt150->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt150->Sumw2();

  h_TrigGJmass_Jetpt190 = new TH1F("h_TrigGJmass_Jetpt190", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt190->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt190->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt190->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt190->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt190->Sumw2();

  //Defining the histogram from filling number of events after various cuts
  const int nbins_qstar = 11;
  const int nbins_bstar = 16;

  TString CutFlowLabels_qstar[nbins_qstar] = {"Total", "HLT", "PrimaryVtx_&_METFilters", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "DPhi", "DEta", "MassCut"};
  TString CutFlowLabels_bstar[nbins_bstar] = {"Total", "HLT", "PrimaryVtx_&_METFilters", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "1BTag", "1BTag_DPhi", "1BTag_DEta", "1BTag_MassCut", "0BTag", "0BTag_DPhi", "0BTag_DEta", "0BTag_MassCut"};

  h_CutFlow_qstar = new TH1F("h_CutFlow_qstar", "Events Passing Various Cuts for qstar", nbins_qstar, 0, nbins_qstar);
  h_CutFlow_qstar->GetYaxis()->SetTitle("Events");         h_CutFlow_qstar->GetYaxis()->CenterTitle();
  for(int i = 0; i < nbins_qstar; i++){
    h_CutFlow_qstar->GetXaxis()->SetBinLabel(i+1, CutFlowLabels_qstar[i]);
  }
  h_CutFlow_bstar = new TH1F("h_CutFlow_bstar", "Events Passing Various Cuts for bstar", nbins_bstar, 0, nbins_bstar);
  h_CutFlow_bstar->GetYaxis()->SetTitle("Events");         h_CutFlow_bstar->GetYaxis()->CenterTitle();
  for(int i = 0; i < nbins_bstar; i++){
    h_CutFlow_bstar->GetXaxis()->SetBinLabel(i+1, CutFlowLabels_bstar[i]);
  }


}
#endif // #ifdef PostAnalyzer_Data_cxx
 
EOF


cat > analysis_${filenameTag}.C <<EOF
#include "PostAnalyzer_Data.C"
#include "TROOT.h"
int main(){
    PostAnalyzer_Data a;
    a.Loop();
    return 0;
}

EOF


####Compilation

g++ -Wno-deprecated analysis_${filenameTag}.C -o ${filenameTag}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`

####Execution
#./${filenameTag}.exe

###Submit jobs

chmod 775 MakeCondorFiles_local.csh

./MakeCondorFiles_local.csh ${filenameTag}

##**************************************************
@ sf = ${sf} + ${r}
@ maxf = ${maxf} + ${r}

end ##end of while loop
##**************************************************
end ##end of for loop##


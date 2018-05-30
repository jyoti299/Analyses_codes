#define PostAnalyzer_Data_cxx
#include "PostAnalyzer_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//To be removed in script
//------------------------
int main(){
  PostAnalyzer_Data t;
  t.Loop();
  return 0;
}
//-----------------------

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

   Cut_Photon_pt = 200.0; // GeV
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 200.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.5;

   Cut_GJInvtMass = 695.0;

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

   //Comment it in script
   //---------------------------
   //Define Output file here
   file = new TFile("PostAnalyzer_Data.root", "RECREATE");
   //---------------------------

   //Uncomment this in script
   /*
   //Define Output file here
   //TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   //file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   file = new TFile(OutputFile+".root", "RECREATE");
   */

   //Define Histograms here
   BookHistograms();

   //Defining CSVv2 bTag Operating Point (LOOSE, MEDIUM, TIGHT OR RESHAPING)  
   std::string CSV_WP = "L"; // required for Tagger (L, M or T)              

   //*********************************************************************************************************//
   //Get Event No., Lumi No. and Run no. for high invariant mass events
   Double_t htmass    = 2000.0;
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
      cout << "metFilters = " << metFilters << endl;


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
      Pass_GJdEta = false;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      GoodVertex = nGoodVtx;
      if(GoodVertex > 0) HasPrimaryVtx = true;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("medium"); //All photons passing loose id, residual spikes and pt > 30

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
      GoodIsoJets = GoodJets(PC, "tight"); // All jets passing dR > 0.5, tight jet id and pt > 30.0 
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) < Cut_Jet_eta);
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
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
         
	if(HasPrimaryVtx && metFilters == 0){
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


   


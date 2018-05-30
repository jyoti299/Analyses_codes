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

   Cut_Photon_pt = 190.0; // GeV
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 190.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 2.5;
   Cut_GJdEta = 2.0;

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
   //HLT.push_back(HLT_Photon175_v);
   //   HLT.push_back(HLT_Photon250_NoHE_v);
   //If some events passing Photon175 and failing 165_HE10 ==> these have more than 10% of their energy in HCAL so these are essentially junk and 
   //should not be considered. So using only 165_HE10 is best option.
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
   std::string CSV_WP = "M"; // required for Tagger (L, M or T)              

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
      Pass_GJdPhi = false;
      Pass_GJdEta = false;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("loose"); //All photons passing loose id, residual spikes and pt > 30

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
      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);

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
      
      h_CutFlow->Fill(0.5);
   
      if(Pass_HLT){
	h_CutFlow->Fill(1.5);
         
	if(HasPrimaryVtx && metFilters == 0){
	  h_CutFlow->Fill(2.5);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow->Fill(3.5);
	  
	    if(PC > -1){
	      h_CutFlow->Fill(4.5);
	      
	      if(JC > -1){
		h_CutFlow->Fill(5.5);
		  
		if(Pass_JetPt){
		  h_CutFlow->Fill(6.5);

		  if(Pass_JetEta){
		    h_CutFlow->Fill(7.5);

		    if(Pass_GJdPhi){
		      h_CutFlow->Fill(8.5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);
		      
			//Photon distributions noBTag_noMasscut
			h_PhotonPt[0][0]               ->Fill((*phoEt)[PC]);
			h_PhotonEta[0][0]              ->Fill((*phoSCEta)[PC]);
			h_PhotonPhi[0][0]              ->Fill((*phoSCPhi)[PC]);
			h_Photon_SigmaIEtaIEta[0][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			h_Photon_R9[0][0]              ->Fill((*phoR9)[PC]);
			h_Photon_HoverE[0][0]          ->Fill((*phoHoverE)[PC]);
			h_Photon_EleVeto[0][0]         ->Fill((*phoEleVeto)[PC]);
			h_Photon_CorrPFChIso[0][0]     ->Fill((*phoPFChIso)[PC]);
			h_Photon_CorrPFNeuIso[0][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			h_Photon_CorrPFPhoIso[0][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
		      
			//Jet distributions noBTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
			h_bJetPt[0][0]                 ->Fill((*jetPt)[JC]);
			h_bJetEta[0][0]                ->Fill((*jetEta)[JC]);
			h_bJetPhi[0][0]                ->Fill((*jetPhi)[JC]);
			h_bJet_NHEF[0][0]              ->Fill((*jetNHF)[JC]);
			h_bJet_NEEF[0][0]              ->Fill((*jetNEF)[JC]);
			h_bJet_NConst[0][0]            ->Fill((*jetNConstituents)[JC]);
			h_bJet_CHEF[0][0]              ->Fill((*jetCHF)[JC]);
			h_bJet_ChMult[0][0]            ->Fill((*jetNCH)[JC]);
			h_bJet_CEEF[0][0]              ->Fill((*jetCEF)[JC]);
		      
			//Photon+Jet distributions noBTag_noMasscut
			h_GbJetInvtMass_VarBin[0][0]   ->Fill(GetInvtMass(PC, JC));
			h_GbJetInvtMass_UnitBin[0][0]  ->Fill(GetInvtMass(PC, JC));
			h_GbJet_dEta[0][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			h_GbJet_dPhi[0][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			h_GbJet_dR[0][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			h_cosThetaStar[0][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));
		      
			//PFMet distributions for noBTag_noMasscut
			h_PFMet[0][0]                  ->Fill(pfMET);
			h_SumPFMet[0][0]               ->Fill(pfMETsumEt);
                        h_MetBySumMET[0][0]            ->Fill(pfMET/pfMETsumEt);
			h_PFMetVsGJmass[0][0]          ->Fill(GetInvtMass(PC, JC), pfMET);
                        h_PFMetOverSumEtVsGJmass[0][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
		      
			//Photon vs Jet dist for noBTag_noMasscut
			h_PhPt_vs_bJetPt[0][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			h_PhEta_vs_bJetEta[0][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
		      
			//CSVv2 discriminator distributions for noBTag_noMasscut
			h_CSVv2Dist[0][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			h_CSVv2_vs_bJetPt[0][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			h_CSVv2_vs_bJetEta[0][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
		      
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
		      
			//Photon and Jet index for noBTag_noMasscut only
			h_PC                          ->Fill(PC);
			h_JC                          ->Fill(JC);
		      
			//--------------------------------------------------------------------------------
			//Getting Photon and Jet distributions after mass cut
			if(Pass_GJInvtMass){

			  h_CutFlow->Fill(10.5);
			
			  //Photon distributions noBTag_Masscut
			  h_PhotonPt[0][1]               ->Fill((*phoEt)[PC]);
			  h_PhotonEta[0][1]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[0][1]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[0][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[0][1]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[0][1]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[0][1]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[0][1]     ->Fill((*phoPFChIso)[PC]);
			  h_Photon_CorrPFNeuIso[0][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[0][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
			
			  //Jet distributions noBTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
			  h_bJetPt[0][1]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[0][1]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[0][1]                ->Fill((*jetPhi)[JC]);
			  h_bJet_NHEF[0][1]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[0][1]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[0][1]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[0][1]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[0][1]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[0][1]              ->Fill((*jetCEF)[JC]);
			
			  //Photon+Jet distributions noBTag_Masscut
			  h_GbJetInvtMass_VarBin[0][1]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[0][1]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[0][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[0][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[0][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[0][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));

			  //PFMet distributions for noBTag_Masscut
			  h_PFMet[0][1]                  ->Fill(pfMET);
			  h_SumPFMet[0][1]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[0][1]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[0][1]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[0][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			  //Photon vs Jet dist for noBTag_Masscut
			  h_PhPt_vs_bJetPt[0][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[0][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
			
			  //CSVv2 discriminator distributions for noBTag_Masscut
			  h_CSVv2Dist[0][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[0][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[0][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);

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
			
			}// Closing if(Pass_GJInvtMass) here so to get btag distributions without mass cut
			//--------------------------------------------------------------------------------

			//Photon and bJet distributions for 1BTag category
			if(Pass_CSVv2bTag){

			  h_CutFlow->Fill(11.5);

			  //Photon distributions 1BTag_noMasscut
			  h_PhotonPt[1][0]               ->Fill((*phoEt)[PC]);
			  h_PhotonEta[1][0]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[1][0]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[1][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[1][0]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[1][0]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[1][0]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[1][0]     ->Fill((*phoPFChIso)[PC]);
			  h_Photon_CorrPFNeuIso[1][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[1][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			  //bJet distributions 1BTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function
			  h_bJetPt[1][0]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[1][0]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[1][0]                ->Fill((*jetPhi)[JC]);
			  h_bJet_NHEF[1][0]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[1][0]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[1][0]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[1][0]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[1][0]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[1][0]              ->Fill((*jetCEF)[JC]);

			  //Photon+bJet distributions 1BTag_noMasscut
			  h_GbJetInvtMass_VarBin[1][0]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[1][0]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[1][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[1][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[1][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[1][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));

			  //PFMet distributions for 1BTag_noMasscut
			  h_PFMet[1][0]                  ->Fill(pfMET);
			  h_SumPFMet[1][0]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[1][0]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[1][0]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[1][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			  //Photon vs Jet dist for 1BTag_noMasscut
			  h_PhPt_vs_bJetPt[1][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[1][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
			
			  //CSVv2 discriminator distributions for 1BTag_noMasscut
			  h_CSVv2Dist[1][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[1][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[1][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);

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

			  //--------------------------------------------------------------------------------
			  //Getting Photon and bJet-1Tag distributions after mass cut
			  if(Pass_GJInvtMass){

			    h_CutFlow->Fill(12.5);

			    //Photon distributions 1BTag_Masscut
			    h_PhotonPt[1][1]               ->Fill((*phoEt)[PC]);
			    h_PhotonEta[1][1]              ->Fill((*phoSCEta)[PC]);
			    h_PhotonPhi[1][1]              ->Fill((*phoSCPhi)[PC]);
			    h_Photon_SigmaIEtaIEta[1][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			    h_Photon_R9[1][1]              ->Fill((*phoR9)[PC]);
			    h_Photon_HoverE[1][1]          ->Fill((*phoHoverE)[PC]);
			    h_Photon_EleVeto[1][1]         ->Fill((*phoEleVeto)[PC]);
			    h_Photon_CorrPFChIso[1][1]     ->Fill((*phoPFChIso)[PC]);
			    h_Photon_CorrPFNeuIso[1][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFPhoIso[1][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			    //bJet distributions 1BTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			    h_bJetPt[1][1]                 ->Fill((*jetPt)[JC]);
			    h_bJetEta[1][1]                ->Fill((*jetEta)[JC]);
			    h_bJetPhi[1][1]                ->Fill((*jetPhi)[JC]);
			    h_bJet_NHEF[1][1]              ->Fill((*jetNHF)[JC]);
			    h_bJet_NEEF[1][1]              ->Fill((*jetNEF)[JC]);
			    h_bJet_NConst[1][1]            ->Fill((*jetNConstituents)[JC]);
			    h_bJet_CHEF[1][1]              ->Fill((*jetCHF)[JC]);
			    h_bJet_ChMult[1][1]            ->Fill((*jetNCH)[JC]);
			    h_bJet_CEEF[1][1]              ->Fill((*jetCEF)[JC]);

			    //Photon+bJet-1Tag distributions 1BTag_Masscut
			    h_GbJetInvtMass_VarBin[1][1]   ->Fill(GetInvtMass(PC, JC));
			    h_GbJetInvtMass_UnitBin[1][1]  ->Fill(GetInvtMass(PC, JC));
			    h_GbJet_dEta[1][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			    h_GbJet_dPhi[1][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_GbJet_dR[1][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_cosThetaStar[1][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));

			    //PFMet distributions for 1BTag_Masscut
			    h_PFMet[1][1]                  ->Fill(pfMET);
			    h_SumPFMet[1][1]               ->Fill(pfMETsumEt);
			    h_MetBySumMET[1][1]            ->Fill(pfMET/pfMETsumEt);
			    h_PFMetVsGJmass[1][1]          ->Fill(GetInvtMass(PC, JC), pfMET);
			    h_PFMetOverSumEtVsGJmass[1][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			    //Photon vs Jet dist for 1BTag_Masscut
			    h_PhPt_vs_bJetPt[1][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			    h_PhEta_vs_bJetEta[1][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);

			    //CSVv2 discriminator distributions for 1BTag_Masscut
			    h_CSVv2Dist[1][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetPt[1][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetEta[1][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);

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

			  }// Closing if(Pass_GJInvtMass) here for 1bTag
			  //--------------------------------------------------------------------------------
			
			}//Closing if(Pass_CSVv2bTag)
		      
			//Photon and bJet distributions for 0BTag category
			else{

			  h_CutFlow->Fill(13.5);

			  //Photon distributions 0BTag_noMasscut
			  h_PhotonPt[2][0]               ->Fill((*phoEt)[PC]);
			  h_PhotonEta[2][0]              ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[2][0]              ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[2][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[2][0]              ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[2][0]          ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[2][0]         ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[2][0]     ->Fill((*phoPFChIso)[PC]);
			  h_Photon_CorrPFNeuIso[2][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[2][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			  //bJet distributions 0BTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			  h_bJetPt[2][0]                 ->Fill((*jetPt)[JC]);
			  h_bJetEta[2][0]                ->Fill((*jetEta)[JC]);
			  h_bJetPhi[2][0]                ->Fill((*jetPhi)[JC]);
			  h_bJet_NHEF[2][0]              ->Fill((*jetNHF)[JC]);
			  h_bJet_NEEF[2][0]              ->Fill((*jetNEF)[JC]);
			  h_bJet_NConst[2][0]            ->Fill((*jetNConstituents)[JC]);
			  h_bJet_CHEF[2][0]              ->Fill((*jetCHF)[JC]);
			  h_bJet_ChMult[2][0]            ->Fill((*jetNCH)[JC]);
			  h_bJet_CEEF[2][0]              ->Fill((*jetCEF)[JC]);

			  //Photon+bJet distributions 0BTag_noMasscut
			  h_GbJetInvtMass_VarBin[2][0]   ->Fill(GetInvtMass(PC, JC));
			  h_GbJetInvtMass_UnitBin[2][0]  ->Fill(GetInvtMass(PC, JC));
			  h_GbJet_dEta[2][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GbJet_dPhi[2][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GbJet_dR[2][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_cosThetaStar[2][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));

			  //PFMet distributions for 0BTag_noMasscut
			  h_PFMet[2][0]                  ->Fill(pfMET);
			  h_SumPFMet[2][0]               ->Fill(pfMETsumEt);
			  h_MetBySumMET[2][0]            ->Fill(pfMET/pfMETsumEt);
			  h_PFMetVsGJmass[2][0]          ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[2][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);
			
			  //Photon vs Jet dist for 0BTag_noMasscut
			  h_PhPt_vs_bJetPt[2][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			  h_PhEta_vs_bJetEta[2][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);

			  //CSVv2 discriminator distributions for 0BTag_noMasscut
			  h_CSVv2Dist[2][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetPt[2][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			  h_CSVv2_vs_bJetEta[2][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);

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
			
			  //--------------------------------------------------------------------------------
			  //Getting Photon and bJet-0Tag distributions after mass cut
			  if(Pass_GJInvtMass){

			    h_CutFlow->Fill(14.5);

			    //Photon distributions 0BTag_Masscut
			    h_PhotonPt[2][1]               ->Fill((*phoEt)[PC]);
			    h_PhotonEta[2][1]              ->Fill((*phoSCEta)[PC]);
			    h_PhotonPhi[2][1]              ->Fill((*phoSCPhi)[PC]);
			    h_Photon_SigmaIEtaIEta[2][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			    h_Photon_R9[2][1]              ->Fill((*phoR9)[PC]);
			    h_Photon_HoverE[2][1]          ->Fill((*phoHoverE)[PC]);
			    h_Photon_EleVeto[2][1]         ->Fill((*phoEleVeto)[PC]);
			    h_Photon_CorrPFChIso[2][1]     ->Fill((*phoPFChIso)[PC]);
			    h_Photon_CorrPFNeuIso[2][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			    h_Photon_CorrPFPhoIso[2][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			    //bJet distributions 0BTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			    h_bJetPt[2][1]                 ->Fill((*jetPt)[JC]);
			    h_bJetEta[2][1]                ->Fill((*jetEta)[JC]);
			    h_bJetPhi[2][1]                ->Fill((*jetPhi)[JC]);
			    h_bJet_NHEF[2][1]              ->Fill((*jetNHF)[JC]);
			    h_bJet_NEEF[2][1]              ->Fill((*jetNEF)[JC]);
			    h_bJet_NConst[2][1]            ->Fill((*jetNConstituents)[JC]);
			    h_bJet_CHEF[2][1]              ->Fill((*jetCHF)[JC]);
			    h_bJet_ChMult[2][1]            ->Fill((*jetNCH)[JC]);
			    h_bJet_CEEF[2][1]              ->Fill((*jetCEF)[JC]);
			  
			    //Photon+bJet distributions 0BTag_Masscut
			    h_GbJetInvtMass_VarBin[2][1]   ->Fill(GetInvtMass(PC, JC));
			    h_GbJetInvtMass_UnitBin[2][1]  ->Fill(GetInvtMass(PC, JC));
			    h_GbJet_dEta[2][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			    h_GbJet_dPhi[2][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_GbJet_dR[2][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));
			    h_cosThetaStar[2][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]));

			    //PFMet distributions for 0BTag_Masscut
			    h_PFMet[2][1]                  ->Fill(pfMET);
			    h_SumPFMet[2][1]               ->Fill(pfMETsumEt);
			    h_MetBySumMET[2][1]            ->Fill(pfMET/pfMETsumEt);
			    h_PFMetVsGJmass[2][1]          ->Fill(GetInvtMass(PC, JC), pfMET);
			    h_PFMetOverSumEtVsGJmass[2][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			    //Photon vs Jet dist for 0BTag_Masscut
			    h_PhPt_vs_bJetPt[2][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC]);
			    h_PhEta_vs_bJetEta[2][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC]);
			  
			    //CSVv2 discriminator distributions for 0BTag_Masscut
			    h_CSVv2Dist[2][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetPt[2][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);
			    h_CSVv2_vs_bJetEta[2][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC]);

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

			  }// Closing if(Pass_GJInvtMass) here for 0BTag_Masscut
			//--------------------------------------------------------------------------------
			}//Closing else(PassCSVv2BTag)		      
		      
		      }//Pass_GJdEta
		    }//Pass_GJdPhi
		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//PC > -1
	  }//PhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
      
      //---------------------------------------------------------------------------------   
      //Trigger Turn-on
      Bool_t PassHLTNum = false;
      Bool_t PassHLTDeno = false;
      Bool_t PassHLT_Pre = false;

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
	  if(PassHLTNum && PassHLT_Pre){
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


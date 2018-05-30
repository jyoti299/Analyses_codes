#define PostAnalyzer_MC_cxx
#include "PostAnalyzer_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//To be removed in script
//-----------------------
int main(){
  PostAnalyzer_MC t;
  t.Loop();
  return 0;
}
//------------------------

void PostAnalyzer_MC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PostAnalyzer_MC.C
//      root> PostAnalyzer_MC t
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
   //Luminosity of Data Compared
   Lumi = 2670.555; //(pb^{-1}), COMPLETE 2015 25ns DATA WITH SILVER JSON AFTER CHANGING REFERENCE FILE IN BIRL CALC (DONE CENTRALLY BY CMS LUMI GRP)

   //To be removed in script
   //-----------------------
   Float_t XS = 5813.0;
   //-----------------------

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
   
   //To be removed in script
   //-----------------------
   //Define Output file here
   file = new TFile("PostAnalyzer_MC.root", "RECREATE");
   //-----------------------

   //Uncomment this in script
   //Define Output file here
   //   TString OutputPath = "${destinationDir}/";
   //   TString OutputFile = "${filenameTag}";
   //  // file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   //   file = new TFile(OutputFile+".root", "RECREATE");

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   PileupReWeighting();

   //Defining CSVc2 bTag Operating Point (LOOSE, MEDIUM, TIGHT OR RESHAPING (for boosted btag discs))
   BTagEntry::OperatingPoint CSV_OP = BTagEntry::OP_MEDIUM; // required for SF calculation
   std::string CSV_WP = "M"; // required for Tagger (L,M or T)

   //Event For loop starts from here
   Long64_t nentries = fChain->GetEntries();
   cout << "<Total entries: " << nentries << endl; 
   Long64_t nbytes = 0, nb = 0;

   /*
   //running for loop to get the total genweight
   Long64_t nb1 = 0;
   Float_t Tot_genWt = 0;
   for (Long64_t bentry=0; bentry<nentries;bentry++) {
     Long64_t kentry = LoadTree(bentry);
     nb1 = fChain->GetEntry(bentry);
     Tot_genWt += genWeight;
   }
   cout << "Tot genWt = " << Tot_genWt << endl;
   //---------------------------------------------
   */

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //For 76X, use nentries in place of totalEvents[${sampleIndex}], however, this can be used for 74X as well. As for these samples, there is no exta cuts at ntuples level. So events in the samples are same as the actual ones.

      //Uncomment this in script  
      //Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/${totalEvents[${sampleIndex}]};

      //To be removed in script                                                                       
      //-----------------------                                                                                                                     
      Lumi_EvtWt = (Lumi*XS)/9956130;//2000069
      //-----------------------                                                                                                 

      PU_EvtWt = PUWeights((*puTrue)[0]); //Since TrueNumofInt is same for an event, so any value of the vector puTrue can be taken. 
                                          //(see definition here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information)
      PreBTag_EvtWt = Lumi_EvtWt * PU_EvtWt;
      
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
      Pass_HLT = true; //Always true for MC
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("loose");  //All photons passing loose id, residual spikes and pt > 30.0

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

      //Primary vertex, trueNumofInt and number of photon and jets for noCut
      h_trueNumofInt       ->Fill((*puTrue)[0]);
      h_goodPV_noWt        ->Fill(GoodVertex);
      h_goodPV_PUWt        ->Fill(GoodVertex, PU_EvtWt);

      h_goodPV_TotalWt[0]              ->Fill(GoodVertex, PreBTag_EvtWt);
      h_nIsoPhotons[0]                 ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
      h_nGoodPhotons[0]                ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of isolated photons with pt > cut and eta < cut 
      for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
	h_IsoPhotonIdxVsPt[0]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
      }
      for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
	h_GoodPhotonIdxVsPt[0]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
      }				     
      h_nJets[0]                       ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
      for(int ij = 0; ij < GoodIsoJets.size(); ij++){
	h_JetIdxVsPt[0]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
      }

      h_CutFlow->Fill(0.5);
      h_CutFlow_PreBTagWts->Fill(0.5, PreBTag_EvtWt);
      h_CutFlow_TotalWts->Fill(0.5, PreBTag_EvtWt);

      if(Pass_HLT){
	h_CutFlow->Fill(1.5);
	h_CutFlow_PreBTagWts->Fill(1.5, PreBTag_EvtWt);
	h_CutFlow_TotalWts->Fill(1.5, PreBTag_EvtWt);

	if(HasPrimaryVtx){
	  h_CutFlow->Fill(2.5);
	  h_CutFlow_PreBTagWts->Fill(2.5, PreBTag_EvtWt);
	  h_CutFlow_TotalWts->Fill(2.5, PreBTag_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow->Fill(3.5);
	    h_CutFlow_PreBTagWts->Fill(3.5, PreBTag_EvtWt);
	    h_CutFlow_TotalWts->Fill(3.5, PreBTag_EvtWt);

	    if(PC > -1){
	      h_CutFlow->Fill(4.5);
	      h_CutFlow_PreBTagWts->Fill(4.5, PreBTag_EvtWt);
	      h_CutFlow_TotalWts->Fill(4.5, PreBTag_EvtWt);

	      if(JC > -1){
		h_CutFlow->Fill(5.5);
		h_CutFlow_PreBTagWts->Fill(5.5, PreBTag_EvtWt);
		h_CutFlow_TotalWts->Fill(5.5, PreBTag_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow->Fill(6.5);
		  h_CutFlow_PreBTagWts->Fill(6.5, PreBTag_EvtWt);
		  h_CutFlow_TotalWts->Fill(6.5, PreBTag_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow->Fill(7.5);
		    h_CutFlow_PreBTagWts->Fill(7.5, PreBTag_EvtWt);
		    h_CutFlow_TotalWts->Fill(7.5, PreBTag_EvtWt);

		    if(Pass_GJdPhi){
		      h_CutFlow->Fill(8.5);
		      h_CutFlow_PreBTagWts->Fill(8.5, PreBTag_EvtWt);
		      h_CutFlow_TotalWts->Fill(8.5, PreBTag_EvtWt);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);
			h_CutFlow_PreBTagWts->Fill(9.5, PreBTag_EvtWt);
			h_CutFlow_TotalWts->Fill(9.5, PreBTag_EvtWt);
		      
			//Photon distributions noBTag_noMasscut
			h_PhotonPt[0][0]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			h_PhotonEta[0][0]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			h_PhotonPhi[0][0]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			h_Photon_SigmaIEtaIEta[0][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			h_Photon_R9[0][0]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			h_Photon_HoverE[0][0]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			h_Photon_EleVeto[0][0]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			h_Photon_CorrPFChIso[0][0]     ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			h_Photon_CorrPFNeuIso[0][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			h_Photon_CorrPFPhoIso[0][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
		      
			//Jet distributions noBTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
			h_bJetPt[0][0]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			h_bJetEta[0][0]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			h_bJetPhi[0][0]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			h_bJet_NHEF[0][0]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			h_bJet_NEEF[0][0]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			h_bJet_NConst[0][0]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			h_bJet_CHEF[0][0]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			h_bJet_ChMult[0][0]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			h_bJet_CEEF[0][0]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
		      
			//Photon+Jet distributions noBTag_noMasscut
			h_GbJetInvtMass_VarBin[0][0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GbJetInvtMass_UnitBin[0][0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GbJet_dEta[0][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			h_GbJet_dPhi[0][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			h_GbJet_dR[0][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			h_cosThetaStar[0][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
		      
			//PFMet distributions for noBTag_noMasscut
			h_PFMet[0][0]                  ->Fill(pfMET, PreBTag_EvtWt);
			h_SumPFMet[0][0]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
                        h_MetBySumMET[0][0]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			h_PFMetVsGJmass[0][0]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
                        h_PFMetOverSumEtVsGJmass[0][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
		      
			//Photon vs Jet dist for noBTag_noMasscut
			h_PhPt_vs_bJetPt[0][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			h_PhEta_vs_bJetEta[0][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
		      
			//CSVv2 discriminator distributions for noBTag_noMasscut
			h_CSVv2Dist[0][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			h_CSVv2_vs_bJetPt[0][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			h_CSVv2_vs_bJetEta[0][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
		      
			//Primary vertex and number of photon and jets for noBTag_noMasscut
			h_goodPV_TotalWt[1]            ->Fill(GoodVertex, PreBTag_EvtWt);
			h_nIsoPhotons[1]               ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			h_nGoodPhotons[1]              ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of iso photons with pt>cut and eta<cut 
			for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			  h_IsoPhotonIdxVsPt[1]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			}
			for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			  h_GoodPhotonIdxVsPt[1]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			}
			h_nJets[1]                     ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			  h_JetIdxVsPt[1]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			}
		      
			//Photon and Jet index for noBTag_noMasscut only
			h_PC                          ->Fill(PC, PreBTag_EvtWt);
			h_JC                          ->Fill(JC, PreBTag_EvtWt);
		      
			//--------------------------------------------------------------------------------
			//Getting Photon and Jet distributions after mass cut
			if(Pass_GJInvtMass){

			  h_CutFlow->Fill(10.5);
			  h_CutFlow_PreBTagWts->Fill(10.5, PreBTag_EvtWt);
			  h_CutFlow_TotalWts->Fill(10.5, PreBTag_EvtWt);
			
			  //Photon distributions noBTag_Masscut
			  h_PhotonPt[0][1]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			  h_PhotonEta[0][1]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			  h_PhotonPhi[0][1]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			  h_Photon_SigmaIEtaIEta[0][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			  h_Photon_R9[0][1]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			  h_Photon_HoverE[0][1]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			  h_Photon_EleVeto[0][1]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFChIso[0][1]     ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFNeuIso[0][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFPhoIso[0][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			
			  //Jet distributions noBTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
			  h_bJetPt[0][1]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			  h_bJetEta[0][1]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			  h_bJetPhi[0][1]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			  h_bJet_NHEF[0][1]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			  h_bJet_NEEF[0][1]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			  h_bJet_NConst[0][1]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			  h_bJet_CHEF[0][1]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			  h_bJet_ChMult[0][1]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			  h_bJet_CEEF[0][1]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			
			  //Photon+Jet distributions noBTag_Masscut
			  h_GbJetInvtMass_VarBin[0][1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJetInvtMass_UnitBin[0][1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJet_dEta[0][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			  h_GbJet_dPhi[0][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_GbJet_dR[0][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_cosThetaStar[0][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);

			  //PFMet distributions for noBTag_Masscut
			  h_PFMet[0][1]                  ->Fill(pfMET, PreBTag_EvtWt);
			  h_SumPFMet[0][1]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			  h_MetBySumMET[0][1]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_PFMetVsGJmass[0][1]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			  h_PFMetOverSumEtVsGJmass[0][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);

			  //Photon vs Jet dist for noBTag_Masscut
			  h_PhPt_vs_bJetPt[0][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			  h_PhEta_vs_bJetEta[0][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
			
			  //CSVv2 discriminator distributions for noBTag_Masscut
			  h_CSVv2Dist[0][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetPt[0][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetEta[0][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  //Primary vertex and number of photon and jets for noBTag_Masscut
			  h_goodPV_TotalWt[2]            ->Fill(GoodVertex, PreBTag_EvtWt);
			  h_nIsoPhotons[2]               ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			  h_nGoodPhotons[2]              ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[2]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[2]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			  }
			  h_nJets[2]                     ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[2]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			  }
						
			}// Closing if(Pass_GJInvtMass) here so to get btag distributions without mass cut
			//--------------------------------------------------------------------------------
		      
			if(Pass_CSVv2bTag){

			  h_CutFlow->Fill(11.5);
			  h_CutFlow_PreBTagWts->Fill(11.5, PreBTag_EvtWt);

			  //Here using PartonFlavor to get the gen level jet information. In general, HadronFlavor should be used.
			  //Use HadronFlavor for 76X. Here no HadronFlavor info in 74X, that's why using PartonFlavor

			  Double_t SF, Wt_1Tag, Wt_0Tag;
			  BTagEntry::JetFlavor JF;
			  std::string sys_type = "central"; //central is required to get scale factors (up and down for uncertainties)
			  if(fabs((*jetPartonID)[JC]) == 5) JF = BTagEntry::FLAV_B; //b
			  if(fabs((*jetPartonID)[JC]) == 4) JF = BTagEntry::FLAV_C; //c
			  if(fabs((*jetPartonID)[JC]) == 1 || fabs((*jetPartonID)[JC]) == 2 || fabs((*jetPartonID)[JC]) == 3 || fabs((*jetPartonID)[JC]) == 21)                                                  JF = BTagEntry::FLAV_UDSG; //u,d,s,g

			  if((*jetPartonID)[JC] == 0) JF = BTagEntry::FLAV_UDSG; //If Parton_ID ==0, Jet has undefined flavor. Means no Jet within 
                                                                                 //dR < 0.4, These jets are mainly soft jets from pileup int.
                                                                                 //On applying, high pt cut, their fraction reduced to < 1%
                                                                                 //So safe to either igonre them or keep in light jet category.
                                                                                 //Read NOTE in the twiki: 
                                                                      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2010OpenData
			                                                         //This problem will not be there for HadronFlavor.
                                                                                 //As hadronflavor = 0 => light jets.

			  SF = CSVv2bTagSF(CSV_OP, JF, sys_type, (*jetPt)[JC], (*jetEta)[JC]);

			  Wt_1Tag = BTagEventWeight(SF, 1); // SF
			  Wt_0Tag = BTagEventWeight(SF, 0); // (1-SF)

			  Total_EvtWt_1tag = PreBTag_EvtWt * Wt_1Tag;
			  Total_EvtWt_0tag = PreBTag_EvtWt * Wt_0Tag;

			  //For PassCSVTag, Both 1 and 0 btag categories will be filled with corresponding weights.
			  //1Btag Category (BTagWt = SF for 1 BTag Category)
			  h_CutFlow_TotalWts->Fill(11.5, Total_EvtWt_1tag);

			  //Photon distributions 1BTag_noMasscut
			  h_PhotonPt[1][0]               ->Fill((*phoEt)[PC], Total_EvtWt_1tag);
			  h_PhotonEta[1][0]              ->Fill((*phoSCEta)[PC], Total_EvtWt_1tag);
			  h_PhotonPhi[1][0]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_1tag);
			  h_Photon_SigmaIEtaIEta[1][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_1tag);
			  h_Photon_R9[1][0]              ->Fill((*phoR9)[PC], Total_EvtWt_1tag);
			  h_Photon_HoverE[1][0]          ->Fill((*phoHoverE)[PC], Total_EvtWt_1tag);
			  h_Photon_EleVeto[1][0]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_1tag);
			  h_Photon_CorrPFChIso[1][0]     ->Fill((*phoPFChIso)[PC], Total_EvtWt_1tag);
			  h_Photon_CorrPFNeuIso[1][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			  h_Photon_CorrPFPhoIso[1][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);

			  //bJet distributions 1BTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function
			  h_bJetPt[1][0]                 ->Fill((*jetPt)[JC], Total_EvtWt_1tag);
			  h_bJetEta[1][0]                ->Fill((*jetEta)[JC], Total_EvtWt_1tag);
			  h_bJetPhi[1][0]                ->Fill((*jetPhi)[JC], Total_EvtWt_1tag);
			  h_bJet_NHEF[1][0]              ->Fill((*jetNHF)[JC], Total_EvtWt_1tag);
			  h_bJet_NEEF[1][0]              ->Fill((*jetNEF)[JC], Total_EvtWt_1tag);
			  h_bJet_NConst[1][0]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_1tag);
			  h_bJet_CHEF[1][0]              ->Fill((*jetCHF)[JC], Total_EvtWt_1tag);
			  h_bJet_ChMult[1][0]            ->Fill((*jetNCH)[JC], Total_EvtWt_1tag);
			  h_bJet_CEEF[1][0]              ->Fill((*jetCEF)[JC], Total_EvtWt_1tag);

			  //Photon+bJet distributions 1BTag_noMasscut
			  h_GbJetInvtMass_VarBin[1][0]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			  h_GbJetInvtMass_UnitBin[1][0]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			  h_GbJet_dEta[1][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
			  h_GbJet_dPhi[1][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			  h_GbJet_dR[1][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			  h_cosThetaStar[1][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);

			  //PFMet distributions for 1BTag_noMasscut
			  h_PFMet[1][0]                  ->Fill(pfMET, Total_EvtWt_1tag);
			  h_SumPFMet[1][0]               ->Fill(pfMETsumEt, Total_EvtWt_1tag);
			  h_MetBySumMET[1][0]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_1tag);
			  h_PFMetVsGJmass[1][0]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_1tag);
			  h_PFMetOverSumEtVsGJmass[1][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_1tag);

			  //Photon vs Jet dist for 1BTag_noMasscut
			  h_PhPt_vs_bJetPt[1][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_1tag);
			  h_PhEta_vs_bJetEta[1][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_1tag);
			
			  //CSVv2 discriminator distributions for 1BTag_noMasscut
			  h_CSVv2Dist[1][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);
			  h_CSVv2_vs_bJetPt[1][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);
			  h_CSVv2_vs_bJetEta[1][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);
			  //Primary vertex and number of photon and jets for 1BTag_noMasscut
			  h_goodPV_TotalWt[3]            ->Fill(GoodVertex, Total_EvtWt_1tag);
			  h_nIsoPhotons[3]               ->Fill(GoodIsoPhotons.size(), Total_EvtWt_1tag);  // Tot # of isolated photons
			  h_nGoodPhotons[3]              ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_1tag); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[3]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_1tag);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[3]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_1tag);
			  }
			  h_nJets[3]                     ->Fill(GoodIsoJets.size(), Total_EvtWt_1tag);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[3]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_1tag);
			  }

			  //0BTag Category (BTagWt = 1-SF for 0Btag category)
			  h_CutFlow_TotalWts->Fill(13.5, Total_EvtWt_0tag);

			  //Photon distributions 0BTag_noMasscut 
			  h_PhotonPt[2][0]               ->Fill((*phoEt)[PC], Total_EvtWt_0tag);
			  h_PhotonEta[2][0]              ->Fill((*phoSCEta)[PC], Total_EvtWt_0tag);
			  h_PhotonPhi[2][0]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_0tag);
			  h_Photon_SigmaIEtaIEta[2][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_0tag);
			  h_Photon_R9[2][0]              ->Fill((*phoR9)[PC], Total_EvtWt_0tag);
			  h_Photon_HoverE[2][0]          ->Fill((*phoHoverE)[PC], Total_EvtWt_0tag);
			  h_Photon_EleVeto[2][0]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_0tag);
			  h_Photon_CorrPFChIso[2][0]     ->Fill((*phoPFChIso)[PC], Total_EvtWt_0tag);
			  h_Photon_CorrPFNeuIso[2][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			  h_Photon_CorrPFPhoIso[2][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			  //bJet distributions 0BTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			  h_bJetPt[2][0]                 ->Fill((*jetPt)[JC], Total_EvtWt_0tag);
			  h_bJetEta[2][0]                ->Fill((*jetEta)[JC], Total_EvtWt_0tag);
			  h_bJetPhi[2][0]                ->Fill((*jetPhi)[JC], Total_EvtWt_0tag);
			  h_bJet_NHEF[2][0]              ->Fill((*jetNHF)[JC], Total_EvtWt_0tag);
			  h_bJet_NEEF[2][0]              ->Fill((*jetNEF)[JC], Total_EvtWt_0tag);
			  h_bJet_NConst[2][0]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_0tag);
			  h_bJet_CHEF[2][0]              ->Fill((*jetCHF)[JC], Total_EvtWt_0tag);
			  h_bJet_ChMult[2][0]            ->Fill((*jetNCH)[JC], Total_EvtWt_0tag);
			  h_bJet_CEEF[2][0]              ->Fill((*jetCEF)[JC], Total_EvtWt_0tag);

			  //Photon+bJet distributions 0BTag_noMasscut
			  h_GbJetInvtMass_VarBin[2][0]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			  h_GbJetInvtMass_UnitBin[2][0]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			  h_GbJet_dEta[2][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
			  h_GbJet_dPhi[2][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			  h_GbJet_dR[2][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			  h_cosThetaStar[2][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);

			  //PFMet distributions for 0BTag_noMasscut
			  h_PFMet[2][0]                  ->Fill(pfMET, Total_EvtWt_0tag);
			  h_SumPFMet[2][0]               ->Fill(pfMETsumEt, Total_EvtWt_0tag);
			  h_MetBySumMET[2][0]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_0tag);
			  h_PFMetVsGJmass[2][0]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_0tag);
			  h_PFMetOverSumEtVsGJmass[2][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_0tag);
			
			  //Photon vs Jet dist for 0BTag_noMasscut
			  h_PhPt_vs_bJetPt[2][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_0tag);
			  h_PhEta_vs_bJetEta[2][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_0tag);

			  //CSVv2 discriminator distributions for 0BTag_noMasscut
			  h_CSVv2Dist[2][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			  h_CSVv2_vs_bJetPt[2][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			  h_CSVv2_vs_bJetEta[2][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			  //Primary vertex and number of photon and jets for 0BTag_noMasscut
			  h_goodPV_TotalWt[5]            ->Fill(GoodVertex, Total_EvtWt_0tag);
			  h_nIsoPhotons[5]               ->Fill(GoodIsoPhotons.size(), Total_EvtWt_0tag);  // Tot # of isolated photons
			  h_nGoodPhotons[5]              ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_0tag); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[5]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_0tag);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[5]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_0tag);
			  }
			  h_nJets[5]                     ->Fill(GoodIsoJets.size(), Total_EvtWt_0tag);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[5]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_0tag);
			  }

			  if(Pass_GJInvtMass){
			    h_CutFlow->Fill(12.5);
			    h_CutFlow_PreBTagWts->Fill(12.5, PreBTag_EvtWt);

			    //1Btag Category (BTagWt = SF for 1 tag category)
			    h_CutFlow_TotalWts->Fill(12.5, Total_EvtWt_1tag);

			    //Photon distributions 1BTag_Masscut
			    h_PhotonPt[1][1]               ->Fill((*phoEt)[PC], Total_EvtWt_1tag);
			    h_PhotonEta[1][1]              ->Fill((*phoSCEta)[PC], Total_EvtWt_1tag);
			    h_PhotonPhi[1][1]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_1tag);
			    h_Photon_SigmaIEtaIEta[1][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_1tag);
			    h_Photon_R9[1][1]              ->Fill((*phoR9)[PC], Total_EvtWt_1tag);
			    h_Photon_HoverE[1][1]          ->Fill((*phoHoverE)[PC], Total_EvtWt_1tag);
			    h_Photon_EleVeto[1][1]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_1tag);
			    h_Photon_CorrPFChIso[1][1]     ->Fill((*phoPFChIso)[PC], Total_EvtWt_1tag);
			    h_Photon_CorrPFNeuIso[1][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			    h_Photon_CorrPFPhoIso[1][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);

			    //bJet distributions 1BTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			    h_bJetPt[1][1]                 ->Fill((*jetPt)[JC], Total_EvtWt_1tag);
			    h_bJetEta[1][1]                ->Fill((*jetEta)[JC], Total_EvtWt_1tag);
			    h_bJetPhi[1][1]                ->Fill((*jetPhi)[JC], Total_EvtWt_1tag);
			    h_bJet_NHEF[1][1]              ->Fill((*jetNHF)[JC], Total_EvtWt_1tag);
			    h_bJet_NEEF[1][1]              ->Fill((*jetNEF)[JC], Total_EvtWt_1tag);
			    h_bJet_NConst[1][1]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_1tag);
			    h_bJet_CHEF[1][1]              ->Fill((*jetCHF)[JC], Total_EvtWt_1tag);
			    h_bJet_ChMult[1][1]            ->Fill((*jetNCH)[JC], Total_EvtWt_1tag);
			    h_bJet_CEEF[1][1]              ->Fill((*jetCEF)[JC], Total_EvtWt_1tag);

			    //Photon+bJet-1Tag distributions 1BTag_Masscut
			    h_GbJetInvtMass_VarBin[1][1]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			    h_GbJetInvtMass_UnitBin[1][1]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			    h_GbJet_dEta[1][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
			    h_GbJet_dPhi[1][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			    h_GbJet_dR[1][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			    h_cosThetaStar[1][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);

			    //PFMet distributions for 1BTag_Masscut
			    h_PFMet[1][1]                  ->Fill(pfMET, Total_EvtWt_1tag);
			    h_SumPFMet[1][1]               ->Fill(pfMETsumEt, Total_EvtWt_1tag);
			    h_MetBySumMET[1][1]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_1tag);
			    h_PFMetVsGJmass[1][1]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_1tag);
			    h_PFMetOverSumEtVsGJmass[1][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_1tag);

			    //Photon vs Jet dist for 1BTag_Masscut
			    h_PhPt_vs_bJetPt[1][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_1tag);
			    h_PhEta_vs_bJetEta[1][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_1tag);

			    //CSVv2 discriminator distributions for 1BTag_Masscut
			    h_CSVv2Dist[1][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);
			    h_CSVv2_vs_bJetPt[1][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);
			    h_CSVv2_vs_bJetEta[1][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_1tag);

			    //Primary vertex and number of photon and jets for 1BTag_Masscut
			    h_goodPV_TotalWt[4]           ->Fill(GoodVertex, Total_EvtWt_1tag);
			    h_nIsoPhotons[4]               ->Fill(GoodIsoPhotons.size(), Total_EvtWt_1tag);  // Tot # of isolated photons
			    h_nGoodPhotons[4]              ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_1tag); // Tot # of iso photons with pt>cut and eta<cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[4]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_1tag);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[4]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_1tag);
			    }
			    h_nJets[4]                     ->Fill(GoodIsoJets.size(), Total_EvtWt_1tag);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[4]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_1tag);
			    }

			    //0BTag Category (BTag Wt = 1-SF for 0tag category)
			    h_CutFlow_TotalWts->Fill(14.5, Total_EvtWt_0tag);

			    //Photon distributions 0BTag_Masscut
			    h_PhotonPt[2][1]               ->Fill((*phoEt)[PC], Total_EvtWt_0tag);
			    h_PhotonEta[2][1]              ->Fill((*phoSCEta)[PC], Total_EvtWt_0tag);
			    h_PhotonPhi[2][1]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_0tag);
			    h_Photon_SigmaIEtaIEta[2][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_0tag);
			    h_Photon_R9[2][1]              ->Fill((*phoR9)[PC], Total_EvtWt_0tag);
			    h_Photon_HoverE[2][1]          ->Fill((*phoHoverE)[PC], Total_EvtWt_0tag);
			    h_Photon_EleVeto[2][1]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_0tag);
			    h_Photon_CorrPFChIso[2][1]     ->Fill((*phoPFChIso)[PC], Total_EvtWt_0tag);
			    h_Photon_CorrPFNeuIso[2][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			    h_Photon_CorrPFPhoIso[2][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);

			    //bJet distributions 0BTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			    h_bJetPt[2][1]                 ->Fill((*jetPt)[JC], Total_EvtWt_0tag);
			    h_bJetEta[2][1]                ->Fill((*jetEta)[JC], Total_EvtWt_0tag);
			    h_bJetPhi[2][1]                ->Fill((*jetPhi)[JC], Total_EvtWt_0tag);
			    h_bJet_NHEF[2][1]              ->Fill((*jetNHF)[JC], Total_EvtWt_0tag);
			    h_bJet_NEEF[2][1]              ->Fill((*jetNEF)[JC], Total_EvtWt_0tag);
			    h_bJet_NConst[2][1]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_0tag);
			    h_bJet_CHEF[2][1]              ->Fill((*jetCHF)[JC], Total_EvtWt_0tag);
			    h_bJet_ChMult[2][1]            ->Fill((*jetNCH)[JC], Total_EvtWt_0tag);
			    h_bJet_CEEF[2][1]              ->Fill((*jetCEF)[JC], Total_EvtWt_0tag);
			  
			    //Photon+bJet distributions 0BTag_Masscut
			    h_GbJetInvtMass_VarBin[2][1]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			    h_GbJetInvtMass_UnitBin[2][1]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			    h_GbJet_dEta[2][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
			    h_GbJet_dPhi[2][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			    h_GbJet_dR[2][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			    h_cosThetaStar[2][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);

			    //PFMet distributions for 0BTag_Masscut
			    h_PFMet[2][1]                  ->Fill(pfMET, Total_EvtWt_0tag);
			    h_SumPFMet[2][1]               ->Fill(pfMETsumEt, Total_EvtWt_0tag);
			    h_MetBySumMET[2][1]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_0tag);
			    h_PFMetVsGJmass[2][1]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_0tag);
			    h_PFMetOverSumEtVsGJmass[2][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_0tag);

			    //Photon vs Jet dist for 0BTag_Masscut
			    h_PhPt_vs_bJetPt[2][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_0tag);
			    h_PhEta_vs_bJetEta[2][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_0tag);
			  
			    //CSVv2 discriminator distributions for 0BTag_Masscut
			    h_CSVv2Dist[2][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			    h_CSVv2_vs_bJetPt[2][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			    h_CSVv2_vs_bJetEta[2][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], Total_EvtWt_0tag);
			    //Primary vertex and number of photon and jets for 0BTag_Masscut
			    h_goodPV_TotalWt[6]           ->Fill(GoodVertex, Total_EvtWt_0tag);
			    h_nIsoPhotons[6]               ->Fill(GoodIsoPhotons.size(), Total_EvtWt_0tag);  // Tot # of isolated photons
			    h_nGoodPhotons[6]              ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_0tag); // Tot # of iso photons with pt>cut and eta<cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[6]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_0tag);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[6]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_0tag);
			    }
			    h_nJets[6]                     ->Fill(GoodIsoJets.size(), Total_EvtWt_0tag);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[6]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_0tag);
			    }

			  } //if(Pass_GJInvtMass)
			} //if(Pass_CSVv2bTag)

			else{ //Events not passing b tag disc
			
			  //BTagWt = 1 for non passing events
			  h_CutFlow->Fill(13.5);
			  h_CutFlow_PreBTagWts->Fill(13.5, PreBTag_EvtWt);
			  h_CutFlow_TotalWts->Fill(13.5, PreBTag_EvtWt);

			  //Photon distributions 0BTag_noMasscut 
			  h_PhotonPt[2][0]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			  h_PhotonEta[2][0]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			  h_PhotonPhi[2][0]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			  h_Photon_SigmaIEtaIEta[2][0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			  h_Photon_R9[2][0]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			  h_Photon_HoverE[2][0]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			  h_Photon_EleVeto[2][0]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFChIso[2][0]     ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFNeuIso[2][0]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFPhoIso[2][0]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  //bJet distributions 0BTag_noMasscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			  h_bJetPt[2][0]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			  h_bJetEta[2][0]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			  h_bJetPhi[2][0]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			  h_bJet_NHEF[2][0]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			  h_bJet_NEEF[2][0]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			  h_bJet_NConst[2][0]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			  h_bJet_CHEF[2][0]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			  h_bJet_ChMult[2][0]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			  h_bJet_CEEF[2][0]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);

			  //Photon+bJet distributions 0BTag_noMasscut
			  h_GbJetInvtMass_VarBin[2][0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJetInvtMass_UnitBin[2][0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJet_dEta[2][0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			  h_GbJet_dPhi[2][0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_GbJet_dR[2][0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_cosThetaStar[2][0]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);

			  //PFMet distributions for 0BTag_noMasscut
			  h_PFMet[2][0]                  ->Fill(pfMET, PreBTag_EvtWt);
			  h_SumPFMet[2][0]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			  h_MetBySumMET[2][0]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_PFMetVsGJmass[2][0]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			  h_PFMetOverSumEtVsGJmass[2][0] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
			
			  //Photon vs Jet dist for 0BTag_noMasscut
			  h_PhPt_vs_bJetPt[2][0]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			  h_PhEta_vs_bJetEta[2][0]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);

			  //CSVv2 discriminator distributions for 0BTag_noMasscut
			  h_CSVv2Dist[2][0]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetPt[2][0]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetEta[2][0]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			  //Primary vertex and number of photon and jets for 0BTag_noMasscut
			  h_goodPV_TotalWt[5]            ->Fill(GoodVertex, PreBTag_EvtWt);
			  h_nIsoPhotons[5]               ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			  h_nGoodPhotons[5]              ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of iso photons with pt>cut and eta<cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[5]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[5]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			  }
			  h_nJets[5]                     ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[5]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			  }

			  if(Pass_GJInvtMass){
			    h_CutFlow->Fill(14.5);
			    h_CutFlow_PreBTagWts->Fill(14.5, PreBTag_EvtWt);
			    h_CutFlow_TotalWts->Fill(14.5, PreBTag_EvtWt);

			    //Photon distributions 0BTag_Masscut
			    h_PhotonPt[2][1]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			    h_PhotonEta[2][1]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			    h_PhotonPhi[2][1]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			    h_Photon_SigmaIEtaIEta[2][1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			    h_Photon_R9[2][1]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			    h_Photon_HoverE[2][1]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			    h_Photon_EleVeto[2][1]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			    h_Photon_CorrPFChIso[2][1]     ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			    h_Photon_CorrPFNeuIso[2][1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			    h_Photon_CorrPFPhoIso[2][1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);

			    //bJet distributions 0BTag_Masscut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms())
			    h_bJetPt[2][1]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			    h_bJetEta[2][1]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			    h_bJetPhi[2][1]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			    h_bJet_NHEF[2][1]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			    h_bJet_NEEF[2][1]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			    h_bJet_NConst[2][1]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			    h_bJet_CHEF[2][1]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			    h_bJet_ChMult[2][1]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			    h_bJet_CEEF[2][1]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			  
			    //Photon+bJet distributions 0BTag_Masscut
			    h_GbJetInvtMass_VarBin[2][1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GbJetInvtMass_UnitBin[2][1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GbJet_dEta[2][1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			    h_GbJet_dPhi[2][1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			    h_GbJet_dR[2][1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			    h_cosThetaStar[2][1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);

			    //PFMet distributions for 0BTag_Masscut
			    h_PFMet[2][1]                  ->Fill(pfMET, PreBTag_EvtWt);
			    h_SumPFMet[2][1]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			    h_MetBySumMET[2][1]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			    h_PFMetVsGJmass[2][1]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			    h_PFMetOverSumEtVsGJmass[2][1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);

			    //Photon vs Jet dist for 0BTag_Masscut
			    h_PhPt_vs_bJetPt[2][1]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			    h_PhEta_vs_bJetEta[2][1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
			  
			    //CSVv2 discriminator distributions for 0BTag_Masscut
			    h_CSVv2Dist[2][1]              ->Fill((*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			    h_CSVv2_vs_bJetPt[2][1]        ->Fill((*jetPt)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);
			    h_CSVv2_vs_bJetEta[2][1]       ->Fill((*jetEta)[JC], (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[JC], PreBTag_EvtWt);

			    //Primary vertex and number of photon and jets for 0BTag_Masscut
			    h_goodPV_TotalWt[6]            ->Fill(GoodVertex, PreBTag_EvtWt);
			    h_nIsoPhotons[6]               ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			    h_nGoodPhotons[6]              ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of iso photons with pt>cut and eta<cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[6]        ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[6]       ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			    }
			    h_nJets[6]                     ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[6]              ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			    }

			  } //if(Pass_GJInvtMass)
			} //else
		      } //if(Pass_GJdEta)
		    }//Pass_GJdPhi
		  } //Pass_JetEta
		} //Pass_JetPt
	      } //JC > -1 
	    }//Pass_PhoPtEta 
	  }//PassPhotonID
	}//HasPrimaryVtx 
      }//Pass_HLT  

      //Mass turn on
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

      //-----------------------------------------------------------------
      //CSVv2 BTag SF Errors

      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(GoodIsoPhotons.size() > 0){
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
		    if(Pass_GJdPhi){
		      if(Pass_GJdEta){
			if(Pass_CSVv2bTag){

			  h_CutFlow_BTagSFerr->Fill(0.5);

			  Double_t SF, SFup, SFdown;
			  Double_t Wt_1Tag_SF, Wt_1Tag_SFup, Wt_1Tag_SFdown;
			  Double_t Wt_0Tag_SF, Wt_0Tag_SFup, Wt_0Tag_SFdown;
			  BTagEntry::JetFlavor JF;
			  std::string sys_type_c = "central"; 
			  std::string sys_type_u = "up"; 
			  std::string sys_type_d = "down"; 

			  if(fabs((*jetPartonID)[JC]) == 5) JF = BTagEntry::FLAV_B; //b
			  if(fabs((*jetPartonID)[JC]) == 4) JF = BTagEntry::FLAV_C; //c
			  if(fabs((*jetPartonID)[JC]) == 1 || fabs((*jetPartonID)[JC]) == 2 || fabs((*jetPartonID)[JC]) == 3 || fabs((*jetPartonID)[JC]) == 21)                                                  JF = BTagEntry::FLAV_UDSG; //u,d,s,g

			  if(fabs((*jetPartonID)[JC]) == 0) JF = BTagEntry::FLAV_UDSG;

			  SF = CSVv2bTagSF(CSV_OP, JF, sys_type_c, (*jetPt)[JC], (*jetEta)[JC]);
			  SFup = CSVv2bTagSF(CSV_OP, JF, sys_type_u, (*jetPt)[JC], (*jetEta)[JC]);
			  SFdown = CSVv2bTagSF(CSV_OP, JF, sys_type_d, (*jetPt)[JC], (*jetEta)[JC]);

			  Wt_1Tag_SF = BTagEventWeight(SF, 1);
			  Wt_1Tag_SFup = BTagEventWeight(SFup, 1);
			  Wt_1Tag_SFdown = BTagEventWeight(SFdown, 1);
			  Wt_0Tag_SF = BTagEventWeight(SF, 0);
			  Wt_0Tag_SFup = BTagEventWeight(SFup, 0);
			  Wt_0Tag_SFdown = BTagEventWeight(SFdown, 0);

			  h_CutFlow_BTagSFerr->Fill(1.5, Wt_1Tag_SF);
			  h_CutFlow_BTagSFerr->Fill(2.5, Wt_1Tag_SFup);
			  h_CutFlow_BTagSFerr->Fill(3.5, Wt_1Tag_SFdown);

			  h_CutFlow_BTagSFerr->Fill(4.5, Wt_0Tag_SF);
			  h_CutFlow_BTagSFerr->Fill(5.5, Wt_0Tag_SFup);
			  h_CutFlow_BTagSFerr->Fill(6.5, Wt_0Tag_SFdown);

			  if(Pass_GJInvtMass){
			    h_CutFlow_BTagSFerr->Fill(7.5);

			    h_CutFlow_BTagSFerr->Fill(8.5, Wt_1Tag_SF);
			    h_CutFlow_BTagSFerr->Fill(9.5, Wt_1Tag_SFup);
			    h_CutFlow_BTagSFerr->Fill(10.5, Wt_1Tag_SFdown);

			    h_CutFlow_BTagSFerr->Fill(11.5, Wt_0Tag_SF);
			    h_CutFlow_BTagSFerr->Fill(12.5, Wt_0Tag_SFup);
			    h_CutFlow_BTagSFerr->Fill(13.5, Wt_0Tag_SFdown);
 
			  }
			} //if(Pass_CSVv2bTag)
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      } //if(Pass_HLT)
      //-----------------------------------------------------------------
 
      }//for jentry
}//Loop()


		      

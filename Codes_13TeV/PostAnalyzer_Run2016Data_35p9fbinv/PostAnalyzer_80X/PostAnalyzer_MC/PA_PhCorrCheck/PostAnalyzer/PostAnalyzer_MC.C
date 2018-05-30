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
   Lumi = 35866.0; //(pb^{-1}) //Re-miniAOD//  (36813.0 pb-1)Rereco-BCDEFG_PromptReco-H  (24487 pb-1)Run2016BCDEFG_PromptReco

   //To be removed in script
   //-----------------------
   Float_t XS = 5813.0;
   //-----------------------

   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 200.0; // GeV                                                                                                             
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 170.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.5;
   
   Cut_GJInvtMass = 700.0;
   
   Cut_PhId = "medium";
   Cut_JetId = "tight";

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
   BTagEntry::OperatingPoint CSV_OP = BTagEntry::OP_LOOSE; // required for SF calculation
   std::string CSV_WP = "L"; // required for Tagger (L,M or T)

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
     cout << "genWt = " << genWeight << endl;
     Tot_genWt += genWeight;
   }
   cout << "Tot genWt = " << Tot_genWt << endl;
   //---------------------------------------------
   */

   //Switches for overlap removal
   Double_t dR_overlap = 0.05;
   Bool_t Remove_QCD = false;
   Bool_t Remove_GJ = false;
   //uncomment in script
   //   if(${QCD}) Remove_QCD = true;
   //   if(${GJ}) Remove_GJ = true;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Uncomment this in script  
      //Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/nentries;

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
      //No dphi cut
      Pass_GJdPhi = true;
      Pass_GJdEta = false;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;
      
      //Running different functions     
      Pass_HLT = true; //Always true for MC
      GoodVertex = nGoodVtx;
      if(GoodVertex > 0) HasPrimaryVtx = true;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons(Cut_PhId);  //All photons passing loose id, residual spikes and pt > 30.0

      GoodIsoBarrelPhotons.clear();
      if(GoodIsoPhotons.size() != 0){
	for(Int_t ip = 0; ip <  GoodIsoPhotons.size(); ip++){
	  if( (*phoCalibEt)[GoodIsoPhotons[ip]] > Cut_Photon_pt && fabs((*phoSCEta)[GoodIsoPhotons[ip]]) < Cut_Photon_eta){
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
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);

      h_CutFlow_qstar->Fill(0.5);
      h_CutFlow_bstar->Fill(0.5);
      h_CutFlowWt_qstar->Fill(0.5, PreBTag_EvtWt);
      h_CutFlowWt_bstar->Fill(0.5, PreBTag_EvtWt);
      h_CutFlowTotalWt_bstar->Fill(0.5, PreBTag_EvtWt);

      //Primary vertex info for noCut only
      h_trueNumofInt       ->Fill((*puTrue)[0]);
      h_goodPV_noWt        ->Fill(GoodVertex);
      h_goodPV_LumiWt      ->Fill(GoodVertex, Lumi_EvtWt);
      h_goodPV_PUWt        ->Fill(GoodVertex, PU_EvtWt);

      //------------------------------------------------------------
      //Photon distributions noCut
      if(nPho > 0 && nJet > 0){
	h_PhotonPt[0]               ->Fill((*phoEt)[0], PreBTag_EvtWt);
	h_PhotonCalibPt[0]          ->Fill((*phoCalibEt)[0], PreBTag_EvtWt);
	h_PhotonEta[0]              ->Fill((*phoSCEta)[0], PreBTag_EvtWt);
	h_PhotonPhi[0]              ->Fill((*phoSCPhi)[0], PreBTag_EvtWt);
	h_Photon_SigmaIEtaIEta[0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[0], PreBTag_EvtWt);
	h_Photon_R9[0]              ->Fill((*phoR9)[0], PreBTag_EvtWt);
	h_Photon_HoverE[0]          ->Fill((*phoHoverE)[0], PreBTag_EvtWt);
	h_Photon_EleVeto[0]         ->Fill((*phoEleVeto)[0], PreBTag_EvtWt);
	h_Photon_CorrPFChIso[0]     ->Fill(TMath::Max(((*phoPFChIso)[0] - rho*EAChargedHadrons((*phoSCEta)[0])), 0.0), PreBTag_EvtWt);
	h_Photon_CorrPFNeuIso[0]    ->Fill(TMath::Max(((*phoPFNeuIso)[0] - rho*EANeutralHadrons((*phoSCEta)[0])), 0.0), PreBTag_EvtWt);
	h_Photon_CorrPFPhoIso[0]    ->Fill(TMath::Max(((*phoPFPhoIso)[0] - rho*EAPhotons((*phoSCEta)[0])), 0.0), PreBTag_EvtWt);
     
	//Jet distributions noCut (Hists named h_bJet to define Jet and bJet dist in one go in BookHistograms() function)
	h_bJetPt[0]                 ->Fill((*jetPt)[0], PreBTag_EvtWt);
	h_bJetEta[0]                ->Fill((*jetEta)[0], PreBTag_EvtWt);
	h_bJetPhi[0]                ->Fill((*jetPhi)[0], PreBTag_EvtWt);
	h_bJet_Mt[0]                ->Fill((*jetMt)[0], PreBTag_EvtWt);
	h_bJet_area[0]              ->Fill((*jetArea)[0], PreBTag_EvtWt);
	h_bJet_Mass[0]              ->Fill((*jetMass)[0], PreBTag_EvtWt);
	h_bJet_NHEF[0]              ->Fill((*jetNHF)[0], PreBTag_EvtWt);
	h_bJet_NEEF[0]              ->Fill((*jetNEF)[0], PreBTag_EvtWt);
        h_bJet_NConst[0]            ->Fill((*jetNConstituents)[0], PreBTag_EvtWt);
	h_bJet_CHEF[0]              ->Fill((*jetCHF)[0], PreBTag_EvtWt);
	h_bJet_ChMult[0]            ->Fill((*jetNCH)[0], PreBTag_EvtWt);
	h_bJet_CEEF[0]              ->Fill((*jetCEF)[0], PreBTag_EvtWt);
	h_bJet_MUF[0]               ->Fill((*jetMUF)[0], PreBTag_EvtWt);
	h_bJet_NNP[0]               ->Fill((*jetNNP)[0], PreBTag_EvtWt);
		      
	//Photon+Jet distributions noCut
	h_GbJetInvtMass_VarBin[0]   ->Fill(GetInvtMass(0, 0), PreBTag_EvtWt);
	h_GbJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(0, 0), PreBTag_EvtWt);
	h_GbJet_dEta[0]             ->Fill(GetdEta((*phoSCEta)[0], (*jetEta)[0]), PreBTag_EvtWt);
	h_GbJet_dPhi[0]             ->Fill(GetdPhi((*phoSCPhi)[0], (*jetPhi)[0]), PreBTag_EvtWt);
	h_GbJet_dR[0]               ->Fill(GetdR((*phoSCEta)[0], (*jetEta)[0], (*phoSCPhi)[0], (*jetPhi)[0]), PreBTag_EvtWt);
	h_cosThetaStar[0]           ->Fill(GetCosThetaStar((*phoSCEta)[0], (*jetEta)[0]), PreBTag_EvtWt);
		      
	//PFMet distributions for noCut
	h_PFMet[0]                  ->Fill(pfMET, PreBTag_EvtWt);
	h_SumPFMet[0]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
	h_MetBySumMET[0]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
	h_PFMetVsGJmass[0]          ->Fill(GetInvtMass(0, 0), pfMET, PreBTag_EvtWt);
	h_PFMetOverSumEtVsGJmass[0] ->Fill(GetInvtMass(0, 0), pfMET/pfMETsumEt, PreBTag_EvtWt);
	h_MetByPhPt[0]              ->Fill(pfMET/(*phoEt)[0], PreBTag_EvtWt);    

	//Photon vs Jet dist for noCut
	h_PhPt_vs_bJetPt[0]         ->Fill((*phoEt)[0], (*jetPt)[0], PreBTag_EvtWt);
	h_PhEta_vs_bJetEta[0]       ->Fill((*phoSCEta)[0], (*jetEta)[0], PreBTag_EvtWt);
		      
	//CSVv2 discriminator distributions for noCut
	h_CSVv2Dist[0]              ->Fill((*jetCSV2BJetTags)[0], PreBTag_EvtWt);
	h_CSVv2_vs_bJetPt[0]        ->Fill((*jetPt)[0], (*jetCSV2BJetTags)[0], PreBTag_EvtWt);
	h_CSVv2_vs_bJetEta[0]       ->Fill((*jetEta)[0], (*jetCSV2BJetTags)[0], PreBTag_EvtWt);

	//Primary vertex and number of photon and jets for noCut
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
      }
      //------------------------------------------------------------
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	h_CutFlow_qstar->Fill(1.5);
	h_CutFlowWt_qstar->Fill(1.5, PreBTag_EvtWt);
         
	if(HasPrimaryVtx){
	  h_CutFlow_qstar->Fill(2.5);
	  h_CutFlowWt_qstar->Fill(2.5, PreBTag_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_qstar->Fill(3.5);
	    h_CutFlowWt_qstar->Fill(3.5, PreBTag_EvtWt);

	    if(PC > -1){
	      h_CutFlow_qstar->Fill(4.5);
	      h_CutFlowWt_qstar->Fill(4.5, PreBTag_EvtWt);

	      //---------------------------------------------------------------------
	      //Overlap removal and gen distributions
	      //QCD/GJ overlap removal (overlap is removed here only as the events is already lost from the loop, so not available at b* level)
	      Int_t PC_match = MatchedPromptGenPhotonToReco(PC); //Matched prompt gen photon with dR < 0.1 and dpt/pt < 0.1   
	      Bool_t rm_QCD = false;
	      Bool_t rm_GJ = false;
	      double minDR = 10.0;

	      if(PC_match > -1){ //true matched photons including overlaps

		for(Int_t ii = 0; ii < nMC; ii++){
		  if(ii == PC_match) continue;
		  if((*mcStatus)[ii] != 22 && (*mcStatus)[ii] != 23) continue;
		  if(fabs((*mcPID)[ii]) > 21) continue;

		  double dR_temp = GetdR((*mcEta)[PC_match], (*mcEta)[ii], (*mcPhi)[PC_match], (*mcPhi)[ii]);
		  if(dR_temp < minDR) minDR = dR_temp;
		}
		if(minDR > dR_overlap)  rm_QCD = true;
		if(minDR < dR_overlap)  rm_GJ = true;
	      }

	      h_matchedPromptGenPhoton_GenParton_dR[0]->Fill(minDR, PreBTag_EvtWt);
	      h_matchedPromptGenPhoton_GenParton_dPToverPT[0]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_match])/(*phoEt)[PC], PreBTag_EvtWt);
	 
	      if(Remove_QCD && rm_QCD) continue; //Removing QCD event
	      if(Remove_GJ && rm_GJ) continue; //Removing GJ event

	      h_matchedPromptGenPhoton_GenParton_dR[1]->Fill(minDR, PreBTag_EvtWt);
	      h_matchedPromptGenPhoton_GenParton_dPToverPT[1]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_match])/(*phoEt)[PC], PreBTag_EvtWt);

	      //dR dists for non prompts
	      Int_t PC_GJmatch = MatchedNonPromptGenPhotonToReco(PC); //Matched non-prompt gen photon with dR < 0.1 and dpt/pt < 0.1
	      double minGJDR = 10.0;

	      if(PC_GJmatch > -1){

		for(Int_t ii = 0; ii < nMC; ii++){
		  if(ii == PC_GJmatch) continue;
		  if((*mcStatus)[ii] != 22 && (*mcStatus)[ii] != 23) continue;
		  if(fabs((*mcPID)[ii]) > 21) continue;

		  double dR_GJtemp = GetdR((*mcEta)[PC_GJmatch], (*mcEta)[ii], (*mcPhi)[PC_GJmatch], (*mcPhi)[ii]);
		  if(dR_GJtemp < minGJDR) minGJDR = dR_GJtemp;
		}
	      }

	      h_matchedNonPromptGenPhoton_GenParton_dR[0]->Fill(minGJDR, PreBTag_EvtWt);
	      h_matchedNonPromptGenPhoton_GenParton_dPToverPT[0]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_GJmatch])/(*phoEt)[PC], PreBTag_EvtWt);

	      h_matchedNonPromptGenPhoton_GenParton_dR[1]->Fill(minGJDR, PreBTag_EvtWt);
	      h_matchedNonPromptGenPhoton_GenParton_dPToverPT[1]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_GJmatch])/(*phoEt)[PC], PreBTag_EvtWt);
	      //---------------------------------------------------------------------

	      h_CutFlow_qstar->Fill(5.5);
	      h_CutFlowWt_qstar->Fill(5.5, PreBTag_EvtWt);

	      if(JC > -1){
		h_CutFlow_qstar->Fill(6.5);
		h_CutFlowWt_qstar->Fill(6.5, PreBTag_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow_qstar->Fill(7.5);
		  h_CutFlowWt_qstar->Fill(7.5, PreBTag_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow_qstar->Fill(8.5);
		    h_CutFlowWt_qstar->Fill(8.5, PreBTag_EvtWt);

 		    if(Pass_GJdPhi){
		      h_CutFlow_qstar->Fill(9.5);
		      h_CutFlowWt_qstar->Fill(9.5, PreBTag_EvtWt);

		      if(Pass_GJdEta){
			h_CutFlow_qstar->Fill(10.5);
			h_CutFlowWt_qstar->Fill(10.5, PreBTag_EvtWt);

			h_matchedPromptGenPhoton_GenParton_dR[2]->Fill(minDR, PreBTag_EvtWt);
			h_matchedPromptGenPhoton_GenParton_dPToverPT[2]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_match])/(*phoEt)[PC], PreBTag_EvtWt);
			h_matchedNonPromptGenPhoton_GenParton_dR[2]->Fill(minGJDR, PreBTag_EvtWt);
			h_matchedNonPromptGenPhoton_GenParton_dPToverPT[2]->Fill(fabs((*phoEt)[PC] - (*mcPt)[PC_GJmatch])/(*phoEt)[PC], PreBTag_EvtWt);

			//Evt Wts
			h_PU_EvtWt->Fill(PU_EvtWt);
			h_PreBTag_EvtWt->Fill(PreBTag_EvtWt);

			//----------------------------------------------------------
			//[1]
			//Photon distributions noBTag_noMasscut
			h_PhotonPt[1]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			h_PhotonCalibPt[1]          ->Fill((*phoCalibEt)[PC], PreBTag_EvtWt);
			h_PhotonEta[1]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			h_PhotonPhi[1]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			h_Photon_SigmaIEtaIEta[1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			h_Photon_R9[1]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			h_Photon_HoverE[1]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			h_Photon_EleVeto[1]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			h_Photon_CorrPFChIso[1]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			h_Photon_CorrPFNeuIso[1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			h_Photon_CorrPFPhoIso[1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
		      
			//Jet distributions noBTag_noMasscut 
			h_bJetPt[1]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			h_bJetEta[1]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			h_bJetPhi[1]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			h_bJet_Mt[1]                ->Fill((*jetMt)[JC], PreBTag_EvtWt);
			h_bJet_area[1]              ->Fill((*jetArea)[JC], PreBTag_EvtWt);
			h_bJet_Mass[1]              ->Fill((*jetMass)[JC], PreBTag_EvtWt);
			h_bJet_NHEF[1]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			h_bJet_NEEF[1]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			h_bJet_NConst[1]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			h_bJet_CHEF[1]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			h_bJet_ChMult[1]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			h_bJet_CEEF[1]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			h_bJet_MUF[1]               ->Fill((*jetMUF)[JC], PreBTag_EvtWt);
			h_bJet_NNP[1]               ->Fill((*jetNNP)[JC], PreBTag_EvtWt);
		      
			//Photon+Jet distributions noBTag_noMasscut 
			h_GbJetInvtMass_VarBin[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GbJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GbJet_dEta[1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			h_GbJet_dPhi[1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			h_GbJet_dR[1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			h_cosThetaStar[1]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
		      
			//PFMet distributions for noBTag_noMasscut 
			h_PFMet[1]                  ->Fill(pfMET, PreBTag_EvtWt);
			h_SumPFMet[1]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			h_MetBySumMET[1]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			h_PFMetVsGJmass[1]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			h_PFMetOverSumEtVsGJmass[1] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
		      	h_MetByPhPt[1]              ->Fill(pfMET/(*phoEt)[PC], PreBTag_EvtWt);    

			//Photon vs Jet dist for noBTag_noMasscut 
			h_PhPt_vs_bJetPt[1]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			h_PhEta_vs_bJetEta[1]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
		      
			//CSVv2 discriminator distributions for noBTag_noMasscut 
			h_CSVv2Dist[1]              ->Fill((*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			h_CSVv2_vs_bJetPt[1]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			h_CSVv2_vs_bJetEta[1]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);

			//Primary vertex and number of photon and jets for noBTag_noMasscut 
			h_goodPV_LumiWt_noMassCut        ->Fill(GoodVertex, Lumi_EvtWt);
			h_goodPV_TotalWt[1]              ->Fill(GoodVertex, PreBTag_EvtWt);
			h_nIsoPhotons[1]                 ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			h_nGoodPhotons[1]                ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of isolated photons with pt > cut and eta < cut 
			for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			  h_IsoPhotonIdxVsPt[1]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			}
			for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			  h_GoodPhotonIdxVsPt[1]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			}				    
			h_nJets[1]                       ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			  h_JetIdxVsPt[1]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			}
			//------------------------------------------------------------

			//Photon and Jet index for noBTag_noMasscut only
			h_PC                          ->Fill(PC);
			h_JC                          ->Fill(JC);
		      
			if(Pass_GJInvtMass){
			  h_CutFlow_qstar->Fill(11.5);
			  h_CutFlowWt_qstar->Fill(11.5, PreBTag_EvtWt);
			  
			  //----------------------------------------------------------
			  //[2]
			  //Photon distributions noBTag_Masscut
			  h_PhotonPt[2]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			  h_PhotonCalibPt[2]          ->Fill((*phoCalibEt)[PC], PreBTag_EvtWt);
			  h_PhotonEta[2]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			  h_PhotonPhi[2]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			  h_Photon_SigmaIEtaIEta[2]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			  h_Photon_R9[2]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			  h_Photon_HoverE[2]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			  h_Photon_EleVeto[2]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFChIso[2]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFNeuIso[2]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFPhoIso[2]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
		      
			  //Jet distributions noBTag_Masscut  
			  h_bJetPt[2]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			  h_bJetEta[2]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			  h_bJetPhi[2]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			  h_bJet_Mt[2]                ->Fill((*jetMt)[JC], PreBTag_EvtWt);
			  h_bJet_area[2]              ->Fill((*jetArea)[JC], PreBTag_EvtWt);
			  h_bJet_Mass[2]              ->Fill((*jetMass)[JC], PreBTag_EvtWt);
			  h_bJet_NHEF[2]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			  h_bJet_NEEF[2]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			  h_bJet_NConst[2]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			  h_bJet_CHEF[2]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			  h_bJet_ChMult[2]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			  h_bJet_CEEF[2]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			  h_bJet_MUF[2]               ->Fill((*jetMUF)[JC], PreBTag_EvtWt);
			  h_bJet_NNP[2]               ->Fill((*jetNNP)[JC], PreBTag_EvtWt);
		      
			  //Photon+Jet distributions noBTag_Masscut 
			  h_GbJetInvtMass_VarBin[2]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJetInvtMass_UnitBin[2]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJet_dEta[2]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			  h_GbJet_dPhi[2]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_GbJet_dR[2]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_cosThetaStar[2]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
		      
			  //PFMet distributions for noBTag_Masscut 
			  h_PFMet[2]                  ->Fill(pfMET, PreBTag_EvtWt);
			  h_SumPFMet[2]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			  h_MetBySumMET[2]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_PFMetVsGJmass[2]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			  h_PFMetOverSumEtVsGJmass[2] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_MetByPhPt[2]              ->Fill(pfMET/(*phoEt)[PC], PreBTag_EvtWt);   

			  //Photon vs Jet dist for noBTag_Masscut 
			  h_PhPt_vs_bJetPt[2]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			  h_PhEta_vs_bJetEta[2]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
		      
			  //CSVv2 discriminator distributions for noBTag_Masscut 
			  h_CSVv2Dist[2]              ->Fill((*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetPt[2]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetEta[2]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);

			  //Primary vertex and number of photon and jets for noBTag_Masscut 
			  h_goodPV_TotalWt[2]              ->Fill(GoodVertex, PreBTag_EvtWt);
			  h_nIsoPhotons[2]                 ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			  h_nGoodPhotons[2]                ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of isolated photons with pt > cut and eta < cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[2]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[2]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			  }				    
			  h_nJets[2]                       ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[2]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			  }
			  //------------------------------------------------------------
			
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
	h_CutFlowWt_bstar->Fill(1.5, PreBTag_EvtWt);
	h_CutFlowTotalWt_bstar->Fill(1.5, PreBTag_EvtWt);
         
	if(HasPrimaryVtx){
	  h_CutFlow_bstar->Fill(2.5);
	  h_CutFlowWt_bstar->Fill(2.5, PreBTag_EvtWt);
	  h_CutFlowTotalWt_bstar->Fill(2.5, PreBTag_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_bstar->Fill(3.5);
	    h_CutFlowWt_bstar->Fill(3.5, PreBTag_EvtWt);
	    h_CutFlowTotalWt_bstar->Fill(3.5, PreBTag_EvtWt);

	    if(PC > -1){
	      h_CutFlow_bstar->Fill(4.5);
	      h_CutFlowWt_bstar->Fill(4.5, PreBTag_EvtWt);
	      h_CutFlowTotalWt_bstar->Fill(4.5, PreBTag_EvtWt);

	      if(JC > -1){
		h_CutFlow_bstar->Fill(5.5);
		h_CutFlowWt_bstar->Fill(5.5, PreBTag_EvtWt);
		h_CutFlowTotalWt_bstar->Fill(5.5, PreBTag_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow_bstar->Fill(6.5);
		  h_CutFlowWt_bstar->Fill(6.5, PreBTag_EvtWt);
		  h_CutFlowTotalWt_bstar->Fill(6.5, PreBTag_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow_bstar->Fill(7.5);
		    h_CutFlowWt_bstar->Fill(7.5, PreBTag_EvtWt);
		    h_CutFlowTotalWt_bstar->Fill(7.5, PreBTag_EvtWt);

		    if(Pass_CSVv2bTag){
		      h_CutFlow_bstar->Fill(8.5);
		      h_CutFlowWt_bstar->Fill(8.5, PreBTag_EvtWt);

		      Double_t SF, Wt_1Tag, Wt_0Tag;
		      Double_t SFerr;
		      BTagEntry::JetFlavor JF;
		      std::string sys_type = "central"; //central is required to get scale factors (up and down for uncertainties) 
		      if(fabs((*jetHadFlvr)[JC]) == 5) JF = BTagEntry::FLAV_B; //b
		      if(fabs((*jetHadFlvr)[JC]) == 4) JF = BTagEntry::FLAV_C; //c
		      if(fabs((*jetHadFlvr)[JC]) == 0) JF = BTagEntry::FLAV_UDSG; //u,d,s,g,undefined

		      SF = CSVv2bTagSF_auto(CSV_OP, JF, sys_type, (*jetPt)[JC], (*jetEta)[JC]);
		      SFerr = CSVv2bTagSF_auto(CSV_OP, JF, "up", (*jetPt)[JC], (*jetEta)[JC]);

		      h_bTag_SF->Fill(SF);

		      if(JF == BTagEntry::FLAV_B){
			h_BTagSF_vs_pt[0]->Fill((*jetPt)[JC], SF, PreBTag_EvtWt);
			h_BTagSFerr_vs_pt[0]->Fill((*jetPt)[JC], SFerr-SF, PreBTag_EvtWt);
		      }
		      if(JF == BTagEntry::FLAV_C){
			h_BTagSF_vs_pt[1]->Fill((*jetPt)[JC], SF, PreBTag_EvtWt);
			h_BTagSFerr_vs_pt[1]->Fill((*jetPt)[JC], SFerr-SF, PreBTag_EvtWt);
		      }
		      if(JF == BTagEntry::FLAV_UDSG){
			h_BTagSF_vs_pt[2]->Fill((*jetPt)[JC], SF, PreBTag_EvtWt);
			h_BTagSFerr_vs_pt[2]->Fill((*jetPt)[JC], SFerr-SF, PreBTag_EvtWt);
		      }

		      Wt_1Tag = BTagEventWeight(SF, 1); // SF
		      Wt_0Tag = BTagEventWeight(SF, 0); // (1-SF)  

		      h_bTag_EvtWt_1Tag->Fill(Wt_1Tag);
		      h_bTag_EvtWt_0Tag->Fill(Wt_0Tag);

		      Total_EvtWt_1tag = PreBTag_EvtWt * Wt_1Tag;
		      Total_EvtWt_0tag = PreBTag_EvtWt * Wt_0Tag;

		      h_CutFlowTotalWt_bstar->Fill(8.5, Total_EvtWt_1tag);
		      h_CutFlowTotalWt_bstar->Fill(12.5, Total_EvtWt_0tag);

		      if(Pass_GJdPhi){
			h_CutFlow_bstar->Fill(9.5);
			h_CutFlowWt_bstar->Fill(9.5, PreBTag_EvtWt);
			h_CutFlowTotalWt_bstar->Fill(9.5, Total_EvtWt_1tag);
			h_CutFlowTotalWt_bstar->Fill(13.5, Total_EvtWt_0tag);

			if(Pass_GJdEta){
			  h_CutFlow_bstar->Fill(10.5);
			  h_CutFlowWt_bstar->Fill(10.5, PreBTag_EvtWt);
			  h_CutFlowTotalWt_bstar->Fill(10.5, Total_EvtWt_1tag);
			  h_CutFlowTotalWt_bstar->Fill(14.5, Total_EvtWt_0tag);

			  //For PassCSVv2Tag, Both 1 and 0 btag categories will be filled with corresponding weights.
			  //1Btag Category (BTagWt = SF for 1 BTag Category)
			  //----------------------------------------------------------
			  //[3]
			  //Photon distributions 1BTag_noMasscut
			  h_PhotonPt[3]               ->Fill((*phoEt)[PC], Total_EvtWt_1tag);
			  h_PhotonCalibPt[3]          ->Fill((*phoCalibEt)[PC], Total_EvtWt_1tag);
			  h_PhotonEta[3]              ->Fill((*phoSCEta)[PC], Total_EvtWt_1tag);
			  h_PhotonPhi[3]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_1tag);
			  h_Photon_SigmaIEtaIEta[3]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_1tag);
			  h_Photon_R9[3]              ->Fill((*phoR9)[PC], Total_EvtWt_1tag);
			  h_Photon_HoverE[3]          ->Fill((*phoHoverE)[PC], Total_EvtWt_1tag);
			  h_Photon_EleVeto[3]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_1tag);
			  h_Photon_CorrPFChIso[3]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			  h_Photon_CorrPFNeuIso[3]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			  h_Photon_CorrPFPhoIso[3]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);		      

			  //Jet distributions 1BTag_noMasscut 
			  h_bJetPt[3]                 ->Fill((*jetPt)[JC], Total_EvtWt_1tag);
			  h_bJetEta[3]                ->Fill((*jetEta)[JC], Total_EvtWt_1tag);
			  h_bJetPhi[3]                ->Fill((*jetPhi)[JC], Total_EvtWt_1tag);
			  h_bJet_Mt[3]                ->Fill((*jetMt)[JC], Total_EvtWt_1tag);
			  h_bJet_area[3]              ->Fill((*jetArea)[JC], Total_EvtWt_1tag);
			  h_bJet_Mass[3]              ->Fill((*jetMass)[JC], Total_EvtWt_1tag);
			  h_bJet_NHEF[3]              ->Fill((*jetNHF)[JC], Total_EvtWt_1tag);
			  h_bJet_NEEF[3]              ->Fill((*jetNEF)[JC], Total_EvtWt_1tag);
			  h_bJet_NConst[3]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_1tag);
			  h_bJet_CHEF[3]              ->Fill((*jetCHF)[JC], Total_EvtWt_1tag);
			  h_bJet_ChMult[3]            ->Fill((*jetNCH)[JC], Total_EvtWt_1tag);
			  h_bJet_CEEF[3]              ->Fill((*jetCEF)[JC], Total_EvtWt_1tag);
			  h_bJet_MUF[3]               ->Fill((*jetMUF)[JC], Total_EvtWt_1tag);
			  h_bJet_NNP[3]               ->Fill((*jetNNP)[JC], Total_EvtWt_1tag);
		      
			  //Photon+Jet distributions 1BTag_noMasscut 
			  h_GbJetInvtMass_VarBin[3]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			  h_GbJetInvtMass_UnitBin[3]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			  h_GbJet_dEta[3]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
			  h_GbJet_dPhi[3]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			  h_GbJet_dR[3]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			  h_cosThetaStar[3]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
		      
			  //PFMet distributions for 1BTag_noMasscut 
			  h_PFMet[3]                  ->Fill(pfMET, Total_EvtWt_1tag);
			  h_SumPFMet[3]               ->Fill(pfMETsumEt, Total_EvtWt_1tag);
			  h_MetBySumMET[3]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_1tag);
			  h_PFMetVsGJmass[3]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_1tag);
			  h_PFMetOverSumEtVsGJmass[3] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_1tag);
			  h_MetByPhPt[3]              ->Fill(pfMET/(*phoEt)[PC], Total_EvtWt_1tag);   

			  //Photon vs Jet dist for 1BTag_noMasscut 
			  h_PhPt_vs_bJetPt[3]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_1tag);
			  h_PhEta_vs_bJetEta[3]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_1tag);
		      
			  //CSVv2 discriminator distributions for 1BTag_noMasscut 
			  h_CSVv2Dist[3]              ->Fill((*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);
			  h_CSVv2_vs_bJetPt[3]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);
			  h_CSVv2_vs_bJetEta[3]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);

			  //Primary vertex and number of photon and jets for 1BTag_noMasscut 
			  h_goodPV_TotalWt[3]              ->Fill(GoodVertex, Total_EvtWt_1tag);
			  h_nIsoPhotons[3]                 ->Fill(GoodIsoPhotons.size(), Total_EvtWt_1tag);  // Tot # of isolated photons
			  h_nGoodPhotons[3]                ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_1tag); // Tot # of isolated photons with pt > cut and eta < cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[3]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_1tag);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[3]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_1tag);
			  }				    
			  h_nJets[3]                       ->Fill(GoodIsoJets.size(), Total_EvtWt_1tag);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[3]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_1tag);
			  }
			  //------------------------------------------------------------

			  //0BTag Category (BTagWt = 1-SF for 0Btag category)
			  //Photon distributions 0BTag_noMasscut
			  h_PhotonPt[5]               ->Fill((*phoEt)[PC], Total_EvtWt_0tag);
			  h_PhotonCalibPt[5]          ->Fill((*phoCalibEt)[PC], Total_EvtWt_0tag);
			  h_PhotonEta[5]              ->Fill((*phoSCEta)[PC], Total_EvtWt_0tag);
			  h_PhotonPhi[5]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_0tag);
			  h_Photon_SigmaIEtaIEta[5]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_0tag);
			  h_Photon_R9[5]              ->Fill((*phoR9)[PC], Total_EvtWt_0tag);
			  h_Photon_HoverE[5]          ->Fill((*phoHoverE)[PC], Total_EvtWt_0tag);
			  h_Photon_EleVeto[5]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_0tag);
			  h_Photon_CorrPFChIso[5]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			  h_Photon_CorrPFNeuIso[5]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			  h_Photon_CorrPFPhoIso[5]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
		      
			  //Jet distributions 0BTag_noMasscut 
			  h_bJetPt[5]                 ->Fill((*jetPt)[JC], Total_EvtWt_0tag);
			  h_bJetEta[5]                ->Fill((*jetEta)[JC], Total_EvtWt_0tag);
			  h_bJetPhi[5]                ->Fill((*jetPhi)[JC], Total_EvtWt_0tag);
			  h_bJet_Mt[5]                ->Fill((*jetMt)[JC], Total_EvtWt_0tag);
			  h_bJet_area[5]              ->Fill((*jetArea)[JC], Total_EvtWt_0tag);
			  h_bJet_Mass[5]              ->Fill((*jetMass)[JC], Total_EvtWt_0tag);
			  h_bJet_NHEF[5]              ->Fill((*jetNHF)[JC], Total_EvtWt_0tag);
			  h_bJet_NEEF[5]              ->Fill((*jetNEF)[JC], Total_EvtWt_0tag);
			  h_bJet_NConst[5]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_0tag);
			  h_bJet_CHEF[5]              ->Fill((*jetCHF)[JC], Total_EvtWt_0tag);
			  h_bJet_ChMult[5]            ->Fill((*jetNCH)[JC], Total_EvtWt_0tag);
			  h_bJet_CEEF[5]              ->Fill((*jetCEF)[JC], Total_EvtWt_0tag);
			  h_bJet_MUF[5]               ->Fill((*jetMUF)[JC], Total_EvtWt_0tag);
			  h_bJet_NNP[5]               ->Fill((*jetNNP)[JC], Total_EvtWt_0tag);
		      
			  //Photon+Jet distributions 0BTag_noMasscut 
			  h_GbJetInvtMass_VarBin[5]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			  h_GbJetInvtMass_UnitBin[5]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			  h_GbJet_dEta[5]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
			  h_GbJet_dPhi[5]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			  h_GbJet_dR[5]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			  h_cosThetaStar[5]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
		      
			  //PFMet distributions for 0BTag_noMasscut 
			  h_PFMet[5]                  ->Fill(pfMET, Total_EvtWt_0tag);
			  h_SumPFMet[5]               ->Fill(pfMETsumEt, Total_EvtWt_0tag);
			  h_MetBySumMET[5]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_0tag);
			  h_PFMetVsGJmass[5]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_0tag);
			  h_PFMetOverSumEtVsGJmass[5] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_0tag);
			  h_MetByPhPt[5]              ->Fill(pfMET/(*phoEt)[PC], Total_EvtWt_0tag);  

			  //Photon vs Jet dist for 0BTag_noMasscut 
			  h_PhPt_vs_bJetPt[5]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_0tag);
			  h_PhEta_vs_bJetEta[5]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_0tag);
		      
			  //CSVv2 discriminator distributions for 0BTag_noMasscut 
			  h_CSVv2Dist[5]              ->Fill((*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);
			  h_CSVv2_vs_bJetPt[5]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);
			  h_CSVv2_vs_bJetEta[5]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);

			  //Primary vertex and number of photon and jets for 0BTag_noMasscut 
			  h_goodPV_TotalWt[5]              ->Fill(GoodVertex, Total_EvtWt_0tag);
			  h_nIsoPhotons[5]                 ->Fill(GoodIsoPhotons.size(), Total_EvtWt_0tag);  // Tot # of isolated photons
			  h_nGoodPhotons[5]                ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_0tag); // Tot # of isolated photons with pt > cut and eta < cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[5]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_0tag);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[5]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_0tag);
			  }				    
			  h_nJets[5]                       ->Fill(GoodIsoJets.size(), Total_EvtWt_0tag);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[5]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_0tag);
			  }
			  //------------------------------------------------------------

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(11.5);
			    h_CutFlowWt_bstar->Fill(11.5, PreBTag_EvtWt);
			    h_CutFlowTotalWt_bstar->Fill(11.5, Total_EvtWt_1tag);
			    h_CutFlowTotalWt_bstar->Fill(15.5, Total_EvtWt_0tag);

				
			    //1Btag Category (BTagWt = SF for 1 BTag Category)
			    //----------------------------------------------------------
			    //[4]
			    //Photon distributions 1BTag_Masscut
			    h_PhotonPt[4]               ->Fill((*phoEt)[PC], Total_EvtWt_1tag);
			    h_PhotonCalibPt[4]          ->Fill((*phoCalibEt)[PC], Total_EvtWt_1tag);
			    h_PhotonEta[4]              ->Fill((*phoSCEta)[PC], Total_EvtWt_1tag);
			    h_PhotonPhi[4]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_1tag);
			    h_Photon_SigmaIEtaIEta[4]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_1tag);
			    h_Photon_R9[4]              ->Fill((*phoR9)[PC], Total_EvtWt_1tag);
			    h_Photon_HoverE[4]          ->Fill((*phoHoverE)[PC], Total_EvtWt_1tag);
			    h_Photon_EleVeto[4]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_1tag);
			    h_Photon_CorrPFChIso[4]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			    h_Photon_CorrPFNeuIso[4]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
			    h_Photon_CorrPFPhoIso[4]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_1tag);
		      
			    //Jet distributions 1BTag_Masscut 
			    h_bJetPt[4]                 ->Fill((*jetPt)[JC], Total_EvtWt_1tag);
			    h_bJetEta[4]                ->Fill((*jetEta)[JC], Total_EvtWt_1tag);
			    h_bJetPhi[4]                ->Fill((*jetPhi)[JC], Total_EvtWt_1tag);
			    h_bJet_Mt[4]                ->Fill((*jetMt)[JC], Total_EvtWt_1tag);
			    h_bJet_area[4]              ->Fill((*jetArea)[JC], Total_EvtWt_1tag);
			    h_bJet_Mass[4]              ->Fill((*jetMass)[JC], Total_EvtWt_1tag);
			    h_bJet_NHEF[4]              ->Fill((*jetNHF)[JC], Total_EvtWt_1tag);
			    h_bJet_NEEF[4]              ->Fill((*jetNEF)[JC], Total_EvtWt_1tag);
			    h_bJet_NConst[4]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_1tag);
			    h_bJet_CHEF[4]              ->Fill((*jetCHF)[JC], Total_EvtWt_1tag);
			    h_bJet_ChMult[4]            ->Fill((*jetNCH)[JC], Total_EvtWt_1tag);
			    h_bJet_CEEF[4]              ->Fill((*jetCEF)[JC], Total_EvtWt_1tag);
			    h_bJet_MUF[4]               ->Fill((*jetMUF)[JC], Total_EvtWt_1tag);
			    h_bJet_NNP[4]               ->Fill((*jetNNP)[JC], Total_EvtWt_1tag);
			    
			    //Photon+Jet distributions 1BTag_Masscut 
			    h_GbJetInvtMass_VarBin[4]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			    h_GbJetInvtMass_UnitBin[4]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			    h_GbJet_dEta[4]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
			    h_GbJet_dPhi[4]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			    h_GbJet_dR[4]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_1tag);
			    h_cosThetaStar[4]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_1tag);
		     
			    //PFMet distributions for 1BTag_Masscut 
			    h_PFMet[4]                  ->Fill(pfMET, Total_EvtWt_1tag);
			    h_SumPFMet[4]               ->Fill(pfMETsumEt, Total_EvtWt_1tag);
			    h_MetBySumMET[4]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_1tag);
			    h_PFMetVsGJmass[4]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_1tag);
			    h_PFMetOverSumEtVsGJmass[4] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_1tag);
			    h_MetByPhPt[4]              ->Fill(pfMET/(*phoEt)[PC], Total_EvtWt_1tag);   

			    //Photon vs Jet dist for 1BTag_Masscut 
			    h_PhPt_vs_bJetPt[4]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_1tag);
			    h_PhEta_vs_bJetEta[4]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_1tag);
		      
			    //CSVv2 discriminator distributions for 1BTag_Masscut 
			    h_CSVv2Dist[4]              ->Fill((*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);
			    h_CSVv2_vs_bJetPt[4]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);
			    h_CSVv2_vs_bJetEta[4]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_1tag);

			    //Primary vertex and number of photon and jets for 1BTag_Masscut 
			    h_goodPV_TotalWt[4]              ->Fill(GoodVertex, Total_EvtWt_1tag);
			    h_nIsoPhotons[4]                 ->Fill(GoodIsoPhotons.size(), Total_EvtWt_1tag);  // Tot # of isolated photons
			    h_nGoodPhotons[4]                ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_1tag); // Tot # of isolated photons with pt > cut and eta < cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[4]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_1tag);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[4]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_1tag);
			    }				   
			    h_nJets[4]                       ->Fill(GoodIsoJets.size(), Total_EvtWt_1tag);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[4]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_1tag);
			    }
			    //------------------------------------------------------------

			    //0BTag Category (BTagWt = 1-SF for 0Btag category)
			    //Photon distributions 0BTag_Masscut
			    h_PhotonPt[6]               ->Fill((*phoEt)[PC], Total_EvtWt_0tag);
			    h_PhotonCalibPt[6]          ->Fill((*phoCalibEt)[PC], Total_EvtWt_0tag);
			    h_PhotonEta[6]              ->Fill((*phoSCEta)[PC], Total_EvtWt_0tag);
			    h_PhotonPhi[6]              ->Fill((*phoSCPhi)[PC], Total_EvtWt_0tag);
			    h_Photon_SigmaIEtaIEta[6]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt_0tag);
			    h_Photon_R9[6]              ->Fill((*phoR9)[PC], Total_EvtWt_0tag);
			    h_Photon_HoverE[6]          ->Fill((*phoHoverE)[PC], Total_EvtWt_0tag);
			    h_Photon_EleVeto[6]         ->Fill((*phoEleVeto)[PC], Total_EvtWt_0tag);
			    h_Photon_CorrPFChIso[6]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			    h_Photon_CorrPFNeuIso[6]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
			    h_Photon_CorrPFPhoIso[6]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt_0tag);
		      
			    //Jet distributions 0BTag_Masscut 
			    h_bJetPt[6]                 ->Fill((*jetPt)[JC], Total_EvtWt_0tag);
			    h_bJetEta[6]                ->Fill((*jetEta)[JC], Total_EvtWt_0tag);
			    h_bJetPhi[6]                ->Fill((*jetPhi)[JC], Total_EvtWt_0tag);
			    h_bJet_Mt[6]                ->Fill((*jetMt)[JC], Total_EvtWt_0tag);
			    h_bJet_area[6]              ->Fill((*jetArea)[JC], Total_EvtWt_0tag);
			    h_bJet_Mass[6]              ->Fill((*jetMass)[JC], Total_EvtWt_0tag);
			    h_bJet_NHEF[6]              ->Fill((*jetNHF)[JC], Total_EvtWt_0tag);
			    h_bJet_NEEF[6]              ->Fill((*jetNEF)[JC], Total_EvtWt_0tag);
			    h_bJet_NConst[6]            ->Fill((*jetNConstituents)[JC], Total_EvtWt_0tag);
			    h_bJet_CHEF[6]              ->Fill((*jetCHF)[JC], Total_EvtWt_0tag);
			    h_bJet_ChMult[6]            ->Fill((*jetNCH)[JC], Total_EvtWt_0tag);
			    h_bJet_CEEF[6]              ->Fill((*jetCEF)[JC], Total_EvtWt_0tag);
			    h_bJet_MUF[6]               ->Fill((*jetMUF)[JC], Total_EvtWt_0tag);
			    h_bJet_NNP[6]               ->Fill((*jetNNP)[JC], Total_EvtWt_0tag);
		      
			    //Photon+Jet distributions 0BTag_Masscut 
			    h_GbJetInvtMass_VarBin[6]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			    h_GbJetInvtMass_UnitBin[6]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);
			    h_GbJet_dEta[6]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
			    h_GbJet_dPhi[6]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			    h_GbJet_dR[6]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt_0tag);
			    h_cosThetaStar[6]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt_0tag);
		      
			    //PFMet distributions for 0BTag_Masscut 
			    h_PFMet[6]                  ->Fill(pfMET, Total_EvtWt_0tag);
			    h_SumPFMet[6]               ->Fill(pfMETsumEt, Total_EvtWt_0tag);
			    h_MetBySumMET[6]            ->Fill(pfMET/pfMETsumEt, Total_EvtWt_0tag);
			    h_PFMetVsGJmass[6]          ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt_0tag);
			    h_PFMetOverSumEtVsGJmass[6] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt_0tag);
			    h_MetByPhPt[6]              ->Fill(pfMET/(*phoEt)[PC], Total_EvtWt_0tag);   

			    //Photon vs Jet dist for 0BTag_Masscut 
			    h_PhPt_vs_bJetPt[6]         ->Fill((*phoEt)[PC], (*jetPt)[JC], Total_EvtWt_0tag);
			    h_PhEta_vs_bJetEta[6]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], Total_EvtWt_0tag);
		      
			    //CSVv2 discriminator distributions for 0BTag_Masscut 
			    h_CSVv2Dist[6]              ->Fill((*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);
			    h_CSVv2_vs_bJetPt[6]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);
			    h_CSVv2_vs_bJetEta[6]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], Total_EvtWt_0tag);

			    //Primary vertex and number of photon and jets for 0BTag_Masscut 
			    h_goodPV_TotalWt[6]              ->Fill(GoodVertex, Total_EvtWt_0tag);
			    h_nIsoPhotons[6]                 ->Fill(GoodIsoPhotons.size(), Total_EvtWt_0tag);  // Tot # of isolated photons
			    h_nGoodPhotons[6]                ->Fill(GoodIsoBarrelPhotons.size(), Total_EvtWt_0tag); // Tot # of isolated photons with pt > cut and eta < cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[6]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, Total_EvtWt_0tag);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[6]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, Total_EvtWt_0tag);
			    }				    
			    h_nJets[6]                       ->Fill(GoodIsoJets.size(), Total_EvtWt_0tag);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[6]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, Total_EvtWt_0tag);
			    }
			    //------------------------------------------------------------

			  }//if(Pass_GJInvtMass) inside if(Pass_CSVv2bTag)
			}//if(Pass_GJdEta) inside if(Pass_CSVv2bTag)
		      }//if(Pass_GJdPhi) inside if(Pass_CSVv2bTag)		       

		    }//if(Pass_CSVv2bTag)		     
		    else{
		      h_CutFlow_bstar->Fill(12.5);
		      h_CutFlowWt_bstar->Fill(12.5, PreBTag_EvtWt);
		      h_CutFlowTotalWt_bstar->Fill(12.5, PreBTag_EvtWt);

		      if(Pass_GJdPhi){
			h_CutFlow_bstar->Fill(13.5);
			h_CutFlowWt_bstar->Fill(13.5, PreBTag_EvtWt);
			h_CutFlowTotalWt_bstar->Fill(13.5, PreBTag_EvtWt);

			if(Pass_GJdEta){
			  h_CutFlow_bstar->Fill(14.5);
			  h_CutFlowWt_bstar->Fill(14.5, PreBTag_EvtWt);
			  h_CutFlowTotalWt_bstar->Fill(14.5, PreBTag_EvtWt);

			  //----------------------------------------------------------
			  //[5]
			  //Photon distributions 0BTag_noMasscut
			  h_PhotonPt[5]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			  h_PhotonCalibPt[5]          ->Fill((*phoCalibEt)[PC], PreBTag_EvtWt);
			  h_PhotonEta[5]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			  h_PhotonPhi[5]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			  h_Photon_SigmaIEtaIEta[5]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			  h_Photon_R9[5]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			  h_Photon_HoverE[5]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			  h_Photon_EleVeto[5]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFChIso[5]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFNeuIso[5]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFPhoIso[5]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
		      
			  //Jet distributions 0BTag_noMasscut  
			  h_bJetPt[5]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			  h_bJetEta[5]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			  h_bJetPhi[5]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			  h_bJet_Mt[5]                ->Fill((*jetMt)[JC], PreBTag_EvtWt);
			  h_bJet_area[5]              ->Fill((*jetArea)[JC], PreBTag_EvtWt);
			  h_bJet_Mass[5]              ->Fill((*jetMass)[JC], PreBTag_EvtWt);
			  h_bJet_NHEF[5]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			  h_bJet_NEEF[5]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			  h_bJet_NConst[5]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			  h_bJet_CHEF[5]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			  h_bJet_ChMult[5]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			  h_bJet_CEEF[5]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			  h_bJet_MUF[5]               ->Fill((*jetMUF)[JC], PreBTag_EvtWt);
			  h_bJet_NNP[5]               ->Fill((*jetNNP)[JC], PreBTag_EvtWt);
		      
			  //Photon+Jet distributions 0BTag_noMasscut 
			  h_GbJetInvtMass_VarBin[5]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJetInvtMass_UnitBin[5]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GbJet_dEta[5]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			  h_GbJet_dPhi[5]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_GbJet_dR[5]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_cosThetaStar[5]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
		      
			  //PFMet distributions for 0BTag_noMasscut 
			  h_PFMet[5]                  ->Fill(pfMET, PreBTag_EvtWt);
			  h_SumPFMet[5]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			  h_MetBySumMET[5]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_PFMetVsGJmass[5]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			  h_PFMetOverSumEtVsGJmass[5] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
			  h_MetByPhPt[5]              ->Fill(pfMET/(*phoEt)[PC], PreBTag_EvtWt);   

			  //Photon vs Jet dist for 0BTag_noMasscut 
			  h_PhPt_vs_bJetPt[5]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			  h_PhEta_vs_bJetEta[5]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
		      
			  //CSVv2 discriminator distributions for 0BTag_noMasscut 
			  h_CSVv2Dist[5]              ->Fill((*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetPt[5]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			  h_CSVv2_vs_bJetEta[5]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);

			  //Primary vertex and number of photon and jets for 0BTag_noMasscut 
			  h_goodPV_TotalWt[5]              ->Fill(GoodVertex, PreBTag_EvtWt);
			  h_nIsoPhotons[5]                 ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			  h_nGoodPhotons[5]                ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of isolated photons with pt > cut and eta < cut 
			  for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			    h_IsoPhotonIdxVsPt[5]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			  }
			  for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			    h_GoodPhotonIdxVsPt[5]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			  }				    
			  h_nJets[5]                       ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			  for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			    h_JetIdxVsPt[5]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			  }
			  //------------------------------------------------------------

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(15.5);
			    h_CutFlowWt_bstar->Fill(15.5, PreBTag_EvtWt);
			    h_CutFlowTotalWt_bstar->Fill(15.5, PreBTag_EvtWt);

			    //----------------------------------------------------------
			    //[6]
			    //Photon distributions 0BTag_Masscut
			    h_PhotonPt[6]               ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			    h_PhotonCalibPt[6]          ->Fill((*phoCalibEt)[PC], PreBTag_EvtWt);
			    h_PhotonEta[6]              ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			    h_PhotonPhi[6]              ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			    h_Photon_SigmaIEtaIEta[6]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			    h_Photon_R9[6]              ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			    h_Photon_HoverE[6]          ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			    h_Photon_EleVeto[6]         ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			    h_Photon_CorrPFChIso[6]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			    h_Photon_CorrPFNeuIso[6]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			    h_Photon_CorrPFPhoIso[6]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
		      
			    //Jet distributions 0BTag_Masscut  
			    h_bJetPt[6]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			    h_bJetEta[6]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			    h_bJetPhi[6]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			    h_bJet_Mt[6]                ->Fill((*jetMt)[JC], PreBTag_EvtWt);
			    h_bJet_area[6]              ->Fill((*jetArea)[JC], PreBTag_EvtWt);
			    h_bJet_Mass[6]              ->Fill((*jetMass)[JC], PreBTag_EvtWt);
			    h_bJet_NHEF[6]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			    h_bJet_NEEF[6]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			    h_bJet_NConst[6]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			    h_bJet_CHEF[6]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			    h_bJet_ChMult[6]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			    h_bJet_CEEF[6]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);
			    h_bJet_MUF[6]               ->Fill((*jetMUF)[JC], PreBTag_EvtWt);
			    h_bJet_NNP[6]               ->Fill((*jetNNP)[JC], PreBTag_EvtWt);
		      
			    //Photon+Jet distributions 0BTag_Masscut 
			    h_GbJetInvtMass_VarBin[6]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GbJetInvtMass_UnitBin[6]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GbJet_dEta[6]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			    h_GbJet_dPhi[6]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			    h_GbJet_dR[6]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			    h_cosThetaStar[6]           ->Fill(GetCosThetaStar((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
		      
			    //PFMet distributions for 0BTag_Masscut 
			    h_PFMet[6]                  ->Fill(pfMET, PreBTag_EvtWt);
			    h_SumPFMet[6]               ->Fill(pfMETsumEt, PreBTag_EvtWt);
			    h_MetBySumMET[6]            ->Fill(pfMET/pfMETsumEt, PreBTag_EvtWt);
			    h_PFMetVsGJmass[6]          ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			    h_PFMetOverSumEtVsGJmass[6] ->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);
			    h_MetByPhPt[6]              ->Fill(pfMET/(*phoEt)[PC], PreBTag_EvtWt);   

			    //Photon vs Jet dist for 0BTag_Masscut 
			    h_PhPt_vs_bJetPt[6]         ->Fill((*phoEt)[PC], (*jetPt)[JC], PreBTag_EvtWt);
			    h_PhEta_vs_bJetEta[6]       ->Fill((*phoSCEta)[PC], (*jetEta)[JC], PreBTag_EvtWt);
		      
			    //CSVv2 discriminator distributions for 0BTag_Masscut 
			    h_CSVv2Dist[6]              ->Fill((*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			    h_CSVv2_vs_bJetPt[6]        ->Fill((*jetPt)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);
			    h_CSVv2_vs_bJetEta[6]       ->Fill((*jetEta)[JC], (*jetCSV2BJetTags)[JC], PreBTag_EvtWt);

			    //Primary vertex and number of photon and jets for 0BTag_Masscut 
			    h_goodPV_TotalWt[6]              ->Fill(GoodVertex, PreBTag_EvtWt);
			    h_nIsoPhotons[6]                 ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);  // Tot # of isolated photons
			    h_nGoodPhotons[6]                ->Fill(GoodIsoBarrelPhotons.size(), PreBTag_EvtWt); // Tot # of isolated photons with pt > cut and eta < cut 
			    for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
			      h_IsoPhotonIdxVsPt[6]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1, PreBTag_EvtWt);
			    }
			    for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
			      h_GoodPhotonIdxVsPt[6]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1, PreBTag_EvtWt);
			    }				    
			    h_nJets[6]                       ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);
			    for(int ij = 0; ij < GoodIsoJets.size(); ij++){
			      h_JetIdxVsPt[6]                ->Fill((*jetPt)[GoodIsoJets[ij]], ij+1, PreBTag_EvtWt);
			    }
			    //-----------------------------------------------------------

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
      /*
      //CSVv2 BTag Efficiencies
      vector<Int_t> goodjets;
      goodjets.clear();
      for(Int_t ii = 0; ii < nJet; ii++){
	if(JetId(ii, "tight")){
	  if(fabs((*jetEta)[ii]) < 2.4){
	    goodjets.push_back(ii);
	  }
	}
      }

      for(Int_t jj = 0; jj < goodjets.size(); jj++){
	Int_t Selectedjet = goodjets[jj];
	if(fabs((*jetHadFlvr)[Selectedjet]) == 5){
	  h_bEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_bEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}
	if(fabs((*jetHadFlvr)[Selectedjet]) == 4){
	  h_cEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_cEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}
	if(fabs((*jetHadFlvr)[Selectedjet]) == 0){
	  h_udsgEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_udsgEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}

	if(CSVv2bTag(Selectedjet, CSV_WP)){

	  Double_t SF_eff, Wt;
          BTagEntry::JetFlavor JF_eff;
	  std::string sys_type_eff = "central"; //central is required to get scale factors (up and down for uncertainties) 
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 5) JF_eff = BTagEntry::FLAV_B; //b
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 4) JF_eff = BTagEntry::FLAV_C; //c
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 0) JF_eff = BTagEntry::FLAV_UDSG; //u,d,s,g,undefined

	  SF_eff = CSVv2bTagSF_auto(CSV_OP, JF_eff, sys_type_eff, (*jetPt)[Selectedjet], (*jetEta)[Selectedjet]);
		    
	  Wt = BTagEventWeight(SF_eff, 1); // SF		     
	  if(JF_eff == BTagEntry::FLAV_B){
	    h_bEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_bEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_bEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_bEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
	  }
	  if(JF_eff == BTagEntry::FLAV_C){
	    h_cEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_cEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_cEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_cEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
  	  }
	  if(JF_eff == BTagEntry::FLAV_UDSG){
	    h_udsgEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_udsgEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_udsgEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_udsgEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
	  }
	}
      }
      //---------------------------
      */
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

			  if(fabs((*jetHadFlvr)[JC]) == 5) JF = BTagEntry::FLAV_B; //b
			  if(fabs((*jetHadFlvr)[JC]) == 4) JF = BTagEntry::FLAV_C; //c
			  if(fabs((*jetHadFlvr)[JC]) == 0) JF = BTagEntry::FLAV_UDSG; //u,d,s,g,undefined

			  SF = CSVv2bTagSF_auto(CSV_OP, JF, sys_type_c, (*jetPt)[JC], (*jetEta)[JC]);
			  SFup = CSVv2bTagSF_auto(CSV_OP, JF, sys_type_u, (*jetPt)[JC], (*jetEta)[JC]);
			  SFdown = CSVv2bTagSF_auto(CSV_OP, JF, sys_type_d, (*jetPt)[JC], (*jetEta)[JC]);

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


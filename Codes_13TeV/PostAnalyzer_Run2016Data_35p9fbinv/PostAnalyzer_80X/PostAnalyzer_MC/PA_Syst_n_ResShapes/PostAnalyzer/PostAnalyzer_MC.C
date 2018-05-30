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
   PileupReWeighting_XSm5();
   PileupReWeighting_XSp5();

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
     Tot_genWt += genWeight;
   }
   cout << "Tot genWt = " << Tot_genWt << endl;
   //---------------------------------------------
   */

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      cout << "<Analyzing entry: " << jentry << endl;
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

      PU_EvtWt_XSm5 = PUWeights_XSm5((*puTrue)[0]);
      Total_EvtWt_XSm5 = Lumi_EvtWt * PU_EvtWt_XSm5;

      PU_EvtWt_XSp5 = PUWeights_XSp5((*puTrue)[0]);
      Total_EvtWt_XSp5 = Lumi_EvtWt * PU_EvtWt_XSp5;

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
	  if( (*phoEt)[GoodIsoPhotons[ip]] > Cut_Photon_pt && fabs((*phoSCEta)[GoodIsoPhotons[ip]]) < Cut_Photon_eta){
	  GoodIsoBarrelPhotons.push_back(GoodIsoPhotons[ip]);
	  }
	}
      }

      if(GoodIsoBarrelPhotons.size() != 0) PC = GoodIsoBarrelPhotons[0];

      GoodIsoJets.clear();
      GoodIsoJets = GoodJets(PC, Cut_JetId); // All jets passing dR > 0.5, tight jet id and pt > 30.0
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

      h_CutFlowWt_qstar_XSm5->Fill(0.5, Total_EvtWt_XSm5);
      h_CutFlowWt_qstar_XSp5->Fill(0.5, Total_EvtWt_XSp5);


      for(int i = 0; i < nJet; i++) cout << (*jetJERSmearing)[i] << endl;

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	h_CutFlow_qstar->Fill(1.5);
	h_CutFlowWt_qstar->Fill(1.5, PreBTag_EvtWt);
	h_CutFlowWt_qstar_XSm5->Fill(1.5, Total_EvtWt_XSm5);
	h_CutFlowWt_qstar_XSp5->Fill(1.5, Total_EvtWt_XSp5);
         
	if(HasPrimaryVtx){
	  h_CutFlow_qstar->Fill(2.5);
	  h_CutFlowWt_qstar->Fill(2.5, PreBTag_EvtWt);
	  h_CutFlowWt_qstar_XSm5->Fill(2.5, Total_EvtWt_XSm5);
	  h_CutFlowWt_qstar_XSp5->Fill(2.5, Total_EvtWt_XSp5);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_qstar->Fill(3.5);
	    h_CutFlowWt_qstar->Fill(3.5, PreBTag_EvtWt);
	    h_CutFlowWt_qstar_XSm5->Fill(3.5, Total_EvtWt_XSm5);
	    h_CutFlowWt_qstar_XSp5->Fill(3.5, Total_EvtWt_XSp5);

	    if(PC > -1){
	      h_CutFlow_qstar->Fill(4.5);
	      h_CutFlowWt_qstar->Fill(4.5, PreBTag_EvtWt);
	      h_CutFlowWt_qstar_XSm5->Fill(4.5, Total_EvtWt_XSm5);
	      h_CutFlowWt_qstar_XSp5->Fill(4.5, Total_EvtWt_XSp5);

	      if(JC > -1){
		h_CutFlow_qstar->Fill(5.5);
		h_CutFlowWt_qstar->Fill(5.5, PreBTag_EvtWt);
		h_CutFlowWt_qstar_XSm5->Fill(5.5, Total_EvtWt_XSm5);
		h_CutFlowWt_qstar_XSp5->Fill(5.5, Total_EvtWt_XSp5);

		//Defining uncert here so that the same functions can be used for both q* and b* loops.
		//JES uncertainty
		PFJet4Vec_JESup.SetPxPyPzE((*jetPx)[JC]*(1+(*jetJECUnc)[JC]), (*jetPy)[JC]*(1+(*jetJECUnc)[JC]), (*jetPz)[JC]*(1+(*jetJECUnc)[JC]), (*jetEn)[JC]*(1+(*jetJECUnc)[JC]));
		PFJet4Vec_JESdown.SetPxPyPzE((*jetPx)[JC]*(1-(*jetJECUnc)[JC]), (*jetPy)[JC]*(1-(*jetJECUnc)[JC]), (*jetPz)[JC]*(1-(*jetJECUnc)[JC]), (*jetEn)[JC]*(1-(*jetJECUnc)[JC]));

		GJ4Vec_JESup.SetPxPyPzE((*phoPx)[PC]+PFJet4Vec_JESup[0], (*phoPy)[PC]+PFJet4Vec_JESup[1], (*phoPz)[PC]+PFJet4Vec_JESup[2], (*phoE)[PC]+PFJet4Vec_JESup[3]);
		GJ4Vec_JESdown.SetPxPyPzE((*phoPx)[PC]+PFJet4Vec_JESdown[0], (*phoPy)[PC]+PFJet4Vec_JESdown[1], (*phoPz)[PC]+PFJet4Vec_JESdown[2], (*phoE)[PC]+PFJet4Vec_JESdown[3]);

		Mass_JESup = pow((GJ4Vec_JESup[3]*GJ4Vec_JESup[3] - GJ4Vec_JESup[0]*GJ4Vec_JESup[0] - GJ4Vec_JESup[1]*GJ4Vec_JESup[1] - GJ4Vec_JESup[2]*GJ4Vec_JESup[2]), 0.5);
		Mass_JESdown = pow((GJ4Vec_JESdown[3]*GJ4Vec_JESdown[3] - GJ4Vec_JESdown[0]*GJ4Vec_JESdown[0] - GJ4Vec_JESdown[1]*GJ4Vec_JESdown[1] - GJ4Vec_JESdown[2]*GJ4Vec_JESdown[2]), 0.5);

		//JER uncertainty
		PFJet4Vec_JER.SetPxPyPzE((*jetPx)[JC]*(*jetJERSmearing)[JC], (*jetPy)[JC]*(*jetJERSmearing)[JC], (*jetPz)[JC]*(*jetJERSmearing)[JC], (*jetEn)[JC]*(*jetJERSmearing)[JC]);

		GJ4Vec_JER.SetPxPyPzE((*phoPx)[PC]+PFJet4Vec_JER[0], (*phoPy)[PC]+PFJet4Vec_JER[1], (*phoPz)[PC]+PFJet4Vec_JER[2], (*phoE)[PC]+PFJet4Vec_JER[3]);

		Mass_JER = pow((GJ4Vec_JER[3]*GJ4Vec_JER[3] - GJ4Vec_JER[0]*GJ4Vec_JER[0] - GJ4Vec_JER[1]*GJ4Vec_JER[1] - GJ4Vec_JER[2]*GJ4Vec_JER[2]), 0.5);

		//PES uncertainty
		Photon4Vec_PESup.SetPxPyPzE((*phoPx)[PC]*(1-(*phoScaleCorrUnc)[PC]), (*phoPy)[PC]*(1-(*phoScaleCorrUnc)[PC]), (*phoPz)[PC]*(1-(*phoScaleCorrUnc)[PC]), (*phoE)[PC]*(1-(*phoScaleCorrUnc)[PC]));
		Photon4Vec_PESdown.SetPxPyPzE((*phoPx)[PC]*(1+(*phoScaleCorrUnc)[PC]), (*phoPy)[PC]*(1+(*phoScaleCorrUnc)[PC]), (*phoPz)[PC]*(1+(*phoScaleCorrUnc)[PC]), (*phoE)[PC]*(1+(*phoScaleCorrUnc)[PC]));

		GJ4Vec_PESup.SetPxPyPzE(Photon4Vec_PESup[0]+(*jetPx)[JC], Photon4Vec_PESup[1]+(*jetPy)[JC], Photon4Vec_PESup[2]+(*jetPz)[JC], Photon4Vec_PESup[3]+(*jetEn)[JC]);
		GJ4Vec_PESdown.SetPxPyPzE(Photon4Vec_PESdown[0]+(*jetPx)[JC], Photon4Vec_PESdown[1]+(*jetPy)[JC], Photon4Vec_PESdown[2]+(*jetPz)[JC], Photon4Vec_PESdown[3]+(*jetEn)[JC]);

		Mass_PESup = pow((GJ4Vec_PESup[3]*GJ4Vec_PESup[3] - GJ4Vec_PESup[0]*GJ4Vec_PESup[0] - GJ4Vec_PESup[1]*GJ4Vec_PESup[1] - GJ4Vec_PESup[2]*GJ4Vec_PESup[2]), 0.5);
		Mass_PESdown = pow((GJ4Vec_PESdown[3]*GJ4Vec_PESdown[3] - GJ4Vec_PESdown[0]*GJ4Vec_PESdown[0] - GJ4Vec_PESdown[1]*GJ4Vec_PESdown[1] - GJ4Vec_PESdown[2]*GJ4Vec_PESdown[2]), 0.5);

		//PER uncertainty
		Photon4Vec_PER.SetPxPyPzE((*phoPx)[PC]*(*phoSmearUnc_nominal)[PC], (*phoPy)[PC]*(*phoSmearUnc_nominal)[PC], (*phoPz)[PC]*(*phoSmearUnc_nominal)[PC], (*phoE)[PC]*(*phoSmearUnc_nominal)[PC]);

		GJ4Vec_PER.SetPxPyPzE(Photon4Vec_PER[0]+(*jetPx)[JC], Photon4Vec_PER[1]+(*jetPy)[JC], Photon4Vec_PER[2]+(*jetPz)[JC], Photon4Vec_PER[3]+(*jetEn)[JC]);

		Mass_PER = pow((GJ4Vec_PER[3]*GJ4Vec_PER[3] - GJ4Vec_PER[0]*GJ4Vec_PER[0] - GJ4Vec_PER[1]*GJ4Vec_PER[1] - GJ4Vec_PER[2]*GJ4Vec_PER[2]), 0.5);

		if(Pass_JetPt){
		  h_CutFlow_qstar->Fill(6.5);
		  h_CutFlowWt_qstar->Fill(6.5, PreBTag_EvtWt);
		  h_CutFlowWt_qstar_XSm5->Fill(6.5, Total_EvtWt_XSm5);
		  h_CutFlowWt_qstar_XSp5->Fill(6.5, Total_EvtWt_XSp5);

		  if(Pass_JetEta){
		    h_CutFlow_qstar->Fill(7.5);
		    h_CutFlowWt_qstar->Fill(7.5, PreBTag_EvtWt);
		    h_CutFlowWt_qstar_XSm5->Fill(7.5, Total_EvtWt_XSm5);
		    h_CutFlowWt_qstar_XSp5->Fill(7.5, Total_EvtWt_XSp5);

 		    if(Pass_GJdPhi){
		      h_CutFlow_qstar->Fill(8.5);
		      h_CutFlowWt_qstar->Fill(8.5, PreBTag_EvtWt);
		      h_CutFlowWt_qstar_XSm5->Fill(8.5, Total_EvtWt_XSm5);
		      h_CutFlowWt_qstar_XSp5->Fill(8.5, Total_EvtWt_XSp5);


		      if(Pass_GJdEta){
			h_CutFlow_qstar->Fill(9.5);
			h_CutFlowWt_qstar->Fill(9.5, PreBTag_EvtWt);
			h_CutFlowWt_qstar_XSm5->Fill(9.5, Total_EvtWt_XSm5);
			h_CutFlowWt_qstar_XSp5->Fill(9.5, Total_EvtWt_XSp5);

			//QCD Overlap removal			
			//uncomment it in script
			//if(${QCD} && IsOverlappedEvent(PC)) continue;

			h_CutFlow_qstar->Fill(10.5);
			h_CutFlowWt_qstar->Fill(10.5, PreBTag_EvtWt);
			h_CutFlowWt_qstar_XSm5->Fill(10.5, Total_EvtWt_XSm5);
			h_CutFlowWt_qstar_XSp5->Fill(10.5, Total_EvtWt_XSp5);

			//JES
			h_profile_JES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_profile_JES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			h_mass_JES[0][0]->Fill(Mass_JESup, PreBTag_EvtWt);
			h_mass_JES[1][0]->Fill(Mass_JESdown, PreBTag_EvtWt);

			h_mass_X_bin1[0]->Fill(GetInvtMass(PC, JC)/1000, PreBTag_EvtWt);

			h_mass_X_bin1_JES[0][0]->Fill(Mass_JESup/1000, PreBTag_EvtWt);
			h_mass_X_bin1_JES[1][0]->Fill(Mass_JESdown/1000, PreBTag_EvtWt);

			//JER
			h_profile_JER[0]->Fill(GetInvtMass(PC, JC), (Mass_JER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_mass_JER[0]->Fill(Mass_JER, PreBTag_EvtWt);
			h_mass_X_bin1_JER[0]->Fill(Mass_JER/1000, PreBTag_EvtWt);

			//PES
			h_profile_PES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_profile_PES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			h_mass_PES[0][0]->Fill(Mass_PESup, PreBTag_EvtWt);
			h_mass_PES[1][0]->Fill(Mass_PESdown, PreBTag_EvtWt);

			h_mass_X_bin1_PES[0][0]->Fill(Mass_PESup/1000, PreBTag_EvtWt);
			h_mass_X_bin1_PES[1][0]->Fill(Mass_PESdown/1000, PreBTag_EvtWt);

			//PER
			h_profile_PER[0]->Fill(GetInvtMass(PC, JC), (Mass_PER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_mass_PER[0]->Fill(Mass_PER, PreBTag_EvtWt);
			h_mass_X_bin1_PER[0]->Fill(Mass_PER/1000, PreBTag_EvtWt);

			//Invt Mass Dist Unit binning                                                    
                        h_GbJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow_qstar->Fill(11.5);
			  h_CutFlowWt_qstar->Fill(11.5, PreBTag_EvtWt);
			  h_CutFlowWt_qstar_XSm5->Fill(11.5, Total_EvtWt_XSm5);
			  h_CutFlowWt_qstar_XSp5->Fill(11.5, Total_EvtWt_XSp5);

			  //Invt Mass Dist Unit binning                                                           
                          h_GbJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

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
		      Double_t SFup, Wt_1Tag_up, Wt_0Tag_up;
		      Double_t SFdown, Wt_1Tag_down, Wt_0Tag_down;
		      BTagEntry::JetFlavor JF;
		      std::string sys_type = "central"; //central is required to get scale factors (up and down for uncertainties) 
		      if(fabs((*jetHadFlvr)[JC]) == 5) JF = BTagEntry::FLAV_B; //b
		      if(fabs((*jetHadFlvr)[JC]) == 4) JF = BTagEntry::FLAV_C; //c
		      if(fabs((*jetHadFlvr)[JC]) == 0) JF = BTagEntry::FLAV_UDSG; //u,d,s,g,undefined

		      SF = CSVv2bTagSF_auto(CSV_OP, JF, sys_type, (*jetPt)[JC], (*jetEta)[JC]);
		      SFup = CSVv2bTagSF_auto(CSV_OP, JF, "up", (*jetPt)[JC], (*jetEta)[JC]);
		      SFdown = CSVv2bTagSF_auto(CSV_OP, JF, "down", (*jetPt)[JC], (*jetEta)[JC]);

		      Wt_1Tag = BTagEventWeight(SF, 1); // SF
		      Wt_0Tag = BTagEventWeight(SF, 0); // (1-SF)  
		      Wt_1Tag_up = BTagEventWeight(SFup, 1); // SF
		      Wt_0Tag_up = BTagEventWeight(SFup, 0); // (1-SF)  
		      Wt_1Tag_down = BTagEventWeight(SFdown, 1); // SF
		      Wt_0Tag_down = BTagEventWeight(SFdown, 0); // (1-SF)  

		      Total_EvtWt_1tag = PreBTag_EvtWt * Wt_1Tag;
		      Total_EvtWt_0tag = PreBTag_EvtWt * Wt_0Tag;
		      Total_EvtWt_1tag_up = PreBTag_EvtWt * Wt_1Tag_up;
		      Total_EvtWt_0tag_up = PreBTag_EvtWt * Wt_0Tag_up;
		      Total_EvtWt_1tag_down = PreBTag_EvtWt * Wt_1Tag_down;
		      Total_EvtWt_0tag_down = PreBTag_EvtWt * Wt_0Tag_down;

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

			  //JES
			  //1BTag
			  h_profile_JES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_profile_JES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			  h_mass_JES[0][1]->Fill(Mass_JESup, Total_EvtWt_1tag);
			  h_mass_JES[1][1]->Fill(Mass_JESdown, Total_EvtWt_1tag);

			  h_mass_X_bin1[1]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_1tag);

			  h_mass_X_bin1_JES[0][1]->Fill(Mass_JESup/1000, Total_EvtWt_1tag);
			  h_mass_X_bin1_JES[1][1]->Fill(Mass_JESdown/1000, Total_EvtWt_1tag);

			  //0BTag
			  h_mass_JES[0][2]->Fill(Mass_JESup, Total_EvtWt_0tag);
			  h_mass_JES[1][2]->Fill(Mass_JESdown, Total_EvtWt_0tag);

			  h_mass_X_bin1[2]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_0tag);

			  h_mass_X_bin1_JES[0][2]->Fill(Mass_JESup/1000, Total_EvtWt_0tag);
			  h_mass_X_bin1_JES[1][2]->Fill(Mass_JESdown/1000, Total_EvtWt_0tag);

			  //JER
			  h_profile_JER[1]->Fill(GetInvtMass(PC, JC), (Mass_JER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_mass_JER[1]->Fill(Mass_JER, Total_EvtWt_1tag);
			  h_mass_JER[2]->Fill(Mass_JER, Total_EvtWt_0tag);

			  h_mass_X_bin1_JER[1]->Fill(Mass_JER/1000, Total_EvtWt_1tag);
			  h_mass_X_bin1_JER[2]->Fill(Mass_JER/1000, Total_EvtWt_0tag);

			  //PES
			  //1BTag
			  h_profile_PES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_profile_PES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			  h_mass_PES[0][1]->Fill(Mass_PESup, Total_EvtWt_1tag);
			  h_mass_PES[1][1]->Fill(Mass_PESdown, Total_EvtWt_1tag);

			  h_mass_X_bin1_PES[0][1]->Fill(Mass_PESup/1000, Total_EvtWt_1tag);
			  h_mass_X_bin1_PES[1][1]->Fill(Mass_PESdown/1000, Total_EvtWt_1tag);

			  //0BTag
			  h_mass_PES[0][2]->Fill(Mass_PESup, Total_EvtWt_0tag);
			  h_mass_PES[1][2]->Fill(Mass_PESdown, Total_EvtWt_0tag);

			  h_mass_X_bin1_PES[0][2]->Fill(Mass_PESup/1000, Total_EvtWt_0tag);
			  h_mass_X_bin1_PES[1][2]->Fill(Mass_PESdown/1000, Total_EvtWt_0tag);

			  //PER
			  h_profile_PER[1]->Fill(GetInvtMass(PC, JC), (Mass_PER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_mass_PER[1]->Fill(Mass_PER, Total_EvtWt_1tag);
			  h_mass_PER[2]->Fill(Mass_PER, Total_EvtWt_0tag);

			  h_mass_X_bin1_PER[1]->Fill(Mass_PER/1000, Total_EvtWt_1tag);
			  h_mass_X_bin1_PER[2]->Fill(Mass_PER/1000, Total_EvtWt_0tag);

			  //BSF
			  h_profile_BSF[0]->Fill(GetInvtMass(PC, JC), (SFup - SF)/SF);
			  h_profile_BSF[1]->Fill(GetInvtMass(PC, JC), (SFdown - SF)/SF);
			 
			  h_mass_X_bin1_BSF[0][0]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_1tag_up);
			  h_mass_X_bin1_BSF[0][1]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_1tag_down);
			  h_mass_X_bin1_BSF[1][0]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_0tag_up);
			  h_mass_X_bin1_BSF[1][1]->Fill(GetInvtMass(PC, JC)/1000, Total_EvtWt_0tag_down);

			  //Invt Mass Dist Unit binning 
			  h_GbJetInvtMass_UnitBin[2]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			  h_GbJetInvtMass_UnitBin[4]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(11.5);
			    h_CutFlowWt_bstar->Fill(11.5, PreBTag_EvtWt);
			    h_CutFlowTotalWt_bstar->Fill(11.5, Total_EvtWt_1tag);
			    h_CutFlowTotalWt_bstar->Fill(15.5, Total_EvtWt_0tag);

			    //Invt Mass Dist Unit binning 
			    h_GbJetInvtMass_UnitBin[3]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_1tag);
			    h_GbJetInvtMass_UnitBin[5]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt_0tag);

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

			  //JES
			  h_profile_JES[0][2]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_profile_JES[1][2]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			  h_mass_JES[0][2]->Fill(Mass_JESup, PreBTag_EvtWt);
			  h_mass_JES[1][2]->Fill(Mass_JESdown, PreBTag_EvtWt);

			  h_mass_X_bin1[2]->Fill(GetInvtMass(PC, JC)/1000, PreBTag_EvtWt);

			  h_mass_X_bin1_JES[0][2]->Fill(Mass_JESup/1000, PreBTag_EvtWt);
			  h_mass_X_bin1_JES[1][2]->Fill(Mass_JESdown/1000, PreBTag_EvtWt);

			  //JER
			  h_profile_JER[2]->Fill(GetInvtMass(PC, JC), (Mass_JER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_mass_JER[2]->Fill(Mass_JER, PreBTag_EvtWt);
			  h_mass_X_bin1_JER[2]->Fill(Mass_JER/1000, PreBTag_EvtWt);

			  //PES
			  h_profile_PES[0][2]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_profile_PES[1][2]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			  h_mass_PES[0][2]->Fill(Mass_PESup, PreBTag_EvtWt);
			  h_mass_PES[1][2]->Fill(Mass_PESdown, PreBTag_EvtWt);

			  h_mass_X_bin1_PES[0][2]->Fill(Mass_PESup/1000, PreBTag_EvtWt);
			  h_mass_X_bin1_PES[1][2]->Fill(Mass_PESdown/1000, PreBTag_EvtWt);

			  //PER
			  h_profile_PER[2]->Fill(GetInvtMass(PC, JC), (Mass_PER - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  h_mass_PER[2]->Fill(Mass_PER, PreBTag_EvtWt);
			  h_mass_X_bin1_PER[2]->Fill(Mass_PER/1000, PreBTag_EvtWt);

			  //BSF                                                            
                          h_mass_X_bin1_BSF[1][0]->Fill(GetInvtMass(PC, JC)/1000, PreBTag_EvtWt);
                          h_mass_X_bin1_BSF[1][1]->Fill(GetInvtMass(PC, JC)/1000, PreBTag_EvtWt);

			  //Invt Mass Dist Unit binning
			  h_GbJetInvtMass_UnitBin[4]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			  if(Pass_GJInvtMass){
			    h_CutFlow_bstar->Fill(15.5);
			    h_CutFlowWt_bstar->Fill(15.5, PreBTag_EvtWt);
			    h_CutFlowTotalWt_bstar->Fill(15.5, PreBTag_EvtWt);

			    //Invt Mass Dist Unit binning
			    h_GbJetInvtMass_UnitBin[5]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

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
 
   }//for jentry

}//Loop()





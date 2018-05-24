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

   Cut_GJInvtMass = 560.0;

   //triggers used in the study (efficiency of these three triggers is same as using all triggers)
   vector<ULong64_t> HLT;
   HLT.clear();
   HLT.push_back(HLT_Photon175_v);
   HLT.push_back(HLT_Photon250_NoHE_v);
   HLT.push_back(HLT_Photon165_HE10_v);

   //triggers for the denominator of trigger turn on
   vector<ULong64_t> HLT_pre;
   HLT_pre.clear();
   HLT_pre.push_back(HLT_Photon75_v);
   HLT_pre.push_back(HLT_Photon90_v);
   HLT_pre.push_back(HLT_Photon120_v);

   //Comment it in script
   //---------------------------
   //Define Output file here
   file = new TFile("PostAnalyzer_Data.root", "RECREATE");
   //---------------------------

   //Uncomment this in script
   /*
   //Define Output file here
   TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   */

   //Define Histograms here
   BookHistograms();

   //   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << "<Total entries" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "Analyzing entry:" << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      GoodVertex = 0;

      //Initially making all bools to false
      Pass_HLT = false;
      HasPrimaryVtx = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);
   
      //Photon Efficiency
      PC_G = -1;
      PC_L = -1;
      PC_M = -1;
      PC_T = -1;
      PC_H = -1;

      JC_G = -1;
      JC_L = -1;
      JC_M = -1;
      JC_T = -1;
      JC_H = -1;

      PC_G = LooserIdforEff();
      PC_L = FirstGoodPhoton("loose");
      PC_M = FirstGoodPhoton("medium");
      PC_T = FirstGoodPhoton("tight");
      PC_H = FirstHighPtIDPhoton();

      JC_G = FirstGoodJet(PC_G);
      JC_L = FirstGoodJet(PC_L);
      JC_M = FirstGoodJet(PC_M);
      JC_T = FirstGoodJet(PC_T);
      JC_H = FirstGoodJet(PC_H);

      h_CutFlow_PhEff->Fill(0.5);

      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC_G > -1){
	    if(fabs((*phoSCEta)[PC_G]) <= 1.4442){//Barrel 
	      h_CutFlow_PhEff->Fill(1.5);

	      h_PassLooserPhId_EB->Fill((*phoEt)[PC_G]);

	      if(JC_G > -1) h_PassLooserInvtMass_EB->Fill(GetInvtMass(PC_G, JC_G));
	      
	      if(PC_L > -1 && PC_G == PC_L){
		h_CutFlow_PhEff->Fill(2.5);

		h_PassPhIdLoose_EB->Fill((*phoEt)[PC_L]);

		if(JC_L > -1) h_PassInvtMassLoose_EB->Fill(GetInvtMass(PC_L, JC_L));
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_CutFlow_PhEff->Fill(3.5);

		h_PassPhIdMedium_EB->Fill((*phoEt)[PC_M]);

		if(JC_M > -1) h_PassInvtMassMedium_EB->Fill(GetInvtMass(PC_M, JC_M));
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_CutFlow_PhEff->Fill(4.5);

		h_PassPhIdTight_EB->Fill((*phoEt)[PC_T]);

		if(JC_T > -1) h_PassInvtMassTight_EB->Fill(GetInvtMass(PC_T, JC_T));
	      }
	    
	      if(PC_H > -1 && PC_G == PC_H){
		h_CutFlow_PhEff->Fill(5.5);

		h_PassPhIdHighPt_EB->Fill((*phoEt)[PC_H]);

		if(JC_H > -1) h_PassInvtMassHighPt_EB->Fill(GetInvtMass(PC_H, JC_H));
	      }
	    }
	    if(fabs((*phoSCEta)[PC_G]) < 2.5 && fabs((*phoSCEta)[PC_G]) >= 1.5666){//Endcap  
	      h_CutFlow_PhEff->Fill(6.5);

	      h_PassLooserPhId_EE->Fill((*phoEt)[PC_G]);

	      if(JC_G > -1) h_PassLooserInvtMass_EE->Fill(GetInvtMass(PC_G, JC_G));

	      if(PC_L > -1 && PC_G == PC_L){
		h_CutFlow_PhEff->Fill(7.5);

		h_PassPhIdLoose_EE->Fill((*phoEt)[PC_L]);

		if(JC_L > -1) h_PassInvtMassLoose_EE->Fill(GetInvtMass(PC_L, JC_L));
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_CutFlow_PhEff->Fill(8.5);

		h_PassPhIdMedium_EE->Fill((*phoEt)[PC_M]);

		if(JC_M > -1) h_PassInvtMassMedium_EE->Fill(GetInvtMass(PC_M, JC_M));
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_CutFlow_PhEff->Fill(9.5);

		h_PassPhIdTight_EE->Fill((*phoEt)[PC_T]);

		if(JC_T > -1) h_PassInvtMassTight_EE->Fill(GetInvtMass(PC_T, JC_T));
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_CutFlow_PhEff->Fill(10.5);

		h_PassPhIdHighPt_EE->Fill((*phoEt)[PC_H]);

		if(JC_H > -1) h_PassInvtMassHighPt_EE->Fill(GetInvtMass(PC_H, JC_H));
	      }
	    }	  
	  }
	}
      }

  


   }//for jentry
}//Loop() 


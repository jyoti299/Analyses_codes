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

   //True if masscut is applied
   IfMassCut = true;

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
   //   cout << "<Total entries" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "Analyzing entry:" << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

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
      Pass_GJInvtMass = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT_Photon165_HE10_v);
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);
      PC = GoodPhoton("medium");
      JC = GoodJet(PC);
      if(PC > -1) Pass_PhoPt = ((*phoEt)[PC] > Cut_Photon_pt);
      if(PC > -1) Pass_PhoEtaEB = (fabs((*phoSCEta)[PC]) <= Cut_Photon_eta);
      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) <= Cut_Jet_eta);
      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             

      //=====================================Standard skelton===========================================
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC > -1){
	    if(Pass_PhoPt){
	      if(Pass_PhoEtaEB){
		if(JC > -1){
		  if(Pass_JetPt){
		    if(Pass_JetEta){
		      if(Pass_GJdPhi){
			if(Pass_GJdEta){
			  if(Pass_GJInvtMass && IfMassCut){
			  }//Pass_GJInvtMass
			}//Pass_GJdEta
		      }//Pass_GJdPhi
		    }//Pass_JetEta
		  }//Pass_JetPt
		}//JC > -1
	      }//Pass_PhoEtaEB
	    }//Pass_PhoPt
	  }//PC > -1
	}//HasPrimaryVtx
      }//Pass_HLT
      //==================================================================================================

      //Extra for advance
      PC_id = -1;
      PC_id = PassPhotonIDNoIso("medium");
      
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC_id > -1){

	    h_PFChIsoVsPt[0]                 ->Fill((*phoEt)[PC_id], (*phoPFChIso)[PC_id]);
	    h_PFChIsoVsRho[0]                ->Fill(rho, (*phoPFChIso)[PC_id]);
	    h_prof_PFChIsoVsPt[0]            ->Fill((*phoEt)[PC_id], (*phoPFChIso)[PC_id]);
	    h_prof_PFChIsoVsRho[0]           ->Fill(rho, (*phoPFChIso)[PC_id]);
	    
	    h_PFNeuIsoVsPt[0][0]             ->Fill((*phoEt)[PC_id], (*phoPFNeuIso)[PC_id]);
	    h_PFNeuIsoVsRho[0][0]            ->Fill(rho, (*phoPFNeuIso)[PC_id]);
	    h_prof_PFNeuIsoVsPt[0][0]        ->Fill((*phoEt)[PC_id], (*phoPFNeuIso)[PC_id]);
	    h_prof_PFNeuIsoVsRho[0][0]       ->Fill(rho, (*phoPFNeuIso)[PC_id]);

	    h_PFPhIsoVsPt[0][0]              ->Fill((*phoEt)[PC_id], (*phoPFPhoIso)[PC_id]);
	    h_PFPhIsoVsRho[0][0]             ->Fill(rho, (*phoPFPhoIso)[PC_id]);
	    h_prof_PFPhIsoVsPt[0][0]         ->Fill((*phoEt)[PC_id], (*phoPFPhoIso)[PC_id]);
	    h_prof_PFPhIsoVsRho[0][0]        ->Fill(rho, (*phoPFPhoIso)[PC_id]);
				
	    h_PFNeuIsoVsPt[0][1]      ->Fill((*phoEt)[PC_id], TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_PFNeuIsoVsRho[0][1]     ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_prof_PFNeuIsoVsPt[0][1] ->Fill((*phoEt)[PC_id], TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_prof_PFNeuIsoVsRho[0][1]->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));

	    h_PFPhIsoVsPt[0][1]       ->Fill((*phoEt)[PC_id], TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_PFPhIsoVsRho[0][1]      ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_prof_PFPhIsoVsPt[0][1]  ->Fill((*phoEt)[PC_id], TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));
	    h_prof_PFPhIsoVsRho[0][1] ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC_id] - rho*EAPhotons((*phoSCEta)[PC_id])), 0.0));

	  }
	}
      }
      
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC > -1){
	    
	    h_PFChIsoVsPt[1]                 ->Fill((*phoEt)[PC], (*phoPFChIso)[PC]);
      	    h_PFChIsoVsRho[1]                ->Fill(rho, (*phoPFChIso)[PC]);
      	    h_prof_PFChIsoVsPt[1]            ->Fill((*phoEt)[PC], (*phoPFChIso)[PC]);
	    h_prof_PFChIsoVsRho[1]           ->Fill(rho, (*phoPFChIso)[PC]);
	    
	    h_PFNeuIsoVsPt[1][0]             ->Fill((*phoEt)[PC], (*phoPFNeuIso)[PC]);
	    h_PFNeuIsoVsRho[1][0]            ->Fill(rho, (*phoPFNeuIso)[PC]);
	    h_prof_PFNeuIsoVsPt[1][0]        ->Fill((*phoEt)[PC], (*phoPFNeuIso)[PC]);
	    h_prof_PFNeuIsoVsRho[1][0]       ->Fill(rho, (*phoPFNeuIso)[PC]);

	    h_PFPhIsoVsPt[1][0]              ->Fill((*phoEt)[PC], (*phoPFPhoIso)[PC]);
	    h_PFPhIsoVsRho[1][0]             ->Fill(rho, (*phoPFPhoIso)[PC]);
	    h_prof_PFPhIsoVsPt[1][0]         ->Fill((*phoEt)[PC], (*phoPFPhoIso)[PC]);
      	    h_prof_PFPhIsoVsRho[1][0]        ->Fill(rho, (*phoPFPhoIso)[PC]);
				
	    h_PFNeuIsoVsPt[1][1]      ->Fill((*phoEt)[PC], TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
	    h_PFNeuIsoVsRho[1][1]     ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
	    h_prof_PFNeuIsoVsPt[1][1] ->Fill((*phoEt)[PC], TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
	    h_prof_PFNeuIsoVsRho[1][1]->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

	    h_PFPhIsoVsPt[1][1]       ->Fill((*phoEt)[PC], TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
	    h_PFPhIsoVsRho[1][1]      ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
            h_prof_PFPhIsoVsPt[1][1]  ->Fill((*phoEt)[PC], TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
     	    h_prof_PFPhIsoVsRho[1][1] ->Fill(rho, TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

          }
        }
      }   

  

   }//for jentry


}//Loop() 





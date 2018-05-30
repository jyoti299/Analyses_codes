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

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.8;

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
   //HLT.push_back(HLT_Photon250_NoHE_v);

   //triggers for the denominator of trigger turn on
   vector<ULong64_t> HLT_deno;
   HLT_deno.clear();
   HLT_deno.push_back(HLT_Mu50_v);

   //Trigger bits for photon matching
   vector<UInt_t> PhoTrigObjs;
   PhoTrigObjs.push_back(hltEG165HE10Filter);
   //PhoTrigObjs.push_back(hltEG250erEtFilter);

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
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      //      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);
      
      h_CutFlow_qstar->Fill(0.5);

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

			if(Pass_GJInvtMass){
			  h_CutFlow_qstar->Fill(10.5);
			  
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

      //---------------------------------------------------------------------------------   
      //Trigger Turn-on
      Bool_t PassHLTNum = false;
      Bool_t PassHLTDeno = false;
      Bool_t PassHLT_Pre = false;
      //      cout << PassHLT_Pre << endl;
      ULong64_t HLT_pre = HLT[0];

      PassHLTNum = PassHLT(HLT);
      PassHLTDeno = PassHLTMu(HLT_deno);
      PassHLT_Pre = PassHLT_Prescale(HLT_pre);

      std::vector<Int_t> TrigPhotons;
      std::vector<Int_t> TrigMuons;

      TrigPhotons.clear();
      for(Int_t i = 0; i < nPho; i++){
	if(CutBasedPhotonID(i, "medium") && fabs((*phoSCEta)[i]) < 2.4 && (*phoEt)[i] > 30.0) TrigPhotons.push_back(i);
      }

      UShort_t muLoose = 0;
      UShort_t muHighPt = 4;

      Int_t MuIdx = -1;
      if(PassHLTDeno){
	for(int k = 0; k < nMu; k++){
          if(MuonId(k, muHighPt) && (*muPt)[k] > 50.0 && (*muEta)[k] < 2.4){
	    MuIdx = k;
	    break;
	  }
	}
      }

      if(MuIdx > -1){
	for(Int_t ii = 0; ii < TrigPhotons.size(); ii++){
	  if(GetdR((*phoSCEta)[TrigPhotons[ii]], (*muEta)[MuIdx], (*phoSCPhi)[TrigPhotons[ii]], (*muPhi)[MuIdx]) > 0.1){
	    h_TrigPhotonPt[0]->Fill((*phoEt)[TrigPhotons[ii]]);
	    if(PassHLTNum){
	      h_TrigPhotonPt[1]->Fill((*phoEt)[TrigPhotons[ii]]);
	      break;
	    }
	  }
	}
      }
      /*
      if(PassHLTDeno){
	for(int k = 0; k < nMu; k++){
	  if(MuonId(k, muHighPt) && (*muPt)[k] > 50.0 && (*muEta)[k] < 2.4){
	    for(Int_t ii = 0; ii < TrigPhotons.size(); ii++){
	      if(GetdR((*phoSCEta)[TrigPhotons[ii]], (*muEta)[k], (*phoSCPhi)[TrigPhotons[ii]], (*muPhi)[k]) > 0.1){
		h_TrigPhotonPt[0]->Fill((*phoEt)[TrigPhotons[ii]]);
		if(PassHLTNum){		 
		  //		  if(PhoTrigObjMatching(ii, PhoTrigObjs)){
		  h_TrigPhotonPt[1]->Fill((*phoEt)[TrigPhotons[ii]]);
		  break;
		    //		  }
		  //            }                         	}
	      }	     
	    }	    
	  }
	  break;
	}
      }
      */
   }//for jentry
}//Loop() 


   


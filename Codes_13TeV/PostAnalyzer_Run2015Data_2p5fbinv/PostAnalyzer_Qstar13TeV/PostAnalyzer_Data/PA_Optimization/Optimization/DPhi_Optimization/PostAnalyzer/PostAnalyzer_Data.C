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

   //List of single photon triggers in data
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
   HLT.push_back(HLT_Photon175_v);
   HLT.push_back(HLT_Photon250_NoHE_v);
   HLT.push_back(HLT_Photon165_HE10_v);

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

      //Defining bools for optimization
      Bool_t Pass_dPhi1p5 = false;
      Bool_t Pass_dPhi2p0 = false;
      Bool_t Pass_dPhi2p5 = false;
      Bool_t Pass_dPhi3p0 = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("loose");

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
      GoodIsoJets = GoodJets(PC);
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) <= Cut_Jet_eta);
      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             

      //Getting bools for optimization
      Pass_dPhi1p5 = ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > 1.5);
      Pass_dPhi2p0 = ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > 2.0);
      Pass_dPhi2p5 = ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > 2.5);
      Pass_dPhi3p0 = ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > 3.0);

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

		    //dphi = 1.5
		    if(Pass_dPhi1p5){
		      h_CutFlow->Fill(8.5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);

			h_GJetInvtMass_bin40_dphi1p5[0]    ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_VarBin_dphi1p5[0]   ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_UnitBin_dphi1p5[0]  ->Fill(GetInvtMass(PC, JC));

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(10.5);

			  h_GJetInvtMass_bin40_dphi1p5[1]    ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_VarBin_dphi1p5[1]   ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_UnitBin_dphi1p5[1]  ->Fill(GetInvtMass(PC, JC));

			}//Pass_GJInvtMass
		      }//Pass_GJdEta
		    }//Pass_DPhi1p5

		    //dphi = 2.0
		    if(Pass_dPhi2p0){
		      h_CutFlow->Fill(11.5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(12.5);

			h_GJetInvtMass_bin40_dphi2p0[0]    ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_VarBin_dphi2p0[0]   ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_UnitBin_dphi2p0[0]  ->Fill(GetInvtMass(PC, JC));

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(13.5);

			  h_GJetInvtMass_bin40_dphi2p0[1]    ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_VarBin_dphi2p0[1]   ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_UnitBin_dphi2p0[1]  ->Fill(GetInvtMass(PC, JC));

			}//Pass_GJInvtMass
		      }//Pass_GJdEta
		    }//Pass_DPhi2p0

		    //dphi = 2.5
		    if(Pass_dPhi2p5){
		      h_CutFlow->Fill(14.5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(15.5);

			h_GJetInvtMass_bin40_dphi2p5[0]    ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_VarBin_dphi2p5[0]   ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_UnitBin_dphi2p5[0]  ->Fill(GetInvtMass(PC, JC));

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(16.5);

			  h_GJetInvtMass_bin40_dphi2p5[1]    ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_VarBin_dphi2p5[1]   ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_UnitBin_dphi2p5[1]  ->Fill(GetInvtMass(PC, JC));

			}//Pass_GJInvtMass
		      }//Pass_GJdEta
		    }//Pass_DPhi2p5

		    //dphi = 3.0
		    if(Pass_dPhi3p0){
		      h_CutFlow->Fill(17.5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(18.5);

			h_GJetInvtMass_bin40_dphi3p0[0]    ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_VarBin_dphi3p0[0]   ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_UnitBin_dphi3p0[0]  ->Fill(GetInvtMass(PC, JC));

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(19.5);

			  h_GJetInvtMass_bin40_dphi3p0[1]    ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_VarBin_dphi3p0[1]   ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_UnitBin_dphi3p0[1]  ->Fill(GetInvtMass(PC, JC));

			}//Pass_GJInvtMass
		      }//Pass_GJdEta
		    }//Pass_DPhi3p0

		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//Pass_PhoPtEta
	  }//PassPhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
   
   

   }//for jentry
}//Loop() 


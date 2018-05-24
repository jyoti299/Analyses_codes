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

   //*********************************************************************************************************//
   //Get Event No., Lumi No. and Run no. for highest invariant mass event
   Double_t htmass = 0.0;
   Int_t RunNo = 0;
   Long64_t EvtNo = 0;
   Int_t LumiNo = 0;
   Double_t htPhoPt = 0.0;
   Double_t htJetPt = 0.0;
   Double_t htPhoEta = 0.0;
   Double_t htJetEta = 0.0;
   Double_t htPhoPhi = 0.0;
   Double_t htJetPhi = 0.0;
   //********************************************************************************************************//

   //   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   //   cout << "<Total entries" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;

   //For Table
   cout << "Event table for events in invariant mass range 1700-2200 GeV/c2" << endl;
   cout << "{\\bf Evt No. } & {\\bf Mass } & {\\bf Photon Pt } & {\\bf Jet Pt } & {\\bf Photon eta } & {\\bf Jet eta } & {\\bf Photon phi } & {\\bf Jet phi } & {\\bf nPhotons (pt > 30) } & {\\bf nPhotons (pt > 190) } & {\\bf nJets } \\\ " << endl;  

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

      //To get nPhotons with leading photon pt > 190 gev and all other > 30 gev.
      Int_t Lead = -1;
      GoodIsoLeadingPhotons.clear();
      if(GoodIsoPhotons.size() != 0){
	for(Int_t ii = 0; ii < GoodIsoPhotons.size(); ii++){
	  if((*phoEt)[GoodIsoPhotons[ii]] > Cut_Photon_pt && fabs((*phoSCEta)[GoodIsoPhotons[ii]]) < Cut_Photon_eta){
	    Lead = ii;
	    break;
	  }
	}
	if(Lead > -1) GoodIsoLeadingPhotons.push_back(GoodIsoPhotons[Lead]);
      	for(Int_t ij = 0; ij < GoodIsoPhotons.size(); ij++){
     	  if(ij == Lead) continue;
	  if( (*phoEt)[GoodIsoPhotons[ij]] > 30 && fabs((*phoSCEta)[GoodIsoPhotons[ij]]) < Cut_Photon_eta){
	    GoodIsoLeadingPhotons.push_back(GoodIsoPhotons[ij]);
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

      if(Pass_HLT){
	if(HasPrimaryVtx && metFilters == 0){
	  if(GoodIsoPhotons.size() > 0){
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
		    if(Pass_GJdPhi){
		      if(Pass_GJdEta){

			//Getting high invariant mass events
			Double_t InvtMass = GetInvtMass(PC, JC);

			if( InvtMass > 1700.0 && InvtMass < 2200.0 ){
			  RunNo = run;
			  EvtNo = event;
			  LumiNo = lumis;			   
			  htPhoPt = (*phoEt)[PC];
			  htJetPt = (*jetPt)[JC];
			  htPhoEta = (*phoSCEta)[PC];
			  htJetEta = (*jetEta)[JC];
			  htPhoPhi = (*phoSCPhi)[PC];
			  htJetPhi = (*jetPhi)[JC];

			  cout << EvtNo << " & "  << InvtMass << "  &  " << htPhoPt << "  &  " << htJetPt << "  &  " << htPhoEta << "  &  " << htJetEta << "  &  " << htPhoPhi << "  &  " << htJetPhi << "  &  " << GoodIsoLeadingPhotons.size() <<  "  &  " << GoodIsoBarrelPhotons.size() << "  &  "  <<  GoodIsoJets.size() << " \\\ " << endl;
 

			}
			    //--------------------------------------
		            
		      }//Pass_GJdEta
		    }//Pass_GJdPhi
		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//Pass_PhoPtEta
	  }//PassPhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
   
   }//for jentry
}//Loop() 


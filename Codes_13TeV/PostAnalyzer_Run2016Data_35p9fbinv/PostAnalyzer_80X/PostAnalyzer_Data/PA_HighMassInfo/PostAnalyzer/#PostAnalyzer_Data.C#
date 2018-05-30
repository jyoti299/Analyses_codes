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

   Cut_Jet_pt = 170.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.5;

   Cut_GJInvtMass = 700.0;

   Cut_PhId = "medium";
   Cut_JetId = "tight";

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
   Double_t htmass    = 3500.0;
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
      //No dphi cut
      Pass_GJdPhi = true;
      Pass_GJdEta = false;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;

      //Running different functions 
      Pass_HLT = PassHLT(HLT);
      GoodVertex = nGoodVtx;
      if(GoodVertex > 0) HasPrimaryVtx = true;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons(Cut_PhId); //All photons passing loose id, residual spikes and pt > 30

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
      GoodIsoJets = GoodJets(Cut_JetId); // All jets passing dR > 0.5, tight jet id and pt > 30.0 
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) < Cut_Jet_eta);
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);
      
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	if(HasPrimaryVtx && metFilters == 1536){
	  if(GoodIsoPhotons.size() > 0){
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
 		    if(Pass_GJdPhi){
		      if(Pass_GJdEta){
			Double_t InvtMass = GetInvtMass(PC, JC);
			if(InvtMass > htmass ){

			  cout << "                       Entry = " << jentry << endl;

			  cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			  cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			  cout <<"                                        Qstar EVENT                                            "<< endl;
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

			  cout <<"                                                                                               "<< endl;
			  cout <<"+++++++++++++++++++++++++++++++++++ EVENT INFO ==> ENTRY = " << jentry << "+++++++++++++++++++++++++++++++"<< endl;
			  cout <<"| InvtMass (GeV/c2) | Ph_Pt (GeV/c) | Jet_Pt (GeV/c) |  Ph_Eta  |  Jet_Eta  |  Ph_Phi  |  Jet_Phi  |  Run_No  |  Event_No  |  Lumi_No  |" << endl;
			  cout << "|   " << InvtMass << "   |   " << htPhoPt << "   |   " << htbJetPt << "   |   " << htPhoEta << "   |   " << htbJetEta << "   |   " << htPhoPhi << "   |   " << htbJetPhi << "   |   " << RunNo << "   |   " << EvtNo << "   |   " << LumiNo << "   |   " << endl;
			  cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			  cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			  cout <<"                                                                                               "<< endl;

			}
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
	if(HasPrimaryVtx && metFilters == 1536){
	  if(GoodIsoPhotons.size() > 0){
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
		    if(Pass_CSVv2bTag){
		      if(Pass_GJdPhi){
			if(Pass_GJdEta){
			  Double_t InvtMass = GetInvtMass(PC, JC);
			  if(InvtMass > htmass ){

			    cout << "                       Entry = " << jentry << endl;

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

			    cout <<"                                                                                               "<< endl;
			    cout <<"+++++++++++++++++++++++++++++++++++ EVENT INFO ==> ENTRY = " << jentry << "+++++++++++++++++++++++++++++++"<< endl;
			    cout <<"| InvtMass (GeV/c2) | Ph_Pt (GeV/c) | Jet_Pt (GeV/c) |  Ph_Eta  |  Jet_Eta  |  Ph_Phi  |  Jet_Phi  |  Run_No  |  Event_No  |  Lumi_No  |" << endl;
			    cout << "|   " << InvtMass << "   |   " << htPhoPt << "   |   " << htbJetPt << "   |   " << htPhoEta << "   |   " << htbJetEta << "   |   " << htPhoPhi << "   |   " << htbJetPhi << "   |   " << RunNo << "   |   " << EvtNo << "   |   " << LumiNo << "   |   " << endl;
			    cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			    cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			    cout <<"                                                                                               "<< endl;

			  }
			}//if(Pass_GJdEta) inside if(Pass_CSVv2bTag)
		      }//if(Pass_GJdPhi) inside if(Pass_CSVv2bTag)		       
		      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    }//if(Pass_CSVv2bTag)		     
		    else{
		      if(Pass_GJdPhi){
			if(Pass_GJdEta){
			  Double_t InvtMass = GetInvtMass(PC, JC);
			  if(InvtMass > htmass ){

			    cout << "                       Entry = " << jentry << endl;

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

			    cout <<"                                                                                               "<< endl;
			    cout <<"+++++++++++++++++++++++++++++++++++ EVENT INFO ==> ENTRY = " << jentry << "+++++++++++++++++++++++++++++++"<< endl;
			    cout <<"| InvtMass (GeV/c2) | Ph_Pt (GeV/c) | Jet_Pt (GeV/c) |  Ph_Eta  |  Jet_Eta  |  Ph_Phi  |  Jet_Phi  |  Run_No  |  Event_No  |  Lumi_No  |" << endl;
			    cout << "|   " << InvtMass << "   |   " << htPhoPt << "   |   " << htbJetPt << "   |   " << htPhoEta << "   |   " << htbJetEta << "   |   " << htPhoPhi << "   |   " << htbJetPhi << "   |   " << RunNo << "   |   " << EvtNo << "   |   " << LumiNo << "   |   " << endl;
			    cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
			    cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
			    cout <<"                                                                                               "<< endl;

			  }
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


   


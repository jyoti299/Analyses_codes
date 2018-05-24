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

   //*********************************************************************************************************//
   //Get Event No., Lumi No. and Run no. for highest invariant mass event
   Double_t htmass = 2500.0;
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
   cout << "<Total entries" << nentries << endl;
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
      Pass_HLT = PassHLT(HLT);
      HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("medium");

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

      h_goodPV[0]     ->Fill(GoodVertex);
      h_nPhotons[0]   ->Fill(GoodIsoPhotons.size());
      h_nJets[0]      ->Fill(GoodIsoJets.size());

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

			h_PC                       ->Fill(PC);
			h_JC                       ->Fill(JC);

			h_PhotonPt[0]              ->Fill((*phoEt)[PC]);
			h_PhotonEta[0]             ->Fill((*phoSCEta)[PC]);
			h_PhotonPhi[0]             ->Fill((*phoSCPhi)[PC]);
			h_Photon_SigmaIEtaIEta[0]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			h_Photon_R9[0]             ->Fill((*phoR9)[PC]);
			h_Photon_HoverE[0]         ->Fill((*phoHoverE)[PC]);
			h_Photon_EleVeto[0]        ->Fill((*phoEleVeto)[PC]);
			h_Photon_CorrPFChIso[0]    ->Fill((*phoPFChIso)[PC]);
			h_Photon_CorrPFNeuIso[0]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			h_Photon_CorrPFPhoIso[0]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			h_JetPt[0]                 ->Fill((*jetPt)[JC]);
			h_JetEta[0]                ->Fill((*jetEta)[JC]);
			h_JetPhi[0]                ->Fill((*jetPhi)[JC]);
			h_Jet_NHEF[0]              ->Fill((*jetNHF)[JC]);
			h_Jet_NEEF[0]              ->Fill((*jetNEF)[JC]);
			h_Jet_NConst[0]            ->Fill((*jetNConstituents)[JC]);
			h_Jet_CHEF[0]              ->Fill((*jetCHF)[JC]);
			h_Jet_ChMult[0]            ->Fill((*jetNCH)[JC]);
			h_Jet_CEEF[0]              ->Fill((*jetCEF)[JC]);

			h_GJetInvtMass_bin40[0]    ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_VarBin[0]   ->Fill(GetInvtMass(PC, JC));
			h_GJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(PC, JC));
			h_GJet_dEta[0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			h_GJet_dPhi[0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			h_GJet_dR[0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));

			h_PFMet[0]                 ->Fill(pfMET);
			h_PFMetVsGJmass[0]         ->Fill(GetInvtMass(PC, JC), pfMET);
                        h_PFMetOverSumEtVsGJmass[0]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			h_goodPV[1]                ->Fill(GoodVertex);
			h_nPhotons[1]              ->Fill(GoodIsoPhotons.size());
			h_nJets[1]                 ->Fill(GoodIsoJets.size());


			if(Pass_GJInvtMass && IfMassCut){
			  h_CutFlow->Fill(10.5);

			  //Getting high invariant mass events
			  Double_t InvtMass = GetInvtMass(PC, JC);

			  if( InvtMass > htmass ){
			    RunNo = run;
			    EvtNo = event;
			    LumiNo = lumis;			    
			    htPhoPt = (*phoEt)[PC];
			    htJetPt = (*jetPt)[JC];
			    htPhoEta = (*phoSCEta)[PC];
			    htJetEta = (*jetEta)[JC];
			    htPhoPhi = (*phoSCPhi)[PC];
			    htJetPhi = (*jetPhi)[JC];

			    cout << "                                                                                " << endl;
			    cout << "********* Event with InvtMass > 2500 GeV, Entry = " << jentry << " ***********" << endl;
			    cout << " | InvtMass (GeV/c2) | Photon Pt (GeV/c) | Jet Pt (GeV/c) | Photon Eta | Jet Eta | Photon Phi | Jet Phi | Run Number | Event Number | Lumi Number | " << endl;
			    cout << " |  " << InvtMass << "  |  " << htPhoPt << "  |  " << htJetPt << "  |  " << htPhoEta << "  |  " << htJetEta << "  |  " << htPhoPhi << "  |  " << htJetPhi << "  |  " << RunNo << "  |  " << EvtNo << "  |  " << LumiNo << "  |  " << endl;
			    cout << "********************************************************************************" << endl;
			    cout << "                                                                                " << endl;

			  }
			    //--------------------------------------

			  h_PhotonPt[1]              ->Fill((*phoEt)[PC]);
			  h_PhotonEta[1]             ->Fill((*phoSCEta)[PC]);
			  h_PhotonPhi[1]             ->Fill((*phoSCPhi)[PC]);
			  h_Photon_SigmaIEtaIEta[1]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
			  h_Photon_R9[1]             ->Fill((*phoR9)[PC]);
			  h_Photon_HoverE[1]         ->Fill((*phoHoverE)[PC]);
			  h_Photon_EleVeto[1]        ->Fill((*phoEleVeto)[PC]);
			  h_Photon_CorrPFChIso[1]    ->Fill((*phoPFChIso)[PC]);
			  h_Photon_CorrPFNeuIso[1]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
			  h_Photon_CorrPFPhoIso[1]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));

			  h_JetPt[1]                 ->Fill((*jetPt)[JC]);
			  h_JetEta[1]                ->Fill((*jetEta)[JC]);
			  h_JetPhi[1]                ->Fill((*jetPhi)[JC]);
			  h_Jet_NHEF[1]              ->Fill((*jetNHF)[JC]);
			  h_Jet_NEEF[1]              ->Fill((*jetNEF)[JC]);
			  h_Jet_NConst[1]            ->Fill((*jetNConstituents)[JC]);
			  h_Jet_CHEF[1]              ->Fill((*jetCHF)[JC]);
			  h_Jet_ChMult[1]            ->Fill((*jetNCH)[JC]);
			  h_Jet_CEEF[1]              ->Fill((*jetCEF)[JC]);

			  h_GJetInvtMass_bin40[1]    ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_VarBin[1]   ->Fill(GetInvtMass(PC, JC));
			  h_GJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC));
			  h_GJet_dEta[1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]));
			  h_GJet_dPhi[1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]));
			  h_GJet_dR[1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]));

			  h_PFMet[1]                 ->Fill(pfMET);
			  h_PFMetVsGJmass[1]         ->Fill(GetInvtMass(PC, JC), pfMET);
			  h_PFMetOverSumEtVsGJmass[1]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt);

			  h_goodPV[2]                ->Fill(GoodVertex);
			  h_nPhotons[2]              ->Fill(GoodIsoPhotons.size());
			  h_nJets[2]                 ->Fill(GoodIsoJets.size());


			  if(pfMET < 500){
			    h_CutFlow->Fill(11.5);

			    if(pfMET < 250){
			      h_CutFlow->Fill(12.5);

			      if(pfMET < 100){
				h_CutFlow->Fill(13.5);

				if(pfMET < 50){
				  h_CutFlow->Fill(14.5);

				}//MET500
			      }//MET250
			    }//MET100
			  }//MET50
			}//Pass_GJInvtMass
		      }//Pass_GJdEta
		    }//Pass_GJdPhi
		  }//Pass_JetEta
		}//Pass_JetPt
	      }//JC > -1
	    }//Pass_PhoPtEta
	  }//PassPhotonID
	}//HasPrimaryVtx
      }//Pass_HLT
   
   
      //Trigger Turn-on
      Bool_t PassHLTNum = false;
      Bool_t PassHLTDeno = false;

      PassHLTNum = PassHLT(HLT);
      PassHLTDeno = PassHLT(HLT_pre);

      if(PC > -1){
	if(fabs((*phoSCEta)[PC]) <= 1.4442){//Barrel
	  if(PassHLTDeno){
	    h_PassHLTDeno_EB->Fill((*phoEt)[PC]);
	    if(PassHLTNum){
	      h_PassHLTNum_EB->Fill((*phoEt)[PC]);
	    }
	  }
	}
	if(fabs((*phoSCEta)[PC]) < 2.5 && fabs((*phoSCEta)[PC]) >= 1.5666){//Endcap 
	  if(PassHLTDeno){
	    h_PassHLTDeno_EE->Fill((*phoEt)[PC]);
            if(PassHLTNum){
              h_PassHLTNum_EE->Fill((*phoEt)[PC]);
            }
          }
	}
      }


      //Checking Trigger Biasing
      vector<ULong64_t> HLT_175;
      vector<ULong64_t> HLT_250;
      HLT_175.clear();
      HLT_250.clear();
      HLT_175.push_back(HLT_Photon175_v);
      HLT_250.push_back(HLT_Photon250_NoHE_v);

      Bool_t Pass_HLT175 = false;
      Bool_t Pass_HLT250 = false;

      Pass_HLT175 = PassHLT(HLT_175);
      Pass_HLT250 = PassHLT(HLT_250);

      h_CutFlow_trig->Fill(0.5);
      if(Pass_HLT250){
	h_CutFlow_trig->Fill(1.5);
	if(Pass_HLT175){
	  h_CutFlow_trig->Fill(2.5);
	}
      }


      //Photon Efficiency
      PC_G = -1;
      PC_L = -1;
      PC_M = -1;
      PC_T = -1;
      PC_H = -1;

      PC_G = LooserIdforEff();
      PC_L = FirstGoodPhoton("loose");
      PC_M = FirstGoodPhoton("medium");
      PC_T = FirstGoodPhoton("tight");
      PC_H = FirstHighPtIDPhoton();

      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC_G > -1){
	    if(fabs((*phoSCEta)[PC_G]) <= 1.4442){//Barrel 

	      h_PassLooserPhId_EB->Fill((*phoEt)[PC_G]);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EB->Fill((*phoEt)[PC_L]);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EB->Fill((*phoEt)[PC_M]);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EB->Fill((*phoEt)[PC_T]);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EB->Fill((*phoEt)[PC_H]);
	      }
	    }
	    if(fabs((*phoSCEta)[PC_G]) < 2.5 && fabs((*phoSCEta)[PC_G]) >= 1.5666){//Endcap  

	      h_PassLooserPhId_EE->Fill((*phoEt)[PC_G]);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EE->Fill((*phoEt)[PC_L]);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EE->Fill((*phoEt)[PC_M]);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EE->Fill((*phoEt)[PC_T]);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EE->Fill((*phoEt)[PC_H]);
	      }
	    }	  
	  }
	}
      }

  


   }//for jentry
}//Loop() 


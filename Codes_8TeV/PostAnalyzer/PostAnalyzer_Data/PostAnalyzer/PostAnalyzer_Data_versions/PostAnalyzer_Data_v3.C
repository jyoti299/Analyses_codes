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
//      Root > .L PostAnalyzer_Data.C
//      Root > PostAnalyzer_Data t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Cut_Photon_pt = 170.0;
   Cut_Photon_eta = 1.4442;
   Cut_Jet_pt = 170.0;
   Cut_Jet_eta = 2.5;
   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 2.0;
   Cut_GJInvtMass = 560.0; 

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
   
   //*******************************************************************************************************//
   //Defining the histogram from filling number of events after various cuts

   const int nbins = 17;
   TString CutFlowLabels[nbins] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag", "PassCSVMBTag", "PassCSVTBTag", "PassPhotonID_NoLICTD"};

   h_CutFlowTable = new TH1F("h_CutFlowTable", "Events Passing Various Cuts", nbins, 0, nbins);
   h_CutFlowTable->GetYaxis()->SetTitle("Events");         h_CutFlowTable->GetYaxis()->CenterTitle();
   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->GetXaxis()->SetBinLabel(i+1, CutFlowLabels[i]);
   }

   Long64_t CutFlowNumber[nbins];
   for(int i = 0; i < nbins; i++){
     CutFlowNumber[i] = 0;
   }
  
   //These variables defined just to check deta and dphi dependence on the b disc efficiency in choosing leading jet as b jet
   Long64_t njets = 0;
   Long64_t ncsvlbjets = 0;
   Long64_t ncsvmbjets = 0;
   Long64_t ncsvtbjets = 0;

   //********************************************************************************************************//

   //Make it true if mass cut to be applied otherwise false
   IfMassCut = false;

   //*********************************************************************************************************//
   //Get Event No., Lumi No. and Run no. for highest invarient mass event
   Double_t htmass = 0.0;
   UInt_t RunNo = 0;
   UInt_t EvtNo = 0;
   UInt_t LumiNo = 0;
   UInt_t BXNo = 0;
   Double_t htPhoPt = 0.0;
   Double_t htBJetPt = 0.0;

   //********************************************************************************************************//

   Long64_t nentries = fChain->GetEntries();
   cout << "no. of entries " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;

   for(Long64_t jentry = 0; jentry < nentries; jentry++){

     //      cout << "++++++++++++++++++Analyzing entry++++++++++++" << jentry << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      PC = -1;
      JC = -1;
      GoodVertex = 0;

      Pass_HLT = false;
      NoScrapingEvt = false;
      HasPrimaryVtx = false;
      Pass_PhoPtcut = false;
      Pass_PhoEtaEBcut = false;
      Pass_JetPtcut = false;
      Pass_JetEtacut = false;
      Pass_GJdPhicut = false;
      Pass_GJdEtacut = false;
      if(IfMassCut == true){
	Pass_GJInvtMasscut = false;
      }else{
	Pass_GJInvtMasscut = true;
      }
      Pass_CSVLBTag = false;
      Pass_CSVMBTag = false;      
      Pass_CSVTBTag = false;

      Pass_HLT = PassHLT();
      NoScrapingEvt = NonScrapingEvt();
      GoodVertex = GoodPrimaryVtx();

      h_goodPV->Fill(GoodVertex);

      if(GoodVertex > 0) HasPrimaryVtx = true;
      PC = GetPhotonPassingAllCuts();
      JC = GetJetPassingIDnMatchedToPhoton(PC);
      if(PC > -1) Pass_PhoPtcut = ((*Photon_pt)[PC] > Cut_Photon_pt);
      if(PC > -1) Pass_PhoEtaEBcut = (fabs((*Photon_SC_eta)[PC]) < Cut_Photon_eta);     
      if(JC > -1) Pass_JetPtcut = ((*PFPatJet_pt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEtacut = (fabs((*PFPatJet_eta)[JC]) < Cut_Jet_eta);
      if(PC > -1 && JC > -1) Pass_GJdPhicut = ((GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEtacut = ((GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC])) < Cut_GJdEta);
      if(IfMassCut == true){
	if(PC > -1 && JC > -1){
	  Pass_GJInvtMasscut = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);
	}
      }
      if(JC > -1) Pass_CSVLBTag = PassCSVLBTag(JC);
      if(JC > -1) Pass_CSVMBTag = PassCSVMBTag(JC);
      if(JC > -1) Pass_CSVTBTag = PassCSVTBTag(JC);

      CutFlowNumber[0]++;

      if(Pass_HLT){
	CutFlowNumber[1]++;

	if(NoScrapingEvt){
	  CutFlowNumber[2]++;

	  if(HasPrimaryVtx){
	    CutFlowNumber[3]++;

	    if(PC > -1){
	      CutFlowNumber[4]++;

	      if(Pass_PhoPtcut){
		CutFlowNumber[5]++;
    
		if(Pass_PhoEtaEBcut){
		  CutFlowNumber[6]++;

		  if(JC > -1){
		    CutFlowNumber[7]++;

		    if(Pass_JetPtcut){
		      CutFlowNumber[8]++;

		      if(Pass_JetEtacut){
			CutFlowNumber[9]++;

			if(Pass_GJdPhicut){
			  CutFlowNumber[10]++;

			  if(Pass_GJdEtacut){
			    CutFlowNumber[11]++;

			    if(Pass_GJInvtMasscut){
			      CutFlowNumber[12]++;
			     
			      h_PhotonPt->Fill((*Photon_pt)[PC]);
			      h_PhotonEta->Fill((*Photon_SC_eta)[PC]);
			      h_PhotonPhi->Fill((*Photon_phi)[PC]);
			      h_PhotonSigmaIEtaIEta->Fill((*Photon_SigmaIEtaIEta)[PC]);
			      h_PhotonSigmaIPhiIPhi->Fill((*Photon_SigmaIPhiIPhi)[PC]);
			      h_Photon_r9->Fill((*Photon_r9)[PC]);
			      h_Photon_SingleTowerHoE->Fill((*Photon_SingleTowerHoE)[PC]);
			      h_Photon_PFIsoCharged03->Fill((*PFIsoCharged03)[PC]);
			      h_Photon_PFIsoNeutral03->Fill((*PFIsoNeutral03)[PC]);
			      h_Photon_PFIsoPhoton03->Fill((*PFIsoPhoton03)[PC]);
			      h_Photon_PFIsoSum03->Fill((*PFIsoSum03)[PC]);

			      h_JetPt->Fill((*PFPatJet_pt)[JC]);
			      h_JetEta->Fill((*PFPatJet_eta)[JC]);
			      h_JetPhi->Fill((*PFPatJet_phi)[JC]);
			      h_Jet_NeutralHadEnergyFrac->Fill((*PFPatJet_NeutralHadEnergyFrac)[JC]);
			      h_Jet_NeutralEmEnergyFrac->Fill((*PFPatJet_NeutralEmEnergyFrac)[JC]);
			      h_Jet_NConstituents->Fill((*PFPatJet_NConstituents)[JC]);
			      h_Jet_ChargedHadEnergyFrac->Fill((*PFPatJet_ChargedHadEnergyFrac)[JC]);
			      h_Jet_ChargedMult->Fill((*PFPatJet_ChargedMult)[JC]);
			      h_Jet_ChargedEmEnergyFrac->Fill((*PFPatJet_ChargedEmEnergyFrac)[JC]);

			      h_GJetInvtMass_binning40->Fill(GetInvtMass(PC, JC));
			      h_GJetInvtMass_VariableBinning->Fill(GetInvtMass(PC, JC));
			      h_GJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			      h_GJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			      h_GJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      h_BJetDiscByCSV->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      h_PC->Fill(PC);
			      h_JC->Fill(JC);

			      if(Pass_CSVLBTag){
			    
				CutFlowNumber[13]++;

				Double_t InvtMass = GetInvtMass(PC, JC);

				if( InvtMass > htmass ){
				  htmass = InvtMass;
				  RunNo = RunNumber;
				  EvtNo = EventNumber;
				  LumiNo = LumiNumber;
				  BXNo = BXNumber;
				  htPhoPt = (*Photon_pt)[PC];
				  htBJetPt = (*PFPatJet_pt)[JC];
				}

				h_CSVL_BJetPt->Fill((*PFPatJet_pt)[JC]);
				h_CSVL_BJetEta->Fill((*PFPatJet_eta)[JC]);
				h_CSVL_BJetPhi->Fill((*PFPatJet_phi)[JC]);

				h_CSVL_GBJetInvtMass_binning40->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetInvtMass_VariableBinning->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVL_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVL_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVL->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      }
			      if(Pass_CSVMBTag){
			    
				CutFlowNumber[14]++;

				h_CSVM_BJetPt->Fill((*PFPatJet_pt)[JC]);
				h_CSVM_BJetEta->Fill((*PFPatJet_eta)[JC]);
				h_CSVM_BJetPhi->Fill((*PFPatJet_phi)[JC]);

				h_CSVM_GBJetInvtMass_binning40->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetInvtMass_VariableBinning->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVM_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVM_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVM->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      }
			      if(Pass_CSVTBTag){
			    
				CutFlowNumber[15]++;

				h_CSVT_BJetPt->Fill((*PFPatJet_pt)[JC]);
				h_CSVT_BJetEta->Fill((*PFPatJet_eta)[JC]);
				h_CSVT_BJetPhi->Fill((*PFPatJet_phi)[JC]);

				h_CSVT_GBJetInvtMass_binning40->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetInvtMass_VariableBinning->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVT_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVT_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVT->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      }//END_OF if(Pass_CSVTBTag)			     		 
			    }//END_OF if(Pass_GJInvtMasscut)
			  }//END_OF if(Pass_GJdEtacut) 
			}//END_OF if(Pass_GJdPhicut) 
		      }//END_OF if(Pass_JetEtacut)
		    }//END_OF if(Pass_JetPtcut)
		  }//END_OF if(JC > -1) 
		}//END_OF if(Pass_PhoEtaEBcut)
	      }//END_OF if(Pass_PhoPtcut)
	    }//END_OF if(PC > -1)
	  }//END_OF if(HasPrimaryVtx)
	}//END_OF if(NoScrapingEvt)
      }//END_OF if(Pass_HLT)
   

      //To Get the efficiency of LICTD Cut
      Int_t PC_NoLICTD = -1;
      Int_t JC_NoLICTD = -1;
      PC_NoLICTD = GetPhotonPassingAllCuts_NoLICTD();
      JC_NoLICTD = GetJetPassingIDnMatchedToPhoton(PC_NoLICTD);
      if(PC_NoLICTD > -1) Pass_PhoPtcut = ((*Photon_pt)[PC_NoLICTD] > Cut_Photon_pt);
      if(PC_NoLICTD > -1) Pass_PhoEtaEBcut = (fabs((*Photon_SC_eta)[PC_NoLICTD]) < Cut_Photon_eta);
      if(JC_NoLICTD > -1) Pass_JetPtcut = ((*PFPatJet_pt)[JC_NoLICTD] > Cut_Jet_pt);
      if(JC_NoLICTD > -1) Pass_JetEtacut = (fabs((*PFPatJet_eta)[JC_NoLICTD]) < Cut_Jet_eta);
      if(PC_NoLICTD > -1 && JC_NoLICTD > -1) Pass_GJdPhicut = ((GetdPhi((*Photon_phi)[PC_NoLICTD], (*PFPatJet_phi)[JC_NoLICTD])) > Cut_GJdPhi);
      if(PC_NoLICTD > -1 && JC_NoLICTD > -1) Pass_GJdEtacut = ((GetdEta((*Photon_SC_eta)[PC_NoLICTD], (*PFPatJet_eta)[JC_NoLICTD])) < Cut_GJdEta);
      if(IfMassCut == true){
	if(PC_NoLICTD > -1 && JC_NoLICTD > -1){
	  Pass_GJInvtMasscut = ((GetInvtMass(PC_NoLICTD, JC_NoLICTD)) > Cut_GJInvtMass);
	}
      }
    
      if(Pass_HLT){
	if(NoScrapingEvt){
	  if(HasPrimaryVtx){
	    if(PC_NoLICTD > -1){
	      if(Pass_PhoPtcut){
		if(Pass_PhoEtaEBcut){
		  if(JC_NoLICTD > -1){
		    if(Pass_JetPtcut){
		      if(Pass_JetEtacut){
			if(Pass_GJdPhicut){
			  if(Pass_GJdEtacut){
			    if(Pass_GJInvtMasscut){
			      CutFlowNumber[16]++;
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      //Getting All_probes and Pass_probes plots for Trigger Efficiency

      Bool_t Pass_HLTPhoton150 = false;
      Bool_t Pass_HLTPhoton75 = false;
      Bool_t Pass_HLTPhoton90 = false;

      Pass_HLTPhoton150 = PassHLT();
      Pass_HLTPhoton75 = PassHLT_Photon75();
      Pass_HLTPhoton90 = PassHLT_Photon90();

      if(Pass_HLTPhoton75 == true || Pass_HLTPhoton90 == true){
	h_PassHLT75OR90_AllProbes->Fill((*Photon_pt)[0]);

	if(Pass_HLTPhoton150 == true){
	  h_PassHLT150_PassProbes->Fill((*Photon_pt)[0]);
	}
      }




      //Getting the number and fraction of Photons, Jets and BJets in each event
      int nPhotons = 0;
      int nJets = 0;
      int nCSVLBJets = 0;
      int nCSVMBJets = 0;
      int nCSVTBJets = 0;
      nPhotons = Photon_n;
      nJets = PFPatJet_n;
      for(int i = 0; i < PFPatJet_n; i++){
	if(PassCSVLBTag(i) == 1){
	  nCSVLBJets++;
	}
	if(PassCSVMBTag(i) == 1){
	  nCSVMBJets++;
	}
	if(PassCSVTBTag(i) == 1){
	  nCSVTBJets++;
	}
      }
      float frac_CSVL = (float)nCSVLBJets/(float)PFPatJet_n;
      float frac_CSVM = (float)nCSVMBJets/(float)PFPatJet_n;
      float frac_CSVT = (float)nCSVTBJets/(float)PFPatJet_n;
      h_nPhotons->Fill(nPhotons);
      h_nJets->Fill(nJets);
      h_nCSVLBJets->Fill(nCSVLBJets);
      h_nCSVMBJets->Fill(nCSVMBJets);
      h_nCSVTBJets->Fill(nCSVTBJets);
      h_CSVL_BJetsFrac->Fill(frac_CSVL);
      h_CSVM_BJetsFrac->Fill(frac_CSVM);
      h_CSVT_BJetsFrac->Fill(frac_CSVT);

      //Filling histogram for PhotonIdx vs PhotonPt
      for(int i = 0; i < Photon_n; i++){
	h_PhotonIdxVsPt->Fill((*Photon_pt)[i], i);
      }
      //Filling histogram for JetIdx vs JetPt
      for(int i = 0; i < PFPatJet_n; i++){
	h_JetIdxVsPt->Fill((*PFPatJet_pt)[i], i);

      }
      //Filling histogram for CSVL-BJetIdx vs CSVL-BJetPt
      for(int i = 0; i < PFPatJet_n; i++){
	if(PassCSVLBTag(i) == 1){
	  h_CSVLBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i);
	}
      }
      //Filling histogram for CSVM-BJetIdx vs CSVM-BJetPt
      for(int i = 0; i < PFPatJet_n; i++){
	if(PassCSVMBTag(i) == 1){
	  h_CSVMBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i);
	}
      }
      //Filling histogram for CSVT-BJetIdx vs CSVT-BJetPt
      for(int i = 0; i < PFPatJet_n; i++){
	if(PassCSVTBTag(i) == 1){
	  h_CSVTBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i);
	}
      }
   
      //Getting no. of leading jets and corresponding B jets without deta and dphi cut
      if(Pass_HLT){
	if(NoScrapingEvt){
	  if(HasPrimaryVtx){
	    if(PC > -1){	    
	      if(Pass_PhoPtcut){	        
		if(Pass_PhoEtaEBcut){	       
		  if(JC > -1){		    
		    if(Pass_JetPtcut){		      
		      if(Pass_JetEtacut){		
			njets++;	
			if(Pass_CSVLBTag){
			  ncsvlbjets++;
			}
			if(Pass_CSVMBTag){
			  ncsvmbjets++;
			}
			if(Pass_CSVTBTag){
			  ncsvtbjets++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      
   }//END_OF for loop 

   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->SetBinContent(i+1, CutFlowNumber[i]);
   }

   //Efficiency of various cuts
   Eff_PassHLT        = (double)CutFlowNumber[1]/(double)CutFlowNumber[0];
   Eff_PassScraping   = (double)CutFlowNumber[2]/(double)CutFlowNumber[1];
   Eff_PassPrimaryVtx = (double)CutFlowNumber[3]/(double)CutFlowNumber[2];
   Eff_PassPhotonID   = (double)CutFlowNumber[4]/(double)CutFlowNumber[3];
   Eff_PassPhotonPt   = (double)CutFlowNumber[5]/(double)CutFlowNumber[4];
   Eff_PassPhotonEta  = (double)CutFlowNumber[6]/(double)CutFlowNumber[5];
   Eff_PassJetID      = (double)CutFlowNumber[7]/(double)CutFlowNumber[6];
   Eff_PassJetPt      = (double)CutFlowNumber[8]/(double)CutFlowNumber[7];
   Eff_PassJetEta     = (double)CutFlowNumber[9]/(double)CutFlowNumber[8];
   Eff_PassDPhi       = (double)CutFlowNumber[10]/(double)CutFlowNumber[9];
   Eff_PassDEta       = (double)CutFlowNumber[11]/(double)CutFlowNumber[10];
   Eff_PassGJInvtMass = (double)CutFlowNumber[12]/(double)CutFlowNumber[11];
   Eff_PassCSVLBTag   = (double)CutFlowNumber[13]/(double)CutFlowNumber[12];
   Eff_PassCSVMBTag   = (double)CutFlowNumber[14]/(double)CutFlowNumber[12];
   Eff_PassCSVTBTag   = (double)CutFlowNumber[15]/(double)CutFlowNumber[12];
   Eff_NodetadphiMasscut_PassCSVL = (double)ncsvlbjets/(double)njets; 
   Eff_NodetadphiMasscut_PassCSVM = (double)ncsvmbjets/(double)njets;
   Eff_NodetadphiMasscut_PassCSVT = (double)ncsvtbjets/(double)njets;
   Eff_LICTD          = (double)CutFlowNumber[12]/(double)CutFlowNumber[16];

   //---------------------------------------------------------------------------------------------------
   cout << "****************************************************************************************************" << endl;
   cout << " Total no. of events = " << CutFlowNumber[0] << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "No. of events passing HLT = " << CutFlowNumber[1] << endl;
   cout << "No. of events passing PrimaryVtx = " << CutFlowNumber[2] << endl;
   cout << "No. of events passing Scraping = " << CutFlowNumber[3] << endl;
   cout << "No. of events passing PhotonID = " << CutFlowNumber[4] << endl;
   cout << "No. of events passing PhotonPt = " << CutFlowNumber[5] << endl;
   cout << "No. of events passing PhotonEta = " << CutFlowNumber[6] << endl;
   cout << "No. of events passing JetID = " << CutFlowNumber[7] << endl;
   cout << "No. of events passing JetPt = " << CutFlowNumber[8] << endl;
   cout << "No. of events passing JetEta = " << CutFlowNumber[9] << endl;
   cout << "No. of events passing DPhiCut = " << CutFlowNumber[10] << endl;
   cout << "No. of events passing DetaCut = " << CutFlowNumber[11] << endl;
   cout << "No. of events passing InvtMass Cut = " << CutFlowNumber[12] << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "No. of events passing CSVL BTag = " << CutFlowNumber[13] << endl;
   cout << "No. of events passing CSVM BTag = " << CutFlowNumber[14] << endl;
   cout << "No. of events passing CSVT BTag = " << CutFlowNumber[15] << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "No. of events having leading jet without deta, dphi and mass cut = " << njets << endl;
   cout << "No. of events passing CSVL BTag without deta, dphi and mass cut = " << ncsvlbjets << endl;
   cout << "No. of events passing CSVM BTag without deta, dphi and mass cut = " << ncsvmbjets << endl;
   cout << "No. of events passing CSVT BTag without deta, dphi and mass cut = " << ncsvtbjets << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing with LICTD cut = " << CutFlowNumber[12] << endl;
   cout << "No. of events passing without LICTD cut = " << CutFlowNumber[16] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "*********************************Efficiency of various cuts******************************************" << endl;
   cout << "Eff_PassHLT                = " << Eff_PassHLT*100 << "%" << endl;
   cout << "Eff_PassScraping           = " << Eff_PassScraping*100 << "%" << endl;
   cout << "Eff_PassPrimaryVtx         = " << Eff_PassPrimaryVtx*100 << "%" << endl;
   cout << "Eff_PassPhotonID           = " << Eff_PassPhotonID*100 << "%" << endl;
   cout << "Eff_PassPhotonPt           = " << Eff_PassPhotonPt*100 << "%" << endl;
   cout << "Eff_PassPhotonEta          = " << Eff_PassPhotonEta*100 << "%" << endl;
   cout << "Eff_PassJetID              = " << Eff_PassJetID*100 << "%" << endl;
   cout << "Eff_PassJetPt              = " << Eff_PassJetPt*100 << "%" << endl;
   cout << "Eff_PassJetEta             = " << Eff_PassJetEta*100 << "%" << endl;
   cout << "Eff_PassDPhi               = " << Eff_PassDPhi*100 << "%" << endl;
   cout << "Eff_PassDEta               = " << Eff_PassDEta*100 << "%" << endl;
   cout << "Eff_PassGJInvtMass         = " << Eff_PassGJInvtMass*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag           = " << Eff_PassCSVLBTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag           = " << Eff_PassCSVMBTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag           = " << Eff_PassCSVTBTag*100 << "%" << endl;
   cout << "Eff_NodetadphiMasscut_PassCSVL = " << Eff_NodetadphiMasscut_PassCSVL*100 << "%" << endl;
   cout << "Eff_NodetadphMasscut_PassCSVM = " << Eff_NodetadphiMasscut_PassCSVM*100 << "%" << endl;
   cout << "Eff_NodetadphiMasscut_PassCSVT = " << Eff_NodetadphiMasscut_PassCSVT*100 << "%" << endl;
   cout << "Eff_LICTD                  = " << Eff_LICTD*100 << "%" << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "********************************HIGHEST MASS EVENT INFORMATION**************************************" << endl;
   cout << "Highest InvtMass of Photon and CSVL-bJet = " << htmass << " GeV/c2" << endl;
   cout << "Highest mass event's Photon Pt = " << htPhoPt << " GeV/c" << endl;
   cout << "Highets mass event's CSVL-bJet Pt = " << htBJetPt << " GeV/c" << endl;
   cout << "Highest mass evnent's Run Number = " << RunNo << endl;
   cout << "Highest mass evnent's Event Number = " << EvtNo << endl;
   cout << "Highest mass evnent's Lumi Number = " << LumiNo << endl;
   cout << "Highest mass evnent's BX Number = " << BXNo << endl;
   cout << "****************************************************************************************************" << endl;
 

 
}


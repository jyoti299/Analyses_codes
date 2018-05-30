#define PostAnalyzer_Data_cxx
#include "PostAnalyzer_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

int main(){
  PostAnalyzer_Data t;
  t.Loop();
  return 0;
}

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
   Cut_Jet_eta = 3.0;
   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 2.0;
   Cut_GJInvtMass = 560.0; //defined but not used yet

   //Define Output file here
   file = new TFile("PostAnalyzer_Data.root", "RECREATE");

   //Define Histograms here
   BookHistograms();
   
   //*******************************************************************************************************//
   //Defining the histogram from filling number of events after various cuts

   const int nbins = 15;
   TString CutFlowLabels[nbins] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassCSVLBTag", "PassCSVMBTag", "PassCSVTBTag"};

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

   Long64_t nentries = fChain->GetEntries();
   cout << "no. of entries " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;

   for(Long64_t jentry = 0; jentry <nentries; jentry++){

     //     cout << "++++++++++++++++++Analyzing entry++++++++++++" << jentry << endl;

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
      //      Pass_GJInvtMasscut = false;
      Pass_CSVLBTag = false;
      Pass_CSVMBTag = false;      
      Pass_CSVTBTag = false;

      Pass_HLT = PassHLT();
      NoScrapingEvt = NonScrapingEvt();
      GoodVertex = GoodPrimaryVtx();
      if(GoodVertex > 0) HasPrimaryVtx = true;
      PC = GetPhotonPassingAllCuts();
      JC = GetJetPassingIDnMatchedToPhoton(PC);
      if(PC > -1) Pass_PhoPtcut = ((*Photon_pt)[PC] > Cut_Photon_pt);
      if(PC > -1) Pass_PhoEtaEBcut = (fabs((*Photon_SC_eta)[PC]) < Cut_Photon_eta);     
      if(JC > -1) Pass_JetPtcut = ((*PFPatJet_pt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEtacut = (fabs((*PFPatJet_eta)[JC]) < Cut_Jet_eta);
      if(PC > -1 && JC > -1) Pass_GJdPhicut = ((GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEtacut = ((GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC])) < Cut_GJdEta);
      //      if(PC > -1 && JC > -1) Pass_GJInvtMasscut = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);
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

			    h_PhotonPt->Fill((*Photon_pt)[PC]);
			    h_PhotonEta->Fill((*Photon_SC_eta)[PC]);
			    h_PhotonPhi->Fill((*Photon_phi)[PC]);

			    h_JetPt->Fill((*PFPatJet_pt)[JC]);
			    h_JetEta->Fill((*PFPatJet_eta)[JC]);
			    h_JetPhi->Fill((*PFPatJet_phi)[JC]);

			    h_GJetInvtMass->Fill(GetInvtMass(PC, JC));
			    h_GJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			    h_GJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			    h_GJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			    h_BJetDiscByCSV->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			    h_PC->Fill(PC);
			    h_JC->Fill(JC);

			    if(Pass_CSVLBTag){
			    
			      CutFlowNumber[12]++;

			      h_CSVL_BJetPt->Fill((*PFPatJet_pt)[JC]);
			      h_CSVL_BJetEta->Fill((*PFPatJet_eta)[JC]);
			      h_CSVL_BJetPhi->Fill((*PFPatJet_phi)[JC]);

			      h_CSVL_GBJetInvtMass->Fill(GetInvtMass(PC, JC));
			      h_CSVL_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			      h_CSVL_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			      h_CSVL_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			    }
			    if(Pass_CSVMBTag){
			    
			      CutFlowNumber[13]++;

			      h_CSVM_BJetPt->Fill((*PFPatJet_pt)[JC]);
			      h_CSVM_BJetEta->Fill((*PFPatJet_eta)[JC]);
			      h_CSVM_BJetPhi->Fill((*PFPatJet_phi)[JC]);

			      h_CSVM_GBJetInvtMass->Fill(GetInvtMass(PC, JC));
			      h_CSVM_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			      h_CSVM_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			      h_CSVM_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			    }
			    if(Pass_CSVTBTag){
			    
			      CutFlowNumber[14]++;

			      h_CSVT_BJetPt->Fill((*PFPatJet_pt)[JC]);
			      h_CSVT_BJetEta->Fill((*PFPatJet_eta)[JC]);
			      h_CSVT_BJetPhi->Fill((*PFPatJet_phi)[JC]);

			      h_CSVT_GBJetInvtMass->Fill(GetInvtMass(PC, JC));
			      h_CSVT_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			      h_CSVT_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			      h_CSVT_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			    }//END_OF if(Pass_CSVTBTag)			     		 
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
     

      //Getting the fraction of BJets in each event
      int nCSVLBJets = 0;
      int nCSVMBJets = 0;
      int nCSVTBJets = 0;
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
      h_CSVL_BJetsFrac->Fill(frac_CSVL);
      h_CSVM_BJetsFrac->Fill(frac_CSVM);
      h_CSVT_BJetsFrac->Fill(frac_CSVT);

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

   Eff_CSVL = 100.0*(double)CutFlowNumber[12]/(double)CutFlowNumber[11];
   Eff_CSVM = 100.0*(double)CutFlowNumber[13]/(double)CutFlowNumber[11];
   Eff_CSVT  =100.0*(double)CutFlowNumber[14]/(double)CutFlowNumber[11];
   Eff_Nodetadphicut_CSVL = 100.0*(double)ncsvlbjets/(double)njets;
   Eff_Nodetadphicut_CSVM = 100.0*(double)ncsvmbjets/(double)njets;
   Eff_Nodetadphicut_CSVT = 100.0*(double)ncsvtbjets/(double)njets;


   cout << " Total no. of events = " << CutFlowNumber[0] << endl;
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
   cout << "No. of events passing CSVL BTag = " << CutFlowNumber[12] << endl;
   cout << "No. of events passing CSVM BTag = " << CutFlowNumber[13] << endl;
   cout << "No. of events passing CSVT BTag = " << CutFlowNumber[14] << endl;
   cout << "-----------------------------------------------------------" << endl;
   cout << "No. of events having leading jet without deta and dphi cut = " << njets << endl;
   cout << "No. of events passing CSVL BTag without deta and dphi cut = " << ncsvlbjets << endl;
   cout << "No. of events passing CSVM BTag without deta and dphi cut = " << ncsvmbjets << endl;
   cout << "No. of events passing CSVT BTag without deta and dphi cut = " << ncsvtbjets << endl;
   cout << "------------------------------------------------------------" << endl;
   cout << "Eff_CSVL = " << Eff_CSVL << endl;
   cout << "Eff_CSVM = " << Eff_CSVM << endl;
   cout << "Eff_CSVT = " << Eff_CSVT << endl;
   cout << "Eff_Nodetadphicut_CSVL = " << Eff_Nodetadphicut_CSVL << endl;
   cout << "Eff_Nodetadphicut_CSVM = " << Eff_Nodetadphicut_CSVM << endl;
   cout << "Eff_Nodetadphicut_CSVT = " << Eff_Nodetadphicut_CSVT << endl;
}

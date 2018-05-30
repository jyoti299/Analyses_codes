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
   Cut_GJInvtMass = 560.0;

   //Define Output file here
   file = new TFile("PostAnalyzer_Data.root", "RECREATE");

   //Define Histograms here
   BookHistograms();
   
   //*******************************************************************************************************//
   //Defining the histogram from filling number of events after various cuts

   const int nbins = 14;
   TString CutFlowLabels[nbins] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassInvtMass", "PassCSVMBTag"};

   h_CutFlowTable = new TH1F("h_CutFlowTable", "Events Passing Various Cuts", nbins, 0, nbins);
   h_CutFlowTable->GetYaxis()->SetTitle("Events");         h_CutFlowTable->GetYaxis()->CenterTitle();
   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->GetXaxis()->SetBinLabel(i+1, CutFlowLabels[i]);
   }

   Long64_t CutFlowNumber[nbins];
   for(int i = 0; i < nbins; i++){
     CutFlowNumber[i] = 0;
   }
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
      Pass_GJInvtMasscut = false;
      Pass_CSVMBTag = false;      
      
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
      if(PC > -1 && JC > -1) Pass_GJInvtMasscut = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);
      if(JC > -1) Pass_CSVMBTag = PassCSVMBTag(JC);

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

			      h_JetPt->Fill((*PFPatJet_pt)[JC]);
			      h_JetEta->Fill((*PFPatJet_eta)[JC]);
			      h_JetPhi->Fill((*PFPatJet_phi)[JC]);

			      h_GJetInvtMass->Fill(GetInvtMass(PC, JC));
			      h_GJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
			      h_GJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
			      h_GJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      h_PC->Fill(PC);
			      h_JC->Fill(JC);
			      
			      if(Pass_CSVMBTag){
				CutFlowNumber[13]++;

				h_BJetPt->Fill((*PFPatJet_pt)[JC]);
				h_BJetEta->Fill((*PFPatJet_eta)[JC]);
				h_BJetPhi->Fill((*PFPatJet_phi)[JC]);

				h_GBJetInvtMass->Fill(GetInvtMass(PC, JC));
				h_GBJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_GBJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_GBJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      } //END_OF if(Pass_CSVMBTag)
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
     

      //Getting the fraction of BJets in each event
      int nBJets = 0;
      for(int i = 0; i < PFPatJet_n; i++){
	if(PassCSVMBTag(i) == 1){
	  nBJets++;
	}
      }
      float frac = (float)nBJets/(float)PFPatJet_n;
      h_BJetsFrac->Fill(frac);
      
   }//END_OF for loop 

   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->SetBinContent(i+1, CutFlowNumber[i]);
   }


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
   cout << "No. of events passing MassCut = " << CutFlowNumber[12] << endl;
   cout << "No. of events passing BTag = " << CutFlowNumber[13] << endl;
}

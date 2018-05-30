#!/bin/bash

for i in Data_Run2012A_22Jan_AOD_NTuples Data_Run2012B_22Jan_AOD_NTuples Data_Run2012C_22Jan_AOD_NTuples 
#for i in Data_Run2012B_22Jan_AOD_NTuples 
do 

#This one needs to change to save results in a new directory
PostAna_this=PostAna_Results_MassCut532_jetEta2p4_dphi1p5_deta2p0_PhIdLoose

eosDir=/eos/uscms/store/user/rocky86
lpcqstarDir=/eos/uscms/store/user/lpcqstar

sourceDir_rocky=/eos/uscms/store/user/rocky86/NTuples/Data/Data22Jan2013ReReco/${i}/
sourceDir_lpcqstar=/eos/uscms/store/user/lpcqstar/NTuples/Data/Data22Jan2013ReReco/${i}/

#if [ ${i} = "Data_Run2012D_22Jan_AOD_NTuples" ]
#then
#sourceDir=${sourceDir_lpcqstar}
#else
sourceDir=${sourceDir_rocky}
#fi

if [ -f PostAnalyzer_Data.C ]; then
echo "++++++++++++++ Deleting PostAnalyzer_Data.C ++++++++++++++"
rm PostAnalyzer_Data.C
fi

if [ -f PostAnalyzer_Data.h ]; then
echo "++++++++++++++ Deleting PostAnalyzer_Data.h ++++++++++++++"
rm PostAnalyzer_Data.h
fi

filenameTag=${i}

myDir=${lpcqstarDir}/Rocky
PostAna_Dir=${myDir}/PostAnalyzer_Results
PostAna_Dir_final=${PostAna_Dir}/${PostAna_this}
outputDir=${PostAna_Dir_final}/Data
destinationDir=${outputDir}
#destinationDir=${outputDir}/Check_Data #used this just to check save files in another directory, dont need it for normal job submission

if [ ! -d ${myDir} ]; then
echo "Making Directory ${myDir}"
mkdir ${myDir}
chmod 775 ${myDir}
fi

if [ ! -d ${PostAna_Dir} ]; then
echo "Making Directory ${PostAna_Dir}"
mkdir ${PostAna_Dir}
chmod 775 ${PostAna_Dir}
fi

if [ ! -d ${PostAna_Dir_final} ]; then
echo "Making Directory ${PostAna_Dir_final}"
mkdir ${PostAna_Dir_final}
chmod 775 ${PostAna_Dir_final}
fi

if [ ! -d ${destinationDir} ]; then
echo "Making Directory ${destinationDir}"
mkdir ${destinationDir}
chmod 775 ${destinationDir}
fi

echo "filenameTag = ${filenameTag}"
echo "destinationDir = ${destinationDir}"
echo "sourceDir = ${sourceDir}"


########## Making PostAnalyzer_Data.C ############
cat > PostAnalyzer_Data.C <<EOF
#define PostAnalyzer_Data_cxx
#include "PostAnalyzer_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
   Cut_Jet_eta = 2.4;
   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 2.0;
   Cut_GJInvtMass = 532.0; 

   //Comment it in script
   //---------------------------
   //Define Output file here
   //file = new TFile("PostAnalyzer_Data.root", "RECREATE");
   //---------------------------

   //Uncomment this in script
   
   //Define Output file here
   TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   

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
  
   //********************************************************************************************************//

   //Make it true if mass cut to be applied otherwise false
   IfMassCut = true;

   //Make true for which Photon Id to apply
   IfLoosePhId = true;
   IfMediumPhId = false;
   IfTightPhId = false;

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
      if(IfLoosePhId == true){
	PC = GetPhotonPassingAllCuts_WithLooseID();
      }
      if(IfMediumPhId == true){
	PC = GetPhotonPassingAllCuts_WithMediumID();
      }
      if(IfTightPhId == true){
	PC = GetPhotonPassingAllCuts_WithTightID();
      }      
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
			      h_GJetInvtMass_UnitBinning->Fill(GetInvtMass(PC, JC));
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

				h_CSVL_BJetPt_1bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVL_BJetEta_1bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVL_BJetPhi_1bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVL_GBJetInvtMass_binning40_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetInvtMass_VariableBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetInvtMass_UnitBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetdEta_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVL_GBJetdPhi_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVL_GBJetdR_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVL->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);
			      
			      }else{

				h_CSVL_BJetPt_0bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVL_BJetEta_0bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVL_BJetPhi_0bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVL_GBJetInvtMass_binning40_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetInvtMass_VariableBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetInvtMass_UnitBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVL_GBJetdEta_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVL_GBJetdPhi_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVL_GBJetdR_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      }
			      if(Pass_CSVMBTag){
			    
				CutFlowNumber[14]++;

				h_CSVM_BJetPt_1bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVM_BJetEta_1bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVM_BJetPhi_1bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVM_GBJetInvtMass_binning40_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetInvtMass_VariableBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetInvtMass_UnitBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetdEta_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVM_GBJetdPhi_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVM_GBJetdR_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVM->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      }else{

				h_CSVM_BJetPt_0bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVM_BJetEta_0bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVM_BJetPhi_0bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVM_GBJetInvtMass_binning40_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetInvtMass_VariableBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetInvtMass_UnitBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVM_GBJetdEta_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVM_GBJetdPhi_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVM_GBJetdR_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      }
			      if(Pass_CSVTBTag){
			    
				CutFlowNumber[15]++;

				h_CSVT_BJetPt_1bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVT_BJetEta_1bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVT_BJetPhi_1bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVT_GBJetInvtMass_binning40_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetInvtMass_VariableBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetInvtMass_UnitBinning_1bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetdEta_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVT_GBJetdPhi_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVT_GBJetdR_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

				h_BJetDiscByCSV_PassingCSVT->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC]);

			      }else{

				h_CSVT_BJetPt_0bTag->Fill((*PFPatJet_pt)[JC]);
				h_CSVT_BJetEta_0bTag->Fill((*PFPatJet_eta)[JC]);
				h_CSVT_BJetPhi_0bTag->Fill((*PFPatJet_phi)[JC]);

				h_CSVT_GBJetInvtMass_binning40_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetInvtMass_VariableBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetInvtMass_UnitBinning_0bTag->Fill(GetInvtMass(PC, JC));
				h_CSVT_GBJetdEta_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]));
				h_CSVT_GBJetdPhi_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]));
				h_CSVT_GBJetdR_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]));

			      }//END_OF else(Pass_CSVTBTag)			     		 
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
      if(IfLoosePhId == true){
	PC_NoLICTD = GetPhotonPassingAllCuts_WithLooseID_NoLICTD();
      }
      if(IfMediumPhId == true){
	PC_NoLICTD = GetPhotonPassingAllCuts_WithMediumID_NoLICTD();
      }
      if(IfTightPhId == true){
	PC_NoLICTD = GetPhotonPassingAllCuts_WithTightID_NoLICTD();
      }
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

      if(PC > -1){
	if(fabs((*Photon_SC_eta)[PC]) < 1.4442){
	  if(Pass_HLTPhoton75 == true || Pass_HLTPhoton90 == true){
	    h_PassHLT75OR90_EBAllProbes->Fill((*Photon_pt)[PC]);

	    if(Pass_HLTPhoton150 == true){
	      h_PassHLT150_EBPassProbes->Fill((*Photon_pt)[PC]);
	    }
	  }
	}

	if(fabs((*Photon_SC_eta)[PC]) < 2.5 && fabs((*Photon_SC_eta)[PC]) > 1.5666){
	  if(Pass_HLTPhoton75 == true || Pass_HLTPhoton90 == true){
	    h_PassHLT75OR90_EEAllProbes->Fill((*Photon_pt)[PC]);

	    if(Pass_HLTPhoton150 == true){
	      h_PassHLT150_EEPassProbes->Fill((*Photon_pt)[PC]);
	    }
	  }
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





EOF



########## Making PostAnalyzer_Data.h ############
cat > PostAnalyzer_Data.h <<EOF
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 28 07:27:12 2014 by ROOT version 5.32/00
// from TChain myEvent/
//////////////////////////////////////////////////////////

#ifndef PostAnalyzer_Data_h
#define PostAnalyzer_Data_h

//ROOT include files
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

//c++ include files
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <map>

using namespace std;
using namespace ROOT;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxIsScrapingEvent = 1;

class PostAnalyzer_Data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //Define Output File
   TFile *file;

   //If need to apply mass cut or not
   Bool_t IfMassCut;

   //For Loose, Medium or Tight Photon Id
   Bool_t IfLoosePhId;
   Bool_t IfMediumPhId;
   Bool_t IfTightPhId ;

   //Self Declared Variables
   Bool_t Pass_HLT;
   Bool_t NoScrapingEvt;
   Bool_t HasPrimaryVtx;
   Bool_t Pass_PhoPtcut;
   Bool_t Pass_PhoEtaEBcut;
   Bool_t Pass_JetPtcut;
   Bool_t Pass_JetEtacut;
   Bool_t Pass_GJdPhicut;
   Bool_t Pass_GJdEtacut;
   Bool_t Pass_GJInvtMasscut;
   Bool_t Pass_CSVLBTag;
   Bool_t Pass_CSVMBTag;
   Bool_t Pass_CSVTBTag;

   Int_t GoodVertex;
   Int_t PC, JC;
  
   Double_t Cut_Vtx_z; //(cm)
   Double_t Cut_Vtx_ndof;
   Double_t Cut_Vtx_rho; //(cm)

   Double_t Cut_Photon_pt; //(GeV)
   Double_t Cut_Photon_eta; 
   Double_t Cut_Jet_pt; //(GeV)
   Double_t Cut_Jet_eta; 
   Double_t Cut_GJdPhi;
   Double_t Cut_GJdEta;
   Double_t Cut_GJInvtMass;

   Double_t Eff_PassHLT, Eff_PassScraping, Eff_PassPrimaryVtx, Eff_PassPhotonID, Eff_PassPhotonPt, Eff_PassPhotonEta;
   Double_t Eff_PassJetID, Eff_PassJetPt, Eff_PassJetEta, Eff_PassDPhi, Eff_PassDEta, Eff_PassGJInvtMass;
   Double_t Eff_PassCSVLBTag, Eff_PassCSVMBTag, Eff_PassCSVTBTag;
   Double_t Eff_LICTD;

   //Declaration of Histograms
   //Historgrams for Photon variables
   TH1F *h_PhotonPt;
   TH1F *h_PhotonEta;
   TH1F *h_PhotonPhi;
   TH1F *h_PhotonSigmaIEtaIEta;
   TH1F *h_PhotonSigmaIPhiIPhi;
   TH1F *h_Photon_r9;
   TH1F *h_Photon_SingleTowerHoE;
   TH1F *h_Photon_PFIsoCharged03;
   TH1F *h_Photon_PFIsoNeutral03;
   TH1F *h_Photon_PFIsoPhoton03;
   TH1F *h_Photon_PFIsoSum03;

   //Histograms for Jet variables
   TH1F *h_JetPt;
   TH1F *h_JetEta;
   TH1F *h_JetPhi;
   TH1F *h_Jet_NeutralHadEnergyFrac;
   TH1F *h_Jet_NeutralEmEnergyFrac;
   TH1F *h_Jet_NConstituents;
   TH1F *h_Jet_ChargedHadEnergyFrac;
   TH1F *h_Jet_ChargedMult;
   TH1F *h_Jet_ChargedEmEnergyFrac;

   //Histograms for Photon-Jet variables
   TH1F *h_GJetInvtMass_binning40;
   TH1F *h_GJetInvtMass_VariableBinning;
   TH1F *h_GJetInvtMass_UnitBinning;
   TH1F *h_GJetdEta;
   TH1F *h_GJetdPhi;
   TH1F *h_GJetdR;

   //Histograms for CSVL BJets (0bTag)
   TH1F *h_CSVL_BJetPt_0bTag;
   TH1F *h_CSVL_BJetEta_0bTag;
   TH1F *h_CSVL_BJetPhi_0bTag;

   //Histograms for Photon-CSVLBJet (0bTag)
   TH1F *h_CSVL_GBJetInvtMass_binning40_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_0bTag;
   TH1F *h_CSVL_GBJetdEta_0bTag;
   TH1F *h_CSVL_GBJetdPhi_0bTag;
   TH1F *h_CSVL_GBJetdR_0bTag;

   //Histograms for CSVL BJets (1bTag)
   TH1F *h_CSVL_BJetPt_1bTag;
   TH1F *h_CSVL_BJetEta_1bTag;
   TH1F *h_CSVL_BJetPhi_1bTag;

   //Histograms for Photon-CSVLBJet (1bTag)
   TH1F *h_CSVL_GBJetInvtMass_binning40_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_1bTag;
   TH1F *h_CSVL_GBJetdEta_1bTag;
   TH1F *h_CSVL_GBJetdPhi_1bTag;
   TH1F *h_CSVL_GBJetdR_1bTag;

   //Histograms for CSVM BJets (0bTag)
   TH1F *h_CSVM_BJetPt_0bTag;
   TH1F *h_CSVM_BJetEta_0bTag;
   TH1F *h_CSVM_BJetPhi_0bTag;

   //Histograms for Photon-CSVMBJet (0bTag)
   TH1F *h_CSVM_GBJetInvtMass_binning40_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_0bTag;
   TH1F *h_CSVM_GBJetdEta_0bTag;
   TH1F *h_CSVM_GBJetdPhi_0bTag;
   TH1F *h_CSVM_GBJetdR_0bTag;

   //Histograms for CSVM BJets (1bTag)
   TH1F *h_CSVM_BJetPt_1bTag;
   TH1F *h_CSVM_BJetEta_1bTag;
   TH1F *h_CSVM_BJetPhi_1bTag;

   //Histograms for Photon-CSVMBJet (1bTag)
   TH1F *h_CSVM_GBJetInvtMass_binning40_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_1bTag;
   TH1F *h_CSVM_GBJetdEta_1bTag;
   TH1F *h_CSVM_GBJetdPhi_1bTag;
   TH1F *h_CSVM_GBJetdR_1bTag;

   //Histograms for CSVT BJets (0bTag)
   TH1F *h_CSVT_BJetPt_0bTag;
   TH1F *h_CSVT_BJetEta_0bTag;
   TH1F *h_CSVT_BJetPhi_0bTag;

   //Histograms for Photon-CSVTBJet (0bTag)
   TH1F *h_CSVT_GBJetInvtMass_binning40_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_0bTag;
   TH1F *h_CSVT_GBJetdEta_0bTag;
   TH1F *h_CSVT_GBJetdPhi_0bTag;
   TH1F *h_CSVT_GBJetdR_0bTag;

   //Histograms for CSVT BJets (1bTag)
   TH1F *h_CSVT_BJetPt_1bTag;
   TH1F *h_CSVT_BJetEta_1bTag;
   TH1F *h_CSVT_BJetPhi_1bTag;

   //Histograms for Photon-CSVTBJet (1bTag)
   TH1F *h_CSVT_GBJetInvtMass_binning40_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_1bTag;
   TH1F *h_CSVT_GBJetdEta_1bTag;
   TH1F *h_CSVT_GBJetdPhi_1bTag;
   TH1F *h_CSVT_GBJetdR_1bTag;

   //Histograms for BJetDisc before and after cut
   TH1F *h_BJetDiscByCSV;
   TH1F *h_BJetDiscByCSV_PassingCSVL;
   TH1F *h_BJetDiscByCSV_PassingCSVM;
   TH1F *h_BJetDiscByCSV_PassingCSVT;

   //Histograms for trigger efficiency
   TH1F *h_PassHLT150_EBPassProbes;
   TH1F *h_PassHLT75OR90_EBAllProbes;
   TH1F *h_PassHLT150_EEPassProbes;
   TH1F *h_PassHLT75OR90_EEAllProbes;

   //Histogram filling goodPV for pileup
   TH1F *h_goodPV;

   TH1F *h_nPhotons;
   TH1F *h_nJets;
   TH1F *h_nCSVLBJets;
   TH1F *h_nCSVMBJets;
   TH1F *h_nCSVTBJets;
   TH1F *h_CSVL_BJetsFrac;
   TH1F *h_CSVM_BJetsFrac;
   TH1F *h_CSVT_BJetsFrac;

   TH2F *h_PhotonIdxVsPt;
   TH2F *h_JetIdxVsPt;
   TH2F *h_CSVLBJetIdxVsPt;
   TH2F *h_CSVMBJetIdxVsPt;
   TH2F *h_CSVTBJetIdxVsPt;

   TH1F *h_PC, *h_JC;

   TH1F *h_CutFlowTable;

   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          LumiNumber;
   UInt_t          BXNumber;
   UInt_t          totalIntensityBeam1;
   UInt_t          totalIntensityBeam2;
   Float_t         avgInsDelLumi;
   Float_t         avgInsDelLumiErr;
   Float_t         avgInsRecLumi;
   Float_t         avgInsRecLumiErr;
   Int_t           Photon_n;
   vector<double>  *Photon_E;
   vector<double>  *Photon_et;
   vector<double>  *Photon_pt;
   vector<double>  *Photon_eta;
   vector<double>  *Photon_phi;
   vector<double>  *Photon_theta;
   vector<double>  *Photon_px;
   vector<double>  *Photon_py;
   vector<double>  *Photon_pz;
   vector<double>  *Photon_vx;
   vector<double>  *Photon_vy;
   vector<double>  *Photon_vz;
   vector<float>   *Photon_r9;
   vector<float>   *Photon_maxEnergyXtal;
   vector<float>   *Photon_e1x5;
   vector<float>   *Photon_e2x5;
   vector<float>   *Photon_e3x3;
   vector<float>   *Photon_e5x5;
   vector<float>   *Photon_r1x5;
   vector<float>   *Photon_r2x5;
   vector<float>   *Photon_SigmaEtaEta;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaEtaPhi;
   vector<float>   *Photon_SigmaIEtaIPhi;
   vector<float>   *Photon_SigmaPhiPhi;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_roundness;
   vector<float>   *Photon_angle;
   vector<float>   *Photon_swissCross;
   vector<float>   *Photon_s9;
   vector<float>   *Photon_e4Overe1;
   vector<float>   *Photon_e6Overe2;
   vector<float>   *Photon_e2Overe9;
   vector<float>   *Photon_rookFraction;
   vector<bool>    *Photon_isEB;
   vector<bool>    *Photon_isEE;
   vector<bool>    *Photon_isEBGap;
   vector<bool>    *Photon_isEEGap;
   vector<bool>    *Photon_isEBEEGap;
   vector<float>   *Photon_ecalRecHitSumEtConeDR03;
   vector<float>   *Photon_hcalTowerSumEtConeDR03;
   vector<float>   *Photon_hcalDepth1TowerSumEtConeDR03;
   vector<float>   *Photon_hcalDepth2TowerSumEtConeDR03;
   vector<float>   *Photon_trkSumPtSolidConeDR03;
   vector<float>   *Photon_trkSumPtHollowConeDR03;
   vector<int>     *Photon_nTrkSolidConeDR03;
   vector<int>     *Photon_nTrkHollowConeDR03;
   vector<float>   *Photon_ecalRecHitSumEtConeDR04;
   vector<float>   *Photon_hcalTowerSumEtConeDR04;
   vector<float>   *Photon_hcalDepth1TowerSumEtConeDR04;
   vector<float>   *Photon_hcalDepth2TowerSumEtConeDR04;
   vector<float>   *Photon_trkSumPtSolidConeDR04;
   vector<float>   *Photon_trkSumPtHollowConeDR04;
   vector<int>     *Photon_nTrkSolidConeDR04;
   vector<int>     *Photon_nTrkHollowConeDR04;
   vector<float>   *Photon_HoE;
   vector<float>   *Photon_SingleTowerHoE;
   vector<bool>    *Photon_hasConvTrk;
   vector<bool>    *Photon_hasPixelSeed;
   vector<bool>    *passedConvSafeElectronVeto;
   vector<int>     *Photon_SC_nOfBasicClusters;
   vector<double>  *Photon_SC_rawEnergy;
   vector<double>  *Photon_SC_preShowerEnergy;
   vector<double>  *Photon_SC_energy;
   vector<double>  *Photon_SC_eta;
   vector<double>  *Photon_SC_phi;
   vector<double>  *Photon_SC_x;
   vector<double>  *Photon_SC_y;
   vector<double>  *Photon_SC_z;
   vector<double>  *Photon_SC_etaWidth;
   vector<double>  *Photon_SC_phiWidth;
   vector<float>   *Photon_mipChi2;
   vector<float>   *Photon_mipTotEnergy;
   vector<float>   *Photon_mipSlope;
   vector<float>   *Photon_mipIntercept;
   vector<int>     *Photon_mipNhitCone;
   vector<bool>    *Photon_mipIsHalo;
   vector<unsigned int> *Photon_nConvTracks;
   vector<bool>    *Photon_isConverted;
   vector<float>   *Photon_pairInvariantMass;
   vector<float>   *Photon_pairCotThetaSeparation;
   vector<float>   *Photon_pairMomentum_x;
   vector<float>   *Photon_pairMomentum_y;
   vector<float>   *Photon_pairMomentum_z;
   vector<float>   *Photon_EoverP;
   vector<float>   *Photon_conv_vx;
   vector<float>   *Photon_conv_vy;
   vector<float>   *Photon_conv_vz;
   vector<float>   *Photon_zOfPrimaryVtxFromTrks;
   vector<float>   *Photon_distOfMinimumApproach;
   vector<float>   *Photon_dPhiTracksAtVtx;
   vector<float>   *Photon_dPhiTracksAtEcal;
   vector<float>   *Photon_dEtaTracksAtEcal;
   vector<int>     *Photon_nCrystals;
   vector<vector<float> > *Photon_xtal_timing;
   vector<vector<float> > *Photon_xtal_timeErr;
   vector<float>   *Photon_avgTimeAllxtals;
   vector<vector<float> > *Photon_xtal_energy;
   vector<vector<int> > *Photon_xtal_EBieta;
   vector<vector<int> > *Photon_xtal_EBiphi;
   vector<vector<int> > *Photon_xtal_EBrecoFlag;
   vector<double>  *PFIsoPhoton03;
   vector<double>  *PFIsoNeutral03;
   vector<double>  *PFIsoCharged03;
   vector<double>  *PFIsoSum03;
   vector<double>  *PFIsoChargedWorstvtx03;
   Int_t           PFPatJet_n;
   vector<double>  *PFPatJet_E;
   vector<double>  *PFPatJet_et;
   vector<double>  *PFPatJet_pt;
   vector<double>  *PFPatJet_eta;
   vector<double>  *PFPatJet_phi;
   vector<double>  *PFPatJet_theta;
   vector<double>  *PFPatJet_px;
   vector<double>  *PFPatJet_py;
   vector<double>  *PFPatJet_pz;
   vector<double>  *PFPatJet_vx;
   vector<double>  *PFPatJet_vy;
   vector<double>  *PFPatJet_vz;
   vector<float>   *PFPatJet_ChargedEmEnergy;
   vector<float>   *PFPatJet_ChargedEmEnergyFrac;
   vector<float>   *PFPatJet_ChargedHadEnergy;
   vector<float>   *PFPatJet_ChargedHadEnergyFrac;
   vector<int>     *PFPatJet_ChargedHadMult;
   vector<int>     *PFPatJet_ChargedMult;
   vector<int>     *PFPatJet_NConstituents;
   vector<float>   *PFPatJet_HFEMEnergy;
   vector<float>   *PFPatJet_HFEMEnergyFrac;
   vector<int>     *PFPatJet_HFEMMult;
   vector<float>   *PFPatJet_HFHadEnergy;
   vector<float>   *PFPatJet_HFHadEnergyFrac;
   vector<int>     *PFPatJet_HFHadMult;
   vector<float>   *PFPatJet_NeutralEmEnergy;
   vector<float>   *PFPatJet_NeutralEmEnergyFrac;
   vector<float>   *PFPatJet_NeutralHadEnergy;
   vector<float>   *PFPatJet_NeutralHadEnergyFrac;
   vector<int>     *PFPatJet_NeutralHadMult;
   vector<int>     *PFPatJet_NeutralMult;
   vector<double>  *PFPatJet_jecUncertainity;
   vector<float>   *PFPatJet_puJetIdCutBased_MVA;
   vector<float>   *PFPatJet_puJetIdSimple_MVA;
   vector<float>   *PFPatJet_puJetIdFull_MVA;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_loose;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_medium;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_tight;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_loose;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_medium;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_tight;
   vector<bool>    *PFPatJet_PassPUJetIdFull_loose;
   vector<bool>    *PFPatJet_PassPUJetIdFull_medium;
   vector<bool>    *PFPatJet_PassPUJetIdFull_tight;
   vector<float>   *PFPatJet_BJetDiscrByTrackCountingHighEff;
   vector<float>   *PFPatJet_BJetDiscrByTrackCountingHighPur;
   vector<float>   *PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;
   vector<float>   *PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;
   vector<float>   *PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;
   vector<float>   *PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA;
   vector<float>   *PFPatJet_BJetDiscrByjetProbabilityBJetTags;
   vector<float>   *PFPatJet_BJetDiscrByjetBProbabilityBJetTags;
   vector<float>   *PFPatJet_BJetDiscrBySoftElectronBJetTags;
   vector<float>   *PFPatJet_BJetDiscrBySoftMuonBJetTags;
   Int_t           PFRecoJet_n;
   vector<double>  *PFRecoJet_E;
   vector<double>  *PFRecoJet_et;
   vector<double>  *PFRecoJet_pt;
   vector<double>  *PFRecoJet_eta;
   vector<double>  *PFRecoJet_phi;
   vector<double>  *PFRecoJet_theta;
   vector<double>  *PFRecoJet_px;
   vector<double>  *PFRecoJet_py;
   vector<double>  *PFRecoJet_pz;
   vector<double>  *PFRecoJet_vx;
   vector<double>  *PFRecoJet_vy;
   vector<double>  *PFRecoJet_vz;
   vector<float>   *PFRecoJet_ChargedEmEnergy;
   vector<float>   *PFRecoJet_ChargedEmEnergyFrac;
   vector<float>   *PFRecoJet_ChargedHadEnergy;
   vector<float>   *PFRecoJet_ChargedHadEnergyFrac;
   vector<int>     *PFRecoJet_ChargedHadMult;
   vector<int>     *PFRecoJet_ChargedMult;
   vector<int>     *PFRecoJet_NConstituents;
   vector<float>   *PFRecoJet_HFEMEnergy;
   vector<float>   *PFRecoJet_HFEMEnergyFrac;
   vector<int>     *PFRecoJet_HFEMMult;
   vector<float>   *PFRecoJet_HFHadEnergy;
   vector<float>   *PFRecoJet_HFHadEnergyFrac;
   vector<int>     *PFRecoJet_HFHadMult;
   vector<float>   *PFRecoJet_NeutralEmEnergy;
   vector<float>   *PFRecoJet_NeutralEmEnergyFrac;
   vector<float>   *PFRecoJet_NeutralHadEnergy;
   vector<float>   *PFRecoJet_NeutralHadEnergyFrac;
   vector<int>     *PFRecoJet_NeutralHadMult;
   vector<int>     *PFRecoJet_NeutralMult;
   vector<double>  *PFRecoJet_jecUncertainity;
   vector<float>   *PFRecoJet_BJetDiscrByTrackCountingHighEff;
   vector<float>   *PFRecoJet_BJetDiscrByTrackCountingHighPur;
   vector<float>   *PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;
   vector<float>   *PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;
   vector<float>   *PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;
   vector<float>   *PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;
   vector<float>   *PFRecoJet_BJetDiscrByjetProbabilityBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrByjetBProbabilityBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrBySoftElectronBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrBySoftMuonBJetTags;
   Int_t           Vertex_n;
   vector<double>  *Vertex_x;
   vector<double>  *Vertex_y;
   vector<double>  *Vertex_z;
   vector<double>  *Vertex_chi2;
   vector<double>  *Vertex_nchi2;
   vector<double>  *Vertex_ndof;
   vector<int>     *Vertex_tracksSize;
   vector<bool>    *Vertex_isFake;
   vector<bool>    *Vertex_isValid;
   vector<double>  *Vertex_d0;
   Int_t           Tracks_n;
   vector<double>  *Track_pt;
   vector<double>  *Track_px;
   vector<double>  *Track_py;
   vector<double>  *Track_pz;
   vector<double>  *Track_vx;
   vector<double>  *Track_vy;
   vector<double>  *Track_vz;
   vector<double>  *Track_eta;
   vector<double>  *Track_phi;
   vector<double>  *Track_theta;
   vector<double>  *Track_chi2;
   Bool_t          IsScrapingEvent_;
   Float_t         Scraping_FractionOfGoodTracks;
   vector<string>  *HLT_Photon_triggers;
   vector<int>     *HLT_Photon_trig_prescales;
   vector<bool>    *HLT_Photon_ifTriggerPassed;
   Int_t           HLT_Photon_nTriggers;
   vector<int>     *HLT_Photon_triggerIndex;
   vector<int>     *HLT_Photon_nFilters;
   vector<string>  *HLT_Photon_FilterNames;
   vector<int>     *HLT_Photon_trigger_FilterStartPosition;
   vector<int>     *HLT_Photon_trigger_FilterEndPosition;
   vector<vector<double> > *HLT_Photon_FilterObjects_pt;
   vector<vector<double> > *HLT_Photon_FilterObjects_eta;
   vector<vector<double> > *HLT_Photon_FilterObjects_phi;
   Float_t         rho;
   Float_t         sigma;
   Float_t         rho25;
   Float_t         sigma25;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_BXNumber;   //!
   TBranch        *b_totalIntensityBeam1;   //!
   TBranch        *b_totalIntensityBeam2;   //!
   TBranch        *b_avgInsDelLumi;   //!
   TBranch        *b_avgInsDelLumiErr;   //!
   TBranch        *b_avgInsRecLumi;   //!
   TBranch        *b_avgInsRecLumiErr;   //!
   TBranch        *b_Photon_n;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_theta;   //!
   TBranch        *b_Photon_px;   //!
   TBranch        *b_Photon_py;   //!
   TBranch        *b_Photon_pz;   //!
   TBranch        *b_Photon_vx;   //!
   TBranch        *b_Photon_vy;   //!
   TBranch        *b_Photon_vz;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_maxEnergyXtal;   //!
   TBranch        *b_Photon_e1x5;   //!
   TBranch        *b_Photon_e2x5;   //!
   TBranch        *b_Photon_e3x3;   //!
   TBranch        *b_Photon_e5x5;   //!
   TBranch        *b_Photon_r1x5;   //!
   TBranch        *b_Photon_r2x5;   //!
   TBranch        *b_Photon_SigmaEtaEta;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_SigmaEtaPhi;   //!
   TBranch        *b_Photon_SigmaIEtaIPhi;   //!
   TBranch        *b_Photon_SigmaPhiPhi;   //!
   TBranch        *b_Photon_SigmaIPhiIPhi;   //!
   TBranch        *b_Photon_roundness;   //!
   TBranch        *b_Photon_angle;   //!
   TBranch        *b_Photon_swissCross;   //!
   TBranch        *b_Photon_s9;   //!
   TBranch        *b_Photon_e4Overe1;   //!
   TBranch        *b_Photon_e6Overe2;   //!
   TBranch        *b_Photon_e2Overe9;   //!
   TBranch        *b_Photon_rookFraction;   //!
   TBranch        *b_Photon_isEB;   //!
   TBranch        *b_Photon_isEE;   //!
   TBranch        *b_Photon_isEBGap;   //!
   TBranch        *b_Photon_isEEGap;   //!
   TBranch        *b_Photon_isEBEEGap;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_nTrkSolidConeDR03;   //!
   TBranch        *b_Photon_nTrkHollowConeDR03;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_nTrkSolidConeDR04;   //!
   TBranch        *b_Photon_nTrkHollowConeDR04;   //!
   TBranch        *b_Photon_HoE;   //!
   TBranch        *b_Photon_SingleTowerHoE;   //!
   TBranch        *b_Photon_hasConvTrk;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_passedConvSafeElectronVeto;   //!
   TBranch        *b_Photon_SC_nOfBasicClusters;   //!
   TBranch        *b_Photon_SC_rawEnergy;   //!
   TBranch        *b_Photon_SC_preShowerEnergy;   //!
   TBranch        *b_Photon_SC_energy;   //!
   TBranch        *b_Photon_SC_eta;   //!
   TBranch        *b_Photon_SC_phi;   //!
   TBranch        *b_Photon_SC_x;   //!
   TBranch        *b_Photon_SC_y;   //!
   TBranch        *b_Photon_SC_z;   //!
   TBranch        *b_Photon_SC_etaWidth;   //!
   TBranch        *b_Photon_SC_phiWidth;   //!
   TBranch        *b_Photon_mipChi2;   //!
   TBranch        *b_Photon_mipTotEnergy;   //!
   TBranch        *b_Photon_mipSlope;   //!
   TBranch        *b_Photon_mipIntercept;   //!
   TBranch        *b_Photon_mipNhitCone;   //!
   TBranch        *b_Photon_mipIsHalo;   //!
   TBranch        *b_Photon_nConvTracks;   //!
   TBranch        *b_Photon_isConverted;   //!
   TBranch        *b_Photon_pairInvariantMass;   //!
   TBranch        *b_Photon_pairCotThetaSeparation;   //!
   TBranch        *b_Photon_pairMomentum_x;   //!
   TBranch        *b_Photon_pairMomentum_y;   //!
   TBranch        *b_Photon_pairMomentum_z;   //!
   TBranch        *b_Photon_EoverP;   //!
   TBranch        *b_Photon_conv_vx;   //!
   TBranch        *b_Photon_conv_vy;   //!
   TBranch        *b_Photon_conv_vz;   //!
   TBranch        *b_Photon_zOfPrimaryVtxFromTrks;   //!
   TBranch        *b_Photon_distOfMinimumApproach;   //!
   TBranch        *b_Photon_dPhiTracksAtVtx;   //!
   TBranch        *b_Photon_dPhiTracksAtEcal;   //!
   TBranch        *b_Photon_dEtaTracksAtEcal;   //!
   TBranch        *b_Photon_nCrystals;   //!
   TBranch        *b_Photon_xtal_timing;   //!
   TBranch        *b_Photon_xtal_timeErr;   //!
   TBranch        *b_Photon_avgTimeAllxtals;   //!
   TBranch        *b_Photon_xtal_energy;   //!
   TBranch        *b_Photon_xtal_EBieta;   //!
   TBranch        *b_Photon_xtal_EBiphi;   //!
   TBranch        *b_Photon_xtal_EBrecoFlag;   //!
   TBranch        *b_PFIsoPhoton03;   //!
   TBranch        *b_PFIsoNeutral03;   //!
   TBranch        *b_PFIsoCharged03;   //!
   TBranch        *b_PFIsoSum03;   //!
   TBranch        *b_PFIsoChargedWorstvtx03;   //!
   TBranch        *b_PFPatJet_n;   //!
   TBranch        *b_PFPatJet_E;   //!
   TBranch        *b_PFPatJet_et;   //!
   TBranch        *b_PFPatJet_pt;   //!
   TBranch        *b_PFPatJet_eta;   //!
   TBranch        *b_PFPatJet_phi;   //!
   TBranch        *b_PFPatJet_theta;   //!
   TBranch        *b_PFPatJet_px;   //!
   TBranch        *b_PFPatJet_py;   //!
   TBranch        *b_PFPatJet_pz;   //!
   TBranch        *b_PFPatJet_vx;   //!
   TBranch        *b_PFPatJet_vy;   //!
   TBranch        *b_PFPatJet_vz;   //!
   TBranch        *b_PFPatJet_ChargedEmEnergy;   //!
   TBranch        *b_PFPatJet_ChargedEmEnergyFrac;   //!
   TBranch        *b_PFPatJet_ChargedHadEnergy;   //!
   TBranch        *b_PFPatJet_ChargedHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_ChargedHadMult;   //!
   TBranch        *b_PFPatJet_ChargedMult;   //!
   TBranch        *b_PFPatJet_NConstituents;   //!
   TBranch        *b_PFPatJet_HFEMEnergy;   //!
   TBranch        *b_PFPatJet_HFEMEnergyFrac;   //!
   TBranch        *b_PFPatJet_HFEMMult;   //!
   TBranch        *b_PFPatJet_HFHadEnergy;   //!
   TBranch        *b_PFPatJet_HFHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_HFHadMult;   //!
   TBranch        *b_PFPatJet_NeutralEmEnergy;   //!
   TBranch        *b_PFPatJet_NeutralEmEnergyFrac;   //!
   TBranch        *b_PFPatJet_NeutralHadEnergy;   //!
   TBranch        *b_PFPatJet_NeutralHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_NeutralHadMult;   //!
   TBranch        *b_PFPatJet_NeutralMult;   //!
   TBranch        *b_PFPatJet_jecUncertainity;   //!
   TBranch        *b_PFPatJet_puJetIdCutBased_MVA;   //!
   TBranch        *b_PFPatJet_puJetIdSimple_MVA;   //!
   TBranch        *b_PFPatJet_puJetIdFull_MVA;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_tight;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_tight;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_tight;   //!
   TBranch        *b_PFPatJet_BJetDiscrByTrackCountingHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrByTrackCountingHighPur;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;   //!
   TBranch        *b_PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA;   //!
   TBranch        *b_PFPatJet_BJetDiscrByjetProbabilityBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrByjetBProbabilityBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySoftElectronBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySoftMuonBJetTags;   //!
   TBranch        *b_PFRecoJet_n;   //!
   TBranch        *b_PFRecoJet_E;   //!
   TBranch        *b_PFRecoJet_et;   //!
   TBranch        *b_PFRecoJet_pt;   //!
   TBranch        *b_PFRecoJet_eta;   //!
   TBranch        *b_PFRecoJet_phi;   //!
   TBranch        *b_PFRecoJet_theta;   //!
   TBranch        *b_PFRecoJet_px;   //!
   TBranch        *b_PFRecoJet_py;   //!
   TBranch        *b_PFRecoJet_pz;   //!
   TBranch        *b_PFRecoJet_vx;   //!
   TBranch        *b_PFRecoJet_vy;   //!
   TBranch        *b_PFRecoJet_vz;   //!
   TBranch        *b_PFRecoJet_ChargedEmEnergy;   //!
   TBranch        *b_PFRecoJet_ChargedEmEnergyFrac;   //!
   TBranch        *b_PFRecoJet_ChargedHadEnergy;   //!
   TBranch        *b_PFRecoJet_ChargedHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_ChargedHadMult;   //!
   TBranch        *b_PFRecoJet_ChargedMult;   //!
   TBranch        *b_PFRecoJet_NConstituents;   //!
   TBranch        *b_PFRecoJet_HFEMEnergy;   //!
   TBranch        *b_PFRecoJet_HFEMEnergyFrac;   //!
   TBranch        *b_PFRecoJet_HFEMMult;   //!
   TBranch        *b_PFRecoJet_HFHadEnergy;   //!
   TBranch        *b_PFRecoJet_HFHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_HFHadMult;   //!
   TBranch        *b_PFRecoJet_NeutralEmEnergy;   //!
   TBranch        *b_PFRecoJet_NeutralEmEnergyFrac;   //!
   TBranch        *b_PFRecoJet_NeutralHadEnergy;   //!
   TBranch        *b_PFRecoJet_NeutralHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_NeutralHadMult;   //!
   TBranch        *b_PFRecoJet_NeutralMult;   //!
   TBranch        *b_PFRecoJet_jecUncertainity;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByTrackCountingHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByTrackCountingHighPur;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByjetProbabilityBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByjetBProbabilityBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySoftElectronBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySoftMuonBJetTags;   //!
   TBranch        *b_Vertex_n;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_nchi2;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_tracksSize;   //!
   TBranch        *b_Vertex_isFake;   //!
   TBranch        *b_Vertex_isValid;   //!
   TBranch        *b_Vertex_d0;   //!
   TBranch        *b_Tracks_n;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_Track_theta;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_IsScrapingEvent_;   //!
   TBranch        *b_Scraping_FractionOfGoodTracks;   //!
   TBranch        *b_HLT_Photon_triggers;   //!
   TBranch        *b_HLT_Photon_trig_prescales;   //!
   TBranch        *b_HLT_Photon_ifTriggerPassed;   //!
   TBranch        *b_HLT_Photon_nTriggers;   //!
   TBranch        *b_HLT_Photon_triggerIndex;   //!
   TBranch        *b_HLT_Photon_nFilters;   //!
   TBranch        *b_HLT_Photon_FilterNames;   //!
   TBranch        *b_HLT_Photon_trigger_FilterStartPosition;   //!
   TBranch        *b_HLT_Photon_trigger_FilterEndPosition;   //!
   TBranch        *b_HLT_Photon_FilterObjects_pt;   //!
   TBranch        *b_HLT_Photon_FilterObjects_eta;   //!
   TBranch        *b_HLT_Photon_FilterObjects_phi;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_sigma;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_sigma25;   //!

   PostAnalyzer_Data(TTree *tree=0);
   virtual ~PostAnalyzer_Data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);  

   //Self Declared Functions
   virtual Bool_t   PassHLT();
   virtual Bool_t   PassHLT_Photon75();
   virtual Bool_t   PassHLT_Photon90();
   virtual Bool_t   NonScrapingEvt();
   virtual Int_t    GoodPrimaryVtx();
   virtual Bool_t   ResSpikes(Int_t);
   virtual Bool_t   ResSpikes_NoLICTD(Int_t);
   virtual Double_t GetLICTD(Int_t);
   virtual Double_t EAChargedHadrons(Double_t);
   virtual Double_t EANeutralHadrons(Double_t);
   virtual Double_t EAPhotons(Double_t);
   virtual Bool_t   TightPhotonIdnPFIso(Int_t);
   virtual Bool_t   MediumPhotonIdnPFIso(Int_t);
   virtual Bool_t   LoosePhotonIdnPFIso(Int_t);
   virtual Bool_t   TightJetId(Int_t);
   virtual Double_t GetdEta(Double_t, Double_t);
   virtual Double_t GetdPhi(Double_t, Double_t);
   virtual Double_t GetdR(Double_t, Double_t, Double_t, Double_t);
   virtual Double_t GetInvtMass(Int_t, Int_t);
   virtual Int_t    GetPhotonPassingAllCuts_WithTightID();
   virtual Int_t    GetPhotonPassingAllCuts_WithMediumID();
   virtual Int_t    GetPhotonPassingAllCuts_WithLooseID();
   virtual Int_t    GetPhotonPassingAllCuts_WithTightID_NoLICTD();
   virtual Int_t    GetPhotonPassingAllCuts_WithMediumID_NoLICTD();
   virtual Int_t    GetPhotonPassingAllCuts_WithLooseID_NoLICTD();
   virtual Int_t    GetJetPassingIDnMatchedToPhoton(Int_t);
   virtual Bool_t   PassCSVLBTag(Int_t);
   virtual Bool_t   PassCSVMBTag(Int_t);
   virtual Bool_t   PassCSVTBTag(Int_t);
   virtual void     BookHistograms();
};

#endif

#ifdef PostAnalyzer_Data_cxx
PostAnalyzer_Data::PostAnalyzer_Data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("myEvent",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("myEvent","");
/*
      //Comment out this part in script
      //      chain->Add("/eos/uscms/store/user/rocky86/NTuples/Data/Data22Jan2013ReReco/Data_Run2012A_22Jan_AOD_NTuples/Data_Run2012A_22Jan_AOD_NTuples_39_1_6Oc.root/myEvent");
      chain->Add("/eos/uscms/store/user/rocky86/NTuples/Data/Data22Jan2013ReReco/Data_Run2012A_22Jan_AOD_NTuples/Data_Run2012A_22Jan_AOD_NTuples_412_1_Mzi.root/myEvent");
*/
      //Uncomment this part in script
      
      ///--------------------Use this part while submitting job in one go for a dataset--------------///
      TString main_path = "${sourceDir}";

      TSystemDirectory sourceDir("sysDir",main_path);
      TList* fileList = sourceDir.GetListOfFiles();
      TIter next(fileList);
      TSystemFile* fileName;

      int fileNumber = 1;
      int maxFiles = -1;

      while ((fileName = (TSystemFile*)next()) && fileNumber > maxFiles){
        if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."){continue;}

        TString FullPathInputFile = (main_path+fileName->GetName());

        //cout << FullPathInputFile << endl;

        chain->Add(FullPathInputFile+"/myEvent");

        fileNumber++;

      }

      cout << "Total files in this set = " << fileNumber - 1 << endl;
      ///-------------------------------------------------------------------------------------------------///
/*
      ///-----------Only change Tstring part with this part for job submission in parts for a dataset-----///
      ifstream dataset;
      dataset.open("${subdataset}", ifstream::in);
      char datafilename[300];

      int a = 0;
      for(Int_t i = ${sf}; i <= ${maxf} && i <= ${Total_files}; i++){      
	dataset >> datafilename;   ////dataset >> datafilename always starts from the 1st line if the file is just opened and will start from the 
                                   //// next to last line if already opened.
	string fname(datafilename);
	string main_path = "${sourceDir}";
	string FullPathInputFile = (main_path+datafilename);
	string myevt = "/myEvent";

	chain->Add((FullPathInputFile+myevt).c_str());
	a++;

      //cout << (FullPathInputFile+myevt).c_str() << endl;

      }
      cout << "Total Files in this job are = " << a << endl;
      ///------------------------------------------------------------------------------------------------///


      */

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PostAnalyzer_Data::~PostAnalyzer_Data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   file->cd();
   file->Write();
   file->Close();
}

Int_t PostAnalyzer_Data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PostAnalyzer_Data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PostAnalyzer_Data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Photon_E = 0;
   Photon_et = 0;
   Photon_pt = 0;
   Photon_eta = 0;
   Photon_phi = 0;
   Photon_theta = 0;
   Photon_px = 0;
   Photon_py = 0;
   Photon_pz = 0;
   Photon_vx = 0;
   Photon_vy = 0;
   Photon_vz = 0;
   Photon_r9 = 0;
   Photon_maxEnergyXtal = 0;
   Photon_e1x5 = 0;
   Photon_e2x5 = 0;
   Photon_e3x3 = 0;
   Photon_e5x5 = 0;
   Photon_r1x5 = 0;
   Photon_r2x5 = 0;
   Photon_SigmaEtaEta = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_SigmaEtaPhi = 0;
   Photon_SigmaIEtaIPhi = 0;
   Photon_SigmaPhiPhi = 0;
   Photon_SigmaIPhiIPhi = 0;
   Photon_roundness = 0;
   Photon_angle = 0;
   Photon_swissCross = 0;
   Photon_s9 = 0;
   Photon_e4Overe1 = 0;
   Photon_e6Overe2 = 0;
   Photon_e2Overe9 = 0;
   Photon_rookFraction = 0;
   Photon_isEB = 0;
   Photon_isEE = 0;
   Photon_isEBGap = 0;
   Photon_isEEGap = 0;
   Photon_isEBEEGap = 0;
   Photon_ecalRecHitSumEtConeDR03 = 0;
   Photon_hcalTowerSumEtConeDR03 = 0;
   Photon_hcalDepth1TowerSumEtConeDR03 = 0;
   Photon_hcalDepth2TowerSumEtConeDR03 = 0;
   Photon_trkSumPtSolidConeDR03 = 0;
   Photon_trkSumPtHollowConeDR03 = 0;
   Photon_nTrkSolidConeDR03 = 0;
   Photon_nTrkHollowConeDR03 = 0;
   Photon_ecalRecHitSumEtConeDR04 = 0;
   Photon_hcalTowerSumEtConeDR04 = 0;
   Photon_hcalDepth1TowerSumEtConeDR04 = 0;
   Photon_hcalDepth2TowerSumEtConeDR04 = 0;
   Photon_trkSumPtSolidConeDR04 = 0;
   Photon_trkSumPtHollowConeDR04 = 0;
   Photon_nTrkSolidConeDR04 = 0;
   Photon_nTrkHollowConeDR04 = 0;
   Photon_HoE = 0;
   Photon_SingleTowerHoE = 0;
   Photon_hasConvTrk = 0;
   Photon_hasPixelSeed = 0;
   passedConvSafeElectronVeto = 0;
   Photon_SC_nOfBasicClusters = 0;
   Photon_SC_rawEnergy = 0;
   Photon_SC_preShowerEnergy = 0;
   Photon_SC_energy = 0;
   Photon_SC_eta = 0;
   Photon_SC_phi = 0;
   Photon_SC_x = 0;
   Photon_SC_y = 0;
   Photon_SC_z = 0;
   Photon_SC_etaWidth = 0;
   Photon_SC_phiWidth = 0;
   Photon_mipChi2 = 0;
   Photon_mipTotEnergy = 0;
   Photon_mipSlope = 0;
   Photon_mipIntercept = 0;
   Photon_mipNhitCone = 0;
   Photon_mipIsHalo = 0;
   Photon_nConvTracks = 0;
   Photon_isConverted = 0;
   Photon_pairInvariantMass = 0;
   Photon_pairCotThetaSeparation = 0;
   Photon_pairMomentum_x = 0;
   Photon_pairMomentum_y = 0;
   Photon_pairMomentum_z = 0;
   Photon_EoverP = 0;
   Photon_conv_vx = 0;
   Photon_conv_vy = 0;
   Photon_conv_vz = 0;
   Photon_zOfPrimaryVtxFromTrks = 0;
   Photon_distOfMinimumApproach = 0;
   Photon_dPhiTracksAtVtx = 0;
   Photon_dPhiTracksAtEcal = 0;
   Photon_dEtaTracksAtEcal = 0;
   Photon_nCrystals = 0;
   Photon_xtal_timing = 0;
   Photon_xtal_timeErr = 0;
   Photon_avgTimeAllxtals = 0;
   Photon_xtal_energy = 0;
   Photon_xtal_EBieta = 0;
   Photon_xtal_EBiphi = 0;
   Photon_xtal_EBrecoFlag = 0;
   PFIsoPhoton03 = 0;
   PFIsoNeutral03 = 0;
   PFIsoCharged03 = 0;
   PFIsoSum03 = 0;
   PFIsoChargedWorstvtx03 = 0;
   PFPatJet_E = 0;
   PFPatJet_et = 0;
   PFPatJet_pt = 0;
   PFPatJet_eta = 0;
   PFPatJet_phi = 0;
   PFPatJet_theta = 0;
   PFPatJet_px = 0;
   PFPatJet_py = 0;
   PFPatJet_pz = 0;
   PFPatJet_vx = 0;
   PFPatJet_vy = 0;
   PFPatJet_vz = 0;
   PFPatJet_ChargedEmEnergy = 0;
   PFPatJet_ChargedEmEnergyFrac = 0;
   PFPatJet_ChargedHadEnergy = 0;
   PFPatJet_ChargedHadEnergyFrac = 0;
   PFPatJet_ChargedHadMult = 0;
   PFPatJet_ChargedMult = 0;
   PFPatJet_NConstituents = 0;
   PFPatJet_HFEMEnergy = 0;
   PFPatJet_HFEMEnergyFrac = 0;
   PFPatJet_HFEMMult = 0;
   PFPatJet_HFHadEnergy = 0;
   PFPatJet_HFHadEnergyFrac = 0;
   PFPatJet_HFHadMult = 0;
   PFPatJet_NeutralEmEnergy = 0;
   PFPatJet_NeutralEmEnergyFrac = 0;
   PFPatJet_NeutralHadEnergy = 0;
   PFPatJet_NeutralHadEnergyFrac = 0;
   PFPatJet_NeutralHadMult = 0;
   PFPatJet_NeutralMult = 0;
   PFPatJet_jecUncertainity = 0;
   PFPatJet_puJetIdCutBased_MVA = 0;
   PFPatJet_puJetIdSimple_MVA = 0;
   PFPatJet_puJetIdFull_MVA = 0;
   PFPatJet_PassPUJetIdCutBased_loose = 0;
   PFPatJet_PassPUJetIdCutBased_medium = 0;
   PFPatJet_PassPUJetIdCutBased_tight = 0;
   PFPatJet_PassPUJetIdSimple_loose = 0;
   PFPatJet_PassPUJetIdSimple_medium = 0;
   PFPatJet_PassPUJetIdSimple_tight = 0;
   PFPatJet_PassPUJetIdFull_loose = 0;
   PFPatJet_PassPUJetIdFull_medium = 0;
   PFPatJet_PassPUJetIdFull_tight = 0;
   PFPatJet_BJetDiscrByTrackCountingHighEff = 0;
   PFPatJet_BJetDiscrByTrackCountingHighPur = 0;
   PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff = 0;
   PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur = 0;
   PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff = 0;
   PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA = 0;
   PFPatJet_BJetDiscrByjetProbabilityBJetTags = 0;
   PFPatJet_BJetDiscrByjetBProbabilityBJetTags = 0;
   PFPatJet_BJetDiscrBySoftElectronBJetTags = 0;
   PFPatJet_BJetDiscrBySoftMuonBJetTags = 0;
   PFRecoJet_E = 0;
   PFRecoJet_et = 0;
   PFRecoJet_pt = 0;
   PFRecoJet_eta = 0;
   PFRecoJet_phi = 0;
   PFRecoJet_theta = 0;
   PFRecoJet_px = 0;
   PFRecoJet_py = 0;
   PFRecoJet_pz = 0;
   PFRecoJet_vx = 0;
   PFRecoJet_vy = 0;
   PFRecoJet_vz = 0;
   PFRecoJet_ChargedEmEnergy = 0;
   PFRecoJet_ChargedEmEnergyFrac = 0;
   PFRecoJet_ChargedHadEnergy = 0;
   PFRecoJet_ChargedHadEnergyFrac = 0;
   PFRecoJet_ChargedHadMult = 0;
   PFRecoJet_ChargedMult = 0;
   PFRecoJet_NConstituents = 0;
   PFRecoJet_HFEMEnergy = 0;
   PFRecoJet_HFEMEnergyFrac = 0;
   PFRecoJet_HFEMMult = 0;
   PFRecoJet_HFHadEnergy = 0;
   PFRecoJet_HFHadEnergyFrac = 0;
   PFRecoJet_HFHadMult = 0;
   PFRecoJet_NeutralEmEnergy = 0;
   PFRecoJet_NeutralEmEnergyFrac = 0;
   PFRecoJet_NeutralHadEnergy = 0;
   PFRecoJet_NeutralHadEnergyFrac = 0;
   PFRecoJet_NeutralHadMult = 0;
   PFRecoJet_NeutralMult = 0;
   PFRecoJet_jecUncertainity = 0;
   PFRecoJet_BJetDiscrByTrackCountingHighEff = 0;
   PFRecoJet_BJetDiscrByTrackCountingHighPur = 0;
   PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff = 0;
   PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur = 0;
   PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff = 0;
   PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA = 0;
   PFRecoJet_BJetDiscrByjetProbabilityBJetTags = 0;
   PFRecoJet_BJetDiscrByjetBProbabilityBJetTags = 0;
   PFRecoJet_BJetDiscrBySoftElectronBJetTags = 0;
   PFRecoJet_BJetDiscrBySoftMuonBJetTags = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   Vertex_chi2 = 0;
   Vertex_nchi2 = 0;
   Vertex_ndof = 0;
   Vertex_tracksSize = 0;
   Vertex_isFake = 0;
   Vertex_isValid = 0;
   Vertex_d0 = 0;
   Track_pt = 0;
   Track_px = 0;
   Track_py = 0;
   Track_pz = 0;
   Track_vx = 0;
   Track_vy = 0;
   Track_vz = 0;
   Track_eta = 0;
   Track_phi = 0;
   Track_theta = 0;
   Track_chi2 = 0;
   HLT_Photon_triggers = 0;
   HLT_Photon_trig_prescales = 0;
   HLT_Photon_ifTriggerPassed = 0;
   HLT_Photon_triggerIndex = 0;
   HLT_Photon_nFilters = 0;
   HLT_Photon_FilterNames = 0;
   HLT_Photon_trigger_FilterStartPosition = 0;
   HLT_Photon_trigger_FilterEndPosition = 0;
   HLT_Photon_FilterObjects_pt = 0;
   HLT_Photon_FilterObjects_eta = 0;
   HLT_Photon_FilterObjects_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("BXNumber", &BXNumber, &b_BXNumber);
   fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
   fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
   fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
   fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
   fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
   fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
   fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
   fChain->SetBranchAddress("Photon_E", &Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_theta", &Photon_theta, &b_Photon_theta);
   fChain->SetBranchAddress("Photon_px", &Photon_px, &b_Photon_px);
   fChain->SetBranchAddress("Photon_py", &Photon_py, &b_Photon_py);
   fChain->SetBranchAddress("Photon_pz", &Photon_pz, &b_Photon_pz);
   fChain->SetBranchAddress("Photon_vx", &Photon_vx, &b_Photon_vx);
   fChain->SetBranchAddress("Photon_vy", &Photon_vy, &b_Photon_vy);
   fChain->SetBranchAddress("Photon_vz", &Photon_vz, &b_Photon_vz);
   fChain->SetBranchAddress("Photon_r9", &Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_maxEnergyXtal", &Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
   fChain->SetBranchAddress("Photon_e1x5", &Photon_e1x5, &b_Photon_e1x5);
   fChain->SetBranchAddress("Photon_e2x5", &Photon_e2x5, &b_Photon_e2x5);
   fChain->SetBranchAddress("Photon_e3x3", &Photon_e3x3, &b_Photon_e3x3);
   fChain->SetBranchAddress("Photon_e5x5", &Photon_e5x5, &b_Photon_e5x5);
   fChain->SetBranchAddress("Photon_r1x5", &Photon_r1x5, &b_Photon_r1x5);
   fChain->SetBranchAddress("Photon_r2x5", &Photon_r2x5, &b_Photon_r2x5);
   fChain->SetBranchAddress("Photon_SigmaEtaEta", &Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_SigmaEtaPhi", &Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
   fChain->SetBranchAddress("Photon_SigmaIEtaIPhi", &Photon_SigmaIEtaIPhi, &b_Photon_SigmaIEtaIPhi);
   fChain->SetBranchAddress("Photon_SigmaPhiPhi", &Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
   fChain->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
   fChain->SetBranchAddress("Photon_roundness", &Photon_roundness, &b_Photon_roundness);
   fChain->SetBranchAddress("Photon_angle", &Photon_angle, &b_Photon_angle);
   fChain->SetBranchAddress("Photon_swissCross", &Photon_swissCross, &b_Photon_swissCross);
   fChain->SetBranchAddress("Photon_s9", &Photon_s9, &b_Photon_s9);
   fChain->SetBranchAddress("Photon_e4Overe1", &Photon_e4Overe1, &b_Photon_e4Overe1);
   fChain->SetBranchAddress("Photon_e6Overe2", &Photon_e6Overe2, &b_Photon_e6Overe2);
   fChain->SetBranchAddress("Photon_e2Overe9", &Photon_e2Overe9, &b_Photon_e2Overe9);
   fChain->SetBranchAddress("Photon_rookFraction", &Photon_rookFraction, &b_Photon_rookFraction);
   fChain->SetBranchAddress("Photon_isEB", &Photon_isEB, &b_Photon_isEB);
   fChain->SetBranchAddress("Photon_isEE", &Photon_isEE, &b_Photon_isEE);
   fChain->SetBranchAddress("Photon_isEBGap", &Photon_isEBGap, &b_Photon_isEBGap);
   fChain->SetBranchAddress("Photon_isEEGap", &Photon_isEEGap, &b_Photon_isEEGap);
   fChain->SetBranchAddress("Photon_isEBEEGap", &Photon_isEBEEGap, &b_Photon_isEBEEGap);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", &Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", &Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", &Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", &Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", &Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", &Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", &Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", &Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", &Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", &Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", &Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", &Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", &Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", &Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", &Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", &Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
   fChain->SetBranchAddress("Photon_HoE", &Photon_HoE, &b_Photon_HoE);
   fChain->SetBranchAddress("Photon_SingleTowerHoE", &Photon_SingleTowerHoE, &b_Photon_SingleTowerHoE);
   fChain->SetBranchAddress("Photon_hasConvTrk", &Photon_hasConvTrk, &b_Photon_hasConvTrk);
   fChain->SetBranchAddress("Photon_hasPixelSeed", &Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   fChain->SetBranchAddress("passedConvSafeElectronVeto", &passedConvSafeElectronVeto, &b_passedConvSafeElectronVeto);
   fChain->SetBranchAddress("Photon_SC_nOfBasicClusters", &Photon_SC_nOfBasicClusters, &b_Photon_SC_nOfBasicClusters);
   fChain->SetBranchAddress("Photon_SC_rawEnergy", &Photon_SC_rawEnergy, &b_Photon_SC_rawEnergy);
   fChain->SetBranchAddress("Photon_SC_preShowerEnergy", &Photon_SC_preShowerEnergy, &b_Photon_SC_preShowerEnergy);
   fChain->SetBranchAddress("Photon_SC_energy", &Photon_SC_energy, &b_Photon_SC_energy);
   fChain->SetBranchAddress("Photon_SC_eta", &Photon_SC_eta, &b_Photon_SC_eta);
   fChain->SetBranchAddress("Photon_SC_phi", &Photon_SC_phi, &b_Photon_SC_phi);
   fChain->SetBranchAddress("Photon_SC_x", &Photon_SC_x, &b_Photon_SC_x);
   fChain->SetBranchAddress("Photon_SC_y", &Photon_SC_y, &b_Photon_SC_y);
   fChain->SetBranchAddress("Photon_SC_z", &Photon_SC_z, &b_Photon_SC_z);
   fChain->SetBranchAddress("Photon_SC_etaWidth", &Photon_SC_etaWidth, &b_Photon_SC_etaWidth);
   fChain->SetBranchAddress("Photon_SC_phiWidth", &Photon_SC_phiWidth, &b_Photon_SC_phiWidth);
   fChain->SetBranchAddress("Photon_mipChi2", &Photon_mipChi2, &b_Photon_mipChi2);
   fChain->SetBranchAddress("Photon_mipTotEnergy", &Photon_mipTotEnergy, &b_Photon_mipTotEnergy);
   fChain->SetBranchAddress("Photon_mipSlope", &Photon_mipSlope, &b_Photon_mipSlope);
   fChain->SetBranchAddress("Photon_mipIntercept", &Photon_mipIntercept, &b_Photon_mipIntercept);
   fChain->SetBranchAddress("Photon_mipNhitCone", &Photon_mipNhitCone, &b_Photon_mipNhitCone);
   fChain->SetBranchAddress("Photon_mipIsHalo", &Photon_mipIsHalo, &b_Photon_mipIsHalo);
   fChain->SetBranchAddress("Photon_nConvTracks", &Photon_nConvTracks, &b_Photon_nConvTracks);
   fChain->SetBranchAddress("Photon_isConverted", &Photon_isConverted, &b_Photon_isConverted);
   fChain->SetBranchAddress("Photon_pairInvariantMass", &Photon_pairInvariantMass, &b_Photon_pairInvariantMass);
   fChain->SetBranchAddress("Photon_pairCotThetaSeparation", &Photon_pairCotThetaSeparation, &b_Photon_pairCotThetaSeparation);
   fChain->SetBranchAddress("Photon_pairMomentum_x", &Photon_pairMomentum_x, &b_Photon_pairMomentum_x);
   fChain->SetBranchAddress("Photon_pairMomentum_y", &Photon_pairMomentum_y, &b_Photon_pairMomentum_y);
   fChain->SetBranchAddress("Photon_pairMomentum_z", &Photon_pairMomentum_z, &b_Photon_pairMomentum_z);
   fChain->SetBranchAddress("Photon_EoverP", &Photon_EoverP, &b_Photon_EoverP);
   fChain->SetBranchAddress("Photon_conv_vx", &Photon_conv_vx, &b_Photon_conv_vx);
   fChain->SetBranchAddress("Photon_conv_vy", &Photon_conv_vy, &b_Photon_conv_vy);
   fChain->SetBranchAddress("Photon_conv_vz", &Photon_conv_vz, &b_Photon_conv_vz);
   fChain->SetBranchAddress("Photon_zOfPrimaryVtxFromTrks", &Photon_zOfPrimaryVtxFromTrks, &b_Photon_zOfPrimaryVtxFromTrks);
   fChain->SetBranchAddress("Photon_distOfMinimumApproach", &Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
   fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", &Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
   fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", &Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
   fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", &Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
   fChain->SetBranchAddress("Photon_nCrystals", &Photon_nCrystals, &b_Photon_nCrystals);
   fChain->SetBranchAddress("Photon_xtal_timing", &Photon_xtal_timing, &b_Photon_xtal_timing);
   fChain->SetBranchAddress("Photon_xtal_timeErr", &Photon_xtal_timeErr, &b_Photon_xtal_timeErr);
   fChain->SetBranchAddress("Photon_avgTimeAllxtals", &Photon_avgTimeAllxtals, &b_Photon_avgTimeAllxtals);
   fChain->SetBranchAddress("Photon_xtal_energy", &Photon_xtal_energy, &b_Photon_xtal_energy);
   fChain->SetBranchAddress("Photon_xtal_EBieta", &Photon_xtal_EBieta, &b_Photon_xtal_EBieta);
   fChain->SetBranchAddress("Photon_xtal_EBiphi", &Photon_xtal_EBiphi, &b_Photon_xtal_EBiphi);
   fChain->SetBranchAddress("Photon_xtal_EBrecoFlag", &Photon_xtal_EBrecoFlag, &b_Photon_xtal_EBrecoFlag);
   fChain->SetBranchAddress("PFIsoPhoton03", &PFIsoPhoton03, &b_PFIsoPhoton03);
   fChain->SetBranchAddress("PFIsoNeutral03", &PFIsoNeutral03, &b_PFIsoNeutral03);
   fChain->SetBranchAddress("PFIsoCharged03", &PFIsoCharged03, &b_PFIsoCharged03);
   fChain->SetBranchAddress("PFIsoSum03", &PFIsoSum03, &b_PFIsoSum03);
   fChain->SetBranchAddress("PFIsoChargedWorstvtx03", &PFIsoChargedWorstvtx03, &b_PFIsoChargedWorstvtx03);
   fChain->SetBranchAddress("PFPatJet_n", &PFPatJet_n, &b_PFPatJet_n);
   fChain->SetBranchAddress("PFPatJet_E", &PFPatJet_E, &b_PFPatJet_E);
   fChain->SetBranchAddress("PFPatJet_et", &PFPatJet_et, &b_PFPatJet_et);
   fChain->SetBranchAddress("PFPatJet_pt", &PFPatJet_pt, &b_PFPatJet_pt);
   fChain->SetBranchAddress("PFPatJet_eta", &PFPatJet_eta, &b_PFPatJet_eta);
   fChain->SetBranchAddress("PFPatJet_phi", &PFPatJet_phi, &b_PFPatJet_phi);
   fChain->SetBranchAddress("PFPatJet_theta", &PFPatJet_theta, &b_PFPatJet_theta);
   fChain->SetBranchAddress("PFPatJet_px", &PFPatJet_px, &b_PFPatJet_px);
   fChain->SetBranchAddress("PFPatJet_py", &PFPatJet_py, &b_PFPatJet_py);
   fChain->SetBranchAddress("PFPatJet_pz", &PFPatJet_pz, &b_PFPatJet_pz);
   fChain->SetBranchAddress("PFPatJet_vx", &PFPatJet_vx, &b_PFPatJet_vx);
   fChain->SetBranchAddress("PFPatJet_vy", &PFPatJet_vy, &b_PFPatJet_vy);
   fChain->SetBranchAddress("PFPatJet_vz", &PFPatJet_vz, &b_PFPatJet_vz);
   fChain->SetBranchAddress("PFPatJet_ChargedEmEnergy", &PFPatJet_ChargedEmEnergy, &b_PFPatJet_ChargedEmEnergy);
   fChain->SetBranchAddress("PFPatJet_ChargedEmEnergyFrac", &PFPatJet_ChargedEmEnergyFrac, &b_PFPatJet_ChargedEmEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_ChargedHadEnergy", &PFPatJet_ChargedHadEnergy, &b_PFPatJet_ChargedHadEnergy);
   fChain->SetBranchAddress("PFPatJet_ChargedHadEnergyFrac", &PFPatJet_ChargedHadEnergyFrac, &b_PFPatJet_ChargedHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_ChargedHadMult", &PFPatJet_ChargedHadMult, &b_PFPatJet_ChargedHadMult);
   fChain->SetBranchAddress("PFPatJet_ChargedMult", &PFPatJet_ChargedMult, &b_PFPatJet_ChargedMult);
   fChain->SetBranchAddress("PFPatJet_NConstituents", &PFPatJet_NConstituents, &b_PFPatJet_NConstituents);
   fChain->SetBranchAddress("PFPatJet_HFEMEnergy", &PFPatJet_HFEMEnergy, &b_PFPatJet_HFEMEnergy);
   fChain->SetBranchAddress("PFPatJet_HFEMEnergyFrac", &PFPatJet_HFEMEnergyFrac, &b_PFPatJet_HFEMEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_HFEMMult", &PFPatJet_HFEMMult, &b_PFPatJet_HFEMMult);
   fChain->SetBranchAddress("PFPatJet_HFHadEnergy", &PFPatJet_HFHadEnergy, &b_PFPatJet_HFHadEnergy);
   fChain->SetBranchAddress("PFPatJet_HFHadEnergyFrac", &PFPatJet_HFHadEnergyFrac, &b_PFPatJet_HFHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_HFHadMult", &PFPatJet_HFHadMult, &b_PFPatJet_HFHadMult);
   fChain->SetBranchAddress("PFPatJet_NeutralEmEnergy", &PFPatJet_NeutralEmEnergy, &b_PFPatJet_NeutralEmEnergy);
   fChain->SetBranchAddress("PFPatJet_NeutralEmEnergyFrac", &PFPatJet_NeutralEmEnergyFrac, &b_PFPatJet_NeutralEmEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_NeutralHadEnergy", &PFPatJet_NeutralHadEnergy, &b_PFPatJet_NeutralHadEnergy);
   fChain->SetBranchAddress("PFPatJet_NeutralHadEnergyFrac", &PFPatJet_NeutralHadEnergyFrac, &b_PFPatJet_NeutralHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_NeutralHadMult", &PFPatJet_NeutralHadMult, &b_PFPatJet_NeutralHadMult);
   fChain->SetBranchAddress("PFPatJet_NeutralMult", &PFPatJet_NeutralMult, &b_PFPatJet_NeutralMult);
   fChain->SetBranchAddress("PFPatJet_jecUncertainity", &PFPatJet_jecUncertainity, &b_PFPatJet_jecUncertainity);
   fChain->SetBranchAddress("PFPatJet_puJetIdCutBased_MVA", &PFPatJet_puJetIdCutBased_MVA, &b_PFPatJet_puJetIdCutBased_MVA);
   fChain->SetBranchAddress("PFPatJet_puJetIdSimple_MVA", &PFPatJet_puJetIdSimple_MVA, &b_PFPatJet_puJetIdSimple_MVA);
   fChain->SetBranchAddress("PFPatJet_puJetIdFull_MVA", &PFPatJet_puJetIdFull_MVA, &b_PFPatJet_puJetIdFull_MVA);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_loose", &PFPatJet_PassPUJetIdCutBased_loose, &b_PFPatJet_PassPUJetIdCutBased_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_medium", &PFPatJet_PassPUJetIdCutBased_medium, &b_PFPatJet_PassPUJetIdCutBased_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_tight", &PFPatJet_PassPUJetIdCutBased_tight, &b_PFPatJet_PassPUJetIdCutBased_tight);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_loose", &PFPatJet_PassPUJetIdSimple_loose, &b_PFPatJet_PassPUJetIdSimple_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_medium", &PFPatJet_PassPUJetIdSimple_medium, &b_PFPatJet_PassPUJetIdSimple_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_tight", &PFPatJet_PassPUJetIdSimple_tight, &b_PFPatJet_PassPUJetIdSimple_tight);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_loose", &PFPatJet_PassPUJetIdFull_loose, &b_PFPatJet_PassPUJetIdFull_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_medium", &PFPatJet_PassPUJetIdFull_medium, &b_PFPatJet_PassPUJetIdFull_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_tight", &PFPatJet_PassPUJetIdFull_tight, &b_PFPatJet_PassPUJetIdFull_tight);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByTrackCountingHighEff", &PFPatJet_BJetDiscrByTrackCountingHighEff, &b_PFPatJet_BJetDiscrByTrackCountingHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByTrackCountingHighPur", &PFPatJet_BJetDiscrByTrackCountingHighPur, &b_PFPatJet_BJetDiscrByTrackCountingHighPur);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff, &b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur, &b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff, &b_PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA, &b_PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByjetProbabilityBJetTags", &PFPatJet_BJetDiscrByjetProbabilityBJetTags, &b_PFPatJet_BJetDiscrByjetProbabilityBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByjetBProbabilityBJetTags", &PFPatJet_BJetDiscrByjetBProbabilityBJetTags, &b_PFPatJet_BJetDiscrByjetBProbabilityBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySoftElectronBJetTags", &PFPatJet_BJetDiscrBySoftElectronBJetTags, &b_PFPatJet_BJetDiscrBySoftElectronBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySoftMuonBJetTags", &PFPatJet_BJetDiscrBySoftMuonBJetTags, &b_PFPatJet_BJetDiscrBySoftMuonBJetTags);
   fChain->SetBranchAddress("PFRecoJet_n", &PFRecoJet_n, &b_PFRecoJet_n);
   fChain->SetBranchAddress("PFRecoJet_E", &PFRecoJet_E, &b_PFRecoJet_E);
   fChain->SetBranchAddress("PFRecoJet_et", &PFRecoJet_et, &b_PFRecoJet_et);
   fChain->SetBranchAddress("PFRecoJet_pt", &PFRecoJet_pt, &b_PFRecoJet_pt);
   fChain->SetBranchAddress("PFRecoJet_eta", &PFRecoJet_eta, &b_PFRecoJet_eta);
   fChain->SetBranchAddress("PFRecoJet_phi", &PFRecoJet_phi, &b_PFRecoJet_phi);
   fChain->SetBranchAddress("PFRecoJet_theta", &PFRecoJet_theta, &b_PFRecoJet_theta);
   fChain->SetBranchAddress("PFRecoJet_px", &PFRecoJet_px, &b_PFRecoJet_px);
   fChain->SetBranchAddress("PFRecoJet_py", &PFRecoJet_py, &b_PFRecoJet_py);
   fChain->SetBranchAddress("PFRecoJet_pz", &PFRecoJet_pz, &b_PFRecoJet_pz);
   fChain->SetBranchAddress("PFRecoJet_vx", &PFRecoJet_vx, &b_PFRecoJet_vx);
   fChain->SetBranchAddress("PFRecoJet_vy", &PFRecoJet_vy, &b_PFRecoJet_vy);
   fChain->SetBranchAddress("PFRecoJet_vz", &PFRecoJet_vz, &b_PFRecoJet_vz);
   fChain->SetBranchAddress("PFRecoJet_ChargedEmEnergy", &PFRecoJet_ChargedEmEnergy, &b_PFRecoJet_ChargedEmEnergy);
   fChain->SetBranchAddress("PFRecoJet_ChargedEmEnergyFrac", &PFRecoJet_ChargedEmEnergyFrac, &b_PFRecoJet_ChargedEmEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadEnergy", &PFRecoJet_ChargedHadEnergy, &b_PFRecoJet_ChargedHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadEnergyFrac", &PFRecoJet_ChargedHadEnergyFrac, &b_PFRecoJet_ChargedHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadMult", &PFRecoJet_ChargedHadMult, &b_PFRecoJet_ChargedHadMult);
   fChain->SetBranchAddress("PFRecoJet_ChargedMult", &PFRecoJet_ChargedMult, &b_PFRecoJet_ChargedMult);
   fChain->SetBranchAddress("PFRecoJet_NConstituents", &PFRecoJet_NConstituents, &b_PFRecoJet_NConstituents);
   fChain->SetBranchAddress("PFRecoJet_HFEMEnergy", &PFRecoJet_HFEMEnergy, &b_PFRecoJet_HFEMEnergy);
   fChain->SetBranchAddress("PFRecoJet_HFEMEnergyFrac", &PFRecoJet_HFEMEnergyFrac, &b_PFRecoJet_HFEMEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_HFEMMult", &PFRecoJet_HFEMMult, &b_PFRecoJet_HFEMMult);
   fChain->SetBranchAddress("PFRecoJet_HFHadEnergy", &PFRecoJet_HFHadEnergy, &b_PFRecoJet_HFHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_HFHadEnergyFrac", &PFRecoJet_HFHadEnergyFrac, &b_PFRecoJet_HFHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_HFHadMult", &PFRecoJet_HFHadMult, &b_PFRecoJet_HFHadMult);
   fChain->SetBranchAddress("PFRecoJet_NeutralEmEnergy", &PFRecoJet_NeutralEmEnergy, &b_PFRecoJet_NeutralEmEnergy);
   fChain->SetBranchAddress("PFRecoJet_NeutralEmEnergyFrac", &PFRecoJet_NeutralEmEnergyFrac, &b_PFRecoJet_NeutralEmEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadEnergy", &PFRecoJet_NeutralHadEnergy, &b_PFRecoJet_NeutralHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadEnergyFrac", &PFRecoJet_NeutralHadEnergyFrac, &b_PFRecoJet_NeutralHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadMult", &PFRecoJet_NeutralHadMult, &b_PFRecoJet_NeutralHadMult);
   fChain->SetBranchAddress("PFRecoJet_NeutralMult", &PFRecoJet_NeutralMult, &b_PFRecoJet_NeutralMult);
   fChain->SetBranchAddress("PFRecoJet_jecUncertainity", &PFRecoJet_jecUncertainity, &b_PFRecoJet_jecUncertainity);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByTrackCountingHighEff", &PFRecoJet_BJetDiscrByTrackCountingHighEff, &b_PFRecoJet_BJetDiscrByTrackCountingHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByTrackCountingHighPur", &PFRecoJet_BJetDiscrByTrackCountingHighPur, &b_PFRecoJet_BJetDiscrByTrackCountingHighPur);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff, &b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur, &b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff, &b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA, &b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByjetProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetProbabilityBJetTags, &b_PFRecoJet_BJetDiscrByjetProbabilityBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByjetBProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetBProbabilityBJetTags, &b_PFRecoJet_BJetDiscrByjetBProbabilityBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySoftElectronBJetTags", &PFRecoJet_BJetDiscrBySoftElectronBJetTags, &b_PFRecoJet_BJetDiscrBySoftElectronBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySoftMuonBJetTags", &PFRecoJet_BJetDiscrBySoftMuonBJetTags, &b_PFRecoJet_BJetDiscrBySoftMuonBJetTags);
   fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex_chi2", &Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_nchi2", &Vertex_nchi2, &b_Vertex_nchi2);
   fChain->SetBranchAddress("Vertex_ndof", &Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_tracksSize", &Vertex_tracksSize, &b_Vertex_tracksSize);
   fChain->SetBranchAddress("Vertex_isFake", &Vertex_isFake, &b_Vertex_isFake);
   fChain->SetBranchAddress("Vertex_isValid", &Vertex_isValid, &b_Vertex_isValid);
   fChain->SetBranchAddress("Vertex_d0", &Vertex_d0, &b_Vertex_d0);
   fChain->SetBranchAddress("Tracks_n", &Tracks_n, &b_Tracks_n);
   fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_px", &Track_px, &b_Track_px);
   fChain->SetBranchAddress("Track_py", &Track_py, &b_Track_py);
   fChain->SetBranchAddress("Track_pz", &Track_pz, &b_Track_pz);
   fChain->SetBranchAddress("Track_vx", &Track_vx, &b_Track_vx);
   fChain->SetBranchAddress("Track_vy", &Track_vy, &b_Track_vy);
   fChain->SetBranchAddress("Track_vz", &Track_vz, &b_Track_vz);
   fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_phi", &Track_phi, &b_Track_phi);
   fChain->SetBranchAddress("Track_theta", &Track_theta, &b_Track_theta);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("IsScrapingEvent_", &IsScrapingEvent_, &b_IsScrapingEvent_);
   fChain->SetBranchAddress("Scraping_FractionOfGoodTracks", &Scraping_FractionOfGoodTracks, &b_Scraping_FractionOfGoodTracks);
   fChain->SetBranchAddress("HLT_Photon_triggers", &HLT_Photon_triggers, &b_HLT_Photon_triggers);
   fChain->SetBranchAddress("HLT_Photon_trig_prescales", &HLT_Photon_trig_prescales, &b_HLT_Photon_trig_prescales);
   fChain->SetBranchAddress("HLT_Photon_ifTriggerPassed", &HLT_Photon_ifTriggerPassed, &b_HLT_Photon_ifTriggerPassed);
   fChain->SetBranchAddress("HLT_Photon_nTriggers", &HLT_Photon_nTriggers, &b_HLT_Photon_nTriggers);
   fChain->SetBranchAddress("HLT_Photon_triggerIndex", &HLT_Photon_triggerIndex, &b_HLT_Photon_triggerIndex);
   fChain->SetBranchAddress("HLT_Photon_nFilters", &HLT_Photon_nFilters, &b_HLT_Photon_nFilters);
   fChain->SetBranchAddress("HLT_Photon_FilterNames", &HLT_Photon_FilterNames, &b_HLT_Photon_FilterNames);
   fChain->SetBranchAddress("HLT_Photon_trigger_FilterStartPosition", &HLT_Photon_trigger_FilterStartPosition, &b_HLT_Photon_trigger_FilterStartPosition);
   fChain->SetBranchAddress("HLT_Photon_trigger_FilterEndPosition", &HLT_Photon_trigger_FilterEndPosition, &b_HLT_Photon_trigger_FilterEndPosition);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_pt", &HLT_Photon_FilterObjects_pt, &b_HLT_Photon_FilterObjects_pt);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_eta", &HLT_Photon_FilterObjects_eta, &b_HLT_Photon_FilterObjects_eta);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_phi", &HLT_Photon_FilterObjects_phi, &b_HLT_Photon_FilterObjects_phi);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("sigma", &sigma, &b_sigma);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("sigma25", &sigma25, &b_sigma25);
   Notify();
}

Bool_t PostAnalyzer_Data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PostAnalyzer_Data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PostAnalyzer_Data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t PostAnalyzer_Data::PassHLT(){
  Bool_t passedHLT = false;
  std::string hlt = "HLT_Photon150";
  int hlt_idx = -1;
  for(unsigned int i = 0; i < HLT_Photon_triggers->size(); i++){
    if((*HLT_Photon_triggers)[i].find(hlt.c_str())!=std::string::npos){
      hlt_idx = i;
      break;
    }
  }
  if((*HLT_Photon_ifTriggerPassed)[hlt_idx] == true && (*HLT_Photon_trig_prescales)[hlt_idx] == 1){
    passedHLT = true;
  }
  return passedHLT;
}

Bool_t PostAnalyzer_Data::PassHLT_Photon75(){
  Bool_t passedHLT = false;
  std::string hlt = "HLT_Photon75_CaloIdVL";
  int hlt_idx = -1;
  for(unsigned int i = 0; i < HLT_Photon_triggers->size(); i++){
    if((*HLT_Photon_triggers)[i].find(hlt.c_str())!=std::string::npos){
      hlt_idx = i;
      break;
    }
  }
  if((*HLT_Photon_ifTriggerPassed)[hlt_idx] == true){
    passedHLT = true;
  }
  return passedHLT;
}

Bool_t PostAnalyzer_Data::PassHLT_Photon90(){
  Bool_t passedHLT = false;
  std::string hlt = "HLT_Photon90_CaloIdVL";
  int hlt_idx = -1;
  for(unsigned int i = 0; i < HLT_Photon_triggers->size(); i++){
    if((*HLT_Photon_triggers)[i].find(hlt.c_str())!=std::string::npos){
      hlt_idx = i;
      break;
    }
  }
  if((*HLT_Photon_ifTriggerPassed)[hlt_idx] == true){
    passedHLT = true;
  }
  return passedHLT;
}

Bool_t PostAnalyzer_Data::NonScrapingEvt(){
  Bool_t scrapingEvt = !IsScrapingEvent_;
  return scrapingEvt;
}

Int_t PostAnalyzer_Data::GoodPrimaryVtx(){
  Int_t GoodVtx = 0;
  for(Int_t i = 0; i < Vertex_n; i++){
    if( (fabs((*Vertex_z)[i])) <= Cut_Vtx_z && 
        (*Vertex_ndof)[i] >= Cut_Vtx_ndof   &&
        !((*Vertex_isFake)[i])              && 
        (fabs((*Vertex_d0)[i])) <= Cut_Vtx_rho ){
      GoodVtx++;
    }
  }
  return GoodVtx;
}

Bool_t PostAnalyzer_Data::ResSpikes(Int_t i){
  Bool_t spikes = false;
  if( (fabs((*Photon_xtal_timing)[i][0])) < 3.0 &&   //(*Photon_xtal_timing)[i][0]) is the time of arrival for ith photon and seed crystal 
      (*Photon_SigmaIEtaIEta)[i] > 0.001        &&
      (*Photon_SigmaIPhiIPhi)[i] > 0.001        &&
      fabs(GetLICTD(i)) < 5.0                   &&   //LICTD is the largest time difference between the seed crystal and the any other crystal
      (*Photon_r9)[i] < 1.0){
    spikes = true;
  }
  return spikes;
}
      
Double_t PostAnalyzer_Data::GetLICTD(Int_t i){

  Double_t SeedTime = -999;
  Double_t SeedE    = -999;
  Int_t CrysIdx     = -1;

  for(Int_t k = 0; k < (*Photon_nCrystals)[i]; k++){
    Float_t CrysE = (*Photon_xtal_energy)[i][k];
    if(CrysE > SeedE){
      SeedE    = CrysE;
      SeedTime = (*Photon_xtal_timing)[i][k];
      CrysIdx  = k;
    }
  }

  Float_t LICTD = 99.0;

  if(fabs(SeedTime) < 3.0){
    LICTD = 0.0;
    Int_t CrysCrys   = -1;
    Int_t CrysThresh = 0;

    for(Int_t k = 0; k < (*Photon_nCrystals)[i]; k++){
      if(CrysIdx == k)continue;
      Float_t CrysE = (*Photon_xtal_energy)[i][k];

      if(CrysE > 1.0){
        CrysThresh++;
        Float_t timeDiff = (*Photon_xtal_timing)[i][CrysIdx] - (*Photon_xtal_timing)[i][k];
        if(fabs(timeDiff) > fabs(LICTD)){
          LICTD    = timeDiff;
          CrysCrys = k;
        }
      }
    }
  }

  return LICTD;
}

Bool_t PostAnalyzer_Data::ResSpikes_NoLICTD(Int_t i){
  Bool_t spikes = false;
  if( (fabs((*Photon_xtal_timing)[i][0])) < 3.0 &&   //(*Photon_xtal_timing)[i][0]) is the time of arrival for ith photon and seed crystal 
      (*Photon_SigmaIEtaIEta)[i] > 0.001        &&
      (*Photon_SigmaIPhiIPhi)[i] > 0.001        &&
      (*Photon_r9)[i] < 1.0){
    spikes = true;
  }
  return spikes;
}



Double_t PostAnalyzer_Data::EAChargedHadrons(Double_t eta){
  Double_t EffArea = 0;
  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.012;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.010;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.014;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.012;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.016;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.020;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.012;

  return EffArea;
}

Double_t PostAnalyzer_Data::EANeutralHadrons(Double_t eta){
  Double_t EffArea = 0;
  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.030;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.057;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.039;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.015;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.024;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.039;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.072;

  return EffArea;
}

Double_t PostAnalyzer_Data::EAPhotons(Double_t eta){
  Double_t EffArea = 0;
  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.148;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.130;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.112;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.216;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.262;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.260;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.266;

  return EffArea;
}

Bool_t PostAnalyzer_Data::TightPhotonIdnPFIso(Int_t i){
  Bool_t ID = false;
  if((fabs((*Photon_SC_eta)[i])) <= 1.4442){ //Barrel
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.011       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 0.7                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 0.4 + 0.04*((*Photon_pt)[i]) &&
      (TMath::Max(((*PFIsoPhoton03)[i] - rho*EAPhotons((*Photon_SC_eta)[i])), 0.0)) < 0.5 + 0.005*((*Photon_pt)[i]);
  }

  if((fabs((*Photon_SC_eta)[i])) > 1.4442){ //Endcap
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.031       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 0.5                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 1.5 + 0.04*((*Photon_pt)[i]) &&
      (TMath::Max(((*PFIsoPhoton03)[i] - rho*EAPhotons((*Photon_SC_eta)[i])), 0.0)) < 1.0 + 0.005*((*Photon_pt)[i]);
  }

  return ID;
}


Bool_t PostAnalyzer_Data::MediumPhotonIdnPFIso(Int_t i){
  Bool_t ID = false;
  if((fabs((*Photon_SC_eta)[i])) <= 1.4442){ //Barrel
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.011       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 1.5                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 1.0 + 0.04*((*Photon_pt)[i]) &&
      (TMath::Max(((*PFIsoPhoton03)[i] - rho*EAPhotons((*Photon_SC_eta)[i])), 0.0)) < 0.7 + 0.005*((*Photon_pt)[i]);
  }

  if((fabs((*Photon_SC_eta)[i])) > 1.4442){ //Endcap
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.033       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 1.2                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 1.5 + 0.04*((*Photon_pt)[i]) &&
      (TMath::Max(((*PFIsoPhoton03)[i] - rho*EAPhotons((*Photon_SC_eta)[i])), 0.0)) < 1.0 + 0.005*((*Photon_pt)[i]);
  }

  return ID;
}

Bool_t PostAnalyzer_Data::LoosePhotonIdnPFIso(Int_t i){
  Bool_t ID = false;
  if((fabs((*Photon_SC_eta)[i])) <= 1.4442){ //Barrel
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.012       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 2.6                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 3.5 + 0.04*((*Photon_pt)[i]) &&
      (TMath::Max(((*PFIsoPhoton03)[i] - rho*EAPhotons((*Photon_SC_eta)[i])), 0.0)) < 1.3 + 0.005*((*Photon_pt)[i]);
  }

  if((fabs((*Photon_SC_eta)[i])) > 1.4442){ //Endcap
    ID = (*passedConvSafeElectronVeto)[i] == 1 &&
      (*Photon_SingleTowerHoE)[i] < 0.05       &&
      (*Photon_SigmaIEtaIEta)[i] < 0.034       &&
      (TMath::Max(((*PFIsoCharged03)[i] - rho*EAChargedHadrons((*Photon_SC_eta)[i])), 0.0)) < 2.3                          &&
      (TMath::Max(((*PFIsoNeutral03)[i] - rho*EANeutralHadrons((*Photon_SC_eta)[i])), 0.0)) < 2.9 + 0.04*((*Photon_pt)[i]);
  }

  return ID;
}


Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithTightID(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes(i);
    Bool_t ID = TightPhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}

Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithMediumID(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes(i);
    Bool_t ID = MediumPhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}

Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithLooseID(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes(i);
    Bool_t ID = LoosePhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}


Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithTightID_NoLICTD(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes_NoLICTD(i);
    Bool_t ID = TightPhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}

Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithMediumID_NoLICTD(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes_NoLICTD(i);
    Bool_t ID = MediumPhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}

Int_t PostAnalyzer_Data::GetPhotonPassingAllCuts_WithLooseID_NoLICTD(){
  Int_t pc = -1;
  for(Int_t i = 0; i < Photon_n; i++){
    Bool_t resSpikes = ResSpikes_NoLICTD(i);
    Bool_t ID = LoosePhotonIdnPFIso(i);
    if(resSpikes && ID){
      pc = i;
      break;
    }
  }
  return pc;
}


Int_t PostAnalyzer_Data::GetJetPassingIDnMatchedToPhoton(Int_t pc){
  Int_t jc = -1;
  for(Int_t i = 0; i < PFPatJet_n; i++){
    Double_t dR = -1.0;
    if(pc > -1) dR = GetdR((*Photon_SC_eta)[pc], (*PFPatJet_eta)[i], (*Photon_phi)[pc], (*PFPatJet_phi)[i]);

    Bool_t ID = TightJetId(i);
    if(dR > 0.5 && ID){
      jc = i;
      break;
    }
  }
  return jc;
}

Bool_t PostAnalyzer_Data::TightJetId(Int_t i){
  Bool_t ID = false;
  if(fabs((*PFPatJet_eta)[i]) <= 2.4){
    ID = (*PFPatJet_NeutralHadEnergyFrac)[i] < 0.90 &&
      (*PFPatJet_NeutralEmEnergyFrac)[i] < 0.90     &&
      (*PFPatJet_NConstituents)[i] > 1              &&
      (*PFPatJet_ChargedHadEnergyFrac)[i] > 0       &&
      (*PFPatJet_ChargedMult)[i] > 0                &&
      (*PFPatJet_ChargedEmEnergyFrac)[i] < 0.99;
  }

  if(fabs((*PFPatJet_eta)[i]) > 2.4){
    ID = (*PFPatJet_NeutralHadEnergyFrac)[i] < 0.90 &&
      (*PFPatJet_NeutralEmEnergyFrac)[i] < 0.90     &&
      (*PFPatJet_NConstituents)[i] > 1;
  }

  return ID;
}

Double_t PostAnalyzer_Data::GetdEta(Double_t eta1, Double_t eta2){

  Double_t dEta = fabs(eta1 - eta2);
  return dEta;
}

Double_t PostAnalyzer_Data::GetdPhi(Double_t phi1, Double_t phi2){

  Double_t dphi = (phi1 - phi2);
  Double_t twoPi = 2.0*(TMath::Pi());

  if(dphi < 0) dphi = - dphi;
  if(dphi >= (twoPi - dphi)) dphi = twoPi - dphi;

  return dphi;
}

Double_t PostAnalyzer_Data::GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

  Double_t dEta = GetdEta(eta1, eta2);
  Double_t dPhi = GetdPhi(phi1, phi2);

  Double_t dR = 0.0;
  dR = sqrt(dEta*dEta + dPhi*dPhi);

  return dR;
}

Double_t PostAnalyzer_Data::GetInvtMass(Int_t ph, Int_t jet){

  Double_t mass = 0.0;

  Double_t E  = (*Photon_E)[ph] + (*PFPatJet_E)[jet];
  Double_t PX = (*Photon_px)[ph] + (*PFPatJet_px)[jet];
  Double_t PY = (*Photon_py)[ph] + (*PFPatJet_py)[jet];
  Double_t PZ = (*Photon_pz)[ph] + (*PFPatJet_pz)[jet];

  mass = sqrt(E*E - PX*PX - PY*PY - PZ*PZ);

  return mass;
}

Bool_t PostAnalyzer_Data::PassCSVLBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.244){
    passBtag = true;
  }
  return passBtag;
}

Bool_t PostAnalyzer_Data::PassCSVMBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.679){
    passBtag = true;
  }
  return passBtag;
}

Bool_t PostAnalyzer_Data::PassCSVTBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.898){
    passBtag = true;
  }
  return passBtag;
}


void PostAnalyzer_Data::BookHistograms(){
  file->cd();
  //Variable mass bins
  int nbins = 95;
  double xbins[96] = {0, 29, 59, 90, 121, 153, 186, 220, 255, 291, 328, 366, 406, 447, 489, 532, 576, 622, 669, 718, 768, 819, 872, 926, 982, 1040, 1099, 1160, 1223, 1288, 1355, 1424, 1495, 1568, 1643, 1720, 1799, 1881, 1965, 2052, 2141, 2233, 2328, 2425, 2525, 2628, 2734, 2843, 2955, 3071, 3190, 3313, 3439, 3569, 3703, 3841, 3983, 4129, 4279, 4434, 4593, 4757, 4926, 5100, 5279, 5463, 5653, 5848, 6049, 6256, 6469, 6688, 6914, 7146, 7385, 7631, 7884, 8145, 8413, 8689, 8973, 9266, 9567, 9877, 10196, 10525, 10863, 11211, 11570, 11939, 12319, 12710, 13113, 13527, 13954, 14393 };

  h_PhotonPt = new TH1F("h_PhotonPt", "Pt Distribution of Photons", 100, 20.0, 2520.0);
  h_PhotonPt->GetYaxis()->SetTitle("Events/25 GeV");         h_PhotonPt->GetYaxis()->CenterTitle();
  h_PhotonPt->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PhotonPt->GetXaxis()->CenterTitle();
  h_PhotonPt->Sumw2();

  h_PhotonEta = new TH1F("h_PhotonEta", "Eta Distribution of Photons", 100, -2.5, 2.5);
  h_PhotonEta->GetYaxis()->SetTitle("Events");               h_PhotonEta->GetYaxis()->CenterTitle();
  h_PhotonEta->GetXaxis()->SetTitle("{#eta}^{#gamma}");      h_PhotonEta->GetXaxis()->CenterTitle();
  h_PhotonEta->Sumw2();

  h_PhotonPhi = new TH1F("h_PhotonPhi", "Phi Distribution of Photons", 200, 0.0, 10.0);
  h_PhotonPhi->GetYaxis()->SetTitle("Events");               h_PhotonPhi->GetYaxis()->CenterTitle();
  h_PhotonPhi->GetXaxis()->SetTitle("{#phi}^{#gamma}");      h_PhotonPhi->GetXaxis()->CenterTitle();
  h_PhotonPhi->Sumw2();

  h_PhotonSigmaIEtaIEta = new TH1F("h_PhotonSigmaIEtaIEta", "SigmaIEtaIEta Distribution of Photons", 1000, 0.0, 0.05);
  h_PhotonSigmaIEtaIEta->GetYaxis()->SetTitle("Events");                h_PhotonSigmaIEtaIEta->GetYaxis()->CenterTitle();
  h_PhotonSigmaIEtaIEta->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");  h_PhotonSigmaIEtaIEta->GetXaxis()->CenterTitle();
  h_PhotonSigmaIEtaIEta->Sumw2();

  h_PhotonSigmaIPhiIPhi = new TH1F("h_PhotonSigmaIPhiIPhi", "SigmaIPhiIPhi Distribution of Photons", 1000,0.0,0.05);
  h_PhotonSigmaIPhiIPhi->GetYaxis()->SetTitle("Events");                h_PhotonSigmaIPhiIPhi->GetYaxis()->CenterTitle();
  h_PhotonSigmaIPhiIPhi->GetXaxis()->SetTitle("#sigma_{i#phi i#phi}");  h_PhotonSigmaIPhiIPhi->GetXaxis()->CenterTitle();
  h_PhotonSigmaIPhiIPhi->Sumw2();

  h_Photon_r9 = new TH1F("h_Photon_r9", "R9 Distribution of Photons", 100, 0.0, 10.0);
  h_Photon_r9->GetYaxis()->SetTitle("Events");         h_Photon_r9->GetYaxis()->CenterTitle();
  h_Photon_r9->GetXaxis()->SetTitle("Photon_r9");      h_Photon_r9->GetXaxis()->CenterTitle();
  h_Photon_r9->Sumw2();

  h_Photon_SingleTowerHoE = new TH1F("h_Photon_SingleTowerHoE", "H/E Distribution of Photons", 100, 0.0, 0.1);
  h_Photon_SingleTowerHoE->GetYaxis()->SetTitle("Events");      h_Photon_SingleTowerHoE->GetYaxis()->CenterTitle();
  h_Photon_SingleTowerHoE->GetXaxis()->SetTitle("Photon_HoE");  h_Photon_SingleTowerHoE->GetXaxis()->CenterTitle();
  h_Photon_SingleTowerHoE->Sumw2();

  h_Photon_PFIsoCharged03 = new TH1F("h_Photon_PFIsoCharged03", "PF Charged Isolation of Photons", 1500, 0.0, 15.0);
  h_Photon_PFIsoCharged03->GetYaxis()->SetTitle("Events");                h_Photon_PFIsoCharged03->GetYaxis()->CenterTitle();
  h_Photon_PFIsoCharged03->GetXaxis()->SetTitle("Photon_PFChargedIso");   h_Photon_PFIsoCharged03->GetXaxis()->CenterTitle();
  h_Photon_PFIsoCharged03->Sumw2();

  h_Photon_PFIsoNeutral03 = new TH1F("h_Photon_PFIsoNeutral03", "PF Neutral Isolation of Photons", 1500, 0.0, 15.0);
  h_Photon_PFIsoNeutral03->GetYaxis()->SetTitle("Events");                h_Photon_PFIsoNeutral03->GetYaxis()->CenterTitle();
  h_Photon_PFIsoNeutral03->GetXaxis()->SetTitle("Photon_PFNeutralIso");   h_Photon_PFIsoNeutral03->GetXaxis()->CenterTitle();
  h_Photon_PFIsoNeutral03->Sumw2();

  h_Photon_PFIsoPhoton03 = new TH1F("h_Photon_PFIsoPhoton03", "PF Photon Isolation of Photons", 1500, 0.0, 15.0);
  h_Photon_PFIsoPhoton03->GetYaxis()->SetTitle("Events");                 h_Photon_PFIsoPhoton03->GetYaxis()->CenterTitle();
  h_Photon_PFIsoPhoton03->GetXaxis()->SetTitle("Photon_PFPhotonIso");     h_Photon_PFIsoPhoton03->GetXaxis()->CenterTitle();
  h_Photon_PFIsoPhoton03->Sumw2();

  h_Photon_PFIsoSum03 = new TH1F("h_Photon_PFIsoSum03", "PF Isolation Sum of Photons", 1500, 0.0, 15.0);
  h_Photon_PFIsoSum03->GetYaxis()->SetTitle("Events");                    h_Photon_PFIsoSum03->GetYaxis()->CenterTitle();
  h_Photon_PFIsoSum03->GetXaxis()->SetTitle("Photon_PFIsoSum");           h_Photon_PFIsoSum03->GetXaxis()->CenterTitle();
  h_Photon_PFIsoSum03->Sumw2();

  h_JetPt = new TH1F("h_JetPt", "Pt Distribution of Jets", 100, 20.0, 2520.0);
  h_JetPt->GetYaxis()->SetTitle("Events/25 GeV");            h_JetPt->GetYaxis()->CenterTitle();
  h_JetPt->GetXaxis()->SetTitle("P_{T}^{Jet}");              h_JetPt->GetXaxis()->CenterTitle();
  h_JetPt->Sumw2();

  h_JetEta = new TH1F("h_JetEta", "Eta Distribution of Jets", 200, -5.0, 5.0);
  h_JetEta->GetYaxis()->SetTitle("Events");                  h_JetEta->GetYaxis()->CenterTitle();
  h_JetEta->GetXaxis()->SetTitle("{#eta}^{Jet}");            h_JetEta->GetXaxis()->CenterTitle();
  h_JetEta->Sumw2();

  h_JetPhi = new TH1F("h_JetPhi", "Phi Distribution of Jets", 200, 0.0, 10.0);
  h_JetPhi->GetYaxis()->SetTitle("Events");                  h_JetPhi->GetYaxis()->CenterTitle();
  h_JetPhi->GetXaxis()->SetTitle("{#phi}^{Jet}");            h_JetPhi->GetXaxis()->CenterTitle();
  h_JetPhi->Sumw2();

  h_Jet_NeutralHadEnergyFrac = new TH1F("h_Jet_NeutralHadEnergyFrac", "Neutral Hadron Energy Fraction of Jets", 100, 0.0, 1.0);
  h_Jet_NeutralHadEnergyFrac->GetYaxis()->SetTitle("Events");                     h_Jet_NeutralHadEnergyFrac->GetYaxis()->CenterTitle();
  h_Jet_NeutralHadEnergyFrac->GetXaxis()->SetTitle("Jet_NeutralHadEnergyFrac");   h_Jet_NeutralHadEnergyFrac->GetXaxis()->CenterTitle();
  h_Jet_NeutralHadEnergyFrac->Sumw2();

  h_Jet_NeutralEmEnergyFrac = new TH1F("h_Jet_NeutralEmEnergyFrac", "Neutral EM Energy Fraction of Jets", 100, 0.0, 1.0);
  h_Jet_NeutralEmEnergyFrac->GetYaxis()->SetTitle("Events");                      h_Jet_NeutralEmEnergyFrac->GetYaxis()->CenterTitle();
  h_Jet_NeutralEmEnergyFrac->GetXaxis()->SetTitle("Jet_NeutralEmEnergyFrac");     h_Jet_NeutralEmEnergyFrac->GetXaxis()->CenterTitle();
  h_Jet_NeutralEmEnergyFrac->Sumw2();

  h_Jet_NConstituents = new TH1F("h_Jet_NConstituents", "No. of Constituents of Jets", 100, 0.0, 100.0);
  h_Jet_NConstituents->GetYaxis()->SetTitle("Events");               h_Jet_NConstituents->GetYaxis()->CenterTitle();
  h_Jet_NConstituents->GetXaxis()->SetTitle("Jet_NConstituents");    h_Jet_NConstituents->GetXaxis()->CenterTitle();
  h_Jet_NConstituents->Sumw2();

  h_Jet_ChargedHadEnergyFrac = new TH1F("h_Jet_ChargedHadEnergyFrac", "Charged Hadron Energy Fraction of Jets", 100, 0.0, 1.0);
  h_Jet_ChargedHadEnergyFrac->GetYaxis()->SetTitle("Events");                    h_Jet_ChargedHadEnergyFrac->GetYaxis()->CenterTitle();
  h_Jet_ChargedHadEnergyFrac->GetXaxis()->SetTitle("Jet_ChargedHadEnergyFrac");  h_Jet_ChargedHadEnergyFrac->GetXaxis()->CenterTitle();
  h_Jet_ChargedHadEnergyFrac->Sumw2();

  h_Jet_ChargedMult = new TH1F("h_Jet_ChargedMult", "Charged Multiplicity of Jets", 100, 0.0, 100.0);
  h_Jet_ChargedMult->GetYaxis()->SetTitle("Events");             h_Jet_ChargedMult->GetYaxis()->CenterTitle();
  h_Jet_ChargedMult->GetXaxis()->SetTitle("Jet_ChargedMult");    h_Jet_ChargedMult->GetXaxis()->CenterTitle();
  h_Jet_ChargedMult->Sumw2();

  h_Jet_ChargedEmEnergyFrac = new TH1F("h_Jet_ChargedEmEnergyFrac", "Charged EM Energy Fraction of Jets", 100, 0.0, 1.0);
  h_Jet_ChargedEmEnergyFrac->GetYaxis()->SetTitle("Events");                    h_Jet_ChargedEmEnergyFrac->GetYaxis()->CenterTitle();
  h_Jet_ChargedEmEnergyFrac->GetXaxis()->SetTitle("Jet_ChargedEmEnergyFrac");   h_Jet_ChargedEmEnergyFrac->GetXaxis()->CenterTitle();
  h_Jet_ChargedEmEnergyFrac->Sumw2();

  h_GJetInvtMass_binning40 = new TH1F("h_GJetInvtMass_binning40", "Invarient Mass Distribution of Gamma+Jet", 125, 0.0, 5000.0);
  h_GJetInvtMass_binning40->GetYaxis()->SetTitle("Events/40 GeV");     h_GJetInvtMass_binning40->GetYaxis()->CenterTitle();
  h_GJetInvtMass_binning40->GetXaxis()->SetTitle("M({#gamma} + Jet)"); h_GJetInvtMass_binning40->GetXaxis()->CenterTitle();
  h_GJetInvtMass_binning40->Sumw2();

  h_GJetInvtMass_VariableBinning = new TH1F("h_GJetInvtMass_VariableBinning", "Invarient Mass Distribution of Gamma+Jet", nbins, xbins);
  h_GJetInvtMass_VariableBinning->GetYaxis()->SetTitle("Events/Variable");     h_GJetInvtMass_VariableBinning->GetYaxis()->CenterTitle();
  h_GJetInvtMass_VariableBinning->GetXaxis()->SetTitle("M({#gamma} + Jet)"); h_GJetInvtMass_VariableBinning->GetXaxis()->CenterTitle();
  h_GJetInvtMass_VariableBinning->Sumw2();

  h_GJetInvtMass_UnitBinning = new TH1F("h_GJetInvtMass_UnitBinning", "Invarient Mass Distribution of Gamma+Jet", 14000, 0, 14000);
  h_GJetInvtMass_UnitBinning->GetYaxis()->SetTitle("Events/UnitBin");      h_GJetInvtMass_UnitBinning->GetYaxis()->CenterTitle();
  h_GJetInvtMass_UnitBinning->GetXaxis()->SetTitle("M({#gamma} + Jet)");   h_GJetInvtMass_UnitBinning->GetXaxis()->CenterTitle();
  h_GJetInvtMass_UnitBinning->Sumw2(); 

  h_GJetdEta = new TH1F("h_GJetdEta", "DeltaEta Distribution Between Gamma & Jet", 120, 0.0, 6.0);
  h_GJetdEta->GetYaxis()->SetTitle("Events");                h_GJetdEta->GetYaxis()->CenterTitle();
  h_GJetdEta->GetXaxis()->SetTitle("#Delta #eta");           h_GJetdEta->GetXaxis()->CenterTitle();
  h_GJetdEta->Sumw2();

  h_GJetdPhi = new TH1F("h_GJetdPhi", "DeltaPhi Distribution Between Gamma & Jet", 64, 0.0, 3.2);
  h_GJetdPhi->GetYaxis()->SetTitle("Events");                h_GJetdPhi->GetYaxis()->CenterTitle();
  h_GJetdPhi->GetXaxis()->SetTitle("#Delta #phi");           h_GJetdPhi->GetXaxis()->CenterTitle();
  h_GJetdPhi->Sumw2();

  h_GJetdR = new TH1F("h_GJetdR", "DeltaR Distribution Between Gamma & Jet", 100, 0.0, 10.0);
  h_GJetdR->GetYaxis()->SetTitle("Events");                  h_GJetdR->GetYaxis()->CenterTitle();
  h_GJetdR->GetXaxis()->SetTitle("#Delta R");                h_GJetdR->GetXaxis()->CenterTitle();
  h_GJetdR->Sumw2();

  h_CSVL_BJetPt_0bTag = new TH1F("h_CSVL_BJetPt_0bTag", "Pt Distribution of Jets failing CSVL (0bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_0bTag->Sumw2();

  h_CSVL_BJetEta_0bTag = new TH1F("h_CSVL_BJetEta_0bTag", "Eta Distribution of Jets failing CSVL (0bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_0bTag->Sumw2();

  h_CSVL_BJetPhi_0bTag = new TH1F("h_CSVL_BJetPhi_0bTag", "Phi Distribution of Jets failing CSVL (0bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_0bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_0bTag", "Invarient Mass Distribution of Gamma+NoCSVL-Jet (0bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_0bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVL-Jet (0bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_0bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVL-Jet (0bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVL_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVL_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_0bTag->Sumw2();

  h_CSVL_GBJetdEta_0bTag = new TH1F("h_CSVL_GBJetdEta_0bTag", "DeltaEta Distribution Between Gamma & NoCSVL-Jet (0bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_0bTag->Sumw2();

  h_CSVL_GBJetdPhi_0bTag = new TH1F("h_CSVL_GBJetdPhi_0bTag", "DeltaPhi Distribution Between Gamma & NoCSVL-Jet (0bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_0bTag->Sumw2();

  h_CSVL_GBJetdR_0bTag = new TH1F("h_CSVL_GBJetdR_0bTag", "DeltaR Distribution Between Gamma & NoCSVL-Jet (0bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_0bTag->Sumw2();

  h_CSVL_BJetPt_1bTag = new TH1F("h_CSVL_BJetPt_1bTag", "Pt Distribution of BJets passing CSVL (1bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_1bTag->Sumw2();

  h_CSVL_BJetEta_1bTag = new TH1F("h_CSVL_BJetEta_1bTag", "Eta Distribution of BJets passing CSVL (1bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_1bTag->Sumw2();

  h_CSVL_BJetPhi_1bTag = new TH1F("h_CSVL_BJetPhi_1bTag", "Phi Distribution of BJets passing CSVL (1bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_1bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_1bTag", "Invarient Mass Distribution of Gamma+CSVL-BJet (1bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_1bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVL-BJet (1bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_1bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVL-BJet (1bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVL_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVL_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_1bTag->Sumw2();

  h_CSVL_GBJetdEta_1bTag = new TH1F("h_CSVL_GBJetdEta_1bTag", "DeltaEta Distribution Between Gamma & CSVL-BJet (1bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_1bTag->Sumw2();

  h_CSVL_GBJetdPhi_1bTag = new TH1F("h_CSVL_GBJetdPhi_1bTag", "DeltaPhi Distribution Between Gamma & CSVL-BJet (1bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_1bTag->Sumw2();

  h_CSVL_GBJetdR_1bTag = new TH1F("h_CSVL_GBJetdR_1bTag", "DeltaR Distribution Between Gamma & CSVL-BJet (1bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_1bTag->Sumw2();

  h_CSVM_BJetPt_0bTag = new TH1F("h_CSVM_BJetPt_0bTag", "Pt Distribution of Jets failing CSVM (0bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_0bTag->Sumw2();

  h_CSVM_BJetEta_0bTag = new TH1F("h_CSVM_BJetEta_0bTag", "Eta Distribution of Jets failing CSVM (0bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_0bTag->Sumw2();

  h_CSVM_BJetPhi_0bTag = new TH1F("h_CSVM_BJetPhi_0bTag", "Phi Distribution of Jets failing CSVM (0bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_0bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_0bTag", "Invarient Mass Distribution of Gamma+NoCSVM-Jet (0bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_0bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVM-Jet (0bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVM_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_0bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVM-Jet (0bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVM_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVM_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_0bTag->Sumw2();

  h_CSVM_GBJetdEta_0bTag = new TH1F("h_CSVM_GBJetdEta_0bTag", "DeltaEta Distribution Between Gamma & NoCSVM-Jet (0bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_0bTag->Sumw2();

  h_CSVM_GBJetdPhi_0bTag = new TH1F("h_CSVM_GBJetdPhi_0bTag", "DeltaPhi Distribution Between Gamma & NoCSVM-Jet (0bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_0bTag->Sumw2();

  h_CSVM_GBJetdR_0bTag = new TH1F("h_CSVM_GBJetdR_0bTag", "DeltaR Distribution Between Gamma & NoCSVM-Jet (0bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_0bTag->Sumw2();

  h_CSVM_BJetPt_1bTag = new TH1F("h_CSVM_BJetPt_1bTag", "Pt Distribution of BJets passing CSVM (1bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_1bTag->Sumw2();

  h_CSVM_BJetEta_1bTag = new TH1F("h_CSVM_BJetEta_1bTag", "Eta Distribution of BJets passing CSVM (1bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_1bTag->Sumw2();

  h_CSVM_BJetPhi_1bTag = new TH1F("h_CSVM_BJetPhi_1bTag", "Phi Distribution of BJets passing CSVM (1bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_1bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_1bTag", "Invarient Mass Distribution of Gamma+CSVM-BJet (1bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_1bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVM-BJet (1bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVM_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_1bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVM-BJet (1bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVM_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVM_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_1bTag->Sumw2();

  h_CSVM_GBJetdEta_1bTag = new TH1F("h_CSVM_GBJetdEta_1bTag", "DeltaEta Distribution Between Gamma & CSVM-BJet (1bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_1bTag->Sumw2();

  h_CSVM_GBJetdPhi_1bTag = new TH1F("h_CSVM_GBJetdPhi_1bTag", "DeltaPhi Distribution Between Gamma & CSVM-BJet (1bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_1bTag->Sumw2();

  h_CSVM_GBJetdR_1bTag = new TH1F("h_CSVM_GBJetdR_1bTag", "DeltaR Distribution Between Gamma & CSVM-BJet (1bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_1bTag->Sumw2();

  h_CSVT_BJetPt_0bTag = new TH1F("h_CSVT_BJetPt_0bTag", "Pt Distribution of Jets failing CSVT (0bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_0bTag->Sumw2();

  h_CSVT_BJetEta_0bTag = new TH1F("h_CSVT_BJetEta_0bTag", "Eta Distribution of Jets failing CSVT (0bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_0bTag->Sumw2();

  h_CSVT_BJetPhi_0bTag = new TH1F("h_CSVT_BJetPhi_0bTag", "Phi Distribution of Jets failing CSVT (0bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_0bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_0bTag", "Invarient Mass Distribution of Gamma+NoCSVT-Jet (0bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_0bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVT-Jet (0bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_0bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_0bTag", "Invarient Mass Dist. of Gamma+NoCSVT-Jet (0bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVT_GBJetInvtMass_UnitBinning_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVT_GBJetInvtMass_UnitBinning_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_0bTag->Sumw2();

  h_CSVT_GBJetdEta_0bTag = new TH1F("h_CSVT_GBJetdEta_0bTag", "DeltaEta Distribution Between Gamma & NoCSVT-Jet (0bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_0bTag->Sumw2();

  h_CSVT_GBJetdPhi_0bTag = new TH1F("h_CSVT_GBJetdPhi_0bTag", "DeltaPhi Distribution Between Gamma & NoCSVT-Jet (0bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_0bTag->Sumw2();

  h_CSVT_GBJetdR_0bTag = new TH1F("h_CSVT_GBJetdR_0bTag", "DeltaR Distribution Between Gamma & NoCSVT-Jet (0bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_0bTag->Sumw2();

  h_CSVT_BJetPt_1bTag = new TH1F("h_CSVT_BJetPt_1bTag", "Pt Distribution of BJets passing CSVT (1bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_1bTag->Sumw2();

  h_CSVT_BJetEta_1bTag = new TH1F("h_CSVT_BJetEta_1bTag", "Eta Distribution of BJets passing CSVT (1bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_1bTag->Sumw2();

  h_CSVT_BJetPhi_1bTag = new TH1F("h_CSVT_BJetPhi_1bTag", "Phi Distribution of BJets passing CSVT (1bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_1bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_1bTag", "Invarient Mass Distribution of Gamma+CSVT-BJet (1bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_1bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVT-BJet (1bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_1bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_1bTag", "Invarient Mass Dist. of Gamma+CSVT-BJet (1bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->SetTitle("Events/UnitBin");     h_CSVT_GBJetInvtMass_UnitBinning_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVT_GBJetInvtMass_UnitBinning_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_1bTag->Sumw2();

  h_CSVT_GBJetdEta_1bTag = new TH1F("h_CSVT_GBJetdEta_1bTag", "DeltaEta Distribution Between Gamma & CSVT-BJet (1bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_1bTag->Sumw2();

  h_CSVT_GBJetdPhi_1bTag = new TH1F("h_CSVT_GBJetdPhi_1bTag", "DeltaPhi Distribution Between Gamma & CSVT-BJet (1bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_1bTag->Sumw2();

  h_CSVT_GBJetdR_1bTag = new TH1F("h_CSVT_GBJetdR_1bTag", "DeltaR Distribution Between Gamma & CSVT-BJet (1bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_1bTag->Sumw2();

  h_BJetDiscByCSV = new TH1F("h_BJetDiscByCSV", "Value of CSV BJet Discriminator for leading jet in each event", 100, 0.0, 5.0);
  h_BJetDiscByCSV->GetYaxis()->SetTitle("Events");                 h_BJetDiscByCSV->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV->GetXaxis()->SetTitle("BJetDisc_CSV");           h_BJetDiscByCSV->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV->Sumw2();

  h_BJetDiscByCSV_PassingCSVL = new TH1F("h_BJetDiscByCSV_PassingCSVL", "Value of CSV BJet Disc for leading jet passing CSVL", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVL->GetYaxis()->SetTitle("Events");               h_BJetDiscByCSV_PassingCSVL->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL->GetXaxis()->SetTitle("BJetDisc_CSVL");         h_BJetDiscByCSV_PassingCSVL->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL->Sumw2();

  h_BJetDiscByCSV_PassingCSVM = new TH1F("h_BJetDiscByCSV_PassingCSVM", "Value of CSV BJet Disc for leading jet passing CSVM", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVM->GetYaxis()->SetTitle("Events");               h_BJetDiscByCSV_PassingCSVM->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM->GetXaxis()->SetTitle("BJetDisc_CSVM");        h_BJetDiscByCSV_PassingCSVM->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM->Sumw2();

  h_BJetDiscByCSV_PassingCSVT = new TH1F("h_BJetDiscByCSV_PassingCSVT", "Value of CSV BJet Disc for leading jet passing CSVT", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVT->GetYaxis()->SetTitle("Events");               h_BJetDiscByCSV_PassingCSVT->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT->GetXaxis()->SetTitle("BJetDisc_CSVT");        h_BJetDiscByCSV_PassingCSVT->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT->Sumw2();

  h_PassHLT150_EBPassProbes = new TH1F("h_PassHLT150_EBPassProbes", "Pt Distribution of Photons for events passing HLT_Photon150 in EB", 300, 0.0, 1500.0);
  h_PassHLT150_EBPassProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT150_EBPassProbes->GetYaxis()->CenterTitle();
  h_PassHLT150_EBPassProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT150_EBPassProbes->GetXaxis()->CenterTitle();
  h_PassHLT150_EBPassProbes->Sumw2();

  h_PassHLT75OR90_EBAllProbes = new TH1F("h_PassHLT75OR90_EBAllProbes", "Pt Distribution of Photons for events passing HLT_Photon75 OR HLT_Photon90 in EB", 300, 0.0, 1500.0);
  h_PassHLT75OR90_EBAllProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT75OR90_EBAllProbes->GetYaxis()->CenterTitle();
  h_PassHLT75OR90_EBAllProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT75OR90_EBAllProbes->GetXaxis()->CenterTitle();
  h_PassHLT75OR90_EBAllProbes->Sumw2();

  h_PassHLT150_EEPassProbes = new TH1F("h_PassHLT150_EEPassProbes", "Pt Distribution of Photons for events passing HLT_Photon150 in EE", 300, 0.0, 1500.0);
  h_PassHLT150_EEPassProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT150_EEPassProbes->GetYaxis()->CenterTitle();
  h_PassHLT150_EEPassProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT150_EEPassProbes->GetXaxis()->CenterTitle();
  h_PassHLT150_EEPassProbes->Sumw2();

  h_PassHLT75OR90_EEAllProbes = new TH1F("h_PassHLT75OR90_EEAllProbes", "Pt Distribution of Photons for events passing HLT_Photon75 OR HLT_Photon90 in EE", 300, 0.0, 1500.0);
  h_PassHLT75OR90_EEAllProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT75OR90_EEAllProbes->GetYaxis()->CenterTitle();
  h_PassHLT75OR90_EEAllProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT75OR90_EEAllProbes->GetXaxis()->CenterTitle();
  h_PassHLT75OR90_EEAllProbes->Sumw2();

  h_goodPV = new TH1F("h_goodPV", "No. of Good Primary Vertices", 60, 0, 60);
  h_goodPV->GetYaxis()->SetTitle("Events");        h_goodPV->GetYaxis()->CenterTitle();
  h_goodPV->GetXaxis()->SetTitle("nPV (Data)");    h_goodPV->GetXaxis()->CenterTitle();
  h_goodPV->Sumw2();

  h_nPhotons = new TH1F("h_nPhotons", "Number of Photons in each event", 400, 0, 200);
  h_nPhotons->GetYaxis()->SetTitle("Events");      h_nPhotons->GetYaxis()->CenterTitle();
  h_nPhotons->GetXaxis()->SetTitle("nPhotons");    h_nPhotons->GetXaxis()->CenterTitle();
  h_nPhotons->Sumw2();

  h_nJets = new TH1F("h_nJets", "Number of PFJets in each event", 400, 0, 200);
  h_nJets->GetYaxis()->SetTitle("Events");         h_nJets->GetYaxis()->CenterTitle();
  h_nJets->GetXaxis()->SetTitle("nJets");          h_nJets->GetXaxis()->CenterTitle();
  h_nJets->Sumw2();

  h_nCSVLBJets = new TH1F("h_nCSVLBJets", "Number of CSVL BJets in each event", 400, 0, 200);
  h_nCSVLBJets->GetYaxis()->SetTitle("Events");      h_nCSVLBJets->GetYaxis()->CenterTitle();
  h_nCSVLBJets->GetXaxis()->SetTitle("nCSVLBJets");  h_nCSVLBJets->GetXaxis()->CenterTitle();
  h_nCSVLBJets->Sumw2();
 
  h_nCSVMBJets = new TH1F("h_nCSVMBJets", "Number of CSVM BJets in each event", 400, 0, 200);
  h_nCSVMBJets->GetYaxis()->SetTitle("Events");      h_nCSVMBJets->GetYaxis()->CenterTitle();
  h_nCSVMBJets->GetXaxis()->SetTitle("nCSVMBJets");  h_nCSVMBJets->GetXaxis()->CenterTitle();
  h_nCSVMBJets->Sumw2();
  
  h_nCSVTBJets = new TH1F("h_nCSVTBJets", "Number of CSVT BJets in each event", 400, 0, 200);
  h_nCSVTBJets->GetYaxis()->SetTitle("Events");      h_nCSVTBJets->GetYaxis()->CenterTitle();
  h_nCSVTBJets->GetXaxis()->SetTitle("nCSVTBJets");  h_nCSVTBJets->GetXaxis()->CenterTitle();
  h_nCSVTBJets->Sumw2(); 

  h_CSVL_BJetsFrac = new TH1F("h_CSVL_BJetsFrac", "Fraction of CSVL passing BJets among all Jets", 100, 0.0, 1.0);
  h_CSVL_BJetsFrac->GetYaxis()->SetTitle("Events");                h_CSVL_BJetsFrac->GetYaxis()->CenterTitle();
  h_CSVL_BJetsFrac->GetXaxis()->SetTitle("Fraction of BJets");     h_CSVL_BJetsFrac->GetXaxis()->CenterTitle();
  h_CSVL_BJetsFrac->Sumw2();

  h_CSVM_BJetsFrac = new TH1F("h_CSVM_BJetsFrac", "Fraction of CSVM passing BJets among all Jets", 100, 0.0, 1.0);
  h_CSVM_BJetsFrac->GetYaxis()->SetTitle("Events");                h_CSVM_BJetsFrac->GetYaxis()->CenterTitle();
  h_CSVM_BJetsFrac->GetXaxis()->SetTitle("Fraction of BJets");     h_CSVM_BJetsFrac->GetXaxis()->CenterTitle();
  h_CSVM_BJetsFrac->Sumw2();

  h_CSVT_BJetsFrac = new TH1F("h_CSVT_BJetsFrac", "Fraction of CSVT passing BJets among all Jets", 100, 0.0, 1.0);
  h_CSVT_BJetsFrac->GetYaxis()->SetTitle("Events");                h_CSVT_BJetsFrac->GetYaxis()->CenterTitle();
  h_CSVT_BJetsFrac->GetXaxis()->SetTitle("Fraction of BJets");     h_CSVT_BJetsFrac->GetXaxis()->CenterTitle();
  h_CSVT_BJetsFrac->Sumw2();

  h_PhotonIdxVsPt = new TH2F("h_PhotonIdxVsPt", "Index of Photon in PhotonColl Vs. Pt of Photon", 100, 20.0, 2520.0, 100, 0, 100);
  h_PhotonIdxVsPt->GetYaxis()->SetTitle("Idx Of Photon");            h_PhotonIdxVsPt->GetYaxis()->CenterTitle();
  h_PhotonIdxVsPt->GetXaxis()->SetTitle("P_{T}^{#gamma}");           h_PhotonIdxVsPt->GetXaxis()->CenterTitle();
  h_PhotonIdxVsPt->Sumw2();

  h_JetIdxVsPt = new TH2F("h_JetIdxVsPt", "Index of Jet in JetColl Vs. Pt of Jet", 100, 20.0, 2520.0, 100, 0, 100);
  h_JetIdxVsPt->GetYaxis()->SetTitle("Idx Of Jet");               h_JetIdxVsPt->GetYaxis()->CenterTitle();
  h_JetIdxVsPt->GetXaxis()->SetTitle("P_{T}^{Jet}");              h_JetIdxVsPt->GetXaxis()->CenterTitle();
  h_JetIdxVsPt->Sumw2();

  h_CSVLBJetIdxVsPt = new TH2F("h_CSVLBJetIdxVsPt", "Index of CSVL-BJet in JetColl Vs. Pt of BJet", 100, 20.0, 2520.0, 100, 0, 100);
  h_CSVLBJetIdxVsPt->GetYaxis()->SetTitle("Idx Of CSVL-BJet");              h_CSVLBJetIdxVsPt->GetYaxis()->CenterTitle();
  h_CSVLBJetIdxVsPt->GetXaxis()->SetTitle("P_{T}^{CSVL-BJet}");             h_CSVLBJetIdxVsPt->GetXaxis()->CenterTitle();
  h_CSVLBJetIdxVsPt->Sumw2();

  h_CSVMBJetIdxVsPt = new TH2F("h_CSVMBJetIdxVsPt", "Index of CSVM-BJet in JetColl Vs. Pt of BJet", 100, 20.0, 2520.0, 100, 0, 100);
  h_CSVMBJetIdxVsPt->GetYaxis()->SetTitle("Idx Of CSVM-BJet");              h_CSVMBJetIdxVsPt->GetYaxis()->CenterTitle();
  h_CSVMBJetIdxVsPt->GetXaxis()->SetTitle("P_{T}^{CSVM-BJet}");             h_CSVMBJetIdxVsPt->GetXaxis()->CenterTitle();
  h_CSVMBJetIdxVsPt->Sumw2();

  h_CSVTBJetIdxVsPt = new TH2F("h_CSVTBJetIdxVsPt", "Index of CSVT-BJet in JetColl Vs. Pt of BJet", 100, 20.0, 2520.0, 100, 0, 100);
  h_CSVTBJetIdxVsPt->GetYaxis()->SetTitle("Idx Of CSVT-BJet");              h_CSVTBJetIdxVsPt->GetYaxis()->CenterTitle();
  h_CSVTBJetIdxVsPt->GetXaxis()->SetTitle("P_{T}^{CSVT-BJet}");             h_CSVTBJetIdxVsPt->GetXaxis()->CenterTitle();
  h_CSVTBJetIdxVsPt->Sumw2();

  h_PC = new TH1F("h_PC", "Photon Candidate", 10, 0, 10);
  h_PC->GetYaxis()->SetTitle("Events");                       h_PC->GetYaxis()->CenterTitle();
  h_PC->GetXaxis()->SetTitle("Position of Photon");           h_PC->GetXaxis()->CenterTitle();
  h_PC->Sumw2();

  h_JC = new TH1F("h_JC", "Jet Candidate", 20, 0, 20);
  h_JC->GetYaxis()->SetTitle("Events");                       h_JC->GetYaxis()->CenterTitle();
  h_JC->GetXaxis()->SetTitle("Position of Jet");              h_JC->GetXaxis()->CenterTitle();
  h_JC->Sumw2();

}





#endif // #ifdef PostAnalyzer_Data_cxx



 
EOF


cat > analysis_${filenameTag}.C <<EOF
#include "PostAnalyzer_Data.C"
#include "TROOT.h"
int main(){
    PostAnalyzer_Data a;
    a.Loop();
    return 0;
}

EOF


####Compilation

g++ -Wno-deprecated analysis_${filenameTag}.C -o ${filenameTag}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`

####Execution
#./${filenameTag}.exe

###Submit jobs

chmod 775 MakeCondorFiles.csh

./MakeCondorFiles.csh ${filenameTag}


#if [ ! -d Junk ]; then
#echo "Making Directory Junk"
#mkdir Junk
#chmod 775 -R Junk
#fi
#
#mv analysis_${filenameTag}.C Junk/analysis_${filenameTag}.C
#mv ${filenameTag}.exe Junk/${filenameTag}.exe
#mv Job_${filenameTag}.csh Junk/Job_${filenameTag}.csh
#mv condor_${filenameTag} Junk/condor_${filenameTag}


done ##end of for loop##


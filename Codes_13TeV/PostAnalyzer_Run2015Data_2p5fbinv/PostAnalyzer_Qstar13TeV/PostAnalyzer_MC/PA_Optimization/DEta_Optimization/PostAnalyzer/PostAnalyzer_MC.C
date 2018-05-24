#define PostAnalyzer_MC_cxx
#include "PostAnalyzer_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//To be removed in script
//-----------------------
int main(){
  PostAnalyzer_MC t;
  t.Loop();
  return 0;
}
//------------------------

void PostAnalyzer_MC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PostAnalyzer_MC.C
//      root> PostAnalyzer_MC t
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

   //Luminosity of Data Compared
   //Lumi = 2197.327; // (pb^{-1}), COMPLETE 2015 25ns DATA WITH GOLDEN JSON
   Lumi = 2502.816; // (pb^{-1}), COMPLETE 2015 25ns DATA WITH SILVER JSON

   //To be removed in script
   //-----------------------
   Float_t XS = 30.12207;
   //-----------------------

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

   //To be removed in script
   //-----------------------
   //Define Output file here
   file = new TFile("PostAnalyzer_MC.root", "RECREATE");
   //-----------------------

   //Uncomment this in script
   //Define Output file here
   //   TString OutputPath = "${destinationDir}/";
   //   TString OutputFile = "${filenameTag}";
   //   file = new TFile(OutputPath+OutputFile+".root", "RECREATE");

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   PileupReWeighting();

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   //   cout << "<Total entries" << nentries << endl; 
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Uncomment this in script  
      //Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/${totalEvents[${sampleIndex}]};

      //To be removed in script                                                                       
      //-----------------------                                                                                                                     
      Lumi_EvtWt = Lumi*117276/3364368;//2000069
      //-----------------------                                                                                                 

      PU_EvtWt = PUWeights((*puTrue)[0]);
      PreBTag_EvtWt = Lumi_EvtWt * PU_EvtWt;
      
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
      Bool_t Pass_dEta1p0 = false;
      Bool_t Pass_dEta1p2 = false;
      Bool_t Pass_dEta1p5 = false;
      Bool_t Pass_dEta2p0 = false;
      Bool_t Pass_dEta2p2 = false;

      //Running different functions     
      Pass_HLT = true;
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
      if(PC > -1 && JC > -1){
	Pass_dEta1p0 = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < 1.0);
	Pass_dEta1p2 = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < 1.2);
	Pass_dEta1p5 = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < 1.5);
	Pass_dEta2p0 = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < 2.0);
	Pass_dEta2p2 = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < 2.2);
      }

      h_CutFlow->Fill(0.5);
      h_CutFlowWithWts->Fill(0.5, PreBTag_EvtWt);

      if(Pass_HLT){
	h_CutFlow->Fill(1.5);
	h_CutFlowWithWts->Fill(1.5, PreBTag_EvtWt);

	if(HasPrimaryVtx){
	  h_CutFlow->Fill(2.5);
	  h_CutFlowWithWts->Fill(2.5, PreBTag_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow->Fill(3.5);
	    h_CutFlowWithWts->Fill(3.5, PreBTag_EvtWt);

	    if(PC > -1){
	      h_CutFlow->Fill(4.5);
	      h_CutFlowWithWts->Fill(4.5, PreBTag_EvtWt);

	      if(JC > -1){
		h_CutFlow->Fill(5.5);
		h_CutFlowWithWts->Fill(5.5, PreBTag_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow->Fill(6.5);
		  h_CutFlowWithWts->Fill(6.5, PreBTag_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow->Fill(7.5);
		    h_CutFlowWithWts->Fill(7.5, PreBTag_EvtWt);

		    if(Pass_GJdPhi){
		      h_CutFlow->Fill(8.5);
		      h_CutFlowWithWts->Fill(8.5, PreBTag_EvtWt);

		      //deta = 1.0
		      if(Pass_dEta1p0){
			h_CutFlow->Fill(9.5);
			h_CutFlowWithWts->Fill(9.5, PreBTag_EvtWt);
			
			//Uncomment in script
			//h_mass_X_bin1_deta1p0            ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			h_GJetInvtMass_bin40_deta1p0[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_VarBin_deta1p0[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_UnitBin_deta1p0[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(10.5);
			  h_CutFlowWithWts->Fill(10.5, PreBTag_EvtWt);

			  h_GJetInvtMass_bin40_deta1p0[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_VarBin_deta1p0[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_UnitBin_deta1p0[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			}//Pass_GJInvtMass
		      }//Pass_GJdEta1p0

		      //deta = 1.2
		      if(Pass_dEta1p2){
			h_CutFlow->Fill(11.5);
			h_CutFlowWithWts->Fill(11.5, PreBTag_EvtWt);

			//Uncomment in script
			//h_mass_X_bin1_deta1p2            ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			h_GJetInvtMass_bin40_deta1p2[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_VarBin_deta1p2[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_UnitBin_deta1p2[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(12.5);
			  h_CutFlowWithWts->Fill(12.5, PreBTag_EvtWt);

			  h_GJetInvtMass_bin40_deta1p2[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_VarBin_deta1p2[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_UnitBin_deta1p2[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			}//Pass_GJInvtMass
		      }//Pass_GJdEta1p2

		      //deta = 1.5
		      if(Pass_dEta1p5){
			h_CutFlow->Fill(13.5);
			h_CutFlowWithWts->Fill(13.5, PreBTag_EvtWt);

			//Uncomment in script
			//h_mass_X_bin1_deta1p5            ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			h_GJetInvtMass_bin40_deta1p5[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_VarBin_deta1p5[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_UnitBin_deta1p5[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(14.5);
			  h_CutFlowWithWts->Fill(14.5, PreBTag_EvtWt);

			  h_GJetInvtMass_bin40_deta1p5[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_VarBin_deta1p5[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_UnitBin_deta1p5[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			}//Pass_GJInvtMass
		      }//Pass_GJdEta1p5

		      //deta = 2.0
		      if(Pass_dEta2p0){
			h_CutFlow->Fill(15.5);
			h_CutFlowWithWts->Fill(15.5, PreBTag_EvtWt);

			//Uncomment in script
			//h_mass_X_bin1_deta2p0            ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			h_GJetInvtMass_bin40_deta2p0[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_VarBin_deta2p0[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_UnitBin_deta2p0[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(16.5);
			  h_CutFlowWithWts->Fill(16.5, PreBTag_EvtWt);

			  h_GJetInvtMass_bin40_deta2p0[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_VarBin_deta2p0[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_UnitBin_deta2p0[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			}//Pass_GJInvtMass
		      }//Pass_GJdEta2p0

		      //deta 2.2
		      if(Pass_dEta2p2){
			h_CutFlow->Fill(17.5);
			h_CutFlowWithWts->Fill(17.5, PreBTag_EvtWt);

			//Uncomment in script
			//h_mass_X_bin1_deta2p2            ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			h_GJetInvtMass_bin40_deta2p2[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_VarBin_deta2p2[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			h_GJetInvtMass_UnitBin_deta2p2[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);

			if(Pass_GJInvtMass){
			  h_CutFlow->Fill(18.5);
			  h_CutFlowWithWts->Fill(18.5, PreBTag_EvtWt);

			  h_GJetInvtMass_bin40_deta2p2[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_VarBin_deta2p2[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJetInvtMass_UnitBin_deta2p2[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			
			}//Pass_GJInvtMass       	      
		      }//Pass_GJdEta2p2
		   
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

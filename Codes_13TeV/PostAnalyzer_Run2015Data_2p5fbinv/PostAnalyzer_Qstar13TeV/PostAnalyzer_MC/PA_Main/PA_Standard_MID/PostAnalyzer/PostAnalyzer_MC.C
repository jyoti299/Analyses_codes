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
      PC_Gen = -1;
      JC_Gen = -1;
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
      Pass_HLT = true;
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
      if(PC > -1) PC_Gen = MatchedGenPhotonToReco(PC); 
      if(JC > -1) JC_Gen = MatchedGenJetToReco(JC);

      h_trueNumofInt       ->Fill((*puTrue)[0]);
      h_goodPV             ->Fill(GoodVertex);
      h_goodPV_LumiWt      ->Fill(GoodVertex, Lumi_EvtWt);
      h_goodPV_PUWt        ->Fill(GoodVertex, PU_EvtWt);
      h_goodPV_TotalWt[0]  ->Fill(GoodVertex, PreBTag_EvtWt);
      h_nPhotons[0]        ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);
      h_nJets[0]           ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);

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

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);
			h_CutFlowWithWts->Fill(9.5, PreBTag_EvtWt);

			  h_PC                       ->Fill(PC, PreBTag_EvtWt);
			  h_JC                       ->Fill(JC, PreBTag_EvtWt);

                          h_PhotonPt[0]              ->Fill((*phoEt)[PC], PreBTag_EvtWt);
                          h_PhotonEta[0]             ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
                          h_PhotonPhi[0]             ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			  h_Photon_SigmaIEtaIEta[0]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			  h_Photon_R9[0]             ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			  h_Photon_HoverE[0]         ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			  h_Photon_EleVeto[0]        ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFChIso[0]    ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			  h_Photon_CorrPFNeuIso[0]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			  h_Photon_CorrPFPhoIso[0]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);

                          h_JetPt[0]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
                          h_JetEta[0]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
                          h_JetPhi[0]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			  h_Jet_NHEF[0]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			  h_Jet_NEEF[0]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			  h_Jet_NConst[0]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			  h_Jet_CHEF[0]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			  h_Jet_ChMult[0]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			  h_Jet_CEEF[0]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);

                          h_GJetInvtMass_bin40[0]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
                          h_GJetInvtMass_VarBin[0]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
                          h_GJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			  h_GJet_dEta[0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			  h_GJet_dPhi[0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			  h_GJet_dR[0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);

			  ///Uncomment in script
			  //h_mass_X_bin1              ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, PreBTag_EvtWt);

			  h_PFMet[0]                 ->Fill(pfMET, PreBTag_EvtWt);
			  h_PFMetVsGJmass[0]         ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			  h_PFMetOverSumEtVsGJmass[0]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);

			  h_goodPV_TotalWt[1]        ->Fill(GoodVertex, PreBTag_EvtWt);
			  h_nPhotons[1]              ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);
			  h_nJets[1]                 ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);

			  if(PC_Gen > -1 && JC_Gen > -1){
			    h_Gen_GJetInvtMass_bin40[0]   ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			    h_Gen_GJetInvtMass_VarBin[0]  ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			    h_Gen_GJetInvtMass_UnitBin[0] ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));

			    h_GJetMassResolution[0]    ->Fill(fabs(GetGenLevelInvtMass(PC_Gen, JC_Gen) - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  }


			  if(Pass_GJInvtMass && IfMassCut){
			    h_CutFlow->Fill(10.5);
			    h_CutFlowWithWts->Fill(10.5, PreBTag_EvtWt);

			    h_PhotonPt[1]              ->Fill((*phoEt)[PC], PreBTag_EvtWt);
			    h_PhotonEta[1]             ->Fill((*phoSCEta)[PC], PreBTag_EvtWt);
			    h_PhotonPhi[1]             ->Fill((*phoSCPhi)[PC], PreBTag_EvtWt);
			    h_Photon_SigmaIEtaIEta[1]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], PreBTag_EvtWt);
			    h_Photon_R9[1]             ->Fill((*phoR9)[PC], PreBTag_EvtWt);
			    h_Photon_HoverE[1]         ->Fill((*phoHoverE)[PC], PreBTag_EvtWt);
			    h_Photon_EleVeto[1]        ->Fill((*phoEleVeto)[PC], PreBTag_EvtWt);
			    h_Photon_CorrPFChIso[1]    ->Fill((*phoPFChIso)[PC], PreBTag_EvtWt);
			    h_Photon_CorrPFNeuIso[1]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);
			    h_Photon_CorrPFPhoIso[1]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), PreBTag_EvtWt);

			    h_JetPt[1]                 ->Fill((*jetPt)[JC], PreBTag_EvtWt);
			    h_JetEta[1]                ->Fill((*jetEta)[JC], PreBTag_EvtWt);
			    h_JetPhi[1]                ->Fill((*jetPhi)[JC], PreBTag_EvtWt);
			    h_Jet_NHEF[1]              ->Fill((*jetNHF)[JC], PreBTag_EvtWt);
			    h_Jet_NEEF[1]              ->Fill((*jetNEF)[JC], PreBTag_EvtWt);
			    h_Jet_NConst[1]            ->Fill((*jetNConstituents)[JC], PreBTag_EvtWt);
			    h_Jet_CHEF[1]              ->Fill((*jetCHF)[JC], PreBTag_EvtWt);
			    h_Jet_ChMult[1]            ->Fill((*jetNCH)[JC], PreBTag_EvtWt);
			    h_Jet_CEEF[1]              ->Fill((*jetCEF)[JC], PreBTag_EvtWt);

			    h_GJetInvtMass_bin40[1]    ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GJetInvtMass_VarBin[1]   ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			    h_GJet_dEta[1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), PreBTag_EvtWt);
			    h_GJet_dPhi[1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);
			    h_GJet_dR[1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), PreBTag_EvtWt);

			    h_PFMet[1]                 ->Fill(pfMET, PreBTag_EvtWt);
			    h_PFMetVsGJmass[1]         ->Fill(GetInvtMass(PC, JC), pfMET, PreBTag_EvtWt);
			    h_PFMetOverSumEtVsGJmass[1]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, PreBTag_EvtWt);

			    h_goodPV_TotalWt[2]        ->Fill(GoodVertex, PreBTag_EvtWt);
			    h_nPhotons[2]              ->Fill(GoodIsoPhotons.size(), PreBTag_EvtWt);
			    h_nJets[2]                 ->Fill(GoodIsoJets.size(), PreBTag_EvtWt);

			    if(PC_Gen > -1 && JC_Gen > -1){
			      h_Gen_GJetInvtMass_bin40[1]   ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			      h_Gen_GJetInvtMass_VarBin[1]  ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			      h_Gen_GJetInvtMass_UnitBin[1] ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));

			      h_GJetMassResolution[1]    ->Fill(fabs(GetGenLevelInvtMass(PC_Gen, JC_Gen) - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			    }


			    if(pfMET < 500){
			      h_CutFlow->Fill(11.5);
			      h_CutFlowWithWts->Fill(11.5, PreBTag_EvtWt);

			      if(pfMET < 250){
				h_CutFlow->Fill(12.5);
				h_CutFlowWithWts->Fill(12.5, PreBTag_EvtWt);;

				if(pfMET < 100){
				  h_CutFlow->Fill(13.5);
				  h_CutFlowWithWts->Fill(13.5, PreBTag_EvtWt);

				  if(pfMET < 50){
				    h_CutFlow->Fill(14.5);
				    h_CutFlowWithWts->Fill(14.5, PreBTag_EvtWt);

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
   
      //Photon Efficiency
      PC_G = -1;
      PC_L = -1;
      PC_M = -1;
      PC_T = -1;
      PC_H = -1;

      PC_G = MatchedRecoPhotonToGen_WithGenIsoCut();
      PC_L = FirstGoodPhoton("loose");
      PC_M = FirstGoodPhoton("medium");
      PC_T = FirstGoodPhoton("tight");
      PC_H = FirstHighPtIDPhoton();

      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC_G > -1){
	    if(fabs((*phoSCEta)[PC_G]) <= 1.4442){//Barrel 

	      h_PassGenIsoMatch_EB->Fill((*phoEt)[PC_G], PreBTag_EvtWt);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EB->Fill((*phoEt)[PC_L], PreBTag_EvtWt);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EB->Fill((*phoEt)[PC_M], PreBTag_EvtWt);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EB->Fill((*phoEt)[PC_T], PreBTag_EvtWt);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EB->Fill((*phoEt)[PC_H], PreBTag_EvtWt);
	      }
	    }
	    if(fabs((*phoSCEta)[PC_G]) < 2.5 && fabs((*phoSCEta)[PC_G]) >= 1.5666){//Endcap  

	      h_PassGenIsoMatch_EE->Fill((*phoEt)[PC_G], PreBTag_EvtWt);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EE->Fill((*phoEt)[PC_L], PreBTag_EvtWt);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EE->Fill((*phoEt)[PC_M], PreBTag_EvtWt);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EE->Fill((*phoEt)[PC_T], PreBTag_EvtWt);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EE->Fill((*phoEt)[PC_H], PreBTag_EvtWt);
	      }
	    }	  
	  }
	}
      }

   }//for jentry
}//Loop()

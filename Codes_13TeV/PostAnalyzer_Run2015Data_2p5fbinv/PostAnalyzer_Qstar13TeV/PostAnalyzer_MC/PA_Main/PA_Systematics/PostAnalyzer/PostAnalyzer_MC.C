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
   Lumi = 2670.555; //(pb^{-1}), COMPLETE 2015 25ns DATA WITH SILVER JSON AFTER CHANGING REFERENCE FILE IN BIRL CALC (DONE CENTRALLY BY CMS LUMI GRP)

   //To be removed in script
   //-----------------------
   Float_t XS = 5813.0;
   //-----------------------

   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 190.0; // GeV                                                                                                             
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 190.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.8;

   Cut_GJInvtMass = 560.0;

   //True if masscut is applied 
   IfMassCut = true;

   //Photon energy correction uncertainty (taking a conservative value for now( this value for 8 TeV was 1.5%, so now taking 2%)
   Photon_pecUnc = 0.02;

   //To be removed in script
   //-----------------------
   //Define Output file here
   file = new TFile("PostAnalyzer_MC.root", "RECREATE");
   //-----------------------

   //Uncomment this in script
   //Define Output file here
   //   TString OutputPath = "${destinationDir}/";
   //   TString OutputFile = "${filenameTag}";
   //  // file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   //   file = new TFile(OutputFile+".root", "RECREATE");

   //Uncomment it in script
   IfGJet = true;
   //IfGJet = ${GJ};

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   PileupReWeighting();
   PileupReWeighting_XSm5();
   PileupReWeighting_XSp5();

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << "<Total entries: " << nentries << endl; 
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
      Lumi_EvtWt = (Lumi*XS)/9956130;//2000069
      //-----------------------                                                                                                 

      PU_EvtWt = PUWeights((*puTrue)[0]);
      PreBTag_EvtWt = Lumi_EvtWt * PU_EvtWt;
      Total_EvtWt = PreBTag_EvtWt;

      PU_EvtWt_XSm5 = PUWeights_XSm5((*puTrue)[0]);
      Total_EvtWt_XSm5 = Lumi_EvtWt * PU_EvtWt_XSm5;

      PU_EvtWt_XSp5 = PUWeights_XSp5((*puTrue)[0]);
      Total_EvtWt_XSp5 = Lumi_EvtWt * PU_EvtWt_XSp5;
      
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
      Pass_HLT = true;
      HasPrimaryVtx = hasGoodVtx;
       //HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);

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


      h_CutFlow->Fill(0.5);
      h_CutFlowWithWts->Fill(0.5, Total_EvtWt);
      h_CutFlowWithWts_XSm5->Fill(0.5, Total_EvtWt_XSm5);
      h_CutFlowWithWts_XSp5->Fill(0.5, Total_EvtWt_XSp5);

      if(Pass_HLT){
	h_CutFlow->Fill(1.5);
	h_CutFlowWithWts->Fill(1.5, Total_EvtWt);
	h_CutFlowWithWts_XSm5->Fill(1.5, Total_EvtWt_XSm5);
	h_CutFlowWithWts_XSp5->Fill(1.5, Total_EvtWt_XSp5);

	if(HasPrimaryVtx){
	  h_CutFlow->Fill(2.5);
	  h_CutFlowWithWts->Fill(2.5, Total_EvtWt);
	  h_CutFlowWithWts_XSm5->Fill(2.5, Total_EvtWt_XSm5);
	  h_CutFlowWithWts_XSp5->Fill(2.5, Total_EvtWt_XSp5);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow->Fill(3.5);
	    h_CutFlowWithWts->Fill(3.5, Total_EvtWt);
	    h_CutFlowWithWts_XSm5->Fill(3.5, Total_EvtWt_XSm5);
	    h_CutFlowWithWts_XSp5->Fill(3.5, Total_EvtWt_XSp5);

	    if(PC > -1){
	      h_CutFlow->Fill(4.5);
	      h_CutFlowWithWts->Fill(4.5, Total_EvtWt);
	      h_CutFlowWithWts_XSm5->Fill(4.5, Total_EvtWt_XSm5);
	      h_CutFlowWithWts_XSp5->Fill(4.5, Total_EvtWt_XSp5);

	      if(JC > -1){
		h_CutFlow->Fill(5.5);
		h_CutFlowWithWts->Fill(5.5, Total_EvtWt);
		h_CutFlowWithWts_XSm5->Fill(5.5, Total_EvtWt_XSm5);
		h_CutFlowWithWts_XSp5->Fill(5.5, Total_EvtWt_XSp5);

		if(Pass_JetPt){
		  h_CutFlow->Fill(6.5);
		  h_CutFlowWithWts->Fill(6.5, Total_EvtWt);
		  h_CutFlowWithWts_XSm5->Fill(6.5, Total_EvtWt_XSm5);
		  h_CutFlowWithWts_XSp5->Fill(6.5, Total_EvtWt_XSp5);

		  if(Pass_JetEta){
		    h_CutFlow->Fill(7.5);
		    h_CutFlowWithWts->Fill(7.5, Total_EvtWt);
		    h_CutFlowWithWts_XSm5->Fill(7.5, Total_EvtWt_XSm5);
		    h_CutFlowWithWts_XSp5->Fill(7.5, Total_EvtWt_XSp5);

		    if(Pass_GJdPhi){
		      h_CutFlow->Fill(8.5);
		      h_CutFlowWithWts->Fill(8.5, Total_EvtWt);
		      h_CutFlowWithWts_XSm5->Fill(8.5, Total_EvtWt_XSm5);
		      h_CutFlowWithWts_XSp5->Fill(8.5, Total_EvtWt_XSp5);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);
			h_CutFlowWithWts->Fill(9.5, Total_EvtWt);
			h_CutFlowWithWts_XSm5->Fill(9.5, Total_EvtWt_XSm5);
			h_CutFlowWithWts_XSp5->Fill(9.5, Total_EvtWt_XSp5);

			//JES
			JetPx = (*jetPt)[JC]*TMath::Cos((*jetPhi)[JC]); 
			JetPy = (*jetPt)[JC]*TMath::Sin((*jetPhi)[JC]); 
		        JetPz = (*jetPt)[JC]*TMath::SinH((*jetEta)[JC]);

			phoPx = (*phoEt)[PC]*TMath::Cos((*phoSCPhi)[PC]);
			phoPy = (*phoEt)[PC]*TMath::Sin((*phoSCPhi)[PC]);
			phoPz = (*phoEt)[PC]*TMath::SinH((*phoSCEta)[PC]);

			PFJet4Vec_JESup.SetPxPyPzE(JetPx*(1+(*jetJECUnc)[JC]), JetPy*(1+(*jetJECUnc)[JC]), JetPz*(1+(*jetJECUnc)[JC]), (*jetEn)[JC]*(1+(*jetJECUnc)[JC]));
			PFJet4Vec_JESdown.SetPxPyPzE(JetPx*(1-(*jetJECUnc)[JC]), JetPy*(1-(*jetJECUnc)[JC]), JetPz*(1-(*jetJECUnc)[JC]), (*jetEn)[JC]*(1-(*jetJECUnc)[JC]));

			GJ4Vec_JESup.SetPxPyPzE(phoPx+PFJet4Vec_JESup[0], phoPy+PFJet4Vec_JESup[1], phoPz+PFJet4Vec_JESup[2], (*phoE)[PC]+PFJet4Vec_JESup[3]);
			GJ4Vec_JESdown.SetPxPyPzE(phoPx+PFJet4Vec_JESdown[0], phoPy+PFJet4Vec_JESdown[1], phoPz+PFJet4Vec_JESdown[2], (*phoE)[PC]+PFJet4Vec_JESdown[3]);

			Mass_JESup = pow((GJ4Vec_JESup[3]*GJ4Vec_JESup[3] - GJ4Vec_JESup[0]*GJ4Vec_JESup[0] - GJ4Vec_JESup[1]*GJ4Vec_JESup[1] - GJ4Vec_JESup[2]*GJ4Vec_JESup[2]), 0.5);
			Mass_JESdown = pow((GJ4Vec_JESdown[3]*GJ4Vec_JESdown[3] - GJ4Vec_JESdown[0]*GJ4Vec_JESdown[0] - GJ4Vec_JESdown[1]*GJ4Vec_JESdown[1] - GJ4Vec_JESdown[2]*GJ4Vec_JESdown[2]), 0.5);

			h_profile_JES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_profile_JES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			h_mass_JES[0][0]->Fill(Mass_JESup, Total_EvtWt);
			h_mass_JES[1][0]->Fill(Mass_JESdown, Total_EvtWt);

			h_mass_X_bin1->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, Total_EvtWt);

			h_mass_X_bin1_JES[0]->Fill(Mass_JESup/${mass_norm[${sampleIndex}]}, Total_EvtWt);
			h_mass_X_bin1_JES[1]->Fill(Mass_JESdown/${mass_norm[${sampleIndex}]}, Total_EvtWt);

			//PES
			Photon4Vec_PESup.SetPxPyPzE(phoPx*(1+Photon_pecUnc), phoPy*(1+Photon_pecUnc), phoPz*(1+Photon_pecUnc), (*phoE)[PC]*(1+Photon_pecUnc));
			Photon4Vec_PESdown.SetPxPyPzE(phoPx*(1-Photon_pecUnc), phoPy*(1-Photon_pecUnc), phoPz*(1-Photon_pecUnc), (*phoE)[PC]*(1-Photon_pecUnc));

			GJ4Vec_PESup.SetPxPyPzE(Photon4Vec_PESup[0]+JetPx, Photon4Vec_PESup[1]+JetPy, Photon4Vec_PESup[2]+JetPz, Photon4Vec_PESup[3]+(*jetEn)[JC]);
			GJ4Vec_PESdown.SetPxPyPzE(Photon4Vec_PESdown[0]+JetPx, Photon4Vec_PESdown[1]+JetPy, Photon4Vec_PESdown[2]+JetPz, Photon4Vec_PESdown[3]+(*jetEn)[JC]);

			Mass_PESup = pow((GJ4Vec_PESup[3]*GJ4Vec_PESup[3] - GJ4Vec_PESup[0]*GJ4Vec_PESup[0] - GJ4Vec_PESup[1]*GJ4Vec_PESup[1] - GJ4Vec_PESup[2]*GJ4Vec_PESup[2]), 0.5);
			Mass_PESdown = pow((GJ4Vec_PESdown[3]*GJ4Vec_PESdown[3] - GJ4Vec_PESdown[0]*GJ4Vec_PESdown[0] - GJ4Vec_PESdown[1]*GJ4Vec_PESdown[1] - GJ4Vec_PESdown[2]*GJ4Vec_PESdown[2]), 0.5);

			h_profile_PES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			h_profile_PES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			h_mass_PES[0][0]->Fill(Mass_PESup, Total_EvtWt);
			h_mass_PES[1][0]->Fill(Mass_PESdown, Total_EvtWt);

			h_mass_X_bin1_PES[0]->Fill(Mass_PESup/${mass_norm[${sampleIndex}]}, Total_EvtWt);
			h_mass_X_bin1_PES[1]->Fill(Mass_PESdown/${mass_norm[${sampleIndex}]}, Total_EvtWt);


			if(Pass_GJInvtMass && IfMassCut){
			    h_CutFlow->Fill(10.5);
			    h_CutFlowWithWts->Fill(10.5, Total_EvtWt);
			    h_CutFlowWithWts_XSm5->Fill(10.5, Total_EvtWt_XSm5);
			    h_CutFlowWithWts_XSp5->Fill(10.5, Total_EvtWt_XSp5);

			    //JES
			    h_profile_JES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			    h_profile_JES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			    h_mass_JES[0][1]->Fill(Mass_JESup, Total_EvtWt);
			    h_mass_JES[1][1]->Fill(Mass_JESdown, Total_EvtWt);

			    //PES
			    h_profile_PES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			    h_profile_PES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			    h_mass_PES[0][1]->Fill(Mass_PESup, Total_EvtWt);
			    h_mass_PES[1][1]->Fill(Mass_PESdown, Total_EvtWt);


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
   
   }//for jentry
}//Loop()


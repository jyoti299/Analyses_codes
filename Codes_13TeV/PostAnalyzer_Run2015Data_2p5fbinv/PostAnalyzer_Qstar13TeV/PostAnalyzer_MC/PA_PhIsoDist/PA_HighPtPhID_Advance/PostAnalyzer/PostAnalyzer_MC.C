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
   Lumi = 2197.327; // (pb^{-1}), COMPLETE 2015 25ns DATA WITH GOLDEN JSON

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
   //   cout << "Total entries: " << nentries << endl; 
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //     cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Uncomment this in script  
      //Lumi_EvtWt = Lumi*(${XS[${sampleIndex}]}/${totalEvents[${sampleIndex}]});

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
      PC = PassHighPtPhotonID();
      JC = GoodJet(PC);
      if(PC > -1) Pass_PhoPt = ((*phoEt)[PC] > Cut_Photon_pt);
      if(PC > -1) Pass_PhoEtaEB = (fabs((*phoSCEta)[PC]) <= Cut_Photon_eta);
      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) <= Cut_Jet_eta);
      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(PC > -1) PC_Gen = MatchedPhoton(PC); 
      if(JC > -1) JC_Gen = MatchedJet(JC);

      //=====================================Standard skelton===========================================
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC > -1){
	    if(Pass_PhoPt){
	      if(Pass_PhoEtaEB){
		if(JC > -1){
		  if(Pass_JetPt){
		    if(Pass_JetEta){
		      if(Pass_GJdPhi){
			if(Pass_GJdEta){
			  if(Pass_GJInvtMass && IfMassCut){
			  }//Pass_GJInvtMass
			}//Pass_GJdEta
		      }//Pass_GJdPhi
		    }//Pass_JetEta
		  }//Pass_JetPt
		}//JC > -1
	      }//Pass_PhoEtaEB
	    }//Pass_PhoPt
	  }//PC > -1
	}//HasPrimaryVtx
      }//Pass_HLT
      //==================================================================================================

      //Extra for advance
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(PC > -1){
	    
	    h_PFChIsoVsPt[0]                 ->Fill((*phoEt)[PC], (*phoPFChIso)[PC], PreBTag_EvtWt);
      	    h_PFChIsoVsRho[0]                ->Fill(rho, (*phoPFChIso)[PC], PreBTag_EvtWt);
      	    h_prof_PFChIsoVsPt[0]            ->Fill((*phoEt)[PC], (*phoPFChIso)[PC], PreBTag_EvtWt);
	    h_prof_PFChIsoVsRho[0]           ->Fill(rho, (*phoPFChIso)[PC], PreBTag_EvtWt);
	    
	    h_PFNeuIsoVsPt[0]             ->Fill((*phoEt)[PC], (*phoPFNeuIso)[PC], PreBTag_EvtWt);
	    h_PFNeuIsoVsRho[0]            ->Fill(rho, (*phoPFNeuIso)[PC], PreBTag_EvtWt);
	    h_prof_PFNeuIsoVsPt[0]        ->Fill((*phoEt)[PC], (*phoPFNeuIso)[PC], PreBTag_EvtWt);
	    h_prof_PFNeuIsoVsRho[0]       ->Fill(rho, (*phoPFNeuIso)[PC], PreBTag_EvtWt);

	    h_PFPhIsoVsPt[0][0]              ->Fill((*phoEt)[PC], (*phoPFPhoIso)[PC], PreBTag_EvtWt);
	    h_PFPhIsoVsRho[0][0]             ->Fill(rho, (*phoPFPhoIso)[PC], PreBTag_EvtWt);
	    h_prof_PFPhIsoVsPt[0][0]         ->Fill((*phoEt)[PC], (*phoPFPhoIso)[PC], PreBTag_EvtWt);
      	    h_prof_PFPhIsoVsRho[0][0]        ->Fill(rho, (*phoPFPhoIso)[PC], PreBTag_EvtWt);
				
	    Double_t Kappa = 0.002; // in units of (1/GeV)
	    Double_t Alpha = 0.0; //in units of GeV
	    Double_t Area = 0.0;

	    if((fabs((*phoSCEta)[PC])) <= 1.4442){
	      Alpha = 1.5;
	      Area = EAPhotons_HighPtID((*phoSCEta)[PC]);
	    }
	    if((fabs((*phoSCEta)[PC])) > 1.566){
	      Alpha = 2.0;
	      Area = EAPhotons_HighPtID((*phoSCEta)[PC]);
	    }

	    h_PFPhIsoVsPt[0][1]       ->Fill((*phoEt)[PC], (Alpha + (*phoPFPhoIso)[PC] - rho*Area - Kappa*(*phoEt)[PC]), PreBTag_EvtWt);
	    h_PFPhIsoVsRho[0][1]      ->Fill(rho, (Alpha + (*phoPFPhoIso)[PC] - rho*Area - Kappa*(*phoEt)[PC]), PreBTag_EvtWt);
            h_prof_PFPhIsoVsPt[0][1]  ->Fill((*phoEt)[PC], (Alpha + (*phoPFPhoIso)[PC] - rho*Area - Kappa*(*phoEt)[PC]), PreBTag_EvtWt);
     	    h_prof_PFPhIsoVsRho[0][1] ->Fill(rho, (Alpha + (*phoPFPhoIso)[PC] - rho*Area - Kappa*(*phoEt)[PC]), PreBTag_EvtWt);

          }
        }
      }   

   }//for jentry
}//Loop()

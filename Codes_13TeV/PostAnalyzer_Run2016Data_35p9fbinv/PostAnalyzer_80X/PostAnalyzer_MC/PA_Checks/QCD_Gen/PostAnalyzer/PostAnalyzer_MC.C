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
   Lumi = 24487.0; //(pb^{-1}) //Run2016BCDEFG_PromptReco

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

   Cut_GJdPhi = 2.5;
   Cut_GJdEta = 2.0;
   
   Cut_GJInvtMass = 695.0;
   
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

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   //   PileupReWeighting();

   //Defining CSVc2 bTag Operating Point (LOOSE, MEDIUM, TIGHT OR RESHAPING (for boosted btag discs))
   BTagEntry::OperatingPoint CSV_OP = BTagEntry::OP_MEDIUM; // required for SF calculation
   std::string CSV_WP = "M"; // required for Tagger (L,M or T)

   //Event For loop starts from here
   Long64_t nentries = fChain->GetEntries();
   cout << "<Total entries: " << nentries << endl; 
   Long64_t nbytes = 0, nb = 0;

   /*
   //running for loop to get the total genweight
   Long64_t nb1 = 0;
   Float_t Tot_genWt = 0;
   for (Long64_t bentry=0; bentry<nentries;bentry++) {
     Long64_t kentry = LoadTree(bentry);
     nb1 = fChain->GetEntry(bentry);
     Tot_genWt += genWeight;
   }
   cout << "Tot genWt = " << Tot_genWt << endl;
   //---------------------------------------------
   */
   Int_t n_tot = 0;
   Int_t GJ_Evt = 0; 
   Int_t QCD_Evt = 0; 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //      cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Uncomment this in script  
      //Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/nentries;

      //To be removed in script                                                                       
      //-----------------------                                                                                                                     
      Lumi_EvtWt = (Lumi*XS)/9956130;//2000069
      //-----------------------                                                                                                 

      PU_EvtWt = 1;//PUWeights((*puTrue)[0]); //Since TrueNumofInt is same for an event, so any value of the vector puTrue can be taken. 
                                          //(see definition here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information)
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
      Pass_GJdPhi = true;
      Pass_GJdEta = true;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = true;
      
      //Running different functions     
      Pass_HLT = true; //Always true for MC
      //  HasPrimaryVtx = GoodPrimaryVtx(GoodVertex);
      HasPrimaryVtx = hasGoodVtx;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons("loose");  //All photons passing loose id, residual spikes and pt > 30.0

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
      GoodIsoJets = GoodJets(PC, "tight"); // All jets passing dR > 0.5, tight jet id and pt > 30.0
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) < Cut_Jet_eta);
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      //      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      //      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);

      if(nPho > 0 && nJet > 0){
	//Filling the histogram wihtout any evts wts and with evts wts. 
	h_PhotonPt[0]->Fill((*phoEt)[0]);
	h_JetPt[0]->Fill((*jetPt)[0]);
	h_PhotonPt[1]->Fill((*phoEt)[0], PreBTag_EvtWt);
	h_JetPt[1]->Fill((*jetPt)[0], PreBTag_EvtWt);
      
	Int_t pc_gen0 = -1;
	pc_gen0 = MatchedGenPhotonToReco(0);
	if(pc_gen0 > -1){
	  h_genPhotonPt[0]->Fill((*mcPt)[pc_gen0]);
	  h_genPhotonPt[1]->Fill((*mcPt)[pc_gen0], PreBTag_EvtWt);
	}

	if((*jetGenJetPt)[0] >= 0){
	  h_genJetPt[0]->Fill((*jetGenJetPt)[0]);
	  h_genJetPt[1]->Fill((*jetGenJetPt)[0], PreBTag_EvtWt);
	}

	if((*jetGenPt)[0] >= 0){
	  h_genPartonPt[0]->Fill((*jetGenPt)[0]);
	  h_genPartonPt[1]->Fill((*jetGenPt)[0], PreBTag_EvtWt);
	}
      }

      h_CutFlow_qstar->Fill(0.5);
      h_CutFlowWt_qstar->Fill(0.5, PreBTag_EvtWt);
      
      /*
     //This is to check the fraction of prompt photons in gj/qcd samples before applying any selections
     //as after applying ID and ISO critieria, most events with non prompt photons will go away and hence checking the prompt 
     //photon fraction after the selection will not give correct answer.

      bool gPromptFound=false;
      for(Int_t i = 0; i < nMC; i++){	
	if(fabs((*mcPID)[i]) == 22 && ((*mcStatusFlag)[i] == 2 || (*mcStatusFlag)[i] == 3 || (*mcStatusFlag)[i] == 6)){
	  gPromptFound=true;
	}
      }
      if(gPromptFound) GJ_Evt++;
      if(!gPromptFound) QCD_Evt++;
      */    
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	if(HasPrimaryVtx){
	  if(GoodIsoPhotons.size() > 0){
	    if(PC > -1){
	      if(JC > -1){
		if(Pass_JetPt){
		  if(Pass_JetEta){
 		    if(Pass_GJdPhi){
		      if(Pass_GJdEta){

			h_PhotonPt[2]->Fill((*phoEt)[PC]);
			h_JetPt[2]->Fill((*jetPt)[JC]);
			h_PhotonPt[3]->Fill((*phoEt)[PC], PreBTag_EvtWt);
			h_JetPt[3]->Fill((*jetPt)[JC], PreBTag_EvtWt);
     
			Int_t pc_gen = -1;
			pc_gen = MatchedGenPhotonToReco(PC);
			if(pc_gen > -1){
			  h_genPhotonPt[2]->Fill((*mcPt)[pc_gen]);
			  h_genPhotonPt[3]->Fill((*mcPt)[pc_gen], PreBTag_EvtWt);
			}

			if((*jetGenJetPt)[JC] >= 0){
			  h_genJetPt[2]->Fill((*jetGenJetPt)[JC]);
			  h_genJetPt[3]->Fill((*jetGenJetPt)[JC], PreBTag_EvtWt);
			}

			if((*jetGenPt)[JC] >= 0){
			  h_genPartonPt[2]->Fill((*jetGenPt)[JC]);
			  h_genPartonPt[3]->Fill((*jetGenPt)[JC], PreBTag_EvtWt);
			}

			//Here checking the prompt photon fraction after the final selection
			n_tot++;	
			bool gPromptFound=false;
			for(Int_t i = 0; i < nMC; i++){ 
			  if(fabs((*mcPID)[i]) == 22 && ((*mcStatusFlag)[i] == 2 || (*mcStatusFlag)[i] == 3 || (*mcStatusFlag)[i] == 6)){
			    gPromptFound=true;
			  }
			}
			if(gPromptFound) GJ_Evt++;
			if(!gPromptFound) QCD_Evt++;
		      		      
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
 
      }//for jentry

   cout << "Total = " << n_tot << endl;
   cout << "GJ_Evt = " << GJ_Evt << endl;
   cout << "QCD_Evt = " << QCD_Evt << endl;

}//Loop()


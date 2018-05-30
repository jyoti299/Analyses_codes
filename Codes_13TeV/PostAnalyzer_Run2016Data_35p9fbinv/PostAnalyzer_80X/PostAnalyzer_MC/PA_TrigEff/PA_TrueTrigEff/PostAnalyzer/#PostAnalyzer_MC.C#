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
   Lumi = 36813.0; //(pb^{-1}) //Rereco-BCDEFG_PromptReco-H  (24487 pb-1)Run2016BCDEFG_PromptReco

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
   Cut_GJdEta = 1.5;
   
   Cut_GJInvtMass = 695.0;
   
   Cut_PhId = "loose";
   Cut_JetId = "tight";

   //triggers used in the study (efficiency of these three triggers is same as using all triggers)
   vector<ULong64_t> HLT;
   HLT.clear();
   HLT.push_back(HLT_Photon165_HE10_v);
   //HLT.push_back(HLT_Photon175_v);

   //triggers for the denominator of trigger turn on
   vector<ULong64_t> HLT_deno;
   HLT_deno.clear();
   HLT_deno.push_back(HLT_Photon75_v);
   HLT_deno.push_back(HLT_Photon90_v);
   HLT_deno.push_back(HLT_Photon120_v);

   //Triggers bits for jet matching
   vector<UInt_t> PhodenoTrigObjs;
   PhodenoTrigObjs.clear();
   PhodenoTrigObjs.push_back(hltEG75HEFilter);
   PhodenoTrigObjs.push_back(hltEG90HEFilter);
   PhodenoTrigObjs.push_back(hltEG120HEFilter);

   //Trigger bits for photon matching
   vector<UInt_t> PhoTrigObjs;
   PhoTrigObjs.push_back(hltEG165HE10Filter);
   //PhoTrigObjs.push_back(hltEG175HEFilter);

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
   PileupReWeighting();

   //Defining CSVc2 bTag Operating Point (LOOSE, MEDIUM, TIGHT OR RESHAPING (for boosted btag discs))
   BTagEntry::OperatingPoint CSV_OP = BTagEntry::OP_LOOSE; // required for SF calculation
   std::string CSV_WP = "L"; // required for Tagger (L,M or T)

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
     cout << "genWt = " << genWeight << endl;
     Tot_genWt += genWeight;
   }
   cout << "Tot genWt = " << Tot_genWt << endl;
   //---------------------------------------------
   */

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //    cout << "<Analyzing entry: " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Uncomment this in script  
      //Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/nentries;

      //To be removed in script                                                                       
      //-----------------------                                                                                                                     
      Lumi_EvtWt = (Lumi*XS)/9956130;//2000069
      //-----------------------                                                                                                 

      PU_EvtWt = PUWeights((*puTrue)[0]); //Since TrueNumofInt is same for an event, so any value of the vector puTrue can be taken. 
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
      //No dphi cut
      Pass_GJdPhi = true;
      Pass_GJdEta = false;
      Pass_CSVv2bTag = false;
      Pass_GJInvtMass = false;
      
      //Running different functions     
      Pass_HLT = true; //Always true for MC
      GoodVertex = nGoodVtx;
      if(GoodVertex > 0) HasPrimaryVtx = true;

      GoodIsoPhotons.clear();
      GoodIsoPhotons = GoodPhotons(Cut_PhId);  //All photons passing loose id, residual spikes and pt > 30.0

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
      GoodIsoJets = GoodJets(Cut_JetId); // All jets passing dR > 0.5, tight jet id and pt > 30.0
      if(GoodIsoJets.size() != 0) JC = GoodIsoJets[0];

      if(JC > -1) Pass_JetPt = ((*jetPt)[JC] > Cut_Jet_pt);
      if(JC > -1) Pass_JetEta = (fabs((*jetEta)[JC]) < Cut_Jet_eta);
      //      if(PC > -1 && JC > -1) Pass_GJdPhi =  ((GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC])) > Cut_GJdPhi);
      if(PC > -1 && JC > -1) Pass_GJdEta = ((GetdEta((*phoSCEta)[PC], (*jetEta)[JC])) < Cut_GJdEta);
      if(PC > -1 && JC > -1) Pass_GJInvtMass = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);             
      if(JC > -1) Pass_CSVv2bTag = CSVv2bTag(JC, CSV_WP);

      h_CutFlow_qstar->Fill(0.5);
      h_CutFlow_bstar->Fill(0.5);
      h_CutFlowWt_qstar->Fill(0.5, PreBTag_EvtWt);
      h_CutFlowWt_bstar->Fill(0.5, PreBTag_EvtWt);
      h_CutFlowTotalWt_bstar->Fill(0.5, PreBTag_EvtWt);

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //                     QSTAR
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(Pass_HLT){
	h_CutFlow_qstar->Fill(1.5);
	h_CutFlowWt_qstar->Fill(1.5, PreBTag_EvtWt);
         
	if(HasPrimaryVtx){
	  h_CutFlow_qstar->Fill(2.5);
	  h_CutFlowWt_qstar->Fill(2.5, PreBTag_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow_qstar->Fill(3.5);
	    h_CutFlowWt_qstar->Fill(3.5, PreBTag_EvtWt);

	    if(PC > -1){
	      h_CutFlow_qstar->Fill(4.5);
	      h_CutFlowWt_qstar->Fill(4.5, PreBTag_EvtWt);

	      h_CutFlow_qstar->Fill(5.5);
	      h_CutFlowWt_qstar->Fill(5.5, PreBTag_EvtWt);

	      if(JC > -1){
		h_CutFlow_qstar->Fill(6.5);
		h_CutFlowWt_qstar->Fill(6.5, PreBTag_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow_qstar->Fill(7.5);
		  h_CutFlowWt_qstar->Fill(7.5, PreBTag_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow_qstar->Fill(8.5);
		    h_CutFlowWt_qstar->Fill(8.5, PreBTag_EvtWt);

 		    if(Pass_GJdPhi){
		      h_CutFlow_qstar->Fill(9.5);
		      h_CutFlowWt_qstar->Fill(9.5, PreBTag_EvtWt);

		      if(Pass_GJdEta){
			h_CutFlow_qstar->Fill(10.5);
			h_CutFlowWt_qstar->Fill(10.5, PreBTag_EvtWt);

		      
			if(Pass_GJInvtMass){
			  h_CutFlow_qstar->Fill(11.5);
			  h_CutFlowWt_qstar->Fill(11.5, PreBTag_EvtWt);
			  
			
			}//if(Pass_GJInvtMass)			
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

      //---------------------------------------------------------------------------------   
      //Trigger Turn-on
      Bool_t PassHLTNum = false;
      Bool_t PassHLTDeno = false;
      Bool_t PassHLT_Pre = false;
      //      cout << PassHLT_Pre << endl;
      ULong64_t HLT_pre = HLT[0];

      PassHLTNum = PassHLT(HLT);
      PassHLTDeno = PassHLT(HLT_deno);
      PassHLT_Pre = PassHLT_Prescale(HLT_pre);

      if(nPho > 0){
	h_TrigPhotonPt[0]->Fill((*phoEt)[0]);
	if(PassHLTNum){
	  h_TrigPhotonPt[1]->Fill((*phoEt)[0]);
	}
      }

   }//for jentry

}//Loop()


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
   Lumi = 35866.0; //(pb^{-1}) //Re-miniAOD//  (36813.0 pb-1)Rereco-BCDEFG_PromptReco-H  (24487 pb-1)Run2016BCDEFG_PromptReco

   //To be removed in script
   //-----------------------
   Float_t XS = 5813.0;
   //-----------------------

   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 200.0; // GeV                                                                                                             
   Cut_Photon_eta = 1.4442;

   Cut_Jet_pt = 170.0; // GeV
   Cut_Jet_eta = 2.4;

   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 1.5;
   
   Cut_GJInvtMass = 700.0;
   
   Cut_PhId = "medium";
   Cut_JetId = "tight";

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

   for (Long64_t jentry=0; jentry<10/*nentries*/;jentry++) {
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

      //CSVv2 BTag Efficiencies
      for(Int_t ii = 0; ii < nJet; ii++){
	if(fabs((*jetHadFlvr)[ii]) == 5){
	  h_bEff_vs_pt[0]->Fill((*jetPt)[ii]);
	  h_bEff_vs_eta[0]->Fill((*jetEta)[ii]);
	  if((*jetCSV2BJetTags)[ii] > 0.5426){
	    h_bEff_vs_pt[1]->Fill((*jetPt)[ii]);
	    h_bEff_vs_eta[1]->Fill((*jetEta)[ii]);
	  }
	}
	if(fabs((*jetHadFlvr)[ii]) == 4){
	  h_cEff_vs_pt[0]->Fill((*jetPt)[ii]);
	  h_cEff_vs_eta[0]->Fill((*jetEta)[ii]);
	  if((*jetCSV2BJetTags)[ii] > 0.5426){
	    h_cEff_vs_pt[1]->Fill((*jetPt)[ii]);
	    h_cEff_vs_eta[1]->Fill((*jetEta)[ii]);
	  }
	}
	if(fabs((*jetHadFlvr)[ii]) == 0){
	  h_udsgEff_vs_pt[0]->Fill((*jetPt)[ii]);
	  h_udsgEff_vs_eta[0]->Fill((*jetEta)[ii]);
	  if((*jetCSV2BJetTags)[ii] > 0.5426){
	    h_udsgEff_vs_pt[1]->Fill((*jetPt)[ii]);
	    h_udsgEff_vs_eta[1]->Fill((*jetEta)[ii]);
	  }
	}
      }




      /*
      vector<Int_t> goodjets;
      goodjets.clear();
      for(Int_t ii = 0; ii < nJet; ii++){
	if(JetId(ii, "tight") && fabs((*jetEta)[ii]) < 2.4) goodjets.push_back(ii);
      }

      Int_t Selectedjet = -1;
      for(Int_t jj = 0; jj < goodjets.size(); jj++){
	Selectedjet = goodjets[jj];
	if(fabs((*jetHadFlvr)[Selectedjet]) == 5){
	  h_bEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_bEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}

	if(fabs((*jetHadFlvr)[Selectedjet]) == 4){
	  h_cEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_cEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}

	if(fabs((*jetHadFlvr)[Selectedjet]) == 0){
	  h_udsgEff_vs_pt[0]->Fill((*jetPt)[Selectedjet]);
	  h_udsgEff_vs_eta[0]->Fill((*jetEta)[Selectedjet]);
	}

	if(CSVv2bTag(Selectedjet, CSV_WP)){
	  Double_t SF_eff = 0;
	  Double_t Wt = 0;
          BTagEntry::JetFlavor JF_eff;
	  std::string sys_type_eff = "central"; //central is required to get scale factors (up and down for uncertainties) 
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 5) JF_eff = BTagEntry::FLAV_B; //b
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 4) JF_eff = BTagEntry::FLAV_C; //c
	  if(fabs((*jetHadFlvr)[Selectedjet]) == 0) JF_eff = BTagEntry::FLAV_UDSG; //u,d,s,g,undefined

	  SF_eff = CSVv2bTagSF_auto(CSV_OP, JF_eff, sys_type_eff, (*jetPt)[Selectedjet], (*jetEta)[Selectedjet]);	

	  Wt = BTagEventWeight(SF_eff, 1); // SF		     
	  if(JF_eff == BTagEntry::FLAV_B){
	    h_bEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_bEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_bEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_bEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
	  }

	  if(JF_eff == BTagEntry::FLAV_C){
	    h_cEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_cEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_cEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_cEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
  	  }

	  if(JF_eff == BTagEntry::FLAV_UDSG){
	    h_udsgEff_vs_pt[1]->Fill((*jetPt)[Selectedjet]);
	    h_udsgEff_vs_pt[2]->Fill((*jetPt)[Selectedjet], Wt);
	    h_udsgEff_vs_eta[1]->Fill((*jetEta)[Selectedjet]);
	    h_udsgEff_vs_eta[2]->Fill((*jetEta)[Selectedjet], Wt);
	  }
	}
      }
      */
      //---------------------------
   
   }//for jentry

}//Loop()


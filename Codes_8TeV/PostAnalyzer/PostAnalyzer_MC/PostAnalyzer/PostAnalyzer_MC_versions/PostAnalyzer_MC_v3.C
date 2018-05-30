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
//      Root > .L PostAnalyzer_MC.C
//      Root > PostAnalyzer_MC t
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

   //Luminosity of Data Compared
   Lumi = 19711.4; // (pb^{-1})

   //To be removed in script
   //-----------------------
   Float_t XS = 30.12207;
   //-----------------------

   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 170.0;
   Cut_Photon_eta = 1.4442;
   Cut_Jet_pt = 170.0;
   Cut_Jet_eta = 2.5;
   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 2.0;
   Cut_GJInvtMass = 560.0; 

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

   //*******************************************************************************************************//
   //Defining the histogram from filling number of events after various cuts

   const int nbins = 25;
   const int nbinsWt = 37;
   TString CutFlowLabels[nbins] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag", "PassCSVMBTag", "PassCSVTBTag", "TrueBJets", "TrueBJetsPassingCSVL", "TrueBJetsPassingCSVM", "TrueBJetsPassingCSVT", "NonBJets", "NonBJetsPassingCSVL", "NonBJetsPassingCSVM", "NonBJetsPassingCSVT", "PassPhotonID_NoLICTD"};

   TString CutFlowLabelsWithWts[nbinsWt] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag_Expected", "PassCSVLBTag_noErr_0bTag", "PassCSVLBTag_noErr_1bTag", "PassCSVLBTag_pErr_0bTag", "PassCSVLBTag_pErr_1bTag", "PassCSVMBTag_Expected", "PassCSVMBTag_noErr_0bTag", "PassCSVMBTag_noErr_1bTag", "PassCSVMBTag_pErr_0bTag", "PassCSVMBTag_pErr_1bTag", "PassCSVTBTag_Expected", "PassCSVTBTag_noErr_0bTag", "PassCSVTBTag_noErr_1bTag", "PassCSVTBTag_pErr_0bTag", "PassCSVTBTag_pErr_1bTag", "TrueBJets", "TrueBJetsPassingCSVL", "TrueBJetsPassingCSVM", "TrueBJetsPassingCSVT", "NonBJets", "NonBJetsPassingCSVL", "NonBJetsPassingCSVM", "NonBJetsPassingCSVT", "PassPhotonID_NoLICTD"};

   h_CutFlowTable = new TH1F("h_CutFlowTable", "Events Passing Various Cuts", nbins, 0, nbins);
   h_CutFlowTable->GetYaxis()->SetTitle("Events");            h_CutFlowTable->GetYaxis()->CenterTitle();
   h_CutFlowTableWithWeights = new TH1F("h_CutFlowTableWithWeights", "Events Passing Various Cuts With Weights", nbinsWt, 0, nbinsWt);
   h_CutFlowTableWithWeights->GetYaxis()->SetTitle("Events"); h_CutFlowTableWithWeights->GetYaxis()->CenterTitle(); 
   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->GetXaxis()->SetBinLabel(i+1, CutFlowLabels[i]);    
   }
   for(int i = 0; i < nbinsWt; i++){
     h_CutFlowTableWithWeights->GetXaxis()->SetBinLabel(i+1, CutFlowLabelsWithWts[i]);
   }

   Long64_t CutFlowNumber[nbins];
   Double_t CutFlowNumberWithWeights[nbinsWt]; //This needs to be taken as double as it will consider evnt wts which are doubles, if we take it to
                                              //Long then if evnt wt is less than 1 then it will take that as 0 upon rounding off.
   for(int i = 0; i < nbins; i++){
     CutFlowNumber[i] = 0;
   }
   for(int i = 0; i < nbinsWt; i++){
     CutFlowNumberWithWeights[i] = 0.0;
   }
    
   //These Variables to check with the old 7TeV Scale Factor recommendation
   Double_t CutFlowNum_OldRec7T_passingCSVL_noErr_0bTag = 0.0;
   Double_t CutFlowNum_OldRec7T_passingCSVL_noErr_1bTag = 0.0;
   Double_t CutFlowNum_OldRec7T_passingCSVL_pErr_0bTag = 0.0;
   Double_t CutFlowNum_OldRec7T_passingCSVL_pErr_1bTag = 0.0;

   //These variables are defined just to check deta and dphi dependence on the b disc efficiency in choosing leading jet as b jet
   Double_t njets = 0;
   Double_t ncsvlbjets = 0;
   Double_t ncsvmbjets = 0;
   Double_t ncsvtbjets = 0;

   //********************************************************************************************************//

   //Make it true if mass cut to be applied otherwise false
   IfMassCut = false;

   //*********************************************************************************************************//

   Long64_t nentries = fChain->GetEntries();
   cout << "no. of entries " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;

   for(Long64_t jentry = 0; jentry < nentries; jentry++){

     //     cout << "++++++++++++++++++Analyzing entry++++++++++++" << jentry << endl;

     //Uncomment this in script
     //     Lumi_EvtWt = Lumi*(${XS[${sampleIndex}]}/${totalEvents[${sampleIndex}]});

     //To be removed in script
     //-----------------------
     Lumi_EvtWt = Lumi*XS/Double_t(nentries);//2000069 
     //-----------------------

     PU_EvtWt = PUWeights(trueNumofInteractions);
     PreBTag_EvtWt = Lumi_EvtWt * PU_EvtWt;

     h_trueNumofInteractions->Fill(trueNumofInteractions);
     h_trueNumofInteractions_withpuWt->Fill(trueNumofInteractions, PU_EvtWt);

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

     PC = -1;
     JC = -1;
     GoodVertex = 0;

     Pass_HLT = true;
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

      //      Pass_HLT = PassHLT();
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
     if(IfMassCut == true){
       if(PC > -1 && JC > -1){
	 Pass_GJInvtMasscut = ((GetInvtMass(PC, JC)) > Cut_GJInvtMass);
       }
     }
     if(JC > -1) Pass_CSVLBTag = PassCSVLBTag(JC);
     if(JC > -1) Pass_CSVMBTag = PassCSVMBTag(JC);
     if(JC > -1) Pass_CSVTBTag = PassCSVTBTag(JC);

     CutFlowNumber[0]++;
     CutFlowNumberWithWeights[0] += PreBTag_EvtWt;

     if(Pass_HLT){
       CutFlowNumber[1]++;
       CutFlowNumberWithWeights[1] += PreBTag_EvtWt;

       if(NoScrapingEvt){
	 CutFlowNumber[2]++;
	 CutFlowNumberWithWeights[2] += PreBTag_EvtWt;

	 if(HasPrimaryVtx){
	   CutFlowNumber[3]++;
	   CutFlowNumberWithWeights[3] += PreBTag_EvtWt;

	   if(PC > -1){
	     CutFlowNumber[4]++;
	     CutFlowNumberWithWeights[4] += PreBTag_EvtWt;

	     if(Pass_PhoPtcut){
	       CutFlowNumber[5]++;
	       CutFlowNumberWithWeights[5] += PreBTag_EvtWt;

	       if(Pass_PhoEtaEBcut){
		 CutFlowNumber[6]++;
		 CutFlowNumberWithWeights[6] += PreBTag_EvtWt;

		 if(JC > -1){
		   CutFlowNumber[7]++;
		   CutFlowNumberWithWeights[7] += PreBTag_EvtWt;

		   if(Pass_JetPtcut){
		     CutFlowNumber[8]++;
		     CutFlowNumberWithWeights[8] += PreBTag_EvtWt;

		     if(Pass_JetEtacut){
		       CutFlowNumber[9]++;
		       CutFlowNumberWithWeights[9] += PreBTag_EvtWt;

		       if(Pass_GJdPhicut){
			 CutFlowNumber[10]++;
			 CutFlowNumberWithWeights[10] += PreBTag_EvtWt;

			 if(Pass_GJdEtacut){
			   CutFlowNumber[11]++;
			   CutFlowNumberWithWeights[11] += PreBTag_EvtWt;
      
			   if(Pass_GJInvtMasscut){
			     CutFlowNumber[12]++;
			     CutFlowNumberWithWeights[12] += PreBTag_EvtWt;
 
			     h_PhotonPt->Fill((*Photon_pt)[PC], PreBTag_EvtWt);
			     h_PhotonEta->Fill((*Photon_SC_eta)[PC], PreBTag_EvtWt);
			     h_PhotonPhi->Fill((*Photon_phi)[PC], PreBTag_EvtWt);
			     h_PhotonSigmaIEtaIEta->Fill((*Photon_SigmaIEtaIEta)[PC], PreBTag_EvtWt);
			     h_PhotonSigmaIPhiIPhi->Fill((*Photon_SigmaIPhiIPhi)[PC], PreBTag_EvtWt);
			     h_Photon_r9->Fill((*Photon_r9)[PC], PreBTag_EvtWt);
			     h_Photon_SingleTowerHoE->Fill((*Photon_SingleTowerHoE)[PC], PreBTag_EvtWt);
			     h_Photon_PFIsoCharged03->Fill((*PFIsoCharged03)[PC], PreBTag_EvtWt);
			     h_Photon_PFIsoNeutral03->Fill((*PFIsoNeutral03)[PC], PreBTag_EvtWt);
			     h_Photon_PFIsoPhoton03->Fill((*PFIsoPhoton03)[PC], PreBTag_EvtWt);
			     h_Photon_PFIsoSum03->Fill((*PFIsoSum03)[PC], PreBTag_EvtWt);

			     h_JetPt->Fill((*PFPatJet_pt)[JC], PreBTag_EvtWt);
			     h_JetEta->Fill((*PFPatJet_eta)[JC], PreBTag_EvtWt);
			     h_JetPhi->Fill((*PFPatJet_phi)[JC], PreBTag_EvtWt);
			     h_Jet_NeutralHadEnergyFrac->Fill((*PFPatJet_NeutralHadEnergyFrac)[JC], PreBTag_EvtWt);
			     h_Jet_NeutralEmEnergyFrac->Fill((*PFPatJet_NeutralEmEnergyFrac)[JC], PreBTag_EvtWt);
			     h_Jet_NConstituents->Fill((*PFPatJet_NConstituents)[JC], PreBTag_EvtWt);
			     h_Jet_ChargedHadEnergyFrac->Fill((*PFPatJet_ChargedHadEnergyFrac)[JC], PreBTag_EvtWt);
			     h_Jet_ChargedMult->Fill((*PFPatJet_ChargedMult)[JC], PreBTag_EvtWt);
			     h_Jet_ChargedEmEnergyFrac->Fill((*PFPatJet_ChargedEmEnergyFrac)[JC], PreBTag_EvtWt);

			     h_GJetInvtMass_binning40->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			     h_GJetInvtMass_VariableBinning->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
			     h_GJetdEta->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), PreBTag_EvtWt);
			     h_GJetdPhi->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), PreBTag_EvtWt);
			     h_GJetdR->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), PreBTag_EvtWt);

			     h_BJetDiscByCSV->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], PreBTag_EvtWt);

			     h_PC->Fill(PC, PreBTag_EvtWt);
			     h_JC->Fill(JC, PreBTag_EvtWt);

			     h_Jet_PartonFlavor->Fill((*PFPatJet_JetPartonFlavor)[JC], PreBTag_EvtWt);

			     if(Pass_CSVLBTag){
			       CutFlowNumber[13]++;
			       CutFlowNumberWithWeights[13] += PreBTag_EvtWt;			

			       double SF_CSVL = BTagScaleFactor_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
			       double SFerr_CSVL = BTagSFerr_CSVL((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       Total_CSVLEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL, 0);
			       Total_CSVLEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL, 1);
			       Total_CSVLEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL+SFerr_CSVL, 0);
			       Total_CSVLEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL+SFerr_CSVL, 1);

			       CutFlowNumberWithWeights[14] += Total_CSVLEvtWt_noErr_0bTag;
			       CutFlowNumberWithWeights[15] += Total_CSVLEvtWt_noErr_1bTag;
			       CutFlowNumberWithWeights[16] += Total_CSVLEvtWt_pErr_0bTag;
			       CutFlowNumberWithWeights[17] += Total_CSVLEvtWt_pErr_1bTag;

			       h_CSVL_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_noErr_0bTag);

			       h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_noErr_0bTag);

			       h_CSVL_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_noErr_1bTag);

			       h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_noErr_1bTag);

			       h_CSVL_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       h_CSVL_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_pErr_1bTag);


			       //To check with SF 7TeV Recommendation
			       double SF_CSVL_OldRec7T = BTagScaleFactor_CSVL_OldRec7T((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
			       double SFerr_CSVL_OldRec7T = BTagSFerr_CSVL_OldRec7T((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       double Total_CSVLEvtWt_noErr_0bTag_OldRec7T  = PreBTag_EvtWt * BTagEventWeight(SF_CSVL_OldRec7T, 0);
			       double Total_CSVLEvtWt_noErr_1bTag_OldRec7T  = PreBTag_EvtWt * BTagEventWeight(SF_CSVL_OldRec7T, 1);
			       double Total_CSVLEvtWt_pErr_0bTag_OldRec7T  = PreBTag_EvtWt * BTagEventWeight(SF_CSVL_OldRec7T+SFerr_CSVL_OldRec7T, 0);
			       double Total_CSVLEvtWt_pErr_1bTag_OldRec7T  = PreBTag_EvtWt * BTagEventWeight(SF_CSVL_OldRec7T+SFerr_CSVL_OldRec7T, 1);

			       CutFlowNum_OldRec7T_passingCSVL_noErr_0bTag += Total_CSVLEvtWt_noErr_0bTag_OldRec7T;
			       CutFlowNum_OldRec7T_passingCSVL_noErr_1bTag += Total_CSVLEvtWt_noErr_1bTag_OldRec7T;
			       CutFlowNum_OldRec7T_passingCSVL_pErr_0bTag += Total_CSVLEvtWt_pErr_0bTag_OldRec7T;
			       CutFlowNum_OldRec7T_passingCSVL_pErr_1bTag += Total_CSVLEvtWt_pErr_1bTag_OldRec7T;

			     }
			     if(Pass_CSVMBTag){
			       CutFlowNumber[14]++;
			       CutFlowNumberWithWeights[18] += PreBTag_EvtWt;

			       double SF_CSVM = BTagScaleFactor_CSVM((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
			       double SFerr_CSVM = BTagSFerr_CSVM((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       Total_CSVMEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM, 0);
			       Total_CSVMEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM, 1);
			       Total_CSVMEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM+SFerr_CSVM, 0);
			       Total_CSVMEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM+SFerr_CSVM, 1);

			       CutFlowNumberWithWeights[19] += Total_CSVMEvtWt_noErr_0bTag;
			       CutFlowNumberWithWeights[20] += Total_CSVMEvtWt_noErr_1bTag;
			       CutFlowNumberWithWeights[21] += Total_CSVMEvtWt_pErr_0bTag;
			       CutFlowNumberWithWeights[22] += Total_CSVMEvtWt_pErr_1bTag;

			       h_CSVM_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_noErr_0bTag);

			       h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_noErr_0bTag);

			       h_CSVM_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_noErr_1bTag);

			       h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_noErr_1bTag);

			       h_CSVM_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_pErr_0bTag);

			       h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_pErr_0bTag);

			       h_CSVM_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       h_CSVM_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_pErr_1bTag);


			     }
			     if(Pass_CSVTBTag){
			       CutFlowNumber[15]++;
			       CutFlowNumberWithWeights[23] += PreBTag_EvtWt;

			       double SF_CSVT = BTagScaleFactor_CSVT((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
			       double SFerr_CSVT = BTagSFerr_CSVT((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       Total_CSVTEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT, 0);
			       Total_CSVTEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT, 1);
			       Total_CSVTEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT+SFerr_CSVT, 0);
			       Total_CSVTEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT+SFerr_CSVT, 1);

			       CutFlowNumberWithWeights[24] += Total_CSVTEvtWt_noErr_0bTag;
			       CutFlowNumberWithWeights[25] += Total_CSVTEvtWt_noErr_1bTag;
			       CutFlowNumberWithWeights[26] += Total_CSVTEvtWt_pErr_0bTag;
			       CutFlowNumberWithWeights[27] += Total_CSVTEvtWt_pErr_1bTag;

			       h_CSVT_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_noErr_1bTag);

			       h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_noErr_1bTag);

			       h_CSVT_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_pErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_pErr_0bTag);

			       h_CSVT_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       h_CSVT_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_pErr_1bTag);


			     }//END_OF if(Pass_CSVTBTag)

			     //----------------------------------------------------------------------------------
			     //This part to get true b tag efficiency and Mistag rate of the three discriminators
			     //----------------------------------------------------------------------------------
			     if((*PFPatJet_JetPartonFlavor)[JC] == 5){
			       CutFlowNumber[16]++;
			       CutFlowNumberWithWeights[28] += PreBTag_EvtWt;

			       h_TrueBJetPt->Fill((*PFPatJet_pt)[JC], PreBTag_EvtWt);
			       h_TrueBJetEta->Fill((*PFPatJet_eta)[JC], PreBTag_EvtWt);

			       if(Pass_CSVLBTag){
				 CutFlowNumber[17]++;
				 CutFlowNumberWithWeights[29] += Total_CSVLEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVL->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVL->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVMBTag){
				 CutFlowNumber[18]++;
				 CutFlowNumberWithWeights[30] += Total_CSVMEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVM->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVM->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVTBTag){
				 CutFlowNumber[19]++;
				 CutFlowNumberWithWeights[31] += Total_CSVTEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVT->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVT->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       }

			     }else{
			       CutFlowNumber[20]++;
			       CutFlowNumberWithWeights[32] += PreBTag_EvtWt;

			       h_NonBJetPt->Fill((*PFPatJet_pt)[JC], PreBTag_EvtWt);
			       h_NonBJetEta->Fill((*PFPatJet_eta)[JC], PreBTag_EvtWt);

			       if(Pass_CSVLBTag){
				 CutFlowNumber[21]++;
				 CutFlowNumberWithWeights[33] += Total_CSVLEvtWt_pErr_1bTag;

				 h_NonBJetPtPassingCSVL->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
				 h_NonBJetEtaPassingCSVL->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVMBTag){
				 CutFlowNumber[22]++;
				 CutFlowNumberWithWeights[34] += Total_CSVMEvtWt_pErr_1bTag;

				 h_NonBJetPtPassingCSVM->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
				 h_NonBJetEtaPassingCSVM->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVTBTag){
				 CutFlowNumber[23]++;
				 CutFlowNumberWithWeights[35] += Total_CSVTEvtWt_pErr_1bTag;

				 h_NonBJetPtPassingCSVT->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_1bTag);
				 h_NonBJetEtaPassingCSVT->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       }
			     }
			     //--------------------------------------------------------------------------------

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
     PC_NoLICTD = GetPhotonPassingAllCuts_NoLICTD();
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
			     CutFlowNumber[24]++;
			     CutFlowNumberWithWeights[36] += PreBTag_EvtWt;
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
     h_nPhotons->Fill(nPhotons, PreBTag_EvtWt);
     h_nJets->Fill(nJets, PreBTag_EvtWt);
     h_nCSVLBJets->Fill(nCSVLBJets, Total_CSVLEvtWt_pErr_1bTag);
     h_nCSVMBJets->Fill(nCSVMBJets, Total_CSVMEvtWt_pErr_1bTag);
     h_nCSVTBJets->Fill(nCSVTBJets, Total_CSVTEvtWt_pErr_1bTag);
     h_CSVL_BJetsFrac->Fill(frac_CSVL, Total_CSVLEvtWt_pErr_1bTag/PreBTag_EvtWt);
     h_CSVM_BJetsFrac->Fill(frac_CSVM, Total_CSVMEvtWt_pErr_1bTag/PreBTag_EvtWt);
     h_CSVT_BJetsFrac->Fill(frac_CSVT, Total_CSVTEvtWt_pErr_1bTag/PreBTag_EvtWt);

     //Filling histogram for PhotonIdx vs PhotonPt
     for(int i = 0; i < Photon_n; i++){     
       h_PhotonIdxVsPt->Fill((*Photon_pt)[i], i, PreBTag_EvtWt);       
     }
     //Filling histogram for JetIdx vs JetPt
     for(int i = 0; i < PFPatJet_n; i++){       
       h_JetIdxVsPt->Fill((*PFPatJet_pt)[i], i, PreBTag_EvtWt);
     }
     //Filling histogram for CSVL-BJetIdx vs CSVL-BJetPt
     for(int i = 0; i < PFPatJet_n; i++){
       if(PassCSVLBTag(i) == 1){
	 h_CSVLBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, Total_CSVLEvtWt_pErr_1bTag);
       }
     }
     //Filling histogram for CSVM-BJetIdx vs CSVM-BJetPt
     for(int i = 0; i < PFPatJet_n; i++){
       if(PassCSVMBTag(i) == 1){
	 h_CSVMBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, Total_CSVMEvtWt_pErr_1bTag);
       }
     }
     //Filling histogram for CSVT-BJetIdx vs CSVT-BJetPt
     for(int i = 0; i < PFPatJet_n; i++){
       if(PassCSVTBTag(i) == 1){
	 h_CSVTBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, Total_CSVTEvtWt_pErr_1bTag);
       }
     }

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
		       njets += PreBTag_EvtWt;
		       if(Pass_CSVLBTag){
			 ncsvlbjets += PreBTag_EvtWt;    
		       }
		       if(Pass_CSVMBTag){
			 ncsvmbjets += PreBTag_EvtWt;
		       }
		       if(Pass_CSVTBTag){
			 ncsvtbjets += PreBTag_EvtWt;
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
   for(int i = 0; i < nbinsWt; i++){
     h_CutFlowTableWithWeights->SetBinContent(i+1, CutFlowNumberWithWeights[i]);
   }


   //Efficiency of various cuts
   Eff_PassHLT                           = CutFlowNumberWithWeights[1]/CutFlowNumberWithWeights[0];
   Eff_PassScraping                      = CutFlowNumberWithWeights[2]/CutFlowNumberWithWeights[1];
   Eff_PassPrimaryVtx                    = CutFlowNumberWithWeights[3]/CutFlowNumberWithWeights[2];
   Eff_PassPhotonID                      = CutFlowNumberWithWeights[4]/CutFlowNumberWithWeights[3];
   Eff_PassPhotonPt                      = CutFlowNumberWithWeights[5]/CutFlowNumberWithWeights[4];
   Eff_PassPhotonEta                     = CutFlowNumberWithWeights[6]/CutFlowNumberWithWeights[5];
   Eff_PassJetID                         = CutFlowNumberWithWeights[7]/CutFlowNumberWithWeights[6];
   Eff_PassJetPt                         = CutFlowNumberWithWeights[8]/CutFlowNumberWithWeights[7];
   Eff_PassJetEta                        = CutFlowNumberWithWeights[9]/CutFlowNumberWithWeights[8];
   Eff_PassDPhi                          = CutFlowNumberWithWeights[10]/CutFlowNumberWithWeights[9];
   Eff_PassDEta                          = CutFlowNumberWithWeights[11]/CutFlowNumberWithWeights[10];
   Eff_PassGJInvtMass                    = CutFlowNumberWithWeights[12]/CutFlowNumberWithWeights[11];
   Eff_PassCSVLBTag_Expected             = CutFlowNumberWithWeights[13]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_noErr_0bTag          = CutFlowNumberWithWeights[14]/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_noErr_1bTag          = CutFlowNumberWithWeights[15]/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_pErr_0bTag           = CutFlowNumberWithWeights[16]/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_pErr_1bTag           = CutFlowNumberWithWeights[17]/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_noErr_0bTag_OldRec7T = CutFlowNum_OldRec7T_passingCSVL_noErr_0bTag/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_noErr_1bTag_OldRec7T = CutFlowNum_OldRec7T_passingCSVL_noErr_1bTag/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_pErr_0bTag_OldRec7T  = CutFlowNum_OldRec7T_passingCSVL_pErr_0bTag/CutFlowNumberWithWeights[13];
   Eff_PassCSVLBTag_pErr_1bTag_OldRec7T  = CutFlowNum_OldRec7T_passingCSVL_pErr_1bTag/CutFlowNumberWithWeights[13];
   Eff_PassCSVMBTag_Expected             = CutFlowNumberWithWeights[18]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_noErr_0bTag          = CutFlowNumberWithWeights[19]/CutFlowNumberWithWeights[18];
   Eff_PassCSVMBTag_noErr_1bTag          = CutFlowNumberWithWeights[20]/CutFlowNumberWithWeights[18];
   Eff_PassCSVMBTag_pErr_0bTag           = CutFlowNumberWithWeights[21]/CutFlowNumberWithWeights[18];
   Eff_PassCSVMBTag_pErr_1bTag           = CutFlowNumberWithWeights[22]/CutFlowNumberWithWeights[18];
   Eff_PassCSVTBTag_Expected             = CutFlowNumberWithWeights[23]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_noErr_0bTag          = CutFlowNumberWithWeights[24]/CutFlowNumberWithWeights[23];
   Eff_PassCSVTBTag_noErr_1bTag          = CutFlowNumberWithWeights[25]/CutFlowNumberWithWeights[23];
   Eff_PassCSVTBTag_pErr_0bTag           = CutFlowNumberWithWeights[26]/CutFlowNumberWithWeights[23];
   Eff_PassCSVTBTag_pErr_1bTag           = CutFlowNumberWithWeights[27]/CutFlowNumberWithWeights[23];
   Eff_TrueBJets                         = CutFlowNumberWithWeights[28]/CutFlowNumberWithWeights[12];
   Eff_TrueBJetsPassingCSVL              = CutFlowNumberWithWeights[29]/CutFlowNumberWithWeights[28];
   Eff_TrueBJetsPassingCSVM              = CutFlowNumberWithWeights[30]/CutFlowNumberWithWeights[28];
   Eff_TrueBJetsPassingCSVT              = CutFlowNumberWithWeights[31]/CutFlowNumberWithWeights[28];
   Eff_NonBJets                          = CutFlowNumberWithWeights[32]/CutFlowNumberWithWeights[12];
   Eff_NonBJetsPassingCSVL               = CutFlowNumberWithWeights[33]/CutFlowNumberWithWeights[32];
   Eff_NonBJetsPassingCSVM               = CutFlowNumberWithWeights[34]/CutFlowNumberWithWeights[32];
   Eff_NonBJetsPassingCSVT               = CutFlowNumberWithWeights[35]/CutFlowNumberWithWeights[32];
   Eff_NodetadphiMasscut_PassCSVL_Expected   = ncsvlbjets/njets;
   Eff_NodetadphiMasscut_PassCSVM_Expected   = ncsvmbjets/njets;
   Eff_NodetadphiMasscut_PassCSVT_Expected   = ncsvtbjets/njets;
   Eff_LICTD                             = CutFlowNumberWithWeights[12]/CutFlowNumberWithWeights[36];

   EffError_PassCSVL_0bTag               = fabs(Eff_PassCSVLBTag_pErr_0bTag - Eff_PassCSVLBTag_noErr_0bTag);
   EffError_PassCSVL_1bTag               = fabs(Eff_PassCSVLBTag_pErr_1bTag - Eff_PassCSVLBTag_noErr_1bTag);
   EffError_PassCSVL_0bTag_OldRec7T      = fabs(Eff_PassCSVLBTag_pErr_0bTag_OldRec7T - Eff_PassCSVLBTag_noErr_0bTag_OldRec7T);
   EffError_PassCSVL_1bTag_OldRec7T      = fabs(Eff_PassCSVLBTag_pErr_1bTag_OldRec7T - Eff_PassCSVLBTag_noErr_1bTag_OldRec7T);
   EffError_PassCSVM_0bTag               = fabs(Eff_PassCSVMBTag_pErr_0bTag - Eff_PassCSVMBTag_noErr_0bTag);
   EffError_PassCSVM_1bTag               = fabs(Eff_PassCSVMBTag_pErr_1bTag - Eff_PassCSVMBTag_noErr_1bTag);
   EffError_PassCSVT_0bTag               = fabs(Eff_PassCSVTBTag_pErr_0bTag - Eff_PassCSVTBTag_noErr_0bTag);
   EffError_PassCSVT_1bTag               = fabs(Eff_PassCSVTBTag_pErr_1bTag - Eff_PassCSVTBTag_noErr_1bTag);


   //-------------------------------------------------------------------------------------------------------------
   cout << "****************************************************************************************************" << endl;
   cout << " Total no. of events = " << CutFlowNumber[0] << "," << CutFlowNumberWithWeights[0] << endl;
   cout << "****************************************************************************************************" << endl;
   cout << "No. of events passing HLT = " << CutFlowNumber[1] << "," << CutFlowNumberWithWeights[1] << endl;
   cout << "No. of events passing PrimaryVtx = " << CutFlowNumber[2] << "," << CutFlowNumberWithWeights[2] << endl;
   cout << "No. of events passing Scraping = " << CutFlowNumber[3] << "," << CutFlowNumberWithWeights[3] << endl;
   cout << "No. of events passing PhotonID = " << CutFlowNumber[4] << "," << CutFlowNumberWithWeights[4] << endl;
   cout << "No. of events passing PhotonPt = " << CutFlowNumber[5] << "," << CutFlowNumberWithWeights[5] << endl;
   cout << "No. of events passing PhotonEta = " << CutFlowNumber[6] << "," << CutFlowNumberWithWeights[6] << endl;
   cout << "No. of events passing JetID = " << CutFlowNumber[7] << "," << CutFlowNumberWithWeights[7] << endl;
   cout << "No. of events passing JetPt = " << CutFlowNumber[8] << "," << CutFlowNumberWithWeights[8] << endl;
   cout << "No. of events passing JetEta = " << CutFlowNumber[9] << "," << CutFlowNumberWithWeights[9] << endl;
   cout << "No. of events passing DPhiCut = " << CutFlowNumber[10] << "," << CutFlowNumberWithWeights[10] << endl;
   cout << "No. of events passing DetaCut = " << CutFlowNumber[11] << "," << CutFlowNumberWithWeights[11] << endl;
   cout << "No. of events passing InvtMass Cut = " << CutFlowNumber[12] << "," << CutFlowNumberWithWeights[12] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVL BTag = " << CutFlowNumber[13] << "," << CutFlowNumberWithWeights[13] << endl;
   cout << "No. of events passing CSVL BTag with noErr_0bTag weight = " << CutFlowNumberWithWeights[14] << endl;
   cout << "No. of events passing CSVL BTag with noErr_1bTag weight = " << CutFlowNumberWithWeights[15] << endl;
   cout << "No. of events passing CSVL BTag with pErr_0bTag weight = " << CutFlowNumberWithWeights[16] << endl;
   cout << "No. of events passing CSVL BTag with pErr_1bTag weight = " << CutFlowNumberWithWeights[17] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVL BTag with noErr_0bTag_OldRec7T weight = " << CutFlowNum_OldRec7T_passingCSVL_noErr_0bTag << endl;
   cout << "No. of events passing CSVL BTag with noErr_1bTag_OldRec7T weight = " << CutFlowNum_OldRec7T_passingCSVL_noErr_1bTag << endl;
   cout << "No. of events passing CSVL BTag with pErr_0bTag_OldRec7T weight = " << CutFlowNum_OldRec7T_passingCSVL_pErr_0bTag << endl;
   cout << "No. of events passing CSVL BTag with pErr_1bTag_OldRec7T weight = " << CutFlowNum_OldRec7T_passingCSVL_pErr_1bTag << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVM BTag = " << CutFlowNumber[14] << "," << CutFlowNumberWithWeights[18] << endl;
   cout << "No. of events passing CSVM BTag with noErr_0bTag weight = " << CutFlowNumberWithWeights[19] << endl;
   cout << "No. of events passing CSVM BTag with noErr_1bTag weight = " << CutFlowNumberWithWeights[20] << endl;
   cout << "No. of events passing CSVM BTag with pErr_0bTag weight = " << CutFlowNumberWithWeights[21] << endl;
   cout << "No. of events passing CSVM BTag with pErr_1bTag weight = " << CutFlowNumberWithWeights[22] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVT BTag = " << CutFlowNumber[15] << "," << CutFlowNumberWithWeights[23] << endl;
   cout << "No. of events passing CSVT BTag with noErr_0bTag weight = " << CutFlowNumberWithWeights[24] << endl;
   cout << "No. of events passing CSVT BTag with noErr_1bTag weight = " << CutFlowNumberWithWeights[25] << endl;
   cout << "No. of events passing CSVT BTag with pErr_0bTag weight = " << CutFlowNumberWithWeights[26] << endl;
   cout << "No. of events passing CSVT BTag with pErr_1bTag weight = " << CutFlowNumberWithWeights[27] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of True BJets = " << CutFlowNumber[16] << "," << CutFlowNumberWithWeights[28] << endl;
   cout << "No. of True BJets passing CSVL = " << CutFlowNumber[17] << "," << CutFlowNumberWithWeights[29] << endl;
   cout << "No. of True BJets passing CSVM = " << CutFlowNumber[18] << "," << CutFlowNumberWithWeights[30] << endl;
   cout << "No. of True BJets passing CSVT = " << CutFlowNumber[19] << "," << CutFlowNumberWithWeights[31] << endl;
   cout << "No.of Non BJets =" << CutFlowNumber[20] << "," << CutFlowNumberWithWeights[32] << endl;
   cout << "No.of Non BJets passing CSVL = " << CutFlowNumber[21] << "," << CutFlowNumberWithWeights[33] <<endl;
   cout << "No. of Non BJets passing CSVM = " << CutFlowNumber[22] << "," << CutFlowNumberWithWeights[34] << endl;
   cout << "No. of Non BJets passing CSVT = " << CutFlowNumber[23] << "," << CutFlowNumberWithWeights[35] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No.of events having leading jet without deta, dphi and mass cut = " << njets << endl;
   cout << "No. of events passing CSVL BTag without deta, dphi and mass cut = " << ncsvlbjets << endl;
   cout << "No. of events passing CSVM BTag without deta, dphi and mass cut = " << ncsvmbjets << endl;
   cout << "No. of events passing CSVT BTag without deta, dphi and mass cut = " << ncsvtbjets << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing with LICTD cut = " << CutFlowNumberWithWeights[12] << endl;
   cout << "No. of events passing without LICTD cut without weight = " << CutFlowNumber[24] << endl;
   cout << "No. of events passing without LICTD cut with weight = " << CutFlowNumberWithWeights[36] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "*********************************Efficiency of various cuts******************************************" << endl;
   cout << "Eff_PassHLT         = " << Eff_PassHLT*100 << "%" << endl;
   cout << "Eff_PassScraping    = " << Eff_PassScraping*100 << "%" << endl;
   cout << "Eff_PassPrimaryVtx  = " << Eff_PassPrimaryVtx*100 << "%" << endl;
   cout << "Eff_PassPhotonID    = " << Eff_PassPhotonID*100 << "%" << endl;
   cout << "Eff_PassPhotonPt    = " << Eff_PassPhotonPt*100 << "%" << endl;
   cout << "Eff_PassPhotonEta   = " << Eff_PassPhotonEta*100 << "%" << endl;
   cout << "Eff_PassJetID       = " << Eff_PassJetID*100 << "%" << endl;
   cout << "Eff_PassJetPt       = " << Eff_PassJetPt*100 << "%" << endl;
   cout << "Eff_PassJetEta      = " << Eff_PassJetEta*100 << "%" << endl;
   cout << "Eff_PassDPhi        = " << Eff_PassDPhi*100 << "%" << endl;
   cout << "Eff_PassDEta        = " << Eff_PassDEta*100 << "%" << endl;
   cout << "Eff_PassGJInvtMass  = " << Eff_PassGJInvtMass*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_Expected                                             = " << Eff_PassCSVLBTag_Expected*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_noErr_0bTag(CSVL_MistagRate_noErr)                   = " << Eff_PassCSVLBTag_noErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_noErr_1bTag(CSVL_TrueEff_noErr)                      = " << Eff_PassCSVLBTag_noErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_pErr_0bTag(CSVL_MistagRate_pErr)                     = " << Eff_PassCSVLBTag_pErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_pErr_1bTag(CSVL_TrueEff_pErr)                        = " << Eff_PassCSVLBTag_pErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_noErr_0bTag_OldRec7T(CSVL_MistagRate_noErr_OldRec7T) = " << Eff_PassCSVLBTag_noErr_0bTag_OldRec7T*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_noErr_1bTag_OldRec7T(CSVL_TrueEff_noErr_OldRec7T)    = " << Eff_PassCSVLBTag_noErr_1bTag_OldRec7T*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_pErr_0bTag_OldRec7T(CSVL_MistagRate_pErr_OldRec7T)   = " << Eff_PassCSVLBTag_pErr_0bTag_OldRec7T*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_pErr_1bTag_OldRec7T(CSVL_TrueEff_pErr_OldRec7T)      = " << Eff_PassCSVLBTag_pErr_1bTag_OldRec7T*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_Expected                                             = " << Eff_PassCSVMBTag_Expected*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_noErr_0bTag(CSVM_MistagRate_noErr)                   = " << Eff_PassCSVMBTag_noErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_noErr_1bTag(CSVM_TrueEff_noErr)                      = " << Eff_PassCSVMBTag_noErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_pErr_0bTag(CSVM_MistagRate_pErr)                     = " << Eff_PassCSVMBTag_pErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_pErr_1bTag(CSVM_TrueEff_pErr)                        = " << Eff_PassCSVMBTag_pErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_Expected                                             = " << Eff_PassCSVTBTag_Expected*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_noErr_0bTag(CSVT_MistagRate_noErr)                   = " << Eff_PassCSVTBTag_noErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_noErr_1bTag(CSVT_TrueEff_noErr)                      = " << Eff_PassCSVTBTag_noErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_pErr_0bTag(CSVT_MistagRate_pErr)                     = " << Eff_PassCSVTBTag_pErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_pErr_1bTag(CSVT_TrueEff_pErr)                        = " << Eff_PassCSVTBTag_pErr_1bTag*100 << "%" << endl;
   cout << "Eff_TrueBJets                                                         = " << Eff_TrueBJets*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVL(CSVL_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVL*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVM(CSVM_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVM*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVT(CSVT_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVT*100 << "%" << endl;
   cout << "Eff_NonBJets                                                          = " << Eff_NonBJets*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVL(CSVL_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVL*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVM(CSVM_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVM*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVT(CSVT_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVT*100 << "%" << endl;
   cout << "Eff_NodetadphiMasscut_PassCSVL_Expected                               = " << Eff_NodetadphiMasscut_PassCSVL_Expected*100 << "%" << endl;
   cout << "Eff_NodetadphiMasscut_PassCSVM_Expected                               = " << Eff_NodetadphiMasscut_PassCSVM_Expected*100 << "%" << endl;
   cout << "Eff_NodetadphiMasscut_PassCSVT_Expected                               = " << Eff_NodetadphiMasscut_PassCSVT_Expected*100 << "%" << endl;
   cout << "Eff_LICTD           = " << Eff_LICTD*100 << "%" << endl;
   cout << "******************************************************************************************************" << endl;
   cout << "****************Efficiency and Eff Error Parameters required for Tagging Eff Calculation**************" << endl;
   cout << "******************************************************************************************************" << endl;
   cout << "Eff_PassCSVLBTag_noErr_0bTag          = " << Eff_PassCSVLBTag_noErr_0bTag << endl;
   cout << "EffError_PassCSVL_0bTag               = " << EffError_PassCSVL_0bTag << endl;
   cout << "Eff_PassCSVLBTag_noErr_1bTag          = " << Eff_PassCSVLBTag_noErr_1bTag << endl;
   cout << "EffError_PassCSVL_1bTag               = " << EffError_PassCSVL_1bTag << endl;
   cout << "Eff_PassCSVLBTag_noErr_0bTag_OldRec7T = " << Eff_PassCSVLBTag_noErr_0bTag_OldRec7T << endl;
   cout << "EffError_PassCSVL_0bTag_OldRec7T      = " << EffError_PassCSVL_0bTag_OldRec7T << endl;
   cout << "Eff_PassCSVLBTag_noErr_1bTag_OldRec7T = " << Eff_PassCSVLBTag_noErr_1bTag_OldRec7T << endl;
   cout << "EffError_PassCSVL_1bTag_OldRec7T      = " << EffError_PassCSVL_1bTag_OldRec7T << endl;
   cout << "Eff_PassCSVMBTag_noErr_0bTag          = " << Eff_PassCSVMBTag_noErr_0bTag << endl;
   cout << "EffError_PassCSVM_0bTag               = " << EffError_PassCSVM_0bTag << endl;
   cout << "Eff_PassCSVMBTag_noErr_1bTag          = " << Eff_PassCSVMBTag_noErr_1bTag << endl;
   cout << "EffError_PassCSVM_1bTag               = " << EffError_PassCSVM_1bTag << endl;
   cout << "Eff_PassCSVTBTag_noErr_0bTag          = " << Eff_PassCSVTBTag_noErr_0bTag << endl;
   cout << "EffError_PassCSVT_0bTag               = " << EffError_PassCSVT_0bTag << endl;
   cout << "Eff_PassCSVTBTag_noErr_1bTag          = " << Eff_PassCSVTBTag_noErr_1bTag << endl;
   cout << "EffError_PassCSVT_1bTag               = " << EffError_PassCSVT_1bTag << endl;
   cout << "******************************************************************************************************" << endl;




}















#!/bin/bash

sampleIndex=0

lpcqstarDir=/eos/uscms/store/user/lpcqstar
MCnTuplesDir=/eos/uscms/store/user/lpcqstar/NTuples/MC
GJ=PhotonJet
QCD=QCDDiJet
Bstar=BstarSignal
EWK=ElectroWeak

#This one needs to change to save results in a new directory
PostAna_this=PostAna_Results_MassCut532_jetEta2p0_dphi1p5_deta2p0_PhIdTight

myDir=${lpcqstarDir}/Rocky
PostAna_Dir=${myDir}/PostAnalyzer_Results
PostAna_Dir_final=${PostAna_Dir}/${PostAna_this}
PostAnaMC_Dir=${PostAna_Dir_final}/MC


#GammaJet (one job of 1400to1800 with 8001 so less totalEvents)
for i in G_Pt_120to170 G_Pt_170to300 G_Pt_300to470 G_Pt_470to800 G_Pt_800to1400 G_Pt_1400to1800 G_Pt_1800toInf
do 
XS=( 108.0068   30.12207   2.138632   0.2119244   0.007077847   4.510327E-5   1.867141E-6 )
totalEvents=( 2000043   2000069   2000130   1975231   1973504   1975015   1939122   )

#for i in G_Pt_1800toInf
#do 
#XS=(1.867141E-6)
#totalEvents=(1939122)

sourceDir=${MCnTuplesDir}/${GJ}/${i}/
outputDir=${PostAnaMC_Dir}/${GJ}


#QCDDiJet (five jobs of 120to170 with 8001 so less totalEvents
#for i in QCDDiJet_Pt_120to170 QCDDiJet_Pt_170to300 QCDDiJet_Pt_300to470 QCDDiJet_Pt_470to600 QCDDiJet_Pt_600to800 QCDDiJet_Pt_800to1000 QCDDiJet_Pt_1000to1400 QCDDiJet_Pt_1400to1800 QCDDiJet_Pt_1800toInf
#do
#XS=( 156293.3   34138.15   1759.549   113.8791   26.9921   3.550036   0.737844   0.03352235   0.001829005 )
#totalEvents=( 5836832   5814398   5978500   3994848   3996864   3998563   1964088   2000062   977586 )

#sourceDir=${MCnTuplesDir}/${QCD}/${i}/
#outputDir=${PostAnaMC_Dir}/${QCD}

#Not submitting jobs for this as ntuples has lots of files goes with 8001
#Electroweak
#for i in EWK_WGToLNuG EWK_WJetsToLNu EWK_ZGInclusive
#do
#XS=( 553.9  37510.0  172.1 )
#totalEvents=( 4802358   18393090   6321549  )

#sourceDir=${MCnTuplesDir}/${EWK}/${i}/
#outputDir=${PostAnaMC_Dir}/${EWK}


#BStarSignal (one job f 500_f1p0 with 8001 so less totalEvents)
#for i in BstarToGJ_M500_f0p1 BstarToGJ_M500_f0p5 BstarToGJ_M500_f1p0 BstarToGJ_M700_f0p1 BstarToGJ_M700_f0p5 BstarToGJ_M700_f1p0 BstarToGJ_M1000_f0p1 BstarToGJ_M1000_f0p5 BstarToGJ_M1000_f1p0 BstarToGJ_M1200_f0p1 BstarToGJ_M1200_f0p5 BstarToGJ_M1200_f1p0 BstarToGJ_M1500_f0p1 BstarToGJ_M1500_f0p5 BstarToGJ_M1500_f1p0 BstarToGJ_M1700_f0p1 BstarToGJ_M1700_f0p5 BstarToGJ_M1700_f1p0 BstarToGJ_M2000_f0p1 BstarToGJ_M2000_f0p5 BstarToGJ_M2000_f1p0 BstarToGJ_M2500_f0p1 BstarToGJ_M2500_f0p5 BstarToGJ_M2500_f1p0 BstarToGJ_M3000_f0p1 BstarToGJ_M3000_f0p5 BstarToGJ_M3000_f1p0 BstarToGJ_M4000_f0p1 BstarToGJ_M4000_f0p5 BstarToGJ_M4000_f1p0
#do

#XS=( 0.01377   0.3415   1.348   0.002143   0.05305   0.2096   0.0002366   0.005841   0.02291   0.00006691   0.001667   0.00652   0.00001231   0.0003057   0.001214   0.000004339   0.0001087   0.0004321   0.000001003   0.00002505   0.0001013   0.0000001001   0.000002543   0.00001072   0.00000001103   0.000000291   0.00000134   0.0000000001531   0.000000005162   0.00000003707 )
#totalEvents=( 60000  60000  59700  120000  60000  60000  120000  120000  60000  120000  119626  120000  120000  120000  120000  160000  120000  120000  160000  160000  118402  198980  160000  160000  197872  198200  160000  196962  197291  198922 )

#sourceDir=${MCnTuplesDir}/${Bstar}/${i}/
#outputDir=${PostAnaMC_Dir}/${Bstar}

#This is required for calculating systematics for signal only
if [ ${sourceDir} = "${MCnTuplesDir}/${Bstar}/${i}/" ]
then
Signal=true
else
Signal=false
fi

if [ -f PostAnalyzer_MC.C ]; then
echo "++++++++++++++ Deleting PostAnalyzer_MC.C ++++++++++++++"
rm PostAnalyzer_MC.C
fi

if [ -f PostAnalyzer_MC.h ]; then
echo "++++++++++++++ Deleting PostAnalyzer_MC.h ++++++++++++++"
rm PostAnalyzer_MC.h
fi

filenameTag=${i}
destinationDir=${outputDir}
#destinationDir=${outputDir}/Check_GJ #used this just to check save files in another directory, dont need it for normal job submission

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

if [ ! -d ${PostAnaMC_Dir} ]; then
echo "Making Directory ${PostAnaMC_Dir}"
mkdir ${PostAnaMC_Dir}
chmod 775 ${PostAnaMC_Dir}
fi

if [ ! -d ${destinationDir} ]; then
echo "Making Directory ${destinationDir}"
mkdir ${destinationDir}
chmod 775 ${destinationDir}
fi

echo "filenameTag = ${filenameTag}"
echo "destinationDir = ${destinationDir}"
echo "sourceDir = ${sourceDir}"


cat >> PostAnalyzer_MC.C <<EOF 
#define PostAnalyzer_MC_cxx
#include "PostAnalyzer_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

   Cut_Vtx_z = 24.0;
   Cut_Vtx_ndof = 4.0;
   Cut_Vtx_rho = 2.0;

   Cut_Photon_pt = 170.0;
   Cut_Photon_eta = 1.4442;
   Cut_Jet_pt = 170.0;
   Cut_Jet_eta = 2.0;
   Cut_GJdPhi = 1.5;
   Cut_GJdEta = 2.0;
   Cut_GJInvtMass = 532.0; 

   //Photon energy correction uncertainty (provided centerally)
   Photon_pecUnc = 0.015;

   //********************************************************************************************************//

   //Make it true if mass cut to be applied otherwise false
   IfMassCut = true;

   //Make it true if running on bstar signal samples for pileup systematics calculation
   //Uncomment it in script
   IfSignal = ${Signal};

   //Make true for which Photon Id to apply
   IfLoosePhId = false;
   IfMediumPhId = false;
   IfTightPhId = true;

   //*********************************************************************************************************//

   //To be removed in script 
   //-----------------------
   //Define Output file here
   //file = new TFile("PostAnalyzer_MC.root", "RECREATE");
   //-----------------------

   //Uncomment this in script
   //Define Output file here
   TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   file = new TFile(OutputPath+OutputFile+".root", "RECREATE");

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   PileupReWeighting();
   PileupReWeighting_XSm5();
   PileupReWeighting_XSp5();

   //*******************************************************************************************************//
   //Defining the histogram from filling number of events after various cuts

   const int nbins = 25;
   const int nbinsWt = 43;
   TString CutFlowLabels[nbins] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag", "PassCSVMBTag", "PassCSVTBTag", "TrueBJets", "TrueBJetsPassingCSVL", "TrueBJetsPassingCSVM", "TrueBJetsPassingCSVT", "NonBJets", "NonBJetsPassingCSVL", "NonBJetsPassingCSVM", "NonBJetsPassingCSVT", "PassPhotonID_NoLICTD"};

   TString CutFlowLabelsWithWts[nbinsWt] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag_Expected", "PassCSVLBTag_noErr_0bTag", "PassCSVLBTag_noErr_1bTag", "PassCSVLBTag_pErr_0bTag", "PassCSVLBTag_pErr_1bTag", "PassCSVLBTag_nErr_0bTag", "PassCSVLBTag_nErr_1bTag", "PassCSVMBTag_Expected", "PassCSVMBTag_noErr_0bTag", "PassCSVMBTag_noErr_1bTag", "PassCSVMBTag_pErr_0bTag", "PassCSVMBTag_pErr_1bTag", "PassCSVMBTag_nErr_0bTag", "PassCSVMBTag_nErr_1bTag", "PassCSVTBTag_Expected", "PassCSVTBTag_noErr_0bTag", "PassCSVTBTag_noErr_1bTag", "PassCSVTBTag_pErr_0bTag", "PassCSVTBTag_pErr_1bTag", "PassCSVTBTag_nErr_0bTag", "PassCSVTBTag_nErr_1bTag", "TrueBJets", "TrueBJetsPassingCSVL", "TrueBJetsPassingCSVM", "TrueBJetsPassingCSVT", "NonBJets", "NonBJetsPassingCSVL", "NonBJetsPassingCSVM", "NonBJetsPassingCSVT", "PassPhotonID_NoLICTD"};

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
    
   //--------------------------------------------------------------------------------------------------
   //This is required for pileup systematics calculations to get the number of events for 5% more and less PU reweighting.
   const int nbn = 15;

   TString CutFlowLabelsWithWts_XSm5[nbn] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag_noErr_0bTag", "PassCSVLBTag_noErr_1bTag"};
   TString CutFlowLabelsWithWts_XSp5[nbn] = {"Total", "PassHLT", "PassScraping", "PassPrimaryVtx", "PassPhotonID", "PassPhotonPt", "PassPhotonEta", "PassJetID", "PassJetPt", "PassJetEta", "PassDPhi", "PassDEta", "PassGJInvtMass", "PassCSVLBTag_noErr_0bTag", "PassCSVLBTag_noErr_1bTag"};

   h_CutFlowTableWithWeights_XSm5 = new TH1F("h_CutFlowTableWithWeights_XSm5", "Events Passing Various Cuts With Weights but with 5% less XS in pileup reweighting", nbn, 0, nbn);
   h_CutFlowTableWithWeights_XSm5->GetYaxis()->SetTitle("Events"); h_CutFlowTableWithWeights_XSm5->GetYaxis()->CenterTitle();
   h_CutFlowTableWithWeights_XSp5 = new TH1F("h_CutFlowTableWithWeights_XSp5", "Events Passing Various Cuts With Weights but with 5% more XS in pileup reweighting", nbn, 0, nbn);
   h_CutFlowTableWithWeights_XSp5->GetYaxis()->SetTitle("Events"); h_CutFlowTableWithWeights_XSp5->GetYaxis()->CenterTitle();

   for(int i = 0; i < nbn; i++){
     h_CutFlowTableWithWeights_XSm5->GetXaxis()->SetBinLabel(i+1, CutFlowLabelsWithWts_XSm5[i]);
   }
   for(int i = 0; i < nbn; i++){
     h_CutFlowTableWithWeights_XSp5->GetXaxis()->SetBinLabel(i+1, CutFlowLabelsWithWts_XSp5[i]);
   }

   Double_t CutFlowNumberWithWeights_XSm5[nbn];
   Double_t CutFlowNumberWithWeights_XSp5[nbn];
   for(int i = 0; i < nbn; i++){
     CutFlowNumberWithWeights_XSm5[i] = 0.0;
   }
   for(int i = 0; i < nbn; i++){
     CutFlowNumberWithWeights_XSp5[i] = 0.0;
   }
   //-----------------------------------------------------------------------------------------------------

   Long64_t nentries = fChain->GetEntries();
   cout << "no. of entries " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;

   for(Long64_t jentry = 0; jentry < nentries; jentry++){
    
     //     cout << "++++++++++++++++++Analyzing entry++++++++++++" << jentry << endl;

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

     //Uncomment this in script
     Lumi_EvtWt = Lumi*(${XS[${sampleIndex}]}/${totalEvents[${sampleIndex}]});

     PU_EvtWt = PUWeights(trueNumofInteractions);
     PreBTag_EvtWt = Lumi_EvtWt * PU_EvtWt;

     h_trueNumofInteractions->Fill(trueNumofInteractions);

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

     h_goodPV->Fill(GoodVertex);
     h_goodPV_LumiWt->Fill(GoodVertex, Lumi_EvtWt);
     h_goodPV_PUWt->Fill(GoodVertex, PU_EvtWt);
     h_goodPV_LumiNPUWt->Fill(GoodVertex, PreBTag_EvtWt);

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


     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
			     h_GJetInvtMass_UnitBinning->Fill(GetInvtMass(PC, JC), PreBTag_EvtWt);
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

			       double SF_CSVL, SFperr_CSVL, SFnerr_CSVL;

			       if((*PFPatJet_JetPartonFlavor)[JC] == 5){ //b and bbar 

				 SF_CSVL = GetSFb_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVL = GetSFberr_CSVL((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);
				 SFnerr_CSVL = GetSFberr_CSVL((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       }			    
		     
			       if((*PFPatJet_JetPartonFlavor)[JC] == 4){ // c and cbar

				 SF_CSVL = GetSFc_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVL = GetSFcerr_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVL = GetSFcerr_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
	
			       }
			     
			       if((*PFPatJet_JetPartonFlavor)[JC] == 1 || (*PFPatJet_JetPartonFlavor)[JC] == 2 || (*PFPatJet_JetPartonFlavor)[JC] == 3 || (*PFPatJet_JetPartonFlavor)[JC] == 21){ //u, d, s and g 

				 SF_CSVL = GetSFLight("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVL = GetSFLight_perr("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVL = GetSFLight_nerr("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);

 			       }

			       Total_CSVLEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL, 0);
			       Total_CSVLEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL, 1);
			       Total_CSVLEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL+SFperr_CSVL, 0);
			       Total_CSVLEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL+SFperr_CSVL, 1);
			       Total_CSVLEvtWt_nErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL-SFnerr_CSVL, 0);
                               Total_CSVLEvtWt_nErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVL-SFnerr_CSVL, 1);

			       CutFlowNumberWithWeights[14] += Total_CSVLEvtWt_noErr_0bTag;
			       CutFlowNumberWithWeights[15] += Total_CSVLEvtWt_noErr_1bTag;
			       CutFlowNumberWithWeights[16] += Total_CSVLEvtWt_pErr_0bTag;
			       CutFlowNumberWithWeights[17] += Total_CSVLEvtWt_pErr_1bTag;
			       CutFlowNumberWithWeights[18] += Total_CSVLEvtWt_nErr_0bTag;
			       CutFlowNumberWithWeights[19] += Total_CSVLEvtWt_nErr_1bTag;
			    
			       h_CSVL_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_noErr_0bTag);

			       h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_noErr_0bTag);

			       h_CSVL_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_noErr_1bTag);

			       h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_noErr_1bTag);

			       h_CSVL_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_1bTag);
			       h_CSVL_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       h_CSVL_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_noErr_1bTag);
			       h_CSVL_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_pErr_1bTag);
			     
			     }
			     else{

			       Total_CSVLEvtWt_noErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVLEvtWt_pErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVLEvtWt_nErr_0bTag = PreBTag_EvtWt * 1;

                               CutFlowNumberWithWeights[14] += Total_CSVLEvtWt_noErr_0bTag;
                               CutFlowNumberWithWeights[16] += Total_CSVLEvtWt_pErr_0bTag;
                               CutFlowNumberWithWeights[18] += Total_CSVLEvtWt_nErr_0bTag;
			     
			       h_CSVL_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_noErr_0bTag);
                               h_CSVL_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_noErr_0bTag);
                               h_CSVL_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_noErr_0bTag);

                               h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
                               h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
                               h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_noErr_0bTag);
                               h_CSVL_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);
			       h_CSVL_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_noErr_0bTag);

			       h_CSVL_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);
			       h_CSVL_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVLEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVLEvtWt_pErr_0bTag);

			       h_CSVL_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_noErr_0bTag);		    
			       h_CSVL_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC],Total_CSVLEvtWt_pErr_0bTag);

			     }
			   
			     if(Pass_CSVMBTag){
			       CutFlowNumber[14]++;
			       CutFlowNumberWithWeights[20] += PreBTag_EvtWt;

			       double SF_CSVM, SFperr_CSVM, SFnerr_CSVM;

			       if((*PFPatJet_JetPartonFlavor)[JC] == 5){ //b and bbar

				 SF_CSVM = GetSFb_CSVM((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVM = GetSFberr_CSVM((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);
				 SFnerr_CSVM = GetSFberr_CSVM((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       }

			       if((*PFPatJet_JetPartonFlavor)[JC] == 4){ // c and cbar

				 SF_CSVM = GetSFc_CSVM((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVM = GetSFcerr_CSVM((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVM = GetSFcerr_CSVM((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
	
			       }

			       if((*PFPatJet_JetPartonFlavor)[JC] == 1 || (*PFPatJet_JetPartonFlavor)[JC] == 2 || (*PFPatJet_JetPartonFlavor)[JC] == 3 || (*PFPatJet_JetPartonFlavor)[JC] == 21){ //u, d, s and g 

				 SF_CSVM = GetSFLight("CSV", "M", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVM = GetSFLight_perr("CSV", "M", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVM = GetSFLight_nerr("CSV", "M", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);

 			       }
			     
			       Total_CSVMEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM, 0);
			       Total_CSVMEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM, 1);
			       Total_CSVMEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM+SFperr_CSVM, 0);
			       Total_CSVMEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM+SFperr_CSVM, 1);
			       Total_CSVMEvtWt_nErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM-SFnerr_CSVM, 0);
                               Total_CSVMEvtWt_nErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVM-SFnerr_CSVM, 1);

			       CutFlowNumberWithWeights[21] += Total_CSVMEvtWt_noErr_0bTag;
			       CutFlowNumberWithWeights[22] += Total_CSVMEvtWt_noErr_1bTag;
			       CutFlowNumberWithWeights[23] += Total_CSVMEvtWt_pErr_0bTag;
			       CutFlowNumberWithWeights[24] += Total_CSVMEvtWt_pErr_1bTag;
			       CutFlowNumberWithWeights[25] += Total_CSVMEvtWt_nErr_0bTag;
			       CutFlowNumberWithWeights[26] += Total_CSVMEvtWt_nErr_1bTag; 

			       h_CSVM_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_noErr_0bTag);
			     
			       h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_noErr_0bTag);

			       h_CSVM_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_noErr_1bTag);

			       h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_noErr_1bTag);

			       h_CSVM_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_pErr_0bTag);

			       h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_pErr_0bTag);

			       h_CSVM_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_1bTag);
			       h_CSVM_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       h_CSVM_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_noErr_1bTag);
			       h_CSVM_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_pErr_1bTag);

			     
			     }
			     else{

			       Total_CSVMEvtWt_noErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVMEvtWt_pErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVMEvtWt_nErr_0bTag = PreBTag_EvtWt * 1;

                               CutFlowNumberWithWeights[21] += Total_CSVMEvtWt_noErr_0bTag;
                               CutFlowNumberWithWeights[23] += Total_CSVMEvtWt_pErr_0bTag;
                               CutFlowNumberWithWeights[25] += Total_CSVMEvtWt_nErr_0bTag;

			       h_CSVM_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_noErr_0bTag);

			       h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);
			       h_CSVM_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_noErr_0bTag);

			       h_CSVM_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVMEvtWt_pErr_0bTag);

			       h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);
			       h_CSVM_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVMEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVMEvtWt_pErr_0bTag);

                               h_CSVM_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_noErr_0bTag);
                               h_CSVM_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVMEvtWt_pErr_0bTag);

			     }
			   
			     if(Pass_CSVTBTag){
			       CutFlowNumber[15]++;
			       CutFlowNumberWithWeights[27] += PreBTag_EvtWt;
			     
			       double SF_CSVT, SFperr_CSVT, SFnerr_CSVT;

			       if((*PFPatJet_JetPartonFlavor)[JC] == 5){ //b and bbar 

				 SF_CSVT = GetSFb_CSVT((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVT = GetSFberr_CSVT((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);
				 SFnerr_CSVT = GetSFberr_CSVT((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

			       }
			     
			       if((*PFPatJet_JetPartonFlavor)[JC] == 4){ // c and cbar

				 SF_CSVT = GetSFc_CSVT((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVT = GetSFcerr_CSVT((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVT = GetSFcerr_CSVT((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
	
			       }
			     
			       if((*PFPatJet_JetPartonFlavor)[JC] == 1 || (*PFPatJet_JetPartonFlavor)[JC] == 2 || (*PFPatJet_JetPartonFlavor)[JC] == 3 || (*PFPatJet_JetPartonFlavor)[JC] == 21){ //u, d, s and g 

				 SF_CSVT = GetSFLight("CSV", "T", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVT = GetSFLight_perr("CSV", "T", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVT = GetSFLight_nerr("CSV", "T", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);

 			       }
			     
			       Total_CSVTEvtWt_noErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT, 0);
			       Total_CSVTEvtWt_noErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT, 1);
			       Total_CSVTEvtWt_pErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT+SFperr_CSVT, 0);
			       Total_CSVTEvtWt_pErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT+SFperr_CSVT, 1);
			       Total_CSVTEvtWt_nErr_0bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT-SFnerr_CSVT, 0);
                               Total_CSVTEvtWt_nErr_1bTag = PreBTag_EvtWt * BTagEventWeight(SF_CSVT-SFnerr_CSVT, 1);			     

			       CutFlowNumberWithWeights[28] += Total_CSVTEvtWt_noErr_0bTag;			     
			       CutFlowNumberWithWeights[29] += Total_CSVTEvtWt_noErr_1bTag;			     
			       CutFlowNumberWithWeights[30] += Total_CSVTEvtWt_pErr_0bTag;			    			     
			       CutFlowNumberWithWeights[31] += Total_CSVTEvtWt_pErr_1bTag;
			       CutFlowNumberWithWeights[32] += Total_CSVTEvtWt_nErr_0bTag;
			       CutFlowNumberWithWeights[33] += Total_CSVTEvtWt_nErr_1bTag;
			     
			       h_CSVT_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_BJetPt_noErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJetEta_noErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJetPhi_noErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_noErr_1bTag);

			       h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdEta_noErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdPhi_noErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_GBJetdR_noErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_noErr_1bTag);

			       h_CSVT_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_pErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_pErr_0bTag);

			       h_CSVT_BJetPt_pErr_1bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_BJetEta_pErr_1bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_BJetPhi_pErr_1bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdEta_pErr_1bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdPhi_pErr_1bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_1bTag);
			       h_CSVT_GBJetdR_pErr_1bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_1bTag);

			       h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       h_CSVT_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJet_PartonFlavor_noErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_noErr_1bTag);
			       h_CSVT_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJet_PartonFlavor_pErr_1bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_pErr_1bTag);
			     
			     }
			     else{

			       Total_CSVTEvtWt_noErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVTEvtWt_pErr_0bTag = PreBTag_EvtWt * 1;
			       Total_CSVTEvtWt_nErr_0bTag = PreBTag_EvtWt * 1;

			       CutFlowNumberWithWeights[28] += Total_CSVTEvtWt_noErr_0bTag;
                               CutFlowNumberWithWeights[30] += Total_CSVTEvtWt_pErr_0bTag;
                               CutFlowNumberWithWeights[32] += Total_CSVTEvtWt_nErr_0bTag;

			       h_CSVT_BJetPt_noErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetEta_noErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_BJetPhi_noErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdEta_noErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdPhi_noErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);
			       h_CSVT_GBJetdR_noErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_noErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_noErr_0bTag);

			       h_CSVT_BJetPt_pErr_0bTag->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetEta_pErr_0bTag->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_BJetPhi_pErr_0bTag->Fill((*PFPatJet_phi)[JC], Total_CSVTEvtWt_pErr_0bTag);

			       h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->Fill(GetInvtMass(PC, JC), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdEta_pErr_0bTag->Fill(GetdEta((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdPhi_pErr_0bTag->Fill(GetdPhi((*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);
			       h_CSVT_GBJetdR_pErr_0bTag->Fill(GetdR((*Photon_SC_eta)[PC], (*PFPatJet_eta)[JC], (*Photon_phi)[PC], (*PFPatJet_phi)[JC]), Total_CSVTEvtWt_pErr_0bTag);

			       h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->Fill((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[JC], Total_CSVTEvtWt_pErr_0bTag);

                               h_CSVT_BJet_PartonFlavor_noErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_noErr_0bTag);
                               h_CSVT_BJet_PartonFlavor_pErr_0bTag->Fill((*PFPatJet_JetPartonFlavor)[JC], Total_CSVTEvtWt_pErr_0bTag);
			     
			     }//END_OF else(Pass_CSVTBTag)

			     //----------------------------------------------------------------------------------
			     //This part to get true b tag efficiency and Mistag rate of the three discriminators
			     //----------------------------------------------------------------------------------
			     if((*PFPatJet_JetPartonFlavor)[JC] == 5){
			       CutFlowNumber[16]++;
			       CutFlowNumberWithWeights[34] += PreBTag_EvtWt;

			       h_TrueBJetPt->Fill((*PFPatJet_pt)[JC], PreBTag_EvtWt);
			       h_TrueBJetEta->Fill((*PFPatJet_eta)[JC], PreBTag_EvtWt);

			       if(Pass_CSVLBTag){
				 CutFlowNumber[17]++;
				 CutFlowNumberWithWeights[35] += Total_CSVLEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVL->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVL->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVMBTag){
				 CutFlowNumber[18]++;
				 CutFlowNumberWithWeights[36] += Total_CSVMEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVM->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVM->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVTBTag){
				 CutFlowNumber[19]++;
				 CutFlowNumberWithWeights[37] += Total_CSVTEvtWt_pErr_1bTag;

				 h_TrueBJetPtPassingCSVT->Fill((*PFPatJet_pt)[JC], Total_CSVTEvtWt_pErr_1bTag);
				 h_TrueBJetEtaPassingCSVT->Fill((*PFPatJet_eta)[JC], Total_CSVTEvtWt_pErr_1bTag);

			       }

			     }else{
			       CutFlowNumber[20]++;
			       CutFlowNumberWithWeights[38] += PreBTag_EvtWt;

			       h_NonBJetPt->Fill((*PFPatJet_pt)[JC], PreBTag_EvtWt);
			       h_NonBJetEta->Fill((*PFPatJet_eta)[JC], PreBTag_EvtWt);

			       if(Pass_CSVLBTag){
				 CutFlowNumber[21]++;
				 CutFlowNumberWithWeights[39] += Total_CSVLEvtWt_pErr_1bTag;

				 h_NonBJetPtPassingCSVL->Fill((*PFPatJet_pt)[JC], Total_CSVLEvtWt_pErr_1bTag);
				 h_NonBJetEtaPassingCSVL->Fill((*PFPatJet_eta)[JC], Total_CSVLEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVMBTag){
				 CutFlowNumber[22]++;
				 CutFlowNumberWithWeights[40] += Total_CSVMEvtWt_pErr_1bTag;

				 h_NonBJetPtPassingCSVM->Fill((*PFPatJet_pt)[JC], Total_CSVMEvtWt_pErr_1bTag);
				 h_NonBJetEtaPassingCSVM->Fill((*PFPatJet_eta)[JC], Total_CSVMEvtWt_pErr_1bTag);

			       }
			       if(Pass_CSVTBTag){
				 CutFlowNumber[23]++;
				 CutFlowNumberWithWeights[41] += Total_CSVTEvtWt_pErr_1bTag;

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
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //---------------------------------------------------------------------------
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
     h_nCSVLBJets->Fill(nCSVLBJets, PreBTag_EvtWt);
     h_nCSVMBJets->Fill(nCSVMBJets, PreBTag_EvtWt);
     h_nCSVTBJets->Fill(nCSVTBJets, PreBTag_EvtWt);
     h_CSVL_BJetsFrac->Fill(frac_CSVL);
     h_CSVM_BJetsFrac->Fill(frac_CSVM);
     h_CSVT_BJetsFrac->Fill(frac_CSVT);

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
	 h_CSVLBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, PreBTag_EvtWt);
       }
     }
     //Filling histogram for CSVM-BJetIdx vs CSVM-BJetPt
     for(int i = 0; i < PFPatJet_n; i++){
       if(PassCSVMBTag(i) == 1){
	 h_CSVMBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, PreBTag_EvtWt);
       }
     }
     //Filling histogram for CSVT-BJetIdx vs CSVT-BJetPt
     for(int i = 0; i < PFPatJet_n; i++){
       if(PassCSVTBTag(i) == 1){
	 h_CSVTBJetIdxVsPt->Fill((*PFPatJet_pt)[i], i, PreBTag_EvtWt);
       }
     }
     //--------------------------------------------------
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //---------------------------------------------------------------
     //For Pileup Uncertainty, JES and PES calculation for systematics
     //---------------------------------------------------------------
     if(IfSignal){

       //----------------------------
       //For Pileup Syst, JES and PES
       //----------------------------
       PU_EvtWt_XSm5 = PUWeights_XSm5(trueNumofInteractions);
       PreBTag_EvtWt_XSm5 = Lumi_EvtWt * PU_EvtWt_XSm5;

       PU_EvtWt_XSp5 = PUWeights_XSp5(trueNumofInteractions);
       PreBTag_EvtWt_XSp5 = Lumi_EvtWt * PU_EvtWt_XSp5;

       CutFlowNumberWithWeights_XSm5[0] += PreBTag_EvtWt_XSm5;
       CutFlowNumberWithWeights_XSp5[0] += PreBTag_EvtWt_XSp5;

       if(Pass_HLT){
	 CutFlowNumberWithWeights_XSm5[1] += PreBTag_EvtWt_XSm5;
	 CutFlowNumberWithWeights_XSp5[1] += PreBTag_EvtWt_XSp5;

	 if(NoScrapingEvt){
	   CutFlowNumberWithWeights_XSm5[2] += PreBTag_EvtWt_XSm5;
	   CutFlowNumberWithWeights_XSp5[2] += PreBTag_EvtWt_XSp5;

	   if(HasPrimaryVtx){
	     CutFlowNumberWithWeights_XSm5[3] += PreBTag_EvtWt_XSm5;
	     CutFlowNumberWithWeights_XSp5[3] += PreBTag_EvtWt_XSp5;

	     if(PC > -1){
	       CutFlowNumberWithWeights_XSm5[4] += PreBTag_EvtWt_XSm5;
	       CutFlowNumberWithWeights_XSp5[4] += PreBTag_EvtWt_XSp5;

	       if(Pass_PhoPtcut){
		 CutFlowNumberWithWeights_XSm5[5] += PreBTag_EvtWt_XSm5;
		 CutFlowNumberWithWeights_XSp5[5] += PreBTag_EvtWt_XSp5;

		 if(Pass_PhoEtaEBcut){
		   CutFlowNumberWithWeights_XSm5[6] += PreBTag_EvtWt_XSm5;
		   CutFlowNumberWithWeights_XSp5[6] += PreBTag_EvtWt_XSp5;

		   if(JC > -1){
		     CutFlowNumberWithWeights_XSm5[7] += PreBTag_EvtWt_XSm5;
		     CutFlowNumberWithWeights_XSp5[7] += PreBTag_EvtWt_XSp5;

		     if(Pass_JetPtcut){
		       CutFlowNumberWithWeights_XSm5[8] += PreBTag_EvtWt_XSm5;
		       CutFlowNumberWithWeights_XSp5[8] += PreBTag_EvtWt_XSp5;

		       if(Pass_JetEtacut){
			 CutFlowNumberWithWeights_XSm5[9] += PreBTag_EvtWt_XSm5;
			 CutFlowNumberWithWeights_XSp5[9] += PreBTag_EvtWt_XSp5;

			 if(Pass_GJdPhicut){
			   CutFlowNumberWithWeights_XSm5[10] += PreBTag_EvtWt_XSm5;
			   CutFlowNumberWithWeights_XSp5[10] += PreBTag_EvtWt_XSp5;

			   if(Pass_GJdEtacut){
			     CutFlowNumberWithWeights_XSm5[11] += PreBTag_EvtWt_XSm5;
			     CutFlowNumberWithWeights_XSp5[11] += PreBTag_EvtWt_XSp5;

			     if(Pass_GJInvtMasscut){
			       CutFlowNumberWithWeights_XSm5[12] += PreBTag_EvtWt_XSm5;
			       CutFlowNumberWithWeights_XSp5[12] += PreBTag_EvtWt_XSp5;

			       //-----------------------------------------------------------------------------------------------------------
			       //----------------
			       //For JES
			       //----------------
			       PFJet4Vec_JESup.SetPxPyPzE((*PFPatJet_px)[JC]*(1+(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_py)[JC]*(1+(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_pz)[JC]*(1+(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_E)[JC]*(1+(*PFPatJet_jecUncertainity)[JC]));
			       PFJet4Vec_JESdown.SetPxPyPzE((*PFPatJet_px)[JC]*(1-(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_py)[JC]*(1-(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_pz)[JC]*(1-(*PFPatJet_jecUncertainity)[JC]), (*PFPatJet_E)[JC]*(1-(*PFPatJet_jecUncertainity)[JC]));
			       GJ4Vec_JESup.SetPxPyPzE((*Photon_px)[PC]+PFJet4Vec_JESup[0], (*Photon_py)[PC]+PFJet4Vec_JESup[1], (*Photon_pz)[PC]+PFJet4Vec_JESup[2], (*Photon_E)[PC]+PFJet4Vec_JESup[3]);
			       GJ4Vec_JESdown.SetPxPyPzE((*Photon_px)[PC]+PFJet4Vec_JESdown[0], (*Photon_py)[PC]+PFJet4Vec_JESdown[1], (*Photon_pz)[PC]+PFJet4Vec_JESdown[2], (*Photon_E)[PC]+PFJet4Vec_JESdown[3]);
			       Mass_JESup = pow((GJ4Vec_JESup[3]*GJ4Vec_JESup[3] - GJ4Vec_JESup[0]*GJ4Vec_JESup[0] - GJ4Vec_JESup[1]*GJ4Vec_JESup[1] - GJ4Vec_JESup[2]*GJ4Vec_JESup[2]), 0.5);
			       Mass_JESdown = pow((GJ4Vec_JESdown[3]*GJ4Vec_JESdown[3] - GJ4Vec_JESdown[0]*GJ4Vec_JESdown[0] - GJ4Vec_JESdown[1]*GJ4Vec_JESdown[1] - GJ4Vec_JESdown[2]*GJ4Vec_JESdown[2]), 0.5);

			       h_profile_JES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			       h_profile_JES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			       h_mass_JES[0][0]->Fill(Mass_JESup, PreBTag_EvtWt);
			       h_mass_JES[1][0]->Fill(Mass_JESdown, PreBTag_EvtWt);

			       //----------------
			       //For PES
			       //----------------
			       Photon4Vec_PESup.SetPxPyPzE((*Photon_px)[PC]*(1+Photon_pecUnc), (*Photon_py)[PC]*(1+Photon_pecUnc), (*Photon_pz)[PC]*(1+Photon_pecUnc), (*Photon_E)[PC]*(1+Photon_pecUnc));
			       Photon4Vec_PESdown.SetPxPyPzE((*Photon_px)[PC]*(1-Photon_pecUnc), (*Photon_py)[PC]*(1-Photon_pecUnc), (*Photon_pz)[PC]*(1-Photon_pecUnc), (*Photon_E)[PC]*(1-Photon_pecUnc));
			       GJ4Vec_PESup.SetPxPyPzE(Photon4Vec_PESup[0]+(*PFPatJet_px)[JC], Photon4Vec_PESup[1]+(*PFPatJet_py)[JC], Photon4Vec_PESup[2]+(*PFPatJet_pz)[JC], Photon4Vec_PESup[3]+(*PFPatJet_E)[JC]);
			       GJ4Vec_PESdown.SetPxPyPzE(Photon4Vec_PESdown[0]+(*PFPatJet_px)[JC], Photon4Vec_PESdown[1]+(*PFPatJet_py)[JC], Photon4Vec_PESdown[2]+(*PFPatJet_pz)[JC], Photon4Vec_PESdown[3]+(*PFPatJet_E)[JC]);
			       Mass_PESup = pow((GJ4Vec_PESup[3]*GJ4Vec_PESup[3] - GJ4Vec_PESup[0]*GJ4Vec_PESup[0] - GJ4Vec_PESup[1]*GJ4Vec_PESup[1] - GJ4Vec_PESup[2]*GJ4Vec_PESup[2]), 0.5);
			       Mass_PESdown = pow((GJ4Vec_PESdown[3]*GJ4Vec_PESdown[3] - GJ4Vec_PESdown[0]*GJ4Vec_PESdown[0] - GJ4Vec_PESdown[1]*GJ4Vec_PESdown[1] - GJ4Vec_PESdown[2]*GJ4Vec_PESdown[2]), 0.5);

			       h_profile_PES[0][0]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			       h_profile_PES[1][0]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

			       h_mass_PES[0][0]->Fill(Mass_PESup, PreBTag_EvtWt);
			       h_mass_PES[1][0]->Fill(Mass_PESdown, PreBTag_EvtWt);
			       //-------------------------------------------------------------------------------------------------------

			       double TotEvtWt_noErr_0btag_XSm5;
			       double TotEvtWt_noErr_0btag_XSp5;
			       double TotEvtWt_noErr_1btag_XSm5;
			       double TotEvtWt_noErr_1btag_XSp5;

			       if(Pass_CSVLBTag){

				 double SF_CSVL, SFperr_CSVL, SFnerr_CSVL;

				 if((*PFPatJet_JetPartonFlavor)[JC] == 5){ //b and bbar

				   SF_CSVL = GetSFb_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				   SFperr_CSVL = GetSFberr_CSVL((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);
				   SFnerr_CSVL = GetSFberr_CSVL((*PFPatJet_pt)[JC],(*PFPatJet_eta)[JC]);

				 }

				 if((*PFPatJet_JetPartonFlavor)[JC] == 4){ // c and cbar

				   SF_CSVL = GetSFc_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				   SFperr_CSVL = GetSFcerr_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				   SFnerr_CSVL = GetSFcerr_CSVL((*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);

				 }

				 if((*PFPatJet_JetPartonFlavor)[JC] == 1 || (*PFPatJet_JetPartonFlavor)[JC] == 2 || (*PFPatJet_JetPartonFlavor)[JC] == 3 || (*PFPatJet_JetPartonFlavor)[JC] == 21){ //u, d, s and g 

				 SF_CSVL = GetSFLight("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFperr_CSVL = GetSFLight_perr("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);
				 SFnerr_CSVL = GetSFLight_nerr("CSV", "L", (*PFPatJet_pt)[JC], (*PFPatJet_eta)[JC]);

 			       }
				 TotEvtWt_noErr_0btag_XSm5 = PreBTag_EvtWt_XSm5 * BTagEventWeight(SF_CSVL, 0);
				 TotEvtWt_noErr_0btag_XSp5 = PreBTag_EvtWt_XSp5 * BTagEventWeight(SF_CSVL, 0);
				 TotEvtWt_noErr_1btag_XSm5 = PreBTag_EvtWt_XSm5 * BTagEventWeight(SF_CSVL, 1);
				 TotEvtWt_noErr_1btag_XSp5 = PreBTag_EvtWt_XSp5 * BTagEventWeight(SF_CSVL, 1);


				 //-------------------------For JES--------------------------------
				 h_profile_JES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_JESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
				 h_profile_JES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_JESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

				 h_mass_JES[0][1]->Fill(Mass_JESup, PreBTag_EvtWt*BTagEventWeight(SF_CSVL, 1));
				 h_mass_JES[1][1]->Fill(Mass_JESdown, PreBTag_EvtWt*BTagEventWeight(SF_CSVL, 1));
				 //----------------------------------------------------------------
				 //.........................For PES--------------------------------
				 h_profile_PES[0][1]->Fill(GetInvtMass(PC, JC), (Mass_PESup - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
				 h_profile_PES[1][1]->Fill(GetInvtMass(PC, JC), (Mass_PESdown - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));

				 h_mass_PES[0][1]->Fill(Mass_PESup, PreBTag_EvtWt*BTagEventWeight(SF_CSVL, 1));
				 h_mass_PES[1][1]->Fill(Mass_PESdown, PreBTag_EvtWt*BTagEventWeight(SF_CSVL, 1));
				 //-----------------------------------------------------------------

			       }
			       else{

				 TotEvtWt_noErr_0btag_XSm5 = PreBTag_EvtWt_XSm5;
				 TotEvtWt_noErr_0btag_XSp5 = PreBTag_EvtWt_XSp5;
				 TotEvtWt_noErr_1btag_XSm5 = 0.0;
				 TotEvtWt_noErr_1btag_XSp5 = 0.0;

			       }

			       CutFlowNumberWithWeights_XSm5[13] += TotEvtWt_noErr_0btag_XSm5;
			       CutFlowNumberWithWeights_XSm5[14] += TotEvtWt_noErr_1btag_XSm5;
			       CutFlowNumberWithWeights_XSp5[13] += TotEvtWt_noErr_0btag_XSp5;
			       CutFlowNumberWithWeights_XSp5[14] += TotEvtWt_noErr_1btag_XSp5;
			     
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

     }//END_OF if(IfSignal){
     //------------------------------------------------
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //------------------------------------
     //To Get the efficiency of LICTD Cut 
     //------------------------------------ 
     Int_t PC_NoLICTD = -1;
     Int_t JC_NoLICTD = -1;
     Pass_PhoPtcut_NoLICTD = false;
     Pass_PhoEtaEBcut_NoLICTD = false;
     Pass_JetPtcut_NoLICTD = false;
     Pass_JetEtacut_NoLICTD = false;
     Pass_GJdPhicut_NoLICTD = false;
     Pass_GJdEtacut_NoLICTD = false;
     if(IfMassCut == true){
       Pass_GJInvtMasscut_NoLICTD = false;
     }else{
       Pass_GJInvtMasscut_NoLICTD = true;
     }

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
     if(PC_NoLICTD > -1) Pass_PhoPtcut_NoLICTD = ((*Photon_pt)[PC_NoLICTD] > Cut_Photon_pt);
     if(PC_NoLICTD > -1) Pass_PhoEtaEBcut_NoLICTD = (fabs((*Photon_SC_eta)[PC_NoLICTD]) < Cut_Photon_eta);
     if(JC_NoLICTD > -1) Pass_JetPtcut_NoLICTD = ((*PFPatJet_pt)[JC_NoLICTD] > Cut_Jet_pt);
     if(JC_NoLICTD > -1) Pass_JetEtacut_NoLICTD = (fabs((*PFPatJet_eta)[JC_NoLICTD]) < Cut_Jet_eta);
     if(PC_NoLICTD > -1 && JC_NoLICTD > -1) Pass_GJdPhicut_NoLICTD = ((GetdPhi((*Photon_phi)[PC_NoLICTD], (*PFPatJet_phi)[JC_NoLICTD])) > Cut_GJdPhi);
     if(PC_NoLICTD > -1 && JC_NoLICTD > -1) Pass_GJdEtacut_NoLICTD = ((GetdEta((*Photon_SC_eta)[PC_NoLICTD], (*PFPatJet_eta)[JC_NoLICTD])) < Cut_GJdEta);
     if(IfMassCut == true){
       if(PC_NoLICTD > -1 && JC_NoLICTD > -1){
	 Pass_GJInvtMasscut_NoLICTD = ((GetInvtMass(PC_NoLICTD, JC_NoLICTD)) > Cut_GJInvtMass);
       }
     }

     if(Pass_HLT){
       if(NoScrapingEvt){
	 if(HasPrimaryVtx){
	   if(PC_NoLICTD > -1){
	     if(Pass_PhoPtcut_NoLICTD){
	       if(Pass_PhoEtaEBcut_NoLICTD){
		 if(JC_NoLICTD > -1){
		   if(Pass_JetPtcut_NoLICTD){
		     if(Pass_JetEtacut_NoLICTD){
		       if(Pass_GJdPhicut_NoLICTD){
			 if(Pass_GJdEtacut_NoLICTD){
			   if(Pass_GJInvtMasscut_NoLICTD){
			     CutFlowNumber[24]++;
			     CutFlowNumberWithWeights[42] += PreBTag_EvtWt;
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
     //---------------------------------------------------------------------------
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


   }//END_OF for(jentry) loop 

   for(int i = 0; i < nbins; i++){
     h_CutFlowTable->SetBinContent(i+1, CutFlowNumber[i]);
   }
   for(int i = 0; i < nbinsWt; i++){
     h_CutFlowTableWithWeights->SetBinContent(i+1, CutFlowNumberWithWeights[i]);
   }

   for(int i = 0; i < nbn; i++){
     h_CutFlowTableWithWeights_XSm5->SetBinContent(i+1, CutFlowNumberWithWeights_XSm5[i]);
     h_CutFlowTableWithWeights_XSp5->SetBinContent(i+1, CutFlowNumberWithWeights_XSp5[i]);
   }


   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
   Eff_PassCSVLBTag_noErr_0bTag          = CutFlowNumberWithWeights[14]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_noErr_1bTag          = CutFlowNumberWithWeights[15]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_pErr_0bTag           = CutFlowNumberWithWeights[16]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_pErr_1bTag           = CutFlowNumberWithWeights[17]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_nErr_0bTag           = CutFlowNumberWithWeights[18]/CutFlowNumberWithWeights[12];
   Eff_PassCSVLBTag_nErr_1bTag           = CutFlowNumberWithWeights[19]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_Expected             = CutFlowNumberWithWeights[20]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_noErr_0bTag          = CutFlowNumberWithWeights[21]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_noErr_1bTag          = CutFlowNumberWithWeights[22]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_pErr_0bTag           = CutFlowNumberWithWeights[23]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_pErr_1bTag           = CutFlowNumberWithWeights[24]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_nErr_0bTag           = CutFlowNumberWithWeights[25]/CutFlowNumberWithWeights[12];
   Eff_PassCSVMBTag_nErr_1bTag           = CutFlowNumberWithWeights[26]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_Expected             = CutFlowNumberWithWeights[27]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_noErr_0bTag          = CutFlowNumberWithWeights[28]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_noErr_1bTag          = CutFlowNumberWithWeights[29]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_pErr_0bTag           = CutFlowNumberWithWeights[30]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_pErr_1bTag           = CutFlowNumberWithWeights[31]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_nErr_0bTag           = CutFlowNumberWithWeights[32]/CutFlowNumberWithWeights[12];
   Eff_PassCSVTBTag_nErr_1bTag           = CutFlowNumberWithWeights[33]/CutFlowNumberWithWeights[12];
   Eff_TrueBJets                         = CutFlowNumberWithWeights[34]/CutFlowNumberWithWeights[12];
   Eff_TrueBJetsPassingCSVL              = CutFlowNumberWithWeights[35]/CutFlowNumberWithWeights[34];
   Eff_TrueBJetsPassingCSVM              = CutFlowNumberWithWeights[36]/CutFlowNumberWithWeights[34];
   Eff_TrueBJetsPassingCSVT              = CutFlowNumberWithWeights[37]/CutFlowNumberWithWeights[34];
   Eff_NonBJets                          = CutFlowNumberWithWeights[38]/CutFlowNumberWithWeights[12];
   Eff_NonBJetsPassingCSVL               = CutFlowNumberWithWeights[39]/CutFlowNumberWithWeights[38];
   Eff_NonBJetsPassingCSVM               = CutFlowNumberWithWeights[40]/CutFlowNumberWithWeights[38];
   Eff_NonBJetsPassingCSVT               = CutFlowNumberWithWeights[41]/CutFlowNumberWithWeights[38];
   Eff_LICTD                             = CutFlowNumberWithWeights[12]/CutFlowNumberWithWeights[42];

   Efferr_PassCSVL_0bTag                = fabs((CutFlowNumberWithWeights[16]-CutFlowNumberWithWeights[14])/CutFlowNumberWithWeights[14]);
   Efferr_PassCSVL_1bTag                = fabs((CutFlowNumberWithWeights[17]-CutFlowNumberWithWeights[15])/CutFlowNumberWithWeights[15]);

   Efferr_PassCSVM_0bTag                = fabs((CutFlowNumberWithWeights[23]-CutFlowNumberWithWeights[21])/CutFlowNumberWithWeights[21]);
   Efferr_PassCSVM_1bTag                = fabs((CutFlowNumberWithWeights[24]-CutFlowNumberWithWeights[22])/CutFlowNumberWithWeights[22]);

   Efferr_PassCSVT_0bTag                = fabs((CutFlowNumberWithWeights[30]-CutFlowNumberWithWeights[28])/CutFlowNumberWithWeights[28]);
   Efferr_PassCSVT_1bTag                = fabs((CutFlowNumberWithWeights[31]-CutFlowNumberWithWeights[29])/CutFlowNumberWithWeights[29]);
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
   cout << "No. of events passing CSVL BTag with nErr_0bTag weight = " << CutFlowNumberWithWeights[18] << endl;
   cout << "No. of events passing CSVL BTag with nErr_1bTag weight = " << CutFlowNumberWithWeights[19] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVM BTag = " << CutFlowNumber[14] << "," << CutFlowNumberWithWeights[20] << endl;
   cout << "No. of events passing CSVM BTag with noErr_0bTag weight = " << CutFlowNumberWithWeights[21] << endl;
   cout << "No. of events passing CSVM BTag with noErr_1bTag weight = " << CutFlowNumberWithWeights[22] << endl;
   cout << "No. of events passing CSVM BTag with pErr_0bTag weight = " << CutFlowNumberWithWeights[23] << endl;
   cout << "No. of events passing CSVM BTag with pErr_1bTag weight = " << CutFlowNumberWithWeights[24] << endl;
   cout << "No. of events passing CSVM BTag with nErr_0bTag weight = " << CutFlowNumberWithWeights[25] << endl;
   cout << "No. of events passing CSVM BTag with nErr_1bTag weight = " << CutFlowNumberWithWeights[26] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing CSVT BTag = " << CutFlowNumber[15] << "," << CutFlowNumberWithWeights[27] << endl;
   cout << "No. of events passing CSVT BTag with noErr_0bTag weight = " << CutFlowNumberWithWeights[28] << endl;
   cout << "No. of events passing CSVT BTag with noErr_1bTag weight = " << CutFlowNumberWithWeights[29] << endl;
   cout << "No. of events passing CSVT BTag with pErr_0bTag weight = " << CutFlowNumberWithWeights[30] << endl;
   cout << "No. of events passing CSVT BTag with pErr_1bTag weight = " << CutFlowNumberWithWeights[31] << endl;
   cout << "No. of events passing CSVT BTag with nErr_0bTag weight = " << CutFlowNumberWithWeights[32] << endl;
   cout << "No. of events passing CSVT BTag with nErr_1bTag weight = " << CutFlowNumberWithWeights[33] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of True BJets = " << CutFlowNumber[16] << "," << CutFlowNumberWithWeights[34] << endl;
   cout << "No. of True BJets passing CSVL = " << CutFlowNumber[17] << "," << CutFlowNumberWithWeights[35] << endl;
   cout << "No. of True BJets passing CSVM = " << CutFlowNumber[18] << "," << CutFlowNumberWithWeights[36] << endl;
   cout << "No. of True BJets passing CSVT = " << CutFlowNumber[19] << "," << CutFlowNumberWithWeights[37] << endl;
   cout << "No.of Non BJets =" << CutFlowNumber[20] << "," << CutFlowNumberWithWeights[38] << endl;
   cout << "No.of Non BJets passing CSVL = " << CutFlowNumber[21] << "," << CutFlowNumberWithWeights[39] <<endl;
   cout << "No. of Non BJets passing CSVM = " << CutFlowNumber[22] << "," << CutFlowNumberWithWeights[40] << endl;
   cout << "No. of Non BJets passing CSVT = " << CutFlowNumber[23] << "," << CutFlowNumberWithWeights[41] << endl;
   cout << "*****************************************************************************************************" << endl;
   cout << "No. of events passing with LICTD cut = " << CutFlowNumberWithWeights[12] << endl;
   cout << "No. of events passing without LICTD cut without weight = " << CutFlowNumber[24] << endl;
   cout << "No. of events passing without LICTD cut with weight = " << CutFlowNumberWithWeights[42] << endl;

   cout << "*****************************************************************************************************" << endl;
   cout <<"*********Number of events for Pileup syst corresponding to 5% less XS***********" << endl;
   cout << "****************************************************************************************************" << endl;
   cout << " Total no. of events = " << CutFlowNumberWithWeights_XSm5[0] << endl;
   cout << "No. of events passing HLT = " << CutFlowNumberWithWeights_XSm5[1] << endl;
   cout << "No. of events passing PrimaryVtx =" << CutFlowNumberWithWeights_XSm5[2] << endl;
   cout << "No. of events passing Scraping = " << CutFlowNumberWithWeights_XSm5[3] << endl;
   cout << "No. of events passing PhotonID = " << CutFlowNumberWithWeights_XSm5[4] << endl;
   cout << "No. of events passing PhotonPt = " << CutFlowNumberWithWeights_XSm5[5] << endl;
   cout << "No. of events passing PhotonEta = " << CutFlowNumberWithWeights_XSm5[6] << endl;
   cout << "No. of events passing JetID = " << CutFlowNumberWithWeights_XSm5[7] << endl;
   cout << "No. of events passing JetPt = "  << CutFlowNumberWithWeights_XSm5[8] << endl;
   cout << "No. of events passing JetEta = " << CutFlowNumberWithWeights_XSm5[9] << endl;
   cout << "No. of events passing DPhiCut = " << CutFlowNumberWithWeights_XSm5[10] << endl;
   cout << "No. of events passing DetaCut = " << CutFlowNumberWithWeights_XSm5[11] << endl;
   cout << "No. of events passing InvtMass = "  << CutFlowNumberWithWeights_XSm5[12] << endl;
   cout << "No. of events passing CSVL_noErr_0bTag = " << CutFlowNumberWithWeights_XSm5[13] << endl;
   cout << "No. of events passing CSVL_noErr_1bTag = " << CutFlowNumberWithWeights_XSm5[14] << endl;
   cout << "****************************************************************************************************" << endl;
   cout <<"*********Number of events for Pileup syst corresponding to 5% more XS***********" << endl;
   cout << "****************************************************************************************************" << endl;
   cout << " Total no. of events = " << CutFlowNumberWithWeights_XSp5[0] << endl;
   cout << "No. of events passing HLT = " << CutFlowNumberWithWeights_XSp5[1] << endl;
   cout << "No. of events passing PrimaryVtx =" << CutFlowNumberWithWeights_XSp5[2] << endl;
   cout << "No. of events passing Scraping = " << CutFlowNumberWithWeights_XSp5[3] << endl;
   cout << "No. of events passing PhotonID = " << CutFlowNumberWithWeights_XSp5[4] << endl;
   cout << "No. of events passing PhotonPt = " << CutFlowNumberWithWeights_XSp5[5] << endl;
   cout << "No. of events passing PhotonEta = " << CutFlowNumberWithWeights_XSp5[6] << endl;
   cout << "No. of events passing JetID = " << CutFlowNumberWithWeights_XSp5[7] << endl;
   cout << "No. of events passing JetPt = "  << CutFlowNumberWithWeights_XSp5[8] << endl;
   cout << "No. of events passing JetEta = " << CutFlowNumberWithWeights_XSp5[9] << endl;
   cout << "No. of events passing DPhiCut = " << CutFlowNumberWithWeights_XSp5[10] << endl;
   cout << "No. of events passing DetaCut = " << CutFlowNumberWithWeights_XSp5[11] << endl;
   cout << "No. of events passing InvtMass = "  << CutFlowNumberWithWeights_XSp5[12] << endl;
   cout << "No. of events passing CSVL_noErr_0bTag = " << CutFlowNumberWithWeights_XSp5[13] << endl;
   cout << "No. of events passing CSVL_noErr_1bTag = " << CutFlowNumberWithWeights_XSp5[14] << endl;
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
   cout << "Eff_PassCSVLBTag_nErr_0bTag(CSVL_MistagRate_nErr)                     = " << Eff_PassCSVLBTag_nErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVLBTag_nErr_1bTag(CSVL_TrueEff_nErr)                        = " << Eff_PassCSVLBTag_nErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_Expected                                             = " << Eff_PassCSVMBTag_Expected*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_noErr_0bTag(CSVM_MistagRate_noErr)                   = " << Eff_PassCSVMBTag_noErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_noErr_1bTag(CSVM_TrueEff_noErr)                      = " << Eff_PassCSVMBTag_noErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_pErr_0bTag(CSVM_MistagRate_pErr)                     = " << Eff_PassCSVMBTag_pErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_pErr_1bTag(CSVM_TrueEff_pErr)                        = " << Eff_PassCSVMBTag_pErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_nErr_0bTag(CSVM_MistagRate_nErr)                     = " << Eff_PassCSVMBTag_nErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVMBTag_nErr_1bTag(CSVM_TrueEff_nErr)                        = " << Eff_PassCSVMBTag_nErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_Expected                                             = " << Eff_PassCSVTBTag_Expected*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_noErr_0bTag(CSVT_MistagRate_noErr)                   = " << Eff_PassCSVTBTag_noErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_noErr_1bTag(CSVT_TrueEff_noErr)                      = " << Eff_PassCSVTBTag_noErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_pErr_0bTag(CSVT_MistagRate_pErr)                     = " << Eff_PassCSVTBTag_pErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_pErr_1bTag(CSVT_TrueEff_pErr)                        = " << Eff_PassCSVTBTag_pErr_1bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_nErr_0bTag(CSVT_MistagRate_nErr)                     = " << Eff_PassCSVTBTag_nErr_0bTag*100 << "%" << endl;
   cout << "Eff_PassCSVTBTag_nErr_1bTag(CSVT_TrueEff_nErr)                        = " << Eff_PassCSVTBTag_nErr_1bTag*100 << "%" << endl;
   cout << "Eff_TrueBJets                                                         = " << Eff_TrueBJets*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVL(CSVL_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVL*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVM(CSVM_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVM*100 << "%" << endl;
   cout << "Eff_TrueBJetsPassingCSVT(CSVT_TrueEff_FromMCTruth)                    = " << Eff_TrueBJetsPassingCSVT*100 << "%" << endl;
   cout << "Eff_NonBJets                                                          = " << Eff_NonBJets*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVL(CSVL_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVL*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVM(CSVM_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVM*100 << "%" << endl;
   cout << "Eff_NonBJetsPassingCSVT(CSVT_MistagRate_FromMCTruth)                  = " << Eff_NonBJetsPassingCSVT*100 << "%" << endl;
   cout << "Eff_LICTD           = " << Eff_LICTD*100 << "%" << endl;
   cout << "******************************************************************************************************" << endl;
   cout << "*********Efficiency and Eff Error Parameters required for Tagging Eff plots and Limit setting*********" << endl;
   cout << "**Relative Error = (Evts_pErr-Evts_noErr)/Evts_noErr or (Eff_pErr-Eff_noErr)/Eff_noErr (This one is wrong , not required)**" << endl;
   cout << "**Absolute Error = (Eff_pErr-Eff_noErr) (Required for Tagging Eff plots and for limit setting**" << endl;
   cout << "******************************************************************************************************" << endl;
   cout << "Eff_PassCSVLBTag_noErr_0bTag    = " << Eff_PassCSVLBTag_noErr_0bTag << endl;
   cout << "Relative Efferr_PassCSVL_0bTag     = " << Efferr_PassCSVL_0bTag << endl;
   cout << "Absolute Efferr_PassCSVL_0bTag     = " << Efferr_PassCSVL_0bTag * Eff_PassCSVLBTag_noErr_0bTag<< endl;
   cout << "Eff_PassCSVLBTag_noErr_1bTag    = " << Eff_PassCSVLBTag_noErr_1bTag << endl;
   cout << "Relative Efferr_PassCSVL_1bTag     = " << Efferr_PassCSVL_1bTag << endl;
   cout << "Absolute Efferr_PassCSVL_1bTag     = " << Efferr_PassCSVL_1bTag * Eff_PassCSVLBTag_noErr_1bTag << endl;
   cout << "                                                                              " << endl;
   cout << "Eff_PassCSVMBTag_noErr_0bTag    = " << Eff_PassCSVMBTag_noErr_0bTag << endl;
   cout << "Relative Efferr_PassCSVM_0bTag     = " << Efferr_PassCSVM_0bTag << endl;
   cout << "Absolute Efferr_PassCSVM_0bTag     = " << Efferr_PassCSVM_0bTag * Eff_PassCSVMBTag_noErr_0bTag<< endl;
   cout << "Eff_PassCSVMBTag_noErr_1bTag    = " << Eff_PassCSVMBTag_noErr_1bTag << endl;
   cout << "Relative Efferr_PassCSVM_1bTag     = " << Efferr_PassCSVM_1bTag << endl;
   cout << "Absolute Efferr_PassCSVM_1bTag     = " << Efferr_PassCSVM_1bTag * Eff_PassCSVMBTag_noErr_1bTag<< endl;
   cout << "                                                                              " << endl;
   cout << "Eff_PassCSVTBTag_noErr_0bTag    = " << Eff_PassCSVTBTag_noErr_0bTag << endl;
   cout << "Relative Efferr_PassCSVT_0bTag     = " << Efferr_PassCSVT_0bTag << endl;
   cout << "Absolute Efferr_PassCSVT_0bTag     = " << Efferr_PassCSVT_0bTag * Eff_PassCSVTBTag_noErr_0bTag << endl;
   cout << "Eff_PassCSVTBTag_noErr_1bTag    = " << Eff_PassCSVTBTag_noErr_1bTag << endl;
   cout << "Relative Efferr_PassCSVT_1bTag     = " << Efferr_PassCSVT_1bTag << endl;
   cout << "Absolute Efferr_PassCSVT_1bTag     = " << Efferr_PassCSVT_1bTag * Eff_PassCSVTBTag_noErr_1bTag << endl;
   cout << "******************************************************************************************************" << endl;
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




}






EOF



cat >> PostAnalyzer_MC.h <<EOF 
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 03:42:46 2014 by ROOT version 5.32/00
// from TChain myEvent/
//////////////////////////////////////////////////////////

#ifndef PostAnalyzer_MC_h
#define PostAnalyzer_MC_h

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
#include <TProfile.h>

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
const Int_t kMaxhasMatchedGenPhoton = 1;
const Int_t kMaxhasMatchedGenParton_ToPFJet = 1;
const Int_t kMaxhasMatchedGenJet_ToPFJet = 1;
const Int_t kMaxIsScrapingEvent = 1;

class PostAnalyzer_MC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //Define Output File
   TFile *file;

   //If need to apply mass cut or not
   Bool_t IfMassCut;

   //For running on signal if pileup syst is needed
   Bool_t IfSignal;

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

   //For LICTD Efficiency
   Bool_t Pass_PhoPtcut_NoLICTD;
   Bool_t Pass_PhoEtaEBcut_NoLICTD;
   Bool_t Pass_JetPtcut_NoLICTD;
   Bool_t Pass_JetEtacut_NoLICTD;
   Bool_t Pass_GJdPhicut_NoLICTD;
   Bool_t Pass_GJdEtacut_NoLICTD;
   Bool_t Pass_GJInvtMasscut_NoLICTD;

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

   //Various Efficiency Variables
   Double_t Eff_PassHLT, Eff_PassScraping, Eff_PassPrimaryVtx, Eff_PassPhotonID, Eff_PassPhotonPt, Eff_PassPhotonEta;
   Double_t Eff_PassJetID, Eff_PassJetPt, Eff_PassJetEta, Eff_PassDPhi, Eff_PassDEta, Eff_PassGJInvtMass;
   Double_t Eff_PassCSVLBTag_Expected, Eff_PassCSVLBTag_noErr_0bTag, Eff_PassCSVLBTag_noErr_1bTag, Eff_PassCSVLBTag_pErr_0bTag, Eff_PassCSVLBTag_pErr_1bTag, Eff_PassCSVLBTag_nErr_0bTag, Eff_PassCSVLBTag_nErr_1bTag;
   Double_t Eff_PassCSVMBTag_Expected, Eff_PassCSVMBTag_noErr_0bTag, Eff_PassCSVMBTag_noErr_1bTag, Eff_PassCSVMBTag_pErr_0bTag, Eff_PassCSVMBTag_pErr_1bTag, Eff_PassCSVMBTag_nErr_0bTag, Eff_PassCSVMBTag_nErr_1bTag;           
   Double_t Eff_PassCSVTBTag_Expected, Eff_PassCSVTBTag_noErr_0bTag, Eff_PassCSVTBTag_noErr_1bTag, Eff_PassCSVTBTag_pErr_0bTag, Eff_PassCSVTBTag_pErr_1bTag, Eff_PassCSVTBTag_nErr_0bTag, Eff_PassCSVTBTag_nErr_1bTag;           
   Double_t Eff_TrueBJets, Eff_TrueBJetsPassingCSVL, Eff_TrueBJetsPassingCSVM, Eff_TrueBJetsPassingCSVT;              
   Double_t Eff_NonBJets, Eff_NonBJetsPassingCSVL, Eff_NonBJetsPassingCSVM, Eff_NonBJetsPassingCSVT;               
   Double_t Eff_LICTD;
   Double_t Efferr_PassCSVL_0bTag, Efferr_PassCSVL_1bTag;
   Double_t Efferr_PassCSVM_0bTag, Efferr_PassCSVM_1bTag;
   Double_t Efferr_PassCSVT_0bTag, Efferr_PassCSVT_1bTag;
   
   //Self Decleared Variables for MC only
   Double_t Lumi;
   Double_t Lumi_EvtWt, PU_EvtWt, PreBTag_EvtWt;
   Double_t PU_EvtWt_XSm5, PreBTag_EvtWt_XSm5, PU_EvtWt_XSp5, PreBTag_EvtWt_XSp5;
   Double_t Total_CSVLEvtWt_noErr_0bTag, Total_CSVLEvtWt_noErr_1bTag, Total_CSVLEvtWt_pErr_0bTag, Total_CSVLEvtWt_pErr_1bTag, Total_CSVLEvtWt_nErr_0bTag, Total_CSVLEvtWt_nErr_1bTag;
   Double_t Total_CSVMEvtWt_noErr_0bTag, Total_CSVMEvtWt_noErr_1bTag, Total_CSVMEvtWt_pErr_0bTag, Total_CSVMEvtWt_pErr_1bTag, Total_CSVMEvtWt_nErr_0bTag, Total_CSVMEvtWt_nErr_1bTag;
   Double_t Total_CSVTEvtWt_noErr_0bTag, Total_CSVTEvtWt_noErr_1bTag, Total_CSVTEvtWt_pErr_0bTag, Total_CSVTEvtWt_pErr_1bTag, Total_CSVTEvtWt_nErr_0bTag, Total_CSVTEvtWt_nErr_1bTag;
   //For JES in Signal requried for syst
   TLorentzVector PFJet4Vec_JESup, PFJet4Vec_JESdown;
   TLorentzVector GJ4Vec_JESup, GJ4Vec_JESdown;
   Double_t Mass_JESup, Mass_JESdown;
   //For PES in Signal required for syst
   TLorentzVector Photon4Vec_PESup, Photon4Vec_PESdown;
   TLorentzVector GJ4Vec_PESup, GJ4Vec_PESdown;
   Double_t Mass_PESup, Mass_PESdown; 
   Double_t Photon_pecUnc;

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

   //Histograms for CSVL BJets
   TH1F *h_CSVL_BJetPt_noErr_0bTag;
   TH1F *h_CSVL_BJetEta_noErr_0bTag;
   TH1F *h_CSVL_BJetPhi_noErr_0bTag;

   TH1F *h_CSVL_BJetPt_noErr_1bTag;
   TH1F *h_CSVL_BJetEta_noErr_1bTag;
   TH1F *h_CSVL_BJetPhi_noErr_1bTag;

   TH1F *h_CSVL_BJetPt_pErr_0bTag;
   TH1F *h_CSVL_BJetEta_pErr_0bTag;
   TH1F *h_CSVL_BJetPhi_pErr_0bTag;

   TH1F *h_CSVL_BJetPt_pErr_1bTag;
   TH1F *h_CSVL_BJetEta_pErr_1bTag;
   TH1F *h_CSVL_BJetPhi_pErr_1bTag;

   //Histograms for Photon-CSVLBJet
   TH1F *h_CSVL_GBJetInvtMass_binning40_noErr_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag;
   TH1F *h_CSVL_GBJetdEta_noErr_0bTag;
   TH1F *h_CSVL_GBJetdPhi_noErr_0bTag;
   TH1F *h_CSVL_GBJetdR_noErr_0bTag;

   TH1F *h_CSVL_GBJetInvtMass_binning40_noErr_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag;
   TH1F *h_CSVL_GBJetdEta_noErr_1bTag;
   TH1F *h_CSVL_GBJetdPhi_noErr_1bTag;
   TH1F *h_CSVL_GBJetdR_noErr_1bTag;

   TH1F *h_CSVL_GBJetInvtMass_binning40_pErr_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag;
   TH1F *h_CSVL_GBJetdEta_pErr_0bTag;
   TH1F *h_CSVL_GBJetdPhi_pErr_0bTag;
   TH1F *h_CSVL_GBJetdR_pErr_0bTag;

   TH1F *h_CSVL_GBJetInvtMass_binning40_pErr_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag;
   TH1F *h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag;
   TH1F *h_CSVL_GBJetdEta_pErr_1bTag;
   TH1F *h_CSVL_GBJetdPhi_pErr_1bTag;
   TH1F *h_CSVL_GBJetdR_pErr_1bTag;

   //Histograms for CSVM BJets
   TH1F *h_CSVM_BJetPt_noErr_0bTag;
   TH1F *h_CSVM_BJetEta_noErr_0bTag;
   TH1F *h_CSVM_BJetPhi_noErr_0bTag;

   TH1F *h_CSVM_BJetPt_noErr_1bTag;
   TH1F *h_CSVM_BJetEta_noErr_1bTag;
   TH1F *h_CSVM_BJetPhi_noErr_1bTag;

   TH1F *h_CSVM_BJetPt_pErr_0bTag;
   TH1F *h_CSVM_BJetEta_pErr_0bTag;
   TH1F *h_CSVM_BJetPhi_pErr_0bTag;

   TH1F *h_CSVM_BJetPt_pErr_1bTag;
   TH1F *h_CSVM_BJetEta_pErr_1bTag;
   TH1F *h_CSVM_BJetPhi_pErr_1bTag;

   //Histograms for Photon-CSVMBJet
   TH1F *h_CSVM_GBJetInvtMass_binning40_noErr_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag;
   TH1F *h_CSVM_GBJetdEta_noErr_0bTag;
   TH1F *h_CSVM_GBJetdPhi_noErr_0bTag;
   TH1F *h_CSVM_GBJetdR_noErr_0bTag;

   TH1F *h_CSVM_GBJetInvtMass_binning40_noErr_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag;
   TH1F *h_CSVM_GBJetdEta_noErr_1bTag;
   TH1F *h_CSVM_GBJetdPhi_noErr_1bTag;
   TH1F *h_CSVM_GBJetdR_noErr_1bTag;

   TH1F *h_CSVM_GBJetInvtMass_binning40_pErr_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag;
   TH1F *h_CSVM_GBJetdEta_pErr_0bTag;
   TH1F *h_CSVM_GBJetdPhi_pErr_0bTag;
   TH1F *h_CSVM_GBJetdR_pErr_0bTag;

   TH1F *h_CSVM_GBJetInvtMass_binning40_pErr_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag;
   TH1F *h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag;
   TH1F *h_CSVM_GBJetdEta_pErr_1bTag;
   TH1F *h_CSVM_GBJetdPhi_pErr_1bTag;
   TH1F *h_CSVM_GBJetdR_pErr_1bTag;

   //Histograms for CSVT BJets
   TH1F *h_CSVT_BJetPt_noErr_0bTag;
   TH1F *h_CSVT_BJetEta_noErr_0bTag;
   TH1F *h_CSVT_BJetPhi_noErr_0bTag;

   TH1F *h_CSVT_BJetPt_noErr_1bTag;
   TH1F *h_CSVT_BJetEta_noErr_1bTag;
   TH1F *h_CSVT_BJetPhi_noErr_1bTag;

   TH1F *h_CSVT_BJetPt_pErr_0bTag;
   TH1F *h_CSVT_BJetEta_pErr_0bTag;
   TH1F *h_CSVT_BJetPhi_pErr_0bTag;

   TH1F *h_CSVT_BJetPt_pErr_1bTag;
   TH1F *h_CSVT_BJetEta_pErr_1bTag;
   TH1F *h_CSVT_BJetPhi_pErr_1bTag;

   //Histograms for Photon-CSVTBJet
   TH1F *h_CSVT_GBJetInvtMass_binning40_noErr_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag;
   TH1F *h_CSVT_GBJetdEta_noErr_0bTag;
   TH1F *h_CSVT_GBJetdPhi_noErr_0bTag;
   TH1F *h_CSVT_GBJetdR_noErr_0bTag;

   TH1F *h_CSVT_GBJetInvtMass_binning40_noErr_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag;
   TH1F *h_CSVT_GBJetdEta_noErr_1bTag;
   TH1F *h_CSVT_GBJetdPhi_noErr_1bTag;
   TH1F *h_CSVT_GBJetdR_noErr_1bTag;

   TH1F *h_CSVT_GBJetInvtMass_binning40_pErr_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag;
   TH1F *h_CSVT_GBJetdEta_pErr_0bTag;
   TH1F *h_CSVT_GBJetdPhi_pErr_0bTag;
   TH1F *h_CSVT_GBJetdR_pErr_0bTag;

   TH1F *h_CSVT_GBJetInvtMass_binning40_pErr_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag;
   TH1F *h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag;
   TH1F *h_CSVT_GBJetdEta_pErr_1bTag;
   TH1F *h_CSVT_GBJetdPhi_pErr_1bTag;
   TH1F *h_CSVT_GBJetdR_pErr_1bTag;

   //Histograms for BJetDisc before and after cut
   TH1F *h_BJetDiscByCSV;
   TH1F *h_BJetDiscByCSV_PassingCSVL_noErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVL_noErr_1bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVL_pErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVL_pErr_1bTag;

   TH1F *h_BJetDiscByCSV_PassingCSVM_noErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVM_noErr_1bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVM_pErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVM_pErr_1bTag;

   TH1F *h_BJetDiscByCSV_PassingCSVT_noErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVT_noErr_1bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVT_pErr_0bTag;
   TH1F *h_BJetDiscByCSV_PassingCSVT_pErr_1bTag;

   //Histograms for Pileup effect
   TH1F *h_trueNumofInteractions;
   TH1F *h_goodPV;
   TH1F *h_goodPV_LumiWt;
   TH1F *h_goodPV_PUWt;
   TH1F *h_goodPV_LumiNPUWt;

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
   TH1F *h_CutFlowTableWithWeights;
   TH1F *h_CutFlowTableWithWeights_XSm5; //for pileup rewting with (XS - 5% of XS) cross section value needed for pileup systematics
   TH1F *h_CutFlowTableWithWeights_XSp5; //for pileup rewting with (XS + 5% of XS) cross section value needed for pileup systematics

   //Histograms defined for MC only
   TH1F *h_DataPUDist;
   TH1F *h_DataPUNormDist;
   TH1F *h_MCPUNormDist;
   TH1F *h_PUScaleFactor;
   TH1F *h_PUScaleFactor_XSm5;
   TH1F *h_PUScaleFactor_XSp5;

   TH1F *h_Jet_PartonFlavor;

   TH1F *h_CSVL_BJet_PartonFlavor_noErr_0bTag;
   TH1F *h_CSVL_BJet_PartonFlavor_noErr_1bTag;
   TH1F *h_CSVL_BJet_PartonFlavor_pErr_0bTag;
   TH1F *h_CSVL_BJet_PartonFlavor_pErr_1bTag;

   TH1F *h_CSVM_BJet_PartonFlavor_noErr_0bTag;
   TH1F *h_CSVM_BJet_PartonFlavor_noErr_1bTag;
   TH1F *h_CSVM_BJet_PartonFlavor_pErr_0bTag;
   TH1F *h_CSVM_BJet_PartonFlavor_pErr_1bTag;

   TH1F *h_CSVT_BJet_PartonFlavor_noErr_0bTag;
   TH1F *h_CSVT_BJet_PartonFlavor_noErr_1bTag;
   TH1F *h_CSVT_BJet_PartonFlavor_pErr_0bTag;
   TH1F *h_CSVT_BJet_PartonFlavor_pErr_1bTag;

   TH1F *h_TrueBJetPt;
   TH1F *h_TrueBJetEta;
   TH1F *h_TrueBJetPtPassingCSVL;
   TH1F *h_TrueBJetEtaPassingCSVL;
   TH1F *h_TrueBJetPtPassingCSVM;
   TH1F *h_TrueBJetEtaPassingCSVM;
   TH1F *h_TrueBJetPtPassingCSVT;
   TH1F *h_TrueBJetEtaPassingCSVT;
   TH1F *h_NonBJetPt;
   TH1F *h_NonBJetEta;
   TH1F *h_NonBJetPtPassingCSVL;
   TH1F *h_NonBJetEtaPassingCSVL;
   TH1F *h_NonBJetPtPassingCSVM;
   TH1F *h_NonBJetEtaPassingCSVM;
   TH1F *h_NonBJetPtPassingCSVT;
   TH1F *h_NonBJetEtaPassingCSVT;

   //For JES
   TProfile *h_profile_JES[2][2];
   TH1F *h_mass_JES[2][2];
   //For PES
   TProfile *h_profile_PES[2][2];
   TH1F *h_mass_PES[2][2];

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
   Int_t           npuVertices;
   Int_t           npuVerticesm1;
   Int_t           npuVerticesp1;
   Int_t           npuVerticespm2;
   Float_t         trueNumofInteractions;
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
   Int_t           GenPhoton_n;
   vector<int>     *GenPhoton_status;
   vector<double>  *GenPhoton_E;
   vector<double>  *GenPhoton_et;
   vector<double>  *GenPhoton_pt;
   vector<double>  *GenPhoton_eta;
   vector<double>  *GenPhoton_phi;
   vector<double>  *GenPhoton_px;
   vector<double>  *GenPhoton_py;
   vector<double>  *GenPhoton_pz;
   vector<double>  *GenPhoton_vx;
   vector<double>  *GenPhoton_vy;
   vector<double>  *GenPhoton_vz;
   vector<int>     *GenPhoton_mother_status;
   vector<int>     *GenPhoton_mother_pdgid;
   vector<double>  *GenPhoton_mother_E;
   vector<double>  *GenPhoton_mother_et;
   vector<double>  *GenPhoton_mother_pt;
   vector<double>  *GenPhoton_mother_eta;
   vector<double>  *GenPhoton_mother_phi;
   vector<double>  *GenPhoton_mother_px;
   vector<double>  *GenPhoton_mother_py;
   vector<double>  *GenPhoton_mother_pz;
   vector<double>  *GenPhoton_mother_vx;
   vector<double>  *GenPhoton_mother_vy;
   vector<double>  *GenPhoton_mother_vz;
   vector<int>     *GenPhoton_grandmother_status;
   vector<int>     *GenPhoton_grandmother_pdgid;
   vector<double>  *GenPhoton_grandmother_E;
   vector<double>  *GenPhoton_grandmother_et;
   vector<double>  *GenPhoton_grandmother_pt;
   vector<double>  *GenPhoton_grandmother_eta;
   vector<double>  *GenPhoton_grandmother_phi;
   vector<double>  *GenPhoton_grandmother_px;
   vector<double>  *GenPhoton_grandmother_py;
   vector<double>  *GenPhoton_grandmother_pz;
   vector<double>  *GenPhoton_grandmother_vx;
   vector<double>  *GenPhoton_grandmother_vy;
   vector<double>  *GenPhoton_grandmother_vz;
   vector<bool>    *hasMatchedGenPhoton_;
   vector<int>     *MatchedGenPhoton_status;
   vector<int>     *MatchedGenPhoton_pdgid;
   vector<double>  *MatchedGenPhoton_E;
   vector<double>  *MatchedGenPhoton_et;
   vector<double>  *MatchedGenPhoton_pt;
   vector<double>  *MatchedGenPhoton_eta;
   vector<double>  *MatchedGenPhoton_phi;
   vector<double>  *MatchedGenPhoton_px;
   vector<double>  *MatchedGenPhoton_py;
   vector<double>  *MatchedGenPhoton_pz;
   vector<double>  *MatchedGenPhoton_vx;
   vector<double>  *MatchedGenPhoton_vy;
   vector<double>  *MatchedGenPhoton_vz;
   vector<int>     *MatchedGenPhoton_nMothers;
   vector<vector<int> > *MatchedGenPhoton_mothersPdgid;
   vector<vector<int> > *MatchedGenPhoton_mothersStatus;
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
   vector<bool>    *hasMatchedGenParton_ToPFJet_;
   vector<bool>    *hasMatchedGenJet_ToPFJet_;
   vector<int>     *MatchedGenPartonToPFJet_pdgId;
   vector<int>     *MatchedGenPartonToPFJet_status;
   vector<double>  *MatchedGenJetToPFJet_E;
   vector<double>  *MatchedGenJetToPFJet_pt;
   vector<double>  *MatchedGenJetToPFJet_eta;
   vector<double>  *MatchedGenJetToPFJet_phi;
   vector<int>     *MatchedGenJetToPFJet_nDaughters;
   vector<int>     *PFPatJet_JetPartonFlavor;
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
   Int_t           GenJet_n;
   vector<double>  *GenJet_E;
   vector<double>  *GenJet_et;
   vector<double>  *GenJet_pt;
   vector<double>  *GenJet_eta;
   vector<double>  *GenJet_phi;
   vector<double>  *GenJet_px;
   vector<double>  *GenJet_py;
   vector<double>  *GenJet_pz;
   vector<double>  *GenJet_vx;
   vector<double>  *GenJet_vy;
   vector<double>  *GenJet_vz;
   vector<int>     *GenJet_nDaughters;
   vector<vector<int> > *GenJet_daughters_status;
   vector<vector<int> > *GenJet_daughters_pdgid;
   vector<vector<double> > *GenJet_daughters_E;
   vector<vector<double> > *GenJet_daughters_pt;
   vector<vector<double> > *GenJet_daughters_eta;
   vector<vector<double> > *GenJet_daughters_phi;
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
   TBranch        *b_npuVertices;   //!
   TBranch        *b_npuVerticesm1;   //!
   TBranch        *b_npuVerticesp1;   //!
   TBranch        *b_npuVerticespm2;   //!
   TBranch        *b_trueNumofInteractions;   //!
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
   TBranch        *b_GenPhoton_n;   //!
   TBranch        *b_GenPhoton_status;   //!
   TBranch        *b_GenPhoton_E;   //!
   TBranch        *b_GenPhoton_et;   //!
   TBranch        *b_GenPhoton_pt;   //!
   TBranch        *b_GenPhoton_eta;   //!
   TBranch        *b_GenPhoton_phi;   //!
   TBranch        *b_GenPhoton_px;   //!
   TBranch        *b_GenPhoton_py;   //!
   TBranch        *b_GenPhoton_pz;   //!
   TBranch        *b_GenPhoton_vx;   //!
   TBranch        *b_GenPhoton_vy;   //!
   TBranch        *b_GenPhoton_vz;   //!
   TBranch        *b_GenPhoton_mother_status;   //!
   TBranch        *b_GenPhoton_mother_pdgid;   //!
   TBranch        *b_GenPhoton_mother_E;   //!
   TBranch        *b_GenPhoton_mother_et;   //!
   TBranch        *b_GenPhoton_mother_pt;   //!
   TBranch        *b_GenPhoton_mother_eta;   //!
   TBranch        *b_GenPhoton_mother_phi;   //!
   TBranch        *b_GenPhoton_mother_px;   //!
   TBranch        *b_GenPhoton_mother_py;   //!
   TBranch        *b_GenPhoton_mother_pz;   //!
   TBranch        *b_GenPhoton_mother_vx;   //!
   TBranch        *b_GenPhoton_mother_vy;   //!
   TBranch        *b_GenPhoton_mother_vz;   //!
   TBranch        *b_GenPhoton_grandmother_status;   //!
   TBranch        *b_GenPhoton_grandmother_pdgid;   //!
   TBranch        *b_GenPhoton_grandmother_E;   //!
   TBranch        *b_GenPhoton_grandmother_et;   //!
   TBranch        *b_GenPhoton_grandmother_pt;   //!
   TBranch        *b_GenPhoton_grandmother_eta;   //!
   TBranch        *b_GenPhoton_grandmother_phi;   //!
   TBranch        *b_GenPhoton_grandmother_px;   //!
   TBranch        *b_GenPhoton_grandmother_py;   //!
   TBranch        *b_GenPhoton_grandmother_pz;   //!
   TBranch        *b_GenPhoton_grandmother_vx;   //!
   TBranch        *b_GenPhoton_grandmother_vy;   //!
   TBranch        *b_GenPhoton_grandmother_vz;   //!
   TBranch        *b_hasMatchedGenPhoton_;   //!
   TBranch        *b_MatchedGenPhoton_status;   //!
   TBranch        *b_MatchedGenPhoton_pdgid;   //!
   TBranch        *b_MatchedGenPhoton_E;   //!
   TBranch        *b_MatchedGenPhoton_et;   //!
   TBranch        *b_MatchedGenPhoton_pt;   //!
   TBranch        *b_MatchedGenPhoton_eta;   //!
   TBranch        *b_MatchedGenPhoton_phi;   //!
   TBranch        *b_MatchedGenPhoton_px;   //!
   TBranch        *b_MatchedGenPhoton_py;   //!
   TBranch        *b_MatchedGenPhoton_pz;   //!
   TBranch        *b_MatchedGenPhoton_vx;   //!
   TBranch        *b_MatchedGenPhoton_vy;   //!
   TBranch        *b_MatchedGenPhoton_vz;   //!
   TBranch        *b_MatchedGenPhoton_nMothers;   //!
   TBranch        *b_MatchedGenPhoton_mothersPdgid;   //!
   TBranch        *b_MatchedGenPhoton_mothersStatus;   //!
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
   TBranch        *b_hasMatchedGenParton_ToPFJet_;   //!
   TBranch        *b_hasMatchedGenJet_ToPFJet_;   //!
   TBranch        *b_MatchedGenPartonToPFJet_pdgId;   //!
   TBranch        *b_MatchedGenPartonToPFJet_status;   //!
   TBranch        *b_MatchedGenJetToPFJet_E;   //!
   TBranch        *b_MatchedGenJetToPFJet_pt;   //!
   TBranch        *b_MatchedGenJetToPFJet_eta;   //!
   TBranch        *b_MatchedGenJetToPFJet_phi;   //!
   TBranch        *b_MatchedGenJetToPFJet_nDaughters;   //!
   TBranch        *b_PFPatJet_JetPartonFlavor;   //!
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
   TBranch        *b_GenJet_n;   //!
   TBranch        *b_GenJet_E;   //!
   TBranch        *b_GenJet_et;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_px;   //!
   TBranch        *b_GenJet_py;   //!
   TBranch        *b_GenJet_pz;   //!
   TBranch        *b_GenJet_vx;   //!
   TBranch        *b_GenJet_vy;   //!
   TBranch        *b_GenJet_vz;   //!
   TBranch        *b_GenJet_nDaughters;   //!
   TBranch        *b_GenJet_daughters_status;   //!
   TBranch        *b_GenJet_daughters_pdgid;   //!
   TBranch        *b_GenJet_daughters_E;   //!
   TBranch        *b_GenJet_daughters_pt;   //!
   TBranch        *b_GenJet_daughters_eta;   //!
   TBranch        *b_GenJet_daughters_phi;   //!
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

   PostAnalyzer_MC(TTree *tree=0);
   virtual ~PostAnalyzer_MC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //Self Declared Functions
   virtual Bool_t   PassHLT();
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
   virtual void     PileupReWeighting();
   virtual Double_t PUWeights(Float_t);
   virtual void     PileupReWeighting_XSm5();
   virtual Double_t PUWeights_XSm5(Float_t);
   virtual void     PileupReWeighting_XSp5();
   virtual Double_t PUWeights_XSp5(Float_t);

   //B-Tagging
   virtual Int_t    Getptbin_for_btag(Double_t);
   virtual Double_t GetSFb_CSVL(Double_t, Double_t);
   virtual Double_t GetSFb_CSVM(Double_t, Double_t);
   virtual Double_t GetSFb_CSVT(Double_t, Double_t);
   virtual Double_t GetSFberr_CSVL(Double_t, Double_t);
   virtual Double_t GetSFberr_CSVM(Double_t, Double_t);
   virtual Double_t GetSFberr_CSVT(Double_t, Double_t);
   virtual Double_t GetSFb_CSVL_OldRec7T(Double_t, Double_t); //For old recommendation based on 7TeV Studies (taken from dijet studies, see 
                                                                 //twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods(eg of Method 1c) 
   virtual Double_t GetSFberr_CSVL_OldRec7T(Double_t, Double_t); //For old recommendation got from 7TeV Studies
   virtual Double_t GetSFc_CSVL(Double_t, Double_t);
   virtual Double_t GetSFc_CSVM(Double_t, Double_t);
   virtual Double_t GetSFc_CSVT(Double_t, Double_t);
   virtual Double_t GetSFcerr_CSVL(Double_t, Double_t);
   virtual Double_t GetSFcerr_CSVM(Double_t, Double_t);
   virtual Double_t GetSFcerr_CSVT(Double_t, Double_t);
   virtual Double_t GetSFLight(TString, TString, Double_t, Double_t);
   virtual Double_t GetSFLight_perr(TString, TString, Double_t, Double_t);
   virtual Double_t GetSFLight_nerr(TString, TString, Double_t, Double_t);
   virtual Double_t SFLight(TString, TString, TString, Double_t, Double_t);
   virtual Double_t BTagEventWeight(Double_t, UInt_t);

   virtual void     BookHistograms();

};

#endif

#ifdef PostAnalyzer_MC_cxx
PostAnalyzer_MC::PostAnalyzer_MC(TTree *tree) : fChain(0) 
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
      //To be commented in script
      //-------------------------
      //      chain->Add("/eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_120to170/AODSIM_G_Pt_120to170_164_1_tts.root/myEvent");
      chain->Add("/eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_170to300/AODSIM_G_Pt_170to300_169_1_pd1.root/myEvent");
      chain->Add("/eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_170to300/AODSIM_G_Pt_170to300_35_1_Xsh.root/myEvent");
      //--------------------------
*/
      //Uncomment this part in script
      //-----------------------------
   
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

      //      cout << FullPathInputFile << endl;

        chain->Add(FullPathInputFile+"/myEvent");

        fileNumber++;

      }
      cout << "Total files in this set = " << fileNumber - 1 << endl; 

      //-----------------------------

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PostAnalyzer_MC::~PostAnalyzer_MC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   file->cd();
   file->Write();
   h_CutFlowTable->Write();
   h_CutFlowTableWithWeights->Write();
   h_CutFlowTableWithWeights_XSm5->Write();
   h_CutFlowTableWithWeights_XSp5->Write();
   file->Close();
}

Int_t PostAnalyzer_MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PostAnalyzer_MC::LoadTree(Long64_t entry)
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

void PostAnalyzer_MC::Init(TTree *tree)
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
   GenPhoton_status = 0;
   GenPhoton_E = 0;
   GenPhoton_et = 0;
   GenPhoton_pt = 0;
   GenPhoton_eta = 0;
   GenPhoton_phi = 0;
   GenPhoton_px = 0;
   GenPhoton_py = 0;
   GenPhoton_pz = 0;
   GenPhoton_vx = 0;
   GenPhoton_vy = 0;
   GenPhoton_vz = 0;
   GenPhoton_mother_status = 0;
   GenPhoton_mother_pdgid = 0;
   GenPhoton_mother_E = 0;
   GenPhoton_mother_et = 0;
   GenPhoton_mother_pt = 0;
   GenPhoton_mother_eta = 0;
   GenPhoton_mother_phi = 0;
   GenPhoton_mother_px = 0;
   GenPhoton_mother_py = 0;
   GenPhoton_mother_pz = 0;
   GenPhoton_mother_vx = 0;
   GenPhoton_mother_vy = 0;
   GenPhoton_mother_vz = 0;
   GenPhoton_grandmother_status = 0;
   GenPhoton_grandmother_pdgid = 0;
   GenPhoton_grandmother_E = 0;
   GenPhoton_grandmother_et = 0;
   GenPhoton_grandmother_pt = 0;
   GenPhoton_grandmother_eta = 0;
   GenPhoton_grandmother_phi = 0;
   GenPhoton_grandmother_px = 0;
   GenPhoton_grandmother_py = 0;
   GenPhoton_grandmother_pz = 0;
   GenPhoton_grandmother_vx = 0;
   GenPhoton_grandmother_vy = 0;
   GenPhoton_grandmother_vz = 0;
   hasMatchedGenPhoton_ = 0;
   MatchedGenPhoton_status = 0;
   MatchedGenPhoton_pdgid = 0;
   MatchedGenPhoton_E = 0;
   MatchedGenPhoton_et = 0;
   MatchedGenPhoton_pt = 0;
   MatchedGenPhoton_eta = 0;
   MatchedGenPhoton_phi = 0;
   MatchedGenPhoton_px = 0;
   MatchedGenPhoton_py = 0;
   MatchedGenPhoton_pz = 0;
   MatchedGenPhoton_vx = 0;
   MatchedGenPhoton_vy = 0;
   MatchedGenPhoton_vz = 0;
   MatchedGenPhoton_nMothers = 0;
   MatchedGenPhoton_mothersPdgid = 0;
   MatchedGenPhoton_mothersStatus = 0;
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
   hasMatchedGenParton_ToPFJet_ = 0;
   hasMatchedGenJet_ToPFJet_ = 0;
   MatchedGenPartonToPFJet_pdgId = 0;
   MatchedGenPartonToPFJet_status = 0;
   MatchedGenJetToPFJet_E = 0;
   MatchedGenJetToPFJet_pt = 0;
   MatchedGenJetToPFJet_eta = 0;
   MatchedGenJetToPFJet_phi = 0;
   MatchedGenJetToPFJet_nDaughters = 0;
   PFPatJet_JetPartonFlavor = 0;
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
   GenJet_E = 0;
   GenJet_et = 0;
   GenJet_pt = 0;
   GenJet_eta = 0;
   GenJet_phi = 0;
   GenJet_px = 0;
   GenJet_py = 0;
   GenJet_pz = 0;
   GenJet_vx = 0;
   GenJet_vy = 0;
   GenJet_vz = 0;
   GenJet_nDaughters = 0;
   GenJet_daughters_status = 0;
   GenJet_daughters_pdgid = 0;
   GenJet_daughters_E = 0;
   GenJet_daughters_pt = 0;
   GenJet_daughters_eta = 0;
   GenJet_daughters_phi = 0;
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
   fChain->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
   fChain->SetBranchAddress("npuVerticesm1", &npuVerticesm1, &b_npuVerticesm1);
   fChain->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
   fChain->SetBranchAddress("npuVerticespm2", &npuVerticespm2, &b_npuVerticespm2);
   fChain->SetBranchAddress("trueNumofInteractions", &trueNumofInteractions, &b_trueNumofInteractions);
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
   fChain->SetBranchAddress("GenPhoton_n", &GenPhoton_n, &b_GenPhoton_n);
   fChain->SetBranchAddress("GenPhoton_status", &GenPhoton_status, &b_GenPhoton_status);
   fChain->SetBranchAddress("GenPhoton_E", &GenPhoton_E, &b_GenPhoton_E);
   fChain->SetBranchAddress("GenPhoton_et", &GenPhoton_et, &b_GenPhoton_et);
   fChain->SetBranchAddress("GenPhoton_pt", &GenPhoton_pt, &b_GenPhoton_pt);
   fChain->SetBranchAddress("GenPhoton_eta", &GenPhoton_eta, &b_GenPhoton_eta);
   fChain->SetBranchAddress("GenPhoton_phi", &GenPhoton_phi, &b_GenPhoton_phi);
   fChain->SetBranchAddress("GenPhoton_px", &GenPhoton_px, &b_GenPhoton_px);
   fChain->SetBranchAddress("GenPhoton_py", &GenPhoton_py, &b_GenPhoton_py);
   fChain->SetBranchAddress("GenPhoton_pz", &GenPhoton_pz, &b_GenPhoton_pz);
   fChain->SetBranchAddress("GenPhoton_vx", &GenPhoton_vx, &b_GenPhoton_vx);
   fChain->SetBranchAddress("GenPhoton_vy", &GenPhoton_vy, &b_GenPhoton_vy);
   fChain->SetBranchAddress("GenPhoton_vz", &GenPhoton_vz, &b_GenPhoton_vz);
   fChain->SetBranchAddress("GenPhoton_mother_status", &GenPhoton_mother_status, &b_GenPhoton_mother_status);
   fChain->SetBranchAddress("GenPhoton_mother_pdgid", &GenPhoton_mother_pdgid, &b_GenPhoton_mother_pdgid);
   fChain->SetBranchAddress("GenPhoton_mother_E", &GenPhoton_mother_E, &b_GenPhoton_mother_E);
   fChain->SetBranchAddress("GenPhoton_mother_et", &GenPhoton_mother_et, &b_GenPhoton_mother_et);
   fChain->SetBranchAddress("GenPhoton_mother_pt", &GenPhoton_mother_pt, &b_GenPhoton_mother_pt);
   fChain->SetBranchAddress("GenPhoton_mother_eta", &GenPhoton_mother_eta, &b_GenPhoton_mother_eta);
   fChain->SetBranchAddress("GenPhoton_mother_phi", &GenPhoton_mother_phi, &b_GenPhoton_mother_phi);
   fChain->SetBranchAddress("GenPhoton_mother_px", &GenPhoton_mother_px, &b_GenPhoton_mother_px);
   fChain->SetBranchAddress("GenPhoton_mother_py", &GenPhoton_mother_py, &b_GenPhoton_mother_py);
   fChain->SetBranchAddress("GenPhoton_mother_pz", &GenPhoton_mother_pz, &b_GenPhoton_mother_pz);
   fChain->SetBranchAddress("GenPhoton_mother_vx", &GenPhoton_mother_vx, &b_GenPhoton_mother_vx);
   fChain->SetBranchAddress("GenPhoton_mother_vy", &GenPhoton_mother_vy, &b_GenPhoton_mother_vy);
   fChain->SetBranchAddress("GenPhoton_mother_vz", &GenPhoton_mother_vz, &b_GenPhoton_mother_vz);
   fChain->SetBranchAddress("GenPhoton_grandmother_status", &GenPhoton_grandmother_status, &b_GenPhoton_grandmother_status);
   fChain->SetBranchAddress("GenPhoton_grandmother_pdgid", &GenPhoton_grandmother_pdgid, &b_GenPhoton_grandmother_pdgid);
   fChain->SetBranchAddress("GenPhoton_grandmother_E", &GenPhoton_grandmother_E, &b_GenPhoton_grandmother_E);
   fChain->SetBranchAddress("GenPhoton_grandmother_et", &GenPhoton_grandmother_et, &b_GenPhoton_grandmother_et);
   fChain->SetBranchAddress("GenPhoton_grandmother_pt", &GenPhoton_grandmother_pt, &b_GenPhoton_grandmother_pt);
   fChain->SetBranchAddress("GenPhoton_grandmother_eta", &GenPhoton_grandmother_eta, &b_GenPhoton_grandmother_eta);
   fChain->SetBranchAddress("GenPhoton_grandmother_phi", &GenPhoton_grandmother_phi, &b_GenPhoton_grandmother_phi);
   fChain->SetBranchAddress("GenPhoton_grandmother_px", &GenPhoton_grandmother_px, &b_GenPhoton_grandmother_px);
   fChain->SetBranchAddress("GenPhoton_grandmother_py", &GenPhoton_grandmother_py, &b_GenPhoton_grandmother_py);
   fChain->SetBranchAddress("GenPhoton_grandmother_pz", &GenPhoton_grandmother_pz, &b_GenPhoton_grandmother_pz);
   fChain->SetBranchAddress("GenPhoton_grandmother_vx", &GenPhoton_grandmother_vx, &b_GenPhoton_grandmother_vx);
   fChain->SetBranchAddress("GenPhoton_grandmother_vy", &GenPhoton_grandmother_vy, &b_GenPhoton_grandmother_vy);
   fChain->SetBranchAddress("GenPhoton_grandmother_vz", &GenPhoton_grandmother_vz, &b_GenPhoton_grandmother_vz);
   fChain->SetBranchAddress("hasMatchedGenPhoton_", &hasMatchedGenPhoton_, &b_hasMatchedGenPhoton_);
   fChain->SetBranchAddress("MatchedGenPhoton_status", &MatchedGenPhoton_status, &b_MatchedGenPhoton_status);
   fChain->SetBranchAddress("MatchedGenPhoton_pdgid", &MatchedGenPhoton_pdgid, &b_MatchedGenPhoton_pdgid);
   fChain->SetBranchAddress("MatchedGenPhoton_E", &MatchedGenPhoton_E, &b_MatchedGenPhoton_E);
   fChain->SetBranchAddress("MatchedGenPhoton_et", &MatchedGenPhoton_et, &b_MatchedGenPhoton_et);
   fChain->SetBranchAddress("MatchedGenPhoton_pt", &MatchedGenPhoton_pt, &b_MatchedGenPhoton_pt);
   fChain->SetBranchAddress("MatchedGenPhoton_eta", &MatchedGenPhoton_eta, &b_MatchedGenPhoton_eta);
   fChain->SetBranchAddress("MatchedGenPhoton_phi", &MatchedGenPhoton_phi, &b_MatchedGenPhoton_phi);
   fChain->SetBranchAddress("MatchedGenPhoton_px", &MatchedGenPhoton_px, &b_MatchedGenPhoton_px);
   fChain->SetBranchAddress("MatchedGenPhoton_py", &MatchedGenPhoton_py, &b_MatchedGenPhoton_py);
   fChain->SetBranchAddress("MatchedGenPhoton_pz", &MatchedGenPhoton_pz, &b_MatchedGenPhoton_pz);
   fChain->SetBranchAddress("MatchedGenPhoton_vx", &MatchedGenPhoton_vx, &b_MatchedGenPhoton_vx);
   fChain->SetBranchAddress("MatchedGenPhoton_vy", &MatchedGenPhoton_vy, &b_MatchedGenPhoton_vy);
   fChain->SetBranchAddress("MatchedGenPhoton_vz", &MatchedGenPhoton_vz, &b_MatchedGenPhoton_vz);
   fChain->SetBranchAddress("MatchedGenPhoton_nMothers", &MatchedGenPhoton_nMothers, &b_MatchedGenPhoton_nMothers);
   fChain->SetBranchAddress("MatchedGenPhoton_mothersPdgid", &MatchedGenPhoton_mothersPdgid, &b_MatchedGenPhoton_mothersPdgid);
   fChain->SetBranchAddress("MatchedGenPhoton_mothersStatus", &MatchedGenPhoton_mothersStatus, &b_MatchedGenPhoton_mothersStatus);
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
   fChain->SetBranchAddress("hasMatchedGenParton_ToPFJet_", &hasMatchedGenParton_ToPFJet_, &b_hasMatchedGenParton_ToPFJet_);
   fChain->SetBranchAddress("hasMatchedGenJet_ToPFJet_", &hasMatchedGenJet_ToPFJet_, &b_hasMatchedGenJet_ToPFJet_);
   fChain->SetBranchAddress("MatchedGenPartonToPFJet_pdgId", &MatchedGenPartonToPFJet_pdgId, &b_MatchedGenPartonToPFJet_pdgId);
   fChain->SetBranchAddress("MatchedGenPartonToPFJet_status", &MatchedGenPartonToPFJet_status, &b_MatchedGenPartonToPFJet_status);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_E", &MatchedGenJetToPFJet_E, &b_MatchedGenJetToPFJet_E);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_pt", &MatchedGenJetToPFJet_pt, &b_MatchedGenJetToPFJet_pt);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_eta", &MatchedGenJetToPFJet_eta, &b_MatchedGenJetToPFJet_eta);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_phi", &MatchedGenJetToPFJet_phi, &b_MatchedGenJetToPFJet_phi);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_nDaughters", &MatchedGenJetToPFJet_nDaughters, &b_MatchedGenJetToPFJet_nDaughters);
   fChain->SetBranchAddress("PFPatJet_JetPartonFlavor", &PFPatJet_JetPartonFlavor, &b_PFPatJet_JetPartonFlavor);
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
   fChain->SetBranchAddress("GenJet_n", &GenJet_n, &b_GenJet_n);
   fChain->SetBranchAddress("GenJet_E", &GenJet_E, &b_GenJet_E);
   fChain->SetBranchAddress("GenJet_et", &GenJet_et, &b_GenJet_et);
   fChain->SetBranchAddress("GenJet_pt", &GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("GenJet_eta", &GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_phi", &GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_px", &GenJet_px, &b_GenJet_px);
   fChain->SetBranchAddress("GenJet_py", &GenJet_py, &b_GenJet_py);
   fChain->SetBranchAddress("GenJet_pz", &GenJet_pz, &b_GenJet_pz);
   fChain->SetBranchAddress("GenJet_vx", &GenJet_vx, &b_GenJet_vx);
   fChain->SetBranchAddress("GenJet_vy", &GenJet_vy, &b_GenJet_vy);
   fChain->SetBranchAddress("GenJet_vz", &GenJet_vz, &b_GenJet_vz);
   fChain->SetBranchAddress("GenJet_nDaughters", &GenJet_nDaughters, &b_GenJet_nDaughters);
   fChain->SetBranchAddress("GenJet_daughters_status", &GenJet_daughters_status, &b_GenJet_daughters_status);
   fChain->SetBranchAddress("GenJet_daughters_pdgid", &GenJet_daughters_pdgid, &b_GenJet_daughters_pdgid);
   fChain->SetBranchAddress("GenJet_daughters_E", &GenJet_daughters_E, &b_GenJet_daughters_E);
   fChain->SetBranchAddress("GenJet_daughters_pt", &GenJet_daughters_pt, &b_GenJet_daughters_pt);
   fChain->SetBranchAddress("GenJet_daughters_eta", &GenJet_daughters_eta, &b_GenJet_daughters_eta);
   fChain->SetBranchAddress("GenJet_daughters_phi", &GenJet_daughters_phi, &b_GenJet_daughters_phi);
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

Bool_t PostAnalyzer_MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PostAnalyzer_MC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PostAnalyzer_MC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t PostAnalyzer_MC::PassHLT(){
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

Bool_t PostAnalyzer_MC::NonScrapingEvt(){
  Bool_t scrapingEvt = !IsScrapingEvent_;
  return scrapingEvt;
}

Int_t PostAnalyzer_MC::GoodPrimaryVtx(){
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

Bool_t PostAnalyzer_MC::ResSpikes(Int_t i){
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

Double_t PostAnalyzer_MC::GetLICTD(Int_t i){

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

Bool_t PostAnalyzer_MC::ResSpikes_NoLICTD(Int_t i){
  Bool_t spikes = false;
  if( (fabs((*Photon_xtal_timing)[i][0])) < 3.0 &&   //(*Photon_xtal_timing)[i][0]) is the time of arrival for ith photon and seed crystal
      (*Photon_SigmaIEtaIEta)[i] > 0.001        &&
      (*Photon_SigmaIPhiIPhi)[i] > 0.001        &&
      (*Photon_r9)[i] < 1.0){
    spikes = true;
  }
  return spikes;
}

Double_t PostAnalyzer_MC::EAChargedHadrons(Double_t eta){
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

Double_t PostAnalyzer_MC::EANeutralHadrons(Double_t eta){
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

Double_t PostAnalyzer_MC::EAPhotons(Double_t eta){
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

Bool_t PostAnalyzer_MC::TightPhotonIdnPFIso(Int_t i){
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

Bool_t PostAnalyzer_MC::MediumPhotonIdnPFIso(Int_t i){
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

Bool_t PostAnalyzer_MC::LoosePhotonIdnPFIso(Int_t i){
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


Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithTightID(){
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

Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithMediumID(){
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

Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithLooseID(){
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


Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithTightID_NoLICTD(){
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

Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithMediumID_NoLICTD(){
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

Int_t PostAnalyzer_MC::GetPhotonPassingAllCuts_WithLooseID_NoLICTD(){
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



Int_t PostAnalyzer_MC::GetJetPassingIDnMatchedToPhoton(Int_t pc){
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

Bool_t PostAnalyzer_MC::TightJetId(Int_t i){
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

Double_t PostAnalyzer_MC::GetdEta(Double_t eta1, Double_t eta2){

  Double_t dEta = fabs(eta1 - eta2);
  return dEta;
}

Double_t PostAnalyzer_MC::GetdPhi(Double_t phi1, Double_t phi2){

  Double_t dphi = (phi1 - phi2);
  Double_t twoPi = 2.0*(TMath::Pi());

  if(dphi < 0) dphi = - dphi;
  if(dphi >= (twoPi - dphi)) dphi = twoPi - dphi;

  return dphi;
}

Double_t PostAnalyzer_MC::GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

  Double_t dEta = GetdEta(eta1, eta2);
  Double_t dPhi = GetdPhi(phi1, phi2);

  Double_t dR = 0.0;
  dR = sqrt(dEta*dEta + dPhi*dPhi);

  return dR;
}

Double_t PostAnalyzer_MC::GetInvtMass(Int_t ph, Int_t jet){

  Double_t mass = 0.0;

  Double_t E  = (*Photon_E)[ph] + (*PFPatJet_E)[jet];
  Double_t PX = (*Photon_px)[ph] + (*PFPatJet_px)[jet];
  Double_t PY = (*Photon_py)[ph] + (*PFPatJet_py)[jet];
  Double_t PZ = (*Photon_pz)[ph] + (*PFPatJet_pz)[jet];

  mass = sqrt(E*E - PX*PX - PY*PY - PZ*PZ);

  return mass;
}

Bool_t PostAnalyzer_MC::PassCSVLBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.244){
    passBtag = true;
  }
  return passBtag;
}

Bool_t PostAnalyzer_MC::PassCSVMBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.679){
    passBtag = true;
  }
  return passBtag;
}

Bool_t PostAnalyzer_MC::PassCSVTBTag(Int_t jc){
  Bool_t passBtag = false;
  if((*PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff)[jc] > 0.898){
    passBtag = true;
  }
  return passBtag;
}

void PostAnalyzer_MC::PileupReWeighting(){

  TFile *fData = TFile::Open("/uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/PostAnalyzer_MC/PileupHistograms/Data_Run2012ABCD_PileUpDist/Data_Run2012ABCD_True60_PileupHistogram.root");
  TH1F *dataPU = (TH1F*)fData->Get("pileup");

  std::vector<float> DataPileUp;
  DataPileUp.clear();
  for(Int_t i = 0; i < 60; i++){
    DataPileUp.push_back(dataPU->GetBinContent(i+1));
  }

  for(Int_t ibin = 0; ibin < 60; ibin++){
    h_DataPUDist->SetBinContent(ibin+1, DataPileUp[ibin]); //This to get in output
    h_DataPUNormDist->SetBinContent(ibin+1, DataPileUp[ibin]); //This to get normalized distribution in output
    h_PUScaleFactor->SetBinContent(ibin+1, DataPileUp[ibin]); //This to use in weight calculation
  }

  Double_t Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06}; //Got this from Pileup_MC_Gen_Scenarios twiki

  TH1F *h_MCWeights = new TH1F("h_MCWeights", "MC PileUp Weights", 60, 0, 60);
  for(Int_t i = 0; i < 60; i++){
    h_MCPUNormDist->SetBinContent(i+1, Summer2012_S10[i]); //This to get in output
    h_MCWeights->SetBinContent(i+1, Summer2012_S10[i]); //This to be used in weights calculation
  }

  h_DataPUNormDist->Scale(1.0/h_DataPUNormDist->Integral());
  h_PUScaleFactor->Scale(1.0/h_PUScaleFactor->Integral());

  h_PUScaleFactor->Divide(h_MCWeights);

}

Double_t PostAnalyzer_MC::PUWeights(Float_t npv){
  Int_t bin = h_PUScaleFactor->GetXaxis()->FindBin( npv );
  return h_PUScaleFactor->GetBinContent( bin );
}

//Note:It can be seen on checking the PU weighting for each of XS, XS -5% and XS +5% that on increasing XS, the number of interactions increases.
//and hence the PU of data will increase for higher number of interactions but decrease for low interactions.
//On doing cout of PUWeights for all the three, it can be seen that upto the tureNumOfInt ~ 25, the PU wt increase in the order ps < 0 < m5 and then 
//reverses for high num of interactions. ANd since PU_Wt = data_PU/MC_PU so it means large pu in data for high XS towards high num of interactions.

void PostAnalyzer_MC::PileupReWeighting_XSm5(){

  TFile *fData = TFile::Open("/uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/PostAnalyzer_MC/PileupHistograms/Data_Run2012ABCD_PileUpDist_XSm5/Data_Run2012ABCD_True60_PileupHistogram_XSm5.root");
  TH1F *dataPU = (TH1F*)fData->Get("pileup");

  std::vector<float> DataPileUp;
  DataPileUp.clear();
  for(Int_t i = 0; i < 60; i++){
    DataPileUp.push_back(dataPU->GetBinContent(i+1));
  }

  for(Int_t ibin = 0; ibin < 60; ibin++){
    h_PUScaleFactor_XSm5->SetBinContent(ibin+1, DataPileUp[ibin]); 
  }

  Double_t Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06}; //Got this from Pileup_MC_Gen_Scenarios twiki

  TH1F *h_MCWeights = new TH1F("h_MCWeights", "MC PileUp Weights", 60, 0, 60);
  for(Int_t i = 0; i < 60; i++){
    h_MCWeights->SetBinContent(i+1, Summer2012_S10[i]); 
  }

  h_PUScaleFactor_XSm5->Scale(1.0/h_PUScaleFactor_XSm5->Integral());

  h_PUScaleFactor_XSm5->Divide(h_MCWeights);

}

Double_t PostAnalyzer_MC::PUWeights_XSm5(Float_t npv){
  Int_t bin = h_PUScaleFactor_XSm5->GetXaxis()->FindBin( npv );
  return h_PUScaleFactor_XSm5->GetBinContent( bin );
}

void PostAnalyzer_MC::PileupReWeighting_XSp5(){

  TFile *fData = TFile::Open("/uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/PostAnalyzer_MC/PileupHistograms/Data_Run2012ABCD_PileUpDist_XSp5/Data_Run2012ABCD_True60_PileupHistogram_XSp5.root");
  TH1F *dataPU = (TH1F*)fData->Get("pileup");

  std::vector<float> DataPileUp;
  DataPileUp.clear();
  for(Int_t i = 0; i < 60; i++){
    DataPileUp.push_back(dataPU->GetBinContent(i+1));
  }

  for(Int_t ibin = 0; ibin < 60; ibin++){
    h_PUScaleFactor_XSp5->SetBinContent(ibin+1, DataPileUp[ibin]); 
  }

  Double_t Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06}; //Got this from Pileup_MC_Gen_Scenarios twiki

  TH1F *h_MCWeights = new TH1F("h_MCWeights", "MC PileUp Weights", 60, 0, 60);
  for(Int_t i = 0; i < 60; i++){
    h_MCWeights->SetBinContent(i+1, Summer2012_S10[i]); 
  }

  h_PUScaleFactor_XSp5->Scale(1.0/h_PUScaleFactor_XSp5->Integral());

  h_PUScaleFactor_XSp5->Divide(h_MCWeights);

}

Double_t PostAnalyzer_MC::PUWeights_XSp5(Float_t npv){
  Int_t bin = h_PUScaleFactor_XSp5->GetXaxis()->FindBin( npv );
  return h_PUScaleFactor_XSp5->GetBinContent( bin );
}


//B-Tagging Functions
Int_t PostAnalyzer_MC::Getptbin_for_btag(Double_t pt){
  if(pt<30) return 0;
  else if(pt<40) return 1;
  else if(pt<50) return 2;
  else if(pt<60) return 3;
  else if(pt<70) return 4;
  else if(pt<80) return 5;
  else if(pt<100) return 6;
  else if(pt<120) return 7;
  else if(pt<160) return 8;
  else if(pt<210) return 9;
  else if(pt<260) return 10;
  else if(pt<320) return 11;
  else if(pt<400) return 12;
  else if(pt<500) return 13;
  else if(pt<600) return 14;
  else return 15;
  
}

Double_t PostAnalyzer_MC::GetSFb_CSVL(Double_t jetPt, Double_t jetEta){
  double pt = jetPt;
  // for scale factor extrapolation
  if(pt<20) pt = 20;
  if(pt>800) pt = 800;

  Double_t SF;
  if(fabs(jetEta) <= 2.4){
    SF = 0.997942*((1.+(0.00923753*pt))/(1.+(0.0096119*pt)));
  }
  else{
    SF = 1.0;
  }
  return SF;
}

Double_t PostAnalyzer_MC::GetSFc_CSVL(Double_t jetPt, Double_t jetEta){

  return GetSFb_CSVL(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFberr_CSVL(Double_t jetPt, Double_t jetEta){

  double pt[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};

  double SFerr[] = {0.033299, 0.0146768, 0.013803, 0.0170145, 0.0166976, 0.0137879, 0.0149072, 0.0153068, 0.0133077, 0.0123737, 0.0157152, 0.0175161, 0.0209241, 0.0278605, 0.0346928, 0.0350099};

  int ptBin = Getptbin_for_btag(jetPt);

  Double_t sfErr;
  if(fabs(jetEta) <= 2.4){
    sfErr = SFerr[ptBin];
    if(jetPt < 20 || jetPt > 800) sfErr = 2*sfErr;
  }
  else{
    sfErr = 0.0;
  }

  return sfErr;

}

Double_t PostAnalyzer_MC::GetSFcerr_CSVL(Double_t jetPt, Double_t jetEta){

  return 2*GetSFberr_CSVL(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFb_CSVM(Double_t jetPt, Double_t jetEta){
  double pt = jetPt;
  // for scale factor extrapolation
  if(pt<20) pt = 20;
  if(pt>800) pt = 800;

  Double_t SF;
  if(fabs(jetEta) <= 2.4){
    SF = (0.938887+(0.00017124*pt))+(-2.76366e-07*(pt*pt));
  }
  else{
    SF = 1.0;
  }
  return SF;
}

Double_t PostAnalyzer_MC::GetSFc_CSVM(Double_t jetPt, Double_t jetEta){

  return GetSFb_CSVM(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFberr_CSVM(Double_t jetPt, Double_t jetEta){

  double pt[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};

  double SFerr[] = {0.0415707, 0.0204209, 0.0223227, 0.0206655, 0.0199325, 0.0174121, 0.0202332, 0.0182446, 0.0159777, 0.0218531, 0.0204688, 0.0265191, 0.0313175, 0.0415417, 0.0740446, 0.0596716};

  int ptBin = Getptbin_for_btag(jetPt);

  Double_t sfErr;
  if(fabs(jetEta) <= 2.4){
    sfErr = SFerr[ptBin];
    if(jetPt < 20 || jetPt > 800) sfErr = 2*sfErr;
  }
  else{
    sfErr = 0.0;
  }

  return sfErr;

}

Double_t PostAnalyzer_MC::GetSFcerr_CSVM(Double_t jetPt, Double_t jetEta){

  return 2*GetSFberr_CSVM(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFb_CSVT(Double_t jetPt, Double_t jetEta){
  double pt = jetPt;
  // for scale factor extrapolation
  if(pt<20) pt = 20;
  if(pt>800) pt = 800;

  Double_t SF;
  if(fabs(jetEta) <= 2.4){
    SF = (0.927563+(1.55479e-05*pt))+(-1.90666e-07*(pt*pt));
  }
  else{
    SF = 1.0;
  }
  return SF;
}

Double_t PostAnalyzer_MC::GetSFc_CSVT(Double_t jetPt, Double_t jetEta){

  return GetSFb_CSVT(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFberr_CSVT(Double_t jetPt, Double_t jetEta){

  double pt[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};

  double SFerr[] = {0.0515703, 0.0264008, 0.0272757, 0.0275565, 0.0248745, 0.0218456, 0.0253845, 0.0239588, 0.0271791, 0.0273912, 0.0379822, 0.0411624, 0.0786307, 0.0866832, 0.0942053, 0.102403};

  int ptBin = Getptbin_for_btag(jetPt);

  Double_t sfErr;
  if(fabs(jetEta) <= 2.4){
    sfErr = SFerr[ptBin];
    if(jetPt < 20 || jetPt > 800) sfErr = 2*sfErr;
  }
  else{
    sfErr = 0.0;
  }

  return sfErr;

}

Double_t PostAnalyzer_MC::GetSFcerr_CSVT(Double_t jetPt, Double_t jetEta){

  return 2*GetSFberr_CSVT(jetPt, jetEta);
}

Double_t PostAnalyzer_MC::GetSFb_CSVL_OldRec7T(Double_t jetPt, Double_t jetEta){
  double pt = jetPt;
  // for scale factor extrapolation
  if(pt<30) pt = 30;
  if(pt>670) pt = 670;

  Double_t SF;
  if(fabs(jetEta) <= 2.4){
    SF = 1.02658*((1.+(0.0195388*pt))/(1.+(0.0209145*pt)));
  }
  else{
    SF = 1.0;
  }
  return SF;
}

Double_t PostAnalyzer_MC::GetSFberr_CSVL_OldRec7T(Double_t jetPt, Double_t jetEta){

  double pt[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};

  double SFerr[] = {0.12, 0.0188743, 0.0161816, 0.0139824, 0.0152644, 0.0161226, 0.0157396, 0.0161619, 0.0168747, 0.0257175, 0.026424, 0.0264928, 0.0315127, 0.030734, 0.0438259};

  int ptBin;
  if(jetPt < 500){
    ptBin = Getptbin_for_btag(jetPt);
  }
  else{
    ptBin = 14;
  }

  Double_t sfErr;
  if(fabs(jetEta) <= 2.4){
    sfErr = SFerr[ptBin];
    if(jetPt > 670) sfErr = 2*sfErr;
  }
  else{
    sfErr = 0.0;
  }

  return sfErr;

}

Double_t PostAnalyzer_MC::GetSFLight(TString tagger, TString TaggerStrength, Double_t JetPt, Double_t JetEta){

  return SFLight("mean", tagger, TaggerStrength, JetPt, JetEta);

}

Double_t PostAnalyzer_MC::GetSFLight_perr(TString tagger, TString TaggerStrength, Double_t JetPt, Double_t JetEta){

  Double_t sferr = 0.0;

  Double_t ptmax;
  if((TaggerStrength == "L" && fabs(JetEta) >= 1.5 && fabs(JetEta) <= 2.4) || (TaggerStrength == "M" && fabs(JetEta) >= 1.6 && fabs(JetEta) <= 2.4)){
    ptmax = 850;
  }
  else{
    ptmax = 1000;
  }

  sferr = SFLight("max", tagger, TaggerStrength, JetPt, JetEta) - SFLight("mean", tagger, TaggerStrength, JetPt, JetEta);

  if(JetPt > ptmax) sferr = 2*sferr;

  return sferr;
}

Double_t PostAnalyzer_MC::GetSFLight_nerr(TString tagger, TString TaggerStrength, Double_t JetPt, Double_t JetEta){

  Double_t sferr = 0.0;

  Double_t ptmax;
  if((TaggerStrength == "L" && fabs(JetEta) >= 1.5 && fabs(JetEta) <= 2.4) || (TaggerStrength == "M" && fabs(JetEta) >= 1.6 && fabs(JetEta) <= 2.4)){
    ptmax = 850;
  }
  else{
    ptmax = 1000;
  }

  sferr = SFLight("mean", tagger, TaggerStrength, JetPt, JetEta) - SFLight("min", tagger, TaggerStrength, JetPt, JetEta);

  if(JetPt > ptmax) sferr = 2*sferr;

  return sferr;
}

Double_t PostAnalyzer_MC::SFLight(TString meanminmax, TString tagger, TString TaggerStrength, Double_t JetPt, Double_t JetEta){

  Double_t ptmax;
  if((TaggerStrength == "L" && fabs(JetEta) >= 1.5 && fabs(JetEta) <= 2.4) || (TaggerStrength == "M" && fabs(JetEta) >= 1.6 && fabs(JetEta) <= 2.4)){
    ptmax = 850;
  }
  else{
    ptmax = 1000;
  }

  Double_t sflight = 1.0;
  Double_t Etamin = 0.0;
  Double_t Etamax = 0.0;

  if(fabs(JetEta) <= 2.4){

    if(TaggerStrength == "L" && fabs(JetEta) >= 0.0 && fabs(JetEta) <= 0.5){
      Etamin = 0.0;
      Etamax = 0.5;
    }
    if(TaggerStrength == "L" && fabs(JetEta) >= 0.5 && fabs(JetEta) <= 1.0){
      Etamin = 0.5;
      Etamax = 1.0;
    }
    if(TaggerStrength == "L" && fabs(JetEta) >= 1.0 && fabs(JetEta) <= 1.5){
      Etamin = 1.0;
      Etamax = 1.5;
    }
    if(TaggerStrength == "L" && fabs(JetEta) >= 1.5 && fabs(JetEta) <= 2.4){
      Etamin = 1.5;
      Etamax = 2.4;
    }
    if(TaggerStrength == "M" && fabs(JetEta) >= 0.0 && fabs(JetEta) <= 0.8){
      Etamin = 0.0;
      Etamax = 0.8;
    }
    if(TaggerStrength == "M" && fabs(JetEta) >= 0.8 && fabs(JetEta) <= 1.6){
      Etamin = 0.8;
      Etamax = 1.6;
    }
    if(TaggerStrength == "M" && fabs(JetEta) >= 1.6 && fabs(JetEta) <= 2.4){
      Etamin = 1.6;
      Etamax = 2.4;
    }
    if(TaggerStrength == "T" && fabs(JetEta) >= 0.0 && fabs(JetEta) <= 2.4){
      Etamin = 0.0;
      Etamax = 2.4;
    }
      
    if(JetPt > ptmax) JetPt = ptmax;

    Double_t x = JetPt;

    if(tagger == "CSV", TaggerStrength == "L", Etamin == 0.0, Etamax == 0.5){
      if(meanminmax == "mean") sflight = ((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "L", Etamin == 0.5, Etamax == 1.0){
      if(meanminmax == "mean") sflight = ((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "L", Etamin == 1.0, Etamax == 1.5){
      if(meanminmax == "mean") sflight = ((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "L", Etamin == 1.5, Etamax == 2.4){
      if(meanminmax == "mean") sflight = ((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "M", Etamin == 0.0, Etamax == 0.8){
      if(meanminmax == "mean") sflight = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "M", Etamin == 0.8, Etamax == 1.6){
      if(meanminmax == "mean") sflight = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "M", Etamin == 1.6, Etamax == 2.4){
      if(meanminmax == "mean") sflight = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
    }
    if(tagger == "CSV", TaggerStrength == "T", Etamin == 0.0, Etamax == 2.4){
      if(meanminmax == "mean") sflight = ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
      if(meanminmax == "min") sflight = ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
      if(meanminmax == "max") sflight = ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));
    }

  }
  else{
    sflight = 1.0;
  }

  return sflight;
}

Double_t PostAnalyzer_MC::BTagEventWeight(Double_t ScaleFactor, UInt_t nBTags){

  if( nBTags > 1 )
    {
      cout << "Only one leading jet is considered. Hence, the number of b-tags cannot exceed 1." << endl;     
    }

  /*
    ##################################################################
    Event weight matrix:
    ------------------------------------------------------------------
    nBTags\b-tagged jets  |    0        1             2
    ------------------------------------------------------------------
      0                   |    1      1-SF      (1-SF1)(1-SF2)
                          |
      1                   |    0       SF    SF1(1-SF2)+(1-SF1)SF2
                          |
      2                   |    0        0           SF1SF2
    ##################################################################
    Here
    nBTags = No. of expected b jets from MC truth information
    b-tagged jets = Actual no. of b tagged jets by the discriminator
  */

  double weight = 0;
  double SF = ScaleFactor;

  for(unsigned int i = 0; i <= 1; ++i)
    {
      if( i != nBTags ) continue;

      weight += pow(SF,i)*pow(1-SF,1-i);
    }

  return weight;
}








void PostAnalyzer_MC::BookHistograms(){
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
  h_GJetInvtMass_UnitBinning->GetYaxis()->SetTitle("Events/UnitBin");     h_GJetInvtMass_UnitBinning->GetYaxis()->CenterTitle();
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

  h_CSVL_BJetPt_noErr_0bTag = new TH1F("h_CSVL_BJetPt_noErr_0bTag", "Pt Distribution of BJets passing CSVL (noErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_noErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_noErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_noErr_0bTag->Sumw2();

  h_CSVL_BJetEta_noErr_0bTag = new TH1F("h_CSVL_BJetEta_noErr_0bTag", "Eta Distribution of BJets passing CSVL (noErr_0bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_noErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_noErr_0bTag->Sumw2();

  h_CSVL_BJetPhi_noErr_0bTag = new TH1F("h_CSVL_BJetPhi_noErr_0bTag", "Phi Distribution of BJets passing CSVL (noErr_0bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_noErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_noErr_0bTag->Sumw2();

  h_CSVL_BJetPt_noErr_1bTag = new TH1F("h_CSVL_BJetPt_noErr_1bTag", "Pt Distribution of BJets passing CSVL (noErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_noErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_noErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_noErr_1bTag->Sumw2();

  h_CSVL_BJetEta_noErr_1bTag = new TH1F("h_CSVL_BJetEta_noErr_1bTag", "Eta Distribution of BJets passing CSVL (noErr_1bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_noErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_noErr_1bTag->Sumw2();

  h_CSVL_BJetPhi_noErr_1bTag = new TH1F("h_CSVL_BJetPhi_noErr_1bTag", "Phi Distribution of BJets passing CSVL (noErr_1bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_noErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_noErr_1bTag->Sumw2();

  h_CSVL_BJetPt_pErr_0bTag = new TH1F("h_CSVL_BJetPt_pErr_0bTag", "Pt Distribution of BJets passing CSVL (pErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_pErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_pErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_pErr_0bTag->Sumw2();

  h_CSVL_BJetEta_pErr_0bTag = new TH1F("h_CSVL_BJetEta_pErr_0bTag", "Eta Distribution of BJets passing CSVL (pErr_0bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_pErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_pErr_0bTag->Sumw2();

  h_CSVL_BJetPhi_pErr_0bTag = new TH1F("h_CSVL_BJetPhi_pErr_0bTag", "Phi Distribution of BJets passing CSVL (pErr_0bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_pErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_pErr_0bTag->Sumw2();

  h_CSVL_BJetPt_pErr_1bTag = new TH1F("h_CSVL_BJetPt_pErr_1bTag", "Pt Distribution of BJets passing CSVL (pErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVL_BJetPt_pErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVL_BJetPt_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPt_pErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVL_BJetPt_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPt_pErr_1bTag->Sumw2();

  h_CSVL_BJetEta_pErr_1bTag = new TH1F("h_CSVL_BJetEta_pErr_1bTag", "Eta Distribution of BJets passing CSVL (pErr_1bTag)", 200, -5.0, 5.0);
  h_CSVL_BJetEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetEta_pErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVL_BJetEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetEta_pErr_1bTag->Sumw2();

  h_CSVL_BJetPhi_pErr_1bTag = new TH1F("h_CSVL_BJetPhi_pErr_1bTag", "Phi Distribution of BJets passing CSVL (pErr_1bTag)", 200, 0.0, 10.0);
  h_CSVL_BJetPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVL_BJetPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJetPhi_pErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVL_BJetPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJetPhi_pErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_noErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_noErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_noErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_0bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_0bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_0bTag->Sumw2();

  h_CSVL_GBJetdEta_noErr_0bTag = new TH1F("h_CSVL_GBJetdEta_noErr_0bTag", "DeltaEta Dist Between Gamma & CSVL-BJet (noErr_0bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_noErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_noErr_0bTag->Sumw2();

  h_CSVL_GBJetdPhi_noErr_0bTag = new TH1F("h_CSVL_GBJetdPhi_noErr_0bTag", "DeltaPhi Dist Between Gamma & CSVL-BJet (noErr_0bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_noErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_noErr_0bTag->Sumw2();

  h_CSVL_GBJetdR_noErr_0bTag = new TH1F("h_CSVL_GBJetdR_noErr_0bTag", "DeltaR Dist Between Gamma & CSVL-BJet (noErr_0bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_noErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_noErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_noErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_noErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_noErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_noErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_1bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_noErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (noErr_1bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_noErr_1bTag->Sumw2();

  h_CSVL_GBJetdEta_noErr_1bTag = new TH1F("h_CSVL_GBJetdEta_noErr_1bTag", "DeltaEta Dist Between Gamma & CSVL-BJet (noErr_1bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_noErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_noErr_1bTag->Sumw2();

  h_CSVL_GBJetdPhi_noErr_1bTag = new TH1F("h_CSVL_GBJetdPhi_noErr_1bTag", "DeltaPhi Dist Between Gamma & CSVL-BJet (noErr_1bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_noErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_noErr_1bTag->Sumw2();

  h_CSVL_GBJetdR_noErr_1bTag = new TH1F("h_CSVL_GBJetdR_noErr_1bTag", "DeltaR Dist Between Gamma & CSVL-BJet (noErr_1bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_noErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_noErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_noErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_pErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_pErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_pErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_0bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_0bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_0bTag->Sumw2();

  h_CSVL_GBJetdEta_pErr_0bTag = new TH1F("h_CSVL_GBJetdEta_pErr_0bTag", "DeltaEta Dist Between Gamma & CSVL-BJet (pErr_0bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_pErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_pErr_0bTag->Sumw2();

  h_CSVL_GBJetdPhi_pErr_0bTag = new TH1F("h_CSVL_GBJetdPhi_pErr_0bTag", "DeltaPhi Dist Between Gamma & CSVL-BJet (pErr_0bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_pErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_pErr_0bTag->Sumw2();

  h_CSVL_GBJetdR_pErr_0bTag = new TH1F("h_CSVL_GBJetdR_pErr_0bTag", "DeltaR Dist Between Gamma & CSVL-BJet (pErr_0bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_pErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_pErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_pErr_0bTag->Sumw2();

  h_CSVL_GBJetInvtMass_binning40_pErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_binning40_pErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_binning40_pErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_1bTag)", nbins, xbins);
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_VariableBinning_pErr_1bTag->Sumw2();

  h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag = new TH1F("h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVL-BJet (pErr_1bTag)", 14000, 0, 14000);
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetInvtMass_UnitBinning_pErr_1bTag->Sumw2();

  h_CSVL_GBJetdEta_pErr_1bTag = new TH1F("h_CSVL_GBJetdEta_pErr_1bTag", "DeltaEta Dist Between Gamma & CSVL-BJet (pErr_1bTag)", 120, 0.0, 6.0);
  h_CSVL_GBJetdEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdEta_pErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVL_GBJetdEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdEta_pErr_1bTag->Sumw2();

  h_CSVL_GBJetdPhi_pErr_1bTag = new TH1F("h_CSVL_GBJetdPhi_pErr_1bTag", "DeltaPhi Dist Between Gamma & CSVL-BJet (pErr_1bTag)", 64, 0.0, 3.2);
  h_CSVL_GBJetdPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVL_GBJetdPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_pErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVL_GBJetdPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdPhi_pErr_1bTag->Sumw2();

  h_CSVL_GBJetdR_pErr_1bTag = new TH1F("h_CSVL_GBJetdR_pErr_1bTag", "DeltaR Dist Between Gamma & CSVL-BJet (pErr_1bTag)", 100, 0.0, 10.0);
  h_CSVL_GBJetdR_pErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVL_GBJetdR_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_GBJetdR_pErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVL_GBJetdR_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_GBJetdR_pErr_1bTag->Sumw2();

  h_CSVM_BJetPt_noErr_0bTag = new TH1F("h_CSVM_BJetPt_noErr_0bTag", "Pt Distribution of BJets passing CSVM (noErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_noErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_noErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_noErr_0bTag->Sumw2();

  h_CSVM_BJetEta_noErr_0bTag = new TH1F("h_CSVM_BJetEta_noErr_0bTag", "Eta Distribution of BJets passing CSVM (noErr_0bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_noErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_noErr_0bTag->Sumw2();

  h_CSVM_BJetPhi_noErr_0bTag = new TH1F("h_CSVM_BJetPhi_noErr_0bTag", "Phi Distribution of BJets passing CSVM (noErr_0bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_noErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_noErr_0bTag->Sumw2();

  h_CSVM_BJetPt_noErr_1bTag = new TH1F("h_CSVM_BJetPt_noErr_1bTag", "Pt Distribution of BJets passing CSVM (noErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_noErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_noErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_noErr_1bTag->Sumw2();

  h_CSVM_BJetEta_noErr_1bTag = new TH1F("h_CSVM_BJetEta_noErr_1bTag", "Eta Distribution of BJets passing CSVM (noErr_1bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_noErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_noErr_1bTag->Sumw2();

  h_CSVM_BJetPhi_noErr_1bTag = new TH1F("h_CSVM_BJetPhi_noErr_1bTag", "Phi Distribution of BJets passing CSVM (noErr_1bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_noErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_noErr_1bTag->Sumw2();

  h_CSVM_BJetPt_pErr_0bTag = new TH1F("h_CSVM_BJetPt_pErr_0bTag", "Pt Distribution of BJets passing CSVM (pErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_pErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_pErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_pErr_0bTag->Sumw2();

  h_CSVM_BJetEta_pErr_0bTag = new TH1F("h_CSVM_BJetEta_pErr_0bTag", "Eta Distribution of BJets passing CSVM (pErr_0bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_pErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_pErr_0bTag->Sumw2();

  h_CSVM_BJetPhi_pErr_0bTag = new TH1F("h_CSVM_BJetPhi_pErr_0bTag", "Phi Distribution of BJets passing CSVM (pErr_0bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_pErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_pErr_0bTag->Sumw2();

  h_CSVM_BJetPt_pErr_1bTag = new TH1F("h_CSVM_BJetPt_pErr_1bTag", "Pt Distribution of BJets passing CSVM (pErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVM_BJetPt_pErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVM_BJetPt_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPt_pErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVM_BJetPt_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPt_pErr_1bTag->Sumw2();

  h_CSVM_BJetEta_pErr_1bTag = new TH1F("h_CSVM_BJetEta_pErr_1bTag", "Eta Distribution of BJets passing CSVM (pErr_1bTag)", 200, -5.0, 5.0);
  h_CSVM_BJetEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetEta_pErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVM_BJetEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetEta_pErr_1bTag->Sumw2();

  h_CSVM_BJetPhi_pErr_1bTag = new TH1F("h_CSVM_BJetPhi_pErr_1bTag", "Phi Distribution of BJets passing CSVM (pErr_1bTag)", 200, 0.0, 10.0);
  h_CSVM_BJetPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVM_BJetPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJetPhi_pErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVM_BJetPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJetPhi_pErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_noErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_noErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_noErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_0bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/Variable");   h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_0bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_0bTag->Sumw2();

  h_CSVM_GBJetdEta_noErr_0bTag = new TH1F("h_CSVM_GBJetdEta_noErr_0bTag", "DeltaEta Dist Between Gamma & CSVM-BJet (noErr_0bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_noErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_noErr_0bTag->Sumw2();

  h_CSVM_GBJetdPhi_noErr_0bTag = new TH1F("h_CSVM_GBJetdPhi_noErr_0bTag", "DeltaPhi Dist Between Gamma & CSVM-BJet (noErr_0bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_noErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_noErr_0bTag->Sumw2();

  h_CSVM_GBJetdR_noErr_0bTag = new TH1F("h_CSVM_GBJetdR_noErr_0bTag", "DeltaR Dist Between Gamma & CSVM-BJet (noErr_0bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_noErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_noErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_noErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_noErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_noErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_noErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_1bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");  h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_noErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (noErr_1bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_noErr_1bTag->Sumw2();

  h_CSVM_GBJetdEta_noErr_1bTag = new TH1F("h_CSVM_GBJetdEta_noErr_1bTag", "DeltaEta Dist Between Gamma & CSVM-BJet (noErr_1bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_noErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_noErr_1bTag->Sumw2();

  h_CSVM_GBJetdPhi_noErr_1bTag = new TH1F("h_CSVM_GBJetdPhi_noErr_1bTag", "DeltaPhi Dist Between Gamma & CSVM-BJet (noErr_1bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_noErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_noErr_1bTag->Sumw2();

  h_CSVM_GBJetdR_noErr_1bTag = new TH1F("h_CSVM_GBJetdR_noErr_1bTag", "DeltaR Dist Between Gamma & CSVM-BJet (noErr_1bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_noErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_noErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_noErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_pErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_pErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_pErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_0bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_0bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_0bTag->Sumw2();

  h_CSVM_GBJetdEta_pErr_0bTag = new TH1F("h_CSVM_GBJetdEta_pErr_0bTag", "DeltaEta Dist Between Gamma & CSVM-BJet (pErr_0bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_pErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_pErr_0bTag->Sumw2();

  h_CSVM_GBJetdPhi_pErr_0bTag = new TH1F("h_CSVM_GBJetdPhi_pErr_0bTag", "DeltaPhi Dist Between Gamma & CSVM-BJet (pErr_0bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_pErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_pErr_0bTag->Sumw2();

  h_CSVM_GBJetdR_pErr_0bTag = new TH1F("h_CSVM_GBJetdR_pErr_0bTag", "DeltaR Dist Between Gamma & CSVM-BJet (pErr_0bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_pErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_pErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_pErr_0bTag->Sumw2();

  h_CSVM_GBJetInvtMass_binning40_pErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_binning40_pErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_binning40_pErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_1bTag)", nbins, xbins);
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_VariableBinning_pErr_1bTag->Sumw2();

  h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag = new TH1F("h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVM-BJet (pErr_1bTag)", 14000, 0, 14000);
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetInvtMass_UnitBinning_pErr_1bTag->Sumw2();

  h_CSVM_GBJetdEta_pErr_1bTag = new TH1F("h_CSVM_GBJetdEta_pErr_1bTag", "DeltaEta Dist Between Gamma & CSVM-BJet (pErr_1bTag)", 120, 0.0, 6.0);
  h_CSVM_GBJetdEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdEta_pErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVM_GBJetdEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdEta_pErr_1bTag->Sumw2();

  h_CSVM_GBJetdPhi_pErr_1bTag = new TH1F("h_CSVM_GBJetdPhi_pErr_1bTag", "DeltaPhi Dist Between Gamma & CSVM-BJet (pErr_1bTag)", 64, 0.0, 3.2);
  h_CSVM_GBJetdPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVM_GBJetdPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_pErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVM_GBJetdPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdPhi_pErr_1bTag->Sumw2();

  h_CSVM_GBJetdR_pErr_1bTag = new TH1F("h_CSVM_GBJetdR_pErr_1bTag", "DeltaR Dist Between Gamma & CSVM-BJet (pErr_1bTag)", 100, 0.0, 10.0);
  h_CSVM_GBJetdR_pErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVM_GBJetdR_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_GBJetdR_pErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVM_GBJetdR_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_GBJetdR_pErr_1bTag->Sumw2();

  h_CSVT_BJetPt_noErr_0bTag = new TH1F("h_CSVT_BJetPt_noErr_0bTag", "Pt Distribution of BJets passing CSVT (noErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_noErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_noErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_noErr_0bTag->Sumw2();

  h_CSVT_BJetEta_noErr_0bTag = new TH1F("h_CSVT_BJetEta_noErr_0bTag", "Eta Distribution of BJets passing CSVT (noErr_0bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_noErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_noErr_0bTag->Sumw2();

  h_CSVT_BJetPhi_noErr_0bTag = new TH1F("h_CSVT_BJetPhi_noErr_0bTag", "Phi Distribution of BJets passing CSVT (noErr_0bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_noErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_noErr_0bTag->Sumw2();

  h_CSVT_BJetPt_noErr_1bTag = new TH1F("h_CSVT_BJetPt_noErr_1bTag", "Pt Distribution of BJets passing CSVT (noErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_noErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_noErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_noErr_1bTag->Sumw2();

  h_CSVT_BJetEta_noErr_1bTag = new TH1F("h_CSVT_BJetEta_noErr_1bTag", "Eta Distribution of BJets passing CSVT (noErr_1bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_noErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_noErr_1bTag->Sumw2();

  h_CSVT_BJetPhi_noErr_1bTag = new TH1F("h_CSVT_BJetPhi_noErr_1bTag", "Phi Distribution of BJets passing CSVT (noErr_1bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_noErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_noErr_1bTag->Sumw2();

  h_CSVT_BJetPt_pErr_0bTag = new TH1F("h_CSVT_BJetPt_pErr_0bTag", "Pt Distribution of BJets passing CSVT (pErr_0bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_pErr_0bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_pErr_0bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_pErr_0bTag->Sumw2();

  h_CSVT_BJetEta_pErr_0bTag = new TH1F("h_CSVT_BJetEta_pErr_0bTag", "Eta Distribution of BJets passing CSVT (pErr_0bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_pErr_0bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_pErr_0bTag->Sumw2();

  h_CSVT_BJetPhi_pErr_0bTag = new TH1F("h_CSVT_BJetPhi_pErr_0bTag", "Phi Distribution of BJets passing CSVT (pErr_0bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_pErr_0bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_pErr_0bTag->Sumw2();

  h_CSVT_BJetPt_pErr_1bTag = new TH1F("h_CSVT_BJetPt_pErr_1bTag", "Pt Distribution of BJets passing CSVT (pErr_1bTag)", 100, 20.0, 2520.0);
  h_CSVT_BJetPt_pErr_1bTag->GetYaxis()->SetTitle("Events/25 GeV");            h_CSVT_BJetPt_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPt_pErr_1bTag->GetXaxis()->SetTitle("P_{T}^{BJet}");             h_CSVT_BJetPt_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPt_pErr_1bTag->Sumw2();

  h_CSVT_BJetEta_pErr_1bTag = new TH1F("h_CSVT_BJetEta_pErr_1bTag", "Eta Distribution of BJets passing CSVT (pErr_1bTag)", 200, -5.0, 5.0);
  h_CSVT_BJetEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetEta_pErr_1bTag->GetXaxis()->SetTitle("{#eta}^{BJet}");          h_CSVT_BJetEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetEta_pErr_1bTag->Sumw2();

  h_CSVT_BJetPhi_pErr_1bTag = new TH1F("h_CSVT_BJetPhi_pErr_1bTag", "Phi Distribution of BJets passing CSVT (pErr_1bTag)", 200, 0.0, 10.0);
  h_CSVT_BJetPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                 h_CSVT_BJetPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJetPhi_pErr_1bTag->GetXaxis()->SetTitle("{#phi}^{BJet}");          h_CSVT_BJetPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJetPhi_pErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_noErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_noErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_noErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_0bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_0bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_0bTag->Sumw2();

  h_CSVT_GBJetdEta_noErr_0bTag = new TH1F("h_CSVT_GBJetdEta_noErr_0bTag", "DeltaEta Dist Between Gamma & CSVT-BJet (noErr_0bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_noErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_noErr_0bTag->Sumw2();

  h_CSVT_GBJetdPhi_noErr_0bTag = new TH1F("h_CSVT_GBJetdPhi_noErr_0bTag", "DeltaPhi Dist Between Gamma & CSVT-BJet (noErr_0bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_noErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_noErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_noErr_0bTag->Sumw2();

  h_CSVT_GBJetdR_noErr_0bTag = new TH1F("h_CSVT_GBJetdR_noErr_0bTag", "DeltaR Dist Between Gamma & CSVT-BJet (noErr_0bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_noErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_noErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_noErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_noErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_noErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_noErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_1bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_noErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (noErr_1bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_noErr_1bTag->Sumw2();

  h_CSVT_GBJetdEta_noErr_1bTag = new TH1F("h_CSVT_GBJetdEta_noErr_1bTag", "DeltaEta Dist Between Gamma & CSVT-BJet (noErr_1bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_noErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_noErr_1bTag->Sumw2();

  h_CSVT_GBJetdPhi_noErr_1bTag = new TH1F("h_CSVT_GBJetdPhi_noErr_1bTag", "DeltaPhi Dist Between Gamma & CSVT-BJet (noErr_1bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_noErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_noErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_noErr_1bTag->Sumw2();

  h_CSVT_GBJetdR_noErr_1bTag = new TH1F("h_CSVT_GBJetdR_noErr_1bTag", "DeltaR Dist Between Gamma & CSVT-BJet (noErr_1bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_noErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_noErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_noErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_pErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_pErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_0bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_pErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_0bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_0bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_0bTag->Sumw2();

  h_CSVT_GBJetdEta_pErr_0bTag = new TH1F("h_CSVT_GBJetdEta_pErr_0bTag", "DeltaEta Dist Between Gamma & CSVT-BJet (pErr_0bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_pErr_0bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_pErr_0bTag->Sumw2();

  h_CSVT_GBJetdPhi_pErr_0bTag = new TH1F("h_CSVT_GBJetdPhi_pErr_0bTag", "DeltaPhi Dist Between Gamma & CSVT-BJet (pErr_0bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_pErr_0bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_pErr_0bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_pErr_0bTag->Sumw2();

  h_CSVT_GBJetdR_pErr_0bTag = new TH1F("h_CSVT_GBJetdR_pErr_0bTag", "DeltaR Dist Between Gamma & CSVT-BJet (pErr_0bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_pErr_0bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_pErr_0bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_pErr_0bTag->Sumw2();

  h_CSVT_GBJetInvtMass_binning40_pErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_binning40_pErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_1bTag)", 125, 0.0, 5000.0);
  h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->SetTitle("Events/40 GeV");     h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_binning40_pErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_1bTag)", nbins, xbins);
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/Variable");     h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)"); h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_VariableBinning_pErr_1bTag->Sumw2();

  h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag = new TH1F("h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag", "Int Mass Dist of Gamma+CSVT-BJet (pErr_1bTag)", 14000, 0, 14000);
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->SetTitle("Events/UnitBin");    h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->SetTitle("M({#gamma} + BJet)");   h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetInvtMass_UnitBinning_pErr_1bTag->Sumw2();

  h_CSVT_GBJetdEta_pErr_1bTag = new TH1F("h_CSVT_GBJetdEta_pErr_1bTag", "DeltaEta Dist Between Gamma & CSVT-BJet (pErr_1bTag)", 120, 0.0, 6.0);
  h_CSVT_GBJetdEta_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdEta_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdEta_pErr_1bTag->GetXaxis()->SetTitle("#Delta #eta");           h_CSVT_GBJetdEta_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdEta_pErr_1bTag->Sumw2();

  h_CSVT_GBJetdPhi_pErr_1bTag = new TH1F("h_CSVT_GBJetdPhi_pErr_1bTag", "DeltaPhi Dist Between Gamma & CSVT-BJet (pErr_1bTag)", 64, 0.0, 3.2);
  h_CSVT_GBJetdPhi_pErr_1bTag->GetYaxis()->SetTitle("Events");                h_CSVT_GBJetdPhi_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_pErr_1bTag->GetXaxis()->SetTitle("#Delta #phi");           h_CSVT_GBJetdPhi_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdPhi_pErr_1bTag->Sumw2();

  h_CSVT_GBJetdR_pErr_1bTag = new TH1F("h_CSVT_GBJetdR_pErr_1bTag", "DeltaR Dist Between Gamma & CSVT-BJet (pErr_1bTag)", 100, 0.0, 10.0);
  h_CSVT_GBJetdR_pErr_1bTag->GetYaxis()->SetTitle("Events");                  h_CSVT_GBJetdR_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_GBJetdR_pErr_1bTag->GetXaxis()->SetTitle("#Delta R");                h_CSVT_GBJetdR_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_GBJetdR_pErr_1bTag->Sumw2();

  h_BJetDiscByCSV = new TH1F("h_BJetDiscByCSV", "Value of CSV BJet Discriminator for leading jet in each event", 100, 0.0, 5.0);
  h_BJetDiscByCSV->GetYaxis()->SetTitle("Events");               h_BJetDiscByCSV->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV->GetXaxis()->SetTitle("BJetDisc_CSV");          h_BJetDiscByCSV->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV->Sumw2();

  h_BJetDiscByCSV_PassingCSVL_noErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVL_noErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVL (noErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVL"); h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_noErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVL_noErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVL_noErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVL (noErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVL"); h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_noErr_1bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVL_pErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVL_pErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVL (pErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVL"); h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_pErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVL_pErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVL_pErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVL (pErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVL"); h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVL_pErr_1bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVM_noErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVM_noErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVM (noErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVM"); h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_noErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVM_noErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVM_noErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVM (noErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVM"); h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_noErr_1bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVM_pErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVM_pErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVM (pErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVM"); h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_pErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVM_pErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVM_pErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVM (pErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVM"); h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVM_pErr_1bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVT_noErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVT_noErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVT (noErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVT"); h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_noErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVT_noErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVT_noErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVT (noErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVT"); h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_noErr_1bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVT_pErr_0bTag = new TH1F("h_BJetDiscByCSV_PassingCSVT_pErr_0bTag", "Value of CSV BJet Disc for leading jet passing CSVT (pErr_0bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->GetXaxis()->SetTitle("BJetDisc_CSVT"); h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_pErr_0bTag->Sumw2();

  h_BJetDiscByCSV_PassingCSVT_pErr_1bTag = new TH1F("h_BJetDiscByCSV_PassingCSVT_pErr_1bTag", "Value of CSV BJet Disc for leading jet passing CSVT (pErr_1bTag)", 100, 0.0, 5.0);
  h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->GetYaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->GetXaxis()->SetTitle("BJetDisc_CSVT"); h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->GetXaxis()->CenterTitle();
  h_BJetDiscByCSV_PassingCSVT_pErr_1bTag->Sumw2();

  h_trueNumofInteractions = new TH1F("h_trueNumofInteractions", "True no. of interactions distribution in MC without any weight", 60, 0, 60);
  h_trueNumofInteractions->GetYaxis()->SetTitle("Events");      h_trueNumofInteractions->GetYaxis()->CenterTitle();
  h_trueNumofInteractions->GetXaxis()->SetTitle("nPV_noWt (MC)");  h_trueNumofInteractions->GetXaxis()->CenterTitle();
  h_trueNumofInteractions->Sumw2();

  h_goodPV = new TH1F("h_goodPV", "No. of Good Primary Vertices", 60, 0, 60);
  h_goodPV->GetYaxis()->SetTitle("Events");        h_goodPV->GetYaxis()->CenterTitle();
  h_goodPV->GetXaxis()->SetTitle("nPV (MC)");    h_goodPV->GetXaxis()->CenterTitle();
  h_goodPV->Sumw2();

  h_goodPV_LumiWt = new TH1F("h_goodPV_LumiWt", "No. of Good Primary Vertices with lumi wt", 60, 0, 60);
  h_goodPV_LumiWt->GetYaxis()->SetTitle("Events");        h_goodPV_LumiWt->GetYaxis()->CenterTitle();
  h_goodPV_LumiWt->GetXaxis()->SetTitle("nPV (MC)");    h_goodPV_LumiWt->GetXaxis()->CenterTitle();
  h_goodPV_LumiWt->Sumw2();

  h_goodPV_PUWt = new TH1F("h_goodPV_PUWt", "No. of Good Primary Vertices with pileup wt", 60, 0, 60);
  h_goodPV_PUWt->GetYaxis()->SetTitle("Events");        h_goodPV_PUWt->GetYaxis()->CenterTitle();
  h_goodPV_PUWt->GetXaxis()->SetTitle("nPV (MC)");    h_goodPV_PUWt->GetXaxis()->CenterTitle();
  h_goodPV_PUWt->Sumw2();

  h_goodPV_LumiNPUWt = new TH1F("h_goodPV_LumiNPUWt", "No. of Good Primary Vertices with lumi and pileup wt", 60, 0, 60);
  h_goodPV_LumiNPUWt->GetYaxis()->SetTitle("Events");        h_goodPV_LumiNPUWt->GetYaxis()->CenterTitle();
  h_goodPV_LumiNPUWt->GetXaxis()->SetTitle("nPV (MC)");    h_goodPV_LumiNPUWt->GetXaxis()->CenterTitle();
  h_goodPV_LumiNPUWt->Sumw2();

  h_nPhotons = new TH1F("h_nPhotons", "Number of Photons in each event", 200, 0, 100);
  h_nPhotons->GetYaxis()->SetTitle("Events");     h_nPhotons->GetYaxis()->CenterTitle();
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

  h_DataPUDist = new TH1F("h_DataPUDist", "Data PileUp Distribution", 60, 0, 60);
  h_DataPUDist->GetYaxis()->SetTitle("Events");               h_DataPUDist->GetYaxis()->CenterTitle();
  h_DataPUDist->GetXaxis()->SetTitle("nPUV");                 h_DataPUDist->GetXaxis()->CenterTitle();
  h_DataPUDist->Sumw2();

  h_DataPUNormDist = new TH1F("h_DataPUNormDist", "Normalized Data PileUp Distribution", 60, 0, 60);
  h_DataPUNormDist->GetYaxis()->SetTitle("Events");           h_DataPUNormDist->GetYaxis()->CenterTitle();
  h_DataPUNormDist->GetXaxis()->SetTitle("nPUV");             h_DataPUNormDist->GetXaxis()->CenterTitle();
  h_DataPUNormDist->Sumw2();
  
  h_MCPUNormDist = new TH1F("h_MCPUNormDist", "Normalized MC PileUp Distribution",60, 0, 60);
  h_MCPUNormDist->GetYaxis()->SetTitle("Events");             h_MCPUNormDist->GetYaxis()->CenterTitle();
  h_MCPUNormDist->GetXaxis()->SetTitle("nPUV");               h_MCPUNormDist->GetXaxis()->CenterTitle();
  h_MCPUNormDist->Sumw2();
  
  h_PUScaleFactor = new TH1F("h_PUScaleFactor", "PileUp Scale Factors Distribution", 60, 0, 60);
  h_PUScaleFactor->GetYaxis()->SetTitle("Scale Factor");      h_PUScaleFactor->GetYaxis()->CenterTitle();
  h_PUScaleFactor->GetXaxis()->SetTitle("nPUV");              h_PUScaleFactor->GetXaxis()->CenterTitle();
  h_PUScaleFactor->Sumw2();

  h_PUScaleFactor_XSm5 = new TH1F("h_PUScaleFactor_XSm5", "PileUp Scale Factors Distribution for 5% less XS", 60, 0, 60);
  h_PUScaleFactor_XSm5->GetYaxis()->SetTitle("Scale Factor");      h_PUScaleFactor_XSm5->GetYaxis()->CenterTitle();
  h_PUScaleFactor_XSm5->GetXaxis()->SetTitle("nPUV");              h_PUScaleFactor_XSm5->GetXaxis()->CenterTitle();
  h_PUScaleFactor_XSm5->Sumw2();

  h_PUScaleFactor_XSp5 = new TH1F("h_PUScaleFactor_XSp5", "PileUp Scale Factors Distribution for 5% more XS", 60, 0, 60);
  h_PUScaleFactor_XSp5->GetYaxis()->SetTitle("Scale Factor");      h_PUScaleFactor_XSp5->GetYaxis()->CenterTitle();
  h_PUScaleFactor_XSp5->GetXaxis()->SetTitle("nPUV");              h_PUScaleFactor_XSp5->GetXaxis()->CenterTitle();
  h_PUScaleFactor_XSp5->Sumw2();

  h_Jet_PartonFlavor = new TH1F("h_Jet_PartonFlavor", "True Parton flavour of Jets", 24, 0, 24);
  h_Jet_PartonFlavor->GetYaxis()->SetTitle("Events");        h_Jet_PartonFlavor->GetYaxis()->CenterTitle();
  h_Jet_PartonFlavor->GetXaxis()->SetTitle("PdgId");         h_Jet_PartonFlavor->GetXaxis()->CenterTitle();
  h_Jet_PartonFlavor->Sumw2();

  h_CSVL_BJet_PartonFlavor_noErr_0bTag = new TH1F("h_CSVL_BJet_PartonFlavor_noErr_0bTag", "True Parton flavour of BJets passing CSVL noErr_0bTag", 24, 0, 24);
  h_CSVL_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVL_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVL_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_noErr_0bTag->Sumw2();

  h_CSVL_BJet_PartonFlavor_noErr_1bTag = new TH1F("h_CSVL_BJet_PartonFlavor_noErr_1bTag", "True Parton flavour of BJets passing CSVL noErr_1bTag", 24, 0, 24);
  h_CSVL_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVL_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVL_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_noErr_1bTag->Sumw2();

  h_CSVL_BJet_PartonFlavor_pErr_0bTag = new TH1F("h_CSVL_BJet_PartonFlavor_pErr_0bTag", "True Parton flavour of BJets passing CSVL pErr_0bTag", 24, 0, 24);
  h_CSVL_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVL_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVL_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_pErr_0bTag->Sumw2();

  h_CSVL_BJet_PartonFlavor_pErr_1bTag = new TH1F("h_CSVL_BJet_PartonFlavor_pErr_1bTag", "True Parton flavour of BJets passing CSVL pErr_1bTag", 24, 0, 24);
  h_CSVL_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVL_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVL_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVL_BJet_PartonFlavor_pErr_1bTag->Sumw2();

  h_CSVM_BJet_PartonFlavor_noErr_0bTag = new TH1F("h_CSVM_BJet_PartonFlavor_noErr_0bTag", "True Parton flavour of BJets passing CSVM noErr_0bTag", 24, 0, 24);
  h_CSVM_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVM_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVM_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_noErr_0bTag->Sumw2();

  h_CSVM_BJet_PartonFlavor_noErr_1bTag = new TH1F("h_CSVM_BJet_PartonFlavor_noErr_1bTag", "True Parton flavour of BJets passing CSVM noErr_1bTag", 24, 0, 24);
  h_CSVM_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVM_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVM_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_noErr_1bTag->Sumw2();

  h_CSVM_BJet_PartonFlavor_pErr_0bTag = new TH1F("h_CSVM_BJet_PartonFlavor_pErr_0bTag", "True Parton flavour of BJets passing CSVM pErr_0bTag", 24, 0, 24);
  h_CSVM_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVM_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVM_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_pErr_0bTag->Sumw2();

  h_CSVM_BJet_PartonFlavor_pErr_1bTag = new TH1F("h_CSVM_BJet_PartonFlavor_pErr_1bTag", "True Parton flavour of BJets passing CSVM pErr_1bTag", 24, 0, 24);
  h_CSVM_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVM_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVM_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVM_BJet_PartonFlavor_pErr_1bTag->Sumw2();

  h_CSVT_BJet_PartonFlavor_noErr_0bTag = new TH1F("h_CSVT_BJet_PartonFlavor_noErr_0bTag", "True Parton flavour of BJets passing CSVT noErr_0bTag", 24, 0, 24);
  h_CSVT_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVT_BJet_PartonFlavor_noErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVT_BJet_PartonFlavor_noErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_noErr_0bTag->Sumw2();

  h_CSVT_BJet_PartonFlavor_noErr_1bTag = new TH1F("h_CSVT_BJet_PartonFlavor_noErr_1bTag", "True Parton flavour of BJets passing CSVT noErr_1bTag", 24, 0, 24);
  h_CSVT_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVT_BJet_PartonFlavor_noErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVT_BJet_PartonFlavor_noErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_noErr_1bTag->Sumw2();

  h_CSVT_BJet_PartonFlavor_pErr_0bTag = new TH1F("h_CSVT_BJet_PartonFlavor_pErr_0bTag", "True Parton flavour of BJets passing CSVT pErr_0bTag", 24, 0, 24);
  h_CSVT_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->SetTitle("Events");        h_CSVT_BJet_PartonFlavor_pErr_0bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->SetTitle("PdgId");         h_CSVT_BJet_PartonFlavor_pErr_0bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_pErr_0bTag->Sumw2();

  h_CSVT_BJet_PartonFlavor_pErr_1bTag = new TH1F("h_CSVT_BJet_PartonFlavor_pErr_1bTag", "True Parton flavour of BJets passing CSVT pErr_1bTag", 24, 0, 24);
  h_CSVT_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->SetTitle("Events");        h_CSVT_BJet_PartonFlavor_pErr_1bTag->GetYaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->SetTitle("PdgId");         h_CSVT_BJet_PartonFlavor_pErr_1bTag->GetXaxis()->CenterTitle();
  h_CSVT_BJet_PartonFlavor_pErr_1bTag->Sumw2();

  h_TrueBJetPt = new TH1F("h_TrueBJetPt", "Pt distribution of True B Jets", 100, 20.0, 2520.0);
  h_TrueBJetPt->GetYaxis()->SetTitle("Events/25 GeV");        h_TrueBJetPt->GetYaxis()->CenterTitle();
  h_TrueBJetPt->GetXaxis()->SetTitle("P_{T}^{BJet}");         h_TrueBJetPt->GetXaxis()->CenterTitle();
  h_TrueBJetPt->Sumw2();

  h_TrueBJetEta = new TH1F("h_TrueBJetEta", "Eta Distribution of True B Jets", 200, -5.0, 5.0);
  h_TrueBJetEta->GetYaxis()->SetTitle("Events");              h_TrueBJetEta->GetYaxis()->CenterTitle();
  h_TrueBJetEta->GetXaxis()->SetTitle("{#eta}^{BJet}");       h_TrueBJetEta->GetXaxis()->CenterTitle();
  h_TrueBJetEta->Sumw2();

  h_TrueBJetPtPassingCSVL = new TH1F("h_TrueBJetPtPassingCSVL", "Pt distribution of True B Jets passing CSVL", 100, 20.0, 2520.0);
  h_TrueBJetPtPassingCSVL->GetYaxis()->SetTitle("Events/25 GeV"); h_TrueBJetPtPassingCSVL->GetYaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVL->GetXaxis()->SetTitle("P_{T}^{BJet}");  h_TrueBJetPtPassingCSVL->GetXaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVL->Sumw2();

  h_TrueBJetEtaPassingCSVL = new TH1F("h_TrueBJetEtaPassingCSVL", "Eta Distribution of True B Jets passing CSVL", 200, -5.0, 5.0);
  h_TrueBJetEtaPassingCSVL->GetYaxis()->SetTitle("Events");        h_TrueBJetEtaPassingCSVL->GetYaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVL->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_TrueBJetEtaPassingCSVL->GetXaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVL->Sumw2();

  h_TrueBJetPtPassingCSVM = new TH1F("h_TrueBJetPtPassingCSVM", "Pt distribution of True B Jets passing CSVM", 100, 20.0, 2520.0);
  h_TrueBJetPtPassingCSVM->GetYaxis()->SetTitle("Events/25 GeV"); h_TrueBJetPtPassingCSVM->GetYaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVM->GetXaxis()->SetTitle("P_{T}^{BJet}");  h_TrueBJetPtPassingCSVM->GetXaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVM->Sumw2();

  h_TrueBJetEtaPassingCSVM = new TH1F("h_TrueBJetEtaPassingCSVM", "Eta Distribution of True B Jets passing CSVM", 200, -5.0, 5.0);
  h_TrueBJetEtaPassingCSVM->GetYaxis()->SetTitle("Events");        h_TrueBJetEtaPassingCSVM->GetYaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVM->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_TrueBJetEtaPassingCSVM->GetXaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVM->Sumw2();

  h_TrueBJetPtPassingCSVT = new TH1F("h_TrueBJetPtPassingCSVT", "Pt distribution of True B Jets passing CSVT", 100, 20.0, 2520.0);
  h_TrueBJetPtPassingCSVT->GetYaxis()->SetTitle("Events/25 GeV");  h_TrueBJetPtPassingCSVT->GetYaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVT->GetXaxis()->SetTitle("P_{T}^{BJet}");    h_TrueBJetPtPassingCSVT->GetXaxis()->CenterTitle();
  h_TrueBJetPtPassingCSVT->Sumw2();

  h_TrueBJetEtaPassingCSVT = new TH1F("h_TrueBJetEtaPassingCSVT", "Eta Distribution of True B Jets passing CSVT", 200, -5.0, 5.0);
  h_TrueBJetEtaPassingCSVT->GetYaxis()->SetTitle("Events");        h_TrueBJetEtaPassingCSVT->GetYaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVT->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_TrueBJetEtaPassingCSVT->GetXaxis()->CenterTitle();
  h_TrueBJetEtaPassingCSVT->Sumw2();

  h_NonBJetPt = new TH1F("h_NonBJetPt", "Pt distribution of Non B Jets", 100, 20.0, 2520.0);
  h_NonBJetPt->GetYaxis()->SetTitle("Events/25 GeV");        h_NonBJetPt->GetYaxis()->CenterTitle();
  h_NonBJetPt->GetXaxis()->SetTitle("P_{T}^{BJet}");         h_NonBJetPt->GetXaxis()->CenterTitle();
  h_NonBJetPt->Sumw2();

  h_NonBJetEta = new TH1F("h_NonBJetEta", "Eta Distribution of Non B Jets", 200, -5.0, 5.0);
  h_NonBJetEta->GetYaxis()->SetTitle("Events");              h_NonBJetEta->GetYaxis()->CenterTitle();
  h_NonBJetEta->GetXaxis()->SetTitle("{#eta}^{BJet}");       h_NonBJetEta->GetXaxis()->CenterTitle();
  h_NonBJetEta->Sumw2();

  h_NonBJetPtPassingCSVL = new TH1F("h_NonBJetPtPassingCSVL", "Pt distribution of Non B Jets passing CSVL", 100, 20.0, 2520.0);
  h_NonBJetPtPassingCSVL->GetYaxis()->SetTitle("Events/25 GeV"); h_NonBJetPtPassingCSVL->GetYaxis()->CenterTitle();
  h_NonBJetPtPassingCSVL->GetXaxis()->SetTitle("P_{T}^{BJet}");  h_NonBJetPtPassingCSVL->GetXaxis()->CenterTitle();
  h_NonBJetPtPassingCSVL->Sumw2();

  h_NonBJetEtaPassingCSVL = new TH1F("h_NonBJetEtaPassingCSVL", "Eta Distribution of Non B Jets passing CSVL", 200, -5.0, 5.0);
  h_NonBJetEtaPassingCSVL->GetYaxis()->SetTitle("Events");        h_NonBJetEtaPassingCSVL->GetYaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVL->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_NonBJetEtaPassingCSVL->GetXaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVL->Sumw2();

  h_NonBJetPtPassingCSVM = new TH1F("h_NonBJetPtPassingCSVM", "Pt distribution of Non B Jets passing CSVM", 100, 20.0, 2520.0);
  h_NonBJetPtPassingCSVM->GetYaxis()->SetTitle("Events/25 GeV"); h_NonBJetPtPassingCSVM->GetYaxis()->CenterTitle();
  h_NonBJetPtPassingCSVM->GetXaxis()->SetTitle("P_{T}^{BJet}");  h_NonBJetPtPassingCSVM->GetXaxis()->CenterTitle();
  h_NonBJetPtPassingCSVM->Sumw2();

  h_NonBJetEtaPassingCSVM = new TH1F("h_NonBJetEtaPassingCSVM", "Eta Distribution of Non B Jets passing CSVM", 200, -5.0, 5.0);
  h_NonBJetEtaPassingCSVM->GetYaxis()->SetTitle("Events");        h_NonBJetEtaPassingCSVM->GetYaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVM->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_NonBJetEtaPassingCSVM->GetXaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVM->Sumw2();

  h_NonBJetPtPassingCSVT = new TH1F("h_NonBJetPtPassingCSVT", "Pt distribution of Non B Jets passing CSVT", 100, 20.0, 2520.0);
  h_NonBJetPtPassingCSVT->GetYaxis()->SetTitle("Events/25 GeV");  h_NonBJetPtPassingCSVT->GetYaxis()->CenterTitle();
  h_NonBJetPtPassingCSVT->GetXaxis()->SetTitle("P_{T}^{BJet}");    h_NonBJetPtPassingCSVT->GetXaxis()->CenterTitle();
  h_NonBJetPtPassingCSVT->Sumw2();

  h_NonBJetEtaPassingCSVT = new TH1F("h_NonBJetEtaPassingCSVT", "Eta Distribution of Non B Jets passing CSVT", 200, -5.0, 5.0);
  h_NonBJetEtaPassingCSVT->GetYaxis()->SetTitle("Events");        h_NonBJetEtaPassingCSVT->GetYaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVT->GetXaxis()->SetTitle("{#eta}^{BJet}"); h_NonBJetEtaPassingCSVT->GetXaxis()->CenterTitle();
  h_NonBJetEtaPassingCSVT->Sumw2();

  char hname[100];
  char htit[100];
  std::string scale[2] = {"up", "down"};
  std::string tag[2] = {"_NobTagging", "_WithbTagging"};
  for (int ii = 0; ii < 2; ii++){
    for(int jj = 0; jj < 2; jj++){

      sprintf(hname, "h_profile_JES%s%s", scale[ii].c_str(), tag[jj].c_str());
      sprintf(htit, "Profile of JES%s vs Mass%s", scale[ii].c_str(), tag[jj].c_str());

      h_profile_JES[ii][jj] = new TProfile(hname, htit, nbins, xbins);


      sprintf(hname, "h_mass_JES%s%s", scale[ii].c_str(), tag[jj].c_str());
      sprintf(htit, "Mass Distribution for JES%s%s", scale[ii].c_str(), tag[jj].c_str());

      h_mass_JES[ii][jj] = new TH1F(hname, htit, nbins, xbins);
      h_mass_JES[ii][jj]->GetYaxis()->SetTitle("Events"); h_mass_JES[ii][jj]->GetYaxis()->CenterTitle();
      h_mass_JES[ii][jj]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)"); h_mass_JES[ii][jj]->GetXaxis()->CenterTitle();
      h_mass_JES[ii][jj]->Sumw2();


      sprintf(hname, "h_profile_PES%s%s", scale[ii].c_str(), tag[jj].c_str());
      sprintf(htit, "Profile of PES%s vs Mass%s", scale[ii].c_str(), tag[jj].c_str());

      h_profile_PES[ii][jj] = new TProfile(hname, htit, nbins, xbins);


      sprintf(hname, "h_mass_PES%s%s", scale[ii].c_str(), tag[jj].c_str());
      sprintf(htit, "Mass Distribution for PES%s%s", scale[ii].c_str(), tag[jj].c_str());

      h_mass_PES[ii][jj] = new TH1F(hname, htit, nbins, xbins);
      h_mass_PES[ii][jj]->GetYaxis()->SetTitle("Events"); h_mass_PES[ii][jj]->GetYaxis()->CenterTitle();
      h_mass_PES[ii][jj]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)"); h_mass_PES[ii][jj]->GetXaxis()->CenterTitle();
      h_mass_PES[ii][jj]->Sumw2();

    }
  }







}




#endif // #ifdef PostAnalyzer_MC_cxx




EOF

cat > analysis_${filenameTag}.C <<EOF
#include "PostAnalyzer_MC.C"
#include "TROOT.h"
int main(){
    PostAnalyzer_MC a;
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


((sampleIndex++))

done ##end of for loop##
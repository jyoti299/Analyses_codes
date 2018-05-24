#!/bin/bash

sampleIndex=0

eosDir=/eos/uscms/store/user/rocky86
lpcqstarDir=/eos/uscms/store/user/lpcqstar
leptonjetsDir=/eos/uscms/store/user/leptonjets
MCnTuplesDir=root://cmseos.fnal.gov//eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC
#MCnTuplesDir=root://cmsxrootd.fnal.gov//store/user/lpcqstar/13TeV/NTuples/MC
GJets=GJets
DiJet=DiJet
Signal_Qstar=Signal_Qstar


#myDir=${eosDir}/13TeV
#PostAna_Dir=${myDir}/PostAnalyzer_Results
#PostAna_Dir_final=${PostAna_Dir}/Standard_LID
#PostAnaMC_Dir=${PostAna_Dir_final}/MC

#Uncorrected xs for gamma jet
#XS=(  23080.0  9110.0  2281.0  273.0  94.5  )

#GammaJet (cross-sections corrected using GenXSanalyzer)
GJ=1
for i in GJets_HT40To100  GJets_HT100To200  GJets_HT200To400  GJets_HT400To600  GJets_HT600ToInf
do 
XS=( 20820.0   9201.0   2308.0   275.2   93.31  )
totalEvents=(  4816232  5027618  10424189  2476770   2550765   )

sourceDir=${MCnTuplesDir}/${GJets}/${i}/
#outputDir=${PostAnaMC_Dir}/${GJets}

#QCDDiJet 
QCD=0
#for i in QCD_Pt_120to170  QCD_Pt_170to300  QCD_Pt_300to470  QCD_Pt_470to600  QCD_Pt_600to800  QCD_Pt_800to1000  QCD_Pt_1000to1400  QCD_Pt_1400to1800  QCD_Pt_1800to2400  QCD_Pt_2400to3200  QCD_Pt_3200toInf
#do
#XS=( 471100.0  117276.0  7823.0  648.2  186.9  32.293  9.4183  0.84265  0.114943  0.00682981  0.000165445  )
#totalEvents=( 3458385  3364368  2935633  1937537  1964128  1937216  1487136  197959  193608  194456  192944  )

#sourceDir=${MCnTuplesDir}/${DiJet}/${i}/
#outputDir=${PostAnaMC_Dir}/${DiJet}

#EWK (cross-sections corrected using GenXSanalyzer)
EWK=0
#for i in WGToLNuG  WJetsToLNu  ZGTo2LG  ZJets
#do
#XS=(377.5  50360  123.7  5813)
#totalEvents=(6099599  72220670  4451319  9956130 )

#sourceDir=${MCnTuplesDir}/${i}/
#outputDir=${PostAnaMC_Dir}/${EWK}

#QstarSignal f1p0 (total events checked with CMSDAS, also corrected in twiki   ) 
Qstar=0
#for i in QstarToGJ_M500_f1p0  QstarToGJ_M1000_f1p0  QstarToGJ_M2000_f1p0  QstarToGJ_M3000_f1p0  QstarToGJ_M4000_f1p0  QstarToGJ_M5000_f1p0  QstarToGJ_M6000_f1p0  QstarToGJ_M7000_f1p0  QstarToGJ_M8000_f1p0  QstarToGJ_M9000_f1p0
#do
#XS=( 3.033E2  1.632E1  5.213E-1  4.272E-2  4.8E-3  5.835E-4  7.076E-5  8.66E-6  1.283E-6  2.985E-7  )
#totalEvents=( 298953  299072  100000  74879  74271  74802  74822  49544  49628  50000  )

#########################################################################################################
#QstarSignal f0p5 (total events checked with CMSDAS, also corrected in twiki)
#for i in QstarToGJ_M500_f0p5  QstarToGJ_M1000_f0p5  QstarToGJ_M2000_f0p5  QstarToGJ_M3000_f0p5  QstarToGJ_M4000_f0p5  QstarToGJ_M5000_f0p5  QstarToGJ_M6000_f0p5  QstarToGJ_M7000_f0p5  QstarToGJ_M8000_f0p5  QstarToGJ_M9000_f0p5
#do
#XS=( 7.378E1  4.129E0  1.328E-1  1.095E-2  1.212E-3  1.437E-4  1.62E-5  1.672E-6  1.647E-7  2.329E-8  )
#totalEvents=( 297624  97970  98740  74232  74857  74896  74328  49888  48200  49403  )

#QstarSignal f0p1 (total events checked with CMSDAS, also corrected in twiki   )
#for i in QstarToGJ_M500_f0p1  QstarToGJ_M1000_f0p1  QstarToGJ_M2000_f0p1  QstarToGJ_M3000_f0p1  QstarToGJ_M4000_f0p1  QstarToGJ_M5000_f0p1  QstarToGJ_M6000_f0p1  QstarToGJ_M7000_f0p1  QstarToGJ_M8000_f0p1  QstarToGJ_M9000_f0p1
#do
#XS=( 2.955E0  1.655E-1  5.315E-3  4.356E-4  4.861E-5  5.715E-6  6.241E-7  5.973E-8  4.515E-9  2.655E-10  )
#totalEvents=( 97208  98290  73128  74286  74836  75000  49918  49679  49316  49782  )

#########################################################################################################

#for i in QstarToGJ_M500_f1p0  QstarToGJ_M1000_f1p0  QstarToGJ_M2000_f1p0  QstarToGJ_M3000_f1p0  QstarToGJ_M4000_f1p0  QstarToGJ_M5000_f1p0 QstarToGJ_M6000_f1p0 QstarToGJ_M7000_f1p0
#do
#XS=( 3.033E2  1.632E1  5.213E-1  4.272E-2  4.8E-3  5.835E-4  7.076E-5  8.66E-6 )
#totalEvents=( 298953  299072  100000  74879  74271  74802  74822  49544 )
#mass_norm=( 500.0  1000.0  2000.0  3000.0  4000.0  5000.0  6000.0  7000.0 )

#sourceDir=${MCnTuplesDir}/${Signal_Qstar}/${i}/
#outputDir=${PostAnaMC_Dir}/${Signal_Qstar}

#This is required for calculating systematics and massX dist for signal only
#if [ ${sourceDir} = "${MCnTuplesDir}/${Signal_Qstar}/${i}/" ]
#then
#Signal=true
#else
#Signal=false
#fi

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


#if [ ! -d ${myDir} ]; then
#echo "Making Directory ${myDir}"
#mkdir ${myDir}
#chmod 775 ${myDir}
#fi

#if [ ! -d ${PostAna_Dir} ]; then
#echo "Making Directory ${PostAna_Dir}"
#mkdir ${PostAna_Dir}
#chmod 775 ${PostAna_Dir}
#fi

#if [ ! -d ${PostAna_Dir_final} ]; then
#echo "Making Directory ${PostAna_Dir_final}"
#mkdir ${PostAna_Dir_final}
#chmod 775 ${PostAna_Dir_final}
#fi

#if [ ! -d ${PostAnaMC_Dir} ]; then
#echo "Making Directory ${PostAnaMC_Dir}"
#mkdir ${PostAnaMC_Dir}
#chmod 775 ${PostAnaMC_Dir}
#fi

#if [ ! -d ${destinationDir} ]; then
#echo "Making Directory ${destinationDir}"
#mkdir ${destinationDir}
#chmod 775 ${destinationDir}
#fi

echo "filenameTag = ${filenameTag}"
#echo "destinationDir = ${destinationDir}"
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

   //Uncomment this in script
   //Define Output file here
   TString OutputPath = "${destinationDir}/";
   TString OutputFile = "${filenameTag}";
   // file = new TFile(OutputPath+OutputFile+".root", "RECREATE");
   file = new TFile(OutputFile+".root", "RECREATE");

   //Uncomment it in script
   IfGJet = ${GJ};

   //Define Histograms here
   BookHistograms();

   //Running function for Pile up reweighting
   PileupReWeighting();

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
      Lumi_EvtWt = (Lumi*(${XS[${sampleIndex}]}))/${totalEvents[${sampleIndex}]};

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
      if(PC > -1) PC_Gen = MatchedGenPhotonToReco(PC); 
      if(JC > -1) JC_Gen = MatchedGenJetToReco(JC);

      //Getting total event wt
      if(IfGJet){
	if(PC > -1){
	  Double_t mcfm_kfac = MCFM_kFact((*phoEt)[PC]);
	  //uncomment in script
	  //Double_t mult = (${XS_LO[${sampleIndex}]})/(${XS[${sampleIndex}]});
	  //comment out in script
	  Double_t mult = 1;
	  Total_EvtWt = PreBTag_EvtWt * mcfm_kfac * mult;
	  if(mcfm_kfac == 0 || mcfm_kfac < 0){
	    Total_EvtWt = PreBTag_EvtWt;
	  }
	}
	else{
	  Total_EvtWt = PreBTag_EvtWt;
	}
      }
      else{
	Total_EvtWt = PreBTag_EvtWt;
      }
      //-----------------------

      h_trueNumofInt       ->Fill((*puTrue)[0]);
      h_goodPV             ->Fill(GoodVertex);
      h_goodPV_LumiWt      ->Fill(GoodVertex, Lumi_EvtWt);
      h_goodPV_PUWt        ->Fill(GoodVertex, PU_EvtWt);
      h_goodPV_TotalWt[0]  ->Fill(GoodVertex, Total_EvtWt);
      h_nPhotons[0]        ->Fill(GoodIsoPhotons.size(), Total_EvtWt);
      h_nJets[0]           ->Fill(GoodIsoJets.size(), Total_EvtWt);

      h_CutFlow->Fill(0.5);
      h_CutFlowWithWts->Fill(0.5, Total_EvtWt);

      if(Pass_HLT){
	h_CutFlow->Fill(1.5);
	h_CutFlowWithWts->Fill(1.5, Total_EvtWt);

	if(HasPrimaryVtx){
	  h_CutFlow->Fill(2.5);
	  h_CutFlowWithWts->Fill(2.5, Total_EvtWt);

	  if(GoodIsoPhotons.size() > 0){
	    h_CutFlow->Fill(3.5);
	    h_CutFlowWithWts->Fill(3.5, Total_EvtWt);

	    if(PC > -1){
	      h_CutFlow->Fill(4.5);
	      h_CutFlowWithWts->Fill(4.5, Total_EvtWt);

	      if(JC > -1){
		h_CutFlow->Fill(5.5);
		h_CutFlowWithWts->Fill(5.5, Total_EvtWt);

		if(Pass_JetPt){
		  h_CutFlow->Fill(6.5);
		  h_CutFlowWithWts->Fill(6.5, Total_EvtWt);

		  if(Pass_JetEta){
		    h_CutFlow->Fill(7.5);
		    h_CutFlowWithWts->Fill(7.5, Total_EvtWt);

		    if(Pass_GJdPhi){
		      h_CutFlow->Fill(8.5);
		      h_CutFlowWithWts->Fill(8.5, Total_EvtWt);

		      if(Pass_GJdEta){
			h_CutFlow->Fill(9.5);
			h_CutFlowWithWts->Fill(9.5, Total_EvtWt);

			  h_PC                       ->Fill(PC, Total_EvtWt);
			  h_JC                       ->Fill(JC, Total_EvtWt);

                          h_PhotonPt[0]              ->Fill((*phoEt)[PC], Total_EvtWt);
                          h_PhotonEta[0]             ->Fill((*phoSCEta)[PC], Total_EvtWt);
                          h_PhotonPhi[0]             ->Fill((*phoSCPhi)[PC], Total_EvtWt);
			  h_Photon_SigmaIEtaIEta[0]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt);
			  h_Photon_R9[0]             ->Fill((*phoR9)[PC], Total_EvtWt);
			  h_Photon_HoverE[0]         ->Fill((*phoHoverE)[PC], Total_EvtWt);
			  h_Photon_EleVeto[0]        ->Fill((*phoEleVeto)[PC], Total_EvtWt);
			  h_Photon_CorrPFChIso[0]    ->Fill((*phoPFChIso)[PC], Total_EvtWt);
			  h_Photon_CorrPFNeuIso[0]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt);
			  h_Photon_CorrPFPhoIso[0]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt);

                          h_JetPt[0]                 ->Fill((*jetPt)[JC], Total_EvtWt);
                          h_JetEta[0]                ->Fill((*jetEta)[JC], Total_EvtWt);
                          h_JetPhi[0]                ->Fill((*jetPhi)[JC], Total_EvtWt);
			  h_Jet_NHEF[0]              ->Fill((*jetNHF)[JC], Total_EvtWt);
			  h_Jet_NEEF[0]              ->Fill((*jetNEF)[JC], Total_EvtWt);
			  h_Jet_NConst[0]            ->Fill((*jetNConstituents)[JC], Total_EvtWt);
			  h_Jet_CHEF[0]              ->Fill((*jetCHF)[JC], Total_EvtWt);
			  h_Jet_ChMult[0]            ->Fill((*jetNCH)[JC], Total_EvtWt);
			  h_Jet_CEEF[0]              ->Fill((*jetCEF)[JC], Total_EvtWt);

                          h_GJetInvtMass_bin40[0]    ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
                          h_GJetInvtMass_VarBin[0]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
                          h_GJetInvtMass_UnitBin[0]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
			  h_GJet_dEta[0]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt);
			  h_GJet_dPhi[0]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt);
			  h_GJet_dR[0]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt);

			  ///Uncomment in script for q*
  			    //h_mass_X_bin1              ->Fill(GetInvtMass(PC, JC)/${mass_norm[${sampleIndex}]}, Total_EvtWt);

			  h_PFMet[0]                 ->Fill(pfMET, Total_EvtWt);
			  h_PFMetVsGJmass[0]         ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt);
			  h_PFMetOverSumEtVsGJmass[0]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt);

			  h_goodPV_TotalWt[1]        ->Fill(GoodVertex, Total_EvtWt);
			  h_nPhotons[1]              ->Fill(GoodIsoPhotons.size(), Total_EvtWt);
			  h_nJets[1]                 ->Fill(GoodIsoJets.size(), Total_EvtWt);

			  if(PC_Gen > -1 && JC_Gen > -1){
			    h_Gen_GJetInvtMass_bin40[0]   ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			    h_Gen_GJetInvtMass_VarBin[0]  ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			    h_Gen_GJetInvtMass_UnitBin[0] ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));

			    h_GJetMassResolution[0]    ->Fill(fabs(GetGenLevelInvtMass(PC_Gen, JC_Gen) - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			  }


			  if(Pass_GJInvtMass && IfMassCut){
			    h_CutFlow->Fill(10.5);
			    h_CutFlowWithWts->Fill(10.5, Total_EvtWt);

			    h_PhotonPt[1]              ->Fill((*phoEt)[PC], Total_EvtWt);
			    h_PhotonEta[1]             ->Fill((*phoSCEta)[PC], Total_EvtWt);
			    h_PhotonPhi[1]             ->Fill((*phoSCPhi)[PC], Total_EvtWt);
			    h_Photon_SigmaIEtaIEta[1]  ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC], Total_EvtWt);
			    h_Photon_R9[1]             ->Fill((*phoR9)[PC], Total_EvtWt);
			    h_Photon_HoverE[1]         ->Fill((*phoHoverE)[PC], Total_EvtWt);
			    h_Photon_EleVeto[1]        ->Fill((*phoEleVeto)[PC], Total_EvtWt);
			    h_Photon_CorrPFChIso[1]    ->Fill((*phoPFChIso)[PC], Total_EvtWt);
			    h_Photon_CorrPFNeuIso[1]   ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0), Total_EvtWt);
			    h_Photon_CorrPFPhoIso[1]   ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0), Total_EvtWt);

			    h_JetPt[1]                 ->Fill((*jetPt)[JC], Total_EvtWt);
			    h_JetEta[1]                ->Fill((*jetEta)[JC], Total_EvtWt);
			    h_JetPhi[1]                ->Fill((*jetPhi)[JC], Total_EvtWt);
			    h_Jet_NHEF[1]              ->Fill((*jetNHF)[JC], Total_EvtWt);
			    h_Jet_NEEF[1]              ->Fill((*jetNEF)[JC], Total_EvtWt);
			    h_Jet_NConst[1]            ->Fill((*jetNConstituents)[JC], Total_EvtWt);
			    h_Jet_CHEF[1]              ->Fill((*jetCHF)[JC], Total_EvtWt);
			    h_Jet_ChMult[1]            ->Fill((*jetNCH)[JC], Total_EvtWt);
			    h_Jet_CEEF[1]              ->Fill((*jetCEF)[JC], Total_EvtWt);

			    h_GJetInvtMass_bin40[1]    ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
			    h_GJetInvtMass_VarBin[1]   ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
			    h_GJetInvtMass_UnitBin[1]  ->Fill(GetInvtMass(PC, JC), Total_EvtWt);
			    h_GJet_dEta[1]             ->Fill(GetdEta((*phoSCEta)[PC], (*jetEta)[JC]), Total_EvtWt);
			    h_GJet_dPhi[1]             ->Fill(GetdPhi((*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt);
			    h_GJet_dR[1]               ->Fill(GetdR((*phoSCEta)[PC], (*jetEta)[JC], (*phoSCPhi)[PC], (*jetPhi)[JC]), Total_EvtWt);

			    h_PFMet[1]                 ->Fill(pfMET, Total_EvtWt);
			    h_PFMetVsGJmass[1]         ->Fill(GetInvtMass(PC, JC), pfMET, Total_EvtWt);
			    h_PFMetOverSumEtVsGJmass[1]->Fill(GetInvtMass(PC, JC), pfMET/pfMETsumEt, Total_EvtWt);

			    h_goodPV_TotalWt[2]        ->Fill(GoodVertex, Total_EvtWt);
			    h_nPhotons[2]              ->Fill(GoodIsoPhotons.size(), Total_EvtWt);
			    h_nJets[2]                 ->Fill(GoodIsoJets.size(), Total_EvtWt);

			    if(PC_Gen > -1 && JC_Gen > -1){
			      h_Gen_GJetInvtMass_bin40[1]   ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			      h_Gen_GJetInvtMass_VarBin[1]  ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));
			      h_Gen_GJetInvtMass_UnitBin[1] ->Fill(GetGenLevelInvtMass(PC_Gen, JC_Gen));

			      h_GJetMassResolution[1]    ->Fill(fabs(GetGenLevelInvtMass(PC_Gen, JC_Gen) - GetInvtMass(PC, JC))/GetInvtMass(PC, JC));
			    }


			    if(pfMET < 500){
			      h_CutFlow->Fill(11.5);
			      h_CutFlowWithWts->Fill(11.5, Total_EvtWt);

			      if(pfMET < 250){
				h_CutFlow->Fill(12.5);
				h_CutFlowWithWts->Fill(12.5, Total_EvtWt);;

				if(pfMET < 100){
				  h_CutFlow->Fill(13.5);
				  h_CutFlowWithWts->Fill(13.5, Total_EvtWt);

				  if(pfMET < 50){
				    h_CutFlow->Fill(14.5);
				    h_CutFlowWithWts->Fill(14.5, Total_EvtWt);

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

	      h_PassGenIsoMatch_EB->Fill((*phoEt)[PC_G], Total_EvtWt);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EB->Fill((*phoEt)[PC_L], Total_EvtWt);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EB->Fill((*phoEt)[PC_M], Total_EvtWt);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EB->Fill((*phoEt)[PC_T], Total_EvtWt);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EB->Fill((*phoEt)[PC_H], Total_EvtWt);
	      }
	    }
	    if(fabs((*phoSCEta)[PC_G]) < 2.5 && fabs((*phoSCEta)[PC_G]) >= 1.5666){//Endcap  

	      h_PassGenIsoMatch_EE->Fill((*phoEt)[PC_G], Total_EvtWt);

	      if(PC_L > -1 && PC_G == PC_L){
		h_PassPhIdLoose_EE->Fill((*phoEt)[PC_L], Total_EvtWt);
	      }

	      if(PC_M > -1 && PC_G == PC_M){
		h_PassPhIdMedium_EE->Fill((*phoEt)[PC_M], Total_EvtWt);
	      }

	      if(PC_T > -1 && PC_G == PC_T){
		h_PassPhIdTight_EE->Fill((*phoEt)[PC_T], Total_EvtWt);
	      }

	      if(PC_H > -1 && PC_G == PC_H){
		h_PassPhIdHighPt_EE->Fill((*phoEt)[PC_H], Total_EvtWt);
	      }
	    }	  
	  }
	}
      }

   }//for jentry
}//Loop()

EOF



cat >> PostAnalyzer_MC.h <<EOF 
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 01:52:07 2015 by ROOT version 6.02/05
// from TChain ggNtuplizer/EventTree/
//////////////////////////////////////////////////////////

#ifndef PostAnalyzer_MC_h
#define PostAnalyzer_MC_h

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
#include "TRFIOFile.h"
#include "TKDE.h" 

// Header file for the classes stored in the TTree if any.
#include "vector"

//c++ include files
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <map>
#include <set>

using namespace std;
using namespace ROOT;

class PostAnalyzer_MC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   //Variables defined by me

   //Outut root file
   TFile *file;

   //Bools
   Bool_t _debug;
   Bool_t IfMassCut;
   Bool_t Pass_HLT;
   Bool_t HasPrimaryVtx;
   Bool_t Pass_PhoPt;
   Bool_t Pass_PhoEtaEB;
   Bool_t Pass_JetPt;
   Bool_t Pass_JetEta;
   Bool_t Pass_GJdPhi;
   Bool_t Pass_GJdEta;
   Bool_t Pass_GJInvtMass;

   Bool_t IfGJet;

   //Ints, Doubles etc.
   Int_t GoodVertex;
   Int_t PC, JC;
   Int_t PC_G, PC_L, PC_M, PC_T, PC_H;
   Int_t PC_Gen, JC_Gen;

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

   //Only MC
   Double_t Lumi;
   Double_t Lumi_EvtWt, PU_EvtWt, PreBTag_EvtWt, Total_EvtWt;

   vector<Int_t> GoodIsoPhotons;
   vector<Int_t> GoodIsoBarrelPhotons;
   vector<Int_t> GoodIsoJets;

   //Histograms
   TH1F *h_PhotonPt[2];
   TH1F *h_PhotonEta[2];
   TH1F *h_PhotonPhi[2];
   TH1F *h_Photon_SigmaIEtaIEta[2];
   TH1F *h_Photon_R9[2];
   TH1F *h_Photon_HoverE[2];
   TH1F *h_Photon_EleVeto[2];
   TH1F *h_Photon_CorrPFChIso[2];
   TH1F *h_Photon_CorrPFNeuIso[2];
   TH1F *h_Photon_CorrPFPhoIso[2];

   TH1F *h_JetPt[2];
   TH1F *h_JetEta[2];
   TH1F *h_JetPhi[2];
   TH1F *h_Jet_NHEF[2];
   TH1F *h_Jet_NEEF[2];
   TH1F *h_Jet_NConst[2];
   TH1F *h_Jet_CHEF[2];
   TH1F *h_Jet_ChMult[2];
   TH1F *h_Jet_CEEF[2];

   TH1F *h_GJetInvtMass_bin40[2];
   TH1F *h_GJetInvtMass_VarBin[2];
   TH1F *h_GJetInvtMass_UnitBin[2];
   TH1F *h_GJet_dEta[2];
   TH1F *h_GJet_dPhi[2];
   TH1F *h_GJet_dR[2];

   TH1F *h_mass_X_bin1; //Invt mass divide by mean mass required for interploation of signal shapes.

   TH1F *h_PFMet[2];
   TH2F *h_PFMetVsGJmass[2];
   TH2F *h_PFMetOverSumEtVsGJmass[2];

   //Gen Info for MC only
   TH1F *h_Gen_GJetInvtMass_bin40[2];
   TH1F *h_Gen_GJetInvtMass_VarBin[2];
   TH1F *h_Gen_GJetInvtMass_UnitBin[2];
   TH1F *h_GJetMassResolution[2];

   //Photon Efficiency
   TH1F *h_PassGenIsoMatch_EB;
   TH1F *h_PassGenIsoMatch_EE;
   TH1F *h_PassPhIdLoose_EB;
   TH1F *h_PassPhIdLoose_EE;
   TH1F *h_PassPhIdMedium_EB;
   TH1F *h_PassPhIdMedium_EE;
   TH1F *h_PassPhIdTight_EB;
   TH1F *h_PassPhIdTight_EE;
   TH1F *h_PassPhIdHighPt_EB;
   TH1F *h_PassPhIdHighPt_EE;

   //Pileup
   TH1F *h_trueNumofInt;
   TH1F *h_goodPV;
   TH1F *h_goodPV_LumiWt;
   TH1F *h_goodPV_PUWt;
   TH1F *h_goodPV_TotalWt[3];

   //Others
   TH1F *h_nPhotons[3];
   TH1F *h_nJets[3];

   TH1F *h_PC;
   TH1F *h_JC;

   //Cut Flow
   TH1F *h_CutFlow;
   TH1F *h_CutFlowWithWts;

   //Pileup Reweighting
   TH1F *h_DataPUNormDist;
   TH1F *h_MCPUNormDist;
   TH1F *h_PUScaleFactor;


   //Variables of TTree::EventTree
   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   vector<int>     *nTrksPV;
   vector<float>   *vtx;
   vector<float>   *vty;
   vector<float>   *vtz;
   vector<float>   *vrho;
   vector<int>     *vndof;
   vector<float>   *vchi2;
   vector<bool>    *isFake;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue; //getTrueNumofInteractions, Same for a BX(both in and out of time), so this vector has all the values same for an event.
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<int>     *mcIndex;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR04;
   Int_t           metFilters;
   Float_t         genMET;
   Float_t         genMETPhi;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1MESUp;
   Float_t         pfMET_T1MESDo;
   Float_t         pfMET_T1EESUp;
   Float_t         pfMET_T1EESDo;
   Float_t         pfMET_T1PESUp;
   Float_t         pfMET_T1PESDo;
   Float_t         pfMET_T1TESUp;
   Float_t         pfMET_T1TESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEn;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoSigmaIEtaIPhi;
   vector<float>   *phoSigmaIPhiIPhi;
   vector<float>   *phoE1x3;
   vector<float>   *phoE2x2;
   vector<float>   *phoE2x5Max;
   vector<float>   *phoE5x5;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phomaxXtalenergyFull5x5;
   vector<float>   *phoseedTimeFull5x5;
   vector<float>   *phomaxXtalenergy;
   vector<float>   *phoseedTime;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE1x3Full5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE2x5MaxFull5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoSeedBCE;
   vector<float>   *phoSeedBCEta;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFChIsoFrix1;
   vector<float>   *phoPFChIsoFrix2;
   vector<float>   *phoPFChIsoFrix3;
   vector<float>   *phoPFChIsoFrix4;
   vector<float>   *phoPFChIsoFrix5;
   vector<float>   *phoPFChIsoFrix6;
   vector<float>   *phoPFChIsoFrix7;
   vector<float>   *phoPFChIsoFrix8;
   vector<float>   *phoPFPhoIsoFrix1;
   vector<float>   *phoPFPhoIsoFrix2;
   vector<float>   *phoPFPhoIsoFrix3;
   vector<float>   *phoPFPhoIsoFrix4;
   vector<float>   *phoPFPhoIsoFrix5;
   vector<float>   *phoPFPhoIsoFrix6;
   vector<float>   *phoPFPhoIsoFrix7;
   vector<float>   *phoPFPhoIsoFrix8;
   vector<float>   *phoPFNeuIsoFrix1;
   vector<float>   *phoPFNeuIsoFrix2;
   vector<float>   *phoPFNeuIsoFrix3;
   vector<float>   *phoPFNeuIsoFrix4;
   vector<float>   *phoPFNeuIsoFrix5;
   vector<float>   *phoPFNeuIsoFrix6;
   vector<float>   *phoPFNeuIsoFrix7;
   vector<float>   *phoPFNeuIsoFrix8;
   vector<float>   *phoEcalRecHitSumEtConeDR03;
   vector<float>   *phohcalDepth1TowerSumEtConeDR03;
   vector<float>   *phohcalDepth2TowerSumEtConeDR03;
   vector<float>   *phohcalTowerSumEtConeDR03;
   vector<float>   *photrkSumPtHollowConeDR03;
   vector<float>   *phoIDMVA;
   vector<int>     *phoFiredSingleTrgs;
   vector<int>     *phoFiredDoubleTrgs;
   vector<unsigned short> *phoIDbit;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleESEn;
   vector<float>   *eleESEnP1;
   vector<float>   *eleESEnP2;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaAtCalo;
   vector<float>   *eleSigmaIEtaIEta;
   vector<float>   *eleSigmaIEtaIPhi;
   vector<float>   *eleSigmaIPhiIPhi;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *eleIDMVANonTrg;
   vector<float>   *eleIDMVATrg;
   vector<float>   *eledEtaseedAtVtx;
   vector<float>   *eleE1x5;
   vector<float>   *eleE2x5;
   vector<float>   *eleE5x5;
   vector<float>   *eleE1x5Full5x5;
   vector<float>   *eleE2x5Full5x5;
   vector<float>   *eleE5x5Full5x5;
   vector<float>   *eleR9Full5x5;
   vector<int>     *eleEcalDrivenSeed;
   vector<float>   *eleDr03EcalRecHitSumEt;
   vector<float>   *eleDr03HcalDepth1TowerSumEt;
   vector<float>   *eleDr03HcalDepth2TowerSumEt;
   vector<float>   *eleDr03HcalTowerSumEt;
   vector<float>   *eleDr03TkSumPt;
   vector<float>   *elecaloEnergy;
   vector<float>   *eleTrkdxy;
   vector<float>   *eleKFHits;
   vector<float>   *eleKFChi2;
   vector<vector<float> > *eleGSFPt;
   vector<vector<float> > *eleGSFEta;
   vector<vector<float> > *eleGSFPhi;
   vector<vector<float> > *eleGSFCharge;
   vector<vector<int> > *eleGSFHits;
   vector<vector<int> > *eleGSFMissHits;
   vector<vector<int> > *eleGSFNHitsMax;
   vector<vector<float> > *eleGSFVtxProb;
   vector<vector<float> > *eleGSFlxyPV;
   vector<vector<float> > *eleGSFlxyBS;
   vector<vector<float> > *eleBCEn;
   vector<vector<float> > *eleBCEta;
   vector<vector<float> > *eleBCPhi;
   vector<vector<float> > *eleBCS25;
   vector<vector<float> > *eleBCS15;
   vector<vector<float> > *eleBCSieie;
   vector<vector<float> > *eleBCSieip;
   vector<vector<float> > *eleBCSipip;
   vector<int>     *eleFiredTrgs;
   vector<unsigned short> *eleIDbit;
   vector<float>   *eleESEnP1Raw;
   vector<float>   *eleESEnP2Raw;
   Int_t           nGSFTrk;
   vector<float>   *gsfPt;
   vector<float>   *gsfEta;
   vector<float>   *gsfPhi;
   Int_t           npfHF;
   vector<float>   *pfHFEn;
   vector<float>   *pfHFECALEn;
   vector<float>   *pfHFHCALEn;
   vector<float>   *pfHFPt;
   vector<float>   *pfHFEta;
   vector<float>   *pfHFPhi;
   vector<float>   *pfHFIso;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<bool>    *muIsLooseID;
   vector<bool>    *muIsMediumID;
   vector<bool>    *muIsTightID;
   vector<bool>    *muIsSoftID;
   vector<bool>    *muIsHighPtID;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<int>     *muFiredTrgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   Int_t           nTau;
   vector<bool>    *pfTausDiscriminationByDecayModeFinding;
   vector<bool>    *pfTausDiscriminationByDecayModeFindingNewDMs;
   vector<bool>    *tauByMVA5LooseElectronRejection;
   vector<bool>    *tauByMVA5MediumElectronRejection;
   vector<bool>    *tauByMVA5TightElectronRejection;
   vector<bool>    *tauByMVA5VTightElectronRejection;
   vector<bool>    *tauByLooseMuonRejection3;
   vector<bool>    *tauByTightMuonRejection3;
   vector<bool>    *tauByLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *tauByMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *tauByTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<bool>    *tauByVLooseIsolationMVA3oldDMwLT;
   vector<bool>    *tauByLooseIsolationMVA3oldDMwLT;
   vector<bool>    *tauByMediumIsolationMVA3oldDMwLT;
   vector<bool>    *tauByTightIsolationMVA3oldDMwLT;
   vector<bool>    *tauByVTightIsolationMVA3oldDMwLT;
   vector<bool>    *tauByVVTightIsolationMVA3oldDMwLT;
   vector<float>   *tauByIsolationMVA3oldDMwLTraw;
   vector<bool>    *tauByLooseIsolationMVA3newDMwLT;
   vector<bool>    *tauByVLooseIsolationMVA3newDMwLT;
   vector<bool>    *tauByMediumIsolationMVA3newDMwLT;
   vector<bool>    *tauByTightIsolationMVA3newDMwLT;
   vector<bool>    *tauByVTightIsolationMVA3newDMwLT;
   vector<bool>    *tauByVVTightIsolationMVA3newDMwLT;
   vector<float>   *tauByIsolationMVA3newDMwLTraw;
   vector<float>   *tauEta;
   vector<float>   *tauPhi;
   vector<float>   *tauPt;
   vector<float>   *tauEt;
   vector<float>   *tauCharge;
   vector<float>   *tauP;
   vector<float>   *tauPx;
   vector<float>   *tauPy;
   vector<float>   *tauPz;
   vector<float>   *tauVz;
   vector<float>   *tauEnergy;
   vector<float>   *tauMass;
   vector<float>   *tauDxy;
   vector<float>   *tauZImpact;
   vector<int>     *tauDecayMode;
   vector<bool>    *tauLeadChargedHadronExists;
   vector<float>   *tauLeadChargedHadronEta;
   vector<float>   *tauLeadChargedHadronPhi;
   vector<float>   *tauLeadChargedHadronPt;
   vector<float>   *tauChargedIsoPtSum;
   vector<float>   *tauNeutralIsoPtSum;
   vector<float>   *tauPuCorrPtSum;
   vector<float>   *tauNumSignalPFChargedHadrCands;
   vector<float>   *tauNumSignalPFNeutrHadrCands;
   vector<float>   *tauNumSignalPFGammaCands;
   vector<float>   *tauNumSignalPFCands;
   vector<float>   *tauNumIsolationPFChargedHadrCands;
   vector<float>   *tauNumIsolationPFNeutrHadrCands;
   vector<float>   *tauNumIsolationPFGammaCands;
   vector<float>   *tauNumIsolationPFCands;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetArea;
   vector<float>   *jetpfCombinedInclusiveSecondaryVertexV2BJetTags;
   vector<float>   *jetJetProbabilityBJetTags;
   vector<float>   *jetpfCombinedMVABJetTags;
   vector<int>     *jetPartonID;
   vector<int>     *jetGenJetIndex;
   vector<float>   *jetGenJetEn;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenEn;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   vector<bool>    *jetPFLooseId;
   vector<float>   *jetPUidFullDiscriminant;
   vector<float>   *jetJECUnc;
   vector<int>     *jetFiredTrgs;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<int>     *jetNCH;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConstituents;
   Int_t           nAK8Jet;
   vector<float>   *AK8JetPt;
   vector<float>   *AK8JetEn;
   vector<float>   *AK8JetRawPt;
   vector<float>   *AK8JetRawEn;
   vector<float>   *AK8JetEta;
   vector<float>   *AK8JetPhi;
   vector<float>   *AK8JetMass;
   vector<float>   *AK8Jet_tau1;
   vector<float>   *AK8Jet_tau2;
   vector<float>   *AK8Jet_tau3;
   vector<float>   *AK8JetCHF;
   vector<float>   *AK8JetNHF;
   vector<float>   *AK8JetCEF;
   vector<float>   *AK8JetNEF;
   vector<int>     *AK8JetNCH;
   vector<int>     *AK8Jetnconstituents;
   vector<bool>    *AK8JetPFLooseId;
   vector<float>   *AK8CHSSoftDropJetMass;
   vector<float>   *AK8JetpfBoostedDSVBTag;
   vector<float>   *AK8JetJECUnc;
   vector<int>     *AK8JetPartonID;
   vector<int>     *AK8JetGenJetIndex;
   vector<float>   *AK8JetGenJetEn;
   vector<float>   *AK8JetGenJetPt;
   vector<float>   *AK8JetGenJetEta;
   vector<float>   *AK8JetGenJetPhi;
   vector<int>     *AK8JetGenPartonID;
   vector<float>   *AK8JetGenEn;
   vector<float>   *AK8JetGenPt;
   vector<float>   *AK8JetGenEta;
   vector<float>   *AK8JetGenPhi;
   vector<int>     *AK8JetGenPartonMomID;
   vector<int>     *nAK8softdropSubjet;
   vector<vector<float> > *AK8softdropSubjetPt;
   vector<vector<float> > *AK8softdropSubjetEta;
   vector<vector<float> > *AK8softdropSubjetPhi;
   vector<vector<float> > *AK8softdropSubjetMass;
   vector<vector<float> > *AK8softdropSubjetE;
   vector<vector<int> > *AK8softdropSubjetCharge;
   vector<vector<int> > *AK8softdropSubjetFlavour;
   vector<vector<float> > *AK8softdropSubjetCSV;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_vrho;   //!
   TBranch        *b_vndof;   //!
   TBranch        *b_vchi2;   //!
   TBranch        *b_isFake;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1MESUp;   //!
   TBranch        *b_pfMET_T1MESDo;   //!
   TBranch        *b_pfMET_T1EESUp;   //!
   TBranch        *b_pfMET_T1EESDo;   //!
   TBranch        *b_pfMET_T1PESUp;   //!
   TBranch        *b_pfMET_T1PESDo;   //!
   TBranch        *b_pfMET_T1TESUp;   //!
   TBranch        *b_pfMET_T1TESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phomaxXtalenergyFull5x5;   //!
   TBranch        *b_phoseedTimeFull5x5;   //!
   TBranch        *b_phomaxXtalenergy;   //!
   TBranch        *b_phoseedTime;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE1x3Full5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE2x5MaxFull5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoSeedBCE;   //!
   TBranch        *b_phoSeedBCEta;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFChIsoFrix1;   //!
   TBranch        *b_phoPFChIsoFrix2;   //!
   TBranch        *b_phoPFChIsoFrix3;   //!
   TBranch        *b_phoPFChIsoFrix4;   //!
   TBranch        *b_phoPFChIsoFrix5;   //!
   TBranch        *b_phoPFChIsoFrix6;   //!
   TBranch        *b_phoPFChIsoFrix7;   //!
   TBranch        *b_phoPFChIsoFrix8;   //!
   TBranch        *b_phoPFPhoIsoFrix1;   //!
   TBranch        *b_phoPFPhoIsoFrix2;   //!
   TBranch        *b_phoPFPhoIsoFrix3;   //!
   TBranch        *b_phoPFPhoIsoFrix4;   //!
   TBranch        *b_phoPFPhoIsoFrix5;   //!
   TBranch        *b_phoPFPhoIsoFrix6;   //!
   TBranch        *b_phoPFPhoIsoFrix7;   //!
   TBranch        *b_phoPFPhoIsoFrix8;   //!
   TBranch        *b_phoPFNeuIsoFrix1;   //!
   TBranch        *b_phoPFNeuIsoFrix2;   //!
   TBranch        *b_phoPFNeuIsoFrix3;   //!
   TBranch        *b_phoPFNeuIsoFrix4;   //!
   TBranch        *b_phoPFNeuIsoFrix5;   //!
   TBranch        *b_phoPFNeuIsoFrix6;   //!
   TBranch        *b_phoPFNeuIsoFrix7;   //!
   TBranch        *b_phoPFNeuIsoFrix8;   //!
   TBranch        *b_phoEcalRecHitSumEtConeDR03;   //!
   TBranch        *b_phohcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_phohcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_phohcalTowerSumEtConeDR03;   //!
   TBranch        *b_photrkSumPtHollowConeDR03;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleESEnP1;   //!
   TBranch        *b_eleESEnP2;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_eleIDMVANonTrg;   //!
   TBranch        *b_eleIDMVATrg;   //!
   TBranch        *b_eledEtaseedAtVtx;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE2x5;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE1x5Full5x5;   //!
   TBranch        *b_eleE2x5Full5x5;   //!
   TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleDr03EcalRecHitSumEt;   //!
   TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_eleDr03HcalTowerSumEt;   //!
   TBranch        *b_eleDr03TkSumPt;   //!
   TBranch        *b_elecaloEnergy;   //!
   TBranch        *b_eleTrkdxy;   //!
   TBranch        *b_eleKFHits;   //!
   TBranch        *b_eleKFChi2;   //!
   TBranch        *b_eleGSFPt;   //!
   TBranch        *b_eleGSFEta;   //!
   TBranch        *b_eleGSFPhi;   //!
   TBranch        *b_eleGSFCharge;   //!
   TBranch        *b_eleGSFHits;   //!
   TBranch        *b_eleGSFMissHits;   //!
   TBranch        *b_eleGSFNHitsMax;   //!
   TBranch        *b_eleGSFVtxProb;   //!
   TBranch        *b_eleGSFlxyPV;   //!
   TBranch        *b_eleGSFlxyBS;   //!
   TBranch        *b_eleBCEn;   //!
   TBranch        *b_eleBCEta;   //!
   TBranch        *b_eleBCPhi;   //!
   TBranch        *b_eleBCS25;   //!
   TBranch        *b_eleBCS15;   //!
   TBranch        *b_eleBCSieie;   //!
   TBranch        *b_eleBCSieip;   //!
   TBranch        *b_eleBCSipip;   //!
   TBranch        *b_eleFiredTrgs;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleESEnP1Raw;   //!
   TBranch        *b_eleESEnP2Raw;   //!
   TBranch        *b_nGSFTrk;   //!
   TBranch        *b_gsfPt;   //!
   TBranch        *b_gsfEta;   //!
   TBranch        *b_gsfPhi;   //!
   TBranch        *b_npfHF;   //!
   TBranch        *b_pfHFEn;   //!
   TBranch        *b_pfHFECALEn;   //!
   TBranch        *b_pfHFHCALEn;   //!
   TBranch        *b_pfHFPt;   //!
   TBranch        *b_pfHFEta;   //!
   TBranch        *b_pfHFPhi;   //!
   TBranch        *b_pfHFIso;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIsLooseID;   //!
   TBranch        *b_muIsMediumID;   //!
   TBranch        *b_muIsTightID;   //!
   TBranch        *b_muIsSoftID;   //!
   TBranch        *b_muIsHighPtID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_pfTausDiscriminationByDecayModeFinding;   //!
   TBranch        *b_pfTausDiscriminationByDecayModeFindingNewDMs;   //!
   TBranch        *b_tauByMVA5LooseElectronRejection;   //!
   TBranch        *b_tauByMVA5MediumElectronRejection;   //!
   TBranch        *b_tauByMVA5TightElectronRejection;   //!
   TBranch        *b_tauByMVA5VTightElectronRejection;   //!
   TBranch        *b_tauByLooseMuonRejection3;   //!
   TBranch        *b_tauByTightMuonRejection3;   //!
   TBranch        *b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauByTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tauByVLooseIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByLooseIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByMediumIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByVTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByVVTightIsolationMVA3oldDMwLT;   //!
   TBranch        *b_tauByIsolationMVA3oldDMwLTraw;   //!
   TBranch        *b_tauByLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByVLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByMediumIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByVTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByVVTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_tauByIsolationMVA3newDMwLTraw;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEt;   //!
   TBranch        *b_tauCharge;   //!
   TBranch        *b_tauP;   //!
   TBranch        *b_tauPx;   //!
   TBranch        *b_tauPy;   //!
   TBranch        *b_tauPz;   //!
   TBranch        *b_tauVz;   //!
   TBranch        *b_tauEnergy;   //!
   TBranch        *b_tauMass;   //!
   TBranch        *b_tauDxy;   //!
   TBranch        *b_tauZImpact;   //!
   TBranch        *b_tauDecayMode;   //!
   TBranch        *b_tauLeadChargedHadronExists;   //!
   TBranch        *b_tauLeadChargedHadronEta;   //!
   TBranch        *b_tauLeadChargedHadronPhi;   //!
   TBranch        *b_tauLeadChargedHadronPt;   //!
   TBranch        *b_tauChargedIsoPtSum;   //!
   TBranch        *b_tauNeutralIsoPtSum;   //!
   TBranch        *b_tauPuCorrPtSum;   //!
   TBranch        *b_tauNumSignalPFChargedHadrCands;   //!
   TBranch        *b_tauNumSignalPFNeutrHadrCands;   //!
   TBranch        *b_tauNumSignalPFGammaCands;   //!
   TBranch        *b_tauNumSignalPFCands;   //!
   TBranch        *b_tauNumIsolationPFChargedHadrCands;   //!
   TBranch        *b_tauNumIsolationPFNeutrHadrCands;   //!
   TBranch        *b_tauNumIsolationPFGammaCands;   //!
   TBranch        *b_tauNumIsolationPFCands;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetpfCombinedMVABJetTags;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetPUidFullDiscriminant;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetFiredTrgs;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_nAK8Jet;   //!
   TBranch        *b_AK8JetPt;   //!
   TBranch        *b_AK8JetEn;   //!
   TBranch        *b_AK8JetRawPt;   //!
   TBranch        *b_AK8JetRawEn;   //!
   TBranch        *b_AK8JetEta;   //!
   TBranch        *b_AK8JetPhi;   //!
   TBranch        *b_AK8JetMass;   //!
   TBranch        *b_AK8Jet_tau1;   //!
   TBranch        *b_AK8Jet_tau2;   //!
   TBranch        *b_AK8Jet_tau3;   //!
   TBranch        *b_AK8JetCHF;   //!
   TBranch        *b_AK8JetNHF;   //!
   TBranch        *b_AK8JetCEF;   //!
   TBranch        *b_AK8JetNEF;   //!
   TBranch        *b_AK8JetNCH;   //!
   TBranch        *b_AK8Jetnconstituents;   //!
   TBranch        *b_AK8JetPFLooseId;   //!
   TBranch        *b_AK8CHSSoftDropJetMass;   //!
   TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
   TBranch        *b_AK8JetJECUnc;   //!
   TBranch        *b_AK8JetPartonID;   //!
   TBranch        *b_AK8JetGenJetIndex;   //!
   TBranch        *b_AK8JetGenJetEn;   //!
   TBranch        *b_AK8JetGenJetPt;   //!
   TBranch        *b_AK8JetGenJetEta;   //!
   TBranch        *b_AK8JetGenJetPhi;   //!
   TBranch        *b_AK8JetGenPartonID;   //!
   TBranch        *b_AK8JetGenEn;   //!
   TBranch        *b_AK8JetGenPt;   //!
   TBranch        *b_AK8JetGenEta;   //!
   TBranch        *b_AK8JetGenPhi;   //!
   TBranch        *b_AK8JetGenPartonMomID;   //!
   TBranch        *b_nAK8softdropSubjet;   //!
   TBranch        *b_AK8softdropSubjetPt;   //!
   TBranch        *b_AK8softdropSubjetEta;   //!
   TBranch        *b_AK8softdropSubjetPhi;   //!
   TBranch        *b_AK8softdropSubjetMass;   //!
   TBranch        *b_AK8softdropSubjetE;   //!
   TBranch        *b_AK8softdropSubjetCharge;   //!
   TBranch        *b_AK8softdropSubjetFlavour;   //!
   TBranch        *b_AK8softdropSubjetCSV;   //!

   PostAnalyzer_MC(TTree *tree=0);
   virtual ~PostAnalyzer_MC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //User-Defined functions
   virtual Bool_t   GoodPrimaryVtx(Int_t &GoodVertex);
   virtual Bool_t   CutBasedPhotonID(Int_t ipho, TString phoWP);
   virtual Double_t EANeutralHadrons(Double_t eta);
   virtual Double_t EAPhotons(Double_t eta);
   virtual Int_t    FirstGoodPhoton(TString phoWP); 
   virtual vector<Int_t> GoodPhotons(TString phoWP); 
   virtual Bool_t   TightJetId(Int_t iJet);
   virtual Int_t    FirstGoodJet(Int_t pc);
   virtual vector<Int_t> GoodJets(Int_t pc);

   virtual Double_t GetdEta(Double_t eta1, Double_t eta2);
   virtual Double_t GetdPhi(Double_t phi1, Double_t phi2);
   virtual Double_t GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
   virtual Double_t GetInvtMass(Int_t ph, Int_t jet);

   virtual Int_t    MatchedGenPhotonToReco(Int_t pc);
   virtual Int_t    MatchedGenJetToReco(Int_t jc);
   virtual Double_t GetGenLevelInvtMass(Int_t pc_gen, Int_t jc_gen);
   virtual Int_t    MatchedRecoPhotonToGen_WithGenIsoCut();

   virtual Double_t MCFM_kFact(Double_t phpt);
   virtual void     PileupReWeighting();
   virtual Double_t PUWeights(Float_t npv);

   virtual void     BookHistograms();
   virtual void     WriteHistograms();

   //Extra Functions defined for high Pt Photon ID
   virtual Bool_t   HighPtPhotonID(Int_t ipho);
   virtual Int_t    FirstHighPtIDPhoton(); 
   virtual vector<Int_t> HighPtIDPhotons();
   virtual Double_t EAPhotons_HighPtID(Double_t eta);


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
      f->GetObject("ggNtuplizer/EventTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ggNtuplizer/EventTree","");
/*  
      //chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/DiJet/AOD_QCD_Pt_170to300_1.root/ggNtuplizer/EventTree"); 
chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/DiJet/QCD_Pt_170to300/AOD_QCD_Pt_170to300_10.root/ggNtuplizer/EventTree");
chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/DiJet/QCD_Pt_170to300/AOD_QCD_Pt_170to300_100.root/ggNtuplizer/EventTree");
chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/DiJet/QCD_Pt_170to300/AOD_QCD_Pt_170to300_101.root/ggNtuplizer/EventTree");
      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/Signal_Qstar/QstarToGJ_M1000_f1p0/AOD_QstarToGJ_M1000_f1p0_3.root/ggNtuplizer/EventTree");
      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/Signal_Qstar/QstarToGJ_M1000_f1p0/AOD_QstarToGJ_M1000_f1p0_5.root/ggNtuplizer/EventTree");
      chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/Signal_Qstar/QstarToGJ_M1000_f1p0/AOD_QstarToGJ_M1000_f1p0_4.root/ggNtuplizer/EventTree");*/

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

        chain->Add(FullPathInputFile+"/ggNtuplizer/EventTree");

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
   nTrksPV = 0;
   vtx = 0;
   vty = 0;
   vtz = 0;
   vrho = 0;
   vndof = 0;
   vchi2 = 0;
   isFake = 0;
   pdf = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoSigmaIEtaIEta = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoE1x3 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE5x5 = 0;
   phoESEffSigmaRR = 0;
   phomaxXtalenergyFull5x5 = 0;
   phoseedTimeFull5x5 = 0;
   phomaxXtalenergy = 0;
   phoseedTime = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE1x3Full5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE2x5MaxFull5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoSeedBCE = 0;
   phoSeedBCEta = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoPFChIsoFrix1 = 0;
   phoPFChIsoFrix2 = 0;
   phoPFChIsoFrix3 = 0;
   phoPFChIsoFrix4 = 0;
   phoPFChIsoFrix5 = 0;
   phoPFChIsoFrix6 = 0;
   phoPFChIsoFrix7 = 0;
   phoPFChIsoFrix8 = 0;
   phoPFPhoIsoFrix1 = 0;
   phoPFPhoIsoFrix2 = 0;
   phoPFPhoIsoFrix3 = 0;
   phoPFPhoIsoFrix4 = 0;
   phoPFPhoIsoFrix5 = 0;
   phoPFPhoIsoFrix6 = 0;
   phoPFPhoIsoFrix7 = 0;
   phoPFPhoIsoFrix8 = 0;
   phoPFNeuIsoFrix1 = 0;
   phoPFNeuIsoFrix2 = 0;
   phoPFNeuIsoFrix3 = 0;
   phoPFNeuIsoFrix4 = 0;
   phoPFNeuIsoFrix5 = 0;
   phoPFNeuIsoFrix6 = 0;
   phoPFNeuIsoFrix7 = 0;
   phoPFNeuIsoFrix8 = 0;
   phoEcalRecHitSumEtConeDR03 = 0;
   phohcalDepth1TowerSumEtConeDR03 = 0;
   phohcalDepth2TowerSumEtConeDR03 = 0;
   phohcalTowerSumEtConeDR03 = 0;
   photrkSumPtHollowConeDR03 = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoIDbit = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIPhi = 0;
   eleSigmaIPhiIPhi = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   eleIDMVANonTrg = 0;
   eleIDMVATrg = 0;
   eledEtaseedAtVtx = 0;
   eleE1x5 = 0;
   eleE2x5 = 0;
   eleE5x5 = 0;
   eleE1x5Full5x5 = 0;
   eleE2x5Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   eleDr03EcalRecHitSumEt = 0;
   eleDr03HcalDepth1TowerSumEt = 0;
   eleDr03HcalDepth2TowerSumEt = 0;
   eleDr03HcalTowerSumEt = 0;
   eleDr03TkSumPt = 0;
   elecaloEnergy = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFPt = 0;
   eleGSFEta = 0;
   eleGSFPhi = 0;
   eleGSFCharge = 0;
   eleGSFHits = 0;
   eleGSFMissHits = 0;
   eleGSFNHitsMax = 0;
   eleGSFVtxProb = 0;
   eleGSFlxyPV = 0;
   eleGSFlxyBS = 0;
   eleBCEn = 0;
   eleBCEta = 0;
   eleBCPhi = 0;
   eleBCS25 = 0;
   eleBCS15 = 0;
   eleBCSieie = 0;
   eleBCSieip = 0;
   eleBCSipip = 0;
   eleFiredTrgs = 0;
   eleIDbit = 0;
   eleESEnP1Raw = 0;
   eleESEnP2Raw = 0;
   gsfPt = 0;
   gsfEta = 0;
   gsfPhi = 0;
   pfHFEn = 0;
   pfHFECALEn = 0;
   pfHFHCALEn = 0;
   pfHFPt = 0;
   pfHFEta = 0;
   pfHFPhi = 0;
   pfHFIso = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIsLooseID = 0;
   muIsMediumID = 0;
   muIsTightID = 0;
   muIsSoftID = 0;
   muIsHighPtID = 0;
   muD0 = 0;
   muDz = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muFiredTrgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   pfTausDiscriminationByDecayModeFinding = 0;
   pfTausDiscriminationByDecayModeFindingNewDMs = 0;
   tauByMVA5LooseElectronRejection = 0;
   tauByMVA5MediumElectronRejection = 0;
   tauByMVA5TightElectronRejection = 0;
   tauByMVA5VTightElectronRejection = 0;
   tauByLooseMuonRejection3 = 0;
   tauByTightMuonRejection3 = 0;
   tauByLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauByMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauByTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   tauByVLooseIsolationMVA3oldDMwLT = 0;
   tauByLooseIsolationMVA3oldDMwLT = 0;
   tauByMediumIsolationMVA3oldDMwLT = 0;
   tauByTightIsolationMVA3oldDMwLT = 0;
   tauByVTightIsolationMVA3oldDMwLT = 0;
   tauByVVTightIsolationMVA3oldDMwLT = 0;
   tauByIsolationMVA3oldDMwLTraw = 0;
   tauByLooseIsolationMVA3newDMwLT = 0;
   tauByVLooseIsolationMVA3newDMwLT = 0;
   tauByMediumIsolationMVA3newDMwLT = 0;
   tauByTightIsolationMVA3newDMwLT = 0;
   tauByVTightIsolationMVA3newDMwLT = 0;
   tauByVVTightIsolationMVA3newDMwLT = 0;
   tauByIsolationMVA3newDMwLTraw = 0;
   tauEta = 0;
   tauPhi = 0;
   tauPt = 0;
   tauEt = 0;
   tauCharge = 0;
   tauP = 0;
   tauPx = 0;
   tauPy = 0;
   tauPz = 0;
   tauVz = 0;
   tauEnergy = 0;
   tauMass = 0;
   tauDxy = 0;
   tauZImpact = 0;
   tauDecayMode = 0;
   tauLeadChargedHadronExists = 0;
   tauLeadChargedHadronEta = 0;
   tauLeadChargedHadronPhi = 0;
   tauLeadChargedHadronPt = 0;
   tauChargedIsoPtSum = 0;
   tauNeutralIsoPtSum = 0;
   tauPuCorrPtSum = 0;
   tauNumSignalPFChargedHadrCands = 0;
   tauNumSignalPFNeutrHadrCands = 0;
   tauNumSignalPFGammaCands = 0;
   tauNumSignalPFCands = 0;
   tauNumIsolationPFChargedHadrCands = 0;
   tauNumIsolationPFNeutrHadrCands = 0;
   tauNumIsolationPFGammaCands = 0;
   tauNumIsolationPFCands = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetArea = 0;
   jetpfCombinedInclusiveSecondaryVertexV2BJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetpfCombinedMVABJetTags = 0;
   jetPartonID = 0;
   jetGenJetIndex = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetPFLooseId = 0;
   jetPUidFullDiscriminant = 0;
   jetJECUnc = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   AK8JetPt = 0;
   AK8JetEn = 0;
   AK8JetRawPt = 0;
   AK8JetRawEn = 0;
   AK8JetEta = 0;
   AK8JetPhi = 0;
   AK8JetMass = 0;
   AK8Jet_tau1 = 0;
   AK8Jet_tau2 = 0;
   AK8Jet_tau3 = 0;
   AK8JetCHF = 0;
   AK8JetNHF = 0;
   AK8JetCEF = 0;
   AK8JetNEF = 0;
   AK8JetNCH = 0;
   AK8Jetnconstituents = 0;
   AK8JetPFLooseId = 0;
   AK8CHSSoftDropJetMass = 0;
   AK8JetpfBoostedDSVBTag = 0;
   AK8JetJECUnc = 0;
   AK8JetPartonID = 0;
   AK8JetGenJetIndex = 0;
   AK8JetGenJetEn = 0;
   AK8JetGenJetPt = 0;
   AK8JetGenJetEta = 0;
   AK8JetGenJetPhi = 0;
   AK8JetGenPartonID = 0;
   AK8JetGenEn = 0;
   AK8JetGenPt = 0;
   AK8JetGenEta = 0;
   AK8JetGenPhi = 0;
   AK8JetGenPartonMomID = 0;
   nAK8softdropSubjet = 0;
   AK8softdropSubjetPt = 0;
   AK8softdropSubjetEta = 0;
   AK8softdropSubjetPhi = 0;
   AK8softdropSubjetMass = 0;
   AK8softdropSubjetE = 0;
   AK8softdropSubjetCharge = 0;
   AK8softdropSubjetFlavour = 0;
   AK8softdropSubjetCSV = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("vrho", &vrho, &b_vrho);
   fChain->SetBranchAddress("vndof", &vndof, &b_vndof);
   fChain->SetBranchAddress("vchi2", &vchi2, &b_vchi2);
   fChain->SetBranchAddress("isFake", &isFake, &b_isFake);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1MESUp", &pfMET_T1MESUp, &b_pfMET_T1MESUp);
   fChain->SetBranchAddress("pfMET_T1MESDo", &pfMET_T1MESDo, &b_pfMET_T1MESDo);
   fChain->SetBranchAddress("pfMET_T1EESUp", &pfMET_T1EESUp, &b_pfMET_T1EESUp);
   fChain->SetBranchAddress("pfMET_T1EESDo", &pfMET_T1EESDo, &b_pfMET_T1EESDo);
   fChain->SetBranchAddress("pfMET_T1PESUp", &pfMET_T1PESUp, &b_pfMET_T1PESUp);
   fChain->SetBranchAddress("pfMET_T1PESDo", &pfMET_T1PESDo, &b_pfMET_T1PESDo);
   fChain->SetBranchAddress("pfMET_T1TESUp", &pfMET_T1TESUp, &b_pfMET_T1TESUp);
   fChain->SetBranchAddress("pfMET_T1TESDo", &pfMET_T1TESDo, &b_pfMET_T1TESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phomaxXtalenergyFull5x5", &phomaxXtalenergyFull5x5, &b_phomaxXtalenergyFull5x5);
   fChain->SetBranchAddress("phoseedTimeFull5x5", &phoseedTimeFull5x5, &b_phoseedTimeFull5x5);
   fChain->SetBranchAddress("phomaxXtalenergy", &phomaxXtalenergy, &b_phomaxXtalenergy);
   fChain->SetBranchAddress("phoseedTime", &phoseedTime, &b_phoseedTime);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE2x5MaxFull5x5", &phoE2x5MaxFull5x5, &b_phoE2x5MaxFull5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoSeedBCE", &phoSeedBCE, &b_phoSeedBCE);
   fChain->SetBranchAddress("phoSeedBCEta", &phoSeedBCEta, &b_phoSeedBCEta);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFChIsoFrix1", &phoPFChIsoFrix1, &b_phoPFChIsoFrix1);
   fChain->SetBranchAddress("phoPFChIsoFrix2", &phoPFChIsoFrix2, &b_phoPFChIsoFrix2);
   fChain->SetBranchAddress("phoPFChIsoFrix3", &phoPFChIsoFrix3, &b_phoPFChIsoFrix3);
   fChain->SetBranchAddress("phoPFChIsoFrix4", &phoPFChIsoFrix4, &b_phoPFChIsoFrix4);
   fChain->SetBranchAddress("phoPFChIsoFrix5", &phoPFChIsoFrix5, &b_phoPFChIsoFrix5);
   fChain->SetBranchAddress("phoPFChIsoFrix6", &phoPFChIsoFrix6, &b_phoPFChIsoFrix6);
   fChain->SetBranchAddress("phoPFChIsoFrix7", &phoPFChIsoFrix7, &b_phoPFChIsoFrix7);
   fChain->SetBranchAddress("phoPFChIsoFrix8", &phoPFChIsoFrix8, &b_phoPFChIsoFrix8);
   fChain->SetBranchAddress("phoPFPhoIsoFrix1", &phoPFPhoIsoFrix1, &b_phoPFPhoIsoFrix1);
   fChain->SetBranchAddress("phoPFPhoIsoFrix2", &phoPFPhoIsoFrix2, &b_phoPFPhoIsoFrix2);
   fChain->SetBranchAddress("phoPFPhoIsoFrix3", &phoPFPhoIsoFrix3, &b_phoPFPhoIsoFrix3);
   fChain->SetBranchAddress("phoPFPhoIsoFrix4", &phoPFPhoIsoFrix4, &b_phoPFPhoIsoFrix4);
   fChain->SetBranchAddress("phoPFPhoIsoFrix5", &phoPFPhoIsoFrix5, &b_phoPFPhoIsoFrix5);
   fChain->SetBranchAddress("phoPFPhoIsoFrix6", &phoPFPhoIsoFrix6, &b_phoPFPhoIsoFrix6);
   fChain->SetBranchAddress("phoPFPhoIsoFrix7", &phoPFPhoIsoFrix7, &b_phoPFPhoIsoFrix7);
   fChain->SetBranchAddress("phoPFPhoIsoFrix8", &phoPFPhoIsoFrix8, &b_phoPFPhoIsoFrix8);
   fChain->SetBranchAddress("phoPFNeuIsoFrix1", &phoPFNeuIsoFrix1, &b_phoPFNeuIsoFrix1);
   fChain->SetBranchAddress("phoPFNeuIsoFrix2", &phoPFNeuIsoFrix2, &b_phoPFNeuIsoFrix2);
   fChain->SetBranchAddress("phoPFNeuIsoFrix3", &phoPFNeuIsoFrix3, &b_phoPFNeuIsoFrix3);
   fChain->SetBranchAddress("phoPFNeuIsoFrix4", &phoPFNeuIsoFrix4, &b_phoPFNeuIsoFrix4);
   fChain->SetBranchAddress("phoPFNeuIsoFrix5", &phoPFNeuIsoFrix5, &b_phoPFNeuIsoFrix5);
   fChain->SetBranchAddress("phoPFNeuIsoFrix6", &phoPFNeuIsoFrix6, &b_phoPFNeuIsoFrix6);
   fChain->SetBranchAddress("phoPFNeuIsoFrix7", &phoPFNeuIsoFrix7, &b_phoPFNeuIsoFrix7);
   fChain->SetBranchAddress("phoPFNeuIsoFrix8", &phoPFNeuIsoFrix8, &b_phoPFNeuIsoFrix8);
   fChain->SetBranchAddress("phoEcalRecHitSumEtConeDR03", &phoEcalRecHitSumEtConeDR03, &b_phoEcalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03, &b_phohcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03, &b_phohcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("phohcalTowerSumEtConeDR03", &phohcalTowerSumEtConeDR03, &b_phohcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("photrkSumPtHollowConeDR03", &photrkSumPtHollowConeDR03, &b_photrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("eleIDMVANonTrg", &eleIDMVANonTrg, &b_eleIDMVANonTrg);
   fChain->SetBranchAddress("eleIDMVATrg", &eleIDMVATrg, &b_eleIDMVATrg);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   fChain->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   fChain->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   fChain->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   fChain->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   fChain->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   fChain->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   fChain->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   fChain->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   fChain->SetBranchAddress("eleBCEn", &eleBCEn, &b_eleBCEn);
   fChain->SetBranchAddress("eleBCEta", &eleBCEta, &b_eleBCEta);
   fChain->SetBranchAddress("eleBCPhi", &eleBCPhi, &b_eleBCPhi);
   fChain->SetBranchAddress("eleBCS25", &eleBCS25, &b_eleBCS25);
   fChain->SetBranchAddress("eleBCS15", &eleBCS15, &b_eleBCS15);
   fChain->SetBranchAddress("eleBCSieie", &eleBCSieie, &b_eleBCSieie);
   fChain->SetBranchAddress("eleBCSieip", &eleBCSieip, &b_eleBCSieip);
   fChain->SetBranchAddress("eleBCSipip", &eleBCSipip, &b_eleBCSipip);
   fChain->SetBranchAddress("eleFiredTrgs", &eleFiredTrgs, &b_eleFiredTrgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleESEnP1Raw", &eleESEnP1Raw, &b_eleESEnP1Raw);
   fChain->SetBranchAddress("eleESEnP2Raw", &eleESEnP2Raw, &b_eleESEnP2Raw);
   fChain->SetBranchAddress("nGSFTrk", &nGSFTrk, &b_nGSFTrk);
   fChain->SetBranchAddress("gsfPt", &gsfPt, &b_gsfPt);
   fChain->SetBranchAddress("gsfEta", &gsfEta, &b_gsfEta);
   fChain->SetBranchAddress("gsfPhi", &gsfPhi, &b_gsfPhi);
   fChain->SetBranchAddress("npfHF", &npfHF, &b_npfHF);
   fChain->SetBranchAddress("pfHFEn", &pfHFEn, &b_pfHFEn);
   fChain->SetBranchAddress("pfHFECALEn", &pfHFECALEn, &b_pfHFECALEn);
   fChain->SetBranchAddress("pfHFHCALEn", &pfHFHCALEn, &b_pfHFHCALEn);
   fChain->SetBranchAddress("pfHFPt", &pfHFPt, &b_pfHFPt);
   fChain->SetBranchAddress("pfHFEta", &pfHFEta, &b_pfHFEta);
   fChain->SetBranchAddress("pfHFPhi", &pfHFPhi, &b_pfHFPhi);
   fChain->SetBranchAddress("pfHFIso", &pfHFIso, &b_pfHFIso);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIsLooseID", &muIsLooseID, &b_muIsLooseID);
   fChain->SetBranchAddress("muIsMediumID", &muIsMediumID, &b_muIsMediumID);
   fChain->SetBranchAddress("muIsTightID", &muIsTightID, &b_muIsTightID);
   fChain->SetBranchAddress("muIsSoftID", &muIsSoftID, &b_muIsSoftID);
   fChain->SetBranchAddress("muIsHighPtID", &muIsHighPtID, &b_muIsHighPtID);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("pfTausDiscriminationByDecayModeFinding", &pfTausDiscriminationByDecayModeFinding, &b_pfTausDiscriminationByDecayModeFinding);
   fChain->SetBranchAddress("pfTausDiscriminationByDecayModeFindingNewDMs", &pfTausDiscriminationByDecayModeFindingNewDMs, &b_pfTausDiscriminationByDecayModeFindingNewDMs);
   fChain->SetBranchAddress("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection, &b_tauByMVA5LooseElectronRejection);
   fChain->SetBranchAddress("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection, &b_tauByMVA5MediumElectronRejection);
   fChain->SetBranchAddress("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection, &b_tauByMVA5TightElectronRejection);
   fChain->SetBranchAddress("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection, &b_tauByMVA5VTightElectronRejection);
   fChain->SetBranchAddress("tauByLooseMuonRejection3", &tauByLooseMuonRejection3, &b_tauByLooseMuonRejection3);
   fChain->SetBranchAddress("tauByTightMuonRejection3", &tauByTightMuonRejection3, &b_tauByTightMuonRejection3);
   fChain->SetBranchAddress("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits, &b_tauByTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tauCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT, &b_tauByVLooseIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT, &b_tauByLooseIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT, &b_tauByMediumIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT, &b_tauByTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT, &b_tauByVTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT, &b_tauByVVTightIsolationMVA3oldDMwLT);
   fChain->SetBranchAddress("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw, &b_tauByIsolationMVA3oldDMwLTraw);
   fChain->SetBranchAddress("tauByLooseIsolationMVA3newDMwLT", &tauByLooseIsolationMVA3newDMwLT, &b_tauByLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByVLooseIsolationMVA3newDMwLT", &tauByVLooseIsolationMVA3newDMwLT, &b_tauByVLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByMediumIsolationMVA3newDMwLT", &tauByMediumIsolationMVA3newDMwLT, &b_tauByMediumIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByTightIsolationMVA3newDMwLT", &tauByTightIsolationMVA3newDMwLT, &b_tauByTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByVTightIsolationMVA3newDMwLT", &tauByVTightIsolationMVA3newDMwLT, &b_tauByVTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByVVTightIsolationMVA3newDMwLT", &tauByVVTightIsolationMVA3newDMwLT, &b_tauByVVTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("tauByIsolationMVA3newDMwLTraw", &tauByIsolationMVA3newDMwLTraw, &b_tauByIsolationMVA3newDMwLTraw);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEt", &tauEt, &b_tauEt);
   fChain->SetBranchAddress("tauCharge", &tauCharge, &b_tauCharge);
   fChain->SetBranchAddress("tauP", &tauP, &b_tauP);
   fChain->SetBranchAddress("tauPx", &tauPx, &b_tauPx);
   fChain->SetBranchAddress("tauPy", &tauPy, &b_tauPy);
   fChain->SetBranchAddress("tauPz", &tauPz, &b_tauPz);
   fChain->SetBranchAddress("tauVz", &tauVz, &b_tauVz);
   fChain->SetBranchAddress("tauEnergy", &tauEnergy, &b_tauEnergy);
   fChain->SetBranchAddress("tauMass", &tauMass, &b_tauMass);
   fChain->SetBranchAddress("tauDxy", &tauDxy, &b_tauDxy);
   fChain->SetBranchAddress("tauZImpact", &tauZImpact, &b_tauZImpact);
   fChain->SetBranchAddress("tauDecayMode", &tauDecayMode, &b_tauDecayMode);
   fChain->SetBranchAddress("tauLeadChargedHadronExists", &tauLeadChargedHadronExists, &b_tauLeadChargedHadronExists);
   fChain->SetBranchAddress("tauLeadChargedHadronEta", &tauLeadChargedHadronEta, &b_tauLeadChargedHadronEta);
   fChain->SetBranchAddress("tauLeadChargedHadronPhi", &tauLeadChargedHadronPhi, &b_tauLeadChargedHadronPhi);
   fChain->SetBranchAddress("tauLeadChargedHadronPt", &tauLeadChargedHadronPt, &b_tauLeadChargedHadronPt);
   fChain->SetBranchAddress("tauChargedIsoPtSum", &tauChargedIsoPtSum, &b_tauChargedIsoPtSum);
   fChain->SetBranchAddress("tauNeutralIsoPtSum", &tauNeutralIsoPtSum, &b_tauNeutralIsoPtSum);
   fChain->SetBranchAddress("tauPuCorrPtSum", &tauPuCorrPtSum, &b_tauPuCorrPtSum);
   fChain->SetBranchAddress("tauNumSignalPFChargedHadrCands", &tauNumSignalPFChargedHadrCands, &b_tauNumSignalPFChargedHadrCands);
   fChain->SetBranchAddress("tauNumSignalPFNeutrHadrCands", &tauNumSignalPFNeutrHadrCands, &b_tauNumSignalPFNeutrHadrCands);
   fChain->SetBranchAddress("tauNumSignalPFGammaCands", &tauNumSignalPFGammaCands, &b_tauNumSignalPFGammaCands);
   fChain->SetBranchAddress("tauNumSignalPFCands", &tauNumSignalPFCands, &b_tauNumSignalPFCands);
   fChain->SetBranchAddress("tauNumIsolationPFChargedHadrCands", &tauNumIsolationPFChargedHadrCands, &b_tauNumIsolationPFChargedHadrCands);
   fChain->SetBranchAddress("tauNumIsolationPFNeutrHadrCands", &tauNumIsolationPFNeutrHadrCands, &b_tauNumIsolationPFNeutrHadrCands);
   fChain->SetBranchAddress("tauNumIsolationPFGammaCands", &tauNumIsolationPFGammaCands, &b_tauNumIsolationPFGammaCands);
   fChain->SetBranchAddress("tauNumIsolationPFCands", &tauNumIsolationPFCands, &b_tauNumIsolationPFCands);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetpfCombinedInclusiveSecondaryVertexV2BJetTags", &jetpfCombinedInclusiveSecondaryVertexV2BJetTags, &b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetpfCombinedMVABJetTags", &jetpfCombinedMVABJetTags, &b_jetpfCombinedMVABJetTags);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetPUidFullDiscriminant", &jetPUidFullDiscriminant, &b_jetPUidFullDiscriminant);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("nAK8Jet", &nAK8Jet, &b_nAK8Jet);
   fChain->SetBranchAddress("AK8JetPt", &AK8JetPt, &b_AK8JetPt);
   fChain->SetBranchAddress("AK8JetEn", &AK8JetEn, &b_AK8JetEn);
   fChain->SetBranchAddress("AK8JetRawPt", &AK8JetRawPt, &b_AK8JetRawPt);
   fChain->SetBranchAddress("AK8JetRawEn", &AK8JetRawEn, &b_AK8JetRawEn);
   fChain->SetBranchAddress("AK8JetEta", &AK8JetEta, &b_AK8JetEta);
   fChain->SetBranchAddress("AK8JetPhi", &AK8JetPhi, &b_AK8JetPhi);
   fChain->SetBranchAddress("AK8JetMass", &AK8JetMass, &b_AK8JetMass);
   fChain->SetBranchAddress("AK8Jet_tau1", &AK8Jet_tau1, &b_AK8Jet_tau1);
   fChain->SetBranchAddress("AK8Jet_tau2", &AK8Jet_tau2, &b_AK8Jet_tau2);
   fChain->SetBranchAddress("AK8Jet_tau3", &AK8Jet_tau3, &b_AK8Jet_tau3);
   fChain->SetBranchAddress("AK8JetCHF", &AK8JetCHF, &b_AK8JetCHF);
   fChain->SetBranchAddress("AK8JetNHF", &AK8JetNHF, &b_AK8JetNHF);
   fChain->SetBranchAddress("AK8JetCEF", &AK8JetCEF, &b_AK8JetCEF);
   fChain->SetBranchAddress("AK8JetNEF", &AK8JetNEF, &b_AK8JetNEF);
   fChain->SetBranchAddress("AK8JetNCH", &AK8JetNCH, &b_AK8JetNCH);
   fChain->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
   fChain->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
   fChain->SetBranchAddress("AK8CHSSoftDropJetMass", &AK8CHSSoftDropJetMass, &b_AK8CHSSoftDropJetMass);
   fChain->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
   fChain->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
   fChain->SetBranchAddress("AK8JetPartonID", &AK8JetPartonID, &b_AK8JetPartonID);
   fChain->SetBranchAddress("AK8JetGenJetIndex", &AK8JetGenJetIndex, &b_AK8JetGenJetIndex);
   fChain->SetBranchAddress("AK8JetGenJetEn", &AK8JetGenJetEn, &b_AK8JetGenJetEn);
   fChain->SetBranchAddress("AK8JetGenJetPt", &AK8JetGenJetPt, &b_AK8JetGenJetPt);
   fChain->SetBranchAddress("AK8JetGenJetEta", &AK8JetGenJetEta, &b_AK8JetGenJetEta);
   fChain->SetBranchAddress("AK8JetGenJetPhi", &AK8JetGenJetPhi, &b_AK8JetGenJetPhi);
   fChain->SetBranchAddress("AK8JetGenPartonID", &AK8JetGenPartonID, &b_AK8JetGenPartonID);
   fChain->SetBranchAddress("AK8JetGenEn", &AK8JetGenEn, &b_AK8JetGenEn);
   fChain->SetBranchAddress("AK8JetGenPt", &AK8JetGenPt, &b_AK8JetGenPt);
   fChain->SetBranchAddress("AK8JetGenEta", &AK8JetGenEta, &b_AK8JetGenEta);
   fChain->SetBranchAddress("AK8JetGenPhi", &AK8JetGenPhi, &b_AK8JetGenPhi);
   fChain->SetBranchAddress("AK8JetGenPartonMomID", &AK8JetGenPartonMomID, &b_AK8JetGenPartonMomID);
   fChain->SetBranchAddress("nAK8softdropSubjet", &nAK8softdropSubjet, &b_nAK8softdropSubjet);
   fChain->SetBranchAddress("AK8softdropSubjetPt", &AK8softdropSubjetPt, &b_AK8softdropSubjetPt);
   fChain->SetBranchAddress("AK8softdropSubjetEta", &AK8softdropSubjetEta, &b_AK8softdropSubjetEta);
   fChain->SetBranchAddress("AK8softdropSubjetPhi", &AK8softdropSubjetPhi, &b_AK8softdropSubjetPhi);
   fChain->SetBranchAddress("AK8softdropSubjetMass", &AK8softdropSubjetMass, &b_AK8softdropSubjetMass);
   fChain->SetBranchAddress("AK8softdropSubjetE", &AK8softdropSubjetE, &b_AK8softdropSubjetE);
   fChain->SetBranchAddress("AK8softdropSubjetCharge", &AK8softdropSubjetCharge, &b_AK8softdropSubjetCharge);
   fChain->SetBranchAddress("AK8softdropSubjetFlavour", &AK8softdropSubjetFlavour, &b_AK8softdropSubjetFlavour);
   fChain->SetBranchAddress("AK8softdropSubjetCSV", &AK8softdropSubjetCSV, &b_AK8softdropSubjetCSV);
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

Bool_t PostAnalyzer_MC::GoodPrimaryVtx(Int_t &GoodVertex){

  Bool_t passVtx = false;
  GoodVertex = 0;

  for(Int_t i=0; i < nVtx; ++i){
    if( (fabs((*vtz)[i])) <= Cut_Vtx_z &&
        (*vndof)[i] >= Cut_Vtx_ndof    &&
        !((*isFake)[i])                &&
        (fabs((*vrho)[i])) <= Cut_Vtx_rho )
      GoodVertex++;
  }
  if(GoodVertex > 0) passVtx = true;

  return passVtx;

}

//Cut Based Photon ID for 25ns (https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns)
Bool_t PostAnalyzer_MC::CutBasedPhotonID(Int_t ipho, TString phoWP){

  Bool_t PhID = false;

  if(phoWP == "loose"){ //loose 
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.05)                     &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0102)     &&
        ((*phoPFChIso)[ipho] < 3.32)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 1.92 + 0.014*(*phoEt)[ipho] + 0.000019*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.81 + 0.0053*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.05)                       &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0274)       &&
      ((*phoPFChIso)[ipho] < 1.97)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 11.86 + 0.0139*(*phoEt)[ipho] + 0.000025*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.83 + 0.0034*(*phoEt)[ipho]);

    }
  }
  if(phoWP == "medium"){ //medium
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.05)                     &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0102)     &&
        ((*phoPFChIso)[ipho] < 1.37)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 1.06 + 0.014*(*phoEt)[ipho] + 0.000019*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
	 ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.28 + 0.0053*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.05)                       &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0268)       &&
      ((*phoPFChIso)[ipho] < 1.10)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 2.69 + 0.0139*(*phoEt)[ipho] + 0.000025*(*phoEt)[ipho]*(*phoEt)[ipho])                                           &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.39 + 0.0034*(*phoEt)[ipho]);

    }
  }        
  if(phoWP == "tight"){ //tight
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
	((*phoHoverE)[ipho] < 0.05)                     &&
	((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0100)     &&
        ((*phoPFChIso)[ipho] < 0.76)                    &&
	((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 0.97 + 0.014*(*phoEt)[ipho] + 0.000019*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
	 ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.08 + 0.0053*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.05)                       &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0268)       &&
      ((*phoPFChIso)[ipho] < 0.56)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 2.09 + 0.0139*(*phoEt)[ipho] + 0.000025*(*phoEt)[ipho]*(*phoEt)[ipho])                                           &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 0.16 + 0.0034*(*phoEt)[ipho]);

    }
  }
  return PhID;
}

//(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_det_AN1)
Double_t PostAnalyzer_MC::EANeutralHadrons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.0599;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.0819;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0696;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0360;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.0360;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.0462;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.0656;

  return EffArea;

}

Double_t PostAnalyzer_MC::EAPhotons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.1271;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1101;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0756;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.1175;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1498;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.1857;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.2183;

  return EffArea;

}

Int_t PostAnalyzer_MC::FirstGoodPhoton(TString phoWP){

  Int_t pc = -1;
  Bool_t ID = false;
  for(Int_t i = 0; i < nPho; i++){
    ID = CutBasedPhotonID(i, phoWP);
    if(ID){
      pc = i;
      break;
    }
  }
  return pc;
}

vector<Int_t> PostAnalyzer_MC::GoodPhotons(TString phoWP){

  vector<Int_t> goodphs;
  goodphs.clear();

  for(Int_t i = 0; i < nPho; i++){
    if(CutBasedPhotonID(i, phoWP)){
      goodphs.push_back(i);
    }
  }
  return goodphs;
}
//Taken from https://indico.cern.ch/event/455258/contribution/1/attachments/1173429/1695369/151020_Egamma_PhotonID.pdf
//https://indico.cern.ch/event/460283/contribution/3/attachments/1182205/1716953/preapproval_diphotons_13TeV.pdf (Values from this diphoton preapp talk)
//Link of AN-http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2015/241
Bool_t PostAnalyzer_MC::HighPtPhotonID(Int_t ipho){

  Double_t kappa = 0.002; // in units of (1/GeV)
  Double_t alpha = 0.0; //in units of GeV
  Double_t area = 0.0;

  Bool_t PhID = false;

  if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel
    alpha = 1.5;
    area = EAPhotons_HighPtID((*phoSCEta)[ipho]);

    PhID = ((*phoHoverE)[ipho] < 0.05)              &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0105)   &&
      ((*phoPFChIso)[ipho] < 5)                     &&
      ((alpha + (*phoPFPhoIso)[ipho] - rho*area - kappa*(*phoEt)[ipho]) < 2.75);

  }
  if((fabs((*phoSCEta)[ipho])) > 1.566){ //Endcap
    alpha = 2.0;
    area = EAPhotons_HighPtID((*phoSCEta)[ipho]);

    PhID = ((*phoHoverE)[ipho] < 0.05)              &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.028)    &&
      ((*phoPFChIso)[ipho] < 5)                     &&
      ((alpha + (*phoPFPhoIso)[ipho] - rho*area - kappa*(*phoEt)[ipho]) < 2.0);

  }

  return PhID;
}

Int_t PostAnalyzer_MC::FirstHighPtIDPhoton(){

  Int_t pc = -1;
  Bool_t ID = false;
  for(Int_t i = 0; i < nPho; i++){
    ID = HighPtPhotonID(i);
    if(ID){
      pc = i;
      break;
    }
  }
  return pc;
}

vector<Int_t> PostAnalyzer_MC::HighPtIDPhotons(){

  vector<Int_t> highptphs;
  highptphs.clear();

  for(Int_t i = 0; i < nPho; i++){
    if(HighPtPhotonID(i)){
      highptphs.push_back(i);
    }
  }
  return highptphs;
}



Double_t PostAnalyzer_MC::EAPhotons_HighPtID(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 0.9    ) EffArea = 0.21;
  if( fabs(eta) >= 0.9   && fabs(eta) < 1.4442 ) EffArea = 0.2;
  if( fabs(eta) >= 1.566 && fabs(eta) < 2.0    ) EffArea = 0.14;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2    ) EffArea = 0.22;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.5    ) EffArea = 0.33;

  return EffArea;

}

//Recommended Tight JetID for 13 TeV (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data)
Bool_t PostAnalyzer_MC::TightJetId(Int_t iJet){

  Bool_t JetID = false;

  if(fabs((*jetEta)[iJet]) <= 2.4){

    JetID = ((*jetNHF)[iJet] < 0.90)  &&
      ((*jetNEF)[iJet] < 0.90)        &&
      ((*jetNConstituents)[iJet] > 1) &&
      ((*jetCHF)[iJet] > 0)           &&
      ((*jetNCH)[iJet] > 0)           &&
      ((*jetCEF)[iJet] < 0.99);

  }
  if(fabs((*jetEta)[iJet]) > 2.4 && fabs((*jetEta)[iJet]) <= 3.0){ 

    JetID = ((*jetNHF)[iJet] < 0.90)  &&
      ((*jetNEF)[iJet] < 0.90)        &&
      ((*jetNConstituents)[iJet] > 1);

  }
  return JetID;
}

Int_t PostAnalyzer_MC::FirstGoodJet(Int_t pc){

  Int_t jc = -1;
  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    Double_t dR = -1.0;
    if(pc > -1) dR = GetdR((*phoSCEta)[pc], (*jetEta)[i], (*phoSCPhi)[pc], (*jetPhi)[i]);
    ID = TightJetId(i);
    if(dR > 0.5 && ID){
      jc = i;
      break;
    }
  }
  return jc;
}

vector<Int_t> PostAnalyzer_MC::GoodJets(Int_t pc){

  vector<Int_t> goodjets;
  goodjets.clear();

  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    Double_t dR = -1.0;
    if(pc > -1) dR = GetdR((*phoSCEta)[pc], (*jetEta)[i], (*phoSCPhi)[pc], (*jetPhi)[i]);
    ID = TightJetId(i);
    if(dR > 0.5 && ID){
      goodjets.push_back(i);
    }
  }
  return goodjets;
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

Double_t PostAnalyzer_MC::GetInvtMass(Int_t pho, Int_t jet){

  Double_t mass = 0.0;

  TLorentzVector Pho;
  Pho.SetPtEtaPhiE((*phoEt)[pho], (*phoSCEta)[pho], (*phoSCPhi)[pho], (*phoE)[pho]);

  TLorentzVector Jet;
  Jet.SetPtEtaPhiE((*jetPt)[jet], (*jetEta)[jet], (*jetPhi)[jet], (*jetEn)[jet] );

  mass = (Pho+Jet).M();

  return mass;
}


Int_t PostAnalyzer_MC::MatchedGenPhotonToReco(Int_t pc){

  Int_t pc_gen = -1;

  for(Int_t i = 0; i < nMC; i++){

    Double_t dR = 99.0;
    if((*mcStatus)[i] == 1 && (*mcPID)[i] == 22){ 

      dR = GetdR((*phoSCEta)[pc], (*mcEta)[i], (*phoSCPhi)[pc], (*mcPhi)[i]);

      if(dR < 0.1){
	pc_gen = i;
	break;
      }
    }
  }

  return pc_gen;
}


Int_t PostAnalyzer_MC::MatchedGenJetToReco(Int_t jc){

  Int_t jc_gen = -1;

  if((*jetGenPartonID)[jc] == 21 || fabs((*jetGenPartonID)[jc]) == 1 || fabs((*jetGenPartonID)[jc]) == 2) jc_gen = jc;
 
  return jc_gen;
}


Double_t PostAnalyzer_MC::GetGenLevelInvtMass(Int_t pc_gen, Int_t jc_gen){

  Double_t mass = 0.0;

  TLorentzVector GenPho;
  GenPho.SetPtEtaPhiE((*mcPt)[pc_gen], (*mcEta)[pc_gen], (*mcPhi)[pc_gen], (*mcE)[pc_gen]);

  TLorentzVector GenJet;
  GenJet.SetPtEtaPhiE((*jetGenJetPt)[jc_gen], (*jetGenJetEta)[jc_gen], (*jetGenJetPhi)[jc_gen], (*jetGenJetEn)[jc_gen]);

  mass = (GenPho+GenJet).M();
  
  return mass;

}


Int_t PostAnalyzer_MC::MatchedRecoPhotonToGen_WithGenIsoCut(){

  Int_t pc = -1;

  for(Int_t ipho = 0; ipho < nPho; ipho++){
    Double_t dR = 99.0;
    for(Int_t i = 0; i < nMC; i++){
      if((*mcPID)[i] == 22 && (*mcCalIsoDR03)[i] < 10){

	dR = GetdR((*phoSCEta)[ipho], (*mcEta)[i], (*phoSCPhi)[ipho], (*mcPhi)[i]);
	if(dR < 0.1){
	  pc = ipho;
	  break;
	}
      }
    }
    if(pc > -1) break;
  }
  return pc;
}

Double_t PostAnalyzer_MC::MCFM_kFact(Double_t phpt){

  TFile *f = TFile::Open("/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/PostAnalyzer_MC/MCFM_kFactor/MCFM_kFact_lordvsVirt_fixedscale.root");
  TH1F *h = (TH1F*)f->Get("h_MCFMnlo_over_MCFMlo_DiffXS_vs_PhotPt");
  Int_t bin = h->GetXaxis()->FindBin(phpt);
  return h->GetBinContent( bin );

}
  //For pileup distribution of MC, we need distribution for trueNumofInteractions. In ggNtuplizer, this distribution has been already saved by the name of histogram 'hPUTrue'. So i have taken one sample of each of MC (GJetsHT100to120, QCD_Pt_300to470, Qstar_M1000_f1p0) and hadded all files and got the 'hPUTrue' histogram from those and saved in PileupHistograms/MC_Run2015D_PileUpDist folder. The root script PU.cc used for this getting the hist is also present in the same folder.
void PostAnalyzer_MC::PileupReWeighting(){

  //uncomment in script
  //    Bool_t Pileup_QCD   = ${QCD};
  //    Bool_t Pileup_GJ    = ${GJ};
  //    Bool_t Pileup_EWK   = ${EWK};
  //    Bool_t Pileup_Qstar = ${Qstar};

    TString puMCfile;
    //uncomment in script
    //    if(Pileup_Qstar) puMCfile = "Qstar_PileupHist";
    //    if(Pileup_GJ)    puMCfile = "GJets_PileupHist";
    //    if(Pileup_EWK)   puMCfile = "GJets_PileupHist";
    //    if(Pileup_QCD)   puMCfile = "DiJet_PileupHist";
    
    //remove this in script
    puMCfile = "DiJet_PileupHist";

    //For Golden Json
    //  TFile *fData = TFile::Open("/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/PostAnalyzer_MC/PileupHistograms/Data_Run2015D_PileUpDist_GoldenJson/DataPileupHist_2015D.root");
    //For Silver Json
  TFile *fData = TFile::Open("/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/PostAnalyzer_MC/PileupHistograms/Data_Run2015D_PileUpDist_SilverJson/DataPileupHist_2015D.root");
  TH1F *dataPU = (TH1F*)fData->Get("pileup");
   
  TFile *fMC = TFile::Open("/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/PostAnalyzer_MC/PileupHistograms/MC_Run2015D_PileUpDist/"+puMCfile+".root");
  TH1F *mcPU = (TH1F*)fMC->Get("hPUTrue");

  mcPU->Rebin(5); //because "hPUTrue" has binning of width 0.2 (1000, 0, 200), that is 1000 bins in the range of 0 to 200.

  std::vector<float> DataPileUp;
  std::vector<float> mcPileUp;
  DataPileUp.clear();
  mcPileUp.clear();
  for(Int_t i = 0; i < 50; i++){
    DataPileUp.push_back(dataPU->GetBinContent(i+1));
    mcPileUp.push_back(mcPU->GetBinContent(i+1));
  }

  TH1F *h_MCWeights = new TH1F("h_MCWeights", "MC PileUp Weights", 50, 0, 50);
  for(Int_t ibin = 0; ibin < 50; ibin++){
    h_DataPUNormDist->SetBinContent(ibin+1, DataPileUp[ibin]); //This to get in output
    h_PUScaleFactor->SetBinContent(ibin+1, DataPileUp[ibin]);
    h_MCPUNormDist->SetBinContent(ibin+1, mcPileUp[ibin]);
    h_MCWeights->SetBinContent(ibin+1, mcPileUp[ibin]);
  }

  h_DataPUNormDist->Scale(1.0/h_DataPUNormDist->Integral());
  h_PUScaleFactor->Scale(1.0/h_PUScaleFactor->Integral());
  h_MCPUNormDist->Scale(1.0/h_MCPUNormDist->Integral());
  h_MCWeights->Scale(1.0/h_MCWeights->Integral());

  h_PUScaleFactor->Divide(h_MCWeights);

}

Double_t PostAnalyzer_MC::PUWeights(Float_t npv){
  Int_t bin = h_PUScaleFactor->GetXaxis()->FindBin( npv );
  return h_PUScaleFactor->GetBinContent( bin );
}

void PostAnalyzer_MC::BookHistograms(){
  file->cd();

    char name[100];
    const Int_t nMassBins = 119;
    const Double_t MassBin[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

    std::string cut[2] = {"noMassCut", "MassCut"};
    for( Int_t hist = 0; hist < 2; ++hist ){ 
      sprintf(name, "h_PhotonPt_%s",cut[hist].c_str());
      h_PhotonPt[hist] = new TH1F(name,"Pt distribution of photons",150,0.0,6000.0);
      h_PhotonPt[hist]->GetYaxis()->SetTitle("Events/40 GeV");           h_PhotonPt[hist]->GetYaxis()->CenterTitle();                                
      h_PhotonPt[hist]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_PhotonPt[hist]->GetXaxis()->CenterTitle();                                
      h_PhotonPt[hist]->Sumw2();                                             

      sprintf(name, "h_PhotonEta_%s",cut[hist].c_str());
      h_PhotonEta[hist] = new TH1F(name,"Eta distribution of photons",100,-2.5,2.5);
      h_PhotonEta[hist]->GetYaxis()->SetTitle("Events");           h_PhotonEta[hist]->GetYaxis()->CenterTitle();                               
      h_PhotonEta[hist]->GetXaxis()->SetTitle("#eta^{#gamma}");    h_PhotonEta[hist]->GetXaxis()->CenterTitle();                               
      h_PhotonEta[hist]->Sumw2();                                             

      sprintf(name, "h_PhotonPhi_%s",cut[hist].c_str());
      h_PhotonPhi[hist] = new TH1F(name,"Phi distribution of photons",100,-3.2,3.2);
      h_PhotonPhi[hist]->GetYaxis()->SetTitle("Events");           h_PhotonPhi[hist]->GetYaxis()->CenterTitle();
      h_PhotonPhi[hist]->GetXaxis()->SetTitle("#phi^{#gamma}");    h_PhotonPhi[hist]->GetXaxis()->CenterTitle();
      h_PhotonPhi[hist]->Sumw2();

      sprintf(name, "h_Photon_SigmaIEtaIEta_%s",cut[hist].c_str());
      h_Photon_SigmaIEtaIEta[hist] = new TH1F(name,"Photon SigmaIetaIeta Distribution",100,0.0,0.05);
      h_Photon_SigmaIEtaIEta[hist]->GetYaxis()->SetTitle("Events");                 h_Photon_SigmaIEtaIEta[hist]->GetYaxis()->CenterTitle();
      h_Photon_SigmaIEtaIEta[hist]->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");    h_Photon_SigmaIEtaIEta[hist]->GetXaxis()->CenterTitle();
      h_Photon_SigmaIEtaIEta[hist]->Sumw2();

      sprintf(name, "h_Photon_R9_%s",cut[hist].c_str());
      h_Photon_R9[hist] = new TH1F(name,"Photon R9 Distribution",100,0.0,10.0);
      h_Photon_R9[hist]->GetYaxis()->SetTitle("Events");           h_Photon_R9[hist]->GetYaxis()->CenterTitle();
      h_Photon_R9[hist]->GetXaxis()->SetTitle("Photon_r9");        h_Photon_R9[hist]->GetXaxis()->CenterTitle();
      h_Photon_R9[hist]->Sumw2();

      sprintf(name, "h_Photon_HoverE_%s",cut[hist].c_str());
      h_Photon_HoverE[hist] = new TH1F(name,"Photon HoverE Distribution",50,0.0,0.1);
      h_Photon_HoverE[hist]->GetYaxis()->SetTitle("Events");           h_Photon_HoverE[hist]->GetYaxis()->CenterTitle();
      h_Photon_HoverE[hist]->GetXaxis()->SetTitle("H/E");              h_Photon_HoverE[hist]->GetXaxis()->CenterTitle();
      h_Photon_HoverE[hist]->Sumw2();

      sprintf(name, "h_Photon_EleVeto_%s",cut[hist].c_str());
      h_Photon_EleVeto[hist] = new TH1F(name,"Photon ElectronVeto",3,0,3);
      h_Photon_EleVeto[hist]->GetYaxis()->SetTitle("Events");           h_Photon_EleVeto[hist]->GetYaxis()->CenterTitle();
      h_Photon_EleVeto[hist]->GetXaxis()->SetTitle("ElectronVeto");    h_Photon_EleVeto[hist]->GetXaxis()->CenterTitle();
      h_Photon_EleVeto[hist]->Sumw2();

      sprintf(name, "h_Photon_CorrPFChIso_%s",cut[hist].c_str());
      h_Photon_CorrPFChIso[hist] = new TH1F(name,"Rho Corrected PF Charged Isolation",625,0,25);
      h_Photon_CorrPFChIso[hist]->GetYaxis()->SetTitle("Events");               h_Photon_CorrPFChIso[hist]->GetYaxis()->CenterTitle();
      h_Photon_CorrPFChIso[hist]->GetXaxis()->SetTitle("CorrPFChIso");          h_Photon_CorrPFChIso[hist]->GetXaxis()->CenterTitle();
      h_Photon_CorrPFChIso[hist]->Sumw2();

      sprintf(name, "h_Photon_CorrPFNeuIso_%s",cut[hist].c_str());
      h_Photon_CorrPFNeuIso[hist] = new TH1F(name,"Rho Corrected PF Neutral Isolation",100,0,300);
      h_Photon_CorrPFNeuIso[hist]->GetYaxis()->SetTitle("Events");               h_Photon_CorrPFNeuIso[hist]->GetYaxis()->CenterTitle();
      h_Photon_CorrPFNeuIso[hist]->GetXaxis()->SetTitle("CorrPFNeuIso");         h_Photon_CorrPFNeuIso[hist]->GetXaxis()->CenterTitle();
      h_Photon_CorrPFNeuIso[hist]->Sumw2();

      sprintf(name, "h_Photon_CorrPFPhoIso_%s",cut[hist].c_str());
      h_Photon_CorrPFPhoIso[hist] = new TH1F(name,"Rho Corrected PF Photon Isolation",125,0,25);
      h_Photon_CorrPFPhoIso[hist]->GetYaxis()->SetTitle("Events");           h_Photon_CorrPFPhoIso[hist]->GetYaxis()->CenterTitle();
      h_Photon_CorrPFPhoIso[hist]->GetXaxis()->SetTitle("CorrPFPhoIso");    h_Photon_CorrPFPhoIso[hist]->GetXaxis()->CenterTitle();
      h_Photon_CorrPFPhoIso[hist]->Sumw2();


      sprintf(name, "h_JetPt_%s",cut[hist].c_str());
      h_JetPt[hist] = new TH1F(name,"Pt distribution of jets",150,0.0,6000.0);
      h_JetPt[hist]->GetYaxis()->SetTitle("Events/40 GeV");           h_JetPt[hist]->GetYaxis()->CenterTitle();                                
      h_JetPt[hist]->GetXaxis()->SetTitle("P_{T}^{Jet} (GeV)");       h_JetPt[hist]->GetXaxis()->CenterTitle();                                
      h_JetPt[hist]->Sumw2();                                             

      sprintf(name, "h_JetEta_%s",cut[hist].c_str());
      h_JetEta[hist] = new TH1F(name,"Eta distribution of jets",120,-3.0,3.0);
      h_JetEta[hist]->GetYaxis()->SetTitle("Events");           h_JetEta[hist]->GetYaxis()->CenterTitle();                               
      h_JetEta[hist]->GetXaxis()->SetTitle("#eta^{Jet}");       h_JetEta[hist]->GetXaxis()->CenterTitle();                               
      h_JetEta[hist]->Sumw2();                                             

      sprintf(name, "h_JetPhi_%s",cut[hist].c_str());
      h_JetPhi[hist] = new TH1F(name,"Phi distribution of jets",100,-3.2,3.2);
      h_JetPhi[hist]->GetYaxis()->SetTitle("Events");           h_JetPhi[hist]->GetYaxis()->CenterTitle();
      h_JetPhi[hist]->GetXaxis()->SetTitle("#phi^{Jet}");       h_JetPhi[hist]->GetXaxis()->CenterTitle();
      h_JetPhi[hist]->Sumw2();

      sprintf(name, "h_Jet_NHEF_%s",cut[hist].c_str());
      h_Jet_NHEF[hist] = new TH1F(name,"Neutral Hadron Energy Fraction",25,0,1);
      h_Jet_NHEF[hist]->GetYaxis()->SetTitle("Events");           h_Jet_NHEF[hist]->GetYaxis()->CenterTitle();
      h_Jet_NHEF[hist]->GetXaxis()->SetTitle("Jet_NHEF");         h_Jet_NHEF[hist]->GetXaxis()->CenterTitle();
      h_Jet_NHEF[hist]->Sumw2();

      sprintf(name, "h_Jet_NEEF_%s",cut[hist].c_str());
      h_Jet_NEEF[hist] = new TH1F(name,"Neutral Em Energy Fraction",25,0,1);
      h_Jet_NEEF[hist]->GetYaxis()->SetTitle("Events");           h_Jet_NEEF[hist]->GetYaxis()->CenterTitle();
      h_Jet_NEEF[hist]->GetXaxis()->SetTitle("Jet_NEEF");         h_Jet_NEEF[hist]->GetXaxis()->CenterTitle();
      h_Jet_NEEF[hist]->Sumw2();

      sprintf(name, "h_Jet_NConst_%s",cut[hist].c_str());
      h_Jet_NConst[hist] = new TH1F(name,"No. of Constituents",50,0,50);
      h_Jet_NConst[hist]->GetYaxis()->SetTitle("Events");           h_Jet_NConst[hist]->GetYaxis()->CenterTitle();
      h_Jet_NConst[hist]->GetXaxis()->SetTitle("Jet_NConst");       h_Jet_NConst[hist]->GetXaxis()->CenterTitle();
      h_Jet_NConst[hist]->Sumw2();

      sprintf(name, "h_Jet_CHEF_%s",cut[hist].c_str());
      h_Jet_CHEF[hist] = new TH1F(name,"Charged Hadron Energy Fraction",25,0,1);
      h_Jet_CHEF[hist]->GetYaxis()->SetTitle("Events");           h_Jet_CHEF[hist]->GetYaxis()->CenterTitle();
      h_Jet_CHEF[hist]->GetXaxis()->SetTitle("Jet_CHEF");         h_Jet_CHEF[hist]->GetXaxis()->CenterTitle();
      h_Jet_CHEF[hist]->Sumw2();

      sprintf(name, "h_Jet_ChMult_%s",cut[hist].c_str());
      h_Jet_ChMult[hist] = new TH1F(name,"Charged Multiplicity",50,0,50);
      h_Jet_ChMult[hist]->GetYaxis()->SetTitle("Events");           h_Jet_ChMult[hist]->GetYaxis()->CenterTitle();
      h_Jet_ChMult[hist]->GetXaxis()->SetTitle("Jet_ChMult");       h_Jet_ChMult[hist]->GetXaxis()->CenterTitle();
      h_Jet_ChMult[hist]->Sumw2();

      sprintf(name, "h_Jet_CEEF_%s",cut[hist].c_str());
      h_Jet_CEEF[hist] = new TH1F(name,"Charged Em Energy Fraction",25,0,1);
      h_Jet_CEEF[hist]->GetYaxis()->SetTitle("Events");           h_Jet_CEEF[hist]->GetYaxis()->CenterTitle();
      h_Jet_CEEF[hist]->GetXaxis()->SetTitle("Jet_CEEF");         h_Jet_CEEF[hist]->GetXaxis()->CenterTitle();
      h_Jet_CEEF[hist]->Sumw2();


      sprintf(name, "h_GJetInvtMass_bin40_%s",cut[hist].c_str());
      h_GJetInvtMass_bin40[hist] = new TH1F(name,"Invt mass of photon+jet", 350, 0.0, 14000.0);
      h_GJetInvtMass_bin40[hist]->GetYaxis()->SetTitle("Events/40 GeV");       h_GJetInvtMass_bin40[hist]->GetYaxis()->CenterTitle();
      h_GJetInvtMass_bin40[hist]->GetXaxis()->SetTitle("Invt Mass M_{#gamma jet} (GeV)");  h_GJetInvtMass_bin40[hist]->GetXaxis()->CenterTitle();
      h_GJetInvtMass_bin40[hist]->Sumw2(); 
     
      sprintf(name, "h_GJetInvtMass_VarBin_%s",cut[hist].c_str());
      h_GJetInvtMass_VarBin[hist] = new TH1F(name, "Invt mass of photon+jet", nMassBins, MassBin);
      h_GJetInvtMass_VarBin[hist]->GetYaxis()->SetTitle("Events/VarBin");       h_GJetInvtMass_VarBin[hist]->GetYaxis()->CenterTitle();
      h_GJetInvtMass_VarBin[hist]->GetXaxis()->SetTitle("Invt Mass M_{#gamma jet} (GeV)");  h_GJetInvtMass_VarBin[hist]->GetXaxis()->CenterTitle();
      h_GJetInvtMass_VarBin[hist]->Sumw2();

      sprintf(name, "h_GJetInvtMass_UnitBin_%s",cut[hist].c_str());
      h_GJetInvtMass_UnitBin[hist] = new TH1F(name,"Invt mass of photon+jet", 14000, 0.0, 14000.0);
      h_GJetInvtMass_UnitBin[hist]->GetYaxis()->SetTitle("Events/UnitBin");       h_GJetInvtMass_UnitBin[hist]->GetYaxis()->CenterTitle();
      h_GJetInvtMass_UnitBin[hist]->GetXaxis()->SetTitle("Invt Mass M_{#gamma jet} (GeV)");  h_GJetInvtMass_UnitBin[hist]->GetXaxis()->CenterTitle();
      h_GJetInvtMass_UnitBin[hist]->Sumw2(); 
     
      sprintf(name, "h_GJet_dEta_%s",cut[hist].c_str());
      h_GJet_dEta[hist] = new TH1F(name,"dEta between photon and jet",120,0,6);
      h_GJet_dEta[hist]->GetYaxis()->SetTitle("Events");           h_GJet_dEta[hist]->GetYaxis()->CenterTitle();
      h_GJet_dEta[hist]->GetXaxis()->SetTitle("#Delta #eta");      h_GJet_dEta[hist]->GetXaxis()->CenterTitle();
      h_GJet_dEta[hist]->Sumw2();

      sprintf(name, "h_GJet_dPhi_%s",cut[hist].c_str());
      h_GJet_dPhi[hist] = new TH1F(name,"dPhi between photon and jet",64,0,3.2);
      h_GJet_dPhi[hist]->GetYaxis()->SetTitle("Events");           h_GJet_dPhi[hist]->GetYaxis()->CenterTitle();
      h_GJet_dPhi[hist]->GetXaxis()->SetTitle("#Delta #phi");      h_GJet_dPhi[hist]->GetXaxis()->CenterTitle();
      h_GJet_dPhi[hist]->Sumw2();

      sprintf(name, "h_GJet_dR_%s",cut[hist].c_str());
      h_GJet_dR[hist] = new TH1F(name,"dR between photon and jet",100,0.0,10.0);
      h_GJet_dR[hist]->GetYaxis()->SetTitle("Events");           h_GJet_dR[hist]->GetYaxis()->CenterTitle();
      h_GJet_dR[hist]->GetXaxis()->SetTitle("#Delta R");         h_GJet_dR[hist]->GetXaxis()->CenterTitle();
      h_GJet_dR[hist]->Sumw2();

      sprintf(name, "h_PFMet_%s",cut[hist].c_str());
      h_PFMet[hist] = new TH1F(name, "PFMet distribution", 150,0.0,6000.0);
      h_PFMet[hist]->GetYaxis()->SetTitle("Events");      h_PFMet[hist]->GetYaxis()->CenterTitle();
      h_PFMet[hist]->GetXaxis()->SetTitle("PFMET (GeV)");   h_PFMet[hist]->GetXaxis()->CenterTitle();
      h_PFMet[hist]->Sumw2();

      sprintf(name, "h_PFMetVsGJmass_%s",cut[hist].c_str());
      h_PFMetVsGJmass[hist] = new TH2F(name, "PFMet Vs GJInvtMass", 1400, 0.0, 14000.0, 150, 0.0, 6000.0);
      h_PFMetVsGJmass[hist]->GetYaxis()->SetTitle("PFMET (GeV)");      h_PFMetVsGJmass[hist]->GetYaxis()->CenterTitle();
      h_PFMetVsGJmass[hist]->GetXaxis()->SetTitle("Invt Mass M_{#gamma jet} (GeV)");   h_PFMetVsGJmass[hist]->GetXaxis()->CenterTitle();
      h_PFMetVsGJmass[hist]->Sumw2();

      sprintf(name, "h_PFMetOverSumEtVsGJmass_%s",cut[hist].c_str());
      h_PFMetOverSumEtVsGJmass[hist] = new TH2F(name, "PFMet/SumEt Vs GJInvtMass", 1400, 0.0, 14000.0, 1000, 0.0, 10.0);
      h_PFMetOverSumEtVsGJmass[hist]->GetYaxis()->SetTitle("PFMET/SumEt (GeV)");      h_PFMetOverSumEtVsGJmass[hist]->GetYaxis()->CenterTitle();
      h_PFMetOverSumEtVsGJmass[hist]->GetXaxis()->SetTitle("Invt Mass M_{#gamma jet} (GeV)");  h_PFMetOverSumEtVsGJmass[hist]->GetXaxis()->CenterTitle();
      h_PFMetOverSumEtVsGJmass[hist]->Sumw2();

      sprintf(name, "h_Gen_GJetInvtMass_bin40_%s",cut[hist].c_str());
      h_Gen_GJetInvtMass_bin40[hist] = new TH1F(name,"GenLevel Invt mass of photon+jet", 350, 0.0, 14000.0);
      h_Gen_GJetInvtMass_bin40[hist]->GetYaxis()->SetTitle("Events/40 GeV");       h_Gen_GJetInvtMass_bin40[hist]->GetYaxis()->CenterTitle();
      h_Gen_GJetInvtMass_bin40[hist]->GetXaxis()->SetTitle("GenLevel Invt Mass M_{#gamma jet} (GeV)");  h_Gen_GJetInvtMass_bin40[hist]->GetXaxis()->CenterTitle();
      h_Gen_GJetInvtMass_bin40[hist]->Sumw2(); 
     
      sprintf(name, "h_Gen_GJetInvtMass_VarBin_%s",cut[hist].c_str());
      h_Gen_GJetInvtMass_VarBin[hist] = new TH1F(name, "GenLevel Invt mass of photon+jet", nMassBins, MassBin);
      h_Gen_GJetInvtMass_VarBin[hist]->GetYaxis()->SetTitle("Events/VarBin");       h_Gen_GJetInvtMass_VarBin[hist]->GetYaxis()->CenterTitle();
      h_Gen_GJetInvtMass_VarBin[hist]->GetXaxis()->SetTitle("GenLevel Invt Mass M_{#gamma jet} (GeV)");  h_Gen_GJetInvtMass_VarBin[hist]->GetXaxis()->CenterTitle();
      h_Gen_GJetInvtMass_VarBin[hist]->Sumw2();

      sprintf(name, "h_Gen_GJetInvtMass_UnitBin_%s",cut[hist].c_str());
      h_Gen_GJetInvtMass_UnitBin[hist] = new TH1F(name,"GenLevel Invt mass of photon+jet", 14000, 0.0, 14000.0);
      h_Gen_GJetInvtMass_UnitBin[hist]->GetYaxis()->SetTitle("Events/UnitBin");       h_Gen_GJetInvtMass_UnitBin[hist]->GetYaxis()->CenterTitle();
      h_Gen_GJetInvtMass_UnitBin[hist]->GetXaxis()->SetTitle("GenLevel Invt Mass M_{#gamma jet} (GeV)");  h_Gen_GJetInvtMass_UnitBin[hist]->GetXaxis()->CenterTitle();
      h_GJetInvtMass_UnitBin[hist]->Sumw2(); 

      sprintf(name, "h_GJetMassResolution_%s",cut[hist].c_str());
      h_GJetMassResolution[hist] = new TH1F(name,"Mass Resolution of photon+jet", 100, 0.0, 1.0);
      h_GJetMassResolution[hist]->GetYaxis()->SetTitle("Events");  h_GJetMassResolution[hist]->GetYaxis()->CenterTitle();
      h_GJetMassResolution[hist]->GetXaxis()->SetTitle("Mass Resolution");  h_GJetMassResolution[hist]->GetXaxis()->CenterTitle();
      h_GJetMassResolution[hist]->Sumw2();


    }

    std::string cut1[3] = {"noCut", "noMassCut", "MassCut"};
    for( Int_t hist = 0; hist < 3; ++hist ){
      sprintf(name, "h_goodPV_TotalWt_%s",cut1[hist].c_str());
      h_goodPV_TotalWt[hist] = new TH1F(name,"No. of Good Primary Vertices", 50, 0, 50);
      h_goodPV_TotalWt[hist]->GetYaxis()->SetTitle("Events");   h_goodPV_TotalWt[hist]->GetYaxis()->CenterTitle();
      h_goodPV_TotalWt[hist]->GetXaxis()->SetTitle("nPV");      h_goodPV_TotalWt[hist]->GetXaxis()->CenterTitle();
      h_goodPV_TotalWt[hist]->Sumw2();

      sprintf(name, "h_nPhotons_%s",cut1[hist].c_str());
      h_nPhotons[hist] = new TH1F(name,"No. of Photons", 400, 0, 200);
      h_nPhotons[hist]->GetYaxis()->SetTitle("Events");      h_nPhotons[hist]->GetYaxis()->CenterTitle();
      h_nPhotons[hist]->GetXaxis()->SetTitle("nPhotons");    h_nPhotons[hist]->GetXaxis()->CenterTitle();
      h_nPhotons[hist]->Sumw2();

      sprintf(name, "h_nJets_%s",cut1[hist].c_str());
      h_nJets[hist] = new TH1F(name,"No. of Jets", 400, 0, 200);
      h_nJets[hist]->GetYaxis()->SetTitle("Events");      h_nJets[hist]->GetYaxis()->CenterTitle();
      h_nJets[hist]->GetXaxis()->SetTitle("nJets");       h_nJets[hist]->GetXaxis()->CenterTitle();
      h_nJets[hist]->Sumw2();


    }

    h_mass_X_bin1 = new TH1F("h_mass_X_bin1","mass of photon+jet",14000,0.0,5.0); 
    h_mass_X_bin1->GetYaxis()->SetTitle("Events");                            h_mass_X_bin1->GetYaxis()->CenterTitle();
    h_mass_X_bin1->GetXaxis()->SetTitle("X = M_{#gamma jet}/M_{q^{#ast}}");   h_mass_X_bin1->GetXaxis()->CenterTitle();
    h_mass_X_bin1->Sumw2();

    //Pileup
    h_trueNumofInt = new TH1F("h_trueNumofInt","True Num of Interactions", 50, 0, 50);
    h_trueNumofInt->GetYaxis()->SetTitle("Events");   h_trueNumofInt->GetYaxis()->CenterTitle();
    h_trueNumofInt->GetXaxis()->SetTitle("nPV");      h_trueNumofInt->GetXaxis()->CenterTitle();
    h_trueNumofInt->Sumw2();

    h_goodPV = new TH1F("h_goodPV","No. of Good Primary Vertices", 50, 0, 50);
    h_goodPV->GetYaxis()->SetTitle("Events");   h_goodPV->GetYaxis()->CenterTitle();
    h_goodPV->GetXaxis()->SetTitle("nPV");      h_goodPV->GetXaxis()->CenterTitle();
    h_goodPV->Sumw2();

    h_goodPV_LumiWt = new TH1F("h_goodPV_LumiWt","No. of Good Primary Vertices", 50, 0, 50);
    h_goodPV_LumiWt->GetYaxis()->SetTitle("Events");   h_goodPV_LumiWt->GetYaxis()->CenterTitle();
    h_goodPV_LumiWt->GetXaxis()->SetTitle("nPV");      h_goodPV_LumiWt->GetXaxis()->CenterTitle();
    h_goodPV_LumiWt->Sumw2();

    h_goodPV_PUWt = new TH1F("h_goodPV_PUWt","No. of Good Primary Vertices", 50, 0, 50);
    h_goodPV_PUWt->GetYaxis()->SetTitle("Events");   h_goodPV_PUWt->GetYaxis()->CenterTitle();
    h_goodPV_PUWt->GetXaxis()->SetTitle("nPV");      h_goodPV_PUWt->GetXaxis()->CenterTitle();
    h_goodPV_PUWt->Sumw2();

    //Photon Id Efficiency
    h_PassGenIsoMatch_EB = new TH1F("h_PassGenIsoMatch_EB", "Pt Distribution of Photons matching with Gen in EB", 1200, 0.0, 6000.0);
    h_PassGenIsoMatch_EB->GetYaxis()->SetTitle("Events/5 GeV");   h_PassGenIsoMatch_EB->GetYaxis()->CenterTitle();
    h_PassGenIsoMatch_EB->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassGenIsoMatch_EB->GetXaxis()->CenterTitle();
    h_PassGenIsoMatch_EB->Sumw2();

    h_PassGenIsoMatch_EE = new TH1F("h_PassGenIsoMatch_EE", "Pt Distribution of Photons matching with Gen in EE", 1200, 0.0, 6000.0);
    h_PassGenIsoMatch_EE->GetYaxis()->SetTitle("Events/5 GeV");   h_PassGenIsoMatch_EE->GetYaxis()->CenterTitle();
    h_PassGenIsoMatch_EE->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassGenIsoMatch_EE->GetXaxis()->CenterTitle();
    h_PassGenIsoMatch_EE->Sumw2();

    h_PassPhIdLoose_EB = new TH1F("h_PassPhIdLoose_EB", "Pt Distribution of Photons passing loose PhId in EB", 1200, 0.0, 6000.0);
    h_PassPhIdLoose_EB->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdLoose_EB->GetYaxis()->CenterTitle();
    h_PassPhIdLoose_EB->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdLoose_EB->GetXaxis()->CenterTitle();
    h_PassPhIdLoose_EB->Sumw2();

    h_PassPhIdLoose_EE = new TH1F("h_PassPhIdLoose_EE", "Pt Distribution of Photons passing loose PhId in EE", 1200, 0.0, 6000.0);
    h_PassPhIdLoose_EE->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdLoose_EE->GetYaxis()->CenterTitle();
    h_PassPhIdLoose_EE->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdLoose_EE->GetXaxis()->CenterTitle();
    h_PassPhIdLoose_EE->Sumw2();

    h_PassPhIdMedium_EB = new TH1F("h_PassPhIdMedium_EB", "Pt Distribution of Photons passing medium PhId in EB", 1200, 0.0, 6000.0);
    h_PassPhIdMedium_EB->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdMedium_EB->GetYaxis()->CenterTitle();
    h_PassPhIdMedium_EB->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdMedium_EB->GetXaxis()->CenterTitle();
    h_PassPhIdMedium_EB->Sumw2();

    h_PassPhIdMedium_EE = new TH1F("h_PassPhIdMedium_EE", "Pt Distribution of Photons passing medium PhId in EE", 1200, 0.0, 6000.0);
    h_PassPhIdMedium_EE->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdMedium_EE->GetYaxis()->CenterTitle();
    h_PassPhIdMedium_EE->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdMedium_EE->GetXaxis()->CenterTitle();
    h_PassPhIdMedium_EE->Sumw2();

    h_PassPhIdTight_EB = new TH1F("h_PassPhIdTight_EB", "Pt Distribution of Photons passing tight PhId in EB", 1200, 0.0, 6000.0);
    h_PassPhIdTight_EB->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdTight_EB->GetYaxis()->CenterTitle();
    h_PassPhIdTight_EB->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdTight_EB->GetXaxis()->CenterTitle();
    h_PassPhIdTight_EB->Sumw2();

    h_PassPhIdTight_EE = new TH1F("h_PassPhIdTight_EE", "Pt Distribution of Photons passing tight PhId in EE", 1200, 0.0, 6000.0);
    h_PassPhIdTight_EE->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdTight_EE->GetYaxis()->CenterTitle();
    h_PassPhIdTight_EE->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdTight_EE->GetXaxis()->CenterTitle();
    h_PassPhIdTight_EE->Sumw2();

    h_PassPhIdHighPt_EB = new TH1F("h_PassPhIdHighPt_EB", "Pt Distribution of Photons passing HighPt PhId in EB", 1200, 0.0, 6000.0);
    h_PassPhIdHighPt_EB->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdHighPt_EB->GetYaxis()->CenterTitle();
    h_PassPhIdHighPt_EB->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdHighPt_EB->GetXaxis()->CenterTitle();
    h_PassPhIdHighPt_EB->Sumw2();

    h_PassPhIdHighPt_EE = new TH1F("h_PassPhIdHighPt_EE", "Pt Distribution of Photons passing HighPt PhId in EE", 1200, 0.0, 6000.0);
    h_PassPhIdHighPt_EE->GetYaxis()->SetTitle("Events/5 GeV");   h_PassPhIdHighPt_EE->GetYaxis()->CenterTitle();
    h_PassPhIdHighPt_EE->GetXaxis()->SetTitle("P_{T}^{#gamma}"); h_PassPhIdHighPt_EE->GetXaxis()->CenterTitle();
    h_PassPhIdHighPt_EE->Sumw2();

    //Position of photon and jet
    h_PC = new TH1F("h_PC", "Photon Candidate", 10, 0, 10);
    h_PC->GetYaxis()->SetTitle("Events");                       h_PC->GetYaxis()->CenterTitle();
    h_PC->GetXaxis()->SetTitle("Position of Photon");           h_PC->GetXaxis()->CenterTitle();
    h_PC->Sumw2();

    h_JC = new TH1F("h_JC", "Jet Candidate", 20, 0, 20);
    h_JC->GetYaxis()->SetTitle("Events");                       h_JC->GetYaxis()->CenterTitle();
    h_JC->GetXaxis()->SetTitle("Position of Jet");              h_JC->GetXaxis()->CenterTitle();
    h_JC->Sumw2();

    //Pileup Reweighting
    h_DataPUNormDist = new TH1F("h_DataPUNormDist", "Normalized Data PileUp Distribution", 50, 0, 50);
    h_DataPUNormDist->GetYaxis()->SetTitle("Events");           h_DataPUNormDist->GetYaxis()->CenterTitle();
    h_DataPUNormDist->GetXaxis()->SetTitle("nPUV");             h_DataPUNormDist->GetXaxis()->CenterTitle();
    h_DataPUNormDist->Sumw2();
  
    h_MCPUNormDist = new TH1F("h_MCPUNormDist", "Normalized MC PileUp Distribution",50, 0, 50);
    h_MCPUNormDist->GetYaxis()->SetTitle("Events");             h_MCPUNormDist->GetYaxis()->CenterTitle();
    h_MCPUNormDist->GetXaxis()->SetTitle("nPUV");               h_MCPUNormDist->GetXaxis()->CenterTitle();
    h_MCPUNormDist->Sumw2();
  
    h_PUScaleFactor = new TH1F("h_PUScaleFactor", "PileUp Scale Factors Distribution", 50, 0, 50);
    h_PUScaleFactor->GetYaxis()->SetTitle("SF");      h_PUScaleFactor->GetYaxis()->CenterTitle();
    h_PUScaleFactor->GetXaxis()->SetTitle("nPUV");              h_PUScaleFactor->GetXaxis()->CenterTitle();
    h_PUScaleFactor->Sumw2();



    //Defining the histogram from filling number of events after various cuts
    const int nbins = 15;
    const int nbinsWt = 15;
    TString CutFlowLabels[nbins] = {"Total", "HLT", "PrimaryVtx", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "DPhi", "DEta", "GJInvtMass", "MET500Cut", "MET250Cut", "MET100Cut", "MET50Cut"};
    TString CutFlowLabelsWithWts[nbinsWt] = {"Total", "HLT", "PrimaryVtx", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "DPhi", "DEta", "GJInvtMass", "MET500Cut", "MET250Cut", "MET100Cut", "MET50Cut"};

    h_CutFlow = new TH1F("h_CutFlow", "Events Passing Various Cuts", nbins, 0, nbins);
    h_CutFlow->GetYaxis()->SetTitle("Events");         h_CutFlow->GetYaxis()->CenterTitle();
    for(int i = 0; i < nbins; i++){
      h_CutFlow->GetXaxis()->SetBinLabel(i+1, CutFlowLabels[i]);
    }

    h_CutFlowWithWts = new TH1F("h_CutFlowWithWts", "Events Passing Various Cuts With Weights", nbinsWt, 0, nbinsWt);
    h_CutFlowWithWts->GetYaxis()->SetTitle("Events");         h_CutFlowWithWts->GetYaxis()->CenterTitle();
    for(int i = 0; i < nbinsWt; i++){
      h_CutFlowWithWts->GetXaxis()->SetBinLabel(i+1, CutFlowLabelsWithWts[i]);
    }


}

void PostAnalyzer_MC::WriteHistograms(){

  for( Int_t hist = 0; hist < 2; ++hist){

    h_PhotonPt[hist]               ->Write();
    h_PhotonEta[hist]              ->Write();
    h_PhotonPhi[hist]              ->Write();
    h_Photon_SigmaIEtaIEta[hist]   ->Write();
    h_Photon_R9[hist]              ->Write();
    h_Photon_HoverE[hist]          ->Write();
    h_Photon_EleVeto[hist]         ->Write();
    h_Photon_CorrPFChIso[hist]     ->Write();
    h_Photon_CorrPFNeuIso[hist]    ->Write();
    h_Photon_CorrPFPhoIso[hist]    ->Write();

    h_JetPt[hist]                  ->Write();
    h_JetEta[hist]                 ->Write();
    h_JetPhi[hist]                 ->Write();
    h_Jet_NHEF[hist]               ->Write();
    h_Jet_NEEF[hist]               ->Write();
    h_Jet_NConst[hist]             ->Write();
    h_Jet_CHEF[hist]               ->Write();
    h_Jet_ChMult[hist]             ->Write();
    h_Jet_CEEF[hist]               ->Write();

    h_GJetInvtMass_bin40[hist]     ->Write();
    h_GJetInvtMass_VarBin[hist]    ->Write();
    h_GJetInvtMass_UnitBin[hist]   ->Write();
    h_GJet_dEta[hist]              ->Write();
    h_GJet_dPhi[hist]              ->Write();
    h_GJet_dR[hist]                ->Write();

    h_Gen_GJetInvtMass_bin40[hist]     ->Write();
    h_Gen_GJetInvtMass_VarBin[hist]    ->Write();
    h_Gen_GJetInvtMass_UnitBin[hist]   ->Write();

    h_GJetMassResolution[hist]         ->Write();

  }

  for( Int_t hist = 0; hist < 3; ++hist){ 

    h_goodPV_TotalWt[hist]         ->Write();
    h_nPhotons[hist]               ->Write();
    h_nJets[hist]                  ->Write();

  }

  h_trueNumofInt              ->Write();
  h_goodPV                    ->Write();
  h_goodPV_LumiWt             ->Write();
  h_goodPV_PUWt               ->Write();


  h_PassPhIdLoose_EB          ->Write();
  h_PassPhIdLoose_EE          ->Write();
  h_PassPhIdMedium_EB         ->Write();
  h_PassPhIdMedium_EE         ->Write();

  h_PC                        ->Write();
  h_JC                        ->Write();

  h_CutFlow                   ->Write();
  h_CutFlowWithWts            ->Write();
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
##./${filenameTag}.exe


###Submit jobs

chmod 775 MakeCondorFiles_local.csh

./MakeCondorFiles_local.csh ${filenameTag}


((sampleIndex++))

done ##end of for loop##
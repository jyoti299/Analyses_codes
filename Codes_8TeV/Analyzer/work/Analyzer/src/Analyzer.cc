//______________________________________________________________________________________________________________________________________________
// -*- C++ -*-                                                                                                                                  |
//                                                                                                                                              |
// Package:    Analyzer                                                                                                                         |
// Class:      Analyzer                                                                                                                         |
//                                                                                                                                              |
//class Analyzer Analyzer.cc work/Analyzer/src/Analyzer.cc                                                                                      |
//                                                                                                                                              |
// Description: [one line class summary]                                                                                                        |
//                                                                                                                                              |
// Implementation:                                                                                                                              |
//     [Notes on implementation                                                                                                                 |
//                                                                                                                                              |
//                                                                                                                                              |
// Original Author:  rocky garg                                                                                                                 |
//         Created:  Tue Feb 18 00:52:08 CST 2014                                                                                               |
// $Id$                                                                                                                                         |
//                                                                                                                                              |
//                                                                                                                                              |
//______________________________________________________________________________________________________________________________________________|

// user include files
#include "TString.h"                                                                                                                       
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>


#include "Analyzer.h"



using namespace std;
using namespace ROOT;
using namespace edm;
using namespace ROOT::Math::VectorUtil ;



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig):
           debug_(iConfig.getUntrackedParameter<int>("debug", 0)),
           OutFile_(iConfig.getUntrackedParameter<std::string>("OutFile", "")), 

	   //bools
	   isItAOD_(iConfig.getUntrackedParameter<bool>("isItAOD", false)),
	   saveLuminosityInfo_(iConfig.getUntrackedParameter<bool>("saveLuminosityInfo", false)),
           runMCPileUp_(iConfig.getUntrackedParameter<bool>("runMCPileUp", false)),
           runPhotons_(iConfig.getUntrackedParameter<bool>("runPhotons", false)),
	   savePhotonCrystalInfo_(iConfig.getUntrackedParameter<bool>("savePhotonCrystalInfo", false)),
	   savePhotonPFIsolation_(iConfig.getUntrackedParameter<bool>("savePhotonPFIsolation", false)),
	   savePhotonRecHitsInfo_(iConfig.getUntrackedParameter<bool>("savePhotonRecHitsInfo", false)),
	   runJets_(iConfig.getUntrackedParameter<bool>("runJets", false)),
	   runCaloPatJets_(iConfig.getUntrackedParameter<bool>("runCaloPatJets", false)),
	   runPFPatJets_(iConfig.getUntrackedParameter<bool>("runPFPatJets", false)),
	   runCaloRecoJets_(iConfig.getUntrackedParameter<bool>("runCaloRecoJets", false)),
	   runPFRecoJets_(iConfig.getUntrackedParameter<bool>("runPFRecoJets", false)),
	   savePUJetIdInfo_PFPat_(iConfig.getUntrackedParameter<bool>("savePUJetIdInfo_PFPat", false)),
	   saveBTaggingInfo_CaloPat_(iConfig.getUntrackedParameter<bool>("saveBTaggingInfo_CaloPat", false)),
	   saveBTaggingInfo_PFPat_(iConfig.getUntrackedParameter<bool>("saveBTaggingInfo_PFPat", false)),
	   saveBTaggingInfo_CaloReco_(iConfig.getUntrackedParameter<bool>("saveBTaggingInfo_CaloReco", false)),
	   saveBTaggingInfo_PFReco_(iConfig.getUntrackedParameter<bool>("saveBTaggingInfo_PFReco", false)),
	   runVertex_(iConfig.getUntrackedParameter<bool>("runVertex", false)),
	   runTracks_(iConfig.getUntrackedParameter<bool>("runTracks", false)),
	   runScraping_(iConfig.getUntrackedParameter<bool>("runScraping", false)),
	   runHLT_(iConfig.getUntrackedParameter<bool>("runHLT", false)),


	   //InputTags        
	   MCpileupLabel_(iConfig.getUntrackedParameter<edm::InputTag>("MCpileupTag")),
	   patPhoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("patPhoTag")),
           recoPhoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("recoPhoTag")),
	   pfCandidateLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateTag")),
	   vertexLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag")),
	   BSLabel_(iConfig.getUntrackedParameter<edm::InputTag>("BSTag")),
	   ConvPhoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("ConvPhoTag")),
	   GsfEleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("GsfEleTag")),
	   BarrelrechitLabel_(iConfig.getUntrackedParameter<edm::InputTag>("BarrelrechitTag")),
           EndcaprechitLabel_(iConfig.getUntrackedParameter<edm::InputTag>("EndcaprechitTag")),
	   CaloPatJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloPatJetTag")),
	   PFPatJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFPatJetTag")),
	   CaloRecoJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloRecoJetTag")),
	   PFRecoJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFRecoJetTag")),
	   //b tagging Input tags
	   CaloReco_TrkCountingHighEffBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_TrkCountingHighEffBJetTag")),
	   CaloReco_TrkCountingHighPurBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_TrkCountingHighPurBJetTag")),
	   CaloReco_SimpleSecVtxHighEffBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_SimpleSecVtxHighEffBJetTag")),
	   CaloReco_SimpleSecVtxHighPurBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_SimpleSecVtxHighPurBJetTag")),
	   CaloReco_CombinedSecVtxBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_CombinedSecVtxBJetTag")),
	   CaloReco_CombinedSecVtxMVABJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_CombinedSecVtxMVABJetTag")),
	   CaloReco_JetProbBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_JetProbBJetTag")),
	   CaloReco_JetBProbBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_JetBProbBJetTag")),
	   CaloReco_SoftElectronBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_SoftElectronBJetTag")),
	   CaloReco_SoftMuonBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("CaloReco_SoftMuonBJetTag")),
	   PFReco_TrkCountingHighEffBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_TrkCountingHighEffBJetTag")),
	   PFReco_TrkCountingHighPurBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_TrkCountingHighPurBJetTag")),
	   PFReco_SimpleSecVtxHighEffBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_SimpleSecVtxHighEffBJetTag")),
	   PFReco_SimpleSecVtxHighPurBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_SimpleSecVtxHighPurBJetTag")),
	   PFReco_CombinedSecVtxBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_CombinedSecVtxBJetTag")),
	   PFReco_CombinedSecVtxMVABJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_CombinedSecVtxMVABJetTag")),
	   PFReco_JetProbBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_JetProbBJetTag")),
	   PFReco_JetBProbBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_JetBProbBJetTag")),
	   PFReco_SoftElectronBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_SoftElectronBJetTag")),
	   PFReco_SoftMuonBJetTagsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFReco_SoftMuonBJetTag")),
	   //Track Tag
	   TracksLabel_(iConfig.getUntrackedParameter<edm::InputTag>("TracksTag")),
	   //triggerEvent tag 
	   triggerEventLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
	   //triggerResults Tag
	   triggerResultsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag")),

           //severity flag for shower shape cleaned
           flagnamesEB(iConfig.getParameter<std::vector<std::string> >("RecHitFlagToBeExcludedEB")),
           flagnamesEE(iConfig.getParameter<std::vector<std::string> >("RecHitFlagToBeExcludedEE")),
           severitynamesEB(iConfig.getParameter<std::vector<std::string> >("RecHitSeverityToBeExcludedEB")),
           severitynamesEE(iConfig.getParameter<std::vector<std::string> >("RecHitSeverityToBeExcludedEE")),

	   //Reco Jet Correction Service
	   CaloRecoJetCorrectionService_(iConfig.getUntrackedParameter<std::string>("CaloRecoJetCorrectionService")),
           PFRecoJetCorrectionService_(iConfig.getUntrackedParameter<std::string>("PFRecoJetCorrectionService"))
{
  nEvents = 0;
 
  //Passing flags to the variables for sigmaIphiIphi etc estimation
  flagsexclEB_= StringToEnumValue<EcalRecHit::Flags>(flagnamesEB);
  flagsexclEE_= StringToEnumValue<EcalRecHit::Flags>(flagnamesEE);
  severitiesexclEB_= StringToEnumValue<EcalSeverityLevel::SeverityLevel>(severitynamesEB);
  severitiesexclEE_= StringToEnumValue<EcalSeverityLevel::SeverityLevel>(severitynamesEE);


  //PFIsolation Initialization for Photons
  phoIsolator03.initializePhotonIsolation(kTRUE); //NOTE: this automatically set all the correct default veto values
  phoIsolator03.setConeSize(0.3);



}


//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob()
{
  f       = new TFile(OutFile_.c_str(),"RECREATE");
  myEvent = new TTree("myEvent","A Tree with Histograms");

  //Declearing Branches

  myEvent->Branch("RunNumber", &RunNumber, "RunNumber/i");
  myEvent->Branch("EventNumber", &EventNumber, "EventNumber/i");
  myEvent->Branch("LumiNumber", &LumiNumber, "LumiNumber/i");
  myEvent->Branch("BXNumber", &BXNumber, "BXNumber/i");
  
  if(saveLuminosityInfo_){
    myEvent->Branch("totalIntensityBeam1", &totalIntensityBeam1, "totalIntensityBeam1/i");
    myEvent->Branch("totalIntensityBeam2", &totalIntensityBeam2, "totalIntensityBeam2/i");
    myEvent->Branch("avgInsDelLumi", &avgInsDelLumi, "avgInsDelLumi/F");
    myEvent->Branch("avgInsDelLumiErr", &avgInsDelLumiErr, "avgInsDelLumiErr/F");
    myEvent->Branch("avgInsRecLumi", &avgInsRecLumi, "avgInsRecLumi/F");
    myEvent->Branch("avgInsRecLumiErr", &avgInsRecLumiErr, "avgInsRecLumiErr/F");
  }

  if(runMCPileUp_){
    myEvent->Branch("npuVertices",&npuVertices,"npuVertices/I");
    myEvent->Branch("npuVerticesm1",&npuVerticesm1,"npuVerticesm1/I");
    myEvent->Branch("npuVerticesp1",&npuVerticesp1,"npuVerticesp1/I");
    myEvent->Branch("npuVerticespm2",&npuVerticespm2,"npuVerticespm2/I");                                                        
    myEvent->Branch("trueNumofInteractions",&trueNumofInteractions,"trueNumofInteractions/F");
  }
  
  if(runPhotons_){

    myEvent->Branch("Photon_n", &Photon_n, "Photon_n/I");

    //Branches for Kinematic Variables
    myEvent->Branch("Photon_E", &Photon_E);
    myEvent->Branch("Photon_et", &Photon_et);
    myEvent->Branch("Photon_pt", &Photon_pt);
    myEvent->Branch("Photon_eta", &Photon_eta);
    myEvent->Branch("Photon_phi", &Photon_phi);
    myEvent->Branch(" Photon_theta", &Photon_theta);
    myEvent->Branch("Photon_px", &Photon_px);
    myEvent->Branch("Photon_py", &Photon_py);
    myEvent->Branch("Photon_pz", &Photon_pz);
    myEvent->Branch("Photon_vx", &Photon_vx);
    myEvent->Branch("Photon_vy", &Photon_vy);
    myEvent->Branch("Photon_vz", &Photon_vz);

    //Branches for Shower Shape Variables
    myEvent->Branch("Photon_r9", &Photon_r9);
    myEvent->Branch("Photon_maxEnergyXtal", &Photon_maxEnergyXtal);
    myEvent->Branch("Photon_e1x5", &Photon_e1x5);
    myEvent->Branch("Photon_e2x5", &Photon_e2x5);
    myEvent->Branch("Photon_e3x3", &Photon_e3x3);
    myEvent->Branch("Photon_e5x5", &Photon_e5x5);
    myEvent->Branch("Photon_r1x5", &Photon_r1x5);
    myEvent->Branch("Photon_r2x5", &Photon_r2x5);
    myEvent->Branch("Photon_SigmaEtaEta", &Photon_SigmaEtaEta);
    myEvent->Branch("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta);
    myEvent->Branch("Photon_SigmaEtaPhi", &Photon_SigmaEtaPhi);
    myEvent->Branch("Photon_SigmaIEtaIPhi", &Photon_SigmaIEtaIPhi);
    myEvent->Branch("Photon_SigmaPhiPhi", &Photon_SigmaPhiPhi);
    myEvent->Branch("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi);
    myEvent->Branch("Photon_roundness", &Photon_roundness);
    myEvent->Branch("Photon_angle", &Photon_angle);
    myEvent->Branch("Photon_swissCross", &Photon_swissCross);
    myEvent->Branch("Photon_s9", &Photon_s9);
    myEvent->Branch("Photon_e4Overe1", &Photon_e4Overe1);
    myEvent->Branch("Photon_e6Overe2", &Photon_e6Overe2);
    myEvent->Branch("Photon_e2Overe9", &Photon_e2Overe9);
    myEvent->Branch("Photon_rookFraction", &Photon_rookFraction);

    //Branches for Fiducial Flags Variables
    myEvent->Branch("Photon_isEB", &Photon_isEB);
    myEvent->Branch("Photon_isEE", &Photon_isEE);
    myEvent->Branch("Photon_isEBGap", &Photon_isEBGap);
    myEvent->Branch("Photon_isEEGap", &Photon_isEEGap);
    myEvent->Branch("Photon_isEBEEGap", &Photon_isEBEEGap);

    //Branches for Detector Isolation Variables
    myEvent->Branch("Photon_ecalRecHitSumEtConeDR03", &Photon_ecalRecHitSumEtConeDR03);
    myEvent->Branch("Photon_hcalTowerSumEtConeDR03", &Photon_hcalTowerSumEtConeDR03);
    myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR03", &Photon_hcalDepth1TowerSumEtConeDR03);
    myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR03", &Photon_hcalDepth2TowerSumEtConeDR03);
    myEvent->Branch("Photon_trkSumPtSolidConeDR03", &Photon_trkSumPtSolidConeDR03);
    myEvent->Branch("Photon_trkSumPtHollowConeDR03", &Photon_trkSumPtHollowConeDR03);
    myEvent->Branch("Photon_nTrkSolidConeDR03", &Photon_nTrkSolidConeDR03);
    myEvent->Branch("Photon_nTrkHollowConeDR03", &Photon_nTrkHollowConeDR03);
    myEvent->Branch("Photon_ecalRecHitSumEtConeDR04", &Photon_ecalRecHitSumEtConeDR04);
    myEvent->Branch("Photon_hcalTowerSumEtConeDR04", &Photon_hcalTowerSumEtConeDR04);
    myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR04", &Photon_hcalDepth1TowerSumEtConeDR04);
    myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR04", &Photon_hcalDepth2TowerSumEtConeDR04);
    myEvent->Branch("Photon_trkSumPtSolidConeDR04", &Photon_trkSumPtSolidConeDR04);
    myEvent->Branch("Photon_trkSumPtHollowConeDR04", &Photon_trkSumPtHollowConeDR04);
    myEvent->Branch("Photon_nTrkSolidConeDR04", &Photon_nTrkSolidConeDR04);
    myEvent->Branch("Photon_nTrkHollowConeDR04", &Photon_nTrkHollowConeDR04);

    //Branches for Photon Identification Variables
    myEvent->Branch("Photon_HoE", &Photon_HoE);
    myEvent->Branch("Photon_SingleTowerHoE", &Photon_SingleTowerHoE);
    myEvent->Branch("Photon_hasConvTrk", &Photon_hasConvTrk);
    myEvent->Branch("Photon_hasPixelSeed", &Photon_hasPixelSeed);
    myEvent->Branch("passedConvSafeElectronVeto", &passedConvSafeElectronVeto);

    //Branches for Photon Super Cluster Variables
    myEvent->Branch("Photon_SC_nOfBasicClusters", &Photon_SC_nOfBasicClusters);
    myEvent->Branch("Photon_SC_rawEnergy", &Photon_SC_rawEnergy);
    myEvent->Branch("Photon_SC_preShowerEnergy", &Photon_SC_preShowerEnergy);
    myEvent->Branch("Photon_SC_energy", &Photon_SC_energy);
    myEvent->Branch("Photon_SC_eta", &Photon_SC_eta);
    myEvent->Branch("Photon_SC_phi", &Photon_SC_phi);
    myEvent->Branch("Photon_SC_x", &Photon_SC_x);
    myEvent->Branch("Photon_SC_y", &Photon_SC_y);
    myEvent->Branch("Photon_SC_z", &Photon_SC_z);
    myEvent->Branch("Photon_SC_etaWidth", &Photon_SC_etaWidth);
    myEvent->Branch("Photon_SC_phiWidth", &Photon_SC_phiWidth);

    //Photon MIP Variables
    myEvent->Branch("Photon_mipChi2", &Photon_mipChi2);
    myEvent->Branch("Photon_mipTotEnergy", &Photon_mipTotEnergy);
    myEvent->Branch("Photon_mipSlope", &Photon_mipSlope);
    myEvent->Branch("Photon_mipIntercept", &Photon_mipIntercept);
    myEvent->Branch("Photon_mipNhitCone", &Photon_mipNhitCone);
    myEvent->Branch("Photon_mipIsHalo", &Photon_mipIsHalo);

    //Branches for Converted Photon Variables
    myEvent->Branch("Photon_nConvTracks", &Photon_nConvTracks);
    myEvent->Branch("Photon_isConverted", &Photon_isConverted);
    myEvent->Branch("Photon_pairInvariantMass", &Photon_pairInvariantMass);
    myEvent->Branch("Photon_pairCotThetaSeparation", &Photon_pairCotThetaSeparation);
    myEvent->Branch("Photon_pairMomentum_x", &Photon_pairMomentum_x);
    myEvent->Branch("Photon_pairMomentum_y", &Photon_pairMomentum_y);
    myEvent->Branch("Photon_pairMomentum_z", &Photon_pairMomentum_z);
    myEvent->Branch("Photon_EoverP", &Photon_EoverP);
    myEvent->Branch("Photon_conv_vx", &Photon_conv_vx);
    myEvent->Branch("Photon_conv_vy", &Photon_conv_vy);
    myEvent->Branch("Photon_conv_vz", &Photon_conv_vz);
    myEvent->Branch("Photon_zOfPrimaryVtxFromTrks", &Photon_zOfPrimaryVtxFromTrks);
    myEvent->Branch("Photon_distOfMinimumApproach", &Photon_distOfMinimumApproach);
    myEvent->Branch("Photon_dPhiTracksAtVtx", &Photon_dPhiTracksAtVtx);
    myEvent->Branch("Photon_dPhiTracksAtEcal", &Photon_dPhiTracksAtEcal);
    myEvent->Branch("Photon_dEtaTracksAtEcal", &Photon_dEtaTracksAtEcal);

    //Branches for Photon crystal variables
    if(savePhotonCrystalInfo_){
      myEvent->Branch("Photon_nCrystals", &Photon_nCrystals);
      myEvent->Branch("Photon_xtal_timing", &Photon_xtal_timing);
      myEvent->Branch("Photon_xtal_timeErr", &Photon_xtal_timeErr);
      myEvent->Branch("Photon_avgTimeAllxtals", &Photon_avgTimeAllxtals);
      myEvent->Branch("Photon_xtal_energy", &Photon_xtal_energy);
      myEvent->Branch("Photon_xtal_EBieta", &Photon_xtal_EBieta);
      myEvent->Branch("Photon_xtal_EBiphi", &Photon_xtal_EBiphi);
      myEvent->Branch("Photon_xtal_EBrecoFlag", &Photon_xtal_EBrecoFlag);
    }

    //Branches for PFIsolation Variables
    if(savePhotonPFIsolation_){
      myEvent->Branch("PFIsoPhoton03", &PFIsoPhoton03);
      myEvent->Branch("PFIsoNeutral03", &PFIsoNeutral03);
      myEvent->Branch("PFIsoCharged03", &PFIsoCharged03);
      myEvent->Branch("PFIsoSum03", &PFIsoSum03);
      myEvent->Branch("PFIsoChargedWorstvtx03", &PFIsoChargedWorstvtx03);
    }

  }

  if(runJets_){
    if(runCaloPatJets_){

      myEvent->Branch("CaloPatJet_n", &CaloPatJet_n, "CaloPatJet_n/I");

      //Kinematic Variables
      myEvent->Branch("CaloPatJet_E", &CaloPatJet_E);
      myEvent->Branch("CaloPatJet_et", &CaloPatJet_et);
      myEvent->Branch("CaloPatJet_pt", &CaloPatJet_pt);
      myEvent->Branch("CaloPatJet_eta", &CaloPatJet_eta);
      myEvent->Branch("CaloPatJet_phi", &CaloPatJet_phi);
      myEvent->Branch("CaloPatJet_theta", &CaloPatJet_theta);
      myEvent->Branch("CaloPatJet_px", &CaloPatJet_px);
      myEvent->Branch("CaloPatJet_py", &CaloPatJet_py);
      myEvent->Branch("CaloPatJet_pz", &CaloPatJet_pz);
      myEvent->Branch("CaloPatJet_vx", &CaloPatJet_vx);
      myEvent->Branch("CaloPatJet_vy", &CaloPatJet_vy);
      myEvent->Branch("CaloPatJet_vz", &CaloPatJet_vz);

      //Id and other Variables for Calo Pat Jets
      myEvent->Branch("CaloPatJet_emEnergyFraction", &CaloPatJet_emEnergyFraction);
      myEvent->Branch("CaloPatJet_emEnergyInEB", &CaloPatJet_emEnergyInEB);
      myEvent->Branch("CaloPatJet_emEnergyInEE", &CaloPatJet_emEnergyInEE);
      myEvent->Branch("CaloPatJet_emEnergyInHF", &CaloPatJet_emEnergyInHF);
      myEvent->Branch("CaloPatJet_HadronicEnergyFraction", &CaloPatJet_HadronicEnergyFraction);
      myEvent->Branch("CaloPatJet_NConstituents", &CaloPatJet_NConstituents);
      myEvent->Branch("CaloPatJet_HadronicEnergyInHB", &CaloPatJet_HadronicEnergyInHB);
      myEvent->Branch("CaloPatJet_HadronicEnergyInHE", &CaloPatJet_HadronicEnergyInHE);
      myEvent->Branch("CaloPatJet_HadronicEnergyInHF", &CaloPatJet_HadronicEnergyInHF);
      myEvent->Branch("CaloPatJet_HadronicEnergyInHO", &CaloPatJet_HadronicEnergyInHO);
      myEvent->Branch("CaloPatJet_maxEnergyInEmTowers", &CaloPatJet_maxEnergyInEmTowers);
      myEvent->Branch("CaloPatJet_maxEnergyInHadTowers", &CaloPatJet_maxEnergyInHadTowers);
      myEvent->Branch("CaloPatJet_nConstituents60E", &CaloPatJet_nConstituents60E);
      myEvent->Branch("CaloPatJet_nConstituents90E", &CaloPatJet_nConstituents90E);
      myEvent->Branch("CaloPatJet_towersArea", &CaloPatJet_towersArea);
      myEvent->Branch("CaloPatJet_nTowers", &CaloPatJet_nTowers);
      myEvent->Branch("CaloPatJet_fEinHottestHPD", &CaloPatJet_fEinHottestHPD);
      myEvent->Branch("CaloPatJet_fEinHottestRBX", &CaloPatJet_fEinHottestRBX);
      myEvent->Branch("CaloPatJet_nHitsCarrying90E", &CaloPatJet_nHitsCarrying90E);
      myEvent->Branch("CaloPatJet_nHitsinTowersCarrying90E", &CaloPatJet_nHitsinTowersCarrying90E);
      myEvent->Branch("CaloPatJet_RHF", &CaloPatJet_RHF);

      //JEC Uncertainity variable
      myEvent->Branch("CaloPatJet_jecUncertainity", &CaloPatJet_jecUncertainity);

      if(saveBTaggingInfo_CaloPat_){
	myEvent->Branch("CaloPatJet_BJetDiscrByTrackCountingHighEff", &CaloPatJet_BJetDiscrByTrackCountingHighEff);
	myEvent->Branch("CaloPatJet_BJetDiscrByTrackCountingHighPur", &CaloPatJet_BJetDiscrByTrackCountingHighPur);
	myEvent->Branch("CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighEff", &CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighEff);
	myEvent->Branch("CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighPur", &CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighPur);
	myEvent->Branch("CaloPatJet_BJetDiscrByCombinedSecondaryVertexHighEff", &CaloPatJet_BJetDiscrByCombinedSecondaryVertexHighEff);
	myEvent->Branch("CaloPatJet_BJetDiscrByCombinedSecondaryVertexMVA", &CaloPatJet_BJetDiscrByCombinedSecondaryVertexMVA);
	myEvent->Branch("CaloPatJet_BJetDiscrByjetProbabilityBJetTags", &CaloPatJet_BJetDiscrByjetProbabilityBJetTags);
	myEvent->Branch("CaloPatJet_BJetDiscrByjetBProbabilityBJetTags", &CaloPatJet_BJetDiscrByjetBProbabilityBJetTags);
	myEvent->Branch("CaloPatJet_BJetDiscrBySoftElectronBJetTags", &CaloPatJet_BJetDiscrBySoftElectronBJetTags);
	myEvent->Branch("CaloPatJet_BJetDiscrBySoftMuonBJetTags", &CaloPatJet_BJetDiscrBySoftMuonBJetTags);
	myEvent->Branch("CaloPatJet_JetPartonFlavor", &CaloPatJet_JetPartonFlavor);
      }

    }

    if(runPFPatJets_){
      
      myEvent->Branch("PFPatJet_n", &PFPatJet_n, "PFPatJet_n/I");

      //Kinematic Variables
      myEvent->Branch("PFPatJet_E", &PFPatJet_E);
      myEvent->Branch("PFPatJet_et", &PFPatJet_et);
      myEvent->Branch("PFPatJet_pt", &PFPatJet_pt);
      myEvent->Branch("PFPatJet_eta", &PFPatJet_eta);
      myEvent->Branch("PFPatJet_phi", &PFPatJet_phi);
      myEvent->Branch("PFPatJet_theta", &PFPatJet_theta);
      myEvent->Branch("PFPatJet_px", &PFPatJet_px);
      myEvent->Branch("PFPatJet_py", &PFPatJet_py);
      myEvent->Branch("PFPatJet_pz", &PFPatJet_pz);
      myEvent->Branch("PFPatJet_vx", &PFPatJet_vx);
      myEvent->Branch("PFPatJet_vy", &PFPatJet_vy);
      myEvent->Branch("PFPatJet_vz", &PFPatJet_vz);

      //Id and other Variables for PF Pat Jets
      myEvent->Branch("PFPatJet_ChargedEmEnergy", &PFPatJet_ChargedEmEnergy);
      myEvent->Branch("PFPatJet_ChargedEmEnergyFrac", &PFPatJet_ChargedEmEnergyFrac);
      myEvent->Branch("PFPatJet_ChargedHadEnergy", &PFPatJet_ChargedHadEnergy);
      myEvent->Branch("PFPatJet_ChargedHadEnergyFrac", &PFPatJet_ChargedHadEnergyFrac);
      myEvent->Branch("PFPatJet_ChargedHadMult", &PFPatJet_ChargedHadMult);
      myEvent->Branch("PFPatJet_ChargedMult", &PFPatJet_ChargedMult);
      myEvent->Branch("PFPatJet_NConstituents", &PFPatJet_NConstituents);
      myEvent->Branch("PFPatJet_HFEMEnergy", &PFPatJet_HFEMEnergy);
      myEvent->Branch("PFPatJet_HFEMEnergyFrac", &PFPatJet_HFEMEnergyFrac);
      myEvent->Branch("PFPatJet_HFEMMult", &PFPatJet_HFEMMult);
      myEvent->Branch("PFPatJet_HFHadEnergy", &PFPatJet_HFHadEnergy);
      myEvent->Branch("PFPatJet_HFHadEnergyFrac", &PFPatJet_HFHadEnergyFrac);
      myEvent->Branch("PFPatJet_HFHadMult", &PFPatJet_HFHadMult);
      myEvent->Branch("PFPatJet_NeutralEmEnergy", &PFPatJet_NeutralEmEnergy);
      myEvent->Branch("PFPatJet_NeutralEmEnergyFrac", &PFPatJet_NeutralEmEnergyFrac);
      myEvent->Branch("PFPatJet_NeutralHadEnergy", &PFPatJet_NeutralHadEnergy);
      myEvent->Branch("PFPatJet_NeutralHadEnergyFrac", &PFPatJet_NeutralHadEnergyFrac);
      myEvent->Branch("PFPatJet_NeutralHadMult", &PFPatJet_NeutralHadMult);
      myEvent->Branch("PFPatJet_NeutralMult", &PFPatJet_NeutralMult);

      //JEC Uncertainity variable
      myEvent->Branch("PFPatJet_jecUncertainity", &PFPatJet_jecUncertainity);

      if(savePUJetIdInfo_PFPat_){
	myEvent->Branch("PFPatJet_puJetIdCutBased_MVA", &PFPatJet_puJetIdCutBased_MVA);
	myEvent->Branch("PFPatJet_puJetIdSimple_MVA", &PFPatJet_puJetIdSimple_MVA);
	myEvent->Branch("PFPatJet_puJetIdFull_MVA", &PFPatJet_puJetIdFull_MVA);
	myEvent->Branch("PFPatJet_PassPUJetIdCutBased_loose", &PFPatJet_PassPUJetIdCutBased_loose);
	myEvent->Branch("PFPatJet_PassPUJetIdCutBased_medium", &PFPatJet_PassPUJetIdCutBased_medium);
	myEvent->Branch("PFPatJet_PassPUJetIdCutBased_tight", &PFPatJet_PassPUJetIdCutBased_tight);
	myEvent->Branch("PFPatJet_PassPUJetIdSimple_loose", &PFPatJet_PassPUJetIdSimple_loose);
	myEvent->Branch("PFPatJet_PassPUJetIdSimple_medium", &PFPatJet_PassPUJetIdSimple_medium);
	myEvent->Branch("PFPatJet_PassPUJetIdSimple_tight", &PFPatJet_PassPUJetIdSimple_tight);
	myEvent->Branch("PFPatJet_PassPUJetIdFull_loose", &PFPatJet_PassPUJetIdFull_loose);
	myEvent->Branch("PFPatJet_PassPUJetIdFull_medium", &PFPatJet_PassPUJetIdFull_medium);
	myEvent->Branch("PFPatJet_PassPUJetIdFull_tight", &PFPatJet_PassPUJetIdFull_tight);
      }

      if(saveBTaggingInfo_PFPat_){
	myEvent->Branch("PFPatJet_BJetDiscrByTrackCountingHighEff", &PFPatJet_BJetDiscrByTrackCountingHighEff);
	myEvent->Branch("PFPatJet_BJetDiscrByTrackCountingHighPur", &PFPatJet_BJetDiscrByTrackCountingHighPur);
	myEvent->Branch("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff);
	myEvent->Branch("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur);
	myEvent->Branch("PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff);
	myEvent->Branch("PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA);
	myEvent->Branch("PFPatJet_BJetDiscrByjetProbabilityBJetTags", &PFPatJet_BJetDiscrByjetProbabilityBJetTags);
	myEvent->Branch("PFPatJet_BJetDiscrByjetBProbabilityBJetTags", &PFPatJet_BJetDiscrByjetBProbabilityBJetTags);
	myEvent->Branch("PFPatJet_BJetDiscrBySoftElectronBJetTags", &PFPatJet_BJetDiscrBySoftElectronBJetTags);
	myEvent->Branch("PFPatJet_BJetDiscrBySoftMuonBJetTags", &PFPatJet_BJetDiscrBySoftMuonBJetTags);
	myEvent->Branch("PFPatJet_JetPartonFlavor", &PFPatJet_JetPartonFlavor);
      }
    }

    if(runCaloRecoJets_){

      myEvent->Branch("CaloRecoJet_n", &CaloRecoJet_n, "CaloRecoJet_n/I");

      //Kinematic Variables
      myEvent->Branch("CaloRecoJet_E", &CaloRecoJet_E);
      myEvent->Branch("CaloRecoJet_et", &CaloRecoJet_et);
      myEvent->Branch("CaloRecoJet_pt", &CaloRecoJet_pt);
      myEvent->Branch("CaloRecoJet_eta", &CaloRecoJet_eta);
      myEvent->Branch("CaloRecoJet_phi", &CaloRecoJet_phi);
      myEvent->Branch("CaloRecoJet_theta", &CaloRecoJet_theta);
      myEvent->Branch("CaloRecoJet_px", &CaloRecoJet_px);
      myEvent->Branch("CaloRecoJet_py", &CaloRecoJet_py);
      myEvent->Branch("CaloRecoJet_pz", &CaloRecoJet_pz);
      myEvent->Branch("CaloRecoJet_vx", &CaloRecoJet_vx);
      myEvent->Branch("CaloRecoJet_vy", &CaloRecoJet_vy);
      myEvent->Branch("CaloRecoJet_vz", &CaloRecoJet_vz);

      //Id and other Variables for Calo Reco Jets
      myEvent->Branch("CaloRecoJet_emEnergyFraction", &CaloRecoJet_emEnergyFraction);
      myEvent->Branch("CaloRecoJet_emEnergyInEB", &CaloRecoJet_emEnergyInEB);
      myEvent->Branch("CaloRecoJet_emEnergyInEE", &CaloRecoJet_emEnergyInEE);
      myEvent->Branch("CaloRecoJet_emEnergyInHF", &CaloRecoJet_emEnergyInHF);
      myEvent->Branch("CaloRecoJet_HadronicEnergyFraction", &CaloRecoJet_HadronicEnergyFraction);
      myEvent->Branch("CaloRecoJet_NConstituents", &CaloRecoJet_NConstituents);
      myEvent->Branch("CaloRecoJet_HadronicEnergyInHB", &CaloRecoJet_HadronicEnergyInHB);
      myEvent->Branch("CaloRecoJet_HadronicEnergyInHE", &CaloRecoJet_HadronicEnergyInHE);
      myEvent->Branch("CaloRecoJet_HadronicEnergyInHF", &CaloRecoJet_HadronicEnergyInHF);
      myEvent->Branch("CaloRecoJet_HadronicEnergyInHO", &CaloRecoJet_HadronicEnergyInHO);
      myEvent->Branch("CaloRecoJet_maxEnergyInEmTowers", &CaloRecoJet_maxEnergyInEmTowers);
      myEvent->Branch("CaloRecoJet_maxEnergyInHadTowers", &CaloRecoJet_maxEnergyInHadTowers);
      myEvent->Branch("CaloRecoJet_nConstituents60E", &CaloRecoJet_nConstituents60E);
      myEvent->Branch("CaloRecoJet_nConstituents90E", &CaloRecoJet_nConstituents90E);
      myEvent->Branch("CaloRecoJet_towersArea", &CaloRecoJet_towersArea);

      //JEC Uncertainity variable
      myEvent->Branch("CaloRecoJet_jecUncertainity", &CaloRecoJet_jecUncertainity);

      if(saveBTaggingInfo_CaloReco_){
	myEvent->Branch("CaloRecoJet_BJetDiscrByTrackCountingHighEff", &CaloRecoJet_BJetDiscrByTrackCountingHighEff);
	myEvent->Branch("CaloRecoJet_BJetDiscrByTrackCountingHighPur", &CaloRecoJet_BJetDiscrByTrackCountingHighPur);
	myEvent->Branch("CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff", &CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff);
	myEvent->Branch("CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur", &CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur);
	myEvent->Branch("CaloRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff", &CaloRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff);
	myEvent->Branch("CaloRecoJet_BJetDiscrByCombinedSecondaryVertexMVA", &CaloRecoJet_BJetDiscrByCombinedSecondaryVertexMVA);
	myEvent->Branch("CaloRecoJet_BJetDiscrByjetProbabilityBJetTags", &CaloRecoJet_BJetDiscrByjetProbabilityBJetTags);
	myEvent->Branch("CaloRecoJet_BJetDiscrByjetBProbabilityBJetTags", &CaloRecoJet_BJetDiscrByjetBProbabilityBJetTags);
	myEvent->Branch("CaloRecoJet_BJetDiscrBySoftElectronBJetTags", &CaloRecoJet_BJetDiscrBySoftElectronBJetTags);
	myEvent->Branch("CaloRecoJet_BJetDiscrBySoftMuonBJetTags", &CaloRecoJet_BJetDiscrBySoftMuonBJetTags);
      } 
    }

    if(runPFRecoJets_){

      myEvent->Branch("PFRecoJet_n", &PFRecoJet_n, "PFRecoJet_n/I");

      //Kinematic Variables
      myEvent->Branch("PFRecoJet_E", &PFRecoJet_E);
      myEvent->Branch("PFRecoJet_et", &PFRecoJet_et);
      myEvent->Branch("PFRecoJet_pt", &PFRecoJet_pt);
      myEvent->Branch("PFRecoJet_eta", &PFRecoJet_eta);
      myEvent->Branch("PFRecoJet_phi", &PFRecoJet_phi);
      myEvent->Branch("PFRecoJet_theta", &PFRecoJet_theta);
      myEvent->Branch("PFRecoJet_px", &PFRecoJet_px);
      myEvent->Branch("PFRecoJet_py", &PFRecoJet_py);
      myEvent->Branch("PFRecoJet_pz", &PFRecoJet_pz);
      myEvent->Branch("PFRecoJet_vx", &PFRecoJet_vx);
      myEvent->Branch("PFRecoJet_vy", &PFRecoJet_vy);
      myEvent->Branch("PFRecoJet_vz", &PFRecoJet_vz);

      //Id and other Variables for PF Reco Jets
      myEvent->Branch("PFRecoJet_ChargedEmEnergy", &PFRecoJet_ChargedEmEnergy);
      myEvent->Branch("PFRecoJet_ChargedEmEnergyFrac", &PFRecoJet_ChargedEmEnergyFrac);
      myEvent->Branch("PFRecoJet_ChargedHadEnergy", &PFRecoJet_ChargedHadEnergy);
      myEvent->Branch("PFRecoJet_ChargedHadEnergyFrac", &PFRecoJet_ChargedHadEnergyFrac);
      myEvent->Branch("PFRecoJet_ChargedHadMult", &PFRecoJet_ChargedHadMult);
      myEvent->Branch("PFRecoJet_ChargedMult", &PFRecoJet_ChargedMult);
      myEvent->Branch("PFRecoJet_NConstituents", &PFRecoJet_NConstituents);
      myEvent->Branch("PFRecoJet_HFEMEnergy", &PFRecoJet_HFEMEnergy);
      myEvent->Branch("PFRecoJet_HFEMEnergyFrac", &PFRecoJet_HFEMEnergyFrac);
      myEvent->Branch("PFRecoJet_HFEMMult", &PFRecoJet_HFEMMult);
      myEvent->Branch("PFRecoJet_HFHadEnergy", &PFRecoJet_HFHadEnergy);
      myEvent->Branch("PFRecoJet_HFHadEnergyFrac", &PFRecoJet_HFHadEnergyFrac);
      myEvent->Branch("PFRecoJet_HFHadMult", &PFRecoJet_HFHadMult);
      myEvent->Branch("PFRecoJet_NeutralEmEnergy", &PFRecoJet_NeutralEmEnergy);
      myEvent->Branch("PFRecoJet_NeutralEmEnergyFrac", &PFRecoJet_NeutralEmEnergyFrac);
      myEvent->Branch(" PFRecoJet_NeutralHadEnergy", & PFRecoJet_NeutralHadEnergy);
      myEvent->Branch("PFRecoJet_NeutralHadEnergyFrac", &PFRecoJet_NeutralHadEnergyFrac);
      myEvent->Branch("PFRecoJet_NeutralHadMult", &PFRecoJet_NeutralHadMult);
      myEvent->Branch("PFRecoJet_NeutralMult", &PFRecoJet_NeutralMult);

      //JEC Uncertainity variable
      myEvent->Branch("PFRecoJet_jecUncertainity", &PFRecoJet_jecUncertainity);

      if(saveBTaggingInfo_PFReco_){
	myEvent->Branch("PFRecoJet_BJetDiscrByTrackCountingHighEff", &PFRecoJet_BJetDiscrByTrackCountingHighEff);
	myEvent->Branch("PFRecoJet_BJetDiscrByTrackCountingHighPur", &PFRecoJet_BJetDiscrByTrackCountingHighPur);
	myEvent->Branch("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff);
	myEvent->Branch("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur);
	myEvent->Branch("PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff);
	myEvent->Branch("PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA);
	myEvent->Branch("PFRecoJet_BJetDiscrByjetProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetProbabilityBJetTags);
	myEvent->Branch("PFRecoJet_BJetDiscrByjetBProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetBProbabilityBJetTags);
	myEvent->Branch("PFRecoJet_BJetDiscrBySoftElectronBJetTags", &PFRecoJet_BJetDiscrBySoftElectronBJetTags);
	myEvent->Branch("PFRecoJet_BJetDiscrBySoftMuonBJetTags", &PFRecoJet_BJetDiscrBySoftMuonBJetTags);
      } 
    }
  } //END_OF if(runJets_)

  if(runVertex_){
    myEvent->Branch("Vertex_n", &Vertex_n, "Vertex_n/I");
    myEvent->Branch("Vertex_x", &Vertex_x);
    myEvent->Branch("Vertex_y", &Vertex_y);
    myEvent->Branch("Vertex_z", &Vertex_z);
    myEvent->Branch("Vertex_chi2", &Vertex_chi2);
    myEvent->Branch("Vertex_nchi2", &Vertex_nchi2);
    myEvent->Branch("Vertex_ndof", &Vertex_ndof);
    myEvent->Branch("Vertex_tracksSize", &Vertex_tracksSize);
    myEvent->Branch("Vertex_isFake", &Vertex_isFake);
    myEvent->Branch("Vertex_isValid", &Vertex_isValid);
    myEvent->Branch("Vertex_d0", &Vertex_d0);
  }

  if(runTracks_){
    myEvent->Branch("Tracks_n", &Tracks_n, "Tracks_n/I");
    myEvent->Branch("Track_pt", &Track_pt);
    myEvent->Branch("Track_px", &Track_px);
    myEvent->Branch("Track_py", &Track_py);
    myEvent->Branch("Track_pz", &Track_pz);
    myEvent->Branch("Track_vx", &Track_vx);
    myEvent->Branch("Track_vy", &Track_vy);
    myEvent->Branch("Track_vz", &Track_vz);
    myEvent->Branch("Track_eta", &Track_eta);
    myEvent->Branch("Track_phi", &Track_phi);
    myEvent->Branch("Track_theta", &Track_theta);
    myEvent->Branch("Track_chi2", &Track_chi2);
  }

  if(runScraping_){
    myEvent->Branch("IsScrapingEvent_", &IsScrapingEvent_, "IsScrapingEvent_/O");
    myEvent->Branch("Scraping_FractionOfGoodTracks", &Scraping_FractionOfGoodTracks, "Scraping_FractionOfGoodTracks/F");
  }

  if(runHLT_){
    myEvent->Branch("HLT_Photon_triggers", &HLT_Photon_triggers);
    myEvent->Branch("HLT_Photon_trig_prescales", &HLT_Photon_trig_prescales);
    myEvent->Branch("HLT_Photon_ifTriggerPassed", &HLT_Photon_ifTriggerPassed);
    myEvent->Branch("HLT_Photon_nTriggers", &HLT_Photon_nTriggers, "HLT_Photon_nTriggers/I");
    myEvent->Branch("HLT_Photon_triggerIndex", &HLT_Photon_triggerIndex);
    myEvent->Branch("HLT_Photon_nFilters", &HLT_Photon_nFilters);
//    myEvent->Branch("HLT_Photon_FilterNames", &HLT_Photon_FilterNames);
    myEvent->Branch("HLT_Photon_FilterObjects_pt", &HLT_Photon_FilterObjects_pt);
    myEvent->Branch("HLT_Photon_FilterObjects_eta", &HLT_Photon_FilterObjects_eta);
    myEvent->Branch("HLT_Photon_FilterObjects_phi", &HLT_Photon_FilterObjects_phi);

  }
}


// ------------ method called when starting to processes a run  ------------
void 
Analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  if(debug_){
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"BEGIN NEW RUN: "<< iRun.run() <<endl;
  }

  bool changed(true);
  if(hltConfig_.init(iRun,iSetup,triggerEventLabel_.process(),changed)){
    //        cout << "Table name = " << hltConfig_.tableName() << endl;
    if(changed){
    }
  }
}


// ------------ method called when starting to processes a luminosity block  ------------
void 
Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called for each event  ------------
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Photon_n = 0;
  CaloPatJet_n = 0;
  PFPatJet_n = 0;
  CaloRecoJet_n = 0;
  PFRecoJet_n = 0;
  Vertex_n = 0;
  Tracks_n = 0;

  nEvents++;
 
  ClearVectors();

  if(debug_) {}cout << "PROCESSING EVENT = " << nEvents << endl;

  RunNumber   = iEvent.id().run();
  EventNumber = iEvent.id().event();
  LumiNumber  = iEvent.id().luminosityBlock();
  BXNumber = iEvent.bunchCrossing();
  isRealData = iEvent.isRealData();

  if(debug_) cout << "FOUND BUG IN: RunNo. = " << RunNumber << ", EventNo. = " << EventNumber << endl;


  //////////////////////////
  //                      //
  //  get LuminosityInfo  //
  //                      //
  //////////////////////////

  if(saveLuminosityInfo_){
    //get ConditionsInLumiBlock
    const edm::LuminosityBlock& iLumi = iEvent.getLuminosityBlock();
    edm::Handle<edm::ConditionsInLumiBlock> condInLumiBlock;
    iLumi.getByLabel("conditionsInEdm", condInLumiBlock);

    //initialize to -1 in case isValid fails
    totalIntensityBeam1 = -1;
    totalIntensityBeam2 = -1;
    if(condInLumiBlock.isValid()){//check pointer condInLumiBlock is not null using general isValid() fun for pointers
      totalIntensityBeam1 = condInLumiBlock->totalIntensityBeam1;
      totalIntensityBeam2 = condInLumiBlock->totalIntensityBeam2;
    }
  
    // get LumiSummary
    edm::Handle<LumiSummary> lumiSummary;
    iLumi.getByLabel("lumiProducer", lumiSummary);

    //initialize to -1 in case isValid fails
    avgInsDelLumi = -1.;
    avgInsDelLumiErr = -1.;
    avgInsRecLumi = -1.;
    avgInsRecLumiErr = -1.;
  
    if(lumiSummary.isValid()){ //check pointer lumiSummary is not null using general isValid() fun for all pointers accessible by a dot (.) always//
      if(lumiSummary->isValid()){ //check if data is valid using isValid() fun defined in LumiSummary.h and is accessible by -> like all other funs//
      
	avgInsDelLumi = lumiSummary->avgInsDelLumi();
	avgInsDelLumiErr = lumiSummary->avgInsDelLumiErr();    //data is valid only if run exists from all sources lumi, trg, hlt//
	avgInsRecLumi = lumiSummary->avgInsRecLumi();
	avgInsRecLumiErr = lumiSummary->avgInsRecLumiErr();

      }
    }
  }

  if(debug_) cout << "FOUND BUG AFTER LUMINOSITY INFO" << endl;


  /////////////////////////////////
  //                             //
  // get PileUp Info For MC Only //
  //                             //
  /////////////////////////////////

  if(runMCPileUp_){
    //pile up info from MC only
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(MCpileupLabel_, PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    npuVertices = 0;
    npuVerticesm1 = 0;
    npuVerticesp1 = 0;
    npuVerticespm2 = 0;
    trueNumofInteractions = 0.0;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI){
      trueNumofInteractions = PVI->getTrueNumInteractions();
      if(PVI->getBunchCrossing()==  0){ npuVertices += PVI->getPU_NumInteractions(); }
      if(PVI->getBunchCrossing()== -1){ npuVerticesm1 += PVI->getPU_NumInteractions(); }
      if(PVI->getBunchCrossing()==  1){ npuVerticesp1 += PVI->getPU_NumInteractions(); }
      if(abs(PVI->getBunchCrossing()) >= 2){ npuVerticespm2 += PVI->getPU_NumInteractions(); }
      if(debug_) cout << "PILEUP INFORMATION: BX, nVtx,: " << PVI->getBunchCrossing()               
                      << " ," << PVI->getPU_NumInteractions() << endl;   
    }
  }


  ////////////////////////////
  //                        //
  // Get Photon Information //
  //                        //
  ////////////////////////////

  if(runPhotons_){

    //Accessing pat::Photon Collection
    edm::Handle<pat::PhotonCollection> patPhoColl;
    iEvent.getByLabel(patPhoLabel_, patPhoColl);
    pat::PhotonCollection::const_iterator patPho_itr;

    Photon_n = patPhoColl->size();

    if(debug_) cout << "No. of Photons for event No. = " << EventNumber << " is " << Photon_n << endl;

    //---------------------------------------------------------------   
    //Accessing various collections required for different variables
    //---------------------------------------------------------------

    //Barrel rechits collection (Required for showerShapes, crystalInfo, roundness, angle, swisscross ....)
    edm::Handle<EcalRecHitCollection> Brechit;
    iEvent.getByLabel(BarrelrechitLabel_, Brechit);
    const EcalRecHitCollection* BarrelRecHits = Brechit.product();

    //Endcap rechits collection (Required for showerShapes, crystalInfo, roundness, angle, swisscross ....)
    edm::Handle<EcalRecHitCollection> Erechit;
    iEvent.getByLabel(EndcaprechitLabel_,Erechit);
    const EcalRecHitCollection* EndcapRecHits = Erechit.product();

    //EcalSeverityLevelAlgoRecords (Required for showerShapes, crystalInfo, roundness, angle, swisscross ....)
    edm::ESHandle<EcalSeverityLevelAlgo> sevlv;                                      
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
    const EcalSeverityLevelAlgo* EcalSevAlgo = sevlv.product();

    //CaloTopology records (Required for showerShapes, crystalInfo, roundness, angle, swisscross ....)
    edm::ESHandle<CaloTopology> theCaloTopo;
    iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
    const CaloTopology* CaloTopology = theCaloTopo.product();

    //CaloGeometry records (Required for showerShapes, crystalInfo, roundness, angle, swisscross ....)
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloGeometry* CaloGeometry = geoHandle.product();	

    //BeamSpot Handle (Required for Electron Veto)
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel(BSLabel_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();

    //Photon Conversion Collection (Required for Electron Veto)
    edm::Handle<reco::ConversionCollection> ConversionColl;
    iEvent.getByLabel(ConvPhoLabel_, ConversionColl);

    //GsfElectron Collection (Required for Electron Veto)
    edm::Handle<reco::GsfElectronCollection> GsfEleColl;
    iEvent.getByLabel(GsfEleLabel_, GsfEleColl);

    if(debug_) cout << "Pat Photon Loop Starts" << endl;

    //Pat::Photon loop starts
    for(patPho_itr = patPhoColl->begin(); patPho_itr != patPhoColl->end(); patPho_itr++){
      const pat::Photon* patPho_ptr = &(*patPho_itr);

      //Kinematic Variables
      Photon_E.push_back(patPho_ptr->energy());
      Photon_et.push_back(patPho_ptr->et());
      Photon_pt.push_back(patPho_ptr->pt());
      Photon_eta.push_back(patPho_ptr->eta());
      Photon_phi.push_back(correct_Phi(patPho_ptr->phi()));
      Photon_theta.push_back(patPho_ptr->theta());
      Photon_px.push_back(patPho_ptr->px());
      Photon_py.push_back(patPho_ptr->py());
      Photon_pz.push_back(patPho_ptr->pz());
      Photon_vx.push_back(patPho_ptr->vx());
      Photon_vy.push_back(patPho_ptr->vy());
      Photon_vz.push_back(patPho_ptr->vz());

      //Shower Shape Variables
      Photon_r9.push_back(patPho_ptr->r9());
      Photon_maxEnergyXtal.push_back(patPho_ptr->maxEnergyXtal());
      Photon_e1x5.push_back(patPho_ptr->e1x5());
      Photon_e2x5.push_back(patPho_ptr->e2x5());
      Photon_e3x3.push_back(patPho_ptr->e3x3());
      Photon_e5x5.push_back(patPho_ptr->e5x5());
      Photon_r1x5.push_back(patPho_ptr->r1x5());
      Photon_r2x5.push_back(patPho_ptr->r2x5());
      Photon_SigmaEtaEta.push_back(patPho_ptr->sigmaEtaEta());
      Photon_SigmaIEtaIEta.push_back(patPho_ptr->sigmaIetaIeta());

      //Getting SigmaIEtaIphi, SigmaIphiIphi, roundness, angle, swissCross etc from EcalClusterTools
      //Taking proper rechits depending upon whether photon belongs to EB or EE     
      const EcalRecHitCollection* EcalRecHits =  ( patPho_ptr->isEB() ? BarrelRecHits : EndcapRecHits );
    
      const reco::BasicCluster& seedCluster = *(patPho_ptr->superCluster()->seed());
      //Accessing Id of seed crystal (required for swissCross)     
      const DetId& Id_max = EcalClusterTools::getMaximum(seedCluster, EcalRecHits).first;
      float ene_maxId = recHitE(Id_max, *(EcalRecHits));
      //Getting sigmaIEtaIPhi etc.
      std::vector<float> phoCov = EcalClusterTools::covariances(seedCluster, EcalRecHits, CaloTopology, CaloGeometry, 
                                                                flagsexclEB_,severitiesexclEB_, EcalSevAlgo, 4.7);
      std::vector<float> phoLocalCov = EcalClusterTools::localCovariances(seedCluster, EcalRecHits, CaloTopology, 
                                                                          flagsexclEB_, severitiesexclEB_, EcalSevAlgo, 4.7);
      Photon_SigmaEtaPhi.push_back(phoCov.at(1));
      Photon_SigmaIEtaIPhi.push_back(phoLocalCov.at(1));
      Photon_SigmaPhiPhi.push_back(phoCov.at(2));
      Photon_SigmaIPhiIPhi.push_back(phoLocalCov.at(2));
      
      //Getting roundness and angle
      vector<float> roundness_and_angle;
      if(patPho_ptr->isEB()){
        roundness_and_angle = EcalClusterTools::roundnessBarrelSuperClusters(*(patPho_ptr->superCluster()), *(EcalRecHits), 1, 0);
	Photon_roundness.push_back(roundness_and_angle.at(0));
	Photon_angle.push_back(roundness_and_angle.at(1));
      }else{
	Photon_roundness.push_back(-99); //For endcap we put some -ve values as the fun. computes roundness and angle for barrel only 
	Photon_angle.push_back(-99);     //As canbe seen in the fun. code, it uses eta-phi for xtals while in endcap xtals are numbered in x-y.
      }

      //Getting swissCross
      float swissCross = -99;
      //checking if the Id_max belong to the rechitCollection
      EcalRecHitCollection::const_iterator itr = EcalRecHits->find(Id_max);
      if(itr != EcalRecHits->end()) swissCross = EcalTools::swissCross(Id_max, *(EcalRecHits), 0.08, true);;
      Photon_swissCross.push_back(swissCross);
      
      if(debug_ && swissCross > 0.95) cout << "This photon candidate is an ECAL spike identified by Swiss Cross algorithm" << endl;

      //getting s9 which is the ratio of seed cry energy and the energy in 3x3.
      float s9 = ene_maxId/(patPho_ptr->e3x3());
      Photon_s9.push_back(s9);

      //calling methods to fill e4e1, d6e2, e2e9, rookFraction.
      float e4e1 = GetE4OverE1(Id_max, *(EcalRecHits));
      float e6e2 = GetE6OverE2(Id_max, *(EcalRecHits));
      float e2e9 = GetE2OverE9(Id_max, *(EcalRecHits));
      float rookFrac = GetrookFraction(Id_max, *(EcalRecHits));

      Photon_e4Overe1.push_back(e4e1);
      Photon_e6Overe2.push_back(e6e2);
      Photon_e2Overe9.push_back(e2e9);
      Photon_rookFraction.push_back(rookFrac);

      if(debug_ && e6e2 < 0.04) cout << "Double Spikes in RECO" << endl;
      if(debug_ && rookFrac < 0.1) cout << " This is Noise, not an actual photon" << endl;

      //Fiducial Flags Variables
      Photon_isEB.push_back(patPho_ptr->isEB());
      Photon_isEE.push_back(patPho_ptr->isEE());
      Photon_isEBGap.push_back(patPho_ptr->isEBGap());
      Photon_isEEGap.push_back(patPho_ptr->isEEGap());
      Photon_isEBEEGap.push_back(patPho_ptr->isEBEEGap());

      //Detector Isolation Variables
      Photon_ecalRecHitSumEtConeDR03.push_back(patPho_ptr->ecalRecHitSumEtConeDR03());
      Photon_hcalTowerSumEtConeDR03.push_back(patPho_ptr->hcalTowerSumEtConeDR03());
      Photon_hcalDepth1TowerSumEtConeDR03.push_back(patPho_ptr->hcalDepth1TowerSumEtConeDR03());
      Photon_hcalDepth2TowerSumEtConeDR03.push_back(patPho_ptr->hcalDepth2TowerSumEtConeDR03());
      Photon_trkSumPtSolidConeDR03.push_back(patPho_ptr->trkSumPtSolidConeDR03());
      Photon_trkSumPtHollowConeDR03.push_back(patPho_ptr->trkSumPtHollowConeDR03());
      Photon_nTrkSolidConeDR03.push_back(patPho_ptr->nTrkSolidConeDR03());
      Photon_nTrkHollowConeDR03.push_back(patPho_ptr->nTrkHollowConeDR03());
      Photon_ecalRecHitSumEtConeDR04.push_back(patPho_ptr->ecalRecHitSumEtConeDR04());
      Photon_hcalTowerSumEtConeDR04.push_back(patPho_ptr->hcalTowerSumEtConeDR04());
      Photon_hcalDepth1TowerSumEtConeDR04.push_back(patPho_ptr->hcalDepth1TowerSumEtConeDR04());
      Photon_hcalDepth2TowerSumEtConeDR04.push_back(patPho_ptr->hcalDepth2TowerSumEtConeDR04());
      Photon_trkSumPtSolidConeDR04.push_back(patPho_ptr->trkSumPtSolidConeDR04());
      Photon_trkSumPtHollowConeDR04.push_back(patPho_ptr->trkSumPtHollowConeDR04());
      Photon_nTrkSolidConeDR04.push_back(patPho_ptr->nTrkSolidConeDR04());
      Photon_nTrkHollowConeDR04.push_back(patPho_ptr->nTrkHollowConeDR04());

      //Photon Identification Variables
      Photon_HoE.push_back(patPho_ptr->hadronicOverEm());
      Photon_SingleTowerHoE.push_back(patPho_ptr->hadTowOverEm());
      Photon_hasConvTrk.push_back(patPho_ptr->hasConversionTracks());
      Photon_hasPixelSeed.push_back(patPho_ptr->hasPixelSeed());
      passedConvSafeElectronVeto.push_back(!ConversionTools::hasMatchedPromptElectron(patPho_ptr->superCluster(), GsfEleColl,  
                                                                                      ConversionColl, beamspot.position()));//return true if a photon.


      //Photon Super Cluster Variables
      Photon_SC_nOfBasicClusters.push_back(patPho_ptr->superCluster()->clustersSize());
      Photon_SC_rawEnergy.push_back(patPho_ptr->superCluster()->rawEnergy());
      Photon_SC_preShowerEnergy.push_back(patPho_ptr->superCluster()->preshowerEnergy());
      Photon_SC_energy.push_back(patPho_ptr->superCluster()->energy());
      Photon_SC_eta.push_back(patPho_ptr->superCluster()->eta());
      Photon_SC_phi.push_back(correct_Phi(patPho_ptr->superCluster()->phi()));
      Photon_SC_x.push_back(patPho_ptr->superCluster()->x());
      Photon_SC_y.push_back(patPho_ptr->superCluster()->y());
      Photon_SC_z.push_back(patPho_ptr->superCluster()->z());
      Photon_SC_etaWidth.push_back(patPho_ptr->superCluster()->etaWidth());
      Photon_SC_phiWidth.push_back(patPho_ptr->superCluster()->phiWidth());

      if(debug_) cout << " The SuperCluster energy for this photon is " << patPho_ptr->superCluster()->energy() << endl;

      //Photon MIP Variables
      Photon_mipChi2.push_back(patPho_ptr->mipChi2());
      Photon_mipTotEnergy.push_back(patPho_ptr->mipTotEnergy());
      Photon_mipSlope.push_back(patPho_ptr->mipSlope());
      Photon_mipIntercept.push_back(patPho_ptr->mipIntercept());
      Photon_mipNhitCone.push_back(patPho_ptr->mipNhitCone());
      Photon_mipIsHalo.push_back(patPho_ptr->mipIsHalo());

      //Accessing Conversions Collection 
      reco::ConversionRefVector ConversionV = patPho_ptr->conversions();
      for(unsigned int iConv=0; iConv<ConversionV.size(); iConv++){
	reco::ConversionRef Conv_ref = ConversionV.at(iConv);
	if(Conv_ref->nTracks() < 2) continue;
	if(Conv_ref->conversionVertex().isValid()){

	  //Converted Photon Variables
	  Photon_nConvTracks.push_back(Conv_ref->nTracks());
	  Photon_isConverted.push_back(Conv_ref->isConverted());
	  Photon_pairInvariantMass.push_back(Conv_ref->pairInvariantMass());
	  Photon_pairCotThetaSeparation.push_back(Conv_ref->pairCotThetaSeparation());
	  Photon_pairMomentum_x.push_back(Conv_ref->pairMomentum().x());
	  Photon_pairMomentum_y.push_back(Conv_ref->pairMomentum().y());
	  Photon_pairMomentum_z.push_back(Conv_ref->pairMomentum().z());
	  Photon_EoverP.push_back(Conv_ref->EoverP());
	  Photon_conv_vx.push_back(Conv_ref->conversionVertex().x());
	  Photon_conv_vy.push_back(Conv_ref->conversionVertex().y());
	  Photon_conv_vz.push_back(Conv_ref->conversionVertex().z());
	  Photon_zOfPrimaryVtxFromTrks.push_back(Conv_ref->zOfPrimaryVertexFromTracks());
	  Photon_distOfMinimumApproach.push_back(Conv_ref->distOfMinimumApproach());
	  Photon_dPhiTracksAtVtx.push_back(Conv_ref->dPhiTracksAtVtx());
	  Photon_dPhiTracksAtEcal.push_back(Conv_ref->dPhiTracksAtEcal());
	  Photon_dEtaTracksAtEcal.push_back(Conv_ref->dEtaTracksAtEcal());
	}
      }

      //-----------------------------------
      //Saving crystal info for each photon
      //-----------------------------------
    
      if(savePhotonCrystalInfo_){

	vector< std::pair<DetId, float> > PhotonDetIds_and_Hits = patPho_ptr->superCluster()->hitsAndFractions();
	vector<CrystalInfo> CrystalInfo_container;
	CrystalInfo_container.clear();

	CrystalInfo xtalInfo;
	float total_time = 0.0;
	float timing_avg = 0.0;
	int ncrystals = 0;

	vector< std::pair<DetId, float> >::const_iterator det_itr;

	for(det_itr = PhotonDetIds_and_Hits.begin(); det_itr != PhotonDetIds_and_Hits.end(); det_itr++){
	if (((*det_itr).first).det() == DetId::Ecal && ((*det_itr).first).subdetId() == EcalBarrel){

	  EcalRecHitCollection::const_iterator j = BarrelRecHits->find(((*det_itr).first));
	  EcalRecHitCollection::const_iterator thishit;
	  if(j != BarrelRecHits->end()) thishit = j;
	  else continue;

	  EBDetId detId  = (EBDetId)((*det_itr).first);
	  xtalInfo.rawId = thishit->id().rawId();
	  xtalInfo.energy = thishit->energy();
	  xtalInfo.time = thishit->time();
	  xtalInfo.timeErr = thishit->timeError();
	  xtalInfo.recoFlag = thishit->recoFlag();
	  xtalInfo.ieta = detId.ieta();
	  xtalInfo.iphi = detId.iphi();
	  if(xtalInfo.energy > 0.1){
	    total_time = total_time + xtalInfo.time;
	    ncrystals++;
	  }
	}
	else if (((*det_itr).first).det() == DetId::Ecal && ((*det_itr).first).subdetId() == EcalEndcap){

	  EcalRecHitCollection::const_iterator j = EndcapRecHits->find(((*det_itr).first));
	  EcalRecHitCollection::const_iterator thishit;
	  if(j != EndcapRecHits->end()) thishit = j;
	  else continue;

	  EEDetId detId  = (EEDetId)((*det_itr).first);
	  xtalInfo.rawId = 999;
	  xtalInfo.energy = thishit->energy();
	  xtalInfo.time = thishit->time();
	  xtalInfo.timeErr = thishit->timeError();
	  xtalInfo.recoFlag = thishit->recoFlag();
	  xtalInfo.ieta = -99;
	  xtalInfo.iphi = -99;
	  if(xtalInfo.energy > 0.1){
	    total_time = total_time + xtalInfo.time;
	    ncrystals++;
	  }
	}

	CrystalInfo_container.push_back(xtalInfo);
	}

	std::sort(CrystalInfo_container.begin(), CrystalInfo_container.end(), EnergySortCriterium());

	if(ncrystals != 0) timing_avg = total_time/ncrystals;
	else timing_avg = -99.;

	Photon_nCrystals.push_back(CrystalInfo_container.size());
	Photon_avgTimeAllxtals.push_back(timing_avg);
	
	// clearing temporary vector used for filling 2D vectors
	temp_Photon_xtal_timing.clear();
	temp_Photon_xtal_timeErr.clear();
	temp_Photon_xtal_energy.clear();
	temp_Photon_xtal_EBieta.clear();
	temp_Photon_xtal_EBiphi.clear();
	temp_Photon_xtal_EBrecoFlag.clear();

	for(unsigned int y=0; y<CrystalInfo_container.size(); y++){
	  temp_Photon_xtal_timing.push_back(CrystalInfo_container[y].time);
	  temp_Photon_xtal_timeErr.push_back(CrystalInfo_container[y].timeErr);
	  temp_Photon_xtal_energy.push_back(CrystalInfo_container[y].energy);
	  temp_Photon_xtal_EBieta.push_back(CrystalInfo_container[y].ieta);
	  temp_Photon_xtal_EBiphi.push_back(CrystalInfo_container[y].iphi);
	  temp_Photon_xtal_EBrecoFlag.push_back(CrystalInfo_container[y].recoFlag);
	}
	Photon_xtal_timing.push_back(temp_Photon_xtal_timing);
	Photon_xtal_timeErr.push_back(temp_Photon_xtal_timeErr);
	Photon_xtal_energy.push_back(temp_Photon_xtal_energy);
	Photon_xtal_EBieta.push_back(temp_Photon_xtal_EBieta);
	Photon_xtal_EBiphi.push_back(temp_Photon_xtal_EBiphi);
	Photon_xtal_EBrecoFlag.push_back(temp_Photon_xtal_EBrecoFlag);
      
      } //END_OF if(savePhotonCrystalInfo_)

    } // END_OF for(patPho_itr = patPhoColl->begin(); patPho_itr != patPhoColl->end(); patPho_itr++)


    //-------------------------------
    // PF ISOLATION PROCEDURE STARTS
    //------------------------------- 

    if(savePhotonPFIsolation_){
      edm::Handle<reco::PhotonCollection> recoPhoColl;
      bool found = iEvent.getByLabel(recoPhoLabel_, recoPhoColl);
      reco::PhotonCollection::const_iterator recoPho_itr;

      if(!found ) {
	std::ostringstream  err;
        err<<" cannot get Photons: "
	   <<recoPhoLabel_<<std::endl;
	edm::LogError("Analyzer")<<err.str();
	throw cms::Exception( "MissingProduct", err.str());
      }

      std::vector<reco::Photon> recoPhoton_container;
      recoPhoton_container.clear();

      // To pt sort reco::Photons in same order as that of pat::Photons which are pt sorted by construction
      for(patPho_itr = patPhoColl->begin(); patPho_itr != patPhoColl->end(); ++patPho_itr){
	const pat::Photon* patPho_ptr = &(*patPho_itr);

	for(recoPho_itr = recoPhoColl->begin(); recoPho_itr !=recoPhoColl->end(); ++recoPho_itr){
	  const reco::Photon* recoPho_ptr = &(*recoPho_itr);
     
	  if(patPho_ptr->et() != recoPho_ptr->et()) continue;

	  recoPhoton_container.push_back(*recoPho_itr);
	}
      }

      //Accessing various collections required for PFIsolation

      // PFCanndidate Collection
      edm::Handle<reco::PFCandidateCollection> pfCandidatesH;
      iEvent.getByLabel(pfCandidateLabel_, pfCandidatesH);
      const reco::PFCandidateCollection thePfColl = *(pfCandidatesH.product());

      // Vertex Collection
      edm::Handle<reco::VertexCollection>  vertexCollection;
      iEvent.getByLabel(vertexLabel_, vertexCollection);

      // Reference to primary vertex or vertex at position 0
      unsigned int ivtx = 0;
      reco::VertexRef myVtxRef(vertexCollection, ivtx);

      //reco::Photon loop starts
      for(unsigned int x = 0; x<recoPhoton_container.size(); x++){
	const reco::Photon* recoPho_ptr = &(recoPhoton_container.at(x));

	phoIsolator03.fGetIsolation(recoPho_ptr, &thePfColl, myVtxRef, vertexCollection); // calling this function will set all the isolation 
                                                                                          // values in corresponding vectors
	PFIsoPhoton03.push_back(phoIsolator03.getIsolationPhoton());
	PFIsoNeutral03.push_back(phoIsolator03.getIsolationNeutral());
	PFIsoCharged03.push_back(phoIsolator03.getIsolationCharged());
	PFIsoSum03.push_back(phoIsolator03.fGetIsolation(recoPho_ptr, &thePfColl, myVtxRef, vertexCollection));

	if(debug_){
	  cout << "PFIsoPhoton03[" << x << "] = " << PFIsoPhoton03.at(x) << endl;
	  cout << "PFIsoNeutral03[" << x << "] = " << PFIsoNeutral03.at(x) << endl;
	  cout << "PFIsoCharged03[" << x << "] = " << PFIsoCharged03.at(x) << endl;
	  cout << "PFIsoSum03[" << x << "] = " << PFIsoSum03.at(x) << endl;
	}

	//To get the worst charged hadron isolation corresponding to a non primary vertex 
	PFIsoChargedWorstvtx03.insert(PFIsoChargedWorstvtx03.begin()+x, 0); //Inserting zero at the xth position of vector.
    
	for(unsigned int iv=0; iv<vertexCollection->size(); iv++){
	  reco::VertexRef thisVtxRef(vertexCollection, iv);
	  phoIsolator03.fGetIsolation(recoPho_ptr, &thePfColl, thisVtxRef, vertexCollection);
	  double thisChargedHadronIso = phoIsolator03.getIsolationCharged();
	  if(thisChargedHadronIso > PFIsoChargedWorstvtx03.at(x)) PFIsoChargedWorstvtx03.at(x) = thisChargedHadronIso;
	}
      } //END_OF for(unsigned int x = 0; x<recoPhoton_container.size(); x++)
    } //END_OF if(savePhotonPFIsolation_) 
  
    if(savePhotonRecHitsInfo_){
 





    } //END_OF if(savePhotonRecHitsInfo_)
  } // END_OF if(runPhotons_)
  

  ////////////////////////////
  //                        //
  //  Get Jets Information  //
  //                        //
  ////////////////////////////

  if(runJets_){

    if(debug_) cout << "INSIDE runJets_ LOOP" << endl;

    if(runCaloPatJets_){

      //Accessing pat::CaloJets Collection
      edm::Handle<pat::JetCollection> CaloPatJetColl;
      iEvent.getByLabel(CaloPatJetLabel_, CaloPatJetColl);
      pat::JetCollection::const_iterator CaloPatJet_itr; 

      //Accessing Jec Uncertainities
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5Calo",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

      for(CaloPatJet_itr = CaloPatJetColl->begin(); CaloPatJet_itr != CaloPatJetColl->end(); CaloPatJet_itr++){
	const pat::Jet* CaloPatJet_ptr = &(*CaloPatJet_itr);
	if(CaloPatJet_ptr->pt() > 30){
	  CaloPatJet_n++;

	  if(debug_) cout << "Calo Pat Jet Pt is = " << CaloPatJet_ptr->pt() << endl;

	  //Kinematic Variables
	  CaloPatJet_E.push_back(CaloPatJet_ptr->energy());
	  CaloPatJet_et.push_back(CaloPatJet_ptr->et());
	  CaloPatJet_pt.push_back(CaloPatJet_ptr->pt());
	  CaloPatJet_eta.push_back(CaloPatJet_ptr->eta());
	  CaloPatJet_phi.push_back(correct_Phi(CaloPatJet_ptr->phi()));
	  CaloPatJet_theta.push_back(CaloPatJet_ptr->theta());
	  CaloPatJet_px.push_back(CaloPatJet_ptr->px());
	  CaloPatJet_py.push_back(CaloPatJet_ptr->py());
	  CaloPatJet_pz.push_back(CaloPatJet_ptr->pz());
	  CaloPatJet_vx.push_back(CaloPatJet_ptr->vx());
	  CaloPatJet_vy.push_back(CaloPatJet_ptr->vy());
	  CaloPatJet_vz.push_back(CaloPatJet_ptr->vz());
	
	  //Id and other Variables for Calo Pat Jets
	  CaloPatJet_emEnergyFraction.push_back(CaloPatJet_ptr->emEnergyFraction());
	  CaloPatJet_emEnergyInEB.push_back(CaloPatJet_ptr->emEnergyInEB());
	  CaloPatJet_emEnergyInEE.push_back(CaloPatJet_ptr->emEnergyInEE());
	  CaloPatJet_emEnergyInHF.push_back(CaloPatJet_ptr->emEnergyInHF());
	  CaloPatJet_HadronicEnergyFraction.push_back(CaloPatJet_ptr->energyFractionHadronic());
	  CaloPatJet_NConstituents.push_back(CaloPatJet_ptr->getCaloConstituents().size());
	  CaloPatJet_HadronicEnergyInHB.push_back(CaloPatJet_ptr->hadEnergyInHB());
	  CaloPatJet_HadronicEnergyInHE.push_back(CaloPatJet_ptr->hadEnergyInHE());
	  CaloPatJet_HadronicEnergyInHF.push_back(CaloPatJet_ptr->hadEnergyInHF());
	  CaloPatJet_HadronicEnergyInHO.push_back(CaloPatJet_ptr->hadEnergyInHO());
	  CaloPatJet_maxEnergyInEmTowers.push_back(CaloPatJet_ptr->maxEInEmTowers());
	  CaloPatJet_maxEnergyInHadTowers.push_back(CaloPatJet_ptr->maxEInHadTowers());
	  CaloPatJet_nConstituents60E.push_back(CaloPatJet_ptr->n60());
	  CaloPatJet_nConstituents90E.push_back(CaloPatJet_ptr->n90());
	  CaloPatJet_towersArea.push_back(CaloPatJet_ptr->towersArea());
	  CaloPatJet_nTowers.push_back(CaloPatJet_ptr->jetID().nECALTowers + CaloPatJet_ptr->jetID().nHCALTowers);
	  CaloPatJet_fEinHottestHPD.push_back(CaloPatJet_ptr->jetID().fHPD);
	  CaloPatJet_fEinHottestRBX.push_back(CaloPatJet_ptr->jetID().fRBX);
	  CaloPatJet_nHitsCarrying90E.push_back(CaloPatJet_ptr->jetID().n90Hits);
	  CaloPatJet_nHitsinTowersCarrying90E.push_back(CaloPatJet_ptr->jetID().hitsInN90);
	  CaloPatJet_RHF.push_back((CaloPatJet_ptr->jetID().fLong - CaloPatJet_ptr->jetID().fShort)/(CaloPatJet_ptr->jetID().fLong + CaloPatJet_ptr->jetID().fShort));

	  //Jec Uncertainities 
	  jecUnc->setJetEta(CaloPatJet_ptr->eta());
	  jecUnc->setJetPt(CaloPatJet_ptr->pt());
	  double unc = jecUnc->getUncertainty(true);
	  CaloPatJet_jecUncertainity.push_back(unc);
	    
	  if(saveBTaggingInfo_CaloPat_){

	    CaloPatJet_BJetDiscrByTrackCountingHighEff.push_back(CaloPatJet_ptr->bDiscriminator("trackCountingHighEffBJetTags"));
	    CaloPatJet_BJetDiscrByTrackCountingHighPur.push_back(CaloPatJet_ptr->bDiscriminator("trackCountingHighPurBJetTags"));
	    CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighEff.push_back(CaloPatJet_ptr->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	    CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighPur.push_back(CaloPatJet_ptr->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	    CaloPatJet_BJetDiscrByCombinedSecondaryVertexHighEff.push_back(CaloPatJet_ptr->bDiscriminator("combinedSecondaryVertexBJetTags"));
	    CaloPatJet_BJetDiscrByCombinedSecondaryVertexMVA.push_back(CaloPatJet_ptr->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	    CaloPatJet_BJetDiscrByjetProbabilityBJetTags.push_back(CaloPatJet_ptr->bDiscriminator("jetProbabilityBJetTags"));
	    CaloPatJet_BJetDiscrByjetBProbabilityBJetTags.push_back(CaloPatJet_ptr->bDiscriminator("jetBProbabilityBJetTags"));
	    CaloPatJet_BJetDiscrBySoftElectronBJetTags.push_back(CaloPatJet_ptr->bDiscriminator("softElectronBJetTags"));
	    CaloPatJet_BJetDiscrBySoftMuonBJetTags.push_back(CaloPatJet_ptr->bDiscriminator("softMuonBJetTags"));
	    CaloPatJet_JetPartonFlavor.push_back(fabs(CaloPatJet_ptr->partonFlavour()));

	  } // END_OF if(saveBTaggingInfo_CaloPat_)	  
	} // END_OF if(CaloPatJet_ptr->pt() > 30)
      } // END_OF for(CaloPatJet_itr = CaloPatJetColl->begin(); CaloPatJet_itr != CaloPatJetColl->end(); CaloPatJet_itr++)

      delete jecUnc;

      if(debug_) cout << "End of CaloPatJet Loop" << endl;

    } //END_OF if(runCaloPatJets_)

    if(runPFPatJets_){

      //Accessing pat::PFJets Collection
      edm::Handle<pat::JetCollection> PFPatJetColl;
      iEvent.getByLabel(PFPatJetLabel_, PFPatJetColl);
      pat::JetCollection::const_iterator PFPatJet_itr; 

      //Accessing Jec Uncertainities
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

      for(PFPatJet_itr = PFPatJetColl->begin(); PFPatJet_itr != PFPatJetColl->end(); PFPatJet_itr++){
	const pat::Jet* PFPatJet_ptr = &(*PFPatJet_itr);
	if(PFPatJet_ptr->pt() > 30){
	  PFPatJet_n++;

	  if(debug_) cout << "PF Pat Jet Pt is = " << PFPatJet_ptr->pt() << endl;

	  //Kinematic Variables
	  PFPatJet_E.push_back(PFPatJet_ptr->energy());
	  PFPatJet_et.push_back(PFPatJet_ptr->et());
	  PFPatJet_pt.push_back(PFPatJet_ptr->pt());
	  PFPatJet_eta.push_back(PFPatJet_ptr->eta());
	  PFPatJet_phi.push_back(correct_Phi(PFPatJet_ptr->phi()));
	  PFPatJet_theta.push_back(PFPatJet_ptr->theta());
	  PFPatJet_px.push_back(PFPatJet_ptr->px());
	  PFPatJet_py.push_back(PFPatJet_ptr->py());
	  PFPatJet_pz.push_back(PFPatJet_ptr->pz());
	  PFPatJet_vx.push_back(PFPatJet_ptr->vx());
	  PFPatJet_vy.push_back(PFPatJet_ptr->vy());
	  PFPatJet_vz.push_back(PFPatJet_ptr->vz());
	 
	  //Id and other Variables for PF Pat Jets
	  PFPatJet_ChargedEmEnergy.push_back(PFPatJet_ptr->chargedEmEnergy());
	  PFPatJet_ChargedEmEnergyFrac.push_back(PFPatJet_ptr->chargedEmEnergyFraction());
	  PFPatJet_ChargedHadEnergy.push_back(PFPatJet_ptr->chargedHadronEnergy());
	  PFPatJet_ChargedHadEnergyFrac.push_back(PFPatJet_ptr->chargedHadronEnergyFraction());
	  PFPatJet_ChargedHadMult.push_back(PFPatJet_ptr->chargedHadronMultiplicity());
	  PFPatJet_ChargedMult.push_back(PFPatJet_ptr->chargedMultiplicity());
	  PFPatJet_NConstituents.push_back(PFPatJet_ptr->getPFConstituents().size());
	  PFPatJet_HFEMEnergy.push_back(PFPatJet_ptr->HFEMEnergy());
	  PFPatJet_HFEMEnergyFrac.push_back(PFPatJet_ptr->HFEMEnergyFraction());
	  PFPatJet_HFEMMult.push_back(PFPatJet_ptr->HFEMMultiplicity());
	  PFPatJet_HFHadEnergy.push_back(PFPatJet_ptr->HFHadronEnergy());
	  PFPatJet_HFHadEnergyFrac.push_back(PFPatJet_ptr->HFHadronEnergyFraction());
	  PFPatJet_HFHadMult.push_back(PFPatJet_ptr->HFHadronMultiplicity());
	  PFPatJet_NeutralEmEnergy.push_back(PFPatJet_ptr->neutralEmEnergy());
	  PFPatJet_NeutralEmEnergyFrac.push_back(PFPatJet_ptr->neutralEmEnergyFraction());
	  PFPatJet_NeutralHadEnergy.push_back(PFPatJet_ptr->neutralHadronEnergy());
	  PFPatJet_NeutralHadEnergyFrac.push_back(PFPatJet_ptr->neutralHadronEnergyFraction());
	  PFPatJet_NeutralHadMult.push_back(PFPatJet_ptr->neutralHadronMultiplicity());
	  PFPatJet_NeutralMult.push_back(PFPatJet_ptr->neutralMultiplicity());

	  //Jec Uncertainities
	  jecUnc->setJetEta(PFPatJet_ptr->eta());
	  jecUnc->setJetPt(PFPatJet_ptr->pt());
	  double unc = jecUnc->getUncertainty(true);
	  PFPatJet_jecUncertainity.push_back(unc);
	  
	  if(saveBTaggingInfo_PFPat_){

	    PFPatJet_BJetDiscrByTrackCountingHighEff.push_back(PFPatJet_ptr->bDiscriminator("trackCountingHighEffBJetTags"));
	    PFPatJet_BJetDiscrByTrackCountingHighPur.push_back(PFPatJet_ptr->bDiscriminator("trackCountingHighPurBJetTags"));
	    PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff.push_back(PFPatJet_ptr->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	    PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur.push_back(PFPatJet_ptr->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	    PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff.push_back(PFPatJet_ptr->bDiscriminator("combinedSecondaryVertexBJetTags"));
	    PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA.push_back(PFPatJet_ptr->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	    PFPatJet_BJetDiscrByjetProbabilityBJetTags.push_back(PFPatJet_ptr->bDiscriminator("jetProbabilityBJetTags"));
	    PFPatJet_BJetDiscrByjetBProbabilityBJetTags.push_back(PFPatJet_ptr->bDiscriminator("jetBProbabilityBJetTags"));
	    PFPatJet_BJetDiscrBySoftElectronBJetTags.push_back(PFPatJet_ptr->bDiscriminator("softElectronBJetTags"));
	    PFPatJet_BJetDiscrBySoftMuonBJetTags.push_back(PFPatJet_ptr->bDiscriminator("softMuonBJetTags"));
	    PFPatJet_JetPartonFlavor.push_back(fabs(PFPatJet_ptr->partonFlavour()));

	  } //END_OF if(saveBTaggingInfo_PFPat_)
	} //END_OF if(PFPatJet_ptr->pt() > 30)
      } //END_OF for(PFPatJet_itr = PFPatJetColl->begin(); PFPatJet_itr != PFPatJetColl->end(); PFPatJet_itr++)

      delete jecUnc;
      
      if(savePUJetIdInfo_PFPat_){
	
	if(debug_) cout << "Inside PU Jet Id for PF Pat Jet loop" << endl;

        Handle<edm::View<pat::Jet> > puJets;  //need jets accessed with view as pujet id requires to match to ref_base only possible with view
	iEvent.getByLabel(PFPatJetLabel_, puJets);

        //pileup based jetId:MVA (Containing the value of discriminent)
	Handle<ValueMap<float> > puJetIdMVA_CutBased;
        iEvent.getByLabel(edm::InputTag("puJetMva","cutbasedDiscriminant"),puJetIdMVA_CutBased);
	Handle<ValueMap<float> > puJetIdMVA_Simple;
        iEvent.getByLabel(edm::InputTag("puJetMva","simpleDiscriminant"),puJetIdMVA_Simple);
	Handle<ValueMap<float> > puJetIdMVA_Full;
        iEvent.getByLabel(edm::InputTag("puJetMva","fullDiscriminant"),puJetIdMVA_Full);
        //ID Flags (Containing the Id of MVA)
	Handle<ValueMap<int> > puJetIdFlag_CutBased;
        iEvent.getByLabel(edm::InputTag("puJetMva","cutbasedId"),puJetIdFlag_CutBased);
	Handle<ValueMap<int> > puJetIdFlag_Simple;
        iEvent.getByLabel(edm::InputTag("puJetMva","simpleId"),puJetIdFlag_Simple);
	Handle<ValueMap<int> > puJetIdFlag_Full;
        iEvent.getByLabel(edm::InputTag("puJetMva","fullId"),puJetIdFlag_Full);

	for(unsigned int i = 0; i < puJets->size(); i++){
	  const pat::Jet puJet = puJets->at(i);
	  if(puJet.pt() > 30){
	    float mva_cutbased = (*puJetIdMVA_CutBased)[puJets->refAt(i)];
	    float mva_simple = (*puJetIdMVA_Simple)[puJets->refAt(i)];
	    float mva_full = (*puJetIdMVA_Full)[puJets->refAt(i)];

	    int idflag_cutbased = (*puJetIdFlag_CutBased)[puJets->refAt(i)];
	    int idflag_simple = (*puJetIdFlag_Simple)[puJets->refAt(i)];
	    int idflag_full = (*puJetIdFlag_Full)[puJets->refAt(i)];

	    PFPatJet_puJetIdCutBased_MVA.push_back(mva_cutbased);
	    PFPatJet_puJetIdSimple_MVA.push_back(mva_simple);
	    PFPatJet_puJetIdFull_MVA.push_back(mva_full);

	    bool idflag_cutbased_loose = false;
	    bool idflag_cutbased_medium = false;
	    bool idflag_cutbased_tight = false;

	    bool idflag_simple_loose = false;
	    bool idflag_simple_medium = false;
	    bool idflag_simple_tight = false;

	    bool idflag_full_loose = false;
	    bool idflag_full_medium = false;
	    bool idflag_full_tight = false;

	    if(PileupJetIdentifier::passJetId(idflag_cutbased, PileupJetIdentifier::kLoose)) idflag_cutbased_loose = true;
            if(PileupJetIdentifier::passJetId(idflag_cutbased, PileupJetIdentifier::kMedium)) idflag_cutbased_medium = true;
            if(PileupJetIdentifier::passJetId(idflag_cutbased, PileupJetIdentifier::kTight)) idflag_cutbased_tight = true;

            if(PileupJetIdentifier::passJetId(idflag_simple, PileupJetIdentifier::kLoose)) idflag_simple_loose = true;
            if(PileupJetIdentifier::passJetId(idflag_simple, PileupJetIdentifier::kMedium)) idflag_simple_medium = true;
            if(PileupJetIdentifier::passJetId(idflag_simple, PileupJetIdentifier::kTight)) idflag_simple_tight = true;

            if(PileupJetIdentifier::passJetId(idflag_full, PileupJetIdentifier::kLoose)) idflag_full_loose = true;
            if(PileupJetIdentifier::passJetId(idflag_full, PileupJetIdentifier::kMedium)) idflag_full_medium = true;
            if(PileupJetIdentifier::passJetId(idflag_full, PileupJetIdentifier::kTight)) idflag_full_tight = true;
	   
	    PFPatJet_PassPUJetIdCutBased_loose.push_back(idflag_cutbased_loose);
	    PFPatJet_PassPUJetIdCutBased_medium.push_back(idflag_cutbased_medium);
	    PFPatJet_PassPUJetIdCutBased_tight.push_back(idflag_cutbased_tight);

	    PFPatJet_PassPUJetIdSimple_loose.push_back(idflag_simple_loose);
	    PFPatJet_PassPUJetIdSimple_medium.push_back(idflag_simple_medium);
	    PFPatJet_PassPUJetIdSimple_tight.push_back(idflag_simple_tight);

	    PFPatJet_PassPUJetIdFull_loose.push_back(idflag_full_loose);
	    PFPatJet_PassPUJetIdFull_medium.push_back(idflag_full_medium);
	    PFPatJet_PassPUJetIdFull_tight.push_back(idflag_full_tight);

	  } //END_OF if(puJet.pt() > 30)
	} //END_OF for(unsigned int i = 0; i < puJets->size(); i++)	
      } //END_OF if(savePUJetIdInfo_PFPat_)
      
      if(debug_) cout << "End of PF Pat Jet Loop" << endl;

    } //END_OF if(runPFPatJets_)

    if(runCaloRecoJets_){

      //Accessing reco::CaloJets Collection
      edm::Handle<reco::CaloJetCollection> CaloRecoJetColl;
      iEvent.getByLabel(CaloRecoJetLabel_, CaloRecoJetColl);
      reco::CaloJetCollection::const_iterator CaloRecoJet_itr; 

      std::vector<reco::CaloJet> CaloRecoJet_container;
      CaloRecoJet_container.clear();

      const JetCorrector* corrector = JetCorrector::getJetCorrector(CaloRecoJetCorrectionService_, iSetup);

      for(CaloRecoJet_itr = CaloRecoJetColl->begin(); CaloRecoJet_itr != CaloRecoJetColl->end(); CaloRecoJet_itr++){
        reco::CaloJet CaloRecoJet_corrected = *CaloRecoJet_itr;
	double correction = corrector->correction(*CaloRecoJet_itr, iEvent, iSetup);
	CaloRecoJet_corrected.scaleEnergy(correction);
	CaloRecoJet_container.push_back(CaloRecoJet_corrected);
      }
    
      if(CaloRecoJet_container.size() > 1){
	std::sort(CaloRecoJet_container.begin(), CaloRecoJet_container.end(), PtSortCriteriumCaloRecoJet());
      }

      if(debug_) cout << "CaloRecoJet_container size = " << CaloRecoJet_container.size() << endl;

      //Accessing Jec Uncertainities
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5Calo",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

      for(unsigned int calorecojet = 0; calorecojet < CaloRecoJet_container.size(); calorecojet++){
	const reco::CaloJetRef CaloRecoJet_ptr(&(CaloRecoJet_container), calorecojet);
	if(CaloRecoJet_ptr->pt() > 30){
	  CaloRecoJet_n++;

	  //Kinematic Variables
	  CaloRecoJet_E.push_back(CaloRecoJet_ptr->energy());
	  CaloRecoJet_et.push_back(CaloRecoJet_ptr->et());
	  CaloRecoJet_pt.push_back(CaloRecoJet_ptr->pt());
	  CaloRecoJet_eta.push_back(CaloRecoJet_ptr->eta());
	  CaloRecoJet_phi.push_back(correct_Phi(CaloRecoJet_ptr->phi()));
	  CaloRecoJet_theta.push_back(CaloRecoJet_ptr->theta());
	  CaloRecoJet_px.push_back(CaloRecoJet_ptr->px());
	  CaloRecoJet_py.push_back(CaloRecoJet_ptr->py());
	  CaloRecoJet_pz.push_back(CaloRecoJet_ptr->pz());
	  CaloRecoJet_vx.push_back(CaloRecoJet_ptr->vx());
	  CaloRecoJet_vy.push_back(CaloRecoJet_ptr->vy());
	  CaloRecoJet_vz.push_back(CaloRecoJet_ptr->vz());

	  //Id and other Variables for Calo Reco Jets
	  CaloRecoJet_emEnergyFraction.push_back(CaloRecoJet_ptr->emEnergyFraction());
	  CaloRecoJet_emEnergyInEB.push_back(CaloRecoJet_ptr->emEnergyInEB());
	  CaloRecoJet_emEnergyInEE.push_back(CaloRecoJet_ptr->emEnergyInEE());
	  CaloRecoJet_emEnergyInHF.push_back(CaloRecoJet_ptr->emEnergyInHF());
	  CaloRecoJet_HadronicEnergyFraction.push_back(CaloRecoJet_ptr->energyFractionHadronic());
	  CaloRecoJet_NConstituents.push_back(CaloRecoJet_ptr->getCaloConstituents().size());
	  CaloRecoJet_HadronicEnergyInHB.push_back(CaloRecoJet_ptr->hadEnergyInHB());
	  CaloRecoJet_HadronicEnergyInHE.push_back(CaloRecoJet_ptr->hadEnergyInHE());
	  CaloRecoJet_HadronicEnergyInHF.push_back(CaloRecoJet_ptr->hadEnergyInHF());
	  CaloRecoJet_HadronicEnergyInHO.push_back(CaloRecoJet_ptr->hadEnergyInHO());
	  CaloRecoJet_maxEnergyInEmTowers.push_back(CaloRecoJet_ptr->maxEInEmTowers());
	  CaloRecoJet_maxEnergyInHadTowers.push_back(CaloRecoJet_ptr->maxEInHadTowers());
	  CaloRecoJet_nConstituents60E.push_back(CaloRecoJet_ptr->n60());
	  CaloRecoJet_nConstituents90E.push_back(CaloRecoJet_ptr->n90());
	  CaloRecoJet_towersArea.push_back(CaloRecoJet_ptr->towersArea());
	  
	  //Jec Uncertainities
	  jecUnc->setJetEta(CaloRecoJet_ptr->eta());
	  jecUnc->setJetPt(CaloRecoJet_ptr->pt());
	  double unc = jecUnc->getUncertainty(true);
	  CaloRecoJet_jecUncertainity.push_back(unc);

	} //END_OF if(CaloRecoJet_ptr->pt() > 30)    
      } //END_OF for(unsigned int calorecojet = 0; calorecojet < CaloRecoJet_container.size(); calorecojet++)

      delete jecUnc;

      if(saveBTaggingInfo_CaloReco_){
	Handle<reco::JetTagCollection> CaloReco_bTagCollTCHE;
	iEvent.getByLabel(CaloReco_TrkCountingHighEffBJetTagsLabel_, CaloReco_bTagCollTCHE);
	const reco::JetTagCollection & CaloReco_bTagColl_TCHE = *(CaloReco_bTagCollTCHE.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollTCHP;
	iEvent.getByLabel(CaloReco_TrkCountingHighPurBJetTagsLabel_, CaloReco_bTagCollTCHP);
	const reco::JetTagCollection & CaloReco_bTagColl_TCHP = *(CaloReco_bTagCollTCHP.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollSSVHE;
	iEvent.getByLabel(CaloReco_SimpleSecVtxHighEffBJetTagsLabel_, CaloReco_bTagCollSSVHE);
	const reco::JetTagCollection & CaloReco_bTagColl_SSVHE = *(CaloReco_bTagCollSSVHE.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollSSVHP;
	iEvent.getByLabel(CaloReco_SimpleSecVtxHighPurBJetTagsLabel_, CaloReco_bTagCollSSVHP);
	const reco::JetTagCollection & CaloReco_bTagColl_SSVHP = *(CaloReco_bTagCollSSVHP.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollCSV;
	iEvent.getByLabel(CaloReco_CombinedSecVtxBJetTagsLabel_, CaloReco_bTagCollCSV);
	const reco::JetTagCollection & CaloReco_bTagColl_CSV = *(CaloReco_bTagCollCSV.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollCSVMVA;
	iEvent.getByLabel(CaloReco_CombinedSecVtxMVABJetTagsLabel_, CaloReco_bTagCollCSVMVA);
	const reco::JetTagCollection & CaloReco_bTagColl_CSVMVA = *(CaloReco_bTagCollCSVMVA.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollJP;
	iEvent.getByLabel(CaloReco_JetProbBJetTagsLabel_, CaloReco_bTagCollJP);
	const reco::JetTagCollection & CaloReco_bTagColl_JP = *(CaloReco_bTagCollJP.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollJBP;
	iEvent.getByLabel(CaloReco_JetBProbBJetTagsLabel_, CaloReco_bTagCollJBP);
	const reco::JetTagCollection & CaloReco_bTagColl_JBP = *(CaloReco_bTagCollJBP.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollSE;
	iEvent.getByLabel(CaloReco_SoftElectronBJetTagsLabel_, CaloReco_bTagCollSE);
	const reco::JetTagCollection & CaloReco_bTagColl_SE = *(CaloReco_bTagCollSE.product());

	Handle<reco::JetTagCollection> CaloReco_bTagCollSM;
	iEvent.getByLabel(CaloReco_SoftMuonBJetTagsLabel_, CaloReco_bTagCollSM);
	const reco::JetTagCollection & CaloReco_bTagColl_SM = *(CaloReco_bTagCollSM.product());

	for(unsigned int calorecojet = 0; calorecojet < CaloRecoJet_container.size(); calorecojet++){
	  const reco::CaloJet & caloJet = CaloRecoJet_container.at(calorecojet);

	  float TCHE = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_TCHE, caloJet); 
	  float TCHP = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_TCHP, caloJet);
	  float SSVHE = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_SSVHE, caloJet);
	  float SSVHP = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_SSVHP, caloJet);
	  float CSV = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_CSV, caloJet);
	  float CSVMVA = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_CSVMVA, caloJet);
	  float JP = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_JP, caloJet);
	  float JBP = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_JBP, caloJet);
	  float SE = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_SE, caloJet);
	  float SM = GetMatchedBTagDiscrValueToCaloJet(CaloReco_bTagColl_SM, caloJet);

	  CaloRecoJet_BJetDiscrByTrackCountingHighEff.push_back(TCHE);
	  CaloRecoJet_BJetDiscrByTrackCountingHighPur.push_back(TCHP);
	  CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff.push_back(SSVHE);
	  CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur.push_back(SSVHP);
	  CaloRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff.push_back(CSV);
	  CaloRecoJet_BJetDiscrByCombinedSecondaryVertexMVA.push_back(CSVMVA);
	  CaloRecoJet_BJetDiscrByjetProbabilityBJetTags.push_back(JP);
	  CaloRecoJet_BJetDiscrByjetBProbabilityBJetTags.push_back(JBP);
	  CaloRecoJet_BJetDiscrBySoftElectronBJetTags.push_back(SE);
	  CaloRecoJet_BJetDiscrBySoftMuonBJetTags.push_back(SM);

	}
      } //END_OF if(saveBTaggingInfo_CaloReco_)

      if(debug_) cout << "End of Calo Reco Jet Loop" << endl;

    } //END_OF if(runCaloRecoJets_)

    if(runPFRecoJets_){

      //Accessing reco::PFJets Collection
      edm::Handle<reco::PFJetCollection> PFRecoJetColl;
      iEvent.getByLabel(PFRecoJetLabel_, PFRecoJetColl);
      reco::PFJetCollection::const_iterator PFRecoJet_itr;

      std::vector<reco::PFJet> PFRecoJet_container;
      PFRecoJet_container.clear();

      const JetCorrector* corrector = JetCorrector::getJetCorrector(PFRecoJetCorrectionService_, iSetup);

      for(PFRecoJet_itr = PFRecoJetColl->begin(); PFRecoJet_itr != PFRecoJetColl->end(); PFRecoJet_itr++){
	reco::PFJet PFRecoJet_corrected = *PFRecoJet_itr;
	double correction = corrector->correction(*PFRecoJet_itr, iEvent, iSetup);
	PFRecoJet_corrected.scaleEnergy(correction);
	PFRecoJet_container.push_back(PFRecoJet_corrected);
      }

      if(PFRecoJet_container.size() > 1){
	std::sort(PFRecoJet_container.begin(), PFRecoJet_container.end(), PtSortCriteriumPFRecoJet());
      }

      if(debug_) cout << "PFRecoJet_container size = " << PFRecoJet_container.size() << endl;

      //Accessing Jec Uncertainities
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

      for(unsigned int pfrecojet = 0; pfrecojet < PFRecoJet_container.size(); pfrecojet++){
	const reco::PFJetRef PFRecoJet_ptr(&(PFRecoJet_container), pfrecojet);
	if(PFRecoJet_ptr->pt() > 30){
	  PFRecoJet_n++;

	  //Kinematic Variables
	  PFRecoJet_E.push_back(PFRecoJet_ptr->energy());
	  PFRecoJet_et.push_back(PFRecoJet_ptr->et());
	  PFRecoJet_pt.push_back(PFRecoJet_ptr->pt());
	  PFRecoJet_eta.push_back(PFRecoJet_ptr->eta());
	  PFRecoJet_phi.push_back(correct_Phi(PFRecoJet_ptr->phi()));
	  PFRecoJet_theta.push_back(PFRecoJet_ptr->theta());
	  PFRecoJet_px.push_back(PFRecoJet_ptr->px());
	  PFRecoJet_py.push_back(PFRecoJet_ptr->py());
	  PFRecoJet_pz.push_back(PFRecoJet_ptr->pz());
	  PFRecoJet_vx.push_back(PFRecoJet_ptr->vx());
	  PFRecoJet_vy.push_back(PFRecoJet_ptr->vy());
	  PFRecoJet_vz.push_back(PFRecoJet_ptr->vz());

	  //Id and other Variables for PF Reco Jets
	  PFRecoJet_ChargedEmEnergy.push_back(PFRecoJet_ptr->chargedEmEnergy());
	  PFRecoJet_ChargedEmEnergyFrac.push_back(PFRecoJet_ptr->chargedEmEnergyFraction());
	  PFRecoJet_ChargedHadEnergy.push_back(PFRecoJet_ptr->chargedHadronEnergy());
	  PFRecoJet_ChargedHadEnergyFrac.push_back(PFRecoJet_ptr->chargedHadronEnergyFraction());
	  PFRecoJet_ChargedHadMult.push_back(PFRecoJet_ptr->chargedHadronMultiplicity());
	  PFRecoJet_ChargedMult.push_back(PFRecoJet_ptr->chargedMultiplicity());
	  PFRecoJet_NConstituents.push_back(PFRecoJet_ptr->getPFConstituents().size());
	  PFRecoJet_HFEMEnergy.push_back(PFRecoJet_ptr->HFEMEnergy());
	  PFRecoJet_HFEMEnergyFrac.push_back(PFRecoJet_ptr->HFEMEnergyFraction());
	  PFRecoJet_HFEMMult.push_back(PFRecoJet_ptr->HFEMMultiplicity());
	  PFRecoJet_HFHadEnergy.push_back(PFRecoJet_ptr->HFHadronEnergy());
	  PFRecoJet_HFHadEnergyFrac.push_back(PFRecoJet_ptr->HFHadronEnergyFraction());
	  PFRecoJet_HFHadMult.push_back(PFRecoJet_ptr->HFHadronMultiplicity());
	  PFRecoJet_NeutralEmEnergy.push_back(PFRecoJet_ptr->neutralEmEnergy());
	  PFRecoJet_NeutralEmEnergyFrac.push_back(PFRecoJet_ptr->neutralEmEnergyFraction());
	  PFRecoJet_NeutralHadEnergy.push_back(PFRecoJet_ptr->neutralHadronEnergy());
	  PFRecoJet_NeutralHadEnergyFrac.push_back(PFRecoJet_ptr->neutralHadronEnergyFraction());
	  PFRecoJet_NeutralHadMult.push_back(PFRecoJet_ptr->neutralHadronMultiplicity());
	  PFRecoJet_NeutralMult.push_back(PFRecoJet_ptr->neutralMultiplicity());

	  //Jec Uncertainities
	  jecUnc->setJetEta(PFRecoJet_ptr->eta());
	  jecUnc->setJetPt(PFRecoJet_ptr->pt());
	  double unc = jecUnc->getUncertainty(true);
	  PFRecoJet_jecUncertainity.push_back(unc);

	} //END_OF if(PFRecoJet_ptr->pt() > 30)
      } //END_OF for(unsigned int pfrecojet = 0; pfrecojet < PFRecoJet_container.size(); pfrecojet++)

      delete jecUnc;

      if(saveBTaggingInfo_PFReco_){
	Handle<reco::JetTagCollection> PFReco_bTagCollTCHE;
	iEvent.getByLabel(PFReco_TrkCountingHighEffBJetTagsLabel_, PFReco_bTagCollTCHE);
	const reco::JetTagCollection & PFReco_bTagColl_TCHE = *(PFReco_bTagCollTCHE.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollTCHP;
	iEvent.getByLabel(PFReco_TrkCountingHighPurBJetTagsLabel_, PFReco_bTagCollTCHP);
	const reco::JetTagCollection & PFReco_bTagColl_TCHP = *(PFReco_bTagCollTCHP.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollSSVHE;
	iEvent.getByLabel(PFReco_SimpleSecVtxHighEffBJetTagsLabel_, PFReco_bTagCollSSVHE);
	const reco::JetTagCollection & PFReco_bTagColl_SSVHE = *(PFReco_bTagCollSSVHE.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollSSVHP;
	iEvent.getByLabel(PFReco_SimpleSecVtxHighPurBJetTagsLabel_, PFReco_bTagCollSSVHP);
	const reco::JetTagCollection & PFReco_bTagColl_SSVHP = *(PFReco_bTagCollSSVHP.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollCSV;
	iEvent.getByLabel(PFReco_CombinedSecVtxBJetTagsLabel_, PFReco_bTagCollCSV);
	const reco::JetTagCollection & PFReco_bTagColl_CSV = *(PFReco_bTagCollCSV.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollCSVMVA;
	iEvent.getByLabel(PFReco_CombinedSecVtxMVABJetTagsLabel_, PFReco_bTagCollCSVMVA);
	const reco::JetTagCollection & PFReco_bTagColl_CSVMVA = *(PFReco_bTagCollCSVMVA.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollJP;
	iEvent.getByLabel(PFReco_JetProbBJetTagsLabel_, PFReco_bTagCollJP);
	const reco::JetTagCollection & PFReco_bTagColl_JP = *(PFReco_bTagCollJP.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollJBP;
	iEvent.getByLabel(PFReco_JetBProbBJetTagsLabel_, PFReco_bTagCollJBP);
	const reco::JetTagCollection & PFReco_bTagColl_JBP = *(PFReco_bTagCollJBP.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollSE;
	iEvent.getByLabel(PFReco_SoftElectronBJetTagsLabel_, PFReco_bTagCollSE);
	const reco::JetTagCollection & PFReco_bTagColl_SE = *(PFReco_bTagCollSE.product());

	Handle<reco::JetTagCollection> PFReco_bTagCollSM;
	iEvent.getByLabel(PFReco_SoftMuonBJetTagsLabel_, PFReco_bTagCollSM);
	const reco::JetTagCollection & PFReco_bTagColl_SM = *(PFReco_bTagCollSM.product());
    
	for(unsigned int pfrecojet = 0; pfrecojet < PFRecoJet_container.size(); pfrecojet++){
	  const reco::PFJet & pfJet = PFRecoJet_container.at(pfrecojet);

	  float TCHE = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_TCHE, pfJet);
	  float TCHP = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_TCHP, pfJet);
	  float SSVHE = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_SSVHE, pfJet);
	  float SSVHP = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_SSVHP, pfJet);
	  float CSV = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_CSV, pfJet);
	  float CSVMVA = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_CSVMVA, pfJet);
	  float JP = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_JP, pfJet);
	  float JBP = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_JBP, pfJet);
	  float SE = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_SE, pfJet);
	  float SM = GetMatchedBTagDiscrValueToPFJet(PFReco_bTagColl_SM, pfJet);

	  PFRecoJet_BJetDiscrByTrackCountingHighEff.push_back(TCHE);
	  PFRecoJet_BJetDiscrByTrackCountingHighPur.push_back(TCHP);
	  PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff.push_back(SSVHE);
	  PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur.push_back(SSVHP);
	  PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff.push_back(CSV);
	  PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA.push_back(CSVMVA);
	  PFRecoJet_BJetDiscrByjetProbabilityBJetTags.push_back(JP);
	  PFRecoJet_BJetDiscrByjetBProbabilityBJetTags.push_back(JBP);
	  PFRecoJet_BJetDiscrBySoftElectronBJetTags.push_back(SE);
	  PFRecoJet_BJetDiscrBySoftMuonBJetTags.push_back(SM);

	}
      } //END_OF if(saveBTaggingInfo_PFReco_)

      if(debug_) cout << "End of PF Reco Jet Loop" << endl;

    } //END_OF if(runPFRecoJets_)

  } //END_OF if(runJets_)


  //////////////////////////////
  //                          //
  //  Get Vertex Information  //
  //                          //
  //////////////////////////////

  if(runVertex_){
    Handle<reco::VertexCollection> recVtxsColl;
    iEvent.getByLabel(vertexLabel_, recVtxsColl);

    reco::VertexCollection::const_iterator Vtx_itr;

    for(Vtx_itr = recVtxsColl->begin(); Vtx_itr != recVtxsColl->end(); Vtx_itr++){
      reco::Vertex myVtx = (*Vtx_itr);
      Vertex_n++;

      Vertex_x.push_back(myVtx.x());
      Vertex_y.push_back(myVtx.y());
      Vertex_z.push_back(myVtx.z());
      Vertex_chi2.push_back(myVtx.chi2());
      Vertex_nchi2.push_back(myVtx.normalizedChi2());
      Vertex_ndof.push_back(myVtx.ndof());
      Vertex_tracksSize.push_back(myVtx.tracksSize());
      Vertex_isFake.push_back(myVtx.isFake());
      Vertex_isValid.push_back(myVtx.isValid());
      Vertex_d0.push_back(myVtx.position().rho());
    }
  }//END_OF if(runVertex_)


  //////////////////////////////
  //                          //
  //  Get Tracks Information  //
  //                          //
  //////////////////////////////

  if(runTracks_){
    Handle<reco::TrackCollection> tracksColl;
    iEvent.getByLabel(TracksLabel_, tracksColl);
    reco::TrackCollection::const_iterator Track_itr;

    std::vector<reco::Track> Track_container;
    Track_container.clear();

    for(Track_itr = tracksColl->begin(); Track_itr != tracksColl->end(); Track_itr++){
      const reco::Track* Track_ptr = &(*Track_itr);
      if(Track_ptr->pt() > 0.5){
	Track_container.push_back(*Track_ptr);
      }
    }

    if(Track_container.size() > 1){
      std::sort(Track_container.begin(), Track_container.end(), PtSortCriteriumTracks());
    }

    for(unsigned int trk = 0; trk < Track_container.size(); trk++){
      const reco::TrackRef Trk_ptr(&(Track_container), trk);
      Tracks_n++;
      Track_pt.push_back(Trk_ptr->pt());
      Track_px.push_back(Trk_ptr->px());
      Track_py.push_back(Trk_ptr->py());
      Track_pz.push_back(Trk_ptr->pz());
      Track_vx.push_back(Trk_ptr->vx());
      Track_vy.push_back(Trk_ptr->vy());
      Track_vz.push_back(Trk_ptr->vz());
      Track_eta.push_back(Trk_ptr->eta());
      Track_phi.push_back(correct_Phi(Trk_ptr->phi()));
      Track_theta.push_back(Trk_ptr->theta());
      Track_chi2.push_back(Trk_ptr->chi2());
    }     
  } //END_OF if(runTracks_)


  ////////////////////////////////
  //                            //
  //  Get Scraping Information  //
  //                            //
  ////////////////////////////////

  if(runScraping_){

    //This code for scraping event is taken from DPGAnalysis/Skims/src/FilterOutScraping.cc (CMSSW_5_3_8)

    IsScrapingEvent_ = true;
    Scraping_FractionOfGoodTracks = 0.0;

    // get GeneralTracks collection
    Handle<reco::TrackCollection> trkRef;
    iEvent.getByLabel(TracksLabel_, trkRef);
    const reco::TrackCollection* trkColl = trkRef.product();

    int num_highpurity=0;
    reco::TrackBase::TrackQuality trackQuality_ = reco::TrackBase::qualityByName("highPurity");

    if(trkColl->size() > 10){

      for(reco::TrackCollection::const_iterator trk_itr = trkColl->begin(); trk_itr != trkColl->end(); trk_itr++){
	if(trk_itr->quality(trackQuality_)) num_highpurity++;
      }

      Scraping_FractionOfGoodTracks = (float)num_highpurity/(float)trkColl->size();
      if(Scraping_FractionOfGoodTracks > 0.2) IsScrapingEvent_ = false;

    }else{
      //if less than 10 Tracks, mark the event as noScrapingEvent
      IsScrapingEvent_ = false;
    }
  } //END_OF if(runScraping_)


  ///////////////////////////////
  //                           //
  //  Get Trigger Information  //
  //                           //
  ///////////////////////////////

  //L1





  //HLT

  if(runHLT_){

    std::string hltName = "HLT_Photon";

    //-----------------------------------------------------
    //Accessing TriggerResults and triggerNames collections
    //-----------------------------------------------------
    edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
    iEvent.getByLabel(triggerEventLabel_,triggerEventHandle_);
    if(!triggerEventHandle_.isValid()){
      cout << "HLT TriggerEvent Handle is not valid______Can't get the product" << endl;
    }
    else{
      edm::Handle<edm::TriggerResults> triggerResults;
      iEvent.getByLabel(triggerResultsLabel_, triggerResults);

      const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
      const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();

      for(unsigned int trig = 0; trig < triggerNames_.size(); trig++){
	if(triggerNames_[trig].find(hltName.c_str())!=std::string::npos){
	  HLT_Photon_triggers.push_back(triggerNames_[trig]);
	}
      }



      for(unsigned int l = 0; l < HLT_Photon_triggers.size(); l++){
	cout << "trigger = " << HLT_Photon_triggers[l] << endl;}



      HLT_Photon_nTriggers = HLT_Photon_triggers.size();


      for(int i = 0; i < HLT_Photon_nTriggers; i++){
	int idx = (int)triggerNames.triggerIndex(HLT_Photon_triggers[i]);
	bool accept = triggerResults->accept(idx);
	int pre = (int)hltConfig_.prescaleValue(iEvent, iSetup, HLT_Photon_triggers[i]);

	HLT_Photon_triggerIndex.push_back(idx);
	HLT_Photon_ifTriggerPassed.push_back(accept);
	HLT_Photon_trig_prescales.push_back(pre);
      }


      for(int o = 0; o < HLT_Photon_nTriggers; o++){
	cout << "trig index = " << HLT_Photon_triggerIndex[o] << endl;
	cout << "accept ? = " << HLT_Photon_ifTriggerPassed[o] << endl;
	cout << "pre = " << HLT_Photon_trig_prescales[o] << endl;}


      cout << "No. of triggers = " << HLT_Photon_nTriggers << "," << HLT_Photon_triggerIndex.size() << "," << HLT_Photon_trig_prescales.size() << endl;



      for(unsigned int i = 0; i < HLT_Photon_triggers.size(); i++){
	vector<std::string> saveTagFilters;
	saveTagFilters.clear();
	if(hltConfig_.triggerIndex(HLT_Photon_triggers[i])<hltConfig_.size()){
	  saveTagFilters = hltConfig_.saveTagsModules(HLT_Photon_triggers[i]);
	}
	HLT_Photon_nFilters.push_back(saveTagFilters.size());
	HLT_Photon_FilterNames.push_back(saveTagFilters);
      }


      for(unsigned int k = 0; k < HLT_Photon_FilterNames.size(); k++){
	for(unsigned int n = 0; n < HLT_Photon_FilterNames[k].size(); n++){
	  cout << " Filter name for [" << k << "][" << n << "] = " << HLT_Photon_FilterNames[k][n] << endl;
	}
      }


      cout << " HLT_Photon_FilterNames.size() = " << HLT_Photon_FilterNames.size() << " and " << HLT_Photon_FilterNames[2].size() << endl;

      for(unsigned int k = 0; k < HLT_Photon_FilterNames.size(); k++){
	for(unsigned int n = 0; n < HLT_Photon_FilterNames[k].size(); n++){
	  vector<double> Obj_pt;
	  vector<double> Obj_eta;
	  vector<double> Obj_phi;
	  Obj_pt.clear();
	  Obj_eta.clear();
	  Obj_phi.clear();
	  trigger::size_type filterIndex = triggerEventHandle_->filterIndex(edm::InputTag(HLT_Photon_FilterNames[k][n], "", triggerEventLabel_.process()));
	  if(filterIndex < triggerEventHandle_->sizeFilters()){
	    const trigger::Keys& trigKeys = triggerEventHandle_->filterKeys(filterIndex);
	    const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle_->getObjects();
	    for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
	      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	      Obj_pt.push_back(obj.pt());
	      Obj_eta.push_back(obj.eta());
	      Obj_phi.push_back(correct_Phi(obj.phi()));
	    }
	  }
	  HLT_Photon_FilterObjects_pt.push_back(Obj_pt);
	  HLT_Photon_FilterObjects_eta.push_back(Obj_eta);
	  HLT_Photon_FilterObjects_phi.push_back(Obj_phi);
	}
      }


      for(unsigned int k = 0; k < HLT_Photon_FilterObjects_pt.size(); k++){
	for(unsigned int n = 0; n < HLT_Photon_FilterObjects_pt[k].size(); n++){
 	  cout << " Obj pt for [" << k << "][" << n << "] = " << HLT_Photon_FilterObjects_pt[k][n] << endl;}}

      cout << " HLT_Photon_FilterObjects_pt.size()" << HLT_Photon_FilterObjects_pt.size()<<endl;


    }
  } //END_OF if(runHLT_)

 
   myEvent->Fill();
   if(debug_) cout << "ANALYZE LOOP DONE" << endl;
}


// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analyzer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  if(isRealData){

    edm::Handle<LumiSummary> lumiSummary;                       ///////////////////////////////////////////////////////////////////////////////////
    iLumi.getByLabel("lumiProducer", lumiSummary);              // InsLumi is lumi in 1 sec and 1 LumiBlock is 93.244 seconds long. So to get    //
                                                                // lumi in a LumiBlock we multiply insLumi by 93.244. This method is implemented //  
    deliveredLumi += lumiSummary->avgInsDelLumi()*93.244;       // in cmssw in RecoLuminosity/LumiProducer/plugins/LumiCalculator.cc.            //
    recordedLumi += deliveredLumi*lumiSummary->liveFrac();      // By using this method in endLuminosityBlock will enable to keep on adding Lumis//
    cout << "deliveredLumi = " << deliveredLumi << ","          // from each lumiblock and will give the total lumi for all the events           //  
         << "recordedLumi = " << recordedLumi << endl;          ///////////////////////////////////////////////////////////////////////////////////
  }                                                        
}                                                          


// ------------ method called when ending the processing of a run  ------------
void 
Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob() 
{
  f->WriteTObject(myEvent);
  delete myEvent;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//
//Destructor
//
Analyzer::~Analyzer()
{
  f->Close();
  delete f;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}




//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);

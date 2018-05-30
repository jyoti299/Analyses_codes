// system include files
//#include <memory>
//#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//Event Info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

//Lumi Info
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

//PileUp Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//PatPhoton
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

//StringToEnumValue required for severity flags
#include "CommonTools/Utils/interface/StringToEnumValue.h"

//EcalClusterTools and Ecaltools needed for cov, localCov, swissCross, roundness, angle etc.
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"

//PFPhotonIsolation
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

//Photon Electron Veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//EcalRecHit Collection
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//EcalSeverityLevelAlgos
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

//CaloTopology
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

//CaloGeometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

//Pat Jet Collection
#include "DataFormats/PatCandidates/interface/Jet.h"

//Reco Jet Collections
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

//Reco Jet correction on the fly
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/JetReco/interface/JetID.h"

//Jet Energy Correction Uncertainity
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//PileUp Jet Id
#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include "CMGTools/External/interface/PileupJetIdAlgo.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

//Trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"



//ROOT include files
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TRint.h"
#include <memory>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class std::vector< std::vector< int > >+;
#pragma link C++ class std::vector< std::vector< std::string > >+;
#pragma link C++ class std::vector< std::vector< TString > >+;
#pragma link C++ class std::vector< std::vector< float > >+;
#pragma link C++ class std::vector< std::vector< bool > >+;
#pragma extra_include "std::vector";
#endif

//Utility Function Prototypes (defined in Utilities.cc)
double correct_Phi(double phi);
float correct_Phi(float phi);
double Theta(double eta);
float Theta(float eta);
double Pl(double P, double Pt);
float Pl(float P, float Pt);
float recHitE( const  DetId id,  const EcalRecHitCollection &recHits );
float recHitE( const DetId id, const EcalRecHitCollection & recHits, int di, int dj );
float recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits );
float recHitE( const DetId id, const EcalRecHitCollection &recHits, bool useTimingInfo );
const std::vector<DetId> neighbours(const DetId& id);



//
// class declaration
//

class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      void ClearVectors();
      float GetE4OverE1(const DetId& id, const EcalRecHitCollection& rhs);
      float GetE6OverE2(const DetId& id, const EcalRecHitCollection& rhs);
      float GetE2OverE9(const DetId& id, const EcalRecHitCollection& rhs);
      float GetE2OverE9(const DetId& id, const EcalRecHitCollection& rhs, float recHitEtThreshold, 
                        float recHitEtThreshold2, bool avoidIeta85, bool KillSecondHit);
      float GetrookFraction(const DetId& id, const EcalRecHitCollection& rhs);
      float GetMatchedBTagDiscrValueToCaloJet(const reco::JetTagCollection & bTags, const reco::CaloJet & jet);
      float GetMatchedBTagDiscrValueToPFJet(const reco::JetTagCollection & bTags, const reco::PFJet & jet);


      //--------------------------------------------------
      // ----------member data ---------------------------
      //--------------------------------------------------

      TFile                                     *f;
      TTree                                     *myEvent;

      std::string                               OutFile_;
      int                                       debug_;
      unsigned int                              RunNumber, EventNumber, LumiNumber, BXNumber;
      int                                       nEvents;
      bool                                      isRealData, isItAOD_;


      //----------- DECLEARING PILEUP VARIABLES-----------
      bool                                      runMCPileUp_;
      int                                       npuVertices, npuVerticesm1, npuVerticesp1, npuVerticespm2;
      float                                     trueNumofInteractions;

      //----------- DECLEARING LUMINOSITY VARIABLES-------
      bool                                      saveLuminosityInfo_;
      unsigned int                              totalIntensityBeam1, totalIntensityBeam2;
      float                                     avgInsDelLumi, avgInsDelLumiErr, avgInsRecLumi, avgInsRecLumiErr; 
      float                                     deliveredLumi, recordedLumi;


      //----------- DECLEARING PHOTONS VARIABLES-----------
      bool                                      runPhotons_;
      bool                                      savePhotonCrystalInfo_;
      bool                                      savePhotonPFIsolation_;
      bool                                      savePhotonRecHitsInfo_;
      int                                       Photon_n;


      //Kinematic Variables
      vector<double>                            Photon_E;
      vector<double>                            Photon_et;
      vector<double>                            Photon_pt;
      vector<double>                            Photon_eta;
      vector<double>                            Photon_phi;
      vector<double>                            Photon_theta;
      vector<double>                            Photon_px;
      vector<double>                            Photon_py;
      vector<double>                            Photon_pz;
      vector<double>                            Photon_vx;
      vector<double>                            Photon_vy;
      vector<double>                            Photon_vz;

      //Shower Shape Variables
      vector<float>                            Photon_r9;
      vector<float>                            Photon_maxEnergyXtal;
      vector<float>                            Photon_e1x5;
      vector<float>                            Photon_e2x5;
      vector<float>                            Photon_e3x3;
      vector<float>                            Photon_e5x5;
      vector<float>                            Photon_r1x5;
      vector<float>                            Photon_r2x5;
      vector<float>                            Photon_SigmaEtaEta;
      vector<float>                            Photon_SigmaIEtaIEta;
      vector<float>                            Photon_SigmaEtaPhi;
      vector<float>                            Photon_SigmaIEtaIPhi;
      vector<float>                            Photon_SigmaPhiPhi;
      vector<float>                            Photon_SigmaIPhiIPhi;
      vector<float>                            Photon_roundness;   
      vector<float>                            Photon_angle;
      vector<float>                            Photon_swissCross;
      vector<float>                            Photon_s9;
      vector<float>                            Photon_e4Overe1;
      vector<float>                            Photon_e6Overe2;
      vector<float>                            Photon_e2Overe9; 
      vector<float>                            Photon_rookFraction;

      //Flags for sigmaIphiIphi etc.
      vector<int>                              flagsexclEB_;
      vector<int>                              flagsexclEE_;
      vector<int>                              severitiesexclEB_;
      vector<int>                              severitiesexclEE_;
      const vector<std::string>                flagnamesEB;
      const vector<std::string>                flagnamesEE;
      const vector<std::string>                severitynamesEB;
      const vector<std::string>                severitynamesEE;

      //Fiducial Flags Variables
      vector<bool>                             Photon_isEB;
      vector<bool>                             Photon_isEE;
      vector<bool>                             Photon_isEBGap;
      vector<bool>                             Photon_isEEGap;
      vector<bool>                             Photon_isEBEEGap;
             
      //Detector Isolation Variables
      vector<float>                            Photon_ecalRecHitSumEtConeDR03;
      vector<float>                            Photon_hcalTowerSumEtConeDR03;
      vector<float>                            Photon_hcalDepth1TowerSumEtConeDR03;
      vector<float>                            Photon_hcalDepth2TowerSumEtConeDR03;
      vector<float>                            Photon_trkSumPtSolidConeDR03;
      vector<float>                            Photon_trkSumPtHollowConeDR03;
      vector<int>                              Photon_nTrkSolidConeDR03;
      vector<int>                              Photon_nTrkHollowConeDR03;
      vector<float>                            Photon_ecalRecHitSumEtConeDR04;
      vector<float>                            Photon_hcalTowerSumEtConeDR04;
      vector<float>                            Photon_hcalDepth1TowerSumEtConeDR04;
      vector<float>                            Photon_hcalDepth2TowerSumEtConeDR04;
      vector<float>                            Photon_trkSumPtSolidConeDR04;
      vector<float>                            Photon_trkSumPtHollowConeDR04;
      vector<int>                              Photon_nTrkSolidConeDR04;
      vector<int>                              Photon_nTrkHollowConeDR04;

      //Photon Identification Variables
      vector<float>                            Photon_HoE;
      vector<float>                            Photon_SingleTowerHoE;
      vector<bool>                             Photon_hasConvTrk;
      vector<bool>                             Photon_hasPixelSeed;
      vector<bool>                             passedConvSafeElectronVeto;  

      //Photon Super Cluster Variables
      vector<int>                              Photon_SC_nOfBasicClusters;
      vector<double>                           Photon_SC_rawEnergy;
      vector<double>                           Photon_SC_preShowerEnergy;
      vector<double>                           Photon_SC_energy;
      vector<double>                           Photon_SC_eta;
      vector<double>                           Photon_SC_phi;
      vector<double>                           Photon_SC_x;
      vector<double>                           Photon_SC_y;
      vector<double>                           Photon_SC_z;
      vector<double>                           Photon_SC_etaWidth;
      vector<double>                           Photon_SC_phiWidth;
    
      //Photon MIP Variables
      vector<float>                            Photon_mipChi2;
      vector<float>                            Photon_mipTotEnergy;
      vector<float>                            Photon_mipSlope;
      vector<float>                            Photon_mipIntercept;
      vector<int>                              Photon_mipNhitCone;
      vector<bool>                             Photon_mipIsHalo;

      //Converted Photon Variables
      vector<unsigned int>                     Photon_nConvTracks;
      vector<bool>                             Photon_isConverted;
      vector<float>                            Photon_pairInvariantMass;
      vector<float>                            Photon_pairCotThetaSeparation;
      vector<float>                            Photon_pairMomentum_x;
      vector<float>                            Photon_pairMomentum_y;
      vector<float>                            Photon_pairMomentum_z;
      vector<float>                            Photon_EoverP;
      vector<float>                            Photon_conv_vx;
      vector<float>                            Photon_conv_vy;
      vector<float>                            Photon_conv_vz;
      vector<float>                            Photon_zOfPrimaryVtxFromTrks;
      vector<float>                            Photon_distOfMinimumApproach;
      vector<float>                            Photon_dPhiTracksAtVtx;
      vector<float>                            Photon_dPhiTracksAtEcal;
      vector<float>                            Photon_dEtaTracksAtEcal;

      //Photon Crystal variables
      vector<int>                              Photon_nCrystals; //No. of xtals in which hits are found for a photon     
      vector<vector<float> >                   Photon_xtal_timing;
      vector<vector<float> >                   Photon_xtal_timeErr;
      vector<float>                            Photon_avgTimeAllxtals;
      vector<vector<float> >                   Photon_xtal_energy;
      vector<vector<int> >                     Photon_xtal_EBieta;
      vector<vector<int> >                     Photon_xtal_EBiphi;
      vector<vector<int> >                     Photon_xtal_EBrecoFlag;
      vector<float>                            temp_Photon_xtal_timing;
      vector<float>                            temp_Photon_xtal_timeErr;
      vector<float>                            temp_Photon_xtal_energy;
      vector<int>                              temp_Photon_xtal_EBieta;
      vector<int>                              temp_Photon_xtal_EBiphi;
      vector<int>                              temp_Photon_xtal_EBrecoFlag;

      //PF Isolation Variables
      PFIsolationEstimator                      phoIsolator03;
      vector<double>                            PFIsoPhoton03;
      vector<double>                            PFIsoNeutral03;
      vector<double>                            PFIsoCharged03;
      vector<double>                            PFIsoSum03;
      vector<double>                            PFIsoChargedWorstvtx03;

      //----------- DECLEARING JETS VARIABLES-----------
    
      bool                                      runJets_;
      bool                                      runCaloPatJets_;
      bool                                      runPFPatJets_;
      bool                                      runCaloRecoJets_;
      bool                                      runPFRecoJets_;
      bool                                      savePUJetIdInfo_PFPat_;
      bool                                      saveBTaggingInfo_CaloPat_;
      bool                                      saveBTaggingInfo_PFPat_;
      bool                                      saveBTaggingInfo_CaloReco_;
      bool                                      saveBTaggingInfo_PFReco_;

      int                                       CaloPatJet_n;
      int                                       PFPatJet_n;
      int                                       CaloRecoJet_n;
      int                                       PFRecoJet_n;
      std::string                               CaloRecoJetCorrectionService_;
      std::string                               PFRecoJetCorrectionService_;

      //Kinematic Variables for Calo Pat Jets
      vector<double>                            CaloPatJet_E;
      vector<double>                            CaloPatJet_et;
      vector<double>                            CaloPatJet_pt;
      vector<double>                            CaloPatJet_eta;
      vector<double>                            CaloPatJet_phi;
      vector<double>                            CaloPatJet_theta;
      vector<double>                            CaloPatJet_px;
      vector<double>                            CaloPatJet_py;
      vector<double>                            CaloPatJet_pz;
      vector<double>                            CaloPatJet_vx;
      vector<double>                            CaloPatJet_vy;
      vector<double>                            CaloPatJet_vz;

      //Kinematic Variables for PF Pat Jets
      vector<double>                            PFPatJet_E;
      vector<double>                            PFPatJet_et;
      vector<double>                            PFPatJet_pt;
      vector<double>                            PFPatJet_eta;
      vector<double>                            PFPatJet_phi;
      vector<double>                            PFPatJet_theta;
      vector<double>                            PFPatJet_px;
      vector<double>                            PFPatJet_py;
      vector<double>                            PFPatJet_pz;
      vector<double>                            PFPatJet_vx;
      vector<double>                            PFPatJet_vy;
      vector<double>                            PFPatJet_vz;

      //Kinematic Variables for Reco Calo Jets
      vector<double>                            CaloRecoJet_E;
      vector<double>                            CaloRecoJet_et;
      vector<double>                            CaloRecoJet_pt;
      vector<double>                            CaloRecoJet_eta;
      vector<double>                            CaloRecoJet_phi;
      vector<double>                            CaloRecoJet_theta;
      vector<double>                            CaloRecoJet_px;
      vector<double>                            CaloRecoJet_py;
      vector<double>                            CaloRecoJet_pz;
      vector<double>                            CaloRecoJet_vx;
      vector<double>                            CaloRecoJet_vy;
      vector<double>                            CaloRecoJet_vz;

      //Kinematic Variables for Reco PF Jets
      vector<double>                            PFRecoJet_E;
      vector<double>                            PFRecoJet_et;
      vector<double>                            PFRecoJet_pt;
      vector<double>                            PFRecoJet_eta;
      vector<double>                            PFRecoJet_phi;
      vector<double>                            PFRecoJet_theta;
      vector<double>                            PFRecoJet_px;
      vector<double>                            PFRecoJet_py;
      vector<double>                            PFRecoJet_pz;
      vector<double>                            PFRecoJet_vx;
      vector<double>                            PFRecoJet_vy;
      vector<double>                            PFRecoJet_vz;

      //Id and other Variables for Calo Pat Jets
      vector<float>                             CaloPatJet_emEnergyFraction;
      vector<float>                             CaloPatJet_emEnergyInEB;
      vector<float>                             CaloPatJet_emEnergyInEE;
      vector<float>                             CaloPatJet_emEnergyInHF;
      vector<float>                             CaloPatJet_HadronicEnergyFraction;
      vector<int>                               CaloPatJet_NConstituents;
      vector<float>                             CaloPatJet_HadronicEnergyInHB;
      vector<float>                             CaloPatJet_HadronicEnergyInHE;
      vector<float>                             CaloPatJet_HadronicEnergyInHF;
      vector<float>                             CaloPatJet_HadronicEnergyInHO;
      vector<float>                             CaloPatJet_maxEnergyInEmTowers;
      vector<float>                             CaloPatJet_maxEnergyInHadTowers;
      vector<int>                               CaloPatJet_nConstituents60E;
      vector<int>                               CaloPatJet_nConstituents90E;
      vector<float>                             CaloPatJet_towersArea;
      vector<int>                               CaloPatJet_nTowers;
      vector<float>                             CaloPatJet_fEinHottestHPD;  
      vector<float>                             CaloPatJet_fEinHottestRBX;
      vector<int>                               CaloPatJet_nHitsCarrying90E;
      vector<int>                               CaloPatJet_nHitsinTowersCarrying90E;
      vector<float>                             CaloPatJet_RHF;

      //Id and other Variables for PF Pat Jets
      vector<float>                             PFPatJet_ChargedEmEnergy;
      vector<float>                             PFPatJet_ChargedEmEnergyFrac;
      vector<float>                             PFPatJet_ChargedHadEnergy;
      vector<float>                             PFPatJet_ChargedHadEnergyFrac;
      vector<int>                               PFPatJet_ChargedHadMult;
      vector<int>                               PFPatJet_ChargedMult;
      vector<int>                               PFPatJet_NConstituents;
      vector<float>                             PFPatJet_HFEMEnergy;
      vector<float>                             PFPatJet_HFEMEnergyFrac;
      vector<int>                               PFPatJet_HFEMMult;
      vector<float>                             PFPatJet_HFHadEnergy;
      vector<float>                             PFPatJet_HFHadEnergyFrac;
      vector<int>                               PFPatJet_HFHadMult;
      vector<float>                             PFPatJet_NeutralEmEnergy;
      vector<float>                             PFPatJet_NeutralEmEnergyFrac;
      vector<float>                             PFPatJet_NeutralHadEnergy;
      vector<float>                             PFPatJet_NeutralHadEnergyFrac;
      vector<int>                               PFPatJet_NeutralHadMult;
      vector<int>                               PFPatJet_NeutralMult;

      //Id and other Variables for Calo Reco Jets
      vector<float>                             CaloRecoJet_emEnergyFraction;
      vector<float>                             CaloRecoJet_emEnergyInEB;
      vector<float>                             CaloRecoJet_emEnergyInEE;
      vector<float>                             CaloRecoJet_emEnergyInHF;
      vector<float>                             CaloRecoJet_HadronicEnergyFraction;
      vector<int>                               CaloRecoJet_NConstituents;
      vector<float>                             CaloRecoJet_HadronicEnergyInHB;
      vector<float>                             CaloRecoJet_HadronicEnergyInHE;
      vector<float>                             CaloRecoJet_HadronicEnergyInHF;
      vector<float>                             CaloRecoJet_HadronicEnergyInHO;
      vector<float>                             CaloRecoJet_maxEnergyInEmTowers;
      vector<float>                             CaloRecoJet_maxEnergyInHadTowers;
      vector<int>                               CaloRecoJet_nConstituents60E;
      vector<int>                               CaloRecoJet_nConstituents90E;
      vector<float>                             CaloRecoJet_towersArea;

      //Id and other Variables for PF Reco Jets
      vector<float>                             PFRecoJet_ChargedEmEnergy;
      vector<float>                             PFRecoJet_ChargedEmEnergyFrac;
      vector<float>                             PFRecoJet_ChargedHadEnergy;
      vector<float>                             PFRecoJet_ChargedHadEnergyFrac;
      vector<int>                               PFRecoJet_ChargedHadMult;
      vector<int>                               PFRecoJet_ChargedMult;
      vector<int>                               PFRecoJet_NConstituents;
      vector<float>                             PFRecoJet_HFEMEnergy;
      vector<float>                             PFRecoJet_HFEMEnergyFrac;
      vector<int>                               PFRecoJet_HFEMMult;
      vector<float>                             PFRecoJet_HFHadEnergy;
      vector<float>                             PFRecoJet_HFHadEnergyFrac;
      vector<int>                               PFRecoJet_HFHadMult;
      vector<float>                             PFRecoJet_NeutralEmEnergy;
      vector<float>                             PFRecoJet_NeutralEmEnergyFrac;
      vector<float>                             PFRecoJet_NeutralHadEnergy;
      vector<float>                             PFRecoJet_NeutralHadEnergyFrac;
      vector<int>                               PFRecoJet_NeutralHadMult;
      vector<int>                               PFRecoJet_NeutralMult;

      //JEC Uncertainity variables
      vector<double>                             CaloPatJet_jecUncertainity;
      vector<double>                             PFPatJet_jecUncertainity;
      vector<double>                             CaloRecoJet_jecUncertainity;
      vector<double>                             PFRecoJet_jecUncertainity;

      //PileUp Jet Id variables for PFPat Jets (PU Jet Id stored only for PF PAT)
      vector<float>                              PFPatJet_puJetIdCutBased_MVA;
      vector<float>                              PFPatJet_puJetIdSimple_MVA;
      vector<float>                              PFPatJet_puJetIdFull_MVA;

      vector<bool>                               PFPatJet_PassPUJetIdCutBased_loose;
      vector<bool>                               PFPatJet_PassPUJetIdCutBased_medium;
      vector<bool>                               PFPatJet_PassPUJetIdCutBased_tight;

      vector<bool>                               PFPatJet_PassPUJetIdSimple_loose;
      vector<bool>                               PFPatJet_PassPUJetIdSimple_medium;
      vector<bool>                               PFPatJet_PassPUJetIdSimple_tight;

      vector<bool>                               PFPatJet_PassPUJetIdFull_loose;
      vector<bool>                               PFPatJet_PassPUJetIdFull_medium;
      vector<bool>                               PFPatJet_PassPUJetIdFull_tight;

      //B-tagging variables for CaloPatJets
      vector<float>                              CaloPatJet_BJetDiscrByTrackCountingHighEff;
      vector<float>                              CaloPatJet_BJetDiscrByTrackCountingHighPur;
      vector<float>                              CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;
      vector<float>                              CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;
      vector<float>                              CaloPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;
      vector<float>                              CaloPatJet_BJetDiscrByCombinedSecondaryVertexMVA;
      vector<float>                              CaloPatJet_BJetDiscrByjetProbabilityBJetTags;
      vector<float>                              CaloPatJet_BJetDiscrByjetBProbabilityBJetTags;
      vector<float>                              CaloPatJet_BJetDiscrBySoftElectronBJetTags;
      vector<float>                              CaloPatJet_BJetDiscrBySoftMuonBJetTags;
      vector<int>                                CaloPatJet_JetPartonFlavor;

      //B-tagging variables for PFPatJets
      vector<float>                              PFPatJet_BJetDiscrByTrackCountingHighEff;
      vector<float>                              PFPatJet_BJetDiscrByTrackCountingHighPur;
      vector<float>                              PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;
      vector<float>                              PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;
      vector<float>                              PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;
      vector<float>                              PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA;
      vector<float>                              PFPatJet_BJetDiscrByjetProbabilityBJetTags;
      vector<float>                              PFPatJet_BJetDiscrByjetBProbabilityBJetTags;
      vector<float>                              PFPatJet_BJetDiscrBySoftElectronBJetTags;
      vector<float>                              PFPatJet_BJetDiscrBySoftMuonBJetTags;
      vector<int>                                PFPatJet_JetPartonFlavor;

      //B-tagging variables for CaloRecoJets
      vector<float>                              CaloRecoJet_BJetDiscrByTrackCountingHighEff;
      vector<float>                              CaloRecoJet_BJetDiscrByTrackCountingHighPur;
      vector<float>                              CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;
      vector<float>                              CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;
      vector<float>                              CaloRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;
      vector<float>                              CaloRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;
      vector<float>                              CaloRecoJet_BJetDiscrByjetProbabilityBJetTags;
      vector<float>                              CaloRecoJet_BJetDiscrByjetBProbabilityBJetTags;
      vector<float>                              CaloRecoJet_BJetDiscrBySoftElectronBJetTags;
      vector<float>                              CaloRecoJet_BJetDiscrBySoftMuonBJetTags;

      //B-tagging variables for PFRecoJet
      vector<float>                              PFRecoJet_BJetDiscrByTrackCountingHighEff;
      vector<float>                              PFRecoJet_BJetDiscrByTrackCountingHighPur;
      vector<float>                              PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;
      vector<float>                              PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;
      vector<float>                              PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;
      vector<float>                              PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;
      vector<float>                              PFRecoJet_BJetDiscrByjetProbabilityBJetTags;
      vector<float>                              PFRecoJet_BJetDiscrByjetBProbabilityBJetTags;
      vector<float>                              PFRecoJet_BJetDiscrBySoftElectronBJetTags;
      vector<float>                              PFRecoJet_BJetDiscrBySoftMuonBJetTags;

      //------------------Vertex Info Variables------------------

      bool                                       runVertex_;
      int                                        Vertex_n;
      vector<double>                             Vertex_x;
      vector<double>                             Vertex_y;
      vector<double>                             Vertex_z;
      vector<double>                             Vertex_chi2;
      vector<double>                             Vertex_nchi2;  //normalized chi2(chi2/ndof)
      vector<double>                             Vertex_ndof;
      vector<int>                                Vertex_tracksSize;
      vector<bool>                               Vertex_isFake;
      vector<bool>                               Vertex_isValid;
      vector<double>                             Vertex_d0;

      //-------------------Track Info Variables------------------

      bool                                      runTracks_;
      int                                       Tracks_n;
      vector<double>                            Track_pt;
      vector<double>                            Track_px;
      vector<double>                            Track_py;
      vector<double>                            Track_pz;
      vector<double>                            Track_vx;
      vector<double>                            Track_vy;
      vector<double>                            Track_vz;
      vector<double>                            Track_eta;
      vector<double>                            Track_phi;
      vector<double>                            Track_theta;
      vector<double>                            Track_chi2;

      //--------------------Scrapping Info-----------------------

      bool                                      runScraping_;
      bool                                      IsScrapingEvent_;
      float                                     Scraping_FractionOfGoodTracks;

      //--------------------Trigger Info-------------------------
      //+++++++++++++++++++++++ HLT +++++++++++++++++++++++++++++                

      bool                                      runHLT_;
      HLTConfigProvider                         hltConfig_;
      vector<std::string>                       HLT_Photon_triggers;
      vector<int>                               HLT_Photon_trig_prescales;
      vector<bool>                              HLT_Photon_ifTriggerPassed;
      int                                       HLT_Photon_nTriggers; 
      vector<int>                               HLT_Photon_triggerIndex;
      vector<int>                               HLT_Photon_nFilters; //No. of saveTag = true filters in each trigger
      std::vector< std::vector< std::string > >              HLT_Photon_FilterNames; //FilterNames for each of trigger
      vector<vector<double> >                   HLT_Photon_FilterObjects_pt; //Obj pt for each filter
      vector<vector<double> >                   HLT_Photon_FilterObjects_eta;
      vector<vector<double> >                   HLT_Photon_FilterObjects_phi;














      //+++++++++++++++++++++++ L1 +++++++++++++++++++++++++++++++
      bool runL1;




























      //EBRecHit Collection variables
      int Photon_EBRecHit_size;
      float EBRecHit_eta[10000];
      float EBRecHit_phi[10000];
      int EBRecHit_ieta[10000];
      int EBRecHit_iphi[10000];
      float EBRecHit_e[10000];
      float EBRecHit_et[10000];
      int EBRecHit_flag[10000];
      float EBRecHit_time[10000];



      //EErecHit variables;
      int EERecHit_size;
      float EERecHit_eta[10000];
      float EERecHit_phi[10000];
      int EERecHit_ieta[10000];
      int EERecHit_iphi[10000];
      float EERecHit_e[10000];
      float EERecHit_et[10000];
      int EERecHit_flag[10000];
      float EERecHit_time[10000];








      int  pho_mGenmompdgId[200][100];
      int  pho_mGenpdgId[200];
      int  pho_nummoth[200];








      //----------- INPUT TAGS------------------------------
      edm::InputTag MCpileupLabel_;
      edm::InputTag patPhoLabel_;
      edm::InputTag recoPhoLabel_;
      edm::InputTag pfCandidateLabel_;
      edm::InputTag vertexLabel_;
      edm::InputTag BSLabel_;
      edm::InputTag ConvPhoLabel_;
      edm::InputTag GsfEleLabel_;
      edm::InputTag BarrelrechitLabel_;
      edm::InputTag EndcaprechitLabel_;
      edm::InputTag CaloPatJetLabel_;
      edm::InputTag PFPatJetLabel_;
      edm::InputTag CaloRecoJetLabel_;
      edm::InputTag PFRecoJetLabel_;
      edm::InputTag CaloReco_TrkCountingHighEffBJetTagsLabel_;
      edm::InputTag CaloReco_TrkCountingHighPurBJetTagsLabel_;
      edm::InputTag CaloReco_SimpleSecVtxHighEffBJetTagsLabel_;
      edm::InputTag CaloReco_SimpleSecVtxHighPurBJetTagsLabel_;
      edm::InputTag CaloReco_CombinedSecVtxBJetTagsLabel_;
      edm::InputTag CaloReco_CombinedSecVtxMVABJetTagsLabel_;
      edm::InputTag CaloReco_JetProbBJetTagsLabel_;
      edm::InputTag CaloReco_JetBProbBJetTagsLabel_;
      edm::InputTag CaloReco_SoftElectronBJetTagsLabel_;
      edm::InputTag CaloReco_SoftMuonBJetTagsLabel_;
      edm::InputTag PFReco_TrkCountingHighEffBJetTagsLabel_;
      edm::InputTag PFReco_TrkCountingHighPurBJetTagsLabel_;
      edm::InputTag PFReco_SimpleSecVtxHighEffBJetTagsLabel_;
      edm::InputTag PFReco_SimpleSecVtxHighPurBJetTagsLabel_;
      edm::InputTag PFReco_CombinedSecVtxBJetTagsLabel_;
      edm::InputTag PFReco_CombinedSecVtxMVABJetTagsLabel_;
      edm::InputTag PFReco_JetProbBJetTagsLabel_;
      edm::InputTag PFReco_JetBProbBJetTagsLabel_;
      edm::InputTag PFReco_SoftElectronBJetTagsLabel_;
      edm::InputTag PFReco_SoftMuonBJetTagsLabel_;
      edm::InputTag TracksLabel_;
      edm::InputTag triggerEventLabel_;
      edm::InputTag triggerResultsLabel_;


};

struct CrystalInfo{
  int rawId;
  int ieta;
  int iphi;
  int ix;
  int iy;
  double energy;
  double time;
  double timeErr;
  int recoFlag;
  CrystalInfo():
    rawId(-999),
    ieta(-999),
    iphi(-999),
    energy(-999),
    time(-999),
    timeErr(-999),
    recoFlag(-999)
  {}
};


void Analyzer::ClearVectors(){

  //Kinematic Variables
  Photon_E.clear();
  Photon_et.clear();
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_theta.clear();
  Photon_px.clear();
  Photon_py.clear();
  Photon_pz.clear();
  Photon_vx.clear();
  Photon_vy.clear();
  Photon_vz.clear();

  //Shower Shape Variables
  Photon_r9.clear();
  Photon_maxEnergyXtal.clear();
  Photon_e1x5.clear();
  Photon_e2x5.clear();
  Photon_e3x3.clear();
  Photon_e5x5.clear();
  Photon_r1x5.clear();
  Photon_r2x5.clear();
  Photon_SigmaEtaEta.clear();
  Photon_SigmaIEtaIEta.clear();
  Photon_SigmaEtaPhi.clear();
  Photon_SigmaIEtaIPhi.clear();
  Photon_SigmaPhiPhi.clear();
  Photon_SigmaIPhiIPhi.clear();
  Photon_roundness.clear();   
  Photon_angle.clear();
  Photon_swissCross.clear();
  Photon_s9.clear();
  Photon_e4Overe1.clear();
  Photon_e6Overe2.clear();
  Photon_e2Overe9.clear(); 
  Photon_rookFraction.clear();


  //Fiducial flags Variables
  Photon_isEB.clear();
  Photon_isEE.clear();
  Photon_isEBGap.clear();
  Photon_isEEGap.clear();
  Photon_isEBEEGap.clear();
 
  //Detector Isolation Variables
  Photon_ecalRecHitSumEtConeDR03.clear();
  Photon_hcalTowerSumEtConeDR03.clear();
  Photon_hcalDepth1TowerSumEtConeDR03.clear();
  Photon_hcalDepth2TowerSumEtConeDR03.clear();
  Photon_trkSumPtSolidConeDR03.clear();
  Photon_trkSumPtHollowConeDR03.clear();
  Photon_nTrkSolidConeDR03.clear();
  Photon_nTrkHollowConeDR03.clear();
  Photon_ecalRecHitSumEtConeDR04.clear();
  Photon_hcalTowerSumEtConeDR04.clear();
  Photon_hcalDepth1TowerSumEtConeDR04.clear();
  Photon_hcalDepth2TowerSumEtConeDR04.clear();
  Photon_trkSumPtSolidConeDR04.clear();
  Photon_trkSumPtHollowConeDR04.clear();
  Photon_nTrkSolidConeDR04.clear();
  Photon_nTrkHollowConeDR04.clear();

  //Photon Identification Variables
  Photon_HoE.clear();
  Photon_SingleTowerHoE.clear();
  Photon_hasConvTrk.clear();
  Photon_hasPixelSeed.clear();
  passedConvSafeElectronVeto.clear();  

  //Photon Super Cluster Variables
  Photon_SC_nOfBasicClusters.clear();
  Photon_SC_rawEnergy.clear();
  Photon_SC_preShowerEnergy.clear();
  Photon_SC_energy.clear();
  Photon_SC_eta.clear();
  Photon_SC_phi.clear();
  Photon_SC_x.clear();
  Photon_SC_y.clear();
  Photon_SC_z.clear();
  Photon_SC_etaWidth.clear();
  Photon_SC_phiWidth.clear();
      
  //Photon MIP Variables
  Photon_mipChi2.clear();
  Photon_mipTotEnergy.clear();
  Photon_mipSlope.clear();
  Photon_mipIntercept.clear();
  Photon_mipNhitCone.clear();
  Photon_mipIsHalo.clear();

  //Converted Photon Variables
  Photon_nConvTracks.clear();
  Photon_isConverted.clear();
  Photon_pairInvariantMass.clear();
  Photon_pairCotThetaSeparation.clear();
  Photon_pairMomentum_x.clear();
  Photon_pairMomentum_y.clear();
  Photon_pairMomentum_z.clear();
  Photon_EoverP.clear();
  Photon_conv_vx.clear();
  Photon_conv_vy.clear();
  Photon_conv_vz.clear();
  Photon_zOfPrimaryVtxFromTrks.clear();
  Photon_distOfMinimumApproach.clear();
  Photon_dPhiTracksAtVtx.clear();
  Photon_dPhiTracksAtEcal.clear();
  Photon_dEtaTracksAtEcal.clear();
                  
  //Photon Crystal Variables
  Photon_nCrystals.clear();     
  Photon_xtal_timing.clear();
  Photon_xtal_timeErr.clear();
  Photon_avgTimeAllxtals.clear();
  Photon_xtal_energy.clear();
  Photon_xtal_EBieta.clear();
  Photon_xtal_EBiphi.clear();
  Photon_xtal_EBrecoFlag.clear();

  //PFIsolation Variables
  PFIsoPhoton03.clear();
  PFIsoNeutral03.clear();
  PFIsoCharged03.clear();
  PFIsoSum03.clear();
  PFIsoChargedWorstvtx03.clear();

  //CLEARING JET VARIABLES

  //Kinematic Variables for Calo Pat Jets
  CaloPatJet_E.clear();
  CaloPatJet_et.clear();
  CaloPatJet_pt.clear();
  CaloPatJet_eta.clear();
  CaloPatJet_phi.clear();
  CaloPatJet_theta.clear();
  CaloPatJet_px.clear();
  CaloPatJet_py.clear();
  CaloPatJet_pz.clear();
  CaloPatJet_vx.clear();
  CaloPatJet_vy.clear();
  CaloPatJet_vz.clear();

  //Kinematic Variables for PF Pat Jets
  PFPatJet_E.clear();
  PFPatJet_et.clear();
  PFPatJet_pt.clear();
  PFPatJet_eta.clear();
  PFPatJet_phi.clear();
  PFPatJet_theta.clear();
  PFPatJet_px.clear();
  PFPatJet_py.clear();
  PFPatJet_pz.clear();
  PFPatJet_vx.clear();
  PFPatJet_vy.clear();
  PFPatJet_vz.clear();

  //Kinematic Variables for Reco Calo Jets
  CaloRecoJet_E.clear();
  CaloRecoJet_et.clear();
  CaloRecoJet_pt.clear();
  CaloRecoJet_eta.clear();
  CaloRecoJet_phi.clear();
  CaloRecoJet_theta.clear();
  CaloRecoJet_px.clear();
  CaloRecoJet_py.clear();
  CaloRecoJet_pz.clear();
  CaloRecoJet_vx.clear();
  CaloRecoJet_vy.clear();
  CaloRecoJet_vz.clear();

  //Kinematic Variables for Reco PF Jets
  PFRecoJet_E.clear();
  PFRecoJet_et.clear();
  PFRecoJet_pt.clear();
  PFRecoJet_eta.clear();
  PFRecoJet_phi.clear();
  PFRecoJet_theta.clear();
  PFRecoJet_px.clear();
  PFRecoJet_py.clear();
  PFRecoJet_pz.clear();
  PFRecoJet_vx.clear();
  PFRecoJet_vy.clear();
  PFRecoJet_vz.clear();

  //Id and other Variables for Calo Pat Jets
  CaloPatJet_emEnergyFraction.clear();
  CaloPatJet_emEnergyInEB.clear();
  CaloPatJet_emEnergyInEE.clear();
  CaloPatJet_emEnergyInHF.clear();
  CaloPatJet_HadronicEnergyFraction.clear();
  CaloPatJet_NConstituents.clear();
  CaloPatJet_HadronicEnergyInHB.clear();
  CaloPatJet_HadronicEnergyInHE.clear();
  CaloPatJet_HadronicEnergyInHF.clear();
  CaloPatJet_HadronicEnergyInHO.clear();
  CaloPatJet_maxEnergyInEmTowers.clear();
  CaloPatJet_maxEnergyInHadTowers.clear();
  CaloPatJet_nConstituents60E.clear();
  CaloPatJet_nConstituents90E.clear();
  CaloPatJet_towersArea.clear();
  CaloPatJet_nTowers.clear();
  CaloPatJet_fEinHottestHPD.clear(); 
  CaloPatJet_fEinHottestRBX.clear();
  CaloPatJet_nHitsCarrying90E.clear();
  CaloPatJet_nHitsinTowersCarrying90E.clear();
  CaloPatJet_RHF.clear();

  //Id and other Variables for PF Pat Jets
  PFPatJet_ChargedEmEnergy.clear();
  PFPatJet_ChargedEmEnergyFrac.clear();
  PFPatJet_ChargedHadEnergy.clear();
  PFPatJet_ChargedHadEnergyFrac.clear();
  PFPatJet_ChargedHadMult.clear();
  PFPatJet_ChargedMult.clear();
  PFPatJet_NConstituents.clear();
  PFPatJet_HFEMEnergy.clear();
  PFPatJet_HFEMEnergyFrac.clear();
  PFPatJet_HFEMMult.clear();
  PFPatJet_HFHadEnergy.clear();
  PFPatJet_HFHadEnergyFrac.clear();
  PFPatJet_HFHadMult.clear();
  PFPatJet_NeutralEmEnergy.clear();
  PFPatJet_NeutralEmEnergyFrac.clear();
  PFPatJet_NeutralHadEnergy.clear();
  PFPatJet_NeutralHadEnergyFrac.clear();
  PFPatJet_NeutralHadMult.clear();
  PFPatJet_NeutralMult.clear();

  //Id and other Variables for Calo Reco Jets
  CaloRecoJet_emEnergyFraction.clear();
  CaloRecoJet_emEnergyInEB.clear();
  CaloRecoJet_emEnergyInEE.clear();
  CaloRecoJet_emEnergyInHF.clear();
  CaloRecoJet_HadronicEnergyFraction.clear();
  CaloRecoJet_NConstituents.clear();
  CaloRecoJet_HadronicEnergyInHB.clear();
  CaloRecoJet_HadronicEnergyInHE.clear();
  CaloRecoJet_HadronicEnergyInHF.clear();
  CaloRecoJet_HadronicEnergyInHO.clear();
  CaloRecoJet_maxEnergyInEmTowers.clear();
  CaloRecoJet_maxEnergyInHadTowers.clear();
  CaloRecoJet_nConstituents60E.clear();
  CaloRecoJet_nConstituents90E.clear();
  CaloRecoJet_towersArea.clear();

  //Id and other Variables for PF Reco Jets
  PFRecoJet_ChargedEmEnergy.clear();
  PFRecoJet_ChargedEmEnergyFrac.clear();
  PFRecoJet_ChargedHadEnergy.clear();
  PFRecoJet_ChargedHadEnergyFrac.clear();
  PFRecoJet_ChargedHadMult.clear();
  PFRecoJet_ChargedMult.clear();
  PFRecoJet_NConstituents.clear();
  PFRecoJet_HFEMEnergy.clear();
  PFRecoJet_HFEMEnergyFrac.clear();
  PFRecoJet_HFEMMult.clear();
  PFRecoJet_HFHadEnergy.clear();
  PFRecoJet_HFHadEnergyFrac.clear();
  PFRecoJet_HFHadMult.clear();
  PFRecoJet_NeutralEmEnergy.clear();
  PFRecoJet_NeutralEmEnergyFrac.clear();
  PFRecoJet_NeutralHadEnergy.clear();
  PFRecoJet_NeutralHadEnergyFrac.clear();
  PFRecoJet_NeutralHadMult.clear();
  PFRecoJet_NeutralMult.clear();

  //JEC Uncertainity variables
  CaloPatJet_jecUncertainity.clear();
  PFPatJet_jecUncertainity.clear();
  CaloRecoJet_jecUncertainity.clear();
  PFRecoJet_jecUncertainity.clear();

  //PileUp Jet Id variables for PFPat Jets
  PFPatJet_puJetIdCutBased_MVA.clear();
  PFPatJet_puJetIdSimple_MVA.clear();
  PFPatJet_puJetIdFull_MVA.clear();
  PFPatJet_PassPUJetIdCutBased_loose.clear();
  PFPatJet_PassPUJetIdCutBased_medium.clear();
  PFPatJet_PassPUJetIdCutBased_tight.clear();
  PFPatJet_PassPUJetIdSimple_loose.clear();
  PFPatJet_PassPUJetIdSimple_medium.clear();
  PFPatJet_PassPUJetIdSimple_tight.clear();
  PFPatJet_PassPUJetIdFull_loose.clear();
  PFPatJet_PassPUJetIdFull_medium.clear();
  PFPatJet_PassPUJetIdFull_tight.clear();

  //B-tagging variables for CaloPatJets
  CaloPatJet_BJetDiscrByTrackCountingHighEff.clear();
  CaloPatJet_BJetDiscrByTrackCountingHighPur.clear();
  CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighEff.clear();
  CaloPatJet_BJetDiscrBySimpleSecondaryVertexHighPur.clear();
  CaloPatJet_BJetDiscrByCombinedSecondaryVertexHighEff.clear();
  CaloPatJet_BJetDiscrByCombinedSecondaryVertexMVA.clear();
  CaloPatJet_BJetDiscrByjetProbabilityBJetTags.clear();
  CaloPatJet_BJetDiscrByjetBProbabilityBJetTags.clear();
  CaloPatJet_BJetDiscrBySoftElectronBJetTags.clear();
  CaloPatJet_BJetDiscrBySoftMuonBJetTags.clear();
  CaloPatJet_JetPartonFlavor.clear();

  //B-tagging variables for PFPatJets
  PFPatJet_BJetDiscrByTrackCountingHighEff.clear();
  PFPatJet_BJetDiscrByTrackCountingHighPur.clear();
  PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff.clear();
  PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur.clear();
  PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff.clear();
  PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA.clear();
  PFPatJet_BJetDiscrByjetProbabilityBJetTags.clear();
  PFPatJet_BJetDiscrByjetBProbabilityBJetTags.clear();
  PFPatJet_BJetDiscrBySoftElectronBJetTags.clear();
  PFPatJet_BJetDiscrBySoftMuonBJetTags.clear();
  PFPatJet_JetPartonFlavor.clear();

  //B-tagging variables for CaloRecoJets
  CaloRecoJet_BJetDiscrByTrackCountingHighEff.clear();
  CaloRecoJet_BJetDiscrByTrackCountingHighPur.clear();
  CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff.clear();
  CaloRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur.clear();
  CaloRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff.clear();
  CaloRecoJet_BJetDiscrByCombinedSecondaryVertexMVA.clear();
  CaloRecoJet_BJetDiscrByjetProbabilityBJetTags.clear();
  CaloRecoJet_BJetDiscrByjetBProbabilityBJetTags.clear();
  CaloRecoJet_BJetDiscrBySoftElectronBJetTags.clear();
  CaloRecoJet_BJetDiscrBySoftMuonBJetTags.clear();

  //B-tagging variables for PFRecoJet
  PFRecoJet_BJetDiscrByTrackCountingHighEff.clear();
  PFRecoJet_BJetDiscrByTrackCountingHighPur.clear();
  PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff.clear();
  PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur.clear();
  PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff.clear();
  PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA.clear();
  PFRecoJet_BJetDiscrByjetProbabilityBJetTags.clear();
  PFRecoJet_BJetDiscrByjetBProbabilityBJetTags.clear();
  PFRecoJet_BJetDiscrBySoftElectronBJetTags.clear();
  PFRecoJet_BJetDiscrBySoftMuonBJetTags.clear();

  //Vertex Info Vaiables
  Vertex_x.clear();
  Vertex_y.clear();
  Vertex_z.clear();
  Vertex_chi2.clear();
  Vertex_nchi2.clear();
  Vertex_ndof.clear();
  Vertex_tracksSize.clear();
  Vertex_isFake.clear();
  Vertex_isValid.clear();
  Vertex_d0.clear();

  //Track Info Vaiables
  Track_pt.clear();
  Track_px.clear();
  Track_py.clear();
  Track_pz.clear();
  Track_vx.clear();
  Track_vy.clear();
  Track_vz.clear();
  Track_eta.clear();
  Track_phi.clear();
  Track_theta.clear();
  Track_chi2.clear();

  //Trigger Info Variables
  HLT_Photon_triggers.clear();
  HLT_Photon_trig_prescales.clear();
  HLT_Photon_ifTriggerPassed.clear();
  HLT_Photon_triggerIndex.clear();
  HLT_Photon_nFilters.clear();
  HLT_Photon_FilterNames.clear();
  HLT_Photon_FilterObjects_pt.clear();
  HLT_Photon_FilterObjects_eta.clear();
  HLT_Photon_FilterObjects_phi.clear();



}

//Taken from RecoLocalCalo/EcalRecAlgos/src/EcalCleaningAlgo.cc
float Analyzer::GetE4OverE1(const DetId& id, const EcalRecHitCollection& rhs){

   float s4 = 0;
   float e1 = recHitE( id, rhs, false );

   if ( e1 == 0 ) return 0;
   const std::vector<DetId>& neighs =  neighbours(id);
   for (size_t i=0; i<neighs.size(); ++i){
     // avoid hits out of time when making s4
     s4+=recHitE(neighs[i],rhs, true);
   }
   return s4 / e1;

}


// Taken from RecoLocalCalo/EcalRecAlgos/src/EcalCleaningAlgo.cc
float Analyzer::GetE6OverE2(const DetId& id, const EcalRecHitCollection& rhs){

   float s4_1 = 0;
   float s4_2 = 0;
   float e1 = recHitE( id, rhs , false );

   float maxene=0;
   DetId maxid;

   if ( e1 == 0 ) return 0;

   const std::vector<DetId>& neighs =  neighbours(id);

   // find the most energetic neighbour ignoring time info
   for (size_t i=0; i<neighs.size(); ++i){
     float ene = recHitE(neighs[i],rhs,false);
     if (ene>maxene) {
       maxene=ene;
       maxid = neighs[i];
     }
   }

   float e2=maxene;

   s4_1 = GetE4OverE1(id,rhs)* e1;
   s4_2 = GetE4OverE1(maxid,rhs)* e2;

   return (s4_1 + s4_2) / (e1+e2) -1. ;

}

float Analyzer::GetE2OverE9(const DetId& id, const EcalRecHitCollection& rhs){
  float e2e9 = 0;
  if ( id.subdetId() == EcalBarrel ) {
    EBDetId ebId( id );
    float e1  = recHitE( id, rhs );
    float e2  = 0;
    float s9  = 0;

    for ( int deta = -1; deta <= +1; deta++ )
      {
	for ( int dphi = -1; dphi <= +1; dphi++ )
	  {
	    float etmp=recHitE( id, rhs, deta, dphi );
	    s9 += etmp;
	    if (etmp > e2 && !(deta==0 && dphi==0)) {
	      e2=etmp;
	    }
	  }
      }
    float s2=e1+e2;
    if (s9!=0) e2e9= s2/s9;
  }
  return e2e9;
}

float GetE2OverE9(const DetId& id, const EcalRecHitCollection& rhs, float recHitEtThreshold = 10.0,
		  float recHitEtThreshold2 = 1.0, bool avoidIeta85 = false, bool KillSecondHit = true)
  // ====================================================================================
  // taken from CMSSW/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc CMSSW_3_9_0_pre5
  // ====================================================================================
{

  // compute e2overe9
  //   | | | |
  //   +-+-+-+
  //   | |1|2|
  //   +-+-+-+
  //   | | | |
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
  //   rechit 1 must have E_t > recHitEtThreshold
  //   rechit 2 must have E_t > recHitEtThreshold2
  //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
  //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
  //   provided that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0

  if ( id.subdetId() == EcalBarrel ) {
    EBDetId ebId( id );
    // avoid recHits at |eta|=85 where one side of the neighbours is missing
    if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;
    // select recHits with Et above recHitEtThreshold
    float e1 = recHitE( id, rhs );
    float ete1=recHitApproxEt( id, rhs );
    // check that rechit E_t is above threshold
    if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;
    if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;

    float e2=-1;
    float ete2=0;
    float s9 = 0;
    // coordinates of 2nd hit relative to central hit
    int e2eta=0;
    int e2phi=0;

    // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
    for ( int deta = -1; deta <= +1; ++deta ) {
      for ( int dphi = -1; dphi <= +1; ++dphi ) {
	// compute 3x3 energy
	float etmp=recHitE( id, rhs, deta, dphi );
	s9 += etmp;
	EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
	float eapproxet=recHitApproxEt( idtmp, rhs );
	// remember 2nd highest energy deposit (above threshold) in 3x3 array
	if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	  e2=etmp;
	  ete2=eapproxet;
	  e2eta=deta;
	  e2phi=dphi;
	}
      }
    }
    if ( e1 == 0 )  return 0;
    // return 0 if 2nd hit is below threshold
    if ( e2 == -1 ) return 0;
    // compute e2/e9 centered around 1st hit
    float e2nd=e1+e2;
    float e2e9=0;

    if (s9!=0) e2e9=e2nd/s9;
    // if central hit has higher energy than 2nd hit
    //  return e2/e9 if 1st hit is above E_t threshold
    if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;
    // if second hit has higher energy than 1st hit
    if ( e2 > e1 ) {
      // return 0 if user does not want to flag 2nd hit, or
      // hits are below E_t thresholds - note here we
      // now assume the 2nd hit to be the leading hit.

      if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
        return 0;
      }
      else {
        // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2
        float s92nd=0;
        float e2nd_prime=0;
        int e2prime_eta=0;
        int e2prime_phi=0;

        EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

        for ( int deta = -1; deta <= +1; ++deta ) {
	  for ( int dphi = -1; dphi <= +1; ++dphi ) {

	    // compute 3x3 energy
	    float etmp=recHitE( secondid, rhs, deta, dphi );
	    s92nd += etmp;

	    if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	      e2nd_prime=etmp;
	      e2prime_eta=deta;
	      e2prime_phi=dphi;
	    }
	  }
        }
        // if highest energy hit around E2 is not the same as the input hit, return 0;
        if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	  {
            return 0;
	  }
        // compute E2/E9 around second hit
        float e2e9_2=0;
        if (s92nd!=0) e2e9_2=e2nd/s92nd;
        //   return the value of E2/E9 calculated around 2nd hit
	return e2e9_2;
      }
    }
  } else if ( id.subdetId() == EcalEndcap ) {
    // only used for EB at the moment
    return 0;
  }
  return 0;
}

float Analyzer::GetrookFraction(const DetId& id, const EcalRecHitCollection& rhs){
 
  float e1 = recHitE( id, rhs , false );

  float maxene=0;
  DetId maxid;

  if ( e1 == 0 ) return 0;

  const std::vector<DetId>& neighs =  neighbours(id);

  // find the most energetic neighbour ignoring time info
  for (size_t i=0; i<neighs.size(); ++i){
    float ene = recHitE(neighs[i],rhs,false);
    if (ene>maxene) {
      maxene=ene;
      maxid = neighs[i];
    }
  }

  float e2=maxene;

  return e2/e1;
}

float Analyzer::GetMatchedBTagDiscrValueToCaloJet(const reco::JetTagCollection & bTags, const reco::CaloJet & jet){
  float discrValue = -999;
  for(reco::JetTagCollection::const_iterator iTag_itr = bTags.begin(); iTag_itr != bTags.end(); iTag_itr++){
    const reco::JetTag* iTag = &(*iTag_itr);
    if(sqrt(pow(iTag->first->eta() - jet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),jet.phi()), 2)) == 0.){
      discrValue = iTag->second;
      break;
    }
  }
  return discrValue;
}

float Analyzer::GetMatchedBTagDiscrValueToPFJet(const reco::JetTagCollection & bTags, const reco::PFJet & jet){
  float discrValue = -999;
  for(reco::JetTagCollection::const_iterator iTag_itr = bTags.begin(); iTag_itr != bTags.end(); iTag_itr++){
    const reco::JetTag* iTag = &(*iTag_itr);
    if(sqrt(pow(iTag->first->eta() - jet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),jet.phi()), 2)) == 0.){
      discrValue = iTag->second;
      break;
    }
  }
  return discrValue;
}


//Class for sorting of Photon Crystal energy
class EnergySortCriterium{
 public:
  bool operator() (CrystalInfo p1,CrystalInfo p2 ){
    return p1.energy > p2.energy;
  }
};

//Class for sorting of CaloRecoJet Pt
class PtSortCriteriumCaloRecoJet{
public:
  bool operator() (reco::CaloJet p1, reco::CaloJet p2){
    return p1.pt() > p2.pt();
  }
};

//Class for sorting of PFRecoJet Pt
class PtSortCriteriumPFRecoJet{
public:
  bool operator() (reco::PFJet p1, reco::PFJet p2){
    return p1.pt() > p2.pt();
  }
};

//Class for sorting of Tracks Pt
class PtSortCriteriumTracks{
public:
  bool operator() (reco::Track p1, reco::Track p2){
    return p1.pt() > p2.pt();
  }
};







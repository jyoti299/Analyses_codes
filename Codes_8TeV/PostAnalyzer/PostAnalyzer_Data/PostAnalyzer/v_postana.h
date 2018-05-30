//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  4 05:07:17 2012 by ROOT version 5.32/00
// from TTree myEvent/a tree with histograms
// found on file: dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data/Run2012A_recover6Aug/AOD_Output_Run2012A_re
//cover6Aug_10_1_0jb.root
//////////////////////////////////////////////////////////

#ifndef PostAnalyzerData_h
#define PostAnalyzerData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TVector.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include <map>
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<bool>+;
#endif

using namespace std;
using namespace ROOT;

// Fixed size dimensions of array or collections stored in the TTree if any.

class PostAnalyzerData {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Variables defined by me

        TFile *f1;

        Bool_t MC, DATA;
        Double_t EvtWeight, OnlyEvtWeight;
        Double_t Lumi;

        std::vector<int> foundPhoton;
        std::vector<int> foundJet;

        Float_t Cvertex_z;
        Float_t Cvertex_ndof;
        Float_t Cvertex_rho;

	Double_t photon_pt_cut;
        Double_t photon_eta_cut;
        Double_t jet_pt_cut;
	Double_t jet_eta_cut;
        Double_t mass_cut;
        Double_t PhotJetDPhi_cut;
	Double_t PhotJetDEta_cut;

        Int_t goodVertex;
        Int_t PC;
        Int_t JC;

        // Variables used for the turn-on curves
        Bool_t passedHLT150;
        Bool_t passed90HLTall;

	// Histos Declarations
        TH1F *h_CutFlowTable;
	TH1F *h_PC, *h_JC;

	TH1F *h_ptPhoton;
        TH1F *h_etaPhoton;
        TH1F *h_ptJet;
        TH1F *h_etaJet;
        TH1F *h_mass_bin25;

        TH1F *h_Photon_SigmaIetaIeta;
        TH1F *h_DR_PhotonJet;
        TH1F *h_dEta;
	TH1F *h_dphi;

        TH1F *h_nPhoton;
        TH1F *h_nJet;

       // Declaration of leaf types
	Int_t           nevents;
        UInt_t          run;
        UInt_t          event;
        UInt_t          luminosityBlock;
	UInt_t          beamCrossing;
        UInt_t          totalIntensityBeam1;
        UInt_t          totalIntensityBeam2;
        Float_t         avgInsDelLumi;
	Float_t         avgInsDelLumiErr;
	Float_t         avgInsRecLumi;
        Float_t         avgInsRecLumiErr;
	Int_t           ntriggers;
	vector<string>  *triggernames;
	vector<int>     *triggerprescales;
        vector<bool>    *ifTriggerpassed;
        vector<float>   *ObjectPt;
        vector<float>   *ObjectEta;
	vector<float>   *ObjectPhi;
        vector<string>  *FilterNames;
        vector<int>     *FilterStartPosition;
	vector<int>     *FilterEndPosition;
        vector<int>     *ObjectStartPosition;
        vector<int>     *ObjectEndPosition;
    	// Vertex variables
	Int_t           Vertex_n;
        Float_t         Vertex_x[200];   //[Vertex_n]
        Float_t         Vertex_y[200];   //[Vertex_n]
	Float_t         Vertex_z[200];   //[Vertex_n]
        Int_t           Vertex_tracksize[200];   //[Vertex_n]
	Int_t           Vertex_ndof[200];   //[Vertex_n]
	Float_t         Vertex_chi2[200];   //[Vertex_n]
	Float_t         Vertex_d0[200];   //[Vertex_n]
        Bool_t          Vertex_isFake[200];   //[Vertex_n]

        //scraping variables
	Bool_t          Scraping_isScrapingEvent;
        Int_t           Scraping_numOfTracks;
        Float_t         Scraping_fractionOfGoodTracks;

	// Track variables
	Int_t           Track_n;
        Float_t         Track_px[1000];   //[Track_n]
	Float_t         Track_py[1000];   //[Track_n]
	Float_t         Track_pz[1000];   //[Track_n]
	Float_t         Track_vx[1000];   //[Track_n]
        Float_t         Track_vy[1000];   //[Track_n]
        Float_t         Track_vz[1000];   //[Track_n]
        Float_t         Track_pt[1000];   //[Track_n]
	Float_t         Track_eta[1000];   //[Track_n]
        Float_t         Track_phi[1000];   //[Track_n]

        //PFJets
        Int_t           pfJet_n;
        Float_t         pfJet_px[200];   //[pfJet_n]
        Float_t         pfJet_py[200];   //[pfJet_n]
        Float_t         pfJet_E[200];   //[pfJet_n]
        Float_t         pfJet_pz[200];   //[pfJet_n]
        Float_t         pfJet_vx[200];   //[pfJet_n]
        Float_t         pfJet_vy[200];   //[pfJet_n]
        Float_t         pfJet_vz[200];   //[pfJet_n]
        Float_t         pfJet_pt[200];   //[pfJet_n]
        Float_t         pfJet_eta[200];   //[pfJet_n]
        Float_t         pfJet_phi[200];   //[pfJet_n]
        Float_t         pfjet_CEF[200];   //[pfJet_n]
        Float_t         pfjet_CHF[200];   //[pfJet_n]
        Float_t         pfjet_NEF[200];   //[pfJet_n]
        Float_t         pfjet_NHF[200];   //[pfJet_n]
        Int_t           pfjet_NCH[200];   //[pfJet_n]
        Float_t         pfjet_HFHAE[200];   //[pfJet_n]
        Float_t         pfjet_HFEME[200];   //[pfJet_n]
        Int_t           pfjet_NConstituents[200];   //[pfJet_n]
        Int_t           pfJet_partonFlavor[200];   //[pfJet_n]
        Int_t           pfJet_partonStatus[200];   //[pfJet_n]

        //PU based Jet Id
        Float_t         pujetIdFull_mva[200];   //[pfJet_n]
        Float_t         pujetIdSimple_mva[200];   //[pfJet_n]
        Float_t         pujetIdCutBased_mva[200];   //[pfJet_n]

        Int_t           pujetIdFull_loose[200];   //[pfJet_n]
        Int_t           pujetIdFull_medium[200];   //[pfJet_n]
        Int_t           pujetIdFull_tight[200];   //[pfJet_n]

	// Int_t           pujetIdSimple_loose[200];   //[pfJet_n]
	// Int_t           pujetIdFull_medium[200];   //[pfJet_n]
        // Int_t           pujetIdFull_tight[200];   //[pfJet_n]

        Int_t           pujetIdSimple_loose[200];   //[pfJet_n]
        Int_t           pujetIdSimple_medium[200];   //[pfJet_n]
        Int_t           pujetIdSimple_tight[200];   //[pfJet_n]

        Int_t           pujetIdCutBased_loose[200];   //[pfJet_n]
        Int_t           pujetIdCutBased_medium[200];   //[pfJet_n]
        Int_t           pujetIdCutBased_tight[200];   //[pfJet_n]

        Float_t         pfjet_TrackCountHiEffBJetTags[200];   //[pfJet_n]
        Float_t         pfjet_TrackCountHiPurBJetTags[200];   //[pfJet_n]
        Float_t         pfjet_SimpleSVHiEffBJetTags[200];   //[pfJet_n]
        Float_t         pfjet_SimpleSVHiPurBJetTags[200];   //[pfJet_n]
        Float_t         pfJet_jecUncer[200];   //[pfJet_n]
        Float_t         pfJet_jecCorr[200];   //[pfJet_n]

        // Some uncorrectd jet information
        Float_t         ucpfJet_px[200];   //[pfJet_n]
        Float_t         ucpfJet_py[200];   //[pfJet_n]
        Float_t         ucpfJet_E[200];   //[pfJet_n]
        Float_t         ucpfJet_pz[200];   //[pfJet_n]
        Float_t         ucpfJet_pt[200];   //[pfJet_n]
        Float_t         ucpfJet_eta[200];   //[pfJet_n]
	Float_t         ucpfJet_phi[200];   //[pfJet_n]

        // Electron Collection
        Int_t           Electron_n;
        Float_t         Electron_px[200];   //[Electron_n]
	Float_t         Electron_py[200];   //[Electron_n]
        Float_t         Electron_pz[200];   //[Electron_n]
        Float_t         Electron_vx[200];   //[Electron_n]
        Float_t         Electron_vy[200];   //[Electron_n]
	Float_t         Electron_vz[200];   //[Electron_n]
        Float_t         Electron_pt[200];   //[Electron_n]
        Float_t         Electron_eta[200];   //[Electron_n]
        Float_t         Electron_phi[200];   //[Electron_n]
	Float_t         Electron_energy[200];   //[Electron_n]
        Float_t         Electron_charge[200];   //[Electron_n]
        Float_t         Electron_trkIso[200];   //[Electron_n]
        Float_t         Electron_ecalIso[200];   //[Electron_n]
	Float_t         Electron_hcalIso[200];   //[Electron_n]
        Float_t         Electron_SigmaIetaIeta[200];   //[Electron_n]
        Float_t         Electron_dEtaIn[200];   //[Electron_n]
        Float_t         Electron_dPhiIn[200];   //[Electron_n]
        Float_t         Electron_HoE[200];   //[Electron_n]
        Float_t         Electron_sc_energy[200];   //[Electron_n]
        Float_t         Electron_sc_eta[200];   //[Electron_n]
	Float_t         Electron_sc_phi[200];   //[Electron_n]

        // Muon Collection
        Int_t           Muon_n;
        Float_t         Muon_px[200];   //[Muon_n]
        Float_t         Muon_py[200];   //[Muon_n]
        Float_t         Muon_pz[200];   //[Muon_n]
        Float_t         Muon_vx[200];   //[Muon_n]
	Float_t         Muon_vy[200];   //[Muon_n]
	Float_t         Muon_vz[200];   //[Muon_n]
        Float_t         Muon_pt[200];   //[Muon_n]
        Float_t         Muon_eta[200];   //[Muon_n]
        Float_t         Muon_phi[200];   //[Muon_n]
	Float_t         Muon_energy[200];   //[Muon_n]
        Float_t         Muon_charge[200];   //[Muon_n]
        Bool_t          Muon_isGlobalMuon[200];   //[Muon_n]
        Bool_t          Muon_isTrackerMuon[200];   //[Muon_n]
	Bool_t          Muon_isStandAloneMuon[200];   //[Muon_n]
        Bool_t          Muon_InnerTrack_isNonnull[200];   //[Muon_n]
        Bool_t          Muon_OuterTrack_isNonnull[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_y[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_InnerPoint_z[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_InnerPoint_px[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_InnerPoint_py[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_pz[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_x[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_y[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_z[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_px[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_py[200];   //[Muon_n]
        Float_t         Muon_OuterTrack_OuterPoint_pz[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_y[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_InnerPoint_z[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_InnerPoint_px[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_InnerPoint_py[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_InnerPoint_pz[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_OuterPoint_x[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_OuterPoint_y[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_z[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_px[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_OuterPoint_py[200];   //[Muon_n]
        Float_t         Muon_InnerTrack_OuterPoint_pz[200];   //[Muon_n]
        Float_t         Muon_trackIso[200];   //[Muon_n]
	Float_t         Muon_ecalIso[200];   //[Muon_n]
        Float_t         Muon_hcalIso[200];   //[Muon_n]
        Float_t         Muon_relIso[200];   //[Muon_n]
        Int_t           Muon_normChi2[200];   //[Muon_n]
	Int_t           Muon_validHits[200];   //[Muon_n]
        Int_t           Muon_tkHits[200];   //[Muon_n]
        Int_t           Muon_pixHits[200];   //[Muon_n]
        Int_t           Muon_numberOfMatches[200];   //[Muon_n]
	Float_t         Muon_OuterPoint_x[200];   //[Muon_n]
        Float_t         Muon_OuterPoint_y[200];   //[Muon_n]
        Float_t         Muon_OuterPoint_z[200];   //[Muon_n]
        Float_t         Muon_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_InnerPoint_y[200];   //[Muon_n]
        Float_t         Muon_InnerPoint_z[200];   //[Muon_n]

        //Cosmic Muons
        Int_t           CosmicMuon_n;
        Float_t         CosmicMuon_px[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_py[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_pz[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_pt[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_eta[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_phi[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_energy[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_charge[200];   //[CosmicMuon_n]
        Bool_t          CosmicMuon_isGlobalMuon[200];   //[CosmicMuon_n]
        Bool_t          CosmicMuon_isTrackerMuon[200];   //[CosmicMuon_n]
        Bool_t          CosmicMuon_isStandAloneMuon[200];   //[CosmicMuon_n]
        Bool_t          CosmicMuon_InnerTrack_isNonnull[200];   //[CosmicMuon_n]
        Bool_t          CosmicMuon_OuterTrack_isNonnull[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_x[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_y[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_z[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_px[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_py[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_x[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_y[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_z[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_px[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_py[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterTrack_OuterPoint_pz[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_InnerPoint_x[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_InnerPoint_y[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_InnerPoint_z[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_InnerPoint_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_py[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_x[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_y[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_z[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_px[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_py[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_InnerTrack_OuterPoint_pz[200];   //[CosmicMuon_n]

        //for  AOD only
        Float_t         CosmicMuon_OuterPoint_x[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterPoint_y[200];   //[CosmicMuon_n]
        Float_t         CosmicMuon_OuterPoint_z[200];   //[CosmicMuon_n]

        // Photon Information
        Int_t           Photon_n;
        Float_t         Photon_E[200];   //[Photon_n]
        Float_t         Photon_pt[200];   //[Photon_n]
        Float_t         Photon_eta[200];   //[Photon_n]
        Float_t         Photon_phi[200];   //[Photon_n]
        Float_t         Photon_theta[200];   //[Photon_n]
        Float_t         Photon_et[200];   //[Photon_n]
        Float_t         Photon_swissCross[200];   //[Photon_n]
        Float_t         Photon_e6e2[200];   //[Photon_n]
        Float_t         Photon_e4e1[200];   //[Photon_n]
        Float_t         Photonr9[200];   //[Photon_n]
        Float_t         Photon_e1x5[200];   //[Photon_n]
        Float_t         Photon_e2x5[200];   //[Photon_n]
        Float_t         Photon_e3x3[200];   //[Photon_n]
        Float_t         Photon_e5x5[200];   //[Photon_n]
        Float_t         Photon_r1x5[200];   //[Photon_n]
        Float_t         Photon_r2x5[200];   //[Photon_n]
	Float_t         Photon_maxEnergyXtal[200];   //[Photon_n]
        Float_t         Photon_SigmaEtaEta[200];   //[Photon_n]
        Float_t         Photon_SigmaIetaIeta[200];   //[Photon_n]
        Float_t         Photon_SigmaEtaPhi[200];   //[Photon_n]
        Float_t         Photon_SigmaIetaIphi[200];   //[Photon_n]
        Float_t         Photon_SigmaPhiPhi[200];   //[Photon_n]
        Float_t         Photon_SigmaIphiIphi[200];   //[Photon_n]
        Float_t         Photon_Roundness[200];   //[Photon_n]
	Float_t         Photon_Angle[200];   //[Photon_n]

        Float_t         Photon_ecalRecHitSumEtConeDR03[200];   //[Photon_n]
        Float_t         Photon_hcalTowerSumEtConeDR03[200];   //[Photon_n]
        Float_t         Photon_trkSumPtSolidConeDR03[200];   //[Photon_n]
	Float_t         Photon_trkSumPtHollowConeDR03[200];   //[Photon_n]
        Int_t           Photon_nTrkSolidConeDR03[200];   //[Photon_n]
        Int_t           Photon_nTrkHollowConeDR03[200];   //[Photon_n]
        Float_t         Photon_hcalDepth1TowerSumEtConeDR03[200];   //[Photon_n]
        Float_t         Photon_hcalDepth2TowerSumEtConeDR03[200];   //[Photon_n]
        Float_t         Photon_ecalRecHitSumEtConeDR04[200];   //[Photon_n]
        Float_t         Photon_hcalTowerSumEtConeDR04[200];   //[Photon_n]
        Float_t         Photon_trkSumPtSolidConeDR04[200];   //[Photon_n]
        Float_t         Photon_trkSumPtHollowConeDR04[200];   //[Photon_n]
        Int_t           Photon_nTrkSolidConeDR04[200];   //[Photon_n]
        Int_t           Photon_nTrkHollowConeDR04[200];   //[Photon_n]
        Float_t         Photon_hcalDepth1TowerSumEtConeDR04[200];   //[Photon_n]
        Float_t         Photon_hcalDepth2TowerSumEtConeDR04[200];   //[Photon_n]

        Bool_t          Photon_hasPixelSeed[200];   //[Photon_n]
        Bool_t          Photon_isEB[200];   //[Photon_n]
        Bool_t          Photon_isEE[200];   //[Photon_n]
        Bool_t          Photon_isEBGap[200];   //[Photon_n]
        Bool_t          Photon_isEEGap[200];   //[Photon_n]
	Bool_t          Photon_isEBEEGap[200];   //[Photon_n]
        Float_t         Photon_e2e9[200];   //[Photon_n]
        Float_t         Photon_HoE[200];   //[Photon_n]
        Float_t         Photon_HoEnew[200];   //[Photon_n]
        Float_t         Photon_px[200];   //[Photon_n]
        Float_t         Photon_py[200];   //[Photon_n]
        Float_t         Photon_pz[200];   //[Photon_n]
        Float_t         Photon_vx[200];   //[Photon_n]
	Float_t         Photon_vy[200];   //[Photon_n]
	Float_t         Photon_vz[200];   //[Photon_n]
        Int_t           Photon_no_of_basic_clusters[200];   //[Photon_n]
        Float_t         Photon_sc_energy[200];   //[Photon_n]
        Float_t         Photon_sc_eta[200];   //[Photon_n]
	Float_t         Photon_sc_phi[200];   //[Photon_n]
        Float_t         Photon_sc_x[200];   //[Photon_n]
        Float_t         Photon_sc_y[200];   //[Photon_n]
        Float_t         Photon_sc_z[200];   //[Photon_n]
        Float_t         Photon_etaWidth[200];   //[Photon_n]
        Float_t         Photon_phiWidth[200];   //[Photon_n]
        Float_t         Photon_sc_et[200];   //[Photon_n]

        // Gen matched photon variables
        Float_t         matchphotonE[200];   //[Photon_n]
        Float_t         matchphotonpt[200];   //[Photon_n]
        Float_t         matchphotoneta[200];   //[Photon_n]
        Float_t         matchphotonphi[200];   //[Photon_n]
	Float_t         matchphotonpx[200];   //[Photon_n]
        Float_t         matchphotonpy[200];   //[Photon_n]
        Float_t         matchphotonpz[200];   //[Photon_n]
        Bool_t          ismatchedphoton[200];   //[Photon_n]

        //Converted photon variables
	Bool_t          Photon_hasConvTrk[200];   //[Photon_n]
        Int_t           Photon_ntracks[200];   //[Photon_n]
        Bool_t          Photon_isconverted[200];   //[Photon_n]
        Float_t         Photon_pairInvmass[200];   //[Photon_n]
        Float_t         Photon_pairCotThetaSeperation[200];   //[Photon_n]
        Float_t         Photon_pairmomentumX[200];   //[Photon_n]
        Float_t         Photon_pairmomentumY[200];   //[Photon_n]
        Float_t         Photon_pairmomentumZ[200];   //[Photon_n]
	Float_t         Photon_EoverP[200];   //[Photon_n]
	Float_t         Photon_ConvVx[200];   //[Photon_n]
        Float_t         Photon_ConvVy[200];   //[Photon_n]
        Float_t         Photon_ConvVz[200];   //[Photon_n]
        Float_t         Photon_ZOfPrimaryVertex[200];   //[Photon_n]
	Float_t         Photon_distOfMinimumApproach[200];   //[Photon_n]
        Float_t         Photon_dPhiTracksAtVtx[200];   //[Photon_n]
        Float_t         Photon_dPhiTracksAtEcal[200];   //[Photon_n]
        Float_t         Photon_dEtaTracksAtEcal[200];   //[Photon_n]

        Int_t           npho;
        Bool_t          Photon_Electronveto[200];   //[npho]
	Float_t         PFiso_Charged03[200];   //[npho]
        Float_t         PFiso_Photon03[200];   //[npho]
        Float_t         PFiso_Neutral03[200];   //[npho]
        Float_t         PFiso_Sum03[200];   //[npho]
        Float_t         PFWorstiso_Charged03[200];   //[npho]

	Int_t           Photon_ncrys[200];   //[Photon_n]
        Float_t         Photon_timing_xtal[200][100];   //[Photon_n]
        Float_t         Photon_timingavg_xtal[200];   //[Photon_n]
        Float_t         Photon_energy_xtal[200][100];   //[Photon_n]
	Int_t           Photon_ieta_xtalEB[200][100];   //[Photon_n]
        Int_t           Photon_iphi_xtalEB[200][100];   //[Photon_n]
	Int_t           Photon_recoFlag_xtalEB[200][100];   //[Photon_n]
        Float_t         Photon_timeError_xtal[200][100];   //[Photon_n]
        Float_t         Photon_s9[200];   //[Photon_n]

        //MIP variable
        Float_t         Photon_mipChi2[200];   //[Photon_n]
        Float_t         Photon_mipTotEnergy[200];   //[Photon_n]
        Float_t         Photon_mipSlope[200];   //[Photon_n]
        Float_t         Photon_mipIntercept[200];   //[Photon_n]
        Int_t           Photon_mipNhitCone[200];   //[Photon_n]
        Bool_t          Photon_mipIsHalo[200];   //[Photon_n]

        // EBrechit variables
        Int_t           EBRecHit_size;
        Float_t         EBRecHit_eta[10000];   //[EBRecHit_size]
        Float_t         EBRecHit_phi[10000];   //[EBRecHit_size]
        Int_t           EBRecHit_ieta[10000];   //[EBRecHit_size]
        Int_t           EBRecHit_iphi[10000];   //[EBRecHit_size]
        Float_t         EBRecHit_e[10000];   //[EBRecHit_size]
        Float_t         EBRecHit_et[10000];   //[EBRecHit_size]
        Int_t           EBRecHit_flag[10000];   //[EBRecHit_size]

        // EErechit variables
        Float_t         EBRecHit_time[10000];   //[EBRecHit_size]
        Int_t           EERecHit_size;
        Float_t         EERecHit_eta[10000];   //[EERecHit_size]
        Float_t         EERecHit_phi[10000];   //[EERecHit_size]
        Int_t           EERecHit_ieta[10000];   //[EERecHit_size]
        Int_t           EERecHit_iphi[10000];   //[EERecHit_size]
        Float_t         EERecHit_e[10000];   //[EERecHit_size]
        Float_t         EERecHit_et[10000];   //[EERecHit_size]
        Int_t           EERecHit_flag[10000];   //[EERecHit_size]
        Float_t         EERecHit_time[10000];   //[EERecHit_size]

        // PFMET variables
        Float_t         PFMetPt[6];
        Float_t         PFMetPx[6];
        Float_t         PFMetPy[6];
	Float_t         PFMetPhi[6];
        Float_t         PFMetSumEt[6];
        Float_t         Delta_phiPF;

        //CaloTowers
        Int_t           CaloTower_n;
	Float_t         CaloTower_eta[5000];   //[CaloTower_n]
        Float_t         CaloTower_phi[5000];   //[CaloTower_n]
        Float_t         CaloTower_E[5000];   //[CaloTower_n]
        Float_t         CaloTower_Et[5000];   //[CaloTower_n]
	Float_t         CaloTower_emEnergy[5000];   //[CaloTower_n]
        Float_t         CaloTower_hadEnergy[5000];   //[CaloTower_n]
        Float_t         CaloTower_p[5000];   //[CaloTower_n]
        Float_t         CaloTower_EMEt[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadEt[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadPhi[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadEta[5000];   //[CaloTower_n]
        Float_t         CaloTower_EMPhi[5000];   //[CaloTower_n]
        Float_t         CaloTower_EMEta[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadX[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadY[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadZ[5000];   //[CaloTower_n]
	Float_t         CaloTower_HE_E[5000];   //[CaloTower_n]
        Float_t         CaloTower_HB_E[5000];   //[CaloTower_n]
        Float_t         CaloTower_EMTime[5000];   //[CaloTower_n]
        Float_t         CaloTower_HadTime[5000];   //[CaloTower_n]

        Float_t         rho;
        Float_t         sigma;

        Float_t         rho25;
        Float_t         sigma25;

        // List of branches
        TBranch        *b_nevents;   //!
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
        TBranch        *b_ntriggers;   //!
        TBranch        *b_triggernames;   //!
        TBranch        *b_triggerprescales;   //!
        TBranch        *b_ifTriggerpassed;   //!
        TBranch        *b_ObjectPt;   //!
        TBranch        *b_ObjectEta;   //!
        TBranch        *b_ObjectPhi;   //!
        TBranch        *b_FilterNames;   //!
        TBranch        *b_FilterStartPosition;   //!
	TBranch        *b_FilterEndPosition;   //!
        TBranch        *b_ObjectStartPosition;   //!
        TBranch        *b_ObjectEndPosition;   //!
        TBranch        *b_Vertex_n;   //!
	TBranch        *b_Vertex_x;   //!
        TBranch        *b_Vertex_y;   //!
        TBranch        *b_Vertex_z;   //!
	TBranch        *b_Vertex_tracksize;   //!
        TBranch        *b_Vertex_ndof;   //!
        TBranch        *b_Vertex_chi2;   //!
	TBranch        *b_Vertex_d0;   //!
        TBranch        *b_Vertex_isFake;   //!
        TBranch        *b_Scraping_isScrapingEvent;   //!
	TBranch        *b_Scraping_numOfTracks;   //!
        TBranch        *b_Scraping_fractionOfGoodTracks;   //!
        TBranch        *b_Track_n;   //!
	TBranch        *b_Track_px;   //!
        TBranch        *b_Track_py;   //!
        TBranch        *b_Track_pz;   //!
        TBranch        *b_Track_vx;   //!
	TBranch        *b_Track_vy;   //!
        TBranch        *b_Track_vz;   //!
        TBranch        *b_Track_pt;   //!
        TBranch        *b_Track_eta;   //!
        TBranch        *b_Track_phi;   //!
        TBranch        *b_pfJet_n;   //!
        TBranch        *b_pfJet_px;   //!
        TBranch        *b_pfJet_py;   //!
        TBranch        *b_pfJet_E;   //!
        TBranch        *b_pfJet_pz;   //!
        TBranch        *b_pfJet_vx;   //!
        TBranch        *b_pfJet_vy;   //!
        TBranch        *b_pfJet_vz;   //!
        TBranch        *b_pfJet_pt;   //!
        TBranch        *b_pfJet_eta;   //!
        TBranch        *b_pfJet_phi;   //!
        TBranch        *b_pfjet_CEF;   //!
        TBranch        *b_pfjet_CHF;   //!
        TBranch        *b_pfjet_NEF;   //!
        TBranch        *b_pfjet_NHF;   //!
        TBranch        *b_pfjet_NCH;   //!
        TBranch        *b_pfjet_HFHAE;   //!
        TBranch        *b_pfjet_HFEME;   //!
        TBranch        *b_pfjet_NConstituents;   //!
        TBranch        *b_pfJet_partonFlavor;   //!
        TBranch        *b_pfJet_partonStatus;   //!
        TBranch        *b_pujetIdFull_mva;   //!
        TBranch        *b_pujetIdSimple_mva;   //!
        TBranch        *b_pujetIdCutBased_mva;   //!
        TBranch        *b_pujetIdFull_loose;   //!
        TBranch        *b_pujetIdFull_medium;   //!
        TBranch        *b_pujetIdFull_tight;   //!
        TBranch        *b_pujetIdSimple_loose;   //!
        TBranch        *b_pujetIdSimple_medium;   //!
        TBranch        *b_pujetIdSimple_tight;   //!
        TBranch        *b_pujetIdCutBased_loose;   //!
        TBranch        *b_pujetIdCutBased_medium;   //!
        TBranch        *b_pujetIdCutBased_tight;   //!
        TBranch        *b_pfjet_TrackCountHiEffBJetTags;   //!
        TBranch        *b_pfjet_TrackCountHiPurBJetTags;   //!
        TBranch        *b_pfjet_SimpleSVHiEffBJetTags;   //!
        TBranch        *b_pfjet_SimpleSVHiPurBJetTags;   //!
        TBranch        *b_pfJet_jecUncer;   //!
        TBranch        *b_pfJet_jecCorr;   //!
        TBranch        *b_ucpfJet_px;   //!
        TBranch        *b_ucpfJet_py;   //!
        TBranch        *b_ucpfJet_E;   //!
        TBranch        *b_ucpfJet_pz;   //!
        TBranch        *b_ucpfJet_pt;   //!
        TBranch        *b_ucpfJet_eta;   //!
        TBranch        *b_ucpfJet_phi;   //!
        TBranch        *b_Electron_n;   //!
        TBranch        *b_Electron_px;   //!
        TBranch        *b_Electron_py;   //!
        TBranch        *b_Electron_pz;   //!
        TBranch        *b_Electron_vx;   //!
        TBranch        *b_Electron_vy;   //!
        TBranch        *b_Electron_vz;   //!
        TBranch        *b_Electron_pt;   //!
        TBranch        *b_Electron_eta;   //!
        TBranch        *b_Electron_phi;   //!
        TBranch        *b_Electron_energy;   //!
        TBranch        *b_Electron_charge;   //!
        TBranch        *b_Electron_trkIso;   //!
        TBranch        *b_Electron_ecalIso;   //!
        TBranch        *b_Electron_hcalIso;   //!
        TBranch        *b_Electron_SigmaIetaIeta;   //!
        TBranch        *b_Electron_dEtaIn;   //!
        TBranch        *b_Electron_dPhiIn;   //!
        TBranch        *b_Electron_HoE;   //!
        TBranch        *b_Electron_sc_energy;   //!
        TBranch        *b_Electron_sc_eta;   //!
        TBranch        *b_Electron_sc_phi;   //!
        TBranch        *b_Muon_n;   //!
        TBranch        *b_Muon_px;   //!
        TBranch        *b_Muon_py;   //!
        TBranch        *b_Muon_pz;   //!
        TBranch        *b_Muon_vx;   //!
        TBranch        *b_Muon_vy;   //!
        TBranch        *b_Muon_vz;   //!
        TBranch        *b_Muon_pt;   //!
        TBranch        *b_Muon_eta;   //!
        TBranch        *b_Muon_phi;   //!
        TBranch        *b_Muon_energy;   //!
        TBranch        *b_Muon_charge;   //!
        TBranch        *b_Muon_isGlobalMuon;   //!
        TBranch        *b_Muon_isTrackerMuon;   //!
        TBranch        *b_Muon_isStandAloneMuon;   //!
        TBranch        *b_Muon_InnerTrack_isNonnull;   //!
        TBranch        *b_Muon_OuterTrack_isNonnull;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_x;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_y;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_z;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_px;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_py;   //!
        TBranch        *b_Muon_OuterTrack_InnerPoint_pz;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_x;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_y;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_z;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_px;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_py;   //!
        TBranch        *b_Muon_OuterTrack_OuterPoint_pz;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_x;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_y;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_z;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_px;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_py;   //!
        TBranch        *b_Muon_InnerTrack_InnerPoint_pz;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_x;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_y;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_z;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_px;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_py;   //!
        TBranch        *b_Muon_InnerTrack_OuterPoint_pz;   //!
        TBranch        *b_Muon_trackIso;   //!
        TBranch        *b_Muon_ecalIso;   //!
        TBranch        *b_Muon_hcalIso;   //!
        TBranch        *b_Muon_relIso;   //!
        TBranch        *b_Muon_normChi2;   //!
        TBranch        *b_Muon_validHits;   //!
        TBranch        *b_Muon_tkHits;   //!
        TBranch        *b_Muon_pixHits;   //!
        TBranch        *b_Muon_numberOfMatches;   //!
        TBranch        *b_Muon_OuterPoint_x;   //!
        TBranch        *b_Muon_OuterPoint_y;   //!
        TBranch        *b_Muon_OuterPoint_z;   //!
        TBranch        *b_Muon_InnerPoint_x;   //!
        TBranch        *b_Muon_InnerPoint_y;   //!
        TBranch        *b_Muon_InnerPoint_z;   //!
        TBranch        *b_CosmicMuon_n;   //!
        TBranch        *b_CosmicMuon_px;   //!
        TBranch        *b_CosmicMuon_py;   //!
        TBranch        *b_CosmicMuon_pz;   //!
        TBranch        *b_CosmicMuon_pt;   //!
        TBranch        *b_CosmicMuon_eta;   //!
        TBranch        *b_CosmicMuon_phi;   //!
        TBranch        *b_CosmicMuon_energy;   //!
        TBranch        *b_CosmicMuon_charge;   //!
        TBranch        *b_CosmicMuon_isGlobalMuon;   //!
        TBranch        *b_CosmicMuon_isTrackerMuon;   //!
        TBranch        *b_CosmicMuon_isStandAloneMuon;   //!
        TBranch        *b_CosmicMuon_InnerTrack_isNonnull;   //!
        TBranch        *b_CosmicMuon_OuterTrack_isNonnull;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_x;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_y;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_z;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_px;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_py;   //!
        TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_pz;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_x;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_y;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_z;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_px;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_py;   //!
        TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_pz;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_x;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_y;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_z;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_px;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_py;   //!
        TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_pz;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_x;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_y;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_z;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_px;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_py;   //!
        TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_pz;   //!
        TBranch        *b_CosmicMuon_OuterPoint_x;   //!
        TBranch        *b_CosmicMuon_OuterPoint_y;   //!
        TBranch        *b_CosmicMuon_OuterPoint_z;   //!
        TBranch        *b_Photon_n;   //!
        TBranch        *b_Photon_E;   //!
        TBranch        *b_Photon_pt;   //!
        TBranch        *b_Photon_eta;   //!
        TBranch        *b_Photon_phi;   //!
        TBranch        *b_Photon_theta;   //!
        TBranch        *b_Photon_et;   //!
        TBranch        *b_Photon_swissCross;   //!
        TBranch        *b_Photon_e6e2;   //!
        TBranch        *b_Photon_e4e1;   //!
        TBranch        *b_Photonr9;   //!
        TBranch        *b_Photon_e1x5;   //!
        TBranch        *b_Photon_e2x5;   //!
        TBranch        *b_Photon_e3x3;   //!
        TBranch        *b_Photon_e5x5;   //!
        TBranch        *b_Photon_r1x5;   //!
        TBranch        *b_Photon_r2x5;   //!
        TBranch        *b_Photon_maxEnergyXtal;   //!
        TBranch        *b_Photon_SigmaEtaEta;   //!
        TBranch        *b_Photon_SigmaIetaIeta;   //!
        TBranch        *b_Photon_SigmaEtaPhi;   //!
        TBranch        *b_Photon_SigmaIetaIphi;   //!
        TBranch        *b_Photon_SigmaPhiPhi;   //!
        TBranch        *b_Photon_SigmaIphiIphi;   //!
        TBranch        *b_Photon_Roundness;   //!
        TBranch        *b_Photon_Angle;   //!
        TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
        TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
        TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
        TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
        TBranch        *b_Photon_nTrkSolidConeDR03;   //!
        TBranch        *b_Photon_nTrkHollowConeDR03;   //!
        TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
        TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
        TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
        TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
        TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
        TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
        TBranch        *b_Photon_nTrkSolidConeDR04;   //!
        TBranch        *b_Photon_nTrkHollowConeDR04;   //!
        TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
        TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
        TBranch        *b_Photon_hasPixelSeed;   //!
        TBranch        *b_Photon_isEB;   //!
        TBranch        *b_Photon_isEE;   //!
        TBranch        *b_Photon_isEBGap;   //!
        TBranch        *b_Photon_isEEGap;   //!
        TBranch        *b_Photon_isEBEEGap;   //!
        TBranch        *b_Photon_e2e9;   //!
        TBranch        *b_Photon_HoE;   //!
        TBranch        *b_Photon_HoEnew;   //!
        TBranch        *b_Photon_px;   //!
        TBranch        *b_Photon_py;   //!
        TBranch        *b_Photon_pz;   //!
        TBranch        *b_Photon_vx;   //!
        TBranch        *b_Photon_vy;   //!
        TBranch        *b_Photon_vz;   //!
        TBranch        *b_Photon_no_of_basic_clusters;   //!
        TBranch        *b_Photon_sc_energy;   //!
        TBranch        *b_Photon_sc_eta;   //!
        TBranch        *b_Photon_sc_phi;   //!
        TBranch        *b_Photon_sc_x;   //!
        TBranch        *b_Photon_sc_y;   //!
        TBranch        *b_Photon_sc_z;   //!
        TBranch        *b_Photon_etaWidth;   //!
        TBranch        *b_Photon_phiWidth;   //!
        TBranch        *b_Photon_sc_et;   //!
        TBranch        *b_matchphotonE;   //!
        TBranch        *b_matchphotonpt;   //!
        TBranch        *b_matchphotoneta;   //!
        TBranch        *b_matchphotonphi;   //!
        TBranch        *b_matchphotonpx;   //!
        TBranch        *b_matchphotonpy;   //!
        TBranch        *b_matchphotonpz;   //!
        TBranch        *b_ismatchedphoton;   //!
        TBranch        *b_Photon_hasConvTrk;   //!
        TBranch        *b_Photon_ntracks;   //!
        TBranch        *b_Photon_isconverted;   //!
        TBranch        *b_Photon_pairInvmass;   //!
        TBranch        *b_Photon_pairCotThetaSeperation;   //!
        TBranch        *b_Photon_pairmomentumX;   //!
        TBranch        *b_Photon_pairmomentumY;   //!
        TBranch        *b_Photon_pairmomentumZ;   //!
        TBranch        *b_Photon_EoverP;   //!
        TBranch        *b_Photon_ConvVx;   //!
        TBranch        *b_Photon_ConvVy;   //!
        TBranch        *b_Photon_ConvVz;   //!
        TBranch        *b_Photon_ZOfPrimaryVertex;   //!
        TBranch        *b_Photon_distOfMinimumApproach;   //!
        TBranch        *b_Photon_dPhiTracksAtVtx;   //!
        TBranch        *b_Photon_dPhiTracksAtEcal;   //!
        TBranch        *b_Photon_dEtaTracksAtEcal;   //!
        TBranch        *b_npho;   //!
        TBranch        *b_Photon_Electronveto;   //!
        TBranch        *b_PFiso_Charged03;   //!
        TBranch        *b_PFiso_Photon03;   //!
        TBranch        *b_PFiso_Neutral03;   //!
        TBranch        *b_PFiso_Sum03;   //!
        TBranch        *b_PFWorstiso_Charged03;   //!
        TBranch        *b_Photon_ncrys;   //!
        TBranch        *b_Photon_timing_xtal;   //!
        TBranch        *b_Photon_timingavg_xtal;   //!
        TBranch        *b_Photon_energy_xtal;   //!
        TBranch        *b_Photon_ieta_xtalEB;   //!
        TBranch        *b_Photon_iphi_xtalEB;   //!
        TBranch        *b_Photon_recoFlag_xtalEB;   //!
        TBranch        *b_Photon_timeError_xtal;   //!
        TBranch        *b_Photon_s9;   //!
        TBranch        *b_Photon_mipChi2;   //!
        TBranch        *b_Photon_mipTotEnergy;   //!
        TBranch        *b_Photon_mipSlope;   //!
        TBranch        *b_Photon_mipIntercept;   //!
        TBranch        *b_Photon_mipNhitCone;   //!
        TBranch        *b_Photon_mipIsHalo;   //!
        TBranch        *b_EBRecHit_size;   //!
        TBranch        *b_EBRecHit_eta;   //!
        TBranch        *b_EBRecHit_phi;   //!
        TBranch        *b_EBRecHit_ieta;   //!
        TBranch        *b_EBRecHit_iphi;   //!
        TBranch        *b_EBRecHit_e;   //!
        TBranch        *b_EBRecHit_et;   //!
        TBranch        *b_EBRecHit_flag;   //!
        TBranch        *b_EBRecHit_time;   //!
        TBranch        *b_EERecHit_size;   //!
        TBranch        *b_EERecHit_eta;   //!
        TBranch        *b_EERecHit_phi;   //!
        TBranch        *b_EERecHit_ieta;   //!
        TBranch        *b_EERecHit_iphi;   //!
        TBranch        *b_EERecHit_e;   //!
        TBranch        *b_EERecHit_et;   //!
        TBranch        *b_EERecHit_flag;   //!
        TBranch        *b_EERecHit_time;   //!
        TBranch        *b_PFMetPt;   //!
        TBranch        *b_PFMetPx;   //!
        TBranch        *b_PFMetPy;   //!
        TBranch        *b_PFMetPhi;   //!
        TBranch        *b_PFMetSumEt;   //!
        TBranch        *b_Delta_phiPF;   //!
        TBranch        *b_CaloTower_n;   //!
        TBranch        *b_CaloTower_eta;   //!
        TBranch        *b_CaloTower_phi;   //!
        TBranch        *b_CaloTower_E;   //!
        TBranch        *b_CaloTower_Et;   //!
        TBranch        *b_CaloTower_emEnergy;   //!
        TBranch        *b_CaloTower_hadEnergy;   //!
        TBranch        *b_CaloTower_p;   //!
        TBranch        *b_CaloTower_EMEt;   //!
        TBranch        *b_CaloTower_HadEt;   //!
        TBranch        *b_CaloTower_HadPhi;   //!
        TBranch        *b_CaloTower_HadEta;   //!
        TBranch        *b_CaloTower_EMPhi;   //!
        TBranch        *b_CaloTower_EMEta;   //!
        TBranch        *b_CaloTower_HadX;   //!
        TBranch        *b_CaloTower_HadY;   //!
        TBranch        *b_CaloTower_HadZ;   //!
        TBranch        *b_CaloTower_HE_E;   //!
        TBranch        *b_CaloTower_HB_E;   //!
        TBranch        *b_CaloTower_EMTime;   //!
        TBranch        *b_CaloTower_HadTime;   //!
        TBranch        *b_rho;   //!
        TBranch        *b_sigma;   //!
        TBranch        *b_rho25;   //!
        TBranch        *b_sigma25;   //!

        PostAnalyzerData();
        virtual ~PostAnalyzerData();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
        virtual Bool_t   PassHLT(Bool_t &passedHLT150, Bool_t &passed90HLTall);
        virtual Bool_t   NonScraping();
        virtual Bool_t   NoSpike(Int_t ipho);
        virtual Bool_t   PrimaryVertex(Int_t &goodVertex);
        virtual Double_t getRapidity(Double_t r_E, Double_t r_Pz);
        virtual Double_t getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
        virtual Double_t getDEta(Double_t eta1, Double_t eta2);
        virtual Double_t getDPhi(Double_t phit1, Double_t phi2);
        virtual Double_t getMass(Int_t pho_i, Int_t jet_i);
        virtual Double_t getLICTD(Int_t i);
        virtual Bool_t   TightPhotonPFIso(Int_t ipho);
        virtual Bool_t   MediumPhotonPFIso(Int_t ipho);
        virtual Bool_t   LoosePhotonPFIso(Int_t ipho);
        virtual Bool_t   TightJetID( Int_t ijet);
        virtual Double_t EAElectroncharged(Double_t eta);
        virtual Double_t EAElectronneutral(Double_t eta);
        virtual Double_t EAElectronphoton(Double_t eta);
        virtual Double_t EXOtightPhoID(Int_t ipho);
        virtual void     BookHistos();
};

#endif

#ifdef PostAnalyzerData_cxx
PostAnalyzerData::PostAnalyzerData(){

    TChain *chain = new TChain("myEvent");
    //add input files here
      
      chain->Add("/eos/uscms/store/user/varun/Datasets8TeV/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_100_1_YzW.root/myEvent");
      chain->Add("/eos/uscms/store/user/varun/Datasets8TeV/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_101_1_ZGH.root/myEvent");
      chain->Add("/eos/uscms/store/user/varun/Datasets8TeV/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_102_1_QTy.root/myEvent");
      chain->Add("/eos/uscms/store/user/varun/Datasets8TeV/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_103_1_nGN.root/myEvent");
      chain->Add("/eos/uscms/store/user/varun/Datasets8TeV/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_104_1_11Q.root/myEvent");
      /*

    chain->Add("root://xrootd.unl.edu//store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_100_1_YzW.root/myEvent");
    chain->Add("root://xrootd.unl.edu//store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_101_1_ZGH.root/myEvent");
    chain->Add("root://xrootd.unl.edu//store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_102_1_QTy.root/myEvent");
    chain->Add("root://xrootd.unl.edu//store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_103_1_nGN.root/myEvent");
    chain->Add("root://xrootd.unl.edu//store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run2012A_22Jan_104_1_11Q.root/myEvent");

      */
    /*     chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run201 \
2A_22Jan_100_1_YzW.root/myEvent");
     chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run201\
2A_22Jan_101_1_ZGH.root/myEvent");
     chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run201\
2A_22Jan_102_1_QTy.root/myEvent");
     chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run201\
2A_22Jan_103_1_nGN.root/myEvent");
     chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/Run2012A_22Jan/AOD_Output_Run201\
2A_22Jan_104_1_11Q.root/myEvent");
    */
Init(chain);
}

PostAnalyzerData::~PostAnalyzerData(){

    if (!fChain) return;
    delete fChain->GetCurrentFile();
    f1->cd();
    f1->Write();
    f1->Close();

}

Int_t PostAnalyzerData::GetEntry(Long64_t entry){
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t PostAnalyzerData::LoadTree(Long64_t entry){
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

void PostAnalyzerData::Init(TTree *tree) {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    triggernames = 0;
    triggerprescales = 0;
    ifTriggerpassed = 0;
    ObjectPt = 0;
    ObjectEta = 0;
    ObjectPhi = 0;
    FilterNames = 0;
    FilterStartPosition = 0;
    FilterEndPosition = 0;
    ObjectStartPosition = 0;
    ObjectEndPosition = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("nevents", &nevents, &b_nevents);
    fChain->SetBranchAddress("run", &run, &b_RunNumber);
    fChain->SetBranchAddress("event", &event, &b_EventNumber);
    fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_LumiNumber);
    fChain->SetBranchAddress("beamCrossing", &beamCrossing, &b_BXNumber);
    fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
    fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
    fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
    fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
    fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
    fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
    fChain->SetBranchAddress("ntriggers", &ntriggers, &b_ntriggers);
    fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
    fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
    fChain->SetBranchAddress("ifTriggerpassed", &ifTriggerpassed, &b_ifTriggerpassed);
    fChain->SetBranchAddress("ObjectPt", &ObjectPt, &b_ObjectPt);
    fChain->SetBranchAddress("ObjectEta", &ObjectEta, &b_ObjectEta);
    fChain->SetBranchAddress("ObjectPhi", &ObjectPhi, &b_ObjectPhi);
    fChain->SetBranchAddress("FilterNames", &FilterNames, &b_FilterNames);
    fChain->SetBranchAddress("FilterStartPosition", &FilterStartPosition, &b_FilterStartPosition);
    fChain->SetBranchAddress("FilterEndPosition", &FilterEndPosition, &b_FilterEndPosition);
    fChain->SetBranchAddress("ObjectStartPosition", &ObjectStartPosition, &b_ObjectStartPosition);
    fChain->SetBranchAddress("ObjectEndPosition", &ObjectEndPosition, &b_ObjectEndPosition);
    fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
    fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
    fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
    fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
    fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize, &b_Vertex_tracksize);
    fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof, &b_Vertex_ndof);
    fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2, &b_Vertex_chi2);
    fChain->SetBranchAddress("Vertex_d0", Vertex_d0, &b_Vertex_d0);
    fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake, &b_Vertex_isFake);
    fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent, &b_Scraping_isScrapingEvent);
    fChain->SetBranchAddress("Scraping_numOfTracks", &Scraping_numOfTracks, &b_Scraping_numOfTracks);
    fChain->SetBranchAddress("Scraping_fractionOfGoodTracks", &Scraping_fractionOfGoodTracks, &b_Scraping_fractionOfGoodTracks);
    fChain->SetBranchAddress("Track_n", &Track_n, &b_Track_n);
    fChain->SetBranchAddress("Track_px", Track_px, &b_Track_px);
    fChain->SetBranchAddress("Track_py", Track_py, &b_Track_py);
    fChain->SetBranchAddress("Track_pz", Track_pz, &b_Track_pz);
    fChain->SetBranchAddress("Track_vx", Track_vx, &b_Track_vx);
    fChain->SetBranchAddress("Track_vy", Track_vy, &b_Track_vy);
    fChain->SetBranchAddress("Track_vz", Track_vz, &b_Track_vz);
    fChain->SetBranchAddress("Track_pt", Track_pt, &b_Track_pt);
    fChain->SetBranchAddress("Track_eta", Track_eta, &b_Track_eta);
    fChain->SetBranchAddress("Track_phi", Track_phi, &b_Track_phi);
    fChain->SetBranchAddress("pfJet_n", &pfJet_n, &b_pfJet_n);
    fChain->SetBranchAddress("pfJet_px", pfJet_px, &b_pfJet_px);
    fChain->SetBranchAddress("pfJet_py", pfJet_py, &b_pfJet_py);
    fChain->SetBranchAddress("pfJet_E", pfJet_E, &b_pfJet_E);
    fChain->SetBranchAddress("pfJet_pz", pfJet_pz, &b_pfJet_pz);
    fChain->SetBranchAddress("pfJet_vx", pfJet_vx, &b_pfJet_vx);
    fChain->SetBranchAddress("pfJet_vy", pfJet_vy, &b_pfJet_vy);
    fChain->SetBranchAddress("pfJet_vz", pfJet_vz, &b_pfJet_vz);
    fChain->SetBranchAddress("pfJet_pt", pfJet_pt, &b_pfJet_pt);
    fChain->SetBranchAddress("pfJet_eta", pfJet_eta, &b_pfJet_eta);
    fChain->SetBranchAddress("pfJet_phi", pfJet_phi, &b_pfJet_phi);
    fChain->SetBranchAddress("pfjet_CEF", pfjet_CEF, &b_pfjet_CEF);
    fChain->SetBranchAddress("pfjet_CHF", pfjet_CHF, &b_pfjet_CHF);
    fChain->SetBranchAddress("pfjet_NEF", pfjet_NEF, &b_pfjet_NEF);
    fChain->SetBranchAddress("pfjet_NHF", pfjet_NHF, &b_pfjet_NHF);
    fChain->SetBranchAddress("pfjet_NCH", pfjet_NCH, &b_pfjet_NCH);
    fChain->SetBranchAddress("pfjet_HFHAE", pfjet_HFHAE, &b_pfjet_HFHAE);
    fChain->SetBranchAddress("pfjet_HFEME", pfjet_HFEME, &b_pfjet_HFEME);
    fChain->SetBranchAddress("pfjet_NConstituents", pfjet_NConstituents, &b_pfjet_NConstituents);
    fChain->SetBranchAddress("pfJet_partonFlavor", pfJet_partonFlavor, &b_pfJet_partonFlavor);
    fChain->SetBranchAddress("pfJet_partonStatus", pfJet_partonStatus, &b_pfJet_partonStatus);
    fChain->SetBranchAddress("pujetIdFull_mva", pujetIdFull_mva, &b_pujetIdFull_mva);
    fChain->SetBranchAddress("pujetIdSimple_mva", pujetIdSimple_mva, &b_pujetIdSimple_mva);
    fChain->SetBranchAddress("pujetIdCutBased_mva", pujetIdCutBased_mva, &b_pujetIdCutBased_mva);
    fChain->SetBranchAddress("pujetIdFull_loose", pujetIdFull_loose, &b_pujetIdFull_loose);
    fChain->SetBranchAddress("pujetIdFull_medium", pujetIdFull_medium, &b_pujetIdFull_medium);
    fChain->SetBranchAddress("pujetIdFull_tight", pujetIdFull_tight, &b_pujetIdFull_tight);
    fChain->SetBranchAddress("pujetIdSimple_loose", pujetIdSimple_loose, &b_pujetIdSimple_loose);
    fChain->SetBranchAddress("pujetIdSimple_medium", pujetIdSimple_medium, &b_pujetIdSimple_medium);
    fChain->SetBranchAddress("pujetIdSimple_tight", pujetIdSimple_tight, &b_pujetIdSimple_tight);
    fChain->SetBranchAddress("pujetIdCutBased_loose", pujetIdCutBased_loose, &b_pujetIdCutBased_loose);
    fChain->SetBranchAddress("pujetIdCutBased_medium", pujetIdCutBased_medium, &b_pujetIdCutBased_medium);
    fChain->SetBranchAddress("pujetIdCutBased_tight", pujetIdCutBased_tight, &b_pujetIdCutBased_tight);
    fChain->SetBranchAddress("pfjet_TrackCountHiEffBJetTags", pfjet_TrackCountHiEffBJetTags, &b_pfjet_TrackCountHiEffBJetTags);
    fChain->SetBranchAddress("pfjet_TrackCountHiPurBJetTags", pfjet_TrackCountHiPurBJetTags, &b_pfjet_TrackCountHiPurBJetTags);
    fChain->SetBranchAddress("pfjet_SimpleSVHiEffBJetTags", pfjet_SimpleSVHiEffBJetTags, &b_pfjet_SimpleSVHiEffBJetTags);
    fChain->SetBranchAddress("pfjet_SimpleSVHiPurBJetTags", pfjet_SimpleSVHiPurBJetTags, &b_pfjet_SimpleSVHiPurBJetTags);
    fChain->SetBranchAddress("pfJet_jecUncer", pfJet_jecUncer, &b_pfJet_jecUncer);
    fChain->SetBranchAddress("pfJet_jecCorr", pfJet_jecCorr, &b_pfJet_jecCorr);
    fChain->SetBranchAddress("ucpfJet_px", ucpfJet_px, &b_ucpfJet_px);
    fChain->SetBranchAddress("ucpfJet_py", ucpfJet_py, &b_ucpfJet_py);
    fChain->SetBranchAddress("ucpfJet_E", ucpfJet_E, &b_ucpfJet_E);
    fChain->SetBranchAddress("ucpfJet_pz", ucpfJet_pz, &b_ucpfJet_pz);
    fChain->SetBranchAddress("ucpfJet_pt", ucpfJet_pt, &b_ucpfJet_pt);
    fChain->SetBranchAddress("ucpfJet_eta", ucpfJet_eta, &b_ucpfJet_eta);
    fChain->SetBranchAddress("ucpfJet_phi", ucpfJet_phi, &b_ucpfJet_phi);
    fChain->SetBranchAddress("Electron_n", &Electron_n, &b_Electron_n);
    fChain->SetBranchAddress("Electron_px", Electron_px, &b_Electron_px);
    fChain->SetBranchAddress("Electron_py", Electron_py, &b_Electron_py);
    fChain->SetBranchAddress("Electron_pz", Electron_pz, &b_Electron_pz);
    fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
    fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
    fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
    fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
    fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
    fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
    fChain->SetBranchAddress("Electron_energy", Electron_energy, &b_Electron_energy);
    fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
    fChain->SetBranchAddress("Electron_trkIso", Electron_trkIso, &b_Electron_trkIso);
    fChain->SetBranchAddress("Electron_ecalIso", Electron_ecalIso, &b_Electron_ecalIso);
    fChain->SetBranchAddress("Electron_hcalIso", Electron_hcalIso, &b_Electron_hcalIso);
    fChain->SetBranchAddress("Electron_SigmaIetaIeta", Electron_SigmaIetaIeta, &b_Electron_SigmaIetaIeta);
    fChain->SetBranchAddress("Electron_dEtaIn", Electron_dEtaIn, &b_Electron_dEtaIn);
    fChain->SetBranchAddress("Electron_dPhiIn", Electron_dPhiIn, &b_Electron_dPhiIn);
    fChain->SetBranchAddress("Electron_HoE", Electron_HoE, &b_Electron_HoE);
    fChain->SetBranchAddress("Electron_sc_energy", Electron_sc_energy, &b_Electron_sc_energy);
    fChain->SetBranchAddress("Electron_sc_eta", Electron_sc_eta, &b_Electron_sc_eta);
    fChain->SetBranchAddress("Electron_sc_phi", Electron_sc_phi, &b_Electron_sc_phi);
    fChain->SetBranchAddress("Muon_n", &Muon_n, &b_Muon_n);
    fChain->SetBranchAddress("Muon_px", Muon_px, &b_Muon_px);
    fChain->SetBranchAddress("Muon_py", Muon_py, &b_Muon_py);
    fChain->SetBranchAddress("Muon_pz", Muon_pz, &b_Muon_pz);
    fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
    fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
    fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
    fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
    fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
    fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
    fChain->SetBranchAddress("Muon_energy", Muon_energy, &b_Muon_energy);
    fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
    fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
    fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
    fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
    fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull, &b_Muon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull, &b_Muon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_x", Muon_OuterTrack_InnerPoint_x, &b_Muon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_y", Muon_OuterTrack_InnerPoint_y, &b_Muon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_z", Muon_OuterTrack_InnerPoint_z, &b_Muon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_px", Muon_OuterTrack_InnerPoint_px, &b_Muon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_py", Muon_OuterTrack_InnerPoint_py, &b_Muon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_pz", Muon_OuterTrack_InnerPoint_pz, &b_Muon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_x", Muon_OuterTrack_OuterPoint_x, &b_Muon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_y", Muon_OuterTrack_OuterPoint_y, &b_Muon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_z", Muon_OuterTrack_OuterPoint_z, &b_Muon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_px", Muon_OuterTrack_OuterPoint_px, &b_Muon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_py", Muon_OuterTrack_OuterPoint_py, &b_Muon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_pz", Muon_OuterTrack_OuterPoint_pz, &b_Muon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_x", Muon_InnerTrack_InnerPoint_x, &b_Muon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_y", Muon_InnerTrack_InnerPoint_y, &b_Muon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_z", Muon_InnerTrack_InnerPoint_z, &b_Muon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_px", Muon_InnerTrack_InnerPoint_px, &b_Muon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_py", Muon_InnerTrack_InnerPoint_py, &b_Muon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_pz", Muon_InnerTrack_InnerPoint_pz, &b_Muon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_x", Muon_InnerTrack_OuterPoint_x, &b_Muon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_y", Muon_InnerTrack_OuterPoint_y, &b_Muon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_z", Muon_InnerTrack_OuterPoint_z, &b_Muon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_px", Muon_InnerTrack_OuterPoint_px, &b_Muon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_py", Muon_InnerTrack_OuterPoint_py, &b_Muon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_pz", Muon_InnerTrack_OuterPoint_pz, &b_Muon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_trackIso", Muon_trackIso, &b_Muon_trackIso);
    fChain->SetBranchAddress("Muon_ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
    fChain->SetBranchAddress("Muon_hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
    fChain->SetBranchAddress("Muon_relIso", Muon_relIso, &b_Muon_relIso);
    fChain->SetBranchAddress("Muon_normChi2", Muon_normChi2, &b_Muon_normChi2);
    fChain->SetBranchAddress("Muon_validHits", Muon_validHits, &b_Muon_validHits);
    fChain->SetBranchAddress("Muon_tkHits", Muon_tkHits, &b_Muon_tkHits);
    fChain->SetBranchAddress("Muon_pixHits", Muon_pixHits, &b_Muon_pixHits);
    fChain->SetBranchAddress("Muon_numberOfMatches", Muon_numberOfMatches, &b_Muon_numberOfMatches);
    fChain->SetBranchAddress("Muon_OuterPoint_x", Muon_OuterPoint_x, &b_Muon_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterPoint_y", Muon_OuterPoint_y, &b_Muon_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterPoint_z", Muon_OuterPoint_z, &b_Muon_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerPoint_x", Muon_InnerPoint_x, &b_Muon_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerPoint_y", Muon_InnerPoint_y, &b_Muon_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerPoint_z", Muon_InnerPoint_z, &b_Muon_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_n", &CosmicMuon_n, &b_CosmicMuon_n);
    fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px, &b_CosmicMuon_px);
    fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py, &b_CosmicMuon_py);
    fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz, &b_CosmicMuon_pz);
    fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt, &b_CosmicMuon_pt);
    fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta, &b_CosmicMuon_eta);
    fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi, &b_CosmicMuon_phi);
    fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy, &b_CosmicMuon_energy);
    fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge, &b_CosmicMuon_charge);
    fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon, &b_CosmicMuon_isGlobalMuon);
    fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon, &b_CosmicMuon_isTrackerMuon);
    fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon, &b_CosmicMuon_isStandAloneMuon);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull, &b_CosmicMuon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull, &b_CosmicMuon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_x", CosmicMuon_OuterTrack_InnerPoint_x, &b_CosmicMuon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_y", CosmicMuon_OuterTrack_InnerPoint_y, &b_CosmicMuon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_z", CosmicMuon_OuterTrack_InnerPoint_z, &b_CosmicMuon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_px", CosmicMuon_OuterTrack_InnerPoint_px, &b_CosmicMuon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_py", CosmicMuon_OuterTrack_InnerPoint_py, &b_CosmicMuon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_pz", CosmicMuon_OuterTrack_InnerPoint_pz, &b_CosmicMuon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_x", CosmicMuon_OuterTrack_OuterPoint_x, &b_CosmicMuon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_y", CosmicMuon_OuterTrack_OuterPoint_y, &b_CosmicMuon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_z", CosmicMuon_OuterTrack_OuterPoint_z, &b_CosmicMuon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_px", CosmicMuon_OuterTrack_OuterPoint_px, &b_CosmicMuon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_py", CosmicMuon_OuterTrack_OuterPoint_py, &b_CosmicMuon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_pz", CosmicMuon_OuterTrack_OuterPoint_pz, &b_CosmicMuon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_x", CosmicMuon_InnerTrack_InnerPoint_x, &b_CosmicMuon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_y", CosmicMuon_InnerTrack_InnerPoint_y, &b_CosmicMuon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_z", CosmicMuon_InnerTrack_InnerPoint_z, &b_CosmicMuon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_px", CosmicMuon_InnerTrack_InnerPoint_px, &b_CosmicMuon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_py", CosmicMuon_InnerTrack_InnerPoint_py, &b_CosmicMuon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_pz", CosmicMuon_InnerTrack_InnerPoint_pz, &b_CosmicMuon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_x", CosmicMuon_InnerTrack_OuterPoint_x, &b_CosmicMuon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_y", CosmicMuon_InnerTrack_OuterPoint_y, &b_CosmicMuon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_z", CosmicMuon_InnerTrack_OuterPoint_z, &b_CosmicMuon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_px", CosmicMuon_InnerTrack_OuterPoint_px, &b_CosmicMuon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_py", CosmicMuon_InnerTrack_OuterPoint_py, &b_CosmicMuon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_pz", CosmicMuon_InnerTrack_OuterPoint_pz, &b_CosmicMuon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_x", CosmicMuon_OuterPoint_x, &b_CosmicMuon_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_y", CosmicMuon_OuterPoint_y, &b_CosmicMuon_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_z", CosmicMuon_OuterPoint_z, &b_CosmicMuon_OuterPoint_z);
    fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
    fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
    fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
    fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
    fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
    fChain->SetBranchAddress("Photon_theta", Photon_theta, &b_Photon_theta);
    fChain->SetBranchAddress("Photon_et", Photon_et, &b_Photon_et);
    fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross, &b_Photon_swissCross);
    fChain->SetBranchAddress("Photon_e6e2", Photon_e6e2, &b_Photon_e6e2);
    fChain->SetBranchAddress("Photon_e4e1", Photon_e4e1, &b_Photon_e4e1);
    fChain->SetBranchAddress("Photonr9", Photonr9, &b_Photonr9);
    fChain->SetBranchAddress("Photon_e1x5", Photon_e1x5, &b_Photon_e1x5);
    fChain->SetBranchAddress("Photon_e2x5", Photon_e2x5, &b_Photon_e2x5);
    fChain->SetBranchAddress("Photon_e3x3", Photon_e3x3, &b_Photon_e3x3);
    fChain->SetBranchAddress("Photon_e5x5", Photon_e5x5, &b_Photon_e5x5);
    fChain->SetBranchAddress("Photon_r1x5", Photon_r1x5, &b_Photon_r1x5);
    fChain->SetBranchAddress("Photon_r2x5", Photon_r2x5, &b_Photon_r2x5);
    fChain->SetBranchAddress("Photon_maxEnergyXtal", Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
    fChain->SetBranchAddress("Photon_SigmaEtaEta", Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
    fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
    fChain->SetBranchAddress("Photon_SigmaEtaPhi", Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
    fChain->SetBranchAddress("Photon_SigmaIetaIphi", Photon_SigmaIetaIphi, &b_Photon_SigmaIetaIphi);
    fChain->SetBranchAddress("Photon_SigmaPhiPhi", Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
    fChain->SetBranchAddress("Photon_SigmaIphiIphi", Photon_SigmaIphiIphi, &b_Photon_SigmaIphiIphi);
    fChain->SetBranchAddress("Photon_Roundness", Photon_Roundness, &b_Photon_Roundness);
    fChain->SetBranchAddress("Photon_Angle", Photon_Angle, &b_Photon_Angle);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
    fChain->SetBranchAddress("Photon_isEB", Photon_isEB, &b_Photon_isEB);
    fChain->SetBranchAddress("Photon_isEE", Photon_isEE, &b_Photon_isEE);
    fChain->SetBranchAddress("Photon_isEBGap", Photon_isEBGap, &b_Photon_isEBGap);
    fChain->SetBranchAddress("Photon_isEEGap", Photon_isEEGap, &b_Photon_isEEGap);
    fChain->SetBranchAddress("Photon_isEBEEGap", Photon_isEBEEGap, &b_Photon_isEBEEGap);
    fChain->SetBranchAddress("Photon_e2e9", Photon_e2e9, &b_Photon_e2e9);
    fChain->SetBranchAddress("Photon_HoE", Photon_HoE, &b_Photon_HoE);
    fChain->SetBranchAddress("Photon_HoEnew", Photon_HoEnew, &b_Photon_HoEnew);
    fChain->SetBranchAddress("Photon_px", Photon_px, &b_Photon_px);
    fChain->SetBranchAddress("Photon_py", Photon_py, &b_Photon_py);
    fChain->SetBranchAddress("Photon_pz", Photon_pz, &b_Photon_pz);
    fChain->SetBranchAddress("Photon_vx", Photon_vx, &b_Photon_vx);
    fChain->SetBranchAddress("Photon_vy", Photon_vy, &b_Photon_vy);
    fChain->SetBranchAddress("Photon_vz", Photon_vz, &b_Photon_vz);
    fChain->SetBranchAddress("Photon_no_of_basic_clusters", Photon_no_of_basic_clusters, &b_Photon_no_of_basic_clusters);
    fChain->SetBranchAddress("Photon_sc_energy", Photon_sc_energy, &b_Photon_sc_energy);
    fChain->SetBranchAddress("Photon_sc_eta", Photon_sc_eta, &b_Photon_sc_eta);
    fChain->SetBranchAddress("Photon_sc_phi", Photon_sc_phi, &b_Photon_sc_phi);
    fChain->SetBranchAddress("Photon_sc_x", Photon_sc_x, &b_Photon_sc_x);
    fChain->SetBranchAddress("Photon_sc_y", Photon_sc_y, &b_Photon_sc_y);
    fChain->SetBranchAddress("Photon_sc_z", Photon_sc_z, &b_Photon_sc_z);
    fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
    fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
    fChain->SetBranchAddress("Photon_sc_et", Photon_sc_et, &b_Photon_sc_et);
    fChain->SetBranchAddress("matchphotonE", matchphotonE, &b_matchphotonE);
    fChain->SetBranchAddress("matchphotonpt", matchphotonpt, &b_matchphotonpt);
    fChain->SetBranchAddress("matchphotoneta", matchphotoneta, &b_matchphotoneta);
    fChain->SetBranchAddress("matchphotonphi", matchphotonphi, &b_matchphotonphi);
    fChain->SetBranchAddress("matchphotonpx", matchphotonpx, &b_matchphotonpx);
    fChain->SetBranchAddress("matchphotonpy", matchphotonpy, &b_matchphotonpy);
    fChain->SetBranchAddress("matchphotonpz", matchphotonpz, &b_matchphotonpz);
    fChain->SetBranchAddress("ismatchedphoton", ismatchedphoton, &b_ismatchedphoton);
    fChain->SetBranchAddress("Photon_hasConvTrk", Photon_hasConvTrk, &b_Photon_hasConvTrk);
    fChain->SetBranchAddress("Photon_ntracks", Photon_ntracks, &b_Photon_ntracks);
    fChain->SetBranchAddress("Photon_isconverted", Photon_isconverted, &b_Photon_isconverted);
    fChain->SetBranchAddress("Photon_pairInvmass", Photon_pairInvmass, &b_Photon_pairInvmass);
    fChain->SetBranchAddress("Photon_pairCotThetaSeperation", Photon_pairCotThetaSeperation, &b_Photon_pairCotThetaSeperation);
    fChain->SetBranchAddress("Photon_pairmomentumX", Photon_pairmomentumX, &b_Photon_pairmomentumX);
    fChain->SetBranchAddress("Photon_pairmomentumY", Photon_pairmomentumY, &b_Photon_pairmomentumY);
    fChain->SetBranchAddress("Photon_pairmomentumZ", Photon_pairmomentumZ, &b_Photon_pairmomentumZ);
    fChain->SetBranchAddress("Photon_EoverP", Photon_EoverP, &b_Photon_EoverP);
    fChain->SetBranchAddress("Photon_ConvVx", Photon_ConvVx, &b_Photon_ConvVx);
    fChain->SetBranchAddress("Photon_ConvVy", Photon_ConvVy, &b_Photon_ConvVy);
    fChain->SetBranchAddress("Photon_ConvVz", Photon_ConvVz, &b_Photon_ConvVz);
    fChain->SetBranchAddress("Photon_ZOfPrimaryVertex", Photon_ZOfPrimaryVertex, &b_Photon_ZOfPrimaryVertex);
    fChain->SetBranchAddress("Photon_distOfMinimumApproach", Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
    fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
    fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
    fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
    fChain->SetBranchAddress("npho", &npho, &b_npho);
    fChain->SetBranchAddress("Photon_Electronveto", Photon_Electronveto, &b_Photon_Electronveto);
    fChain->SetBranchAddress("PFiso_Charged03", PFiso_Charged03, &b_PFiso_Charged03);
    fChain->SetBranchAddress("PFiso_Photon03", PFiso_Photon03, &b_PFiso_Photon03);
    fChain->SetBranchAddress("PFiso_Neutral03", PFiso_Neutral03, &b_PFiso_Neutral03);
    fChain->SetBranchAddress("PFiso_Sum03", PFiso_Sum03, &b_PFiso_Sum03);
    fChain->SetBranchAddress("PFWorstiso_Charged03", PFWorstiso_Charged03, &b_PFWorstiso_Charged03);
    fChain->SetBranchAddress("Photon_ncrys", Photon_ncrys, &b_Photon_ncrys);
    fChain->SetBranchAddress("Photon_timing_xtal", Photon_timing_xtal, &b_Photon_timing_xtal);
    fChain->SetBranchAddress("Photon_timingavg_xtal", Photon_timingavg_xtal, &b_Photon_timingavg_xtal);
    fChain->SetBranchAddress("Photon_energy_xtal", Photon_energy_xtal, &b_Photon_energy_xtal);
    fChain->SetBranchAddress("Photon_ieta_xtalEB", Photon_ieta_xtalEB, &b_Photon_ieta_xtalEB);
    fChain->SetBranchAddress("Photon_iphi_xtalEB", Photon_iphi_xtalEB, &b_Photon_iphi_xtalEB);
    fChain->SetBranchAddress("Photon_recoFlag_xtalEB", Photon_recoFlag_xtalEB, &b_Photon_recoFlag_xtalEB);
    fChain->SetBranchAddress("Photon_timeError_xtal", Photon_timeError_xtal, &b_Photon_timeError_xtal);
    fChain->SetBranchAddress("Photon_s9", Photon_s9, &b_Photon_s9);
    fChain->SetBranchAddress("Photon_mipChi2", Photon_mipChi2, &b_Photon_mipChi2);
    fChain->SetBranchAddress("Photon_mipTotEnergy", Photon_mipTotEnergy, &b_Photon_mipTotEnergy);
    fChain->SetBranchAddress("Photon_mipSlope", Photon_mipSlope, &b_Photon_mipSlope);
    fChain->SetBranchAddress("Photon_mipIntercept", Photon_mipIntercept, &b_Photon_mipIntercept);
    fChain->SetBranchAddress("Photon_mipNhitCone", Photon_mipNhitCone, &b_Photon_mipNhitCone);
    fChain->SetBranchAddress("Photon_mipIsHalo", Photon_mipIsHalo, &b_Photon_mipIsHalo);
    fChain->SetBranchAddress("EBRecHit_size", &EBRecHit_size, &b_EBRecHit_size);
    fChain->SetBranchAddress("EBRecHit_eta", EBRecHit_eta, &b_EBRecHit_eta);
    fChain->SetBranchAddress("EBRecHit_phi", EBRecHit_phi, &b_EBRecHit_phi);
    fChain->SetBranchAddress("EBRecHit_ieta", EBRecHit_ieta, &b_EBRecHit_ieta);
    fChain->SetBranchAddress("EBRecHit_iphi", EBRecHit_iphi, &b_EBRecHit_iphi);
    fChain->SetBranchAddress("EBRecHit_e", EBRecHit_e, &b_EBRecHit_e);
    fChain->SetBranchAddress("EBRecHit_et", EBRecHit_et, &b_EBRecHit_et);
    fChain->SetBranchAddress("EBRecHit_flag", EBRecHit_flag, &b_EBRecHit_flag);
    fChain->SetBranchAddress("EBRecHit_time", EBRecHit_time, &b_EBRecHit_time);
    fChain->SetBranchAddress("EERecHit_size", &EERecHit_size, &b_EERecHit_size);
    fChain->SetBranchAddress("EERecHit_eta", EERecHit_eta, &b_EERecHit_eta);
    fChain->SetBranchAddress("EERecHit_phi", EERecHit_phi, &b_EERecHit_phi);
    fChain->SetBranchAddress("EERecHit_ieta", EERecHit_ieta, &b_EERecHit_ieta);
    fChain->SetBranchAddress("EERecHit_iphi", EERecHit_iphi, &b_EERecHit_iphi);
    fChain->SetBranchAddress("EERecHit_e", EERecHit_e, &b_EERecHit_e);
    fChain->SetBranchAddress("EERecHit_et", EERecHit_et, &b_EERecHit_et);
    fChain->SetBranchAddress("EERecHit_flag", EERecHit_flag, &b_EERecHit_flag);
    fChain->SetBranchAddress("EERecHit_time", EERecHit_time, &b_EERecHit_time);
    fChain->SetBranchAddress("PFMetPt", PFMetPt, &b_PFMetPt);
    fChain->SetBranchAddress("PFMetPx", PFMetPx, &b_PFMetPx);
    fChain->SetBranchAddress("PFMetPy", PFMetPy, &b_PFMetPy);
    fChain->SetBranchAddress("PFMetPhi", PFMetPhi, &b_PFMetPhi);
    fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt, &b_PFMetSumEt);
    fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF, &b_Delta_phiPF);
    fChain->SetBranchAddress("CaloTower_n", &CaloTower_n, &b_CaloTower_n);
    fChain->SetBranchAddress("CaloTower_eta", CaloTower_eta, &b_CaloTower_eta);
    fChain->SetBranchAddress("CaloTower_phi", CaloTower_phi, &b_CaloTower_phi);
    fChain->SetBranchAddress("CaloTower_E", CaloTower_E, &b_CaloTower_E);
    fChain->SetBranchAddress("CaloTower_Et", CaloTower_Et, &b_CaloTower_Et);
    fChain->SetBranchAddress("CaloTower_emEnergy", CaloTower_emEnergy, &b_CaloTower_emEnergy);
    fChain->SetBranchAddress("CaloTower_hadEnergy", CaloTower_hadEnergy, &b_CaloTower_hadEnergy);
    fChain->SetBranchAddress("CaloTower_p", CaloTower_p, &b_CaloTower_p);
    fChain->SetBranchAddress("CaloTower_EMEt", CaloTower_EMEt, &b_CaloTower_EMEt);
    fChain->SetBranchAddress("CaloTower_HadEt", CaloTower_HadEt, &b_CaloTower_HadEt);
    fChain->SetBranchAddress("CaloTower_HadPhi", CaloTower_HadPhi, &b_CaloTower_HadPhi);
    fChain->SetBranchAddress("CaloTower_HadEta", CaloTower_HadEta, &b_CaloTower_HadEta);
    fChain->SetBranchAddress("CaloTower_EMPhi", CaloTower_EMPhi, &b_CaloTower_EMPhi);
    fChain->SetBranchAddress("CaloTower_EMEta", CaloTower_EMEta, &b_CaloTower_EMEta);
    fChain->SetBranchAddress("CaloTower_HadX", CaloTower_HadX, &b_CaloTower_HadX);
    fChain->SetBranchAddress("CaloTower_HadY", CaloTower_HadY, &b_CaloTower_HadY);
    fChain->SetBranchAddress("CaloTower_HadZ", CaloTower_HadZ, &b_CaloTower_HadZ);
    fChain->SetBranchAddress("CaloTower_HE_E", CaloTower_HE_E, &b_CaloTower_HE_E);
    fChain->SetBranchAddress("CaloTower_HB_E", CaloTower_HB_E, &b_CaloTower_HB_E);
    fChain->SetBranchAddress("CaloTower_EMTime", CaloTower_EMTime, &b_CaloTower_EMTime);
    fChain->SetBranchAddress("CaloTower_HadTime", CaloTower_HadTime, &b_CaloTower_HadTime);
    fChain->SetBranchAddress("rho", &rho, &b_rho);
    fChain->SetBranchAddress("sigma", &sigma, &b_sigma);
    fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
    fChain->SetBranchAddress("sigma25", &sigma25, &b_sigma25);
    Notify();
}

Bool_t PostAnalyzerData::Notify(){
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void PostAnalyzerData::Show(Long64_t entry){
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t PostAnalyzerData::Cut(Long64_t entry){
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}


// HLT for data
Bool_t PostAnalyzerData::PassHLT(Bool_t &passedHLT150, Bool_t &passed90HLTall){
    Bool_t passedHLT = false;
    Int_t HLTcounter = 0;
    const int SelectedTriggers = 10;

    TString HLTPhotonTriggers[SelectedTriggers]={
        "HLT_Photon150_v1",
        "HLT_Photon150_v2",
        "HLT_Photon150_v3",
        "HLT_Photon150_v4",
        "HLT_Photon150_v5",
        "HLT_Photon150_v6",
        "HLT_Photon150_v7",
        "HLT_Photon150_v8",
	"HLT_Photon150_v9",
        "HLT_Photon150_v10"
    };

    for(Int_t p = 0; p < (*triggernames).size(); ++p){
        for(Int_t k = 0; k < SelectedTriggers; ++k){
            if(HLTPhotonTriggers[k] == ((*triggernames)[p])){
		if( (*triggerprescales)[p] == 1 && ((*ifTriggerpassed)[p]) ){
                    HLTcounter++;
                }
            }
        }
    }
    if(HLTcounter > 0) passedHLT = true;
    return passedHLT;
}

//---------------------
// Scraping variable
Bool_t PostAnalyzerData::NonScraping(){
    Bool_t noScraping = !Scraping_isScrapingEvent ;
    return noScraping;
}

//--------------------
// Spike Cuts
Bool_t PostAnalyzerData::NoSpike(Int_t ipho){
    Bool_t passSpike = false;

    if(  fabs(getLICTD(ipho))                 < 5.0    &&
            fabs(Photon_timing_xtal[ipho][0]) < 3.0    &&
            Photon_SigmaIetaIeta[ipho]        > 0.001  &&
            Photon_SigmaIphiIphi[ipho]        > 0.001  &&
            Photonr9[ipho]                    < 1.0){
        passSpike = true ;
    }

    return passSpike;

}
 
//----------------
//Primary Vertex
Bool_t PostAnalyzerData::PrimaryVertex(Int_t &goodVertex){
    Bool_t VertexAccepted = false;
    goodVertex=0;

    for(Int_t i=0; i<Vertex_n && i<200;++i){
        if(  fabs(Vertex_z[i]) <= Cvertex_z    &&
                Vertex_ndof[i]    >= Cvertex_ndof &&
                !Vertex_isFake[i]                 &&
                fabs(Vertex_d0[i])<= Cvertex_rho)

            goodVertex++ ;
    }
    if(goodVertex > 0)VertexAccepted=true;

    return VertexAccepted;
}

//-----------------
// Compute Rapidity
Double_t PostAnalyzerData::getRapidity(Double_t r_E, Double_t r_Pz){
    Double_t rapdty  =  0.5*log( (r_E + r_Pz)/(r_E - r_Pz) );

    return rapdty ;
}

//-----------------
// Compute DeltaR
Double_t PostAnalyzerData::getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

    Double_t dPhi  = fabs(phi1 - phi2);
    Double_t twopi = 2.0*(TMath::Pi());
    Double_t dEta = fabs(eta1 - eta2);
    Double_t DR=0.0;

    if(dPhi < 0) dPhi = - dPhi;
    if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;

    DR = pow((dPhi*dPhi + dEta*dEta),0.5);
    return DR;
}

//-----------------
//Compute Delta Eta
Double_t PostAnalyzerData::getDEta(Double_t eta1, Double_t eta2){

    Double_t dEta = fabs(eta1 - eta2);
    return dEta;
}

//------------------
//Compute Delta Phi
Double_t PostAnalyzerData::getDPhi(Double_t phi1, Double_t phi2){

    Double_t dPhi  = fabs(phi1 - phi2);
    Double_t twopi = 2.0*(TMath::Pi());

    if(dPhi < 0) dPhi = - dPhi;
    if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;

    return dPhi;
}

//----------------
//Compute Invariant mass for Gamma and Jet
Double_t PostAnalyzerData::getMass(Int_t pho_i, Int_t jet_i){

    Double_t mass=0.0;

    Double_t E  = Photon_E[pho_i] + pfJet_E[jet_i];
    Double_t PX = Photon_px[pho_i] + pfJet_px[jet_i];
    Double_t PY = Photon_py[pho_i] + pfJet_py[jet_i];
    Double_t PZ = Photon_pz[pho_i] + pfJet_pz[jet_i];

    mass = pow((E*E - PX*PX - PY*PY - PZ*PZ),0.5);

    return mass;
}

//-----------------
// Compute LICTD LargeIntraClusterTimeDifference
Double_t PostAnalyzerData::getLICTD(Int_t i){

    Double_t SeedTime = -999;
    Double_t SeedE    = -999;

    Int_t CrysIdx     = -1;

    for(Int_t k=0; k<Photon_ncrys[i]&&k<100;++k){
	Float_t CrysE = Photon_energy_xtal[i][k];
        if(CrysE > SeedE){
            SeedE    = CrysE ;
            SeedTime = Photon_timing_xtal[i][k];
            CrysIdx  = k;
        }
    }

    Float_t LICTD = 99.0;

    if(fabs(SeedTime) < 3.0){
	LICTD = 0.0;
        Int_t CrysCrys   = -1;
        Int_t CrysThresh = 0;

	for(Int_t k=0;k<Photon_ncrys[i]&&k<100;++k){
            if(CrysIdx == k)continue;
            Float_t CrysE = Photon_energy_xtal[i][k];

            if(CrysE > 1.0){
		CrysThresh++;
		Float_t timeDiff = Photon_timing_xtal[i][CrysIdx] - Photon_timing_xtal[i][k];
                if(fabs(timeDiff) > fabs(LICTD)){
                    LICTD    = timeDiff;
                    CrysCrys = k;
		}
            }
        }
    }

    return LICTD;
}


//Tight Photon PF ISo for 2012 to select photon candidate
Bool_t PostAnalyzerData::TightPhotonPFIso(Int_t ipho){

    Bool_t tightID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
        tightID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.011)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 0.7)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 0.4+0.04*Photon_pt[ipho])       &&
            (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 0.5+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
        tightID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.031)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 0.5)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.5+0.04*Photon_pt[ipho])       &&
            (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);

    }
    return tightID;
}

Bool_t PostAnalyzerData::MediumPhotonPFIso(Int_t ipho){

    Bool_t mediumID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
        mediumID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.011)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 1.5)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.0+0.04*Photon_pt[ipho])       &&
            (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 0.7+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
        mediumID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.033)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 1.2)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.5+0.04*Photon_pt[ipho])       &&
            (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);

    }
    return mediumID;
}

Bool_t PostAnalyzerData::LoosePhotonPFIso(Int_t ipho){

    Bool_t looseID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
        looseID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.012)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 2.6)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 3.5+0.04*Photon_pt[ipho])       &&
            (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.3+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
        looseID = (Photon_HoEnew[ipho]   < 0.05)          &&
            (Photon_SigmaIetaIeta[ipho]  < 0.034)         &&
            (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 2.3)                            &&
            (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 2.9+0.04*Photon_pt[ipho])       &&
            //      (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
            (Photon_Electronveto[ipho]  == 1);

    }
    return looseID;
}


// Tight Jet ID and jet seperation from photon to select an isolated jet.
Bool_t PostAnalyzerData::TightJetID( Int_t ijet){

    Bool_t ID = false;
    if(fabs(pfJet_eta[ijet]) <= 2.4){
        ID = (pfjet_NEF[ijet]           < 0.90) &&
            (pfjet_NHF[ijet]            < 0.90) &&
            (pfjet_NConstituents[ijet]  > 1)    &&
            (pfjet_CEF[ijet]            < 0.99) &&
            (pfjet_CHF[ijet]            > 0)    &&
            (pfjet_NCH[ijet]            > 0);
    }

    if(fabs(pfJet_eta[ijet]) > 2.4){
	ID = (pfjet_NEF[ijet]           < 0.90) &&
            (pfjet_NHF[ijet]            < 0.90) &&
            (pfjet_NConstituents[ijet]  > 1);
    }
    return ID;
}

// Effective area to be needed in PF Iso for photon ID
Double_t PostAnalyzerData::EAElectroncharged(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.012;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.010;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.014;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.012;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.016;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.020;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.012;

    return EffectiveArea;
}

Double_t PostAnalyzerData::EAElectronneutral(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.030;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.057;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.039;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.015;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.024;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.039;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.072;

    return EffectiveArea;
}

Double_t PostAnalyzerData::EAElectronphoton(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.148;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.130;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.112;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.216;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.262;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.260;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.266;

    return EffectiveArea;
}

Double_t PostAnalyzerData::EXOtightPhoID(Int_t ipho){

    Bool_t ID = false;

    if(fabs(Photon_sc_eta[ipho]) < 1.4442 ){
        ID = (  Photon_HoE[ipho] < 0.05                                                         &&
                Photon_ecalRecHitSumEtConeDR04[ipho] < 4.2 + 0.006*Photon_pt[ipho]  + 0.1830*rho25 &&
        	Photon_hcalTowerSumEtConeDR04[ipho]  < 2.2 + 0.0025*Photon_pt[ipho] + 0.0620*rho25 &&
                Photon_trkSumPtHollowConeDR04[ipho]  < 2.0 + 0.001*Photon_pt[ipho]  + 0.0167*rho25 &&
                Photon_SigmaIetaIeta[ipho]           < 0.011                                    &&
        	(!Photon_hasPixelSeed[ipho])
             );
    }
    if(fabs(Photon_sc_eta[ipho]) > 1.4442 ){
        ID = (  Photon_HoE[ipho] < 0.05                                                      &&
                Photon_ecalRecHitSumEtConeDR04[ipho] < 4.2+0.006*Photon_pt[ipho]  + 0.090*rho25 &&
                Photon_hcalTowerSumEtConeDR04[ipho]  < 2.2+0.0025*Photon_pt[ipho] + 0.180*rho25 &&
                Photon_trkSumPtHollowConeDR04[ipho]  < 2.0+0.001*Photon_pt[ipho]  + 0.032*rho25 &&
                Photon_SigmaIetaIeta[ipho]           < 0.030                                 &&
                (!Photon_hasPixelSeed[ipho])
             );
    }

    return ID;
}

void PostAnalyzerData::BookHistos(){

  // Define Histograms here +++++++++++++++++++++++++++++
  f1->cd();

  h_ptPhoton  = new TH1F("h_ptPhoton_final","pt of photon",100,20.0,2520.0);  // As we are putting a cut of 170.
  h_ptPhoton->GetYaxis()->SetTitle("Events/25 GeV");           h_ptPhoton->GetYaxis()->CenterTitle();
  h_ptPhoton->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton->GetXaxis()->CenterTitle();
  h_ptPhoton->Sumw2();

  h_ptJet = new TH1F("h_ptJet_final","pt of jet      ",100,20.0,2520.0); // As we are putting a cut of 170.
  h_ptJet->GetYaxis()->SetTitle("Events/25 GeV");       h_ptJet->GetYaxis()->CenterTitle();
  h_ptJet->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet->GetXaxis()->CenterTitle();
  h_ptJet->Sumw2();

  h_mass_bin25 = new TH1F("h_mass_bin25_final","mass of photon+jet",160,0.0,4000.0);
  h_mass_bin25->GetYaxis()->SetTitle("Events/25 GeV");               h_mass_bin25->GetYaxis()->CenterTitle();
  h_mass_bin25->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_bin25->GetXaxis()->CenterTitle();
  h_mass_bin25->Sumw2();

  h_etaPhoton =new TH1F("h_etaPhoton_final","eta of photon",100,-2.5,2.5);
  h_etaPhoton->GetYaxis()->SetTitle("Events");        h_etaPhoton->GetYaxis()->CenterTitle();
  h_etaPhoton->GetXaxis()->SetTitle("#eta^{#gamma}"); h_etaPhoton->GetXaxis()->CenterTitle();
  h_etaPhoton->Sumw2();

  h_etaJet = new TH1F("h_etaJet_final","eta of jet ",200,-5.0,5.0);
  h_etaJet->GetYaxis()->SetTitle("Events");      h_etaJet->GetYaxis()->CenterTitle();
  h_etaJet->GetXaxis()->SetTitle("#eta^{jet}");  h_etaJet->GetXaxis()->CenterTitle();
  h_etaJet->Sumw2();
 
  h_Photon_SigmaIetaIeta = new TH1F("h_Photon_SigmaIetaIeta_final","Photon SigmaIetaIeta",100,0.0,0.05);
  h_Photon_SigmaIetaIeta->GetYaxis()->SetTitle("Events");                 h_Photon_SigmaIetaIeta->GetYaxis()->CenterTitle();
  h_Photon_SigmaIetaIeta->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");   h_Photon_SigmaIetaIeta->GetXaxis()->CenterTitle();
  h_Photon_SigmaIetaIeta->Sumw2();

  h_DR_PhotonJet = new TH1F("h_DR_PhotonJet_final","DeltaR between photon n jet",100,0.0,10.0);
  h_DR_PhotonJet->GetYaxis()->SetTitle("Events");        h_DR_PhotonJet->GetYaxis()->CenterTitle();
  h_DR_PhotonJet->GetXaxis()->SetTitle("#Delta R");      h_DR_PhotonJet->GetXaxis()->CenterTitle();
  h_DR_PhotonJet->Sumw2();

  h_dEta = new TH1F("h_dEta_final","dEta of photon+jet",120,0,6);
  h_dEta->GetYaxis()->SetTitle("Events");        h_dEta->GetYaxis()->CenterTitle();
  h_dEta->GetXaxis()->SetTitle("#Delta #eta");   h_dEta->GetXaxis()->CenterTitle();
  h_dEta->Sumw2();

  h_dphi = new TH1F("h_dphi_final","dphi of photon+jet",64,0,3.2);
  h_dphi->GetYaxis()->SetTitle("Events");        h_dphi->GetYaxis()->CenterTitle();
  h_dphi->GetXaxis()->SetTitle("#Delta #phi");   h_dphi->GetXaxis()->CenterTitle();
  h_dphi->Sumw2();

  h_nPhoton = new TH1F("h_nPhoton_final","no of Photons",10,0,10);

  h_nJet = new TH1F("h_nJet_final","no of Jets",20,0,20);

  h_PC = new TH1F ("h_PC","Photon Candidate", 10, 0, 10);

  h_JC = new TH1F ("h_JC","Photon Candidate", 20, 0, 20);

}


#endif // #ifdef PostAnalyzerData_cxx



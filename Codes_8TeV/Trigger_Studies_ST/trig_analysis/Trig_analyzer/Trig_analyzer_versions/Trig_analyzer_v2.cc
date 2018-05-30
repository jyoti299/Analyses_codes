// -*- C++ -*-
//
// Package:    Trig_analyzer
// Class:      Trig_analyzer
// 
/**\class Trig_analyzer Trig_analyzer.cc trig_analysis/Trig_analyzer/src/Trig_analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  rocky garg
//         Created:  Sat Mar  1 00:10:14 CST 2014
// $Id$
//
//


//User Defined Files
#include <iostream>
#include "Trig_analyzer.h"



using namespace std;
using namespace ROOT;
using namespace edm;


//
// constants, enums and typedefs 
//

//
// static data member definitions
//


//Class for sorting of variables
class PtSortCriteriumElectron{
public:
  bool operator() (reco::GsfElectron p1, reco::GsfElectron p2){
    return p1.pt() > p2.pt();
  }
};




//
// constructors 
//

Trig_analyzer::Trig_analyzer(const edm::ParameterSet& iConfig):
  filterName_(iConfig.getUntrackedParameter<std::string>("filterName", "hltPhoton135HEFilter")),
  triggerPathName_(iConfig.getUntrackedParameter<std::string>("triggerPathName", "")),
  gsfelectron_(iConfig.getUntrackedParameter<edm::InputTag>("gsfelectronTag")),
  triggerEvent_(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
  triggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag")),
  ConvCollLabel_(iConfig.getUntrackedParameter<edm::InputTag>("ConvCollTag")),
  BSLabel_(iConfig.getUntrackedParameter<edm::InputTag>("BSTag")),
  primaryVertexLabel_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertexTag")),
  pfCandidateLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateTag")),
  rhoLabel25_(iConfig.getUntrackedParameter<edm::InputTag>("rhoInputTag25")),
  EAtargetToken_(iConfig.getUntrackedParameter<std::string>("EAtarget")),
  pTOfflineProbeCut_(iConfig.getUntrackedParameter<double>("pTOfflineProbeCut", 0.0))
{
  
  nEvents = 0;

  //PFIsolation Initialization for Electrons
  eleIsolator03.initializeElectronIsolation(kTRUE); //NOTE: this automatically set all the correct default veto values
  eleIsolator03.setConeSize(0.3);



}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
Trig_analyzer::beginJob()
{
  f = new TFile("Trig_analyzer.root", "recreate");
  t = new TTree("T", "Trigger_tree");
 
  //Declearing Branches for trigger Objects
  t->Branch("TrigObj_Pt", &trigObj_pt);
  t->Branch("TrigObj_P", &trigObj_p);
  t->Branch("TrigObj_Px", &trigObj_px);
  t->Branch("TrigObj_Py", &trigObj_py);
  t->Branch("TrigObj_Pz", &trigObj_pz);
  t->Branch("TrigObj_Energy", &trigObj_energy);
  t->Branch("TrigObj_Et", &trigObj_et);
  t->Branch("TrigObj_Eta", &trigObj_eta);
  t->Branch("TrigObj_Phi", &trigObj_phi);
  t->Branch("TrigObj_Mass", &trigObj_mass);

  //Declearing Branches for Electron Variables
  t->Branch("Electron_n", &Electron_n, "Electron_n/I");
  //Kinematic variables (no cuts)
  t->Branch("Electron_nocut_E", &Electron_nocut_E);
  t->Branch("Electron_nocut_et", &Electron_nocut_et);
  t->Branch("Electron_nocut_pt", &Electron_nocut_pt);
  t->Branch("Electron_nocut_eta", &Electron_nocut_eta);
  t->Branch("Electron_nocut_phi", &Electron_nocut_phi);
  t->Branch("Electron_nocut_theta", &Electron_nocut_theta);
  t->Branch("Electron_nocut_px", &Electron_nocut_px);
  t->Branch("Electron_nocut_py", &Electron_nocut_py);
  t->Branch("Electron_nocut_pz", &Electron_nocut_pz);
  t->Branch("Electron_nocut_charge", &Electron_nocut_charge);
  //Id variables (no cuts)
  t->Branch("Electron_nocut_sigmaIEtaIEta", &Electron_nocut_sigmaIEtaIEta);
  t->Branch("Electron_nocut_dEtaIn", &Electron_nocut_dEtaIn);
  t->Branch("Electron_nocut_dPhiIn", &Electron_nocut_dPhiIn);
  t->Branch("Electron_nocut_HoE", &Electron_nocut_HoE);
  t->Branch("Electron_nocut_ooemoop", &Electron_nocut_ooemoop);
  t->Branch("Electron_nocut_eopIn", &Electron_nocut_eopIn);
  t->Branch("Electron_nocut_fBrem", &Electron_nocut_fBrem);
  //Impact Parameter Variables (no cuts)
  t->Branch("Electron_nocut_d0vtx", &Electron_nocut_d0vtx);
  t->Branch("Electron_nocut_dzvtx", &Electron_nocut_dzvtx);
  //Conversion Variables (no cuts)
  t->Branch("Electron_nocut_vtxFitConversion", &Electron_nocut_vtxFitConversion);
  t->Branch("Electron_nocut_misHits", &Electron_nocut_misHits);
  //Detector Isolation Variables (no cuts)
  t->Branch("Electron_nocut_trackIso03", &Electron_nocut_trackIso03);
  t->Branch("Electron_nocut_ecalIso03", &Electron_nocut_ecalIso03);
  t->Branch("Electron_nocut_hcalIso03", &Electron_nocut_hcalIso03);
  //PFIsolation Variables (no cuts)
  t->Branch("Electron_nocut_PFIsoCh03", &Electron_nocut_PFIsoCh03);
  t->Branch("Electron_nocut_PFIsoNh03", &Electron_nocut_PFIsoNh03);
  t->Branch("Electron_nocut_PFIsoEM03", &Electron_nocut_PFIsoEM03);
  t->Branch("rho25", &rho25, "rho25/F"); 
  t->Branch("Electron_nocut_PFIsoSum_rhoCorr", &Electron_nocut_PFIsoSum_rhoCorr);

  //Kinematic Variables (with cuts and mass matching)
  t->Branch("Electron_allProbes_E", &Electron_allProbes_E);
  t->Branch("Electron_allProbes_et", &Electron_allProbes_et);
  t->Branch("Electron_allProbes_pt", &Electron_allProbes_pt);
  t->Branch("Electron_allProbes_eta", &Electron_allProbes_eta);
  t->Branch("Electron_allProbes_phi", &Electron_allProbes_phi);
  t->Branch("Electron_allProbes_theta", &Electron_allProbes_theta);
  t->Branch("Electron_allProbes_px", &Electron_allProbes_px);
  t->Branch("Electron_allProbes_py", &Electron_allProbes_py);
  t->Branch("Electron_allProbes_pz", &Electron_allProbes_pz);
  t->Branch("Electron_allProbes_charge", &Electron_allProbes_charge);
  //Id Variables (with cuts and mass matching)
  t->Branch("Electron_allProbes_sigmaIEtaIEta", &Electron_allProbes_sigmaIEtaIEta);
  t->Branch("Electron_allProbes_dEtaIn", &Electron_allProbes_dEtaIn);
  t->Branch("Electron_allProbes_dPhiIn", &Electron_allProbes_dPhiIn);
  t->Branch("Electron_allProbes_HoE", &Electron_allProbes_HoE);
  t->Branch("Electron_allProbes_ooemoop", &Electron_allProbes_ooemoop);
  t->Branch("Electron_allProbes_eopIn", &Electron_allProbes_eopIn);
  t->Branch("Electron_allProbes_fBrem", &Electron_allProbes_fBrem);
  //Impact Parameter Variables (with cuts and mass matching)
  t->Branch("Electron_allProbes_d0vtx", &Electron_allProbes_d0vtx);
  t->Branch("Electron_allProbes_dzvtx", &Electron_allProbes_dzvtx);


  //Declearing Histograms
  h_deltaR = new TH1F("h_deltaR", "DeltaR between RECO and HLT Object", 100, 0.0, 10.0);
  h_deltaR->GetYaxis()->SetTitle("\# of Objects");        h_deltaR->GetYaxis()->CenterTitle();
  h_deltaR->GetXaxis()->SetTitle("#Delta R");             h_deltaR->GetXaxis()->CenterTitle();
  //h_deltaR->Sumw2();

  //Here change only nbins and xbins[nbins] value if need different no. of bins and different final bin.
  //And change i*5 to i*x if required to change bin size from 5 to x.
  int nbins = 76;//31;
  float xbins[nbins+1];
  for(int i = 0; i < nbins; i++){
    xbins[i] = i*2;
  }
  xbins[nbins] = 200;
 
  h_allProbes = new TH1F("h_allProbes", "Pt distribution of all Probes", nbins, xbins);
  h_allProbes->GetYaxis()->SetTitle("Objects/2 GeV");        h_allProbes->GetYaxis()->CenterTitle();
  h_allProbes->GetXaxis()->SetTitle("P_{T}(GeV)");             h_allProbes->GetXaxis()->CenterTitle();
  //  h_allProbes->Sumw2();

  h_passProbes = new TH1F("h_passProbes", "Pt distribution of pass Probes", nbins, xbins);
  h_passProbes->GetYaxis()->SetTitle("Objects/2 GeV");          h_passProbes->GetYaxis()->CenterTitle();
  h_passProbes->GetXaxis()->SetTitle("P_{T}(GeV)");             h_passProbes->GetXaxis()->CenterTitle();
  //h_passProbes->Sumw2();


  h_efficiency = new TH1F("h_efficiency", "Efficiency vs. Pt of reco objects", nbins, xbins);
  h_efficiency->GetYaxis()->SetTitle("Efficiency/2 GeV");          h_efficiency->GetYaxis()->CenterTitle();
  h_efficiency->GetXaxis()->SetTitle("P_{T}(GeV)");          h_efficiency->GetXaxis()->CenterTitle();
  //h_efficiency->Sumw2();

  //Getting correct EAtarget for Effective area calculation
  if( EAtargetToken_ == "EleEANoCorr")
    EAtarget  =ElectronEffectiveArea::kEleEANoCorr;
  else if( EAtargetToken_ == "EleEAData2011")
    EAtarget  =ElectronEffectiveArea::kEleEAData2011;
  else if( EAtargetToken_ == "EleEASummer11MC")
    EAtarget  =ElectronEffectiveArea::kEleEASummer11MC;
  else if( EAtargetToken_ == "EleEAFall11MC")
    EAtarget  =ElectronEffectiveArea::kEleEAFall11MC;
  else if( EAtargetToken_ == "EleEAData2012")
    EAtarget  =ElectronEffectiveArea::kEleEAData2012;
  else
    EAtarget  =ElectronEffectiveArea::kEleEAData2012;

}

// ------------ method called when starting to processes a run  ------------
void 
Trig_analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if(hltConfig_.init(iRun,iSetup,"HLT",changed)){
    cout << "Table name = " << hltConfig_.tableName() << endl;
    if(changed){
    }
  }
}

// ------------ method called for each event  ------------
void
Trig_analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
  nEvents++;
  cout << "PROCESSING EVENT = " << nEvents << endl;

  Electron_n = 0;  

  Clear();
  
  //-----------------------------------------------------
  //Accessing TriggerResults and triggerNames collections
  //-----------------------------------------------------
  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  iEvent.getByLabel(triggerEvent_,triggerEventHandle);

  if (!triggerEventHandle.isValid()) {
  std::cout << "Error! Can't get the product "<< triggerEvent_.label() << std::endl;
  } else{
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResults_, triggerResults);

    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
    const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();

    //---------------------------------------------------------------------------------
    //This part just saving some trigger object variables passing a filter filterName_.
    //This has nothing to do with trigger efficiency calculation.
    //---------------------------------------------------------------------------------
        
    trigger::size_type filterIndex = triggerEventHandle->filterIndex(edm::InputTag(filterName_, "", triggerEvent_.process()));
    if(filterIndex < triggerEventHandle->sizeFilters()){
      const trigger::Keys& trigKeys = triggerEventHandle->filterKeys(filterIndex);
      const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle->getObjects();
      for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//cout << obj.id() << endl;

	trigObj_pt.push_back(obj.pt());
	trigObj_p.push_back(obj.p());
	trigObj_px.push_back(obj.px());
	trigObj_py.push_back(obj.py());
	trigObj_pz.push_back(obj.pz());
	trigObj_energy.push_back(obj.energy());
	trigObj_et.push_back(obj.et());
	trigObj_eta.push_back(obj.eta());
	trigObj_phi.push_back(obj.phi());
	trigObj_mass.push_back(obj.mass());
      } 
    }
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //-----------------------------------------------------------
    //Getting the correct version name of triggerPath triggerPathName_
    //-----------------------------------------------------------
    int triggerPathName_idx = -1;
    for(unsigned int i=0; i<triggerNames_.size(); i++){
      if(triggerNames_[i].find(triggerPathName_.c_str())!=std::string::npos){
	triggerPathName_idx = i;
	break;
      }
    }
    triggerPathName_ = triggerNames_[triggerPathName_idx];
      
    //-------------------------------------------------------------
    //Getting the list of filters with savetags = true in tighPath_
    //-------------------------------------------------------------
    if(hltConfig_.triggerIndex(triggerPathName_)<hltConfig_.size()){
      filtersOfTriggerPath_ = hltConfig_.saveTagsModules(triggerPathName_);
    }
    
    //----------------------------------------------
    //Initializing the container for trigger objects
    //----------------------------------------------
    std::vector<trigger::TriggerObject> ObjectsOfTriggerPath_;
    ObjectsOfTriggerPath_.clear();

    //------------------------------------------------------------------------------------------
    //Getting all the objects passing filters of triggerPathName_ and filling into the container
    //------------------------------------------------------------------------------------------
    for(size_t filterNr=0; filterNr<filtersOfTriggerPath_.size();filterNr++){
      trigger::size_type filterIndex = triggerEventHandle->filterIndex(edm::InputTag(filtersOfTriggerPath_[filterNr], "", triggerEvent_.process()));
      if(filterIndex < triggerEventHandle->sizeFilters()){
        const trigger::Keys& trigKeys = triggerEventHandle->filterKeys(filterIndex);
        const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle->getObjects();
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          ObjectsOfTriggerPath_.push_back(obj);
        }
      }
    }
    
    //Accessing electron collection
    edm::Handle<reco::GsfElectronCollection> gsfEleColl;
    iEvent.getByLabel(gsfelectron_, gsfEleColl);
    
    //Accessing conversion collection
    edm::Handle<reco::ConversionCollection> ConversionColl;
    iEvent.getByLabel(ConvCollLabel_, ConversionColl);

    //Accessing PFCandidate collection
    edm::Handle<reco::PFCandidateCollection> pfCandidatesH;
    iEvent.getByLabel(pfCandidateLabel_, pfCandidatesH);
    const reco::PFCandidateCollection thePfColl = *(pfCandidatesH.product());

    //Accessing beam spot handle
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel(BSLabel_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();

    //Accessing vertex collection
    edm::Handle<reco::VertexCollection> vtxColl;
    iEvent.getByLabel(primaryVertexLabel_, vtxColl);

    //getting rho handle for isolation
    edm::Handle<double> rhoHandle25;
    iEvent.getByLabel(rhoLabel25_, rhoHandle25);
    rho25=0.;
    if(rhoHandle25.isValid()) {
      rho25= *(rhoHandle25.product());
    }

  

    //Getting Sorted electron collection
    std::vector<reco::GsfElectron> gsfEle_container;
    gsfEle_container.clear(); 
    for(reco::GsfElectronCollection::const_iterator ele_itr = gsfEleColl->begin(); ele_itr != gsfEleColl->end(); ele_itr++){
    gsfEle_container.push_back(*ele_itr);
    }
    if(gsfEle_container.size() > 1){
      std::sort(gsfEle_container.begin(), gsfEle_container.end(), PtSortCriteriumElectron());
    }

    Electron_n = gsfEle_container.size();

    for(unsigned int ele = 0; ele < gsfEle_container.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(gsfEle_container), ele);

      //Kintematic Varables
      Electron_nocut_E.push_back(ele_ptr->energy());
      Electron_nocut_et.push_back(ele_ptr->et());
      Electron_nocut_pt.push_back(ele_ptr->pt());
      Electron_nocut_eta.push_back(ele_ptr->eta());
      Electron_nocut_phi.push_back(correct_Phi(ele_ptr->phi()));
      Electron_nocut_theta.push_back(ele_ptr->theta());
      Electron_nocut_px.push_back(ele_ptr->px());
      Electron_nocut_py.push_back(ele_ptr->py());
      Electron_nocut_pz.push_back(ele_ptr->pz());
      Electron_nocut_charge.push_back(ele_ptr->charge());

      //Id variables
      Electron_nocut_sigmaIEtaIEta.push_back(ele_ptr->sigmaIetaIeta());
      Electron_nocut_dEtaIn.push_back(ele_ptr->deltaEtaSuperClusterTrackAtVtx());
      Electron_nocut_dPhiIn.push_back(ele_ptr->deltaPhiSuperClusterTrackAtVtx());
      Electron_nocut_HoE.push_back(ele_ptr->hadronicOverEm());
      Electron_nocut_ooemoop.push_back(1.0/ele_ptr->ecalEnergy() - ele_ptr->eSuperClusterOverP()/ele_ptr->ecalEnergy());
      Electron_nocut_eopIn.push_back(ele_ptr->eSuperClusterOverP());
      Electron_nocut_fBrem.push_back(ele_ptr->fbrem());

      //Impact Parameter Variables
      float d0vtx         = 0.0;
      float dzvtx         = 0.0;
      if (vtxColl->size() > 0) {
	reco::VertexRef vtx(vtxColl, 0);
        d0vtx = ele_ptr->gsfTrack()->dxy(vtx->position());
        dzvtx = ele_ptr->gsfTrack()->dz(vtx->position());
      } else {
        d0vtx = ele_ptr->gsfTrack()->dxy();
        dzvtx = ele_ptr->gsfTrack()->dz();
      }
      Electron_nocut_d0vtx.push_back(d0vtx);
      Electron_nocut_dzvtx.push_back(dzvtx);

      //Conversion Variables
      Electron_nocut_vtxFitConversion.push_back(ConversionTools::hasMatchedConversion(*(ele_ptr), ConversionColl, beamspot.position()));
      Electron_nocut_misHits.push_back(ele_ptr->gsfTrack()->trackerExpectedHitsInner().numberOfHits());

      //Detector Isolation Variables
      Electron_nocut_trackIso03.push_back(ele_ptr->dr03TkSumPt());
      Electron_nocut_ecalIso03.push_back(ele_ptr->dr03EcalRecHitSumEt());
      Electron_nocut_hcalIso03.push_back(ele_ptr->dr03HcalTowerSumEt());

      //PFIsolation Variables
      unsigned int ivtx = 0;
      reco::VertexRef myVtxRef(vtxColl, ivtx);

      const reco::GsfElectron* electron = &(gsfEle_container.at(ele));

      eleIsolator03.fGetIsolation(electron, &thePfColl, myVtxRef, vtxColl);

      float PFIsoCh03 = eleIsolator03.getIsolationCharged();
      float PFIsoNh03 = eleIsolator03.getIsolationNeutral();
      float PFIsoEM03 = eleIsolator03.getIsolationPhoton();

      Electron_nocut_PFIsoCh03.push_back(PFIsoCh03);
      Electron_nocut_PFIsoNh03.push_back(PFIsoNh03);
      Electron_nocut_PFIsoEM03.push_back(PFIsoEM03);

      float Eff_area = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, 
                                                                       ele_ptr->superCluster()->eta(), EAtarget);
      double rhoPrime = std::max(rho25, 0.0);
      float PFIsoEM_Nh_03 = std::max((PFIsoNh03 + PFIsoEM03 - rhoPrime * Eff_area), 0.0);
      float PFIso_total = PFIsoCh03 + PFIsoEM_Nh_03;

      Electron_nocut_PFIsoSum_rhoCorr.push_back(PFIso_total); 

    } //END_OF for(unsigned int ele = 0; ele < gsfEle_container.size(); ele++)


    std::vector<reco::GsfElectron> Positron_tightId_tagContainer;
    Positron_tightId_tagContainer.clear();

    std::vector<reco::GsfElectron> Electron_looseId_probeContainer;
    Electron_looseId_probeContainer.clear();

    for(unsigned int i = 0; i < gsfEle_container.size(); i++){
      if(gsfEle_container[i].pdgId() == -11){ //For Positrons
	if(gsfEle_container[i].isEB()){
	  if(fabs(Electron_nocut_dEtaIn[i])                          < 0.004 &&
             fabs(Electron_nocut_dPhiIn[i])                          < 0.03  &&
             Electron_nocut_sigmaIEtaIEta[i]                         < 0.01  &&
             Electron_nocut_HoE[i]                                   < 0.12  &&
             fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
             fabs(Electron_nocut_dzvtx[i])                           < 0.1   &&
             fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
             Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.10  &&
             !(Electron_nocut_vtxFitConversion[i])                           &&
             Electron_nocut_misHits[i]                               <= 0
	     ){
	    Positron_tightId_tagContainer.push_back(gsfEle_container[i]);
	  }
	}
	else{ //EE
	  if(Electron_nocut_pt[i] > 20.0){
	    if(fabs(Electron_nocut_dEtaIn[i])                          < 0.005 &&
               fabs(Electron_nocut_dPhiIn[i])                          < 0.02  &&
               Electron_nocut_sigmaIEtaIEta[i]                         < 0.03  &&
               Electron_nocut_HoE[i]                                   < 0.10  &&
               fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
               fabs(Electron_nocut_dzvtx[i])                           < 0.1   &&
               fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
               Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.10  &&
               !(Electron_nocut_vtxFitConversion[i])                           &&
               Electron_nocut_misHits[i]                               <= 0
	       ){
	      Positron_tightId_tagContainer.push_back(gsfEle_container[i]);
	    }
	  }
	  else if(Electron_nocut_pt[i] < 20.0){
	    if(fabs(Electron_nocut_dEtaIn[i])                          < 0.005 &&
               fabs(Electron_nocut_dPhiIn[i])                          < 0.02  &&
               Electron_nocut_sigmaIEtaIEta[i]                         < 0.03  &&
               Electron_nocut_HoE[i]                                   < 0.10  &&
               fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
               fabs(Electron_nocut_dzvtx[i])                           < 0.1   &&
               fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
               Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.07  &&
               !(Electron_nocut_vtxFitConversion[i])                           &&
               Electron_nocut_misHits[i]                               <= 0
	       ){
	      Positron_tightId_tagContainer.push_back(gsfEle_container[i]);
	    }
	  }
	} //END_OF else{
      } //END_OF if(gsfEle_container[i].pdgId() == -11)
      else if(gsfEle_container[i].pdgId() == 11){ //For Electrons
	if(gsfEle_container[i].isEB()){
	  if(fabs(Electron_nocut_dEtaIn[i])                          < 0.007 &&
             fabs(Electron_nocut_dPhiIn[i])                          < 0.15  &&
             Electron_nocut_sigmaIEtaIEta[i]                         < 0.01  &&
             Electron_nocut_HoE[i]                                   < 0.12  &&
             fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
             fabs(Electron_nocut_dzvtx[i])                           < 0.2   &&
             fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
             Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.15  &&
             !(Electron_nocut_vtxFitConversion[i])                           &&
             Electron_nocut_misHits[i]                               <= 1
	     ){
	    Electron_looseId_probeContainer.push_back(gsfEle_container[i]);
	  }
	}
	else{ //EE
	  if(Electron_nocut_pt[i] > 20.0){
	    if(fabs(Electron_nocut_dEtaIn[i])                          < 0.009 &&
               fabs(Electron_nocut_dPhiIn[i])                          < 0.10  &&
               Electron_nocut_sigmaIEtaIEta[i]                         < 0.03  &&
               Electron_nocut_HoE[i]                                   < 0.10  &&
               fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
               fabs(Electron_nocut_dzvtx[i])                           < 0.2   &&
               fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
               Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.15  &&
               !(Electron_nocut_vtxFitConversion[i])                           &&
               Electron_nocut_misHits[i]                               <= 1
	       ){
	      Electron_looseId_probeContainer.push_back(gsfEle_container[i]);
	    }
	  }
	  else if(Electron_nocut_pt[i] < 20.0){
	    if(fabs(Electron_nocut_dEtaIn[i])                          < 0.009 &&
               fabs(Electron_nocut_dPhiIn[i])                          < 0.10  &&
               Electron_nocut_sigmaIEtaIEta[i]                         < 0.03  &&
               Electron_nocut_HoE[i]                                   < 0.10  &&
               fabs(Electron_nocut_d0vtx[i])                           < 0.02  &&
               fabs(Electron_nocut_dzvtx[i])                           < 0.2   &&
               fabs(Electron_nocut_ooemoop[i])                         < 0.05  &&
               Electron_nocut_PFIsoSum_rhoCorr[i]/Electron_nocut_pt[i] < 0.10  &&
               !(Electron_nocut_vtxFitConversion[i])                           &&
               Electron_nocut_misHits[i]                               <= 1
               ){
              Electron_looseId_probeContainer.push_back(gsfEle_container[i]);
            }
          }          
	} //END_OF else{
      } //END_OF else if(gsfEle_container[i].pdgId() == 11)
    } //END_OF for(unsigned int i = 0; i < gsfEle_container.size(); i++
    

    std::vector<reco::GsfElectron> Electron_allProbes_container;
    Electron_allProbes_container.clear();

    double Invarient_m = 0.0;
    for(unsigned int x = 0; x < Positron_tightId_tagContainer.size(); x++){
      for(unsigned int y = 0; y < Electron_looseId_probeContainer.size(); y++){
        double E = Positron_tightId_tagContainer[x].energy() + Electron_looseId_probeContainer[y].energy();
        double px = Positron_tightId_tagContainer[x].px() + Electron_looseId_probeContainer[y].px();
        double py = Positron_tightId_tagContainer[x].py() + Electron_looseId_probeContainer[y].py();
	double pz = Positron_tightId_tagContainer[x].pz() + Electron_looseId_probeContainer[y].pz();
	Invarient_m = sqrt(E*E - px*px - py*py - pz*pz);
	if(Invarient_m >= 70 && Invarient_m <= 120){
	  Electron_allProbes_container.push_back(Electron_looseId_probeContainer[y]);
	  break;
	}
      }
    }
    
    for(unsigned int ele = 0; ele < Electron_allProbes_container.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(Electron_allProbes_container), ele);

      //Kintematic Varables (with cuts and mass matching)
      Electron_allProbes_E.push_back(ele_ptr->energy());
      Electron_allProbes_et.push_back(ele_ptr->et());
      Electron_allProbes_pt.push_back(ele_ptr->pt());
      Electron_allProbes_eta.push_back(ele_ptr->eta());
      Electron_allProbes_phi.push_back(correct_Phi(ele_ptr->phi()));
      Electron_allProbes_theta.push_back(ele_ptr->theta());
      Electron_allProbes_px.push_back(ele_ptr->px());
      Electron_allProbes_py.push_back(ele_ptr->py());
      Electron_allProbes_pz.push_back(ele_ptr->pz());
      Electron_allProbes_charge.push_back(ele_ptr->charge());

      //Id variables
      Electron_allProbes_sigmaIEtaIEta.push_back(ele_ptr->sigmaIetaIeta());
      Electron_allProbes_dEtaIn.push_back(ele_ptr->deltaEtaSuperClusterTrackAtVtx());
      Electron_allProbes_dPhiIn.push_back(ele_ptr->deltaPhiSuperClusterTrackAtVtx());
      Electron_allProbes_HoE.push_back(ele_ptr->hadronicOverEm());
      Electron_allProbes_ooemoop.push_back(1.0/ele_ptr->ecalEnergy() - ele_ptr->eSuperClusterOverP()/ele_ptr->ecalEnergy());
      Electron_allProbes_eopIn.push_back(ele_ptr->eSuperClusterOverP());
      Electron_allProbes_fBrem.push_back(ele_ptr->fbrem());

      //Impact Parameter Variables
      float d0vtx         = 0.0;
      float dzvtx         = 0.0;
      if (vtxColl->size() > 0) {
	reco::VertexRef vtx(vtxColl, 0);
        d0vtx = ele_ptr->gsfTrack()->dxy(vtx->position());
        dzvtx = ele_ptr->gsfTrack()->dz(vtx->position());
      } else {
        d0vtx = ele_ptr->gsfTrack()->dxy();
        dzvtx = ele_ptr->gsfTrack()->dz();
      }
      Electron_allProbes_d0vtx.push_back(d0vtx);
      Electron_allProbes_dzvtx.push_back(dzvtx);
    }


    for(unsigned int ele = 0; ele < Electron_allProbes_container.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(Electron_allProbes_container), ele);
      
      if(ele_ptr->pt() > pTOfflineProbeCut_){
      h_allProbes->Fill(ele_ptr->pt());

      //This loop to get deltaR distribution only
      for(unsigned int x = 0; x < ObjectsOfTriggerPath_.size(); x++){
        const trigger::TriggerObject* trig_ptr = &(ObjectsOfTriggerPath_.at(x));
	double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	h_deltaR->Fill(deltaR);
      }

      const trigger::TriggerObject* trigPtr_min = &(ObjectsOfTriggerPath_.at(0));
      double min_dR = getDeltaR(ele_ptr->eta(), trigPtr_min->eta(), ele_ptr->phi(), trigPtr_min->phi());
    
      for(unsigned int x = 1; x < ObjectsOfTriggerPath_.size(); x++){
	const trigger::TriggerObject* trig_ptr = &(ObjectsOfTriggerPath_.at(x));
	double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	if(deltaR < min_dR){
	  min_dR = deltaR;
	}
      }

      if(min_dR < 0.1){
 	  h_passProbes->Fill(ele_ptr->pt());
	}
      } //END_OF if(ele_ptr->pt() > pTOfflineProbeCut_)
    } //END_OF for(unsigned int ele = 0; ele < Electron_allProbes_container.size(); ele++)


    //Getting efficiency plot using Divide method of TH1.
    h_efficiency->Divide(h_passProbes, h_allProbes, 1., 1., "");
    
  }
  

  

  t->Fill();

}

// ------------ method called when ending the processing of a run  ------------
void 
Trig_analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Trig_analyzer::endJob() 
{
  f->cd();
  h_allProbes->Write();
  h_deltaR->Write();
  h_passProbes->Write();
  h_efficiency->Write();

  f->WriteTObject(t);
  delete t;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Trig_analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Trig_analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Trig_analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//
// Destructors
//

Trig_analyzer::~Trig_analyzer()
{
   f->Close();
   delete f;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}



//define this as a plug-in
DEFINE_FWK_MODULE(Trig_analyzer);

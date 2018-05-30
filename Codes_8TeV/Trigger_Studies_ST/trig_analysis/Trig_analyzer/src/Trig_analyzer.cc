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
  OutputFile_(iConfig.getUntrackedParameter<std::string>("OutputFile", "")),
  tnp_trigger_(iConfig.getUntrackedParameter<std::string>("tnp_trigger", "")),
  tnp_leadinglegfilter_(iConfig.getUntrackedParameter<std::string>("tnp_leadinglegfilter", "")),
  tnp_trailinglegfilter_(iConfig.getUntrackedParameter<std::string>("tnp_trailinglegfilter", "")),
  Ele27_WP80_(iConfig.getUntrackedParameter<std::string>("Ele27_WP80", "")),
  Ele27_WP80_final_filter_(iConfig.getUntrackedParameter<std::string>("Ele27_WP80_final_filter", "")),
  gsfelectron_(iConfig.getUntrackedParameter<edm::InputTag>("gsfelectronTag")),
  triggerEvent_(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
  triggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag")),

  pTOfflineProbeCut_(iConfig.getUntrackedParameter<double>("pTOfflineProbeCut", 0.0))
{
  nEvents = 0;
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
Trig_analyzer::beginJob()
{
  f = new TFile(OutputFile_.c_str(), "recreate");

  //Declearing Histograms
  h_deltaR_EB = new TH1F("h_deltaR_EB", "DeltaR between RECO and HLT Object in Barrel", 100, 0.0, 10.0);
  h_deltaR_EB->GetYaxis()->SetTitle("\# of Objects");        h_deltaR_EB->GetYaxis()->CenterTitle();
  h_deltaR_EB->GetXaxis()->SetTitle("#Delta R");             h_deltaR_EB->GetXaxis()->CenterTitle();
  h_deltaR_EB->Sumw2();

  h_deltaR_EE = new TH1F("h_deltaR_EE", "DeltaR between RECO and HLT Object in Endcap", 100, 0.0, 10.0);
  h_deltaR_EE->GetYaxis()->SetTitle("\# of Objects");        h_deltaR_EE->GetYaxis()->CenterTitle();
  h_deltaR_EE->GetXaxis()->SetTitle("#Delta R");             h_deltaR_EE->GetXaxis()->CenterTitle();
  h_deltaR_EE->Sumw2();

  //Here change only nbins and xbins[nbins] value if need different no. of bins and different final bin.
  //And change i*5 to i*x if required to change bin size from 5 to x.
  int nbins = 76;//31;
  float xbins[nbins+1];
  for(int i = 0; i < nbins; i++){
    xbins[i] = i*2;
  }
  xbins[nbins] = 200;
 
  h_allProbes_EB = new TH1F("h_allProbes_EB", "Pt distribution of all Probes in Barrel Region", nbins, xbins);
  h_allProbes_EB->GetYaxis()->SetTitle("Objects/2 GeV");        h_allProbes_EB->GetYaxis()->CenterTitle();
  h_allProbes_EB->GetXaxis()->SetTitle("P_{T}(GeV)");           h_allProbes_EB->GetXaxis()->CenterTitle();
  h_allProbes_EB->Sumw2();
 
  h_allProbes_EE = new TH1F("h_allProbes_EE", "Pt distribution of all Probes in Endcap Region", nbins, xbins);
  h_allProbes_EE->GetYaxis()->SetTitle("Objects/2 GeV");        h_allProbes_EE->GetYaxis()->CenterTitle();
  h_allProbes_EE->GetXaxis()->SetTitle("P_{T}(GeV)");           h_allProbes_EE->GetXaxis()->CenterTitle();
  h_allProbes_EE->Sumw2();

  h_passProbes_EB = new TH1F("h_passProbes_EB", "Pt distribution of pass Probes in Barrel Region", nbins, xbins);
  h_passProbes_EB->GetYaxis()->SetTitle("Objects/2 GeV");       h_passProbes_EB->GetYaxis()->CenterTitle();
  h_passProbes_EB->GetXaxis()->SetTitle("P_{T}(GeV)");          h_passProbes_EB->GetXaxis()->CenterTitle();
  h_passProbes_EB->Sumw2();

  h_passProbes_EE = new TH1F("h_passProbes_EE", "Pt distribution of pass Probes in Endcap Region", nbins, xbins);
  h_passProbes_EE->GetYaxis()->SetTitle("Objects/2 GeV");       h_passProbes_EE->GetYaxis()->CenterTitle();
  h_passProbes_EE->GetXaxis()->SetTitle("P_{T}(GeV)");          h_passProbes_EE->GetXaxis()->CenterTitle();
  h_passProbes_EE->Sumw2();

  h_efficiency_EB = new TH1F("h_efficiency_EB", "Efficiency vs. Pt of reco objects in Barrel Region", nbins, xbins);
  h_efficiency_EB->GetYaxis()->SetTitle("Efficiency/2 GeV");    h_efficiency_EB->GetYaxis()->CenterTitle();
  h_efficiency_EB->GetXaxis()->SetTitle("P_{T}(GeV)");          h_efficiency_EB->GetXaxis()->CenterTitle();
  h_efficiency_EB->Sumw2();

  h_efficiency_EE = new TH1F("h_efficiency_EE", "Efficiency vs. Pt of reco objects in Endcap Region", nbins, xbins);
  h_efficiency_EE->GetYaxis()->SetTitle("Efficiency/2 GeV");    h_efficiency_EE->GetYaxis()->CenterTitle();
  h_efficiency_EE->GetXaxis()->SetTitle("P_{T}(GeV)");          h_efficiency_EE->GetXaxis()->CenterTitle();
  h_efficiency_EE->Sumw2();

}


// ------------ method called when starting to processes a run  ------------
void 
Trig_analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if(hltConfig_.init(iRun,iSetup,triggerEvent_.process(),changed)){
    //        cout << "Table name = " << hltConfig_.tableName() << endl;
    if(changed){
    }
  }
}

// ------------ method called for each event  ------------
void
Trig_analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  nEvents++;
  //    cout << "PROCESSING EVENT = " << nEvents << endl;

  Electron_n = 0;  

  Clear();
  
  //-----------------------------------------------------
  //Accessing TriggerResults and triggerNames collections
  //-----------------------------------------------------
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  iEvent.getByLabel(triggerEvent_,triggerEventHandle_);
  if(!triggerEventHandle_.isValid()){
    cout << "HLT TriggerEvent Handle is not valid______Can't get the product" << endl;
  }
  else{
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResults_, triggerResults);

    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
    const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();

    //------------------------------------------
    //Getting the correct version of tnp trigger
    //------------------------------------------
    int tnp_trigger_idx = -1;
    for(unsigned int i=0; i<triggerNames_.size(); i++){
      if(triggerNames_[i].find(tnp_trigger_.c_str())!=std::string::npos){
	tnp_trigger_idx = i;
	break;
      }
    }
    tnp_trigger_ = triggerNames_[tnp_trigger_idx];

    //--------------------------------------------------------------
    //Getting list of all the savetags = true filters of tnp_trigger
    //--------------------------------------------------------------
    if(hltConfig_.triggerIndex(tnp_trigger_)<hltConfig_.size()){
      tnp_filters_ = hltConfig_.saveTagsModules(tnp_trigger_);
    }    

    //---------------------------------------------------------------------
    //Accessing trigger object collection passing the tnp_leadingleg_filter
    //---------------------------------------------------------------------
    std::vector<trigger::TriggerObject> Ele17_triggerObjs_;
    Ele17_triggerObjs_.clear();

    trigger::size_type leadingleg_filterIdx = triggerEventHandle_->filterIndex(edm::InputTag(tnp_leadinglegfilter_, "", triggerEvent_.process()));
    if(leadingleg_filterIdx < triggerEventHandle_->sizeFilters()){

      const trigger::Keys& trigKeys = triggerEventHandle_->filterKeys(leadingleg_filterIdx);
      const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle_->getObjects();
      for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	Ele17_triggerObjs_.push_back(obj);
      }
    }

    //   if(Ele17_triggerObjs_.size() > 0)  cout << "Ele17_triggerObjs_ size is  = " << Ele17_triggerObjs_.size() << endl;

    //----------------------------------------------------------------------
    //Accessing trigger object collection passing the tnp_trailingleg_filter
    //----------------------------------------------------------------------
    std::vector<trigger::TriggerObject> Ele8_triggerObjs_;
    Ele8_triggerObjs_.clear();

    trigger::size_type trailingleg_filterIdx = triggerEventHandle_->filterIndex(edm::InputTag(tnp_trailinglegfilter_, "", triggerEvent_.process()));
    if(trailingleg_filterIdx < triggerEventHandle_->sizeFilters()){

      const trigger::Keys& trigKeys = triggerEventHandle_->filterKeys(trailingleg_filterIdx);
      const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle_->getObjects();
      for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
        const trigger::TriggerObject& obj = trigObjColl[*keyIt];
        Ele8_triggerObjs_.push_back(obj);
      }
    }

    //      if(Ele8_triggerObjs_.size() > 0) cout << "Ele8_triggerObjs size is = " << Ele8_triggerObjs_.size() << endl;
 
    //-------------------------------------------------
    //Getting the correct version of Ele27_WP80 trigger
    //-------------------------------------------------
    int Ele27_WP80_idx = -1;
    for(unsigned int i=0; i<triggerNames_.size(); i++){
      if(triggerNames_[i].find(Ele27_WP80_.c_str())!=std::string::npos){
        Ele27_WP80_idx = i;
        break;
      }
    }
    Ele27_WP80_ = triggerNames_[Ele27_WP80_idx];

    //-----------------------------------------------------------------
    //Getting list of all saveTags = true filters of Ele27_WP80 trigger
    //-----------------------------------------------------------------
    if(hltConfig_.triggerIndex(Ele27_WP80_)<hltConfig_.size()){
      filtersOf_Ele27_WP80_ = hltConfig_.saveTagsModules(Ele27_WP80_);
    }

    //------------------------------------------------------------------------
    //Accessing trigger object collection passing filter of Ele27_WP80 trigger
    //------------------------------------------------------------------------
    std::vector<trigger::TriggerObject> Ele27_WP80_triggerObjs_;
    Ele27_WP80_triggerObjs_.clear();

    trigger::size_type Ele27_WP80_filterIdx = triggerEventHandle_->filterIndex(edm::InputTag(Ele27_WP80_final_filter_, "", triggerEvent_.process()));
    if(Ele27_WP80_filterIdx < triggerEventHandle_->sizeFilters()){

	const trigger::Keys& trigKeys = triggerEventHandle_->filterKeys(Ele27_WP80_filterIdx);
	const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle_->getObjects();
	for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
	  const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	  Ele27_WP80_triggerObjs_.push_back(obj);
	}
    }

    //        if(Ele27_WP80_triggerObjs_.size() > 0) cout << "Ele27_WP80_triggerObjs_ size  is = " << Ele27_WP80_triggerObjs_.size() << endl;
 
    //Accessing electron collection
    edm::Handle<reco::GsfElectronCollection> gsfEleColl;
    iEvent.getByLabel(gsfelectron_, gsfEleColl);

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
    
    //---------------------------------------------------------------------------------------------
    //Getting tag collection of gsfElectrons matching with leading leg Ele17 of tnp and pt > 27 GeV 
    //---------------------------------------------------------------------------------------------
    std::vector<reco::GsfElectron> Tag_container_; 
    Tag_container_.clear();

    std::vector<reco::GsfElectron> restOfgsfEle_container_;
    restOfgsfEle_container_.clear();

    for(unsigned int ele = 0; ele < gsfEle_container.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(gsfEle_container), ele);

      double ele17min_dR = 99.0;
     
      if(Ele17_triggerObjs_.size() > 0){
	const trigger::TriggerObject* ele17ptr_min = &(Ele17_triggerObjs_.at(0));
	ele17min_dR = getDeltaR(ele_ptr->eta(), ele17ptr_min->eta(), ele_ptr->phi(), ele17ptr_min->phi());
     
	if(Ele17_triggerObjs_.size() > 1){
	  for(unsigned int x = 1; x < Ele17_triggerObjs_.size(); x++){
	    const trigger::TriggerObject* ele17_ptr = &(Ele17_triggerObjs_.at(x));
	    double deltaR = getDeltaR(ele_ptr->eta(), ele17_ptr->eta(), ele_ptr->phi(), ele17_ptr->phi());
	    if(deltaR < ele17min_dR){
	      ele17min_dR = deltaR;
	    }
	  }
	}
      }
    
      if(ele17min_dR < 0.1 && ele_ptr->pt() > 27){
	Tag_container_.push_back(*ele_ptr);
      }else{
	restOfgsfEle_container_.push_back(*ele_ptr);
      }
    }
   
    //        if(Tag_container_.size() > 0) cout << "Tag container size = " << Tag_container_.size() << endl;

    //---------------------------------------------------------------------------------------------------------
    //Getting probe collection of gsfElectrons matching with trailing leg of tnp from the rest of gsf electrons
    //--------------------------------------------------------------------------------------------------------- 
    std::vector<reco::GsfElectron> NoIdProbe_container_; //NoId since trailing electron in tnp is chosen without any id cut
    NoIdProbe_container_.clear();

    for(unsigned int ele = 0; ele < restOfgsfEle_container_.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(restOfgsfEle_container_), ele);

      double ele8min_dR = 99.0;

      if(Ele8_triggerObjs_.size() > 0){
	const trigger::TriggerObject* ele8ptr_min = &(Ele8_triggerObjs_.at(0));
	ele8min_dR = getDeltaR(ele_ptr->eta(), ele8ptr_min->eta(), ele_ptr->phi(), ele8ptr_min->phi());
      
	if(Ele8_triggerObjs_.size() > 1){
	  for(unsigned int x = 1; x < Ele8_triggerObjs_.size(); x++){
	    const trigger::TriggerObject* ele8_ptr = &(Ele8_triggerObjs_.at(x));
	    double deltaR = getDeltaR(ele_ptr->eta(), ele8_ptr->eta(), ele_ptr->phi(), ele8_ptr->phi());
	    if(deltaR < ele8min_dR){
	      ele8min_dR = deltaR;
	    }
	  }
	}
      }
      
      if(ele8min_dR < 0.1){
	NoIdProbe_container_.push_back(*ele_ptr);
      }
    }

    //        if(NoIdProbe_container_.size() > 0) cout << "NoIdProbe container size =" << NoIdProbe_container_.size() << endl;

    std::vector<reco::GsfElectron> NoMatchingProbe_container_;
    NoMatchingProbe_container_.clear();

    for(unsigned int x = 0; x < NoIdProbe_container_.size(); x++){
      const reco::GsfElectronRef probe_ptr(&(NoIdProbe_container_), x);

      float sigmaIetaIeta = probe_ptr->sigmaIetaIeta();
      float dEtaIn = probe_ptr->deltaEtaSuperClusterTrackAtVtx();
      float dPhiIn = probe_ptr->deltaPhiSuperClusterTrackAtVtx();
      float HoE = probe_ptr->hadronicOverEm();
      float trackIso03 = probe_ptr->dr03TkSumPt();
      float ecalIso03 = probe_ptr->dr03EcalRecHitSumEt();
      float hcalIso03 = probe_ptr->dr03HcalTowerSumEt();

      if(probe_ptr->isEB()){
	if(sigmaIetaIeta                  <   0.01      &&
           fabs(dEtaIn)                   <   0.007     &&
           fabs(dPhiIn)                   <   0.15      &&
           HoE                            <   0.12      &&
           trackIso03/probe_ptr->pt()     <   0.2       &&
           ecalIso03/probe_ptr->pt()      <   0.2       &&
           hcalIso03/probe_ptr->pt()      <   0.2      
	   ){
	  NoMatchingProbe_container_.push_back(*probe_ptr);
	}
      }else{
	if(sigmaIetaIeta                  <   0.03      &&
           fabs(dEtaIn)                   <   0.009     &&
           fabs(dPhiIn)                   <   0.10      &&
           HoE                            <   0.10      &&
           trackIso03/probe_ptr->pt()     <   0.2       &&
           ecalIso03/probe_ptr->pt()      <   0.2       &&
           hcalIso03/probe_ptr->pt()      <   0.2
	   ){
	  NoMatchingProbe_container_.push_back(*probe_ptr);
	}
      }
    }

    //        if(NoMatchingProbe_container_.size() > 0) cout << "NoMatchingProbe container size =" << NoMatchingProbe_container_.size() << endl;

    std::vector<reco::GsfElectron> AllProbes_container_;
    AllProbes_container_.clear();

    double Invarient_m = 0.0;
    for(unsigned int x = 0; x < NoMatchingProbe_container_.size(); x++){
      for(unsigned int y = 0; y < Tag_container_.size(); y++){
	
        double E = NoMatchingProbe_container_[x].energy() + Tag_container_[y].energy();
        double px = NoMatchingProbe_container_[x].px() + Tag_container_[y].px();
        double py = NoMatchingProbe_container_[x].py() + Tag_container_[y].py();
	double pz = NoMatchingProbe_container_[x].pz() + Tag_container_[y].pz();

	Invarient_m = sqrt(E*E - px*px - py*py - pz*pz);

	if(Invarient_m >= 70 && Invarient_m <= 120){
	  AllProbes_container_.push_back(NoMatchingProbe_container_[x]);
	  break;
	}
      }
    }

    //    if(AllProbes_container_.size() > 0) cout << "All Probes container size =" << AllProbes_container_.size() << endl;

    //    int count_all_eb = 0;
    //    int count_all_ee = 0;
    //    int count_pass_eb = 0;
    //    int count_pass_ee = 0;

    for(unsigned int ele = 0; ele < AllProbes_container_.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(AllProbes_container_), ele);

      if(ele_ptr->isEB()){
        h_allProbes_EB->Fill(ele_ptr->pt());
	if(ele_ptr->pt() > 27){
	  h_passProbes_EB->Fill(ele_ptr->pt());
	}
      }

      if(ele_ptr->isEE()){
        h_allProbes_EE->Fill(ele_ptr->pt());
        if(ele_ptr->pt() > 27){
          h_passProbes_EE->Fill(ele_ptr->pt());
	}
      }
    }

    /*
    for(unsigned int ele = 0; ele < AllProbes_container_.size(); ele++){
      const reco::GsfElectronRef ele_ptr(&(AllProbes_container_), ele);
      
      if(ele_ptr->isEB()){
	h_allProbes_EB->Fill(ele_ptr->pt());
	//	count_all_eb++;

	//This loop to get deltaR distribution only
	if(Ele27_WP80_triggerObjs_.size() > 0){
	  for(unsigned int x = 0; x < Ele27_WP80_triggerObjs_.size(); x++){
	    const trigger::TriggerObject* trig_ptr = &(Ele27_WP80_triggerObjs_.at(x));
	    double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	    h_deltaR_EB->Fill(deltaR);
	  }
	}

	double ele27_wp80min_dR = 99.0;

	if(Ele27_WP80_triggerObjs_.size() > 0){
	  const trigger::TriggerObject* ele27_wp80ptr_min = &(Ele27_WP80_triggerObjs_.at(0));
	  ele27_wp80min_dR = getDeltaR(ele_ptr->eta(), ele27_wp80ptr_min->eta(), ele_ptr->phi(), ele27_wp80ptr_min->phi());
    
	  if(Ele27_WP80_triggerObjs_.size() > 1){
	    for(unsigned int x = 1; x < Ele27_WP80_triggerObjs_.size(); x++){
	      const trigger::TriggerObject* ele27_wp80ptr = &(Ele27_WP80_triggerObjs_.at(x));
	      double deltaR = getDeltaR(ele_ptr->eta(), ele27_wp80ptr->eta(), ele_ptr->phi(), ele27_wp80ptr->phi());
	      if(deltaR < ele27_wp80min_dR){
		ele27_wp80min_dR = deltaR;
	      }
	    }
	  }
	}

	if(ele27_wp80min_dR < 0.1){
	  //	  count_pass_eb++;
	  h_passProbes_EB->Fill(ele_ptr->pt());
	}
      }

      if(ele_ptr->isEE()){
	h_allProbes_EE->Fill(ele_ptr->pt());
	//	count_all_ee++;

	//This loop to get deltaR distribution only
	if(Ele27_WP80_triggerObjs_.size() > 0){
	  for(unsigned int x = 0; x < Ele27_WP80_triggerObjs_.size(); x++){
	    const trigger::TriggerObject* trig_ptr = &(Ele27_WP80_triggerObjs_.at(x));
	    double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	    h_deltaR_EE->Fill(deltaR);
	  }
	}

	double ele27_wp80min_dR = 99.0;

	if(Ele27_WP80_triggerObjs_.size() > 0){
	  const trigger::TriggerObject* ele27_wp80ptr_min = &(Ele27_WP80_triggerObjs_.at(0));
	  ele27_wp80min_dR = getDeltaR(ele_ptr->eta(), ele27_wp80ptr_min->eta(), ele_ptr->phi(), ele27_wp80ptr_min->phi());
    
	  if(Ele27_WP80_triggerObjs_.size() > 1){
	    for(unsigned int x = 1; x < Ele27_WP80_triggerObjs_.size(); x++){
	      const trigger::TriggerObject* ele27_wp80ptr = &(Ele27_WP80_triggerObjs_.at(x));
	      double deltaR = getDeltaR(ele_ptr->eta(), ele27_wp80ptr->eta(), ele_ptr->phi(), ele27_wp80ptr->phi());
	      if(deltaR < ele27_wp80min_dR){
		ele27_wp80min_dR = deltaR;
	      }
	    }
	  }
	}

	if(ele27_wp80min_dR < 0.1){
	  //	  count_pass_ee++;
	  h_passProbes_EE->Fill(ele_ptr->pt());
	}
      }
    }
    */

    //  if(count_all_eb > 0) cout << "All probe container in barrel size = " << count_all_eb << endl;
    //  if(count_all_ee > 0) cout << "All probe container in endcap size = " << count_all_ee << endl;
    //  if(count_pass_eb > 0) cout << "Pass probe container in barrel size = " << count_pass_eb << endl;
    //  if(count_pass_ee > 0) cout << "Pass probe container in endcap size = " << count_pass_ee << endl;
 
  
  // Getting efficiency plot using Divide method of TH1.
  h_efficiency_EB->Divide(h_passProbes_EB, h_allProbes_EB, 1., 1., "");
  h_efficiency_EE->Divide(h_passProbes_EE, h_allProbes_EE, 1., 1., "");
  
  }
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
  h_allProbes_EB->Write();
  h_deltaR_EB->Write();
  h_passProbes_EB->Write();
  h_efficiency_EB->Write();
  h_allProbes_EE->Write();
  h_deltaR_EE->Write();
  h_passProbes_EE->Write();
  h_efficiency_EE->Write();
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

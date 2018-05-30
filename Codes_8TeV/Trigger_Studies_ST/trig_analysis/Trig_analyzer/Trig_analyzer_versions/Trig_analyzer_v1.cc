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

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// constants, enums and typedefs 
//

//
// static data member definitions
//


//
// constructors 
//

Trig_analyzer::Trig_analyzer(const edm::ParameterSet& iConfig)
{
  tightPath_            = iConfig.getUntrackedParameter<std::string>("tightPath", "");
  gsfelectron_          = iConfig.getUntrackedParameter<edm::InputTag>("gsfelectronTag");
  triggerEvent_         = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag");
  triggerResults_       = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag");
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
Trig_analyzer::beginJob()
{
  f = new TFile("Trig_analyzer.root", "recreate");
  
  //Declearing Histograms
  h_deltaR = new TH1F("h_deltaR", "DeltaR between RECO and HLT Object", 100, 0.0, 10.0);
  h_deltaR->GetYaxis()->SetTitle("\# of Objects");        h_deltaR->GetYaxis()->CenterTitle();
  h_deltaR->GetXaxis()->SetTitle("#Delta R");             h_deltaR->GetXaxis()->CenterTitle();
  //h_deltaR->Sumw2();

  h_noCut = new TH1F("h_noCut", "Pt distribution of reco objects with no trigger Cut", 100, 0.0, 500.0);
  h_noCut->GetYaxis()->SetTitle("Objects/5 GeV");        h_noCut->GetYaxis()->CenterTitle();
  h_noCut->GetXaxis()->SetTitle("P_{T}(GeV)");             h_noCut->GetXaxis()->CenterTitle();
  //h_noCut->Sumw2();

  h_tightPath = new TH1F("h_tightPath", "Pt distribution of reco objects with a tight trigger Cut", 100, 0.0, 500.0);
  h_tightPath->GetYaxis()->SetTitle("Objects/5 GeV");          h_tightPath->GetYaxis()->CenterTitle();
  h_tightPath->GetXaxis()->SetTitle("P_{T}(GeV)");             h_tightPath->GetXaxis()->CenterTitle();
  //h_tightPath->Sumw2();

  h_efficiency = new TH1F("h_efficiency", "Efficiency vs. Pt of reco objects", 100, 0.0, 500.0);
  h_efficiency->GetYaxis()->SetTitle("Efficiency");          h_efficiency->GetYaxis()->CenterTitle();
  h_efficiency->GetXaxis()->SetTitle("P_{T}(GeV)");             h_efficiency->GetXaxis()->CenterTitle();
  //h_efficiency->Sumw2();

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
  
  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  iEvent.getByLabel(triggerEvent_,triggerEventHandle);

  if (!triggerEventHandle.isValid()) {
  std::cout << "Error! Can't get the product "<< triggerEvent_.label() << std::endl;
  } else{
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResults_, triggerResults);

    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
    const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
   
    int tightPath_idx = -1;
    for(unsigned int i=0; i<triggerNames_.size(); i++){
      if(triggerNames_[i].find(tightPath_.c_str())!=std::string::npos){
	tightPath_idx = i;
	break;
      }
    }
    tightPath_ = triggerNames_[tightPath_idx];
   
    if(hltConfig_.triggerIndex(tightPath_)<hltConfig_.size()){
	filtersOfTightPath_ = hltConfig_.saveTagsModules(tightPath_);
    }
    
    std::vector<trigger::TriggerObject> ObjectsOfTightPath_;
    ObjectsOfTightPath_.clear();

    for(size_t filterNr=0; filterNr<filtersOfTightPath_.size();filterNr++){
      trigger::size_type filterIndex = triggerEventHandle->filterIndex(edm::InputTag(filtersOfTightPath_[filterNr], "", triggerEvent_.process()));
      if(filterIndex < triggerEventHandle->sizeFilters()){
        const trigger::Keys& trigKeys = triggerEventHandle->filterKeys(filterIndex);
        const trigger::TriggerObjectCollection & trigObjColl = triggerEventHandle->getObjects();
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); keyIt++){
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          ObjectsOfTightPath_.push_back(obj);
        }
      }
    }
    
    edm::Handle<reco::GsfElectronCollection> gsfEleColl;
    iEvent.getByLabel(gsfelectron_, gsfEleColl);
    reco::GsfElectronCollection::const_iterator ele_itr;
    
    for(ele_itr=gsfEleColl->begin(); ele_itr!=gsfEleColl->end(); ele_itr++){
      const reco::GsfElectron* ele_ptr = &(*ele_itr);

      h_noCut->Fill(ele_ptr->pt());

      //This loop to get deltaR distribution only
      for(unsigned int x = 0; x < ObjectsOfTightPath_.size(); x++){
        const trigger::TriggerObject* trig_ptr = &(ObjectsOfTightPath_.at(x));
	double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	h_deltaR->Fill(deltaR);
      }

      const trigger::TriggerObject* trigPtr_min = &(ObjectsOfTightPath_.at(0));
      double min_dR = getDeltaR(ele_ptr->eta(), trigPtr_min->eta(), ele_ptr->phi(), trigPtr_min->phi());
    
      for(unsigned int x = 1; x < ObjectsOfTightPath_.size(); x++){
	const trigger::TriggerObject* trig_ptr = &(ObjectsOfTightPath_.at(x));
	double deltaR = getDeltaR(ele_ptr->eta(), trig_ptr->eta(), ele_ptr->phi(), trig_ptr->phi());
	if(deltaR < min_dR){
	  min_dR = deltaR;
	}
      }

      if(min_dR < 0.1){
 	  h_tightPath->Fill(ele_ptr->pt());
	}
    }
    //Getting efficiency plot using Divide method of TH1.
    h_efficiency->Divide(h_tightPath, h_noCut, 1., 1., "");
    /*
    TGraphAsymmErrors *h_efficiency = new TGraphAsymmErrors();
    h_efficiency->SetName("h_efficiency");
    h_efficiency->BayesDivide(h_tightPath,h_noCut);
    h_efficiency->GetXaxis()->SetTitle("P_{T} (GeV)");
    h_efficiency->GetXaxis()->CenterTitle();
    h_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_efficiency->GetYaxis()->CenterTitle();
    */
 
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
  h_noCut->Write();
  h_deltaR->Write();
  h_tightPath->Write();
  h_efficiency->Write();
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

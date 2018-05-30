#include "ADDmonophoton/Skimmer/interface/MonoPhotonSkimmer.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include <vector>
#include <TMath.h>
#include <iostream>
using namespace std;

MonoPhotonSkimmer::MonoPhotonSkimmer(const edm::ParameterSet& iConfig){

    //now do what ever initialization is needed
    _phoTag = iConfig.getParameter<edm::InputTag>("phoTag");

    _hadoveremEB = iConfig.getParameter<double>("hadoveremEB");
    _scHighPtThreshEB = iConfig.getParameter<double>("scHighPtThreshEB");
    _photonEtaThreshEB =  iConfig.getParameter<double>("photonEtaThreshEB");

}


MonoPhotonSkimmer::~MonoPhotonSkimmer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool MonoPhotonSkimmer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<reco::PhotonCollection> photonColl;
    iEvent.getByLabel(_phoTag, photonColl);
    const reco::PhotonCollection *photons = photonColl.product();  

    //Iterate over photon collection.
    std::vector<reco::Photon> PreselPhotons;
    PreselPhotons.clear();

    bool decision=false;

    reco::PhotonCollection::const_iterator pho;
    for (pho = (*photons).begin(); pho!= (*photons).end(); pho++)
    {  
	if(pho->pt() >  _scHighPtThreshEB && fabs(pho->eta())< _photonEtaThreshEB && pho->hadTowOverEm() < _hadoveremEB )
	{ 
	    PreselPhotons.push_back(*pho);
	}

    }//Loop over pat::Photons

    //atleat 1 photon
    if(PreselPhotons.size() > 0 )decision=true;
    if(PreselPhotons.size() <= 0)decision=false;

    return decision;
}

// ------------ method called once each job just before starting event loop  ------------
void MonoPhotonSkimmer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MonoPhotonSkimmer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoPhotonSkimmer);

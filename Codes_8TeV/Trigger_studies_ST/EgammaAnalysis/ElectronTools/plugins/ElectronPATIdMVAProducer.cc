// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//
// class declaration
//

using namespace std;
using namespace reco;
class ElectronPATIdMVAProducer : public edm::EDProducer {
	public:
		explicit ElectronPATIdMVAProducer(const edm::ParameterSet&);
		virtual ~ElectronPATIdMVAProducer();

	private:
		virtual void produce(edm::Event&, const edm::EventSetup&);

		// ----------member data ---------------------------
  bool verbose_;
  edm::InputTag electronTag_;
  double _Rho;
  edm::InputTag eventrho;
  string method_;
  vector<string> mvaWeightFiles_;
  
  
  EGammaMvaEleEstimator* mvaID_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronPATIdMVAProducer::ElectronPATIdMVAProducer(const edm::ParameterSet& iConfig) {
        verbose_ = iConfig.getUntrackedParameter<bool>("verbose", false);
	electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
	method_ = iConfig.getParameter<string>("method");
	std::vector<string> fpMvaWeightFiles = iConfig.getParameter<std::vector<std::string> >("mvaWeightFile");
	eventrho = iConfig.getParameter<edm::InputTag>("Rho");

        produces<edm::ValueMap<float> >();

        mvaID_ = new EGammaMvaEleEstimator();
 
        EGammaMvaEleEstimator::MVAType type_;
        
	 
	type_ = EGammaMvaEleEstimator::kTrigNoIP;
	 
	

        bool manualCat_ = true;

	string path_mvaWeightFileEleID;
	for(unsigned ifile=0 ; ifile < fpMvaWeightFiles.size() ; ++ifile) {
	  path_mvaWeightFileEleID = edm::FileInPath ( fpMvaWeightFiles[ifile].c_str() ).fullPath();
	  mvaWeightFiles_.push_back(path_mvaWeightFileEleID);
	}
	
        mvaID_->initialize(method_, type_, manualCat_, mvaWeightFiles_);

}


ElectronPATIdMVAProducer::~ElectronPATIdMVAProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void ElectronPATIdMVAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

        std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>() );



      

       

	Handle<pat::ElectronCollection> egCollection;
	iEvent.getByLabel(electronTag_,egCollection);
        const pat::ElectronCollection egCandidates = (*egCollection.product());

	_Rho=0;
	edm::Handle<double> rhoPtr;
	//const edm::InputTag eventrho("kt6PFJets", "rho");
	iEvent.getByLabel(eventrho,rhoPtr);
	_Rho=*rhoPtr;

        std::vector<float> values;
        values.reserve(egCollection->size());
   
        for ( pat::ElectronCollection::const_iterator egIter = egCandidates.begin(); egIter != egCandidates.end(); ++egIter) {

          double mvaVal = -999999;
	  
	  
	  
	  
	  mvaVal = mvaID_->mvaValue( *egIter, _Rho, verbose_);
	  
	  
	  values.push_back( mvaVal ); 
	}

        edm::ValueMap<float>::Filler filler(*out);
        filler.insert(egCollection, values.begin(), values.end() );
	filler.fill();

	iEvent.put(out);

   
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronPATIdMVAProducer);

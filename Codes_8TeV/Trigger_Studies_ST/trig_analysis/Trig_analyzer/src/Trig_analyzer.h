// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

//Root Header Files
#include "TROOT.h"
#include "TObject.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TMath.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"


#include <cmath>
using namespace edm;
using namespace std;
using namespace ROOT;
using namespace trigger;

//
// class declaration
//

class Trig_analyzer : public edm::EDAnalyzer {
   public:
      explicit Trig_analyzer(const edm::ParameterSet&);
      ~Trig_analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void Clear();
      double getDeltaR(double, double, double, double);
      double correct_Phi(double phi);
      float correct_Phi(float phi);

      // ----------member data ---------------------------
      TFile                      *f;
      std::string                OutputFile_;


      // Defining Histograms for Turns-Ons
      TH1F                      *h_allProbes_EB, *h_passProbes_EB, *h_efficiency_EB;
      TH1F                      *h_deltaR_EB;
      TH1F                      *h_allProbes_EE, *h_passProbes_EE, *h_efficiency_EE;
      TH1F                      *h_deltaR_EE;


      HLTConfigProvider          hltConfig_;

      std::string                tnp_trigger_;
      std::string                tnp_leadinglegfilter_;
      std::string                tnp_trailinglegfilter_;
      vector<std::string>        tnp_filters_;
      std::string                Ele27_WP80_;
      std::string                Ele27_WP80_final_filter_;
      vector<std::string>        filtersOf_Ele27_WP80_;
 
      int                                       nEvents;

      //Electron Variables
      int                                       Electron_n;

      //Pt cut while getting efficiency
      double                                    pTOfflineProbeCut_;

      //InputTags
      edm::InputTag              gsfelectron_;
      edm::InputTag              triggerEvent_;
      edm::InputTag              triggerResults_; 
};

void Trig_analyzer::Clear(){
  tnp_filters_.clear();
  filtersOf_Ele27_WP80_.clear();
}

//Compute DeltaR
double Trig_analyzer::getDeltaR(double eta1, double eta2, double phi1, double phi2){
  double dEta = fabs(eta1 - eta2);
  double dPhi  = fabs(phi1 - phi2);
  double twopi = 2.0*(TMath::Pi());
  if(dPhi < 0) dPhi = - dPhi;
  if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;
  double dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

//Compute correct phi
double Trig_analyzer::correct_Phi(double phi){
  return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}

float Trig_analyzer::correct_Phi(float phi){
  return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}



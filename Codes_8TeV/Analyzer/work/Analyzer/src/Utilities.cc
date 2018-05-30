#include "TMath.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


//Changing phi<0 to (2*pi + phi) to make phi varying from 0 to 360 instead of -180 to +180. 
double correct_Phi(double phi){
  return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}

float correct_Phi(float phi){
  return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}



//Evaluating theta from eta using the formula eta = -log(tan(theta/2)) => theta =2*atan(exp(-eta)) where atan is tan inverse.
double Theta(double eta){
  double theta = 2. * atan(exp(-eta));
  return theta;
}

float Theta(float eta){
  float theta = 2. * atan(exp(-eta));
  return theta;
}



//Getting longitudinal momentum from total and transverse momentum.
double Pl(double P, double Pt){
  double pl = sqrt(pow(P,2)-pow(Pt,2));
  return pl;
}

float Pl(float P, float Pt){
  float pl = sqrt(pow(P,2)-pow(Pt,2));
  return pl;
}



//Functions to calculate the recHit Energy for the DetId given as an argument from the rechitCollection also given as argument.
//Calculating energy of given DetId. Taken from RecoEcal/EgammaCoreTools/src/EcalTools.cc
float recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}



//Calculting energy of DetId corresponding to di, dj. Taken from RecoEcal/EgammaCoreTools/src/EcalTools.cc
float recHitE( const DetId id, const EcalRecHitCollection & recHits, int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY

  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );

  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}



//Calculating Et using the formula E = Et cos(ieta) = Et cosh(eta). Taken from RecoEcal/EgammaCoreTools/src/EcalTools.cc
float recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits ){
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}



//Calculating energy for DetId if greater than a threshold and also using timing info. Taken from RecoLocalCalo/EcalRecAlgos/src/EcalCleaningAlgo.cc
float recHitE( const DetId id, const EcalRecHitCollection &recHits, bool useTimingInfo )
{
  //These values are taken from:RecoLocalCalo/EcalRecAlgos/python/ecalCleaningAlgo.py
  double e4e1Threshold_barrel_  = 0.080;
  double e4e1Threshold_endcap_  = 0.300;
  double ignoreOutOfTimeThresh_ = 1e9;

  if ( id.rawId() == 0 ) return 0;

  float threshold = e4e1Threshold_barrel_;
  if ( id.subdetId() == EcalEndcap) threshold = e4e1Threshold_endcap_;

  EcalRecHitCollection::const_iterator it = recHits.find( id );
  if ( it != recHits.end() ){
    float ene= (*it).energy();

    // ignore out of time in EB when making e4e1 if so configured
    if (useTimingInfo){

      if (id.subdetId()==EcalBarrel &&
          it->checkFlag(EcalRecHit::kOutOfTime)
          && ene>ignoreOutOfTimeThresh_) return 0;
    }

    // ignore hits below threshold
    if (ene < threshold) return 0;

    // else return the energy of this hit
    return ene;
  }
  return 0;
}



// four neighbours in the swiss cross around id. Taken from RecoLocalCalo/EcalRecAlgos/src/EcalCleaningAlgo.cc
const std::vector<DetId> neighbours(const DetId& id){

  std::vector<DetId> ret;

  if ( id.subdetId() == EcalBarrel) {

    ret.push_back( EBDetId::offsetBy( id,  1, 0 ));
    ret.push_back( EBDetId::offsetBy( id, -1, 0 ));
    ret.push_back( EBDetId::offsetBy( id,  0, 1 ));
    ret.push_back( EBDetId::offsetBy( id,  0,-1 ));
  }
  // nobody understands what polymorphism is for, sgrunt !
  else  if (id.subdetId() == EcalEndcap) {
    ret.push_back( EEDetId::offsetBy( id,  1, 0 ));
    ret.push_back( EEDetId::offsetBy( id, -1, 0 ));
    ret.push_back( EEDetId::offsetBy( id,  0, 1 ));
    ret.push_back( EEDetId::offsetBy( id,  0,-1 ));
  }
  return ret;
}




#include <iostream>
#include <fstream>
#include <TMath.h>
#include <vector>

using namespace std;
using namespace ROOT;

int main(){

  ifstream arrayFile;

  arrayFile.open("array.txt");
  char ch[100];

  vector<double> expectedUpperBounds;

  while(arrayFile >> ch){  
 
    double num = atof(ch);
    if (num < 10.0){
      expectedUpperBounds.push_back(num);
      //cout << num << endl;
    }
  }

  arrayFile.close();

  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  unsigned int nit=expectedUpperBounds.size();

  cout << "array size = " << nit << endl;

  // sort the vector with limits                                                                                                                     
  std::sort(expectedUpperBounds.begin(), expectedUpperBounds.end());

  // median for the expected limit                                                                                                                   
  median = TMath::Median(nit, &expectedUpperBounds[0]);

  // quantiles for the expected limit bands                                                                                                           
  double prob[4]; // array with quantile boundaries                                                                                                  
  prob[0] = 0.021;
  prob[1] = 0.159;
  prob[2] = 0.841;
  prob[3] = 0.979;

  // array for the results                                                                                                                            
  double quantiles[4];
  TMath::Quantiles(nit, 4, &expectedUpperBounds[0], quantiles, prob); // evaluate quantiles                                                           
  onesigma.first=quantiles[1];
  onesigma.second=quantiles[2];
  twosigma.first=quantiles[0];
  twosigma.second=quantiles[3];

  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}

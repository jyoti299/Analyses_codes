#ifndef __BINNED_DATA_HH__
#define __BINNED_DATA_HH__

#include "fit.hh"


class TH1D;

TH1D* getData(const std::vector<std::string>& filenames, const char* histname, int nbins, double* bins);

TH1D* getSignalCDF(const char* filename1, const char* histname1, const char* filename2, const char* histname2, const double BR, const double eff_h, const double eff_l, const std::string& postfix = "");

#endif

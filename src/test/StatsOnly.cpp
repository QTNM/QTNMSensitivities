/*
  StatsOnly.cpp

  Make some plots of necessary exposures 
*/

#include "Experiment.h"

#include "Formatting.h"

#include "src/include/CommonNumbers.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace sens;

int main(int argc, char *argv[])
{
  TString outputFile{argv[1]};
  TFile *fout{new TFile(outputFile, "RECREATE")};

  // Background rate
  const double bkgRate{1e-6}; 
  
  // Set up scale
  const double nAtomsMax{3e20};
  const double nAtomsMin{1e14};
  const int nAtomPnts{500};
  const double liveTime{YEAR};

  double logAtomDiff{(log10(nAtomsMax) - log10(nAtomsMin)) / double(nAtomPnts - 1)};
  TGraph *gr90CLStats{new TGraph()};
  SetGraphAttr(gr90CLStats);
  gr90CLStats->SetTitle("; Exposure [atoms #times years]; 90% CL m_{#beta} limit [eV]");

  for (int n{0}; n < nAtomPnts; n++)
  {
    double nAtoms{nAtomsMin * pow(10, double(n) * logAtomDiff)};
    double r{nAtoms * BF1EV / TLIFETIME};
    Experiment exp(r, bkgRate, liveTime, {});

    double cl90{exp.Compute90PcCL()};
    gr90CLStats->SetPoint(n, nAtoms, cl90);
  }

  fout->cd();
  gr90CLStats->Write("gr90CLStats");

  fout->Close();
  delete fout;
  return 0;
}
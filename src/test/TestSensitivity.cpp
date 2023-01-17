/*
  TestSensitivity.cpp

  Test out the basic sensitivity calculations
*/

#include "Experiment.h"

#include "Formatting.h"

#include "src/include/CommonNumbers.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"

using namespace sens;

int main(int argc, char *argv[])
{
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "RECREATE");

  const double bkg{1e-6}; // Background events per second per eV
  const double numberDensity{1e19}; // Number density per m^3
  const double volume{10}; // metres^3
  const double efficiency{0.075};
  const double sigRate{efficiency * numberDensity * volume * BF1EV / TLIFETIME};

  

  const double tMin{3600};
  const double tMax{3 * YEAR}; 
  const int nPnts{1000};
  TGraph *gr90CL = new TGraph();
  gr90CL->SetTitle("Neutrino mass limit; Live time [s]; 90% CL mass limit [eV]");
  SetGraphAttr(gr90CL);
  TGraph *grSigmaMSq = new TGraph();
  grSigmaMSq->SetTitle("Neutrino mass limit; Live time [s]; #sigma_{m_{#beta}^{2}} [eV^{2}]");
  SetGraphAttr(grSigmaMSq);

  for (int n{0}; n < nPnts; n++)
  {
    // Log distribute points
    double logDiff{(log10(tMax) - log10(tMin)) / double(nPnts - 1)};
    double time{tMin * pow(10, double(n) * logDiff)};
    Experiment exp(sigRate, bkg, time);
    gr90CL->SetPoint(n, time, exp.Compute90PcCL());
    grSigmaMSq->SetPoint(n, time, exp.ComputeSigmaMSq());
  }

  fout->cd();
  gr90CL->Write("gr90CL");
  grSigmaMSq->Write("grSigmaMSq");

  fout->Close();
  delete fout;
  return 0;
}
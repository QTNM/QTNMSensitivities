/*
  TestSensitivity.cpp

  Test out the basic sensitivity calculations
*/

#include "Experiment.h"
#include "SystematicEffect.h"

#include "Formatting.h"

#include "src/include/CommonNumbers.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"

#include <iostream>
#include <vector>
#include <memory>

using namespace sens;

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "RECREATE");

  const double bkg{1e-6}; // Background events per second per eV
  const double numberDensity{1e19}; // Number density per m^3
  const double volume{10}; // metres^3
  const double efficiency{0.075};
  const double sigRate{efficiency * numberDensity * volume * BF1EV / TLIFETIME};
  const double eOpt{sqrt(bkg / sigRate)};

  cout << "Sig rate =  " << sigRate << "\tBkg rate = " << bkg << endl;
  cout << "E_opt = " << eOpt << " eV\n";

  // Test out some of the systematic effects to check they're realistic values
  const double bField{1};
  const double gasNumberDensity{1e18};
  const double gasTemperature{30};
  std::vector<std::unique_ptr<SystematicEffect>> systs;
  systs.emplace_back(new MolecularTritium());
  systs.emplace_back(new BFieldUncertainty(1e-6, bField));
  systs.emplace_back(new GasScattering(gasNumberDensity, bField, true));
  systs.emplace_back(new ThermalBroadening(gasTemperature, true));
  for (auto const& s : systs)
  {
    cout << s->GetName() << "\tWidth = " << s->GetWidth() << " eV\t Uncertainty = " << s->GetUnc() << " eV\n";
  }

  const double tMin{10};
  const double tMax{100 * YEAR}; 
  const int nPnts{1000};
  TGraph *gr90CL = new TGraph();
  gr90CL->SetTitle("Statistics only; Live time [s]; 90% CL mass limit [eV]");
  SetGraphAttr(gr90CL);
  TGraph *gr90CLMol = new TGraph();
  gr90CLMol->SetTitle("Molecular tritium; Live time [s]; 90% CL mass limit [eV]");
  SetGraphAttr(gr90CLMol);
  gr90CLMol->SetLineColor(kRed);

  // std::vector<std::unique_ptr<SystematicEffect>> systsMol;
  std::vector<SystematicEffect*> systsMol;
  systsMol.push_back(new MolecularTritium());
  for (int n{0}; n < nPnts; n++)
  {
    // Log distribute points
    double logDiff{(log10(tMax) - log10(tMin)) / double(nPnts - 1)};
    double time{tMin * pow(10, double(n) * logDiff)};
    Experiment expStats(sigRate, bkg, time, {});
    Experiment expMol(sigRate, bkg, time, systsMol);
    gr90CL->SetPoint(n, time, expStats.Compute90PcCL());
    gr90CLMol->SetPoint(n, time, expMol.Compute90PcCL());
  }

  fout->cd();
  gr90CL->Write("gr90CL");
  gr90CLMol->Write("gr90CLMol");

  fout->Close();
  delete fout;
  return 0;
}
/*
  TestSensitivity.cpp

  Test out the basic sensitivity calculations
*/

#include <iostream>
#include <memory>
#include <vector>

#include "Experiment.h"
#include "Formatting.h"
#include "SystematicEffect.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "src/include/CommonNumbers.h"

using namespace sens;

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "RECREATE");

  const double bkg{1e-6};            // Background events per second per eV
  const double numberDensity{1e19};  // Number density per m^3
  const double volume{10};           // metres^3
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
  for (auto const &s : systs) {
    cout << s->GetName() << "\tWidth = " << s->GetWidth()
         << " eV\t Uncertainty = " << s->GetUnc() << " eV\n";
  }

  const double tMin{10};
  const double tMax{100 * YEAR};
  const int nPnts{1000};
  TGraph *gr90CL = new TGraph();
  gr90CL->SetTitle("Statistics only; Live time [s]; 90% CL mass limit [eV]");
  SetGraphAttr(gr90CL);
  TGraph *gr90CLMol = new TGraph();
  gr90CLMol->SetTitle(
      "Molecular tritium; Live time [s]; 90% CL mass limit [eV]");
  SetGraphAttr(gr90CLMol);
  gr90CLMol->SetLineColor(kRed);

  // std::vector<std::unique_ptr<SystematicEffect>> systsMol;
  std::vector<SystematicEffect> systsMol;
  systsMol.push_back(MolecularTritium());
  for (int n{0}; n < nPnts; n++) {
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

  cout << "\n";

  // Aim to reproduce Project 8 sensitivities
  auto *mg = new TMultiGraph(
      "mg",
      "arXiv:1309.7093; Effective volume [m^{3}]; 90% CL mass limit [eV]");
  const double bGraph{0.98};  // Tesla
  const double gasTempMol{30};
  const double gasTempAtom{1};
  const double liveTime{3e7};  // seconds

  // 1: Molecular tritium, density of 3 x 10^11 cm^-3
  const double nGraph1{3e11 * 1e6};
  std::vector<SystematicEffect> systs1;
  systs1.push_back(MolecularTritium());
  systs1.push_back(BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs1.push_back(GasScattering(nGraph1, bGraph, false));
  systs1.push_back(ThermalBroadening(gasTempMol, false));
  TGraph *gr1 = new TGraph();
  SetGraphAttr(gr1);
  gr1->SetTitle("T_{2}, 3 #times 10^{11} cm^{-3}");
  gr1->SetLineColor(kBlue);

  // 2: Molecular tritium, density of 3 x 10^12 cm^-3
  const double nGraph2{3e12 * 1e6};
  std::vector<SystematicEffect> systs2;
  systs2.push_back(MolecularTritium());
  systs2.push_back(BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs2.push_back(GasScattering(nGraph2, bGraph, false));
  systs2.push_back(ThermalBroadening(gasTempMol, false));
  TGraph *gr2 = new TGraph();
  SetGraphAttr(gr2);
  gr2->SetTitle("T_{2}, 3 #times 10^{12} cm^{-3}");
  gr2->SetLineColor(kBlue);
  gr2->SetLineStyle(2);

  // 3: Molecular tritium, density of 3 x 10^13 cm^-3
  const double nGraph3{3e13 * 1e6};
  std::vector<SystematicEffect> systs3;
  systs3.push_back(MolecularTritium());
  systs3.push_back(BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs3.push_back(GasScattering(nGraph3, bGraph, false));
  systs3.push_back(ThermalBroadening(gasTempMol, false));
  TGraph *gr3 = new TGraph();
  SetGraphAttr(gr3);
  gr3->SetTitle("T_{2}, 3 #times 10^{13} cm^{-3}");
  gr3->SetLineColor(kBlue);
  gr3->SetLineStyle(10);

  // 4: Atomic tritium, density of 1 x 10^12 cm^-3
  const double nGraph4{1e12 * 1e6};
  std::vector<SystematicEffect> systs4;
  systs4.push_back(BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs4.push_back(GasScattering(nGraph4, bGraph, true));
  systs4.push_back(ThermalBroadening(gasTempAtom, true));
  TGraph *gr4 = new TGraph();
  SetGraphAttr(gr4);
  gr4->SetTitle("T, 1 #times 10^{12} cm^{-3}");
  gr4->SetLineColor(kRed);

  const double vEffMin{1e-6};
  const double vEffMax{1e3};
  const double nVPnts{1000};
  double logVDiff{(log10(vEffMax) - log10(vEffMin)) / double(nVPnts - 1)};
  for (int n{0}; n < nVPnts; n++) {
    double v{vEffMin * pow(10, double(n) * logVDiff)};
    double r1{nGraph1 * v * BF1EV / TLIFETIME};
    double r2{nGraph2 * v * BF1EV / TLIFETIME};
    double r3{nGraph3 * v * BF1EV / TLIFETIME};
    double r4{nGraph4 * v * BF1EV / TLIFETIME};
    Experiment exp1(r1, bkg, liveTime, systs1);
    Experiment exp2(r2, bkg, liveTime, systs2);
    Experiment exp3(r3, bkg, liveTime, systs3);
    Experiment exp4(r4, bkg, liveTime, systs4);
    gr1->SetPoint(n, v, exp1.Compute90PcCL());
    gr2->SetPoint(n, v, exp2.Compute90PcCL());
    gr3->SetPoint(n, v, exp3.Compute90PcCL());
    gr4->SetPoint(n, v, exp4.Compute90PcCL());
  }
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);

  fout->cd();
  gr1->Write("gr1");
  gr2->Write("gr2");
  gr3->Write("gr3");
  gr4->Write("gr4");
  mg->Write();

  fout->Close();
  delete fout;
  return 0;
}
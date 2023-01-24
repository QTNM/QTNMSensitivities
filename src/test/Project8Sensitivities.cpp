/*
  Project8Sensitivities.cpp

  Aim to reproduce the some Project 8 sensitivity plots
  These plots are found in arXiv:1309.7093
*/

#include "Experiment.h"
#include "SystematicEffect.h"

#include "Formatting.h"

#include "src/include/CommonNumbers.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMultiGraph.h"

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

  const double bkgRate{1e-6};
  auto *mg90CL = new TMultiGraph("mg90CL", "arXiv:1309.7093; Effective volume [m^{3}]; 90% CL mass limit [eV]");
  auto *mgSigma = new TMultiGraph("mgSigma", "arXiv:1309.7093; Effective volume [m^{3}]; #sigma_{m_{#beta}^{2}} [eV^{2}]");
  const double bGraph{0.98};   // Tesla
  const double gasTempMol{30}; // K
  const double gasTempAtom{1}; // K
  const double liveTime{3e7};  // seconds

  // 1: Molecular tritium, density of 3 x 10^11 cm^-3
  const double nAtom1{3e11 * 1e6};
  const double nMol1{nAtom1 / 2};
  std::vector<SystematicEffect *> systs1;
  systs1.push_back(new MolecularTritium());
  systs1.push_back(new BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs1.push_back(new GasScattering(nMol1, bGraph, false));
  systs1.push_back(new ThermalBroadening(gasTempMol, false));
  TGraph *gr90CL1 = new TGraph();
  SetGraphAttr(gr90CL1);
  gr90CL1->SetTitle("T_{2}, 3 #times 10^{11} cm^{-3}");
  gr90CL1->SetLineColor(kBlue);
  TGraph *grSigma1 = new TGraph();
  SetGraphAttr(grSigma1);
  grSigma1->SetTitle("T_{2}, 3 #times 10^{11} cm^{-3}");
  grSigma1->SetLineColor(kBlue);

  // 2: Molecular tritium, density of 3 x 10^12 cm^-3
  const double nAtom2{3e12 * 1e6};
  const double nMol2{nAtom2 / 2};
  std::vector<SystematicEffect *> systs2;
  systs2.push_back(new MolecularTritium());
  systs2.push_back(new BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs2.push_back(new GasScattering(nMol2, bGraph, false));
  systs2.push_back(new ThermalBroadening(gasTempMol, false));
  TGraph *gr90CL2 = new TGraph();
  SetGraphAttr(gr90CL2);
  gr90CL2->SetTitle("T_{2}, 3 #times 10^{12} cm^{-3}");
  gr90CL2->SetLineColor(kBlue);
  gr90CL2->SetLineStyle(2);
  TGraph *grSigma2 = new TGraph();
  SetGraphAttr(grSigma2);
  grSigma2->SetTitle("T_{2}, 3 #times 10^{12} cm^{-3}");
  grSigma2->SetLineColor(kBlue);
  grSigma2->SetLineStyle(2);

  // 3: Molecular tritium, density of 3 x 10^13 cm^-3
  const double nAtom3{3e13 * 1e6};
  const double nMol3{nAtom3 / 2};
  std::vector<SystematicEffect *> systs3;
  systs3.push_back(new MolecularTritium());
  systs3.push_back(new BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs3.push_back(new GasScattering(nMol3, bGraph, false));
  systs3.push_back(new ThermalBroadening(gasTempMol, false));
  TGraph *gr90CL3 = new TGraph();
  SetGraphAttr(gr90CL3);
  gr90CL3->SetTitle("T_{2}, 3 #times 10^{13} cm^{-3}");
  gr90CL3->SetLineColor(kBlue);
  gr90CL3->SetLineStyle(10);
  TGraph *grSigma3 = new TGraph();
  SetGraphAttr(grSigma3);
  grSigma3->SetTitle("T_{2}, 3 #times 10^{13} cm^{-3}");
  grSigma3->SetLineColor(kBlue);
  grSigma3->SetLineStyle(10);

  // 4: Atomic tritium, density of 1 x 10^12 cm^-3
  const double nMol4{1e12 * 1e6};
  const double nAtom4{nMol4};
  std::vector<SystematicEffect *> systs4;
  systs4.push_back(new BFieldUncertainty(0.1 * 1e-6 * bGraph, bGraph));
  systs4.push_back(new GasScattering(nMol4, bGraph, true));
  systs4.push_back(new ThermalBroadening(gasTempAtom, true));
  TGraph *gr90CL4 = new TGraph();
  SetGraphAttr(gr90CL4);
  gr90CL4->SetTitle("T, 1 #times 10^{12} cm^{-3}; Effective volume [m^{3}];  90% CL mass limit [eV]");
  gr90CL4->SetLineColor(kRed);
  TGraph *grSigma4 = new TGraph();
  SetGraphAttr(grSigma4);
  grSigma4->SetTitle("T, 1 #times 10^{12} cm^{-3}; Effective volume [m^{3}]; #sigma_{m_{#beta}^{2}} [eV^{2}]");
  grSigma4->SetLineColor(kRed);

  // Print the systematic widths
  cout << "Molecular tritium, n = 3 x 10^11 cm^-3\n";
  double sigmaMSyst1{0};
  for (auto const &s : systs1)
  {
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\n";
    sigmaMSyst1 += pow(s->GetWidth(), 4) * pow(s->GetUnc() / s->GetWidth(), 2);
  }
  sigmaMSyst1 = 4 * sqrt(sigmaMSyst1);
  cout << "sigmaM = " << sqrt(sigmaMSyst1) << "eV\n";

  cout << "\nMolecular tritium, n = 3 x 10^12 cm^-3\n";
  double sigmaMSyst2{0};
  for (auto const &s : systs2)
  {
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\n";
    sigmaMSyst2 += pow(s->GetWidth(), 4) * pow(s->GetUnc() / s->GetWidth(), 2);
  }
  sigmaMSyst2 = 4 * sqrt(sigmaMSyst2);
  cout << "sigmaM = " << sqrt(sigmaMSyst2) << "eV\n";

  cout << "\nMolecular tritium, n = 3 x 10^13 cm^-3\n";
  double sigmaMSyst3{0};
  for (auto const &s : systs3)
  {
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\n";
    sigmaMSyst3 += pow(s->GetWidth(), 4) * pow(s->GetUnc() / s->GetWidth(), 2);
  }
  sigmaMSyst3 = 4 * sqrt(sigmaMSyst3);
  cout << "sigmaM = " << sqrt(sigmaMSyst3) << "eV\n";

  cout << "\nAtomic tritium, n = 1 x 10^12 cm^-3\n";
  double sigmaMSyst4{0};
  for (auto const &s : systs4)
  {
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\n";
    sigmaMSyst4 += pow(s->GetWidth(), 4) * pow(s->GetUnc() / s->GetWidth(), 2);
  }
  sigmaMSyst4 = 4 * sqrt(sigmaMSyst4);
  cout << "sigmaM = " << sqrt(sigmaMSyst4) << "eV\n";

  const double vEffMin{1e-6};
  const double vEffMax{1e3};
  const double nVPnts{1000};
  double logVDiff{(log10(vEffMax) - log10(vEffMin)) / double(nVPnts - 1)};
  for (int n{0}; n < nVPnts; n++)
  {
    double v{vEffMin * pow(10, double(n) * logVDiff)};
    double r1{nAtom1 * v * BF1EV / TLIFETIME};
    double r2{nAtom2 * v * BF1EV / TLIFETIME};
    double r3{nAtom3 * v * BF1EV / TLIFETIME};
    double r4{nAtom4 * v * BF1EV / TLIFETIME};
    Experiment exp1(r1, bkgRate, liveTime, systs1);
    Experiment exp2(r2, bkgRate, liveTime, systs2);
    Experiment exp3(r3, bkgRate, liveTime, systs3);
    Experiment exp4(r4, bkgRate, liveTime, systs4);
    gr90CL1->SetPoint(n, v, exp1.Compute90PcCL());
    gr90CL2->SetPoint(n, v, exp2.Compute90PcCL());
    gr90CL3->SetPoint(n, v, exp3.Compute90PcCL());
    gr90CL4->SetPoint(n, v, exp4.Compute90PcCL());
    grSigma1->SetPoint(n, v, exp1.ComputeSigmaMSq());
    grSigma2->SetPoint(n, v, exp2.ComputeSigmaMSq());
    grSigma3->SetPoint(n, v, exp3.ComputeSigmaMSq());
    grSigma4->SetPoint(n, v, exp4.ComputeSigmaMSq());
  }
  gr90CL1->GetXaxis()->SetRangeUser(1e-6, 1e3);
  gr90CL2->GetXaxis()->SetRangeUser(1e-6, 1e3);
  gr90CL3->GetXaxis()->SetRangeUser(1e-6, 1e3);
  gr90CL4->GetXaxis()->SetRangeUser(1e-6, 1e3);
  grSigma1->GetXaxis()->SetRangeUser(1e-6, 1e3);
  grSigma2->GetXaxis()->SetRangeUser(1e-6, 1e3);
  grSigma3->GetXaxis()->SetRangeUser(1e-6, 1e3);
  grSigma4->GetXaxis()->SetRangeUser(1e-6, 1e3);

  gr90CL1->GetYaxis()->SetRangeUser(0.02, 10);
  gr90CL2->GetYaxis()->SetRangeUser(0.02, 10);
  gr90CL3->GetYaxis()->SetRangeUser(0.02, 10);
  gr90CL4->GetYaxis()->SetRangeUser(0.02, 10);
  grSigma1->GetYaxis()->SetRangeUser(0.0003, 100);
  grSigma2->GetYaxis()->SetRangeUser(0.0003, 100);
  grSigma3->GetYaxis()->SetRangeUser(0.0003, 100);
  grSigma4->GetYaxis()->SetRangeUser(0.0003, 100);

  fout->cd();
  mg90CL->Add(gr90CL1);
  mg90CL->Add(gr90CL2);
  mg90CL->Add(gr90CL3);
  mg90CL->Add(gr90CL4);
  gr90CL1->Write("gr90CL1");
  gr90CL2->Write("gr90CL2");
  gr90CL3->Write("gr90CL3");
  gr90CL4->Write("gr90CL4");
  mg90CL->Write();

  mgSigma->Add(grSigma1);
  mgSigma->Add(grSigma2);
  mgSigma->Add(grSigma3);
  mgSigma->Add(grSigma4);
  grSigma1->Write("grSigma1");
  grSigma2->Write("grSigma2");
  grSigma3->Write("grSigma3");
  grSigma4->Write("grSigma4");
  mgSigma->Write();

  // Now look at reproducing the plot from arXiv:2102.00594
  const double energyResPPM{3e-6}; // energy resolution = 3 ppm RMS 
  const double energyRes{energyResPPM * 18.6e3}; // eV
  // Experiment A: Molecular tritium, n = 2 x 10^17 m^-3
  const double nAtomsA{2e17};
  const double nMolA{nAtomsA / 2};
  std::vector<SystematicEffect *> systsA;
  systsA.push_back(new MolecularTritium());
  systsA.push_back(new SystematicEffect(energyRes, energyRes * 0.01, "Energy Res. Syst"));
  systsA.push_back(new GasScattering(nMolA, bGraph, false));
  systsA.push_back(new ThermalBroadening(gasTempMol, false));
  TGraph *gr90CLA = new TGraph();
  SetGraphAttr(gr90CLA);
  gr90CLA->SetLineColor(kBlue);
  gr90CLA->SetTitle("arXiv:2102.00594; Volume #times Efficiency #times Time [m^{3} yr]; 90% CL m_{#beta} [eV]");

  // Experiment B: Atomic tritium, n = 7 x 10^17 m^-3
  const double nAtomsB{7e17};
  const double nMolB{nAtomsB};
  std::vector<SystematicEffect *> systsB;
  systsB.push_back(new SystematicEffect(energyRes, energyRes * 0.01, "Energy Res. Syst"));
  systsB.push_back(new GasScattering(nMolB, bGraph, true));
  systsB.push_back(new ThermalBroadening(gasTempAtom, true));
  TGraph *gr90CLB = new TGraph();
  SetGraphAttr(gr90CLB);
  gr90CLB->SetLineColor(kRed);
  gr90CLB->SetTitle("T, 7 #times 10^{17} m^{-3}; Volume #times Efficiency #times Time [m^{3} yr]; 90% CL m_{#beta} [eV]");

  cout << "============================\n";
  cout << "arXiv:2102.00594\n";
  cout << "============================\n";

  cout << "\nMolecular tritium\n";
  for (auto const & s : systsA)
  {
    double mSigma{sqrt(4 * s->GetUnc() * s->GetWidth())};
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\t sigma_m = " << mSigma << " eV\n";
  }

  cout << "\nAtomic tritium\n";
  for (auto const & s : systsB)
  {
    double mSigma{sqrt(4 * s->GetUnc() * s->GetWidth())};
    cout << s->GetName() << "\t" << s->GetWidth() << " eV\t" << s->GetUnc() << " eV\t sigma_m = " << mSigma << " eV\n";
  }

  // Exposures are given in m^3 yr
  const double exposureMin{1e-6};
  const double exposureMax{1e3};
  double logExpDiff{(log10(exposureMax) - log10(exposureMin)) / double(nVPnts - 1)};
  for (int n{0}; n < nVPnts; n++)
  {
    double exp{exposureMin * pow(10, double(n) * logExpDiff)};
    double timeYear{1};       // Time in years
    double v{exp / timeYear}; // Volume in m^3
    double rA{nAtomsA * v * BF1EV / TLIFETIME};
    double rB{nAtomsB * v * BF1EV / TLIFETIME};
    Experiment expA(rA, bkgRate, YEAR, systsA);
    Experiment expB(rB, bkgRate, YEAR, systsB);
    gr90CLA->SetPoint(n, exp, expA.Compute90PcCL());
    gr90CLB->SetPoint(n, exp, expB.Compute90PcCL());
  }
  gr90CLA->GetXaxis()->SetRangeUser(1e-6, 1e3);
  gr90CLB->GetXaxis()->SetRangeUser(1e-6, 1e3);
  gr90CLA->GetYaxis()->SetRangeUser(0.02, 10);
  gr90CLB->GetYaxis()->SetRangeUser(0.02, 10);
  fout->cd();
  gr90CLA->Write("gr90CLA");
  gr90CLB->Write("gr90CLB");

  fout->Close();
  delete fout;
  return 0;
}
/*
  WhitePaperPlots.cpp

  03/08/2023
*/

#include <unistd.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Experiment.h"
#include "Formatting.h"
#include "SystematicEffect.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "src/include/CommonNumbers.h"
#include "src/include/FundamentalConstants.h"

using namespace sens;

int main(int argc, char* argv[]) {
  int opt{};
  std::string plotFileStr{};
  while ((opt = getopt(argc, argv, ":o:")) != -1) {
    switch (opt) {
      case 'o':
        plotFileStr = optarg;
        break;

      case ':':
        std::cout << "Option requires an argument\n";
        break;

      case '?':
        std::cout << "Unknown option\n";
        break;
    }
  }

  TString plotFile{plotFileStr};
  TFile fout(plotFile, "recreate");

  // Basic setup
  const double nAtom1{1e12 * 1e6};
  const double eKE{18.6e3};  // eV
  const double bField{0.7};  // Tesla
  const double gamma{1 + eKE * QE / (EMASS * CLIGHT * CLIGHT)};
  const double cycFreq{(1 / (2 * M_PI)) * QE * bField /
                       (EMASS + eKE * QE / pow(CLIGHT, 2))};
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";
  std::cout << "Gamma = " << gamma << std::endl;
  const double deltaB{0.1e-6 * bField};  // Tesla
  const double systWidth{0.01};
  const double deltaFOverF{0.1e-6};
  const double deltaEOverE{deltaFOverF * gamma / (gamma - 1)};
  const double deltaE{deltaEOverE * eKE};

  // Spread in atom speeds
  const double sigmaV{20};  // m/s
  // Atom temperature
  const double effT{T_MASS_DA * DALTON * sigmaV * sigmaV / (3 * KB)};  // K
  std::cout << "Effective atom temperature = " << effT << " K\n";
  // Define thermal broadening systematic
  ThermalBroadening thermalSyst(effT, true);
  // Magnetic field systematics
  BFieldUncertainty bSyst100ppt(1e-7 * bField, bField);
  BFieldUncertainty bSyst1ppm(1e-6 * bField, bField);
  // Gas scattering systematic
  GasScattering scatSyst(nAtom1, bField, true);
  // Energy resolution systematics
  SystematicEffect eResSyst100ppt(eKE * 0.1e-6, eKE * 0.1e-6 * 0.01,
                                  "100ppt energy resolution");
  SystematicEffect eResSyst1ppm(eKE * 1e-6, eKE * 1e-6 * 0.01,
                                "1ppm energy resolution");
  SystematicEffect eResSyst10ppm(eKE * 10e-6, eKE * 10e-6 * 0.01,
                                 "10ppm energy resolution");
  SystematicEffect eResSyst30ppm(eKE * 30e-6, eKE * 30e-6 * 0.01,
                                 "30ppm energy resolution");

  ////////////////////// Check optimum density //////////////////////
  const double expVolumeDensCheck{10};      // m^3
  const double runTimeDensCheck{2 * YEAR};  // seconds
  const double rhoMin{1e15};                // m^-3
  const double rhoMax{1e20};                // m^-3
  const uint nRhoPnts{500};
  const double logRhoDiff{(log10(rhoMax) - log10(rhoMin)) /
                          double(nRhoPnts - 1)};
  // Create graph
  auto grRhoVar500mT = new TGraph();
  SetGraphAttr(grRhoVar500mT);
  grRhoVar500mT->SetLineWidth(3);
  grRhoVar500mT->SetTitle(
      "0.5 T; Number density [cm^{-3}]; 90% CL on m_{#beta} [eV]");
  grRhoVar500mT->SetLineColor(kCyan + 1);
  auto grRhoVar700mT = new TGraph();
  SetGraphAttr(grRhoVar700mT);
  grRhoVar700mT->SetLineWidth(3);
  grRhoVar700mT->SetTitle(
      "0.7 T; Number density [cm^{-3}]; 90% CL on m_{#beta} [eV]");
  grRhoVar700mT->SetLineColor(kBlack);
  auto grRhoVar1000mT = new TGraph();
  SetGraphAttr(grRhoVar1000mT);
  grRhoVar1000mT->SetLineWidth(3);
  grRhoVar1000mT->SetTitle(
      "1.0 T; Number density [cm^{-3}]; 90% CL on m_{#beta} [eV]");
  grRhoVar1000mT->SetLineColor(kOrange + 1);

  // Loop through different densities
  for (int i{0}; i < nRhoPnts; i++) {
    double rho{rhoMin * pow(10, double(i) * logRhoDiff)};
    // Create scattering systematic
    GasScattering tempScatSyst500mT(rho, 0.5, true);
    GasScattering tempScatSyst700mT(rho, 0.7, true);
    GasScattering tempScatSyst1000mT(rho, 1.0, true);

    std::vector<SystematicEffect> systsDensCheck500mT{
        eResSyst1ppm, thermalSyst, bSyst100ppt, tempScatSyst500mT};
    std::vector<SystematicEffect> systsDensCheck700mT{
        eResSyst1ppm, thermalSyst, bSyst100ppt, tempScatSyst700mT};
    std::vector<SystematicEffect> systsDensCheck1000mT{
        eResSyst1ppm, thermalSyst, bSyst100ppt, tempScatSyst1000mT};

    double r{rho * expVolumeDensCheck * BF1EV / TLIFETIME};
    double bkgDensCheck{1e-3 * r};
    Experiment expDensCheck500mT(r - bkgDensCheck, bkgDensCheck,
                                 runTimeDensCheck, systsDensCheck500mT);
    grRhoVar500mT->SetPoint(i, rho / 1e6, expDensCheck500mT.Compute90PcCL());
    Experiment expDensCheck700mT(r - bkgDensCheck, bkgDensCheck,
                                 runTimeDensCheck, systsDensCheck700mT);
    grRhoVar700mT->SetPoint(i, rho / 1e6, expDensCheck700mT.Compute90PcCL());
    Experiment expDensCheck1000mT(r - bkgDensCheck, bkgDensCheck,
                                  runTimeDensCheck, systsDensCheck1000mT);
    grRhoVar1000mT->SetPoint(i, rho / 1e6, expDensCheck1000mT.Compute90PcCL());
  }

  fout.cd();
  grRhoVar500mT->Write("grRhoVar500mT");
  grRhoVar700mT->Write("grRhoVar700mT");
  grRhoVar1000mT->Write("grRhoVar1000mT");

  // Now plot a few different options
  // Keep 0.7 T as field choice in all cases and don't change the temperature of
  // the atoms

  // Sample 1: N = 10^13 cm^-3
  // 1 ppm field uniformity and 10 ppm energy resolution
  const double n1{1e19};  // m^-3
  std::vector<SystematicEffect> systs1{thermalSyst};
  systs1.push_back(bSyst1ppm);
  systs1.push_back(eResSyst10ppm);
  systs1.push_back(GasScattering(n1, bField, true));
  auto gr90CL_1 = new TGraph();
  SetGraphAttr(gr90CL_1);
  gr90CL_1->SetLineWidth(3);
  gr90CL_1->SetLineColor(kRed);
  gr90CL_1->SetTitle(
      "N = 10^{13} cm^{-3}, #sigma_{B} = 1 ppm, #sigma_{E} = 10 ppm");
  gr90CL_1->GetXaxis()->SetTitle("Exposure [m^{3} yr]");
  gr90CL_1->GetYaxis()->SetTitle("90% CL on m_{#beta} [eV]");

  auto grSigma_1 = new TGraph();
  SetGraphAttr(grSigma_1);
  grSigma_1->SetLineWidth(3);
  grSigma_1->SetLineColor(kRed);
  grSigma_1->SetTitle(
      "N = 10^{13} cm^{-3}, #sigma_{B} = 1 ppm, #sigma_{E} = 10 ppm");
  grSigma_1->GetXaxis()->SetTitle(
      "Efficiency #times Volume #times Time [m^{3} yr]");
  grSigma_1->GetYaxis()->SetTitle("#sigma_{m^{2}_{#beta}} [eV^{2}]");

  // Sample 2: N = 2 x 10^12 cm^-3
  // 1 ppm field uniformity and 10 ppm energy resolution
  const double n2{2e18};  // m^-3
  std::vector<SystematicEffect> systs2{thermalSyst};
  systs2.push_back(bSyst1ppm);
  systs2.push_back(eResSyst10ppm);
  systs2.push_back(GasScattering(n2, bField, true));
  auto gr90CL_2 = new TGraph();
  SetGraphAttr(gr90CL_2);
  gr90CL_2->SetLineWidth(3);
  gr90CL_2->SetLineColor(kCyan + 1);
  gr90CL_2->SetTitle(
      "N = 2 #times 10^{12} cm^{-3}, #sigma_{B} = 1 ppm, #sigma_{E} = 10 ppm");
  gr90CL_2->GetXaxis()->SetTitle(
      "Efficiency #times Volume #times Time [m^{3} yr]");
  gr90CL_2->GetYaxis()->SetTitle("90% CL on m_{#beta} [eV]");

  auto grSigma_2 = new TGraph();
  SetGraphAttr(grSigma_2);
  grSigma_2->SetLineWidth(3);
  grSigma_2->SetLineColor(kCyan + 1);
  grSigma_2->SetTitle(
      "N = 2 #times 10^{12} cm^{-3}, #sigma_{B} = 1 ppm, #sigma_{E} = 10 ppm");
  grSigma_2->GetXaxis()->SetTitle(
      "Efficiency #times Volume #times Time [m^{3} yr]");
  grSigma_2->GetYaxis()->SetTitle("#sigma_{m^{2}_{#beta}} [eV^{2}]");

  // Sample 3: N = 2 x 10^12 cm^-3
  // 0.1 ppm field uniformity and 0.1 ppm energy resolution
  const double n3{2e18};  // m^-3
  std::vector<SystematicEffect> systs3{thermalSyst};
  systs3.push_back(bSyst100ppt);
  systs3.push_back(eResSyst100ppt);
  systs3.push_back(GasScattering(n3, bField, true));
  auto gr90CL_3 = new TGraph();
  SetGraphAttr(gr90CL_3);
  gr90CL_3->SetLineWidth(3);
  gr90CL_3->SetLineColor(kOrange + 1);
  gr90CL_3->SetTitle(
      "N = 2 #times 10^{12} cm^{-3}, #sigma_{B} = 0.1 ppm, #sigma_{E} = 0.1 "
      "ppm");
  gr90CL_3->GetXaxis()->SetTitle(
      "Efficiency #times Volume #times Time [m^{3} yr]");
  gr90CL_3->GetYaxis()->SetTitle("90% CL on m_{#beta} [eV^{2}]");

  auto grSigma_3 = new TGraph();
  SetGraphAttr(grSigma_3);
  grSigma_3->SetLineWidth(3);
  grSigma_3->SetLineColor(kOrange + 1);
  grSigma_3->SetTitle(
      "N = 2 #times 10^{12} cm^{-3}, #sigma_{B} = 1 ppm, #sigma_{E} = 10 ppm");
  grSigma_3->GetXaxis()->SetTitle(
      "Efficiency #times Volume #times Time [m^{3} yr]");
  grSigma_3->GetYaxis()->SetTitle("#sigma_{m^{2}_{#beta}} [eV^{2}]");

  const uint nExpPnts{500};
  const double expMin{1e-5};
  const double expMax{1e4};
  const double logExpDiff{(log10(expMax) - log10(expMin)) /
                          double(nExpPnts - 1)};
  const double v{10};  // m^3
  for (int i{0}; i < nExpPnts; i++) {
    const double exposure{expMin * pow(10, double(i) * logExpDiff)};
    const double timeYears{exposure / v};
    const double time{timeYears * YEAR};
    double r1{n1 * v * BF1EV / TLIFETIME};
    double b1{1e-4 * r1};
    double r2{n2 * v * BF1EV / TLIFETIME};
    double b2{1e-4 * r2};
    double r3{n3 * v * BF1EV / TLIFETIME};
    double b3{1e-4 * r3};
    Experiment exp1(r1 - b1, b1 + 1e-6, time, systs1);
    Experiment exp2(r2 - b2, b2 + 1e-6, time, systs2);
    Experiment exp3(r3 - b3, b3 + 1e-6, time, systs3);
    gr90CL_1->SetPoint(i, exposure, exp1.Compute90PcCL());
    gr90CL_2->SetPoint(i, exposure, exp2.Compute90PcCL());
    gr90CL_3->SetPoint(i, exposure, exp3.Compute90PcCL());
    grSigma_1->SetPoint(i, exposure, exp1.ComputeSigmaMSq());
    grSigma_2->SetPoint(i, exposure, exp2.ComputeSigmaMSq());
    grSigma_3->SetPoint(i, exposure, exp3.ComputeSigmaMSq());
  }

  fout.cd();
  gr90CL_1->Write("gr90CL_1");
  gr90CL_2->Write("gr90CL_2");
  gr90CL_3->Write("gr90CL_3");
  grSigma_1->Write("grSigma_1");
  grSigma_2->Write("grSigma_2");
  grSigma_3->Write("grSigma_3");

  fout.Close();
  return 0;
}
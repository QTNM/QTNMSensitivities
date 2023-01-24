/*
  DopplerBroadening.cpp

  Plot the Doppler broadening as a function of temperature
*/

#include "SystematicEffect.h"

#include "Formatting.h"

#include "src/include/CommonNumbers.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"

#include <iostream>
#include <vector>

using namespace sens;

int main(int argc, char *argv[])
{
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "RECREATE");

  const double TMax{30};     // Kelvin
  const double TMin{100e-3}; // Kelvin
  const int nPnts{300};
  double logTDiff{(log10(TMax) - log10(TMin)) / double(nPnts - 1)};

  TGraph *grSigmaEAtom = new TGraph();
  SetGraphAttr(grSigmaEAtom);
  grSigmaEAtom->SetTitle("Atomic tritium; Atom/molecule T [K]; #sigma_{E} [eV]");
  grSigmaEAtom->SetLineColor(kRed);
  TGraph *grSigmaEMol  = new TGraph();
  SetGraphAttr(grSigmaEMol);
  grSigmaEMol->SetTitle("Molecular tritium; Atom/molecule T [K]; #sigma_{E} [eV]");
  grSigmaEMol->SetLineColor(kBlue);

  for (int n{0}; n < nPnts; n++)
  {
    double T{TMin * pow(10, double(n) * logTDiff)};
    ThermalBroadening thermAtom(T, true);
    grSigmaEAtom->SetPoint(n, T, thermAtom.GetWidth());
    ThermalBroadening thermMol(T, false);
    grSigmaEMol->SetPoint(n, T, thermMol.GetWidth());
  }
  fout->cd();
  grSigmaEAtom->Write("grSigmaEAtom");
  grSigmaEMol->Write("grSigmaEMol");

  fout->Close();
  delete fout;
  return 0;
}
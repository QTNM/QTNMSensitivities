// SystematicEffect.cpp

#include "SystematicEffect.h"
#include "src/include/CommonNumbers.h"
#include "src/include/FundamentalConstants.h"

#include <math.h>

void sens::SystematicEffect::SetRMSWidth(double rms)
{
  rmsWidth = rms;
}

void sens::SystematicEffect::SetWidthUnc(double unc)
{
  widthUnc = unc;
}

void sens::SystematicEffect::SetSystName(std::string name)
{
  systName = name;
}

sens::MolecularTritium::MolecularTritium()
{
  SetSystName("Molecular Tritium Syst");
  // Broadening of endpoint from molecular tritium in eV
  // From arXiv:1502.03497[nucl-ex]
  SetRMSWidth(0.4);
  SetWidthUnc(GetWidth() * 0.01);
}

sens::BFieldUncertainty::BFieldUncertainty(double fieldVar, double bAvg)
{
  SetSystName("Magnetic Field Syst");
  // Need to account for how uncertainty in frequency relates to energy
  SetRMSWidth((GAMMA_END / (GAMMA_END - 1)) * 18.6e3 * fieldVar / bAvg);
  SetWidthUnc(GetWidth() * 0.01);
}

sens::GasScattering::GasScattering(double n, double b, bool isAtomicTritium)
{
  SetSystName("Gas Scattering Syst");
  // Set the correct scattering cross section
  double xsec{isAtomicTritium ? 9e-23 : 3.4e-22};
  double f_c{(1 / (2 * M_PI)) * QE * b / (EMASS + 18.6e3 * QE / (CLIGHT * CLIGHT))};
  SetRMSWidth((18.6e3 / (GAMMA_END - 1)) * BETA_END * CLIGHT * xsec * n / (2 * M_PI * f_c));
  SetWidthUnc(GetWidth() * 0.01);
}

sens::ThermalBroadening::ThermalBroadening(double T, bool isAtomicTritium)
{
  SetSystName("Thermal Broadening Syst");
  double mass{isAtomicTritium ? DALTON * T_MASS_DA : DALTON * T2_MASS_DA};
  // RMS frequency change
  double sigmaFOverF{sqrt(KB * T / (mass * CLIGHT * CLIGHT))};
  SetRMSWidth(sigmaFOverF * 18.6e3 * GAMMA_END / (GAMMA_END - 1));
  SetWidthUnc(GetWidth() * 0.01);
}
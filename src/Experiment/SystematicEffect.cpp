// SystematicEffect.cpp

#include "SystematicEffect.h"
#include "src/include/CommonNumbers.h"
#include "src/include/FundamentalConstants.h"

#include <math.h>

sens::MolecularTritium::MolecularTritium()
{
  // Broadening of endpoint from molecular tritium in eV
  // From arXiv:1502.03497[nucl-ex]
  rmsWidth = 0.4;
  widthUnc = rmsWidth * 0.01;
}

sens::BFieldUncertainty::BFieldUncertainty(double fieldVar, double bAvg)
{
  // Need to account for how uncertainty in frequency relates to energy
  rmsWidth = (GAMMA_END / (GAMMA_END - 1)) * 18.6e3 * fieldVar / bAvg;
  widthUnc = rmsWidth * 0.01;
}

sens::GasScattering::GasScattering(double n, double b, bool isAtomicTritium)
{
  // Set the correct scattering cross section
  double xsec{isAtomicTritium ? 9e-15 : 3.4e-14};
  double f_c{(1 / (2 * M_PI)) * QE * b / (EMASS + 18.6e3 * QE / (CLIGHT * CLIGHT))};
  rmsWidth = (18.6e3 / (GAMMA_END - 1)) * BETA_END * CLIGHT * xsec * n / (2 * M_PI * f_c);
  widthUnc = rmsWidth * 0.01;
}

sens::ThermalBroadening::ThermalBroadening(double T, bool isAtomicTritium)
{
  double mass{isAtomicTritium ? DALTON * T_MASS_DA : DALTON * T2_MASS_DA};
  // Calculate average velocity of atom or molecule
  double vAvg{sqrt(3 * KB * T / mass)};
  // Now calculate the energy difference for an endpoint electron
  double beta1{(BETA_END * CLIGHT + vAvg) / CLIGHT};
  double deltaE{((1 / sqrt(1 - beta1 * beta1)) - (1 / sqrt(1 - BETA_END * BETA_END))) * EMASS * CLIGHT * CLIGHT / QE};
  rmsWidth = deltaE;
  widthUnc = rmsWidth * 0.01;
}
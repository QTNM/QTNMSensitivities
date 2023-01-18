// SystematicEffect.cpp

#include "SystematicEffect.h"
#include "src/include/CommonNumbers.h"

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
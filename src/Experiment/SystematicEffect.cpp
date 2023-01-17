// SystematicEffect.cpp

#include "SystematicEffect.h"

sens::MolecularTritium::MolecularTritium()
{ 
  // Broadening of endpoint from molecular tritium in eV
  // From arXiv:1502.03497[nucl-ex]
  rmsWidth = 0.4;
  widthUnc = rmsWidth * 0.01;
}
/*
  FundamentalConstants.h

  Several constants of nature it's useful to keep track of
*/

#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

namespace sens
{
  // Electron rest mass in kilograms
  inline constexpr double EMASS{9.1093837015e-31};

  // Permeability of free space
  inline constexpr double EPSILON0{8.8541878128e-12};

  // Permeability of free space
  inline constexpr double MU0{1.25663706212e-6};

  // Elementary charge in coulombs
  inline constexpr double QE{1.602176634e-19};

  // Speed of light in a vacuum
  inline constexpr double CLIGHT{299792458};

  // Atomic mass unit/Dalton in kg
  inline constexpr double DALTON{1.6605390666e-27};
}

#endif
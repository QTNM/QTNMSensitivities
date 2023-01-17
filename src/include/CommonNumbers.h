/*
  CommonNumbers.h

  These aren't fundamental constants but they are frequently used
*/

#include <cmath>

namespace sens
{
  // One year in seconds
  inline constexpr double YEAR{31536000};

  // Branching fraction to last eV of tritium spectrum 
  inline constexpr double BF1EV{2e-13};

  // Mean tritium lifetime in seconds
  inline constexpr double TLIFETIME{560975924};
}
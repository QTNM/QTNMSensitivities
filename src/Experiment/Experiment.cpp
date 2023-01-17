/*
  Experiment.cpp
*/

#include "Experiment.h"

#include <cmath>

double sens::Experiment::ComputeSigmaMSq()
{
  // Calculate the optimal stats only energy interval
  double eOpt{sqrt(bkgRate / sigRate)};
  double premult{2 / (3 * sigRate * liveTime)};
  return premult * sqrt(sigRate * liveTime * eOpt + bkgRate * liveTime / eOpt);
}

double sens::Experiment::Compute90PcCL()
{
  return sqrt(1.28 * ComputeSigmaMSq());
}
/*
  Experiment.cpp
*/

#include "Experiment.h"

#include <cmath>

double sens::Experiment::ComputeSigmaMSq(double liveTime)
{
  // Calculate the optimal stats only energy interval
  double eOpt{sqrt(bkgRate / sigRate)};
  double premult{2 / (3 * sigRate * liveTime)};
  return premult * sqrt(sigRate * liveTime * eOpt + bkgRate * liveTime / eOpt);
}
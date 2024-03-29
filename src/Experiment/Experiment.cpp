/*
  Experiment.cpp
*/

#include "Experiment.h"

#include <cmath>

#include "SystematicEffect.h"

sens::Experiment::Experiment(double r, double b, double t,
                             std::vector<SystematicEffect> systs)
    : sigRate(r), bkgRate(b), liveTime(t), systematics(systs) {}

sens::Experiment::~Experiment() { systematics.clear(); }

double sens::Experiment::ComputeSigmaMSq() {
  // Calculate the optimal energy interval
  double deltaE{bkgRate / sigRate};
  double systUnc{0};
  for (auto &s : systematics) {
    deltaE += 8 * log(2) * s.GetWidth() * s.GetWidth();
    systUnc += pow(s.GetWidth(), 4) * pow(s.GetUnc() / s.GetWidth(), 2);
  }
  deltaE = sqrt(deltaE);
  systUnc = 4 * sqrt(systUnc);

  // Statistical uncertainty
  double statUnc{
      2 / (3 * sigRate * liveTime) *
      sqrt(sigRate * liveTime * deltaE + bkgRate * liveTime / deltaE)};

  return sqrt(systUnc * systUnc + statUnc * statUnc);
}

double sens::Experiment::Compute90PcCL() {
  return sqrt(1.28 * ComputeSigmaMSq());
}
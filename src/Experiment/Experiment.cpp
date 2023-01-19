/*
  Experiment.cpp
*/

#include "Experiment.h"
#include "SystematicEffect.h"

#include <cmath>

sens::Experiment::Experiment(double r, double b, double t,
                             /*std::vector<std::unique_ptr<SystematicEffect>> systs)*/
                             std::vector<SystematicEffect*> systs)
{
  sigRate = r;
  bkgRate = b;
  liveTime = t;
  for (auto const &s : systs)
  {
    SystematicEffect systCopy(s->GetWidth(), s->GetUnc(), s->GetName());
    systematics.emplace_back(std::make_unique<SystematicEffect>(systCopy));
  }
}

sens::Experiment::~Experiment()
{
  systematics.clear();
}

double sens::Experiment::ComputeSigmaMSq()
{
  // Calculate the optimal energy interval
  double deltaE{bkgRate / sigRate};
  double systUnc{0};
  for (auto const &s : systematics)
  {
    deltaE += 8 * log(2) * s->GetWidth() * s->GetWidth();
    systUnc += pow(s->GetWidth(), 4) * pow(s->GetUnc() / s->GetWidth(), 2);
  }
  deltaE = sqrt(deltaE);
  systUnc = 4 * sqrt(systUnc);

  // Statistical uncertainty
  double statUnc{2 / (3 * sigRate * liveTime) * sqrt(sigRate * liveTime * deltaE + bkgRate * liveTime / deltaE)};

  return sqrt(systUnc * systUnc + statUnc * statUnc);
}

double sens::Experiment::Compute90PcCL()
{
  return sqrt(1.28 * ComputeSigmaMSq());
}
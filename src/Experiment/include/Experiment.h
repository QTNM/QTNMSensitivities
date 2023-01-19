#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "SystematicEffect.h"

#include <vector>
#include <memory>

namespace sens
{
  class Experiment
  {
  private:
    double sigRate;
    double bkgRate;
    double liveTime;
    std::vector<std::unique_ptr<SystematicEffect>> systematics;

  public:
    /// @brief Parametrised constructor
    /// @param r Signal rate in the last eV of the spectrum (counts s^-1)
    /// @param b Signal rate in the last eV of the spectrum (counts eV^-1 s^-1)
    /// @param t Live time of the experiment in seconds
    /// @param systs Vector of pointers to systematic effects
    Experiment(double r, double b, double t,
               std::vector<SystematicEffect*> systs);

    /// @brief Destructor
    ~Experiment();

    /// @brief Computes the statistical uncertainty on the square of the nu mass
    /// @return Statistical uncertainty on effective nu mass in eV^2
    double ComputeSigmaMSq();

    /// @brief Computes the 90% CL on the effective nu mass
    /// @return 90% CL on effective nu mass in eV
    double Compute90PcCL();
  };
}

#endif
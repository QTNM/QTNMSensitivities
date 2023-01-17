#ifndef EXPERIMENT_H
#define EXPERIMENT_H

namespace sens
{
  class Experiment
  {
  private:
    double sigRate;
    double bkgRate;
    double liveTime;
    
  public:
    /// @brief Parametrised constructor
    /// @param r Signal rate in the last eV of the spectrum (counts s^-1)
    /// @param b Signal rate in the last eV of the spectrum (counts eV^-1 s^-1)
    /// @param t Live time of the experiment in seconds
    Experiment(double r, double b, double t) : sigRate(r), bkgRate(b), liveTime(t) {}

    /// @brief Computes the statistical uncertainty on the square of the nu mass
    /// @return Statistical uncertainty on effective nu mass in eV^2
    double ComputeSigmaMSq();

    /// @brief Computes the 90% CL on the effective nu mass
    /// @return 90% CL on effective nu mass in eV
    double Compute90PcCL();
  };
}

#endif
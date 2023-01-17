#ifndef EXPERIMENT_H
#define EXPERIMENT_H

namespace sens
{
  class Experiment
  {
  private:
    double sigRate;
    double bkgRate;
    
  public:
    /// @brief Parametrised constructor
    /// @param r Signal rate in the last eV of the spectrum (counts s^-1)
    /// @param b Signal rate in the last eV of the spectrum (counts eV^-1 s^-1)
    Experiment(double r, double b) : sigRate(r), bkgRate(b) {}

    /// @brief Computes the statistical uncertainty on the square of the nu mass
    /// @param liveTime Uptime of experiment in seconds
    /// @return Statistical uncertainty on effective nu mass in eV^2
    double ComputeSigmaMSq(double liveTime);

    /// @brief Computes the 90% CL on the effective nu mass
    /// @param liveTime Uptime of experiment in seconds
    /// @return 90% CL on effective nu mass in eV
    double Compute90PcCL(double liveTime);
  };
}

#endif
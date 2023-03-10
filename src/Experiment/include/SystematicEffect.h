/*
  SystematicEffect.h

  Basic representation of a systematic that is just an RMS with and an
  uncertainty on that RMS
*/

#ifndef SYSTEMATIC_EFFECT_H
#define SYSTEMATIC_EFFECT_H

#include <string>

namespace sens
{
  class SystematicEffect
  {
  private:
    double rmsWidth; // RMS width in eV
    double widthUnc; // Uncertainity on RMS in eV
    std::string systName; // Systematic name

  protected:
    /// @brief RMS setter
    /// @param rms 
    void SetRMSWidth(double rms);

    /// @brief Width uncertainty setter
    /// @param unc 
    void SetWidthUnc(double unc);

    /// @brief Systematic name setter
    /// @param name 
    void SetSystName(std::string name);

  public:
    /// @brief Parametrised constructor
    /// @param width RMS width in eV
    /// @param unc Uncertainty on RMS in eV
    SystematicEffect(double width = 0, double unc = 0, std::string name = "")
        : rmsWidth(width), widthUnc(unc), systName(name) {}

    /// @brief Get the RMS width
    /// @return The RMS width in eV
    double GetWidth() { return rmsWidth; }

    /// @brief Get the uncertainty on the width
    /// @return Uncertainty on RMS in eV
    double GetUnc() { return widthUnc; }

    /// @brief Get the systematic name
    /// @return A string containing the systematic name
    std::string GetName() { return systName; }
  };

  /// Systematic representing effect of molecular tritium
  class MolecularTritium : public SystematicEffect
  {
  public:
    /// @brief Constructor
    MolecularTritium();
  };

  /// Systematic representing lack of uncertainty about magnetic field
  class BFieldUncertainty : public SystematicEffect
  {
  public:
    /// @brief Parametrised constructor
    /// @param fieldVar RMS variation in magnetic field in tesla
    /// @param bAvg Average magnetic field (in Tesla)
    BFieldUncertainty(double fieldVar, double bAvg);
  };

  /// Systematic uncertainty representing scattering from gas atoms/molecules
  class GasScattering : public SystematicEffect
  {
  public:
    /// @brief Parametrised constructor
    /// @param n Number density of atoms/molecules in m^-3
    /// @param b Magnetic field in tesla
    /// @param isAtomicTritium Boolean indicating if atomic tritium is used
    GasScattering(double n, double b, bool isAtomicTritium);
  };

  /// Systematic uncertainty representing doppler broadening from thermal energy
  class ThermalBroadening : public SystematicEffect
  {
  public:
    /// @brief Parametrised constructor
    /// @param T Temperature of atoms/molecules in kelvin
    /// @param isAtomicTritium Boolean indicating if atomic tritium is used
    ThermalBroadening(double T, bool isAtomicTritium);
  };
}

#endif
/*
  SystematicEffect.h

  Basic representation of a systematic that is just an RMS with and an
  uncertainty on that RMS
*/

namespace sens
{
  class SystematicEffect
  {
  private:
    double rmsWidth; // RMS width in eV
    double widthUnc; // Uncertainity on RMS in eV

  public:
    /// @brief Parametrised constructor
    /// @param width RMS width in eV
    /// @param unc Uncertainty on RMS in eV
    SystematicEffect(double width, double unc) : rmsWidth(width), widthUnc(unc) {}

    /// @brief Get the RMS width
    /// @return The RMS width in eV
    double GetWidth() { return rmsWidth; }

    /// @brief Get the uncertainty on the width
    /// @return Uncertainty on RMS in eV
    double GetUnc() { return widthUnc; }
  };
}
/*
  Formatting.h
*/

#include "TGraph.h"

namespace sens
{
  /// @brief 
  /// @param gr 
  void SetGraphAttr(TGraph *gr);

  /// @brief Formats TGraphs to with desired label sizes
  /// @param gr Unique pointer to TGraph
  void SetGraphAttr(std::unique_ptr<TGraph> &gr);
}
/// Formatting.cpp

#include "Formatting.h"

#include "TAxis.h"

void sens::SetGraphAttr(TGraph *gr)
{
  gr->SetLineWidth(2);
  gr->SetLineColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
}

void sens::SetGraphAttr(std::unique_ptr<TGraph> &gr)
{
  gr->SetLineWidth(2);
  gr->SetLineColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
}

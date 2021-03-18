#include "TROOT.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include <iostream>

//TF1 * fsb = NULL;
//#TF1 * fsr = NULL;

//Double_t function_divide(Double_t *x, Double_t *par)
//{
//  const Double_t xx = x[0];
//  return fsr->Eval(xx)/fsb->Eval(xx);
//}

Double_t landauCustom(Double_t *x, Double_t *par) {
  return TMath::Landau(*x,par[0],par[1]);
}

TF1 * landauFit(TH1D *hist) {
  TF1 *lfit = new TF1("testlfit","TMath::Landau(x,[0],[1])",900,5000);
  //lfit->SetParameter(0,hist->GetMean());
  lfit->SetParameter(0,1592);
  //lfit->SetParameter(1,hist->GetRMS());
  lfit->SetParameter(1,180);
  hist->Fit("testlfit","R");
  TF1 * testfit = hist->GetFunction("testlfit");
  std::cout<<lfit->GetParameters()<<std::endl;
  return (testfit);
}

TH1D *histDoubler(TH1D *hist){
  hist->Scale(2);
  return hist;
}



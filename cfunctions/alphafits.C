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

// this fucntion take 6 parameters 0..2 describe landau1
// 3..5 describe landau2
// returns landau1(x) / landau2(x)
Double_t landauRatio(Double_t *x, Double_t *par) {
  double numer = landauCustom(x,par);
  double denom = landauCustom(x,&par[3]);
  if (denom==0) return 1.0;  // project against divide by 0
  return numer/denom;	  
}

// implement a regular C++ function for maximum flexibility
Double_t landauModel(Double_t *X, Double_t *par){
  double x =X[0];  // get ring of pointer/array notation
  return par[0] * TMath::Landau(x,par[1],par[2]);
}

TF1 * landauFit(TH1D *hist) {
  TF1 *lfit = new TF1("testl",landauModel,900,5000,3);//options: low range, high range, num param

  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double rms = hist->GetRMS();
  lfit->SetParameter(0,amp);
  lfit->SetParameter(1,max);
  lfit->SetParameter(2,rms);
  hist->Fit("testl","LR+");
  hist->Draw();
  TF1* fitout = hist->GetFunction("testl");
  return fitout;
}

TF1 * landauFitExplicit(TH1D *hist) {
  TF1 *lfit = new TF1("testlfit","[0]*TMath::Landau(x,[1],[2])",500,5000);
  //lfit->SetParameter(0,hist->GetMean());
  lfit->SetParameter(0,30);
  //lfit->SetParameter(1,hist->GetRMS());
  lfit->SetParameter(1,1592);
  lfit->SetParameter(1,250);
  hist->Fit("testlfit","R0");
  TF1 * testfit = hist->GetFunction("testlfit");
  std::cout<<lfit->GetParameters()<<std::endl;
  return (testfit);
}

TH1D *histDoubler(TH1D *hist){
  hist->Scale(2);
  return hist;
}

TF1 * alphaRatioMaker(TH1D *histsr,TH1D *histsb) {
  TF1 *srfit = new TF1("srfit","[0]*TMath::Landau(x,[1],[2])",500,5000);
  TF1 *sbfit = new TF1("sbfit","[0]*TMath::Landau(x,[1],[2])",500,5000);
  srfit->SetParameter(0,30);
  srfit->SetParameter(1,1592);
  srfit->SetParameter(1,250);
  histsr->Fit("srfit","R0");
  TF1 * testfitsr = histsr->GetFunction("srfit");
  sbfit->SetParameter(0,34);
  sbfit->SetParameter(1,1598);
  sbfit->SetParameter(1,250);
  histsb->Fit("sbfit","R0");
  TF1 * testfitsb = histsb->GetFunction("sbfit");
  TF1 * alpha = new TF1("alpha","srfit/sbfit");
  return alpha;
}



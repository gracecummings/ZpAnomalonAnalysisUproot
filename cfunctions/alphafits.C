#include "TROOT.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include <iostream>
#include <cstring>

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

// implement a regular C++ function for maximum flexibility
Double_t landauModel(Double_t *X, Double_t *par){
  double x =X[0];  // get ring of pointer/array notation
  return par[0] * TMath::Landau(x,par[1],par[2]);
}

// this fucntion take 6 parameters 0..2 describe landau1
// 3..5 describe landau2
// returns landau1(x) / landau2(x)
Double_t landauRatio(Double_t *x, Double_t *par) {
  double numer = landauModel(x,par);
  double denom = landauModel(x,&par[3]);
  if (denom==0) return 1.0;  // project against divide by 0
  return numer/denom;	  
}

TF1 * landauFit(TH1D *hist, char *name) {
  TF1 *lfit = new TF1(name,landauModel,900,5000,3);//options: low range, high range, num param
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double rms = hist->GetRMS();
  lfit->SetParameter(0,amp);
  lfit->SetParameter(1,max);
  lfit->SetParameter(2,rms);
  double samp = lfit->GetParameter(0);
  hist->Fit(name,"LR0+");
  //hist->Draw();
  TF1* fitout = hist->GetFunction(name);
  double acamp = fitout->GetParameter(0);
  //std::cout<<" guessed amp: "<<amp<<std::endl;
  //std::cout<<"     set amp: "<<samp<<std::endl;
  //std::cout<<"post fit amp: "<<acamp<<std::endl;
  return fitout;
}

TF1 * alphaRatioMaker(TH1D *hsb, TH1D *hsr){
  string sb = "sbl";
  string sr = "srl";
  int len = sb.length();
  char sbl[len+1];
  char srl[len+1];
  strcpy(sbl,sb.c_str());
  strcpy(srl,sr.c_str());
  TF1 *sbfit= landauFit(hsb,sbl);
  TF1 *srfit= landauFit(hsr,srl);
  Double_t sboff = sbfit->GetParameter(0);
  Double_t sbmpv = sbfit->GetParameter(1);
  Double_t sbsig = sbfit->GetParameter(2);
  Double_t sroff = srfit->GetParameter(0);
  Double_t srmpv = srfit->GetParameter(1);
  Double_t srsig = srfit->GetParameter(2);
  TF1 *alpha = new TF1("alpha",landauRatio,900,5000,6);
  alpha->SetParameter(0,sroff);
  alpha->SetParameter(1,srmpv);
  alpha->SetParameter(2,srsig);
  alpha->SetParameter(3,sboff);
  alpha->SetParameter(4,sbmpv);
  alpha->SetParameter(5,sbsig);
 
  //Double_t landau1mpv = alpha->GetParameter(1);
  //std::cout<<"Set Most Probable Value: "<<srmpv<<std::endl;
  //std::cout<<"Read Most Probable Valu: "<<landau1mpv<<std::endl;
  //alpha->Draw();
  return alpha;
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






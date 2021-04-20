#include "TROOT.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <cstring>

Double_t landauCustom(Double_t *x, Double_t *par) {
  return TMath::Landau(*x,par[0],par[1]);
}

// implement a regular C++ function for maximum flexibility
Double_t landauModel(Double_t *X, Double_t *par){
  double x =X[0];  // get ring of pointer/array notation
  return par[0] * TMath::Landau(x,par[1],par[2]);
}

Double_t guessDecayConstant(TH1D *hist,double max) {
  double halfmax  = max/2;
  int halflifebin = hist->FindLastBinAbove(halfmax);
  Double_t halflife = hist->GetBinCenter(halflifebin);
  Double_t guess = TMath::Log(0.5)/halflife;
  return guess;
}

Double_t expModel(Double_t *X, Double_t *par){
  double x = X[0];
  //Double_t arg= 0;
  //if (par[2] != 0) arg = (x[0] - par[1]/par[2]);
  Double_t fitval = par[0]*TMath::Exp(par[1]*x);
  return fitval;
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

Double_t expRatio(Double_t *x, Double_t *par) {
  double numer = expModel(x,par);
  double denom = expModel(x,&par[2]);
  if (denom==0) return -1.0;
  return numer/denom;
}

//TF1 * landauFit(TH1D *hist, char *name) {
TF1 * landauFit(TH1D *hist,TString name, TString opt="LR0+") {
  TF1 *lfit = new TF1(name,landauModel,900,5000,3);//options: low range, high range, num param
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double rms = hist->GetRMS();
  lfit->SetParameter(0,amp);
  lfit->SetParameter(1,max);
  lfit->SetParameter(2,rms);
  double samp = lfit->GetParameter(0);
  //hist->Fit(name,"LR0+");
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  double acamp = fitout->GetParameter(0);
  //std::cout<<" guessed amp: "<<amp<<std::endl;
  //std::cout<<"     set amp: "<<samp<<std::endl;
  //std::cout<<"post fit amp: "<<acamp<<std::endl;
  return fitout;
}

TF1 * alphaRatioMakerLandau(TH1D *hsb, TH1D *hsr){
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

TF1 * expFit(TH1D *hist, TString name, TString opt="LR0+") {
  TF1 *expfit = new TF1(name,expModel,1500,3000,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambda = guessDecayConstant(hist,amp);
  expfit->SetParameter(0,amp);
  expfit->SetParameter(1,lambda);
  //double samp = lfit->GetParameter(0);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
  //double acamp = fitout->GetParameter(0);
}

TH1D * expFitErrBands(TH1D *hist, TString name, TString opt="LR0+") {
  TH1D *empty = new TH1D("empty","Attempt at error bands",100,500,5000);
  TF1 *expfit = new TF1(name,expModel,1500,3000,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambda = guessDecayConstant(hist,amp);
  expfit->SetParameter(0,amp);
  expfit->SetParameter(1,lambda);
  //double samp = lfit->GetParameter(0);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  //above is the old fit
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(empty);
  
  return empty;
  //double acamp = fitout->GetParameter(0);
}


TF1 * alphaRatioMakerExp(TH1D *hsb, TH1D *hsr){
  string sb = "sbl";
  string sr = "srl";
  int len = sb.length();
  char sbl[len+1];
  char srl[len+1];
  strcpy(sbl,sb.c_str());
  strcpy(srl,sr.c_str());
  TF1 *sbfit= expFit(hsb,sbl,"R0+");
  TF1 *srfit= expFit(hsr,srl,"R0+");
  Double_t sbamp = sbfit->GetParameter(0);
  Double_t sblambda = sbfit->GetParameter(1);
  Double_t sramp = srfit->GetParameter(0);
  Double_t srlambda = srfit->GetParameter(1);
  TF1 *alpha = new TF1("alpha",expRatio,1500,3000,4);
  alpha->SetParameter(0,sramp);
  alpha->SetParameter(1,srlambda);
  alpha->SetParameter(2,sbamp);
  alpha->SetParameter(3,sblambda);
  return alpha;
}

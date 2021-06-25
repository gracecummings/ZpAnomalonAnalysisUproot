#include "Math/WrappedTF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <cstring>
#include <cmath>

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
  //Double_t fitval = par[0]*TMath::Exp(par[1]*x);//old way, works great
  //Double_t fitval = TMath::Exp(TMath::Log(par[0])+par[1]*x);//same as above
  Double_t fitval = TMath::Exp(par[0]+par[1]*x);//new way, but keeping params closer
  return fitval;
}

Double_t expNModel(Double_t *X, Double_t *par){
  double x = X[0];
  Double_t fitval = par[0]*TMath::Exp(par[1]*x+par[2]/x);
  return fitval;
}

Double_t expRatio(Double_t *x, Double_t *par) {
  double numer = expModel(x,par);
  double denom = expModel(x,&par[2]);
  if (denom==0) return -1.0;
  return numer/denom;
}

Double_t expNRatio(Double_t *x, Double_t *par) {
  double numer = expNModel(x,par);
  double denom = expNModel(x,&par[2]);
  if (denom==0) return -1.0;
  return numer/denom;
}


TF1 * expFit(TH1D *hist, TString name, TString opt="R0+",int lowr=1500, int highr=3000) {
  TF1 *expfit = new TF1(name,expModel,lowr,highr,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambda = guessDecayConstant(hist,amp);
  expfit->SetParameter(0,TMath::Log(amp));//just use amp for old way
  //expfit->SetParameter(0,amp);//just use amp for old way
  expfit->SetParameter(1,lambda);
  //double samp = lfit->GetParameter(0);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
  //double acamp = fitout->GetParameter(0);
}

TF1 * expNFit(TH1D *hist, TString name, TString opt="R0+",int lowr=1500, int highr=5000) {
  TF1 *expNfit = new TF1(name,expNModel,1500,5000,3);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambda = guessDecayConstant(hist,amp);
  expNfit->SetParameter(0,amp);
  expNfit->SetParameter(1,lambda);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

// here we can do error propagation for an arbitrary function
// matrix and vector notation is used
// return error on function evaluated at x, sigma_f(x;par)
double GetError(ROOT::Math::WrappedTF1 &f, Double_t x, Double_t *pars, TMatrixD &COV, Int_t npar){
  Double_t *grad=new Double_t[npar];
  f.ParameterGradient(x,pars,grad);
  TVectorD g(npar,grad);
  double var = g*(COV*g);
  delete[] grad;
  return TMath::Sqrt(var);
}

TH1D * expFitErrBands(TH1D *hist, TString name, TString opt="R0+",Int_t nsig=2,int lowr=1500, int highr=3000) {
  const Int_t nPars=2;
  Double_t pars[nPars], grad[nPars];
  int nbins  = hist->GetNbinsX();
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambda = guessDecayConstant(hist,amp);
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  TH1D *errhist = new TH1D("errhist","Histogram with 2 sigma band for fit",nbins,histlowed,histhied);
  TF1 *expfit = new TF1(name,expModel,lowr,highr,nPars);

  expfit->SetParameter(0,TMath::Log(amp));
  expfit->SetParameter(1,lambda);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  fitout->GetParameters(pars);

  ROOT::Math::WrappedTF1 wfit(*fitout);

  //std::cout<<"Bin Width "<<binwidth<<std::endl;
  //std::cout<<"Bin number "<<hist->GetNbinsX()<<std::endl;
  //std::cout<<"First bin low edge "<<hist->GetBinLowEdge(1)<<std::endl;
  //std::cout<<"Last  bin low edge "<<hist->GetBinLowEdge(hist->GetNbinsX())+binwidth<<std::endl;

  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if ((fitplace - fitplacer) < 0.5) {
      fitstartbin = fitplacer+1;
    }
  else {
    fitstartbin = fitplacer;
  }
    
  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  int finalbin = (highr-histlowed)/binwidth;
  for (int ib = fitstartbin;ib<=finalbin;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    errhist->SetBinContent(ib,(*fitout)(x));
    errhist->SetBinError(ib,nsig*sigma);
  }

  //debug errors to show values are the same
  //for (int ib=12;ib<=finalbin;ib++){
  //double x = hist->GetBinCenter(ib);
  // double fitval = (*fitout)(x);
  //double sigma = GetError(wfit,x,pars,*COV,nPars);
  //std::cout<<"sigma on the fit point is "<<2*sigma<<std::endl;
  //std::cout<<"             fit point is "<<fitval<<std::endl;
  //std::cout<<"          err hist val is "<<errhist->GetBinContent(ib)<<std::endl;
  //std::cout<<"          err hist err is "<<errhist->GetBinError(ib)<<std::endl;
  //}

  return errhist;//returns a histogram with the fit values as the bins and the 1 sigma band the error
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
  TF1 *alpha = new TF1("alpha",expRatio,1500,5000,4);
  alpha->SetParameter(0,sramp);
  alpha->SetParameter(1,srlambda);
  alpha->SetParameter(2,sbamp);
  alpha->SetParameter(3,sblambda);
  return alpha;
}


///////Development for normalization//

Double_t poly5Model(Double_t *X, Double_t *par){
  double x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0)+par[5]*std::pow(x,5);
  return fitval;
}

TF1 * poly5Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5fit = new TF1(name,poly5Model,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

///////Unused Functions/////

/*
// this fucntion take 6 parameters 0..2 describe landau1
// 3..5 describe landau2
// returns landau1(x) / landau2(x)
Double_t landauRatio(Double_t *x, Double_t *par) {
  double numer = landauModel(x,par);
  double denom = landauModel(x,&par[3]);
  if (denom==0) return 1.0;  // project against divide by 0
  return numer/denom;	  
}

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
*/


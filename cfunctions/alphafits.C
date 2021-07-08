
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
  Double_t x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0)+par[5]*std::pow(x,5);
  return fitval;
}

TF1 * poly5Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5fit = new TF1(name,poly5Model,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

Double_t gausModel(Double_t *X, Double_t *par){
  double x = X[0];
  Double_t fitval = par[0]*TMath::Gaus(x,par[1],par[2]);
  return fitval;
}

TF1 * gausFit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400,int mguess=90,int sigguess = 5) {
  TF1 *gausfit = new TF1(name,gausModel,lowr,highr,5);
  gausfit->SetParameter(1,mguess);
  gausfit->SetParameter(2,sigguess);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

Double_t gaus2Model2(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[5]*(par[0]*TMath::Gaus(x,par[1],par[2])+(1-par[0])*TMath::Gaus(x,par[3],par[4]));
  return fitval;
}

TF1 * gaus2Fit2(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

Double_t gaus2Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]*TMath::Gaus(x,par[1],par[2])+(par[5])*TMath::Gaus(x,par[3],par[4]);
  return fitval;
}

TF1 * gaus2Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *gaus2fit = new TF1(name,gaus2Model,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

Double_t ErfExp(Double_t x, Double_t c, Double_t offset, Double_t width){
    if(width<1e-2)width=1e-2;
    if (c==0)c=-1e-7;
	return TMath::Exp(c*x)*(1.+TMath::Erf((x-offset)/width))/2. ;
}

Double_t gausErfExpModel(Double_t *X, Double_t *par){
  double x = X[0];
  Double_t fitval = par[0]*TMath::Gaus(x,par[1],par[2])+par[3]*ErfExp(x,par[4],par[5],par[6]);
  return fitval;
}

Double_t gausPoly1Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]*TMath::Gaus(x,par[1],par[2])+par[3]*x+par[4];
  return fitval;
}


TF1 * gausErfExpFit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400,int mguess=90,int sigguess = 5) {
  TF1 *gausErffit = new TF1(name,gausErfExpModel,lowr,highr,7);
  gausErffit->SetParameter(1,mguess);
  gausErffit->SetParameter(2,sigguess);
  //gausErffit->SetParameter(3,30);
  //gausErffit->SetParameter(4,30);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * gausPoly1Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400,int mguess=90,int sigguess = 5) {
  TF1 *gausPoly1fit = new TF1(name,gausPoly1Model,lowr,highr,5);
  gausPoly1fit->SetParameter(1,mguess);
  gausPoly1fit->SetParameter(2,sigguess);
  gausPoly1fit->SetParameter(5,hist->GetBinContent(6));
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

Double_t totalBkgModel(Double_t *X, Double_t *par){
  //Double_t x = X[0];
  //Double_t fitval = par[17]*poly5Model(X,&par[0])+par[18]*gaus2Model(X,&par[6])+par[19]*gausPoly1Model(X,&par[12]);
  Double_t fitval = par[17]*poly5Model(X,&par[0])+gaus2Model2(X,&par[6])+par[18]*gausPoly1Model(X,&par[12]);
  return fitval;
}

Double_t totalBkgModelBlind(Double_t *X, Double_t *par){
  if (X[0] > 70. && X[0] < 150.) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t fitval = poly5Model(X,&par[0])+gaus2Model(X,&par[6])+gausPoly1Model(X,&par[12]);
  return fitval;
}

vector<TF1 *> totalFit(TH1D *hist, TH1D *dyhist, TH1D *tthist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *dyfit = poly5Fit(dyhist,"dyl","R0+",30,250);
  TF1 *ttfit = gaus2Fit2(tthist,"ttl","R0+",30,400);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","R0+",30,250);

  Double_t par[17];
  dyfit->GetParameters(&par[0]);
  ttfit->GetParameters(&par[6]);
  vvfit->GetParameters(&par[12]);

  TF1 *totalfit = new TF1("totalfit",totalBkgModel,30,250,19);//17 without norms
  totalfit->FixParameter(0,par[0]);
  totalfit->FixParameter(1,par[1]);
  totalfit->FixParameter(2,par[2]);
  totalfit->FixParameter(3,par[3]);
  totalfit->FixParameter(4,par[4]);
  totalfit->FixParameter(5,par[5]);
  totalfit->FixParameter(6,par[6]);
  totalfit->FixParameter(7,par[7]);
  totalfit->FixParameter(8,par[8]);
  totalfit->FixParameter(9,par[9]);
  totalfit->FixParameter(10,par[10]);
  totalfit->FixParameter(11,par[11]);//norm of ttbar
  totalfit->FixParameter(12,par[12]);
  totalfit->FixParameter(13,par[13]);
  totalfit->FixParameter(14,par[14]);
  totalfit->FixParameter(15,par[15]);
  totalfit->FixParameter(16,par[16]);
  //totalfit->SetParameters(17,1);
  //totalfit->SetParameters(18,1);
  //totalfit->SetParameters(19,1);
  hist->Fit("totalfit","R0L+");
  TF1 *totmcfit = hist->GetFunction("totalfit");

  Double_t parData[17];
  totmcfit->GetParameters(parData);
  //std::cout<<"Z peak guess "<<parData[13]<<std::endl;
  //std::cout<<"t peak guess "<<parData[9]<<std::endl;

  //Now do the extrapolation
  TF1 *sbdatfit = new TF1("sbdatfit",totalBkgModelBlind,30,250,17);
  sbdatfit->SetParameters(parData);
  sbdatfit->FixParameter(1,parData[1]);
  sbdatfit->FixParameter(2,parData[2]);
  sbdatfit->FixParameter(3,parData[3]);
  sbdatfit->FixParameter(4,parData[4]);
  sbdatfit->FixParameter(5,parData[5]);
  sbdatfit->FixParameter(7,parData[7]);
  sbdatfit->FixParameter(8,parData[8]);
  sbdatfit->FixParameter(9,parData[9]);
  sbdatfit->FixParameter(10,parData[10]);
  sbdatfit->FixParameter(12,parData[12]);
  sbdatfit->FixParameter(13,parData[13]);
  sbdatfit->FixParameter(14,parData[14]);
  sbdatfit->FixParameter(15,parData[15]);

  dathist->Fit("sbdatfit","LR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfit");

  Double_t parExtrap[17];
  totsbdatfit->GetParameters(parExtrap);
  TF1 *totextrap = new TF1("totextrap",totalBkgModel,30,250,17);
  totextrap->FixParameter(0,parExtrap[0]);
  totextrap->FixParameter(1,parExtrap[1]);
  totextrap->FixParameter(2,parExtrap[2]);
  totextrap->FixParameter(3,parExtrap[3]);
  totextrap->FixParameter(4,parExtrap[4]);
  totextrap->FixParameter(5,parExtrap[5]);
  totextrap->FixParameter(6,parExtrap[6]);
  totextrap->FixParameter(7,parExtrap[7]);
  totextrap->FixParameter(8,parExtrap[8]);
  totextrap->FixParameter(9,parExtrap[9]);
  totextrap->FixParameter(10,parExtrap[10]);
  totextrap->FixParameter(11,parExtrap[11]);
  totextrap->FixParameter(12,parExtrap[12]);
  totextrap->FixParameter(13,parExtrap[13]);
  totextrap->FixParameter(14,parExtrap[14]);
  totextrap->FixParameter(15,parExtrap[15]);
  totextrap->FixParameter(16,parExtrap[16]);
  dathist->Fit("totextrap","QR0+");
  TF1 *totnorm = dathist->GetFunction("totextrap");

  //Now seperate fits for visualization
  //broken
  TF1 *flsb = new TF1("flsb",totalBkgModelBlind,30,70,17);
  flsb->FixParameter(0,parExtrap[0]);
  flsb->FixParameter(1,parExtrap[1]);
  flsb->FixParameter(2,parExtrap[2]);
  flsb->FixParameter(3,parExtrap[3]);
  flsb->FixParameter(4,parExtrap[4]);
  flsb->FixParameter(5,parExtrap[5]);
  flsb->FixParameter(6,parExtrap[6]);
  flsb->FixParameter(7,parExtrap[7]);
  flsb->FixParameter(8,parExtrap[8]);
  flsb->FixParameter(9,parExtrap[9]);
  flsb->FixParameter(10,parExtrap[10]);
  flsb->FixParameter(11,parExtrap[11]);
  flsb->FixParameter(12,parExtrap[12]);
  flsb->FixParameter(13,parExtrap[13]);
  flsb->FixParameter(14,parExtrap[14]);
  flsb->FixParameter(15,parExtrap[15]);
  flsb->FixParameter(16,parExtrap[16]);
  dathist->Fit("flsb","QR0+");
  TF1 *lsbdatfit = dathist->GetFunction("flsb");

  TF1 *fhsb = new TF1("fhsb",totalBkgModelBlind,150,250,17);
  fhsb->FixParameter(0,parExtrap[0]);
  fhsb->FixParameter(1,parExtrap[1]);
  fhsb->FixParameter(2,parExtrap[2]);
  fhsb->FixParameter(3,parExtrap[3]);
  fhsb->FixParameter(4,parExtrap[4]);
  fhsb->FixParameter(5,parExtrap[5]);
  fhsb->FixParameter(6,parExtrap[6]);
  fhsb->FixParameter(7,parExtrap[7]);
  fhsb->FixParameter(8,parExtrap[8]);
  fhsb->FixParameter(9,parExtrap[9]);
  fhsb->FixParameter(10,parExtrap[10]);
  fhsb->FixParameter(11,parExtrap[11]);
  fhsb->FixParameter(12,parExtrap[12]);
  fhsb->FixParameter(13,parExtrap[13]);
  fhsb->FixParameter(14,parExtrap[14]);
  fhsb->FixParameter(15,parExtrap[15]);
  fhsb->FixParameter(16,parExtrap[16]);
  dathist->Fit("fhsb","QR0+");
  TF1 *hsbdatfit = dathist->GetFunction("fhsb");


  
  std::vector<TF1*> fitvector;
  fitvector.push_back(totmcfit);
  fitvector.push_back(totsbdatfit);
  fitvector.push_back(totnorm);
  fitvector.push_back(lsbdatfit);
  fitvector.push_back(hsbdatfit);
  return fitvector;
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


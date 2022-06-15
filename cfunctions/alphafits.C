#include "TFitResult.h"
#include "Math/WrappedTF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector.h"
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

Double_t guessDecayConstantOG(TH1D *hist,double max) {
  double halfmax  = max/2;
  int halflifebin = hist->FindLastBinAbove(halfmax);
  Double_t halflife = hist->GetBinCenter(halflifebin);
  Double_t guess = TMath::Log(0.5)/halflife;
  return guess;
}

Double_t guessDecayConstant(TH1D *hist,double max) {
  int maxbin = hist->GetMaximumBin();
  int intbin = maxbin+4;
  int xval = hist->GetBinCenter(intbin);
  Double_t bincont = hist->GetBinContent(intbin);
  Double_t guess = (TMath::Log(bincont)-TMath::Log(max))/xval;
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

Double_t expRatioMultiply(Double_t *x, Double_t *par) {
  double numer = expModel(x,par);
  double denom = expModel(x,&par[2]);
  double multiplier = expModel(x,&par[4]);
  if (denom==0) return -1.0;
  return numer/denom*multiplier;
}

TF1 * expFit(TH1D *hist, TString name, TString opt="R0+",int lowr=1500, int highr=3000) {
  TF1 *expfit = new TF1(name,expModel,lowr,highr,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambdaOG = guessDecayConstantOG(hist,amp);
  double_t lambda = guessDecayConstant(hist,amp);

  std::cout<<" Max guess> "<<amp<<std::endl;
  std::cout<<" lamda guess OG > "<<lambdaOG<<std::endl;
  std::cout<<" lamda guess> "<<lambda<<std::endl;
  expfit->SetParameter(0,TMath::Log(amp));//just use amp for old way
  //expfit->SetParameter(0,amp);//just use amp for old way
  expfit->SetParameter(1,lambda);
  //double samp = lfit->GetParameter(0);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
  //double acamp = fitout->GetParameter(0);
}

TF1 * expFitSetParsAndErrs(TString name,TVector pars, int lowr,int highr){
  TF1 *expfit = new TF1(name,expModel,lowr,highr,2);
  expfit->SetParameter(0,pars[0]);
  expfit->SetParameter(1,pars[1]);
  return expfit;
}

TVectorD expFitDecorrParamsShiftedUp(TH1D *hist, TString name, TString opt="R0+",int lowr=1500, int highr=3000) {
  //Do the basic fit
  TF1 *expfit = new TF1(name,expModel,lowr,highr,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambdaOG = guessDecayConstantOG(hist,amp);
  double_t lambda = guessDecayConstant(hist,amp);
  expfit->SetParameter(0,TMath::Log(amp));//just use amp for old way
  expfit->SetParameter(1,lambda);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  //Make the stuff to hold the shifted params
  const Int_t nPars = 2;
  Double_t pars[nPars], grad[nPars];
  TVectorD oriparams = TVectorD(nPars);
  TVectorD vdecorrparams = TVectorD(nPars);
  TVectorD vdecorrerrs = TVectorD(nPars);
  TVectorD vecout = TVectorD(nPars);

  //shift the params
  fitout->GetParameters(pars);
  for (int i = 0;i<nPars;i++){
    oriparams[i] = pars[i];
  }
  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  TMatrixD eigvec = COV->EigenVectors(vdecorrerrs);
  eigvec.Transpose(eigvec);
  vdecorrparams = eigvec*oriparams;
  vecout = eigvec.Invert()*(vdecorrparams+vdecorrerrs);
  return vecout;
}

TVectorD expFitDecorrParamsShiftedDown(TH1D *hist, TString name, TString opt="R0+",int lowr=1500, int highr=3000) {
  //Do the basic fit
  TF1 *expfit = new TF1(name,expModel,lowr,highr,2);
  int binmax = hist->GetMaximumBin();
  double max = hist->GetXaxis()->GetBinCenter(binmax);
  double amp = hist->GetMaximum();
  double_t lambdaOG = guessDecayConstantOG(hist,amp);
  double_t lambda = guessDecayConstant(hist,amp);
  expfit->SetParameter(0,TMath::Log(amp));//just use amp for old way
  expfit->SetParameter(1,lambda);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  //Make the stuff to hold the shifted params
  const Int_t nPars = 2;
  Double_t pars[nPars], grad[nPars];
  TVectorD oriparams = TVectorD(nPars);
  TVectorD vdecorrparams = TVectorD(nPars);
  TVectorD vdecorrerrs = TVectorD(nPars);
  TVectorD vecout = TVectorD(nPars);

  //shift the params
  fitout->GetParameters(pars);
  for (int i = 0;i<nPars;i++){
    oriparams[i] = pars[i];
  }
  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  TMatrixD eigvec = COV->EigenVectors(vdecorrerrs);
  eigvec.Transpose(eigvec);
  vdecorrparams = eigvec*oriparams;
  vecout = eigvec.Invert()*(vdecorrparams-vdecorrerrs);
  return vecout;
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
  //std::cout<<"the variance "<<var<<std::endl;
  
  delete[] grad;
 return TMath::Sqrt(var);
}

double GetErrorDecorrelated(ROOT::Math::WrappedTF1 &f, Double_t x, Double_t *pars, TMatrixD &COV, Int_t npar){
  Double_t *grad=new Double_t[npar];
  f.ParameterGradient(x,pars,grad);
  TVectorD g(npar,grad);
  double var = g*(COV*g);
  //std::cout<<"the variance "<<var<<std::endl;
  
  delete[] grad;
 return TMath::Sqrt(var);
}

double GetErrorGEC(ROOT::Math::WrappedTF1 &f, Double_t x, Double_t *pars, TMatrixDSym &COV, Int_t npar){
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

  //std::cout<<"hist low edge  "<<histlowed<<std::endl;
  //std::cout<<"fit low edge   "<<lowr<<std::endl;
  
  float fitplace = (lowr-histlowed)/binwidth;
  //std::cout<<"fit place "<<fitplace<<std::endl;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  //example: hist low edge = 0, fit low = 1500, binwidth = 200
  //fitplace = 7.5, we actually know that 1500 is in bin 8
  //we want the fit histogram to start at bin 8
  //fit placer, when rounded, gives us 8
  if (std::abs(fitplace - fitplacer) < 0.5) {
      fitstartbin = fitplacer+1;
      //std::cout<<"calulated diff "<<(fitplace - fitplacer)<<std::endl;
    }
  else {
    fitstartbin = fitplacer;
  }
  //std::cout<<"rounded fit placer "<<fitplacer<<std::endl;
  //std::cout<<"fit start bin "<<fitstartbin<<std::endl;
  //std::cout<<"Value in starting bin "<<hist->GetBinContent(fitplacer)<<std::endl;
  //std::cout<<"Filling the UNC Bin histogram "<<std::endl;
    
  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  int finalbin = (highr-histlowed)/binwidth;
  for (int ib = fitstartbin;ib<=finalbin;ib++){
    //std::cout<<"Bin number "<<ib<<std::endl;
    double x = hist->GetBinCenter(ib);
    //std::cout<<"Bin Center "<<x<<std::endl;
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    //std::cout<<"Fit value  "<<(*fitout)(x)<<std::endl;
    //std::cout<<"Hist value  "<<hist->GetBinContent(ib)<<std::endl;
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


TF1 * alphaRatioMakerExp(TH1D *hsb, TH1D *hsr,int lowrsb=1500, int highrsb=3000,int lowrsr=1500, int highrsr=3000){
  string sb = "sbl";
  string sr = "srl";
  int len = sb.length();
  char sbl[len+1];
  char srl[len+1];
  strcpy(sbl,sb.c_str());
  strcpy(srl,sr.c_str());
  TF1 *sbfit= expFit(hsb,sbl,"R0+",lowrsb,highrsb);//added ranges
  TF1 *srfit= expFit(hsr,srl,"R0+",lowrsr,highrsr);
  Double_t sbamp = sbfit->GetParameter(0);
  Double_t sblambda = sbfit->GetParameter(1);
  Double_t sramp = srfit->GetParameter(0);
  Double_t srlambda = srfit->GetParameter(1);
  TF1 *alpha = new TF1("alpha",expRatio,1500,5000,4);//was to 5000 previously
  alpha->SetParameter(0,sramp);
  alpha->SetParameter(1,srlambda);
  alpha->SetParameter(2,sbamp);
  alpha->SetParameter(3,sblambda);
  return alpha;
}


TF1 * alphaRatioMakerExpExternalParams(TF1 *sbfit,TF1 *srfit, TString name) {
  Double_t sbamp = sbfit->GetParameter(0);
  Double_t sblambda = sbfit->GetParameter(1);
  Double_t sramp = srfit->GetParameter(0);
  Double_t srlambda = srfit->GetParameter(1);
  TF1 *alpha = new TF1(name,expRatio,1500,5000,4);//was to 5000 previously
  alpha->SetParameter(0,sramp);
  alpha->SetParameter(1,srlambda);
  alpha->SetParameter(2,sbamp);
  alpha->SetParameter(3,sblambda);
  return alpha;
}

//TF1 *datafitforsubtraction;
//TF1 *ttsbfitforsubtraction;
//TF1 *vvsbfitforsubtraction;

Double_t subtractedFitsForDataSB(Double_t *x,Double_t *par){
  Double_t data = expModel(x,par);
  Double_t tt   = expModel(x,&par[2]);
  Double_t vv   = expModel(x,&par[4]);
  return data-tt-vv;
}

TF1 * subtractionFromFits(TF1 *dataf, TF1 *ttf, TF1 *vvf,TString name){
  Double_t par[6];
  dataf->GetParameters(&par[0]);
  ttf->GetParameters(&par[2]);
  vvf->GetParameters(&par[4]);
  TF1 *subtracteddatasb = new TF1(name,subtractedFitsForDataSB,1500,3000,6);
  subtracteddatasb->SetParameter(0,par[0]);
  subtracteddatasb->SetParameter(1,par[1]);
  subtracteddatasb->SetParameter(2,par[2]);
  subtracteddatasb->SetParameter(3,par[3]);
  subtracteddatasb->SetParameter(4,par[4]);
  subtracteddatasb->SetParameter(5,par[5]);
  return subtracteddatasb;
}

Double_t alphaExtrapBuiltFromFunctions(Double_t *x,Double_t *par){
  Double_t subdat = subtractedFitsForDataSB(x,par);
  Double_t alpha  = expRatio(x,&par[6]);
  return subdat*alpha;
}

TF1 * alphaExtrapolationFromFits(TF1* subsb,TF1* alphar,TString name){
  Double_t par[10];
  subsb->GetParameters(&par[0]);
  alphar->GetParameters(&par[6]);
  TF1 *alphaextraptf1 = new TF1(name,alphaExtrapBuiltFromFunctions,1500,3000,10);
  //alphaextraptf1->SetParameters(par);
  alphaextraptf1->SetParameter(0,par[0]);
  alphaextraptf1->SetParameter(1,par[1]);
  alphaextraptf1->SetParameter(2,par[2]);
  alphaextraptf1->SetParameter(3,par[3]);
  alphaextraptf1->SetParameter(4,par[4]);
  alphaextraptf1->SetParameter(5,par[5]);
  alphaextraptf1->SetParameter(6,par[6]);
  alphaextraptf1->SetParameter(7,par[7]);
  alphaextraptf1->SetParameter(8,par[8]);
  alphaextraptf1->SetParameter(9,par[9]);
  return alphaextraptf1;
}

TF1 * alphaExtrapolation(TH1D *hsb, TH1D *hsr, TH1D *hdatsb,int lowrsb=1500, int highrsb=3000,int lowrsr=1500, int highrsr=3000,int lowrdat = 1500,int highrdat = 3000){
  string sb = "sbl";
  string sr = "srl";
  string dt = "dtl";
  int len = sb.length();
  char sbl[len+1];
  char srl[len+1];
  char dtl[len+1];
  strcpy(sbl,sb.c_str());
  strcpy(srl,sr.c_str());
  strcpy(dtl,dt.c_str());
  TF1 *sbfit= expFit(hsb,sbl,"R0+",lowrsb,highrsb);//added ranges
  TF1 *srfit= expFit(hsr,srl,"R0+",lowrsr,highrsr);
  //TF1 *sbdatfit= expFit(hdatsb,dtl,"R0+",lowrdat,highrdat);
  TF1 *sbdatfit= expFit(hdatsb,dtl,"R0+",lowrdat,3000);
  //TF1 *sbdatfit= expFit(hdatsb,dtl,"R0+",1500,2500);//jecup
  Double_t sbamp = sbfit->GetParameter(0);
  Double_t sblambda = sbfit->GetParameter(1);
  Double_t sramp = srfit->GetParameter(0);
  Double_t srlambda = srfit->GetParameter(1);
  Double_t sbdatamp = sbdatfit->GetParameter(0);
  Double_t sbdatlambda = sbdatfit->GetParameter(1);
  TF1 *srextrap = new TF1("srextrap",expRatioMultiply,1500,5000,6);
  srextrap->SetParameter(0,sramp);
  srextrap->SetParameter(1,srlambda);
  srextrap->SetParameter(2,sbamp);
  srextrap->SetParameter(3,sblambda);
  srextrap->SetParameter(4,sbdatamp);
  srextrap->SetParameter(5,sbdatlambda);
  return srextrap;
}

TH1D * alphaExtrapolationHist(TH1D *hsb, TH1D *hsr, TH1D *hdatsb,int rebindiv=1,int lowrsb=1500, int highrsb=3000,int lowrsr=1500, int highrsr=3000,int lowrdat = 1500,int highrdat = 3000){
  string sb = "sbl";
  string sr = "srl";
  string dt = "dtl";
  int len = sb.length();
  char sbl[len+1];
  char srl[len+1];
  char dtl[len+1];
  strcpy(sbl,sb.c_str());
  strcpy(srl,sr.c_str());
  strcpy(dtl,dt.c_str());
  float nbins     = hsb->GetNbinsX();
  float binwidth  = hsb->GetBinWidth(1);
  float histlowed = hsb->GetBinLowEdge(1);
  float histhied  = hsb->GetBinLowEdge(nbins)+binwidth;

  
  //TF1 *sbfit= expFit(hsb,sbl,"R0+",1500,5000);
  //TF1 *srfit= expFit(hsr,srl,"R0+",1500,4000);
  //TF1 *sbdatfit= expFit(hdatsb,dtl,"R0+",1500,3000);
  TF1 *sbfit= expFit(hsb,sbl,"R0+",lowrsb,highrsb);//added ranges
  TF1 *srfit= expFit(hsr,srl,"R0+",lowrsr,highrsr);

  ///limits might be weird because of the limits of the others
  TF1 *sbdatfit= expFit(hdatsb,dtl,"R0+",lowrdat,3000);///THIS IS THE BUG
  //This is a fit to a histogram that is made form subtracting fits from fits
  Double_t sbamp = sbfit->GetParameter(0);
  Double_t sblambda = sbfit->GetParameter(1);
  Double_t sramp = srfit->GetParameter(0);
  Double_t srlambda = srfit->GetParameter(1);
  Double_t sbdatamp = sbdatfit->GetParameter(0);
  Double_t sbdatlambda = sbdatfit->GetParameter(1);
  //limits!!! NEEDs attention
  TF1 *srextrap = new TF1("srextrap",expRatioMultiply,1500,5000,6);
  srextrap->SetParameter(0,sramp);
  srextrap->SetParameter(1,srlambda);
  srextrap->SetParameter(2,sbamp);
  srextrap->SetParameter(3,sblambda);
  srextrap->SetParameter(4,sbdatamp);
  srextrap->SetParameter(5,sbdatlambda);

  //Make the histogram
  TH1D *extrphist = new TH1D("extrphist","Histogram with DY extrap to SR",nbins/rebindiv,histlowed,histhied);
  for (int ib = 0; ib <= (nbins/rebindiv); ib++){
    double x = extrphist->GetBinCenter(ib);
    double val = (*srextrap)(x);
    extrphist->SetBinContent(ib,val);
  }

  return extrphist;
  
}


///////Development for normalization//

Double_t poly5Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0)+par[5]*std::pow(x,5);
  return fitval;
}

Double_t poly5Model4(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[1]*(x-par[0])+par[2]*std::pow(x-par[0],2.0)+par[3]*std::pow(x-par[0],3.0)+par[4]*std::pow(x-par[0],4.0)+par[5]*std::pow(x-par[0],5);
  return fitval;
}

Double_t poly5Model3(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]*x+par[1]*std::pow(x,2.0)+par[2]*std::pow(x,3.0)+par[3]*std::pow(x,4.0)+par[4]*std::pow(x,5);
  return fitval;
}

Double_t poly5Model5(Double_t *X, Double_t *par){
  Double_t x = X[0]-250;
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3)+par[4]*std::pow(x,4)+par[5]*std::pow(x,5);
  return fitval;
}

Double_t poly6Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0)+par[5]*std::pow(x,6)+par[5]*std::pow(x,6);
  return fitval;
}

Double_t poly4Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0);
  return fitval;
}

Double_t poly3Model(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0);
  return fitval;
}

TF1 * poly5Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5fit = new TF1(name,poly5Model,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly5mod3Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5mod3fit = new TF1(name,poly5Model3,lowr,highr,5);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly5mod5Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly5mod5FitSetParsAndErrs(TString name, TVector pars,int lowr=30, int highr=400) {
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  poly5mod5fit->SetParameter(0,pars[0]);
  poly5mod5fit->SetParameter(1,pars[1]);
  poly5mod5fit->SetParameter(2,pars[2]);
  poly5mod5fit->SetParameter(3,pars[3]);
  poly5mod5fit->SetParameter(4,pars[4]);
  poly5mod5fit->SetParameter(5,pars[5]);
  return poly5mod5fit;
}

TH1D * poly5mod5FitErrBands(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  TH1D *errhist = new TH1D("errhist","Histogram with 1 sigma band for fit",nbins,histlowed,histhied);

  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
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
    errhist->SetBinError(ib,sigma);
  }
  return errhist;
}

TGraph *errup;
TGraph *errdwn;

Double_t upfunc(Double_t *x,Double_t *t){return errup->Eval(x[0]);}
Double_t dwnfunc(Double_t *x,Double_t *t){return errdwn->Eval(x[0]);}

vector<TF1 *> poly5mod5FitErrFunctionsGraphs(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  //TH1D *errhistup = new TH1D("errhistup","Histogram with bin value as nomval+1sigma",nbins,histlowed,histhied);
  //TH1D *errhistdwn = new TH1D("errhistdwn","Histogram with bin value as nomval+1sigma",nbins,histlowed,histhied);

  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
    fitstartbin = fitplacer+1;
  }
  else {
    fitstartbin = fitplacer;
  }
  int finalbin = (highr-histlowed)/binwidth;

  errup = new TGraph(finalbin+1);
  errdwn = new TGraph(finalbin+1);
  
  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  for (int ib = fitstartbin;ib<=finalbin+1;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    errup->SetPoint(ib,x,(*fitout)(x)+sigma);
    errdwn->SetPoint(ib,x,(*fitout)(x)-sigma);
  }

  TF1 *fitoutup  = new TF1("fitoutup",upfunc,lowr,highr,0);
  TF1 *fitoutdwn = new TF1("fitoutdwn",dwnfunc,lowr,highr,0);

  //std::cout<<"DEBUG DEBUG DEBUG"<<std::endl;
  //for (int ib = fitstartbin;ib<=finalbin+1;ib++){
  //double x = hist->GetBinCenter(ib);
    //double sigma = GetError(wfit,x,pars,*COV,nPars);
    //double upgx;
    //double upgy;
    //errup->GetPoint(ib,upgx,upgy);
    //std::cout<<"The histogram up value y: "<<upgy<<std::endl;
    //std::cout<<"The       fit up value: "<<(*fitoutup)(x)<<std::endl;
    //std::cout<<"The histogram dn value: "<<errdwn->GetPoint(ib)<<std::endl;
    //std::cout<<"The       fit dn value: "<<(*fitoutdwn)(x)<<std::endl;
  //}

  //std::cout<<fitoutup<<std::endl;
  //std::cout<<fitoutdwn<<std::endl;

  //TFormula* upform = fitoutup->GetFormula();
  //std::cout<<"The formula param 0: "<<upform->GetParameter(0)<<std::endl;
  vector<TF1 *> vecout;
  vecout.push_back(fitoutup);
  vecout.push_back(fitoutdwn);
  return vecout;
}

vector<TF1 *> poly5mod5FitErrFunctionsHists(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  TH1D *errhistup = new TH1D("errhistup","Histogram with bin value as nomval+1sigma",nbins,histlowed,histhied);
  TH1D *errhistdwn = new TH1D("errhistdwn","Histogram with bin value as nomval+1sigma",nbins,histlowed,histhied);

  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
    fitstartbin = fitplacer+1;
  }
  else {
    fitstartbin = fitplacer;
  }
  int finalbin = (highr-histlowed)/binwidth;

  //errup = new TGraph(finalbin+1);
  //errdwn = new TGraph(finalbin+1);
  
  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  for (int ib = fitstartbin;ib<=finalbin+1;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    errhistup->SetBinContent(ib,x,(*fitout)(x)+sigma);
    errhistdwn->SetBinContent(ib,x,(*fitout)(x)-sigma);
  }

  TF1 *poly5mod5fitup = new TF1("fitforupuncs",poly5Model5,lowr,highr,6);
  TF1 *poly5mod5fitdn = new TF1("fitfordnuncs",poly5Model5,lowr,highr,6);
  errhistup->Fit("fitforupuncs","R0+");
  errhistdwn->Fit("fitfordnuncs","R0+");
  TF1 *fitoutup = errhistup->GetFunction("fitforupuncs");
  TF1 *fitoutdwn = errhistdwn->GetFunction("fitfordnuncs");

  std::cout<<"DEBUG DEBUG DEBUG"<<std::endl;
  for (int ib = fitstartbin;ib<=finalbin+1;ib++){
    double x = hist->GetBinCenter(ib);
    //double sigma = GetError(wfit,x,pars,*COV,nPars);
    std::cout<<"The histogram up value: "<<errhistup->GetBinContent(ib)<<std::endl;
    std::cout<<"The       fit up value: "<<(*fitoutup)(x)<<std::endl;
    std::cout<<"The histogram dn value: "<<errhistdwn->GetBinContent(ib)<<std::endl;
    std::cout<<"The       fit dn value: "<<(*fitoutdwn)(x)<<std::endl;
  }
    

  std::cout<<fitoutup<<std::endl;
  std::cout<<fitoutdwn<<std::endl;
  vector<TF1 *> vecout;
  vecout.push_back(fitoutup);
  vecout.push_back(fitoutdwn);
  return vecout;
}

TVectorD poly5mod5FitDecorrelatedUncs(TH1D *hist,TString name,TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  TVectorD vecout;
  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  //std::cout<<"The standard covariance matrix"<<std::endl;
  //COV->Print();
  TMatrixD eigvec = COV->EigenVectors(vecout);
  //std::cout<<"The eigenvators of the covariance matrix"<<std::endl;
  //eigvec.Print();
  //std::cout<<"The eigenvalues of the covariance matrix"<<std::endl;
  //vecout.Print();
  TVectorD vdecorrerrs(nPars);
  for (int i = 0;i<nPars;++i){
    vdecorrerrs[i] = sqrt(vecout[i]);
  }
  //vdecorrerrs.Print();
  return vdecorrerrs;
}

TVectorD poly5mod5FitDecorrelatedParams(TH1D *hist,TString name,TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  TVectorD oriparams = TVectorD(6);
  TVectorD vdecorrparams = TVectorD(6);
  TVectorD vecout;

  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  fitout->GetParameters(pars);
  for (int i = 0;i<nPars;i++){
    oriparams[i] = pars[i];
  }
  //std::cout<<"The original params"<<std::endl;
  //oriparams.Print();

  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  //std::cout<<"The standard covariance matrix"<<std::endl;
  //COV->Print();
  TMatrixD eigvec = COV->EigenVectors(vecout);
  //std::cout<<"The eigenvators of the covariance matrix"<<std::endl;
  //eigvec.Print();
  eigvec.Transpose(eigvec);
  //std::cout<<"The eigenvators of the covariance matrix transposed, hopefully"<<std::endl;
  //eigvec.Print();
  vdecorrparams = eigvec*oriparams;
  //vdecorrparams.Print() ;
  return vdecorrparams;
}



TVectorD poly5mod5FitDecorrParamsShifted(TH1D *hist,TString name,int multi,TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *poly5mod5fit = new TF1(name,poly5Model5,lowr,highr,6);
  TVectorD oriparams = TVectorD(6);
  TVectorD vdecorrparams = TVectorD(6);
  TVectorD vdecorrerrs = TVectorD(6);
  TVectorD vecout = TVectorD(6);

  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  fitout->GetParameters(pars);
  for (int i = 0;i<nPars;i++){
    oriparams[i] = pars[i];
  }
  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  TMatrixD eigvec = COV->EigenVectors(vdecorrerrs);
  eigvec.Transpose(eigvec);
  vdecorrparams = eigvec*oriparams;
  vecout = eigvec.Invert()*(vdecorrparams+vdecorrerrs);
  return vecout;
}

TF1 * poly5mod4Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5mod4fit = new TF1(name,poly5Model4,lowr,highr,6);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly4Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly4fit = new TF1(name,poly4Model,lowr,highr,5);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly3Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly3fit = new TF1(name,poly3Model,lowr,highr,4);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * poly6Fit(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly6fit = new TF1(name,poly6Model,lowr,highr,7);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}


Double_t poly5Model2(Double_t *X, Double_t *par){
  Double_t x = X[0];
  Double_t fitval = par[6]*(par[0]+par[1]*x+par[2]*std::pow(x,2.0)+par[3]*std::pow(x,3.0)+par[4]*std::pow(x,4.0)+par[5]*std::pow(x,5));
  return fitval;
}

TF1 * poly5Fit2(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  TF1 *poly5fit = new TF1(name,poly5Model,lowr,highr,7);
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

TVectorD gaus2Fit2DecorrParamsShifted(TH1D *hist,TString name,int multi,TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TVectorD oriparams = TVectorD(6);
  TVectorD vdecorrparams = TVectorD(6);
  TVectorD vdecorrerrs = TVectorD(6);
  TVectorD vecout = TVectorD(6);

  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  fitout->GetParameters(pars);
  for (int i = 0;i<nPars;i++){
    oriparams[i] = pars[i];
  }
  ROOT::Math::WrappedTF1 wfit(*fitout);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  TMatrixD eigvec = COV->EigenVectors(vdecorrerrs);
  eigvec.Transpose(eigvec);
  vdecorrparams = eigvec*oriparams;
  vecout = eigvec.Invert()*(vdecorrparams+vdecorrerrs);
  return vecout;
}


TGraph *erruptt;
TGraph *errdwntt;

Double_t upfunctt(Double_t *x,Double_t *t){return erruptt->Eval(x[0]);}
Double_t dwnfunctt(Double_t *x,Double_t *t){return errdwntt->Eval(x[0]);}

vector<TF1 *> gaus2Fit2ErrFunctionsGraphs(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;
  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
    fitstartbin = fitplacer+1;
  }
  else {
    fitstartbin = fitplacer;
  }
  int finalbin = (highr-histlowed)/binwidth;

  erruptt = new TGraph(finalbin+1);
  errdwntt = new TGraph(finalbin+1);
  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
  for (int ib = fitstartbin;ib<=finalbin+1;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    erruptt->SetPoint(ib,x,(*fitout)(x)+sigma);
    errdwntt->SetPoint(ib,x,(*fitout)(x)-sigma);
  }

  TF1 *fitoutup  = new TF1("fitoutup",upfunctt,lowr,highr,0);
  TF1 *fitoutdwn = new TF1("fitoutdwn",dwnfunctt,lowr,highr,0);

  /*
  std::cout<<"DEBUG DEBUG DEBUG TTBAR ENVELOPE FUNCTIONS"<<std::endl;
  for (int ib = fitstartbin;ib<=finalbin+1;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetError(wfit,x,pars,*COV,nPars);
    double upgx;
    double upgy;
    errup->GetPoint(ib,upgx,upgy);
    std::cout<<"The histogram up value y: "<<upgy<<std::endl;
    std::cout<<"The       fit up value: "<<(*fitoutup)(x)<<std::endl;
  }
  */
  vector<TF1 *> vecout;
  vecout.push_back(fitoutup);
  vecout.push_back(fitoutdwn);
  return vecout;
}


  

TH1D * gaus2Fit2ErrBands(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  TH1D *errhist = new TH1D("errhist","Histogram with 1 sigma band for fit",nbins,histlowed,histhied);

  fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
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
    errhist->SetBinError(ib,sigma);
  }
  return errhist;
}

TH1D * gaus2Fit2ErrBands2(TH1D *hist, TString name, TString opt="R0+",int lowr=30, int highr=400) {
  const Int_t nPars = 6;
  Double_t pars[nPars], grad[nPars];
  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(1,120);
  gaus2fit->SetParameter(2,40);
  gaus2fit->SetParameter(3,175);
  gaus2fit->SetParameter(4,10);
  TFitResultPtr result = hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);
  //TFitResultPtr res = hist->Fit(

  int nbins  = hist->GetNbinsX();
  float binwidth  = hist->GetBinWidth(1);
  float histlowed = hist->GetBinLowEdge(1);
  float histhied  = hist->GetBinLowEdge(nbins)+binwidth;

  TH1D *errhist = new TH1D("errhist","Histogram with 1 sigma band for fit",nbins,histlowed,histhied);
  TMatrixDSym cov = result->GetCovarianceMatrix();
  
  //fitout->GetParameters(pars);
  ROOT::Math::WrappedTF1 wfit(*fitout);
  float fitplace = (lowr-histlowed)/binwidth;
  float fitplacer = round((lowr-histlowed)/binwidth);
  float fitstartbin;
  if (std::abs(fitplace - fitplacer) < 0.5) {
    fitstartbin = fitplacer+1;
  }
  else {
    fitstartbin = fitplacer;
  }

  int finalbin = (highr-histlowed)/binwidth;
  for (int ib = fitstartbin;ib<=finalbin;ib++){
    double x = hist->GetBinCenter(ib);
    double sigma = GetErrorGEC(wfit,x,pars,cov,nPars);
    errhist->SetBinContent(ib,(*fitout)(x));
    errhist->SetBinError(ib,sigma);
  }
  //*/
  return errhist;
}

TF1 * gaus2Fit2SetParsAndErrs(TString name,TVector pars,int lowr=30, int highr=400) {
  TF1 *gaus2fit = new TF1(name,gaus2Model2,lowr,highr,6);
  gaus2fit->SetParameter(0,pars[0]);
  gaus2fit->SetParameter(1,pars[1]);
  gaus2fit->SetParameter(2,pars[2]);
  gaus2fit->SetParameter(3,pars[3]);
  gaus2fit->SetParameter(4,pars[4]);
  gaus2fit->SetParameter(5,pars[5]);
  //hist->Fit(name,opt);
  //TF1* fitout = hist->GetFunction(name);

  return gaus2fit;
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
  Double_t fitval = par[0]*(TMath::Gaus(x,par[1],par[2])+par[3]*x+par[4]);
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
  //gausPoly1fit->SetParameter(5,hist->GetBinContent(6));
  hist->Fit(name,opt);
  TF1* fitout = hist->GetFunction(name);

  return fitout;
}

TF1 * gausPoly1FitSetParsAndErrs(TString name,TVector pars,int lowr=30, int highr=400) {
  TF1 *gausPoly1fit = new TF1(name,gausPoly1Model,lowr,highr,5);
  gausPoly1fit->SetParameter(0,pars[0]);
  gausPoly1fit->SetParameter(1,pars[1]);
  gausPoly1fit->SetParameter(2,pars[2]);
  gausPoly1fit->SetParameter(3,pars[3]);
  gausPoly1fit->SetParameter(4,pars[4]);
  //gausPoly1fit->SetParameter(5,pars[5]);
  return gausPoly1fit;
}

Double_t totalBkgModel(Double_t *X, Double_t *par){
  Double_t fitval = par[17]*poly5Model5(X,&par[0])+gaus2Model2(X,&par[6])+gausPoly1Model(X,&par[12]);
  return fitval;
}

int lowsbl = 30;
int lowsbh = 70;
int highsbl = 150;
int highsbh = 250;
Bool_t validationregion;
int validationlsbh = 55.;

Double_t totalBkgModelBlind(Double_t *X, Double_t *par){
  if (validationregion && X[0] > validationlsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  if (!validationregion && X[0] > lowsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t fitval = par[17]*poly5Model5(X,&par[0])+gaus2Model2(X,&par[6])+gausPoly1Model(X,&par[12]);
  return fitval;
}

TF1* dyfunctionInTotalEnvelope;
TF1* ttfunctionInTotalEnvelope;

Double_t totalBkgModelTTEnvelope(Double_t *X, Double_t *par){
  Double_t fitval = par[11]*poly5Model5(X,&par[0])+ttfunctionInTotalEnvelope->Eval(*X)+gausPoly1Model(X,&par[6]);
  return fitval;
}

Double_t totalBkgModelDYEnvelope(Double_t *X, Double_t *par){
  //Double_t x = X[0];
  Double_t fitval = par[11]*dyfunctionInTotalEnvelope->Eval(*X)+gaus2Model2(X,&par[0])+gausPoly1Model(X,&par[6]);
  return fitval;
}

Double_t totalBkgModelBlindTTEnvelope(Double_t *X, Double_t *par){
  if (validationregion && X[0] > validationlsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  if (!validationregion && X[0] > lowsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t fitval = par[11]*poly5Model5(X,&par[0])+ttfunctionInTotalEnvelope->Eval(*X)+gausPoly1Model(X,&par[6]);
  return fitval;
}

Double_t totalBkgModelBlindDYEnvelope(Double_t *X, Double_t *par){
  if (validationregion && X[0] > validationlsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  if (!validationregion && X[0] > lowsbh && X[0] < highsbl) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t fitval = par[11]*dyfunctionInTotalEnvelope->Eval(*X)+gaus2Model2(X,&par[0])+gausPoly1Model(X,&par[6]);
  return fitval;
}

vector<TF1 *> totalFitTTEnvelope(TH1D *hist, TF1 *ttfit, TH1D *dyhist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",int lsbl=30,int hsbh=250) {
  TF1 *dyfit = poly5mod5Fit(dyhist,"dyl","QER0+",lsbl,225);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","QER0+",lsbl,hsbh);
  Double_t parIn[12];
  dyfit->GetParameters(&parIn[0]);
  vvfit->GetParameters(&parIn[6]);
  ttfunctionInTotalEnvelope = ttfit;

  std::cout<<"    You have the max fixed at "<<225<<std::endl;
  //std::cout<<"    Doing total bkg fit "<<std::endl;
  TF1 *totalfit = new TF1("totalfitttenv",totalBkgModelTTEnvelope,lsbl,225,12);//17 without norms
  totalfit->FixParameter(0,parIn[0]);
  totalfit->FixParameter(1,parIn[1]);
  totalfit->FixParameter(2,parIn[2]);
  totalfit->FixParameter(3,parIn[3]);
  totalfit->FixParameter(4,parIn[4]);
  totalfit->FixParameter(5,parIn[5]);
  totalfit->FixParameter(6,parIn[6]);
  totalfit->FixParameter(7,parIn[7]);
  totalfit->FixParameter(8,parIn[8]);
  totalfit->FixParameter(9,parIn[9]);
  totalfit->FixParameter(10,parIn[10]);
  //totalfit->FixParameter(11,parIn[11]);
  hist->Fit("totalfitttenv","ER0L+");
  TF1 *totmcfit = hist->GetFunction("totalfitttenv");
  
  //Now do the extrapolation
  //std::cout<<"    Doing sb data fit "<<std::endl;
  TF1 *sbdatfit = new TF1("sbdatfitttshiftenv",totalBkgModelBlindTTEnvelope,lsbl,225,12);
  sbdatfit->FixParameter(0,parIn[0]);
  sbdatfit->FixParameter(1,parIn[1]);
  sbdatfit->FixParameter(2,parIn[2]);
  sbdatfit->FixParameter(3,parIn[3]);
  sbdatfit->FixParameter(4,parIn[4]);
  sbdatfit->FixParameter(5,parIn[5]);
  sbdatfit->FixParameter(6,parIn[6]);
  sbdatfit->FixParameter(7,parIn[7]);
  sbdatfit->FixParameter(8,parIn[8]);
  sbdatfit->FixParameter(9,parIn[9]);
  sbdatfit->FixParameter(10,parIn[10]);
  //sbdatfit->FixParameter(11,parIn[11]);
  dathist->Fit("sbdatfitttshiftenv","ELR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfitttshiftenv");

  //Double_t parExtrap[18];
  //totsbdatfit->GetParameters(parExtrap);
  TF1 *totextrap = new TF1("totextrapttenv",totalBkgModelTTEnvelope,lsbl,225,12);
  totextrap->FixParameter(0,parIn[0]);
  totextrap->FixParameter(1,parIn[1]);
  totextrap->FixParameter(2,parIn[2]);
  totextrap->FixParameter(3,parIn[3]);
  totextrap->FixParameter(4,parIn[4]);
  totextrap->FixParameter(5,parIn[5]);
  totextrap->FixParameter(6,parIn[6]);
  totextrap->FixParameter(7,parIn[7]);
  totextrap->FixParameter(8,parIn[8]);
  totextrap->FixParameter(9,parIn[9]);
  totextrap->FixParameter(10,parIn[10]);
  totextrap->FixParameter(11,totsbdatfit->GetParameter(11));
  //totextrap->FixParameter(12,totsbdatfit->GetParameter(12));

  std::vector<TF1*> fitvector;
  fitvector.push_back(totmcfit);
  fitvector.push_back(totsbdatfit);
  fitvector.push_back(totextrap);
  return fitvector;
}

vector<TF1 *> totalFitDYEnvelope(TH1D *hist, TF1 *dyfit, TH1D *tthist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",int lsbl=30,int hsbh=250) {
  
  //std::cout<<"  TT fit    "<<std::endl;
  TF1 *ttfit = gaus2Fit2(tthist,"ttl","QER0+",lsbl,400);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","QER0+",lsbl,hsbh);
  Double_t parIn[12];
  ttfit->GetParameters(&parIn[0]);
  vvfit->GetParameters(&parIn[6]);
  dyfunctionInTotalEnvelope = dyfit;

  std::cout<<"    You have the max fixed at "<<225<<std::endl;
  //std::cout<<"    Doing total bkg fit "<<std::endl;
  TF1 *totalfit = new TF1("totalfitdyenv",totalBkgModelDYEnvelope,lsbl,225,12);//17 without norms
  totalfit->FixParameter(0,parIn[0]);
  totalfit->FixParameter(1,parIn[1]);
  totalfit->FixParameter(2,parIn[2]);
  totalfit->FixParameter(3,parIn[3]);
  totalfit->FixParameter(4,parIn[4]);
  totalfit->FixParameter(5,parIn[5]);
  totalfit->FixParameter(6,parIn[6]);
  totalfit->FixParameter(7,parIn[7]);
  totalfit->FixParameter(8,parIn[8]);
  totalfit->FixParameter(9,parIn[9]);
  totalfit->FixParameter(10,parIn[10]);
  //totalfit->FixParameter(11,parIn[11]);
  hist->Fit("totalfitdyenv","ER0L+");
  TF1 *totmcfit = hist->GetFunction("totalfitdyenv");
  
  //Now do the extrapolation
  //std::cout<<"    Doing sb data fit "<<std::endl;
  TF1 *sbdatfit = new TF1("sbdatfitdyshiftenv",totalBkgModelBlindDYEnvelope,lsbl,225,12);
  sbdatfit->FixParameter(0,parIn[0]);
  sbdatfit->FixParameter(1,parIn[1]);
  sbdatfit->FixParameter(2,parIn[2]);
  sbdatfit->FixParameter(3,parIn[3]);
  sbdatfit->FixParameter(4,parIn[4]);
  sbdatfit->FixParameter(5,parIn[5]);
  sbdatfit->FixParameter(6,parIn[6]);
  sbdatfit->FixParameter(7,parIn[7]);
  sbdatfit->FixParameter(8,parIn[8]);
  sbdatfit->FixParameter(9,parIn[9]);
  sbdatfit->FixParameter(10,parIn[10]);
  //sbdatfit->FixParameter(11,parIn[11]);
  dathist->Fit("sbdatfitdyshiftenv","ELR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfitdyshiftenv");

  //Double_t parExtrap[18];
  //totsbdatfit->GetParameters(parExtrap);
  TF1 *totextrap = new TF1("totextrapdyenv",totalBkgModelDYEnvelope,lsbl,225,12);
  totextrap->FixParameter(0,parIn[0]);
  totextrap->FixParameter(1,parIn[1]);
  totextrap->FixParameter(2,parIn[2]);
  totextrap->FixParameter(3,parIn[3]);
  totextrap->FixParameter(4,parIn[4]);
  totextrap->FixParameter(5,parIn[5]);
  totextrap->FixParameter(6,parIn[6]);
  totextrap->FixParameter(7,parIn[7]);
  totextrap->FixParameter(8,parIn[8]);
  totextrap->FixParameter(9,parIn[9]);
  totextrap->FixParameter(10,parIn[10]);
  totextrap->FixParameter(11,totsbdatfit->GetParameter(11));
  //totextrap->FixParameter(12,totsbdatfit->GetParameter(12));

  std::vector<TF1*> fitvector;
  fitvector.push_back(totmcfit);
  fitvector.push_back(totsbdatfit);
  fitvector.push_back(totextrap);
  return fitvector;
}

//*/
vector<TF1 *> totalFitDYTemplateVaried(TH1D *hist, TF1 dyfit, TH1D *tthist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",int lsbl=30,int hsbh=250) {
  //std::cout<<"+++++Doing Norm fits DY Shifted++++++"<<std::endl;
  //std::cout<<"  TT fit    "<<std::endl;
  TF1 *ttfit = gaus2Fit2(tthist,"ttl","QER0+",lsbl,400);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","QER0+",lsbl,hsbh);
  Double_t par[17];
  dyfit.GetParameters(&par[0]);
  ttfit->GetParameters(&par[6]);
  vvfit->GetParameters(&par[12]);

  std::cout<<"    You have the max fixed at "<<225<<std::endl;
  //std::cout<<"    Doing total bkg fit "<<std::endl;
  TF1 *totalfit = new TF1("totalfitdy",totalBkgModel,lsbl,225,18);//17 without norms
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
  totalfit->FixParameter(12,par[12]);//norm of dib
  totalfit->FixParameter(13,par[13]);
  totalfit->FixParameter(14,par[14]);
  totalfit->FixParameter(15,par[15]);
  totalfit->FixParameter(16,par[16]);
  hist->Fit("totalfitdy","QER0L+");
  TF1 *totmcfit = hist->GetFunction("totalfitdy");
  
  //Now do the extrapolation
  //std::cout<<"    Doing sb data fit "<<std::endl;
  TF1 *sbdatfit = new TF1("sbdatfitdyshift",totalBkgModelBlind,lsbl,225,18);
  sbdatfit->FixParameter(0,par[0]);
  sbdatfit->FixParameter(1,par[1]);
  sbdatfit->FixParameter(2,par[2]);
  sbdatfit->FixParameter(3,par[3]);
  sbdatfit->FixParameter(4,par[4]);
  sbdatfit->FixParameter(5,par[5]);
  sbdatfit->FixParameter(6,par[6]);
  sbdatfit->FixParameter(7,par[7]);
  sbdatfit->FixParameter(8,par[8]);
  sbdatfit->FixParameter(9,par[9]);
  sbdatfit->FixParameter(10,par[10]);
  sbdatfit->FixParameter(11,par[11]);//norm of ttbar
  sbdatfit->FixParameter(12,par[12]);//norm of dib
  sbdatfit->FixParameter(13,par[13]);
  sbdatfit->FixParameter(14,par[14]);
  sbdatfit->FixParameter(15,par[15]);
  sbdatfit->FixParameter(16,par[16]);
  dathist->Fit("sbdatfitdyshift","QELR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfitdyshift");

  //Double_t parExtrap[18];
  //totsbdatfit->GetParameters(parExtrap);
  TF1 *totextrap = new TF1("totextrap",totalBkgModel,lsbl,225,18);
  totextrap->FixParameter(0,par[0]);
  totextrap->FixParameter(1,par[1]);
  totextrap->FixParameter(2,par[2]);
  totextrap->FixParameter(3,par[3]);
  totextrap->FixParameter(4,par[4]);
  totextrap->FixParameter(5,par[5]);
  totextrap->FixParameter(6,par[6]);
  totextrap->FixParameter(7,par[7]);
  totextrap->FixParameter(8,par[8]);
  totextrap->FixParameter(9,par[9]);
  totextrap->FixParameter(10,par[10]);
  totextrap->FixParameter(11,par[11]);
  totextrap->FixParameter(12,par[12]);
  totextrap->FixParameter(13,par[13]);
  totextrap->FixParameter(14,par[14]);
  totextrap->FixParameter(15,par[15]);
  totextrap->FixParameter(16,par[16]);
  totextrap->FixParameter(17,totsbdatfit->GetParameter(17));

  std::vector<TF1*> fitvector;
  fitvector.push_back(totmcfit);
  fitvector.push_back(totsbdatfit);
  fitvector.push_back(totextrap);
  return fitvector;
}

vector<TF1 *> totalFitTTTemplateVaried(TH1D *hist, TF1 ttfit, TH1D *dyhist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",int lsbl=30,int hsbh=250) {
  TF1 *dyfit = poly5mod5Fit(dyhist,"dyl","QER0+",lsbl,225);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","QER0+",lsbl,hsbh);
  Double_t par[17];
  dyfit->GetParameters(&par[0]);
  ttfit.GetParameters(&par[6]);
  vvfit->GetParameters(&par[12]);
  std::cout<<"    You have the max fixed at "<<225<<std::endl;
  //std::cout<<"    Doing total bkg fit "<<std::endl;
  TF1 *totalfit = new TF1("totalfittt",totalBkgModel,lsbl,225,18);//17 without norms
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
  totalfit->FixParameter(12,par[12]);//norm of dib
  totalfit->FixParameter(13,par[13]);
  totalfit->FixParameter(14,par[14]);
  totalfit->FixParameter(15,par[15]);
  totalfit->FixParameter(16,par[16]);
  hist->Fit("totalfittt","QER0L+");
  TF1 *totmcfit = hist->GetFunction("totalfittt");
  
  //Now do the extrapolation
  //std::cout<<"    Doing sb data fit "<<std::endl;
  TF1 *sbdatfit = new TF1("sbdatfitttshift",totalBkgModelBlind,lsbl,225,18);
  sbdatfit->FixParameter(0,par[0]);
  sbdatfit->FixParameter(1,par[1]);
  sbdatfit->FixParameter(2,par[2]);
  sbdatfit->FixParameter(3,par[3]);
  sbdatfit->FixParameter(4,par[4]);
  sbdatfit->FixParameter(5,par[5]);
  sbdatfit->FixParameter(6,par[6]);
  sbdatfit->FixParameter(7,par[7]);
  sbdatfit->FixParameter(8,par[8]);
  sbdatfit->FixParameter(9,par[9]);
  sbdatfit->FixParameter(10,par[10]);
  sbdatfit->FixParameter(11,par[11]);//norm of ttbar
  sbdatfit->FixParameter(12,par[12]);//norm of dib
  sbdatfit->FixParameter(13,par[13]);
  sbdatfit->FixParameter(14,par[14]);
  sbdatfit->FixParameter(15,par[15]);
  sbdatfit->FixParameter(16,par[16]);
  dathist->Fit("sbdatfitttshift","QELR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfitttshift");

  //Double_t parExtrap[18];
  //totsbdatfit->GetParameters(parExtrap);
  TF1 *totextrap = new TF1("totextrap",totalBkgModel,lsbl,225,18);
  totextrap->FixParameter(0,par[0]);
  totextrap->FixParameter(1,par[1]);
  totextrap->FixParameter(2,par[2]);
  totextrap->FixParameter(3,par[3]);
  totextrap->FixParameter(4,par[4]);
  totextrap->FixParameter(5,par[5]);
  totextrap->FixParameter(6,par[6]);
  totextrap->FixParameter(7,par[7]);
  totextrap->FixParameter(8,par[8]);
  totextrap->FixParameter(9,par[9]);
  totextrap->FixParameter(10,par[10]);
  totextrap->FixParameter(11,par[11]);
  totextrap->FixParameter(12,par[12]);
  totextrap->FixParameter(13,par[13]);
  totextrap->FixParameter(14,par[14]);
  totextrap->FixParameter(15,par[15]);
  totextrap->FixParameter(16,par[16]);
  totextrap->FixParameter(17,totsbdatfit->GetParameter(17));

  std::vector<TF1*> fitvector;
  fitvector.push_back(totmcfit);
  fitvector.push_back(totsbdatfit);
  fitvector.push_back(totextrap);
  return fitvector;
}

vector<TF1 *> totalFit(TH1D *hist, TH1D *dyhist, TH1D *tthist, TH1D *vvhist, TH1D *dathist, TString opt="R0+",Bool_t vr= false,int lsbl=lowsbl,int lsbh = lowsbh,int hsbl=highsbl,int hsbh=225) {
  //TF1 *dyfit = poly5Fit(dyhist,"dyl","QR0+",lsbl,hsbh);
  TF1 *dyfit = poly5mod5Fit(dyhist,"dyl","QER0+",lsbl,225);
  TF1 *ttfit = gaus2Fit2(tthist,"ttl","QER0+",lsbl,400);
  TF1 *vvfit = gausPoly1Fit(vvhist,"vvl","QER0+",lsbl,hsbh);

  Double_t par[17];
  dyfit->GetParameters(&par[0]);
  ttfit->GetParameters(&par[6]);
  vvfit->GetParameters(&par[12]);

  std::cout<<"You have the max fixed at "<<225<<std::endl;
  std::cout<<"Doing total MC background fit with fixed templates"<<std::endl;
  TF1 *totalfit = new TF1("totalfit",totalBkgModel,lsbl,225,18);//17 without norms
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
  totalfit->FixParameter(12,par[12]);//norm of dib
  totalfit->FixParameter(13,par[13]);
  totalfit->FixParameter(14,par[14]);
  totalfit->FixParameter(15,par[15]);
  totalfit->FixParameter(16,par[16]);
  hist->Fit("totalfit","QER0L+");
  TF1 *totmcfit = hist->GetFunction("totalfit");

  Double_t parData[17];
  totmcfit->GetParameters(parData);

  if (vr){
    std::cout<<"Using validation studies SB"<<std::endl;
    lsbh = validationlsbh;
    validationregion = true;
  }
  
  //Now do the extrapolation
  std::cout<<"Doing sideband fit with fixed templates"<<std::endl;
  TF1 *sbdatfit = new TF1("sbdatfit",totalBkgModelBlind,lsbl,hsbh,18);
  sbdatfit->SetParameters(parData);
  sbdatfit->FixParameter(0,parData[0]);
  sbdatfit->FixParameter(1,parData[1]);
  sbdatfit->FixParameter(2,parData[2]);
  sbdatfit->FixParameter(3,parData[3]);
  sbdatfit->FixParameter(4,parData[4]);
  sbdatfit->FixParameter(5,parData[5]);
  sbdatfit->FixParameter(6,parData[6]);
  sbdatfit->FixParameter(7,parData[7]);
  sbdatfit->FixParameter(8,parData[8]);
  sbdatfit->FixParameter(9,parData[9]);
  sbdatfit->FixParameter(10,parData[10]);
  sbdatfit->FixParameter(11,parData[11]);
  sbdatfit->FixParameter(12,parData[12]);
  sbdatfit->FixParameter(13,parData[13]);
  sbdatfit->FixParameter(14,parData[14]);
  sbdatfit->FixParameter(15,parData[15]);
  sbdatfit->FixParameter(16,parData[16]);

  dathist->Fit("sbdatfit","QELR0+");
  TF1 *totsbdatfit = dathist->GetFunction("sbdatfit");

  Double_t parExtrap[18];
  totsbdatfit->GetParameters(parExtrap);
  std::cout<<"remaking the extrapolation fit with new params"<<std::endl;
  TF1 *totextrap = new TF1("totextrap",totalBkgModel,lsbl,hsbh,18);
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
  totextrap->FixParameter(17,parExtrap[17]);
  dathist->Fit("totextrap","QER0+");
  TF1 *totnorm = dathist->GetFunction("totextrap");

  //Now seperate fits for visualization
  TF1 *flsb = new TF1("flsb",totalBkgModelBlind,lsbl,lsbh,18);
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
  flsb->FixParameter(17,parExtrap[17]);
  dathist->Fit("flsb","QR0+");
  TF1 *lsbdatfit = dathist->GetFunction("flsb");

  TF1 *fhsb = new TF1("fhsb",totalBkgModelBlind,hsbl,hsbh,18);
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
  fhsb->FixParameter(17,parExtrap[17]);
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


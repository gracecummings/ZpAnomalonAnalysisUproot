import tdrstyle
import CMS_lumi
import ROOT
import glob
import os
import gecorg as go
import numpy as np
import pandas as pd
import configparser

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

def setLogAxis(pad,islog):
    if islog:
        pad.SetLogy()

def plotMsd(pad,hist,islog=False,logmin=0.1,isData=False):
    maxi = hist.GetMaximum()
    mr   = round(maxi,0)
    histmax = mr+mr*0.30
    histmin = 0
    if islog:
        histmax = mr*10
        histmin = logmin
    hist.SetMaximum(histmax)
    hist.SetMinimum(histmin)
    hist.SetMarkerStyle(8)
    hist.SetMarkerSize(0.5)
    if isData:
        hist.SetMarkerColor(ROOT.kBlack)
        drawopts = "SAMEE2"
    else:
        hist.SetMarkerColor(ROOT.kBlue)
        drawopts = "E1"
    xax = hist.GetXaxis()
    yax = hist.GetYaxis()
    xax.SetTitle("m_{j}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    
    hist.Draw(drawopts)

ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
ROOT.gSystem.Load("cfunctions/alphafits_C")

if __name__=='__main__':

    #will replace with command line options
    path    = 'analysis_output_ZpAnomalon/2021-06-22_alphaMethStuff/'
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'

    bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    data = go.run2(path,zptcut,hptcut,metcut,btagwp)

    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_h_sd')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    empty9 = empty.Clone()

    htrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","tr","h_h_sd")
    htrtt = bkgs.getAddedHist(empty6,"TT","tr","h_h_sd")
    htrzz = bkgs.getAddedHist(empty7,"ZZTo2L2Q","tr","h_h_sd")
    htrwz = bkgs.getAddedHist(empty8,"WZTo2L2Q","tr","h_h_sd")
    htrvv = htrzz.Clone()
    htrvv.Add(htrwz)

    hdatsb = data.getAddedHist(empty9,"sb","h_h_sd")

    #makes some fits
    dyfit = ROOT.poly5Fit(htrdy,"dyl","R0+",30,250)
    ttfit = ROOT.gaus2Fit(htrtt,"ttl","R0+",30,400)
    vvfit = ROOT.gausErfFit(htrvv,"vvl","R0+",30,400,90,5)

    #label
    dyleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    dyleg.AddEntry(htrdy,"DY","ep")
    dyleg.AddEntry(dyfit,"5th deg poly fit","l")
    dyleg.SetBorderSize(0)
    ttleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    ttleg.AddEntry(htrtt,"TT","ep")
    ttleg.AddEntry(ttfit,"2 Gaussian fit","l")
    ttleg.SetBorderSize(0)
    vvleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    vvleg.AddEntry(htrtt,"VV","ep")
    vvleg.AddEntry(ttfit,"ErfExpGaus Fit","l")
    vvleg.SetBorderSize(0)
    


    #make some output
    tc = ROOT.TCanvas("tc","shapes",1100,400)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,1.0)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,1.0)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,1.0)


    tc.cd()
    p11.Draw()
    p11.cd()
    plotMsd(p11,htrdy)
    CMS_lumi.CMS_lumi(p11,4,13)
    dyfit.Draw('same')
    dyleg.Draw()
    p11.Update()
    
    tc.cd()
    p12.Draw()
    p12.cd()
    plotMsd(p12,htrtt)
    CMS_lumi.CMS_lumi(p12,4,13)
    ttfit.Draw("same")
    ttleg.Draw()
    p12.Update()

    tc.cd()
    p13.Draw()
    p13.cd()
    plotMsd(p13,htrvv)
    CMS_lumi.CMS_lumi(p13,4,13)
    vvfit.Draw("same")
    vvleg.Draw()
    p13.Update()

    tc.cd()

    normshapes = go.makeOutFile('Run2_2017_2018','norm_shapes','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(normshapes)

    #ttbarhist = go.makeOutFile('Run2_2017_2018','ttbar_hist','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    #rootfile = ROOT.TFile(ttbarhist,"recreate")
    #htrtt.Write()
    #rootfile.Close()

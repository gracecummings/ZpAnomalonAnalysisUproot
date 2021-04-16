import argparse
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

def plotMzp(pad,hist,histmax,textForPave):
    hist.SetMaximum(histmax)
    hist.GetXaxis().SetTitle("M_{Z'}")
    hist.GetYaxis().SetTitle("Events / 45 GeV")
    hist.Draw()

if __name__=='__main__':

    #make the output
    tc = ROOT.TCanvas("tc","alphafits",1100,400)
    #tc.Divide(3,1)
    p1 = ROOT.TPad("p1","sb",0,0,0.33,1.0)
    p2 = ROOT.TPad("p2","sr",0.33,0,0.66,1.0)
    p3 = ROOT.TPad("p3","alpha",0.66,0,1.0,1.0)
    
    #will replace with command line options
    #bkg_dir = 'analysis_output_ZpAnomalon/2021-03-31/'
    bkg_dir = 'analysis_output_ZpAnomalon/2021-04-13/'
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    
    bkgs = go.backgrounds(bkg_dir,zptcut,hptcut,metcut,btagwp)

    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    hsb = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
    
    tf2 = ROOT.TFile(bkgs.f17dyjetsr[0])
    empty2 = tf2.Get('h_zp_jigm')
    empty2.Reset("ICESM")
    hsr = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_zp_jigm")

    ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("cfunctions/alphafits_C")
    
    #Draw
    histmax = 17.
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)

    p1.Draw()
    p1.cd()
    plotMzp(p1,hsb,histmax,'sideband')
    CMS_lumi.CMS_lumi(p1,4,13)
    p1.Update()
    sbfit = ROOT.expFit(hsb,"sbl","QR0+")
    sbfit.Draw("SAME")

    label = ROOT.TPaveText(2500,histmax/2,4500,histmax/2+histmax*.2,"NB")
    label.AddText("DY MC Sideband")
    label.AddText("30 < m_{hcand,SD} < 70")
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label.Draw()
    p1.Update()
     
    tc.cd()
    p2.Draw()
    p2.cd()
    plotMzp(p2,hsr,4,"foo")
    srfit = ROOT.expFit(hsr,"srl","R0+")
    srfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p2,4,13)
    p2.Update()

    label2 = ROOT.TPaveText(2500,4/2,4500,4/2+4*.2,"NB")
    label2.AddText("DY MC Signal Region")
    label2.AddText("110 <= m_{hcand,SD} < 150")#higgs mass
    #label2.AddText("70 < m_{hcand,SD} < 110")
    label2.SetFillColor(0)
    label2.Draw()
    p2.Update()
    
    tc.cd()
    p3.Draw()
    p3.cd()
    alpha = ROOT.alphaRatioMakerExp(hsb,hsr)
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()
    #CMS_lumi.CMS_lumi(p3,4,13)
    #p3.Update()

    
    figname = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figname)

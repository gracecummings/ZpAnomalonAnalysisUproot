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
    #tc = ROOT.TCanvas("tc","alphafits",1100,400)
    #p1 = ROOT.TPad("p1","sb",0,0,0.33,1.0)
    #p2 = ROOT.TPad("p2","sr",0.33,0,0.66,1.0)
    #p3 = ROOT.TPad("p3","alpha",0.66,0,1.0,1.0)

    #tc2 = ROOT.TCanvas("tc2","ttvvshapes",1100,400)
    #p21 = ROOT.TPad("p21","ttbar",0,0,0.33,1.0)
    #p22 = ROOT.TPad("p22","vv",0.33,0,0.66,1.0)
    #p23 = ROOT.TPad("p23","fittodata",0.66,0,1.0,1.0)

    #make the output
    tc = ROOT.TCanvas("tc","alphafits",1100,800)
    p11 = ROOT.TPad("p11","ttbar",0,0,0.33,.5)
    p12 = ROOT.TPad("p12","vv",0.33,0,0.66,.5)
    p13 = ROOT.TPad("p13","fittodata",0.66,0,1.0,.5)
    p21 = ROOT.TPad("p21","ttbar",0,.5,0.33,1.0)
    p22 = ROOT.TPad("p22","vv",0.33,.5,0.66,1.0)
    p23 = ROOT.TPad("p23","fittodata",0.66,.5,1.0,1.0)
    
    #will replace with command line options
    bkg_dir = 'analysis_output_ZpAnomalon/2021-04-16/'
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    
    bkgs = go.backgrounds(bkg_dir,zptcut,hptcut,metcut,btagwp)

    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    
    hsbdy = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
    hsrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_zp_jigm")
    hsbtt = bkgs.getAddedHist(empty3,"TT","sb","h_zp_jigm")
    hsbzz = bkgs.getAddedHist(empty4,"ZZTo2L2Q","sb","h_zp_jigm")
    hsbwz = bkgs.getAddedHist(empty5,"WZTo2L2Q","sb","h_zp_jigm")
    hsbvv = hsbzz.Add(hsbwz)

    ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("cfunctions/alphafits_C")
    
    #Draw
    histmax = 15.
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)

    tc.cd()
    p11.Draw()
    p11.cd()
    plotMzp(p11,hsbdy,histmax,'sideband')
    CMS_lumi.CMS_lumi(p11,4,13)
    p11.Update()
    sbfit = ROOT.expFit(hsbdy,"sbl","QR0+")
    sbfit.Draw("SAME")
    label = ROOT.TPaveText(.5,.5,.9,.7,"NBNDC")
    label.AddText("DY MC SB")
    label.AddText("30 < m_{hcand,SD} < 70")
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label.Draw()
    p11.Update()
     
    tc.cd()
    p12.Draw()
    p12.cd()
    plotMzp(p12,hsrdy,6,"foo")
    srfit = ROOT.expFit(hsrdy,"srl","QR0+")
    srfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p12,4,13)
    p12.Update()

    label2 = ROOT.TPaveText(.5,.5,.9,.7,"NBNDC")
    label2.AddText("DY MC SR")
    label2.AddText("110 <= m_{hcand,SD} < 150")#higgs mass
    label2.AddText("70 < m_{hcand,SD} < 110")
    label2.SetFillColor(0)
    label2.Draw()
    p12.Update()
    
    tc.cd()
    p13.Draw()
    p13.cd()
    alpha = ROOT.alphaRatioMakerExp(hsbdy,hsrdy)
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()

    tc.cd()
    p21.Draw()
    p21.cd()
    plotMzp(p21,hsbtt,5,"fff")
    CMS_lumi.CMS_lumi(p21,4,13)
    p21.Update()

    tc.cd()
    p22.Draw()
    p22.cd()
    plotMzp(p22,hsbzz,5,"fff")
    CMS_lumi.CMS_lumi(p22,4,13)
    p22.Update()
    figalpha = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    #figshapes = go.makeOutFile('Run2_2017_2018','alpha_shapes','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figalpha)
    #xtc2.SaveAs(figshapes)

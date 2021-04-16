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

def setLogAxis(pad,islog):
    if islog:
        pad.SetLogy()

def plotMzp(pad,hist,islog=False,logmin=0.1):
    maxi = hist.GetMaximum()
    mr   = round(maxi,0)
    histmax = mr+mr*0.30
    histmin = 0
    if islog:
        histmax = mr*100
        histmin = logmin
    hist.SetMaximum(histmax)
    hist.SetMinimum(histmin)
    hist.SetMarkerStyle(8)
    hist.SetMarkerSize(0.5)
    hist.SetMarkerColor(ROOT.kBlue)
    xax = hist.GetXaxis()
    yax = hist.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 45 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    
    hist.Draw("E1")

def pavedLabel(text):
    label = ROOT.TPaveText(.5,.5,.9,.7,"NBNDC")
    label.AddText(text)
    label.SetFillColor(0)
    label.Draw()

    
if __name__=='__main__':

    #make the output
    #tc = ROOT.TCanvas("tc","alphafits",1100,400)
    #p1 = ROOT.TPad("p1","datsb",0,0,0.5,1.0)
    #p2 = ROOT.TPad("p2","alpha",0.5,0,1.0,1.0)

    tc = ROOT.TCanvas("tc","shapes",1100,800)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,.5)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,.5)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,.5)
    p21 = ROOT.TPad("p21","dysb",0,.5,0.33,1.0)
    p22 = ROOT.TPad("p22","ttsb",0.33,.5,0.66,1.0)
    p23 = ROOT.TPad("p23","vvsb",0.66,.5,1.0,1.0)
        
    
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
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    
    hsbdy = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
    hsrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_zp_jigm")
    hsbtt = bkgs.getAddedHist(empty3,"TT","sb","h_zp_jigm")
    hsrtt = bkgs.getAddedHist(empty6,"TT","sr","h_zp_jigm")
    hsbzz = bkgs.getAddedHist(empty4,"ZZTo2L2Q","sb","h_zp_jigm")
    hsrzz = bkgs.getAddedHist(empty7,"ZZTo2L2Q","sr","h_zp_jigm")
    hsbwz = bkgs.getAddedHist(empty5,"WZTo2L2Q","sb","h_zp_jigm")
    hsrwz = bkgs.getAddedHist(empty8,"WZTo2L2Q","sr","h_zp_jigm")
    hsbvv = hsbzz.Clone()
    hsbvv.Add(hsbwz)
    hsrvv = hsrzz.Clone()
    hsrvv.Add(hsrwz)

    ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("cfunctions/alphafits_C")
    
    #Draw
    histmax = 15.
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)
    islog = True

    tc.cd()
    p21.Draw()
    p21.cd()
    plotMzp(p21,hsbdy)
    CMS_lumi.CMS_lumi(p21,4,13)
    p21.Update()
    sbfit = ROOT.expFit(hsbdy,"sbl","QR0+")
    sbfit.Draw("SAME")
    label = ROOT.TPaveText(.5,.5,.9,.7,"NBNDC")
    label.AddText("DY MC SB")
    label.AddText("30 < m_{hcand,SD} < 70")
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label.Draw()
    p21.Update()
     
    tc.cd()
    p11.Draw()
    p11.cd()
    plotMzp(p11,hsrdy)
    srfit = ROOT.expFit(hsrdy,"srl","QR0+")
    srfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p11,4,13)
    p11.Update()

    label2 = ROOT.TPaveText(.5,.5,.9,.7,"NBNDC")
    label2.AddText("DY MC SR")
    label2.AddText("110 <= m_{hcand,SD} < 150")#higgs mass
    label2.AddText("70 < m_{hcand,SD} < 110")
    label2.SetFillColor(0)
    label2.Draw()
    p11.Update()
    
    #tc.cd()
    #p22.Draw()
    #p22.cd()
    #alpha = ROOT.alphaRatioMakerExp(hsbdy,hsrdy)
    #alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    #alpha.GetXaxis().SetTitle("M_{Z'}")
    #alpha.Draw()

    tc.cd()
    p22.Draw()
    p22.cd()
    setLogAxis(p22,islog)
    plotMzp(p22,hsbtt,islog)
    sbttfit = ROOT.expFit(hsbtt,"sbl","R0+")
    sbttfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p22,4,13)
    p22.Update()
    label3 = ROOT.TPaveText(.5,.65,.9,.7,"NBNDC")
    label3.AddText("TT MC SB")
    label3.SetFillColor(0)
    label3.Draw()
    p22.Update()

    tc.cd()
    p12.Draw()
    p12.cd()
    setLogAxis(p12,islog)
    plotMzp(p12,hsrtt,islog)
    srttfit = ROOT.expFit(hsrtt,"sbl","R0+")
    srttfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p12,4,13)
    p12.Update()
    label4 = ROOT.TPaveText(.5,.65,.9,.7,"NBNDC")
    label4.AddText("TT MC SR")
    label4.SetFillColor(0)
    label4.Draw()
    p12.Update()

    tc.cd()
    p23.Draw()
    p23.cd()
    setLogAxis(p23,islog)
    plotMzp(p23,hsbvv,islog,0.001)
    sbvvfit = ROOT.expFit(hsbvv,"sbl","R0+")
    sbvvfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p23,4,13)
    p23.Update()
    label5 = ROOT.TPaveText(.5,.65,.9,.7,"NBNDC")
    label5.AddText("VV MC SB")
    label5.SetFillColor(0)
    label5.Draw()
    p23.Update()

    tc.cd()
    p13.Draw()
    p13.cd()
    setLogAxis(p13,islog)
    plotMzp(p13,hsrvv,islog,0.001)
    srvvfit = ROOT.expFit(hsrvv,"sbl","R0+")
    srvvfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p13,4,13)
    p13.Update()
    label6 = ROOT.TPaveText(.5,.65,.9,.7,"NBNDC")
    label6.AddText("VV MC SR")
    label6.SetFillColor(0)
    label6.Draw()

    
    figalpha = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    #figshapes = go.makeOutFile('Run2_2017_2018','alpha_shapes','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figalpha)
    #xtc2.SaveAs(figshapes)

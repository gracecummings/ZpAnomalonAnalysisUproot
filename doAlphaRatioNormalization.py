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

    #hsbdy = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_h_sd")
    #hsrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_h_sd")
    #hsbtt = bkgs.getAddedHist(empty3,"TT","sb","h_h_sd")
    #hsrtt = bkgs.getAddedHist(empty6,"TT","sr","h_h_sd")
    #hsbzz = bkgs.getAddedHist(empty4,"ZZTo2L2Q","sb","h_h_sd")
    #hsrzz = bkgs.getAddedHist(empty7,"ZZTo2L2Q","sr","h_h_sd")
    #hsbwz = bkgs.getAddedHist(empty5,"WZTo2L2Q","sb","h_h_sd")
    #hsrwz = bkgs.getAddedHist(empty8,"WZTo2L2Q","sr","h_h_sd")

    htrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","tr","h_h_sd")
    htrtt = bkgs.getAddedHist(empty6,"TT","tr","h_h_sd")
    htrzz = bkgs.getAddedHist(empty7,"ZZTo2L2Q","tr","h_h_sd")
    htrwz = bkgs.getAddedHist(empty8,"WZTo2L2Q","tr","h_h_sd")
    htrvv = htrzz.Clone()
    htrvv.Add(htrwz)
    
    hdatsb = data.getAddedHist(empty9,"sb","h_h_sd")

    #make some output
    tc = ROOT.TCanvas("tc","shapes",1100,400)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,1.0)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,1.0)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,1.0)


    tc.cd()
    p11.Draw()
    p11.cd()
    htrdy.Draw()
    tc.cd()
    p12.Draw()
    p12.cd()
    htrtt.Draw()
    tc.cd()
    p13.Draw()
    p13.cd()
    htrvv.Draw()
    tc.cd()

    normshapes = go.makeOutFile('Run2_2017_2018','norm_shapes','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(normshapes)

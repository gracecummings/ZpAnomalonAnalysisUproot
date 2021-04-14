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

def makeAddedHist(s17,s18,errs17,errs18,xspairs,hsb):
    s17.sort(key = go.orderDY)
    s18.sort(key = go.orderDY)
    errs17.sort(key = go.orderDY)
    errs18.sort(key = go.orderDY)
    bkgdfs = []
    
    for i,f in enumerate(s17):
        #make the hist
        tf = ROOT.TFile(f)
        numevents = float(str(tf.Get('hnevents').GetString()))
        xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
        scale = go.findScale(numevents,xs,41.53)
        h = tf.Get('h_zp_jigm')
        h.Scale(scale)
        hsb.Add(h)

        #calc hist errors
        df = pd.read_pickle(errs17[i])
        sdf = df*scale
        sqrddf = sdf**2
        bkgdfs.append(sqrddf)

    for i,f in enumerate(s18):
        #make the hist
        tf = ROOT.TFile(f)
        numevents = float(str(tf.Get('hnevents').GetString()))
        xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
        scale = go.findScale(numevents,xs,59.74)
        h = tf.Get('h_zp_jigm')
        h.Scale(scale)
        hsb.Add(h)

        #calc hist errors
        df = pd.read_pickle(errs17[i])
        sdf = df*scale
        sqrddf = sdf**2
        bkgdfs.append(sqrddf)

    uncsqdDYJetsdf = sum(bkgdfs)
    uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

    for ibin in range(hsb.GetNbinsX()+1):
        if ibin == 0:
            continue
        else:
            binerr = uncDYJetsdf['h_zp_jigm'][ibin-1]
            hsb.SetBinError(ibin,binerr)
            
    return hsb

def plotMzp(pad,hist,histmax,textForPave):
    hist.SetMaximum(histmax)
    hist.GetXaxis().SetTitle("M_{Z'}")
    hist.GetYaxis().SetTitle("Events / 45 GeV")
    hist.Draw()


#def calcAddedHistErrors(samp,xspairs):
#    samp.sort(key = go.

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
    hist = bkgs.getAddedHist("TT","sb","h_zp_jigm")
    print(hist)

    #scale and add together MC
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
    config.read_file(fp)
    xspairs = config.items("DYJetsToLL")

    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    hsb = tf1.Get('h_zp_jigm')
    hsb.Reset("ICESM")#creates an empty hist with same structure
    hsb = makeAddedHist(bkgs.f17dyjetsb,bkgs.a18dyjetsb,bkgs.f17dyjetsberrs,bkgs.a18dyjetsberrs,xspairs,hsb)
    
    tf2 = ROOT.TFile(bkgs.f17dyjetsr[0])
    hsr = tf2.Get('h_zp_jigm')
    hsr.Reset("ICESM")
    hsr = makeAddedHist(bkgs.f17dyjetsr,bkgs.a18dyjetsr,bkgs.f17dyjetsrerrs,bkgs.a18dyjetsrerrs,xspairs,hsr)
    
    ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("cfunctions/alphafits_C")
    hsrt = hsr.Clone()
    hsbt = hsb.Clone()
    
    #Draw
    histmax = 17.
    ROOT.gStyle.SetOptFit(1011)
    ROOT.gStyle.SetOptStat(0)

    p1.Draw()
    p1.cd()
    plotMzp(p1,hsb,histmax,'sideband')
    tc.Update()
    sbfit = ROOT.expFit(hsbt,"sbl","R0+")
    sbfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p1,4,13)
    p1.Update()

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
    plotMzp(p2,hsr,10,"foo")
    srfit = ROOT.expFit(hsrt,"srl","R0+")
    srfit.Draw("SAME")
    CMS_lumi.CMS_lumi(p2,4,13)
    p2.Update()

    label2 = ROOT.TPaveText(2500,10/2,4500,10/2+10*.2,"NB")
    label2.AddText("DY MC Signal Region")
    label2.AddText("110 <= m_{hcand,SD} < 150")#higgs mass
    #label2.AddText("70 < m_{hcand,SD} < 110")
    label2.SetFillColor(0)
    label2.Draw()
    p2.Update()
    
    tc.cd()
    p3.Draw()
    p3.cd()
    alpha = ROOT.alphaRatioMakerExp(hsbt,hsrt)
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()
    #CMS_lumi.CMS_lumi(p3,4,13)
    #p3.Update()

    
    figname = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figname)

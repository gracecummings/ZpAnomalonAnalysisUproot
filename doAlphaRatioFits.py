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

def plotMzp(pad,hist,islog=False,logmin=0.1,isData=False):
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
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 90 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    
    hist.Draw(drawopts)
    
if __name__=='__main__':

    #make the output
    tc = ROOT.TCanvas("tc","shapes",1100,800)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,.5)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,.5)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,.5)
    p21 = ROOT.TPad("p21","dysb",0,.5,0.33,1.0)
    p22 = ROOT.TPad("p22","ttsb",0.33,.5,0.66,1.0)
    p23 = ROOT.TPad("p23","vvsb",0.66,.5,1.0,1.0)

        
    
    #will replace with command line options
    path    = 'analysis_output_ZpAnomalon/2021-06-22_alphaMethStuff/'
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    dynorm  = 1.
    validation = False
    rstr = "signalblind"
    if validation:
        rstr = "validationblind"
        dynorm = np.load(path+'/Run2_2017_2018_dynormalization_validationblind_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.npy')[0]
    else:
        dynorm = np.load(path+'Run2_2017_2018_dynormalization_signalblind_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.npy')[0]

    print("Using the DY normalization factor: ",dynorm)
    
    bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    data = go.run2(path,zptcut,hptcut,metcut,btagwp)
    
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
    empty9 = empty.Clone()

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

    hdatsb = data.getAddedHist(empty9,"sb","h_zp_jigm")
    hdatsbsub = hdatsb.Clone()

    #Apply the normalization
    hsbdy.Scale(dynorm)
    hsrdy.Scale(dynorm)
    
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
    #hsbdy2 = hsbdy.Clone()
    print("============do dy sb fit====================")
    sbfit = ROOT.expFit(hsbdy,"sbl","R0+",1500,5000)
    uncbands = ROOT.expFitErrBands(hsbdy,"sbl","QR0+",2,1500,5000)
    plotMzp(p21,hsbdy)
    CMS_lumi.CMS_lumi(p21,4,13)
    p21.Update()
    uncbands.SetStats(ROOT.kFALSE)
    uncbands.SetFillColorAlpha(2,0.35)
    uncbands.SetMarkerStyle(8)
    uncbands.SetMarkerSize(0)
    sbfit.SetMarkerStyle(8)
    sbfit.SetMarkerSize(0)
    uncbands.Draw("e4same")#err bands are 2 sigma bands
    #uncbands.Draw()#err bands are 2 sigma bands
    sbfit.Draw("SAME")
    hsbdy.Draw("SAME")

    l21 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l21.AddEntry(hsbdy,"DY Jets SB MC","ep")
    l21.AddEntry(sbfit,"2 Param Exp fit","l")
    l21.AddEntry(uncbands,"2 $\sigma$ uncertainty","f")
    l21.SetBorderSize(0)
    l21.Draw()
    
    label = ROOT.TPaveText(.5,.4,.9,.5,"NBNDC")
    label.AddText("30 < m_{hcand,SD} < 70")
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label.Draw()
    p21.Update()
     
    tc.cd()
    p11.Draw()
    p11.cd()
    plotMzp(p11,hsrdy)
    print("===================doign dy sr fit===================")
    srdyunc = ROOT.expFitErrBands(hsrdy,"sbl","QR0+",2,1500,4000)
    srdyunc.SetFillColor(2)
    srdyunc.SetMarkerSize(0)
    srfit = ROOT.expFit(hsrdy,"srl","R0+",1500,4000)#be aware, diff range from err and sb
    CMS_lumi.CMS_lumi(p11,4,13)
    srdyunc.Draw("e4,same,c")
    srfit.Draw("same")
    hsrdy.Draw("SAME")
    p11.Update()

    l11 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l11.AddEntry(hsrdy,"DY Jets SR MC","ep")
    l11.AddEntry(srfit,"2 Param Exp fit","l")
    l11.AddEntry(srdyunc,"\(2 \sigma\ uncertainty\)","f")
    l11.SetBorderSize(0)
    l11.Draw()
    
    label2 = ROOT.TPaveText(.5,.4,.9,.5,"NBNDC")
    label2.AddText("110 <= m_{hcand,SD} < 150")#higgs mass
    label2.AddText("70 < m_{hcand,SD} < 110")
    label2.SetFillColor(0)
    label2.Draw()
    p11.Update()
    

    tc.cd()
    p22.Draw()
    p22.cd()
    setLogAxis(p22,islog)
    plotMzp(p22,hsbtt,islog)
    sbttfit = ROOT.expFit(hsbtt,"sbl","QR0+",1500,5000)
    sbttunc = ROOT.expFitErrBands(hsbtt,"sbl","QR0+",2,1500,5000)
    sbttunc.SetFillColor(2)
    sbttunc.SetMarkerSize(0)
    CMS_lumi.CMS_lumi(p22,4,13)
    sbttunc.Draw("e3,same,c")
    sbttfit.Draw("SAME")
    hsbtt.Draw("SAME")
    p22.Update()
    l22 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l22.AddEntry(hsbtt,"ttbar SB MC","ep")
    l22.AddEntry(sbttfit,"2 Param Exp fit","l")
    l22.AddEntry(sbttunc,"1 $\sigma$ uncertainty","f")
    l22.SetBorderSize(0)
    l22.Draw()
    p22.Update()


    tc.cd()
    p12.Draw()
    p12.cd()
    sbdatfit = ROOT.expFit(hdatsb,"sbl","QR0+",1500,5000)
    sbdatunc = ROOT.expFitErrBands(hdatsb,"sbl","QR0+",2,1500,5000)
    sbdatunc.SetFillColor(2)
    sbdatunc.SetMarkerSize(0)
    l111 = ROOT.TLegend(0.55,0.4,0.9,0.8)
    l111.AddEntry(hdatsb,"Data SB - full","ep")
    l111.AddEntry(sbdatfit,"2 Param Exp fit","l")
    l111.AddEntry(sbdatunc,"2 $\sigma$ uncertainty","f")
    l111.SetBorderSize(0)

    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    
    hsbkg = ROOT.THStack("hsbkg","")
    hsbdyc = hsbdy.Clone()
    hsbdyc.SetFillColor(bkgcols[0])
    hsbdyc.SetLineColor(bkgcols[0])
    hsbttc = hsbtt.Clone()
    hsbttc.SetFillColor(bkgcols[1])
    hsbttc.SetLineColor(bkgcols[1])
    hsbwzc = hsbwz.Clone() 
    hsbwzc.SetFillColor(bkgcols[2])
    hsbwzc.SetLineColor(bkgcols[2])
    hsbzzc = hsbwz.Clone() 
    hsbzzc.SetFillColor(bkgcols[3])
    hsbzzc.SetLineColor(bkgcols[3])
    hsbkg.Add(hsbzzc)
    hsbkg.Add(hsbwzc)
    hsbkg.Add(hsbttc)
    hsbkg.Add(hsbdyc)
    l111.AddEntry(hsbdyc,"DYJetsToLL","f")
    l111.AddEntry(hsbttc,"TT","f")
    l111.AddEntry(hsbwzc,"WZ","f")
    l111.AddEntry(hsbwzc,"ZZ","f")
    hsbkg.SetMaximum(50.)
    hsbkg.SetMinimum(0.)

    hsbkg.Draw("HIST")
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 45 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    plotMzp(p12,hdatsb,isData=True)
    sbdatunc.Draw("e3,same,c")
    sbdatfit.Draw("SAME")
    hdatsb.Draw("SAME,E1")
    CMS_lumi.CMS_lumi(p12,4,13)
    p12.Update()
    l111.Draw()
    
    tc.cd()
    p23.Draw()
    p23.cd()
    setLogAxis(p23,islog)
    plotMzp(p23,hsbvv,islog,0.001)
    sbvvfit = ROOT.expFit(hsbvv,"sbl","QR0+",1500,5000)
    sbvvunc = ROOT.expFitErrBands(hsbvv,"sbl","QR0+",1,1500,5000)
    sbvvunc.SetFillColor(2)
    sbvvunc.SetMarkerSize(0)
    CMS_lumi.CMS_lumi(p23,4,13)
    sbvvunc.Draw("e4,same,c")
    sbvvfit.Draw("SAME")
    hsbvv.Draw("SAME")

    p23.Update()
    l23 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l23.AddEntry(hsbvv,"VV SB MC","ep")
    l23.AddEntry(sbvvfit,"2 Param Exp fit","l")
    l23.AddEntry(sbvvunc,"1 $\sigma$ uncertainty","f")
    l23.SetBorderSize(0)
    l23.Draw()
    p23.Update()


    tc.cd()
    p13.Draw()
    p13.cd()
    print("=========do alpha ratio fit==================")
    alpha = ROOT.alphaRatioMakerExp(hsbdy,hsrdy)
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()

    
    figshapes = go.makeOutFile('Run2_2017_2018','alpha_shapes_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    datavis = go.makeOutFile('Run2_2017_2018','alpha_sub_tester_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figshapes)

    #Subtracted Background canvas
    tc2 = ROOT.TCanvas("tc","ratio",1200,500)
    pd112 = ROOT.TPad("pd112","datsb",0,0,.33,1)
    pd122 = ROOT.TPad("pd122","alpha",.33,0,.67,1)
    pd132 = ROOT.TPad("pd122","alpha",.67,0,1,1)

    hdatsbsub.Add(sbttunc,-1)
    hdatsbsub.Add(sbvvunc,-1)

    #print("Integral of data sub sideband ",hdatsbsub.Integral())
    #print("Integral of data sideband     ",hdatsb.Integral())
    
    tc2.cd()
    pd112.Draw()
    pd112.cd()
    print("=================doing sb data fit==============")
    sbdatsubfit = ROOT.expFit(hdatsbsub,"sbl","R0+",1500,3000)
    sbdatsubunc = ROOT.expFitErrBands(hdatsbsub,"sbl","QR0+",2,1500,3000)
    sbdatsubunc.SetFillColor(2)
    sbdatsubunc.SetMarkerSize(0)
    #sbfit.SetFillColor(bkgcols[0])
    #sbfit.SetLineColor(bkgcols[0])
    #sbfit.SetFillStyle(1001)
    hsbdy.SetFillColor(bkgcols[0])
    hsbdy.SetLineColor(bkgcols[0])
    plotMzp(pd112,hdatsbsub,isData=True)
    #sbfit.Draw("SAMEC")
    CMS_lumi.CMS_lumi(pd112,4,13)
    hsbdy.Draw("HISTSAME")
    sbdatsubunc.Draw("e3,same,c")
    sbdatsubfit.Draw("SAME")
    hdatsbsub.GetXaxis().SetRangeUser(1500,5000)
    hdatsbsub.Draw("SAME")

    lstack = ROOT.TLegend(0.40,0.6,0.93,0.8)
    lstack.AddEntry(hdatsbsub,"Data SB, VV tt subtracted","ep")
    lstack.AddEntry(sbdatsubfit,"2 Param Exp fit","l")
    lstack.AddEntry(sbdatsubunc,"2 $\sigma$ uncertainty","f")
    #lstack.AddEntry(sbfit,"DY MC Fit","f")
    lstack.AddEntry(hsbdy,"DY SB MC","f")
    lstack.SetBorderSize(0)
    lstack.Draw()

    tc2.cd()
    pd122.Draw()
    pd122.cd()
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()

    tc2.cd()
    pd132.Draw()
    pd132.cd()
    hsrdy.SetFillColor(bkgcols[0])
    hsrdy.SetLineColor(bkgcols[0])
    hdatsbsub1 = hdatsbsub.Clone()#check errors
    print("=========doing sb extrapolation fit==================")
    extrap  = ROOT.alphaExtrapolation(hsbdy,hsrdy,hdatsbsub1)
    hsrdy.GetXaxis().SetRangeUser(1500,5000)
    hsrdy.GetYaxis().SetRangeUser(0,8)
    hsrdy.Draw("HIST")
    extrap.Draw("SAME")
    CMS_lumi.CMS_lumi(pd132,4,13)

    lstack1 = ROOT.TLegend(0.30,0.6,0.93,0.8)
    lstack1.AddEntry(extrap,"DY Predict, alpha*(data SB fit)","p")
    lstack1.AddEntry(hsrdy,"DY MC SR","f")
    lstack1.SetBorderSize(0)
    lstack1.Draw()

    tc2.SaveAs(datavis)

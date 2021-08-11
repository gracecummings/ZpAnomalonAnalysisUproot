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
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
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

    #ranges in question
    lsb = [30,70]
    hsb = [150,400]
    br  = [70,150]
    sr  = [110,150]
    vr  = [55,70]
    
    validation = False
    rstr = "signalblind"
    if validation:
        rstr = "validationblind"

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

    #Gather basics histograms
    hdatsb = data.getAddedHist(empty9,"sb","h_h_sd")
    htrdy  = bkgs.getAddedHist(empty2,"DYJetsToLL","tr","h_h_sd")
    htrtt  = bkgs.getAddedHist(empty6,"TT","tr","h_h_sd")
    htrzz  = bkgs.getAddedHist(empty7,"ZZTo2L2Q","tr","h_h_sd")
    htrwz  = bkgs.getAddedHist(empty8,"WZTo2L2Q","tr","h_h_sd")
    htrvv  = htrzz.Clone()
    htrvv.Add(htrwz)

    #Make overall stackplot, og way
    hsbkg = ROOT.THStack("hsbkg","")
    bkgfiles17 = [bkgs.bkgs["DYJetsToLL"][17]["tr"][0],
                  bkgs.bkgs["TT"][17]["tr"][0],
                  bkgs.bkgs["WZTo2L2Q"][17]["tr"][0],
                  bkgs.bkgs["ZZTo2L2Q"][17]["tr"][0]
    ]
    bkgfiles18 = [bkgs.bkgs["DYJetsToLL"][18]["tr"][0],
                  bkgs.bkgs["TT"][18]["tr"][0],
                  bkgs.bkgs["WZTo2L2Q"][18]["tr"][0],
                  bkgs.bkgs["ZZTo2L2Q"][18]["tr"][0]
    ]

    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    info17 = go.prepBkg(bkgfiles17,bkgnames,bkgcols,"xsects_2017.ini",41.53)
    info18 = go.prepBkg(bkgfiles18,bkgnames,bkgcols,"xsects_2017.ini",59.74)
    stackleg = ROOT.TLegend(0.55,0.40,0.93,0.8)
    go.stackBkgMultiYear(info17,info18,'h_h_sd',hsbkg,stackleg,50,0)


    #makes some fits
    dyfit = ROOT.poly5Fit(htrdy,"dyl","QR0+",30,250)
    #ttfit = ROOT.gaus2Fit(htrtt,"ttl","QR0+",30,400)
    ttfit = ROOT.gaus2Fit2(htrtt,"ttl","QR0+",30,400)
    vvfit = ROOT.gausPoly1Fit(htrvv,"vvl","QR0+",30,250,90,5)
    normfits = ROOT.totalFit(hsbkg.GetStack().Last(),htrdy,htrtt,htrvv,hdatsb,"R0+",validation)
    bkgfit = normfits[0]#fit to un-normalized MC
    sbdatfit = normfits[1]#Fit with blinded region, bad for visuals
    totnormfit = normfits[2]#Uses params of fit above to draw full line
    lsbdatfit = normfits[3]#low sideband of fit above
    hsbdatfit = normfits[4]#high sideband of fit 2 above
    dynormprefit = bkgfit.GetParameters()[17] 
    dynormpostfit = totnormfit.GetParameters()[17]

    #Save the normalization
    normfilename = go.makeOutFile('Run2_2017_2018','dynormalization_'+rstr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    normfile     = open(normfilename,'wb')
    np.save(normfile,np.array([dynormpostfit]))
    normfile.close()
    
    #Define some fit styles
    bkgfit.SetLineStyle(6)
    totnormfit.SetLineColor(ROOT.kBlue)
    totnormfit.SetLineStyle(9)

    #Stack for DY Scaling
    hsbkgnorm = ROOT.THStack("hsbkgnorm","")
    htrdyclone = htrdy.Clone()
    htrdyclone.SetFillColor(bkgcols[0])
    htrdyclone.SetLineColor(bkgcols[0])
    htrdyclone.Scale(dynormpostfit)
    htrtt.SetFillColor(bkgcols[1])
    htrtt.SetLineColor(bkgcols[1])
    htrwz.SetFillColor(bkgcols[2])
    htrwz.SetLineColor(bkgcols[2])
    htrzz.SetFillColor(bkgcols[3])
    htrzz.SetLineColor(bkgcols[3])
    hsbkgnorm.Add(htrzz)
    hsbkgnorm.Add(htrwz)
    hsbkgnorm.Add(htrtt)
    hsbkgnorm.Add(htrdyclone)
    hsbkgnorm.SetMaximum(50.)
    hsbkgnorm.SetMinimum(0.0)

    #labels
    dyleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    dyleg.AddEntry(htrdy,"DY","ep")
    dyleg.AddEntry(dyfit,"5th deg poly fit","l")
    dyleg.SetBorderSize(0)
    ttleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    ttleg.AddEntry(htrtt,"TT","ep")
    ttleg.AddEntry(ttfit,"2 Gaussian fit","l")
    ttleg.SetBorderSize(0)
    vvleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    vvleg.AddEntry(htrvv,"VV","ep")
    vvleg.AddEntry(vvfit,"GausPol1 Fit","l")
    vvleg.SetBorderSize(0)
    stackleg.AddEntry(bkgfit,"Bkg MC fit Prefit - Likelihood","l")
    stackleg.SetBorderSize(0)
    normcomplabel = ROOT.TPaveText(0.17,0.8,0.52,0.9,"NBNDC")
    normcomplabel.AddText("Prefit DY Norm:  {0}".format(round(dynormprefit,4)))
    normcomplabel.AddText("Postfit DY Norm: {0}".format(round(dynormpostfit,4)))
    normcomplabel.SetFillColor(0)
    unnormlabel = ROOT.TPaveText(0.6,0.2,0.93,0.35,"NBNDC")
    unnormlabel.AddText("MC w/out SB data norm")
    unnormlabel.SetFillColor(0)
    normlabel = ROOT.TPaveText(0.6,0.2,0.93,0.35,"NBNDC")
    normlabel.AddText("DY MC with SB data norm")
    normlabel.SetFillColor(0)
    brl = ROOT.TLine(br[0],0,br[0],40.)
    brh = ROOT.TLine(br[1],0,br[1],40.)
    vrl = ROOT.TLine(vr[0],0,vr[0],40.)
    srl = ROOT.TLine(sr[0],0,sr[0],40.)
    srlabel = ROOT.TPaveText(125,35,135,40,"NB")
    srlabel.SetFillColor(0)
    srlabel.AddText("SR")
    vrlabel = ROOT.TPaveText(57.5,35,67.5,40,"NB")
    vrlabel.SetFillColor(0)
    vrlabel.AddText("VR")
    zrlabel = ROOT.TPaveText(85,35,95,40,"NB")
    zrlabel.SetFillColor(0)
    zrlabel.AddText("ZR")
    
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
    
    normshapes = go.makeOutFile('Run2_2017_2018','norm_shapes_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(normshapes)

    tc1 = ROOT.TCanvas("tc1","stacked",1500,600)
    pd11 = ROOT.TPad("pd11","bkgonly",0,0,0.5,1.0)
    pd12 = ROOT.TPad("pd12","datfit",0.5,0.0,1.0,1.0)
    #pd12 = ROOT.TPad("pd12","datfit",0.0,0.0,1.0,1.0)

    pd11.Draw()
    pd11.cd()
    hsbkg.Draw('HIST')
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pd11,4,13)
    bkgfit.Draw('SAME')
    totnormfit.Draw("SAME")
    lsbdatfit.Draw("SAME")
    hsbdatfit.Draw("SAME")
    hdatsb.SetMarkerStyle(8)
    hdatsb.SetMarkerSize(0.5)
    hdatsb.SetMarkerColor(ROOT.kBlack)
    hdatsb.Draw("SAME")
    stackleg.AddEntry(sbdatfit,"Fit to Data SB","l")
    stackleg.AddEntry(totnormfit,"Data SB Extrap - Likelihood","l")
    stackleg.AddEntry(sbdatfit,"Data SB","ep")
    stackleg.Draw()
    unnormlabel.Draw()
    brl.Draw()
    brh.Draw()
    srl.Draw()
    srlabel.Draw()
    zrlabel.Draw()
    if validation:
        vrl.Draw()
        vrlabel.Draw()
    pd11.Update()
    tc1.cd()
    tc1.Update()

    pd12.Draw()
    pd12.cd()
    hsbkgnorm.Draw('HIST')
    xax = hsbkgnorm.GetXaxis()
    yax = hsbkgnorm.GetYaxis()
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pd12,4,13)
    totnormfit.Draw("SAME")
    lsbdatfit.Draw("SAME")
    hsbdatfit.Draw("SAME")
    hdatsb.SetMarkerStyle(8)
    hdatsb.SetMarkerSize(0.5)
    hdatsb.SetMarkerColor(ROOT.kBlack)
    hdatsb.Draw("SAME")
    bkgfit.Draw("SAME")
    stackleg.Draw()
    normcomplabel.Draw()
    normlabel.Draw()
    brl.Draw()
    brh.Draw()
    srl.Draw()
    srlabel.Draw()
    zrlabel.Draw()
    if validation:
        vrl.Draw()    
        vrlabel.Draw()
    p12.Update()
    tc1.cd()
    tc1.Update()


    stackedfit = go.makeOutFile('Run2_2017_2018','norm_stackfit_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(stackedfit)

    
    #ttbarhist = go.makeOutFile('Run2_2017_2018','ttbar_hist','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    #rootfile = ROOT.TFile(ttbarhist,"recreate")
    #htrtt.Write()
    #rootfile.Close()

import argparse
import ROOT
import glob
import os
import gecorg as go
import numpy as np
import pandas as pd
import configparser

def makeAddedHist(s17,s18,xspairs,hsb):
    s17.sort(key = go.orderDY)
    s18.sort(key = go.orderDY)
    
    print("infunction")
    for i,f in enumerate(s17):
        tf = ROOT.TFile(f)
        numevents = float(str(tf.Get('hnevents').GetString()))
        xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
        scale = go.findScale(numevents,xs,41.53)
        h = tf.Get('h_zp_jigm')
        h.Scale(scale)
        hsb.Add(h)

    for i,f in enumerate(s18):
        tf = ROOT.TFile(f)
        numevents = float(str(tf.Get('hnevents').GetString()))
        xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
        scale = go.findScale(numevents,xs,59.74)
        h = tf.Get('h_zp_jigm')
        h.Scale(scale)
        hsb.Add(h)

    return hsb

if __name__=='__main__':

    #make the output
    tc = ROOT.TCanvas("tc","alphafits",1100,400)
    tc.Divide(3,1)
    
    #will replace with command line options
    bkg_dir = 'analysis_output_ZpAnomalon/2021-03-16/'
    zptcut  = '200.0'
    hptcut  = '300.0'
    metcut  = '300.0'
    btagwp  = '0.8'
    
    #gather the MC files
    f17dyjetsb = glob.glob(str(bkg_dir)+'/Fall17.DYJetsToLL_M-50_HT*_upout_DeepMassDecorrelTagZHbbvsQCD_sideband_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    f17dyjetsr = glob.glob(str(bkg_dir)+'/Fall17.DYJetsToLL_M-50_HT*_upout_DeepMassDecorrelTagZHbbvsQCD_signalr_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    a18dyjetsb = glob.glob(str(bkg_dir)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_DeepMassDecorrelTagZHbbvsQCD_sideband_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    a18dyjetsr = glob.glob(str(bkg_dir)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_DeepMassDecorrelTagZHbbvsQCD_signalr_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    
    #scale and add together MC
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
    config.read_file(fp)
    xspairs = config.items("DYJetsToLL")

    tf1 = ROOT.TFile(f17dyjetsb[0])
    hsb = tf1.Get('h_zp_jigm')
    hsb.Reset("ICESM")#creates an empty hist with same structure
    hsb = makeAddedHist(f17dyjetsb,a18dyjetsb,xspairs,hsb)
    
    tf2 = ROOT.TFile(f17dyjetsr[0])
    hsr = tf2.Get('h_zp_jigm')
    hsr.Reset("ICESM")
    hsr = makeAddedHist(f17dyjetsr,a18dyjetsr,xspairs,hsr)
    
    #crysb = ROOT.TF1("crysb",ROOT.Math.CrystalBall())
    
    
    #Draw
    ROOT.gStyle.SetOptFit(1011)
    tc.cd(1)
    hsb.Fit("landau","WLR+","",900,5000)
    lsb = hsb.GetFunction("landau")
    tc.cd(2)
    hsr.Fit("landau","WLR+","",900,5000)
    lsr = hsr.GetFunction("landau")

    fsb = lsb.Clone("fsb")
    fsr = lsr.Clone("fsr")
    alpha = ROOT.TF1("alpha","fsr/fsb",900,500)

    #ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    #ROOT.gSystem.Load("cfunctions/alphafits_C")
    #alpha = ROOT.TF1("alpha",ROOT.function_divide,900,5000,0)
    #print(type(alpha))
    
    #tc.cd(3)
    #alpha.Draw()
    #fsr.Draw("SAME")
    #alpha.Draw("SAME")


    
    figname = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figname)

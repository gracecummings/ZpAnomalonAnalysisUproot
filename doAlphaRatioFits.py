import argparse
import ROOT
import glob
import os
import gecorg as go
import numpy as np
import pandas as pd
import configparser

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

#def calcAddedHistErrors(samp,xspairs):
#    samp.sort(key = go.

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

    #gather errors
    f17dyjetsberrs = glob.glob(str(bkg_dir)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_DeepMassDecorrelTagZHbbvsQCD_sideband_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    f17dyjetsrerrs = glob.glob(str(bkg_dir)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_DeepMassDecorrelTagZHbbvsQCD_signalr_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    a18dyjetsberrs = glob.glob(str(bkg_dir)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_DeepMassDecorrelTagZHbbvsQCD_sideband_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    a18dyjetsrerrs = glob.glob(str(bkg_dir)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_DeepMassDecorrelTagZHbbvsQCD_signalr_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    
    #scale and add together MC
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
    config.read_file(fp)
    xspairs = config.items("DYJetsToLL")

    tf1 = ROOT.TFile(f17dyjetsb[0])
    hsb = tf1.Get('h_zp_jigm')
    hsb.Reset("ICESM")#creates an empty hist with same structure
    hsb = makeAddedHist(f17dyjetsb,a18dyjetsb,f17dyjetsberrs,a18dyjetsberrs,xspairs,hsb)
    
    tf2 = ROOT.TFile(f17dyjetsr[0])
    hsr = tf2.Get('h_zp_jigm')
    hsr.Reset("ICESM")
    hsr = makeAddedHist(f17dyjetsr,a18dyjetsr,f17dyjetsrerrs,a18dyjetsrerrs,xspairs,hsr)
    
    ROOT.gSystem.CompileMacro("cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("cfunctions/alphafits_C")
    hsrt = hsr.Clone()
    hsbt = hsb.Clone()
    
    #Draw
    ROOT.gStyle.SetOptFit(1011)
    tc.cd(1)
    sbfit = ROOT.landauFit(hsbt,"sbl")
    hsb.Draw()
    sbfit.Draw("SAME")
    tc.cd(2)
    srfit = ROOT.landauFit(hsrt,"srl")
    hsr.Draw()
    srfit.Draw("SAME")
    #srfit.Draw()
    tc.cd(3)
    alpha = ROOT.alphaRatioMaker(hsbt,hsrt)
    #print(alpha)
    alpha.Draw()

    figname = go.makeOutFile('Run2_2017_2018','alpha_fits','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figname)

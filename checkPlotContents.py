import ROOT
import gecorg as go
import numpy as np
import pandas as pd

if __name__=='__main__':

    #will replace with command line options
    #path    = 'analysis_output_ZpAnomalon/2021-05-18/'#reclustering
    path    = 'analysis_output_ZpAnomalon/2021-05-30/'#no reclustering
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'

    bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    sigs = go.signal(path,zptcut,hptcut,metcut,btagwp,10,101.27)

    #tf1 = ROOT.TFile(sigs.sigsr[0])
    #empty = tf1.Get('h_zp_jigm')
    #empty.Reset("ICESM")#creates an empty hist with same structure
    pmap = np.zeros((len(sigs.sigsr),9))
    for i,sig in enumerate(sigs.sigsr):
        mzp,mnd,mns = go.massPoints(sig.replace(path+"ZpAnomalonHZ_UFO-",""))
        zpnddiff = float(mzp)-2*float(mnd)
        ndnsdiff = float(mnd)-float(mns)
        tf = ROOT.TFile(sig)
        hdr = tf.Get('h_dr_zh')
        yieldl = hdr.Integral(0,4)
        yieldh = hdr.Integral(5,hdr.GetNbinsX())

        pmap[i,:] = [float(mzp),float(mnd),float(mns),zpnddiff,ndnsdiff,yieldl,yieldh,round(yieldl/25000*100,2),round(yieldh/25000*100,2)]

    print("mZp   mnd   mns   zpnddiff     ndnsdiff    overlap    nonoverlap    overlapeff    nonovereff")
    print(pmap)

    #for row in pmap:

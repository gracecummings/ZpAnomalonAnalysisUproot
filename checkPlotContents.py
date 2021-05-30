import ROOT
import gecorg as go
import numpy as np
import pandas as pd

if __name__=='__main__':

    #will replace with command line options
    #pathrec    = 'analysis_output_ZpAnomalon/2021-05-18/'#reclustering sr only
    pathrec    = 'analysis_output_ZpAnomalon/2021-05-31/'#reclustering total r
    #pathori    = 'analysis_output_ZpAnomalon/2021-05-30/'#no reclustering sr only
    pathori    = 'analysis_output_ZpAnomalon/2021-05-31_fullregion_noreclustering/'#no reclustering totalr
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'

    #bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    sigsrec = go.signal(pathrec,zptcut,hptcut,metcut,btagwp,10,101.27)
    sigsori = go.signal(pathori,zptcut,hptcut,metcut,btagwp,10,101.27)

    #Get deltaR yeilds
    #sigs = sigsori
    #pmap = np.zeros((len(sigs.sigsr),9))
    #for i,sig in enumerate(sigs.sigsr):
    #    mzp,mnd,mns = go.massPoints(sig.replace(path+"ZpAnomalonHZ_UFO-",""))
    #    zpnddiff = float(mzp)-2*float(mnd)
    #    ndnsdiff = float(mnd)-float(mns)
    #    tf = ROOT.TFile(sig)
    #    hdr = tf.Get('h_dr_zh')
    #    yieldl = hdr.Integral(0,4)
    #    yieldh = hdr.Integral(5,hdr.GetNbinsX())

    #    pmap[i,:] = [float(mzp),float(mnd),float(mns),zpnddiff,ndnsdiff,yieldl,yieldh,round(yieldl/25000*100,2),round(yieldh/25000*100,2)]

    #print("mZp   mnd   mns   zpnddiff     ndnsdiff    overlap    nonoverlap    overlapeff    nonovereff")
    #print(pmap)

    #for i,sig in enumerate(sigsrec.sigsr):
    for i,sig in enumerate(sigsrec.sigfl):
        #print(sig)
        #print(sigsori.sigsr[i])
        mzp,mnd,mns = go.massPoints(sig.replace(pathrec+"ZpAnomalonHZ_UFO-",""))
        tfrec = ROOT.TFile(sig)
        #tfori = ROOT.TFile(sigsori.sigsr[i])####
        tfori = ROOT.TFile(sigsori.sigfl[i])####
        hsdmrec = tfrec.Get('h_h_sd')
        hsdmori = tfori.Get('h_h_sd')
        hptrec  = tfrec.Get('h_h_pt')
        hptori  = tfori.Get('h_h_pt')
        maxsd = hsdmrec.GetMaximum()
        maxpt = hptrec.GetMaximum()
        hsdmori.SetMaximum(maxsd+50)
        hptori.SetMaximum(maxpt+50)
        hsdmrec.SetLineColor(ROOT.kRed)
        hptrec.SetLineColor(ROOT.kRed)
        hsdmori.SetStats(0)
        hptori.SetStats(0)
        hsdmrec.SetStats(0)
        hptrec.SetStats(0)
        
        tc = ROOT.TCanvas("tc","comp",1200,600)
        legsd = ROOT.TLegend(0.50,0.50,0.88,0.88)
        legpt = ROOT.TLegend(0.50,0.50,0.88,0.88)
        
        legsd.AddEntry(hsdmori,"No reclustering")
        legsd.AddEntry(hptrec,"w/ reclustering")
        legpt.AddEntry(hptori,"No reclustering")
        legpt.AddEntry(hptrec,"w/ reclustering")
        
        tc.Divide(2,1)
        tc.cd(1)
        hsdmori.Draw("HIST")
        hsdmrec.Draw("HISTSAME")
        legsd.Draw()
        tc.cd(2)
        hptori.Draw("HIST")
        hptrec.Draw("HISTSAME")
        legpt.Draw()
        tc.cd()
        tc.Update()

        pngname = go.makeOutFile("Zp"+str(mzp)+"_ND"+str(mnd)+"_NS"+str(mns),"recluscomp_totalr",'.png',zptcut,hptcut,metcut,btagwp)
        tc.SaveAs(pngname)
        

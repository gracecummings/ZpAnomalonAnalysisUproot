import ROOT
import gecorg as go
import numpy as np
import pandas as pd
import glob

if __name__=='__main__':

    #will replace with command line options
    pathrec    = 'analysis_output_ZpAnomalon/2021-06-07_reclusteredJets/'#reclustering sr only
    #pathrec    = 'analysis_output_ZpAnomalon/2021-05-31_totalr_reclustering/'#reclustering total r
    pathori    = 'analysis_output_ZpAnomalon/2021-06-07_nonreclusteredJets/'#no reclustering sr only
    #pathori    = 'analysis_output_ZpAnomalon/2021-05-31_fullregion_noreclustering/'#no reclustering totalr
    pathgen    = 'analysis_output_ZpAnomalon/2021-06-07_genLevelInfo/'
    #pathrec = 'analysis_output_ZpAnomalon/2021-06-01_totalr_reclustered_0sdm_0btag/'
    #pathori = 'analysis_output_ZpAnomalon/2021-06-01_totalr_orignal_0sdm_0btag/'
    
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    #btagwp  = '0.0'

    #bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    sigsrec = go.signal(pathrec,zptcut,hptcut,metcut,btagwp,10,101.27)
    sigsori = go.signal(pathori,zptcut,hptcut,metcut,btagwp,10,101.27)
    siggen  = glob.glob(pathgen+'ZpAnomalonHZ*upout_GENONLY*.root')
    
    #Get deltaR yeilds
    sigs = sigsrec
    path = pathrec
    pmap = np.zeros((len(sigs.sigsr),5))
    #for i,sig in enumerate(sigs.sigsr):
    for i,sig in enumerate(siggen):
        #mzp,mnd,mns = go.massPoints(sig.replace(path+"ZpAnomalonHZ_UFO-",""))
        part = sig.replace(pathgen+"ZpAnomalonHZ_UFO-","").split("_upout")[0]
        mzp,mnd,mns = go.massPoints(part)
        #print(mzp,mnd,mns)
        zpnddiff = float(mzp)-2*float(mnd)
        ndnsdiff = float(mnd)-float(mns)
        tf = ROOT.TFile(sig)
        hdr = tf.Get('h_dr_zh')
        yieldl = hdr.Integral(0,4)
        yieldh = hdr.Integral(5,hdr.GetNbinsX())

        pmap[i,:] = [float(mzp),float(mnd),float(mns),yieldl,yieldh]

    #print("mZp   mnd   mns   zpnddiff     ndnsdiff    overlap    nonoverlap    overlapeff    nonovereff")
    print(pmap)

    #for i,sig in enumerate(sigsrec.sigsr):
    #for i,sig in enumerate(sigsrec.sigfl):
        #print(sig)
        #print(sigsori.sigsr[i])
     #   mzp,mnd,mns = go.massPoints(sig.replace(pathrec+"ZpAnomalonHZ_UFO-",""))
     #   tfrec = ROOT.TFile(sig)
     #   tfori = ROOT.TFile(sigsori.sigsr[i])####
     #   #tfori = ROOT.TFile(sigsori.sigfl[i])####
     #   hsdmrec = tfrec.Get('h_h_sd')
     #   hsdmori = tfori.Get('h_h_sd')
     #   hptrec  = tfrec.Get('h_h_pt')
     #   hptori  = tfori.Get('h_h_pt')
     #   maxsd = hsdmrec.GetMaximum()
     #   maxpt = hptrec.GetMaximum()
     #   hsdmori.SetMaximum(maxsd+50)
     #   hptori.SetMaximum(maxpt+50)
     #   hsdmrec.SetLineColor(ROOT.kRed)
     #   hptrec.SetLineColor(ROOT.kRed)
     #   hsdmori.SetStats(0)
     #   hptori.SetStats(0)
     #   hsdmrec.SetStats(0)
     #   hptrec.SetStats(0)
     #   hsdmori.GetXaxis().SetTitle("h cand soft drop mass")
     #   hptori.GetXaxis().SetTitle("h cand p_{T}")
     #   hsdmori.GetYaxis().SetTitle("Events")
     #   hptori.GetYaxis().SetTitle("Events")#

#        hsdmrec.GetXaxis().SetTitle("h cand soft drop mass")
#        hptrec.GetXaxis().SetTitle("h cand p_{T}")
#        hsdmrec.GetYaxis().SetTitle("Events")
#        hptrec.GetYaxis().SetTitle("Events")##

#        tc = ROOT.TCanvas("tc","comp",1200,600)
#        legsd = ROOT.TLegend(0.50,0.75,0.88,0.88)
#        legpt = ROOT.TLegend(0.50,0.75,0.88,0.88)
#        
#        legsd.AddEntry(hsdmori,"No reclustering")
#        legsd.AddEntry(hptrec,"w/ reclustering")
#        legpt.AddEntry(hptori,"No reclustering")
#        legpt.AddEntry(hptrec,"w/ reclustering")
 #       
 #       tc.Divide(2,1)
 #       tc.cd(1)
 #       hsdmori.Draw("HIST")
 #       hsdmrec.Draw("HISTSAME")
 #       legsd.Draw()
 #       tc.cd(2)
 #       hptori.Draw("HIST")
 #       hptrec.Draw("HISTSAME")
 #       legpt.Draw()
 #       tc.cd()
 #       tc.Update()#

#        pngname = go.makeOutFile("Zp"+str(mzp)+"_ND"+str(mnd)+"_NS"+str(mns),"recluscomp_totalr",'.png',zptcut,hptcut,metcut,btagwp)
 #       tc.SaveAs(pngname)
        

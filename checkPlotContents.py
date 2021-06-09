import ROOT
import gecorg as go
import numpy as np
import pandas as pd
import glob

def scaleFinder(pair):
    #pair = [hsdmrec,hsdmori]
    maxs = [x.GetMaximum() for x in pair]
    maxidx  = maxs.index(max(maxs))
    minidx  = maxs.index(min(maxs))
    scale = maxs[maxidx]/maxs[minidx]
    return scale,maxidx,max(maxs)
    
if __name__=='__main__':

    #will replace with command line options
    pathrec    = 'analysis_output_ZpAnomalon/2021-06-08_reclusteredJets_skimlevel/'
    #pathrec    = 'analysis_output_ZpAnomalon/2021-05-31_totalr_reclustering/'
    pathori    = 'analysis_output_ZpAnomalon/2021-06-08_nonreclusteredJets_skimlevel/'
    #pathori    = 'analysis_output_ZpAnomalon/2021-05-31_fullregion_noreclustering/'#no reclustering totalr
    pathgen    = 'analysis_output_ZpAnomalon/2021-06-07_genLevelInfo/'
    #pathrec = 'analysis_output_ZpAnomalon/2021-06-01_totalr_reclustered_0sdm_0btag/'
    #pathori = 'analysis_output_ZpAnomalon/2021-06-01_totalr_orignal_0sdm_0btag/'
    
    zptcut  = '150.0'
    hptcut  = '0.0'
    metcut  = '0.0'
    btagwp  = '0.8'


    #bkgs = go.backgrounds(path,zptcut,hptcut,metcut,btagwp)
    sigsrec = go.signal(pathrec,zptcut,hptcut,metcut,btagwp,10,101.27)
    sigsori = go.signal(pathori,zptcut,hptcut,metcut,btagwp,10,101.27)
    #siggen  = glob.glob(pathgen+'ZpAnomalonHZ*upout_GENONLY*.root')
    
    #Get deltaR yeilds
    #sigs = sigsori
    #path = pathori
    #pmap = np.zeros((len(sigs.sigsr),5))
    #for i,sig in enumerate(sigs.sigsr):
    #for i,sig in enumerate(siggen):
    #    mzp,mnd,mns = go.massPoints(sig.replace(path+"ZpAnomalonHZ_UFO-",""))
        #part = sig.replace(pathgen+"ZpAnomalonHZ_UFO-","").split("_upout")[0]
        #mzp,mnd,mns = go.massPoints(part)
        #print(mzp,mnd,mns)
    #    zpnddiff = float(mzp)-2*float(mnd)
    #    ndnsdiff = float(mnd)-float(mns)
    #    tf = ROOT.TFile(sig)
    #    hdr = tf.Get('h_dr_zh')
    #    yieldl = hdr.Integral(0,4)
    #    yieldh = hdr.Integral(5,hdr.GetNbinsX())#

    #    pmap[i,:] = [float(mzp),float(mnd),float(mns),yieldl,yieldh]

    #print("mZp   mnd   mns   zpnddiff     ndnsdiff    overlap    nonoverlap    overlapeff    nonovereff")
    #print(pmap)

    for i,sig in enumerate(sigsrec.sigsr):
    #for i,sig in enumerate(sigsrec.sigfl):
        print(sig)
        print(sigsori.sigsr[i])
        mzp,mnd,mns = go.massPoints(sig.replace(pathrec+"ZpAnomalonHZ_UFO-",""))
        tfrec = ROOT.TFile(sig)
        tfori = ROOT.TFile(sigsori.sigsr[i])####
        #tfori = ROOT.TFile(sigsori.sigfl[i])####
        hsdmrec = tfrec.Get('h_h_sd')
        hsdmori = tfori.Get('h_h_sd')
        hptrec  = tfrec.Get('h_h_pt')
        hptori  = tfori.Get('h_h_pt')
        hphirec  = tfrec.Get('h_h_phiw')
        hphiori  = tfori.Get('h_h_phiw')
        hetarec  = tfrec.Get('h_h_eta')
        hetaori  = tfori.Get('h_h_eta')

        #Unscaled verions
        #maxsd = max(hsdmrec.GetMaximum(),hsdmori.GetMaximum())
        #maxpt = max(hptrec.GetMaximum(),hptori.GetMaximum())
        #maxphi = max(hphirec.GetMaximum(),hphiori.GetMaximum())
        #maxeta = max(hetarec.GetMaximum(),hetaori.GetMaximum())
        #hsdmori.SetMaximum(maxsd+50)
        #hptori.SetMaximum(maxpt+50)
        #hphiori.SetMaximum(maxphi+50)
        #hetaori.SetMaximum(maxeta+50)

        hsdmrec.Scale(1/hsdmrec.Integral())
        hptrec.Scale(1/hptrec.Integral())
        hphirec.Scale(1/hphirec.Integral())
        hetarec.Scale(1/hetarec.Integral())
        hsdmori.Scale(1/hsdmori.Integral())
        hptori.Scale(1/hptori.Integral())
        hphiori.Scale(1/hphiori.Integral())
        hetaori.Scale(1/hetaori.Integral())

        #maxlambda = lambda x : x.GetMaximum()
        sdpair = [hsdmrec,hsdmori]
        sdmaxs  = [x.GetMaximum() for x in sdpair]
        maxidx  = sdmaxs.index(max(sdmaxs))
        minidx  = sdmaxs.index(min(sdmaxs))
        scale = sdmaxs[maxidx]/sdmaxs[minidx]
        sdpair[minidx].Scale(scale)

        ptpair = [hptrec,hptori]
        ptmaxs  = [x.GetMaximum() for x in ptpair]
        maxidx  = ptmaxs.index(max(ptmaxs))
        minidx  = ptmaxs.index(min(ptmaxs))
        scale = ptmaxs[maxidx]/ptmaxs[minidx]
        ptpair[minidx].Scale(scale)

        phipair = [hphirec,hphiori]
        phimaxs  = [x.GetMaximum() for x in phipair]
        maxidx  = phimaxs.index(max(phimaxs))
        minidx  = phimaxs.index(min(phimaxs))
        scale = phimaxs[maxidx]/phimaxs[minidx]
        phipair[minidx].Scale(scale)

        etapair = [hetarec,hetaori]
        etamaxs  = [x.GetMaximum() for x in etapair]
        maxidx  = etamaxs.index(max(etamaxs))
        minidx  = etamaxs.index(min(etamaxs))
        scale = etamaxs[maxidx]/etamaxs[minidx]
        etapair[minidx].Scale(scale)

        hsdmori.SetMaximum(max(sdmaxs)+max(sdmaxs)*.25)
        hptori.SetMaximum(max(ptmaxs)+max(ptmaxs)*.25)
        hphiori.SetMaximum(max(phimaxs)+max(phimaxs)*.25)
        hetaori.SetMaximum(max(etamaxs)+max(etamaxs)*.25)

        hsdmrec.SetLineColor(ROOT.kRed)
        hptrec.SetLineColor(ROOT.kRed)
        hphirec.SetLineColor(ROOT.kRed)
        hetarec.SetLineColor(ROOT.kRed)
        hsdmori.SetStats(0)
        hptori.SetStats(0)
        hphiori.SetStats(0)
        hetaori.SetStats(0)
        hsdmrec.SetStats(0)
        hptrec.SetStats(0)
        hphirec.SetStats(0)
        hetarec.SetStats(0)
        hsdmori.GetXaxis().SetTitle("h cand soft drop mass")
        hptori.GetXaxis().SetTitle("h cand \(p_{T}\)")
        hphiori.GetXaxis().SetTitle("h cand \(\phi\)")
        hetaori.GetXaxis().SetTitle("h cand \(\eta\)")
        hsdmori.GetYaxis().SetTitle("Events")
        hptori.GetYaxis().SetTitle("Events")
        hphiori.GetYaxis().SetTitle("Events")
        hetaori.GetYaxis().SetTitle("Events")

        hsdmrec.GetXaxis().SetTitle("h cand soft drop mass")
        hptrec.GetXaxis().SetTitle("h cand \(p_{T}\)")
        hphirec.GetXaxis().SetTitle("h cand \(\phi\)")
        hetarec.GetXaxis().SetTitle("h cand \(\eta\)")
        hsdmrec.GetYaxis().SetTitle("Events")
        hptrec.GetYaxis().SetTitle("Events")##
        hphirec.GetYaxis().SetTitle("Events")
        hetarec.GetYaxis().SetTitle("Events")

        tc = ROOT.TCanvas("tc","comp",1200,1200)
        legsd = ROOT.TLegend(0.50,0.75,0.88,0.88)
        legpt = ROOT.TLegend(0.50,0.75,0.88,0.88)
        legphi = ROOT.TLegend(0.50,0.75,0.88,0.88)
        legeta = ROOT.TLegend(0.50,0.75,0.88,0.88)
        
        legsd.AddEntry(hsdmori,"No reclustering")
        legsd.AddEntry(hptrec,"w/ reclustering")
        legpt.AddEntry(hptori,"No reclustering")
        legpt.AddEntry(hptrec,"w/ reclustering")
        legphi.AddEntry(hsdmori,"No reclustering")
        legphi.AddEntry(hptrec,"w/ reclustering")
        legeta.AddEntry(hptori,"No reclustering")
        legeta.AddEntry(hptrec,"w/ reclustering")
        
        tc.Divide(2,2)
        tc.cd(1)
        hsdmori.Draw("HIST")
        hsdmrec.Draw("HISTSAME")
        legsd.Draw()
        tc.cd(2)
        hptori.Draw("HIST")
        hptrec.Draw("HISTSAME")
        legpt.Draw()
        tc.cd(3)
        hphiori.Draw("HIST")
        hphirec.Draw("HISTSAME")
        legphi.Draw()
        tc.cd(4)
        hetaori.Draw("HIST")
        hetarec.Draw("HISTSAME")
        legeta.Draw()
        tc.cd()
        tc.Update()#

        pngname = go.makeOutFile("Zp"+str(mzp)+"_ND"+str(mnd)+"_NS"+str(mns),"recluscomp_totalr",'.png',zptcut,hptcut,metcut,btagwp)
        tc.SaveAs(pngname)
        

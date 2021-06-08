import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg as go


parser = argparse.ArgumentParser()

def boostUnc(values,weights,nbins,binstart,binstop):
    boosth = bh.Histogram(bh.axis.Regular(bins=nbins,start=binstart,stop=binstop),storage=bh.storage.Weight())
    boosth.fill(values,weight=weights)
    boostvar = boosth.view().variance
    boosterr = np.sqrt(boostvar)
    return boosterr

def wrapPhi(phi):
    if phi < 0:
        wphi = -1*phi
    else:
        wphi = phi
    return wphi

def wrapDeltaPhi(dphi):
    if dphi > 3.14159:
        dp = 2*3.14159 - dphi
    else:
        dp = dphi
    return dp
 

def deltaPhi(v1phi,v2phi):
    dp = abs(v1phi-v2phi)
    wdp = dp.map(wrapDeltaPhi)
    #This returns a df with the same number of events as input
    #If a cut is introduced, need to carry weight column
    return wdp

def deltaR(v1phi,v2phi,v1eta,v2eta):
    dR = ((v2phi-v1phi)**2+(v2eta-v1eta)**2)**(1/2)
    return dR

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100.0,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250.0,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50.0,help = "pT cut on met")
    parser.add_argument("-sdm","--sdmCut",type=float,default = 10.0,help = "lowest soft drop mass cut")
    parser.add_argument("-date","--date",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-c","--comboregion",type=bool,help = "do you want combined SR and SB?")
    args = parser.parse_args()

    samp   = args.sample
    sdmcut = args.sdmCut
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    metcut = args.metPtCut
    btaggr = args.btagger
    btagwp = args.btagWP
    sr     = args.signalregion
    comb   = args.comboregion


    inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*.root')
    print("    Doing selections on:")
    print(inputfiles[:1])
    stype,year = go.sampleType(samp)
    print(stype)
    if sr and stype != 0:
        print("    using signal region selections")
    elif comb and stype != 0:
        print("    using full region selections")
    else:
        print("    using sideband selections")
    
    branches = [#b'ZCandidate_*',
                b'ghCandidate_*',
                b'gzCandidate_*',
                #b'METclean',
                #b'METPhiclean',
                #b'ZPrime_mass_est',
                #b'ND_mass_est',
                #b'NS_mass_est',
                #b'event_weight',
                #b'LMuCandidate_*',
                #b'sLMuCandidate_*',
    ]

    
    #events = up3.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
    #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
    #print(events)
    
    #for b in events:
    #    fdf = b
        #print(type(b))
        #print(b.keys())
        #print(b)
        #print("Doing SD mass lower cut of :" ,sdmcut)
        #do some cuts
        #sddf   = b[b['hCandidate_sd'] > sdmcut]
        #metdf  = sddf[sddf['METclean'] > metcut]
        #zptdf  = metdf[metdf['ZCandidate_pt'] > zptcut]
        #hptdf  = zptdf[zptdf['hCandidate_pt'] > hptcut]
        #btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]
        
        #srup   = btdf[btdf['hCandidate_sd'] > 70.]
        #bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
        #srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
        #lowsb  = btdf[btdf['hCandidate_sd'] <= 70.]
        #highsb = btdf[btdf['hCandidate_sd'] >= 150.]
        #sbdf   = pd.concat([lowsb,highsb])

    tree = up3.open(inputfiles[0])["PreSelection;1"]
    fdf = tree.pandas.df(branches=branches)
        
    #This will have to come out of the loop if true iteration is added
    region = "sideband"
    if stype != 0:
        if sr:
            #fdf = srdf
            region = "signalr"
        elif comb:
            #fdf = btdf
            region = "totalr"
        else:
            #fdf = sbdf
            region = "sideband"
    else:
        fdf = sbdf

    print("    number of passing events ",len(fdf))
    #print(fdf)
    #print("number of btag passing events ",len(btdf))

    #calculated quantities
    deltaphizhdf   = deltaPhi(fdf['gzCandidate_phi'],fdf['ghCandidate_phi'])
    deltaRzhdf     = deltaR(fdf['gzCandidate_phi'],fdf['ghCandidate_phi'],fdf['gzCandidate_eta'],fdf['ghCandidate_eta'])

    #lets make some histograms.
    rootfilename  = go.makeOutFile(samp,'upout_GENONLY_'+region+'_'+btaggr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))#need to update for btagger
    rootOutFile   = up3.recreate(rootfilename,compression = None)

    #print(fdf['gzCandidate_pt'])
    
    rootOutFile["h_z_pt"]       = np.histogram(fdf['gzCandidate_pt'],bins=80,range=(0,800))
    rootOutFile["h_z_phi"]      = np.histogram(fdf['gzCandidate_phi'],bins=100,range=(-3.14159,3.14159))
    rootOutFile["h_z_phiw"]     = np.histogram(fdf['gzCandidate_phi'].map(wrapPhi),bins=100,range=(0,3.14159))#wrapped version of phi
    rootOutFile["h_z_eta"]      = np.histogram(fdf['gzCandidate_eta'],bins=100,range=(-5,5))
    rootOutFile["h_z_m"]        = np.histogram(fdf['gzCandidate_m'],bins=100,range=(40,140))
    rootOutFile["h_h_pt"]       = np.histogram(fdf['ghCandidate_pt'],bins=40,range=(200,1200))
    rootOutFile["h_h_phi"]      = np.histogram(fdf['ghCandidate_phi'],bins=100,range=(-3.14159,3.14159))
    rootOutFile["h_h_phiw"]     = np.histogram(fdf['ghCandidate_phi'].map(wrapPhi),bins=100,range=(0,3.14159))#wrapped version of phi
    rootOutFile["h_h_eta"]      = np.histogram(fdf['ghCandidate_eta'],bins=100,range=(-5,5))
    rootOutFile["h_h_m"]        = np.histogram(fdf['ghCandidate_m'],bins=80,range=(0,400))
    rootOutFile["h_dphi_zh"]    = np.histogram(deltaphizhdf,bins=100,range=(0,3.14159))
    rootOutFile["h_dr_zh"]      = np.histogram(deltaRzhdf,bins=30,range=(0,6))
    
    #Book Keeping
    f = up3.open(inputfiles[0])
    #np.save(npOutFile,np.array([f['hnorigevnts'].values[0]]))
    rootOutFile["hnevents"]      = str(f['hnorigevnts'].values[0])

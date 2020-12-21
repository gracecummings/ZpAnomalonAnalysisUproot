import uproot4 as up4
import uproot as up3
import pandas as pd
import numpy as np
import argparse
import glob
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-s","--isSig", type=bool,help="is this a signal sample?")
    parser.add_argument("-d","--isData", type=bool,help="is this data")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50,help = "pT cut on met")
    args = parser.parse_args()

    samp   = args.sample
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    metcut = args.metPtCut
    btaggr = args.btagger
    btagwp = args.btagWP
    isSig  = args.isSig
    isData = args.isData

    inputfiles = glob.glob('../RestFrames/analysis_output_ZpAnomalon/2020-12-16/'+samp+'*_topiary.root')

    branches = [b'ZCandidate_*',
                b'hCandidate_*',
                b'METclean',
                b'METPhiclean',
                b'ZPrime_mass_est',
                b'ND_mass_est',
                b'NS_mass_est',
    ]

    #Example of plots
    outFile = up3.recreate("upout.root",compression = None)
    
    #events = up3.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
    events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)#pandas dies with ZCandidates


    #print(jets.ls))
    
    for b in events:
        #print(type(b))
        print(b.keys())
        #print(b)

        #do some cuts
        hetadf = b[np.abs(b['hCandidate_eta']) < 2.4]
        metdf  = hetadf[hetadf['METclean'] > metcut]
        zptdf  = metdf[metdf['ZCandidate_pt'] > zptcut]
        hptdf  = zptdf[zptdf['hCandidate_pt'] > hptcut]
        #hbtgdf = hptdf[hptdf['hCandidate_'+btaggr] > btagwp]
        #print(hptdf)
        #fdf is always the last dataframe
        fdf = hptdf

    #lets make some histograms.
    outFile["h_z_pt"]    = np.histogram(fdf['ZCandidate_pt'],bins=80,range=(0,800))
    outFile["h_z_phi"]   = np.histogram(fdf['ZCandidate_phi'],bins=100,range=(0,3.14159))#needs to fit range
    outFile["h_z_eta"]   = np.histogram(fdf['ZCandidate_eta'],bins=100,range=(-5,5))
    outFile["h_z_m"]     = np.histogram(fdf['ZCandidate_m'],bins=40,range=(70,110))
    outFile["h_h_pt"]    = np.histogram(fdf['hCandidate_pt'],bins=40,range=(200,1200))
    outFile["h_h_phi"]   = np.histogram(fdf['hCandidate_phi'],bins=100,range=(0,3.14159))#needs to fit range
    outFile["h_h_eta"]   = np.histogram(fdf['hCandidate_eta'],bins=100,range=(-5,5))
    outFile["h_h_m"]     = np.histogram(fdf['hCandidate_m'],bins=80,range=(0,400))
    outFile["h_h_sd"]    = np.histogram(fdf['hCandidate_sd'],bins=80,range=(0,4000))
    outFile["h_met"]     = np.histogram(fdf['METclean'],bins=78,range=(50,2000))
    outFile["h_met_phi"] = np.histogram(fdf['METPhiclean'],bins=100,range=(0,3.14159))#needs to fit range
    outFile["h_zp_jigm"] = np.histogram(fdf['ZPrime_mass_est'],bins=100,range=(500,5000))
    outFile["h_nd_jigm"] = np.histogram(fdf['ND_mass_est'],bins=130,range=(0,1300))
    outFile["h_ns_jigm"] = np.histogram(fdf['NS_mass_est'],bins=200,range=(0,500))

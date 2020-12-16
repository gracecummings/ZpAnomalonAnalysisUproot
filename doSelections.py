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
    parser.add_argument("-s","--isSig", type=bool,help="is this a signal sample?")
    parser.add_argument("-d","--isData", type=bool,help="is this data")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250,help = "pT cut on h")
    args = parser.parse_args()

    samp   = args.sample
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    btagwp = args.btagWP
    isSig  = args.isSig
    isData = args.isData

    inputfiles = glob.glob('../RestFrames/analysis_output_ZpAnomalon/2020-12-16/'+samp+'*_topiary.root')

    branches = [b'ZCandidate_pt',
                b'ZCandidate_eta',
                b'ZCandidate_phi',
                b'ZCandidate_m',
                b'hCandidate_pt',
                b'hCandidate_eta',
                b'hCandidate_phi',
                b'hCandidate_m',
                b'hCandidate_sd',
                b'METclean',
    ]

    #events = up3.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
    events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)#pandas dies with ZCandidates


    #print(jets.ls))
    
    for b in events:
        print(type(b))
        print(b.keys())
        print(b)



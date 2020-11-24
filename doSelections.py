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

    inputfiles = glob.glob('../dataHandling/MCskims_redo/skimsoct28/'+samp+'*.root')

    branches = [b'JetsAK8Clean_ID',
                b'JetsAK8Clean_softDropMass',
                b'JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD',
                #b'JetsAK8Clean_DeepMassDecorrelTagbbvsLight',
                #b'JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD',
                #b'JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD',
                #b'JetsAK8Clean_DeepMassDecorrelTagZvsQCD',
                #b'JetsAK8Clean_DeepTagHbbvsQCD',
                #b'JetsAK8Clean_DeepTagZbbvsQCD',
                #b'JetsAK8Clean_DeepTagZvsQCD',
                #b'JetsAK8Clean_doubleBDiscriminator',
                #b'JetsAK8Clean_jecFactor',
                #b'JetsAK8Clean_jerFactor',
                #b'JetsAK8Clean_origIndex',
                #b'JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb',
                b'ZCandidates',]

    df = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)

    #print(jets.ls))
    
    for d in df:
        #print("I gave it something valid?")
        print(type(d))
        print(d.keys())
    

import ROOT
import glob
import argparse
import os
import sys
import gecorg as go
from datetime import date

parser = argparse.ArgumentParser()

if __name__=="__main__":
    parser.add_argument("-s","--sample",help="sample name")
    parser.add_argument("-y","--year",help="analysis year")
    args = parser.parse_args()
    samp = args.sample
    year = args.year
    samptype = -1

    samptype = go.sampleType(samp)
    if samptype < 0:
        print "You have a problem, we do not undertand the sample coding"

    origevnts = 0
    #if not isSig:
    if samptype != 1:
        inChain = ROOT.TChain("PreSelection")
        inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        for f in inputs:
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
    else:
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inChain.Add("../dataHandling/"+year+"/"+samp+"*.root")
        origevnts = inChain.GetEntries()

    outFile = go.makeOutFile(samp,'topiary')

    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    ROOT.gSystem.Load('TreeMakerTopiary_C')

    print "Making topiary of ",samp
    print "     Sample type ",samptype
    print "     Events in TChain: ",inChain.GetEntries()
    print ("     Original data set had {0} events.").format(origevnts)
    print "    Saving topiary in ",outFile

    
    topiary = ROOT.TreeMakerTopiary(inChain)
    topiary.Loop(outFile,origevnts,samptype)


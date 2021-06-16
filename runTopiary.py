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
    #year = args.year
    samptype = -1

    samptype,year = go.sampleType(samp)
    year = str(year)
    #samptype = 1
    #checkedyear = 17
    if samptype < 0:
        print("You have a problem, we do not undertand the sample coding")

    origevnts = 0
    
    if samptype != 1:
        inChain = ROOT.TChain("PreSelection")
        inputs  = glob.glob("../dataHandling/20"+year+"/"+samp+"*.root")
        
        for f in inputs:
            #print(f)
            #if "_chunk108__" in f:
            #    continue
            #if "_chunk100_" in f:
            #    continue
            #if "chunk478" in f:
            #    continue
            #if "chunk75" in f:
            #    continue
            #if "chunk483" in f:
            #    continue
            #if "chunk375" in f:
            #    continue
            #if "chunk39" in f:
            #    continue
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
    else:
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inChain.Add("../dataHandling/20"+year+"/"+samp+"*.root")
        origevnts = inChain.GetEntries()

    outFile = go.makeOutFile(samp,'topiary_mumu','.root','0.0','250.0','0.0','0.0')#Needs to become dynamic with cuts
    print( "Making topiary of ",samp)
    print("     Sample type ",samptype)
    print("     Sample Year ",year)
    print("     Events in TChain: ",inChain.GetEntries())
    print(("     Original data set had {0} events in type.").format(origevnts))
    print("    Saving topiary in ",outFile)


    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","kfc")
    ROOT.gSystem.Load('TreeMakerTopiary_C')
    topiary = ROOT.TreeMakerTopiary(inChain,samptype,int(year))
    topiary.Loop(outFile,origevnts,samptype,int(year))




    



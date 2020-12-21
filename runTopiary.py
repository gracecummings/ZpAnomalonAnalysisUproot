import ROOT
import glob
import argparse
import os
import sys
from datetime import date

parser = argparse.ArgumentParser()

if __name__=="__main__":
    parser.add_argument("-s","--sample",help="sample name")
    parser.add_argument("-y","--year",help="analysis year")
    #parser.add_argument("-A","--anomalon",type=bool,help="is this signal?")
    args = parser.parse_args()
    samp = args.sample
    year = args.year
    #isSig = args.anomalon
    samptype = -1

    #Make code for type of sample
    if "Run" in samp:
        samptype = 0
    if "ZpAnomalon" in samp:
        samptype = 1
        #isSig = True
    if "DYJetsToLL" in samp:
         samptype = 2
    if "TTTo" in samp:
        samptype = 3
    if "WZTo" in samp:
        samptype = 4
    if "ZZTo" in samp:
        samptype = 5
    else:
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

    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+samp+"_topiary.root"
        
    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    ROOT.gSystem.Load('TreeMakerTopiary_C')

    print "Making topiary of ",samp
    print ("    Original data set had {0} events, trimmed tree is smaller.").format(origevnts)
    print "    Saving topiary in ",outFile

    
    topiary = ROOT.TreeMakerTopiary(inChain)
    topiary.Loop(outFile,origevnts,samptype)


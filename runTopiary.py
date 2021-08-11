import ROOT
import glob
import argparse
import os
import sys
#import gecorg as go
from datetime import date

def channelEncoding(string):
    #channel code comes from a 3 bit binary
    #there are 3 possible canidates in an event
    #z->mumu, z->ee, z->emu
    #the number in decimal is this flag
    #ex. Z->mumu, no others: 100, or 4
    if "mumu" == string:
        channel = 4
        message = "Using the Z->mumu selections"
    elif "ee" == string:
        channel = 2
        message = "Using the Z->ee selections"
    elif "emu" == string:
        channel = 1
        message = "Usiing the ttbar background emu selections"
    else:
        channel = -1
    return channel,message

def sampleType(sampstring):
    #Make numerical code for type of sample
    if "Run" in sampstring:
        samptype = 0
    elif "ZpAnomalon" in sampstring:
        samptype = 1
        year = 17
    elif "DYJetsToLL" in sampstring:
         samptype = 2
    elif "TTTo" in sampstring:
        samptype = 3
    elif "WZTo" in sampstring:
        samptype = 4
    elif "ZZTo" in sampstring:
        samptype = 5
    else:
        samptype = -1

    if "2018" in sampstring:
        year = 18
    if "Autumn18" in sampstring:
        year = 18
    if "2017" in sampstring:
        year = 17
    if "Fall17" in sampstring:
        year = 17
    if "2016" in sampstring:
        year = 16

    return samptype,year

def makeOutFile(sampstring,descrip,ftype,zptcut,hptcut,metcut,btagwp):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+descrip+"_Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+ftype
    return outFile


parser = argparse.ArgumentParser()

if __name__=="__main__":
    parser.add_argument("-s","--sample",help="sample name")
    parser.add_argument("-c","--channel",help="string for channel: mumu,ee, or emu")
    args = parser.parse_args() 
    samp = args.sample
    samptype = -1

    #Check what you are working with
    samptype,checkedyear = sampleType(samp)
    channel,channelprint = channelEncoding(args.channel)
    year = "20"+str(checkedyear)

    #Do not start if you are not prepared
    if samptype < 0:
        print "You have a problem, we do not undertand the sample coding"
        sys.exit()
    if channel <= 0:
        print "You have a problem, no channel given"
        sys.exit()
    origevnts = 0

    #Prepare your TChain
    if samptype != 1 and ".root" not in samp:
        #for non-signal samples
        inChain = ROOT.TChain("PreSelection")
        inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        for f in inputs:
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
    elif ".root" in samp:
        #for debug, and ntuple in working directory
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inChain.Add(samp)
    else:
        #for signal
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inputs = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        inChain.Add("../dataHandling/"+year+"/"+samp+"*.root")
        origevnts = inChain.GetEntries()

    outFile = makeOutFile(samp,'topiary_'+args.channel,'.root','0.0','250.0','0.0','0.0')#Needs to become dynamic with cuts
    print "Making topiary of ",samp
    print "     Sample type ",samptype
    print "     Sample Year ",year
    print "    ",channelprint
    print "     Events in TChain: ",inChain.GetEntries()
    print ("     Original data set had {0} events in type.").format(origevnts)
    print "    Saving topiary in ",outFile


    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    #ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","kfc")
    ROOT.gSystem.Load('TreeMakerTopiary_C')
    topiary = ROOT.TreeMakerTopiary(inChain,samptype,checkedyear,channel)
    topiary.Loop(outFile,origevnts,samptype,checkedyear,channel)




    



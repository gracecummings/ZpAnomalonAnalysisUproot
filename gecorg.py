import os
import sys
from datetime import date

def sampleType(sampstring):
    #Make numerical code for type of sample
    if "Run" in sampstring:
        samptype = 0
    elif "ZpAnomalon" in sampstring:
        samptype = 1
        #isSig = True
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
    return samptype

def makeOutFile(sampstring,suffix):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+suffix+".root"
    return outFile



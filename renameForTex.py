import argparse
import glob
import gecorg as go
import subprocess

parser = argparse.ArgumentParser()

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
    args = parser.parse_args()

    samp   = args.sample
    sdmcut = args.sdmCut
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    metcut = args.metPtCut
    btaggr = args.btagger
    btagwp = args.btagWP

    inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.png')
    for f in inputfiles:
        protoname = f.split("_Zpt")[0]
        name = protoname.split(args.date+"/")[-1]
        newfilename  = go.makeOutFile(name,'tex','.png',str(int(zptcut)),str(int(hptcut)),str(int(metcut)),'loose')
        subprocess.run(["mv",f,newfilename])

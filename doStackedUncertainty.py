import numpy as np
import pandas as pd
import gecorg_py3 as go
import configparser
import glob


if __name__=='__main__':
    #will need to add the stuff for scale factor calc and so on, that will need to be propagated
    zptcut = 200.0
    hptcut = 250.0
    metcut = 250.0
    lumi   = 41.52
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcounts = go.gatherBkg('analysis_output_ZpAnomalon/2021-01-12','totalevents',zptcut,hptcut,metcut)
    bkgerrfs  = go.gatherBkg('analysis_output_ZpAnomalon/2021-01-12','selected_errors',zptcut,hptcut,metcut)

    #find luminosity scaling for each background
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('xsects_2017.ini')
    config.readfp(fp)
    bkgdfs = {}
    for name in bkgnames:
        if name == "DYJetsToLL":
            bkgcounts[name].sort(key=go.orderFall17DY)
            bkgerrfs[name].sort(key=go.orderFall17DY)
        bkgxspairs = config.items(name)
        bkgdfs[name] = []
        for b,bkgbin in enumerate(bkgerrfs[name]):
            binname = bkgxspairs[b][0]
            xs = bkgxspairs[b][1]
            df = pd.read_pickle(bkgbin)
            origevnts = np.load(bkgcounts[name][b])[0]
            scale = go.findScale(origevnts,float(xs),lumi)
            sdf = df*scale
            sqrddf = sdf**2
            bkgdfs[name].append(sqrddf)


    #Now here we have a dictionary, with a key for each squared background error
    #From here we can calculate the uncertainty on each background individually
    #We also can do the uncertainty of all backgrounds together
    
    uncsqdDYJetsdf = sum(bkgdfs["DYJetsToLL"])
    uncsqdTTdf     = sum(bkgdfs["TT"])
    uncsqdWZdf     = sum(bkgdfs["WZTo2L2Q"])
    uncsqdZZdf     = sum(bkgdfs["ZZTo2L2Q"])
    uncsqdAlldf    = uncsqdDYJetsdf+uncsqdTTdf+uncsqdWZdf+uncsqdZZdf
    uncAlldf       = uncsqdAlldf**(1/2)

    totalUncFile = go.makeOutFile('Fall17.AllZpAnomalonBkgs','unc','.pkl',str(zptcut),str(hptcut),str(metcut))
    uncAlldf.to_pickle(totalUncFile)


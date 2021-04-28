import argparse
import numpy as np
import pandas as pd
import gecorg as go
import configparser
import glob

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-y","--year", type=float,help = "year of samples eg. 2017 -> 17")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-date","--date", type=str,help = "date folder with files to use")
    parser.add_argument("-r","--region",help="region of phase space: totalr,sideband, or signalr")
    args = parser.parse_args()

    zptcut = args.zptcut
    hptcut = args.hptcut
    metcut = args.metcut
    btagwp = args.btagwp
    year   = args.year
    reg    = args.region
    
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcounts = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'totalevents_'+reg,zptcut,hptcut,metcut,btagwp,year)
    bkgerrfs  = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'selected_errors_'+reg,zptcut,hptcut,metcut,btagwp,year)

    print("calculaing statistical uncertainties")    
    #find luminosity scaling for each background
    config = configparser.RawConfigParser()
    config.optionxform = str
    if year == 17:
        fp = open('xsects_2017.ini')
        mcprefix  = 'Fall17'
        datprefix = 'Run2017'
        lumi = 41.53
    if year == 18:
        fp = open('xsects_2018.ini')
        mcprefix  = 'Autumn18'
        datprefix = 'Run2018'
        lumi = 59.74

        
    config.read_file(fp)
    bkgdfs = {}
    for n,name in enumerate(bkgnames):
        if name == "DYJetsToLL":
            bkgcounts[n].sort(key=go.orderDY)
            bkgerrfs[n].sort(key=go.orderDY)
        bkgxspairs = config.items(name)
        bkgdfs[name] = []
        for b,bkgbin in enumerate(bkgerrfs[n]):
            binname = bkgxspairs[b][0]
            xs = bkgxspairs[b][1]
            df = pd.read_pickle(bkgbin)
            origevnts = np.load(bkgcounts[n][b])[0]
            scale = go.findScale(origevnts,float(xs),lumi)
            sdf = df*scale
            sqrddf = sdf**2
            bkgdfs[name].append(sqrddf)

    dy   = False
    tt   = False
    zz   = False
    wz   = False
    data = False
    
    #Now here we have a dictionary, with
    #a key for each squared background error
    #From here we can calculate the uncertainty on each background individually
    #We also can do the uncertainty of all backgrounds together

    if len(bkgdfs["DYJetsToLL"]) > 0:
        dy = True
        uncsqdDYJetsdf = sum(bkgdfs["DYJetsToLL"])
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)
        saveDYJetsunc  = go.makeOutFile(mcprefix+'.DYJetsToLL','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        uncDYJetsdf.to_pickle(saveDYJetsunc)

    if len(bkgdfs["TT"]) > 0:
        tt = True
        uncsqdTTdf     = sum(bkgdfs["TT"])
        uncTTdf        = uncsqdTTdf**(1/2)
        saveTTunc      = go.makeOutFile(mcprefix+'.TT','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        uncTTdf.to_pickle(saveTTunc)

    if len(bkgdfs["WZTo2L2Q"]) > 0:
        wz = True
        uncsqdWZdf     = sum(bkgdfs["WZTo2L2Q"])
        uncWZdf        = uncsqdWZdf**(1/2)
        saveWZ2L2Qunc  = go.makeOutFile(mcprefix+'.WZ2L2Q','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        uncWZdf.to_pickle(saveWZ2L2Qunc)

    if len(bkgdfs["ZZTo2L2Q"]) > 0:
        zz = True
        uncsqdZZdf     = sum(bkgdfs["ZZTo2L2Q"])
        uncZZdf        = uncsqdZZdf**(1/2)
        saveZZ2L2Qunc  = go.makeOutFile(mcprefix+'.ZZ2L2Q','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        uncZZdf.to_pickle(saveZZ2L2Qunc)
    
    #Combine all for stacks
    if (dy and tt and wz and zz):
        print("All background samples accounted for, doing total uncertainty stacking for ",year)
        uncsqdAlldf    = uncsqdDYJetsdf+uncsqdTTdf+uncsqdWZdf+uncsqdZZdf
        uncAlldf       = uncsqdAlldf**(1/2)
        totalUncFile = go.makeOutFile(mcprefix+'.AllZpAnomalonBkgs','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        uncAlldf.to_pickle(totalUncFile)

    daterrfs = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+datprefix+'*selected_errors_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
    datadfs = []
    if len(daterrfs) > 0:
        for d  in daterrfs:
            babyname = d.split('.SingleMuon')[0]
            name     = babyname.split('/')[-1]
            df = pd.read_pickle(d)
            sqrddf = df**2
            #datadfs[name] = sqrddf
            datadfs.append(sqrddf)

        datuncsum = sum(datadfs)
        datuncall = datuncsum**(1/2)

        datFileName = go.makeOutFile(datprefix+'.AllZpAnomalonData','unc_'+reg,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        datuncall.to_pickle(datFileName)



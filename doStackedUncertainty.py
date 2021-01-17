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

    #will only work once move stacker to python3
    #totalUncFile = go.makeOutFile('Fall17.AllZpAnomalonBkgs','unc','.pkl',str(zptcut),str(hptcut),str(metcut))
    #uncAlldf.to_pickle(totalUncFile)

    #saving as a numpy zip file hard coded right now because no mistake goes unpunished.
    #stacker just needs to be moved to python three to fix this mess
    #THIS IS A HACK AND I HATE IT
    fileName = go.makeOutFile('Fall17.AllZpAnomalonBkgs','unc','.zpz',str(zptcut),str(hptcut),str(metcut))
    npF = open(fileName,'wb')
    np.savez(npF,
             h_z_pt  = uncAlldf['h_z_pt'].values,
             h_z_eta = uncAlldf['h_z_eta'].values,
             h_z_m   = uncAlldf['h_z_m'].values,
             h_h_pt  = uncAlldf['h_h_pt'].values,
             h_h_eta = uncAlldf['h_h_eta'].values,
             h_h_m   = uncAlldf['h_h_m'].values,
             h_h_sd  = uncAlldf['h_h_sd'].values,
             h_met   = uncAlldf['h_met'].values,
             h_zp_jigm = uncAlldf['h_zp_jigm'].values,
             h_nd_jigm = uncAlldf['h_nd_jigm'].values,
             h_ns_jigm = uncAlldf['h_ns_jigm'].values
             )


    #Now, do the same for data. This can clearly be combined into a function of somekind,
    #but we need results now!

    #datcountfs = glob.glob('analysis_output_ZpAnomalon/2021-01-14/Run2017*totalevents_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.npy')
    daterrfs = glob.glob('analysis_output_ZpAnomalon/2021-01-14/Run2017*selected_errors_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.pkl')
    datadfs = []
    for d  in daterrfs:
        babyname = d.split('.SingleMuon')[0]
        name     = babyname.split('/')[-1]
        df = pd.read_pickle(d)
        sqrddf = df**2
        #datadfs[name] = sqrddf
        datadfs.append(sqrddf)

    datuncsum = sum(datadfs)
    datuncall = datuncsum**(1/2)

    datFileName = go.makeOutFile('Fall17.AllZpAnomalonData','unc','.zpz',str(zptcut),str(hptcut),str(metcut))

    npdatF = open(datFileName,'wb')
    np.savez(npdatF,
             h_z_pt  = datuncall['h_z_pt'].values,
             h_z_eta = datuncall['h_z_eta'].values,
             h_z_m   = datuncall['h_z_m'].values,
             h_h_pt  = datuncall['h_h_pt'].values,
             h_h_eta = datuncall['h_h_eta'].values,
             h_h_m   = datuncall['h_h_m'].values,
             h_h_sd  = datuncall['h_h_sd'].values,
             h_met   = datuncall['h_met'].values,
             h_zp_jigm = datuncall['h_zp_jigm'].values,
             h_nd_jigm = datuncall['h_nd_jigm'].values,
             h_ns_jigm = datuncall['h_ns_jigm'].values
             )


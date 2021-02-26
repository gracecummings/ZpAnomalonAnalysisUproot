import argparse
import numpy as np
import pandas as pd
import gecorg_py3 as go
import configparser
import glob

def saveUncertainties(uncdf,filename):
    npF = open(filename,'wb')
    np.savez(npF,
             h_z_pt  = uncdf['h_z_pt'].values,
             h_z_eta = uncdf['h_z_eta'].values,
             h_z_m   = uncdf['h_z_m'].values,
             h_h_pt  = uncdf['h_h_pt'].values,
             h_h_eta = uncdf['h_h_eta'].values,
             h_h_m   = uncdf['h_h_m'].values,
             h_h_sd  = uncdf['h_h_sd'].values,
             h_met   = uncdf['h_met'].values,
             h_zp_jigm = uncdf['h_zp_jigm'].values,
             h_nd_jigm = uncdf['h_nd_jigm'].values,
             h_ns_jigm = uncdf['h_ns_jigm'].values,
             h_btag    = uncdf['h_btag'].values
             )


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-L","--lumi", type=float,default = 41.53, help = "integrated luminosity for scale in fb^-1")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-y","--year", type=float,help = "year of samples eg. 2017 -> 17")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-date","--date", type=str,help = "date folder with files to use")
    args = parser.parse_args()

    zptcut = args.zptcut
    hptcut = args.hptcut
    metcut = args.metcut
    btagwp = args.btagwp
    lumi   = args.lumi
    year   = args.year
    
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcounts = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'totalevents',zptcut,hptcut,metcut,btagwp,year)
    bkgerrfs  = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'selected_errors',zptcut,hptcut,metcut,btagwp,year)

    print("Calculating statistical uncertainties for stacked, weighted background and stacked data")    
    #find luminosity scaling for each background
    config = configparser.RawConfigParser()
    config.optionxform = str
    if year == 17:
        fp = open('xsects_2017.ini')
        mcprefix  = 'Fall17'
        datprefix = 'Run2017'
    if year == 18:
        fp = open('xsects_2018.ini')
        mcprefix  = 'Aumtumn18'
        datprefix = 'Run2018'
    config.readfp(fp)
    bkgdfs = {}
    for name in bkgnames:
        if name == "DYJetsToLL":
            bkgcounts[name].sort(key=go.orderDY)
            bkgerrfs[name].sort(key=go.orderDY)
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
    uncDYJetsdf    = uncsqdDYJetsdf**(1/2)
    saveDYJetsunc  = go.makeOutFile(mcprefix+'.DYJetsToLL','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveDYJetsuncCsv  = go.makeOutFile(mcprefix+'.DYJetsToLL','unc','.csv',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveUncertainties(uncDYJetsdf,saveDYJetsunc)
    uncDYJetsdf.to_csv()
    
    uncsqdTTdf     = sum(bkgdfs["TT"])
    uncTTdf        = uncsqdTTdf**(1/2)
    saveTTunc      = go.makeOutFile(mcprefix+'.TT','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveUncertainties(uncTTdf,saveTTunc)
    
    uncsqdWZdf     = sum(bkgdfs["WZTo2L2Q"])
    uncWZdf        = uncsqdWZdf**(1/2)
    saveWZ2L2Qunc  = go.makeOutFile(mcprefix+'.WZ2L2Q','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveUncertainties(uncWZdf,saveWZ2L2Qunc)
    
    uncsqdZZdf     = sum(bkgdfs["ZZTo2L2Q"])
    uncZZdf        = uncsqdZZdf**(1/2)
    saveZZ2L2Qunc  = go.makeOutFile(mcprefix+'.ZZ2L2Q','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveUncertainties(uncZZdf,saveWZ2L2Qunc)
    
    #Combine all for stacks
    uncsqdAlldf    = uncsqdDYJetsdf+uncsqdTTdf+uncsqdWZdf+uncsqdZZdf
    uncAlldf       = uncsqdAlldf**(1/2)
    saveAllBkgunc  = go.makeOutFile(mcprefix+'.AllZpAnomalonBkgs','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    saveUncertainties(uncAlldf,saveAllBkgunc)

    
    

    #will only work once move stacker to python3
    #totalUncFile = go.makeOutFile('Fall17.AllZpAnomalonBkgs','unc','.pkl',str(zptcut),str(hptcut),str(metcut))
    #uncAlldf.to_pickle(totalUncFile)

    #saving as a numpy zip file hard coded right now because no mistake goes unpunished.
    #stacker just needs to be moved to python three to fix this mess
    #THIS IS A HACK AND I HATE IT
    #fileName = go.makeOutFile(mcprefix+'.AllZpAnomalonBkgs','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    #npF = open(fileName,'wb')
    #np.savez(npF,
    #         h_z_pt  = uncAlldf['h_z_pt'].values,
    #         h_z_eta = uncAlldf['h_z_eta'].values,
    #         h_z_m   = uncAlldf['h_z_m'].values,
    #         h_h_pt  = uncAlldf['h_h_pt'].values,
    #         h_h_eta = uncAlldf['h_h_eta'].values,
    #         h_h_m   = uncAlldf['h_h_m'].values,
    #         h_h_sd  = uncAlldf['h_h_sd'].values,
    #         h_met   = uncAlldf['h_met'].values,
    #         h_zp_jigm = uncAlldf['h_zp_jigm'].values,
    #         h_nd_jigm = uncAlldf['h_nd_jigm'].values,
    #         h_ns_jigm = uncAlldf['h_ns_jigm'].values,
    #         h_btag    = uncAlldf['h_btag'].values
    #         )


    #Now, do the same for data. This can clearly be combined into a function of somekind,
    #but we need results now!

    datcountfs = glob.glob('analysis_output_ZpAnomalon/2021-01-14/Run2017*totalevents_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.npy')
    daterrfs = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+datprefix+'*selected_errors*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')

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

    datFileName = go.makeOutFile(datprefix+'.AllZpAnomalonData','unc','.npz',str(zptcut),str(hptcut),str(metcut),str(btagwp))

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
             h_ns_jigm = datuncall['h_ns_jigm'].values,
             h_btag    = datuncall['h_btag'].values
             )


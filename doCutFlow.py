import subprocess
import gecorg as go
import argparse
import ROOT
import numpy as np

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-L","--lumi", type=float,default = 41.53, help = "integrated luminosity for scale in fb^-1")
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-date","--date",help="date folder with plots to stack")
    parser.add_argument("-y","--year", type=float,help = "year of samples eg. 2017 -> 17")
    args = parser.parse_args()

    #Get command line parameters
    lumi          = args.lumi
    sig_xsec      = args.xsec
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut
    btagwp        = args.btagwp
    year          = args.year

    bkgfiles = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'upout',zptcut,hptcut,metcut,btagwp,year)
    #bkguncs  = np.load('analysis_output_ZpAnomalon/'+args.date+'/Fall17.AllZpAnomalonBkgs_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    bkginfo  = go.prepBkg(bkgfiles,bkgnames,bkgcols,'xsects_2017.ini',lumi,"yes")#gathers xs scales

    metstr = "Percent Passing \(MET > "+str(metcut)+"\)"
    zptstr = "Percent Passing \(Z p_{T} > "+str(zptcut)+"\)"
    hptstr = "Percent Passing Fat Jet \(p_{T} > "+str(hptcut)+"\)"
    btgstr = "Percent Passing {0} btag WP".format(str(btagwp))
    sbstr  = "Percent in side band"
    
    
    cfdict = {metstr:{},zptstr:{},hptstr:{},btgstr:{},sbstr:{}}
    totorig = 0
    totymet = 0
    totyzpt = 0
    totyhpt = 0
    totybtg = 0
    totysb  = 0
    
    for bkg in bkginfo:
        origevnts = 0
        yieldmet  = 0
        yieldzpt  = 0
        yieldhpt  = 0
        yieldbtag = 0
        yieldsb   = 0
        #cfdict
        for bkgbin in bkg["binlist"]:
            scale         = bkgbin["scale"]
            evntstot      = float(str(bkgbin["tfile"].Get('hnevents').GetString()))*scale
            evntspassmet  = float(str(bkgbin["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
            evntspasszpt  = float(str(bkgbin["tfile"].Get('hnevents_pZ').GetString()))*scale
            evntspasshpt  = float(str(bkgbin["tfile"].Get('hnevents_ph').GetString()))*scale
            evntspassbtag = float(str(bkgbin["tfile"].Get('hnevents_btag').GetString()))*scale
            evntspasssb   = float(str(bkgbin["tfile"].Get('hnevents_sb').GetString()))*scale
            origevnts += evntstot
            yieldmet  += evntspassmet
            yieldzpt  += evntspasszpt
            yieldhpt  += evntspasshpt
            yieldbtag += evntspassbtag
            yieldsb   += evntspasssb

        totorig += origevnts
        totymet += yieldmet
        totyzpt += yieldzpt
        totyhpt += yieldhpt
        totybtg += yieldbtag
        totysb  += yieldsb

        cfdict[metstr][bkg["name"]] = round(yieldmet/origevnts*100,4)
        cfdict[zptstr][bkg["name"]] = round(yieldzpt/origevnts*100,4)
        cfdict[hptstr][bkg["name"]] = round(yieldhpt/origevnts*100,4)
        cfdict[btgstr][bkg["name"]] = round(yieldbtag/origevnts*100,4)
        cfdict[sbstr][bkg["name"]] = round(yieldsb/origevnts*100,4)


    cfdict[metstr]["totbkg"] = round(totymet/totorig*100,4)
    cfdict[zptstr]["totbkg"] = round(totyzpt/totorig*100,4)
    cfdict[hptstr]["totbkg"] = round(totyhpt/totorig*100,4)
    cfdict[btgstr]["totbkg"] = round(totybtg/totorig*100,4)
    cfdict[sbstr]["totbkg"]  = round(totysb/totorig*100,4)

    #gets cuts in order of loosest to tightest based on yield
    cutordered = sorted(cfdict,key= lambda x : -1*cfdict[x]["totbkg"])

    cutFlowTableTex = go.makeOutFile("ZpAnomalon_mumu_bkg","cutflow",'.tex',str(int(zptcut)),str(int(hptcut)),str(int(metcut)),str(btagwp).split('.')[-1]+'E-10')
    cftab = open(cutFlowTableTex,"w")
    cftab.write(r'\begin{table}[htbp]')
    cftab.write('\n')
    cftab.write(r'\begin{center}')
    cftab.write('\n')
    cftab.write(r'\begin{tabular}{l | c | c | c | c | c | c }')#need to make dynamic
    cftab.write('\n')
    cftab.write("\hline\hline\n")
    cftab.write(r'cut description & ZZ & WZ & ttbar & DY+jets & Total \\')
    cftab.write('\n')
    cftab.write("\hline\n")
    for key in cutordered:
        #tabline = r''+key+' & '+str(cfdict[key]['ZZTo2L2Q'])+' & '+str(cfdict[key]['WZTo2L2Q'])+' & '+str(cfdict[key]['TT'])+' & '+str(cfdict[key]['DYJetsToLL'])+' & '+str(cfdict[key]['totbkg'])+' \\'
        tabline = key+' & '+str(cfdict[key]['ZZTo2L2Q'])+' & '+str(cfdict[key]['WZTo2L2Q'])+' & '+str(cfdict[key]['TT'])+' & '+str(cfdict[key]['DYJetsToLL'])+' & '+str(cfdict[key]['totbkg'])+' \\'
        cftab.write(tabline)
        cftab.write('\\')
        cftab.write('\n')
    cftab.write("\end{tabular}\n")
    cftab.write("\caption{Cut Flow based on MC weighted yields with luminosity scaling for 2017.}\n")
    cftab.write("\label{tab:cutflow}\n")
    cftab.write("\end{center}\n")
    cftab.write("\end{table}\n")
        
    

    

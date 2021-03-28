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

    #bkgupout17 = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'upout',zptcut,hptcut,metcut,btagwp,17)
    #bkgtopia17 = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'topiary',0.0,250.0,0.0,0.0,17)
    bkgupout17 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-26','upout',zptcut,hptcut,metcut,btagwp,17)
    bkgtopia17 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-26','topiary',0.0,250.0,0.0,0.0,17)
    #bkgupout18 = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'upout',zptcut,hptcut,metcut,btagwp,18)
    #bkgtopia18 = go.gatherBkg('analysis_output_ZpAnomalon/'+args.date,'topiary',0.0,250.0,0.0,0.0,18)
    
    #bkguncs  = np.load('analysis_output_ZpAnomalon/'+args.date+'/Fall17.AllZpAnomalonBkgs_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')
    
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    bkginfo17  = go.prepBkg(bkgupout17,bkgnames,bkgcols,'xsects_2017.ini',41.53,"yes")#gathers xs scales

    skimstr  = "Events with \( > 0\) fat jets and \( Z(\ell^{+}\ell^{-}) cand \)"
    trigstr  = "Events Pass trigger"
    zrecostr = "Events with \( 70 < m_{\ell\ell} < 110 \)"
    hrecostr = "Events Pass $H$ reco "
    metstr = "Events Passing \(MET > "+str(metcut)+"\)"
    zptstr = "Events Passing \(Z p_{T} > "+str(zptcut)+"\)"
    hptstr = "Events Passing Fat Jet \(p_{T} > "+str(hptcut)+"\)"
    btgstr = "Events Passing {0} btag WP".format(str(btagwp))
    sbstr  = "Events in side band"
    
    
    cfdict = {metstr:{},zptstr:{},hptstr:{},btgstr:{},sbstr:{},skimstr:{},trigstr:{},zrecostr:{},hrecostr:{}}
    totorig = 0
    totskim = 0
    tottrig = 0
    totzrec = 0
    tothrec = 0
    totymet = 0
    totyzpt = 0
    totyhpt = 0
    totybtg = 0
    totysb  = 0
    
    for bkg in bkginfo17:
        toplist = []
        for top in bkgtopia17:
            if bkg["name"] in top[0]:
                toplist = top
        origevnts = 0
        yieldskim = 0
        yieldtrig = 0
        yieldzrec = 0
        yieldhrec = 0
        yieldmet  = 0
        yieldzpt  = 0
        yieldhpt  = 0
        yieldbtag = 0
        yieldsb   = 0
        #cfdict
        for b,bkgbin in enumerate(bkg["binlist"]):
            reltop = 'str'
            for bkgtopia in toplist:
                if bkgbin["binname"] in bkgtopia:
                    reltop = bkgtopia
            topiary       = ROOT.TFile(reltop)        
            scale         = bkgbin["scale"]
            evntstot      = float(str(bkgbin["tfile"].Get('hnevents').GetString()))*scale
            evntspassmet  = float(str(bkgbin["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
            evntspasszpt  = float(str(bkgbin["tfile"].Get('hnevents_pZ').GetString()))*scale
            evntspasshpt  = float(str(bkgbin["tfile"].Get('hnevents_ph').GetString()))*scale
            evntspassbtag = float(str(bkgbin["tfile"].Get('hnevents_btag').GetString()))*scale
            evntspasssb   = float(str(bkgbin["tfile"].Get('hnevents_sb').GetString()))*scale
            evntsinskim   = topiary.Get('hnskimed').GetBinContent(1)*scale
            evntspasstrig = topiary.Get('htrigpass').GetBinContent(1)*scale
            evntspasszrec = topiary.Get('hZpass').GetBinContent(1)*scale
            evntspasshrec = topiary.Get('hHpass').GetBinContent(1)*scale
            
            origevnts += evntstot
            yieldskim += evntsinskim
            yieldtrig += evntspasstrig
            yieldzrec += evntspasszrec
            yieldhrec += evntspasshrec
            yieldmet  += evntspassmet
            yieldzpt  += evntspasszpt
            yieldhpt  += evntspasshpt
            yieldbtag += evntspassbtag
            yieldsb   += evntspasssb

        totorig += origevnts
        totskim += yieldskim
        tottrig += yieldtrig
        totzrec += yieldzrec
        tothrec += yieldhrec
        totymet += yieldmet
        totyzpt += yieldzpt
        totyhpt += yieldhpt
        totybtg += yieldbtag
        totysb  += yieldsb

        #cfdict[skimstr][bkg["name"]] = round(yieldskim/origevnts*100,4)
        #cfdict[trigstr][bkg["name"]] = round(yieldtrig/origevnts*100,4)
        #cfdict[zrecostr][bkg["name"]] = round(yieldzrec/origevnts*100,4)
        #cfdict[hrecostr][bkg["name"]] = round(yieldhrec/origevnts*100,4)
        #cfdict[metstr][bkg["name"]] = round(yieldmet/origevnts*100,4)
        #cfdict[zptstr][bkg["name"]] = round(yieldzpt/origevnts*100,4)
        #cfdict[hptstr][bkg["name"]] = round(yieldhpt/origevnts*100,4)
        #cfdict[btgstr][bkg["name"]] = round(yieldbtag/origevnts*100,4)
        #cfdict[sbstr][bkg["name"]] = round(yieldsb/origevnts*100,4)

        cfdict[skimstr][bkg["name"]] = round(yieldskim,2)
        cfdict[trigstr][bkg["name"]] = round(yieldtrig,2)
        cfdict[zrecostr][bkg["name"]] = round(yieldzrec,2)
        cfdict[hrecostr][bkg["name"]] = round(yieldhrec,2)
        cfdict[metstr][bkg["name"]] = round(yieldmet,2)
        cfdict[zptstr][bkg["name"]] = round(yieldzpt,2)
        cfdict[hptstr][bkg["name"]] = round(yieldhpt,2)
        cfdict[btgstr][bkg["name"]] = round(yieldbtag,2)
        cfdict[sbstr][bkg["name"]] = round(yieldsb,2)
        

    #cfdict[skimstr]["totbkg"] = round(totskim/origevnts*100,2)
    #cfdict[trigstr]["totbkg"] = round(tottrig/origevnts*100,2)
    #cfdict[zrecostr]["totbkg"] = round(totzrec/origevnts*100,2)
    #cfdict[hrecostr]["totbkg"] = round(tothrec/origevnts*100,2)
    #cfdict[metstr]["totbkg"] = round(totymet/totorig*100,2)
    #cfdict[zptstr]["totbkg"] = round(totyzpt/totorig*100,2)
    #cfdict[hptstr]["totbkg"] = round(totyhpt/totorig*100,2)
    #cfdict[btgstr]["totbkg"] = round(totybtg/totorig*100,2)
    #cfdict[sbstr]["totbkg"]  = round(totysb/totorig*100,2)

    cfdict[skimstr]["totbkg"] = round(totskim,2)
    cfdict[trigstr]["totbkg"] = round(tottrig,2)
    cfdict[zrecostr]["totbkg"] = round(totzrec,2)
    cfdict[hrecostr]["totbkg"] = round(tothrec,2)
    cfdict[metstr]["totbkg"] = round(totymet,2)
    cfdict[zptstr]["totbkg"] = round(totyzpt,2)
    cfdict[hptstr]["totbkg"] = round(totyhpt,2)
    cfdict[btgstr]["totbkg"] = round(totybtg,2)
    cfdict[sbstr]["totbkg"]  = round(totysb,2)

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
        
    

    

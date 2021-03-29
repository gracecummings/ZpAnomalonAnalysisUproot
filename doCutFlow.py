import subprocess
import gecorg as go
import argparse
import ROOT
import numpy as np

def gatherYields(topiaryf,bkgbinf):
    topiary       = ROOT.TFile(topiaryf)        
    scale         = bkgbinf["scale"]
    evntstot      = float(str(bkgbinf["tfile"].Get('hnevents').GetString()))*scale
    evntspassmet  = float(str(bkgbinf["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
    evntspasszpt  = float(str(bkgbinf["tfile"].Get('hnevents_pZ').GetString()))*scale
    evntspasshpt  = float(str(bkgbinf["tfile"].Get('hnevents_ph').GetString()))*scale
    evntspassbtag = float(str(bkgbinf["tfile"].Get('hnevents_btag').GetString()))*scale
    evntspasssb   = float(str(bkgbinf["tfile"].Get('hnevents_sb').GetString()))*scale
    evntsinskim   = topiary.Get('hnskimed').GetBinContent(1)*scale
    evntspasstrig = topiary.Get('htrigpass').GetBinContent(1)*scale
    evntspasszrec = topiary.Get('hZpass').GetBinContent(1)*scale
    evntspasshrec = topiary.Get('hHpass').GetBinContent(1)*scale

    return evntstot,evntsinskim,evntspasstrig,evntspasszrec,evntspasshrec,evntspassmet,evntspasszpt,evntspasshpt,evntspassbtag,evntspasssb

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
    bkgupout17 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-29','upout',zptcut,hptcut,metcut,btagwp,17)
    bkgtopia17 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-26','topiary',0.0,250.0,0.0,0.0,17)
    bkgupout18 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-29','upout',zptcut,hptcut,metcut,btagwp,18)
    bkgtopia18 = go.gatherBkg('analysis_output_ZpAnomalon/2021-03-28','topiary',0.0,250.0,0.0,0.0,18)
    
    #bkguncs  = np.load('analysis_output_ZpAnomalon/'+args.date+'/Fall17.AllZpAnomalonBkgs_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')
    
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    bkginfo17  = go.prepBkg(bkgupout17,bkgnames,bkgcols,'xsects_2017.ini',41.53,"yes")#gathers xs scales
    bkginfo18  = go.prepBkg(bkgupout18,bkgnames,bkgcols,'xsects_2017.ini',59.74,"yes")

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
    cf17dict = {metstr:{},zptstr:{},hptstr:{},btgstr:{},sbstr:{},skimstr:{},trigstr:{},zrecostr:{},hrecostr:{}}
    cf18dict = {metstr:{},zptstr:{},hptstr:{},btgstr:{},sbstr:{},skimstr:{},trigstr:{},zrecostr:{},hrecostr:{}}
    
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
    
    totorig17 = 0
    totskim17 = 0
    tottrig17 = 0
    totzrec17 = 0
    tothrec17 = 0
    totymet17 = 0
    totyzpt17 = 0
    totyhpt17 = 0
    totybtg17 = 0
    totysb17  = 0
    
    totorig18 = 0
    totskim18 = 0
    tottrig18 = 0
    totzrec18 = 0
    tothrec18 = 0
    totymet18 = 0
    totyzpt18 = 0
    totyhpt18 = 0
    totybtg18 = 0
    totysb18  = 0

    for bkg in bkginfo17:
        toplist17 = []
        toplist18 = []
        uplist18  = []
        for top in bkgtopia17:
            if bkg["name"] in top[0]:
                toplist17 = top
        for top in bkgtopia18:
            if bkg["name"] in top[0]:
                toplist18 = top
        for up in bkginfo18:
            if bkg["name"] == up["name"]:
                uplist18 = up

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

        origevnts17 = 0
        yieldskim17 = 0
        yieldtrig17 = 0
        yieldzrec17 = 0
        yieldhrec17 = 0
        yieldmet17  = 0
        yieldzpt17  = 0
        yieldhpt17  = 0
        yieldbtag17 = 0
        yieldsb17   = 0
        
        origevnts18 = 0
        yieldskim18 = 0
        yieldtrig18 = 0
        yieldzrec18 = 0
        yieldhrec18 = 0
        yieldmet18  = 0
        yieldzpt18  = 0
        yieldhpt18  = 0
        yieldbtag18 = 0
        yieldsb18   = 0

        for b,bkgbin in enumerate(bkg["binlist"]):
            reltop17 = 'str'
            reltop18 = 'str'
            relup18  = {}
            for bkgtopia in toplist17:
                if bkgbin["binname"] in bkgtopia:
                    reltop17 = bkgtopia
            for bkgtopia in toplist18:
                if bkgbin["binname"] in bkgtopia:
                    reltop18 = bkgtopia
            for bkgbin18 in uplist18["binlist"]:
                if bkgbin["binname"] in bkgbin18["binname"]:
                    relup18 = bkgbin18
            
            topiary       = ROOT.TFile(reltop17)        
            scale         = bkgbin["scale"]
            evntstot17      = float(str(bkgbin["tfile"].Get('hnevents').GetString()))*scale
            evntspassmet17  = float(str(bkgbin["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
            evntspasszpt17  = float(str(bkgbin["tfile"].Get('hnevents_pZ').GetString()))*scale
            evntspasshpt17  = float(str(bkgbin["tfile"].Get('hnevents_ph').GetString()))*scale
            evntspassbtag17 = float(str(bkgbin["tfile"].Get('hnevents_btag').GetString()))*scale
            evntspasssb17   = float(str(bkgbin["tfile"].Get('hnevents_sb').GetString()))*scale
            evntsinskim17   = topiary.Get('hnskimed').GetBinContent(1)*scale
            evntspasstrig17 = topiary.Get('htrigpass').GetBinContent(1)*scale
            evntspasszrec17 = topiary.Get('hZpass').GetBinContent(1)*scale
            evntspasshrec17 = topiary.Get('hHpass').GetBinContent(1)*scale

            etot18,eskim18,etrig18,ezrec18,ehrec18,emet18,ezpt18,ehpt18,ebtag18,esb18 = gatherYields(reltop18,relup18)

            origevnts += (evntstot17+etot18)
            yieldskim += (evntsinskim17+eskim18)
            yieldtrig += (evntspasstrig17+etrig18)
            yieldzrec += (evntspasszrec17+ezrec18)
            yieldhrec += (evntspasshrec17+ehrec18)
            yieldmet  += (evntspassmet17+emet18)
            yieldzpt  += (evntspasszpt17+ezpt18)
            yieldhpt  += (evntspasshpt17+ehpt18)
            yieldbtag += (evntspassbtag17+ebtag18)
            yieldsb   += (evntspasssb17+esb18)

            origevnts17 += evntstot17
            yieldskim17 += evntsinskim17
            yieldtrig17 += evntspasstrig17
            yieldzrec17 += evntspasszrec17
            yieldhrec17 += evntspasshrec17
            yieldmet17  += evntspassmet17
            yieldzpt17  += evntspasszpt17
            yieldhpt17  += evntspasshpt17
            yieldbtag17 += evntspassbtag17
            yieldsb17   += evntspasssb17

            origevnts18 += etot18
            yieldskim18 += eskim18
            yieldtrig18 += etrig18
            yieldzrec18 += ezrec18
            yieldhrec18 += ehrec18
            yieldmet18  += emet18
            yieldzpt18  += ezpt18
            yieldhpt18  += ehpt18
            yieldbtag18 += ebtag18
            yieldsb18   += esb18


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

        totorig17 += origevnts17
        totskim17 += yieldskim17
        tottrig17 += yieldtrig17
        totzrec17 += yieldzrec17
        tothrec17 += yieldhrec17
        totymet17 += yieldmet17
        totyzpt17 += yieldzpt17
        totyhpt17 += yieldhpt17
        totybtg17 += yieldbtag17
        totysb17  += yieldsb17

        totorig18 += origevnts18
        totskim18 += yieldskim18
        tottrig18 += yieldtrig18
        totzrec18 += yieldzrec18
        tothrec18 += yieldhrec18
        totymet18 += yieldmet18
        totyzpt18 += yieldzpt18
        totyhpt18 += yieldhpt18
        totybtg18 += yieldbtag18
        totysb18  += yieldsb18

        cfdict[skimstr][bkg["name"]] = round(yieldskim,2)
        cfdict[trigstr][bkg["name"]] = round(yieldtrig,2)
        cfdict[zrecostr][bkg["name"]] = round(yieldzrec,2)
        cfdict[hrecostr][bkg["name"]] = round(yieldhrec,2)
        cfdict[metstr][bkg["name"]] = round(yieldmet,2)
        cfdict[zptstr][bkg["name"]] = round(yieldzpt,2)
        cfdict[hptstr][bkg["name"]] = round(yieldhpt,2)
        cfdict[btgstr][bkg["name"]] = round(yieldbtag,2)
        cfdict[sbstr][bkg["name"]] = round(yieldsb,2)

        cf17dict[skimstr][bkg["name"]] = round(yieldskim17,2)
        cf17dict[trigstr][bkg["name"]] = round(yieldtrig17,2)
        cf17dict[zrecostr][bkg["name"]] = round(yieldzrec17,2)
        cf17dict[hrecostr][bkg["name"]] = round(yieldhrec17,2)
        cf17dict[metstr][bkg["name"]] = round(yieldmet17,2)
        cf17dict[zptstr][bkg["name"]] = round(yieldzpt17,2)
        cf17dict[hptstr][bkg["name"]] = round(yieldhpt17,2)
        cf17dict[btgstr][bkg["name"]] = round(yieldbtag17,2)
        cf17dict[sbstr][bkg["name"]] = round(yieldsb17,2)

        cf18dict[skimstr][bkg["name"]] = round(yieldskim18,2)
        cf18dict[trigstr][bkg["name"]] = round(yieldtrig18,2)
        cf18dict[zrecostr][bkg["name"]] = round(yieldzrec18,2)
        cf18dict[hrecostr][bkg["name"]] = round(yieldhrec18,2)
        cf18dict[metstr][bkg["name"]] = round(yieldmet18,2)
        cf18dict[zptstr][bkg["name"]] = round(yieldzpt18,2)
        cf18dict[hptstr][bkg["name"]] = round(yieldhpt18,2)
        cf18dict[btgstr][bkg["name"]] = round(yieldbtag18,2)
        cf18dict[sbstr][bkg["name"]] = round(yieldsb18,2)



    cfdict[skimstr]["totbkg"] = round(totskim,2)
    cfdict[trigstr]["totbkg"] = round(tottrig,2)
    cfdict[zrecostr]["totbkg"] = round(totzrec,2)
    cfdict[hrecostr]["totbkg"] = round(tothrec,2)
    cfdict[metstr]["totbkg"] = round(totymet,2)
    cfdict[zptstr]["totbkg"] = round(totyzpt,2)
    cfdict[hptstr]["totbkg"] = round(totyhpt,2)
    cfdict[btgstr]["totbkg"] = round(totybtg,2)
    cfdict[sbstr]["totbkg"]  = round(totysb,2)

    cf17dict[skimstr]["totbkg"] = round(totskim17,2)
    cf17dict[trigstr]["totbkg"] = round(tottrig17,2)
    cf17dict[zrecostr]["totbkg"] = round(totzrec17,2)
    cf17dict[hrecostr]["totbkg"] = round(tothrec17,2)
    cf17dict[metstr]["totbkg"] = round(totymet17,2)
    cf17dict[zptstr]["totbkg"] = round(totyzpt17,2)
    cf17dict[hptstr]["totbkg"] = round(totyhpt17,2)
    cf17dict[btgstr]["totbkg"] = round(totybtg17,2)
    cf17dict[sbstr]["totbkg"]  = round(totysb17,2)

    cf18dict[skimstr]["totbkg"] = round(totskim18,2)
    cf18dict[trigstr]["totbkg"] = round(tottrig18,2)
    cf18dict[zrecostr]["totbkg"] = round(totzrec18,2)
    cf18dict[hrecostr]["totbkg"] = round(tothrec18,2)
    cf18dict[metstr]["totbkg"] = round(totymet18,2)
    cf18dict[zptstr]["totbkg"] = round(totyzpt18,2)
    cf18dict[hptstr]["totbkg"] = round(totyhpt18,2)
    cf18dict[btgstr]["totbkg"] = round(totybtg18,2)
    cf18dict[sbstr]["totbkg"]  = round(totysb18,2)

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
    cftab.write("\caption{Cut Flow based on MC weighted yields with luminosity scaling for 2017 and 2018. Cut description for skims only relates to 2018, 2017 skims excluded the fat jet requirement.}\n")
    cftab.write("\label{tab:cutflow}\n")
    cftab.write("\end{center}\n")
    cftab.write("\end{table}\n")
    cftab.close()
        
    print("Cut       |      2017 yield      |     2018 yield  | Total ")
    print("total e   | {0} | {1} | {2} ".format(totorig17,totorig18,totorig))
    print("tot skim  | {0} | {1} | {2} ".format(totskim17,totskim18,totskim))
    print("tot trig  | {0} | {1} | {2} ".format(tottrig17,tottrig18,tottrig))
    print("tot zreco | {0} | {1} | {2} ".format(totzrec17,totzrec18,totzrec))
    print("tot hreco | {0} | {1} | {2} ".format(tothrec17,tothrec18,tothrec))
    print("tot met   | {0} | {1} | {2} ".format(totymet17,totymet18,totymet))
    print("tot zpt   | {0} | {1} | {2} ".format(totyzpt17,totyzpt18,totyzpt))
    print("tot hpt   | {0} | {1} | {2} ".format(totyhpt17,totyhpt18,totyhpt))
    print("tot btag  | {0} | {1} | {2} ".format(totybtg17,totybtg18,totybtg))
    print("tot sb    | {0} | {1} | {2} ".format(totysb17,totysb18,totysb))

    cutgrid = open("cutgrid.txt","a")
    cutgrid.write("{0} {1} {2} {3} {4} {5} {6}\n".format(str(zptcut),str(metcut),str(cfdict[sbstr]["ZZTo2L2Q"]),str(cfdict[sbstr]["WZTo2L2Q"]),str(cfdict[sbstr]["TT"]),str(cfdict[sbstr]["DYJetsToLL"]),str(cfdict[sbstr]["totbkg"])))

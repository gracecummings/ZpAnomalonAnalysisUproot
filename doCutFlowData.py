import subprocess
import gecorg as go
import argparse
import ROOT
import numpy as np
import glob

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
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-date","--date",help="date folder with plots to stack")
    args = parser.parse_args()

    #Get command line parameters
    sig_xsec      = args.xsec
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut
    btagwp        = args.btagwp

    dat17top = glob.glob("analysis_output_ZpAnomalon/2021-03-26/Run2017*_topiary_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root")
    dat18top = glob.glob("analysis_output_ZpAnomalon/2021-03-28/Run2018*_topiary_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root")

    dat17up = glob.glob("analysis_output_ZpAnomalon/2021-04-06/Run2017*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut"+str(zptcut)+"_Hptcut"+str(hptcut)+"_metcut"+str(metcut)+"_btagwp"+str(btagwp)+".root")
    dat18up = glob.glob("analysis_output_ZpAnomalon/2021-04-06/Run2018*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut"+str(zptcut)+"_Hptcut"+str(hptcut)+"_metcut"+str(metcut)+"_btagwp"+str(btagwp)+".root")

    tops = [dat17top,dat18top]
    outs = [dat17up,dat18up]
    
    skimstr  = "Events with \( > 0\) fat jets and \( Z(\ell^{+}\ell^{-}) cand \)"
    trigstr  = "Events Pass trigger"
    zrecostr = "Events with \( 70 < m_{\ell\ell} < 110 \)"
    hrecostr = "Events Pass $H$ reco "
    sbstr  = "Events in side band "
    metstr = "Events Passing \(MET > "+str(metcut)+"\)"
    zptstr = "Events Passing \(Z p_{T} > "+str(zptcut)+"\)"
    hptstr = "Events Passing Fat Jet \(p_{T} > "+str(hptcut)+"\)"
    btgstr = "Events Passing {0} btag WP".format(str(btagwp))

    
    
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

    tot17orig = 0
    tot17skim = 0
    tot17trig = 0
    tot17zrec = 0
    tot17hrec = 0
    tot17ymet = 0
    tot17yzpt = 0
    tot17yhpt = 0
    tot17ybtg = 0
    tot17ysb  = 0

    tot18orig = 0
    tot18skim = 0
    tot18trig = 0
    tot18zrec = 0
    tot18hrec = 0
    tot18ymet = 0
    tot18yzpt = 0
    tot18yhpt = 0
    tot18ybtg = 0
    tot18ysb  = 0

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
    orig17evnts = 0
    yield17skim = 0
    yield17trig = 0
    yield17zrec = 0
    yield17hrec = 0
    yield17sb   = 0
    yield17met  = 0
    yield17zpt  = 0
    yield17hpt  = 0
    yield17btag = 0
    orig18evnts = 0
    yield18skim = 0
    yield18trig = 0
    yield18zrec = 0
    yield18hrec = 0
    yield18sb   = 0
    yield18met  = 0
    yield18zpt  = 0
    yield18hpt  = 0
    yield18btag = 0

    
    for d,dat in enumerate(outs):
        for e,era in enumerate(dat):
            topiary       = ROOT.TFile(tops[d][e])        
            evntstot      = float(str(ROOT.TFile(era).Get('hnevents').GetString()))
            evntspassmet  = float(str(ROOT.TFile(era).Get('hnevents_pMET').GetString()))
            evntspasszpt  = float(str(ROOT.TFile(era).Get('hnevents_pZ').GetString()))
            evntspasshpt  = float(str(ROOT.TFile(era).Get('hnevents_ph').GetString()))
            evntspassbtag = float(str(ROOT.TFile(era).Get('hnevents_btag').GetString()))
            evntspasssb   = float(str(ROOT.TFile(era).Get('hnevents_sb').GetString()))
            evntsinskim   = topiary.Get('hnskimed').GetBinContent(1)
            evntspasstrig = topiary.Get('htrigpass').GetBinContent(1)
            evntspasszrec = topiary.Get('hZpass').GetBinContent(1)
            evntspasshrec = topiary.Get('hHpass').GetBinContent(1)

            origevnts += evntstot
            yieldskim += evntsinskim
            yieldtrig += evntspasstrig
            yieldzrec += evntspasszrec
            yieldhrec += evntspasshrec
            yieldsb   += evntspasssb
            yieldmet  += evntspassmet
            yieldzpt  += evntspasszpt
            yieldhpt  += evntspasshpt
            yieldbtag += evntspassbtag
            
            if "Run2017" in era:
                orig17evnts += evntstot
                yield17skim += evntsinskim
                yield17trig += evntspasstrig
                yield17zrec += evntspasszrec
                yield17hrec += evntspasshrec
                yield17sb   += evntspasssb
                yield17met  += evntspassmet
                yield17zpt  += evntspasszpt
                yield17hpt  += evntspasshpt
                yield17btag += evntspassbtag

            if "Run2018" in era:
                orig18evnts += evntstot
                yield18skim += evntsinskim
                yield18trig += evntspasstrig
                yield18zrec += evntspasszrec
                yield18hrec += evntspasshrec
                yield18sb   += evntspasssb
                yield18met  += evntspassmet
                yield18zpt  += evntspasszpt
                yield18hpt  += evntspasshpt
                yield18btag += evntspassbtag


    cfdict[skimstr]["Run 2"] = round(yieldskim,2)
    cfdict[trigstr]["Run 2"] = round(yieldtrig,2)
    cfdict[zrecostr]["Run 2"] = round(yieldzrec,2)
    cfdict[hrecostr]["Run 2"] = round(yieldhrec,2)
    cfdict[metstr]["Run 2"] = round(yieldmet,2)
    cfdict[zptstr]["Run 2"] = round(yieldzpt,2)
    cfdict[hptstr]["Run 2"] = round(yieldhpt,2)
    cfdict[btgstr]["Run 2"] = round(yieldbtag,2)
    cfdict[sbstr]["Run 2"] = round(yieldsb,2)

    cfdict[skimstr]["2017"] = round(yield17skim,2)
    cfdict[trigstr]["2017"] = round(yield17trig,2)
    cfdict[zrecostr]["2017"] = round(yield17zrec,2)
    cfdict[hrecostr]["2017"] = round(yield17hrec,2)
    cfdict[metstr]["2017"] = round(yield17met,2)
    cfdict[zptstr]["2017"] = round(yield17zpt,2)
    cfdict[hptstr]["2017"] = round(yield17hpt,2)
    cfdict[btgstr]["2017"] = round(yield17btag,2)
    cfdict[sbstr]["2017"] = round(yield17sb,2)

    cfdict[skimstr]["2018"] = round(yield18skim,2)
    cfdict[trigstr]["2018"] = round(yield18trig,2)
    cfdict[zrecostr]["2018"] = round(yield18zrec,2)
    cfdict[hrecostr]["2018"] = round(yield18hrec,2)
    cfdict[metstr]["2018"] = round(yield18met,2)
    cfdict[zptstr]["2018"] = round(yield18zpt,2)
    cfdict[hptstr]["2018"] = round(yield18hpt,2)
    cfdict[btgstr]["2018"] = round(yield18btag,2)
    cfdict[sbstr]["2018"] = round(yield18sb,2)
        
    cutFlowTableTex = go.makeOutFile("ZpAnomalon_mumu_data","cutflow",'.tex',str(int(zptcut)),str(int(hptcut)),str(int(metcut)),str(btagwp).split('.')[-1]+'E-10')
    cftab = open(cutFlowTableTex,"w")
    cftab.write(r'\begin{table}[htbp]')
    cftab.write('\n')
    cftab.write(r'\begin{center}')
    cftab.write('\n')
    cftab.write(r'\begin{tabular}{l | c | c | c }')#need to make dynamic
    cftab.write('\n')
    cftab.write("\hline\hline\n")
    cftab.write(r'cut description & 2017 & 2018 & Combined\\')
    cftab.write('\n')
    cftab.write("\hline\n")
    cftab.write(skimstr+' & '+str(cfdict[skimstr]['2017'])+' & '+str(cfdict[skimstr]['2018'])+' & '+str(cfdict[skimstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(trigstr+' & '+str(cfdict[trigstr]['2017'])+' & '+str(cfdict[trigstr]['2018'])+' & '+str(cfdict[trigstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(zrecostr+' & '+str(cfdict[zrecostr]['2017'])+' & '+str(cfdict[zrecostr]['2018'])+' & '+str(cfdict[zrecostr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(hrecostr+' & '+str(cfdict[hrecostr]['2017'])+' & '+str(cfdict[hrecostr]['2018'])+' & '+str(cfdict[hrecostr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write("\hline\n")
    cftab.write(sbstr+' & '+str(cfdict[sbstr]['2017'])+' & '+str(cfdict[sbstr]['2018'])+' & '+str(cfdict[sbstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(metstr+' & '+str(cfdict[metstr]['2017'])+' & '+str(cfdict[metstr]['2018'])+' & '+str(cfdict[metstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(zptstr+' & '+str(cfdict[zptstr]['2017'])+' & '+str(cfdict[zptstr]['2018'])+' & '+str(cfdict[zptstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(hptstr+' & '+str(cfdict[hptstr]['2017'])+' & '+str(cfdict[hptstr]['2018'])+' & '+str(cfdict[hptstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(btgstr+' & '+str(cfdict[btgstr]['2017'])+' & '+str(cfdict[btgstr]['2018'])+' & '+str(cfdict[btgstr]['Run 2']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write("\end{tabular}\n")
    cftab.write("\caption{Data cutflow in sideband region.}\n")
    cftab.write("\label{tab:dateff}\n")
    cftab.write("\end{center}\n")
    cftab.write("\end{table}\n")
    cftab.close()

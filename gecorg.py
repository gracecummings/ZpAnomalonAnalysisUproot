import os
import sys
import glob
import ROOT
import configparser
import pandas as pd
from datetime import date

def sampleType(sampstring):
    #Make numerical code for type of sample
    if "Run" in sampstring:
        samptype = 0
    elif "ZpAnomalon" in sampstring:
        samptype = 1
        year = 17
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

    if "2018" in sampstring:
        year = 18
    if "Autumn18" in sampstring:
        year = 18
    if "2017" in sampstring:
        year = 17
    if "Fall17" in sampstring:
        year = 17
    if "2016" in sampstring:
        year = 16

    return samptype,year

def makeOutFile(sampstring,descrip,ftype,zptcut,hptcut,metcut,btagwp):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+descrip+"_Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+ftype
    return outFile

def orderDY(histFile):
    s1 = histFile.split("to")[0]
    s2 = s1.split("HT-")[1]
    return int(s2)       

def orderFall17TT(histFile):#NROKEN AT THE MOMENT
    s1 = histFile.split("Events")[0]
    s2 = s1.split("_")[-1]
    return int(s2)

def massPoints(nameSig):#do this for full name
    s1  = nameSig.split("-")
    #print s1
    mzp = int(s1[0].split("Zp")[1])
    mnd = int(s1[1].split("ND")[1])
    ms1 = s1[2].split("_upout")[0]
    mns = int(ms1.split("NS")[1])
    return mzp,mnd,mns

def nameSignal(histFile):
    s1 = histFile.split("ZpAnomalonHZ_UFO-")[1]
    s2 = s1.split("_upout")[0]
    return s2

def findScale(prodnum,lumi,xsec):
    expecnum = xsec*lumi
    scalefac = expecnum/prodnum
    return  scalefac

def colsFromPalette(samplist,palname):
    collist = []
    ROOT.gStyle.SetPalette(palname)
    cols = ROOT.TColor.GetPalette()
    colsnum = cols.GetSize()
    for i in range(len(samplist)):
        collist.append(cols.At(0+i*int(colsnum/len(samplist))))
    collist.reverse()
    return collist

def gatherBkg(bkg_dir,descrip,zptcut,hptcut,metcut,btagwp,year):
    if year == 18:
        mcprefix = 'Autumn18'
    if year == 17:
        mcprefix = 'Fall17'
        
    DYJetsToLL = glob.glob(str(bkg_dir)+'/'+mcprefix+'.DYJetsToLL_M-50_HT*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    TT         = glob.glob(str(bkg_dir)+'/'+mcprefix+'.TTT*_'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                         
    WZTo2L2Q   = glob.glob(str(bkg_dir)+'/'+mcprefix+'.WZTo2L2Q*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                    
    ZZTo2L2Q   = glob.glob(str(bkg_dir)+'/'+mcprefix+'.ZZTo2L2Q*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                    
    bkgfiles   = [DYJetsToLL,TT,WZTo2L2Q,ZZTo2L2Q]
    return bkgfiles

def prepSig(sigfiles,sig_colors,sig_xsec,lumi):
    sig_info = []
    for s,sig in enumerate(sigfiles):
        sig_dict = {}                                                                                   
        sig_dict["tfile"] = ROOT.TFile(sig)                                                             
        sig_samplesize    = str(sig_dict["tfile"].Get('hnevents').GetString())
        sig_dict["scale"] = findScale(float(sig_samplesize),sig_xsec,lumi)
        sig_dict["name"]  = nameSignal(sig)
        mzp,mnd,mns       = massPoints(sig_dict["name"])
        sig_dict["mzp"]   = mzp#
        sig_dict["mnd"]   = mnd
        sig_dict["mns"]   = mns
        sig_info.append(sig_dict)                                                                    
    #Sort Signals by ND mass, then by Zp mass                                                           
    sig_info = sorted(sig_info,key = lambda sig: (sig["mnd"],sig["mzp"],sig["mns"]))                    
    for s,sig in enumerate(sig_info):                                                                   
        sig["color"] = sig_colors[s]

    return sig_info

def prepBkg(bkgfiles,bkgnames,bkg_colors,ini_file,lumi,flag="yes"):
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open(ini_file)
    config.read_file(fp)
    bkg_info = []
    for b,bkg in enumerate(bkgfiles):
        bkg_binsum   = {}
        bkg_binlist  = []
        bkg_channel  = bkgnames[b]
        bkg_expyield = 0
        # bkg xs from .ini file
        bkgbin_xs_pairs = config.items(bkg_channel)
        if bkg_channel == "DYJetsToLL":
            #orders smallest HT to largest
            bkg.sort(key = orderDY)
        elif bkg_channel == "TT":
            #sorts in alphabetical order 
            bkg.sort()                                                     
        else:
            if flag == "no":
                break
        #loop through each process bin or categrory
        for s,bkgbin in enumerate(bkg):
            bkgbin_dict = {}
            bkgbin_dict["binname"] = bkgbin_xs_pairs[s][0]
            bkgbin_dict["tfile"]   = ROOT.TFile(bkgbin)
            bkgbin_sampsize        = str(bkgbin_dict["tfile"].Get('hnevents').GetString())
            bkgbin_xs              = float(bkgbin_xs_pairs[s][1].split()[0])*1000#Into Femtobarn
            bkgbin_dict["scale"]   = findScale(float(bkgbin_sampsize),bkgbin_xs,lumi)
            bkgbin_dict["color"]   = bkg_colors[b]
            #get the number of passing events
            bkgbin_yield           = float(str(bkgbin_dict["tfile"].Get('hnevents_pZ').GetString()))#last cut hist ##goes back to btag
            bkg_expyield          += bkgbin_yield*bkgbin_dict["scale"]
            bkg_binlist.append(bkgbin_dict)
            bkg_binsum["expyield"] = bkg_expyield

        bkg_binsum["binlist"] = bkg_binlist
        bkg_binsum["name"]    = bkg_channel
        bkg_info.append(bkg_binsum)

    #Sort the backgrounds from the smallest yields to largest        
    bkg_info = sorted(bkg_info, key = lambda bkg:bkg["expyield"])
    return bkg_info
#
def stackBkg(bkg_info,hist_to_stack,hsbkg,legend,stack_max,stack_min):
    for bkg in bkg_info:
        for b,bkgbin in enumerate(bkg["binlist"]):
            hbkg = bkgbin["tfile"].Get(hist_to_stack)
            hbkg.SetStats(0)
            hbkg.Scale(bkgbin["scale"])
            hbkg.SetFillColor(bkgbin["color"])
            hbkg.SetLineColor(bkgbin["color"])
            hbkg.SetMaximum(stack_max)
            hbkg.SetMinimum(stack_min)
            hsbkg.Add(hbkg)
            hsbkg.Draw("HIST")
            hsbkg.SetMaximum(stack_max)
            hsbkg.SetMinimum(stack_min)
            if b == len(bkg["binlist"])-1:
                legend.AddEntry(hbkg,bkg["name"],"f")
                
def stackBkgMultiYear(bkginfo0,bkginfo1,hist_to_stack,hsbkg,legend,stack_max,stack_min):
    for i,bkg in enumerate(bkginfo0):
        #print(bkg["name"])
        #print(bkginfo1[i]["name"])
        for b,bkgbin in enumerate(bkg["binlist"]):
            #print(bkgbin["binname"])
            #print(bkginfo1[i]["binlist"][b]["binname"])
            hbkg0 = bkgbin["tfile"].Get(hist_to_stack)
            hbkg0.SetStats(0)
            hbkg0.Scale(bkgbin["scale"])
            hbkg0.SetFillColor(bkgbin["color"])
            hbkg0.SetLineColor(bkgbin["color"])
            hbkg0.SetMaximum(stack_max)
            hbkg0.SetMinimum(stack_min)
            hbkg1 = bkginfo1[i]["binlist"][b]["tfile"].Get(hist_to_stack)
            hbkg1.SetStats(0)
            hbkg1.Scale(bkginfo1[i]["binlist"][b]["scale"])
            hbkg1.SetFillColor(bkginfo1[i]["binlist"][b]["color"])
            hbkg1.SetLineColor(bkginfo1[i]["binlist"][b]["color"])
            hbkg1.SetMaximum(stack_max)
            hbkg1.SetMinimum(stack_min)

            hsbkg.Add(hbkg0)
            hsbkg.Add(hbkg1)
            hsbkg.Draw("HIST")
            hsbkg.SetMaximum(stack_max)
            hsbkg.SetMinimum(stack_min)
            
            if b == len(bkg["binlist"])-1:
                legend.AddEntry(hbkg1,bkg["name"],"f")

def saveNpUncertainties(uncdf,filename):
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

class backgrounds:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp

        #gather background MC files
        self.f17dyjetsb = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsr = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsb = glob.glob(str(path)+'/Fall17.TTT*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsr = glob.glob(str(path)+'/Fall17.TTT*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsb = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsr = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsb = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsr = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsb = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsr = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsb = glob.glob(str(path)+'/Autumn18.TTT*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsr = glob.glob(str(path)+'/Autumn18.TTT*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsb = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsr = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsb = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsr = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')


        #gather errors
        self.f17dyjetsberrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsrerrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsberrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsrerrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsberrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsrerrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsberrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsrerrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        self.a18dyjetsberrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsrerrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsberrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsrerrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsberrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsrerrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsberrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_sideband_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsrerrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_signalr_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        self.bkgs = {"DYJetsToLL":
                {18:
                 {"sb":[self.a18dyjetsb,self.a18dyjetsberrs],
                  "sr":[self.a18dyjetsr,self.a18dyjetsrerrs]},
                 17:
                 {"sb":[self.f17dyjetsb,self.f17dyjetsberrs],
                  "sr":[self.f17dyjetsr,self.f17dyjetsrerrs]},
                },
                "TT":
                {18:
                 {"sb":[self.a18ttsb,self.a18ttsberrs],
                  "sr":[self.a18ttsr,self.a18ttsrerrs]},
                 17:
                 {"sb":[self.f17ttsb,self.f17ttsberrs],
                  "sr":[self.f17ttsr,self.f17ttsrerrs]},
                },
                "WZTo2L2Q":
                {18:
                 {"sb":[self.a18wzsb,self.a18wzsberrs],
                  "sr":[self.a18wzsr,self.a18wzsrerrs]},
                 17:
                 {"sb":[self.f17wzsb,self.f17wzsberrs],
                  "sr":[self.f17wzsr,self.f17wzsrerrs]},
                },
                "ZZTo2L2Q":
                {18:
                 {"sb":[self.a18zzsb,self.a18zzsberrs],
                  "sr":[self.a18zzsr,self.a18zzsrerrs]},
                 17:
                 {"sb":[self.f17zzsb,self.f17zzsberrs],
                  "sr":[self.f17zzsr,self.f17zzsrerrs]},
                }
        }

        self.config = configparser.RawConfigParser()
        self.config.optionxform = str
        fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
        self.config.read_file(fp)

        
    def getAddedHist(self,hist,samp,region,hname,years = [17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        bkgdfs  = []
        holder = []
        
        for year in years:
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TTTo" in samp:
                files.sort()
                errs.sort()
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                tf = ROOT.TFile(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
                scale = findScale(numevents,xs,lumi)
                h = tf.Get(hname)
                h.Scale(scale)
                hist.Add(h)
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)

                debugstring = name+" "+str(xs)+" "+str(scale)
                holder.append(debugstring)

        uncsqdDYJetsdf = sum(bkgdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = uncDYJetsdf['h_zp_jigm'][ibin-1]
                hist.SetBinError(ibin,binerr)

        return hist,holder

    #Add something that does this nicely within this class
    #no craziness
    #def getStackofBkgs(self,hstk,region,hname,legend,years = [17,18]):
    #    bkgfiles17 = [self.bkgs["DYJetsToLL"][17][region][0],
    #                  self.bkgs["TT"][17][region][0],
    #                  self.bkgs["WZTo2L2Q"][17][region][0],
    #                  self.bkgs["ZZTo2L2Q"][17][region][0]
    #                  ]
    #    bkgfiles18 = [self.bkgs["DYJetsToLL"][18][region][0],
    #                  self.bkgs["TT"][18][region][0],
    #                  self.bkgs["WZTo2L2Q"][18][region][0],
    #                  self.bkgs["ZZTo2L2Q"][18][region][0]
    #                  ]
    #    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    #    bkgcols  = colsFromPalette(bkgnames,ROOT.kLake)

    #    info17 = prepBkg(bkgfiles17,bkgnames,bkgcols,"xsects_2017.ini",41.53)
    #    info18 = prepBkg(bkgfiles18,bkgnames,bkgcols,"xsects_2017.ini",59.74)
    #    stackBkgMultiYear(info17,info18,hname,hstk,legend,18,0)


class run2:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp

        #gather data files
        self.run17sb = glob.glob(str(path)+'/Run2017*upout_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18sb = glob.glob(str(path)+'/Run2018*upout_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17sr = glob.glob(str(path)+'/Run2017*upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18sr = glob.glob(str(path)+'/Run2018*upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

        #gather errors
        self.run17sberrs = glob.glob(str(path)+'/Run2017*selected_errors_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18sberrs = glob.glob(str(path)+'/Run2018*selected_errors_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17srerrs = glob.glob(str(path)+'/Run2017*selected_errors_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18srerrs = glob.glob(str(path)+'/Run2018*selected_errors_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')

        self.data = {18:
                      {"sb":[self.run18sb,self.run18sberrs],
                       "sr":[self.run18sr,self.run18srerrs],
                       },
                      17:
                      {"sb":[self.run17sb,self.run17sberrs],
                       "sr":[self.run17sr,self.run17srerrs],
                       }
                      }

    def getAddedHist(self,hist,region,hname,years = [17,18]):
        data = self.data
        datadfs = []

        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile(f)
                h = tf.Get(hname)
                hist.Add(h)

                #calc errs
                df = pd.read_pickle(errs[i])
                sqrddf = df**2
                datadfs.append(sqrddf)

        uncsqddf = sum(datadfs)
        uncdf   = uncsqddf**(1/2)
        
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = uncdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)

        return hist


class signal:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,xs,lumi):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        
        #gather signal plots
        self.sigsr = glob.glob(str(path)+'/Zp*_upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.sigsb = glob.glob(str(path)+'/Zp*_upout_sideband*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.sigfl = glob.glob(str(path)+'/Zp*_upout_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        
        sig_colors = colsFromPalette(self.sigsr,ROOT.kCMYK)
        
        #prep signals
        self.prepsigsr = prepSig(self.sigsr,sig_colors,xs,lumi)

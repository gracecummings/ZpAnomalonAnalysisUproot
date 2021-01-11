import argparse
import ROOT
import glob
import os
import gecorg
import numpy as np
from datetime import date
from ROOT import kOrange, kViolet, kCyan, kGreen, kPink, kAzure, kMagenta, kBlue, kBird
from math import sqrt

if __name__=='__main__':
    #build module objects
    parser = argparse.ArgumentParser()

    #Define parser imputs
    parser.add_argument("-L","--lumi", type=float,default = 137, help = "integrated luminosity for scale in fb^-1")
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    args = parser.parse_args()

    #Get command line parameters
    lumi          = args.lumi
    sig_xsec      = args.xsec
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut

    #samples
    bkgfiles = gecorg.gatherBkg('analysis_output_ZpAnomalon/2021-01-05/','upout',zptcut,hptcut,metcut)
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    sigfiles = glob.glob('analysis_output_ZpAnomalon/2021-01-05/Zp*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.root')
    datfiles = glob.glob('analysis_output_ZpAnomalon/2021-01-05/Run2017*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.root')
    

    #Prep signals
    sig_colors = gecorg.colsFromPalette(sigfiles,ROOT.kCMYK)
    sig_info   = gecorg.prepSig(sigfiles,sig_colors,sig_xsec,lumi)

    #Prep backgrounds
    bkg_colors = gecorg.colsFromPalette(bkgnames,ROOT.kLake)
    bkg_info   = gecorg.prepBkg(bkgfiles,bkgnames,bkg_colors,'xsects_2017.ini',lumi)

    #Prep Data
    dat_info = [ROOT.TFile(dat) for dat in datfiles]

    #some beauty stuff
    max_plot = 1000.
    min_plot = 0.
    titles = {
        "h_z_pt":"Z p_{T}",
        "h_z_eta":"Z \eta",
        "h_z_m":"m_{Z}",
        "h_h_pt":"Higgs p_{T}",
        "h_h_eta":"Higgs \eta",
        "h_h_m":"m_{h}",
        "h_h_sd":"Higgs Soft Drop Mass",
        "h_met":"p_{T}^{miss}",
        "h_zp_jigm":"Jigsaw Mass Estimator Z'",
        "h_nd_jigm":"Jigsaw Mass Estimator ND",
        "h_ns_jigm":"Jigsaw Mass Estimator NS",
        "h_weights":"event weights",
    }


    #Make the stacks
    keys = dat_info[0].GetListOfKeys()
    for i,key in enumerate(keys):
        hname = key.GetName()
        if 'hnevents' in hname:
            break
        
        leg = ROOT.TLegend(0.45,0.55,0.90,0.88)

        #bkg stack
        hsbkg = ROOT.THStack('hsbkg','')
        gecorg.stackBkg(bkg_info,hname,hsbkg,leg,max_plot,min_plot)

        #data hist
        hsdat = dat_info[0].Get(hname)
        hsdat.SetStats(0)
        hsdat.SetMaximum(max_plot)
        hsdat.SetMinimum(min_plot)
        hsdat.SetMarkerStyle(8)
        leg.AddEntry(hsdat,"Data")
        for d,datf in enumerate(dat_info[1:]):
            hdat = datf.Get(hname)
            hdat.SetStats(0)
            hdat.SetMaximum(max_plot)
            hdat.SetMinimum(min_plot)
            hsdat.Add(hdat)
        hsdat.SetBinErrorOption(1)

        
        #Make a multigraph
        mg = ROOT.TMultiGraph()
                
        #Prep the pads
        tc = ROOT.TCanvas("tc",hname,600,800)
        p1 = ROOT.TPad("p1","stack_"+hname,0,0.4,1.0,1.0)
        #p1.SetLogy()
        #p1.SetBottomMargin(0)
        p1.SetLeftMargin(0.15)
        p1.SetRightMargin(0.05)
        p2 = ROOT.TPad("p2","signif_"+hname,0,0.0,1.0,0.4)
        #p2.SetTopMargin(0)
        p2.SetRightMargin(.05)
        p2.SetLeftMargin(0.15)
        p2.SetBottomMargin(0.2)
        
        #Prepare first pad for stack
        p1.Draw()
        p1.cd()

        #Draw the stack
        hsbkg.Draw("HIST")#add PFC for palette drawing
        hsbkg.GetXaxis().SetTitle(titles[hname])
        hsbkg.GetXaxis().SetTitleSize(0.05)
        hsbkg.GetYaxis().SetTitle("Events")
        hsbkg.GetYaxis().SetTitleSize(0.05)
        hsdat.Draw("HISTSAMEPE")
        tc.Modified()
        
        #Add the signal plots
        max_max = 0
        for p,masspoint in enumerate(sig_info):
            hsig = masspoint["tfile"].Get(hname)
            hsig.SetStats(0)
            hsig.Scale(masspoint["scale"])
            hsig.SetLineColor(masspoint["color"])
            hsig.SetLineWidth(2)
            hsig.SetMaximum(max_plot)
            hsig.SetMinimum(min_plot)
            hsig.Draw("HISTSAME")
            leg.AddEntry(hsig,"Zp"+str(masspoint["mzp"])+" ND"+str(masspoint["mnd"])+" NS"+str(masspoint["mns"])+" "+str(sig_xsec/1000)+" pb","l")

        #Now the ratio calculation
        hsumb = hsbkg.GetStack().Last()
        binlist = np.zeros(hsumb.GetNbinsX())
        fill    = np.zeros(hsumb.GetNbinsX())
        ratiolist = np.zeros(hsumb.GetNbinsX())
        pulllist = np.zeros(hsumb.GetNbinsX())
        rerrlist = np.zeros(hsumb.GetNbinsX())

        #Load in the uncertainties
        #uncf = np.load()
        
        for ibin in range(hsumb.GetNbinsX()):
            bincen = hsumb.GetBinCenter(ibin)
            bkgmc  = hsumb.GetBinContent(ibin)
            data   = hsdat.GetBinContent(ibin)
            datunc = hsdat.GetBinError(ibin)
            binlist[ibin] = bincen
            if bkgmc != 0:
                ratiolist[ibin] = data/bkgmc
                #pulllist[ibin]  = (data-bkgmc)/datunc
                rerrlist[ibin] = datunc/bkgmc
            if bkgmc == 0:
                ratiolist[ibin] = -1
                rerrlist[ibin] = 0
        
        #remove underflow bin
        ratiolist = np.delete(ratiolist,0)
        binlist   = np.delete(binlist,0)
        rerrlist  = np.delete(rerrlist,0)

        #Build the graphs
        tg = ROOT.TGraphErrors((hsumb.GetNbinsX()-1),binlist,ratiolist,fill,rerrlist)
        tg.SetTitle("")
        tg.SetMarkerStyle(8)
        mg.Add(tg)

        #Make the second pad with the ratio plot
        tc.cd()
        p2.Draw()
        p2.cd()
        l = ROOT.TLine(binlist[0],1,binlist[-1]+hsumb.GetBinWidth(1),1)
        #Now, the beauty aspects
        temp_max = np.amax(ratiolist)
        if temp_max > max_max:
            max_max = temp_max

        ratio_max = 2.0
        if max_max < 2.0:
            ratio_max = 2.0
        if max_max > 5:
            ratio_max = 5.0
        else:
            ratio_max = max_max
        
        mg.SetTitle("")
        #x axis
        mg.GetXaxis().SetTitle("bin center")
        mg.GetXaxis().SetTitleSize(0.07)
        mg.GetXaxis().SetLabelSize(0.05)
        mg.GetXaxis().SetLimits(binlist[0],binlist[-1]+hsumb.GetBinWidth(1))
        #y axis
        mg.GetYaxis().SetTitle("data/MC")
        mg.GetYaxis().SetTitleSize(0.07)
        mg.GetYaxis().SetTitleOffset(.7)
        mg.GetYaxis().SetLabelSize(0.05)
        mg.SetMinimum(0.0)
        #mg.SetMaximum(ratio_max)
        mg.SetMaximum(1.5)
        l.Draw()
        mg.Draw("AP")
        
        #Go back to previous pad so next kinematic plots draw
        tc.cd()
        p1.cd()
        
        #Draw the legent with everything added
        if not "50rig" in hname:
            leg.SetBorderSize(0)
            leg.Draw()

        #Draw the legent with everything added
        leg.SetBorderSize(0)
        leg.Draw()

        #Save the plot
        pngname = gecorg.makeOutFile(hname,'ratio','.png',str(zptcut),str(hptcut),str(metcut))
        tc.SaveAs(pngname)

import argparse
import ROOT
import glob
#import os
import gecorg
import numpy as np
#from datetime import date
#from ROOT import kOrange, kViolet, kCyan, kGreen, kPink, kAzure, kMagenta, kBlue, kBird
#from math import sqrt

if __name__=='__main__':
    #build module objects
    parser = argparse.ArgumentParser()

    #Define parser imputs
    parser.add_argument("-L","--lumi", type=float,default = 137, help = "integrated luminosity for scale in fb^-1")
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-p","--plot", type=str,help = "desired plot name to make the optimization for")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-date","--date",help="date folder with plots to stack")
    args = parser.parse_args()

    #Get command line parameters
    lumi          = args.lumi
    sig_xsec      = args.xsec
    released_plot = 'h_btag'
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut
    btaggers      = ['DeepMassDecorrelTagHbbvsQCD',#move to command line option?
                     'DeepMassDecorrelTagZbbvsQCD',
                     'DeepMassDecorrelTagZHbbvsQCD',
                     'pfMassIndependentDeepDoubleBvLJetTagsProbHbb',
                     ]
    sigsamp = 'Zp1200-ND175-NS1'
    plotmax       = 100000000000.0
    
    #Prep the plot
    mg = ROOT.TMultiGraph()
    gleg = ROOT.TLegend(0.45,0.1,0.90,0.3)
    tc = ROOT.TCanvas("tc","btagrco",600,600)
    

    btagcols = gecorg.colsFromPalette(btaggers,ROOT.kColorPrintableOnGrey)
    
    for i,btag in enumerate(btaggers):
        #Samples
        bkgfiles = gecorg.gatherBkg('analysis_output_ZpAnomalon/'+args.date+'/','upout_'+btag,zptcut,hptcut,metcut)
        bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
        sigfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/ZpAnomalonHZ_UFO-'+sigsamp+'_upout_'+btag+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'.root')

        #Prep signals
        sig_colors = gecorg.colsFromPalette(sigfiles,ROOT.kCMYK)
        sig_info   = gecorg.prepSig(sigfiles,sig_colors,sig_xsec,lumi)

        #Prep backgrounds
        bkg_colors = gecorg.colsFromPalette(bkgnames,ROOT.kLake)
        bkg_info   = gecorg.prepBkg(bkgfiles,bkgnames,bkg_colors,'xsects_2017.ini',lumi)

        #Make the stacked plot
        
        hname = released_plot
        hsbkg = ROOT.THStack('hsbkg','')
        leg = ROOT.TLegend(0.45,0.55,0.90,0.88)
        gecorg.stackBkg(bkg_info,released_plot,hsbkg,leg,plotmax,0.0)

        #Signal Plot
        hsig = sig_info[0]["tfile"].Get(released_plot)
        hsig.Scale(sig_info[0]["scale"])

        #LUT with titles
        titles = {
            "h_z_pt":"Z pT",
            "h_h_pt":"Higgs pT",
            "h_h_sd":"Higgs Soft Drop Mass",
            "h_met":"pT miss",
            "h_zp_jigm":"Jigsaw Mass Estimator Z'",
            "h_nd_jigm":"Jigsaw Mass Estimator ND",
            "h_btag":btag,
        }
        
        #Now the efficiency calculation
        hsum    = hsbkg.GetStack().Last()
        totbkg  = hsum.Integral()
        totsig  = hsig.Integral()
        cutlist = np.zeros(hsum.GetNbinsX())
        bkgeffs = np.zeros(hsum.GetNbinsX())
        sigeffs = np.zeros(hsum.GetNbinsX())

        if 'ZHbbvsQCD' in btag:
            print btag
            #print(totbkg)
            #print(totsig)
            for ibin in range(hsum.GetNbinsX()):
                bkg_pass = 0
                sig_pass = 0
                for b in range(ibin,hsum.GetNbinsX()+1):
                    #print(hsum.GetBinContent(b))
                    bkg_pass += hsum.GetBinContent(b)
                    sig_pass += hsig.GetBinContent(b)

                #print(bkg_pass)
                #print(sig_pass)
                bkgeffs[ibin] = bkg_pass/totbkg
                sigeffs[ibin] = sig_pass/totsig
                cutlist[ibin] = hsum.GetBinLowEdge(ibin)

                print "Bin low edge ",hsum.GetBinLowEdge(ibin)
                print "  bkg bin content ",hsum.GetBinContent(ibin)
                print "  sig bin content ",hsig.GetBinContent(ibin)
                print "        total bkg ",totbkg
                print "    bkg from edge ",bkg_pass
                print "        total sig ",totsig
                print "    sig from edge ",sig_pass
                #print(hsum.GetBinLowEdge(ibin))


        #print(cutlist,bkgeffs,sigeffs)
        #remove underflow bin
        bkgeffs = np.delete(bkgeffs,0)
        sigeffs = np.delete(sigeffs,0)
        cutlist = np.delete(cutlist,0)

        tg = ROOT.TGraph(hsum.GetNbinsX()-1,sigeffs,bkgeffs)
        tg.SetTitle(btag)
        tg.SetLineWidth(2)
        tg.SetLineColor(btagcols[i])
        mg.Add(tg)

        gleg.AddEntry(tg,btag)

    tc.cd()
    tc.SetLogy()
    mg.Draw("AL")
    mg.GetXaxis().SetTitle("signal efficiency")
    mg.GetXaxis().SetLimits(0,1.0)
    mg.GetYaxis().SetTitle("background efficiency")
    mg.GetYaxis().SetTitleOffset(1.1)
    mg.GetYaxis().SetLabelSize(0.025)
    mg.SetMinimum(0.0001)
    mg.SetMaximum(1.0)
    gleg.Draw()
    tc.Update()
    tc.Draw()

    #Save the plot
    pngname = gecorg.makeOutFile(hname,'btagger_roc_'+sigsamp,'.png',str(zptcut),str(hptcut),str(metcut))
    tc.SaveAs(pngname)

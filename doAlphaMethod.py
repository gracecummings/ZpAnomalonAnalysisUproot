import argparse
import ROOT
import glob
import os
import gecorg_py2 as gecorg
import numpy as np
#from datetime import date
from ROOT import kOrange, kViolet, kCyan, kGreen, kPink, kAzure, kMagenta, kBlue, kBird
from math import sqrt

if __name__=='__main__':
    #build module objects
    parser = argparse.ArgumentParser()

    #Define parser imputs
    parser.add_argument("-L","--lumi", type=float,default = 41.53, help = "integrated luminosity for scale in fb^-1")
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    #parser.add_argument("-p","--plot", type=str,help = "desired plot name to make the optimization for")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-y","--year", type=int,help = "2 digit indicator of year")
    parser.add_argument("-date","--date",help="date folder with plots to stack")
    args = parser.parse_args()

    #Get command line parameters
    lumi          = args.lumi
    sig_xsec      = args.xsec
    #released_plot = args.plot
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut
    btagwp        = args.btagwp
    year          = args.year
    plotmax       = 15.0
    
    #Samples
    #this is a hack right now, need to decide how to do sideband and signal better
    bkgfilessb = gecorg.gatherBkg('analysis_output_ZpAnomalon/'+args.date+'/sideband_only/','upout',zptcut,hptcut,metcut,btagwp,year)
    bkgfilessr = gecorg.gatherBkg('analysis_output_ZpAnomalon/'+args.date+'/signalregion_only','upout',zptcut,hptcut,metcut,btagwp,year)
    datfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/sideband_only/Run2017*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')#not changed for new naming yet

    #These are the uncertainties
    bkguncssb  = np.load('analysis_output_ZpAnomalon/'+args.date+'/sideband_only/Fall17.AllZpAnomalonBkgs_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')
    bkguncssr  = np.load('analysis_output_ZpAnomalon/'+args.date+'/signalregion_only/Fall17.AllZpAnomalonBkgs_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')

    datuncs  = np.load('analysis_output_ZpAnomalon/'+args.date+'/sideband_only/Fall17.AllZpAnomalonData_unc_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npz')

    
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]

    #Prep backgrounds and data
    bkg_colors = gecorg.colsFromPalette(bkgnames,ROOT.kLake)
    bkg_infosb = gecorg.prepBkg(bkgfilessb,bkgnames,bkg_colors,'xsects_2017.ini',lumi)
    bkg_infosr = gecorg.prepBkg(bkgfilessr,bkgnames,bkg_colors,'xsects_2017.ini',lumi)
    dat_info   = [ROOT.TFile(dat) for dat in datfiles]
    
    #Make the stacked plot
    hname = 'h_zp_jigm'
    #print bkguncssb[hname]
    #print bkguncssr[hname]
    legsb = ROOT.TLegend(0.45,0.55,0.90,0.88)
    legsr = ROOT.TLegend(0.45,0.55,0.90,0.88)

    hsbkgsb = ROOT.THStack('hsbkgsb','')
    hsbkgsr = ROOT.THStack('hsbkgsr','')
    gecorg.stackBkg(bkg_infosb,hname,hsbkgsb,legsb,plotmax,0.0)
    gecorg.stackBkg(bkg_infosr,hname,hsbkgsr,legsr,plotmax,0.0)


    #LUT with titles
    titles = {
        "h_z_pt":"Z pT",
        "h_h_pt":"Higgs pT",
        "h_h_sd":"Higgs Soft Drop Mass",
        "h_met":"pT miss",
        "h_zp_jigm":"Jigsaw Mass Estimator Z'",
        "h_nd_jigm":"Jigsaw Mass Estimator ND",
        }

    #And the data stack
    hsdat = dat_info[0].Get(hname)
    hsdat.SetStats(hsdat)
    hsdat.SetMaximum(plotmax)
    hsdat.SetMinimum(0.0)
    hsdat.SetMarkerStyle(8)
    hsdat.SetBinErrorOption(1)
    legsb.AddEntry(hsdat,"Data")
    for d,datf in enumerate(dat_info[1:]):
        hdat = datf.Get(hname)
        hdat.SetStats(0)
        hdat.SetMaximum(plotmax)
        hdat.SetMinimum(0.0)
        hsdat.Add(hdat)

    #Define the Bkg Estimate Hist
    hbest = ROOT.TH1F("hbest","Bkg Est",100,500,5000)
    hbest.SetMarkerStyle(8)
    hbest.SetMarkerColor(4)
    hbest.SetLineColor(4)
        
    #Alpha Canvas
    tc = ROOT.TCanvas("tc",hname,1200,800)
    #Pad with SR MC and alpha ratio extrapolation 
    p1 = ROOT.TPad("p1","stack_"+hname,0.0,0.4,0.5,1.0)
    #p1.SetLogy()
    #p1.SetBottomMargin(0)
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    #Pad with the ratio of extrapolation to MC
    p2 = ROOT.TPad("p2","alpha ratio",0.0,0.0,0.5,0.4)
    #p2.SetTopMargin(0)
    p2.SetRightMargin(.05)
    p2.SetLeftMargin(0.15)
    p2.SetBottomMargin(0.2)
    #Pad with alpha ratios
    p3 = ROOT.TPad("p3","alpha_"+hname,0.5,0.4,1.0,1.0)
    #p1.SetBottomMargin(0)
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)

    #Alpha Ratio
    #Access Sum bkgmc stacks
    hsumsb = hsbkgsb.GetStack().Last()
    hsumsr = hsbkgsr.GetStack().Last()
    hsumsb.SetBinErrorOption(1)
    hsumsr.SetBinErrorOption(1)

    #print "max of sideband ",hsumsb.GetBinCenter(hsumsb.GetMaximumBin())
    #print "max of signal region ",hsumsr.GetBinCenter(hsumsr.GetMaximumBin())

    print "events in sideband: ",hsumsb.Integral()
    print "events in signal r: ",hsumsr.Integral()
    #Test of Normalization
    #hsumsb.Scale(1/hsumsb.Integral())
    #hsumsr.Scale(1/hsumsr.Integral())
    
    #initilize
    max_max    = 0.0
    binlist    = np.zeros(hsumsb.GetNbinsX())
    alphalist  = np.zeros(hsumsb.GetNbinsX())
    alphaerrl  = np.zeros(hsumsb.GetNbinsX())
    ratiolist  = np.zeros(hsumsb.GetNbinsX())
    rerrlist   = np.zeros(hsumsb.GetNbinsX())
    fill       = np.zeros(hsumsb.GetNbinsX())

    #Calculate bin-by-bin alpha ratio and extrapolations
    for ibin in range(hsumsb.GetNbinsX()):
        bincen = hsumsb.GetBinCenter(ibin)
        bkgsb  = hsumsb.GetBinContent(ibin)
        bkgsr  = hsumsr.GetBinContent(ibin)
        #sbunc  = hsumsb.GetBinError(ibin)
        #srunc  = hsumsr.GetBinError(ibin)
        data   = hsdat.GetBinContent(ibin)
        #datunc = hsdat.GetBinError(ibin)
        if bkgsb !=0:
            alpha  = bkgsr/bkgsb
        else:
            alpha = -1
        if ibin != 0:
            hsdat.SetBinError(ibin,datuncs[hname][ibin-1])
            datunc = datuncs[hname][ibin-1]
            sbunc  = bkguncssb[hname][ibin-1]
            srunc  = bkguncssr[hname][ibin-1]
            #if data != 0:
            bkgest = alpha*data
            if bkgsb != 0 and bkgsr != 0:
                alphaunc = alpha*sqrt((sbunc/bkgsb)**2+(srunc/bkgsr)**2)
                if alpha != 0:
                    estunc   = bkgest*sqrt((alphaunc/alpha)**2+(datunc/data)**2)
                else:
                    estunc = 0.0
                ratio    = bkgest/bkgsr
                if bkgest != 0:
                    ratiounc = ratio*sqrt((estunc/bkgest)**2+(srunc/bkgsr)**2)
                else:
                    ratiounc = 0
            else:
                alphaunc = 0#alpha*sqrt((sbunc/bkgsb)**2+(srunc/bkgsr)**2)
                if alpha != 0:
                    estunc   = bkgest*sqrt((alphaunc/alpha)**2)
                else:
                    estunc = 0
                ratio    = -1
                ratiounc = 0
        else:
            alphaunc = -1
            ratio = -1
            ratiounc = -1
            bkgest = -1
            estunc = -1

        #fill lists for plotting
        alphalist[ibin] = alpha
        binlist[ibin]   = bincen
        alphaerrl[ibin] = alphaunc
        ratiolist[ibin]  = ratio
        rerrlist[ibin]   = ratiounc 
        hbest.SetBinContent(ibin,bkgest)
        hbest.SetBinError(ibin,estunc)

    #remove underflow bin
    alphalist  = np.delete(alphalist,0)
    binlist    = np.delete(binlist,0)
    alphaerrl  = np.delete(alphaerrl,0)
    ratiolist  = np.delete(ratiolist,0)
    rerrlist  = np.delete(rerrlist,0)
    
    #Build the graphs
    ta = ROOT.TGraphErrors(hsumsb.GetNbinsX()-1,binlist,alphalist,fill,alphaerrl)
    ta.SetTitle("")
    ta.SetMarkerStyle(8)
    ta.SetMarkerColor(4)
    ta.SetLineColor(4)

     
    #beauty aspects
    ta.SetTitle("")
    ta.GetXaxis().SetTitle("RJR Contraboost Zp Mass Estimator")
    ta.GetXaxis().SetTitleSize(0.05)
    ta.GetXaxis().SetLabelSize(0.035)
    ta.GetXaxis().SetLimits(binlist[0],binlist[-1]+hsumsb.GetBinWidth(1))
    #y axis
    ta.GetYaxis().SetTitle("alpha ratio")
    ta.GetYaxis().SetTitleSize(0.05)
    ta.GetYaxis().SetTitleOffset(.8)
    ta.GetYaxis().SetLabelSize(0.035)
    ta.SetMinimum(0.5)
    ta.SetMaximum(1.5)

    tr = ROOT.TGraphErrors(hsumsb.GetNbinsX()-1,binlist,ratiolist,fill,rerrlist)
    tr.SetTitle("")
    tr.SetMarkerStyle(8)
    tr.SetMarkerColor(4)
    tr.SetLineColor(4)
     
    #beauty aspects
    tr.SetTitle("")
    tr.GetXaxis().SetTitle("RJR Contraboost Zp Mass Estimator")
    tr.GetXaxis().SetTitleSize(0.07)
    tr.GetXaxis().SetLabelSize(0.05)
    tr.GetXaxis().SetLimits(binlist[0],binlist[-1]+hsumsb.GetBinWidth(1))
    #y axis
    tr.GetYaxis().SetTitle("(estimated bkg)/MC")
    tr.GetYaxis().SetTitleSize(0.07)
    tr.GetYaxis().SetTitleOffset(.7)
    tr.GetYaxis().SetLabelSize(0.05)
    tr.SetMinimum(0)
    tr.SetMaximum(1.5)
    l = ROOT.TLine(binlist[0],1,binlist[-1]+hsumsb.GetBinWidth(1),1)

    #let's Draw Stuff
    #Pad with SR Bkg MC
    p1.Draw()
    p1.cd()
    hsbkgsr.Draw("HIST")#add PFC for palette drawing
    #hsumsr.Draw("E3")
    hsbkgsr.GetXaxis().SetTitle("RJR Contraboost Zp Mass Estimator")
    hsbkgsr.GetXaxis().SetTitleSize(0.05)
    hsbkgsr.GetYaxis().SetTitle("Events")
    hsbkgsr.GetYaxis().SetTitleSize(0.05)
    hbest.Draw("HISTSAMEPE1")
    legsr.AddEntry(hbest,"Bkg Est")
    legsr.SetBorderSize(0)
    legsr.Draw()
    tc.Modified()

    #Make the pad with the ratio plot
    tc.cd()
    p2.Draw()
    p2.cd()
    tr.Draw("AP")
    l.Draw()
    
    tc.cd()
    p3.Draw()
    p3.cd()
    ta.Draw("AP")
    
    tc.Update()

    #save the canvas
    pngname = gecorg.makeOutFile(hname,'alpha','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(pngname)

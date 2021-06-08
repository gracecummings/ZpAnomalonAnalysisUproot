import subprocess
from datetime import date

if __name__=='__main__':

    #steps to run
    #assumes you have run the whole thing at the start of the day
    #steps = {"selections":True,"uncs":True,"ratios":True,"opts":True}
    steps = {"topiary":False,"selections":False,"uncs":False,"ratios":True,"opts":False,"cutflow":False}
    
    #cut list, Zpt, Hpt, met,btagger,btagwp
    cutlist = [#['0.0','300.0','0.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['150.0','300.0','0.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['0.0','300.0','200.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ['150.0','300.0','200.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ]

    lumi = "41.53"

    #year specifics, year as 2 digit end, integrated luminosity
    eras = [["17","41.53"],
            ["18","59.74"]
            ]

    plots = ['h_z_pt']#['h_h_pt','h_z_pt','h_met','h_nd_jigm','h_zp_jigm','h_h_sd','h_btag']

    #2021-06-07 for reclustering, leading muons, and gen info
    #topiary sample list: dateforfolder, samplename
    samplelist = [#['2021-06-07_reclusteredJets','Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                  #['2021-06-07_reclusteredJets','Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-03-26','Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                  #['2021-03-26','Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-03-26','Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-03-26','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-06-03','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-06-02_nonreclus_topiary','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-03-26','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-06-03','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-06-02_nonreclus_topiary','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-03-26','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-06-03','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-06-02_nonreclus_topiary','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  #['2021-03-26','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  #['2021-06-03','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],\
                  #['2021-06-02_nonreclus_topiary','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  #['2021-03-26','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  #['2021-06-03','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  #['2021-06-07_reclusteredJets','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  #['2021-06-07_reclusteredJets_alternate_reclustopiary_nominalselections','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  #['2021-06-02_nonreclus_topiary','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'], 
                  #['2021-06-04','ZpAnomalonHZ_UFO-Zp3000-ND500-NS1'],
                  #['2021-06-04','Run2017B-31Mar2018-v1.SingleMuon'],
                  #['2021-06-04','Run2017C-31Mar2018-v1.SingleMuon'],
                  #['2021-06-04','Run2017D-31Mar2018-v1.SingleMuon'],
                  #['2021-06-04','Run2017E-31Mar2018-v1.SingleMuon'],
                  #['2021-06-04','Run2017F-31Mar2018-v1.SingleMuon'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets_alternate_reclustopiary_nominalselections','Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-05-29','Run2018C-17Sep2018-v1.SingleMuon'],
                  #['2021-05-29','Run2018B-17Sep2018-v1.SingleMuon'],
                  #['2021-05-29','Run2018A-17Sep2018-v1.SingleMuon'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-06-07_reclusteredJets','Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-05-29','Autumn18.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                  #['2021-05-29','Autumn18.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-05-29','Autumn18.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  #['2021-03-26','Run2017B-31Mar2018-v1.SingleMuon'],
                  #['2021-03-26','Run2017C-31Mar2018-v1.SingleMuon'],
                  #['2021-03-26','Run2017D-31Mar2018-v1.SingleMuon'],
                  #['2021-03-26','Run2017E-31Mar2018-v1.SingleMuon'],
                  #['2021-03-26','Run2017F-31Mar2018-v1.SingleMuon'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2021-06-07_reclusteredJets','Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-03-26','Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  #['2021-04-27','Run2018C-17Sep2018-v1.SingleMuon'],
                  #['2021-04-27','Run2018B-17Sep2018-v1.SingleMuon'],
                  #['2021-04-27','Run2018A-17Sep2018-v1.SingleMuon'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  #['2021-04-27','Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia'],
                  ]

    if steps["topiary"]:
        for samp in samplelist:
            subprocess.run(["python","runTopiary.py","-s",samp[1],"-y","2017"])
    
    for cut in cutlist:
        print("Doing ZpT cut {0}, HpT cut {1}, MET cut {2}, btag wp {3}".format(cut[0],cut[1],cut[2],cut[4]))
        #do selections
        if steps["selections"]:
            for samp in samplelist:
                subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0]])
                subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-sr","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-c","True"])

        for era in eras:
            print("   Beginning plottng and analysis for year {0}, with a luminosity of {1}".format(era[0],era[1]))
            if steps["uncs"]:
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","totalr"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","sideband"])
                subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_reclusteredJets","-y",era[0],"-r","sideband"])
                subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_reclusteredJets","-y",era[0],"-r","signalr"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_nonreclusteredJets","-y",era[0],"-r","sideband"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_nonreclusteredJets","-y",era[0],"-r","signalr"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","signalr"])

        #stack all  
        if steps["ratios"]:
            #subprocess.run(["python","stackAll.py","-L",era[1],"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today())])
            #subprocess.run(["python","stackAll.py","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-r","totalr"])
            subprocess.run(["python","stackAll.py","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_reclusteredJets","-r","signalr"])
            subprocess.run(["python","stackAll.py","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_reclusteredJets","-r","sideband"])
            #subprocess.run(["python","stackAll.py","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-06-07_reclusteredJets","-r","signalr"])

        #Optimization Plots
        if steps["opts"]:
            for plot in plots:
                subprocess.run(["python","stackForOptimization.py","-L",era[1],"-x","10.0","-p",plot,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])
                #subprocess.run(["python","stackForOptimization.py","-L",lumi,"-x","100.0","-p",plot,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",'2021-02-03'])

        if steps["cutflow"]:
            print("Creating cutflow table")
            #subprocess.run(["python","doCutFlow.py","-L",era[1],"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])
            subprocess.run(["python","doCutFlowSig.py","-L",era[1],"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])
            

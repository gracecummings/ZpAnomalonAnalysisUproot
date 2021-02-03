import subprocess
from datetime import date

if __name__=='__main__':

    #steps to run
    #assumes you have run the whole thing at the start of the day
    #steps = {"selections":True,"uncs":True,"ratios":True,"opts":True}
    steps = {"selections":False,"uncs":False,"ratios":True,"opts":False,"btagcomp":False}
    
    #cut list, Zpt, Hpt, met,btagger,btagwp
    cutlist = [['200.0','300.0','300.0','DeepMassDecorrelTagZHbbvsQCD','0.0'],
               ['200.0','300.0','300.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ['200.0','300.0','300.0','DeepMassDecorrelTagZHbbvsQCD','0.97'],
               #['200.0','300.0','0.0','DeepMassDecorrelTagZHbbvsQCD','0.97'],
               #['200.0','300.0','200.0','DeepMassDecorrelTagZHbbvsQCD','0.97'],
               ]

    lumi = "41.53"

    plots = ['h_h_pt','h_z_pt','h_met','h_nd_jigm','h_zp_jigm','h_h_sd']

    #topiary sample list: dateforfolder, samplename
    samplelist = [['2020-12-29','Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ['2020-12-29','Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                  ['2020-12-29','Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ['2020-12-29','Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp2000-ND175-NS1'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  ['2020-12-29','ZpAnomalonHZ_UFO-Zp3000-ND500-NS1'],
                  ['2021-01-05','Run2017B-31Mar2018-v1.SingleMuon'],
                  ['2021-01-05','Run2017C-31Mar2018-v1.SingleMuon'],
                  ['2021-01-05','Run2017D-31Mar2018-v1.SingleMuon'],
                  ['2021-01-05','Run2017E-31Mar2018-v1.SingleMuon'],
                  ['2021-01-05','Run2017F-31Mar2018-v1.SingleMuon'],
                  ]

    
    for cut in cutlist:

        #do selections
        if steps["selections"]:
            for samp in samplelist:
                subprocess.run(["python3","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0]])

        #do uncertainites
        if steps["uncs"]:
            subprocess.run(["python3","doStackedUncertainty.py","-L",lumi,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today())])

        #stack all
        if steps["ratios"]:
            subprocess.run(["python","stackAll.py","-L",lumi,"-x","100.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today())])

        #Optimization Plots
        if steps["opts"]:
            for plot in plots:
                subprocess.run(["python","stackForOptimization.py","-L",lumi,"-x","100.0","-p",plot,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today())])

    #btagger comps
    #if steps["btagcomp"]:
    #    print("got btagger")
            

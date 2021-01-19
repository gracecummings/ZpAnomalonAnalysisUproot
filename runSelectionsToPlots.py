import subprocess
from datetime import date

if __name__=='__main__':

    #cut list, Zpt, Hpt, met
    cutlist = [['200.0','250.0','0.0'],
               #[200.0,250.0,300.0],
               #[200.0,0.0,300.0],
               #[200,300.0,300.0],
               #[0.0,300.0,300.0]
               ]

    lumi = "41.53"

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
        for samp in samplelist:
            subprocess.run(["python3","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-date",samp[0]])

        #do uncertainites
        subprocess.run(["python3","doStackedUncertainty.py","-L",lumi,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-date",str(date.today())])

        #stack all
        subprocess.run(["python","stackAll.py","-L",lumi,"-x","100.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-date",str(date.today())])
        

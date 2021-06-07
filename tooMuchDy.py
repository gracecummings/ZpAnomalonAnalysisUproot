import gecorg as go
import ROOT

goodpath  = 'analysis_output_ZpAnomalon/2021-04-28/'
badpath   = 'analysis_output_ZpAnomalon/2021-06-07_reclusteredJets'

zptcut  = '150.0'
hptcut  = '300.0'
metcut  = '200.0'
btagwp  = '0.8'

goodbkgs = go.backgrounds(goodpath,zptcut,hptcut,metcut,btagwp)
badbkgs  = go.backgrounds(badpath,zptcut,hptcut,metcut,btagwp)                                          
g17dysb = goodbkgs.f17dyjetsb
g17dysr = goodbkgs.f17dyjetsr
g18dysb = goodbkgs.a18dyjetsb
g18dysr = goodbkgs.a18dyjetsr

tf1 = ROOT.TFile(goodbkgs.f17dyjetsb[0])
empty = tf1.Get('h_zp_jigm')
empty.Reset("ICESM")

tf2 = ROOT.TFile(badbkgs.f17dyjetsb[0])
empty1 = tf2.Get('h_zp_jigm')
empty1.Reset("ICESM")


goodhist,goodinfo = goodbkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
badhist,badinfo = badbkgs.getAddedHist(empty1,"DYJetsToLL","sb","h_zp_jigm")

print(len(goodinfo))
print(len(badinfo))

for i,inf in enumerate(goodinfo):
    #print(badinfo[i])
    ginfos = inf.split(" ")
    binfos = badinfo[i].split(" ")

    print(ginfos[0],binfos[0])



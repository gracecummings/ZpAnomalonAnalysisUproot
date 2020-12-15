import ROOT
import glob

if __name__=="__main__":
    f = 'dyjets_1200to2500_skim.root'
    outFile = 'test_topiary.root'
    inChain = ROOT.TChain("PreSelection")
    inChain.Add(f)
    tf = ROOT.TFile.Open(f)
    origevnts = tf.Get("hnevents").GetBinContent(1)
    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    ROOT.gSystem.Load('TreeMakerTopiary_C')
    topiary = ROOT.TreeMakerTopiary(inChain)
    topiary.Loop(outFile,origevnts)
    #inputfs = glob.glob(path)
    #origevnts = [ROOT.TFile.Open(f).Get("hnevents").GetBinContent(1) for f in inputfs]
    #totorig = sum(origevnts)

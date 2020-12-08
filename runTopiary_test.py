import ROOT

if __name__=="__main__":
    f = 'dyjets_1200to2500_skim.root'
    outFile = 'test_topiary.root'
    inChain = ROOT.TChain("PreSelection")
    inChain.Add(f)
    ROOT.gSystem.CompileMacro("TreeMakerTopiary.C","g0ck")
    ROOT.gSystem.Load('TreeMakerTopiary_C')
    topiary = ROOT.TreeMakerTopiary(inChain)
    topiary.Loop(outFile)

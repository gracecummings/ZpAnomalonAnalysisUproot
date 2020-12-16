#define TreeMakerTopiary_cxx
#include "TreeMakerTopiary.h"
#include "RestFrames/RestFrames.hh"
#include <TLeaf.h>
#include <iostream>

RestFrames::RFKey ensure_autoload(1);
using namespace RestFrames;

void TreeMakerTopiary::Loop(std::string outputFileName, float totalOriginalEvents)
{
//   In a ROOT session, you can do:
//      root> .L TreeMakerTopiary.C
//      root> TreeMakerTopiary t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("JetsAK8Clean*",1);
   fChain->SetBranchStatus("Muons*",1);
   fChain->SetBranchStatus("METclean*",1);
   fChain->SetBranchStatus("SelectedMuons*",1);
   fChain->SetBranchStatus("ZCandidates",1);


   //Initialize Stuff
   TLorentzVector hCandidate;
   TLorentzVector ZCandidate;
   double ZCandidate_pt;
   double ZCandidate_phi;
   double ZCandidate_eta;
   double ZCandidate_m;
   double hCandidate_pt;
   double hCandidate_phi;
   double hCandidate_eta;
   double hCandidate_m;
   double hCandidate_sd;
   double hCandidate_dmdhbbvqcd;
   double hCandidate_dmdzbbvqcd;
   double hCandidate_dmdzhbbvqcd;
   double hCandidate_middb;
   double mEstZp;
   double mEstND;
   double mEstNS;
   
   //Define the skimmed skim  output file and tree
   TFile* trimFile = new TFile(outputFileName.c_str(),"recreate");
   TTree* trimTree = fChain->CloneTree(0);
   TH1F*  hnskimed = new TH1F("hnskimed","number of events at skim level",1,0,1);
   TH1F*  hnorigevnts = new TH1F("hnorigevnts","original number of events, preskim",1,0,1);
   TBranch *hCand     = trimTree->Branch("hCandidate","TLorentzVector",&hCandidate);
   TBranch *hCand_pt  = trimTree->Branch("hCandidate_pt",&hCandidate_pt,"hCandidate_pt/D");
   TBranch *hCand_phi = trimTree->Branch("hCandidate_phi",&hCandidate_phi,"hCandidate_phi/D");
   TBranch *hCand_eta = trimTree->Branch("hCandidate_eta",&hCandidate_eta,"hCandidate_eta/D");
   TBranch *hCand_m   = trimTree->Branch("hCandidate_m",&hCandidate_m,"hCandidate_m/D");
   TBranch *hCand_sd  = trimTree->Branch("hCandidate_sd",&hCandidate_sd,"hCandidate_sd/D");
   TBranch *hCand_dmdhbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagHbbvsQCD",&hCandidate_dmdhbbvqcd,"hCandidate_dmdhbbvqcd/D");
   TBranch *hCand_dmdzbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagZbbvsQCD",&hCandidate_dmdzbbvqcd,"hCandidate_dmdzbbvqcd/D");
   TBranch *hCand_dmdzhbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagZHbbvsQCD",&hCandidate_dmdzhbbvqcd,"hCandidate_dmdzhbbvqcd/D");
   TBranch *hCand_middb  = trimTree->Branch("hCandidate__pfMassIndependentDeepDoubleBvLJetTagsProbHbb",&hCandidate_middb,"hCandidate_middb/D");
   TBranch *ZCand     = trimTree->Branch("ZCandidate","TLorentzVector",&ZCandidate);
   TBranch *ZCand_pt  = trimTree->Branch("ZCandidate_pt",&ZCandidate_pt,"ZCandidate_pt/D");
   TBranch *ZCand_phi = trimTree->Branch("ZCandidate_phi",&ZCandidate_phi,"ZCandidate_phi/D");
   TBranch *ZCand_eta = trimTree->Branch("ZCandidate_eta",&ZCandidate_eta,"ZCandidate_eta/D");
   TBranch *ZCand_m   = trimTree->Branch("ZCandidate_m",&ZCandidate_m,"ZCandidate_m/D");
   TBranch *ZpMest    = trimTree->Branch("ZPrime_mass_est",&mEstZp,"mEstZp/D");
   TBranch *NDMest    = trimTree->Branch("ND_mass_est",&mEstND,"mEstND/D");
   TBranch *NSMest    = trimTree->Branch("NS_mass_est",&mEstNS,"mEstNS/D");
   hnskimed->SetBinContent(1,nentries);
   hnorigevnts->SetBinContent(1,totalOriginalEvents);

   float zmwinlow = 70.;
   float zmwinhi  = 110.;
   float hptcut   = 250.;

   //Recursive Jigsaw Part
   LabRecoFrame         LAB("LAB","LAB");
   DecayRecoFrame       Zp("Zp","Z'");
   DecayRecoFrame       ND("ND","N_{D}");
   DecayRecoFrame       NDbar("NDbar","N_{Dbar}");
   VisibleRecoFrame     Z("Z","Z");
   InvisibleRecoFrame   NS("NS","N_{S}");
   VisibleRecoFrame     h("h","h");
   InvisibleRecoFrame   NSbar("NSbar","Z_{Sbar}");

   LAB.SetChildFrame(Zp);
   Zp.AddChildFrame(ND);
   Zp.AddChildFrame(NDbar);
   ND.AddChildFrame(Z);
   ND.AddChildFrame(NS);
   NDbar.AddChildFrame(h);
   NDbar.AddChildFrame(NSbar);

   LAB.InitializeTree();
   
   // Invisible Group
   InvisibleGroup INV("INV","NS NS Jigsaws");
   INV.AddFrame(NS);
   INV.AddFrame(NSbar);
   
   // Set NS NS~ mass equal to Z h mass
   SetMassInvJigsaw NSNSM("NSNSM", "M_{NSNS} = m_{Zh}");
   INV.AddJigsaw(NSNSM);
   
   SetRapidityInvJigsaw NSNSR("NSNSR", "#eta_{NSNS} = #eta_{ZH}");
   INV.AddJigsaw(NSNSR);
   NSNSR.AddVisibleFrames(LAB.GetListVisibleFrames());
   
   //MinMassesSqInvJigsaw MinMND("MinMND","min M_{D}, M_{ND}= M_{NDbar}",2);
   ContraBoostInvJigsaw MinMND("MinMND","min M_{ND}, M_{ND}= M_{NDbar}");
   INV.AddJigsaw(MinMND);
   MinMND.AddVisibleFrame(Z, 0);
   MinMND.AddVisibleFrame(h, 1);
   MinMND.AddInvisibleFrame(NS, 0);
   MinMND.AddInvisibleFrame(NSbar, 1);

   LAB.InitializeAnalysis();
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      //Define some boos to signal events to write
      bool passZ   = false;
      bool passh   = false;
      bool passMET = false;

      //A counter, for my sanity
      if (jentry%20000 == 0) {
      	std::cout<<"    analyzing event "<<jentry<<std::endl;
      }

      //Z Candidate Build
      unsigned int nZs = ZCandidates->size();
      TLorentzVector theZ;
      double baseZdiff = 99999;
      if (nZs > 0) {
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidates->begin(); zit != ZCandidates->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() >= zmwinlow) && (zit->M() <= zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Phi(),zit->Eta(),zit->M());
	    passZ = true;
	  }
	}
      }

      //Higgs Candidate Build
      unsigned long nfat = JetsAK8Clean->size();
      TLorentzVector theh;
      TLorentzVector fat;
      double basehdiff = 99999;
      double hsd;
      double fsd;
      double hdmdhbbvqcd;
      double hdmdzbbvqcd;
      double hdmdzhbbvqcd;
      double hmiddb;
      
      if (nfat > 0) {
	for (unsigned long i =0; i < nfat; ++i) {
	  fat = JetsAK8Clean->at(i);
	  fsd = JetsAK8Clean_softDropMass->at(i);
	  double masshdiff = std::abs(125.18 - fsd);
	  if ((masshdiff < basehdiff) && (fat.Pt() > hptcut)) {
	    basehdiff = masshdiff;
	    theh = fat;
	    hsd = fsd;
	    hdmdhbbvqcd  = JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD->at(i);
	    hdmdzbbvqcd  = JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD->at(i);
	    hdmdzhbbvqcd = JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD->at(i);
	    hmiddb = JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i);
	    passh = true;
	  }
	}
      }

      //MET
      double ptmiss     = fChain->GetLeaf("METclean")->GetValue(0);
      double ptmiss_phi = fChain->GetLeaf("METPhiclean")->GetValue();
      double ptmiss_px  = ptmiss*std::cos(ptmiss_phi);
      double ptmiss_py  = ptmiss*std::sin(ptmiss_phi);
      TVector3 met3     = TVector3(ptmiss_px,ptmiss_py,0.0);
      
      //recursive jigsaw
      LAB.ClearEvent();
      INV.SetLabFrameThreeVector(met3);
      Z.SetLabFrameFourVector(theZ);
      h.SetLabFrameFourVector(theh);
      LAB.AnalyzeEvent();
      mEstZp = Zp.GetMass();
      mEstND = ND.GetMass();
      mEstNS = NS.GetMass();
      
      if (passZ) {
	ZCandidate = theZ;
	ZCandidate_pt  = theZ.Pt();
	ZCandidate_phi = theZ.Phi();
	ZCandidate_eta = theZ.Eta();
	ZCandidate_m   = theZ.M();
      }

      if (passh) {
	hCandidate = theh;
	hCandidate_pt  = theh.Pt();
	hCandidate_phi = theh.Phi();
	hCandidate_eta = theh.Eta();
	hCandidate_m   = theh.M();
	hCandidate_sd  = hsd;
	hCandidate_dmdhbbvqcd = hdmdhbbvqcd;
	hCandidate_dmdzbbvqcd = hdmdzbbvqcd;
	hCandidate_dmdzhbbvqcd = hdmdzhbbvqcd;
	hCandidate_middb = hmiddb;
      }

      
      //debug
      //if (jentry == 2000) {
      //break;
      //}

      //Fill the Tree
      if (Cut(ientry) < 0) continue;
      if (passZ && passh) {
	trimTree->Fill();
	}
      
   }

   trimFile->Write();
   trimFile->Close();   
   //std::cout<<"trimmed to "<< passEvents <<" events"<<std::endl;
   std::cout<<"Completed your topiary garden, hopefully your tastes have not changed."<<std::endl;

}

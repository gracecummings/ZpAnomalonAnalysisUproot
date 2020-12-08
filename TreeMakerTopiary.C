#define TreeMakerTopiary_cxx
#include "TreeMakerTopiary.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void TreeMakerTopiary::Loop(std::string outputFileName)
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
   
   //Define the skimmed skim  output file and tree
   TFile* trimFile = new TFile(outputFileName.c_str(),"recreate");
   TTree* trimTree = fChain->CloneTree(0);
   TH1F*  hnskimed = new TH1F("hnskimed","first skimmed entires entries",1,0,1);
   TBranch *hCand     = trimTree->Branch("hCandidate",&hCandidate,"hCandidate");
   TBranch *ZCand     = trimTree->Branch("ZCandidate",&ZCandidate,"ZCandidate");
   TBranch *ZCand_pt  = trimTree->Branch("ZCandidate_pt",&ZCandidate_pt,"ZCandidate_pt");
   TBranch *ZCand_phi = trimTree->Branch("ZCandidate_phi",&ZCandidate_phi,"ZCandidate_phi");
   TBranch *ZCand_eta = trimTree->Branch("ZCandidate_eta",&ZCandidate_eta,"ZCandidate_eta");
   TBranch *ZCand_m   = trimTree->Branch("ZCandidate_m",&ZCandidate_eta,"ZCandidate_m");
   hnskimed->SetBinContent(1,nentries);


   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      //Define some boos to signal events to write
      bool passZ   = false;
      bool passh   = false;
      bool passMET = false;

      //A counter, for my sanity
      //if (jentry%20000 == 0) {
      //	std::cout<<"skimming event "<<jentry<<std::endl;
      //}

      //Z Candidate Check
      unsigned int nZs = ZCandidates->size();
      TLorentzVector theZ;
      double baseZdiff = 99999;
      if (nZs > 0) {
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidates->begin(); zit != ZCandidates->end(); ++zit) {
	  double massZdiff = std::abs(91.18 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() >= 70.) && (zit->M() <= 110.)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Phi(),zit->Eta(),zit->M());
	    std::cout<<"The Zit  pt "<<zit->Pt()<<std::endl;
	    std::cout<<"   The Z pt "<<theZ.Pt()<<std::endl;
	    passZ = true;
	  }
	}
      }
	
      if (passZ) {
	//ZCandidate = theZ;
	ZCandidate_pt  = theZ.Pt();
	ZCandidate_phi = theZ.Phi();
	ZCandidate_eta = theZ.Eta();
	ZCandidate_m   = theZ.M();
	//ZCand->Fill();
	ZCand_pt->Fill();
	ZCand_phi->Fill();
	ZCand_eta->Fill();
	ZCand_m->Fill();
      }
      
      //debug
      if (jentry == 2000) {
	break;
      }

      //Fill the Tree
      if (Cut(ientry) < 0) continue;
      if (passZ) {
	trimTree->Fill();
      }
      
   }

   trimFile->Write();
   trimFile->Close();   
   //std::cout<<"trimmed to "<< passEvents <<" events"<<std::endl;
   std::cout<<"DONE"<<std::endl;

}

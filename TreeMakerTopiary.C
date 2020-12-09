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
   double hCandidate_pt;
   double hCandidate_phi;
   double hCandidate_eta;
   double hCandidate_m;
   
   //Define the skimmed skim  output file and tree
   TFile* trimFile = new TFile(outputFileName.c_str(),"recreate");
   TTree* trimTree = fChain->CloneTree(0);
   TH1F*  hnskimed = new TH1F("hnskimed","first skimmed entires entries",1,0,1);
   TBranch *hCand     = trimTree->Branch("hCandidate","TLorentzVector",&hCandidate);
   TBranch *hCand_pt  = trimTree->Branch("hCandidate_pt",&hCandidate_pt,"hCandidate_pt/D");
   TBranch *hCand_phi = trimTree->Branch("hCandidate_phi",&hCandidate_phi,"hCandidate_phi/D");
   TBranch *hCand_eta = trimTree->Branch("hCandidate_eta",&hCandidate_eta,"hCandidate_eta/D");
   TBranch *hCand_m   = trimTree->Branch("hCandidate_m",&hCandidate_m,"hCandidate_m/D");
   TBranch *ZCand     = trimTree->Branch("ZCandidate","TLorentzVector",&ZCandidate);
   TBranch *ZCand_pt  = trimTree->Branch("ZCandidate_pt",&ZCandidate_pt,"ZCandidate_pt/D");
   TBranch *ZCand_phi = trimTree->Branch("ZCandidate_phi",&ZCandidate_phi,"ZCandidate_phi/D");
   TBranch *ZCand_eta = trimTree->Branch("ZCandidate_eta",&ZCandidate_eta,"ZCandidate_eta/D");
   TBranch *ZCand_m   = trimTree->Branch("ZCandidate_m",&ZCandidate_m,"ZCandidate_m/D");
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
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() >= 70.) && (zit->M() <= 110.)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Phi(),zit->Eta(),zit->M());
	    //std::cout<<"The Zit  pt "<<zit->Pt()<<std::endl;
	    //std::cout<<"   The Z pt "<<theZ.Pt()<<std::endl;
	    passZ = true;
	  }
	}
      }

      //Higgs Candidates
      unsigned int nfat = JetsAK8Clean->size();
      TLorentzVector theh;
      double basehdiff = 99999;
      if (nfat > 0) {
	std::vector<TLorentzVector>::iterator fit;
	for (fit = JetsAK8Clean->begin(); fit != JetsAK8Clean->end(); ++fit) {
	  double masshdiff = std::abs(125.18 - fit->M());
	  if ((masshdiff < basehdiff) && (fit->Pt() > 250.)) {
	    basehdiff = masshdiff;
	    theh.SetPtEtaPhiM(fit->Pt(),fit->Phi(),fit->Eta(),fit->M());
	    //std::cout<<"The fit  pt "<<fit->Pt()<<std::endl;
	    //std::cout<<"   The h pt "<<theh.Pt()<<std::endl;
	    passh = true;
	  }
	}
      }
      
	
      if (passZ) {
	ZCandidate = theZ;
	ZCandidate_pt  = theZ.Pt();
	ZCand_pt->Fill();
	ZCandidate_phi = theZ.Phi();
	ZCandidate_eta = theZ.Eta();
	ZCandidate_m   = theZ.M();
      }

      if (passh) {
	hCandidate = theh;
	hCandidate_pt  = theh.Pt();
	//std::cout<<"   The h pt "<<hCandidate_pt<<std::endl;
	hCandidate_phi = theh.Phi();
	hCandidate_eta = theh.Eta();
	hCandidate_m   = theh.M();
      }
      //debug
      if (jentry == 200) {
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

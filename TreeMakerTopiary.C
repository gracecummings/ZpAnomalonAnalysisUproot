#define TreeMakerTopiary_cxx
#include "TreeMakerTopiary.h"
#include "RestFrames/RestFrames.hh"
#include <TLeaf.h>
#include <iostream>

RestFrames::RFKey ensure_autoload(1);
using namespace RestFrames;

void TreeMakerTopiary::Loop(std::string outputFileName, float totalOriginalEvents, int sampleType,int year)
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("TriggerPass",1);
   fChain->SetBranchStatus("JetsAK8Clean*",1);
   fChain->SetBranchStatus("Muons*",1);
   fChain->SetBranchStatus("METclean*",1);
   fChain->SetBranchStatus("METPhiclean",1);
   fChain->SetBranchStatus("SelectedMuons*",1);
   fChain->SetBranchStatus("ZCandidates",1);
   fChain->SetBranchStatus("eeBadScFilter",1);

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
   float  evntw = 1;

      //Define the skimmed skim  output file and tree
   TFile* trimFile = new TFile(outputFileName.c_str(),"recreate");
   TTree* trimTree = fChain->CloneTree(0);
   TH1F*  hnskimed = new TH1F("hnskimed","number of events at skim level",1,0,1);
   TH1F*  hnorigevnts = new TH1F("hnorigevnts","original number of events, preskim",1,0,1);
   TH1F*  htrigpass= new TH1F("htrigpass","trigger pass",1,0,1);
   TH1F*  hZpass= new TH1F("hZpass","Z pass",1,0,1);
   TH1F*  hHpass= new TH1F("hHpass","h pass",1,0,1);
   TH1F*  hpass= new TH1F("hpass","passing all req",1,0,1);
   TBranch *hCand     = trimTree->Branch("hCandidate","TLorentzVector",&hCandidate);
   TBranch *hCand_pt  = trimTree->Branch("hCandidate_pt",&hCandidate_pt,"hCandidate_pt/D");
   TBranch *hCand_phi = trimTree->Branch("hCandidate_phi",&hCandidate_phi,"hCandidate_phi/D");
   TBranch *hCand_eta = trimTree->Branch("hCandidate_eta",&hCandidate_eta,"hCandidate_eta/D");
   TBranch *hCand_m   = trimTree->Branch("hCandidate_m",&hCandidate_m,"hCandidate_m/D");
   TBranch *hCand_sd  = trimTree->Branch("hCandidate_sd",&hCandidate_sd,"hCandidate_sd/D");
   TBranch *hCand_dmdhbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagHbbvsQCD",&hCandidate_dmdhbbvqcd,"hCandidate_dmdhbbvqcd/D");
   TBranch *hCand_dmdzbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagZbbvsQCD",&hCandidate_dmdzbbvqcd,"hCandidate_dmdzbbvqcd/D");
   TBranch *hCand_dmdzhbbvqcd = trimTree->Branch("hCandidate_DeepMassDecorrelTagZHbbvsQCD",&hCandidate_dmdzhbbvqcd,"hCandidate_dmdzhbbvqcd/D");
   TBranch *hCand_middb       = trimTree->Branch("hCandidate_pfMassIndependentDeepDoubleBvLJetTagsProbHbb",&hCandidate_middb,"hCandidate_middb/D");
   TBranch *ZCand     = trimTree->Branch("ZCandidate","TLorentzVector",&ZCandidate);
   TBranch *ZCand_pt  = trimTree->Branch("ZCandidate_pt",&ZCandidate_pt,"ZCandidate_pt/D");
   TBranch *ZCand_phi = trimTree->Branch("ZCandidate_phi",&ZCandidate_phi,"ZCandidate_phi/D");
   TBranch *ZCand_eta = trimTree->Branch("ZCandidate_eta",&ZCandidate_eta,"ZCandidate_eta/D");
   TBranch *ZCand_m   = trimTree->Branch("ZCandidate_m",&ZCandidate_m,"ZCandidate_m/D");
   TBranch *ZpMest    = trimTree->Branch("ZPrime_mass_est",&mEstZp,"mEstZp/D");
   TBranch *NDMest    = trimTree->Branch("ND_mass_est",&mEstND,"mEstND/D");
   TBranch *NSMest    = trimTree->Branch("NS_mass_est",&mEstNS,"mEstNS/D");
   TBranch *evntweight = trimTree->Branch("event_weight",&evntw,"evntw/F");
   hnskimed->SetBinContent(1,nentries);
   hnorigevnts->SetBinContent(1,totalOriginalEvents);

   float zmwinlow = 70.;
   float zmwinhi  = 110.;
   float hptcut   = 250.;

   //Recursive Jigsaw Part
   LabRecoFrame         LABcontra("LABcontra","LABcontra");
   DecayRecoFrame       Zp("Zp","Z'");
   DecayRecoFrame       ND("ND","N_{D}");
   DecayRecoFrame       NDbar("NDbar","N_{Dbar}");
   VisibleRecoFrame     Z("Z","Z");
   InvisibleRecoFrame   NS("NS","N_{S}");
   VisibleRecoFrame     h("h","h");
   InvisibleRecoFrame   NSbar("NSbar","Z_{Sbar}");

   LABcontra.SetChildFrame(Zp);
   Zp.AddChildFrame(ND);
   Zp.AddChildFrame(NDbar);
   ND.AddChildFrame(Z);
   ND.AddChildFrame(NS);
   NDbar.AddChildFrame(h);
   NDbar.AddChildFrame(NSbar);

   LABcontra.InitializeTree();
   
   // Invisible Group
   InvisibleGroup INVcontra("INVcontra","NS NS Jigsaws");
   INVcontra.AddFrame(NS);
   INVcontra.AddFrame(NSbar);
   
   // Set NS NS~ mass equal to Z h mass
   SetMassInvJigsaw NSNSM("NSNSM", "M_{NSNS} = m_{Zh}");
   INVcontra.AddJigsaw(NSNSM);
   
   SetRapidityInvJigsaw NSNSR("NSNSR", "#eta_{NSNS} = #eta_{ZH}");
   INVcontra.AddJigsaw(NSNSR);
   NSNSR.AddVisibleFrames(LABcontra.GetListVisibleFrames());
   
   //MinMassesSqInvJigsaw MinMND("MinMND","min M_{D}, M_{ND}= M_{NDbar}",2);
   ContraBoostInvJigsaw MinMND("MinMND","min M_{ND}, M_{ND}= M_{NDbar}");
   INVcontra.AddJigsaw(MinMND);
   MinMND.AddVisibleFrame(Z, 0);
   MinMND.AddVisibleFrame(h, 1);
   MinMND.AddInvisibleFrame(NS, 0);
   MinMND.AddInvisibleFrame(NSbar, 1);

   LABcontra.InitializeAnalysis();

   //Trigger Stuff
   string trgtit;
   string delim  = ",";
   string ourtrg;
   int trgidx    = -1;
   int trgval;
   TFile * fthen;
   TFile * fnow;


   //counters
   int counttrigpass = 0;
   int countzpass    = 0;
   int counthpass    = 0;
   int countmetpass  = 0;
   int countpass     = 0;



   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Define some boos to signal events to write
      bool passZ    = false;
      bool passh    = false;
      bool passMET  = false;
      bool passTrig = false;
      bool passFil  = false;
          
      //A counter, for my sanity
      if (jentry%25000 == 0) {
      	std::cout<<"    analyzing event "<<jentry<<std::endl;
      }

      //Trigger decisions
      size_t pos = 0;
      string token;
      if (year == 18) {
	std::cout<<"Find you 2018 triggers, moron"<<std::endl;
	break;
      }
      if (year == 17) {
	ourtrg = "HLT_Mu50_v";
      }
      if (year == 16) {
	std::cout<<"Find you 2016 triggers, moron"<<std::endl;
	break;
      }
      if (jentry == 0) {
	trgtit = fChain->GetBranch("TriggerPass")->GetTitle();
	fthen = fChain->GetCurrentFile();
	while ((pos = trgtit.find(delim)) != std::string::npos && token != ourtrg) {
	  token = trgtit.substr(0,pos);
	  trgidx += 1;
	  trgtit.erase(0,pos+delim.length());
	}
      }
      else {
	fnow = fChain->GetCurrentFile();
	if (fnow != fthen) {
	  std::cout<< "New file in TChain" <<std::endl;
	  trgtit = fChain->GetBranch("TriggerPass")->GetTitle();
	  trgidx = -1;
	  fthen = fnow;
	  while ((pos = trgtit.find(delim)) != std::string::npos && token != ourtrg) {
	    token = trgtit.substr(0,pos);
	    trgidx += 1;
	    trgtit.erase(0,pos+delim.length());
	  }
	}
      }
      trgval = TriggerPass->at(trgidx);
      if (trgval == 1) {
	passTrig = true;
	counttrigpass += 1;
      }


      //eeBadScFilter
      if (sampleType == 0) {
	if (fChain->GetLeaf("eeBadScFilter")->GetValue() == 1.0) {
	  passFil = true;
	}
      }
	

      //Z Candidate Build
      unsigned int nZs = ZCandidates->size();
      TLorentzVector theZ;
      double baseZdiff = 99999;
      if (nZs > 0) {
	//std::cout<<"Greater than 0 Zs"<<std::endl;
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidates->begin(); zit != ZCandidates->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  //if ((massZdiff < baseZdiff) && (zit->M() >= zmwinlow) && (zit->M() <= zmwinhi)) {
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    //std::cout<<"    massdiff: "<<massZdiff<<std::endl;
	    //std::cout<<"        mass: "<<zit->M()<<std::endl;
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	  }
	}
      }

      //Higgs Candidate Build
      unsigned long nfat = JetsAK8Clean->size();
      TLorentzVector theh;
      theh.SetPtEtaPhiM(0.0,0.0,0.0,0.0);;
      TLorentzVector fat;
      double basehdiff = 99999;
      double basept = 0;
      double hsd = 0;
      double fsd = 0;
      double hdmdhbbvqcd = 0;
      double hdmdzbbvqcd = 0;
      double hdmdzhbbvqcd = 0;
      double hmiddb = 0;
      bool   fid = 0;

      if (nfat > 0) {
	for (unsigned long i =0; i < nfat; ++i) {
	  fat = JetsAK8Clean->at(i);
	  fsd = JetsAK8Clean_softDropMass->at(i);
	  fid = JetsAK8Clean_ID->at(i);
	  double masshdiff = std::abs(125.18 - fsd);
	  if ((masshdiff < basehdiff) && (fat.Pt() > hptcut) && fid && std::abs(fat.Eta()) < 2.4 && (fsd > 10)) {
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
      LABcontra.ClearEvent();
      INVcontra.SetLabFrameThreeVector(met3);
      Z.SetLabFrameFourVector(theZ);
      h.SetLabFrameFourVector(theh);
      LABcontra.AnalyzeEvent();
      mEstZp = Zp.GetMass();
      mEstND = ND.GetMass();
      mEstNS = NS.GetMass();

      if (passZ && passTrig) {
	countzpass +=1 ;
      }

      if (passh && passZ && passTrig) {
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
	ZCandidate = theZ;
	ZCandidate_pt  = theZ.Pt();
	ZCandidate_phi = theZ.Phi();
	ZCandidate_eta = theZ.Eta();
	ZCandidate_m   = theZ.M();
	counthpass += 1;
      }
      
      //Fill the Tree
      evntw = 1.;
      if (Cut(ientry) < 0) continue;
      if (passZ && passh && passTrig && sampleType !=0) {
	trimTree->Fill();
	countpass += 1;
	}
      if (passZ && passh && passTrig && sampleType == 0 && passFil) {
	trimTree->Fill();
	countpass += 1;
      }

   }

   std::cout<<"Passing Trigger req: "<<counttrigpass<<std::endl;
   std::cout<<"Passing Z  req:      "<<countzpass<<std::endl;
   std::cout<<"Passing h  req:      "<<counthpass<<std::endl;
   std::cout<<"Passing    req:      "<<countpass<<std::endl;


   htrigpass->SetBinContent(1,counttrigpass);
   hZpass->SetBinContent(1,countzpass);
   hHpass->SetBinContent(1,counthpass);
   hpass->SetBinContent(1,countpass);
   
   trimFile->Write();
   trimFile->Close();   
   //std::cout<<"trimmed to "<< passEvents <<" events"<<std::endl;
   std::cout<<"Completed your topiary garden, hopefully your tastes have not changed."<<std::endl;

   
}

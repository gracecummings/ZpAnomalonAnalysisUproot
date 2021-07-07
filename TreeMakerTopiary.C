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
   fChain->SetBranchStatus("ZCandidates*",1);
   fChain->SetBranchStatus("SelectedElectrons*",1);
   fChain->SetBranchStatus("eeBadScFilter",1);
   fChain->SetBranchStatus("ZCandidatesMuMu",1);
   fChain->SetBranchStatus("ZCandidatesEE",1);
   fChain->SetBranchStatus("ZCandidatesEU",1);

   if (sampleType != 0){
     fChain->SetBranchStatus("GenParticles*",1);
   }


   TFile qcdnnloFile("../DYCorrection/lindert_qcd_nnlo_sf.root","READ");
   TH1D *hqcdnnlosf  = (TH1D*)qcdnnloFile.Get("eej");
   hqcdnnlosf->SetDirectory(0);
   qcdnnloFile.Close();
   TFile ewknloFile("../DYCorrection/merged_kfactors_zjets.root","READ");
   TH1F *hewknlosf = (TH1F*)ewknloFile.Get("kfactor_monojet_ewk");
   hewknlosf->SetDirectory(0);
   ewknloFile.Close();


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
   double  evntw;
   TLorentzVector LMuCandidate;
   double LMuCandidate_pt;
   double LMuCandidate_phi;
   double LMuCandidate_eta;
   double LMuCandidate_m;
   TLorentzVector sLMuCandidate;
   double sLMuCandidate_pt;
   double sLMuCandidate_phi;
   double sLMuCandidate_eta;
   double sLMuCandidate_m;
   TLorentzVector LEleCandidate;
   double LEleCandidate_pt;
   double LEleCandidate_phi;
   double LEleCandidate_eta;
   double LEleCandidate_m;
   TLorentzVector sLEleCandidate;
   double sLEleCandidate_pt;
   double sLEleCandidate_phi;
   double sLEleCandidate_eta;
   double sLEleCandidate_m;
   double ghCandidate_pt;
   double ghCandidate_phi;
   double ghCandidate_eta;
   double ghCandidate_m;
   double gzCandidate_pt;
   double gzCandidate_phi;
   double gzCandidate_eta;
   double gzCandidate_m;
   double channelflag;

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
   TBranch *evntweight = trimTree->Branch("event_weight",&evntw,"evntw/D");
   TBranch *channelf   = trimTree->Branch("channel_flag",&channelflag,"channelflag/D");
   TBranch *LMuCand     = trimTree->Branch("LMuCandidate","TLorentzVector",&LMuCandidate);
   TBranch *LMuCand_pt  = trimTree->Branch("LMuCandidate_pt",&LMuCandidate_pt,"LMuCandidate_pt/D");
   TBranch *LMuCand_phi = trimTree->Branch("LMuCandidate_phi",&LMuCandidate_phi,"LMuCandidate_phi/D");
   TBranch *LMuCand_eta = trimTree->Branch("LMuCandidate_eta",&LMuCandidate_eta,"LMuCandidate_eta/D");
   TBranch *sLMuCand     = trimTree->Branch("sLMuCandidate","TLorentzVector",&sLMuCandidate);
   TBranch *sLMuCand_pt  = trimTree->Branch("sLMuCandidate_pt",&sLMuCandidate_pt,"sLMuCandidate_pt/D");
   TBranch *sLMuCand_phi = trimTree->Branch("sLMuCandidate_phi",&sLMuCandidate_phi,"sLMuCandidate_phi/D");
   TBranch *sLMuCand_eta = trimTree->Branch("sLMuCandidate_eta",&sLMuCandidate_eta,"sLMuCandidate_eta/D");
   TBranch *LEleCand     = trimTree->Branch("LEleCandidate","TLorentzVector",&LEleCandidate);
   TBranch *LEleCand_pt  = trimTree->Branch("LEleCandidate_pt",&LEleCandidate_pt,"LEleCandidate_pt/D");
   TBranch *LEleCand_phi = trimTree->Branch("LEleCandidate_phi",&LEleCandidate_phi,"LEleCandidate_phi/D");
   TBranch *LEleCand_eta = trimTree->Branch("LEleCandidate_eta",&LEleCandidate_eta,"LEleCandidate_eta/D");
   TBranch *sLEleCand     = trimTree->Branch("sLEleCandidate","TLorentzVector",&sLEleCandidate);
   TBranch *sLEleCand_pt  = trimTree->Branch("sLEleCandidate_pt",&sLEleCandidate_pt,"sLEleCandidate_pt/D");
   TBranch *sLEleCand_phi = trimTree->Branch("sLEleCandidate_phi",&sLEleCandidate_phi,"sLEleCandidate_phi/D");
   TBranch *sLEleCand_eta = trimTree->Branch("sLEleCandidate_eta",&sLEleCandidate_eta,"sLEleCandidate_eta/D");
   TBranch *ghCand_pt  = trimTree->Branch("ghCandidate_pt",&ghCandidate_pt,"ghCandidate_pt/D");
   TBranch *ghCand_phi = trimTree->Branch("ghCandidate_phi",&ghCandidate_phi,"ghCandidate_phi/D");
   TBranch *ghCand_eta = trimTree->Branch("ghCandidate_eta",&ghCandidate_eta,"ghCandidate_eta/D");
   TBranch *ghCand_m   = trimTree->Branch("ghCandidate_m",&ghCandidate_m,"ghCandidate_m/D");
   TBranch *gzCand_pt  = trimTree->Branch("gzCandidate_pt",&gzCandidate_pt,"gzCandidate_pt/D");
   TBranch *gzCand_phi = trimTree->Branch("gzCandidate_phi",&gzCandidate_phi,"gzCandidate_phi/D");
   TBranch *gzCand_eta = trimTree->Branch("gzCandidate_eta",&gzCandidate_eta,"gzCandidate_eta/D");
   TBranch *gzCand_m   = trimTree->Branch("gzCandidate_m",&gzCandidate_m,"gzCandidate_m/D");

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
   TFile * fthen = 0;
   TFile * fnow = 0;

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
      bool mumuchan = false;
      double channel = -1.0;
          
      //A counter, for my sanity
      if (jentry%25000 == 0) {
      	std::cout<<"    analyzing event "<<jentry<<std::endl;
      }

      //debug
      //if (jentry == 200) {
      //break;
      //}

      //Trigger decisions
      size_t pos = 0;
      string token;
      if (year == 18) {
	ourtrg = "HLT_Mu55_v";
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

      //DY+Jets k-factors+GenParticleInfo
      float evntw_hold = 1.;
      TLorentzVector theGenZ;
      TLorentzVector theGenH;
      if (sampleType != 0) {//Not Data
	int gpid;
	unsigned long ngen = GenParticles->size();
	for (unsigned long i = 0; i < ngen; ++i) {
	  int gpid = GenParticles_PdgId->at(i);
	  if (gpid == 23) {
	    theGenZ = GenParticles->at(i);
	  }
	  if (gpid == 25) {
	    theGenH = GenParticles->at(i);
	  }
	}

	//Calculate the k-factor
	float qcdnlosf   = 1;
	double qcdnnlosf = 1;
	float ewknlosf   = 1;
	if (sampleType == 2) {//DY+Jets
	  qcdnlosf = 1.423*exp(-0.002257*theGenZ.Pt())+0.451;
	  if (qcdnlosf <= 0.0) {
	    qcdnlosf = 1.;
	  }
	  int zptbinqcd  = hqcdnnlosf->FindBin(theGenZ.Pt());
	  qcdnnlosf = hqcdnnlosf->GetBinContent(zptbinqcd);
	  if (qcdnnlosf <= 0.0) {
	    qcdnnlosf = 1.;
	  }
	  int zptbinewk = hewknlosf->FindBin(theGenZ.Pt());
	  ewknlosf = hewknlosf->GetBinContent(zptbinewk);
	  if (ewknlosf <= 0.0) {
	    ewknlosf = 1.;
	  }
	  evntw_hold = ewknlosf*qcdnnlosf*qcdnlosf;
	}
      }
      
      //Z exploration
      unsigned int nselmu = SelectedMuons->size();
      unsigned int nselel = SelectedElectrons->size();
      unsigned int nZmumu = ZCandidatesMuMu->size();
      unsigned int nZee = ZCandidatesEE->size();
      unsigned int nZeu = ZCandidatesEU->size();
      TLorentzVector leadmu;
      TLorentzVector subleadmu;
      double muptmax = 0;
      unsigned int nZs = ZCandidates->size();
      TLorentzVector theZ;
      double baseZdiff = 99999;
      //Channel Flags
      ///*
      if (nZmumu > 0 && nZee == 0 && nZeu == 0){
	//in binary 100, 4 in decimal
	channel = 4.;//4 in decimal
	mumuchan = true;
	std::vector<TLorentzVector>::iterator muit;
	for (muit = SelectedMuons->begin(); muit != SelectedMuons->end();++muit) { //might need to move after Z
	  if (muit->Pt() > muptmax) {
	    muptmax = muit->Pt();
	    leadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());
	  }
	  else {
	    subleadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());;
	 }
	}
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidatesMuMu->begin(); zit != ZCandidatesMuMu->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	  }
	}
      }
      if (nZmumu > 0 && nZee > 0 && nZeu == 0) {
	//110 in binary, 6 in decimal
	channel = 6.;//
      }
      if (nZmumu > 0 && nZee == 0 && nZeu > 0) {
	//101 in binary, 5 in decimal
	channel = 5.;//
      }
      if (nZmumu > 0 && nZee > 0 && nZeu > 0) {
	//111 in binary, 7 in decimal
	channel = 7.;
      }
      if (nZmumu == 0 && nZee > 0 && nZeu == 0) {
	//010
	channel = 2.;
      }
      if (nZmumu == 0 && nZee > 0 && nZeu > 0) {
	//011
	channel = 3.;
      }
      if (nZmumu == 0 && nZee == 0 && nZeu > 0) {
	//001
	channel = 1.;
      }
      //*/

      //Z Candidate Build
      //For old ntuples
      /*
      if (nselmu > 0 && nselel == 0) {
      	mumuchan = true;
	std::vector<TLorentzVector>::iterator muit;
	for (muit = SelectedMuons->begin(); muit != SelectedMuons->end();++muit) { //might need to move after Z
	  if (muit->Pt() > muptmax) {
	    muptmax = muit->Pt();
	    leadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());
	  }
	  else {
	    subleadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());;
	 }
	}
      }
      if (nZs > 0) {
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidates->begin(); zit != ZCandidates->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	  }
	}
      }
      */
      
      //Higgs Candidate Build
      //unsigned long nfat = JetsAK8>size();
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

      //reclustered jets
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

      //unreclustered jets
      /*
      if (nfat > 0) {
      for (unsigned long i =0; i < nfat; ++i) {
      	fat = JetsAK8->at(i);
        fsd = JetsAK8_softDropMass->at(i);
        fid = JetsAK8_ID->at(i);
        double masshdiff = std::abs(125.18 - fsd);
        if ((masshdiff < basehdiff) && (fat.Pt() > hptcut) && fid && std::abs(fat.Eta()) < 2.4 && (fsd > 10)) {
	  basehdiff = masshdiff;
          theh = fat;
          hsd = fsd;
          hdmdhbbvqcd  = JetsAK8_DeepMassDecorrelTagHbbvsQCD->at(i);
          hdmdzbbvqcd  = JetsAK8_DeepMassDecorrelTagZbbvsQCD->at(i);
          hdmdzhbbvqcd = JetsAK8_DeepMassDecorrelTagZHbbvsQCD->at(i);
          hmiddb = JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i);
          passh = true;
        }
      }
      }
      */

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

      if (passZ && passTrig && mumuchan) {
	countzpass +=1 ;
      }

      if (passh && passZ && passTrig && mumuchan) {
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
	evntw          = evntw_hold;
	LMuCandidate = leadmu;
	LMuCandidate_pt  = leadmu.Pt();
	LMuCandidate_phi = leadmu.Phi();
	LMuCandidate_eta = leadmu.Eta();
	LMuCandidate_m   = leadmu.M();
	sLMuCandidate = subleadmu;
	sLMuCandidate_pt  = subleadmu.Pt();
	sLMuCandidate_phi = subleadmu.Phi();
	sLMuCandidate_eta = subleadmu.Eta();
	sLMuCandidate_m   = subleadmu.M();
	ghCandidate_pt  = theGenH.Pt();
	ghCandidate_phi = theGenH.Phi();
	ghCandidate_eta = theGenH.Eta();
	ghCandidate_m   = theGenH.M();
	gzCandidate_pt  = theGenZ.Pt();
	gzCandidate_phi = theGenZ.Phi();
	gzCandidate_eta = theGenZ.Eta();
	gzCandidate_m   = theGenZ.M();
	channelflag = channel;
	counthpass += 1;
      }
      
      //Fill the Tree
      if (Cut(ientry) < 0) continue;
      if (passZ && passh && passTrig && sampleType !=0 && mumuchan) {
	trimTree->Fill();
	countpass += 1;
	}
      if (passZ && passh && passTrig && sampleType == 0 && passFil && mumuchan) {
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

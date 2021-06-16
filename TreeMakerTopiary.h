//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  5 12:43:52 2021 by ROOT version 6.14/04
// from TTree PreSelection/PreSelection
// found on file: ../dataHandling/2017/Run2017F-31Mar2018-v1.SingleMuon_0_chunk0_oct08_skim_2640267MiniAODEvents.root
//////////////////////////////////////////////////////////

#ifndef TreeMakerTopiary_h
#define TreeMakerTopiary_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include "vector"

class TreeMakerTopiary {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          RunNum;
   UInt_t          LumiBlockNum;
   ULong64_t       EvtNum;
   Bool_t          BadChargedCandidateFilter;
   Bool_t          BadPFMuonFilter;
   Int_t           BTags;
   Int_t           BTagsclean;
   Int_t           BTagsDeepCSV;
   Int_t           BTagsDeepCSVclean;
   Double_t        CaloMET;
   Double_t        CaloMETPhi;
   Int_t           CSCTightHaloFilter;
   Double_t        DeltaPhi1;
   Double_t        DeltaPhi1clean;
   Double_t        DeltaPhi2;
   Double_t        DeltaPhi2clean;
   Double_t        DeltaPhi3;
   Double_t        DeltaPhi3clean;
   Double_t        DeltaPhi4;
   Double_t        DeltaPhi4clean;
   Int_t           ecalBadCalibFilter;
   Bool_t          ecalBadCalibReducedExtraFilter;
   Bool_t          ecalBadCalibReducedFilter;
   Int_t           EcalDeadCellBoundaryEnergyFilter;
   Int_t           EcalDeadCellTriggerPrimitiveFilter;
   Int_t           eeBadScFilter;
   vector<TLorentzVector> *Electrons;
   vector<int>     *Electrons_charge;
   vector<bool>    *Electrons_mediumID;
   vector<double>  *Electrons_MTW;
   vector<bool>    *Electrons_passIso;
   vector<bool>    *Electrons_tightID;
   Double_t        fixedGridRhoFastjetAll;
   Int_t           globalSuperTightHalo2016Filter;
   Int_t           globalTightHalo2016Filter;
   Bool_t          hasGenPromptPhoton;
   Int_t           HBHEIsoNoiseFilter;
   Int_t           HBHENoiseFilter;
   vector<TLorentzVector> *HLTMuonObjects;
   Double_t        HT;
   Double_t        HT5;
   Double_t        HT5clean;
   Double_t        HTclean;
   Double_t        HTOnline;
   Int_t           isoElectronTracks;
   Int_t           isoElectronTracksclean;
   Int_t           isoMuonTracks;
   Int_t           isoMuonTracksclean;
   Int_t           isoPionTracks;
   Int_t           isoPionTracksclean;
   Bool_t          JetID;
   Bool_t          JetIDAK8;
   Bool_t          JetIDAK8Clean;
   Bool_t          JetIDclean;
   vector<TLorentzVector> *Jets;
   vector<double>  *Jets_axismajor;
   vector<double>  *Jets_axisminor;
   vector<double>  *Jets_bDiscriminatorCSV;
   vector<double>  *Jets_bJetTagDeepCSVBvsAll;
   vector<double>  *Jets_chargedEmEnergyFraction;
   vector<double>  *Jets_chargedHadronEnergyFraction;
   vector<int>     *Jets_chargedHadronMultiplicity;
   vector<int>     *Jets_chargedMultiplicity;
   vector<double>  *Jets_electronEnergyFraction;
   vector<int>     *Jets_electronMultiplicity;
   vector<int>     *Jets_hadronFlavor;
   vector<double>  *Jets_hfEMEnergyFraction;
   vector<double>  *Jets_hfHadronEnergyFraction;
   vector<bool>    *Jets_HTMask;
   vector<bool>    *Jets_ID;
   vector<double>  *Jets_jecFactor;
   vector<double>  *Jets_jecUnc;
   vector<bool>    *Jets_LeptonMask;
   vector<bool>    *Jets_MHTMask;
   vector<int>     *Jets_multiplicity;
   vector<double>  *Jets_muonEnergyFraction;
   vector<int>     *Jets_muonMultiplicity;
   vector<double>  *Jets_neutralEmEnergyFraction;
   vector<double>  *Jets_neutralHadronEnergyFraction;
   vector<int>     *Jets_neutralHadronMultiplicity;
   vector<int>     *Jets_neutralMultiplicity;
   vector<int>     *Jets_partonFlavor;
   vector<double>  *Jets_photonEnergyFraction;
   vector<int>     *Jets_photonMultiplicity;
   vector<double>  *Jets_ptD;
   vector<double>  *Jets_qgLikelihood;
   vector<TLorentzVector> *JetsAK8;
   vector<double>  *JetsAK8_axismajor;
   vector<double>  *JetsAK8_axisminor;
   vector<double>  *JetsAK8_DeepMassDecorrelTagbbvsLight;
   vector<double>  *JetsAK8_DeepMassDecorrelTagHbbvsQCD;
   vector<double>  *JetsAK8_DeepMassDecorrelTagTvsQCD;
   vector<double>  *JetsAK8_DeepMassDecorrelTagWvsQCD;
   vector<double>  *JetsAK8_DeepMassDecorrelTagZbbvsQCD;
   vector<double>  *JetsAK8_DeepMassDecorrelTagZHbbvsQCD;
   vector<double>  *JetsAK8_DeepMassDecorrelTagZvsQCD;
   vector<double>  *JetsAK8_DeepTagHbbvsQCD;
   vector<double>  *JetsAK8_DeepTagTvsQCD;
   vector<double>  *JetsAK8_DeepTagWvsQCD;
   vector<double>  *JetsAK8_DeepTagZbbvsQCD;
   vector<double>  *JetsAK8_DeepTagZvsQCD;
   vector<double>  *JetsAK8_doubleBDiscriminator;
   vector<double>  *JetsAK8_ecfN2b1;
   vector<double>  *JetsAK8_ecfN2b2;
   vector<double>  *JetsAK8_ecfN3b1;
   vector<double>  *JetsAK8_ecfN3b2;
   vector<double>  *JetsAK8_girth;
   vector<bool>    *JetsAK8_ID;
   vector<double>  *JetsAK8_jecFactor;
   vector<double>  *JetsAK8_jecUnc;
   vector<int>     *JetsAK8_multiplicity;
   vector<double>  *JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   vector<double>  *JetsAK8_ptD;
   vector<double>  *JetsAK8_softDropMass;
   vector<vector<TLorentzVector> > *JetsAK8_subjets;
   vector<TLorentzVector> *JetsAK8Clean;
   vector<double>  *JetsAK8Clean_axismajor;
   vector<double>  *JetsAK8Clean_axisminor;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagbbvsLight;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagTvsQCD;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagWvsQCD;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD;
   vector<double>  *JetsAK8Clean_DeepMassDecorrelTagZvsQCD;
   vector<double>  *JetsAK8Clean_DeepTagHbbvsQCD;
   vector<double>  *JetsAK8Clean_DeepTagTvsQCD;
   vector<double>  *JetsAK8Clean_DeepTagWvsQCD;
   vector<double>  *JetsAK8Clean_DeepTagZbbvsQCD;
   vector<double>  *JetsAK8Clean_DeepTagZvsQCD;
   vector<double>  *JetsAK8Clean_doubleBDiscriminator;
   vector<double>  *JetsAK8Clean_girth;
   vector<bool>    *JetsAK8Clean_ID;
   vector<double>  *JetsAK8Clean_jecFactor;
   vector<int>     *JetsAK8Clean_multiplicity;
   vector<double>  *JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   vector<double>  *JetsAK8Clean_ptD;
   vector<double>  *JetsAK8Clean_softDropMass;
   vector<vector<TLorentzVector> > *JetsAK8Clean_subjets;
   vector<TLorentzVector> *Jetsclean;
   vector<double>  *Jetsclean_bDiscriminatorCSV;
   vector<double>  *Jetsclean_bJetTagDeepCSVBvsAll;
   vector<int>     *Jetsclean_hadronFlavor;
   vector<bool>    *Jetsclean_HTMask;
   vector<bool>    *Jetsclean_ID;
   vector<bool>    *Jetsclean_MHTMask;
   vector<int>     *Jetsclean_partonFlavor;
   Double_t        MET;
   Double_t        METclean;
   Double_t        METcleanOrig;
   Double_t        METOrig;
   Double_t        METPhi;
   Double_t        METPhiclean;
   Double_t        METPhicleanOrig;
   Double_t        METPhiOrig;
   Double_t        METSignificance;
   Double_t        METSignificanceclean;
   Double_t        MHT;
   Double_t        MHTclean;
   Double_t        MHTOnline;
   Double_t        MHTPhi;
   Double_t        MHTPhiclean;
   vector<TLorentzVector> *Muons;
   vector<int>     *Muons_charge;
   vector<bool>    *Muons_mediumID;
   vector<double>  *Muons_MTW;
   vector<bool>    *Muons_passIso;
   vector<bool>    *Muons_tightID;
   Int_t           nAllVertices;
   Int_t           NElectrons;
   Int_t           NJets;
   Int_t           NJetsclean;
   Int_t           NMuons;
   Double_t        NonPrefiringProb;
   Double_t        NonPrefiringProbDown;
   Double_t        NonPrefiringProbUp;
   Int_t           NVtx;
   Double_t        PFCaloMETRatio;
   vector<TLorentzVector> *Photons;
   vector<bool>    *Photons_electronFakes;
   vector<bool>    *Photons_fullID;
   vector<double>  *Photons_genMatched;
   vector<double>  *Photons_hadTowOverEM;
   vector<bool>    *Photons_hasPixelSeed;
   vector<double>  *Photons_isEB;
   vector<bool>    *Photons_nonPrompt;
   vector<double>  *Photons_passElectronVeto;
   vector<double>  *Photons_pfChargedIso;
   vector<double>  *Photons_pfChargedIsoRhoCorr;
   vector<double>  *Photons_pfGammaIso;
   vector<double>  *Photons_pfGammaIsoRhoCorr;
   vector<double>  *Photons_pfNeutralIso;
   vector<double>  *Photons_pfNeutralIsoRhoCorr;
   vector<double>  *Photons_sigmaIetaIeta;
   Double_t        PrescaleWeightHT;
   Int_t           PrimaryVertexFilter;
   vector<TLorentzVector> *SelectedElectrons;
   vector<TLorentzVector> *SelectedMuons;
   vector<TLorentzVector> *TAPElectronTracks;
   vector<double>  *TAPElectronTracks_dxypv;
   vector<bool>    *TAPElectronTracks_leptonMatch;
   vector<double>  *TAPElectronTracks_mT;
   vector<double>  *TAPElectronTracks_pfRelIso03chg;
   vector<double>  *TAPElectronTracks_trkiso;
   vector<TLorentzVector> *TAPMuonTracks;
   vector<double>  *TAPMuonTracks_dxypv;
   vector<bool>    *TAPMuonTracks_leptonMatch;
   vector<double>  *TAPMuonTracks_mT;
   vector<double>  *TAPMuonTracks_pfRelIso03chg;
   vector<double>  *TAPMuonTracks_trkiso;
   vector<TLorentzVector> *TAPPionTracks;
   vector<double>  *TAPPionTracks_dxypv;
   vector<bool>    *TAPPionTracks_leptonMatch;
   vector<double>  *TAPPionTracks_mT;
   vector<double>  *TAPPionTracks_pfRelIso03chg;
   vector<double>  *TAPPionTracks_trkiso;
   vector<int>     *TriggerPass;
   vector<int>     *TriggerPrescales;
   vector<int>     *TriggerVersion;
   vector<TLorentzVector> *ZCandidates;
   vector<TLorentzVector> *ZCandidatesMuMu;
   vector<TLorentzVector> *ZCandidatesEE;
   vector<TLorentzVector> *ZCandidatesEU;

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_BadChargedCandidateFilter;   //!
   TBranch        *b_BadPFMuonFilter;   //!
   TBranch        *b_BTags;   //!
   TBranch        *b_BTagsclean;   //!
   TBranch        *b_BTagsDeepCSV;   //!
   TBranch        *b_BTagsDeepCSVclean;   //!
   TBranch        *b_CaloMET;   //!
   TBranch        *b_CaloMETPhi;   //!
   TBranch        *b_CSCTightHaloFilter;   //!
   TBranch        *b_DeltaPhi1;   //!
   TBranch        *b_DeltaPhi1clean;   //!
   TBranch        *b_DeltaPhi2;   //!
   TBranch        *b_DeltaPhi2clean;   //!
   TBranch        *b_DeltaPhi3;   //!
   TBranch        *b_DeltaPhi3clean;   //!
   TBranch        *b_DeltaPhi4;   //!
   TBranch        *b_DeltaPhi4clean;   //!
   TBranch        *b_ecalBadCalibFilter;   //!
   TBranch        *b_ecalBadCalibReducedExtraFilter;   //!
   TBranch        *b_ecalBadCalibReducedFilter;   //!
   TBranch        *b_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_eeBadScFilter;   //!
   TBranch        *b_Electrons;   //!
   TBranch        *b_Electrons_charge;   //!
   TBranch        *b_Electrons_mediumID;   //!
   TBranch        *b_Electrons_MTW;   //!
   TBranch        *b_Electrons_passIso;   //!
   TBranch        *b_Electrons_tightID;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_globalTightHalo2016Filter;   //!
   TBranch        *b_hasGenPromptPhoton;   //!
   TBranch        *b_HBHEIsoNoiseFilter;   //!
   TBranch        *b_HBHENoiseFilter;   //!
   TBranch        *b_HLTMuonObjects;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_HT5;   //!
   TBranch        *b_HT5clean;   //!
   TBranch        *b_HTclean;   //!
   TBranch        *b_HTOnline;   //!
   TBranch        *b_isoElectronTracks;   //!
   TBranch        *b_isoElectronTracksclean;   //!
   TBranch        *b_isoMuonTracks;   //!
   TBranch        *b_isoMuonTracksclean;   //!
   TBranch        *b_isoPionTracks;   //!
   TBranch        *b_isoPionTracksclean;   //!
   TBranch        *b_JetID;   //!
   TBranch        *b_JetIDAK8;   //!
   TBranch        *b_JetIDAK8Clean;   //!
   TBranch        *b_JetIDclean;   //!
   TBranch        *b_Jets;   //!
   TBranch        *b_Jets_axismajor;   //!
   TBranch        *b_Jets_axisminor;   //!
   TBranch        *b_Jets_bDiscriminatorCSV;   //!
   TBranch        *b_Jets_bJetTagDeepCSVBvsAll;   //!
   TBranch        *b_Jets_chargedEmEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronMultiplicity;   //!
   TBranch        *b_Jets_chargedMultiplicity;   //!
   TBranch        *b_Jets_electronEnergyFraction;   //!
   TBranch        *b_Jets_electronMultiplicity;   //!
   TBranch        *b_Jets_hadronFlavor;   //!
   TBranch        *b_Jets_hfEMEnergyFraction;   //!
   TBranch        *b_Jets_hfHadronEnergyFraction;   //!
   TBranch        *b_Jets_HTMask;   //!
   TBranch        *b_Jets_ID;   //!
   TBranch        *b_Jets_jecFactor;   //!
   TBranch        *b_Jets_jecUnc;   //!
   TBranch        *b_Jets_LeptonMask;   //!
   TBranch        *b_Jets_MHTMask;   //!
   TBranch        *b_Jets_multiplicity;   //!
   TBranch        *b_Jets_muonEnergyFraction;   //!
   TBranch        *b_Jets_muonMultiplicity;   //!
   TBranch        *b_Jets_neutralEmEnergyFraction;   //!
   TBranch        *b_Jets_neutralHadronEnergyFraction;   //!
   TBranch        *b_Jets_neutralHadronMultiplicity;   //!
   TBranch        *b_Jets_neutralMultiplicity;   //!
   TBranch        *b_Jets_partonFlavor;   //!
   TBranch        *b_Jets_photonEnergyFraction;   //!
   TBranch        *b_Jets_photonMultiplicity;   //!
   TBranch        *b_Jets_ptD;   //!
   TBranch        *b_Jets_qgLikelihood;   //!
   TBranch        *b_JetsAK8;   //!
   TBranch        *b_JetsAK8_axismajor;   //!
   TBranch        *b_JetsAK8_axisminor;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagbbvsLight;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagTvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagWvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagTvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagWvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagZvsQCD;   //!
   TBranch        *b_JetsAK8_doubleBDiscriminator;   //!
   TBranch        *b_JetsAK8_ecfN2b1;   //!
   TBranch        *b_JetsAK8_ecfN2b2;   //!
   TBranch        *b_JetsAK8_ecfN3b1;   //!
   TBranch        *b_JetsAK8_ecfN3b2;   //!
   TBranch        *b_JetsAK8_girth;   //!
   TBranch        *b_JetsAK8_ID;   //!
   TBranch        *b_JetsAK8_jecFactor;   //!
   TBranch        *b_JetsAK8_jecUnc;   //!
   TBranch        *b_JetsAK8_multiplicity;   //!
   TBranch        *b_JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;   //!
   TBranch        *b_JetsAK8_ptD;   //!
   TBranch        *b_JetsAK8_softDropMass;   //!
   TBranch        *b_JetsAK8_subjets;   //!
   TBranch        *b_JetsAK8Clean;   //!
   TBranch        *b_JetsAK8Clean_axismajor;   //!
   TBranch        *b_JetsAK8Clean_axisminor;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagbbvsLight;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagTvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagWvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepMassDecorrelTagZvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepTagTvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepTagWvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8Clean_DeepTagZvsQCD;   //!
   TBranch        *b_JetsAK8Clean_doubleBDiscriminator;   //!
   TBranch        *b_JetsAK8Clean_girth;   //!
   TBranch        *b_JetsAK8Clean_ID;   //!
   TBranch        *b_JetsAK8Clean_jecFactor;   //!
   TBranch        *b_JetsAK8Clean_multiplicity;   //!
   TBranch        *b_JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;   //!
   TBranch        *b_JetsAK8Clean_ptD;   //!
   TBranch        *b_JetsAK8Clean_softDropMass;   //!
   TBranch        *b_JetsAK8Clean_subjets;   //!
   TBranch        *b_Jetsclean;   //!
   TBranch        *b_Jetsclean_bDiscriminatorCSV;   //!
   TBranch        *b_Jetsclean_bJetTagDeepCSVBvsAll;   //!
   TBranch        *b_Jetsclean_hadronFlavor;   //!
   TBranch        *b_Jetsclean_HTMask;   //!
   TBranch        *b_Jetsclean_ID;   //!
   TBranch        *b_Jetsclean_MHTMask;   //!
   TBranch        *b_Jetsclean_partonFlavor;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METclean;   //!
   TBranch        *b_METcleanOrig;   //!
   TBranch        *b_METOrig;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METPhiclean;   //!
   TBranch        *b_METPhicleanOrig;   //!
   TBranch        *b_METPhiOrig;   //!
   TBranch        *b_METSignificance;   //!
   TBranch        *b_METSignificanceclean;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTclean;   //!
   TBranch        *b_MHTOnline;   //!
   TBranch        *b_MHTPhi;   //!
   TBranch        *b_MHTPhiclean;   //!
   TBranch        *b_Muons;   //!
   TBranch        *b_Muons_charge;   //!
   TBranch        *b_Muons_mediumID;   //!
   TBranch        *b_Muons_MTW;   //!
   TBranch        *b_Muons_passIso;   //!
   TBranch        *b_Muons_tightID;   //!
   TBranch        *b_nAllVertices;   //!
   TBranch        *b_NElectrons;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_NJetsclean;   //!
   TBranch        *b_NMuons;   //!
   TBranch        *b_NonPrefiringProb;   //!
   TBranch        *b_NonPrefiringProbDown;   //!
   TBranch        *b_NonPrefiringProbUp;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_PFCaloMETRatio;   //!
   TBranch        *b_Photons;   //!
   TBranch        *b_Photons_electronFakes;   //!
   TBranch        *b_Photons_fullID;   //!
   TBranch        *b_Photons_genMatched;   //!
   TBranch        *b_Photons_hadTowOverEM;   //!
   TBranch        *b_Photons_hasPixelSeed;   //!
   TBranch        *b_Photons_isEB;   //!
   TBranch        *b_Photons_nonPrompt;   //!
   TBranch        *b_Photons_passElectronVeto;   //!
   TBranch        *b_Photons_pfChargedIso;   //!
   TBranch        *b_Photons_pfChargedIsoRhoCorr;   //!
   TBranch        *b_Photons_pfGammaIso;   //!
   TBranch        *b_Photons_pfGammaIsoRhoCorr;   //!
   TBranch        *b_Photons_pfNeutralIso;   //!
   TBranch        *b_Photons_pfNeutralIsoRhoCorr;   //!
   TBranch        *b_Photons_sigmaIetaIeta;   //!
   TBranch        *b_PrescaleWeightHT;   //!
   TBranch        *b_PrimaryVertexFilter;   //!
   TBranch        *b_SelectedElectrons;   //!
   TBranch        *b_SelectedMuons;   //!
   TBranch        *b_TAPElectronTracks;   //!
   TBranch        *b_TAPElectronTracks_dxypv;   //!
   TBranch        *b_TAPElectronTracks_leptonMatch;   //!
   TBranch        *b_TAPElectronTracks_mT;   //!
   TBranch        *b_TAPElectronTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPElectronTracks_trkiso;   //!
   TBranch        *b_TAPMuonTracks;   //!
   TBranch        *b_TAPMuonTracks_dxypv;   //!
   TBranch        *b_TAPMuonTracks_leptonMatch;   //!
   TBranch        *b_TAPMuonTracks_mT;   //!
   TBranch        *b_TAPMuonTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPMuonTracks_trkiso;   //!
   TBranch        *b_TAPPionTracks;   //!
   TBranch        *b_TAPPionTracks_dxypv;   //!
   TBranch        *b_TAPPionTracks_leptonMatch;   //!
   TBranch        *b_TAPPionTracks_mT;   //!
   TBranch        *b_TAPPionTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPPionTracks_trkiso;   //!
   TBranch        *b_TriggerPass;   //!
   TBranch        *b_TriggerPrescales;   //!
   TBranch        *b_TriggerVersion;   //!
   TBranch        *b_ZCandidates;   //!
   TBranch        *b_ZCandidatesMuMu;   //!
   TBranch        *b_ZCandidatesEE;   //!
   TBranch        *b_ZCandidatesEU;   //!

   //MC Only Variables
   vector<TLorentzVector> *GenElectrons;
   Double_t        GenHT;
   vector<TLorentzVector> *GenJets;
   vector<TLorentzVector> *GenJetsAK8;
   vector<int>     *GenJetsAK8_multiplicity;
   vector<double>  *GenJetsAK8_softDropMass;
   Double_t        GenMET;
   Double_t        GenMETPhi;
   vector<TLorentzVector> *GenMuons;
   vector<TLorentzVector> *GenParticles;
   vector<int>     *GenParticles_ParentId;
   vector<int>     *GenParticles_ParentIdx;
   vector<int>     *GenParticles_PdgId;
   vector<int>     *GenParticles_Status;
   vector<TLorentzVector> *GenTaus;
   vector<bool>    *GenTaus_had;
   Double_t        HT5JECdown;
   Double_t        HT5JECup;
   Double_t        HT5JERdown;
   Double_t        HT5JERup;
   Double_t        HTJECdown;
   Double_t        HTJECup;
   Double_t        HTJERdown;
   Double_t        HTJERup;
   Bool_t          JetIDAK8JECdown;
   Bool_t          JetIDAK8JECup;
   Bool_t          JetIDAK8JERdown;
   Bool_t          JetIDAK8JERup;
   Bool_t          JetIDJECdown;
   Bool_t          JetIDJECup;
   Bool_t          JetIDJERdown;
   Bool_t          JetIDJERup;
   vector<double>  *Jets_jerFactor;
   vector<double>  *Jets_jerFactorDown;
   vector<double>  *Jets_jerFactorUp;
   vector<int>     *Jets_origIndex;
   vector<double>  *JetsAK8_jerFactor;
   vector<double>  *JetsAK8_jerFactorDown;
   vector<double>  *JetsAK8_jerFactorUp;
   vector<int>     *JetsAK8Clean_origIndex;
   vector<double>  *JetsAK8JECup_jerFactor;
   vector<int>     *JetsAK8JECup_origIndex;
   vector<int>     *JetsAK8JERdown_origIndex;
   vector<int>     *JetsAK8JERup_origIndex;
   vector<double>  *JetsJECdown_jerFactor;
   vector<int>     *JetsJECdown_origIndex;
   vector<double>  *JetsJECup_jerFactor;
   vector<int>     *JetsJECup_origIndex;
   vector<int>     *JetsJERdown_origIndex;
   vector<int>     *JetsJERup_origIndex;
   vector<double>  *METDown;
   vector<double>  *METOrigDown;
   vector<double>  *METOrigUp;
   vector<double>  *METPhiDown;
   vector<double>  *METPhiOrigDown;
   vector<double>  *METPhiOrigUp;
   vector<double>  *METPhiUp;
   vector<double>  *METUp;
   vector<float>   *PDFweights;
   vector<float>   *PSweights;
   vector<float>   *ScaleWeights;

   //Definte the Branches
   TBranch        *b_GenElectrons;   //!
   TBranch        *b_GenHT;   //!
   TBranch        *b_GenJets;   //!
   TBranch        *b_GenJetsAK8;   //!
   TBranch        *b_GenJetsAK8_multiplicity;   //!
   TBranch        *b_GenJetsAK8_softDropMass;   
   TBranch        *b_GenMET;  
   TBranch        *b_GenMETPhi;   
   TBranch        *b_GenMuons;   //!
   TBranch        *b_GenParticles;   //!
   TBranch        *b_GenParticles_ParentId;   //!
   TBranch        *b_GenParticles_ParentIdx;   //!
   TBranch        *b_GenParticles_PdgId;   //!
   TBranch        *b_GenParticles_Status;   //!
   TBranch        *b_GenTaus;   //!
   TBranch        *b_GenTaus_had;   //!
   TBranch        *b_HT5JECdown;   //!
   TBranch        *b_HT5JECup;   //!
   TBranch        *b_HT5JERdown;   //!
   TBranch        *b_HT5JERup;   //!
   TBranch        *b_HTJECdown;   //!
   TBranch        *b_HTJECup;   //!
   TBranch        *b_HTJERdown;   //!
   TBranch        *b_HTJERup;   //!
   TBranch        *b_JetIDAK8JECdown;   //!
   TBranch        *b_JetIDAK8JECup;   //!
   TBranch        *b_JetIDAK8JERdown;   //!
   TBranch        *b_JetIDAK8JERup;   //!
   TBranch        *b_JetIDJECdown;   //!
   TBranch        *b_JetIDJECup;   //!
   TBranch        *b_JetIDJERdown;   //!
   TBranch        *b_JetIDJERup;   //!
   TBranch        *b_Jets_jerFactor;   //!
   TBranch        *b_Jets_jerFactorDown;   //!
   TBranch        *b_Jets_jerFactorUp;   //!
   TBranch        *b_Jets_origIndex;
   TBranch        *b_JetsAK8_jerFactor;
   TBranch        *b_JetsAK8_jerFactorDown;
   TBranch        *b_JetsAK8_jerFactorUp;
   TBranch        *b_JetsAK8Clean_origIndex;
   TBranch        *b_JetsAK8JECup_jerFactor;
   TBranch        *b_JetsAK8JECup_origIndex;
   TBranch        *b_JetsAK8JERdown_origIndex;
   TBranch        *b_JetsAK8JERup_origIndex;
   TBranch        *b_JetsJECdown_jerFactor;
   TBranch        *b_JetsJECdown_origIndex;
   TBranch        *b_JetsJECup_jerFactor;
   TBranch        *b_JetsJECup_origIndex;
   TBranch        *b_JetsJERdown_origIndex;
   TBranch        *b_JetsJERup_origIndex;
   TBranch        *b_METDown;
   TBranch        *b_METOrigDown;
   TBranch        *b_METOrigUp;
   TBranch        *b_METPhiDown;
   TBranch        *b_METPhiOrigDown;
   TBranch        *b_METPhiOrigUp;
   TBranch        *b_METPhiUp;
   TBranch        *b_METUp;
   TBranch        *b_PDFweights;
   TBranch        *b_PSweights;
   TBranch        *b_ScaleWeights;


   TreeMakerTopiary(TTree *tree=0,int sampt=0,int year=0);
   virtual ~TreeMakerTopiary();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, int sampt,int year);
   virtual void     Loop(std::string entry,float entry1,int entry2,int year);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeMakerTopiary_cxx
TreeMakerTopiary::TreeMakerTopiary(TTree *tree,int sampt,int year) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../dataHandling/2017/Run2017F-31Mar2018-v1.SingleMuon_0_chunk0_oct08_skim_2640267MiniAODEvents.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../dataHandling/2017/Run2017F-31Mar2018-v1.SingleMuon_0_chunk0_oct08_skim_2640267MiniAODEvents.root");
      }
      f->GetObject("PreSelection",tree);

   }
   Init(tree,sampt,year);
}

TreeMakerTopiary::~TreeMakerTopiary()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeMakerTopiary::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeMakerTopiary::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeMakerTopiary::Init(TTree *tree, int sampt, int year)
{
   //std::cout<<"The sample type is "<<sampt<<std::endl;//this works
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Electrons = 0;
   Electrons_charge = 0;
   Electrons_mediumID = 0;
   Electrons_MTW = 0;
   Electrons_passIso = 0;
   Electrons_tightID = 0;
   HLTMuonObjects = 0;
   Jets = 0;
   Jets_axismajor = 0;
   Jets_axisminor = 0;
   Jets_bDiscriminatorCSV = 0;
   Jets_bJetTagDeepCSVBvsAll = 0;
   Jets_chargedEmEnergyFraction = 0;
   Jets_chargedHadronEnergyFraction = 0;
   Jets_chargedHadronMultiplicity = 0;
   Jets_chargedMultiplicity = 0;
   Jets_electronEnergyFraction = 0;
   Jets_electronMultiplicity = 0;
   Jets_hadronFlavor = 0;
   Jets_hfEMEnergyFraction = 0;
   Jets_hfHadronEnergyFraction = 0;
   Jets_HTMask = 0;
   Jets_ID = 0;
   Jets_jecFactor = 0;
   Jets_jecUnc = 0;
   Jets_LeptonMask = 0;
   Jets_MHTMask = 0;
   Jets_multiplicity = 0;
   Jets_muonEnergyFraction = 0;
   Jets_muonMultiplicity = 0;
   Jets_neutralEmEnergyFraction = 0;
   Jets_neutralHadronEnergyFraction = 0;
   Jets_neutralHadronMultiplicity = 0;
   Jets_neutralMultiplicity = 0;
   Jets_partonFlavor = 0;
   Jets_photonEnergyFraction = 0;
   Jets_photonMultiplicity = 0;
   Jets_ptD = 0;
   Jets_qgLikelihood = 0;
   JetsAK8 = 0;
   JetsAK8_axismajor = 0;
   JetsAK8_axisminor = 0;
   JetsAK8_DeepMassDecorrelTagbbvsLight = 0;
   JetsAK8_DeepMassDecorrelTagHbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagTvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagWvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZHbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZvsQCD = 0;
   JetsAK8_DeepTagHbbvsQCD = 0;
   JetsAK8_DeepTagTvsQCD = 0;
   JetsAK8_DeepTagWvsQCD = 0;
   JetsAK8_DeepTagZbbvsQCD = 0;
   JetsAK8_DeepTagZvsQCD = 0;
   JetsAK8_doubleBDiscriminator = 0;
   JetsAK8_ecfN2b1 = 0;
   JetsAK8_ecfN2b2 = 0;
   JetsAK8_ecfN3b1 = 0;
   JetsAK8_ecfN3b2 = 0;
   JetsAK8_girth = 0;
   JetsAK8_ID = 0;
   JetsAK8_jecFactor = 0;
   JetsAK8_jecUnc = 0;
   JetsAK8_multiplicity = 0;
   JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb = 0;
   JetsAK8_ptD = 0;
   JetsAK8_softDropMass = 0;
   JetsAK8_subjets = 0;
   JetsAK8Clean = 0;
   JetsAK8Clean_axismajor = 0;
   JetsAK8Clean_axisminor = 0;
   JetsAK8Clean_DeepMassDecorrelTagbbvsLight = 0;
   JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD = 0;
   JetsAK8Clean_DeepMassDecorrelTagTvsQCD = 0;
   JetsAK8Clean_DeepMassDecorrelTagWvsQCD = 0;
   JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD = 0;
   JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD = 0;
   JetsAK8Clean_DeepMassDecorrelTagZvsQCD = 0;
   JetsAK8Clean_DeepTagHbbvsQCD = 0;
   JetsAK8Clean_DeepTagTvsQCD = 0;
   JetsAK8Clean_DeepTagWvsQCD = 0;
   JetsAK8Clean_DeepTagZbbvsQCD = 0;
   JetsAK8Clean_DeepTagZvsQCD = 0;
   JetsAK8Clean_doubleBDiscriminator = 0;
   JetsAK8Clean_girth = 0;
   JetsAK8Clean_ID = 0;
   JetsAK8Clean_jecFactor = 0;
   JetsAK8Clean_multiplicity = 0;
   JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb = 0;
   JetsAK8Clean_ptD = 0;
   JetsAK8Clean_softDropMass = 0;
   JetsAK8Clean_subjets = 0;
   Jetsclean = 0;
   Jetsclean_bDiscriminatorCSV = 0;
   Jetsclean_bJetTagDeepCSVBvsAll = 0;
   Jetsclean_hadronFlavor = 0;
   Jetsclean_HTMask = 0;
   Jetsclean_ID = 0;
   Jetsclean_MHTMask = 0;
   Jetsclean_partonFlavor = 0;
   Muons = 0;
   Muons_charge = 0;
   Muons_mediumID = 0;
   Muons_MTW = 0;
   Muons_passIso = 0;
   Muons_tightID = 0;
   Photons = 0;
   Photons_electronFakes = 0;
   Photons_fullID = 0;
   Photons_genMatched = 0;
   Photons_hadTowOverEM = 0;
   Photons_hasPixelSeed = 0;
   Photons_isEB = 0;
   Photons_nonPrompt = 0;
   Photons_passElectronVeto = 0;
   Photons_pfChargedIso = 0;
   Photons_pfChargedIsoRhoCorr = 0;
   Photons_pfGammaIso = 0;
   Photons_pfGammaIsoRhoCorr = 0;
   Photons_pfNeutralIso = 0;
   Photons_pfNeutralIsoRhoCorr = 0;
   Photons_sigmaIetaIeta = 0;
   SelectedElectrons = 0;
   SelectedMuons = 0;
   TAPElectronTracks = 0;
   TAPElectronTracks_dxypv = 0;
   TAPElectronTracks_leptonMatch = 0;
   TAPElectronTracks_mT = 0;
   TAPElectronTracks_pfRelIso03chg = 0;
   TAPElectronTracks_trkiso = 0;
   TAPMuonTracks = 0;
   TAPMuonTracks_dxypv = 0;
   TAPMuonTracks_leptonMatch = 0;
   TAPMuonTracks_mT = 0;
   TAPMuonTracks_pfRelIso03chg = 0;
   TAPMuonTracks_trkiso = 0;
   TAPPionTracks = 0;
   TAPPionTracks_dxypv = 0;
   TAPPionTracks_leptonMatch = 0;
   TAPPionTracks_mT = 0;
   TAPPionTracks_pfRelIso03chg = 0;
   TAPPionTracks_trkiso = 0;
   TriggerPass = 0;
   TriggerPrescales = 0;
   TriggerVersion = 0;
   ZCandidates = 0;
   ZCandidatesMuMu = 0;
   ZCandidatesEE = 0;
   ZCandidatesEU= 0;
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("BadChargedCandidateFilter", &BadChargedCandidateFilter, &b_BadChargedCandidateFilter);
   fChain->SetBranchAddress("BadPFMuonFilter", &BadPFMuonFilter, &b_BadPFMuonFilter);
   fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
   fChain->SetBranchAddress("BTagsclean", &BTagsclean, &b_BTagsclean);
   fChain->SetBranchAddress("BTagsDeepCSV", &BTagsDeepCSV, &b_BTagsDeepCSV);
   fChain->SetBranchAddress("BTagsDeepCSVclean", &BTagsDeepCSVclean, &b_BTagsDeepCSVclean);
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("CaloMETPhi", &CaloMETPhi, &b_CaloMETPhi);
   fChain->SetBranchAddress("CSCTightHaloFilter", &CSCTightHaloFilter, &b_CSCTightHaloFilter);
   fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
   fChain->SetBranchAddress("DeltaPhi1clean", &DeltaPhi1clean, &b_DeltaPhi1clean);
   fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
   fChain->SetBranchAddress("DeltaPhi2clean", &DeltaPhi2clean, &b_DeltaPhi2clean);
   fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
   fChain->SetBranchAddress("DeltaPhi3clean", &DeltaPhi3clean, &b_DeltaPhi3clean);
   fChain->SetBranchAddress("DeltaPhi4", &DeltaPhi4, &b_DeltaPhi4);
   fChain->SetBranchAddress("DeltaPhi4clean", &DeltaPhi4clean, &b_DeltaPhi4clean);
   fChain->SetBranchAddress("ecalBadCalibFilter", &ecalBadCalibFilter, &b_ecalBadCalibFilter);
   fChain->SetBranchAddress("ecalBadCalibReducedExtraFilter", &ecalBadCalibReducedExtraFilter, &b_ecalBadCalibReducedExtraFilter);
   fChain->SetBranchAddress("ecalBadCalibReducedFilter", &ecalBadCalibReducedFilter, &b_ecalBadCalibReducedFilter);
   fChain->SetBranchAddress("EcalDeadCellBoundaryEnergyFilter", &EcalDeadCellBoundaryEnergyFilter, &b_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("eeBadScFilter", &eeBadScFilter, &b_eeBadScFilter);
   fChain->SetBranchAddress("Electrons", &Electrons, &b_Electrons);
   fChain->SetBranchAddress("Electrons_charge", &Electrons_charge, &b_Electrons_charge);
   fChain->SetBranchAddress("Electrons_mediumID", &Electrons_mediumID, &b_Electrons_mediumID);
   fChain->SetBranchAddress("Electrons_MTW", &Electrons_MTW, &b_Electrons_MTW);
   fChain->SetBranchAddress("Electrons_passIso", &Electrons_passIso, &b_Electrons_passIso);
   fChain->SetBranchAddress("Electrons_tightID", &Electrons_tightID, &b_Electrons_tightID);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("globalSuperTightHalo2016Filter", &globalSuperTightHalo2016Filter, &b_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("globalTightHalo2016Filter", &globalTightHalo2016Filter, &b_globalTightHalo2016Filter);
   fChain->SetBranchAddress("hasGenPromptPhoton", &hasGenPromptPhoton, &b_hasGenPromptPhoton);
   fChain->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter, &b_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter, &b_HBHENoiseFilter);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("HT5", &HT5, &b_HT5);
   fChain->SetBranchAddress("HT5clean", &HT5clean, &b_HT5clean);
   fChain->SetBranchAddress("HTclean", &HTclean, &b_HTclean);
   fChain->SetBranchAddress("isoElectronTracks", &isoElectronTracks, &b_isoElectronTracks);
   fChain->SetBranchAddress("isoElectronTracksclean", &isoElectronTracksclean, &b_isoElectronTracksclean);
   fChain->SetBranchAddress("isoMuonTracks", &isoMuonTracks, &b_isoMuonTracks);
   fChain->SetBranchAddress("isoMuonTracksclean", &isoMuonTracksclean, &b_isoMuonTracksclean);
   fChain->SetBranchAddress("isoPionTracks", &isoPionTracks, &b_isoPionTracks);
   fChain->SetBranchAddress("isoPionTracksclean", &isoPionTracksclean, &b_isoPionTracksclean);
   fChain->SetBranchAddress("JetID", &JetID, &b_JetID);
   fChain->SetBranchAddress("JetIDAK8", &JetIDAK8, &b_JetIDAK8);
   fChain->SetBranchAddress("JetIDAK8Clean", &JetIDAK8Clean, &b_JetIDAK8Clean);
   fChain->SetBranchAddress("JetIDclean", &JetIDclean, &b_JetIDclean);
   fChain->SetBranchAddress("Jets", &Jets, &b_Jets);
   fChain->SetBranchAddress("Jets_axismajor", &Jets_axismajor, &b_Jets_axismajor);
   fChain->SetBranchAddress("Jets_axisminor", &Jets_axisminor, &b_Jets_axisminor);
   fChain->SetBranchAddress("Jets_bDiscriminatorCSV", &Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVBvsAll", &Jets_bJetTagDeepCSVBvsAll, &b_Jets_bJetTagDeepCSVBvsAll);
   fChain->SetBranchAddress("Jets_chargedEmEnergyFraction", &Jets_chargedEmEnergyFraction, &b_Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronEnergyFraction", &Jets_chargedHadronEnergyFraction, &b_Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronMultiplicity", &Jets_chargedHadronMultiplicity, &b_Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jets_chargedMultiplicity", &Jets_chargedMultiplicity, &b_Jets_chargedMultiplicity);
   fChain->SetBranchAddress("Jets_electronEnergyFraction", &Jets_electronEnergyFraction, &b_Jets_electronEnergyFraction);
   fChain->SetBranchAddress("Jets_electronMultiplicity", &Jets_electronMultiplicity, &b_Jets_electronMultiplicity);
   fChain->SetBranchAddress("Jets_hadronFlavor", &Jets_hadronFlavor, &b_Jets_hadronFlavor);
   fChain->SetBranchAddress("Jets_hfEMEnergyFraction", &Jets_hfEMEnergyFraction, &b_Jets_hfEMEnergyFraction);
   fChain->SetBranchAddress("Jets_hfHadronEnergyFraction", &Jets_hfHadronEnergyFraction, &b_Jets_hfHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_HTMask", &Jets_HTMask, &b_Jets_HTMask);
   fChain->SetBranchAddress("Jets_ID", &Jets_ID, &b_Jets_ID);
   fChain->SetBranchAddress("Jets_jecFactor", &Jets_jecFactor, &b_Jets_jecFactor);
   fChain->SetBranchAddress("Jets_jecUnc", &Jets_jecUnc, &b_Jets_jecUnc);
   fChain->SetBranchAddress("Jets_LeptonMask", &Jets_LeptonMask, &b_Jets_LeptonMask);
   fChain->SetBranchAddress("Jets_MHTMask", &Jets_MHTMask, &b_Jets_MHTMask);
   fChain->SetBranchAddress("Jets_multiplicity", &Jets_multiplicity, &b_Jets_multiplicity);
   fChain->SetBranchAddress("Jets_muonEnergyFraction", &Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
   fChain->SetBranchAddress("Jets_muonMultiplicity", &Jets_muonMultiplicity, &b_Jets_muonMultiplicity);
   fChain->SetBranchAddress("Jets_neutralEmEnergyFraction", &Jets_neutralEmEnergyFraction, &b_Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronEnergyFraction", &Jets_neutralHadronEnergyFraction, &b_Jets_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronMultiplicity", &Jets_neutralHadronMultiplicity, &b_Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jets_neutralMultiplicity", &Jets_neutralMultiplicity, &b_Jets_neutralMultiplicity);
   fChain->SetBranchAddress("Jets_partonFlavor", &Jets_partonFlavor, &b_Jets_partonFlavor);
   fChain->SetBranchAddress("Jets_photonEnergyFraction", &Jets_photonEnergyFraction, &b_Jets_photonEnergyFraction);
   fChain->SetBranchAddress("Jets_photonMultiplicity", &Jets_photonMultiplicity, &b_Jets_photonMultiplicity);
   fChain->SetBranchAddress("Jets_ptD", &Jets_ptD, &b_Jets_ptD);
   fChain->SetBranchAddress("Jets_qgLikelihood", &Jets_qgLikelihood, &b_Jets_qgLikelihood);
   fChain->SetBranchAddress("JetsAK8", &JetsAK8, &b_JetsAK8);
   fChain->SetBranchAddress("JetsAK8_axismajor", &JetsAK8_axismajor, &b_JetsAK8_axismajor);
   fChain->SetBranchAddress("JetsAK8_axisminor", &JetsAK8_axisminor, &b_JetsAK8_axisminor);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagbbvsLight", &JetsAK8_DeepMassDecorrelTagbbvsLight, &b_JetsAK8_DeepMassDecorrelTagbbvsLight);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagHbbvsQCD", &JetsAK8_DeepMassDecorrelTagHbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagTvsQCD", &JetsAK8_DeepMassDecorrelTagTvsQCD, &b_JetsAK8_DeepMassDecorrelTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagWvsQCD", &JetsAK8_DeepMassDecorrelTagWvsQCD, &b_JetsAK8_DeepMassDecorrelTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZbbvsQCD", &JetsAK8_DeepMassDecorrelTagZbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZHbbvsQCD", &JetsAK8_DeepMassDecorrelTagZHbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagZHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZvsQCD", &JetsAK8_DeepMassDecorrelTagZvsQCD, &b_JetsAK8_DeepMassDecorrelTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagHbbvsQCD", &JetsAK8_DeepTagHbbvsQCD, &b_JetsAK8_DeepTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagTvsQCD", &JetsAK8_DeepTagTvsQCD, &b_JetsAK8_DeepTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagWvsQCD", &JetsAK8_DeepTagWvsQCD, &b_JetsAK8_DeepTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagZbbvsQCD", &JetsAK8_DeepTagZbbvsQCD, &b_JetsAK8_DeepTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagZvsQCD", &JetsAK8_DeepTagZvsQCD, &b_JetsAK8_DeepTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8_doubleBDiscriminator", &JetsAK8_doubleBDiscriminator, &b_JetsAK8_doubleBDiscriminator);
   fChain->SetBranchAddress("JetsAK8_ecfN2b1", &JetsAK8_ecfN2b1, &b_JetsAK8_ecfN2b1);
   fChain->SetBranchAddress("JetsAK8_ecfN2b2", &JetsAK8_ecfN2b2, &b_JetsAK8_ecfN2b2);
   fChain->SetBranchAddress("JetsAK8_ecfN3b1", &JetsAK8_ecfN3b1, &b_JetsAK8_ecfN3b1);
   fChain->SetBranchAddress("JetsAK8_ecfN3b2", &JetsAK8_ecfN3b2, &b_JetsAK8_ecfN3b2);
   fChain->SetBranchAddress("JetsAK8_girth", &JetsAK8_girth, &b_JetsAK8_girth);
   fChain->SetBranchAddress("JetsAK8_ID", &JetsAK8_ID, &b_JetsAK8_ID);
   fChain->SetBranchAddress("JetsAK8_jecFactor", &JetsAK8_jecFactor, &b_JetsAK8_jecFactor);
   fChain->SetBranchAddress("JetsAK8_jecUnc", &JetsAK8_jecUnc, &b_JetsAK8_jecUnc);
   fChain->SetBranchAddress("JetsAK8_multiplicity", &JetsAK8_multiplicity, &b_JetsAK8_multiplicity);
   fChain->SetBranchAddress("JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb", &JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb, &b_JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb);
   fChain->SetBranchAddress("JetsAK8_ptD", &JetsAK8_ptD, &b_JetsAK8_ptD);
   fChain->SetBranchAddress("JetsAK8_softDropMass", &JetsAK8_softDropMass, &b_JetsAK8_softDropMass);
   fChain->SetBranchAddress("JetsAK8_subjets", &JetsAK8_subjets, &b_JetsAK8_subjets);
   fChain->SetBranchAddress("JetsAK8Clean", &JetsAK8Clean, &b_JetsAK8Clean);
   fChain->SetBranchAddress("JetsAK8Clean_axismajor", &JetsAK8Clean_axismajor, &b_JetsAK8Clean_axismajor);
   fChain->SetBranchAddress("JetsAK8Clean_axisminor", &JetsAK8Clean_axisminor, &b_JetsAK8Clean_axisminor);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagbbvsLight", &JetsAK8Clean_DeepMassDecorrelTagbbvsLight, &b_JetsAK8Clean_DeepMassDecorrelTagbbvsLight);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD", &JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagTvsQCD", &JetsAK8Clean_DeepMassDecorrelTagTvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagWvsQCD", &JetsAK8Clean_DeepMassDecorrelTagWvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD", &JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD", &JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepMassDecorrelTagZvsQCD", &JetsAK8Clean_DeepMassDecorrelTagZvsQCD, &b_JetsAK8Clean_DeepMassDecorrelTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepTagHbbvsQCD", &JetsAK8Clean_DeepTagHbbvsQCD, &b_JetsAK8Clean_DeepTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepTagTvsQCD", &JetsAK8Clean_DeepTagTvsQCD, &b_JetsAK8Clean_DeepTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepTagWvsQCD", &JetsAK8Clean_DeepTagWvsQCD, &b_JetsAK8Clean_DeepTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepTagZbbvsQCD", &JetsAK8Clean_DeepTagZbbvsQCD, &b_JetsAK8Clean_DeepTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_DeepTagZvsQCD", &JetsAK8Clean_DeepTagZvsQCD, &b_JetsAK8Clean_DeepTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8Clean_doubleBDiscriminator", &JetsAK8Clean_doubleBDiscriminator, &b_JetsAK8Clean_doubleBDiscriminator);
   fChain->SetBranchAddress("JetsAK8Clean_girth", &JetsAK8Clean_girth, &b_JetsAK8Clean_girth);
   fChain->SetBranchAddress("JetsAK8Clean_ID", &JetsAK8Clean_ID, &b_JetsAK8Clean_ID);
   fChain->SetBranchAddress("JetsAK8Clean_jecFactor", &JetsAK8Clean_jecFactor, &b_JetsAK8Clean_jecFactor);
   fChain->SetBranchAddress("JetsAK8Clean_multiplicity", &JetsAK8Clean_multiplicity, &b_JetsAK8Clean_multiplicity);
   fChain->SetBranchAddress("JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb", &JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb, &b_JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb);
   fChain->SetBranchAddress("JetsAK8Clean_ptD", &JetsAK8Clean_ptD, &b_JetsAK8Clean_ptD);
   fChain->SetBranchAddress("JetsAK8Clean_softDropMass", &JetsAK8Clean_softDropMass, &b_JetsAK8Clean_softDropMass);
   fChain->SetBranchAddress("JetsAK8Clean_subjets", &JetsAK8Clean_subjets, &b_JetsAK8Clean_subjets);
   fChain->SetBranchAddress("Jetsclean", &Jetsclean, &b_Jetsclean);
   fChain->SetBranchAddress("Jetsclean_bDiscriminatorCSV", &Jetsclean_bDiscriminatorCSV, &b_Jetsclean_bDiscriminatorCSV);
   fChain->SetBranchAddress("Jetsclean_bJetTagDeepCSVBvsAll", &Jetsclean_bJetTagDeepCSVBvsAll, &b_Jetsclean_bJetTagDeepCSVBvsAll);
   fChain->SetBranchAddress("Jetsclean_hadronFlavor", &Jetsclean_hadronFlavor, &b_Jetsclean_hadronFlavor);
   fChain->SetBranchAddress("Jetsclean_HTMask", &Jetsclean_HTMask, &b_Jetsclean_HTMask);
   fChain->SetBranchAddress("Jetsclean_ID", &Jetsclean_ID, &b_Jetsclean_ID);
   fChain->SetBranchAddress("Jetsclean_MHTMask", &Jetsclean_MHTMask, &b_Jetsclean_MHTMask);
   fChain->SetBranchAddress("Jetsclean_partonFlavor", &Jetsclean_partonFlavor, &b_Jetsclean_partonFlavor);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METclean", &METclean, &b_METclean);
   fChain->SetBranchAddress("METcleanOrig", &METcleanOrig, &b_METcleanOrig);
   fChain->SetBranchAddress("METOrig", &METOrig, &b_METOrig);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METPhiclean", &METPhiclean, &b_METPhiclean);
   fChain->SetBranchAddress("METPhicleanOrig", &METPhicleanOrig, &b_METPhicleanOrig);
   fChain->SetBranchAddress("METPhiOrig", &METPhiOrig, &b_METPhiOrig);
   fChain->SetBranchAddress("METSignificance", &METSignificance, &b_METSignificance);
   fChain->SetBranchAddress("METSignificanceclean", &METSignificanceclean, &b_METSignificanceclean);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTclean", &MHTclean, &b_MHTclean);
   fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
   fChain->SetBranchAddress("MHTPhiclean", &MHTPhiclean, &b_MHTPhiclean);
   fChain->SetBranchAddress("Muons", &Muons, &b_Muons);
   fChain->SetBranchAddress("Muons_charge", &Muons_charge, &b_Muons_charge);
   fChain->SetBranchAddress("Muons_mediumID", &Muons_mediumID, &b_Muons_mediumID);
   fChain->SetBranchAddress("Muons_MTW", &Muons_MTW, &b_Muons_MTW);
   fChain->SetBranchAddress("Muons_passIso", &Muons_passIso, &b_Muons_passIso);
   fChain->SetBranchAddress("Muons_tightID", &Muons_tightID, &b_Muons_tightID);
   fChain->SetBranchAddress("nAllVertices", &nAllVertices, &b_nAllVertices);
   fChain->SetBranchAddress("NElectrons", &NElectrons, &b_NElectrons);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NJetsclean", &NJetsclean, &b_NJetsclean);
   fChain->SetBranchAddress("NMuons", &NMuons, &b_NMuons);
   fChain->SetBranchAddress("NonPrefiringProb", &NonPrefiringProb, &b_NonPrefiringProb);
   fChain->SetBranchAddress("NonPrefiringProbDown", &NonPrefiringProbDown, &b_NonPrefiringProbDown);
   fChain->SetBranchAddress("NonPrefiringProbUp", &NonPrefiringProbUp, &b_NonPrefiringProbUp);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("PFCaloMETRatio", &PFCaloMETRatio, &b_PFCaloMETRatio);
   fChain->SetBranchAddress("Photons", &Photons, &b_Photons);
   fChain->SetBranchAddress("Photons_electronFakes", &Photons_electronFakes, &b_Photons_electronFakes);
   fChain->SetBranchAddress("Photons_fullID", &Photons_fullID, &b_Photons_fullID);
   fChain->SetBranchAddress("Photons_genMatched", &Photons_genMatched, &b_Photons_genMatched);
   fChain->SetBranchAddress("Photons_hadTowOverEM", &Photons_hadTowOverEM, &b_Photons_hadTowOverEM);
   fChain->SetBranchAddress("Photons_hasPixelSeed", &Photons_hasPixelSeed, &b_Photons_hasPixelSeed);
   fChain->SetBranchAddress("Photons_isEB", &Photons_isEB, &b_Photons_isEB);
   fChain->SetBranchAddress("Photons_nonPrompt", &Photons_nonPrompt, &b_Photons_nonPrompt);
   fChain->SetBranchAddress("Photons_passElectronVeto", &Photons_passElectronVeto, &b_Photons_passElectronVeto);
   fChain->SetBranchAddress("Photons_pfChargedIso", &Photons_pfChargedIso, &b_Photons_pfChargedIso);
   fChain->SetBranchAddress("Photons_pfChargedIsoRhoCorr", &Photons_pfChargedIsoRhoCorr, &b_Photons_pfChargedIsoRhoCorr);
   fChain->SetBranchAddress("Photons_pfGammaIso", &Photons_pfGammaIso, &b_Photons_pfGammaIso);
   fChain->SetBranchAddress("Photons_pfGammaIsoRhoCorr", &Photons_pfGammaIsoRhoCorr, &b_Photons_pfGammaIsoRhoCorr);
   fChain->SetBranchAddress("Photons_pfNeutralIso", &Photons_pfNeutralIso, &b_Photons_pfNeutralIso);
   fChain->SetBranchAddress("Photons_pfNeutralIsoRhoCorr", &Photons_pfNeutralIsoRhoCorr, &b_Photons_pfNeutralIsoRhoCorr);
   fChain->SetBranchAddress("Photons_sigmaIetaIeta", &Photons_sigmaIetaIeta, &b_Photons_sigmaIetaIeta);
   fChain->SetBranchAddress("PrimaryVertexFilter", &PrimaryVertexFilter, &b_PrimaryVertexFilter);
   fChain->SetBranchAddress("SelectedElectrons", &SelectedElectrons, &b_SelectedElectrons);
   fChain->SetBranchAddress("SelectedMuons", &SelectedMuons, &b_SelectedMuons);
   fChain->SetBranchAddress("TAPElectronTracks", &TAPElectronTracks, &b_TAPElectronTracks);
   fChain->SetBranchAddress("TAPElectronTracks_dxypv", &TAPElectronTracks_dxypv, &b_TAPElectronTracks_dxypv);
   fChain->SetBranchAddress("TAPElectronTracks_leptonMatch", &TAPElectronTracks_leptonMatch, &b_TAPElectronTracks_leptonMatch);
   fChain->SetBranchAddress("TAPElectronTracks_mT", &TAPElectronTracks_mT, &b_TAPElectronTracks_mT);
   fChain->SetBranchAddress("TAPElectronTracks_pfRelIso03chg", &TAPElectronTracks_pfRelIso03chg, &b_TAPElectronTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPElectronTracks_trkiso", &TAPElectronTracks_trkiso, &b_TAPElectronTracks_trkiso);
   fChain->SetBranchAddress("TAPMuonTracks", &TAPMuonTracks, &b_TAPMuonTracks);
   fChain->SetBranchAddress("TAPMuonTracks_dxypv", &TAPMuonTracks_dxypv, &b_TAPMuonTracks_dxypv);
   fChain->SetBranchAddress("TAPMuonTracks_leptonMatch", &TAPMuonTracks_leptonMatch, &b_TAPMuonTracks_leptonMatch);
   fChain->SetBranchAddress("TAPMuonTracks_mT", &TAPMuonTracks_mT, &b_TAPMuonTracks_mT);
   fChain->SetBranchAddress("TAPMuonTracks_pfRelIso03chg", &TAPMuonTracks_pfRelIso03chg, &b_TAPMuonTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPMuonTracks_trkiso", &TAPMuonTracks_trkiso, &b_TAPMuonTracks_trkiso);
   fChain->SetBranchAddress("TAPPionTracks", &TAPPionTracks, &b_TAPPionTracks);
   fChain->SetBranchAddress("TAPPionTracks_dxypv", &TAPPionTracks_dxypv, &b_TAPPionTracks_dxypv);
   fChain->SetBranchAddress("TAPPionTracks_leptonMatch", &TAPPionTracks_leptonMatch, &b_TAPPionTracks_leptonMatch);
   fChain->SetBranchAddress("TAPPionTracks_mT", &TAPPionTracks_mT, &b_TAPPionTracks_mT);
   fChain->SetBranchAddress("TAPPionTracks_pfRelIso03chg", &TAPPionTracks_pfRelIso03chg, &b_TAPPionTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPPionTracks_trkiso", &TAPPionTracks_trkiso, &b_TAPPionTracks_trkiso);
   fChain->SetBranchAddress("TriggerPass", &TriggerPass, &b_TriggerPass);
   fChain->SetBranchAddress("TriggerPrescales", &TriggerPrescales, &b_TriggerPrescales);
   fChain->SetBranchAddress("TriggerVersion", &TriggerVersion, &b_TriggerVersion);
   fChain->SetBranchAddress("ZCandidates", &ZCandidates, &b_ZCandidates);
   fChain->SetBranchAddress("ZCandidatesMuMu", &ZCandidatesMuMu, &b_ZCandidatesMuMu);
   fChain->SetBranchAddress("ZCandidatesEE", &ZCandidatesEE, &b_ZCandidatesEE);
   fChain->SetBranchAddress("ZCandidatesEU", &ZCandidatesEU, &b_ZCandidatesEU);

   if (sampt == 0) {
     fChain->SetBranchAddress("HLTMuonObjects", &HLTMuonObjects, &b_HLTMuonObjects);
     fChain->SetBranchAddress("HTOnline", &HTOnline, &b_HTOnline);
     fChain->SetBranchAddress("MHTOnline", &MHTOnline, &b_MHTOnline);
     fChain->SetBranchAddress("PrescaleWeightHT", &PrescaleWeightHT, &b_PrescaleWeightHT);
   }

   if (sampt != 0) {
     //Initialize the MC only branches
     GenElectrons = 0;
     GenHT = 0;
     GenJets = 0;
     GenJetsAK8 = 0;
     GenJetsAK8_multiplicity = 0;
     GenJetsAK8_softDropMass = 0;
     GenMET = 0;
     GenMETPhi = 0;
     GenMuons = 0;
     GenParticles = 0;
     GenParticles_ParentId = 0;
     GenParticles_ParentIdx = 0;
     GenParticles_PdgId = 0;
     GenParticles_Status = 0;
     GenTaus = 0;
     GenTaus_had = 0;
     HT5JECdown = 0;
     HT5JECup = 0;
     HT5JERdown = 0;
     HT5JERup = 0;
     HTJECdown = 0;
     HTJECup = 0;
     HTJERdown = 0;
     HTJERup = 0;
     JetIDAK8JECdown = 0;
     JetIDAK8JECup = 0;
     JetIDAK8JERdown = 0;
     JetIDAK8JERup = 0;
     JetIDJECdown = 0;
     JetIDJECup = 0;
     JetIDJERdown = 0;
     JetIDJERup = 0;
     Jets_jerFactor = 0;
     Jets_jerFactorDown = 0;
     Jets_jerFactorUp = 0;
     Jets_origIndex = 0;
     JetsAK8_jerFactor = 0;
     JetsAK8_jerFactorDown = 0;
     JetsAK8_jerFactorUp = 0;
     JetsAK8Clean_origIndex = 0;
     JetsAK8JECup_jerFactor = 0;
     JetsAK8JECup_origIndex = 0;
     JetsAK8JERdown_origIndex = 0;
     JetsAK8JERup_origIndex = 0;
     JetsJECdown_jerFactor = 0;
     JetsJECdown_origIndex = 0;
     JetsJECup_jerFactor = 0;
     JetsJECup_origIndex = 0;
     JetsJERdown_origIndex = 0;
     JetsJERup_origIndex = 0;
     METDown = 0;
     METOrigDown = 0;
     METOrigUp = 0;
     METPhiDown = 0;
     METPhiOrigDown = 0;
     METPhiOrigUp = 0;
     METPhiUp = 0;
     METUp = 0;
     PDFweights = 0;
     PSweights = 0;
     ScaleWeights = 0;

     //Branches
     //Too many branches confuse root, so only maintaining ones we will
     // a) need
     // b) are big enough to need to turn off
     fChain->SetBranchAddress("GenElectrons", &GenElectrons, &b_GenElectrons);
     fChain->SetBranchAddress("GenHT", &GenHT, &b_GenHT);
     fChain->SetBranchAddress("GenJets", &GenJets, &b_GenJets);
     fChain->SetBranchAddress("GenJetsAK8", &GenJetsAK8, &b_GenJetsAK8);
     fChain->SetBranchAddress("GenJetsAK8_multiplicity", &GenJetsAK8_multiplicity, &b_GenJetsAK8_multiplicity);
     fChain->SetBranchAddress("GenJetsAK8_softDropMass", &GenJetsAK8_softDropMass, &b_GenJetsAK8_softDropMass);
     fChain->SetBranchAddress("GenMET", &GenMET, &b_GenMET);
     fChain->SetBranchAddress("GenMETPhi", &GenMETPhi, &b_GenMETPhi);
     fChain->SetBranchAddress("GenMuons", &GenMuons, &b_GenMuons);
     fChain->SetBranchAddress("GenParticles", &GenParticles, &b_GenParticles);
     fChain->SetBranchAddress("GenParticles_ParentId", &GenParticles_ParentId, &b_GenParticles_ParentId);
     fChain->SetBranchAddress("GenParticles_ParentIdx", &GenParticles_ParentIdx, &b_GenParticles_ParentIdx);
     fChain->SetBranchAddress("GenParticles_PdgId", &GenParticles_PdgId, &b_GenParticles_PdgId);
     fChain->SetBranchAddress("GenParticles_Status", &GenParticles_Status, &b_GenParticles_Status);
     fChain->SetBranchAddress("GenTaus", &GenTaus, &b_GenTaus);
     fChain->SetBranchAddress("GenTaus_had", &GenTaus_had, &b_GenTaus_had);
     fChain->SetBranchAddress("HT5JECdown", &HT5JECdown, &b_HT5JECdown);
     fChain->SetBranchAddress("HT5JECup", &HT5JECup, &b_HT5JECup);
     fChain->SetBranchAddress("HT5JERdown", &HT5JERdown, &b_HT5JERdown);
     fChain->SetBranchAddress("HT5JERup", &HT5JERup, &b_HT5JERup);
     fChain->SetBranchAddress("HTJECdown", &HTJECdown, &b_HTJECdown);
     fChain->SetBranchAddress("HTJECup", &HTJECup, &b_HTJECup);
     fChain->SetBranchAddress("HTJERdown", &HTJERdown, &b_HTJERdown);
     fChain->SetBranchAddress("HTJERup", &HTJERup, &b_HTJERup);
     fChain->SetBranchAddress("JetIDAK8JECdown", &JetIDAK8JECdown, &b_JetIDAK8JECdown);
     fChain->SetBranchAddress("JetIDAK8JECup", &JetIDAK8JECup, &b_JetIDAK8JECup);
     fChain->SetBranchAddress("JetIDAK8JERdown", &JetIDAK8JERdown, &b_JetIDAK8JERdown);
     fChain->SetBranchAddress("JetIDAK8JERup", &JetIDAK8JERup, &b_JetIDAK8JERup);
     fChain->SetBranchAddress("JetIDJECdown", &JetIDJECdown, &b_JetIDJECdown);
     fChain->SetBranchAddress("JetIDJECup", &JetIDJECup, &b_JetIDJECup);
     fChain->SetBranchAddress("JetIDJERdown", &JetIDJERdown, &b_JetIDJERdown);
     fChain->SetBranchAddress("JetIDJERup", &JetIDJERup, &b_JetIDJERup);
     fChain->SetBranchAddress("Jets_jerFactor", &Jets_jerFactor, &b_Jets_jerFactor);
     fChain->SetBranchAddress("Jets_jerFactorDown", &Jets_jerFactorDown, &b_Jets_jerFactorDown);
     fChain->SetBranchAddress("Jets_jerFactorUp", &Jets_jerFactorUp, &b_Jets_jerFactorUp);
     fChain->SetBranchAddress("Jets_origIndex", &Jets_origIndex, &b_Jets_origIndex);
     fChain->SetBranchAddress("JetsAK8_jerFactor", &JetsAK8_jerFactor, &b_JetsAK8_jerFactor);
     fChain->SetBranchAddress("JetsAK8_jerFactorDown", &JetsAK8_jerFactorDown, &b_JetsAK8_jerFactorDown);
     fChain->SetBranchAddress("JetsAK8_jerFactorUp", &JetsAK8_jerFactorUp, &b_JetsAK8_jerFactorUp);
     fChain->SetBranchAddress("JetsAK8Clean_origIndex", &JetsAK8Clean_origIndex, &b_JetsAK8Clean_origIndex);
     fChain->SetBranchAddress("JetsAK8JECup_jerFactor", &JetsAK8JECup_jerFactor, &b_JetsAK8JECup_jerFactor);
     fChain->SetBranchAddress("JetsAK8JECup_origIndex", &JetsAK8JECup_origIndex, &b_JetsAK8JECup_origIndex);
     fChain->SetBranchAddress("JetsAK8JERdown_origIndex", &JetsAK8JERdown_origIndex, &b_JetsAK8JERdown_origIndex);
     fChain->SetBranchAddress("JetsAK8JERup_origIndex", &JetsAK8JERup_origIndex, &b_JetsAK8JERup_origIndex);
     fChain->SetBranchAddress("JetsJECdown_jerFactor", &JetsJECdown_jerFactor, &b_JetsJECdown_jerFactor);
     fChain->SetBranchAddress("JetsJECdown_origIndex", &JetsJECdown_origIndex, &b_JetsJECdown_origIndex);
     fChain->SetBranchAddress("JetsJECup_jerFactor", &JetsJECup_jerFactor, &b_JetsJECup_jerFactor);
     fChain->SetBranchAddress("JetsJECup_origIndex", &JetsJECup_origIndex, &b_JetsJECup_origIndex);
     fChain->SetBranchAddress("JetsJERdown_origIndex", &JetsJERdown_origIndex, &b_JetsJERdown_origIndex);
     fChain->SetBranchAddress("JetsJERup_origIndex", &JetsJERup_origIndex, &b_JetsJERup_origIndex);
     fChain->SetBranchAddress("METDown", &METDown, &b_METDown);
     fChain->SetBranchAddress("METOrigDown", &METOrigDown, &b_METOrigDown);
     fChain->SetBranchAddress("METOrigUp", &METOrigUp, &b_METOrigUp);
     fChain->SetBranchAddress("METPhiDown", &METPhiDown, &b_METPhiDown);
     fChain->SetBranchAddress("METPhiOrigDown", &METPhiOrigDown, &b_METPhiOrigDown);
     fChain->SetBranchAddress("METPhiOrigUp", &METPhiOrigUp, &b_METPhiOrigUp);
     fChain->SetBranchAddress("METPhiUp", &METPhiUp, &b_METPhiUp);
     fChain->SetBranchAddress("METUp", &METUp, &b_METUp);
     fChain->SetBranchAddress("PDFweights", &PDFweights, &b_PDFweights);
     fChain->SetBranchAddress("PSweights", &PSweights, &b_PSweights);
     fChain->SetBranchAddress("ScaleWeights", &ScaleWeights, &b_ScaleWeights);
   }
   
   Notify();
}

Bool_t TreeMakerTopiary::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeMakerTopiary::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeMakerTopiary::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeMakerTopiary_cxx

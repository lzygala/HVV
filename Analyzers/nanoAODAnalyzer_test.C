#include<TApplication.h>
#include<TFile.h>
#include<TMath.h>
#include<TMinuit.h>
#include<TROOT.h>
#include<TSystem.h>
#include<TTree.h>
#include<TVector2.h>
#include<Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TLegend.h>
#include<TLine.h>
#include<TH2F.h>
#include<TPaveText.h>
#include<TPaveStats.h>
#include<TStyle.h>

#include<algorithm>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iostream>
#include<mutex>
#include<string>
#include<thread>
#include<vector>

#include<sys/stat.h>
#include<errno.h>

#ifdef __MAKECINT__ 
#pragma link C++ class vector<vector<double> >+;  
#pragma link C++ class vector<Float_t >+; 
#pragma link C++ class vector<vector<int> >+;  
#pragma link C++ class vector<vector<bool> >+; 
#pragma link C++ class vector<vector<map<int,Float_t>> >+; 
#endif

using namespace std;

/* Variables */


    TH1F* h_CaloMET_phi;
    TH1F* h_CaloMET_pt;
    TH1F* h_CaloMET_sumEt;

    TH1F* h_Electron_eta;
    TH1F* h_Electron_hoe;
    TH1F* h_Electron_mass;
    TH1F* h_Electron_phi;
    TH1F* h_Electron_pt;
    TH1F* h_Electron_r9;
    TH1F* h_Electron_sieie;

    TH1F* h_Muon_eta;
    TH1F* h_Muon_mass;
    TH1F* h_Muon_phi;
    TH1F* h_Muon_pt;

    TH1F* h_Jet_eta;
    TH1F* h_Jet_mass;
    TH1F* h_Jet_phi;
    TH1F* h_Jet_pt;
    TH1F* h_nJets;

    TH1F* h_Gen_Jet_eta;
    TH1F* h_Gen_Jet_mass;
    TH1F* h_Gen_Jet_phi;
    TH1F* h_Gen_Jet_pt;
    TH1F* h_nGenJets;

    TH1F* h_Gen_Part_eta;
    TH1F* h_Gen_Part_mass;
    TH1F* h_Gen_Part_phi;
    TH1F* h_Gen_Part_pt;
    TH1F* h_nGenParts;

    TH1F* h_fatJet_eta;
    TH1F* h_fatJet_mass;
    TH1F* h_fatJet_phi;
    TH1F* h_fatJet_pt;
    TH1F* h_nFatJets;

    TH1F* h_Leading_Jet_eta;
    TH1F* h_Leading_Jet_mass;
    TH1F* h_Leading_Jet_phi;
    TH1F* h_Leading_Jet_pt;

    TH1F* h_Trailing_Jet_eta;
    TH1F* h_Trailing_Jet_mass;
    TH1F* h_Trailing_Jet_phi;
    TH1F* h_Trailing_Jet_pt;

    TH1F* h_Jet_invariantMass;
    TH1F* h_Jet_etaSeparation;

    TH1F* h_Leading_FatJet_eta;
    TH1F* h_Leading_FatJet_mass;
    TH1F* h_Leading_FatJet_phi;
    TH1F* h_Leading_FatJet_pt;
    TH1F* h_Leading_FatJet_deepTag_H;
    TH1F* h_Leading_FatJet_deepTag_WvsQCD;
    TH1F* h_Leading_FatJet_deepTagMD_bbvsLight;

    TH1F* h_Matched_Leading_Jet_eta;
    TH1F* h_Matched_Leading_Jet_mass;
    TH1F* h_Matched_Leading_Jet_phi;
    TH1F* h_Matched_Leading_Jet_pt;

    TH1F* h_Matched_Trailing_Jet_eta;
    TH1F* h_Matched_Trailing_Jet_mass;
    TH1F* h_Matched_Trailing_Jet_phi;
    TH1F* h_Matched_Trailing_Jet_pt;

    TH1F* h_Matched_Jet_invariantMass;
    TH1F* h_Matched_Jet_etaSeparation;

    TH1F* h_LHE_Leading_Jet_eta;
    TH1F* h_LHE_Leading_Jet_mass;
    TH1F* h_LHE_Leading_Jet_phi;
    TH1F* h_LHE_Leading_Jet_pt;

    TH1F* h_LHE_Trailing_Jet_eta;
    TH1F* h_LHE_Trailing_Jet_mass;
    TH1F* h_LHE_Trailing_Jet_phi;
    TH1F* h_LHE_Trailing_Jet_pt;

    TH1F* h_LHE_Jet_invariantMass;
    TH1F* h_LHE_Jet_etaSeparation;

    //TH1F* h_LHE_nJets;


/* Tree Values */
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         HTXS_Higgs_pt;
   Float_t         HTXS_Higgs_y;
   Int_t           HTXS_stage1_1_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_cat_pTjet30GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet30GeV;
   Int_t           HTXS_stage_0;
   Int_t           HTXS_stage_1_pTjet25;
   Int_t           HTXS_stage_1_pTjet30;
   UChar_t         HTXS_njets25;
   UChar_t         HTXS_njets30;
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nCorrT1METJet;
   Float_t         CorrT1METJet_area[21];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[21];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[21];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[21];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[21];   //[nCorrT1METJet]
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[6];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[6];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[6];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[6];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[6];   //[nElectron]
   Float_t         Electron_dxy[6];   //[nElectron]
   Float_t         Electron_dxyErr[6];   //[nElectron]
   Float_t         Electron_dz[6];   //[nElectron]
   Float_t         Electron_dzErr[6];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[6];   //[nElectron]
   Float_t         Electron_energyErr[6];   //[nElectron]
   Float_t         Electron_eta[6];   //[nElectron]
   Float_t         Electron_hoe[6];   //[nElectron]
   Float_t         Electron_ip3d[6];   //[nElectron]
   Float_t         Electron_jetPtRelv2[6];   //[nElectron]
   Float_t         Electron_jetRelIso[6];   //[nElectron]
   Float_t         Electron_mass[6];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[6];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[6];   //[nElectron]
   Float_t         Electron_mvaFall17V1Iso[6];   //[nElectron]
   Float_t         Electron_mvaFall17V1noIso[6];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[6];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[6];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[6];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[6];   //[nElectron]
   Float_t         Electron_phi[6];   //[nElectron]
   Float_t         Electron_pt[6];   //[nElectron]
   Float_t         Electron_r9[6];   //[nElectron]
   Float_t         Electron_sieie[6];   //[nElectron]
   Float_t         Electron_sip3d[6];   //[nElectron]
   Float_t         Electron_mvaTTH[6];   //[nElectron]
   Int_t           Electron_charge[6];   //[nElectron]
   Int_t           Electron_cutBased[6];   //[nElectron]
   Int_t           Electron_cutBased_Fall17_V1[6];   //[nElectron]
   Int_t           Electron_jetIdx[6];   //[nElectron]
   Int_t           Electron_pdgId[6];   //[nElectron]
   Int_t           Electron_photonIdx[6];   //[nElectron]
   Int_t           Electron_tightCharge[6];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[6];   //[nElectron]
   Bool_t          Electron_convVeto[6];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[6];   //[nElectron]
   Bool_t          Electron_isPFcand[6];   //[nElectron]
   UChar_t         Electron_lostHits[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP80[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP90[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WPL[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP80[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP90[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WPL[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[6];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[6];   //[nElectron]
   UChar_t         Electron_seedGain[6];   //[nElectron]

   UInt_t          nFatJet;
   Float_t         FatJet_area[6];   //[nFatJet]
   Float_t         FatJet_btagCMVA[6];   //[nFatJet]
   Float_t         FatJet_btagCSVV2[6];   //[nFatJet]
   Float_t         FatJet_btagDDBvL[6];   //[nFatJet]
   Float_t         FatJet_btagDDCvB[6];   //[nFatJet]
   Float_t         FatJet_btagDDCvL[6];   //[nFatJet]
   Float_t         FatJet_btagDeepB[6];   //[nFatJet]
   Float_t         FatJet_btagHbb[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_H4qvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_HbbvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_TvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_WvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHccvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZbbvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_bbvsLight[6];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ccvsLight[6];   //[nFatJet]
   Float_t         FatJet_deepTag_H[6];   //[nFatJet]
   Float_t         FatJet_deepTag_QCD[6];   //[nFatJet]
   Float_t         FatJet_deepTag_QCDothers[6];   //[nFatJet]
   Float_t         FatJet_deepTag_TvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTag_WvsQCD[6];   //[nFatJet]
   Float_t         FatJet_deepTag_ZvsQCD[6];   //[nFatJet]
   Float_t         FatJet_eta[6];   //[nFatJet]
   Float_t         FatJet_mass[6];   //[nFatJet]
   Float_t         FatJet_msoftdrop[6];   //[nFatJet]
   Float_t         FatJet_n2b1[6];   //[nFatJet]
   Float_t         FatJet_n3b1[6];   //[nFatJet]
   Float_t         FatJet_phi[6];   //[nFatJet]
   Float_t         FatJet_pt[6];   //[nFatJet]
   Float_t         FatJet_rawFactor[6];   //[nFatJet]
   Float_t         FatJet_tau1[6];   //[nFatJet]
   Float_t         FatJet_tau2[6];   //[nFatJet]
   Float_t         FatJet_tau3[6];   //[nFatJet]
   Float_t         FatJet_tau4[6];   //[nFatJet]
   Int_t           FatJet_jetId[6];   //[nFatJet]
   Int_t           FatJet_subJetIdx1[6];   //[nFatJet]
   Int_t           FatJet_subJetIdx2[6];   //[nFatJet]

   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_pt[8];   //[nGenJetAK8]
   UInt_t          nGenJet;
   Float_t         GenJet_eta[21];   //[nGenJet]
   Float_t         GenJet_mass[21];   //[nGenJet]
   Float_t         GenJet_phi[21];   //[nGenJet]
   Float_t         GenJet_pt[21];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[155];   //[nGenPart]
   Float_t         GenPart_mass[155];   //[nGenPart]
   Float_t         GenPart_phi[155];   //[nGenPart]
   Float_t         GenPart_pt[155];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[155];   //[nGenPart]
   Int_t           GenPart_pdgId[155];   //[nGenPart]
   Int_t           GenPart_status[155];   //[nGenPart]
   Int_t           GenPart_statusFlags[155];   //[nGenPart]
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_mass[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_phi[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_pt[16];   //[nSubGenJetAK8]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   UInt_t          nGenVisTau;
   Float_t         GenVisTau_eta[3];   //[nGenVisTau]
   Float_t         GenVisTau_mass[3];   //[nGenVisTau]
   Float_t         GenVisTau_phi[3];   //[nGenVisTau]
   Float_t         GenVisTau_pt[3];   //[nGenVisTau]
   Int_t           GenVisTau_charge[3];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[3];   //[nGenVisTau]
   Int_t           GenVisTau_status[3];   //[nGenVisTau]
   Float_t         genWeight;
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nLHEPdfWeight;
   Float_t         LHEPdfWeight[1];   //[nLHEPdfWeight]
   UInt_t          nLHEReweightingWeight;
   Float_t         LHEReweightingWeight[1];   //[nLHEReweightingWeight]
   UInt_t          nLHEScaleWeight;
   Float_t         LHEScaleWeight[1];   //[nLHEScaleWeight]
   UInt_t          nPSWeight;
   Float_t         PSWeight[4];   //[nPSWeight]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[4];   //[nIsoTrack]
   Float_t         IsoTrack_dz[4];   //[nIsoTrack]
   Float_t         IsoTrack_eta[4];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[4];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[4];   //[nIsoTrack]
   Float_t         IsoTrack_phi[4];   //[nIsoTrack]
   Float_t         IsoTrack_pt[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[4];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[4];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[4];   //[nIsoTrack]

   UInt_t          nJet;
   Float_t         Jet_area[23];   //[nJet]
   Float_t         Jet_btagCMVA[23];   //[nJet]
   Float_t         Jet_btagCSVV2[23];   //[nJet]
   Float_t         Jet_btagDeepB[23];   //[nJet]
   Float_t         Jet_btagDeepC[23];   //[nJet]
   Float_t         Jet_btagDeepFlavB[23];   //[nJet]
   Float_t         Jet_btagDeepFlavC[23];   //[nJet]
   Float_t         Jet_chEmEF[23];   //[nJet]
   Float_t         Jet_chHEF[23];   //[nJet]
   Float_t         Jet_eta[23];   //[nJet]
   Float_t         Jet_jercCHF[23];   //[nJet]
   Float_t         Jet_jercCHPUF[23];   //[nJet]
   Float_t         Jet_mass[23];   //[nJet]
   Float_t         Jet_muEF[23];   //[nJet]
   Float_t         Jet_muonSubtrFactor[23];   //[nJet]
   Float_t         Jet_neEmEF[23];   //[nJet]
   Float_t         Jet_neHEF[23];   //[nJet]
   Float_t         Jet_phi[23];   //[nJet]
   Float_t         Jet_pt[23];   //[nJet]
   Float_t         Jet_qgl[23];   //[nJet]
   Float_t         Jet_rawFactor[23];   //[nJet]
   Float_t         Jet_bRegCorr[23];   //[nJet]
   Float_t         Jet_bRegRes[23];   //[nJet]
   Int_t           Jet_electronIdx1[23];   //[nJet]
   Int_t           Jet_electronIdx2[23];   //[nJet]
   Int_t           Jet_jetId[23];   //[nJet]
   Int_t           Jet_muonIdx1[23];   //[nJet]
   Int_t           Jet_muonIdx2[23];   //[nJet]
   Int_t           Jet_nConstituents[23];   //[nJet]
   Int_t           Jet_nElectrons[23];   //[nJet]
   Int_t           Jet_nMuons[23];   //[nJet]
   Int_t           Jet_puId[23];   //[nJet]

   Float_t         LHE_HT;
   Float_t         LHE_HTIncoming;
   Float_t         LHE_Vpt;
   UChar_t         LHE_Njets;
   UChar_t         LHE_Nb;
   UChar_t         LHE_Nc;
   UChar_t         LHE_Nuds;
   UChar_t         LHE_Nglu;
   UChar_t         LHE_NpNLO;
   UChar_t         LHE_NpLO;
   UInt_t          nLHEPart;
   Float_t         LHEPart_pt[5];   //[nLHEPart]
   Float_t         LHEPart_eta[5];   //[nLHEPart]
   Float_t         LHEPart_phi[5];   //[nLHEPart]
   Float_t         LHEPart_mass[5];   //[nLHEPart]
   Int_t           LHEPart_pdgId[5];   //[nLHEPart]
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   UInt_t          nMuon;
   Float_t         Muon_dxy[8];   //[nMuon]
   Float_t         Muon_dxyErr[8];   //[nMuon]
   Float_t         Muon_dz[8];   //[nMuon]
   Float_t         Muon_dzErr[8];   //[nMuon]
   Float_t         Muon_eta[8];   //[nMuon]
   Float_t         Muon_ip3d[8];   //[nMuon]
   Float_t         Muon_jetPtRelv2[8];   //[nMuon]
   Float_t         Muon_jetRelIso[8];   //[nMuon]
   Float_t         Muon_mass[8];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[8];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[8];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[8];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[8];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[8];   //[nMuon]
   Float_t         Muon_phi[8];   //[nMuon]
   Float_t         Muon_pt[8];   //[nMuon]
   Float_t         Muon_ptErr[8];   //[nMuon]
   Float_t         Muon_segmentComp[8];   //[nMuon]
   Float_t         Muon_sip3d[8];   //[nMuon]
   Float_t         Muon_softMva[8];   //[nMuon]
   Float_t         Muon_tkRelIso[8];   //[nMuon]
   Float_t         Muon_tunepRelPt[8];   //[nMuon]
   Float_t         Muon_mvaLowPt[8];   //[nMuon]
   Float_t         Muon_mvaTTH[8];   //[nMuon]
   Int_t           Muon_charge[8];   //[nMuon]
   Int_t           Muon_jetIdx[8];   //[nMuon]
   Int_t           Muon_nStations[8];   //[nMuon]
   Int_t           Muon_nTrackerLayers[8];   //[nMuon]
   Int_t           Muon_pdgId[8];   //[nMuon]
   Int_t           Muon_tightCharge[8];   //[nMuon]
   UChar_t         Muon_highPtId[8];   //[nMuon]
   Bool_t          Muon_inTimeMuon[8];   //[nMuon]
   Bool_t          Muon_isGlobal[8];   //[nMuon]
   Bool_t          Muon_isPFcand[8];   //[nMuon]
   Bool_t          Muon_isTracker[8];   //[nMuon]
   Bool_t          Muon_looseId[8];   //[nMuon]
   Bool_t          Muon_mediumId[8];   //[nMuon]
   Bool_t          Muon_mediumPromptId[8];   //[nMuon]
   UChar_t         Muon_miniIsoId[8];   //[nMuon]
   UChar_t         Muon_multiIsoId[8];   //[nMuon]
   UChar_t         Muon_mvaId[8];   //[nMuon]
   UChar_t         Muon_pfIsoId[8];   //[nMuon]
   UChar_t         Muon_puppiIsoId[8];   //[nMuon]
   Bool_t          Muon_softId[8];   //[nMuon]
   Bool_t          Muon_softMvaId[8];   //[nMuon]
   Bool_t          Muon_tightId[8];   //[nMuon]
   UChar_t         Muon_tkIsoId[8];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[8];   //[nMuon]
   UInt_t          nPhoton;
   Float_t         Photon_energyErr[8];   //[nPhoton]
   Float_t         Photon_eta[8];   //[nPhoton]
   Float_t         Photon_hoe[8];   //[nPhoton]
   Float_t         Photon_mass[8];   //[nPhoton]
   Float_t         Photon_mvaID[8];   //[nPhoton]
   Float_t         Photon_mvaIDV1[8];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all[8];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg[8];   //[nPhoton]
   Float_t         Photon_phi[8];   //[nPhoton]
   Float_t         Photon_pt[8];   //[nPhoton]
   Float_t         Photon_r9[8];   //[nPhoton]
   Float_t         Photon_sieie[8];   //[nPhoton]
   Int_t           Photon_charge[8];   //[nPhoton]
   Int_t           Photon_cutBasedBitmap[8];   //[nPhoton]
   Int_t           Photon_cutBasedV1Bitmap[8];   //[nPhoton]
   Int_t           Photon_electronIdx[8];   //[nPhoton]
   Int_t           Photon_jetIdx[8];   //[nPhoton]
   Int_t           Photon_pdgId[8];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[8];   //[nPhoton]
   Bool_t          Photon_electronVeto[8];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[8];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[8];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[8];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[8];   //[nPhoton]
   Bool_t          Photon_pixelSeed[8];   //[nPhoton]
   UChar_t         Photon_seedGain[8];   //[nPhoton]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[3];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[3];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[3];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[3];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[3];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[3];   //[nGenDressedLepton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCMVA[12];   //[nSubJet]
   Float_t         SubJet_btagCSVV2[12];   //[nSubJet]
   Float_t         SubJet_btagDeepB[12];   //[nSubJet]
   Float_t         SubJet_eta[12];   //[nSubJet]
   Float_t         SubJet_mass[12];   //[nSubJet]
   Float_t         SubJet_n2b1[12];   //[nSubJet]
   Float_t         SubJet_n3b1[12];   //[nSubJet]
   Float_t         SubJet_phi[12];   //[nSubJet]
   Float_t         SubJet_pt[12];   //[nSubJet]
   Float_t         SubJet_rawFactor[12];   //[nSubJet]
   Float_t         SubJet_tau1[12];   //[nSubJet]
   Float_t         SubJet_tau2[12];   //[nSubJet]
   Float_t         SubJet_tau3[12];   //[nSubJet]
   Float_t         SubJet_tau4[12];   //[nSubJet]
   UInt_t          nTau;
   Float_t         Tau_chargedIso[6];   //[nTau]
   Float_t         Tau_dxy[6];   //[nTau]
   Float_t         Tau_dz[6];   //[nTau]
   Float_t         Tau_eta[6];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[6];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[6];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[6];   //[nTau]
   Float_t         Tau_mass[6];   //[nTau]
   Float_t         Tau_neutralIso[6];   //[nTau]
   Float_t         Tau_phi[6];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[6];   //[nTau]
   Float_t         Tau_pt[6];   //[nTau]
   Float_t         Tau_puCorr[6];   //[nTau]
   Float_t         Tau_rawAntiEle[6];   //[nTau]
   Float_t         Tau_rawAntiEle2018[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2VSe[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2VSjet[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2VSmu[6];   //[nTau]
   Float_t         Tau_rawIso[6];   //[nTau]
   Float_t         Tau_rawIsodR03[6];   //[nTau]
   Float_t         Tau_rawMVAnewDM2017v2[6];   //[nTau]
   Float_t         Tau_rawMVAoldDM[6];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v1[6];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v2[6];   //[nTau]
   Float_t         Tau_rawMVAoldDMdR032017v2[6];   //[nTau]
   Int_t           Tau_charge[6];   //[nTau]
   Int_t           Tau_decayMode[6];   //[nTau]
   Int_t           Tau_jetIdx[6];   //[nTau]
   Int_t           Tau_rawAntiEleCat[6];   //[nTau]
   Int_t           Tau_rawAntiEleCat2018[6];   //[nTau]
   UChar_t         Tau_idAntiEle[6];   //[nTau]
   UChar_t         Tau_idAntiEle2018[6];   //[nTau]
   UChar_t         Tau_idAntiMu[6];   //[nTau]
   Bool_t          Tau_idDecayMode[6];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2VSe[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2VSjet[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2VSmu[6];   //[nTau]
   UChar_t         Tau_idMVAnewDM2017v2[6];   //[nTau]
   UChar_t         Tau_idMVAoldDM[6];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v1[6];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v2[6];   //[nTau]
   UChar_t         Tau_idMVAoldDMdR032017v2[6];   //[nTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[38];   //[nTrigObj]
   Float_t         TrigObj_eta[38];   //[nTrigObj]
   Float_t         TrigObj_phi[38];   //[nTrigObj]
   Float_t         TrigObj_l1pt[38];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[38];   //[nTrigObj]
   Float_t         TrigObj_l2pt[38];   //[nTrigObj]
   Int_t           TrigObj_id[38];   //[nTrigObj]
   Int_t           TrigObj_l1iso[38];   //[nTrigObj]
   Int_t           TrigObj_l1charge[38];   //[nTrigObj]
   Int_t           TrigObj_filterBits[38];   //[nTrigObj]
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[15];   //[nSV]
   Float_t         SV_dlenSig[15];   //[nSV]
   Float_t         SV_pAngle[15];   //[nSV]
   Int_t           Electron_genPartIdx[6];   //[nElectron]
   UChar_t         Electron_genPartFlav[6];   //[nElectron]
   Int_t           GenJetAK8_partonFlavour[8];   //[nGenJetAK8]
   UChar_t         GenJetAK8_hadronFlavour[8];   //[nGenJetAK8]
   Int_t           GenJet_partonFlavour[21];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[21];   //[nGenJet]
   Int_t           Jet_genJetIdx[23];   //[nJet]
   Int_t           Jet_hadronFlavour[23];   //[nJet]
   Int_t           Jet_partonFlavour[23];   //[nJet]
   Int_t           Muon_genPartIdx[8];   //[nMuon]
   UChar_t         Muon_genPartFlav[8];   //[nMuon]
   Int_t           Photon_genPartIdx[8];   //[nPhoton]
   UChar_t         Photon_genPartFlav[8];   //[nPhoton]
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[6];   //[nElectron]
   UChar_t         Jet_cleanmask[23];   //[nJet]
   UChar_t         Muon_cleanmask[8];   //[nMuon]
   UChar_t         Photon_cleanmask[8];   //[nPhoton]
   UChar_t         Tau_cleanmask[6];   //[nTau]
   Float_t         SV_chi2[15];   //[nSV]
   Float_t         SV_eta[15];   //[nSV]
   Float_t         SV_mass[15];   //[nSV]
   Float_t         SV_ndof[15];   //[nSV]
   Float_t         SV_phi[15];   //[nSV]
   Float_t         SV_pt[15];   //[nSV]
   Float_t         SV_x[15];   //[nSV]
   Float_t         SV_y[15];   //[nSV]
   Float_t         SV_z[15];   //[nSV]
   Int_t           Tau_genPartIdx[6];   //[nTau]
   UChar_t         Tau_genPartFlav[6];   //[nTau]
   Bool_t          L1simulation_step;
   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_AK8PFJet360_TrimMass30;
   Bool_t          HLT_AK8PFJet380_TrimMass30;
   Bool_t          HLT_AK8PFJet400_TrimMass30;
   Bool_t          HLT_AK8PFJet420_TrimMass30;
   Bool_t          HLT_AK8PFHT750_TrimMass50;
   Bool_t          HLT_AK8PFHT800_TrimMass50;
   Bool_t          HLT_AK8PFHT850_TrimMass50;
   Bool_t          HLT_AK8PFHT900_TrimMass50;
   Bool_t          HLT_CaloJet500_NoJetID;
   Bool_t          HLT_CaloJet550_NoJetID;
   Bool_t          HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
   Bool_t          HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Ele27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_DoubleMu4_3_Bs;
   Bool_t          HLT_DoubleMu4_3_Jpsi;
   Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
   Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          HLT_DoubleMu3_TkMu_DsTau3Mu;
   Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Bool_t          HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   Bool_t          HLT_Mu3_PFJet40;
   Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          HLT_Mu7p5_Track2_Jpsi;
   Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
   Bool_t          HLT_Mu7p5_Track7_Jpsi;
   Bool_t          HLT_Mu7p5_Track2_Upsilon;
   Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
   Bool_t          HLT_Mu7p5_Track7_Upsilon;
   Bool_t          HLT_Mu3_L1SingleMu5orSingleMu7;
   Bool_t          HLT_DoublePhoton33_CaloIdL;
   Bool_t          HLT_DoublePhoton70;
   Bool_t          HLT_DoublePhoton85;
   Bool_t          HLT_Ele20_WPTight_Gsf;
   Bool_t          HLT_Ele15_WPLoose_Gsf;
   Bool_t          HLT_Ele17_WPLoose_Gsf;
   Bool_t          HLT_Ele20_WPLoose_Gsf;
   Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf;
   Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_Ele28_WPTight_Gsf;
   Bool_t          HLT_Ele30_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_HT450_Beamspot;
   Bool_t          HLT_HT300_Beamspot;
   Bool_t          HLT_ZeroBias_Beamspot;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1;
   Bool_t          HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1;
   Bool_t          HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu30;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t          HLT_L1SingleMu18;
   Bool_t          HLT_L1SingleMu25;
   Bool_t          HLT_L2Mu10;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu50;
   Bool_t          HLT_L2Mu23NoVtx_2Cha;
   Bool_t          HLT_L2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;
   Bool_t          HLT_DoubleL2Mu50;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu25_TkMu0_Onia;
   Bool_t          HLT_Mu30_TkMu0_Psi;
   Bool_t          HLT_Mu30_TkMu0_Upsilon;
   Bool_t          HLT_Mu20_TkMu0_Phi;
   Bool_t          HLT_Mu25_TkMu0_Phi;
   Bool_t          HLT_Mu12;
   Bool_t          HLT_Mu15;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_OldMu100;
   Bool_t          HLT_TkMu100;
   Bool_t          HLT_DiPFJetAve40;
   Bool_t          HLT_DiPFJetAve60;
   Bool_t          HLT_DiPFJetAve80;
   Bool_t          HLT_DiPFJetAve140;
   Bool_t          HLT_DiPFJetAve200;
   Bool_t          HLT_DiPFJetAve260;
   Bool_t          HLT_DiPFJetAve320;
   Bool_t          HLT_DiPFJetAve400;
   Bool_t          HLT_DiPFJetAve500;
   Bool_t          HLT_DiPFJetAve60_HFJEC;
   Bool_t          HLT_DiPFJetAve80_HFJEC;
   Bool_t          HLT_DiPFJetAve100_HFJEC;
   Bool_t          HLT_DiPFJetAve160_HFJEC;
   Bool_t          HLT_DiPFJetAve220_HFJEC;
   Bool_t          HLT_DiPFJetAve300_HFJEC;
   Bool_t          HLT_AK8PFJet15;
   Bool_t          HLT_AK8PFJet25;
   Bool_t          HLT_AK8PFJet40;
   Bool_t          HLT_AK8PFJet60;
   Bool_t          HLT_AK8PFJet80;
   Bool_t          HLT_AK8PFJet140;
   Bool_t          HLT_AK8PFJet200;
   Bool_t          HLT_AK8PFJet260;
   Bool_t          HLT_AK8PFJet320;
   Bool_t          HLT_AK8PFJet400;
   Bool_t          HLT_AK8PFJet450;
   Bool_t          HLT_AK8PFJet500;
   Bool_t          HLT_AK8PFJet550;
   Bool_t          HLT_PFJet15;
   Bool_t          HLT_PFJet25;
   Bool_t          HLT_PFJet40;
   Bool_t          HLT_PFJet60;
   Bool_t          HLT_PFJet80;
   Bool_t          HLT_PFJet140;
   Bool_t          HLT_PFJet200;
   Bool_t          HLT_PFJet260;
   Bool_t          HLT_PFJet320;
   Bool_t          HLT_PFJet400;
   Bool_t          HLT_PFJet450;
   Bool_t          HLT_PFJet500;
   Bool_t          HLT_PFJet550;
   Bool_t          HLT_PFJetFwd15;
   Bool_t          HLT_PFJetFwd25;
   Bool_t          HLT_PFJetFwd40;
   Bool_t          HLT_PFJetFwd60;
   Bool_t          HLT_PFJetFwd80;
   Bool_t          HLT_PFJetFwd140;
   Bool_t          HLT_PFJetFwd200;
   Bool_t          HLT_PFJetFwd260;
   Bool_t          HLT_PFJetFwd320;
   Bool_t          HLT_PFJetFwd400;
   Bool_t          HLT_PFJetFwd450;
   Bool_t          HLT_PFJetFwd500;
   Bool_t          HLT_AK8PFJetFwd15;
   Bool_t          HLT_AK8PFJetFwd25;
   Bool_t          HLT_AK8PFJetFwd40;
   Bool_t          HLT_AK8PFJetFwd60;
   Bool_t          HLT_AK8PFJetFwd80;
   Bool_t          HLT_AK8PFJetFwd140;
   Bool_t          HLT_AK8PFJetFwd200;
   Bool_t          HLT_AK8PFJetFwd260;
   Bool_t          HLT_AK8PFJetFwd320;
   Bool_t          HLT_AK8PFJetFwd400;
   Bool_t          HLT_AK8PFJetFwd450;
   Bool_t          HLT_AK8PFJetFwd500;
   Bool_t          HLT_PFHT180;
   Bool_t          HLT_PFHT250;
   Bool_t          HLT_PFHT370;
   Bool_t          HLT_PFHT430;
   Bool_t          HLT_PFHT510;
   Bool_t          HLT_PFHT590;
   Bool_t          HLT_PFHT680;
   Bool_t          HLT_PFHT780;
   Bool_t          HLT_PFHT890;
   Bool_t          HLT_PFHT1050;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFHT700_PFMET95_PFMHT95_IDTight;
   Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t          HLT_PFHT800_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne110_PFMHT110_IDTight;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight;
   Bool_t          HLT_PFMETTypeOne130_PFMHT130_IDTight;
   Bool_t          HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_L1ETMHadSeeds;
   Bool_t          HLT_CaloMHT90;
   Bool_t          HLT_CaloMET80_NotCleaned;
   Bool_t          HLT_CaloMET90_NotCleaned;
   Bool_t          HLT_CaloMET100_NotCleaned;
   Bool_t          HLT_CaloMET110_NotCleaned;
   Bool_t          HLT_CaloMET250_NotCleaned;
   Bool_t          HLT_CaloMET70_HBHECleaned;
   Bool_t          HLT_CaloMET80_HBHECleaned;
   Bool_t          HLT_CaloMET90_HBHECleaned;
   Bool_t          HLT_CaloMET100_HBHECleaned;
   Bool_t          HLT_CaloMET250_HBHECleaned;
   Bool_t          HLT_CaloMET300_HBHECleaned;
   Bool_t          HLT_CaloMET350_HBHECleaned;
   Bool_t          HLT_PFMET200_NotCleaned;
   Bool_t          HLT_PFMET200_HBHECleaned;
   Bool_t          HLT_PFMET250_HBHECleaned;
   Bool_t          HLT_PFMET300_HBHECleaned;
   Bool_t          HLT_PFMET200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_MET105_IsoTrk50;
   Bool_t          HLT_MET120_IsoTrk50;
   Bool_t          HLT_SingleJet30_Mu12_SinglePFJet40;
   Bool_t          HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets40_CaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets100_CaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets200_CaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets350_CaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          HLT_Photon300_NoHE;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t          HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK8Jet170_DoubleMu5;
   Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet20_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK4DiJet40_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK4DiJet70_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK4DiJet110_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK4DiJet170_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK4Jet300_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK8DiJet170_Mu5_noalgo;
   Bool_t          HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo;
   Bool_t          HLT_BTagMu_AK8Jet300_Mu5_noalgo;
   Bool_t          HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu12_DoublePhoton20;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;
   Bool_t          HLT_Photon20;
   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon100EB_TightID_TightIso;
   Bool_t          HLT_Photon110EB_TightID_TightIso;
   Bool_t          HLT_Photon120EB_TightID_TightIso;
   Bool_t          HLT_Photon100EBHE10;
   Bool_t          HLT_Photon100EEHE10;
   Bool_t          HLT_Photon100EE_TightID_TightIso;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_CaloIdL_PFHT700;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLT_Photon35_TwoProngs35;
   Bool_t          HLT_IsoMu24_TwoProngs35;
   Bool_t          HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t          HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5NoOS;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t          HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5M;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5R;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t          HLT_Dimuon0_LowMass;
   Bool_t          HLT_Dimuon0_LowMass_L1_4;
   Bool_t          HLT_Dimuon0_LowMass_L1_4R;
   Bool_t          HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_L1_TM0;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8_DZ;
   Bool_t          HLT_TripleMu_10_5_5_DZ;
   Bool_t          HLT_TripleMu_12_10_5;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t          HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t          HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t          HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t          HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t          HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t          HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t          HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t          HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t          HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;
   Bool_t          HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;
   Bool_t          HLT_DoubleMu33NoFiltersNoVtxDisplaced;
   Bool_t          HLT_DoubleMu40NoFiltersNoVtxDisplaced;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;
   Bool_t          HLT_HT425;
   Bool_t          HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT500_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet60_DisplacedTrack;
   Bool_t          HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET110;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET120;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET130;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET110;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET120;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET130;
   Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t          HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_Photon23;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;
   Bool_t          HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Bool_t          HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon12_Upsilon_y1p4;
   Bool_t          HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon18_PsiPrime;
   Bool_t          HLT_Dimuon25_Jpsi;
   Bool_t          HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t          HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t          HLT_Dimuon24_Phi_noCorrL1;
   Bool_t          HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t          HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleIsoMu20_eta2p1;
   Bool_t          HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t          HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;
   Bool_t          HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele145_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele200_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele250_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele300_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40;
   Bool_t          HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94;
   Bool_t          HLT_PFHT400_SixPFJet32;
   Bool_t          HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59;
   Bool_t          HLT_PFHT450_SixPFJet36;
   Bool_t          HLT_PFHT350;
   Bool_t          HLT_PFHT350MinPFJet15;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
   Bool_t          HLT_ECALHT800;
   Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          HLT_Physics;
   Bool_t          HLT_Physics_part0;
   Bool_t          HLT_Physics_part1;
   Bool_t          HLT_Physics_part2;
   Bool_t          HLT_Physics_part3;
   Bool_t          HLT_Physics_part4;
   Bool_t          HLT_Physics_part5;
   Bool_t          HLT_Physics_part6;
   Bool_t          HLT_Physics_part7;
   Bool_t          HLT_Random;
   Bool_t          HLT_ZeroBias;
   Bool_t          HLT_ZeroBias_Alignment;
   Bool_t          HLT_ZeroBias_part0;
   Bool_t          HLT_ZeroBias_part1;
   Bool_t          HLT_ZeroBias_part2;
   Bool_t          HLT_ZeroBias_part3;
   Bool_t          HLT_ZeroBias_part4;
   Bool_t          HLT_ZeroBias_part5;
   Bool_t          HLT_ZeroBias_part6;
   Bool_t          HLT_ZeroBias_part7;
   Bool_t          HLT_AK4CaloJet30;
   Bool_t          HLT_AK4CaloJet40;
   Bool_t          HLT_AK4CaloJet50;
   Bool_t          HLT_AK4CaloJet80;
   Bool_t          HLT_AK4CaloJet100;
   Bool_t          HLT_AK4CaloJet120;
   Bool_t          HLT_AK4PFJet30;
   Bool_t          HLT_AK4PFJet50;
   Bool_t          HLT_AK4PFJet80;
   Bool_t          HLT_AK4PFJet100;
   Bool_t          HLT_AK4PFJet120;
   Bool_t          HLT_SinglePhoton10_Eta3p1ForPPRef;
   Bool_t          HLT_SinglePhoton20_Eta3p1ForPPRef;
   Bool_t          HLT_SinglePhoton30_Eta3p1ForPPRef;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   Bool_t          HLT_EcalCalibration;
   Bool_t          HLT_HcalCalibration;
   Bool_t          HLT_L1UnpairedBunchBptxMinus;
   Bool_t          HLT_L1UnpairedBunchBptxPlus;
   Bool_t          HLT_L1NotBptxOR;
   Bool_t          HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          HLT_CDC_L2cosmic_5_er1p0;
   Bool_t          HLT_CDC_L2cosmic_5p5_er1p0;
   Bool_t          HLT_HcalNZS;
   Bool_t          HLT_HcalPhiSym;
   Bool_t          HLT_HcalIsolatedbunch;
   Bool_t          HLT_IsoTrackHB;
   Bool_t          HLT_IsoTrackHE;
   Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          HLT_ZeroBias_IsolatedBunches;
   Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t          HLT_Rsq0p35;
   Bool_t          HLT_Rsq0p40;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200_4jet;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200_4jet;
   Bool_t          HLT_IsoMu27_MET90;
   Bool_t          HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1;
   Bool_t          HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1;
   Bool_t          HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_Mu18_Mu9_SameSign;
   Bool_t          HLT_Mu18_Mu9_SameSign_DZ;
   Bool_t          HLT_Mu18_Mu9;
   Bool_t          HLT_Mu18_Mu9_DZ;
   Bool_t          HLT_Mu20_Mu10_SameSign;
   Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
   Bool_t          HLT_Mu20_Mu10;
   Bool_t          HLT_Mu20_Mu10_DZ;
   Bool_t          HLT_Mu23_Mu12_SameSign;
   Bool_t          HLT_Mu23_Mu12_SameSign_DZ;
   Bool_t          HLT_Mu23_Mu12;
   Bool_t          HLT_Mu23_Mu12_DZ;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8_DCA;
   Bool_t          HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          HLT_QuadPFJet98_83_71_15;
   Bool_t          HLT_QuadPFJet103_88_75_15;
   Bool_t          HLT_QuadPFJet105_88_76_15;
   Bool_t          HLT_QuadPFJet111_90_80_15;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto;
   Bool_t          HLT_Mu12_IP6_part0;
   Bool_t          HLT_Mu12_IP6_part1;
   Bool_t          HLT_Mu12_IP6_part2;
   Bool_t          HLT_Mu12_IP6_part3;
   Bool_t          HLT_Mu12_IP6_part4;
   Bool_t          HLT_Mu9_IP5_part0;
   Bool_t          HLT_Mu9_IP5_part1;
   Bool_t          HLT_Mu9_IP5_part2;
   Bool_t          HLT_Mu9_IP5_part3;
   Bool_t          HLT_Mu9_IP5_part4;
   Bool_t          HLT_Mu7_IP4_part0;
   Bool_t          HLT_Mu7_IP4_part1;
   Bool_t          HLT_Mu7_IP4_part2;
   Bool_t          HLT_Mu7_IP4_part3;
   Bool_t          HLT_Mu7_IP4_part4;
   Bool_t          HLT_Mu9_IP4_part0;
   Bool_t          HLT_Mu9_IP4_part1;
   Bool_t          HLT_Mu9_IP4_part2;
   Bool_t          HLT_Mu9_IP4_part3;
   Bool_t          HLT_Mu9_IP4_part4;
   Bool_t          HLT_Mu8_IP5_part0;
   Bool_t          HLT_Mu8_IP5_part1;
   Bool_t          HLT_Mu8_IP5_part2;
   Bool_t          HLT_Mu8_IP5_part3;
   Bool_t          HLT_Mu8_IP5_part4;
   Bool_t          HLT_Mu8_IP6_part0;
   Bool_t          HLT_Mu8_IP6_part1;
   Bool_t          HLT_Mu8_IP6_part2;
   Bool_t          HLT_Mu8_IP6_part3;
   Bool_t          HLT_Mu8_IP6_part4;
   Bool_t          HLT_Mu9_IP6_part0;
   Bool_t          HLT_Mu9_IP6_part1;
   Bool_t          HLT_Mu9_IP6_part2;
   Bool_t          HLT_Mu9_IP6_part3;
   Bool_t          HLT_Mu9_IP6_part4;
   Bool_t          HLT_Mu8_IP3_part0;
   Bool_t          HLT_Mu8_IP3_part1;
   Bool_t          HLT_Mu8_IP3_part2;
   Bool_t          HLT_Mu8_IP3_part3;
   Bool_t          HLT_Mu8_IP3_part4;
   Bool_t          HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          HLT_TrkMu6NoFiltersNoVtx;
   Bool_t          HLT_TrkMu16NoFiltersNoVtx;
   Bool_t          HLT_DoubleTrkMu_16_6_NoFiltersNoVtx;
   Bool_t          HLTriggerFinalPath;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          L1Reco_step;
   Bool_t          L1_AlwaysTrue;
   Bool_t          L1_BPTX_AND_Ref1_VME;
   Bool_t          L1_BPTX_AND_Ref3_VME;
   Bool_t          L1_BPTX_AND_Ref4_VME;
   Bool_t          L1_BPTX_BeamGas_B1_VME;
   Bool_t          L1_BPTX_BeamGas_B2_VME;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          L1_BPTX_NotOR_VME;
   Bool_t          L1_BPTX_OR_Ref3_VME;
   Bool_t          L1_BPTX_OR_Ref4_VME;
   Bool_t          L1_BPTX_RefAND_VME;
   Bool_t          L1_BptxMinus;
   Bool_t          L1_BptxOR;
   Bool_t          L1_BptxPlus;
   Bool_t          L1_BptxXOR;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          L1_DoubleEG8er2p5_HTT260er;
   Bool_t          L1_DoubleEG8er2p5_HTT280er;
   Bool_t          L1_DoubleEG8er2p5_HTT300er;
   Bool_t          L1_DoubleEG8er2p5_HTT320er;
   Bool_t          L1_DoubleEG8er2p5_HTT340er;
   Bool_t          L1_DoubleEG_15_10_er2p5;
   Bool_t          L1_DoubleEG_20_10_er2p5;
   Bool_t          L1_DoubleEG_22_10_er2p5;
   Bool_t          L1_DoubleEG_25_12_er2p5;
   Bool_t          L1_DoubleEG_25_14_er2p5;
   Bool_t          L1_DoubleEG_27_14_er2p5;
   Bool_t          L1_DoubleEG_LooseIso20_10_er2p5;
   Bool_t          L1_DoubleEG_LooseIso22_10_er2p5;
   Bool_t          L1_DoubleEG_LooseIso22_12_er2p5;
   Bool_t          L1_DoubleEG_LooseIso25_12_er2p5;
   Bool_t          L1_DoubleIsoTau32er2p1;
   Bool_t          L1_DoubleIsoTau34er2p1;
   Bool_t          L1_DoubleIsoTau36er2p1;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet100er2p5;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet120er2p5;
   Bool_t          L1_DoubleJet150er2p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5;
   Bool_t          L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp;
   Bool_t          L1_DoubleJet40er2p5;
   Bool_t          L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28;
   Bool_t          L1_DoubleJet_120_45_DoubleJet45_Mass_Min620;
   Bool_t          L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28;
   Bool_t          L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ;
   Bool_t          L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp;
   Bool_t          L1_DoubleJet_80_30_Mass_Min420_Mu8;
   Bool_t          L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleLooseIsoEG22er2p1;
   Bool_t          L1_DoubleLooseIsoEG24er2p1;
   Bool_t          L1_DoubleMu0;
   Bool_t          L1_DoubleMu0_Mass_Min1;
   Bool_t          L1_DoubleMu0_OQ;
   Bool_t          L1_DoubleMu0_SQ;
   Bool_t          L1_DoubleMu0_SQ_OS;
   Bool_t          L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2p0_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu10_SQ;
   Bool_t          L1_DoubleMu18er2p1;
   Bool_t          L1_DoubleMu3_OS_DoubleEG7p5Upsilon;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_HTT60er;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;
   Bool_t          L1_DoubleMu3_SQ_HTT220er;
   Bool_t          L1_DoubleMu3_SQ_HTT240er;
   Bool_t          L1_DoubleMu3_SQ_HTT260er;
   Bool_t          L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t          L1_DoubleMu4_SQ_EG9er2p5;
   Bool_t          L1_DoubleMu4_SQ_OS;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5_SQ_OS;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu5Upsilon_OS_DoubleEG3;
   Bool_t          L1_DoubleMu5_SQ_EG9er2p5;
   Bool_t          L1_DoubleMu9_SQ;
   Bool_t          L1_DoubleMu_12_5;
   Bool_t          L1_DoubleMu_15_5_SQ;
   Bool_t          L1_DoubleMu_15_7;
   Bool_t          L1_DoubleMu_15_7_Mass_Min1;
   Bool_t          L1_DoubleMu_15_7_SQ;
   Bool_t          L1_DoubleTau70er2p1;
   Bool_t          L1_ETM120;
   Bool_t          L1_ETM150;
   Bool_t          L1_ETMHF100;
   Bool_t          L1_ETMHF100_HTT60er;
   Bool_t          L1_ETMHF110;
   Bool_t          L1_ETMHF110_HTT60er;
   Bool_t          L1_ETMHF110_HTT60er_NotSecondBunchInTrain;
   Bool_t          L1_ETMHF120;
   Bool_t          L1_ETMHF120_HTT60er;
   Bool_t          L1_ETMHF120_NotSecondBunchInTrain;
   Bool_t          L1_ETMHF130;
   Bool_t          L1_ETMHF130_HTT60er;
   Bool_t          L1_ETMHF140;
   Bool_t          L1_ETMHF150;
   Bool_t          L1_ETMHF90_HTT60er;
   Bool_t          L1_ETT1200;
   Bool_t          L1_ETT1600;
   Bool_t          L1_ETT2000;
   Bool_t          L1_FirstBunchAfterTrain;
   Bool_t          L1_FirstBunchBeforeTrain;
   Bool_t          L1_FirstBunchInTrain;
   Bool_t          L1_FirstCollisionInOrbit;
   Bool_t          L1_FirstCollisionInTrain;
   Bool_t          L1_HCAL_LaserMon_Trig;
   Bool_t          L1_HCAL_LaserMon_Veto;
   Bool_t          L1_HTT120er;
   Bool_t          L1_HTT160er;
   Bool_t          L1_HTT200er;
   Bool_t          L1_HTT255er;
   Bool_t          L1_HTT280er;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p4;
   Bool_t          L1_HTT320er;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p4;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;
   Bool_t          L1_HTT360er;
   Bool_t          L1_HTT400er;
   Bool_t          L1_HTT450er;
   Bool_t          L1_IsoEG32er2p5_Mt40;
   Bool_t          L1_IsoEG32er2p5_Mt44;
   Bool_t          L1_IsoEG32er2p5_Mt48;
   Bool_t          L1_IsoTau40er2p1_ETMHF100;
   Bool_t          L1_IsoTau40er2p1_ETMHF110;
   Bool_t          L1_IsoTau40er2p1_ETMHF120;
   Bool_t          L1_IsoTau40er2p1_ETMHF90;
   Bool_t          L1_IsolatedBunch;
   Bool_t          L1_LastBunchInTrain;
   Bool_t          L1_LastCollisionInTrain;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_LooseIsoEG30er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          L1_Mu18er2p1_Tau24er2p1;
   Bool_t          L1_Mu18er2p1_Tau26er2p1;
   Bool_t          L1_Mu20_EG10er2p5;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          L1_Mu22er2p1_Tau70er2p1;
   Bool_t          L1_Mu3_Jet120er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet120er2p5_dR_Max0p8;
   Bool_t          L1_Mu3_Jet16er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet30er2p5;
   Bool_t          L1_Mu3_Jet35er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet60er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet80er2p5_dR_Max0p4;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF40;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF50;
   Bool_t          L1_Mu5_EG23er2p5;
   Bool_t          L1_Mu5_LooseIsoEG20er2p5;
   Bool_t          L1_Mu6_DoubleEG10er2p5;
   Bool_t          L1_Mu6_DoubleEG12er2p5;
   Bool_t          L1_Mu6_DoubleEG15er2p5;
   Bool_t          L1_Mu6_DoubleEG17er2p5;
   Bool_t          L1_Mu6_HTT240er;
   Bool_t          L1_Mu6_HTT250er;
   Bool_t          L1_Mu7_EG23er2p5;
   Bool_t          L1_Mu7_LooseIsoEG20er2p5;
   Bool_t          L1_Mu7_LooseIsoEG23er2p5;
   Bool_t          L1_NotBptxOR;
   Bool_t          L1_QuadJet36er2p5_IsoTau52er2p1;
   Bool_t          L1_QuadJet60er2p5;
   Bool_t          L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;
   Bool_t          L1_QuadMu0;
   Bool_t          L1_QuadMu0_OQ;
   Bool_t          L1_QuadMu0_SQ;
   Bool_t          L1_SecondBunchInTrain;
   Bool_t          L1_SecondLastBunchInTrain;
   Bool_t          L1_SingleEG10er2p5;
   Bool_t          L1_SingleEG15er2p5;
   Bool_t          L1_SingleEG26er2p5;
   Bool_t          L1_SingleEG34er2p5;
   Bool_t          L1_SingleEG36er2p5;
   Bool_t          L1_SingleEG38er2p5;
   Bool_t          L1_SingleEG40er2p5;
   Bool_t          L1_SingleEG42er2p5;
   Bool_t          L1_SingleEG45er2p5;
   Bool_t          L1_SingleEG50;
   Bool_t          L1_SingleEG60;
   Bool_t          L1_SingleEG8er2p5;
   Bool_t          L1_SingleIsoEG24er1p5;
   Bool_t          L1_SingleIsoEG24er2p1;
   Bool_t          L1_SingleIsoEG26er1p5;
   Bool_t          L1_SingleIsoEG26er2p1;
   Bool_t          L1_SingleIsoEG26er2p5;
   Bool_t          L1_SingleIsoEG28er1p5;
   Bool_t          L1_SingleIsoEG28er2p1;
   Bool_t          L1_SingleIsoEG28er2p5;
   Bool_t          L1_SingleIsoEG30er2p1;
   Bool_t          L1_SingleIsoEG30er2p5;
   Bool_t          L1_SingleIsoEG32er2p1;
   Bool_t          L1_SingleIsoEG32er2p5;
   Bool_t          L1_SingleIsoEG34er2p5;
   Bool_t          L1_SingleJet10erHE;
   Bool_t          L1_SingleJet120;
   Bool_t          L1_SingleJet120_FWD3p0;
   Bool_t          L1_SingleJet120er2p5;
   Bool_t          L1_SingleJet12erHE;
   Bool_t          L1_SingleJet140er2p5;
   Bool_t          L1_SingleJet140er2p5_ETMHF80;
   Bool_t          L1_SingleJet140er2p5_ETMHF90;
   Bool_t          L1_SingleJet160er2p5;
   Bool_t          L1_SingleJet180;
   Bool_t          L1_SingleJet180er2p5;
   Bool_t          L1_SingleJet200;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet35;
   Bool_t          L1_SingleJet35_FWD3p0;
   Bool_t          L1_SingleJet35er2p5;
   Bool_t          L1_SingleJet43er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet46er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet60;
   Bool_t          L1_SingleJet60_FWD3p0;
   Bool_t          L1_SingleJet60er2p5;
   Bool_t          L1_SingleJet8erHE;
   Bool_t          L1_SingleJet90;
   Bool_t          L1_SingleJet90_FWD3p0;
   Bool_t          L1_SingleJet90er2p5;
   Bool_t          L1_SingleLooseIsoEG28er1p5;
   Bool_t          L1_SingleLooseIsoEG30er1p5;
   Bool_t          L1_SingleMu0_BMTF;
   Bool_t          L1_SingleMu0_DQ;
   Bool_t          L1_SingleMu0_EMTF;
   Bool_t          L1_SingleMu0_OMTF;
   Bool_t          L1_SingleMu10er1p5;
   Bool_t          L1_SingleMu12_DQ_BMTF;
   Bool_t          L1_SingleMu12_DQ_EMTF;
   Bool_t          L1_SingleMu12_DQ_OMTF;
   Bool_t          L1_SingleMu12er1p5;
   Bool_t          L1_SingleMu14er1p5;
   Bool_t          L1_SingleMu15_DQ;
   Bool_t          L1_SingleMu16er1p5;
   Bool_t          L1_SingleMu18;
   Bool_t          L1_SingleMu18er1p5;
   Bool_t          L1_SingleMu20;
   Bool_t          L1_SingleMu22;
   Bool_t          L1_SingleMu22_BMTF;
   Bool_t          L1_SingleMu22_EMTF;
   Bool_t          L1_SingleMu22_OMTF;
   Bool_t          L1_SingleMu25;
   Bool_t          L1_SingleMu3;
   Bool_t          L1_SingleMu5;
   Bool_t          L1_SingleMu6er1p5;
   Bool_t          L1_SingleMu7;
   Bool_t          L1_SingleMu7_DQ;
   Bool_t          L1_SingleMu7er1p5;
   Bool_t          L1_SingleMu8er1p5;
   Bool_t          L1_SingleMu9er1p5;
   Bool_t          L1_SingleMuCosmics;
   Bool_t          L1_SingleMuCosmics_BMTF;
   Bool_t          L1_SingleMuCosmics_EMTF;
   Bool_t          L1_SingleMuCosmics_OMTF;
   Bool_t          L1_SingleMuOpen;
   Bool_t          L1_SingleMuOpen_NotBptxOR;
   Bool_t          L1_SingleMuOpen_er1p1_NotBptxOR_3BX;
   Bool_t          L1_SingleMuOpen_er1p4_NotBptxOR_3BX;
   Bool_t          L1_SingleTau120er2p1;
   Bool_t          L1_SingleTau130er2p1;
   Bool_t          L1_TOTEM_1;
   Bool_t          L1_TOTEM_2;
   Bool_t          L1_TOTEM_3;
   Bool_t          L1_TOTEM_4;
   Bool_t          L1_TripleEG16er2p5;
   Bool_t          L1_TripleEG_16_12_8_er2p5;
   Bool_t          L1_TripleEG_16_15_8_er2p5;
   Bool_t          L1_TripleEG_18_17_8_er2p5;
   Bool_t          L1_TripleEG_18_18_12_er2p5;
   Bool_t          L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;
   Bool_t          L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;
   Bool_t          L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;
   Bool_t          L1_TripleMu0;
   Bool_t          L1_TripleMu0_OQ;
   Bool_t          L1_TripleMu0_SQ;
   Bool_t          L1_TripleMu3;
   Bool_t          L1_TripleMu3_SQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5_3_3;
   Bool_t          L1_TripleMu_5_3_3_SQ;
   Bool_t          L1_TripleMu_5_3p5_2p5;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_5_3;
   Bool_t          L1_UnpairedBunchBptxMinus;
   Bool_t          L1_UnpairedBunchBptxPlus;
   Bool_t          L1_ZeroBias;
   Bool_t          L1_ZeroBias_copy;
   Bool_t          L1_UnprefireableEvent;


    TTree* EventTree;
    Int_t EvMax;

string cur_time(){
//returns the current time, for logging
    std::time_t tt = std::time(NULL);
    std::string s = std::ctime(&tt);
    return s.substr(0, s.size()-1);
}

void setOutput(string path){
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1){
        if( errno == EEXIST ) {
        // already exists
            return;
        } else {
        // something else
            std::cout << "cannot create folder error:" << strerror(errno) << std::endl;
        }
    }
    return;
}

void swap(int &a, int &b){
    int c = a;
    a = b;
    b = c;
}

void EventLoop(){
//Loops over all events and fills histograms
    int eventLoopMax = EvMax;
    float jetPT_cut = 50.0;
    float jetInvMass_cut = 400.0;
    float jetEtaSep_cut = 2.0;

    for(int iev=0; iev<eventLoopMax;++iev){
        EventTree->GetEvent(iev);
        if(iev % 50 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;

        //plot Gen Level info
        h_nGenJets->Fill(nGenJet);
        for(int iJet=0; iJet<nGenJet; iJet++){
            h_Gen_Jet_eta ->Fill(GenJet_eta[iJet]);
            h_Gen_Jet_mass ->Fill(GenJet_mass[iJet]);
            h_Gen_Jet_phi ->Fill(GenJet_phi[iJet]);
            h_Gen_Jet_pt ->Fill(GenJet_pt[iJet]);

        }
        h_nGenParts->Fill(nGenPart);
        for(int iJet=0; iJet<nGenPart; iJet++){
            h_Gen_Part_eta ->Fill(GenPart_eta[iJet]);
            h_Gen_Part_mass ->Fill(GenPart_mass[iJet]);
            h_Gen_Part_phi ->Fill(GenPart_phi[iJet]);
            h_Gen_Part_pt ->Fill(GenPart_pt[iJet]);

        }


        //match LHE jets - dr^2 = deta^2 + dphi^2
        int lheJet_1 = -1, lheJet_2 = -1;
        for(int iLHE=0; iLHE<nLHEPart; iLHE++){
            if(fabs(LHEPart_pdgId[iLHE]) > 6) continue;
            if(lheJet_1 == -1) lheJet_1 = iLHE;
            else lheJet_2 = iLHE;
        }

        bool LHE_pass_Selection = true;

        if(LHEPart_pt[lheJet_2] > LHEPart_pt[lheJet_1]) swap(lheJet_2, lheJet_1);
        if(LHEPart_pt[lheJet_2] < jetPT_cut || LHEPart_pt[lheJet_1] < jetPT_cut) LHE_pass_Selection = false;
        if(fabs(LHEPart_eta[lheJet_1] - LHEPart_eta[lheJet_2]) < jetEtaSep_cut) LHE_pass_Selection = false;

        ROOT::Math::PtEtaPhiMVector lheJet_1_p4(LHEPart_pt[lheJet_1], LHEPart_eta[lheJet_1], LHEPart_phi[lheJet_1], LHEPart_mass[lheJet_1]);
        ROOT::Math::PtEtaPhiMVector lheJet_2_p4(LHEPart_pt[lheJet_2], LHEPart_eta[lheJet_2], LHEPart_phi[lheJet_2], LHEPart_mass[lheJet_2]);
        ROOT::Math::PtEtaPhiMVector jetSum_lhe = lheJet_1_p4+lheJet_2_p4;
        double invMass_lhe = sqrt(jetSum_lhe.Dot(jetSum_lhe));

        if(invMass_lhe < jetInvMass_cut) LHE_pass_Selection = false;

        if(LHE_pass_Selection == true){

            h_LHE_Leading_Jet_eta ->Fill(LHEPart_eta[lheJet_1]);
            h_LHE_Leading_Jet_mass ->Fill(LHEPart_mass[lheJet_1]);
            h_LHE_Leading_Jet_phi ->Fill(LHEPart_phi[lheJet_1]);
            h_LHE_Leading_Jet_pt ->Fill(LHEPart_pt[lheJet_1]);

            h_LHE_Trailing_Jet_eta ->Fill(LHEPart_eta[lheJet_2]);
            h_LHE_Trailing_Jet_mass ->Fill(LHEPart_mass[lheJet_2]);
            h_LHE_Trailing_Jet_phi ->Fill(LHEPart_phi[lheJet_2]);
            h_LHE_Trailing_Jet_pt ->Fill(LHEPart_pt[lheJet_2]);

            h_LHE_Jet_invariantMass->Fill(invMass_lhe);
            h_LHE_Jet_etaSeparation->Fill(fabs(LHEPart_eta[lheJet_1] - LHEPart_eta[lheJet_2]));
        


            //match jets
            vector<int> matchedJet_candidates_1, matchedJet_candidates_2;
            float dR_cut = 0.4;
            for(int iJet=0; iJet<nJet; iJet++){
                float dPhi_1 = LHEPart_phi[lheJet_1] - Jet_phi[iJet];
                float dPhi_2 = LHEPart_phi[lheJet_2] - Jet_phi[iJet];

                float dEta_1 = LHEPart_eta[lheJet_1] - Jet_eta[iJet];
                float dEta_2 = LHEPart_eta[lheJet_2] - Jet_eta[iJet];

                float dR_1 = sqrt( pow(dPhi_1, 2) + pow(dEta_1, 2) );
                float dR_2 = sqrt( pow(dPhi_2, 2) + pow(dEta_2, 2) );

                if(dR_1 < dR_cut) matchedJet_candidates_1.push_back(iJet);
                if(dR_2 < dR_cut) matchedJet_candidates_2.push_back(iJet);
            }
            
            if(matchedJet_candidates_1.size() != 0 && matchedJet_candidates_2.size() != 0) {
                int matchedJet_1 = matchedJet_candidates_1.at(0);
                int matchedJet_2 = matchedJet_candidates_2.at(0);

                if(matchedJet_candidates_1.size() > 1){
                    for(int iJet=1; iJet<matchedJet_candidates_1.size(); iJet++){
                        if(Jet_pt[matchedJet_candidates_1.at(iJet)] > Jet_pt[matchedJet_1]) matchedJet_1 = matchedJet_candidates_1.at(iJet);
                    }
                }
                if(matchedJet_candidates_2.size() > 1){
                    for(int iJet=1; iJet<matchedJet_candidates_2.size(); iJet++){
                        if(Jet_pt[matchedJet_candidates_2.at(iJet)] > Jet_pt[matchedJet_2]) matchedJet_2 = matchedJet_candidates_2.at(iJet);
                    }
                }

                if(Jet_pt[matchedJet_2] > Jet_pt[matchedJet_1]) swap(matchedJet_2, matchedJet_1);

                ROOT::Math::PtEtaPhiMVector matchedJet_1_p4(Jet_pt[matchedJet_1], Jet_eta[matchedJet_1], Jet_phi[matchedJet_1], Jet_mass[matchedJet_1]);
                ROOT::Math::PtEtaPhiMVector matchedJet_2_p4(Jet_pt[matchedJet_2], Jet_eta[matchedJet_2], Jet_phi[matchedJet_2], Jet_mass[matchedJet_2]);
                ROOT::Math::PtEtaPhiMVector jetSum_match = matchedJet_1_p4+matchedJet_2_p4;
                double invMass_matched = sqrt(jetSum_match.Dot(jetSum_match));


                h_Matched_Leading_Jet_eta ->Fill(Jet_eta[matchedJet_1]);
                h_Matched_Leading_Jet_mass ->Fill(Jet_mass[matchedJet_1]);
                h_Matched_Leading_Jet_phi ->Fill(Jet_phi[matchedJet_1]);
                h_Matched_Leading_Jet_pt ->Fill(Jet_pt[matchedJet_1]);

                h_Matched_Trailing_Jet_eta ->Fill(Jet_eta[matchedJet_2]);
                h_Matched_Trailing_Jet_mass ->Fill(Jet_mass[matchedJet_2]);
                h_Matched_Trailing_Jet_phi ->Fill(Jet_phi[matchedJet_2]);
                h_Matched_Trailing_Jet_pt ->Fill(Jet_pt[matchedJet_2]);

                h_Matched_Jet_invariantMass->Fill(invMass_matched);
                h_Matched_Jet_etaSeparation->Fill(fabs(Jet_eta[matchedJet_1] - Jet_eta[matchedJet_2]));
            }
        }

        

        //find jet pairs
        //requirements: dEta > 2, m_jj > 400, jet_pT > 50
        vector<std::pair<int,int>> jet_pairs;
        for(int iJet_1=0; iJet_1<nJet; iJet_1++){
            if(Jet_pt[iJet_1] < jetPT_cut) continue;

            for(int iJet_2=iJet_1+1; iJet_2<nJet; iJet_2++){
                if(Jet_pt[iJet_2] < jetPT_cut) continue;
                if(fabs(Jet_eta[iJet_1] - Jet_eta[iJet_2]) < jetEtaSep_cut) continue;

                ROOT::Math::PtEtaPhiMVector iJet_1_p4(Jet_pt[iJet_1], Jet_eta[iJet_1], Jet_phi[iJet_1], Jet_mass[iJet_1]);
                ROOT::Math::PtEtaPhiMVector iJet_2_p4(Jet_pt[iJet_2], Jet_eta[iJet_2], Jet_phi[iJet_2], Jet_mass[iJet_2]);
                ROOT::Math::PtEtaPhiMVector jetSum = iJet_1_p4+iJet_2_p4;
                double invMass = sqrt(jetSum.Dot(jetSum));

                if(invMass < jetInvMass_cut) continue;

                if(std::find(jet_pairs.begin(), jet_pairs.end(), std::make_pair(iJet_1, iJet_2)) != jet_pairs.end() ||
                   std::find(jet_pairs.begin(), jet_pairs.end(), std::make_pair(iJet_2, iJet_1)) != jet_pairs.end())
                    continue;

                jet_pairs.push_back(std::make_pair(iJet_1, iJet_2));

            }
        }
        if(jet_pairs.size() == 0) continue;

        int leadJet_1 = jet_pairs.at(0).first;
        int leadJet_2 = jet_pairs.at(0).second;
        if(jet_pairs.size() > 1){

            for(int i=1; i<jet_pairs.size(); i++){
                int tempJet_1 = jet_pairs.at(i).first;
                int tempJet_2 = jet_pairs.at(i).second;

                ROOT::Math::PtEtaPhiMVector leadJet_1_p4(Jet_pt[leadJet_1], Jet_eta[leadJet_1], Jet_phi[leadJet_1], Jet_mass[leadJet_1]);
                ROOT::Math::PtEtaPhiMVector leadJet_2_p4(Jet_pt[leadJet_2], Jet_eta[leadJet_2], Jet_phi[leadJet_2], Jet_mass[leadJet_2]);
                ROOT::Math::PtEtaPhiMVector leadJetSum = leadJet_1_p4+leadJet_2_p4;
                double invMass_lead = sqrt(leadJetSum.Dot(leadJetSum));

                ROOT::Math::PtEtaPhiMVector tempJet_1_p4(Jet_pt[tempJet_1], Jet_eta[tempJet_1], Jet_phi[tempJet_1], Jet_mass[tempJet_1]);
                ROOT::Math::PtEtaPhiMVector tempJet_2_p4(Jet_pt[tempJet_2], Jet_eta[tempJet_2], Jet_phi[tempJet_2], Jet_mass[tempJet_2]);
                ROOT::Math::PtEtaPhiMVector tempJetSum = tempJet_1_p4+tempJet_2_p4;
                double invMass_temp= sqrt(leadJetSum.Dot(tempJetSum));

                if(invMass_temp > invMass_lead ){
                    leadJet_1 = tempJet_1;
                    leadJet_2 = tempJet_2;
                }
            }
        }
        /*int leadJet_1 = 0; int leadJet_2 = -1;

        for(int iJet=1; iJet<nJet; iJet++){
            if(Jet_pt[iJet] > Jet_pt[leadJet_1]){
                leadJet_2 = leadJet_1;
                leadJet_1 = iJet;
            }
            if(leadJet_2 != -1){
                if(Jet_pt[iJet] > Jet_pt[leadJet_2] && Jet_pt[iJet] < Jet_pt[leadJet_1]){
                    leadJet_2 = iJet;
                }
            }
            else leadJet_2 = iJet;

        }*/
        if(leadJet_2 == -1) continue;

        if(Jet_pt[leadJet_2] > Jet_pt[leadJet_1]) swap(leadJet_2,leadJet_1);

        ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[leadJet_1], Jet_eta[leadJet_1], Jet_phi[leadJet_1], Jet_mass[leadJet_1]);
        ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[leadJet_2], Jet_eta[leadJet_2], Jet_phi[leadJet_2], Jet_mass[leadJet_2]);
        ROOT::Math::PtEtaPhiMVector jetSum = leadJet_p4+trailJet_p4;
        double invMass = sqrt(jetSum.Dot(jetSum));

        //if(fabs(Jet_eta[leadJetIdx]) < 2 || fabs(Jet_eta[trailJetIdx]) < 2) continue;

        //find highest pt fatJet
        int leadFatJetIdx = -1;
        for(int iFatJet=0; iFatJet<nFatJet; iFatJet++){
            if(fabs(FatJet_eta[iFatJet]) < 2){
                if(leadFatJetIdx == -1) leadFatJetIdx = iFatJet;
                else if(FatJet_pt[iFatJet] > FatJet_pt[leadFatJetIdx]) leadFatJetIdx = iFatJet;
            }
        }


        h_CaloMET_phi->Fill(CaloMET_phi);
        h_CaloMET_pt ->Fill(CaloMET_pt);
        h_CaloMET_sumEt ->Fill(CaloMET_sumEt);


        for(unsigned int iElectron=0; iElectron<nElectron; iElectron++){

            h_Electron_eta ->Fill(Electron_eta[iElectron]);
            h_Electron_hoe ->Fill(Electron_hoe[iElectron]); 
            h_Electron_mass ->Fill(Electron_mass[iElectron]);
            h_Electron_phi ->Fill(Electron_phi[iElectron]);
            h_Electron_pt ->Fill(Electron_pt[iElectron]);
            h_Electron_r9 ->Fill(Electron_r9[iElectron]);
            h_Electron_sieie ->Fill(Electron_sieie[iElectron]);
        }

        for(int iMuon=0; iMuon<nMuon; iMuon++){

            h_Muon_eta  ->Fill(Muon_eta[iMuon]);
            h_Muon_mass ->Fill(Muon_mass[iMuon]);
            h_Muon_phi ->Fill(Muon_phi[iMuon]);
            h_Muon_pt ->Fill(Muon_pt[iMuon]);
        }

        h_nJets->Fill(nJet);
        for(int iJet=0; iJet<nJet; iJet++){
            h_Jet_eta ->Fill(Jet_eta[iJet]);
            h_Jet_mass ->Fill(Jet_mass[iJet]);
            h_Jet_phi ->Fill(Jet_phi[iJet]);
            h_Jet_pt ->Fill(Jet_pt[iJet]);

        }

        h_nFatJets -> Fill(nFatJet);
        for(int iFatJet=0; iFatJet<nFatJet; iFatJet++){
            h_fatJet_eta ->Fill(FatJet_eta[iFatJet]);
            h_fatJet_mass ->Fill(FatJet_mass[iFatJet]);
            h_fatJet_phi ->Fill(FatJet_phi[iFatJet]);
            h_fatJet_pt ->Fill(FatJet_pt[iFatJet]);

        }

        h_Leading_Jet_eta ->Fill(Jet_eta[leadJet_1]);
        h_Leading_Jet_mass ->Fill(Jet_mass[leadJet_1]);
        h_Leading_Jet_phi ->Fill(Jet_phi[leadJet_1]);
        h_Leading_Jet_pt ->Fill(Jet_pt[leadJet_1]);

        h_Trailing_Jet_eta ->Fill(Jet_eta[leadJet_2]);
        h_Trailing_Jet_mass ->Fill(Jet_mass[leadJet_2]);
        h_Trailing_Jet_phi ->Fill(Jet_phi[leadJet_2]);
        h_Trailing_Jet_pt ->Fill(Jet_pt[leadJet_2]);

        h_Jet_invariantMass->Fill(invMass);
        h_Jet_etaSeparation->Fill(fabs(Jet_eta[leadJet_1] - Jet_eta[leadJet_2]));


        if(leadFatJetIdx != -1){
            h_Leading_FatJet_eta ->Fill(FatJet_eta[leadFatJetIdx]);
            h_Leading_FatJet_mass ->Fill(FatJet_mass[leadFatJetIdx]);
            h_Leading_FatJet_phi ->Fill(FatJet_phi[leadFatJetIdx]);
            h_Leading_FatJet_pt ->Fill(FatJet_pt[leadFatJetIdx]);
            h_Leading_FatJet_deepTagMD_bbvsLight->Fill(FatJet_deepTagMD_bbvsLight[leadFatJetIdx]);
            h_Leading_FatJet_deepTag_H->Fill(FatJet_deepTag_H[leadFatJetIdx]);
            h_Leading_FatJet_deepTag_WvsQCD->Fill(FatJet_deepTag_WvsQCD[leadFatJetIdx]);
        }

    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;
}

void InitHistograms(){

    h_CaloMET_phi = new TH1F("h_CaloMET_phi","h_CaloMET_phi", 200, -3.5, 3.5);
    h_CaloMET_pt = new TH1F("h_CaloMET_pt","h_CaloMET_pt", 200, 0, 200);
    h_CaloMET_sumEt = new TH1F("h_CaloMET_sumEt","h_CaloMET_sumEt", 500, 0, 1000);

    h_Electron_eta = new TH1F("h_Electron_eta","h_Electron_eta", 200, -5.0, 5.0);
    h_Electron_hoe = new TH1F("h_Electron_hoe","h_Electron_hoe", 200, 0, 5.0);
    h_Electron_mass = new TH1F("h_Electron_mass","h_Electron_mass", 200, 0, 0.15);
    h_Electron_phi = new TH1F("h_Electron_phi","h_Electron_phi", 200, -3.5, 3.5);
    h_Electron_pt = new TH1F("h_Electron_pt","h_Electron_pt", 200, 0, 200);
    h_Electron_r9 = new TH1F("h_Electron_r9","h_Electron_r9", 200, 0, 4.0);
    h_Electron_sieie = new TH1F("h_Electron_sieie","h_Electron_sieie", 200, 0, 0.2);

    h_Muon_eta = new TH1F("Muon_eta","Muon_eta", 200, -5.0, 5.0);
    h_Muon_mass = new TH1F("h_Muon_mass","h_Muon_mass", 200, 0, 0.15);
    h_Muon_phi = new TH1F("h_Muon_phi","h_Muon_phi", 200, -3.5, 3.5);
    h_Muon_pt = new TH1F("h_Muon_pt","h_Muon_pt", 200, 0, 200);

    h_Jet_eta = new TH1F("h_Jet_eta","h_Jet_eta", 200, -5.0, 5.0);
    h_Jet_mass = new TH1F("h_Jet_mass","h_Jet_mass", 200, 0, 200);
    h_Jet_phi = new TH1F("h_Jet_phi","h_Jet_phi", 200, -3.5, 3.5);
    h_Jet_pt = new TH1F("h_Jet_pt","h_Jet_pt", 200, 0, 200);
    h_nJets = new TH1F("h_nJets","h_nJets", 50, 0.0, 50);

    h_fatJet_eta = new TH1F("h_fatJet_eta","h_fatJet_eta", 200, -5.0, 5.0);
    h_fatJet_mass = new TH1F("h_fatJet_mass","h_fatJet_mass", 200, 0, 200);
    h_fatJet_phi = new TH1F("h_fatJet_phi","h_fatJet_phi", 200, -3.5, 3.5);
    h_fatJet_pt = new TH1F("h_fatJet_pt","h_fatJet_pt", 100, 0, 2000);
    h_nFatJets = new TH1F("h_nFatJets","h_nFatJets", 50, 0.0, 50);

    h_Leading_Jet_eta = new TH1F("h_Leading_Jet_eta","h_Leading_Jet_eta", 200, -5.0, 5.0);
    h_Leading_Jet_mass = new TH1F("h_Leading_Jet_mass","h_Leading_Jet_mass", 200, 0, 200);
    h_Leading_Jet_phi = new TH1F("h_Leading_Jet_phi","h_Leading_Jet_phi", 200, -3.5, 3.5);
    h_Leading_Jet_pt = new TH1F("h_Leading_Jet_pt","h_Leading_Jet_pt", 100, 0, 2000);

    h_Trailing_Jet_eta = new TH1F("h_Trailing_Jet_eta","h_Trailing_Jet_eta", 200, -5.0, 5.0);
    h_Trailing_Jet_mass = new TH1F("h_Trailing_Jet_mass","h_Trailing_Jet_mass", 200, 0, 200);
    h_Trailing_Jet_phi = new TH1F("h_Trailing_Jet_phi","h_Trailing_Jet_phi", 200, -3.5, 3.5);
    h_Trailing_Jet_pt = new TH1F("h_Trailing_Jet_pt","h_Trailing_Jet_pt", 100, 0, 2000);

    h_Jet_invariantMass = new TH1F("h_Jet_invariantMass","h_Jet_invariantMass", 100, 0, 4000);
    h_Jet_etaSeparation = new TH1F("h_Jet_etaSeparation","h_Jet_etaSeparation", 100, 0, 10);
    
    h_Leading_FatJet_eta = new TH1F("h_Leading_FatJet_eta","h_Leading_FatJet_eta", 200, -5.0, 5.0);
    h_Leading_FatJet_mass = new TH1F("h_Leading_FatJet_mass","h_Leading_FatJet_mass", 200, 0, 200);
    h_Leading_FatJet_phi = new TH1F("h_Leading_FatJet_phi","h_Leading_FatJet_phi", 200, -3.5, 3.5);
    h_Leading_FatJet_pt = new TH1F("h_Leading_FatJet_pt","h_Leading_FatJet_pt", 100, 0, 2000);
    h_Leading_FatJet_deepTag_H = new TH1F("h_Leading_FatJet_deepTag_H","h_Leading_FatJet_deepTag_H", 100, 0, 1);
    h_Leading_FatJet_deepTag_WvsQCD = new TH1F("h_Leading_FatJet_deepTag_WvsQCD","h_Leading_FatJet_deepTag_WvsQCD", 100, 0, 1);
    h_Leading_FatJet_deepTagMD_bbvsLight = new TH1F("h_Leading_FatJet_deepTagMD_bbvsLight","h_Leading_FatJet_deepTagMD_bbvsLight", 100, 0, 1);

    h_Gen_Jet_eta = new TH1F("h_Gen_Jet_eta","h_Gen_Jet_eta", 200, -5.0, 5.0);
    h_Gen_Jet_mass = new TH1F("h_Gen_Jet_mass","h_Gen_Jet_mass", 200, 0, 200);
    h_Gen_Jet_phi = new TH1F("h_Gen_Jet_phi","h_Gen_Jet_phi", 200, -3.5, 3.5);
    h_Gen_Jet_pt = new TH1F("h_Gen_Jet_pt","h_Gen_Jet_pt", 200, 0, 200);
    h_nGenJets = new TH1F("h_nGenJets","h_nGenJets", 50, 0.0, 50);

    h_Gen_Part_eta = new TH1F("h_Gen_Part_eta","h_Gen_Part_eta", 200, -5.0, 5.0);
    h_Gen_Part_mass = new TH1F("h_Gen_Part_mass","h_Gen_Part_mass", 200, 0, 200);
    h_Gen_Part_phi = new TH1F("h_Gen_Part_phi","h_Gen_Part_phi", 200, -3.5, 3.5);
    h_Gen_Part_pt = new TH1F("h_Gen_Part_pt","h_Gen_Part_pt", 200, 0, 200);
    h_nGenParts = new TH1F("h_nGenParts","h_nGenParts", 50, 0.0, 50);

    h_Matched_Leading_Jet_eta = new TH1F("h_Matched_Leading_Jet_eta","h_Matched_Leading_Jet_eta", 200, -5.0, 5.0);
    h_Matched_Leading_Jet_mass = new TH1F("h_Matched_Leading_Jet_mass","h_Matched_Leading_Jet_mass", 200, 0, 200);
    h_Matched_Leading_Jet_phi = new TH1F("h_Matched_Leading_Jet_phi","h_Matched_Leading_Jet_phi", 200, -3.5, 3.5);
    h_Matched_Leading_Jet_pt = new TH1F("h_Matched_Leading_Jet_pt","h_Matched_Leading_Jet_pt", 100, 0, 2000);

    h_Matched_Trailing_Jet_eta = new TH1F("h_Matched_Trailing_Jet_eta","h_Matched_Trailing_Jet_eta", 200, -5.0, 5.0);
    h_Matched_Trailing_Jet_mass = new TH1F("h_Matched_Trailing_Jet_mass","h_Matched_Trailing_Jet_mass", 200, 0, 200);
    h_Matched_Trailing_Jet_phi = new TH1F("h_Matched_Trailing_Jet_phi","h_Matched_Trailing_Jet_phi", 200, -3.5, 3.5);
    h_Matched_Trailing_Jet_pt = new TH1F("h_Matched_Trailing_Jet_pt","h_Matched_Trailing_Jet_pt", 100, 0, 2000);

    h_Matched_Jet_invariantMass = new TH1F("h_Matched_Jet_invariantMass","h_Matched_Jet_invariantMass", 100, 0, 4000);
    h_Matched_Jet_etaSeparation = new TH1F("h_Matched_Jet_etaSeparation","h_Matched_Jet_etaSeparation", 100, 0, 10);
    
    h_LHE_Leading_Jet_eta = new TH1F("h_LHE_Leading_Jet_eta","h_LHE_Leading_Jet_eta", 200, -5.0, 5.0);
    h_LHE_Leading_Jet_mass = new TH1F("h_LHE_Leading_Jet_mass","h_LHE_Leading_Jet_mass", 200, 0, 200);
    h_LHE_Leading_Jet_phi = new TH1F("h_LHE_Leading_Jet_phi","h_LHE_Leading_Jet_phi", 200, -3.5, 3.5);
    h_LHE_Leading_Jet_pt = new TH1F("h_LHE_Leading_Jet_pt","h_LHE_Leading_Jet_pt", 100, 0, 2000);

    h_LHE_Trailing_Jet_eta = new TH1F("h_LHE_Trailing_Jet_eta","h_LHE_Trailing_Jet_eta", 200, -5.0, 5.0);
    h_LHE_Trailing_Jet_mass = new TH1F("h_LHE_Trailing_Jet_mass","h_LHE_Trailing_Jet_mass", 200, 0, 200);
    h_LHE_Trailing_Jet_phi = new TH1F("h_LHE_Trailing_Jet_phi","h_LHE_Trailing_Jet_phi", 200, -3.5, 3.5);
    h_LHE_Trailing_Jet_pt = new TH1F("h_LHE_Trailing_Jet_pt","h_LHE_Trailing_Jet_pt", 100, 0, 2000);

    h_LHE_Jet_invariantMass = new TH1F("h_LHE_Jet_invariantMass","h_LHE_Jet_invariantMass", 100, 0, 4000);
    h_LHE_Jet_etaSeparation = new TH1F("h_LHE_Jet_etaSeparation","h_LHE_Jet_etaSeparation", 100, 0, 10);
    
}

void DrawPlot(TH1F* plot, string name){

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    string outdir = "../Output/ttHTobb/";
    setOutput(outdir);

    plot->Scale(1.0 / plot->Integral());
    plot->SetFillColorAlpha(kGreen + 1,0.35);
    plot->SetTitle(name.c_str());
    plot->Draw("hist");
    c1->SaveAs((outdir+name+".png").c_str());
    c1->SaveAs((outdir+name+".pdf").c_str());

    c1->Delete();
}

void DrawPlot(TH1F* h_1, TH1F* h_2, string h_1_label, string h_2_label, string name){

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    string outdir = "../Output/ttHTobb/";
    setOutput(outdir);

    h_1->Scale(1.0 / h_1->Integral());
    h_2->Scale(1.0 / h_2->Integral());

    std::vector<float> maxima;
    maxima.resize(2);
    maxima[0] = h_1->GetMaximum(); 
    maxima[1] = h_2->GetMaximum();
    std::sort(maxima.begin(),maxima.end()); 
    h_1->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   
    h_1->SetTitle(name.c_str());

    h_1->SetFillColorAlpha(kRed + 1,0.35);
    h_2->SetFillColorAlpha(kBlue + 1,0.35);

   TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();
   

    TLegend* legend = new TLegend(0.79, 0.5, 0.99, 0.66);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04); 
    legend -> AddEntry(h_1,h_1_label.c_str(),"F");
    legend -> AddEntry(h_2,h_2_label.c_str(),"F");
    h_1->Draw("hist");
    gPad -> Update();
    st_old= (TPaveStats*)(h_1->GetListOfFunctions()->FindObject("stats"));
    st_old->SetX1NDC(0.82); //new x start position
    st_old->SetX2NDC(0.99); //new x end position
    st_old->SetY1NDC(0.82); //new y start position
    st_old->SetY2NDC(0.94); //new y end position
    st_old->SetTextColor(kRed+1);
    st_old->Draw("sames");

    h_2->Draw("hist,sames");
    gPad -> Update();
    st_new= (TPaveStats*)(h_2->GetListOfFunctions()->FindObject("stats"));
    st_new->SetX1NDC(0.82); //new x start position
    st_new->SetX2NDC(0.99); //new x end position
    st_new->SetY1NDC(0.68); //new y start position
    st_new->SetY2NDC(0.80); //new y end position
    st_new->SetTextColor(kBlue+1);
    st_new->Draw("sames");
    legend -> Draw("same");

    //h_1->Draw("hist");
    //h_2->Draw("hist, same");

    c1->SaveAs((outdir+name+".png").c_str());
    c1->SaveAs((outdir+name+".pdf").c_str());

    c1->Delete();
}

void DrawPlot(TH1F* h_1, TH1F* h_2, TH1F* h_3, string h_1_label, string h_2_label, string h_3_label, string name){

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    string outdir = "../Output/ttHTobb/";
    setOutput(outdir);

    h_1->Scale(1.0 / h_1->Integral());
    h_2->Scale(1.0 / h_2->Integral());
    h_3->Scale(1.0 / h_3->Integral());

    std::vector<float> maxima;
    maxima.resize(3);
    maxima[0] = h_1->GetMaximum(); 
    maxima[1] = h_2->GetMaximum();
    maxima[2] = h_3->GetMaximum();
    std::sort(maxima.begin(),maxima.end()); 
    h_1->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   
    h_1->SetTitle(name.c_str());
    
    h_1->SetFillColorAlpha(kRed + 1,0.35);
    h_1->SetLineColor(kRed + 1);
    h_2->SetFillColorAlpha(kBlue + 1,0.35);
    h_2->SetLineColor(kBlue + 1);
    h_3->SetFillColorAlpha(kGreen + 1,0.35);
    h_3->SetLineColor(kGreen + 1);

    TPaveStats* st_old = new TPaveStats();
    TPaveStats* st_new = new TPaveStats();
    TPaveStats* st_3 = new TPaveStats();
    TPaveStats* st_ratio = new TPaveStats();
        
    //TCanvas* c = new TCanvas();
    
    TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
    TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
    
    cUp->SetBottomMargin(0.01); 
    cDown->SetTopMargin(0.01); 
    cDown->SetBottomMargin(0.2); 
        
    cUp->Draw();
    cDown->Draw();
        
    cUp->cd();
   

    TLegend* legend = new TLegend(0.65, 0.78, 0.8, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04); 
    legend -> AddEntry(h_1,h_1_label.c_str(),"F");
    legend -> AddEntry(h_2,h_2_label.c_str(),"F");
    legend -> AddEntry(h_3,h_3_label.c_str(),"F");
    h_1->Draw("hist");
    gPad -> Update();
    st_old= (TPaveStats*)(h_1->GetListOfFunctions()->FindObject("stats"));
    st_old->SetX1NDC(0.82); //new x start position
    st_old->SetX2NDC(0.99); //new x end position
    st_old->SetY1NDC(0.82); //new y start position
    st_old->SetY2NDC(0.94); //new y end position
    st_old->SetTextColor(kRed+1);
    st_old->Draw("sames");

    h_2->Draw("hist,sames");
    gPad -> Update();
    st_new= (TPaveStats*)(h_2->GetListOfFunctions()->FindObject("stats"));
    st_new->SetX1NDC(0.82); //new x start position
    st_new->SetX2NDC(0.99); //new x end position
    st_new->SetY1NDC(0.68); //new y start position
    st_new->SetY2NDC(0.80); //new y end position
    st_new->SetTextColor(kBlue+1);
    st_new->Draw("sames");

    h_3->Draw("hist,sames");
    gPad -> Update();
    st_3= (TPaveStats*)(h_3->GetListOfFunctions()->FindObject("stats"));
    st_3->SetX1NDC(0.82); //new x start position
    st_3->SetX2NDC(0.99); //new x end position
    st_3->SetY1NDC(0.54); //new y start position
    st_3->SetY2NDC(0.66); //new y end position
    st_3->SetTextColor(kGreen+1);
    st_3->Draw("sames");

    legend -> Draw("same");


   cDown->cd();
    
   TH1F* histo_ratio=(TH1F*)h_2->Clone("histo_ratio");
   TH1F* histo_ratio_3=(TH1F*)h_3->Clone("histo_ratio_3");
   if(histo_ratio->GetSumw2N()<=0) histo_ratio->Sumw2();
   if(histo_ratio_3->GetSumw2N()<=0) histo_ratio_3->Sumw2();
   histo_ratio->Divide(h_1);
   histo_ratio_3->Divide(h_1);
    
   //histo_ratio -> GetXaxis() -> SetTitle(x_label.c_str());
   histo_ratio -> GetYaxis() -> SetTitle("Jet/LHE Jet");
   histo_ratio -> SetMaximum(2);
   histo_ratio -> SetMinimum(0);
   histo_ratio -> SetMarkerColor(kBlue);
   histo_ratio -> SetLineColor(kBlue);
   histo_ratio -> SetMarkerSize(0.5);
   histo_ratio ->SetTitle("");
   histo_ratio ->SetStats(0);
   histo_ratio_3 ->SetStats(0);

   histo_ratio_3 -> SetMarkerColor(kGreen);
   histo_ratio_3 -> SetLineColor(kGreen);
   histo_ratio_3 -> SetMarkerSize(0.5);
   histo_ratio_3 ->SetTitle("");
   histo_ratio -> GetXaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetXaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleOffset(0.7);
   histo_ratio -> Draw("e");
   TF1* f_const = new TF1("f_1", "[0]",histo_ratio->GetBinCenter(1)-histo_ratio->GetBinWidth(1)/2, histo_ratio->GetBinCenter(histo_ratio->GetNbinsX())+histo_ratio->GetBinWidth(histo_ratio->GetNbinsX())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   f_const -> SetLineWidth(2);
   f_const -> Draw("same");
   histo_ratio -> Draw("e,sames");
   histo_ratio_3 -> Draw("e,sames");
    




    //h_1->Draw("hist");
    //h_2->Draw("hist, same");

    c1->SaveAs((outdir+name+".png").c_str());
    c1->SaveAs((outdir+name+".pdf").c_str());

    c1->Delete();
}

void SaveHistograms(){
    
    DrawPlot(h_CaloMET_phi, "h_CaloMET_phi");
    DrawPlot(h_CaloMET_pt, "h_CaloMET_pt");
    DrawPlot(h_CaloMET_sumEt, "h_CaloMET_sumEt");

    DrawPlot(h_Electron_eta,"h_Electron_eta");
    DrawPlot(h_Electron_hoe,"h_Electron_hoe");
    DrawPlot(h_Electron_mass,"h_Electron_mass");
    DrawPlot(h_Electron_phi,"h_Electron_phi");
    DrawPlot(h_Electron_pt,"h_Electron_pt");
    DrawPlot(h_Electron_r9,"h_Electron_r9");
    DrawPlot(h_Electron_sieie,"h_Electron_sieie");

    DrawPlot(h_Muon_eta,"Muon_eta");
    DrawPlot(h_Muon_mass,"h_Muon_mass");
    DrawPlot(h_Muon_phi,"h_Muon_phi");
    DrawPlot(h_Muon_pt,"h_Muon_pt");

    DrawPlot(h_Jet_eta,"h_Jet_eta");
    DrawPlot(h_Jet_mass,"h_Jet_mass");
    DrawPlot(h_Jet_phi,"h_Jet_phi");
    DrawPlot(h_Jet_pt,"h_Jet_pt");
    DrawPlot(h_nJets,"h_nJets");

    DrawPlot(h_Gen_Jet_eta,"h_Gen_Jet_eta");
    DrawPlot(h_Gen_Jet_mass,"h_Gen_Jet_mass");
    DrawPlot(h_Gen_Jet_phi,"h_Gen_Jet_phi");
    DrawPlot(h_Gen_Jet_pt,"h_Gen_Jet_pt");
    DrawPlot(h_nGenJets,"h_nGenJets");

    DrawPlot(h_Gen_Part_eta,"h_Gen_Part_eta");
    DrawPlot(h_Gen_Part_mass,"h_Gen_Part_mass");
    DrawPlot(h_Gen_Part_phi,"h_Gen_Part_phi");
    DrawPlot(h_Gen_Part_pt,"h_Gen_Part_pt");
    DrawPlot(h_nGenParts,"h_nGenParts");

    DrawPlot(h_fatJet_eta,"h_fatJet_eta");
    DrawPlot(h_fatJet_mass,"h_fatJet_mass");
    DrawPlot(h_fatJet_phi,"h_fatJet_phi");
    DrawPlot(h_fatJet_pt,"h_fatJet_pt");
    DrawPlot(h_nFatJets,"h_nFatJets");

    DrawPlot(h_Leading_Jet_eta,"h_Leading_Jet_eta");
    DrawPlot(h_Leading_Jet_mass,"h_Leading_Jet_mass");
    DrawPlot(h_Leading_Jet_phi,"h_Leading_Jet_phi");
    DrawPlot(h_Leading_Jet_pt,"h_Leading_Jet_pt");

    DrawPlot(h_Trailing_Jet_eta,"h_Trailing_Jet_eta");
    DrawPlot(h_Trailing_Jet_mass,"h_Trailing_Jet_mass");
    DrawPlot(h_Trailing_Jet_phi,"h_Trailing_Jet_phi");
    DrawPlot(h_Trailing_Jet_pt,"h_Trailing_Jet_pt");

    DrawPlot(h_Jet_invariantMass,"h_Jet_invariantMass");    
    DrawPlot(h_Jet_etaSeparation,"h_Jet_etaSeparation");

    DrawPlot(h_Leading_FatJet_eta,"h_Leading_FatJet_eta");
    DrawPlot(h_Leading_FatJet_mass,"h_Leading_FatJet_mass");
    DrawPlot(h_Leading_FatJet_phi,"h_Leading_FatJet_phi");
    DrawPlot(h_Leading_FatJet_pt,"h_Leading_FatJet_pt");
    DrawPlot(h_Leading_FatJet_deepTag_H,"h_Leading_FatJet_deepTag_H");   
    DrawPlot(h_Leading_FatJet_deepTag_WvsQCD,"h_Leading_FatJet_deepTag_WvsQCD"); 
    DrawPlot(h_Leading_FatJet_deepTagMD_bbvsLight,"h_Leading_FatJet_deepTagMD_bbvsLight");

    DrawPlot(h_Matched_Leading_Jet_eta,"h_Matched_Leading_Jet_eta");
    DrawPlot(h_Matched_Leading_Jet_mass,"h_Matched_Leading_Jet_mass");
    DrawPlot(h_Matched_Leading_Jet_phi,"h_Matched_Leading_Jet_phi");
    DrawPlot(h_Matched_Leading_Jet_pt,"h_Matched_Leading_Jet_pt");

    DrawPlot(h_Matched_Trailing_Jet_eta,"h_Matched_Trailing_Jet_eta");
    DrawPlot(h_Matched_Trailing_Jet_mass,"h_Matched_Trailing_Jet_mass");
    DrawPlot(h_Matched_Trailing_Jet_phi,"h_Matched_Trailing_Jet_phi");
    DrawPlot(h_Matched_Trailing_Jet_pt,"h_Matched_Trailing_Jet_pt");

    DrawPlot(h_Matched_Jet_invariantMass,"h_Matched_Jet_invariantMass");    
    DrawPlot(h_Matched_Jet_etaSeparation,"h_Matched_Jet_etaSeparation");

    DrawPlot(h_LHE_Leading_Jet_eta,"h_LHE_Leading_Jet_eta");
    DrawPlot(h_LHE_Leading_Jet_mass,"h_LHE_Leading_Jet_mass");
    DrawPlot(h_LHE_Leading_Jet_phi,"h_LHE_Leading_Jet_phi");
    DrawPlot(h_LHE_Leading_Jet_pt,"h_LHE_Leading_Jet_pt");

    DrawPlot(h_LHE_Trailing_Jet_eta,"h_LHE_Trailing_Jet_eta");
    DrawPlot(h_LHE_Trailing_Jet_mass,"h_LHE_Trailing_Jet_mass");
    DrawPlot(h_LHE_Trailing_Jet_phi,"h_LHE_Trailing_Jet_phi");
    DrawPlot(h_LHE_Trailing_Jet_pt,"h_LHE_Trailing_Jet_pt");

    DrawPlot(h_LHE_Jet_invariantMass,"h_LHE_Jet_invariantMass");    
    DrawPlot(h_LHE_Jet_etaSeparation,"h_LHE_Jet_etaSeparation");

    DrawPlot(h_LHE_Leading_Jet_eta, h_Matched_Leading_Jet_eta, h_Leading_Jet_eta, "LHE Jet", "Matched Jet", "Jet", "h_Leading_Jet_eta_compare");
    DrawPlot(h_LHE_Leading_Jet_mass, h_Matched_Leading_Jet_mass, h_Leading_Jet_mass, "LHE Jet", "Matched Jet", "Jet", "h_Leading_Jet_mass_compare");
    DrawPlot(h_LHE_Leading_Jet_phi, h_Matched_Leading_Jet_phi, h_Leading_Jet_phi, "LHE Jet", "Matched Jet", "Jet", "h_Leading_Jet_phi_compare");
    DrawPlot(h_LHE_Leading_Jet_pt, h_Matched_Leading_Jet_pt, h_Leading_Jet_pt, "LHE Jet", "Matched Jet", "Jet", "h_Leading_Jet_pt_compare");

    DrawPlot(h_LHE_Trailing_Jet_eta, h_Matched_Trailing_Jet_eta, h_Trailing_Jet_eta, "LHE Jet", "Matched Jet", "Jet", "h_Trailing_Jet_eta_compare");
    DrawPlot(h_LHE_Trailing_Jet_mass, h_Matched_Trailing_Jet_mass, h_Trailing_Jet_mass, "LHE Jet", "Matched Jet", "Jet", "h_Trailing_Jet_mass_compare");
    DrawPlot(h_LHE_Trailing_Jet_phi, h_Matched_Trailing_Jet_phi, h_Trailing_Jet_phi, "LHE Jet", "Matched Jet", "Jet", "h_Trailing_Jet_phi_compare");
    DrawPlot(h_LHE_Trailing_Jet_pt, h_Matched_Trailing_Jet_pt, h_Trailing_Jet_pt, "LHE Jet", "Matched Jet", "Jet", "h_Trailing_Jet_pt_compare");

    DrawPlot(h_LHE_Jet_invariantMass, h_Matched_Jet_invariantMass, h_Jet_invariantMass, "LHE Jet", "Matched Jet", "Jet", "h_Jet_invariantMass_compare");
    DrawPlot(h_LHE_Jet_etaSeparation, h_Matched_Jet_etaSeparation, h_Jet_etaSeparation, "LHE Jet", "Matched Jet", "Jet", "h_Jet_etaSeparation_compare");
 
}

void InitTree(TString infileName,float weight=1){

    cout<<cur_time()<<"\tProcessing Tree...\n";
    TFile* infile = new TFile(infileName);
    EventTree=(TTree*)gDirectory->Get("Events");

    EventTree->SetBranchAddress("run", &run);
    EventTree->SetBranchAddress("luminosityBlock", &luminosityBlock);
    EventTree->SetBranchAddress("event", &event);
    EventTree->SetBranchAddress("HTXS_Higgs_pt", &HTXS_Higgs_pt);
    EventTree->SetBranchAddress("HTXS_Higgs_y", &HTXS_Higgs_y);
    EventTree->SetBranchAddress("HTXS_stage1_1_cat_pTjet25GeV", &HTXS_stage1_1_cat_pTjet25GeV);
    EventTree->SetBranchAddress("HTXS_stage1_1_cat_pTjet30GeV", &HTXS_stage1_1_cat_pTjet30GeV);
    EventTree->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet25GeV", &HTXS_stage1_1_fine_cat_pTjet25GeV);
    EventTree->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet30GeV", &HTXS_stage1_1_fine_cat_pTjet30GeV);
    EventTree->SetBranchAddress("HTXS_stage_0", &HTXS_stage_0);
    EventTree->SetBranchAddress("HTXS_stage_1_pTjet25", &HTXS_stage_1_pTjet25);
    EventTree->SetBranchAddress("HTXS_stage_1_pTjet30", &HTXS_stage_1_pTjet30);
    EventTree->SetBranchAddress("HTXS_njets25", &HTXS_njets25);
    EventTree->SetBranchAddress("HTXS_njets30", &HTXS_njets30);
    EventTree->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2);
    EventTree->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB);
    EventTree->SetBranchAddress("CaloMET_phi", &CaloMET_phi);
    EventTree->SetBranchAddress("CaloMET_pt", &CaloMET_pt);
    EventTree->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt);
    EventTree->SetBranchAddress("ChsMET_phi", &ChsMET_phi);
    EventTree->SetBranchAddress("ChsMET_pt", &ChsMET_pt);
    EventTree->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt);
    EventTree->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet);
    EventTree->SetBranchAddress("CorrT1METJet_area", &CorrT1METJet_area);
    EventTree->SetBranchAddress("CorrT1METJet_eta", &CorrT1METJet_eta);
    EventTree->SetBranchAddress("CorrT1METJet_muonSubtrFactor", &CorrT1METJet_muonSubtrFactor);
    EventTree->SetBranchAddress("CorrT1METJet_phi", &CorrT1METJet_phi);
    EventTree->SetBranchAddress("CorrT1METJet_rawPt", &CorrT1METJet_rawPt);
    EventTree->SetBranchAddress("nElectron", &nElectron);
    EventTree->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC);
    EventTree->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt);
    EventTree->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt);
    EventTree->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt);
    EventTree->SetBranchAddress("Electron_dr03TkSumPtHEEP", &Electron_dr03TkSumPtHEEP);
    EventTree->SetBranchAddress("Electron_dxy", &Electron_dxy);
    EventTree->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr);
    EventTree->SetBranchAddress("Electron_dz", &Electron_dz);
    EventTree->SetBranchAddress("Electron_dzErr", &Electron_dzErr);
    EventTree->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv);
    EventTree->SetBranchAddress("Electron_energyErr", &Electron_energyErr);
    EventTree->SetBranchAddress("Electron_eta", &Electron_eta);
    EventTree->SetBranchAddress("Electron_hoe", &Electron_hoe);
    EventTree->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    EventTree->SetBranchAddress("Electron_jetPtRelv2", &Electron_jetPtRelv2);
    EventTree->SetBranchAddress("Electron_jetRelIso", &Electron_jetRelIso);
    EventTree->SetBranchAddress("Electron_mass", &Electron_mass);
    EventTree->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all);
    EventTree->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg);
    EventTree->SetBranchAddress("Electron_mvaFall17V1Iso", &Electron_mvaFall17V1Iso);
    EventTree->SetBranchAddress("Electron_mvaFall17V1noIso", &Electron_mvaFall17V1noIso);
    EventTree->SetBranchAddress("Electron_mvaFall17V2Iso", &Electron_mvaFall17V2Iso);
    EventTree->SetBranchAddress("Electron_mvaFall17V2noIso", &Electron_mvaFall17V2noIso);
    EventTree->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all);
    EventTree->SetBranchAddress("Electron_pfRelIso03_chg", &Electron_pfRelIso03_chg);
    EventTree->SetBranchAddress("Electron_phi", &Electron_phi);
    EventTree->SetBranchAddress("Electron_pt", &Electron_pt);
    EventTree->SetBranchAddress("Electron_r9", &Electron_r9);
    EventTree->SetBranchAddress("Electron_sieie", &Electron_sieie);
    EventTree->SetBranchAddress("Electron_sip3d", &Electron_sip3d);
    EventTree->SetBranchAddress("Electron_mvaTTH", &Electron_mvaTTH);
    EventTree->SetBranchAddress("Electron_charge", &Electron_charge);
    EventTree->SetBranchAddress("Electron_cutBased", &Electron_cutBased);
    EventTree->SetBranchAddress("Electron_cutBased_Fall17_V1", &Electron_cutBased_Fall17_V1);
    EventTree->SetBranchAddress("Electron_jetIdx", &Electron_jetIdx);
    EventTree->SetBranchAddress("Electron_pdgId", &Electron_pdgId);
    EventTree->SetBranchAddress("Electron_photonIdx", &Electron_photonIdx);
    EventTree->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge);
    EventTree->SetBranchAddress("Electron_vidNestedWPBitmap", &Electron_vidNestedWPBitmap);
    EventTree->SetBranchAddress("Electron_convVeto", &Electron_convVeto);
    EventTree->SetBranchAddress("Electron_cutBased_HEEP", &Electron_cutBased_HEEP);
    EventTree->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand);
    EventTree->SetBranchAddress("Electron_lostHits", &Electron_lostHits);
    EventTree->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", &Electron_mvaFall17V1Iso_WP80);
    EventTree->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", &Electron_mvaFall17V1Iso_WP90);
    EventTree->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", &Electron_mvaFall17V1Iso_WPL);
    EventTree->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", &Electron_mvaFall17V1noIso_WP80);
    EventTree->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", &Electron_mvaFall17V1noIso_WP90);
    EventTree->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", &Electron_mvaFall17V1noIso_WPL);
    EventTree->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", &Electron_mvaFall17V2Iso_WP80);
    EventTree->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90);
    EventTree->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", &Electron_mvaFall17V2Iso_WPL);
    EventTree->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", &Electron_mvaFall17V2noIso_WP80);
    EventTree->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", &Electron_mvaFall17V2noIso_WP90);
    EventTree->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", &Electron_mvaFall17V2noIso_WPL);
    EventTree->SetBranchAddress("Electron_seedGain", &Electron_seedGain);
    EventTree->SetBranchAddress("nFatJet", &nFatJet);
    EventTree->SetBranchAddress("FatJet_area", &FatJet_area);
    EventTree->SetBranchAddress("FatJet_btagCMVA", &FatJet_btagCMVA);
    EventTree->SetBranchAddress("FatJet_btagCSVV2", &FatJet_btagCSVV2);
    EventTree->SetBranchAddress("FatJet_btagDDBvL", &FatJet_btagDDBvL);
    EventTree->SetBranchAddress("FatJet_btagDDCvB", &FatJet_btagDDCvB);
    EventTree->SetBranchAddress("FatJet_btagDDCvL", &FatJet_btagDDCvL);
    EventTree->SetBranchAddress("FatJet_btagDeepB", &FatJet_btagDeepB);
    EventTree->SetBranchAddress("FatJet_btagHbb", &FatJet_btagHbb);
    EventTree->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", &FatJet_deepTagMD_H4qvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", &FatJet_deepTagMD_HbbvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_TvsQCD", &FatJet_deepTagMD_TvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_WvsQCD", &FatJet_deepTagMD_WvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", &FatJet_deepTagMD_ZHbbvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", &FatJet_deepTagMD_ZHccvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", &FatJet_deepTagMD_ZbbvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", &FatJet_deepTagMD_ZvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTagMD_bbvsLight", &FatJet_deepTagMD_bbvsLight);
    EventTree->SetBranchAddress("FatJet_deepTagMD_ccvsLight", &FatJet_deepTagMD_ccvsLight);
    EventTree->SetBranchAddress("FatJet_deepTag_H", &FatJet_deepTag_H);
    EventTree->SetBranchAddress("FatJet_deepTag_QCD", &FatJet_deepTag_QCD);
    EventTree->SetBranchAddress("FatJet_deepTag_QCDothers", &FatJet_deepTag_QCDothers);
    EventTree->SetBranchAddress("FatJet_deepTag_TvsQCD", &FatJet_deepTag_TvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTag_WvsQCD", &FatJet_deepTag_WvsQCD);
    EventTree->SetBranchAddress("FatJet_deepTag_ZvsQCD", &FatJet_deepTag_ZvsQCD);
    EventTree->SetBranchAddress("FatJet_eta", &FatJet_eta);
    EventTree->SetBranchAddress("FatJet_mass", &FatJet_mass);
    EventTree->SetBranchAddress("FatJet_msoftdrop", &FatJet_msoftdrop);
    EventTree->SetBranchAddress("FatJet_n2b1", &FatJet_n2b1);
    EventTree->SetBranchAddress("FatJet_n3b1", &FatJet_n3b1);
    EventTree->SetBranchAddress("FatJet_phi", &FatJet_phi);
    EventTree->SetBranchAddress("FatJet_pt", &FatJet_pt);
    EventTree->SetBranchAddress("FatJet_rawFactor", &FatJet_rawFactor);
    EventTree->SetBranchAddress("FatJet_tau1", &FatJet_tau1);
    EventTree->SetBranchAddress("FatJet_tau2", &FatJet_tau2);
    EventTree->SetBranchAddress("FatJet_tau3", &FatJet_tau3);
    EventTree->SetBranchAddress("FatJet_tau4", &FatJet_tau4);
    EventTree->SetBranchAddress("FatJet_jetId", &FatJet_jetId);
    EventTree->SetBranchAddress("FatJet_subJetIdx1", &FatJet_subJetIdx1);
    EventTree->SetBranchAddress("FatJet_subJetIdx2", &FatJet_subJetIdx2);
    EventTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    EventTree->SetBranchAddress("GenJetAK8_eta", &GenJetAK8_eta);
    EventTree->SetBranchAddress("GenJetAK8_mass", &GenJetAK8_mass);
    EventTree->SetBranchAddress("GenJetAK8_phi", &GenJetAK8_phi);
    EventTree->SetBranchAddress("GenJetAK8_pt", &GenJetAK8_pt);
    EventTree->SetBranchAddress("nGenJet", &nGenJet);
    EventTree->SetBranchAddress("GenJet_eta", &GenJet_eta);
    EventTree->SetBranchAddress("GenJet_mass", &GenJet_mass);
    EventTree->SetBranchAddress("GenJet_phi", &GenJet_phi);
    EventTree->SetBranchAddress("GenJet_pt", &GenJet_pt);
    EventTree->SetBranchAddress("nGenPart", &nGenPart);
    EventTree->SetBranchAddress("GenPart_eta", &GenPart_eta);
    EventTree->SetBranchAddress("GenPart_mass", &GenPart_mass);
    EventTree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    EventTree->SetBranchAddress("GenPart_pt", &GenPart_pt);
    EventTree->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    EventTree->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    EventTree->SetBranchAddress("GenPart_status", &GenPart_status);
    EventTree->SetBranchAddress("GenPart_statusFlags", &GenPart_statusFlags);
    EventTree->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    EventTree->SetBranchAddress("SubGenJetAK8_eta", &SubGenJetAK8_eta);
    EventTree->SetBranchAddress("SubGenJetAK8_mass", &SubGenJetAK8_mass);
    EventTree->SetBranchAddress("SubGenJetAK8_phi", &SubGenJetAK8_phi);
    EventTree->SetBranchAddress("SubGenJetAK8_pt", &SubGenJetAK8_pt);
    EventTree->SetBranchAddress("Generator_binvar", &Generator_binvar);
    EventTree->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF);
    EventTree->SetBranchAddress("Generator_weight", &Generator_weight);
    EventTree->SetBranchAddress("Generator_x1", &Generator_x1);
    EventTree->SetBranchAddress("Generator_x2", &Generator_x2);
    EventTree->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1);
    EventTree->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2);
    EventTree->SetBranchAddress("Generator_id1", &Generator_id1);
    EventTree->SetBranchAddress("Generator_id2", &Generator_id2);
    EventTree->SetBranchAddress("nGenVisTau", &nGenVisTau);
    EventTree->SetBranchAddress("GenVisTau_eta", &GenVisTau_eta);
    EventTree->SetBranchAddress("GenVisTau_mass", &GenVisTau_mass);
    EventTree->SetBranchAddress("GenVisTau_phi", &GenVisTau_phi);
    EventTree->SetBranchAddress("GenVisTau_pt", &GenVisTau_pt);
    EventTree->SetBranchAddress("GenVisTau_charge", &GenVisTau_charge);
    EventTree->SetBranchAddress("GenVisTau_genPartIdxMother", &GenVisTau_genPartIdxMother);
    EventTree->SetBranchAddress("GenVisTau_status", &GenVisTau_status);
    EventTree->SetBranchAddress("genWeight", &genWeight);
    EventTree->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP);
    EventTree->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight);
    EventTree->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
    EventTree->SetBranchAddress("nLHEReweightingWeight", &nLHEReweightingWeight);
    EventTree->SetBranchAddress("LHEReweightingWeight", &LHEReweightingWeight);
    EventTree->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight);
    EventTree->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
    EventTree->SetBranchAddress("nPSWeight", &nPSWeight);
    EventTree->SetBranchAddress("PSWeight", &PSWeight);
    EventTree->SetBranchAddress("nIsoTrack", &nIsoTrack);
    EventTree->SetBranchAddress("IsoTrack_dxy", &IsoTrack_dxy);
    EventTree->SetBranchAddress("IsoTrack_dz", &IsoTrack_dz);
    EventTree->SetBranchAddress("IsoTrack_eta", &IsoTrack_eta);
    EventTree->SetBranchAddress("IsoTrack_pfRelIso03_all", &IsoTrack_pfRelIso03_all);
    EventTree->SetBranchAddress("IsoTrack_pfRelIso03_chg", &IsoTrack_pfRelIso03_chg);
    EventTree->SetBranchAddress("IsoTrack_phi", &IsoTrack_phi);
    EventTree->SetBranchAddress("IsoTrack_pt", &IsoTrack_pt);
    EventTree->SetBranchAddress("IsoTrack_miniPFRelIso_all", &IsoTrack_miniPFRelIso_all);
    EventTree->SetBranchAddress("IsoTrack_miniPFRelIso_chg", &IsoTrack_miniPFRelIso_chg);
    EventTree->SetBranchAddress("IsoTrack_fromPV", &IsoTrack_fromPV);
    EventTree->SetBranchAddress("IsoTrack_pdgId", &IsoTrack_pdgId);
    EventTree->SetBranchAddress("IsoTrack_isHighPurityTrack", &IsoTrack_isHighPurityTrack);
    EventTree->SetBranchAddress("IsoTrack_isPFcand", &IsoTrack_isPFcand);
    EventTree->SetBranchAddress("IsoTrack_isFromLostTrack", &IsoTrack_isFromLostTrack);
    EventTree->SetBranchAddress("nJet", &nJet);
    EventTree->SetBranchAddress("Jet_area", &Jet_area);
    EventTree->SetBranchAddress("Jet_btagCMVA", &Jet_btagCMVA);
    EventTree->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2);
    EventTree->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);
    EventTree->SetBranchAddress("Jet_btagDeepC", &Jet_btagDeepC);
    EventTree->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    EventTree->SetBranchAddress("Jet_btagDeepFlavC", &Jet_btagDeepFlavC);
    EventTree->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF);
    EventTree->SetBranchAddress("Jet_chHEF", &Jet_chHEF);
    EventTree->SetBranchAddress("Jet_eta", &Jet_eta);
    EventTree->SetBranchAddress("Jet_jercCHF", &Jet_jercCHF);
    EventTree->SetBranchAddress("Jet_jercCHPUF", &Jet_jercCHPUF);
    EventTree->SetBranchAddress("Jet_mass", &Jet_mass);
    EventTree->SetBranchAddress("Jet_muEF", &Jet_muEF);
    EventTree->SetBranchAddress("Jet_muonSubtrFactor", &Jet_muonSubtrFactor);
    EventTree->SetBranchAddress("Jet_neEmEF", &Jet_neEmEF);
    EventTree->SetBranchAddress("Jet_neHEF", &Jet_neHEF);
    EventTree->SetBranchAddress("Jet_phi", &Jet_phi);
    EventTree->SetBranchAddress("Jet_pt", &Jet_pt);
    EventTree->SetBranchAddress("Jet_qgl", &Jet_qgl);
    EventTree->SetBranchAddress("Jet_rawFactor", &Jet_rawFactor);
    EventTree->SetBranchAddress("Jet_bRegCorr", &Jet_bRegCorr);
    EventTree->SetBranchAddress("Jet_bRegRes", &Jet_bRegRes);
    EventTree->SetBranchAddress("Jet_electronIdx1", &Jet_electronIdx1);
    EventTree->SetBranchAddress("Jet_electronIdx2", &Jet_electronIdx2);
    EventTree->SetBranchAddress("Jet_jetId", &Jet_jetId);
    EventTree->SetBranchAddress("Jet_muonIdx1", &Jet_muonIdx1);
    EventTree->SetBranchAddress("Jet_muonIdx2", &Jet_muonIdx2);
    EventTree->SetBranchAddress("Jet_nConstituents", &Jet_nConstituents);
    EventTree->SetBranchAddress("Jet_nElectrons", &Jet_nElectrons);
    EventTree->SetBranchAddress("Jet_nMuons", &Jet_nMuons);
    EventTree->SetBranchAddress("Jet_puId", &Jet_puId);
    EventTree->SetBranchAddress("LHE_HT", &LHE_HT);
    EventTree->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming);
    EventTree->SetBranchAddress("LHE_Vpt", &LHE_Vpt);
    EventTree->SetBranchAddress("LHE_Njets", &LHE_Njets);
    EventTree->SetBranchAddress("LHE_Nb", &LHE_Nb);
    EventTree->SetBranchAddress("LHE_Nc", &LHE_Nc);
    EventTree->SetBranchAddress("LHE_Nuds", &LHE_Nuds);
    EventTree->SetBranchAddress("LHE_Nglu", &LHE_Nglu);
    EventTree->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO);
    EventTree->SetBranchAddress("LHE_NpLO", &LHE_NpLO);
    EventTree->SetBranchAddress("nLHEPart", &nLHEPart);
    EventTree->SetBranchAddress("LHEPart_pt", &LHEPart_pt);
    EventTree->SetBranchAddress("LHEPart_eta", &LHEPart_eta);
    EventTree->SetBranchAddress("LHEPart_phi", &LHEPart_phi);
    EventTree->SetBranchAddress("LHEPart_mass", &LHEPart_mass);
    EventTree->SetBranchAddress("LHEPart_pdgId", &LHEPart_pdgId);
    EventTree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    EventTree->SetBranchAddress("GenMET_pt", &GenMET_pt);
    EventTree->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX);
    EventTree->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY);
    EventTree->SetBranchAddress("MET_covXX", &MET_covXX);
    EventTree->SetBranchAddress("MET_covXY", &MET_covXY);
    EventTree->SetBranchAddress("MET_covYY", &MET_covYY);
    EventTree->SetBranchAddress("MET_phi", &MET_phi);
    EventTree->SetBranchAddress("MET_pt", &MET_pt);
    EventTree->SetBranchAddress("MET_significance", &MET_significance);
    EventTree->SetBranchAddress("MET_sumEt", &MET_sumEt);
    EventTree->SetBranchAddress("nMuon", &nMuon);
    EventTree->SetBranchAddress("Muon_dxy", &Muon_dxy);
    EventTree->SetBranchAddress("Muon_dxyErr", &Muon_dxyErr);
    EventTree->SetBranchAddress("Muon_dz", &Muon_dz);
    EventTree->SetBranchAddress("Muon_dzErr", &Muon_dzErr);
    EventTree->SetBranchAddress("Muon_eta", &Muon_eta);
    EventTree->SetBranchAddress("Muon_ip3d", &Muon_ip3d);
    EventTree->SetBranchAddress("Muon_jetPtRelv2", &Muon_jetPtRelv2);
    EventTree->SetBranchAddress("Muon_jetRelIso", &Muon_jetRelIso);
    EventTree->SetBranchAddress("Muon_mass", &Muon_mass);
    EventTree->SetBranchAddress("Muon_miniPFRelIso_all", &Muon_miniPFRelIso_all);
    EventTree->SetBranchAddress("Muon_miniPFRelIso_chg", &Muon_miniPFRelIso_chg);
    EventTree->SetBranchAddress("Muon_pfRelIso03_all", &Muon_pfRelIso03_all);
    EventTree->SetBranchAddress("Muon_pfRelIso03_chg", &Muon_pfRelIso03_chg);
    EventTree->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    EventTree->SetBranchAddress("Muon_phi", &Muon_phi);
    EventTree->SetBranchAddress("Muon_pt", &Muon_pt);
    EventTree->SetBranchAddress("Muon_ptErr", &Muon_ptErr);
    EventTree->SetBranchAddress("Muon_segmentComp", &Muon_segmentComp);
    EventTree->SetBranchAddress("Muon_sip3d", &Muon_sip3d);
    EventTree->SetBranchAddress("Muon_softMva", &Muon_softMva);
    EventTree->SetBranchAddress("Muon_tkRelIso", &Muon_tkRelIso);
    EventTree->SetBranchAddress("Muon_tunepRelPt", &Muon_tunepRelPt);
    EventTree->SetBranchAddress("Muon_mvaLowPt", &Muon_mvaLowPt);
    EventTree->SetBranchAddress("Muon_mvaTTH", &Muon_mvaTTH);
    EventTree->SetBranchAddress("Muon_charge", &Muon_charge);
    EventTree->SetBranchAddress("Muon_jetIdx", &Muon_jetIdx);
    EventTree->SetBranchAddress("Muon_nStations", &Muon_nStations);
    EventTree->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
    EventTree->SetBranchAddress("Muon_pdgId", &Muon_pdgId);
    EventTree->SetBranchAddress("Muon_tightCharge", &Muon_tightCharge);
    EventTree->SetBranchAddress("Muon_highPtId", &Muon_highPtId);
    EventTree->SetBranchAddress("Muon_inTimeMuon", &Muon_inTimeMuon);
    EventTree->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal);
    EventTree->SetBranchAddress("Muon_isPFcand", &Muon_isPFcand);
    EventTree->SetBranchAddress("Muon_isTracker", &Muon_isTracker);
    EventTree->SetBranchAddress("Muon_looseId", &Muon_looseId);
    EventTree->SetBranchAddress("Muon_mediumId", &Muon_mediumId);
    EventTree->SetBranchAddress("Muon_mediumPromptId", &Muon_mediumPromptId);
    EventTree->SetBranchAddress("Muon_miniIsoId", &Muon_miniIsoId);
    EventTree->SetBranchAddress("Muon_multiIsoId", &Muon_multiIsoId);
    EventTree->SetBranchAddress("Muon_mvaId", &Muon_mvaId);
    EventTree->SetBranchAddress("Muon_pfIsoId", &Muon_pfIsoId);
    EventTree->SetBranchAddress("Muon_puppiIsoId", &Muon_puppiIsoId);
    EventTree->SetBranchAddress("Muon_softId", &Muon_softId);
    EventTree->SetBranchAddress("Muon_softMvaId", &Muon_softMvaId);
    EventTree->SetBranchAddress("Muon_tightId", &Muon_tightId);
    EventTree->SetBranchAddress("Muon_tkIsoId", &Muon_tkIsoId);
    EventTree->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
    EventTree->SetBranchAddress("nPhoton", &nPhoton);
    EventTree->SetBranchAddress("Photon_energyErr", &Photon_energyErr);
    EventTree->SetBranchAddress("Photon_eta", &Photon_eta);
    EventTree->SetBranchAddress("Photon_hoe", &Photon_hoe);
    EventTree->SetBranchAddress("Photon_mass", &Photon_mass);
    EventTree->SetBranchAddress("Photon_mvaID", &Photon_mvaID);
    EventTree->SetBranchAddress("Photon_mvaIDV1", &Photon_mvaIDV1);
    EventTree->SetBranchAddress("Photon_pfRelIso03_all", &Photon_pfRelIso03_all);
    EventTree->SetBranchAddress("Photon_pfRelIso03_chg", &Photon_pfRelIso03_chg);
    EventTree->SetBranchAddress("Photon_phi", &Photon_phi);
    EventTree->SetBranchAddress("Photon_pt", &Photon_pt);
    EventTree->SetBranchAddress("Photon_r9", &Photon_r9);
    EventTree->SetBranchAddress("Photon_sieie", &Photon_sieie);
    EventTree->SetBranchAddress("Photon_charge", &Photon_charge);
    EventTree->SetBranchAddress("Photon_cutBasedBitmap", &Photon_cutBasedBitmap);
    EventTree->SetBranchAddress("Photon_cutBasedV1Bitmap", &Photon_cutBasedV1Bitmap);
    EventTree->SetBranchAddress("Photon_electronIdx", &Photon_electronIdx);
    EventTree->SetBranchAddress("Photon_jetIdx", &Photon_jetIdx);
    EventTree->SetBranchAddress("Photon_pdgId", &Photon_pdgId);
    EventTree->SetBranchAddress("Photon_vidNestedWPBitmap", &Photon_vidNestedWPBitmap);
    EventTree->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto);
    EventTree->SetBranchAddress("Photon_isScEtaEB", &Photon_isScEtaEB);
    EventTree->SetBranchAddress("Photon_isScEtaEE", &Photon_isScEtaEE);
    EventTree->SetBranchAddress("Photon_mvaID_WP80", &Photon_mvaID_WP80);
    EventTree->SetBranchAddress("Photon_mvaID_WP90", &Photon_mvaID_WP90);
    EventTree->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed);
    EventTree->SetBranchAddress("Photon_seedGain", &Photon_seedGain);
    EventTree->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);
    EventTree->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity);
    EventTree->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity);
    EventTree->SetBranchAddress("Pileup_nPU", &Pileup_nPU);
    EventTree->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT);
    EventTree->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT);
    EventTree->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi);
    EventTree->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    EventTree->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt);
    EventTree->SetBranchAddress("RawMET_phi", &RawMET_phi);
    EventTree->SetBranchAddress("RawMET_pt", &RawMET_pt);
    EventTree->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt);
    EventTree->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll);
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral);
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo);
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp);
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral);
    EventTree->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton);
    EventTree->SetBranchAddress("GenDressedLepton_eta", &GenDressedLepton_eta);
    EventTree->SetBranchAddress("GenDressedLepton_mass", &GenDressedLepton_mass);
    EventTree->SetBranchAddress("GenDressedLepton_phi", &GenDressedLepton_phi);
    EventTree->SetBranchAddress("GenDressedLepton_pt", &GenDressedLepton_pt);
    EventTree->SetBranchAddress("GenDressedLepton_pdgId", &GenDressedLepton_pdgId);
    EventTree->SetBranchAddress("GenDressedLepton_hasTauAnc", &GenDressedLepton_hasTauAnc);
    EventTree->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet);
    EventTree->SetBranchAddress("SoftActivityJet_eta", &SoftActivityJet_eta);
    EventTree->SetBranchAddress("SoftActivityJet_phi", &SoftActivityJet_phi);
    EventTree->SetBranchAddress("SoftActivityJet_pt", &SoftActivityJet_pt);
    EventTree->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT);
    EventTree->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10);
    EventTree->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2);
    EventTree->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5);
    EventTree->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10);
    EventTree->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2);
    EventTree->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5);
    EventTree->SetBranchAddress("nSubJet", &nSubJet);
    EventTree->SetBranchAddress("SubJet_btagCMVA", &SubJet_btagCMVA);
    EventTree->SetBranchAddress("SubJet_btagCSVV2", &SubJet_btagCSVV2);
    EventTree->SetBranchAddress("SubJet_btagDeepB", &SubJet_btagDeepB);
    EventTree->SetBranchAddress("SubJet_eta", &SubJet_eta);
    EventTree->SetBranchAddress("SubJet_mass", &SubJet_mass);
    EventTree->SetBranchAddress("SubJet_n2b1", &SubJet_n2b1);
    EventTree->SetBranchAddress("SubJet_n3b1", &SubJet_n3b1);
    EventTree->SetBranchAddress("SubJet_phi", &SubJet_phi);
    EventTree->SetBranchAddress("SubJet_pt", &SubJet_pt);
    EventTree->SetBranchAddress("SubJet_rawFactor", &SubJet_rawFactor);
    EventTree->SetBranchAddress("SubJet_tau1", &SubJet_tau1);
    EventTree->SetBranchAddress("SubJet_tau2", &SubJet_tau2);
    EventTree->SetBranchAddress("SubJet_tau3", &SubJet_tau3);
    EventTree->SetBranchAddress("SubJet_tau4", &SubJet_tau4);
    EventTree->SetBranchAddress("nTau", &nTau);
    EventTree->SetBranchAddress("Tau_chargedIso", &Tau_chargedIso);
    EventTree->SetBranchAddress("Tau_dxy", &Tau_dxy);
    EventTree->SetBranchAddress("Tau_dz", &Tau_dz);
    EventTree->SetBranchAddress("Tau_eta", &Tau_eta);
    EventTree->SetBranchAddress("Tau_leadTkDeltaEta", &Tau_leadTkDeltaEta);
    EventTree->SetBranchAddress("Tau_leadTkDeltaPhi", &Tau_leadTkDeltaPhi);
    EventTree->SetBranchAddress("Tau_leadTkPtOverTauPt", &Tau_leadTkPtOverTauPt);
    EventTree->SetBranchAddress("Tau_mass", &Tau_mass);
    EventTree->SetBranchAddress("Tau_neutralIso", &Tau_neutralIso);
    EventTree->SetBranchAddress("Tau_phi", &Tau_phi);
    EventTree->SetBranchAddress("Tau_photonsOutsideSignalCone", &Tau_photonsOutsideSignalCone);
    EventTree->SetBranchAddress("Tau_pt", &Tau_pt);
    EventTree->SetBranchAddress("Tau_puCorr", &Tau_puCorr);
    EventTree->SetBranchAddress("Tau_rawAntiEle", &Tau_rawAntiEle);
    EventTree->SetBranchAddress("Tau_rawAntiEle2018", &Tau_rawAntiEle2018);
    EventTree->SetBranchAddress("Tau_rawDeepTau2017v2VSe", &Tau_rawDeepTau2017v2VSe);
    EventTree->SetBranchAddress("Tau_rawDeepTau2017v2VSjet", &Tau_rawDeepTau2017v2VSjet);
    EventTree->SetBranchAddress("Tau_rawDeepTau2017v2VSmu", &Tau_rawDeepTau2017v2VSmu);
    EventTree->SetBranchAddress("Tau_rawIso", &Tau_rawIso);
    EventTree->SetBranchAddress("Tau_rawIsodR03", &Tau_rawIsodR03);
    EventTree->SetBranchAddress("Tau_rawMVAnewDM2017v2", &Tau_rawMVAnewDM2017v2);
    EventTree->SetBranchAddress("Tau_rawMVAoldDM", &Tau_rawMVAoldDM);
    EventTree->SetBranchAddress("Tau_rawMVAoldDM2017v1", &Tau_rawMVAoldDM2017v1);
    EventTree->SetBranchAddress("Tau_rawMVAoldDM2017v2", &Tau_rawMVAoldDM2017v2);
    EventTree->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", &Tau_rawMVAoldDMdR032017v2);
    EventTree->SetBranchAddress("Tau_charge", &Tau_charge);
    EventTree->SetBranchAddress("Tau_decayMode", &Tau_decayMode);
    EventTree->SetBranchAddress("Tau_jetIdx", &Tau_jetIdx);
    EventTree->SetBranchAddress("Tau_rawAntiEleCat", &Tau_rawAntiEleCat);
    EventTree->SetBranchAddress("Tau_rawAntiEleCat2018", &Tau_rawAntiEleCat2018);
    EventTree->SetBranchAddress("Tau_idAntiEle", &Tau_idAntiEle);
    EventTree->SetBranchAddress("Tau_idAntiEle2018", &Tau_idAntiEle2018);
    EventTree->SetBranchAddress("Tau_idAntiMu", &Tau_idAntiMu);
    EventTree->SetBranchAddress("Tau_idDecayMode", &Tau_idDecayMode);
    EventTree->SetBranchAddress("Tau_idDecayModeNewDMs", &Tau_idDecayModeNewDMs);
    EventTree->SetBranchAddress("Tau_idDeepTau2017v2VSe", &Tau_idDeepTau2017v2VSe);
    EventTree->SetBranchAddress("Tau_idDeepTau2017v2VSjet", &Tau_idDeepTau2017v2VSjet);
    EventTree->SetBranchAddress("Tau_idDeepTau2017v2VSmu", &Tau_idDeepTau2017v2VSmu);
    EventTree->SetBranchAddress("Tau_idMVAnewDM2017v2", &Tau_idMVAnewDM2017v2);
    EventTree->SetBranchAddress("Tau_idMVAoldDM", &Tau_idMVAoldDM);
    EventTree->SetBranchAddress("Tau_idMVAoldDM2017v1", &Tau_idMVAoldDM2017v1);
    EventTree->SetBranchAddress("Tau_idMVAoldDM2017v2", &Tau_idMVAoldDM2017v2);
    EventTree->SetBranchAddress("Tau_idMVAoldDMdR032017v2", &Tau_idMVAoldDMdR032017v2);
    EventTree->SetBranchAddress("TkMET_phi", &TkMET_phi);
    EventTree->SetBranchAddress("TkMET_pt", &TkMET_pt);
    EventTree->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt);
    EventTree->SetBranchAddress("nTrigObj", &nTrigObj);
    EventTree->SetBranchAddress("TrigObj_pt", &TrigObj_pt);
    EventTree->SetBranchAddress("TrigObj_eta", &TrigObj_eta);
    EventTree->SetBranchAddress("TrigObj_phi", &TrigObj_phi);
    EventTree->SetBranchAddress("TrigObj_l1pt", &TrigObj_l1pt);
    EventTree->SetBranchAddress("TrigObj_l1pt_2", &TrigObj_l1pt_2);
    EventTree->SetBranchAddress("TrigObj_l2pt", &TrigObj_l2pt);
    EventTree->SetBranchAddress("TrigObj_id", &TrigObj_id);
    EventTree->SetBranchAddress("TrigObj_l1iso", &TrigObj_l1iso);
    EventTree->SetBranchAddress("TrigObj_l1charge", &TrigObj_l1charge);
    EventTree->SetBranchAddress("TrigObj_filterBits", &TrigObj_filterBits);
    EventTree->SetBranchAddress("genTtbarId", &genTtbarId);
    EventTree->SetBranchAddress("nOtherPV", &nOtherPV);
    EventTree->SetBranchAddress("OtherPV_z", &OtherPV_z);
    EventTree->SetBranchAddress("PV_ndof", &PV_ndof);
    EventTree->SetBranchAddress("PV_x", &PV_x);
    EventTree->SetBranchAddress("PV_y", &PV_y);
    EventTree->SetBranchAddress("PV_z", &PV_z);
    EventTree->SetBranchAddress("PV_chi2", &PV_chi2);
    EventTree->SetBranchAddress("PV_score", &PV_score);
    EventTree->SetBranchAddress("PV_npvs", &PV_npvs);
    EventTree->SetBranchAddress("PV_npvsGood", &PV_npvsGood);
    EventTree->SetBranchAddress("nSV", &nSV);
    EventTree->SetBranchAddress("SV_dlen", &SV_dlen);
    EventTree->SetBranchAddress("SV_dlenSig", &SV_dlenSig);
    EventTree->SetBranchAddress("SV_pAngle", &SV_pAngle);
    EventTree->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    EventTree->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
    EventTree->SetBranchAddress("GenJetAK8_partonFlavour", &GenJetAK8_partonFlavour);
    EventTree->SetBranchAddress("GenJetAK8_hadronFlavour", &GenJetAK8_hadronFlavour);
    EventTree->SetBranchAddress("GenJet_partonFlavour", &GenJet_partonFlavour);
    EventTree->SetBranchAddress("GenJet_hadronFlavour", &GenJet_hadronFlavour);
    EventTree->SetBranchAddress("Jet_genJetIdx", &Jet_genJetIdx);
    EventTree->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
    EventTree->SetBranchAddress("Jet_partonFlavour", &Jet_partonFlavour);
    EventTree->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
    EventTree->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
    EventTree->SetBranchAddress("Photon_genPartIdx", &Photon_genPartIdx);
    EventTree->SetBranchAddress("Photon_genPartFlav", &Photon_genPartFlav);
    EventTree->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi);
    EventTree->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt);
    EventTree->SetBranchAddress("Electron_cleanmask", &Electron_cleanmask);
    EventTree->SetBranchAddress("Jet_cleanmask", &Jet_cleanmask);
    EventTree->SetBranchAddress("Muon_cleanmask", &Muon_cleanmask);
    EventTree->SetBranchAddress("Photon_cleanmask", &Photon_cleanmask);
    EventTree->SetBranchAddress("Tau_cleanmask", &Tau_cleanmask);
    EventTree->SetBranchAddress("SV_chi2", &SV_chi2);
    EventTree->SetBranchAddress("SV_eta", &SV_eta);
    EventTree->SetBranchAddress("SV_mass", &SV_mass);
    EventTree->SetBranchAddress("SV_ndof", &SV_ndof);
    EventTree->SetBranchAddress("SV_phi", &SV_phi);
    EventTree->SetBranchAddress("SV_pt", &SV_pt);
    EventTree->SetBranchAddress("SV_x", &SV_x);
    EventTree->SetBranchAddress("SV_y", &SV_y);
    EventTree->SetBranchAddress("SV_z", &SV_z);
    EventTree->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx);
    EventTree->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav);
    EventTree->SetBranchAddress("L1simulation_step", &L1simulation_step);
    EventTree->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath);
    EventTree->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50);
    EventTree->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID);
    EventTree->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID);
    EventTree->SetBranchAddress("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon);
    EventTree->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon);
    EventTree->SetBranchAddress("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", &HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon);
    EventTree->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
    EventTree->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
    EventTree->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27);
    EventTree->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs);
    EventTree->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi);
    EventTree->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu);
    EventTree->SetBranchAddress("HLT_DoubleMu3_TkMu_DsTau3Mu", &HLT_DoubleMu3_TkMu_DsTau3Mu);
    EventTree->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350);
    EventTree->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40);
    EventTree->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu3_L1SingleMu5orSingleMu7", &HLT_Mu3_L1SingleMu5orSingleMu7);
    EventTree->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL);
    EventTree->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70);
    EventTree->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85);
    EventTree->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele15_WPLoose_Gsf", &HLT_Ele15_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele17_WPLoose_Gsf", &HLT_Ele17_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
    EventTree->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT);
    EventTree->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_HT450_Beamspot", &HLT_HT450_Beamspot);
    EventTree->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot);
    EventTree->SetBranchAddress("HLT_ZeroBias_Beamspot", &HLT_ZeroBias_Beamspot);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1);
    EventTree->SetBranchAddress("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20);
    EventTree->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1);
    EventTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
    EventTree->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30);
    EventTree->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX);
    EventTree->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18);
    EventTree->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25);
    EventTree->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10);
    EventTree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX);
    EventTree->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L2Mu50", &HLT_L2Mu50);
    EventTree->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha", &HLT_L2Mu23NoVtx_2Cha);
    EventTree->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha_CosmicSeed", &HLT_L2Mu23NoVtx_2Cha_CosmicSeed);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha", &HLT_DoubleL2Mu23NoVtx_2Cha);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched", &HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha", &HLT_DoubleL2Mu25NoVtx_2Cha);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched", &HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched);
    EventTree->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
    EventTree->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
    EventTree->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
    EventTree->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia);
    EventTree->SetBranchAddress("HLT_Mu30_TkMu0_Psi", &HLT_Mu30_TkMu0_Psi);
    EventTree->SetBranchAddress("HLT_Mu30_TkMu0_Upsilon", &HLT_Mu30_TkMu0_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi);
    EventTree->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi);
    EventTree->SetBranchAddress("HLT_Mu12", &HLT_Mu12);
    EventTree->SetBranchAddress("HLT_Mu15", &HLT_Mu15);
    EventTree->SetBranchAddress("HLT_Mu20", &HLT_Mu20);
    EventTree->SetBranchAddress("HLT_Mu27", &HLT_Mu27);
    EventTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
    EventTree->SetBranchAddress("HLT_Mu55", &HLT_Mu55);
    EventTree->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100);
    EventTree->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100);
    EventTree->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40);
    EventTree->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60);
    EventTree->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80);
    EventTree->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140);
    EventTree->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200);
    EventTree->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260);
    EventTree->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320);
    EventTree->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400);
    EventTree->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500);
    EventTree->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC);
    EventTree->SetBranchAddress("HLT_AK8PFJet15", &HLT_AK8PFJet15);
    EventTree->SetBranchAddress("HLT_AK8PFJet25", &HLT_AK8PFJet25);
    EventTree->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40);
    EventTree->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60);
    EventTree->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80);
    EventTree->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140);
    EventTree->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200);
    EventTree->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260);
    EventTree->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320);
    EventTree->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400);
    EventTree->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450);
    EventTree->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500);
    EventTree->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550);
    EventTree->SetBranchAddress("HLT_PFJet15", &HLT_PFJet15);
    EventTree->SetBranchAddress("HLT_PFJet25", &HLT_PFJet25);
    EventTree->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40);
    EventTree->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60);
    EventTree->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80);
    EventTree->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140);
    EventTree->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200);
    EventTree->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260);
    EventTree->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320);
    EventTree->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400);
    EventTree->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450);
    EventTree->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500);
    EventTree->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550);
    EventTree->SetBranchAddress("HLT_PFJetFwd15", &HLT_PFJetFwd15);
    EventTree->SetBranchAddress("HLT_PFJetFwd25", &HLT_PFJetFwd25);
    EventTree->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40);
    EventTree->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60);
    EventTree->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80);
    EventTree->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140);
    EventTree->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200);
    EventTree->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260);
    EventTree->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320);
    EventTree->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400);
    EventTree->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450);
    EventTree->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd15", &HLT_AK8PFJetFwd15);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd25", &HLT_AK8PFJetFwd25);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450);
    EventTree->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500);
    EventTree->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180);
    EventTree->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250);
    EventTree->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370);
    EventTree->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430);
    EventTree->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510);
    EventTree->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590);
    EventTree->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680);
    EventTree->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780);
    EventTree->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890);
    EventTree->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050);
    EventTree->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight);
    EventTree->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight);
    EventTree->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight);
    EventTree->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight", &HLT_PFHT700_PFMET95_PFMHT95_IDTight);
    EventTree->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight);
    EventTree->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight", &HLT_PFHT800_PFMET85_PFMHT85_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1);
    EventTree->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1);
    EventTree->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1);
    EventTree->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1);
    EventTree->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1);
    EventTree->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight);
    EventTree->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds);
    EventTree->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90);
    EventTree->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned);
    EventTree->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned);
    EventTree->SetBranchAddress("HLT_CaloMET100_NotCleaned", &HLT_CaloMET100_NotCleaned);
    EventTree->SetBranchAddress("HLT_CaloMET110_NotCleaned", &HLT_CaloMET110_NotCleaned);
    EventTree->SetBranchAddress("HLT_CaloMET250_NotCleaned", &HLT_CaloMET250_NotCleaned);
    EventTree->SetBranchAddress("HLT_CaloMET70_HBHECleaned", &HLT_CaloMET70_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET80_HBHECleaned", &HLT_CaloMET80_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET90_HBHECleaned", &HLT_CaloMET90_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET100_HBHECleaned", &HLT_CaloMET100_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET250_HBHECleaned", &HLT_CaloMET250_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET300_HBHECleaned", &HLT_CaloMET300_HBHECleaned);
    EventTree->SetBranchAddress("HLT_CaloMET350_HBHECleaned", &HLT_CaloMET350_HBHECleaned);
    EventTree->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned);
    EventTree->SetBranchAddress("HLT_PFMET200_HBHECleaned", &HLT_PFMET200_HBHECleaned);
    EventTree->SetBranchAddress("HLT_PFMET250_HBHECleaned", &HLT_PFMET250_HBHECleaned);
    EventTree->SetBranchAddress("HLT_PFMET300_HBHECleaned", &HLT_PFMET300_HBHECleaned);
    EventTree->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned", &HLT_PFMET200_HBHE_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50);
    EventTree->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50);
    EventTree->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets40_CaloBTagDeepCSV_p71", &HLT_DoublePFJets40_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets100_CaloBTagDeepCSV_p71", &HLT_DoublePFJets100_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets200_CaloBTagDeepCSV_p71", &HLT_DoublePFJets200_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets350_CaloBTagDeepCSV_p71", &HLT_DoublePFJets350_CaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
    EventTree->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
    EventTree->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
    EventTree->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5", &HLT_BTagMu_AK8Jet170_DoubleMu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5_noalgo", &HLT_BTagMu_AK4DiJet20_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5_noalgo", &HLT_BTagMu_AK4DiJet40_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5_noalgo", &HLT_BTagMu_AK4DiJet70_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5_noalgo", &HLT_BTagMu_AK4DiJet110_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5_noalgo", &HLT_BTagMu_AK4DiJet170_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5_noalgo", &HLT_BTagMu_AK4Jet300_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5_noalgo", &HLT_BTagMu_AK8DiJet170_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo", &HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5_noalgo", &HLT_BTagMu_AK8Jet300_Mu5_noalgo);
    EventTree->SetBranchAddress("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu12_DoublePhoton20", &HLT_Mu12_DoublePhoton20);
    EventTree->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2", &HLT_TriplePhoton_20_20_20_CaloIdLV2);
    EventTree->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL);
    EventTree->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2", &HLT_TriplePhoton_30_30_10_CaloIdLV2);
    EventTree->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL);
    EventTree->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL);
    EventTree->SetBranchAddress("HLT_Photon20", &HLT_Photon20);
    EventTree->SetBranchAddress("HLT_Photon33", &HLT_Photon33);
    EventTree->SetBranchAddress("HLT_Photon50", &HLT_Photon50);
    EventTree->SetBranchAddress("HLT_Photon75", &HLT_Photon75);
    EventTree->SetBranchAddress("HLT_Photon90", &HLT_Photon90);
    EventTree->SetBranchAddress("HLT_Photon120", &HLT_Photon120);
    EventTree->SetBranchAddress("HLT_Photon150", &HLT_Photon150);
    EventTree->SetBranchAddress("HLT_Photon175", &HLT_Photon175);
    EventTree->SetBranchAddress("HLT_Photon200", &HLT_Photon200);
    EventTree->SetBranchAddress("HLT_Photon100EB_TightID_TightIso", &HLT_Photon100EB_TightID_TightIso);
    EventTree->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso);
    EventTree->SetBranchAddress("HLT_Photon120EB_TightID_TightIso", &HLT_Photon120EB_TightID_TightIso);
    EventTree->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10);
    EventTree->SetBranchAddress("HLT_Photon100EEHE10", &HLT_Photon100EEHE10);
    EventTree->SetBranchAddress("HLT_Photon100EE_TightID_TightIso", &HLT_Photon100EE_TightID_TightIso);
    EventTree->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3);
    EventTree->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700", &HLT_Photon90_CaloIdL_PFHT700);
    EventTree->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
    EventTree->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
    EventTree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
    EventTree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
    EventTree->SetBranchAddress("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35);
    EventTree->SetBranchAddress("HLT_IsoMu24_TwoProngs35", &HLT_IsoMu24_TwoProngs35);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R);
    EventTree->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
    EventTree->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DZ", &HLT_TripleMu_5_3_3_Mass3p8_DZ);
    EventTree->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ);
    EventTree->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5);
    EventTree->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
    EventTree->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
    EventTree->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
    EventTree->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
    EventTree->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
    EventTree->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
    EventTree->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
    EventTree->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
    EventTree->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing);
    EventTree->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL);
    EventTree->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtxDisplaced", &HLT_DoubleMu33NoFiltersNoVtxDisplaced);
    EventTree->SetBranchAddress("HLT_DoubleMu40NoFiltersNoVtxDisplaced", &HLT_DoubleMu40NoFiltersNoVtxDisplaced);
    EventTree->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4", &HLT_DoubleMu20_7_Mass0to30_L1_DM4);
    EventTree->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &HLT_DoubleMu20_7_Mass0to30_L1_DM4EG);
    EventTree->SetBranchAddress("HLT_HT425", &HLT_HT425);
    EventTree->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT500_DisplacedDijet40_DisplacedTrack", &HLT_HT500_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive);
    EventTree->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive);
    EventTree->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110);
    EventTree->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120);
    EventTree->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130);
    EventTree->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110);
    EventTree->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120);
    EventTree->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130);
    EventTree->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
    EventTree->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
    EventTree->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55);
    EventTree->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450);
    EventTree->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600);
    EventTree->SetBranchAddress("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
    EventTree->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450);
    EventTree->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight);
    EventTree->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight);
    EventTree->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
    EventTree->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls);
    EventTree->SetBranchAddress("HLT_Dimuon12_Upsilon_y1p4", &HLT_Dimuon12_Upsilon_y1p4);
    EventTree->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls);
    EventTree->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime);
    EventTree->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi);
    EventTree->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1);
    EventTree->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1);
    EventTree->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1);
    EventTree->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1);
    EventTree->SetBranchAddress("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", &HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8);
    EventTree->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
    EventTree->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1);
    EventTree->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_Mu8", &HLT_Mu8);
    EventTree->SetBranchAddress("HLT_Mu17", &HLT_Mu17);
    EventTree->SetBranchAddress("HLT_Mu19", &HLT_Mu19);
    EventTree->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId);
    EventTree->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
    EventTree->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele145_CaloIdVT_GsfTrkIdT", &HLT_Ele145_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele200_CaloIdVT_GsfTrkIdT", &HLT_Ele200_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5);
    EventTree->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40);
    EventTree->SetBranchAddress("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94", &HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94);
    EventTree->SetBranchAddress("HLT_PFHT400_SixPFJet32", &HLT_PFHT400_SixPFJet32);
    EventTree->SetBranchAddress("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59", &HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59);
    EventTree->SetBranchAddress("HLT_PFHT450_SixPFJet36", &HLT_PFHT450_SixPFJet36);
    EventTree->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350);
    EventTree->SetBranchAddress("HLT_PFHT350MinPFJet15", &HLT_PFHT350MinPFJet15);
    EventTree->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL", &HLT_Photon60_R9Id90_CaloIdL_IsoL);
    EventTree->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL);
    EventTree->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
    EventTree->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800);
    EventTree->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70);
    EventTree->SetBranchAddress("HLT_Physics", &HLT_Physics);
    EventTree->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0);
    EventTree->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1);
    EventTree->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2);
    EventTree->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3);
    EventTree->SetBranchAddress("HLT_Physics_part4", &HLT_Physics_part4);
    EventTree->SetBranchAddress("HLT_Physics_part5", &HLT_Physics_part5);
    EventTree->SetBranchAddress("HLT_Physics_part6", &HLT_Physics_part6);
    EventTree->SetBranchAddress("HLT_Physics_part7", &HLT_Physics_part7);
    EventTree->SetBranchAddress("HLT_Random", &HLT_Random);
    EventTree->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias);
    EventTree->SetBranchAddress("HLT_ZeroBias_Alignment", &HLT_ZeroBias_Alignment);
    EventTree->SetBranchAddress("HLT_ZeroBias_part0", &HLT_ZeroBias_part0);
    EventTree->SetBranchAddress("HLT_ZeroBias_part1", &HLT_ZeroBias_part1);
    EventTree->SetBranchAddress("HLT_ZeroBias_part2", &HLT_ZeroBias_part2);
    EventTree->SetBranchAddress("HLT_ZeroBias_part3", &HLT_ZeroBias_part3);
    EventTree->SetBranchAddress("HLT_ZeroBias_part4", &HLT_ZeroBias_part4);
    EventTree->SetBranchAddress("HLT_ZeroBias_part5", &HLT_ZeroBias_part5);
    EventTree->SetBranchAddress("HLT_ZeroBias_part6", &HLT_ZeroBias_part6);
    EventTree->SetBranchAddress("HLT_ZeroBias_part7", &HLT_ZeroBias_part7);
    EventTree->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30);
    EventTree->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40);
    EventTree->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50);
    EventTree->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80);
    EventTree->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100);
    EventTree->SetBranchAddress("HLT_AK4CaloJet120", &HLT_AK4CaloJet120);
    EventTree->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30);
    EventTree->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50);
    EventTree->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80);
    EventTree->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100);
    EventTree->SetBranchAddress("HLT_AK4PFJet120", &HLT_AK4PFJet120);
    EventTree->SetBranchAddress("HLT_SinglePhoton10_Eta3p1ForPPRef", &HLT_SinglePhoton10_Eta3p1ForPPRef);
    EventTree->SetBranchAddress("HLT_SinglePhoton20_Eta3p1ForPPRef", &HLT_SinglePhoton20_Eta3p1ForPPRef);
    EventTree->SetBranchAddress("HLT_SinglePhoton30_Eta3p1ForPPRef", &HLT_SinglePhoton30_Eta3p1ForPPRef);
    EventTree->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose);
    EventTree->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose);
    EventTree->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration);
    EventTree->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration);
    EventTree->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus);
    EventTree->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus);
    EventTree->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR);
    EventTree->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
    EventTree->SetBranchAddress("HLT_CDC_L2cosmic_5_er1p0", &HLT_CDC_L2cosmic_5_er1p0);
    EventTree->SetBranchAddress("HLT_CDC_L2cosmic_5p5_er1p0", &HLT_CDC_L2cosmic_5p5_er1p0);
    EventTree->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS);
    EventTree->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym);
    EventTree->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch);
    EventTree->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB);
    EventTree->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap);
    EventTree->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain);
    EventTree->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain);
    EventTree->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1);
    EventTree->SetBranchAddress("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1);
    EventTree->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_Rsq0p35", &HLT_Rsq0p35);
    EventTree->SetBranchAddress("HLT_Rsq0p40", &HLT_Rsq0p40);
    EventTree->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200", &HLT_RsqMR300_Rsq0p09_MR200);
    EventTree->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200", &HLT_RsqMR320_Rsq0p09_MR200);
    EventTree->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200_4jet", &HLT_RsqMR300_Rsq0p09_MR200_4jet);
    EventTree->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200_4jet", &HLT_RsqMR320_Rsq0p09_MR200_4jet);
    EventTree->SetBranchAddress("HLT_IsoMu27_MET90", &HLT_IsoMu27_MET90);
    EventTree->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3);
    EventTree->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
    EventTree->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign);
    EventTree->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ);
    EventTree->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9);
    EventTree->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ);
    EventTree->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign);
    EventTree->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ);
    EventTree->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12);
    EventTree->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ);
    EventTree->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
    EventTree->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
    EventTree->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
    EventTree->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DCA", &HLT_TripleMu_5_3_3_Mass3p8_DCA);
    EventTree->SetBranchAddress("HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
    EventTree->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
    EventTree->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
    EventTree->SetBranchAddress("HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2);
    EventTree->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2);
    EventTree->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2);
    EventTree->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2);
    EventTree->SetBranchAddress("HLT_QuadPFJet98_83_71_15", &HLT_QuadPFJet98_83_71_15);
    EventTree->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15);
    EventTree->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15);
    EventTree->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15);
    EventTree->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17);
    EventTree->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1);
    EventTree->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02);
    EventTree->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2);
    EventTree->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4);
    EventTree->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55);
    EventTree->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto);
    EventTree->SetBranchAddress("HLT_Mu12_IP6_part0", &HLT_Mu12_IP6_part0);
    EventTree->SetBranchAddress("HLT_Mu12_IP6_part1", &HLT_Mu12_IP6_part1);
    EventTree->SetBranchAddress("HLT_Mu12_IP6_part2", &HLT_Mu12_IP6_part2);
    EventTree->SetBranchAddress("HLT_Mu12_IP6_part3", &HLT_Mu12_IP6_part3);
    EventTree->SetBranchAddress("HLT_Mu12_IP6_part4", &HLT_Mu12_IP6_part4);
    EventTree->SetBranchAddress("HLT_Mu9_IP5_part0", &HLT_Mu9_IP5_part0);
    EventTree->SetBranchAddress("HLT_Mu9_IP5_part1", &HLT_Mu9_IP5_part1);
    EventTree->SetBranchAddress("HLT_Mu9_IP5_part2", &HLT_Mu9_IP5_part2);
    EventTree->SetBranchAddress("HLT_Mu9_IP5_part3", &HLT_Mu9_IP5_part3);
    EventTree->SetBranchAddress("HLT_Mu9_IP5_part4", &HLT_Mu9_IP5_part4);
    EventTree->SetBranchAddress("HLT_Mu7_IP4_part0", &HLT_Mu7_IP4_part0);
    EventTree->SetBranchAddress("HLT_Mu7_IP4_part1", &HLT_Mu7_IP4_part1);
    EventTree->SetBranchAddress("HLT_Mu7_IP4_part2", &HLT_Mu7_IP4_part2);
    EventTree->SetBranchAddress("HLT_Mu7_IP4_part3", &HLT_Mu7_IP4_part3);
    EventTree->SetBranchAddress("HLT_Mu7_IP4_part4", &HLT_Mu7_IP4_part4);
    EventTree->SetBranchAddress("HLT_Mu9_IP4_part0", &HLT_Mu9_IP4_part0);
    EventTree->SetBranchAddress("HLT_Mu9_IP4_part1", &HLT_Mu9_IP4_part1);
    EventTree->SetBranchAddress("HLT_Mu9_IP4_part2", &HLT_Mu9_IP4_part2);
    EventTree->SetBranchAddress("HLT_Mu9_IP4_part3", &HLT_Mu9_IP4_part3);
    EventTree->SetBranchAddress("HLT_Mu9_IP4_part4", &HLT_Mu9_IP4_part4);
    EventTree->SetBranchAddress("HLT_Mu8_IP5_part0", &HLT_Mu8_IP5_part0);
    EventTree->SetBranchAddress("HLT_Mu8_IP5_part1", &HLT_Mu8_IP5_part1);
    EventTree->SetBranchAddress("HLT_Mu8_IP5_part2", &HLT_Mu8_IP5_part2);
    EventTree->SetBranchAddress("HLT_Mu8_IP5_part3", &HLT_Mu8_IP5_part3);
    EventTree->SetBranchAddress("HLT_Mu8_IP5_part4", &HLT_Mu8_IP5_part4);
    EventTree->SetBranchAddress("HLT_Mu8_IP6_part0", &HLT_Mu8_IP6_part0);
    EventTree->SetBranchAddress("HLT_Mu8_IP6_part1", &HLT_Mu8_IP6_part1);
    EventTree->SetBranchAddress("HLT_Mu8_IP6_part2", &HLT_Mu8_IP6_part2);
    EventTree->SetBranchAddress("HLT_Mu8_IP6_part3", &HLT_Mu8_IP6_part3);
    EventTree->SetBranchAddress("HLT_Mu8_IP6_part4", &HLT_Mu8_IP6_part4);
    EventTree->SetBranchAddress("HLT_Mu9_IP6_part0", &HLT_Mu9_IP6_part0);
    EventTree->SetBranchAddress("HLT_Mu9_IP6_part1", &HLT_Mu9_IP6_part1);
    EventTree->SetBranchAddress("HLT_Mu9_IP6_part2", &HLT_Mu9_IP6_part2);
    EventTree->SetBranchAddress("HLT_Mu9_IP6_part3", &HLT_Mu9_IP6_part3);
    EventTree->SetBranchAddress("HLT_Mu9_IP6_part4", &HLT_Mu9_IP6_part4);
    EventTree->SetBranchAddress("HLT_Mu8_IP3_part0", &HLT_Mu8_IP3_part0);
    EventTree->SetBranchAddress("HLT_Mu8_IP3_part1", &HLT_Mu8_IP3_part1);
    EventTree->SetBranchAddress("HLT_Mu8_IP3_part2", &HLT_Mu8_IP3_part2);
    EventTree->SetBranchAddress("HLT_Mu8_IP3_part3", &HLT_Mu8_IP3_part3);
    EventTree->SetBranchAddress("HLT_Mu8_IP3_part4", &HLT_Mu8_IP3_part4);
    EventTree->SetBranchAddress("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
    EventTree->SetBranchAddress("HLT_TrkMu6NoFiltersNoVtx", &HLT_TrkMu6NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_TrkMu16NoFiltersNoVtx", &HLT_TrkMu16NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_DoubleTrkMu_16_6_NoFiltersNoVtx", &HLT_DoubleTrkMu_16_6_NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath);
    EventTree->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    EventTree->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
    EventTree->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter);
    EventTree->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter);
    EventTree->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter);
    EventTree->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
    EventTree->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
    EventTree->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter);
    EventTree->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter);
    EventTree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    EventTree->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter);
    EventTree->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
    EventTree->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
    EventTree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    EventTree->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter);
    EventTree->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters);
    EventTree->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter);
    EventTree->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter);
    EventTree->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
    EventTree->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    EventTree->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter);
    EventTree->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter);
    EventTree->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X);
    EventTree->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X);
    EventTree->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters);
    EventTree->SetBranchAddress("Flag_METFilters", &Flag_METFilters);
    EventTree->SetBranchAddress("L1Reco_step", &L1Reco_step);
    EventTree->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue);
    EventTree->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME);
    EventTree->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME);
    EventTree->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME);
    EventTree->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME);
    EventTree->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME);
    EventTree->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME);
    EventTree->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME);
    EventTree->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME);
    EventTree->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME);
    EventTree->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME);
    EventTree->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME);
    EventTree->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus);
    EventTree->SetBranchAddress("L1_BptxOR", &L1_BptxOR);
    EventTree->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus);
    EventTree->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR);
    EventTree->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
    EventTree->SetBranchAddress("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er);
    EventTree->SetBranchAddress("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er);
    EventTree->SetBranchAddress("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er);
    EventTree->SetBranchAddress("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er);
    EventTree->SetBranchAddress("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er);
    EventTree->SetBranchAddress("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5);
    EventTree->SetBranchAddress("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5);
    EventTree->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1);
    EventTree->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1);
    EventTree->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1);
    EventTree->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6);
    EventTree->SetBranchAddress("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5);
    EventTree->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6);
    EventTree->SetBranchAddress("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5);
    EventTree->SetBranchAddress("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5);
    EventTree->SetBranchAddress("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp);
    EventTree->SetBranchAddress("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5);
    EventTree->SetBranchAddress("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
    EventTree->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
    EventTree->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
    EventTree->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28);
    EventTree->SetBranchAddress("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620);
    EventTree->SetBranchAddress("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28);
    EventTree->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ);
    EventTree->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp);
    EventTree->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8);
    EventTree->SetBranchAddress("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
    EventTree->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1);
    EventTree->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1);
    EventTree->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0);
    EventTree->SetBranchAddress("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1);
    EventTree->SetBranchAddress("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ);
    EventTree->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ);
    EventTree->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS);
    EventTree->SetBranchAddress("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8);
    EventTree->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
    EventTree->SetBranchAddress("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ);
    EventTree->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS);
    EventTree->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
    EventTree->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4);
    EventTree->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4);
    EventTree->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4);
    EventTree->SetBranchAddress("L1_DoubleMu10_SQ", &L1_DoubleMu10_SQ);
    EventTree->SetBranchAddress("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1);
    EventTree->SetBranchAddress("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er);
    EventTree->SetBranchAddress("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er);
    EventTree->SetBranchAddress("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8);
    EventTree->SetBranchAddress("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5);
    EventTree->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS);
    EventTree->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2);
    EventTree->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS);
    EventTree->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
    EventTree->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS);
    EventTree->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18);
    EventTree->SetBranchAddress("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3);
    EventTree->SetBranchAddress("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5);
    EventTree->SetBranchAddress("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ);
    EventTree->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5);
    EventTree->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ);
    EventTree->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7);
    EventTree->SetBranchAddress("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1);
    EventTree->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ);
    EventTree->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1);
    EventTree->SetBranchAddress("L1_ETM120", &L1_ETM120);
    EventTree->SetBranchAddress("L1_ETM150", &L1_ETM150);
    EventTree->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100);
    EventTree->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er);
    EventTree->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110);
    EventTree->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er);
    EventTree->SetBranchAddress("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain);
    EventTree->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120);
    EventTree->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er);
    EventTree->SetBranchAddress("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain);
    EventTree->SetBranchAddress("L1_ETMHF130", &L1_ETMHF130);
    EventTree->SetBranchAddress("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er);
    EventTree->SetBranchAddress("L1_ETMHF140", &L1_ETMHF140);
    EventTree->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150);
    EventTree->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er);
    EventTree->SetBranchAddress("L1_ETT1200", &L1_ETT1200);
    EventTree->SetBranchAddress("L1_ETT1600", &L1_ETT1600);
    EventTree->SetBranchAddress("L1_ETT2000", &L1_ETT2000);
    EventTree->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain);
    EventTree->SetBranchAddress("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain);
    EventTree->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain);
    EventTree->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit);
    EventTree->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain);
    EventTree->SetBranchAddress("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig);
    EventTree->SetBranchAddress("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto);
    EventTree->SetBranchAddress("L1_HTT120er", &L1_HTT120er);
    EventTree->SetBranchAddress("L1_HTT160er", &L1_HTT160er);
    EventTree->SetBranchAddress("L1_HTT200er", &L1_HTT200er);
    EventTree->SetBranchAddress("L1_HTT255er", &L1_HTT255er);
    EventTree->SetBranchAddress("L1_HTT280er", &L1_HTT280er);
    EventTree->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p4", &L1_HTT280er_QuadJet_70_55_40_35_er2p4);
    EventTree->SetBranchAddress("L1_HTT320er", &L1_HTT320er);
    EventTree->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4);
    EventTree->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3);
    EventTree->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3);
    EventTree->SetBranchAddress("L1_HTT360er", &L1_HTT360er);
    EventTree->SetBranchAddress("L1_HTT400er", &L1_HTT400er);
    EventTree->SetBranchAddress("L1_HTT450er", &L1_HTT450er);
    EventTree->SetBranchAddress("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40);
    EventTree->SetBranchAddress("L1_IsoEG32er2p5_Mt44", &L1_IsoEG32er2p5_Mt44);
    EventTree->SetBranchAddress("L1_IsoEG32er2p5_Mt48", &L1_IsoEG32er2p5_Mt48);
    EventTree->SetBranchAddress("L1_IsoTau40er2p1_ETMHF100", &L1_IsoTau40er2p1_ETMHF100);
    EventTree->SetBranchAddress("L1_IsoTau40er2p1_ETMHF110", &L1_IsoTau40er2p1_ETMHF110);
    EventTree->SetBranchAddress("L1_IsoTau40er2p1_ETMHF120", &L1_IsoTau40er2p1_ETMHF120);
    EventTree->SetBranchAddress("L1_IsoTau40er2p1_ETMHF90", &L1_IsoTau40er2p1_ETMHF90);
    EventTree->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch);
    EventTree->SetBranchAddress("L1_LastBunchInTrain", &L1_LastBunchInTrain);
    EventTree->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain);
    EventTree->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
    EventTree->SetBranchAddress("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3);
    EventTree->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er);
    EventTree->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
    EventTree->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er);
    EventTree->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3);
    EventTree->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er);
    EventTree->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3);
    EventTree->SetBranchAddress("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er);
    EventTree->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3);
    EventTree->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND);
    EventTree->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
    EventTree->SetBranchAddress("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6);
    EventTree->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
    EventTree->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1);
    EventTree->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1);
    EventTree->SetBranchAddress("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5);
    EventTree->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1);
    EventTree->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1);
    EventTree->SetBranchAddress("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1);
    EventTree->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1);
    EventTree->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1);
    EventTree->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4);
    EventTree->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8);
    EventTree->SetBranchAddress("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4);
    EventTree->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5);
    EventTree->SetBranchAddress("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4);
    EventTree->SetBranchAddress("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4);
    EventTree->SetBranchAddress("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4);
    EventTree->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40);
    EventTree->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50);
    EventTree->SetBranchAddress("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5);
    EventTree->SetBranchAddress("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5);
    EventTree->SetBranchAddress("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5);
    EventTree->SetBranchAddress("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5);
    EventTree->SetBranchAddress("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5);
    EventTree->SetBranchAddress("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5);
    EventTree->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er);
    EventTree->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er);
    EventTree->SetBranchAddress("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5);
    EventTree->SetBranchAddress("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5);
    EventTree->SetBranchAddress("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5);
    EventTree->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR);
    EventTree->SetBranchAddress("L1_QuadJet36er2p5_IsoTau52er2p1", &L1_QuadJet36er2p5_IsoTau52er2p1);
    EventTree->SetBranchAddress("L1_QuadJet60er2p5", &L1_QuadJet60er2p5);
    EventTree->SetBranchAddress("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0);
    EventTree->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0);
    EventTree->SetBranchAddress("L1_QuadMu0_OQ", &L1_QuadMu0_OQ);
    EventTree->SetBranchAddress("L1_QuadMu0_SQ", &L1_QuadMu0_SQ);
    EventTree->SetBranchAddress("L1_SecondBunchInTrain", &L1_SecondBunchInTrain);
    EventTree->SetBranchAddress("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain);
    EventTree->SetBranchAddress("L1_SingleEG10er2p5", &L1_SingleEG10er2p5);
    EventTree->SetBranchAddress("L1_SingleEG15er2p5", &L1_SingleEG15er2p5);
    EventTree->SetBranchAddress("L1_SingleEG26er2p5", &L1_SingleEG26er2p5);
    EventTree->SetBranchAddress("L1_SingleEG34er2p5", &L1_SingleEG34er2p5);
    EventTree->SetBranchAddress("L1_SingleEG36er2p5", &L1_SingleEG36er2p5);
    EventTree->SetBranchAddress("L1_SingleEG38er2p5", &L1_SingleEG38er2p5);
    EventTree->SetBranchAddress("L1_SingleEG40er2p5", &L1_SingleEG40er2p5);
    EventTree->SetBranchAddress("L1_SingleEG42er2p5", &L1_SingleEG42er2p5);
    EventTree->SetBranchAddress("L1_SingleEG45er2p5", &L1_SingleEG45er2p5);
    EventTree->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50);
    EventTree->SetBranchAddress("L1_SingleEG60", &L1_SingleEG60);
    EventTree->SetBranchAddress("L1_SingleEG8er2p5", &L1_SingleEG8er2p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1);
    EventTree->SetBranchAddress("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1);
    EventTree->SetBranchAddress("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1);
    EventTree->SetBranchAddress("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1);
    EventTree->SetBranchAddress("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1);
    EventTree->SetBranchAddress("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5);
    EventTree->SetBranchAddress("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5);
    EventTree->SetBranchAddress("L1_SingleJet10erHE", &L1_SingleJet10erHE);
    EventTree->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120);
    EventTree->SetBranchAddress("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0);
    EventTree->SetBranchAddress("L1_SingleJet120er2p5", &L1_SingleJet120er2p5);
    EventTree->SetBranchAddress("L1_SingleJet12erHE", &L1_SingleJet12erHE);
    EventTree->SetBranchAddress("L1_SingleJet140er2p5", &L1_SingleJet140er2p5);
    EventTree->SetBranchAddress("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80);
    EventTree->SetBranchAddress("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90);
    EventTree->SetBranchAddress("L1_SingleJet160er2p5", &L1_SingleJet160er2p5);
    EventTree->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180);
    EventTree->SetBranchAddress("L1_SingleJet180er2p5", &L1_SingleJet180er2p5);
    EventTree->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200);
    EventTree->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR);
    EventTree->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX);
    EventTree->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35);
    EventTree->SetBranchAddress("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0);
    EventTree->SetBranchAddress("L1_SingleJet35er2p5", &L1_SingleJet35er2p5);
    EventTree->SetBranchAddress("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX);
    EventTree->SetBranchAddress("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX);
    EventTree->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60);
    EventTree->SetBranchAddress("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0);
    EventTree->SetBranchAddress("L1_SingleJet60er2p5", &L1_SingleJet60er2p5);
    EventTree->SetBranchAddress("L1_SingleJet8erHE", &L1_SingleJet8erHE);
    EventTree->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90);
    EventTree->SetBranchAddress("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0);
    EventTree->SetBranchAddress("L1_SingleJet90er2p5", &L1_SingleJet90er2p5);
    EventTree->SetBranchAddress("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5);
    EventTree->SetBranchAddress("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5);
    EventTree->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF);
    EventTree->SetBranchAddress("L1_SingleMu0_DQ", &L1_SingleMu0_DQ);
    EventTree->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF);
    EventTree->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF);
    EventTree->SetBranchAddress("L1_SingleMu10er1p5", &L1_SingleMu10er1p5);
    EventTree->SetBranchAddress("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF);
    EventTree->SetBranchAddress("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF);
    EventTree->SetBranchAddress("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF);
    EventTree->SetBranchAddress("L1_SingleMu12er1p5", &L1_SingleMu12er1p5);
    EventTree->SetBranchAddress("L1_SingleMu14er1p5", &L1_SingleMu14er1p5);
    EventTree->SetBranchAddress("L1_SingleMu15_DQ", &L1_SingleMu15_DQ);
    EventTree->SetBranchAddress("L1_SingleMu16er1p5", &L1_SingleMu16er1p5);
    EventTree->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18);
    EventTree->SetBranchAddress("L1_SingleMu18er1p5", &L1_SingleMu18er1p5);
    EventTree->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20);
    EventTree->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22);
    EventTree->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF);
    EventTree->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF);
    EventTree->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF);
    EventTree->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25);
    EventTree->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3);
    EventTree->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5);
    EventTree->SetBranchAddress("L1_SingleMu6er1p5", &L1_SingleMu6er1p5);
    EventTree->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7);
    EventTree->SetBranchAddress("L1_SingleMu7_DQ", &L1_SingleMu7_DQ);
    EventTree->SetBranchAddress("L1_SingleMu7er1p5", &L1_SingleMu7er1p5);
    EventTree->SetBranchAddress("L1_SingleMu8er1p5", &L1_SingleMu8er1p5);
    EventTree->SetBranchAddress("L1_SingleMu9er1p5", &L1_SingleMu9er1p5);
    EventTree->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics);
    EventTree->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF);
    EventTree->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF);
    EventTree->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF);
    EventTree->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen);
    EventTree->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR);
    EventTree->SetBranchAddress("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX);
    EventTree->SetBranchAddress("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX);
    EventTree->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1);
    EventTree->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1);
    EventTree->SetBranchAddress("L1_TOTEM_1", &L1_TOTEM_1);
    EventTree->SetBranchAddress("L1_TOTEM_2", &L1_TOTEM_2);
    EventTree->SetBranchAddress("L1_TOTEM_3", &L1_TOTEM_3);
    EventTree->SetBranchAddress("L1_TOTEM_4", &L1_TOTEM_4);
    EventTree->SetBranchAddress("L1_TripleEG16er2p5", &L1_TripleEG16er2p5);
    EventTree->SetBranchAddress("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5);
    EventTree->SetBranchAddress("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5);
    EventTree->SetBranchAddress("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5);
    EventTree->SetBranchAddress("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5);
    EventTree->SetBranchAddress("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5);
    EventTree->SetBranchAddress("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5);
    EventTree->SetBranchAddress("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5);
    EventTree->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0);
    EventTree->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ);
    EventTree->SetBranchAddress("L1_TripleMu0_SQ", &L1_TripleMu0_SQ);
    EventTree->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3);
    EventTree->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ);
    EventTree->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ);
    EventTree->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
    EventTree->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
    EventTree->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3);
    EventTree->SetBranchAddress("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ);
    EventTree->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5);
    EventTree->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
    EventTree->SetBranchAddress("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
    EventTree->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
    EventTree->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3);
    EventTree->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus);
    EventTree->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus);
    EventTree->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias);
    EventTree->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy);
    EventTree->SetBranchAddress("L1_UnprefireableEvent", &L1_UnprefireableEvent);

    EvMax=EventTree->GetEntries();

    cout<<cur_time()<<"\tTree successfully processed!\n";

    EventLoop();
    infile->Close();
    //score_infile->Close();

    return;
}

void nanoAODAnalyzer_test(){
//main program

    //string inputFile = "root://eoscms.cern.ch://eos/cms/store/user/lzygala/HVV/ttHTobb/4A5BE1BE-DA13-E811-BEB0-AC1F6B1AEFFC.root";

    InitHistograms();
    //InitTree(inputFile.c_str()); 

    const char* inDir[3];
    char* dir_xroot[3];
    char* dir[3];
    void* dirp[3];

    inDir[0] = "/eos/cms/store/user/lzygala/HVV/NanoAOD/pp_hw+w-jj_SMNP0_EFT_VBFCut/chw_0.9999000";
    dir_xroot[0] = "root://eoscms.cern.ch///eos/cms/store/user/lzygala/HVV/NanoAOD/pp_hw+w-jj_SMNP0_EFT_VBFCut/chw_0.9999000";
   

    for(int j=0; j<1; j++){


        dir[j] = gSystem->ExpandPathName(inDir[j]);
        dirp[j] = gSystem->OpenDirectory(dir[j]);

        const char* ext = ".root";

        const char* entry;
        const char* filename[10000];
        Int_t n = 0;
        TString str;

        while((entry = (char*)gSystem->GetDirEntry(dirp[j]))) {
            str = entry;
            if(str.EndsWith(ext))
                filename[n++] = gSystem->ConcatFileName(dir_xroot[j], entry);
        }

        for (Int_t i = 0; i < n; i++){ //n
            Printf("\tfile -> %s", filename[i]);
            Printf("\n%s\tfile -> %i / %i", (cur_time()).c_str(), i, n);
            InitTree(filename[i]);
            //EventLoop();
        }

        

    }

    SaveHistograms();
}
#include<TApplication.h>
#include<TFile.h>
#include<TMath.h>
#include<TMinuit.h>
#include<TROOT.h>
#include<TSystem.h>
#include<TTree.h>
#include<TVector2.h>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TLegend.h>
#include<TLine.h>
#include<TH2F.h>
#include<TPaveText.h>
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
    int nSamples;
    vector<string> sampleAddresses;
    vector<string> sampleNames;


    vector<TH1F*> h_CaloMET_phi;
    vector<TH1F*> h_CaloMET_pt;
    vector<TH1F*> h_CaloMET_sumEt;

    vector<TH1F*> h_Electron_eta;
    vector<TH1F*> h_Electron_hoe;
    vector<TH1F*> h_Electron_mass;
    vector<TH1F*> h_Electron_phi;
    vector<TH1F*> h_Electron_pt;
    vector<TH1F*> h_Electron_r9;
    vector<TH1F*> h_Electron_sieie;

    vector<TH1F*> h_Muon_eta;
    vector<TH1F*> h_Muon_mass;
    vector<TH1F*> h_Muon_phi;
    vector<TH1F*> h_Muon_pt;

    vector<TH1F*> h_Jet_eta;
    vector<TH1F*> h_Jet_mass;
    vector<TH1F*> h_Jet_phi;
    vector<TH1F*> h_Jet_pt;


/* Tree Values */
    UInt_t          run;
    UInt_t          luminosityBlock;
    ULong64_t       event;
    Float_t         CaloMET_phi;
    Float_t         CaloMET_pt;
    Float_t         CaloMET_sumEt;
    UInt_t          nElectron;
    Float_t         Electron_deltaEtaSC[11];   //[nElectron]
    Float_t         Electron_dr03EcalRecHitSumEt[11];   //[nElectron]
    Float_t         Electron_dr03HcalDepth1TowerSumEt[11];   //[nElectron]
    Float_t         Electron_dr03TkSumPt[11];   //[nElectron]
    Float_t         Electron_dxy[11];   //[nElectron]
    Float_t         Electron_dxyErr[11];   //[nElectron]
    Float_t         Electron_dz[11];   //[nElectron]
    Float_t         Electron_dzErr[11];   //[nElectron]
    Float_t         Electron_eCorr[11];   //[nElectron]
    Float_t         Electron_eInvMinusPInv[11];   //[nElectron]
    Float_t         Electron_energyErr[11];   //[nElectron]
    Float_t         Electron_eta[11];   //[nElectron]
    Float_t         Electron_hoe[11];   //[nElectron]
    Float_t         Electron_ip3d[11];   //[nElectron]
    Float_t         Electron_mass[11];   //[nElectron]
    Float_t         Electron_miniPFRelIso_all[11];   //[nElectron]
    Float_t         Electron_miniPFRelIso_chg[11];   //[nElectron]
    Float_t         Electron_mvaSpring16GP[11];   //[nElectron]
    Float_t         Electron_mvaSpring16HZZ[11];   //[nElectron]
    Float_t         Electron_pfRelIso03_all[11];   //[nElectron]
    Float_t         Electron_pfRelIso03_chg[11];   //[nElectron]
    Float_t         Electron_phi[11];   //[nElectron]
    Float_t         Electron_pt[11];   //[nElectron]
    Float_t         Electron_r9[11];   //[nElectron]
    Float_t         Electron_sieie[11];   //[nElectron]
    Float_t         Electron_sip3d[11];   //[nElectron]
    Float_t         Electron_mvaTTH[11];   //[nElectron]
    Int_t           Electron_charge[11];   //[nElectron]
    Int_t           Electron_cutBased[11];   //[nElectron]
    Int_t           Electron_cutBased_HLTPreSel[11];   //[nElectron]
    Int_t           Electron_jetIdx[11];   //[nElectron]
    Int_t           Electron_pdgId[11];   //[nElectron]
    Int_t           Electron_photonIdx[11];   //[nElectron]
    Int_t           Electron_tightCharge[11];   //[nElectron]
    Int_t           Electron_vidNestedWPBitmap[11];   //[nElectron]
    Bool_t          Electron_convVeto[11];   //[nElectron]
    Bool_t          Electron_cutBased_HEEP[11];   //[nElectron]
    Bool_t          Electron_isPFcand[11];   //[nElectron]
    UChar_t         Electron_lostHits[11];   //[nElectron]
    Bool_t          Electron_mvaSpring16GP_WP80[11];   //[nElectron]
    Bool_t          Electron_mvaSpring16GP_WP90[11];   //[nElectron]
    Bool_t          Electron_mvaSpring16HZZ_WPL[11];   //[nElectron]
    UInt_t          nFatJet;
    Float_t         FatJet_area[8];   //[nFatJet]
    Float_t         FatJet_btagCMVA[8];   //[nFatJet]
    Float_t         FatJet_btagCSVV2[8];   //[nFatJet]
    Float_t         FatJet_btagDeepB[8];   //[nFatJet]
    Float_t         FatJet_btagHbb[8];   //[nFatJet]
    Float_t         FatJet_eta[8];   //[nFatJet]
    Float_t         FatJet_mass[8];   //[nFatJet]
    Float_t         FatJet_msoftdrop[8];   //[nFatJet]
    Float_t         FatJet_msoftdrop_chs[8];   //[nFatJet]
    Float_t         FatJet_n2b1[8];   //[nFatJet]
    Float_t         FatJet_n3b1[8];   //[nFatJet]
    Float_t         FatJet_phi[8];   //[nFatJet]
    Float_t         FatJet_pt[8];   //[nFatJet]
    Float_t         FatJet_tau1[8];   //[nFatJet]
    Float_t         FatJet_tau2[8];   //[nFatJet]
    Float_t         FatJet_tau3[8];   //[nFatJet]
    Float_t         FatJet_tau4[8];   //[nFatJet]
    Int_t           FatJet_jetId[8];   //[nFatJet]
    Int_t           FatJet_subJetIdx1[8];   //[nFatJet]
    Int_t           FatJet_subJetIdx2[8];   //[nFatJet]
    UInt_t          nGenJetAK8;
    Float_t         GenJetAK8_eta[8];   //[nGenJetAK8]
    Float_t         GenJetAK8_mass[8];   //[nGenJetAK8]
    Float_t         GenJetAK8_phi[8];   //[nGenJetAK8]
    Float_t         GenJetAK8_pt[8];   //[nGenJetAK8]
    UInt_t          nGenJet;
    Float_t         GenJet_eta[28];   //[nGenJet]
    Float_t         GenJet_mass[28];   //[nGenJet]
    Float_t         GenJet_phi[28];   //[nGenJet]
    Float_t         GenJet_pt[28];   //[nGenJet]
    UInt_t          nGenPart;
    Float_t         GenPart_eta[234];   //[nGenPart]
    Float_t         GenPart_mass[234];   //[nGenPart]
    Float_t         GenPart_phi[234];   //[nGenPart]
    Float_t         GenPart_pt[234];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[234];   //[nGenPart]
    Int_t           GenPart_pdgId[234];   //[nGenPart]
    Int_t           GenPart_status[234];   //[nGenPart]
    Int_t           GenPart_statusFlags[234];   //[nGenPart]
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
    Float_t         GenVisTau_eta[4];   //[nGenVisTau]
    Float_t         GenVisTau_mass[4];   //[nGenVisTau]
    Float_t         GenVisTau_phi[4];   //[nGenVisTau]
    Float_t         GenVisTau_pt[4];   //[nGenVisTau]
    Int_t           GenVisTau_charge[4];   //[nGenVisTau]
    Int_t           GenVisTau_genPartIdxMother[4];   //[nGenVisTau]
    Int_t           GenVisTau_status[4];   //[nGenVisTau]
    Float_t         genWeight;
    Float_t         LHEWeight_originalXWGTUP;
    UInt_t          nLHEPdfWeight;
    Float_t         LHEPdfWeight[100];   //[nLHEPdfWeight]
    UInt_t          nLHEScaleWeight;
    Float_t         LHEScaleWeight[9];   //[nLHEScaleWeight]
    UInt_t          nJet;
    Float_t         Jet_area[48];   //[nJet]
    Float_t         Jet_btagCMVA[48];   //[nJet]
    Float_t         Jet_btagCSVV2[48];   //[nJet]
    Float_t         Jet_btagDeepB[48];   //[nJet]
    Float_t         Jet_btagDeepC[48];   //[nJet]
    Float_t         Jet_chEmEF[48];   //[nJet]
    Float_t         Jet_chHEF[48];   //[nJet]
    Float_t         Jet_eta[48];   //[nJet]
    Float_t         Jet_mass[48];   //[nJet]
    Float_t         Jet_neEmEF[48];   //[nJet]
    Float_t         Jet_neHEF[48];   //[nJet]
    Float_t         Jet_phi[48];   //[nJet]
    Float_t         Jet_pt[48];   //[nJet]
    Float_t         Jet_qgl[48];   //[nJet]
    Float_t         Jet_rawFactor[48];   //[nJet]
    Float_t         Jet_bReg[48];   //[nJet]
    Int_t           Jet_electronIdx1[48];   //[nJet]
    Int_t           Jet_electronIdx2[48];   //[nJet]
    Int_t           Jet_jetId[48];   //[nJet]
    Int_t           Jet_muonIdx1[48];   //[nJet]
    Int_t           Jet_muonIdx2[48];   //[nJet]
    Int_t           Jet_nConstituents[48];   //[nJet]
    Int_t           Jet_nElectrons[48];   //[nJet]
    Int_t           Jet_nMuons[48];   //[nJet]
    Int_t           Jet_puId[48];   //[nJet]
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
    Float_t         Muon_dxy[9];   //[nMuon]
    Float_t         Muon_dxyErr[9];   //[nMuon]
    Float_t         Muon_dz[9];   //[nMuon]
    Float_t         Muon_dzErr[9];   //[nMuon]
    Float_t         Muon_eta[9];   //[nMuon]
    Float_t         Muon_ip3d[9];   //[nMuon]
    Float_t         Muon_mass[9];   //[nMuon]
    Float_t         Muon_miniPFRelIso_all[9];   //[nMuon]
    Float_t         Muon_miniPFRelIso_chg[9];   //[nMuon]
    Float_t         Muon_pfRelIso03_all[9];   //[nMuon]
    Float_t         Muon_pfRelIso03_chg[9];   //[nMuon]
    Float_t         Muon_pfRelIso04_all[9];   //[nMuon]
    Float_t         Muon_phi[9];   //[nMuon]
    Float_t         Muon_pt[9];   //[nMuon]
    Float_t         Muon_ptErr[9];   //[nMuon]
    Float_t         Muon_segmentComp[9];   //[nMuon]
    Float_t         Muon_sip3d[9];   //[nMuon]
    Float_t         Muon_mvaTTH[9];   //[nMuon]
    Int_t           Muon_charge[9];   //[nMuon]
    Int_t           Muon_jetIdx[9];   //[nMuon]
    Int_t           Muon_nStations[9];   //[nMuon]
    Int_t           Muon_nTrackerLayers[9];   //[nMuon]
    Int_t           Muon_pdgId[9];   //[nMuon]
    Int_t           Muon_tightCharge[9];   //[nMuon]
    UChar_t         Muon_highPtId[9];   //[nMuon]
    Bool_t          Muon_isPFcand[9];   //[nMuon]
    Bool_t          Muon_mediumId[9];   //[nMuon]
    Bool_t          Muon_softId[9];   //[nMuon]
    Bool_t          Muon_tightId[9];   //[nMuon]
    UInt_t          nPhoton;
    Float_t         Photon_eCorr[13];   //[nPhoton]
    Float_t         Photon_energyErr[13];   //[nPhoton]
    Float_t         Photon_eta[13];   //[nPhoton]
    Float_t         Photon_hoe[13];   //[nPhoton]
    Float_t         Photon_mass[13];   //[nPhoton]
    Float_t         Photon_mvaID[13];   //[nPhoton]
    Float_t         Photon_pfRelIso03_all[13];   //[nPhoton]
    Float_t         Photon_pfRelIso03_chg[13];   //[nPhoton]
    Float_t         Photon_phi[13];   //[nPhoton]
    Float_t         Photon_pt[13];   //[nPhoton]
    Float_t         Photon_r9[13];   //[nPhoton]
    Float_t         Photon_sieie[13];   //[nPhoton]
    Int_t           Photon_charge[13];   //[nPhoton]
    Int_t           Photon_cutBased[13];   //[nPhoton]
    Int_t           Photon_electronIdx[13];   //[nPhoton]
    Int_t           Photon_jetIdx[13];   //[nPhoton]
    Int_t           Photon_pdgId[13];   //[nPhoton]
    Int_t           Photon_vidNestedWPBitmap[13];   //[nPhoton]
    Bool_t          Photon_electronVeto[13];   //[nPhoton]
    Bool_t          Photon_mvaID_WP80[13];   //[nPhoton]
    Bool_t          Photon_mvaID_WP90[13];   //[nPhoton]
    Bool_t          Photon_pixelSeed[13];   //[nPhoton]
    Float_t         Pileup_nTrueInt;
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
    Float_t         fixedGridRhoFastjetCentralCalo;
    Float_t         fixedGridRhoFastjetCentralNeutral;
    UInt_t          nGenDressedLepton;
    Float_t         GenDressedLepton_eta[4];   //[nGenDressedLepton]
    Float_t         GenDressedLepton_mass[4];   //[nGenDressedLepton]
    Float_t         GenDressedLepton_phi[4];   //[nGenDressedLepton]
    Float_t         GenDressedLepton_pt[4];   //[nGenDressedLepton]
    Int_t           GenDressedLepton_pdgId[4];   //[nGenDressedLepton]
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
    Float_t         SubJet_tau1[12];   //[nSubJet]
    Float_t         SubJet_tau2[12];   //[nSubJet]
    Float_t         SubJet_tau3[12];   //[nSubJet]
    Float_t         SubJet_tau4[12];   //[nSubJet]
    UInt_t          nTau;
    Float_t         Tau_chargedIso[7];   //[nTau]
    Float_t         Tau_dxy[7];   //[nTau]
    Float_t         Tau_dz[7];   //[nTau]
    Float_t         Tau_eta[7];   //[nTau]
    Float_t         Tau_footprintCorr[7];   //[nTau]
    Float_t         Tau_leadTkDeltaEta[7];   //[nTau]
    Float_t         Tau_leadTkDeltaPhi[7];   //[nTau]
    Float_t         Tau_leadTkPtOverTauPt[7];   //[nTau]
    Float_t         Tau_mass[7];   //[nTau]
    Float_t         Tau_neutralIso[7];   //[nTau]
    Float_t         Tau_phi[7];   //[nTau]
    Float_t         Tau_photonsOutsideSignalCone[7];   //[nTau]
    Float_t         Tau_pt[7];   //[nTau]
    Float_t         Tau_puCorr[7];   //[nTau]
    Float_t         Tau_rawAntiEle[7];   //[nTau]
    Float_t         Tau_rawIso[7];   //[nTau]
    Float_t         Tau_rawMVAnewDM[7];   //[nTau]
    Float_t         Tau_rawMVAoldDM[7];   //[nTau]
    Float_t         Tau_rawMVAoldDMdR03[7];   //[nTau]
    Int_t           Tau_charge[7];   //[nTau]
    Int_t           Tau_decayMode[7];   //[nTau]
    Int_t           Tau_jetIdx[7];   //[nTau]
    Int_t           Tau_rawAntiEleCat[7];   //[nTau]
    UChar_t         Tau_idAntiEle[7];   //[nTau]
    UChar_t         Tau_idAntiMu[7];   //[nTau]
    Bool_t          Tau_idDecayMode[7];   //[nTau]
    Bool_t          Tau_idDecayModeNewDMs[7];   //[nTau]
    UChar_t         Tau_idMVAnewDM[7];   //[nTau]
    UChar_t         Tau_idMVAoldDM[7];   //[nTau]
    UChar_t         Tau_idMVAoldDMdR03[7];   //[nTau]
    Float_t         TkMET_phi;
    Float_t         TkMET_pt;
    Float_t         TkMET_sumEt;
    UInt_t          nTrigObj;
    Float_t         TrigObj_pt[81];   //[nTrigObj]
    Float_t         TrigObj_eta[81];   //[nTrigObj]
    Float_t         TrigObj_phi[81];   //[nTrigObj]
    Float_t         TrigObj_l1pt[81];   //[nTrigObj]
    Float_t         TrigObj_l1pt_2[81];   //[nTrigObj]
    Float_t         TrigObj_l2pt[81];   //[nTrigObj]
    Int_t           TrigObj_id[81];   //[nTrigObj]
    Int_t           TrigObj_l1iso[81];   //[nTrigObj]
    Int_t           TrigObj_l1charge[81];   //[nTrigObj]
    Int_t           TrigObj_filterBits[81];   //[nTrigObj]
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
    Float_t         SV_dlen[18];   //[nSV]
    Float_t         SV_dlenSig[18];   //[nSV]
    Float_t         SV_pAngle[18];   //[nSV]
    Int_t           Electron_genPartIdx[11];   //[nElectron]
    UChar_t         Electron_genPartFlav[11];   //[nElectron]
    Int_t           GenJetAK8_partonFlavour[8];   //[nGenJetAK8]
    UChar_t         GenJetAK8_hadronFlavour[8];   //[nGenJetAK8]
    Int_t           GenJet_partonFlavour[28];   //[nGenJet]
    UChar_t         GenJet_hadronFlavour[28];   //[nGenJet]
    Int_t           Jet_genJetIdx[48];   //[nJet]
    Int_t           Jet_hadronFlavour[48];   //[nJet]
    Int_t           Jet_partonFlavour[48];   //[nJet]
    Int_t           Muon_genPartIdx[9];   //[nMuon]
    UChar_t         Muon_genPartFlav[9];   //[nMuon]
    Int_t           Photon_genPartIdx[13];   //[nPhoton]
    UChar_t         Photon_genPartFlav[13];   //[nPhoton]
    Float_t         MET_fiducialGenPhi;
    Float_t         MET_fiducialGenPt;
    UChar_t         Electron_cleanmask[11];   //[nElectron]
    UChar_t         Jet_cleanmask[48];   //[nJet]
    UChar_t         Muon_cleanmask[9];   //[nMuon]
    UChar_t         Photon_cleanmask[13];   //[nPhoton]
    UChar_t         Tau_cleanmask[7];   //[nTau]
    Float_t         SV_chi2[18];   //[nSV]
    Float_t         SV_eta[18];   //[nSV]
    Float_t         SV_mass[18];   //[nSV]
    Float_t         SV_ndof[18];   //[nSV]
    Float_t         SV_phi[18];   //[nSV]
    Float_t         SV_pt[18];   //[nSV]
    Float_t         SV_x[18];   //[nSV]
    Float_t         SV_y[18];   //[nSV]
    Float_t         SV_z[18];   //[nSV]
    Int_t           Tau_genPartIdx[7];   //[nTau]
    UChar_t         Tau_genPartFlav[7];   //[nTau]
    Bool_t          L1simulation_step;
    Bool_t          HLTriggerFirstPath;
    Bool_t          HLT_AK8PFJet360_TrimMass30;
    Bool_t          HLT_AK8PFJet400_TrimMass30;
    Bool_t          HLT_AK8PFHT750_TrimMass50;
    Bool_t          HLT_AK8PFHT800_TrimMass50;
    Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;
    Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;
    Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;
    Bool_t          HLT_AK8DiPFJet300_200_TrimMass30;
    Bool_t          HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
    Bool_t          HLT_AK8PFHT650_TrimR0p1PT0p03Mass50;
    Bool_t          HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20;
    Bool_t          HLT_AK8DiPFJet280_200_TrimMass30;
    Bool_t          HLT_AK8DiPFJet250_200_TrimMass30;
    Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;
    Bool_t          HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;
    Bool_t          HLT_CaloJet260;
    Bool_t          HLT_CaloJet500_NoJetID;
    Bool_t          HLT_Dimuon13_PsiPrime;
    Bool_t          HLT_Dimuon13_Upsilon;
    Bool_t          HLT_Dimuon20_Jpsi;
    Bool_t          HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_DoubleEle25_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_DoubleEle33_CaloIdL;
    Bool_t          HLT_DoubleEle33_CaloIdL_MW;
    Bool_t          HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
    Bool_t          HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1;
    Bool_t          HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1;
    Bool_t          HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg;
    Bool_t          HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1;
    Bool_t          HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_DoubleMu33NoFiltersNoVtx;
    Bool_t          HLT_DoubleMu38NoFiltersNoVtx;
    Bool_t          HLT_DoubleMu23NoFiltersNoVtxDisplaced;
    Bool_t          HLT_DoubleMu28NoFiltersNoVtxDisplaced;
    Bool_t          HLT_DoubleMu0;
    Bool_t          HLT_DoubleMu4_3_Bs;
    Bool_t          HLT_DoubleMu4_3_Jpsi_Displaced;
    Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
    Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
    Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
    Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
    Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
    Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
    Bool_t          HLT_Mu7p5_Track2_Jpsi;
    Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
    Bool_t          HLT_Mu7p5_Track7_Jpsi;
    Bool_t          HLT_Mu7p5_Track2_Upsilon;
    Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
    Bool_t          HLT_Mu7p5_Track7_Upsilon;
    Bool_t          HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing;
    Bool_t          HLT_Dimuon0er16_Jpsi_NoVertexing;
    Bool_t          HLT_Dimuon6_Jpsi_NoVertexing;
    Bool_t          HLT_Photon150;
    Bool_t          HLT_Photon90_CaloIdL_HT300;
    Bool_t          HLT_HT250_CaloMET70;
    Bool_t          HLT_DoublePhoton60;
    Bool_t          HLT_DoublePhoton85;
    Bool_t          HLT_Ele17_Ele8_Gsf;
    Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28;
    Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29;
    Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_Ele23_WPLoose_Gsf;
    Bool_t          HLT_Ele23_WPLoose_Gsf_WHbbBoost;
    Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20;
    Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30;
    Bool_t          HLT_Ele25_WPTight_Gsf;
    Bool_t          HLT_Ele25_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele25_eta2p1_WPTight_Gsf;
    Bool_t          HLT_Ele27_WPLoose_Gsf;
    Bool_t          HLT_Ele27_WPLoose_Gsf_WHbbBoost;
    Bool_t          HLT_Ele27_WPTight_Gsf;
    Bool_t          HLT_Ele27_WPTight_Gsf_L1JetTauSeeded;
    Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_Ele27_eta2p1_WPTight_Gsf;
    Bool_t          HLT_Ele30_WPTight_Gsf;
    Bool_t          HLT_Ele30_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf;
    Bool_t          HLT_Ele32_WPTight_Gsf;
    Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf;
    Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_Ele32_eta2p1_WPTight_Gsf;
    Bool_t          HLT_Ele35_WPLoose_Gsf;
    Bool_t          HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;
    Bool_t          HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_Ele45_WPLoose_Gsf;
    Bool_t          HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded;
    Bool_t          HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    Bool_t          HLT_Ele105_CaloIdVT_GsfTrkIdT;
    Bool_t          HLT_Ele30WP60_SC4_Mass55;
    Bool_t          HLT_Ele30WP60_Ele8_Mass55;
    Bool_t          HLT_HT200;
    Bool_t          HLT_HT275;
    Bool_t          HLT_HT325;
    Bool_t          HLT_HT425;
    Bool_t          HLT_HT575;
    Bool_t          HLT_HT410to430;
    Bool_t          HLT_HT430to450;
    Bool_t          HLT_HT450to470;
    Bool_t          HLT_HT470to500;
    Bool_t          HLT_HT500to550;
    Bool_t          HLT_HT550to650;
    Bool_t          HLT_HT650;
    Bool_t          HLT_Mu16_eta2p1_MET30;
    Bool_t          HLT_IsoMu16_eta2p1_MET30;
    Bool_t          HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1;
    Bool_t          HLT_IsoMu17_eta2p1;
    Bool_t          HLT_IsoMu17_eta2p1_LooseIsoPFTau20;
    Bool_t          HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_DoubleIsoMu17_eta2p1;
    Bool_t          HLT_DoubleIsoMu17_eta2p1_noDzCut;
    Bool_t          HLT_IsoMu18;
    Bool_t          HLT_IsoMu19_eta2p1_LooseIsoPFTau20;
    Bool_t          HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20;
    Bool_t          HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu20;
    Bool_t          HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1;
    Bool_t          HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1;
    Bool_t          HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;
    Bool_t          HLT_IsoMu22;
    Bool_t          HLT_IsoMu22_eta2p1;
    Bool_t          HLT_IsoMu24;
    Bool_t          HLT_IsoMu27;
    Bool_t          HLT_IsoTkMu18;
    Bool_t          HLT_IsoTkMu20;
    Bool_t          HLT_IsoTkMu22;
    Bool_t          HLT_IsoTkMu22_eta2p1;
    Bool_t          HLT_IsoTkMu24;
    Bool_t          HLT_IsoTkMu27;
    Bool_t          HLT_JetE30_NoBPTX3BX;
    Bool_t          HLT_JetE30_NoBPTX;
    Bool_t          HLT_JetE50_NoBPTX3BX;
    Bool_t          HLT_JetE70_NoBPTX3BX;
    Bool_t          HLT_L1SingleMu18;
    Bool_t          HLT_L2Mu10;
    Bool_t          HLT_L1SingleMuOpen;
    Bool_t          HLT_L1SingleMuOpen_DT;
    Bool_t          HLT_L2DoubleMu23_NoVertex;
    Bool_t          HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10;
    Bool_t          HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10;
    Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
    Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
    Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
    Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
    Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1;
    Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80;
    Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90;
    Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110;
    Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120;
    Bool_t          HLT_PFTau120_eta2p1;
    Bool_t          HLT_PFTau140_eta2p1;
    Bool_t          HLT_VLooseIsoPFTau120_Trk50_eta2p1;
    Bool_t          HLT_VLooseIsoPFTau140_Trk50_eta2p1;
    Bool_t          HLT_Mu17_Mu8;
    Bool_t          HLT_Mu17_Mu8_DZ;
    Bool_t          HLT_Mu17_Mu8_SameSign;
    Bool_t          HLT_Mu17_Mu8_SameSign_DZ;
    Bool_t          HLT_Mu20_Mu10;
    Bool_t          HLT_Mu20_Mu10_DZ;
    Bool_t          HLT_Mu20_Mu10_SameSign;
    Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
    Bool_t          HLT_Mu17_TkMu8_DZ;
    Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
    Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
    Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
    Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
    Bool_t          HLT_Mu25_TkMu0_dEta18_Onia;
    Bool_t          HLT_Mu27_TkMu8;
    Bool_t          HLT_Mu30_TkMu11;
    Bool_t          HLT_Mu30_eta2p1_PFJet150_PFJet50;
    Bool_t          HLT_Mu40_TkMu11;
    Bool_t          HLT_Mu40_eta2p1_PFJet200_PFJet50;
    Bool_t          HLT_Mu20;
    Bool_t          HLT_TkMu17;
    Bool_t          HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
    Bool_t          HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
    Bool_t          HLT_TkMu20;
    Bool_t          HLT_Mu24_eta2p1;
    Bool_t          HLT_TkMu24_eta2p1;
    Bool_t          HLT_Mu27;
    Bool_t          HLT_TkMu27;
    Bool_t          HLT_Mu45_eta2p1;
    Bool_t          HLT_Mu50;
    Bool_t          HLT_TkMu50;
    Bool_t          HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL;
    Bool_t          HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL;
    Bool_t          HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL;
    Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL;
    Bool_t          HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL;
    Bool_t          HLT_DoubleMu18NoFiltersNoVtx;
    Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight;
    Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose;
    Bool_t          HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose;
    Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight;
    Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose;
    Bool_t          HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose;
    Bool_t          HLT_Mu28NoFiltersNoVtx_CentralCaloJet40;
    Bool_t          HLT_PFHT300_PFMET100;
    Bool_t          HLT_PFHT300_PFMET110;
    Bool_t          HLT_PFHT550_4JetPt50;
    Bool_t          HLT_PFHT650_4JetPt50;
    Bool_t          HLT_PFHT750_4JetPt50;
    Bool_t          HLT_PFHT750_4JetPt70;
    Bool_t          HLT_PFHT750_4JetPt80;
    Bool_t          HLT_PFHT800_4JetPt50;
    Bool_t          HLT_PFHT850_4JetPt50;
    Bool_t          HLT_PFJet15_NoCaloMatched;
    Bool_t          HLT_PFJet25_NoCaloMatched;
    Bool_t          HLT_DiPFJet15_NoCaloMatched;
    Bool_t          HLT_DiPFJet25_NoCaloMatched;
    Bool_t          HLT_DiPFJet15_FBEta3_NoCaloMatched;
    Bool_t          HLT_DiPFJet25_FBEta3_NoCaloMatched;
    Bool_t          HLT_DiPFJetAve15_HFJEC;
    Bool_t          HLT_DiPFJetAve25_HFJEC;
    Bool_t          HLT_DiPFJetAve35_HFJEC;
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
    Bool_t          HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140;
    Bool_t          HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80;
    Bool_t          HLT_DiCentralPFJet170;
    Bool_t          HLT_SingleCentralPFJet170_CFMax0p1;
    Bool_t          HLT_DiCentralPFJet170_CFMax0p1;
    Bool_t          HLT_DiCentralPFJet220_CFMax0p3;
    Bool_t          HLT_DiCentralPFJet330_CFMax0p5;
    Bool_t          HLT_DiCentralPFJet430;
    Bool_t          HLT_PFHT125;
    Bool_t          HLT_PFHT200;
    Bool_t          HLT_PFHT250;
    Bool_t          HLT_PFHT300;
    Bool_t          HLT_PFHT350;
    Bool_t          HLT_PFHT400;
    Bool_t          HLT_PFHT475;
    Bool_t          HLT_PFHT600;
    Bool_t          HLT_PFHT650;
    Bool_t          HLT_PFHT800;
    Bool_t          HLT_PFHT900;
    Bool_t          HLT_PFHT200_PFAlphaT0p51;
    Bool_t          HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57;
    Bool_t          HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63;
    Bool_t          HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55;
    Bool_t          HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58;
    Bool_t          HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53;
    Bool_t          HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54;
    Bool_t          HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52;
    Bool_t          HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53;
    Bool_t          HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51;
    Bool_t          HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52;
    Bool_t          HLT_MET60_IsoTrk35_Loose;
    Bool_t          HLT_MET75_IsoTrk50;
    Bool_t          HLT_MET90_IsoTrk50;
    Bool_t          HLT_PFMET120_BTagCSV_p067;
    Bool_t          HLT_PFMET120_Mu5;
    Bool_t          HLT_PFMET170_NotCleaned;
    Bool_t          HLT_PFMET170_NoiseCleaned;
    Bool_t          HLT_PFMET170_HBHECleaned;
    Bool_t          HLT_PFMET170_JetIdCleaned;
    Bool_t          HLT_PFMET170_BeamHaloCleaned;
    Bool_t          HLT_PFMET170_HBHE_BeamHaloCleaned;
    Bool_t          HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned;
    Bool_t          HLT_PFMET90_PFMHT90_IDTight;
    Bool_t          HLT_PFMET100_PFMHT100_IDTight;
    Bool_t          HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned;
    Bool_t          HLT_PFMET110_PFMHT110_IDTight;
    Bool_t          HLT_PFMET120_PFMHT120_IDTight;
    Bool_t          HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067;
    Bool_t          HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight;
    Bool_t          HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200;
    Bool_t          HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460;
    Bool_t          HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240;
    Bool_t          HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500;
    Bool_t          HLT_QuadPFJet_VBF;
    Bool_t          HLT_L1_TripleJet_VBF;
    Bool_t          HLT_QuadJet45_TripleBTagCSV_p087;
    Bool_t          HLT_QuadJet45_DoubleBTagCSV_p087;
    Bool_t          HLT_DoubleJet90_Double30_TripleBTagCSV_p087;
    Bool_t          HLT_DoubleJet90_Double30_DoubleBTagCSV_p087;
    Bool_t          HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160;
    Bool_t          HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6;
    Bool_t          HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172;
    Bool_t          HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6;
    Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p026;
    Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p014;
    Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350;
    Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350;
    Bool_t          HLT_Photon135_PFMET100;
    Bool_t          HLT_Photon20_CaloIdVL_IsoL;
    Bool_t          HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Photon250_NoHE;
    Bool_t          HLT_Photon300_NoHE;
    Bool_t          HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60;
    Bool_t          HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15;
    Bool_t          HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40;
    Bool_t          HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF;
    Bool_t          HLT_Mu8_TrkIsoVVL;
    Bool_t          HLT_Mu17_TrkIsoVVL;
    Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
    Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
    Bool_t          HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
    Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
    Bool_t          HLT_BTagMu_DiJet20_Mu5;
    Bool_t          HLT_BTagMu_DiJet40_Mu5;
    Bool_t          HLT_BTagMu_DiJet70_Mu5;
    Bool_t          HLT_BTagMu_DiJet110_Mu5;
    Bool_t          HLT_BTagMu_DiJet170_Mu5;
    Bool_t          HLT_BTagMu_Jet300_Mu5;
    Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
    Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded;
    Bool_t          HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
    Bool_t          HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
    Bool_t          HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
    Bool_t          HLT_Mu12_Photon25_CaloIdL;
    Bool_t          HLT_Mu12_Photon25_CaloIdL_L1ISO;
    Bool_t          HLT_Mu12_Photon25_CaloIdL_L1OR;
    Bool_t          HLT_Mu17_Photon22_CaloIdL_L1ISO;
    Bool_t          HLT_Mu17_Photon30_CaloIdL_L1ISO;
    Bool_t          HLT_Mu17_Photon35_CaloIdL_L1ISO;
    Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
    Bool_t          HLT_TripleMu_5_3_3;
    Bool_t          HLT_TripleMu_12_10_5;
    Bool_t          HLT_Mu3er_PFHT140_PFMET125;
    Bool_t          HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067;
    Bool_t          HLT_Mu6_PFHT200_PFMET100;
    Bool_t          HLT_Mu14er_PFMET100;
    Bool_t          HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Ele17_CaloIdL_GsfTrkIdVL;
    Bool_t          HLT_Ele17_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL;
    Bool_t          HLT_PFHT650_WideJetMJJ900DEtaJJ1p5;
    Bool_t          HLT_PFHT650_WideJetMJJ950DEtaJJ1p5;
    Bool_t          HLT_Photon22;
    Bool_t          HLT_Photon30;
    Bool_t          HLT_Photon36;
    Bool_t          HLT_Photon50;
    Bool_t          HLT_Photon75;
    Bool_t          HLT_Photon90;
    Bool_t          HLT_Photon120;
    Bool_t          HLT_Photon175;
    Bool_t          HLT_Photon165_HE10;
    Bool_t          HLT_Photon22_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon30_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon36_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
    Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
    Bool_t          HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
    Bool_t          HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70;
    Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;
    Bool_t          HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55;
    Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;
    Bool_t          HLT_Dimuon0_Jpsi_Muon;
    Bool_t          HLT_Dimuon0_Upsilon_Muon;
    Bool_t          HLT_QuadMuon0_Dimuon0_Jpsi;
    Bool_t          HLT_QuadMuon0_Dimuon0_Upsilon;
    Bool_t          HLT_Rsq0p25_Calo;
    Bool_t          HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo;
    Bool_t          HLT_RsqMR240_Rsq0p09_MR200_Calo;
    Bool_t          HLT_Rsq0p25;
    Bool_t          HLT_Rsq0p30;
    Bool_t          HLT_RsqMR240_Rsq0p09_MR200;
    Bool_t          HLT_RsqMR240_Rsq0p09_MR200_4jet;
    Bool_t          HLT_RsqMR270_Rsq0p09_MR200;
    Bool_t          HLT_RsqMR270_Rsq0p09_MR200_4jet;
    Bool_t          HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200;
    Bool_t          HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
    Bool_t          HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
    Bool_t          HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
    Bool_t          HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
    Bool_t          HLT_HT200_DisplacedDijet40_DisplacedTrack;
    Bool_t          HLT_HT250_DisplacedDijet40_DisplacedTrack;
    Bool_t          HLT_HT350_DisplacedDijet40_DisplacedTrack;
    Bool_t          HLT_HT350_DisplacedDijet80_DisplacedTrack;
    Bool_t          HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack;
    Bool_t          HLT_HT350_DisplacedDijet40_Inclusive;
    Bool_t          HLT_HT400_DisplacedDijet40_Inclusive;
    Bool_t          HLT_HT500_DisplacedDijet40_Inclusive;
    Bool_t          HLT_HT550_DisplacedDijet40_Inclusive;
    Bool_t          HLT_HT550_DisplacedDijet80_Inclusive;
    Bool_t          HLT_HT650_DisplacedDijet80_Inclusive;
    Bool_t          HLT_HT750_DisplacedDijet80_Inclusive;
    Bool_t          HLT_VBF_DisplacedJet40_DisplacedTrack;
    Bool_t          HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5;
    Bool_t          HLT_VBF_DisplacedJet40_TightID_DisplacedTrack;
    Bool_t          HLT_VBF_DisplacedJet40_Hadronic;
    Bool_t          HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack;
    Bool_t          HLT_VBF_DisplacedJet40_TightID_Hadronic;
    Bool_t          HLT_VBF_DisplacedJet40_VTightID_Hadronic;
    Bool_t          HLT_VBF_DisplacedJet40_VVTightID_Hadronic;
    Bool_t          HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack;
    Bool_t          HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack;
    Bool_t          HLT_PFMETNoMu90_PFMHTNoMu90_IDTight;
    Bool_t          HLT_PFMETNoMu100_PFMHTNoMu100_IDTight;
    Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
    Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
    Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight;
    Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight;
    Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
    Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
    Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf_HT200;
    Bool_t          HLT_Photon90_CaloIdL_PFHT500;
    Bool_t          HLT_DoubleMu8_Mass8_PFHT250;
    Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250;
    Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250;
    Bool_t          HLT_DoubleMu8_Mass8_PFHT300;
    Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t          HLT_Mu10_CentralPFJet30_BTagCSV_p13;
    Bool_t          HLT_DoubleMu3_PFMET50;
    Bool_t          HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13;
    Bool_t          HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400;
    Bool_t          HLT_Ele15_IsoVVVL_PFHT350_PFMET50;
    Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
    Bool_t          HLT_Ele15_IsoVVVL_PFHT350;
    Bool_t          HLT_Ele15_IsoVVVL_PFHT400_PFMET50;
    Bool_t          HLT_Ele15_IsoVVVL_PFHT400;
    Bool_t          HLT_Ele50_IsoVVVL_PFHT400;
    Bool_t          HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
    Bool_t          HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
    Bool_t          HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400;
    Bool_t          HLT_Mu15_IsoVVVL_PFHT350_PFMET50;
    Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
    Bool_t          HLT_Mu15_IsoVVVL_PFHT350;
    Bool_t          HLT_Mu15_IsoVVVL_PFHT400_PFMET50;
    Bool_t          HLT_Mu15_IsoVVVL_PFHT400;
    Bool_t          HLT_Mu50_IsoVVVL_PFHT400;
    Bool_t          HLT_Dimuon16_Jpsi;
    Bool_t          HLT_Dimuon10_Jpsi_Barrel;
    Bool_t          HLT_Dimuon8_PsiPrime_Barrel;
    Bool_t          HLT_Dimuon8_Upsilon_Barrel;
    Bool_t          HLT_Dimuon0_Phi_Barrel;
    Bool_t          HLT_Mu16_TkMu0_dEta18_Onia;
    Bool_t          HLT_Mu16_TkMu0_dEta18_Phi;
    Bool_t          HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx;
    Bool_t          HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
    Bool_t          HLT_Mu8;
    Bool_t          HLT_Mu17;
    Bool_t          HLT_Mu3_PFJet40;
    Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
    Bool_t          HLT_Ele12_CaloIdM_TrackIdM_PFJet30;
    Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
    Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
    Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140;
    Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
    Bool_t          HLT_PFHT400_SixJet30_DoubleBTagCSV_p056;
    Bool_t          HLT_PFHT450_SixJet40_BTagCSV_p056;
    Bool_t          HLT_PFHT400_SixJet30;
    Bool_t          HLT_PFHT450_SixJet40;
    Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
    Bool_t          HLT_Mu55;
    Bool_t          HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15;
    Bool_t          HLT_Photon90_CaloIdL_PFHT600;
    Bool_t          HLT_PixelTracks_Multiplicity60ForEndOfFill;
    Bool_t          HLT_PixelTracks_Multiplicity85ForEndOfFill;
    Bool_t          HLT_PixelTracks_Multiplicity110ForEndOfFill;
    Bool_t          HLT_PixelTracks_Multiplicity135ForEndOfFill;
    Bool_t          HLT_PixelTracks_Multiplicity160ForEndOfFill;
    Bool_t          HLT_FullTracks_Multiplicity80;
    Bool_t          HLT_FullTracks_Multiplicity100;
    Bool_t          HLT_FullTracks_Multiplicity130;
    Bool_t          HLT_FullTracks_Multiplicity150;
    Bool_t          HLT_ECALHT800;
    Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
    Bool_t          HLT_Photon125;
    Bool_t          HLT_MET100;
    Bool_t          HLT_MET150;
    Bool_t          HLT_MET200;
    Bool_t          HLT_Ele27_HighEta_Ele20_Mass55;
    Bool_t          HLT_L1FatEvents;
    Bool_t          HLT_Physics;
    Bool_t          HLT_L1FatEvents_part0;
    Bool_t          HLT_L1FatEvents_part1;
    Bool_t          HLT_L1FatEvents_part2;
    Bool_t          HLT_L1FatEvents_part3;
    Bool_t          HLT_Random;
    Bool_t          HLT_ZeroBias;
    Bool_t          HLT_AK4CaloJet30;
    Bool_t          HLT_AK4CaloJet40;
    Bool_t          HLT_AK4CaloJet50;
    Bool_t          HLT_AK4CaloJet80;
    Bool_t          HLT_AK4CaloJet100;
    Bool_t          HLT_AK4PFJet30;
    Bool_t          HLT_AK4PFJet50;
    Bool_t          HLT_AK4PFJet80;
    Bool_t          HLT_AK4PFJet100;
    Bool_t          HLT_HISinglePhoton10;
    Bool_t          HLT_HISinglePhoton15;
    Bool_t          HLT_HISinglePhoton20;
    Bool_t          HLT_HISinglePhoton40;
    Bool_t          HLT_HISinglePhoton60;
    Bool_t          HLT_EcalCalibration;
    Bool_t          HLT_HcalCalibration;
    Bool_t          HLT_GlobalRunHPDNoise;
    Bool_t          HLT_L1BptxMinus;
    Bool_t          HLT_L1BptxPlus;
    Bool_t          HLT_L1NotBptxOR;
    Bool_t          HLT_L1BeamGasMinus;
    Bool_t          HLT_L1BeamGasPlus;
    Bool_t          HLT_L1BptxXOR;
    Bool_t          HLT_L1MinimumBiasHF_OR;
    Bool_t          HLT_L1MinimumBiasHF_AND;
    Bool_t          HLT_HcalNZS;
    Bool_t          HLT_HcalPhiSym;
    Bool_t          HLT_HcalIsolatedbunch;
    Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
    Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap_copy;
    Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS;
    Bool_t          HLT_ZeroBias_IsolatedBunches;
    Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
    Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
    Bool_t          HLT_Photon500;
    Bool_t          HLT_Photon600;
    Bool_t          HLT_Mu300;
    Bool_t          HLT_Mu350;
    Bool_t          HLT_MET250;
    Bool_t          HLT_MET300;
    Bool_t          HLT_MET600;
    Bool_t          HLT_MET700;
    Bool_t          HLT_PFMET300;
    Bool_t          HLT_PFMET400;
    Bool_t          HLT_PFMET500;
    Bool_t          HLT_PFMET600;
    Bool_t          HLT_Ele250_CaloIdVT_GsfTrkIdT;
    Bool_t          HLT_Ele300_CaloIdVT_GsfTrkIdT;
    Bool_t          HLT_HT2000;
    Bool_t          HLT_HT2500;
    Bool_t          HLT_IsoTrackHE;
    Bool_t          HLT_IsoTrackHB;
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
    Bool_t          Flag_goodVertices;
    Bool_t          Flag_eeBadScFilter;
    Bool_t          Flag_ecalLaserCorrFilter;
    Bool_t          Flag_trkPOGFilters;
    Bool_t          Flag_chargedHadronTrackResolutionFilter;
    Bool_t          Flag_muonBadTrackFilter;
    Bool_t          Flag_trkPOG_manystripclus53X;
    Bool_t          Flag_trkPOG_toomanystripclus53X;
    Bool_t          Flag_trkPOG_logErrorTooManyClusters;
    Bool_t          Flag_METFilters;


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

void EventLoop(int sampleIdx){
//Loops over all events and fills histograms
    int eventLoopMax = EvMax;

    for(int iev=0; iev<eventLoopMax;++iev){
        EventTree->GetEvent(iev);
        if(iev % 10000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;


        h_CaloMET_phi.at(sampleIdx)->Fill(CaloMET_phi);
        h_CaloMET_pt.at(sampleIdx)->Fill(CaloMET_pt);
        h_CaloMET_sumEt.at(sampleIdx)->Fill(CaloMET_sumEt);

        for(int iElectron=0; iElectron<nElectron; iElectron++){

            h_Electron_eta.at(sampleIdx)->Fill(Electron_eta[iElectron]);
            h_Electron_hoe.at(sampleIdx)->Fill(Electron_hoe[iElectron]); 
            h_Electron_mass.at(sampleIdx)->Fill(Electron_mass[iElectron]);
            h_Electron_phi.at(sampleIdx)->Fill(Electron_phi[iElectron]);
            h_Electron_pt.at(sampleIdx)->Fill(Electron_pt[iElectron]);
            h_Electron_r9.at(sampleIdx)->Fill(Electron_r9[iElectron]);
            h_Electron_sieie.at(sampleIdx)->Fill(Electron_sieie[iElectron]);
        }
        for(int iMuon=0; iMuon<nMuon; iMuon++){

            h_Muon_eta .at(sampleIdx)->Fill(Muon_eta[iMuon]);
            h_Muon_mass.at(sampleIdx)->Fill(Muon_mass[iMuon]);
            h_Muon_phi.at(sampleIdx)->Fill(Muon_phi[iMuon]);
            h_Muon_pt.at(sampleIdx)->Fill(Muon_pt[iMuon]);
        }
        for(int iJet=0; iJet<nJet; iJet++){
            h_Jet_eta.at(sampleIdx)->Fill(Jet_eta[iJet]);
            h_Jet_mass.at(sampleIdx)->Fill(Jet_mass[iJet]);
            h_Jet_phi.at(sampleIdx)->Fill(Jet_phi[iJet]);
            h_Jet_pt.at(sampleIdx)->Fill(Jet_pt[iJet]);
        }


    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;
}

void InitHistograms(){

    for(int i=0; i<nSamples; i++){

        h_CaloMET_phi[i] = new TH1F(("h_CaloMET_phi_"+to_string(i)).c_str(), ("h_CaloMET_phi_"+to_string(i)).c_str(), 200, -3.5, 3.5);
        h_CaloMET_pt[i] = new TH1F(("h_CaloMET_pt_"+to_string(i)).c_str(), ("h_CaloMET_pt_"+to_string(i)).c_str(), 200, 0, 200);
        h_CaloMET_sumEt[i] = new TH1F(("h_CaloMET_sumEt_"+to_string(i)).c_str(), ("h_CaloMET_sumEt_"+to_string(i)).c_str(), 500, 0, 1000);

        h_Electron_eta[i] = new TH1F(("h_Electron_eta_"+to_string(i)).c_str(), ("h_Electron_eta_"+to_string(i)).c_str(), 200, -5.0, 5.0);
        h_Electron_hoe[i] = new TH1F(("h_Electron_hoe_"+to_string(i)).c_str(), ("h_Electron_hoe_"+to_string(i)).c_str(), 200, 0, 5.0);
        h_Electron_mass[i] = new TH1F(("h_Electron_mass_"+to_string(i)).c_str(), ("h_Electron_mass_"+to_string(i)).c_str(), 200, 0, 0.15);
        h_Electron_phi[i] = new TH1F(("h_Electron_phi_"+to_string(i)).c_str(), ("h_Electron_phi_"+to_string(i)).c_str(), 200, -3.5, 3.5);
        h_Electron_pt[i] = new TH1F(("h_Electron_pt_"+to_string(i)).c_str(), ("h_Electron_pt_"+to_string(i)).c_str(), 200, 0, 200);
        h_Electron_r9[i] = new TH1F(("h_Electron_r9_"+to_string(i)).c_str(), ("h_Electron_r9_"+to_string(i)).c_str(), 200, 0, 4.0);
        h_Electron_sieie[i] = new TH1F(("h_Electron_sieie_"+to_string(i)).c_str(), ("h_Electron_sieie_"+to_string(i)).c_str(), 200, 0, 0.2);

        h_Muon_eta[i] = new TH1F(("h_Muon_eta_"+to_string(i)).c_str(), ("h_Muon_eta_"+to_string(i)).c_str(), 200, -5.0, 5.0);
        h_Muon_mass[i] = new TH1F(("h_Muon_mass_"+to_string(i)).c_str(), ("h_Muon_mass_"+to_string(i)).c_str(), 200, 0, 0.15);
        h_Muon_phi[i] = new TH1F(("h_Muon_phi_"+to_string(i)).c_str(), ("h_Muon_phi_"+to_string(i)).c_str(), 200, -3.5, 3.5);
        h_Muon_pt[i] = new TH1F(("h_Muon_pt_"+to_string(i)).c_str(), ("h_Muon_pt_"+to_string(i)).c_str(), 200, 0, 200);

        h_Jet_eta[i] = new TH1F(("h_Jet_eta_"+to_string(i)).c_str(), ("h_Jet_eta_"+to_string(i)).c_str(), 200, -5.0, 5.0);
        h_Jet_mass[i] = new TH1F(("h_Jet_mass_"+to_string(i)).c_str(), ("h_Jet_mass_"+to_string(i)).c_str(), 200, 0, 200);
        h_Jet_phi[i] = new TH1F(("h_Jet_phi_"+to_string(i)).c_str(), ("h_Jet_phi_"+to_string(i)).c_str(), 200, -3.5, 3.5);
        h_Jet_pt[i] = new TH1F(("h_Jet_pt_"+to_string(i)).c_str(), ("h_Jet_pt_"+to_string(i)).c_str(), 200, 0, 200);
    }
    
}

void DrawPlot(TH1F* plot, string name){

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    string outdir = "../Output/ttHTobb/";
    setOutput(outdir);

    plot->Scale(1.0 / plot->Integral());
    plot->SetFillColorAlpha(kGreen + 1,0.35);

    plot->Draw("hist");
    c1->SaveAs((outdir+name+".png").c_str());
    c1->SaveAs((outdir+name+".pdf").c_str());
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
 
}

void InitTree(TString infileName, int sampleIdx){

    cout<<cur_time()<<"\tProcessing Tree...\n";
    TFile* infile = new TFile(infileName);
    EventTree=(TTree*)gDirectory->Get("Events");


    EventTree->SetBranchAddress("run", &run);
    EventTree->SetBranchAddress("luminosityBlock", &luminosityBlock);
    EventTree->SetBranchAddress("event", &event);
    EventTree->SetBranchAddress("CaloMET_phi", &CaloMET_phi);
    EventTree->SetBranchAddress("CaloMET_pt", &CaloMET_pt);
    EventTree->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt);
    EventTree->SetBranchAddress("nElectron", &nElectron);
    EventTree->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC);
    EventTree->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt);
    EventTree->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt);
    EventTree->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt);
    EventTree->SetBranchAddress("Electron_dxy", &Electron_dxy);
    EventTree->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr);
    EventTree->SetBranchAddress("Electron_dz", &Electron_dz);
    EventTree->SetBranchAddress("Electron_dzErr", &Electron_dzErr);
    EventTree->SetBranchAddress("Electron_eCorr", &Electron_eCorr);
    EventTree->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv);
    EventTree->SetBranchAddress("Electron_energyErr", &Electron_energyErr);
    EventTree->SetBranchAddress("Electron_eta", &Electron_eta);
    EventTree->SetBranchAddress("Electron_hoe", &Electron_hoe);
    EventTree->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    EventTree->SetBranchAddress("Electron_mass", &Electron_mass);
    EventTree->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all);
    EventTree->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg);
    EventTree->SetBranchAddress("Electron_mvaSpring16GP", &Electron_mvaSpring16GP);
    EventTree->SetBranchAddress("Electron_mvaSpring16HZZ", &Electron_mvaSpring16HZZ);
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
    EventTree->SetBranchAddress("Electron_cutBased_HLTPreSel", &Electron_cutBased_HLTPreSel);
    EventTree->SetBranchAddress("Electron_jetIdx", &Electron_jetIdx);
    EventTree->SetBranchAddress("Electron_pdgId", &Electron_pdgId);
    EventTree->SetBranchAddress("Electron_photonIdx", &Electron_photonIdx);
    EventTree->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge);
    EventTree->SetBranchAddress("Electron_vidNestedWPBitmap", &Electron_vidNestedWPBitmap);
    EventTree->SetBranchAddress("Electron_convVeto", &Electron_convVeto);
    EventTree->SetBranchAddress("Electron_cutBased_HEEP", &Electron_cutBased_HEEP);
    EventTree->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand);
    EventTree->SetBranchAddress("Electron_lostHits", &Electron_lostHits);
    EventTree->SetBranchAddress("Electron_mvaSpring16GP_WP80", &Electron_mvaSpring16GP_WP80);
    EventTree->SetBranchAddress("Electron_mvaSpring16GP_WP90", &Electron_mvaSpring16GP_WP90);
    EventTree->SetBranchAddress("Electron_mvaSpring16HZZ_WPL", &Electron_mvaSpring16HZZ_WPL);
    EventTree->SetBranchAddress("nFatJet", &nFatJet);
    EventTree->SetBranchAddress("FatJet_area", &FatJet_area);
    EventTree->SetBranchAddress("FatJet_btagCMVA", &FatJet_btagCMVA);
    EventTree->SetBranchAddress("FatJet_btagCSVV2", &FatJet_btagCSVV2);
    EventTree->SetBranchAddress("FatJet_btagDeepB", &FatJet_btagDeepB);
    EventTree->SetBranchAddress("FatJet_btagHbb", &FatJet_btagHbb);
    EventTree->SetBranchAddress("FatJet_eta", &FatJet_eta);
    EventTree->SetBranchAddress("FatJet_mass", &FatJet_mass);
    EventTree->SetBranchAddress("FatJet_msoftdrop", &FatJet_msoftdrop);
    EventTree->SetBranchAddress("FatJet_msoftdrop_chs", &FatJet_msoftdrop_chs);
    EventTree->SetBranchAddress("FatJet_n2b1", &FatJet_n2b1);
    EventTree->SetBranchAddress("FatJet_n3b1", &FatJet_n3b1);
    EventTree->SetBranchAddress("FatJet_phi", &FatJet_phi);
    EventTree->SetBranchAddress("FatJet_pt", &FatJet_pt);
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
    EventTree->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight);
    EventTree->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
    EventTree->SetBranchAddress("nJet", &nJet);
    EventTree->SetBranchAddress("Jet_area", &Jet_area);
    EventTree->SetBranchAddress("Jet_btagCMVA", &Jet_btagCMVA);
    EventTree->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2);
    EventTree->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);
    EventTree->SetBranchAddress("Jet_btagDeepC", &Jet_btagDeepC);
    EventTree->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF);
    EventTree->SetBranchAddress("Jet_chHEF", &Jet_chHEF);
    EventTree->SetBranchAddress("Jet_eta", &Jet_eta);
    EventTree->SetBranchAddress("Jet_mass", &Jet_mass);
    EventTree->SetBranchAddress("Jet_neEmEF", &Jet_neEmEF);
    EventTree->SetBranchAddress("Jet_neHEF", &Jet_neHEF);
    EventTree->SetBranchAddress("Jet_phi", &Jet_phi);
    EventTree->SetBranchAddress("Jet_pt", &Jet_pt);
    EventTree->SetBranchAddress("Jet_qgl", &Jet_qgl);
    EventTree->SetBranchAddress("Jet_rawFactor", &Jet_rawFactor);
    EventTree->SetBranchAddress("Jet_bReg", &Jet_bReg);
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
    EventTree->SetBranchAddress("Muon_mvaTTH", &Muon_mvaTTH);
    EventTree->SetBranchAddress("Muon_charge", &Muon_charge);
    EventTree->SetBranchAddress("Muon_jetIdx", &Muon_jetIdx);
    EventTree->SetBranchAddress("Muon_nStations", &Muon_nStations);
    EventTree->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
    EventTree->SetBranchAddress("Muon_pdgId", &Muon_pdgId);
    EventTree->SetBranchAddress("Muon_tightCharge", &Muon_tightCharge);
    EventTree->SetBranchAddress("Muon_highPtId", &Muon_highPtId);
    EventTree->SetBranchAddress("Muon_isPFcand", &Muon_isPFcand);
    EventTree->SetBranchAddress("Muon_mediumId", &Muon_mediumId);
    EventTree->SetBranchAddress("Muon_softId", &Muon_softId);
    EventTree->SetBranchAddress("Muon_tightId", &Muon_tightId);
    EventTree->SetBranchAddress("nPhoton", &nPhoton);
    EventTree->SetBranchAddress("Photon_eCorr", &Photon_eCorr);
    EventTree->SetBranchAddress("Photon_energyErr", &Photon_energyErr);
    EventTree->SetBranchAddress("Photon_eta", &Photon_eta);
    EventTree->SetBranchAddress("Photon_hoe", &Photon_hoe);
    EventTree->SetBranchAddress("Photon_mass", &Photon_mass);
    EventTree->SetBranchAddress("Photon_mvaID", &Photon_mvaID);
    EventTree->SetBranchAddress("Photon_pfRelIso03_all", &Photon_pfRelIso03_all);
    EventTree->SetBranchAddress("Photon_pfRelIso03_chg", &Photon_pfRelIso03_chg);
    EventTree->SetBranchAddress("Photon_phi", &Photon_phi);
    EventTree->SetBranchAddress("Photon_pt", &Photon_pt);
    EventTree->SetBranchAddress("Photon_r9", &Photon_r9);
    EventTree->SetBranchAddress("Photon_sieie", &Photon_sieie);
    EventTree->SetBranchAddress("Photon_charge", &Photon_charge);
    EventTree->SetBranchAddress("Photon_cutBased", &Photon_cutBased);
    EventTree->SetBranchAddress("Photon_electronIdx", &Photon_electronIdx);
    EventTree->SetBranchAddress("Photon_jetIdx", &Photon_jetIdx);
    EventTree->SetBranchAddress("Photon_pdgId", &Photon_pdgId);
    EventTree->SetBranchAddress("Photon_vidNestedWPBitmap", &Photon_vidNestedWPBitmap);
    EventTree->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto);
    EventTree->SetBranchAddress("Photon_mvaID_WP80", &Photon_mvaID_WP80);
    EventTree->SetBranchAddress("Photon_mvaID_WP90", &Photon_mvaID_WP90);
    EventTree->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed);
    EventTree->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);
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
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo);
    EventTree->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral);
    EventTree->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton);
    EventTree->SetBranchAddress("GenDressedLepton_eta", &GenDressedLepton_eta);
    EventTree->SetBranchAddress("GenDressedLepton_mass", &GenDressedLepton_mass);
    EventTree->SetBranchAddress("GenDressedLepton_phi", &GenDressedLepton_phi);
    EventTree->SetBranchAddress("GenDressedLepton_pt", &GenDressedLepton_pt);
    EventTree->SetBranchAddress("GenDressedLepton_pdgId", &GenDressedLepton_pdgId);
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
    EventTree->SetBranchAddress("SubJet_tau1", &SubJet_tau1);
    EventTree->SetBranchAddress("SubJet_tau2", &SubJet_tau2);
    EventTree->SetBranchAddress("SubJet_tau3", &SubJet_tau3);
    EventTree->SetBranchAddress("SubJet_tau4", &SubJet_tau4);
    EventTree->SetBranchAddress("nTau", &nTau);
    EventTree->SetBranchAddress("Tau_chargedIso", &Tau_chargedIso);
    EventTree->SetBranchAddress("Tau_dxy", &Tau_dxy);
    EventTree->SetBranchAddress("Tau_dz", &Tau_dz);
    EventTree->SetBranchAddress("Tau_eta", &Tau_eta);
    EventTree->SetBranchAddress("Tau_footprintCorr", &Tau_footprintCorr);
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
    EventTree->SetBranchAddress("Tau_rawIso", &Tau_rawIso);
    EventTree->SetBranchAddress("Tau_rawMVAnewDM", &Tau_rawMVAnewDM);
    EventTree->SetBranchAddress("Tau_rawMVAoldDM", &Tau_rawMVAoldDM);
    EventTree->SetBranchAddress("Tau_rawMVAoldDMdR03", &Tau_rawMVAoldDMdR03);
    EventTree->SetBranchAddress("Tau_charge", &Tau_charge);
    EventTree->SetBranchAddress("Tau_decayMode", &Tau_decayMode);
    EventTree->SetBranchAddress("Tau_jetIdx", &Tau_jetIdx);
    EventTree->SetBranchAddress("Tau_rawAntiEleCat", &Tau_rawAntiEleCat);
    EventTree->SetBranchAddress("Tau_idAntiEle", &Tau_idAntiEle);
    EventTree->SetBranchAddress("Tau_idAntiMu", &Tau_idAntiMu);
    EventTree->SetBranchAddress("Tau_idDecayMode", &Tau_idDecayMode);
    EventTree->SetBranchAddress("Tau_idDecayModeNewDMs", &Tau_idDecayModeNewDMs);
    EventTree->SetBranchAddress("Tau_idMVAnewDM", &Tau_idMVAnewDM);
    EventTree->SetBranchAddress("Tau_idMVAoldDM", &Tau_idMVAoldDM);
    EventTree->SetBranchAddress("Tau_idMVAoldDMdR03", &Tau_idMVAoldDMdR03);
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
    EventTree->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30", &HLT_AK8DiPFJet300_200_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50);
    EventTree->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30", &HLT_AK8DiPFJet280_200_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30", &HLT_AK8DiPFJet250_200_TrimMass30);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20);
    EventTree->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20);
    EventTree->SetBranchAddress("HLT_CaloJet260", &HLT_CaloJet260);
    EventTree->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID);
    EventTree->SetBranchAddress("HLT_Dimuon13_PsiPrime", &HLT_Dimuon13_PsiPrime);
    EventTree->SetBranchAddress("HLT_Dimuon13_Upsilon", &HLT_Dimuon13_Upsilon);
    EventTree->SetBranchAddress("HLT_Dimuon20_Jpsi", &HLT_Dimuon20_Jpsi);
    EventTree->SetBranchAddress("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf", &HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_DoubleEle25_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle25_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL);
    EventTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW);
    EventTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW);
    EventTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1);
    EventTree->SetBranchAddress("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtx", &HLT_DoubleMu33NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_DoubleMu38NoFiltersNoVtx", &HLT_DoubleMu38NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_DoubleMu23NoFiltersNoVtxDisplaced", &HLT_DoubleMu23NoFiltersNoVtxDisplaced);
    EventTree->SetBranchAddress("HLT_DoubleMu28NoFiltersNoVtxDisplaced", &HLT_DoubleMu28NoFiltersNoVtxDisplaced);
    EventTree->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0);
    EventTree->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs);
    EventTree->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
    EventTree->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu);
    EventTree->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced);
    EventTree->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon);
    EventTree->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon);
    EventTree->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing);
    EventTree->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoVertexing);
    EventTree->SetBranchAddress("HLT_Dimuon6_Jpsi_NoVertexing", &HLT_Dimuon6_Jpsi_NoVertexing);
    EventTree->SetBranchAddress("HLT_Photon150", &HLT_Photon150);
    EventTree->SetBranchAddress("HLT_Photon90_CaloIdL_HT300", &HLT_Photon90_CaloIdL_HT300);
    EventTree->SetBranchAddress("HLT_HT250_CaloMET70", &HLT_HT250_CaloMET70);
    EventTree->SetBranchAddress("HLT_DoublePhoton60", &HLT_DoublePhoton60);
    EventTree->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85);
    EventTree->SetBranchAddress("HLT_Ele17_Ele8_Gsf", &HLT_Ele17_Ele8_Gsf);
    EventTree->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28", &HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28);
    EventTree->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29);
    EventTree->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf", &HLT_Ele22_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_WHbbBoost", &HLT_Ele23_WPLoose_Gsf_WHbbBoost);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf", &HLT_Ele24_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30);
    EventTree->SetBranchAddress("HLT_Ele25_WPTight_Gsf", &HLT_Ele25_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele25_eta2p1_WPLoose_Gsf", &HLT_Ele25_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf", &HLT_Ele25_eta2p1_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele27_WPLoose_Gsf", &HLT_Ele27_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_WHbbBoost", &HLT_Ele27_WPLoose_Gsf_WHbbBoost);
    EventTree->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded", &HLT_Ele27_WPTight_Gsf_L1JetTauSeeded);
    EventTree->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf", &HLT_Ele27_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf", &HLT_Ele27_eta2p1_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele30_eta2p1_WPLoose_Gsf", &HLT_Ele30_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf", &HLT_Ele30_eta2p1_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf", &HLT_Ele32_eta2p1_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf);
    EventTree->SetBranchAddress("HLT_Ele35_WPLoose_Gsf", &HLT_Ele35_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50);
    EventTree->SetBranchAddress("HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_Ele45_WPLoose_Gsf", &HLT_Ele45_WPLoose_Gsf);
    EventTree->SetBranchAddress("HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", &HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded);
    EventTree->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50);
    EventTree->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT", &HLT_Ele105_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele30WP60_SC4_Mass55", &HLT_Ele30WP60_SC4_Mass55);
    EventTree->SetBranchAddress("HLT_Ele30WP60_Ele8_Mass55", &HLT_Ele30WP60_Ele8_Mass55);
    EventTree->SetBranchAddress("HLT_HT200", &HLT_HT200);
    EventTree->SetBranchAddress("HLT_HT275", &HLT_HT275);
    EventTree->SetBranchAddress("HLT_HT325", &HLT_HT325);
    EventTree->SetBranchAddress("HLT_HT425", &HLT_HT425);
    EventTree->SetBranchAddress("HLT_HT575", &HLT_HT575);
    EventTree->SetBranchAddress("HLT_HT410to430", &HLT_HT410to430);
    EventTree->SetBranchAddress("HLT_HT430to450", &HLT_HT430to450);
    EventTree->SetBranchAddress("HLT_HT450to470", &HLT_HT450to470);
    EventTree->SetBranchAddress("HLT_HT470to500", &HLT_HT470to500);
    EventTree->SetBranchAddress("HLT_HT500to550", &HLT_HT500to550);
    EventTree->SetBranchAddress("HLT_HT550to650", &HLT_HT550to650);
    EventTree->SetBranchAddress("HLT_HT650", &HLT_HT650);
    EventTree->SetBranchAddress("HLT_Mu16_eta2p1_MET30", &HLT_Mu16_eta2p1_MET30);
    EventTree->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30", &HLT_IsoMu16_eta2p1_MET30);
    EventTree->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1", &HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1);
    EventTree->SetBranchAddress("HLT_IsoMu17_eta2p1", &HLT_IsoMu17_eta2p1);
    EventTree->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20);
    EventTree->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1", &HLT_DoubleIsoMu17_eta2p1);
    EventTree->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1_noDzCut", &HLT_DoubleIsoMu17_eta2p1_noDzCut);
    EventTree->SetBranchAddress("HLT_IsoMu18", &HLT_IsoMu18);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20);
    EventTree->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1);
    EventTree->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
    EventTree->SetBranchAddress("HLT_IsoMu22", &HLT_IsoMu22);
    EventTree->SetBranchAddress("HLT_IsoMu22_eta2p1", &HLT_IsoMu22_eta2p1);
    EventTree->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    EventTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
    EventTree->SetBranchAddress("HLT_IsoTkMu18", &HLT_IsoTkMu18);
    EventTree->SetBranchAddress("HLT_IsoTkMu20", &HLT_IsoTkMu20);
    EventTree->SetBranchAddress("HLT_IsoTkMu22", &HLT_IsoTkMu22);
    EventTree->SetBranchAddress("HLT_IsoTkMu22_eta2p1", &HLT_IsoTkMu22_eta2p1);
    EventTree->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24);
    EventTree->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27);
    EventTree->SetBranchAddress("HLT_JetE30_NoBPTX3BX", &HLT_JetE30_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_JetE30_NoBPTX", &HLT_JetE30_NoBPTX);
    EventTree->SetBranchAddress("HLT_JetE50_NoBPTX3BX", &HLT_JetE50_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_JetE70_NoBPTX3BX", &HLT_JetE70_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18);
    EventTree->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10);
    EventTree->SetBranchAddress("HLT_L1SingleMuOpen", &HLT_L1SingleMuOpen);
    EventTree->SetBranchAddress("HLT_L1SingleMuOpen_DT", &HLT_L1SingleMuOpen_DT);
    EventTree->SetBranchAddress("HLT_L2DoubleMu23_NoVertex", &HLT_L2DoubleMu23_NoVertex);
    EventTree->SetBranchAddress("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10);
    EventTree->SetBranchAddress("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10);
    EventTree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX);
    EventTree->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
    EventTree->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1", &HLT_LooseIsoPFTau50_Trk30_eta2p1);
    EventTree->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80);
    EventTree->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90);
    EventTree->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110);
    EventTree->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120);
    EventTree->SetBranchAddress("HLT_PFTau120_eta2p1", &HLT_PFTau120_eta2p1);
    EventTree->SetBranchAddress("HLT_PFTau140_eta2p1", &HLT_PFTau140_eta2p1);
    EventTree->SetBranchAddress("HLT_VLooseIsoPFTau120_Trk50_eta2p1", &HLT_VLooseIsoPFTau120_Trk50_eta2p1);
    EventTree->SetBranchAddress("HLT_VLooseIsoPFTau140_Trk50_eta2p1", &HLT_VLooseIsoPFTau140_Trk50_eta2p1);
    EventTree->SetBranchAddress("HLT_Mu17_Mu8", &HLT_Mu17_Mu8);
    EventTree->SetBranchAddress("HLT_Mu17_Mu8_DZ", &HLT_Mu17_Mu8_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_Mu8_SameSign", &HLT_Mu17_Mu8_SameSign);
    EventTree->SetBranchAddress("HLT_Mu17_Mu8_SameSign_DZ", &HLT_Mu17_Mu8_SameSign_DZ);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign);
    EventTree->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_TkMu8_DZ", &HLT_Mu17_TkMu8_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu25_TkMu0_dEta18_Onia", &HLT_Mu25_TkMu0_dEta18_Onia);
    EventTree->SetBranchAddress("HLT_Mu27_TkMu8", &HLT_Mu27_TkMu8);
    EventTree->SetBranchAddress("HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11);
    EventTree->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50", &HLT_Mu30_eta2p1_PFJet150_PFJet50);
    EventTree->SetBranchAddress("HLT_Mu40_TkMu11", &HLT_Mu40_TkMu11);
    EventTree->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50", &HLT_Mu40_eta2p1_PFJet200_PFJet50);
    EventTree->SetBranchAddress("HLT_Mu20", &HLT_Mu20);
    EventTree->SetBranchAddress("HLT_TkMu17", &HLT_TkMu17);
    EventTree->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    EventTree->SetBranchAddress("HLT_TkMu20", &HLT_TkMu20);
    EventTree->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1);
    EventTree->SetBranchAddress("HLT_TkMu24_eta2p1", &HLT_TkMu24_eta2p1);
    EventTree->SetBranchAddress("HLT_Mu27", &HLT_Mu27);
    EventTree->SetBranchAddress("HLT_TkMu27", &HLT_TkMu27);
    EventTree->SetBranchAddress("HLT_Mu45_eta2p1", &HLT_Mu45_eta2p1);
    EventTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
    EventTree->SetBranchAddress("HLT_TkMu50", &HLT_TkMu50);
    EventTree->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL", &HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL", &HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL", &HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", &HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL);
    EventTree->SetBranchAddress("HLT_DoubleMu18NoFiltersNoVtx", &HLT_DoubleMu18NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight);
    EventTree->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose);
    EventTree->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose", &HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose);
    EventTree->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight);
    EventTree->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose);
    EventTree->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose);
    EventTree->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_CentralCaloJet40", &HLT_Mu28NoFiltersNoVtx_CentralCaloJet40);
    EventTree->SetBranchAddress("HLT_PFHT300_PFMET100", &HLT_PFHT300_PFMET100);
    EventTree->SetBranchAddress("HLT_PFHT300_PFMET110", &HLT_PFHT300_PFMET110);
    EventTree->SetBranchAddress("HLT_PFHT550_4JetPt50", &HLT_PFHT550_4JetPt50);
    EventTree->SetBranchAddress("HLT_PFHT650_4JetPt50", &HLT_PFHT650_4JetPt50);
    EventTree->SetBranchAddress("HLT_PFHT750_4JetPt50", &HLT_PFHT750_4JetPt50);
    EventTree->SetBranchAddress("HLT_PFHT750_4JetPt70", &HLT_PFHT750_4JetPt70);
    EventTree->SetBranchAddress("HLT_PFHT750_4JetPt80", &HLT_PFHT750_4JetPt80);
    EventTree->SetBranchAddress("HLT_PFHT800_4JetPt50", &HLT_PFHT800_4JetPt50);
    EventTree->SetBranchAddress("HLT_PFHT850_4JetPt50", &HLT_PFHT850_4JetPt50);
    EventTree->SetBranchAddress("HLT_PFJet15_NoCaloMatched", &HLT_PFJet15_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_PFJet25_NoCaloMatched", &HLT_PFJet25_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched);
    EventTree->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC);
    EventTree->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC);
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
    EventTree->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140);
    EventTree->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80);
    EventTree->SetBranchAddress("HLT_DiCentralPFJet170", &HLT_DiCentralPFJet170);
    EventTree->SetBranchAddress("HLT_SingleCentralPFJet170_CFMax0p1", &HLT_SingleCentralPFJet170_CFMax0p1);
    EventTree->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &HLT_DiCentralPFJet170_CFMax0p1);
    EventTree->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &HLT_DiCentralPFJet220_CFMax0p3);
    EventTree->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &HLT_DiCentralPFJet330_CFMax0p5);
    EventTree->SetBranchAddress("HLT_DiCentralPFJet430", &HLT_DiCentralPFJet430);
    EventTree->SetBranchAddress("HLT_PFHT125", &HLT_PFHT125);
    EventTree->SetBranchAddress("HLT_PFHT200", &HLT_PFHT200);
    EventTree->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250);
    EventTree->SetBranchAddress("HLT_PFHT300", &HLT_PFHT300);
    EventTree->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350);
    EventTree->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400);
    EventTree->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475);
    EventTree->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600);
    EventTree->SetBranchAddress("HLT_PFHT650", &HLT_PFHT650);
    EventTree->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800);
    EventTree->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900);
    EventTree->SetBranchAddress("HLT_PFHT200_PFAlphaT0p51", &HLT_PFHT200_PFAlphaT0p51);
    EventTree->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57);
    EventTree->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63);
    EventTree->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55);
    EventTree->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58);
    EventTree->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53);
    EventTree->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54);
    EventTree->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52);
    EventTree->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53);
    EventTree->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51);
    EventTree->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52);
    EventTree->SetBranchAddress("HLT_MET60_IsoTrk35_Loose", &HLT_MET60_IsoTrk35_Loose);
    EventTree->SetBranchAddress("HLT_MET75_IsoTrk50", &HLT_MET75_IsoTrk50);
    EventTree->SetBranchAddress("HLT_MET90_IsoTrk50", &HLT_MET90_IsoTrk50);
    EventTree->SetBranchAddress("HLT_PFMET120_BTagCSV_p067", &HLT_PFMET120_BTagCSV_p067);
    EventTree->SetBranchAddress("HLT_PFMET120_Mu5", &HLT_PFMET120_Mu5);
    EventTree->SetBranchAddress("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned);
    EventTree->SetBranchAddress("HLT_PFMET170_NoiseCleaned", &HLT_PFMET170_NoiseCleaned);
    EventTree->SetBranchAddress("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned);
    EventTree->SetBranchAddress("HLT_PFMET170_JetIdCleaned", &HLT_PFMET170_JetIdCleaned);
    EventTree->SetBranchAddress("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_PFMET170_HBHE_BeamHaloCleaned", &HLT_PFMET170_HBHE_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_PFMET90_PFMHT90_IDTight", &HLT_PFMET90_PFMHT90_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight", &HLT_PFMET100_PFMHT100_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned", &HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned);
    EventTree->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight);
    EventTree->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight);
    EventTree->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067);
    EventTree->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight);
    EventTree->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200);
    EventTree->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460);
    EventTree->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240);
    EventTree->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500);
    EventTree->SetBranchAddress("HLT_QuadPFJet_VBF", &HLT_QuadPFJet_VBF);
    EventTree->SetBranchAddress("HLT_L1_TripleJet_VBF", &HLT_L1_TripleJet_VBF);
    EventTree->SetBranchAddress("HLT_QuadJet45_TripleBTagCSV_p087", &HLT_QuadJet45_TripleBTagCSV_p087);
    EventTree->SetBranchAddress("HLT_QuadJet45_DoubleBTagCSV_p087", &HLT_QuadJet45_DoubleBTagCSV_p087);
    EventTree->SetBranchAddress("HLT_DoubleJet90_Double30_TripleBTagCSV_p087", &HLT_DoubleJet90_Double30_TripleBTagCSV_p087);
    EventTree->SetBranchAddress("HLT_DoubleJet90_Double30_DoubleBTagCSV_p087", &HLT_DoubleJet90_Double30_DoubleBTagCSV_p087);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160", &HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6", &HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6);
    EventTree->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172", &HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172);
    EventTree->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6", &HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026", &HLT_DoubleJetsC100_SingleBTagCSV_p026);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014", &HLT_DoubleJetsC100_SingleBTagCSV_p014);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350);
    EventTree->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350);
    EventTree->SetBranchAddress("HLT_Photon135_PFMET100", &HLT_Photon135_PFMET100);
    EventTree->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL);
    EventTree->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE);
    EventTree->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE);
    EventTree->SetBranchAddress("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60", &HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60);
    EventTree->SetBranchAddress("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15);
    EventTree->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40);
    EventTree->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL);
    EventTree->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
    EventTree->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5", &HLT_BTagMu_DiJet20_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_DiJet40_Mu5", &HLT_BTagMu_DiJet40_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_DiJet70_Mu5", &HLT_BTagMu_DiJet70_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_DiJet110_Mu5", &HLT_BTagMu_DiJet110_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_DiJet170_Mu5", &HLT_BTagMu_DiJet170_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_Jet300_Mu5", &HLT_BTagMu_Jet300_Mu5);
    EventTree->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5);
    EventTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded);
    EventTree->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    EventTree->SetBranchAddress("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL", &HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL", &HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL", &HLT_Mu12_Photon25_CaloIdL);
    EventTree->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1ISO", &HLT_Mu12_Photon25_CaloIdL_L1ISO);
    EventTree->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1OR", &HLT_Mu12_Photon25_CaloIdL_L1OR);
    EventTree->SetBranchAddress("HLT_Mu17_Photon22_CaloIdL_L1ISO", &HLT_Mu17_Photon22_CaloIdL_L1ISO);
    EventTree->SetBranchAddress("HLT_Mu17_Photon30_CaloIdL_L1ISO", &HLT_Mu17_Photon30_CaloIdL_L1ISO);
    EventTree->SetBranchAddress("HLT_Mu17_Photon35_CaloIdL_L1ISO", &HLT_Mu17_Photon35_CaloIdL_L1ISO);
    EventTree->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
    EventTree->SetBranchAddress("HLT_TripleMu_5_3_3", &HLT_TripleMu_5_3_3);
    EventTree->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5);
    EventTree->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125", &HLT_Mu3er_PFHT140_PFMET125);
    EventTree->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067", &HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067);
    EventTree->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100", &HLT_Mu6_PFHT200_PFMET100);
    EventTree->SetBranchAddress("HLT_Mu14er_PFMET100", &HLT_Mu14er_PFMET100);
    EventTree->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Ele17_CaloIdL_GsfTrkIdVL", &HLT_Ele17_CaloIdL_GsfTrkIdVL);
    EventTree->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL);
    EventTree->SetBranchAddress("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ900DEtaJJ1p5);
    EventTree->SetBranchAddress("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ950DEtaJJ1p5);
    EventTree->SetBranchAddress("HLT_Photon22", &HLT_Photon22);
    EventTree->SetBranchAddress("HLT_Photon30", &HLT_Photon30);
    EventTree->SetBranchAddress("HLT_Photon36", &HLT_Photon36);
    EventTree->SetBranchAddress("HLT_Photon50", &HLT_Photon50);
    EventTree->SetBranchAddress("HLT_Photon75", &HLT_Photon75);
    EventTree->SetBranchAddress("HLT_Photon90", &HLT_Photon90);
    EventTree->SetBranchAddress("HLT_Photon120", &HLT_Photon120);
    EventTree->SetBranchAddress("HLT_Photon175", &HLT_Photon175);
    EventTree->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10);
    EventTree->SetBranchAddress("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM);
    EventTree->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
    EventTree->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70);
    EventTree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
    EventTree->SetBranchAddress("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55", &HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55);
    EventTree->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
    EventTree->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon", &HLT_Dimuon0_Jpsi_Muon);
    EventTree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon", &HLT_Dimuon0_Upsilon_Muon);
    EventTree->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Jpsi", &HLT_QuadMuon0_Dimuon0_Jpsi);
    EventTree->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Upsilon", &HLT_QuadMuon0_Dimuon0_Upsilon);
    EventTree->SetBranchAddress("HLT_Rsq0p25_Calo", &HLT_Rsq0p25_Calo);
    EventTree->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo", &HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo);
    EventTree->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_Calo", &HLT_RsqMR240_Rsq0p09_MR200_Calo);
    EventTree->SetBranchAddress("HLT_Rsq0p25", &HLT_Rsq0p25);
    EventTree->SetBranchAddress("HLT_Rsq0p30", &HLT_Rsq0p30);
    EventTree->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200", &HLT_RsqMR240_Rsq0p09_MR200);
    EventTree->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet", &HLT_RsqMR240_Rsq0p09_MR200_4jet);
    EventTree->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200", &HLT_RsqMR270_Rsq0p09_MR200);
    EventTree->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet", &HLT_RsqMR270_Rsq0p09_MR200_4jet);
    EventTree->SetBranchAddress("HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200);
    EventTree->SetBranchAddress("HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
    EventTree->SetBranchAddress("HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
    EventTree->SetBranchAddress("HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
    EventTree->SetBranchAddress("HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
    EventTree->SetBranchAddress("HLT_HT200_DisplacedDijet40_DisplacedTrack", &HLT_HT200_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT250_DisplacedDijet40_DisplacedTrack", &HLT_HT250_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT350_DisplacedDijet40_DisplacedTrack", &HLT_HT350_DisplacedDijet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT350_DisplacedDijet80_DisplacedTrack", &HLT_HT350_DisplacedDijet80_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack", &HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_HT350_DisplacedDijet40_Inclusive", &HLT_HT350_DisplacedDijet40_Inclusive);
    EventTree->SetBranchAddress("HLT_HT400_DisplacedDijet40_Inclusive", &HLT_HT400_DisplacedDijet40_Inclusive);
    EventTree->SetBranchAddress("HLT_HT500_DisplacedDijet40_Inclusive", &HLT_HT500_DisplacedDijet40_Inclusive);
    EventTree->SetBranchAddress("HLT_HT550_DisplacedDijet40_Inclusive", &HLT_HT550_DisplacedDijet40_Inclusive);
    EventTree->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive);
    EventTree->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive);
    EventTree->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack", &HLT_VBF_DisplacedJet40_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5", &HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_TightID_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic", &HLT_VBF_DisplacedJet40_Hadronic);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack", &HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_Hadronic", &HLT_VBF_DisplacedJet40_TightID_Hadronic);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_Hadronic", &HLT_VBF_DisplacedJet40_VTightID_Hadronic);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_Hadronic", &HLT_VBF_DisplacedJet40_VVTightID_Hadronic);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack);
    EventTree->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_PFMETNoMu90_PFMHTNoMu90_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
    EventTree->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
    EventTree->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
    EventTree->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200);
    EventTree->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500);
    EventTree->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT250", &HLT_DoubleMu8_Mass8_PFHT250);
    EventTree->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250);
    EventTree->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250);
    EventTree->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT300", &HLT_DoubleMu8_Mass8_PFHT300);
    EventTree->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300);
    EventTree->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300);
    EventTree->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV_p13", &HLT_Mu10_CentralPFJet30_BTagCSV_p13);
    EventTree->SetBranchAddress("HLT_DoubleMu3_PFMET50", &HLT_DoubleMu3_PFMET50);
    EventTree->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350", &HLT_Ele15_IsoVVVL_PFHT350);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400_PFMET50", &HLT_Ele15_IsoVVVL_PFHT400_PFMET50);
    EventTree->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400", &HLT_Ele15_IsoVVVL_PFHT400);
    EventTree->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT400", &HLT_Ele50_IsoVVVL_PFHT400);
    EventTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
    EventTree->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350", &HLT_Mu15_IsoVVVL_PFHT350);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400_PFMET50", &HLT_Mu15_IsoVVVL_PFHT400_PFMET50);
    EventTree->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400", &HLT_Mu15_IsoVVVL_PFHT400);
    EventTree->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT400", &HLT_Mu50_IsoVVVL_PFHT400);
    EventTree->SetBranchAddress("HLT_Dimuon16_Jpsi", &HLT_Dimuon16_Jpsi);
    EventTree->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", &HLT_Dimuon10_Jpsi_Barrel);
    EventTree->SetBranchAddress("HLT_Dimuon8_PsiPrime_Barrel", &HLT_Dimuon8_PsiPrime_Barrel);
    EventTree->SetBranchAddress("HLT_Dimuon8_Upsilon_Barrel", &HLT_Dimuon8_Upsilon_Barrel);
    EventTree->SetBranchAddress("HLT_Dimuon0_Phi_Barrel", &HLT_Dimuon0_Phi_Barrel);
    EventTree->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Onia", &HLT_Mu16_TkMu0_dEta18_Onia);
    EventTree->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Phi", &HLT_Mu16_TkMu0_dEta18_Phi);
    EventTree->SetBranchAddress("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
    EventTree->SetBranchAddress("HLT_Mu8", &HLT_Mu8);
    EventTree->SetBranchAddress("HLT_Mu17", &HLT_Mu17);
    EventTree->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40);
    EventTree->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
    EventTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140);
    EventTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
    EventTree->SetBranchAddress("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", &HLT_PFHT400_SixJet30_DoubleBTagCSV_p056);
    EventTree->SetBranchAddress("HLT_PFHT450_SixJet40_BTagCSV_p056", &HLT_PFHT450_SixJet40_BTagCSV_p056);
    EventTree->SetBranchAddress("HLT_PFHT400_SixJet30", &HLT_PFHT400_SixJet30);
    EventTree->SetBranchAddress("HLT_PFHT450_SixJet40", &HLT_PFHT450_SixJet40);
    EventTree->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Mu55", &HLT_Mu55);
    EventTree->SetBranchAddress("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15);
    EventTree->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT600", &HLT_Photon90_CaloIdL_PFHT600);
    EventTree->SetBranchAddress("HLT_PixelTracks_Multiplicity60ForEndOfFill", &HLT_PixelTracks_Multiplicity60ForEndOfFill);
    EventTree->SetBranchAddress("HLT_PixelTracks_Multiplicity85ForEndOfFill", &HLT_PixelTracks_Multiplicity85ForEndOfFill);
    EventTree->SetBranchAddress("HLT_PixelTracks_Multiplicity110ForEndOfFill", &HLT_PixelTracks_Multiplicity110ForEndOfFill);
    EventTree->SetBranchAddress("HLT_PixelTracks_Multiplicity135ForEndOfFill", &HLT_PixelTracks_Multiplicity135ForEndOfFill);
    EventTree->SetBranchAddress("HLT_PixelTracks_Multiplicity160ForEndOfFill", &HLT_PixelTracks_Multiplicity160ForEndOfFill);
    EventTree->SetBranchAddress("HLT_FullTracks_Multiplicity80", &HLT_FullTracks_Multiplicity80);
    EventTree->SetBranchAddress("HLT_FullTracks_Multiplicity100", &HLT_FullTracks_Multiplicity100);
    EventTree->SetBranchAddress("HLT_FullTracks_Multiplicity130", &HLT_FullTracks_Multiplicity130);
    EventTree->SetBranchAddress("HLT_FullTracks_Multiplicity150", &HLT_FullTracks_Multiplicity150);
    EventTree->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800);
    EventTree->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70);
    EventTree->SetBranchAddress("HLT_Photon125", &HLT_Photon125);
    EventTree->SetBranchAddress("HLT_MET100", &HLT_MET100);
    EventTree->SetBranchAddress("HLT_MET150", &HLT_MET150);
    EventTree->SetBranchAddress("HLT_MET200", &HLT_MET200);
    EventTree->SetBranchAddress("HLT_Ele27_HighEta_Ele20_Mass55", &HLT_Ele27_HighEta_Ele20_Mass55);
    EventTree->SetBranchAddress("HLT_L1FatEvents", &HLT_L1FatEvents);
    EventTree->SetBranchAddress("HLT_Physics", &HLT_Physics);
    EventTree->SetBranchAddress("HLT_L1FatEvents_part0", &HLT_L1FatEvents_part0);
    EventTree->SetBranchAddress("HLT_L1FatEvents_part1", &HLT_L1FatEvents_part1);
    EventTree->SetBranchAddress("HLT_L1FatEvents_part2", &HLT_L1FatEvents_part2);
    EventTree->SetBranchAddress("HLT_L1FatEvents_part3", &HLT_L1FatEvents_part3);
    EventTree->SetBranchAddress("HLT_Random", &HLT_Random);
    EventTree->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias);
    EventTree->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30);
    EventTree->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40);
    EventTree->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50);
    EventTree->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80);
    EventTree->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100);
    EventTree->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30);
    EventTree->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50);
    EventTree->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80);
    EventTree->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100);
    EventTree->SetBranchAddress("HLT_HISinglePhoton10", &HLT_HISinglePhoton10);
    EventTree->SetBranchAddress("HLT_HISinglePhoton15", &HLT_HISinglePhoton15);
    EventTree->SetBranchAddress("HLT_HISinglePhoton20", &HLT_HISinglePhoton20);
    EventTree->SetBranchAddress("HLT_HISinglePhoton40", &HLT_HISinglePhoton40);
    EventTree->SetBranchAddress("HLT_HISinglePhoton60", &HLT_HISinglePhoton60);
    EventTree->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration);
    EventTree->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration);
    EventTree->SetBranchAddress("HLT_GlobalRunHPDNoise", &HLT_GlobalRunHPDNoise);
    EventTree->SetBranchAddress("HLT_L1BptxMinus", &HLT_L1BptxMinus);
    EventTree->SetBranchAddress("HLT_L1BptxPlus", &HLT_L1BptxPlus);
    EventTree->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR);
    EventTree->SetBranchAddress("HLT_L1BeamGasMinus", &HLT_L1BeamGasMinus);
    EventTree->SetBranchAddress("HLT_L1BeamGasPlus", &HLT_L1BeamGasPlus);
    EventTree->SetBranchAddress("HLT_L1BptxXOR", &HLT_L1BptxXOR);
    EventTree->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR);
    EventTree->SetBranchAddress("HLT_L1MinimumBiasHF_AND", &HLT_L1MinimumBiasHF_AND);
    EventTree->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS);
    EventTree->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym);
    EventTree->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_copy", &HLT_ZeroBias_FirstCollisionAfterAbortGap_copy);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS", &HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS);
    EventTree->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain);
    EventTree->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain);
    EventTree->SetBranchAddress("HLT_Photon500", &HLT_Photon500);
    EventTree->SetBranchAddress("HLT_Photon600", &HLT_Photon600);
    EventTree->SetBranchAddress("HLT_Mu300", &HLT_Mu300);
    EventTree->SetBranchAddress("HLT_Mu350", &HLT_Mu350);
    EventTree->SetBranchAddress("HLT_MET250", &HLT_MET250);
    EventTree->SetBranchAddress("HLT_MET300", &HLT_MET300);
    EventTree->SetBranchAddress("HLT_MET600", &HLT_MET600);
    EventTree->SetBranchAddress("HLT_MET700", &HLT_MET700);
    EventTree->SetBranchAddress("HLT_PFMET300", &HLT_PFMET300);
    EventTree->SetBranchAddress("HLT_PFMET400", &HLT_PFMET400);
    EventTree->SetBranchAddress("HLT_PFMET500", &HLT_PFMET500);
    EventTree->SetBranchAddress("HLT_PFMET600", &HLT_PFMET600);
    EventTree->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT);
    EventTree->SetBranchAddress("HLT_HT2000", &HLT_HT2000);
    EventTree->SetBranchAddress("HLT_HT2500", &HLT_HT2500);
    EventTree->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE);
    EventTree->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB);
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
    EventTree->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
    EventTree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    EventTree->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter);
    EventTree->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters);
    EventTree->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter);
    EventTree->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter);
    EventTree->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X);
    EventTree->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X);
    EventTree->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters);
    EventTree->SetBranchAddress("Flag_METFilters", &Flag_METFilters);
    EvMax=EventTree->GetEntries();

    cout<<cur_time()<<"\tTree successfully processed!\n";

    EventLoop(sampleIdx);
    infile->Close();
    //score_infile->Close();

    return;
}

void RunAnalysis(){

   const char* inDir[3];
   char* dir_xroot[3];
   char* dir[3];
   void* dirp[3];

   inDir[0] = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-20to40_OptPFRHT_OptMustLocalParams";
   dir_xroot[0] = "root://eoscms.cern.ch///eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-20to40_OptPFRHT_OptMustLocalParams";
   inDir[1] = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-20toInf_OptPFRHT_OptMustLocalParams";
   dir_xroot[1] = "root://eoscms.cern.ch///eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-20toInf_OptPFRHT_OptMustLocalParams";
   inDir[2] = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-40toInf_OptPFRHT_OptMustLocalParams";
   dir_xroot[2] = "root://eoscms.cern.ch///eos/cms/store/group/dpg_ecal/alca_ecalcalib/Clustering/lzygala/GJets_sample_Pt-40toInf_OptPFRHT_OptMustLocalParams";
   

    for(int j=0; j<nSamples; j++){

        //if points to single file
        if(0 == sampleAddresses.at(j).compare(sampleAddresses.at(j).length() - 5, 5, ".root")){
            InitTree(sampleAddresses.at(j), j);
        }

        //else points to directory
        else{
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
                //Printf("\tfile -> %s", filename[i]);
                Printf("\n%s\tfile -> %i / %i", (cur_time()).c_str(), i, n);
                InitTree(filename[i],weights[j]);
                //EventLoop();
            }

        }

    }
}

void fullAnalyzer(string signalFile, string signalName, string backgroundFiles[], string backgroundNames[]){
//main program
    nSamples = 1 + sizeof(backgroundFiles);
    sampleAddresses.push_back(signalFile);
    sampleNames.push_back(signalName);

    copy(&backgroundFiles[0], &backgroundFiles[sizeof(backgroundFiles)], back_inserter(sampleAddresses));
    copy(&backgroundNames[0], &backgroundNames[sizeof(backgroundNames)], back_inserter(sampleNames));


    InitHistograms();
    InitTree(inputFile.c_str()); 
    SaveHistograms();
}
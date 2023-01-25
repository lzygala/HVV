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
#include<TTreeReader.h>
#include<TTreeReaderArray.h>
#include<TTreeReaderValue.h>

#include<algorithm>
#include<chrono>
#include<ctime>
#include<fstream>
#include<initializer_list>
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

ofstream outfile;

void EventLoop(TString infileName){

    bool debug = false;

    cout<<"\tProcessing Tree...\n";
    cout<<"\tFile Opened...\n";
    TFile* infile = TFile::Open(infileName);

    TTreeReader myReader("Events", infile);

    TTreeReaderValue<UInt_t> run(myReader, "run");
    TTreeReaderValue<UInt_t> luminosityBlock(myReader, "luminosityBlock");
    TTreeReaderValue<ULong64_t> event(myReader, "event");
    TTreeReaderValue<Float_t> CaloMET_phi(myReader, "CaloMET_phi");
    TTreeReaderValue<Float_t> CaloMET_pt(myReader, "CaloMET_pt");
    TTreeReaderValue<Float_t> CaloMET_sumEt(myReader, "CaloMET_sumEt");
    TTreeReaderValue<Float_t> ChsMET_phi(myReader, "ChsMET_phi");
    TTreeReaderValue<Float_t> ChsMET_pt(myReader, "ChsMET_pt");
    TTreeReaderValue<Float_t> ChsMET_sumEt(myReader, "ChsMET_sumEt");
    TTreeReaderValue<UInt_t> nCorrT1METJet(myReader, "nCorrT1METJet");
    TTreeReaderArray<Float_t> CorrT1METJet_area(myReader, "CorrT1METJet_area");
    TTreeReaderArray<Float_t> CorrT1METJet_eta(myReader, "CorrT1METJet_eta");
    TTreeReaderArray<Float_t> CorrT1METJet_muonSubtrFactor(myReader, "CorrT1METJet_muonSubtrFactor");
    TTreeReaderArray<Float_t> CorrT1METJet_phi(myReader, "CorrT1METJet_phi");
    TTreeReaderArray<Float_t> CorrT1METJet_rawPt(myReader, "CorrT1METJet_rawPt");
    TTreeReaderValue<Float_t> TkMET_phi(myReader, "TkMET_phi");
    TTreeReaderValue<Float_t> TkMET_pt(myReader, "TkMET_pt");
    TTreeReaderValue<Float_t> TkMET_sumEt(myReader, "TkMET_sumEt");

    TTreeReaderValue<UInt_t> nGenPart(myReader, "nGenPart");
    TTreeReaderArray<Float_t> GenPart_eta(myReader, "GenPart_eta");
    TTreeReaderArray<Float_t> GenPart_mass(myReader, "GenPart_mass");
    TTreeReaderArray<Float_t> GenPart_phi(myReader, "GenPart_phi");
    TTreeReaderArray<Float_t> GenPart_pt(myReader, "GenPart_pt");
    TTreeReaderArray<Int_t> GenPart_genPartIdxMother(myReader, "GenPart_genPartIdxMother");
    TTreeReaderArray<Int_t> GenPart_pdgId(myReader, "GenPart_pdgId");
    TTreeReaderArray<Int_t> GenPart_status(myReader, "GenPart_status");
    TTreeReaderArray<Int_t> GenPart_statusFlags(myReader, "GenPart_statusFlags");
    TTreeReaderValue<Float_t> genWeight(myReader, "genWeight");
    TTreeReaderValue<UInt_t> nPSWeight(myReader, "nPSWeight");
    TTreeReaderArray<Float_t> PSWeight(myReader, "PSWeight");
    //if(sampleIdx < 2){
       // TTreeReaderValue<Float_t> LHE_HT(myReader, "LHE_HT");
       // TTreeReaderValue<Float_t> LHE_HTIncoming(myReader, "LHE_HTIncoming");
        //TTreeReaderValue<Float_t> LHE_Vpt(myReader, "LHE_Vpt");
        //TTreeReaderValue<UChar_t> LHE_Njets(myReader, "LHE_Njets");
        //TTreeReaderValue<UChar_t> LHE_Nb(myReader, "LHE_Nb");
        //TTreeReaderValue<UChar_t> LHE_Nc(myReader, "LHE_Nc");
        //TTreeReaderValue<UChar_t> LHE_Nuds(myReader, "LHE_Nuds");
        //TTreeReaderValue<UChar_t> LHE_Nglu(myReader, "LHE_Nglu");
        //TTreeReaderValue<UChar_t> LHE_NpNLO(myReader, "LHE_NpNLO");
       // TTreeReaderValue<UChar_t> LHE_NpLO(myReader, "LHE_NpLO");
        TTreeReaderValue<UInt_t> nLHEPart(myReader, "nLHEPart");
        TTreeReaderArray<Float_t> LHEPart_pt(myReader, "LHEPart_pt");
        TTreeReaderArray<Float_t> LHEPart_eta(myReader, "LHEPart_eta");
        TTreeReaderArray<Float_t> LHEPart_phi(myReader, "LHEPart_phi");
        TTreeReaderArray<Float_t> LHEPart_mass(myReader, "LHEPart_mass");
        TTreeReaderArray<Int_t> LHEPart_pdgId(myReader, "LHEPart_pdgId");
        TTreeReaderValue<Float_t> LHEWeight_originalXWGTUP(myReader, "LHEWeight_originalXWGTUP");
        TTreeReaderValue<UInt_t> nLHEPdfWeight(myReader, "nLHEPdfWeight");
        TTreeReaderArray<Float_t> LHEPdfWeight(myReader, "LHEPdfWeight");
        TTreeReaderValue<UInt_t> nLHEReweightingWeight(myReader, "nLHEReweightingWeight");
        TTreeReaderArray<Float_t> LHEReweightingWeight(myReader, "LHEReweightingWeight");
        TTreeReaderValue<UInt_t> nLHEScaleWeight(myReader, "nLHEScaleWeight");
    //}

    TTreeReaderValue<Float_t> GenMET_phi(myReader, "GenMET_phi");
    TTreeReaderValue<Float_t> GenMET_pt(myReader, "GenMET_pt");
    TTreeReaderValue<Float_t> MET_phi(myReader, "MET_phi");
    TTreeReaderValue<Float_t> MET_pt(myReader, "MET_pt");
    TTreeReaderValue<Float_t> MET_sumEt(myReader, "MET_sumEt");
    TTreeReaderValue<Float_t> PuppiMET_phi(myReader, "PuppiMET_phi");
    TTreeReaderValue<Float_t> PuppiMET_pt(myReader, "PuppiMET_pt");
    TTreeReaderValue<Float_t> PuppiMET_sumEt(myReader, "PuppiMET_sumEt");
    TTreeReaderValue<Float_t> RawMET_phi(myReader, "RawMET_phi");
    TTreeReaderValue<Float_t> RawMET_pt(myReader, "RawMET_pt");
    TTreeReaderValue<Float_t> RawMET_sumEt(myReader, "RawMET_sumEt");


    //int eventLoopMax = EvMax;
    float jetPT_cut = 50.0;
    float jetInvMass_cut = 400.0;
    float jetEtaSep_cut = 2.0;

    Long64_t eventLoopMax = myReader.GetEntries();
    int iev = 0;
    int iLHEPassed = 0;

    //BEGIN EVENT LOOP
    while (myReader.Next()) {
        if(++iev % 10000 == 0) cout<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;
        if(iev > 10) break;

        outfile<<"----------------------------------------------------------------------------------------------------------------\n";
        outfile<<"EVENT #"<<iev<<endl;

        for(int iGenPart=0; iGenPart<*nGenPart; iGenPart++){
            outfile<<"GenPart Idx="<<iGenPart<<"\tpdgid="<<GenPart_pdgId[iGenPart]<<"\tstatus="<<GenPart_status[iGenPart]<<"\tmotherIdx="<<GenPart_genPartIdxMother[iGenPart]<<"\teta="<<GenPart_eta[iGenPart]<<"\tphi="<<GenPart_phi[iGenPart]<<"\tpt="<<GenPart_pt[iGenPart]<<endl;

        }

        outfile<<"----------------------------------------------------------------------------------------------------------------\n\n\n";
    }
}

void GenLevel_info(std::string inputFile){

    outfile.open("./events_kappa_c2v_1.txt", std::ofstream::out | std::ofstream::trunc);

    EventLoop(inputFile);

    outfile.close();
}
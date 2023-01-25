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

using namespace std;

float pi=3.1415927;
float twopi= 2*pi;

int nFailedGenMatchedDecay = 0, 
    nPassGenMatchedDecay_FailHbb = 0,  
    nPassGenMatchedDecay_FailZ= 0,   
    nPassGenMatchedDecay_FailHbb_FailZ= 0, 
    nPassGenMatchedDecay_Total= 0,         
    nPOMatched = 0, 
    nNotPOMatched = 0, 
    nFailPOMatched_FailVBF = 0, 
    nFailPOMatched_FailH = 0, 
    nFailPOMatched_FailZ = 0, 
    nFailPOMatched_FailWp = 0, 
    nFailPOMatched_FailWm = 0;

TTree* Events_basic;

Int_t Higgs_LHEPart_idx;
Int_t Wp_LHEPart_idx;  
Int_t Wm_LHEPart_idx;
Int_t Z_LHEPart_idx;
Int_t VBFLead_LHEPart_idx;
Int_t VBFTrail_LHEPart_idx;

Bool_t GenMatching_Initial_FullMatch;
Bool_t GenMatching_Decayed_FullMatch;
Bool_t POMatching_FullMatch;

Int_t Higgs_GenPart_idx;
Bool_t Higgs_DecayTobbbar;
Int_t Higgs_GenPart_decay_type;  //0=qq, 
Int_t Higgs_GenPart_b_idx;
Int_t Higgs_GenPart_bbar_idx;

Int_t VBFLead_GenPart_idx;
Int_t VBFTrail_GenPart_idx;

Int_t Wp_GenPart_idx;
Int_t Wp_GenPart_decay_type;  //0=qq, 1=enu, 2=munu, 3=taunu
Int_t Wp_GenPart_daughterLead_idx;
Int_t Wp_GenPart_daughterLead_pdgid;
Int_t Wp_GenPart_daughterTrail_idx;
Int_t Wp_GenPart_daughterTrail_pdgid;

Int_t Wm_GenPart_idx;
Int_t Wm_GenPart_decay_type;   //0=qq, 1=enu, 2=munu, 3=taunu
Int_t Wm_GenPart_daughterLead_idx;
Int_t Wm_GenPart_daughterLead_pdgid;
Int_t Wm_GenPart_daughterTrail_idx;
Int_t Wm_GenPart_daughterTrail_pdgid;

Int_t Z_GenPart_idx;
Int_t Z_GenPart_decay_type;    //0=qq, 1=nunu, 2=ee, 3=mumu, 4=tautau
Int_t Z_GenPart_daughterLead_idx;
Int_t Z_GenPart_daughterLead_pdgid;
Int_t Z_GenPart_daughterTrail_idx;
Int_t Z_GenPart_daughterTrail_pdgid;




Int_t Higgs_Matched_AK8_idx;
Int_t Higgs_b_Matched_AK4_idx;
Int_t Higgs_bbar_Matched_AK4_idx;

Int_t VBFLead_Matched_AK4_idx;
Int_t VBFTrail_Matched_AK4_idx;

Int_t Wp_Matched_AK8_idx;
Int_t Wp_Matched_Electron_idx;
Int_t Wp_Matched_Muon_idx;

Int_t Wm_Matched_AK8_idx;
Int_t Wm_Matched_Electron_idx;
Int_t Wm_Matched_Muon_idx;

Int_t Z_Matched_AK8_idx;
Int_t Z_Matched_Electron_1_idx;
Int_t Z_Matched_Electron_2_idx;
Int_t Z_Matched_Muon_1_idx;
Int_t Z_Matched_Muon_2_idx;



Float_t VBF_LHEPart_Mjj;
Float_t VBF_LHEPart_dEta;
Float_t VBF_LHEPart_dR;

Float_t Higgs_GenPart_bb_Mjj;
Float_t Higgs_GenPart_bb_dEta;
Float_t Higgs_GenPart_bb_dR;

Float_t VBF_GenPart_Mjj;
Float_t VBF_GenPart_dEta;
Float_t VBF_GenPart_dR;

Float_t Higgs_bb_Matched_AK4_Mjj;
Float_t Higgs_bb_Matched_AK4_dEta;
Float_t Higgs_bb_Matched_AK4_dR;

Float_t VBF_Matched_AK4_Mjj;
Float_t VBF_Matched_AK4_dEta;
Float_t VBF_Matched_AK4_dR;


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

template<int M, template<typename> class F = std::greater>
struct TupleCompare
{
    template<typename T>
    bool operator()(T const &t1, T const &t2)
    {
        return F<typename tuple_element<M, T>::type>()(std::get<M>(t1), std::get<M>(t2));
    }
};

template <typename T>
bool is_in(const T& v, std::initializer_list<T> lst)
{
    return std::find(std::begin(lst), std::end(lst), v) != std::end(lst);
}

float getDR(float phi_1, float eta_1, float phi_2, float eta_2){
    float dPhi = fabs(phi_1 - phi_2);
    if(dPhi > pi) dPhi -= twopi;

    float dEta = eta_1 - eta_2;

    float dR = sqrt( pow(dPhi, 2) + pow(dEta, 2) );

    return dR;

}

vector<std::tuple<int,int>> getDecayProducts(int genIdx, vector<int> ttra_pdgId, vector<int> ttra_genPartIdxMother){
    //genIdx = index of an initial particle, want to go through chain and find the first decay
    bool foundLastInChain = false; //directly before decay happens
    int lastGenIdx=-1;
    int currIdx=genIdx;
    while(!foundLastInChain){

        for(int ii_genPart=currIdx; ii_genPart<ttra_pdgId.size(); ii_genPart++){
            if(ttra_genPartIdxMother[ii_genPart] == currIdx && ttra_pdgId[ii_genPart] == ttra_pdgId[currIdx]){
                currIdx = ii_genPart;
                break;
            }
            else if(ttra_genPartIdxMother[ii_genPart] == currIdx){
                foundLastInChain = true;
                break;
            }
        }
    }

    vector<std::tuple<int,int>> directDecayProducts;
    for(int ii_genPart=genIdx; ii_genPart<ttra_pdgId.size(); ii_genPart++){
        if(ttra_genPartIdxMother[ii_genPart] == currIdx)
            directDecayProducts.push_back(std::make_tuple(ii_genPart,ttra_pdgId[ii_genPart]));
    }

    return directDecayProducts;
}

int getFinalState(int genIdx, vector<int> ttra_pdgId, vector<int> ttra_genPartIdxMother, vector<int> ttra_genPartStatus){
    //genIdx = index of an initial particle, want to go through chain and find the first decay
    bool foundLastInChain = false; //directly before decay happens
    int FSIdx=-1;
    int currIdx=genIdx;
    while(!foundLastInChain){

        for(int ii_genPart=currIdx; ii_genPart<ttra_pdgId.size(); ii_genPart++){
            if(ttra_genPartIdxMother[ii_genPart] == currIdx && ttra_pdgId[ii_genPart] == ttra_pdgId[currIdx]){
                currIdx = ii_genPart;
                break;
            }
            if(ii_genPart == ttra_pdgId.size()-1){
                foundLastInChain = true;
                break;
            }
        }
    }

    if(ttra_genPartStatus[currIdx] == 1) FSIdx = currIdx;

    return FSIdx;
}


int getDecayType(int d1, int d2){
    if(fabs(d1) < 9 && fabs(d2) < 9) return 0;
    if((fabs(d1) == 11 && fabs(d2) == 12) || (fabs(d1) == 12 && fabs(d2) == 11)) return 1; //e nu
    if((fabs(d1) == 13 && fabs(d2) == 14) || (fabs(d1) == 14 && fabs(d2) == 13)) return 2; // mu nu
    if((fabs(d1) == 15 && fabs(d2) == 16) || (fabs(d1) == 16 && fabs(d2) == 15)) return 3; // tau nu
    if((fabs(d1) == 11 && fabs(d2) == 11)) return 4; //e nu
    if((fabs(d1) == 13 && fabs(d2) == 13)) return 5; //e nu
    if((fabs(d1) == 15 && fabs(d2) == 15)) return 6; //e nu
    if((fabs(d1) == 12 && fabs(d2) == 12) || (fabs(d1) == 14 && fabs(d2) == 14) || (fabs(d1) == 16 && fabs(d2) == 16)) return 7; //e nu

    return -1;
}

void AddBranches(){
    //Candidate Branches in Events Tree
    Events_basic->Branch("Higgs_LHEPart_idx", &Higgs_LHEPart_idx, "Higgs_LHEPart_idx/I");
    Events_basic->Branch("Wp_LHEPart_idx", &Wp_LHEPart_idx, "Wp_LHEPart_idx/I");  
    Events_basic->Branch("Wm_LHEPart_idx", &Wm_LHEPart_idx, "Wm_LHEPart_idx/I");
    Events_basic->Branch("Z_LHEPart_idx", &Z_LHEPart_idx, "Z_LHEPart_idx/I");
    Events_basic->Branch("VBFLead_LHEPart_idx", &VBFLead_LHEPart_idx, "VBFLead_LHEPart_idx/I");
    Events_basic->Branch("VBFTrail_LHEPart_idx", &VBFTrail_LHEPart_idx, "VBFTrail_LHEPart_idx/I");


    Events_basic->Branch("GenMatching_Initial_FullMatch", &GenMatching_Initial_FullMatch, "GenMatching_Initial_FullMatch/B");
    Events_basic->Branch("GenMatching_Decayed_FullMatch", &GenMatching_Decayed_FullMatch, "GenMatching_Decayed_FullMatch/B");
    Events_basic->Branch("POMatching_FullMatch", &POMatching_FullMatch, "POMatching_FullMatch/B");

    Events_basic->Branch("Higgs_GenPart_idx", &Higgs_GenPart_idx, "Higgs_GenPart_idx/I");
    Events_basic->Branch("Higgs_DecayTobbbar", &Higgs_DecayTobbbar, "Higgs_DecayTobbbar/B");
    Events_basic->Branch("Higgs_GenPart_decay_type", &Higgs_GenPart_decay_type, "Higgs_GenPart_decay_type/I");
    Events_basic->Branch("Higgs_GenPart_b_idx", &Higgs_GenPart_b_idx, "Higgs_GenPart_b_idx/I");
    Events_basic->Branch("Higgs_GenPart_bbar_idx", &Higgs_GenPart_bbar_idx, "Higgs_GenPart_bbar_idx/I");

    Events_basic->Branch("VBFLead_GenPart_idx", &VBFLead_GenPart_idx, "VBFLead_GenPart_idx/I");
    Events_basic->Branch("VBFTrail_GenPart_idx", &VBFTrail_GenPart_idx, "VBFTrail_GenPart_idx/I");

    Events_basic->Branch("Wp_GenPart_idx", &Wp_GenPart_idx, "Wp_GenPart_idx/I");
    Events_basic->Branch("Wp_GenPart_decay_type", &Wp_GenPart_decay_type, "Wp_GenPart_decay_type/I");   //0=qq, 1=enu, 2=munu, 3=taunu
    Events_basic->Branch("Wp_GenPart_daughterLead_idx", &Wp_GenPart_daughterLead_idx, "Wp_GenPart_daughterLead_idx/I");
    Events_basic->Branch("Wp_GenPart_daughterLead_pdgid", &Wp_GenPart_daughterLead_pdgid, "Wp_GenPart_daughterLead_pdgid/I");
    Events_basic->Branch("Wp_GenPart_daughterTrail_idx", &Wp_GenPart_daughterTrail_idx, "Wp_GenPart_daughterTrail_idx/I");
    Events_basic->Branch("Wp_GenPart_daughterTrail_pdgid", &Wp_GenPart_daughterTrail_pdgid, "Wp_GenPart_daughterTrail_pdgid/I");

    Events_basic->Branch("Wm_GenPart_idx", &Wm_GenPart_idx, "Wm_GenPart_idx/I");
    Events_basic->Branch("Wm_GenPart_decay_type", &Wm_GenPart_decay_type, "Wm_GenPart_decay_type/I");   //0=qq, 1=enu, 2=munu, 3=taunu
    Events_basic->Branch("Wm_GenPart_daughterLead_idx", &Wm_GenPart_daughterLead_idx, "Wm_GenPart_daughterLead_idx/I");
    Events_basic->Branch("Wm_GenPart_daughterLead_pdgid", &Wm_GenPart_daughterLead_pdgid, "Wm_GenPart_daughterLead_pdgid/I");
    Events_basic->Branch("Wm_GenPart_daughterTrail_idx", &Wm_GenPart_daughterTrail_idx, "Wm_GenPart_daughterTrail_idx/I");
    Events_basic->Branch("Wm_GenPart_daughterTrail_pdgid", &Wm_GenPart_daughterTrail_pdgid, "Wm_GenPart_daughterTrail_pdgid/I");

    Events_basic->Branch("Z_GenPart_idx", &Z_GenPart_idx, "Z_GenPart_idx/I");
    Events_basic->Branch("Z_GenPart_decay_type", &Z_GenPart_decay_type, "Z_GenPart_decay_type/I");    //0=qq, 1=nunu, 2=ee, 3=mumu, 4=tautau
    Events_basic->Branch("Z_GenPart_daughterLead_idx", &Z_GenPart_daughterLead_idx, "Z_GenPart_daughterLead_idx/I");
    Events_basic->Branch("Z_GenPart_daughterLead_pdgid", &Z_GenPart_daughterLead_pdgid, "Z_GenPart_daughterLead_pdgid/I");
    Events_basic->Branch("Z_GenPart_daughterTrail_idx", &Z_GenPart_daughterTrail_idx, "Z_GenPart_daughterTrail_idx/I");
    Events_basic->Branch("Z_GenPart_daughterTrail_pdgid", &Z_GenPart_daughterTrail_pdgid, "Z_GenPart_daughterTrail_pdgid/I");




    Events_basic->Branch("Higgs_Matched_AK8_idx", &Higgs_Matched_AK8_idx, "Higgs_Matched_AK8_idx/I");
    Events_basic->Branch("Higgs_b_Matched_AK4_idx", &Higgs_b_Matched_AK4_idx, "Higgs_b_Matched_AK4_idx/I");
    Events_basic->Branch("Higgs_bbar_Matched_AK4_idx", &Higgs_bbar_Matched_AK4_idx, "Higgs_bbar_Matched_AK4_idx/I");

    Events_basic->Branch("VBFLead_Matched_AK4_idx", &VBFLead_Matched_AK4_idx, "VBFLead_Matched_AK4_idx/I");
    Events_basic->Branch("VBFTrail_Matched_AK4_idx", &VBFTrail_Matched_AK4_idx, "VBFTrail_Matched_AK4_idx/I");

    Events_basic->Branch("Wp_Matched_AK8_idx", &Wp_Matched_AK8_idx, "Wp_Matched_AK8_idx/I");
    Events_basic->Branch("Wp_Matched_Electron_idx", &Wp_Matched_Electron_idx, "Wp_Matched_Electron_idx/I");
    Events_basic->Branch("Wp_Matched_Muon_idx", &Wp_Matched_Muon_idx, "Wp_Matched_Muon_idx/I");

    Events_basic->Branch("Wm_Matched_AK8_idx", &Wm_Matched_AK8_idx, "Wm_Matched_AK8_idx/I");
    Events_basic->Branch("Wm_Matched_Electron_idx", &Wm_Matched_Electron_idx, "Wm_Matched_Electron_idx/I");
    Events_basic->Branch("Wm_Matched_Muon_idx", &Wm_Matched_Muon_idx, "Wm_Matched_Muon_idx/I");



    Events_basic->Branch("Z_Matched_AK8_idx", &Z_Matched_AK8_idx, "Z_Matched_AK8_idx/I");
    Events_basic->Branch("Z_Matched_Electron_1_idx", &Z_Matched_Electron_1_idx, "Z_Matched_Electron_1_idx/I");
    Events_basic->Branch("Z_Matched_Electron_2_idx", &Z_Matched_Electron_2_idx, "Z_Matched_Electron_2_idx/I");
    Events_basic->Branch("Z_Matched_Muon_1_idx", &Z_Matched_Muon_1_idx, "Z_Matched_Muon_1_idx/I");
    Events_basic->Branch("Z_Matched_Muon_2_idx", &Z_Matched_Muon_2_idx, "Z_Matched_Muon_2_idx/I");




    Events_basic->Branch("VBF_LHEPart_Mjj", &VBF_LHEPart_Mjj, "VBF_LHEPart_Mjj/F");
    Events_basic->Branch("VBF_LHEPart_dEta", &VBF_LHEPart_dEta, "VBF_LHEPart_dEta/F");
    Events_basic->Branch("VBF_LHEPart_dR", &VBF_LHEPart_dR, "VBF_LHEPart_dR/F");

    Events_basic->Branch("Higgs_GenPart_bb_Mjj", &Higgs_GenPart_bb_Mjj, "Higgs_GenPart_bb_Mjj/F");
    Events_basic->Branch("Higgs_GenPart_bb_dEta", &Higgs_GenPart_bb_dEta, "Higgs_GenPart_bb_dEta/F");
    Events_basic->Branch("Higgs_GenPart_bb_dR", &Higgs_GenPart_bb_dR, "Higgs_GenPart_bb_dR/F");


    Events_basic->Branch("VBF_GenPart_Mjj", &VBF_GenPart_Mjj, "VBF_GenPart_Mjj/F");
    Events_basic->Branch("VBF_GenPart_dEta", &VBF_GenPart_dEta, "VBF_GenPart_dEta/F");
    Events_basic->Branch("VBF_GenPart_dR", &VBF_GenPart_dR, "VBF_GenPart_dR/F");



    Events_basic->Branch("Higgs_bb_Matched_AK4_Mjj", &Higgs_bb_Matched_AK4_Mjj, "Higgs_bb_Matched_AK4_Mjj/F");
    Events_basic->Branch("Higgs_bb_Matched_AK4_dEta", &Higgs_bb_Matched_AK4_dEta, "Higgs_bb_Matched_AK4_dEta/F");
    Events_basic->Branch("Higgs_bb_Matched_AK4_dR", &Higgs_bb_Matched_AK4_dR, "Higgs_bb_Matched_AK4_dR/F");

    Events_basic->Branch("VBF_Matched_AK4_Mjj", &VBF_Matched_AK4_Mjj, "VBF_Matched_AK4_Mjj/F");
    Events_basic->Branch("VBF_Matched_AK4_dEta", &VBF_Matched_AK4_dEta, "VBF_Matched_AK4_dEta/F");
    Events_basic->Branch("VBF_Matched_AK4_dR", &VBF_Matched_AK4_dR, "VBF_Matched_AK4_dR/F");

}

void EventLoop(TString infileName){

    bool debug = false;

    cout<<cur_time()<<"\tProcessing Tree...\n";
    cout<<cur_time()<<"\tFile Opened...\n";
    TFile* infile = TFile::Open(infileName);
    TTree* Events = (TTree*)infile->Get("Events");

    TFile *outfile = new TFile("TTree_tmp.root", "RECREATE");
    outfile->cd();

    Events_basic = Events->CloneTree(0);
    AddBranches();

    TTree* Events_Passed = Events_basic->CloneTree(0);
    Events_Passed->SetObject("Events_Passed","Events_Passed");


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

    TTreeReaderValue<UInt_t> nElectron(myReader, "nElectron");
    TTreeReaderArray<Float_t> Electron_deltaEtaSC(myReader, "Electron_deltaEtaSC");
    TTreeReaderArray<Float_t> Electron_dr03EcalRecHitSumEt(myReader, "Electron_dr03EcalRecHitSumEt");
    TTreeReaderArray<Float_t> Electron_dr03HcalDepth1TowerSumEt(myReader, "Electron_dr03HcalDepth1TowerSumEt");
    TTreeReaderArray<Float_t> Electron_dr03TkSumPt(myReader, "Electron_dr03TkSumPt");
    TTreeReaderArray<Float_t> Electron_dr03TkSumPtHEEP(myReader, "Electron_dr03TkSumPtHEEP");
    TTreeReaderArray<Float_t> Electron_dxy(myReader, "Electron_dxy");
    TTreeReaderArray<Float_t> Electron_dxyErr(myReader, "Electron_dxyErr");
    TTreeReaderArray<Float_t> Electron_dz(myReader, "Electron_dz");
    TTreeReaderArray<Float_t> Electron_dzErr(myReader, "Electron_dzErr");
    TTreeReaderArray<Float_t> Electron_eInvMinusPInv(myReader, "Electron_eInvMinusPInv");
    TTreeReaderArray<Float_t> Electron_energyErr(myReader, "Electron_energyErr");
    TTreeReaderArray<Float_t> Electron_eta(myReader, "Electron_eta");
    TTreeReaderArray<Float_t> Electron_hoe(myReader, "Electron_hoe");
    TTreeReaderArray<Float_t> Electron_ip3d(myReader, "Electron_ip3d");
    TTreeReaderArray<Float_t> Electron_jetPtRelv2(myReader, "Electron_jetPtRelv2");
    TTreeReaderArray<Float_t> Electron_jetRelIso(myReader, "Electron_jetRelIso");
    TTreeReaderArray<Float_t> Electron_mass(myReader, "Electron_mass");
    TTreeReaderArray<Float_t> Electron_phi(myReader, "Electron_phi");
    TTreeReaderArray<Float_t> Electron_pt(myReader, "Electron_pt");
    TTreeReaderArray<Float_t> Electron_r9(myReader, "Electron_r9");
    TTreeReaderArray<Float_t> Electron_sieie(myReader, "Electron_sieie");
    TTreeReaderArray<Float_t> Electron_sip3d(myReader, "Electron_sip3d");
    TTreeReaderArray<Float_t> Electron_mvaTTH(myReader, "Electron_mvaTTH");
    TTreeReaderArray<Int_t> Electron_charge(myReader, "Electron_charge");
    TTreeReaderArray<Int_t> Electron_cutBased(myReader, "Electron_cutBased");
    TTreeReaderArray<Int_t> Electron_cutBased_Fall17_V1(myReader, "Electron_cutBased_Fall17_V1");
    TTreeReaderArray<Int_t> Electron_jetIdx(myReader, "Electron_jetIdx");
    TTreeReaderArray<Int_t> Electron_pdgId(myReader, "Electron_pdgId");
    TTreeReaderArray<Int_t> Electron_photonIdx(myReader, "Electron_photonIdx");
    TTreeReaderArray<Int_t> Electron_tightCharge(myReader, "Electron_tightCharge");
    TTreeReaderArray<Int_t> Electron_vidNestedWPBitmap(myReader, "Electron_vidNestedWPBitmap");
    TTreeReaderArray<Bool_t> Electron_convVeto(myReader, "Electron_convVeto");
    TTreeReaderArray<Bool_t> Electron_cutBased_HEEP(myReader, "Electron_cutBased_HEEP");
    TTreeReaderArray<Bool_t> Electron_isPFcand(myReader, "Electron_isPFcand");


    TTreeReaderValue<UInt_t> nFatJet(myReader, "nFatJet");
    TTreeReaderArray<Float_t> FatJet_area(myReader, "FatJet_area");
    TTreeReaderArray<Float_t> FatJet_eta(myReader, "FatJet_eta");
    TTreeReaderArray<Float_t> FatJet_mass(myReader, "FatJet_mass");
    TTreeReaderArray<Float_t> FatJet_msoftdrop(myReader, "FatJet_msoftdrop");
    TTreeReaderArray<Float_t> FatJet_phi(myReader, "FatJet_phi");
    TTreeReaderArray<Float_t> FatJet_pt(myReader, "FatJet_pt");
    TTreeReaderArray<Int_t> FatJet_jetId(myReader, "FatJet_jetId");
    TTreeReaderArray<Int_t> FatJet_subJetIdx1(myReader, "FatJet_subJetIdx1");
    TTreeReaderArray<Int_t> FatJet_subJetIdx2(myReader, "FatJet_subJetIdx2");

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

    TTreeReaderValue<UInt_t> nJet(myReader, "nJet");
    TTreeReaderArray<Float_t> Jet_area(myReader, "Jet_area");
    TTreeReaderArray<Float_t> Jet_eta(myReader, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_mass(myReader, "Jet_mass");
    TTreeReaderArray<Float_t> Jet_phi(myReader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_pt(myReader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_qgl(myReader, "Jet_qgl");
    TTreeReaderArray<Float_t> Jet_rawFactor(myReader, "Jet_rawFactor");
    TTreeReaderArray<Float_t> Jet_bRegCorr(myReader, "Jet_bRegCorr");
    TTreeReaderArray<Float_t> Jet_bRegRes(myReader, "Jet_bRegRes");
    TTreeReaderArray<Int_t> Jet_electronIdx1(myReader, "Jet_electronIdx1");
    TTreeReaderArray<Int_t> Jet_electronIdx2(myReader, "Jet_electronIdx2");
    TTreeReaderArray<Int_t> Jet_jetId(myReader, "Jet_jetId");
    TTreeReaderArray<Int_t> Jet_puId(myReader, "Jet_puId");

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

    TTreeReaderValue<UInt_t> nMuon(myReader, "nMuon");
    TTreeReaderArray<Float_t> Muon_eta(myReader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_mass(myReader, "Muon_mass");
    TTreeReaderArray<Float_t> Muon_phi(myReader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_pt(myReader, "Muon_pt");
    TTreeReaderArray<Int_t> Muon_charge(myReader, "Muon_charge");
    TTreeReaderArray<Int_t> Muon_jetIdx(myReader, "Muon_jetIdx");
    TTreeReaderArray<Int_t> Muon_pdgId(myReader, "Muon_pdgId");

    TTreeReaderValue<UInt_t> nTau(myReader, "nTau");
    TTreeReaderArray<Float_t> Tau_eta(myReader, "Tau_eta");
    TTreeReaderArray<Float_t> Tau_mass(myReader, "Tau_mass");
    TTreeReaderArray<Float_t> Tau_phi(myReader, "Tau_phi");
    TTreeReaderArray<Float_t> Tau_pt(myReader, "Tau_pt");

    TTreeReaderArray<Int_t> Electron_genPartIdx(myReader, "Electron_genPartIdx");
    TTreeReaderArray<UChar_t> Electron_genPartFlav(myReader, "Electron_genPartFlav");
    TTreeReaderArray<Int_t> Muon_genPartIdx(myReader, "Muon_genPartIdx");
    TTreeReaderArray<UChar_t> Muon_genPartFlav(myReader, "Muon_genPartFlav");
    TTreeReaderArray<Int_t> Tau_genPartIdx(myReader, "Tau_genPartIdx");
    TTreeReaderArray<UChar_t> Tau_genPartFlav(myReader, "Tau_genPartFlav");

    TTreeReaderArray<Float_t> FatJet_particleNet_HbbvsQCD(myReader, "FatJet_particleNet_HbbvsQCD");
    TTreeReaderArray<Float_t> FatJet_particleNet_WvsQCD(myReader, "FatJet_particleNet_WvsQCD");

    //int eventLoopMax = EvMax;
    float jetPT_cut = 50.0;
    float jetInvMass_cut = 400.0;
    float jetEtaSep_cut = 2.0;

    Long64_t eventLoopMax = myReader.GetEntries();
    int iev = 0;
    int iLHEPassed = 0;
    float eventScale = 1.0;

    vector<float> nLeptons_vec = {0.0,0.0,0.0};

    //BEGIN EVENT LOOP
    while (myReader.Next()) {
        //if(iev == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;
        //if(iev == 1) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;
        if(++iev % 10000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;
        
        //if(iev > 50000) break;
        //if(iev>5) break;


        vector<int> vec_GenPart_pdgId;
        vector<int> vec_GenPart_genPartIdxMother;
        vector<int> vec_GenPart_status;


        /*  -----------------  SETUP  -----------------  */
        GenMatching_Initial_FullMatch = false;
        GenMatching_Decayed_FullMatch = false;
        POMatching_FullMatch = false;


        Wm_GenPart_decay_type = -1; //1 = Leptonically, 0 = Hadronically
        Wp_GenPart_decay_type = -1; //1 = Leptonically, 0 = Hadronically
        Higgs_GenPart_decay_type = -1;

        Higgs_GenPart_idx = -1;//2;
        Wm_GenPart_idx = -1;//4;
        Wp_GenPart_idx = -1;//3;
        VBFLead_GenPart_idx = -1;//5;
        VBFTrail_GenPart_idx = -1;//6;
        Z_GenPart_idx = -1;

        Higgs_GenPart_b_idx = -1; 
        Higgs_GenPart_bbar_idx = -1;

        Z_GenPart_daughterLead_idx = -1;
        Z_GenPart_daughterTrail_idx = -1;
        Z_GenPart_daughterLead_pdgid = -999;
        Z_GenPart_daughterTrail_pdgid = -999;

        Wp_GenPart_daughterLead_idx = -1;
        Wp_GenPart_daughterTrail_idx = -1;
        Wp_GenPart_daughterLead_pdgid = -999;
        Wp_GenPart_daughterTrail_pdgid = -999;

        Wm_GenPart_daughterLead_idx = -1;
        Wm_GenPart_daughterTrail_idx = -1;
        Wm_GenPart_daughterLead_pdgid = -999;
        Wm_GenPart_daughterTrail_pdgid = -999;

        Higgs_Matched_AK8_idx = -1;
        Higgs_b_Matched_AK4_idx = -1;
        Higgs_bbar_Matched_AK4_idx = -1;

        VBFLead_Matched_AK4_idx = -1;
        VBFTrail_Matched_AK4_idx = -1;


        Wp_Matched_AK8_idx = -1;
        Wp_Matched_Electron_idx = -1;
        Wp_Matched_Muon_idx = -1;

        Wm_Matched_AK8_idx = -1;
        Wm_Matched_Electron_idx = -1;
        Wm_Matched_Muon_idx = -1;

        Z_Matched_AK8_idx = -1;
        Z_Matched_Electron_1_idx = -1;
        Z_Matched_Electron_2_idx = -1;
        Z_Matched_Muon_1_idx = -1;
        Z_Matched_Muon_2_idx = -1;


        VBF_LHEPart_Mjj = -1;
        VBF_LHEPart_dEta = -1;
        VBF_LHEPart_dR = -1;

        Higgs_GenPart_bb_Mjj = -1;
        Higgs_GenPart_bb_dEta = -1;
        Higgs_GenPart_bb_dR = -1;

        VBF_GenPart_Mjj = -1;
        VBF_GenPart_dEta = -1;
        VBF_GenPart_dR = -1;

        Higgs_bb_Matched_AK4_Mjj = -1;
        Higgs_bb_Matched_AK4_dEta = -1;
        Higgs_bb_Matched_AK4_dR = -1;

        VBF_Matched_AK4_Mjj = -1;
        VBF_Matched_AK4_dEta = -1;
        VBF_Matched_AK4_dR = -1;

        //cout<<"AAAA"<<endl;

        for(int iJet=0; iJet<*nGenPart; iJet++){
            vec_GenPart_pdgId.push_back(GenPart_pdgId[iJet]);
            vec_GenPart_genPartIdxMother.push_back(GenPart_genPartIdxMother[iJet]);
            vec_GenPart_status.push_back(GenPart_status[iJet]);

            if(debug){
                if(GenPart_pdgId[iJet] == -24) cout<<"Found Wm:\tIdx="<<iJet<<"\tStatus="<<GenPart_status[iJet]<<endl;
                if(GenPart_pdgId[iJet] == 24) cout<<"Found Wp:\tIdx="<<iJet<<"\tStatus="<<GenPart_status[iJet]<<endl;

                int motherIdx = GenPart_genPartIdxMother[iJet];
                if(motherIdx != iJet){
                    if(GenPart_pdgId[motherIdx] == -24) cout<<"Found Wm daughter:\tIdx="<<iJet<<"\tStatus="<<GenPart_status[iJet]<<endl;
                    if(GenPart_pdgId[motherIdx] == 24) cout<<"Found Wp daughter:\tIdx="<<iJet<<"\tStatus="<<GenPart_status[iJet]<<endl;
                }

                //if(GenPart_status[iJet] == 1) cout<<"Found Final State Particle:\tIdx="<<iJet<<"\tStatus="<<GenPart_status[iJet]<<"\tpdgid="<<GenPart_pdgId[iJet]<<"\tMother pdgid="<<GenPart_pdgId[motherIdx]<<endl;

            }
        }


        /*  -----------------  LHE  -----------------  */



        /*  -----------------  GEN  -----------------  */
        
        //cout<<"CCC"<<endl;

        for(int iJet=0; iJet<8; iJet++){
            //cout<<"\tgenpStatus="<<GenPart_status[iJet]<<"\tgenpPdgid="<<GenPart_pdgId[iJet]<<endl;
            if(GenPart_genPartIdxMother[iJet] != 0) continue;
            if(GenPart_pdgId[iJet] == 25) Higgs_GenPart_idx = iJet;
            if(GenPart_pdgId[iJet] == 23) Z_GenPart_idx = iJet;
            if(GenPart_pdgId[iJet] == 24) Wp_GenPart_idx = iJet;
            if(GenPart_pdgId[iJet] == -24) Wm_GenPart_idx = iJet;
            if(fabs(GenPart_pdgId[iJet]) <7 && VBFLead_GenPart_idx == -1) VBFLead_GenPart_idx = iJet;
            if(fabs(GenPart_pdgId[iJet]) <7 && VBFLead_GenPart_idx != -1 && VBFLead_GenPart_idx != iJet) VBFTrail_GenPart_idx = iJet;
        }
        //cout<<"\tVBF1="<<VBFLead_GenPart_idx<<"\tVBF2="<<VBFTrail_GenPart_idx<<"\tH="<<Higgs_GenPart_idx 
        //<<"\tZ="<<Z_GenPart_idx<<"\tWp="<<Wp_GenPart_idx<<"\tWm="<<Wm_GenPart_idx<<endl;

        //get decay products
        vector<std::tuple<int,int>> vec_Higgs_GenPart_DPidxs;
        vector<std::tuple<int,int>> vec_Z_GenPart_DPidxs;
        vector<std::tuple<int,int>> vec_Wp_GenPart_DPidxs;
        vector<std::tuple<int,int>> vec_Wm_GenPart_DPidxs;

        /*cout<<"EVENT #"<<iev<<endl;
        cout<<"\tH: "<<Higgs_GenPart_idx<<endl;
        cout<<"\tgetDecayProducts: ";*/

        if(Higgs_GenPart_idx != -1) vec_Higgs_GenPart_DPidxs = getDecayProducts(Higgs_GenPart_idx, vec_GenPart_pdgId, vec_GenPart_genPartIdxMother);
        if(Z_GenPart_idx != -1) vec_Z_GenPart_DPidxs = getDecayProducts(Z_GenPart_idx, vec_GenPart_pdgId, vec_GenPart_genPartIdxMother);
        if(Wp_GenPart_idx != -1) vec_Wp_GenPart_DPidxs = getDecayProducts(Wp_GenPart_idx, vec_GenPart_pdgId, vec_GenPart_genPartIdxMother);
        if(Wm_GenPart_idx != -1) vec_Wm_GenPart_DPidxs = getDecayProducts(Wm_GenPart_idx, vec_GenPart_pdgId, vec_GenPart_genPartIdxMother);

        /*
        for(int i=0; i< vec_Higgs_GenPart_DPidxs.size(); i++)
            cout<<get<0>(vec_Higgs_GenPart_DPidxs[i])<<", ";
        cout<<endl;*/
        //check decay modes (check if all decay products are quarks)
        if(!vec_Higgs_GenPart_DPidxs.empty() && std::all_of( std::begin(vec_Higgs_GenPart_DPidxs), std::end(vec_Higgs_GenPart_DPidxs), []( std::tuple<int,int> x){ return fabs(get<1>(x))<9; } )) 
            Higgs_GenPart_decay_type = 0;
        if(!vec_Z_GenPart_DPidxs.empty() && std::all_of( std::begin(vec_Z_GenPart_DPidxs), std::end(vec_Z_GenPart_DPidxs), []( std::tuple<int,int> x){ return fabs(get<1>(x))<9; } )) 
            Z_GenPart_decay_type = 0;
        if(!vec_Wp_GenPart_DPidxs.empty() && std::all_of( std::begin(vec_Wp_GenPart_DPidxs), std::end(vec_Wp_GenPart_DPidxs), []( std::tuple<int,int> x){ return fabs(get<1>(x))<9; } )) 
            Wp_GenPart_decay_type = 0;
        if(!vec_Wm_GenPart_DPidxs.empty() && std::all_of( std::begin(vec_Wm_GenPart_DPidxs), std::end(vec_Wm_GenPart_DPidxs), []( std::tuple<int,int> x){ return fabs(get<1>(x))<9; } )) 
            Wm_GenPart_decay_type = 0;

        /*cout<<"H_dau: "<<vec_Higgs_GenPart_DPidxs.size()<<"\tZ_dau: "<<vec_Z_GenPart_DPidxs.size()
        <<"\tWp_dau: "<<vec_Wp_GenPart_DPidxs.size()
        <<"\tWm_dau: "<<vec_Wm_GenPart_DPidxs.size()<<endl;*/


        //cout<<"BBB"<<endl;

        /*for(int iGenPart=0; iGenPart<*nGenPart; iGenPart++){
            cout<<"GenPart Idx="<<iGenPart<<"\tpdgid="<<GenPart_pdgId[iGenPart]<<"\tstatus="<<GenPart_status[iGenPart]<<"\tmotherIdx="<<GenPart_genPartIdxMother[iGenPart]<<"\teta="<<GenPart_eta[iGenPart]<<"\tphi="<<GenPart_phi[iGenPart]<<"\tpt="<<GenPart_pt[iGenPart]<<endl;

        }*/

        //get final state particles (when not quarks)
        vector<int> vec_Higgs_GenPart_FSidxs;
        vector<int> vec_Z_GenPart_FSidxs;
        vector<int> vec_Wp_GenPart_FSidxs;
        vector<int> vec_Wm_GenPart_FSidxs;

        if(Higgs_GenPart_idx != -1 && Higgs_GenPart_decay_type != 0){
            for(int ii=0; ii<vec_Higgs_GenPart_DPidxs.size(); ii++){
               int fs_tmp = getFinalState(get<0>(vec_Higgs_GenPart_DPidxs[ii]), vec_GenPart_pdgId, vec_GenPart_genPartIdxMother, vec_GenPart_status);
               if(fs_tmp != -1) vec_Higgs_GenPart_FSidxs.push_back(fs_tmp);
            }
        }

        /*cout<<"\tgetFS: ";
        for(int i=0; i< vec_Higgs_GenPart_FSidxs.size(); i++)
            cout<<vec_Higgs_GenPart_FSidxs[i]<<", ";
        cout<<endl;*/
        if(Z_GenPart_idx != -1 && Z_GenPart_decay_type != 0){
            for(int ii=0; ii<vec_Z_GenPart_DPidxs.size(); ii++){
               int fs_tmp = getFinalState(get<0>(vec_Z_GenPart_DPidxs[ii]), vec_GenPart_pdgId, vec_GenPart_genPartIdxMother, vec_GenPart_status);
               if(fs_tmp != -1) vec_Z_GenPart_FSidxs.push_back(fs_tmp);
            }
        }
        if(Wp_GenPart_idx != -1 && Wp_GenPart_decay_type != 0){
            for(int ii=0; ii<vec_Wp_GenPart_DPidxs.size(); ii++){
               int fs_tmp = getFinalState(get<0>(vec_Wp_GenPart_DPidxs[ii]), vec_GenPart_pdgId, vec_GenPart_genPartIdxMother, vec_GenPart_status);
               if(fs_tmp != -1) vec_Wp_GenPart_FSidxs.push_back(fs_tmp);
            }
        }
        if(Wm_GenPart_idx != -1 && Wm_GenPart_decay_type != 0){
            for(int ii=0; ii<vec_Wm_GenPart_DPidxs.size(); ii++){
               int fs_tmp = getFinalState(get<0>(vec_Wm_GenPart_DPidxs[ii]), vec_GenPart_pdgId, vec_GenPart_genPartIdxMother, vec_GenPart_status);
               if(fs_tmp != -1) vec_Wm_GenPart_FSidxs.push_back(fs_tmp);
            }
        }
        //cout<<"BBBB"<<endl;
        //cout<<"XXX"<<endl;

        //Higgs b quarks 
        //find GenPart b quarks - pdgid = +-5

        int Higgs_GenPart_daughter_1 = -1,  Higgs_GenPart_daughter_2 = -1;

        if(Higgs_GenPart_idx != -1){
            if(Higgs_GenPart_decay_type == 0 && vec_Higgs_GenPart_DPidxs.size()>1){
                if(std::all_of( std::begin(vec_Higgs_GenPart_DPidxs), std::end(vec_Higgs_GenPart_DPidxs), []( std::tuple<int,int> x){ return fabs(get<1>(x))==5; } )){

                    Higgs_DecayTobbbar = true;
                    Higgs_GenPart_b_idx = get<0>(vec_Higgs_GenPart_DPidxs.at(0));
                    Higgs_GenPart_bbar_idx = get<0>(vec_Higgs_GenPart_DPidxs.at(1));

                    if(GenPart_pdgId[Higgs_GenPart_b_idx] == -5) swap(Higgs_GenPart_b_idx, Higgs_GenPart_bbar_idx);


                    Higgs_GenPart_bb_dR = getDR(GenPart_phi[Higgs_GenPart_b_idx], GenPart_eta[Higgs_GenPart_b_idx], GenPart_phi[Higgs_GenPart_bbar_idx], GenPart_eta[Higgs_GenPart_bbar_idx]);
                    Higgs_GenPart_bb_dEta = fabs(GenPart_eta[Higgs_GenPart_b_idx] - GenPart_eta[Higgs_GenPart_bbar_idx]);
                    
                    ROOT::Math::PtEtaPhiMVector genHb_p4(GenPart_pt[Higgs_GenPart_b_idx], GenPart_eta[Higgs_GenPart_b_idx], GenPart_phi[Higgs_GenPart_b_idx], GenPart_mass[Higgs_GenPart_b_idx]);
                    ROOT::Math::PtEtaPhiMVector genHbbar_p4(GenPart_pt[Higgs_GenPart_bbar_idx], GenPart_eta[Higgs_GenPart_bbar_idx], GenPart_phi[Higgs_GenPart_bbar_idx], GenPart_mass[Higgs_GenPart_bbar_idx]);
                    ROOT::Math::PtEtaPhiMVector genHbb_jetSum = genHb_p4+genHbbar_p4;
                    Higgs_GenPart_bb_Mjj = sqrt(genHbb_jetSum.Dot(genHbb_jetSum));

                }
                Higgs_GenPart_daughter_1 = get<0>(vec_Higgs_GenPart_DPidxs.at(0));
                Higgs_GenPart_daughter_2 = get<0>(vec_Higgs_GenPart_DPidxs.at(1));
            }
            else if(vec_Higgs_GenPart_FSidxs.size() > 1){
                Higgs_GenPart_daughter_1 = vec_Higgs_GenPart_FSidxs.at(0);
                Higgs_GenPart_daughter_2 = vec_Higgs_GenPart_FSidxs.at(1);
                Higgs_GenPart_decay_type = getDecayType(GenPart_pdgId[Higgs_GenPart_daughter_1], GenPart_pdgId[Higgs_GenPart_daughter_2]);
            }
        }

            /*cout<<"\tD1: "<<Higgs_GenPart_daughter_1<<"\tD2: "<<Higgs_GenPart_daughter_2<<"\tB: "<<Higgs_GenPart_b_idx<<"\tBBar: "<<Higgs_GenPart_bbar_idx<<endl;
            cout<<"\tDecayType: "<<Higgs_GenPart_decay_type<<"\tDecaytobbar: "<<Higgs_DecayTobbbar<<endl;
            for(int iGenPart=0; iGenPart<*nGenPart; iGenPart++){
                cout<<"\tGenPart Idx="<<iGenPart<<"\tpdgid="<<GenPart_pdgId[iGenPart]<<"\tstatus="<<GenPart_status[iGenPart]<<"\tmotherIdx="<<GenPart_genPartIdxMother[iGenPart]<<"\teta="<<GenPart_eta[iGenPart]<<"\tphi="<<GenPart_phi[iGenPart]<<"\tpt="<<GenPart_pt[iGenPart]<<endl;

            }*/

        //cout<<"YYYYYY"<<endl;

        if(Z_GenPart_idx != -1){
            if(Z_GenPart_decay_type == 0 && vec_Z_GenPart_DPidxs.size()>1){
                Z_GenPart_daughterLead_idx = get<0>(vec_Z_GenPart_DPidxs.at(0));
                Z_GenPart_daughterTrail_idx = get<0>(vec_Z_GenPart_DPidxs.at(1));
                Z_GenPart_daughterLead_pdgid = GenPart_pdgId[Z_GenPart_daughterLead_idx];
                Z_GenPart_daughterTrail_pdgid = GenPart_pdgId[Z_GenPart_daughterTrail_idx];
            }
            else if(vec_Z_GenPart_FSidxs.size()>1){
                Z_GenPart_daughterLead_idx = vec_Z_GenPart_FSidxs.at(0);
                Z_GenPart_daughterTrail_idx = vec_Z_GenPart_FSidxs.at(1);
                Z_GenPart_daughterLead_pdgid = GenPart_pdgId[Z_GenPart_daughterLead_idx];
                Z_GenPart_daughterTrail_pdgid = GenPart_pdgId[Z_GenPart_daughterTrail_idx];
                Z_GenPart_decay_type = getDecayType(GenPart_pdgId[Z_GenPart_daughterLead_idx], GenPart_pdgId[Z_GenPart_daughterTrail_idx]);
            }
        }

        if(Wp_GenPart_idx != -1){
            if(Wp_GenPart_decay_type == 0 && vec_Wp_GenPart_DPidxs.size()>1){
                Wp_GenPart_daughterLead_idx = get<0>(vec_Wp_GenPart_DPidxs.at(0));
                Wp_GenPart_daughterTrail_idx = get<0>(vec_Wp_GenPart_DPidxs.at(1));
                Wp_GenPart_daughterLead_pdgid = GenPart_pdgId[Wp_GenPart_daughterLead_idx];
                Wp_GenPart_daughterTrail_pdgid = GenPart_pdgId[Wp_GenPart_daughterTrail_idx];
            }
            else if(vec_Wp_GenPart_FSidxs.size()>1){
                Wp_GenPart_daughterLead_idx = vec_Wp_GenPart_FSidxs.at(0);
                Wp_GenPart_daughterTrail_idx = vec_Wp_GenPart_FSidxs.at(1);
                Wp_GenPart_daughterLead_pdgid = GenPart_pdgId[Wp_GenPart_daughterLead_idx];
                Wp_GenPart_daughterTrail_pdgid = GenPart_pdgId[Wp_GenPart_daughterTrail_idx];
                Wp_GenPart_decay_type = getDecayType(GenPart_pdgId[Wp_GenPart_daughterLead_idx], GenPart_pdgId[Wp_GenPart_daughterTrail_idx]);
            }
        }

        if(Wm_GenPart_idx != -1){
            if(Wm_GenPart_decay_type == 0 && vec_Wm_GenPart_DPidxs.size()>1){
                Wm_GenPart_daughterLead_idx = get<0>(vec_Wm_GenPart_DPidxs.at(0));
                Wm_GenPart_daughterTrail_idx = get<0>(vec_Wm_GenPart_DPidxs.at(1));
                Wm_GenPart_daughterLead_pdgid = GenPart_pdgId[Wm_GenPart_daughterLead_idx];
                Wm_GenPart_daughterTrail_pdgid = GenPart_pdgId[Wm_GenPart_daughterTrail_idx];
            }
            else if(vec_Wm_GenPart_FSidxs.size()>1){
                Wm_GenPart_daughterLead_idx = vec_Wm_GenPart_FSidxs.at(0);
                Wm_GenPart_daughterTrail_idx = vec_Wm_GenPart_FSidxs.at(1);
                Wm_GenPart_daughterLead_pdgid = GenPart_pdgId[Wm_GenPart_daughterLead_idx];
                Wm_GenPart_daughterTrail_pdgid = GenPart_pdgId[Wm_GenPart_daughterTrail_idx];
                Wm_GenPart_decay_type = getDecayType(GenPart_pdgId[Wm_GenPart_daughterLead_idx], GenPart_pdgId[Wm_GenPart_daughterTrail_idx]);
            }
        }


        //cout<<"ZZZ"<<endl;
        //cout<<"CCCC"<<endl;
        //Figure out how to tell if gen matching is complete and went well
        //if not, skip object matching
        //also skip object matching if Higgs does not decay to bb

        int nMatchedV = ((Z_GenPart_idx > -1) ? 1:0) + ((Wp_GenPart_idx > -1) ? 1:0) + ((Wm_GenPart_idx > -1) ? 1:0);
        int nMatchedVDecays = ((Z_GenPart_daughterLead_idx > -1 && Z_GenPart_daughterTrail_idx > -1) ? 1:0) 
                            + ((Wp_GenPart_daughterLead_idx > -1 && Wp_GenPart_daughterTrail_idx > -1) ? 1:0) 
                            + ((Wm_GenPart_daughterLead_idx > -1 && Wm_GenPart_daughterTrail_idx > -1) ? 1:0);

        if(nMatchedV == 2 && VBFLead_GenPart_idx > -1 && VBFTrail_GenPart_idx > -1 && Higgs_GenPart_idx > -1){
            GenMatching_Initial_FullMatch = true;

            if(nMatchedVDecays == 2) GenMatching_Decayed_FullMatch = true;
        }

        //if(nMatchedV < 2){cout<<"Less than 2 matched V"<<endl;
        //cout<< "\tZ: "<<Z_GenPart_idx<< " Wp: "<<Wp_GenPart_idx<< " Wm: "<<Wm_GenPart_idx<<endl;}
        /*if(nMatchedV == 3){
            
            cout<<"EVENT #"<<iev<<endl;
            cout<<"\t3 matched V"<<endl;
            cout<< "\tZ: "<<Z_GenPart_idx<< " Wp: "<<Wp_GenPart_idx<< " Wm: "<<Wm_GenPart_idx<<endl;
            for(int iGenPart=0; iGenPart<15; iGenPart++){
                cout<<"\tGenPart Idx="<<iGenPart<<"\tpdgid="<<GenPart_pdgId[iGenPart]<<"\tstatus="<<GenPart_status[iGenPart]<<"\tmotherIdx="<<GenPart_genPartIdxMother[iGenPart]<<"\teta="<<GenPart_eta[iGenPart]<<"\tphi="<<GenPart_phi[iGenPart]<<"\tpt="<<GenPart_pt[iGenPart]<<endl;

            }
        }*/

        //cout<<"HHH"<<endl;

        if(!Higgs_DecayTobbbar)nPassGenMatchedDecay_FailHbb++;
        if(Z_GenPart_decay_type > 0)nPassGenMatchedDecay_FailZ++;
        if(!Higgs_DecayTobbbar || Z_GenPart_decay_type > 0){nPassGenMatchedDecay_FailHbb_FailZ++; continue;}

        if(!GenMatching_Decayed_FullMatch){
            nFailedGenMatchedDecay++;
             continue;
        }

        //cout<<"GGG"<<endl;
        VBF_GenPart_dR = getDR(GenPart_phi[VBFLead_GenPart_idx], GenPart_eta[VBFLead_GenPart_idx], GenPart_phi[VBFTrail_GenPart_idx], GenPart_eta[VBFTrail_GenPart_idx]);
        VBF_GenPart_dEta = fabs(GenPart_eta[VBFLead_GenPart_idx] - GenPart_eta[VBFTrail_GenPart_idx]);
        
        ROOT::Math::PtEtaPhiMVector genvbfL_p4(GenPart_pt[VBFLead_GenPart_idx], GenPart_eta[VBFLead_GenPart_idx], GenPart_phi[VBFLead_GenPart_idx], GenPart_mass[VBFLead_GenPart_idx]);
        ROOT::Math::PtEtaPhiMVector genvbfT_p4(GenPart_pt[VBFTrail_GenPart_idx], GenPart_eta[VBFTrail_GenPart_idx], GenPart_phi[VBFTrail_GenPart_idx], GenPart_mass[VBFTrail_GenPart_idx]);
        ROOT::Math::PtEtaPhiMVector genvbf_jetSum = genvbfL_p4+genvbfT_p4;
        VBF_GenPart_Mjj = sqrt(genvbf_jetSum.Dot(genvbf_jetSum));

        //if(Z_GenPart_decay_type > 0){
           //nPassGenMatchedDecay_FailHbb_FailZ++;

            //Events->GetEntry(myReader.GetCurrentEntry());
            //Events_basic -> Fill();
            //continue;
        //}

        /*  -----------------  MATCHED RECONSTRUCTED  -----------------  */
        nPassGenMatchedDecay_Total++;

        //cout<<"HELLO"<<endl;


        //cout<<"AAA"<<endl;
        //match vbf jets
        vector< std::tuple< int, float > > matchedJet_candidates_1, matchedJet_candidates_2;
        float dR_cut = 0.4;
        for(int iJet=0; iJet<*nJet; iJet++){

            float dR_1 = getDR(GenPart_phi[VBFLead_GenPart_idx], GenPart_eta[VBFLead_GenPart_idx], Jet_phi[iJet], Jet_eta[iJet]);
            float dR_2 = getDR(GenPart_phi[VBFTrail_GenPart_idx], GenPart_eta[VBFTrail_GenPart_idx], Jet_phi[iJet], Jet_eta[iJet]);

            if(dR_1 < dR_cut) matchedJet_candidates_1.push_back(make_tuple(iJet,Jet_pt[iJet]));
            if(dR_2 < dR_cut) matchedJet_candidates_2.push_back(make_tuple(iJet,Jet_pt[iJet]));
        }
        sort(matchedJet_candidates_1.begin(), matchedJet_candidates_1.end(), TupleCompare<1>());
        if(matchedJet_candidates_1.size() > 0) VBFLead_Matched_AK4_idx = get<0>(matchedJet_candidates_1.at(0));

        sort(matchedJet_candidates_2.begin(), matchedJet_candidates_2.end(), TupleCompare<1>());
        if(matchedJet_candidates_2.size() > 0) VBFTrail_Matched_AK4_idx = get<0>(matchedJet_candidates_2.at(0));
        
        if(VBFLead_Matched_AK4_idx > -1 && VBFTrail_Matched_AK4_idx > -1)
            if(Jet_pt[VBFTrail_Matched_AK4_idx] > Jet_pt[VBFLead_Matched_AK4_idx]) 
                swap(VBFLead_Matched_AK4_idx, VBFTrail_Matched_AK4_idx);




        if(Z_GenPart_idx != -1){
            //Only care when Z decays to qq
            if(Z_GenPart_decay_type == 0){
                vector< std::tuple< int, float > > matchedFJet_cands;
                for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
                    float dR = getDR(GenPart_phi[Z_GenPart_idx], GenPart_eta[Z_GenPart_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                    if(dR < dR_cut) matchedFJet_cands.push_back(std::make_tuple(iFatJet, FatJet_pt[iFatJet]));
                }
                sort(matchedFJet_cands.begin(), matchedFJet_cands.end(), TupleCompare<1>());
                if(matchedFJet_cands.size() > 0) Z_Matched_AK8_idx = get<0>(matchedFJet_cands.at(0));
            }
        }
        if(Wp_GenPart_idx != -1){
            if(Wp_GenPart_decay_type == 0){ //FatJet
                vector< std::tuple< int, float > > matchedFJet_cands;
                for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
                    float dR = getDR(GenPart_phi[Wp_GenPart_idx], GenPart_eta[Wp_GenPart_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                    if(dR < dR_cut) matchedFJet_cands.push_back(std::make_tuple(iFatJet, FatJet_pt[iFatJet]));
                }
                sort(matchedFJet_cands.begin(), matchedFJet_cands.end(), TupleCompare<1>());
                if(matchedFJet_cands.size() > 0) Wp_Matched_AK8_idx = get<0>(matchedFJet_cands.at(0));
            }
            if(Wp_GenPart_decay_type == 1){ //Electron
                for(int iEle=0; iEle<*nElectron; iEle++){
                    int genIdx = Electron_genPartIdx[iEle];
                    if(genIdx == -1) continue;
                    if(genIdx == Wp_GenPart_daughterLead_idx || genIdx == Wp_GenPart_daughterTrail_idx){
                        Wp_Matched_Electron_idx = genIdx;
                        break;
                    }

                }
            }
            if(Wp_GenPart_decay_type == 2){ //Muon
                for(int iMu=0; iMu<*nMuon; iMu++){
                    int genIdx = Muon_genPartIdx[iMu];
                    if(genIdx == -1) continue;
                    if(genIdx == Wp_GenPart_daughterLead_idx || genIdx == Wp_GenPart_daughterTrail_idx){
                        Wp_Matched_Muon_idx = genIdx;
                        break;
                    }
                }
            }
        }


        if(Wm_GenPart_idx != -1){
            if(Wm_GenPart_decay_type == 0){ //FatJet
                vector< std::tuple< int, float > > matchedFJet_cands;
                for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
                    float dR = getDR(GenPart_phi[Wm_GenPart_idx], GenPart_eta[Wm_GenPart_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                    if(dR < dR_cut) matchedFJet_cands.push_back(std::make_tuple(iFatJet, FatJet_pt[iFatJet]));
                }
                sort(matchedFJet_cands.begin(), matchedFJet_cands.end(), TupleCompare<1>());
                if(matchedFJet_cands.size() > 0) Wm_Matched_AK8_idx = get<0>(matchedFJet_cands.at(0));
            }
            if(Wm_GenPart_decay_type == 1){ //Electron
                for(int iEle=0; iEle<*nElectron; iEle++){
                    int genIdx = Electron_genPartIdx[iEle];
                    if(genIdx == -1) continue;
                    if(genIdx == Wm_GenPart_daughterLead_idx || genIdx == Wm_GenPart_daughterTrail_idx){
                        Wm_Matched_Electron_idx = genIdx;
                        break;
                    }

                }
            }
            if(Wm_GenPart_decay_type == 2){ //Muon
                for(int iMu=0; iMu<*nMuon; iMu++){
                    int genIdx = Muon_genPartIdx[iMu];
                    if(genIdx == -1) continue;
                    if(genIdx == Wm_GenPart_daughterLead_idx || genIdx == Wm_GenPart_daughterTrail_idx){
                        Wm_Matched_Muon_idx = genIdx;
                        break;
                    }
                }
            }
        }


        //Match Higgs
        vector< std::tuple< int, float > > Matched_H_candidates_fatJets;
        for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
            float dR_1 = getDR(GenPart_phi[Higgs_GenPart_idx], GenPart_eta[Higgs_GenPart_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
            if(dR_1 < 0.4) Matched_H_candidates_fatJets.push_back(make_tuple(iFatJet, FatJet_pt[iFatJet]));
        }
        sort(Matched_H_candidates_fatJets.begin(), Matched_H_candidates_fatJets.end(), TupleCompare<1>());
        if(Matched_H_candidates_fatJets.size() > 0) Higgs_Matched_AK8_idx = get<0>(Matched_H_candidates_fatJets.at(0));

        vector< std::tuple< int, float > > Matched_b_candidates_Jets;
        vector< std::tuple< int, float > > Matched_bbar_candidates_Jets;
        if(Higgs_GenPart_b_idx != -1 && Higgs_GenPart_bbar_idx != -1){
            for(int iJet=0; iJet<*nJet; iJet++){
                float dR_1 = getDR(GenPart_phi[Higgs_GenPart_b_idx], GenPart_eta[Higgs_GenPart_b_idx], Jet_phi[iJet], Jet_eta[iJet]);
                float dR_2 = getDR(GenPart_phi[Higgs_GenPart_bbar_idx], GenPart_eta[Higgs_GenPart_bbar_idx], Jet_phi[iJet], Jet_eta[iJet]);
                if(dR_1 < 0.4) Matched_b_candidates_Jets.push_back(make_tuple(iJet, Jet_pt[iJet]));
                if(dR_2 < 0.4) Matched_bbar_candidates_Jets.push_back(make_tuple(iJet, Jet_pt[iJet]));
            }
        }
        sort(Matched_b_candidates_Jets.begin(), Matched_b_candidates_Jets.end(), TupleCompare<1>());
        if(Matched_b_candidates_Jets.size() > 0) Higgs_b_Matched_AK4_idx = get<0>(Matched_b_candidates_Jets.at(0));
        sort(Matched_bbar_candidates_Jets.begin(), Matched_bbar_candidates_Jets.end(), TupleCompare<1>());
        if(Matched_bbar_candidates_Jets.size() > 0) Higgs_bbar_Matched_AK4_idx = get<0>(Matched_bbar_candidates_Jets.at(0));
        





        bool VBF_PO_Matched = VBFLead_Matched_AK4_idx > -1 && VBFTrail_Matched_AK4_idx > -1;
        bool Higgs_PO_Matched = (Higgs_Matched_AK8_idx > -1) || (Higgs_b_Matched_AK4_idx > -1 && Higgs_bbar_Matched_AK4_idx > -1);
        bool Z_PO_Matched = Z_Matched_AK8_idx > -1;
        bool Wp_PO_Matched = Wp_Matched_AK8_idx > -1 || Wp_Matched_Electron_idx > -1 || Wp_Matched_Muon_idx > -1;
        bool Wm_PO_Matched = Wm_Matched_AK8_idx > -1 || Wm_Matched_Electron_idx > -1 || Wm_Matched_Muon_idx > -1;
        
        int nMatchedV_PO = (Z_PO_Matched ? 1:0) + (Wp_PO_Matched ? 1:0) + (Wm_PO_Matched ? 1:0);

        if(nMatchedV_PO == 2 && VBF_PO_Matched && Higgs_PO_Matched) POMatching_FullMatch = true;

        if(VBF_PO_Matched){

            VBF_Matched_AK4_dR = getDR(Jet_phi[VBFLead_Matched_AK4_idx], Jet_eta[VBFLead_Matched_AK4_idx], Jet_phi[VBFTrail_Matched_AK4_idx], Jet_eta[VBFTrail_Matched_AK4_idx]);
            VBF_Matched_AK4_dEta = fabs(Jet_eta[VBFLead_Matched_AK4_idx] - Jet_eta[VBFTrail_Matched_AK4_idx]);
            
            ROOT::Math::PtEtaPhiMVector povbfL_p4(Jet_pt[VBFLead_Matched_AK4_idx], Jet_eta[VBFLead_Matched_AK4_idx], Jet_phi[VBFLead_Matched_AK4_idx], Jet_mass[VBFLead_Matched_AK4_idx]);
            ROOT::Math::PtEtaPhiMVector povbfT_p4(Jet_pt[VBFTrail_Matched_AK4_idx], Jet_eta[VBFTrail_Matched_AK4_idx], Jet_phi[VBFTrail_Matched_AK4_idx], Jet_mass[VBFTrail_Matched_AK4_idx]);
            ROOT::Math::PtEtaPhiMVector povbf_jetSum = povbfL_p4+povbfT_p4;
            VBF_Matched_AK4_Mjj = sqrt(povbf_jetSum.Dot(povbf_jetSum));
        }
        if(Higgs_b_Matched_AK4_idx > -1 && Higgs_bbar_Matched_AK4_idx > -1){

            Higgs_bb_Matched_AK4_dR = getDR(Jet_phi[Higgs_b_Matched_AK4_idx], Jet_eta[Higgs_b_Matched_AK4_idx], Jet_phi[Higgs_bbar_Matched_AK4_idx], Jet_eta[Higgs_bbar_Matched_AK4_idx]);
            Higgs_bb_Matched_AK4_dEta = fabs(Jet_eta[Higgs_b_Matched_AK4_idx] - Jet_eta[Higgs_bbar_Matched_AK4_idx]);
            
            ROOT::Math::PtEtaPhiMVector pohb_p4(Jet_pt[Higgs_b_Matched_AK4_idx], Jet_eta[Higgs_b_Matched_AK4_idx], Jet_phi[Higgs_b_Matched_AK4_idx], Jet_mass[Higgs_b_Matched_AK4_idx]);
            ROOT::Math::PtEtaPhiMVector pohbbar_p4(Jet_pt[Higgs_bbar_Matched_AK4_idx], Jet_eta[Higgs_bbar_Matched_AK4_idx], Jet_phi[Higgs_bbar_Matched_AK4_idx], Jet_mass[Higgs_bbar_Matched_AK4_idx]);
            ROOT::Math::PtEtaPhiMVector poh_jetSum = pohb_p4+pohbbar_p4;
            Higgs_bb_Matched_AK4_Mjj = sqrt(poh_jetSum.Dot(poh_jetSum));
        }

        if(POMatching_FullMatch){
            //cout<<"A"<<endl;
            nPOMatched++;
        }
        else{
            //cout<<"B"<<endl;
            nNotPOMatched++;
        }

        Events->GetEntry(myReader.GetCurrentEntry());
        Events_basic->Fill();

        if(!VBF_PO_Matched) nFailPOMatched_FailVBF++;
        if(!Higgs_PO_Matched) nFailPOMatched_FailH++;
        if(Z_GenPart_idx > -1 && !Z_PO_Matched) nFailPOMatched_FailZ++;
        if(Wp_GenPart_idx > -1 && !Wp_PO_Matched) nFailPOMatched_FailWp++;
        if(Wm_GenPart_idx > -1 && !Wm_PO_Matched) nFailPOMatched_FailWm++;
                    
    }
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;


    outfile->cd();
    Events_basic->Write();

    outfile->Close();

    cout<<"nFailedGenMatchedDecay: "<<nFailedGenMatchedDecay
        <<"\nnPassGenMatchedDecay_FailHbb: "<<nPassGenMatchedDecay_FailHbb
        <<"\nnPassGenMatchedDecay_FailZ: "<<nPassGenMatchedDecay_FailZ
        <<"\nnPassGenMatchedDecay_FailHbb_FailZ: "<<nPassGenMatchedDecay_FailHbb_FailZ
        <<"\nnPassGenMatchedDecay_Total: "<<nPassGenMatchedDecay_Total
        <<"\nnNotPOMatched: "<<nNotPOMatched
        <<"\nnPOMatched: "<<nPOMatched
        <<"\nnFailPOMatched_FailVBF: "<<nFailPOMatched_FailVBF
        <<"\nnFailPOMatched_FailH: "<<nFailPOMatched_FailH
        <<"\nnFailPOMatched_FailZ: "<<nFailPOMatched_FailZ
        <<"\nnFailPOMatched_FailWp: "<<nFailPOMatched_FailWp
        <<"\nnFailPOMatched_FailWm: "<<nFailPOMatched_FailWm
        <<endl;

}




void GenLevelAnalyzer_TTreeMaker(std::string inputFile){
    //EventLoop("root://eoscms.cern.ch///eos/cms/store/user/lzygala/HVV/NanoAOD/ucsd/RunIISummer20UL18_VBSWZH_incl_C2V_0_Azure_v1.root");
    EventLoop(inputFile);
}
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
    bool use_new_lepton_id = true;

    //OPTIONS
    bool option_HiggsSelectFirst = true;
    bool option_HiggsDeepTagSort = true;

    float pi=3.14159265358979323846;
    float twopi= 2*pi;

    float deepFlavJetTTHCut = 0.0;

    //int eventLoopMax = EvMax;
    float jetPT_trail_cut = 30.0;
    float jetPT_lead_cut = 30.0;
    float jetInvMass_cut = 500.0;
    float jetEtaSep_cut = 2.0;
    float fatJetPTCut_Higgs = 200.0;
    double fatJetDeepTagCut_Higgs = 0.8;
    double fatJetParticleNetCut_WSL_cut = 0.0;
    double fatJetParticleNetCut_Higgs_cut = 0.0;


    float fatJetPTCut = 200;
    int isEle = 0, isMuon = 1, isTau = 2;
    float leptonPTCut = 10;
    float leptonPTCut_lead = 20;

    double lt_l_cut = 0.0;
    double st_l_cut = 0.0;
    double lt_sl_cut = 0.0;
    double st_sl_cut = 0.0;
    
    //deeptag workingpoints 2018 
    int level_N = 2; //0=0.5%, 1=1.0%, 2=2.5%, 3=5.0%
    float deepTag_WvsQCD[4] = { 0.961, 0.918, 0.762, 0.458};
    float deepTagMD_WvsQCD[4] = { 0.806, 0.704, 0.479, 0.245 };

    int nEvents_cuts_leptonicChannel[10] = {0,0,0,0,0,0,0,0,0,0};
    int nEvents_cuts_semileptonicChannel[10] = {0,0,0,0,0,0,0,0,0,0};
    int nEvents_cuts_hadronicChannel[10] = {0,0,0,0,0,0,0,0,0,0};

    float nEvents_weighted_cuts_leptonicChannel[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    float nEvents_weighted_cuts_semileptonicChannel[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    float nEvents_weighted_cuts_hadronicChannel[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


    Int_t RawEventNumber_TotalEvents = 0.0;
    Int_t RawEventNumber_PreSelect_Jets = 0.0;
    Int_t RawEventNumber_PreSelect_FatJets = 0.0;
    Int_t RawEventNumber_PreSelect_LeptonicChannel = 0.0;
    Int_t RawEventNumber_PreSelect_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_Trigger_LeptonicChannel = 0.0;
    Int_t RawEventNumber_Trigger_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_Leptons_LeptonicChannel = 0.0;
    Int_t RawEventNumber_Leptons_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_HiggsFatJet_LeptonicChannel = 0.0;
    Int_t RawEventNumber_HiggsFatJet_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_BJetVeto_LeptonicChannel = 0.0;
    Int_t RawEventNumber_BJetVeto_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_WFatJet_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_VBFJets_LeptonicChannel = 0.0;
    Int_t RawEventNumber_VBFJets_SemiLeptonicChannel = 0.0;
    Int_t RawEventNumber_FinalPassed = 0.0;
    Int_t RawEventNumber_Passed_Leptonic = 0.0;
    Int_t RawEventNumber_Passed_SemiLeptonic = 0.0;
    Int_t RawEventNumber_Rejected = 0.0;
    Int_t RawEventNumber_Passed_BothChannels = 0.0;

    Float_t WeightedEventNumber_TotalEvents = 0.0;
    Float_t WeightedEventNumber_PreSelect_Jets = 0.0;
    Float_t WeightedEventNumber_PreSelect_FatJets = 0.0;
    Float_t WeightedEventNumber_PreSelect_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_PreSelect_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_Trigger_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_Trigger_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_Leptons_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_Leptons_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_HiggsFatJet_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_HiggsFatJet_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_BJetVeto_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_BJetVeto_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_WFatJet_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_VBFJets_LeptonicChannel = 0.0;
    Float_t WeightedEventNumber_VBFJets_SemiLeptonicChannel = 0.0;
    Float_t WeightedEventNumber_FinalPassed = 0.0;
    Float_t WeightedEventNumber_Passed_Leptonic = 0.0;
    Float_t WeightedEventNumber_Passed_SemiLeptonic = 0.0;
    Float_t WeightedEventNumber_Rejected = 0.0;
    Float_t WeightedEventNumber_Passed_BothChannels = 0.0;

    float test_total = 0.0;
    float test_sem = 0.0;

    //TTree Variables
    TTree* Events_basic;
    TTree* Events_allPassed;
    TTree* Events_leptonic;
    TTree* Events_semileptonic;
    TTree* CutFlow_tree;

    Int_t EventType_Leptonic = -1.0;
    Int_t EventType_SemiLeptonic = -1.0;
    Float_t Xsec_genWeight = -1.0;
    Float_t XS = -1.0;
    Float_t SumWeights = -1.0;
    Char_t Sample_name[30];
    Char_t Year[5];

    Int_t CandidateHiggs_FatJet_idx = -1.0;
    Float_t CandidateHiggs_FatJet_pt = -1.0;
    Float_t CandidateHiggs_FatJet_eta = -1.0;
    Float_t CandidateHiggs_FatJet_phi = -1.0;
    Float_t CandidateHiggs_FatJet_particlenetScore = -1.0;
    Float_t CandidateHiggs_FatJet_mass = -1.0;
    Float_t CandidateHiggs_FatJet_msoftdrop = -1.0;

    Int_t CandidateVBF_Lead_Jet_idx = -1.0;
    Float_t CandidateVBF_Lead_Jet_pt = -1.0;
    Float_t CandidateVBF_Lead_Jet_eta = -1.0;
    Float_t CandidateVBF_Lead_Jet_phi = -1.0;

    Int_t CandidateVBF_Trail_Jet_idx = -1.0;
    Float_t CandidateVBF_Trail_Jet_pt = -1.0;
    Float_t CandidateVBF_Trail_Jet_eta = -1.0;
    Float_t CandidateVBF_Trail_Jet_phi = -1.0;

    Float_t CandidateVBF_Jet_invMass = -1.0;
    Float_t CandidateVBF_Jet_etaSep = -1.0;

    Int_t CandidateW_SemiLeptonic_FatJet_idx = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_pt = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_particlenetScore = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_eta = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_phi = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_msoftdrop = -1.0;
    Float_t CandidateW_SemiLeptonic_FatJet_mass = -1.0;
    Int_t CandidateW_SemiLeptonic_FatJet_flavor = -1.0;

    Int_t Candidate_SemiLeptonic_Lepton_idx = -1.0;
    Int_t Candidate_SemiLeptonic_Lepton_type = -1.0;
    Float_t Candidate_SemiLeptonic_Lepton_pt = -1.0;
    Float_t Candidate_SemiLeptonic_Lepton_eta = -1.0;
    Float_t Candidate_SemiLeptonic_Lepton_phi = -1.0;

    Int_t Candidate_Leptonic_Lead_Lepton_idx = -1.0;
    Int_t Candidate_Leptonic_Lead_Lepton_type = -1.0;
    Float_t Candidate_Leptonic_Lead_Lepton_pt = -1.0;
    Float_t Candidate_Leptonic_Lead_Lepton_eta = -1.0;
    Float_t Candidate_Leptonic_Lead_Lepton_phi = -1.0;

    Int_t Candidate_Leptonic_Trail_Lepton_idx = -1.0;
    Int_t Candidate_Leptonic_Trail_Lepton_type = -1.0;
    Float_t Candidate_Leptonic_Trail_Lepton_pt = -1.0;
    Float_t Candidate_Leptonic_Trail_Lepton_eta = -1.0;
    Float_t Candidate_Leptonic_Trail_Lepton_phi = -1.0;

    Float_t Candidate_Leptonic_Lepton_InvMass = -1.0;

    Float_t Candidate_Leptonic_ST = -1.0;
    Float_t Candidate_SemiLeptonic_ST = -1.0;
    Float_t Candidate_Leptonic_LT = -1.0;
    Float_t Candidate_SemiLeptonic_LT = -1.0;


    Float_t Candidate_Leptonic_MT = -1.0;
    Float_t Candidate_SemiLeptonic_MT = -1.0;

    Float_t Candidate_Leptonic_RpT = -1.0;
    Float_t Candidate_SemiLeptonic_RpT = -1.0;

    Int_t nGoodLeptons = 0.0;
    Int_t nGoodElectrons = 0.0;
    Int_t nGoodMuons = 0.0;

    //TH1F* h_LHE_nJets;


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

double getDR(double phi_1, double eta_1, double phi_2, double eta_2){
    double dPhi = fabs(phi_1 - phi_2);
    if(dPhi > pi) dPhi -= twopi;

    double dEta = eta_1 - eta_2;

    double dR = sqrt( pow(dPhi, 2) + pow(dEta, 2) );

    return dR;

}

void AddBranches(){
    //Candidate Branches in Events Tree
    Events_basic->Branch("EventType_Leptonic", &EventType_Leptonic, "EventType_Leptonic/I");
    Events_basic->Branch("EventType_SemiLeptonic", &EventType_SemiLeptonic, "EventType_SemiLeptonic/I");
    Events_basic->Branch("Xsec_genWeight", &Xsec_genWeight, "Xsec_genWeight/F");
    Events_basic->Branch("XS", &XS, "XS/F");
    Events_basic->Branch("SumWeights", &SumWeights, "SumWeights/F");
    Events_basic->Branch("Sample_name", &Sample_name, "Sample_name/C");
    Events_basic->Branch("Year", &Year, "Year/C");

    Events_basic->Branch("CandidateHiggs_FatJet_idx", &CandidateHiggs_FatJet_idx, "CandidateHiggs_FatJet_idx/I");
    Events_basic->Branch("CandidateHiggs_FatJet_pt", &CandidateHiggs_FatJet_pt, "CandidateHiggs_FatJet_pt/F");
    Events_basic->Branch("CandidateHiggs_FatJet_eta", &CandidateHiggs_FatJet_eta, "CandidateHiggs_FatJet_eta/F");
    Events_basic->Branch("CandidateHiggs_FatJet_phi", &CandidateHiggs_FatJet_phi, "CandidateHiggs_FatJet_phi/F");
    Events_basic->Branch("CandidateHiggs_FatJet_particlenetScore", &CandidateHiggs_FatJet_particlenetScore, "CandidateHiggs_FatJet_particlenetScore/F");
    Events_basic->Branch("CandidateHiggs_FatJet_mass", &CandidateHiggs_FatJet_mass, "CandidateHiggs_FatJet_mass/F");
    Events_basic->Branch("CandidateHiggs_FatJet_msoftdrop", &CandidateHiggs_FatJet_msoftdrop, "CandidateHiggs_FatJet_msoftdrop/F");

    Events_basic->Branch("CandidateVBF_Lead_Jet_idx", &CandidateVBF_Lead_Jet_idx, "CandidateVBF_Lead_Jet_idx/I");
    Events_basic->Branch("CandidateVBF_Lead_Jet_pt", &CandidateVBF_Lead_Jet_pt, "CandidateVBF_Lead_Jet_pt/F");
    Events_basic->Branch("CandidateVBF_Lead_Jet_eta", &CandidateVBF_Lead_Jet_eta, "CandidateVBF_Lead_Jet_eta/F");
    Events_basic->Branch("CandidateVBF_Lead_Jet_phi", &CandidateVBF_Lead_Jet_phi, "CandidateVBF_Lead_Jet_phi/F");

    Events_basic->Branch("CandidateVBF_Trail_Jet_idx", &CandidateVBF_Trail_Jet_idx, "CandidateVBF_Trail_Jet_idx/I");
    Events_basic->Branch("CandidateVBF_Trail_Jet_pt", &CandidateVBF_Trail_Jet_pt, "CandidateVBF_Trail_Jet_pt/F");
    Events_basic->Branch("CandidateVBF_Trail_Jet_eta", &CandidateVBF_Trail_Jet_eta, "CandidateVBF_Trail_Jet_eta/F");
    Events_basic->Branch("CandidateVBF_Trail_Jet_phi", &CandidateVBF_Trail_Jet_phi, "CandidateVBF_Trail_Jet_phi/F");

    Events_basic->Branch("CandidateVBF_Jet_invMass", &CandidateVBF_Jet_invMass, "CandidateVBF_Jet_invMass/F");
    Events_basic->Branch("CandidateVBF_Jet_etaSep", &CandidateVBF_Jet_etaSep, "CandidateVBF_Jet_etaSep/F");

    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_idx", &CandidateW_SemiLeptonic_FatJet_idx, "CandidateW_SemiLeptonic_FatJet_idx/I");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_particlenetScore", &CandidateW_SemiLeptonic_FatJet_particlenetScore, "CandidateW_SemiLeptonic_FatJet_particlenetScore/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_pt", &CandidateW_SemiLeptonic_FatJet_pt, "CandidateW_SemiLeptonic_FatJet_pt/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_eta", &CandidateW_SemiLeptonic_FatJet_eta, "CandidateW_SemiLeptonic_FatJet_eta/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_phi", &CandidateW_SemiLeptonic_FatJet_phi, "CandidateW_SemiLeptonic_FatJet_phi/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_mass", &CandidateW_SemiLeptonic_FatJet_mass, "CandidateW_SemiLeptonic_FatJet_mass/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_msoftdrop", &CandidateW_SemiLeptonic_FatJet_msoftdrop, "CandidateW_SemiLeptonic_FatJet_msoftdrop/F");
    Events_basic->Branch("CandidateW_SemiLeptonic_FatJet_flavor", &CandidateW_SemiLeptonic_FatJet_flavor, "CandidateW_SemiLeptonic_FatJet_flavor/I");

    Events_basic->Branch("Candidate_SemiLeptonic_Lepton_idx", &Candidate_SemiLeptonic_Lepton_idx, "Candidate_SemiLeptonic_Lepton_idx/I");
    Events_basic->Branch("Candidate_SemiLeptonic_Lepton_type", &Candidate_SemiLeptonic_Lepton_type, "Candidate_SemiLeptonic_Lepton_type/I");
    Events_basic->Branch("Candidate_SemiLeptonic_Lepton_pt", &Candidate_SemiLeptonic_Lepton_pt, "Candidate_SemiLeptonic_Lepton_pt/F");
    Events_basic->Branch("Candidate_SemiLeptonic_Lepton_eta", &Candidate_SemiLeptonic_Lepton_eta, "Candidate_SemiLeptonic_Lepton_eta/F");
    Events_basic->Branch("Candidate_SemiLeptonic_Lepton_phi", &Candidate_SemiLeptonic_Lepton_phi, "Candidate_SemiLeptonic_Lepton_phi/F");

    Events_basic->Branch("Candidate_Leptonic_Lead_Lepton_idx", &Candidate_Leptonic_Lead_Lepton_idx, "Candidate_Leptonic_Lead_Lepton_idx/I");
    Events_basic->Branch("Candidate_Leptonic_Lead_Lepton_type", &Candidate_Leptonic_Lead_Lepton_type, "Candidate_Leptonic_Lead_Lepton_type/I");
    Events_basic->Branch("Candidate_Leptonic_Lead_Lepton_pt", &Candidate_Leptonic_Lead_Lepton_pt, "Candidate_Leptonic_Lead_Lepton_pt/F");
    Events_basic->Branch("Candidate_Leptonic_Lead_Lepton_eta", &Candidate_Leptonic_Lead_Lepton_eta, "Candidate_Leptonic_Lead_Lepton_eta/F");
    Events_basic->Branch("Candidate_Leptonic_Lead_Lepton_phi", &Candidate_Leptonic_Lead_Lepton_phi, "Candidate_Leptonic_Lead_Lepton_phi/F");

    Events_basic->Branch("Candidate_Leptonic_Trail_Lepton_idx", &Candidate_Leptonic_Trail_Lepton_idx, "Candidate_Leptonic_Trail_Lepton_idx/I");
    Events_basic->Branch("Candidate_Leptonic_Trail_Lepton_type", &Candidate_Leptonic_Trail_Lepton_type, "Candidate_Leptonic_Trail_Lepton_type/I");
    Events_basic->Branch("Candidate_Leptonic_Trail_Lepton_pt", &Candidate_Leptonic_Trail_Lepton_pt, "Candidate_Leptonic_Trail_Lepton_pt/F");
    Events_basic->Branch("Candidate_Leptonic_Trail_Lepton_eta", &Candidate_Leptonic_Trail_Lepton_eta, "Candidate_Leptonic_Trail_Lepton_eta/F");
    Events_basic->Branch("Candidate_Leptonic_Trail_Lepton_phi", &Candidate_Leptonic_Trail_Lepton_phi, "Candidate_Leptonic_Trail_Lepton_phi/F");

    Events_basic->Branch("Candidate_Leptonic_Lepton_InvMass", &Candidate_Leptonic_Lepton_InvMass, "Candidate_Leptonic_Lepton_InvMass/F");

    Events_basic->Branch("Candidate_Leptonic_ST", &Candidate_Leptonic_ST, "Candidate_Leptonic_ST/F");
    Events_basic->Branch("Candidate_SemiLeptonic_ST", &Candidate_SemiLeptonic_ST, "Candidate_SemiLeptonic_ST/F");
    Events_basic->Branch("Candidate_Leptonic_LT", &Candidate_Leptonic_LT, "Candidate_Leptonic_LT/F");
    Events_basic->Branch("Candidate_SemiLeptonic_LT", &Candidate_SemiLeptonic_LT, "Candidate_SemiLeptonic_LT/F");
    Events_basic->Branch("Candidate_Leptonic_MT", &Candidate_Leptonic_MT, "Candidate_Leptonic_MT/F");
    Events_basic->Branch("Candidate_SemiLeptonic_MT", &Candidate_SemiLeptonic_MT, "Candidate_SemiLeptonic_MT/F");
    Events_basic->Branch("Candidate_Leptonic_RpT", &Candidate_Leptonic_RpT, "Candidate_Leptonic_RpT/F");
    Events_basic->Branch("Candidate_SemiLeptonic_RpT", &Candidate_SemiLeptonic_RpT, "Candidate_SemiLeptonic_RpT/F");

    //CutFlow Branches
    CutFlow_tree->Branch("RawEventNumber_TotalEvents", &RawEventNumber_TotalEvents, "RawEventNumber_TotalEvents/I");
    CutFlow_tree->Branch("RawEventNumber_PreSelect_Jets", &RawEventNumber_PreSelect_Jets, "RawEventNumber_PreSelect_Jets/I");
    CutFlow_tree->Branch("RawEventNumber_PreSelect_FatJets", &RawEventNumber_PreSelect_FatJets, "RawEventNumber_PreSelect_FatJets/I");
    CutFlow_tree->Branch("RawEventNumber_PreSelect_LeptonicChannel", &RawEventNumber_PreSelect_LeptonicChannel, "RawEventNumber_PreSelect_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_PreSelect_SemiLeptonicChannel", &RawEventNumber_PreSelect_SemiLeptonicChannel, "RawEventNumber_PreSelect_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_Leptons_LeptonicChannel", &RawEventNumber_Leptons_LeptonicChannel, "RawEventNumber_Leptons_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_Leptons_SemiLeptonicChannel", &RawEventNumber_Leptons_SemiLeptonicChannel, "RawEventNumber_Leptons_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_Trigger_LeptonicChannel", &RawEventNumber_Trigger_LeptonicChannel, "RawEventNumber_Trigger_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_Trigger_SemiLeptonicChannel", &RawEventNumber_Trigger_SemiLeptonicChannel, "RawEventNumber_Trigger_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_HiggsFatJet_LeptonicChannel", &RawEventNumber_HiggsFatJet_LeptonicChannel, "RawEventNumber_HiggsFatJet_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_HiggsFatJet_SemiLeptonicChannel", &RawEventNumber_HiggsFatJet_SemiLeptonicChannel, "RawEventNumber_HiggsFatJet_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_WFatJet_SemiLeptonicChannel", &RawEventNumber_WFatJet_SemiLeptonicChannel, "RawEventNumber_WFatJet_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_BJetVeto_LeptonicChannel", &RawEventNumber_BJetVeto_LeptonicChannel, "RawEventNumber_BJetVeto_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_BJetVeto_SemiLeptonicChannel", &RawEventNumber_BJetVeto_SemiLeptonicChannel, "RawEventNumber_BJetVeto_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_VBFJets_LeptonicChannel", &RawEventNumber_VBFJets_LeptonicChannel, "RawEventNumber_VBFJets_LeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_VBFJets_SemiLeptonicChannel", &RawEventNumber_VBFJets_SemiLeptonicChannel, "RawEventNumber_VBFJets_SemiLeptonicChannel/I");
    CutFlow_tree->Branch("RawEventNumber_FinalPassed", &RawEventNumber_FinalPassed, "RawEventNumber_FinalPassed/I");
    CutFlow_tree->Branch("RawEventNumber_Passed_Leptonic", &RawEventNumber_Passed_Leptonic, "RawEventNumber_Passed_Leptonic/I");
    CutFlow_tree->Branch("RawEventNumber_Passed_SemiLeptonic", &RawEventNumber_Passed_SemiLeptonic, "RawEventNumber_Passed_SemiLeptonic/I");
    CutFlow_tree->Branch("RawEventNumber_Rejected", &RawEventNumber_Rejected, "RawEventNumber_Rejected/I");
    CutFlow_tree->Branch("RawEventNumber_Passed_BothChannels", &RawEventNumber_Passed_BothChannels, "RawEventNumber_Passed_BothChannels/I");

    CutFlow_tree->Branch("WeightedEventNumber_TotalEvents", &WeightedEventNumber_TotalEvents, "WeightedEventNumber_TotalEvents/F");
    CutFlow_tree->Branch("WeightedEventNumber_PreSelect_Jets", &WeightedEventNumber_PreSelect_Jets, "WeightedEventNumber_PreSelect_Jets/F");
    CutFlow_tree->Branch("WeightedEventNumber_PreSelect_FatJets", &WeightedEventNumber_PreSelect_FatJets, "WeightedEventNumber_PreSelect_FatJets/F");
    CutFlow_tree->Branch("WeightedEventNumber_PreSelect_LeptonicChannel", &WeightedEventNumber_PreSelect_LeptonicChannel, "WeightedEventNumber_PreSelect_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_PreSelect_SemiLeptonicChannel", &WeightedEventNumber_PreSelect_SemiLeptonicChannel, "WeightedEventNumber_PreSelect_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_Trigger_LeptonicChannel", &WeightedEventNumber_Trigger_LeptonicChannel, "WeightedEventNumber_Trigger_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_Trigger_SemiLeptonicChannel", &WeightedEventNumber_Trigger_SemiLeptonicChannel, "WeightedEventNumber_Trigger_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_Leptons_LeptonicChannel", &WeightedEventNumber_Leptons_LeptonicChannel, "WeightedEventNumber_Leptons_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_Leptons_SemiLeptonicChannel", &WeightedEventNumber_Leptons_SemiLeptonicChannel, "WeightedEventNumber_Leptons_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_HiggsFatJet_LeptonicChannel", &WeightedEventNumber_HiggsFatJet_LeptonicChannel, "WeightedEventNumber_HiggsFatJet_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_HiggsFatJet_SemiLeptonicChannel", &WeightedEventNumber_HiggsFatJet_SemiLeptonicChannel, "WeightedEventNumber_HiggsFatJet_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_WFatJet_SemiLeptonicChannel", &WeightedEventNumber_WFatJet_SemiLeptonicChannel, "WeightedEventNumber_WFatJet_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_BJetVeto_LeptonicChannel", &WeightedEventNumber_BJetVeto_LeptonicChannel, "WeightedEventNumber_BJetVeto_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_BJetVeto_SemiLeptonicChannel", &WeightedEventNumber_BJetVeto_SemiLeptonicChannel, "WeightedEventNumber_BJetVeto_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_VBFJets_LeptonicChannel", &WeightedEventNumber_VBFJets_LeptonicChannel, "WeightedEventNumber_VBFJets_LeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_VBFJets_SemiLeptonicChannel", &WeightedEventNumber_VBFJets_SemiLeptonicChannel, "WeightedEventNumber_VBFJets_SemiLeptonicChannel/F");
    CutFlow_tree->Branch("WeightedEventNumber_FinalPassed", &WeightedEventNumber_FinalPassed, "WeightedEventNumber_FinalPassed/F");
    CutFlow_tree->Branch("WeightedEventNumber_Passed_Leptonic", &WeightedEventNumber_Passed_Leptonic, "WeightedEventNumber_Passed_Leptonic/F");
    CutFlow_tree->Branch("WeightedEventNumber_Passed_SemiLeptonic", &WeightedEventNumber_Passed_SemiLeptonic, "WeightedEventNumber_Passed_SemiLeptonic/F");
    CutFlow_tree->Branch("WeightedEventNumber_Rejected", &WeightedEventNumber_Rejected, "WeightedEventNumber_Rejected/F");
    CutFlow_tree->Branch("WeightedEventNumber_Passed_BothChannels", &WeightedEventNumber_Passed_BothChannels, "WeightedEventNumber_Passed_BothChannels/F");

    return;
}

void EventLoop_background(TString infileName, string year, float xsec, float sum_weights, float lumi){
//Loops over all events and fills histograms
    bool debug = false;

    cout<<cur_time()<<"\tProcessing Tree...\n";
    cout<<cur_time()<<"\tFile Opened...\n";
    TFile* infile = TFile::Open(infileName);
    TTree* Events = (TTree*)infile->Get("Events");

    cout<<cur_time()<<"\tCreating Output File"<<endl;
    TFile *outfile = new TFile("TTree_tmp.root", "RECREATE");
    outfile->cd();

    CutFlow_tree = new TTree("CutFlow","CutFlow");

    Events_basic = Events->CloneTree(0);
    AddBranches();

    //Split Events into 3 differenct branches
    Events_allPassed = Events_basic->CloneTree(0);
    Events_allPassed->SetObject("Events_allPassed","Events_allPassed");
    Events_leptonic = Events_basic->CloneTree(0);
    Events_leptonic->SetObject("Events_leptonic","Events_leptonic");
    Events_semileptonic= Events_basic->CloneTree(0);
    Events_semileptonic->SetObject("Events_semileptonic","Events_semileptonic");


    TTreeReader myReader("Events", infile);

    TTreeReaderValue<UInt_t> run(myReader, "run");
    TTreeReaderValue<UInt_t> luminosityBlock(myReader, "luminosityBlock");
    TTreeReaderValue<ULong64_t> event(myReader, "event");
    //TTreeReaderValue<Float_t> btagWeight_CSVV2(myReader, "btagWeight_CSVV2");
    //TTreeReaderValue<Float_t> btagWeight_DeepCSVB(myReader, "btagWeight_DeepCSVB");

    TTreeReaderValue<UInt_t> nElectron(myReader, "nElectron");
    TTreeReaderArray<Float_t> Electron_deltaEtaSC(myReader, "Electron_deltaEtaSC");
    //TTreeReaderArray<Float_t> Electron_dr03EcalRecHitSumEt(myReader, "Electron_dr03EcalRecHitSumEt");
    //TTreeReaderArray<Float_t> Electron_dr03HcalDepth1TowerSumEt(myReader, "Electron_dr03HcalDepth1TowerSumEt");
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
    TTreeReaderArray<Float_t> Electron_miniPFRelIso_all(myReader, "Electron_miniPFRelIso_all");
    TTreeReaderArray<Float_t> Electron_miniPFRelIso_chg(myReader, "Electron_miniPFRelIso_chg");
    TTreeReaderArray<Float_t> Electron_mvaFall17V2Iso(myReader, "Electron_mvaFall17V2Iso");
    TTreeReaderArray<Float_t> Electron_mvaFall17V2noIso(myReader, "Electron_mvaFall17V2noIso");
    TTreeReaderArray<Float_t> Electron_pfRelIso03_all(myReader, "Electron_pfRelIso03_all");
    TTreeReaderArray<Float_t> Electron_pfRelIso03_chg(myReader, "Electron_pfRelIso03_chg");
    TTreeReaderArray<Float_t> Electron_phi(myReader, "Electron_phi");
    TTreeReaderArray<Float_t> Electron_pt(myReader, "Electron_pt");
    TTreeReaderArray<Float_t> Electron_r9(myReader, "Electron_r9");
    TTreeReaderArray<Float_t> Electron_sieie(myReader, "Electron_sieie");
    TTreeReaderArray<Float_t> Electron_sip3d(myReader, "Electron_sip3d");
    TTreeReaderArray<Float_t> Electron_mvaTTH(myReader, "Electron_mvaTTH");
    //TTreeReaderArray<Float_t> Electron_mvaTTHUL(myReader, "Electron_mvaTTHUL");
    TTreeReaderArray<Int_t> Electron_charge(myReader, "Electron_charge");
    TTreeReaderArray<Int_t> Electron_cutBased(myReader, "Electron_cutBased");
    TTreeReaderArray<Int_t> Electron_jetIdx(myReader, "Electron_jetIdx");
    TTreeReaderArray<Int_t> Electron_pdgId(myReader, "Electron_pdgId");
    TTreeReaderArray<Int_t> Electron_photonIdx(myReader, "Electron_photonIdx");
    TTreeReaderArray<Int_t> Electron_tightCharge(myReader, "Electron_tightCharge");
    TTreeReaderArray<Int_t> Electron_vidNestedWPBitmap(myReader, "Electron_vidNestedWPBitmap");
    TTreeReaderArray<Bool_t> Electron_convVeto(myReader, "Electron_convVeto");
    TTreeReaderArray<Bool_t> Electron_cutBased_HEEP(myReader, "Electron_cutBased_HEEP");
    TTreeReaderArray<Bool_t> Electron_isPFcand(myReader, "Electron_isPFcand");
    TTreeReaderArray<UChar_t> Electron_lostHits(myReader, "Electron_lostHits");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP80(myReader, "Electron_mvaFall17V2Iso_WP80");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP90(myReader, "Electron_mvaFall17V2Iso_WP90");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WPL(myReader, "Electron_mvaFall17V2Iso_WPL");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP80(myReader, "Electron_mvaFall17V2noIso_WP80");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP90(myReader, "Electron_mvaFall17V2noIso_WP90");
    TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WPL(myReader, "Electron_mvaFall17V2noIso_WPL");
    TTreeReaderArray<UChar_t> Electron_seedGain(myReader, "Electron_seedGain");

    TTreeReaderValue<UInt_t> nFatJet(myReader, "nFatJet");
    TTreeReaderArray<Float_t> FatJet_area(myReader, "FatJet_area");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_H4qvsQCD(myReader, "FatJet_deepTagMD_H4qvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_HbbvsQCD(myReader, "FatJet_deepTagMD_HbbvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_TvsQCD(myReader, "FatJet_deepTagMD_TvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_WvsQCD(myReader, "FatJet_deepTagMD_WvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHbbvsQCD(myReader, "FatJet_deepTagMD_ZHbbvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHccvsQCD(myReader, "FatJet_deepTagMD_ZHccvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_ZbbvsQCD(myReader, "FatJet_deepTagMD_ZbbvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_ZvsQCD(myReader, "FatJet_deepTagMD_ZvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_bbvsLight(myReader, "FatJet_deepTagMD_bbvsLight");
    TTreeReaderArray<Float_t> FatJet_deepTagMD_ccvsLight(myReader, "FatJet_deepTagMD_ccvsLight");
    TTreeReaderArray<Float_t> FatJet_deepTag_H(myReader, "FatJet_deepTag_H");
    TTreeReaderArray<Float_t> FatJet_deepTag_QCD(myReader, "FatJet_deepTag_QCD");
    TTreeReaderArray<Float_t> FatJet_deepTag_QCDothers(myReader, "FatJet_deepTag_QCDothers");
    TTreeReaderArray<Float_t> FatJet_deepTag_TvsQCD(myReader, "FatJet_deepTag_TvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTag_WvsQCD(myReader, "FatJet_deepTag_WvsQCD");
    TTreeReaderArray<Float_t> FatJet_deepTag_ZvsQCD(myReader, "FatJet_deepTag_ZvsQCD");
    TTreeReaderArray<Float_t> FatJet_eta(myReader, "FatJet_eta");
    TTreeReaderArray<Float_t> FatJet_mass(myReader, "FatJet_mass");
    TTreeReaderArray<Float_t> FatJet_msoftdrop(myReader, "FatJet_msoftdrop");
    TTreeReaderArray<Float_t> FatJet_n2b1(myReader, "FatJet_n2b1");
    TTreeReaderArray<Float_t> FatJet_n3b1(myReader, "FatJet_n3b1");
    TTreeReaderArray<Float_t> FatJet_phi(myReader, "FatJet_phi");
    TTreeReaderArray<Float_t> FatJet_pt(myReader, "FatJet_pt");
    TTreeReaderArray<Float_t> FatJet_rawFactor(myReader, "FatJet_rawFactor");
    TTreeReaderArray<Float_t> FatJet_tau1(myReader, "FatJet_tau1");
    TTreeReaderArray<Float_t> FatJet_tau2(myReader, "FatJet_tau2");
    TTreeReaderArray<Float_t> FatJet_tau3(myReader, "FatJet_tau3");
    TTreeReaderArray<Float_t> FatJet_tau4(myReader, "FatJet_tau4");
    TTreeReaderArray<Int_t> FatJet_jetId(myReader, "FatJet_jetId");
    TTreeReaderArray<Int_t> FatJet_subJetIdx1(myReader, "FatJet_subJetIdx1");
    TTreeReaderArray<Int_t> FatJet_subJetIdx2(myReader, "FatJet_subJetIdx2");

    TTreeReaderArray<Float_t> FatJet_particleNet_HbbvsQCD(myReader, "FatJet_particleNet_HbbvsQCD");
    TTreeReaderArray<Float_t> FatJet_particleNet_WvsQCD(myReader, "FatJet_particleNet_WvsQCD");
    TTreeReaderArray<Float_t> FatJet_particleNet_ZvsQCD(myReader, "FatJet_particleNet_ZvsQCD");

    TTreeReaderValue<Float_t> genWeight(myReader, "genWeight");
    TTreeReaderValue<UInt_t> nPSWeight(myReader, "nPSWeight");
    TTreeReaderArray<Float_t> PSWeight(myReader, "PSWeight");
    //TTreeReaderValue<UInt_t> nIsoTrack(myReader, "nIsoTrack");

    TTreeReaderValue<UInt_t> nJet(myReader, "nJet");
    TTreeReaderArray<Float_t> Jet_area(myReader, "Jet_area");
    TTreeReaderArray<Float_t> Jet_btagCSVV2(myReader, "Jet_btagCSVV2");
    TTreeReaderArray<Float_t> Jet_btagDeepB(myReader, "Jet_btagDeepB");
    TTreeReaderArray<Float_t> Jet_btagDeepFlavB(myReader, "Jet_btagDeepFlavB");
    TTreeReaderArray<Float_t> Jet_chEmEF(myReader, "Jet_chEmEF");
    TTreeReaderArray<Float_t> Jet_chHEF(myReader, "Jet_chHEF");
    TTreeReaderArray<Float_t> Jet_eta(myReader, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_mass(myReader, "Jet_mass");
    TTreeReaderArray<Float_t> Jet_muEF(myReader, "Jet_muEF");
    TTreeReaderArray<Float_t> Jet_muonSubtrFactor(myReader, "Jet_muonSubtrFactor");
    TTreeReaderArray<Float_t> Jet_neEmEF(myReader, "Jet_neEmEF");
    TTreeReaderArray<Float_t> Jet_neHEF(myReader, "Jet_neHEF");
    TTreeReaderArray<Float_t> Jet_phi(myReader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_pt(myReader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_qgl(myReader, "Jet_qgl");
    TTreeReaderArray<Float_t> Jet_rawFactor(myReader, "Jet_rawFactor");
    TTreeReaderArray<Float_t> Jet_bRegCorr(myReader, "Jet_bRegCorr");
    TTreeReaderArray<Float_t> Jet_bRegRes(myReader, "Jet_bRegRes");
    TTreeReaderArray<Int_t> Jet_electronIdx1(myReader, "Jet_electronIdx1");
    TTreeReaderArray<Int_t> Jet_electronIdx2(myReader, "Jet_electronIdx2");
    TTreeReaderArray<Int_t> Jet_jetId(myReader, "Jet_jetId");
    TTreeReaderArray<Int_t> Jet_muonIdx1(myReader, "Jet_muonIdx1");
    TTreeReaderArray<Int_t> Jet_muonIdx2(myReader, "Jet_muonIdx2");
    TTreeReaderArray<Int_t> Jet_nElectrons(myReader, "Jet_nElectrons");
    TTreeReaderArray<Int_t> Jet_nMuons(myReader, "Jet_nMuons");
    TTreeReaderArray<Int_t> Jet_puId(myReader, "Jet_puId");

    TTreeReaderValue<UInt_t> nMuon(myReader, "nMuon");
    TTreeReaderArray<Float_t> Muon_dxy(myReader, "Muon_dxy");
    TTreeReaderArray<Float_t> Muon_dxyErr(myReader, "Muon_dxyErr");
    TTreeReaderArray<Float_t> Muon_dz(myReader, "Muon_dz");
    TTreeReaderArray<Float_t> Muon_dzErr(myReader, "Muon_dzErr");
    TTreeReaderArray<Float_t> Muon_eta(myReader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_ip3d(myReader, "Muon_ip3d");
    TTreeReaderArray<Float_t> Muon_jetPtRelv2(myReader, "Muon_jetPtRelv2");
    TTreeReaderArray<Float_t> Muon_jetRelIso(myReader, "Muon_jetRelIso");
    TTreeReaderArray<Float_t> Muon_mass(myReader, "Muon_mass");
    TTreeReaderArray<Float_t> Muon_miniPFRelIso_all(myReader, "Muon_miniPFRelIso_all");
    TTreeReaderArray<Float_t> Muon_miniPFRelIso_chg(myReader, "Muon_miniPFRelIso_chg");
    TTreeReaderArray<Float_t> Muon_pfRelIso03_all(myReader, "Muon_pfRelIso03_all");
    TTreeReaderArray<Float_t> Muon_pfRelIso03_chg(myReader, "Muon_pfRelIso03_chg");
    TTreeReaderArray<Float_t> Muon_pfRelIso04_all(myReader, "Muon_pfRelIso04_all");
    TTreeReaderArray<Float_t> Muon_phi(myReader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_pt(myReader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_ptErr(myReader, "Muon_ptErr");
    TTreeReaderArray<Float_t> Muon_segmentComp(myReader, "Muon_segmentComp");
    TTreeReaderArray<Float_t> Muon_sip3d(myReader, "Muon_sip3d");
    TTreeReaderArray<Float_t> Muon_softMva(myReader, "Muon_softMva");
    TTreeReaderArray<Float_t> Muon_tkRelIso(myReader, "Muon_tkRelIso");
    TTreeReaderArray<Float_t> Muon_tunepRelPt(myReader, "Muon_tunepRelPt");
    TTreeReaderArray<Float_t> Muon_mvaLowPt(myReader, "Muon_mvaLowPt");
    TTreeReaderArray<Float_t> Muon_mvaTTH(myReader, "Muon_mvaTTH");
    //TTreeReaderArray<Float_t> Muon_mvaTTHUL(myReader, "Muon_mvaTTHUL");
    TTreeReaderArray<Int_t> Muon_charge(myReader, "Muon_charge");
    TTreeReaderArray<Int_t> Muon_jetIdx(myReader, "Muon_jetIdx");
    TTreeReaderArray<Int_t> Muon_nStations(myReader, "Muon_nStations");
    TTreeReaderArray<Int_t> Muon_nTrackerLayers(myReader, "Muon_nTrackerLayers");
    TTreeReaderArray<Int_t> Muon_pdgId(myReader, "Muon_pdgId");
    TTreeReaderArray<Int_t> Muon_tightCharge(myReader, "Muon_tightCharge");
    TTreeReaderArray<UChar_t> Muon_highPtId(myReader, "Muon_highPtId");
    TTreeReaderArray<Bool_t> Muon_inTimeMuon(myReader, "Muon_inTimeMuon");
    TTreeReaderArray<Bool_t> Muon_isGlobal(myReader, "Muon_isGlobal");
    TTreeReaderArray<Bool_t> Muon_isPFcand(myReader, "Muon_isPFcand");
    TTreeReaderArray<Bool_t> Muon_isTracker(myReader, "Muon_isTracker");
    TTreeReaderArray<Bool_t> Muon_looseId(myReader, "Muon_looseId");
    TTreeReaderArray<Bool_t> Muon_mediumId(myReader, "Muon_mediumId");
    TTreeReaderArray<Bool_t> Muon_mediumPromptId(myReader, "Muon_mediumPromptId");
    TTreeReaderArray<UChar_t> Muon_miniIsoId(myReader, "Muon_miniIsoId");
    TTreeReaderArray<UChar_t> Muon_multiIsoId(myReader, "Muon_multiIsoId");
    TTreeReaderArray<UChar_t> Muon_mvaId(myReader, "Muon_mvaId");
    TTreeReaderArray<UChar_t> Muon_pfIsoId(myReader, "Muon_pfIsoId");
    TTreeReaderArray<UChar_t> Muon_puppiIsoId(myReader, "Muon_puppiIsoId");
    TTreeReaderArray<Bool_t> Muon_softId(myReader, "Muon_softId");
    TTreeReaderArray<Bool_t> Muon_softMvaId(myReader, "Muon_softMvaId");
    TTreeReaderArray<Bool_t> Muon_tightId(myReader, "Muon_tightId");
    TTreeReaderArray<UChar_t> Muon_tkIsoId(myReader, "Muon_tkIsoId");
    TTreeReaderArray<Bool_t> Muon_triggerIdLoose(myReader, "Muon_triggerIdLoose");

    //TTreeReaderValue<Float_t> PuppiMET_phi(myReader, "PuppiMET_phi");
    //TTreeReaderValue<Float_t> PuppiMET_pt(myReader, "PuppiMET_pt");
    TTreeReaderValue<Float_t> MET_pt(myReader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(myReader, "MET_phi");
    //TTreeReaderValue<Float_t> MET_eta(myReader, "MET_eta");
    TTreeReaderValue<Float_t> MET_sumEt(myReader, "MET_sumEt");
    //TTreeReaderValue<Float_t> PuppiMET_sumEt(myReader, "PuppiMET_sumEt");
    //TTreeReaderValue<Bool_t> HLT_Ele25_eta2p1_WPTight_Gsf(myReader, "HLT_Ele25_eta2p1_WPTight_Gsf");
    //TTreeReaderValue<Bool_t> HLT_IsoTkMu24(myReader, "HLT_IsoTkMu24");
    
    TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(myReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
    TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(myReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL");
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(myReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(myReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(myReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(myReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL");

    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8(myReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8(myReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8(myReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8");
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8(myReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8");

    TTreeReaderValue<Bool_t> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(myReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    TTreeReaderValue<Bool_t> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL(myReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
    TTreeReaderValue<Bool_t> HLT_DiEle27_WPTightCaloOnly_L1DoubleEG(myReader, "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG");
    TTreeReaderValue<Bool_t> HLT_DoubleEle33_CaloIdL_MW(myReader, "HLT_DoubleEle33_CaloIdL_MW");
    TTreeReaderValue<Bool_t> HLT_DoubleEle25_CaloIdL_MW(myReader, "HLT_DoubleEle25_CaloIdL_MW");
    TTreeReaderValue<Bool_t> HLT_DoubleEle27_CaloIdL_MW(myReader, "HLT_DoubleEle27_CaloIdL_MW");
    TTreeReaderValue<Bool_t> HLT_DoublePhoton70(myReader, "HLT_DoublePhoton70");

    TTreeReaderValue<Bool_t> HLT_IsoMu24(myReader, "HLT_IsoMu24");
    TTreeReaderValue<Bool_t> HLT_IsoMu27(myReader, "HLT_IsoMu27");
    TTreeReaderValue<Bool_t> HLT_IsoMu30(myReader, "HLT_IsoMu30");
    TTreeReaderValue<Bool_t> HLT_Mu50(myReader, "HLT_Mu50");

    TTreeReaderValue<Bool_t> HLT_Ele115_CaloIdVT_GsfTrkIdT(myReader, "HLT_Ele115_CaloIdVT_GsfTrkIdT");
    TTreeReaderValue<Bool_t> HLT_Ele27_WPTight_Gsf(myReader, "HLT_Ele27_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele28_WPTight_Gsf(myReader, "HLT_Ele28_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf(myReader, "HLT_Ele32_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele35_WPTight_Gsf(myReader, "HLT_Ele35_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele38_WPTight_Gsf(myReader, "HLT_Ele38_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele40_WPTight_Gsf(myReader, "HLT_Ele40_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf_L1DoubleEG(myReader, "HLT_Ele32_WPTight_Gsf_L1DoubleEG");
    TTreeReaderValue<Bool_t> HLT_Photon200(myReader, "HLT_Photon200");

    TTreeReaderValue<Bool_t> HLT_AK8PFJet500(myReader, "HLT_AK8PFJet500");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet360_TrimMass30(myReader, "HLT_AK8PFJet360_TrimMass30");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet380_TrimMass30(myReader, "HLT_AK8PFJet380_TrimMass30");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet400_TrimMass30(myReader, "HLT_AK8PFJet400_TrimMass30");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet420_TrimMass30(myReader, "HLT_AK8PFJet420_TrimMass30");
    TTreeReaderValue<Bool_t> HLT_AK8PFHT750_TrimMass50(myReader, "HLT_AK8PFHT750_TrimMass50");
    TTreeReaderValue<Bool_t> HLT_AK8PFHT800_TrimMass50(myReader, "HLT_AK8PFHT800_TrimMass50");
    TTreeReaderValue<Bool_t> HLT_AK8PFHT850_TrimMass50(myReader, "HLT_AK8PFHT850_TrimMass50");
    TTreeReaderValue<Bool_t> HLT_AK8PFHT900_TrimMass50(myReader, "HLT_AK8PFHT900_TrimMass50");
    TTreeReaderValue<Bool_t> HLT_PFHT1050(myReader, "HLT_PFHT1050");

    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2(myReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4(myReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4");
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02(myReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02");


    
    
    
    TTreeReaderValue<Bool_t> Flag_goodVertices(myReader, "Flag_goodVertices");
    TTreeReaderValue<Bool_t> Flag_HBHENoiseFilter(myReader, "Flag_HBHENoiseFilter");
    TTreeReaderValue<Bool_t> Flag_HBHENoiseIsoFilter(myReader, "Flag_HBHENoiseIsoFilter");
    TTreeReaderValue<Bool_t> Flag_EcalDeadCellTriggerPrimitiveFilter(myReader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
    TTreeReaderValue<Bool_t> Flag_BadPFMuonFilter(myReader, "Flag_BadPFMuonFilter");
    TTreeReaderValue<Bool_t> Flag_ecalBadCalibFilter(myReader, "Flag_ecalBadCalibFilter");



    Long64_t eventLoopMax = myReader.GetEntries();
    int iev = 0;
    int iLHEPassed = 0;

    //BEGIN EVENT LOOP
    while (myReader.Next()) {
        if(++iev % 100000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<"\t"<<*genWeight<<endl;
        //if(iev > 40000) break;
        //if(iev>5) break;
        float eventScale = *genWeight;

        Xsec_genWeight = 1000.0 * lumi * xsec / sum_weights;

        RawEventNumber_TotalEvents++;
        WeightedEventNumber_TotalEvents+=eventScale;
        test_total +=eventScale;

        EventType_Leptonic = -1.0;
        EventType_SemiLeptonic = -1.0;

        CandidateHiggs_FatJet_idx = -1.0;
        CandidateHiggs_FatJet_particlenetScore = -1.0;
        CandidateHiggs_FatJet_pt = -1.0;
        CandidateHiggs_FatJet_eta = -1.0;
        CandidateHiggs_FatJet_phi = -1.0;
        CandidateHiggs_FatJet_mass = -1.0;
        CandidateHiggs_FatJet_msoftdrop = -1.0;

        CandidateVBF_Lead_Jet_idx = -1.0;
        CandidateVBF_Lead_Jet_pt = -1.0;
        CandidateVBF_Lead_Jet_eta = -1.0;
        CandidateVBF_Lead_Jet_phi = -1.0;

        CandidateVBF_Trail_Jet_idx = -1.0;
        CandidateVBF_Trail_Jet_pt = -1.0;
        CandidateVBF_Trail_Jet_eta = -1.0;
        CandidateVBF_Trail_Jet_phi = -1.0;

        CandidateVBF_Jet_invMass = -1.0;
        CandidateVBF_Jet_etaSep = -1.0;

        CandidateW_SemiLeptonic_FatJet_idx = -1.0;
        CandidateW_SemiLeptonic_FatJet_particlenetScore = -1.0;
        CandidateW_SemiLeptonic_FatJet_pt = -1.0;
        CandidateW_SemiLeptonic_FatJet_eta = -1.0;
        CandidateW_SemiLeptonic_FatJet_phi = -1.0;
        CandidateW_SemiLeptonic_FatJet_mass = -1.0;
        CandidateW_SemiLeptonic_FatJet_msoftdrop = -1.0;

        Candidate_SemiLeptonic_Lepton_idx = -1.0;
        Candidate_SemiLeptonic_Lepton_type = -1.0;
        Candidate_SemiLeptonic_Lepton_pt = -1.0;
        Candidate_SemiLeptonic_Lepton_eta = -1.0;
        Candidate_SemiLeptonic_Lepton_phi = -1.0;

        Candidate_Leptonic_Lead_Lepton_idx = -1.0;
        Candidate_Leptonic_Lead_Lepton_type = -1.0;
        Candidate_Leptonic_Lead_Lepton_pt = -1.0;
        Candidate_Leptonic_Lead_Lepton_eta = -1.0;
        Candidate_Leptonic_Lead_Lepton_phi = -1.0;

        Candidate_Leptonic_Trail_Lepton_idx = -1.0;
        Candidate_Leptonic_Trail_Lepton_type = -1.0;
        Candidate_Leptonic_Trail_Lepton_pt = -1.0;
        Candidate_Leptonic_Trail_Lepton_eta = -1.0;
        Candidate_Leptonic_Trail_Lepton_phi = -1.0;

        Candidate_Leptonic_Lepton_InvMass = -1.0;

        
        //if(++RawEventNumber_TotalEvents > 10000000) continue;

        int nLeptons_base = *nMuon + *nElectron;
        int nFatJet_base = *nFatJet;
        int nJet_base = *nJet;


        int nGoodElectrons = 0, nGoodMuons = 0, nGoodLeptons = 0;
        int nVetoElectrons = 0, nVetoMuons = 0, nVetoLeptons = 0;
        vector<int> vec_GoodElectrons, vec_GoodMuons;
        vector<int> vec_VetoElectrons, vec_VetoMuons;
        //cout<<"\n\nEVENT "<<iev<<endl;
        for(int i = 0; i < *nElectron; i++){
            //cout<<"Electron: "<<i;
            bool eleIso_pt = Electron_pt[i] > 10;
            bool eleIso_pt_veto = Electron_pt[i] > 7;
            bool eleIso_eta = fabs(Electron_eta[i] + Electron_deltaEtaSC[i]) < 2.5;
            bool eleIso_dxy = fabs(Electron_dxy[i]) < 0.05;
            bool eleIso_dz = fabs(Electron_dz[i]) < 0.1;
            bool eleIso_sip3d = fabs(Electron_sip3d[i]) < 8; 
            bool eleIso_Ie = Electron_miniPFRelIso_all[i] < (0.4);

            bool eleIso_sieie = false;
            if(fabs(Electron_eta[i]) < 1.479) eleIso_sieie = Electron_sieie[i] < 0.011;
            else eleIso_sieie = Electron_sieie[i] < 0.03;

            bool eleIso_hoe = Electron_hoe[i] < 0.1;
            bool eleIso_eInvMinusPInv = Electron_eInvMinusPInv[i] > -0.04;
            bool eleIso_convVeto = Electron_convVeto[i];
            bool eleIso_losthits = int(Electron_lostHits[i]) == 0;
            bool eleIso_losthits_veto = int(Electron_lostHits[i]) <= 1;
            bool eleIso_mvaFall17V2Iso_WPL = Electron_mvaFall17V2Iso_WPL[i];
            bool eleIso_mvaTTH = Electron_mvaTTH[i] > 0.8;

            bool eleIso_mvaNearbyJet = true;
            double minDR = 999.0;
            int nearJet = Electron_jetIdx[i];//-1;
            if(nearJet >= 0) { 
                //eleIso_mvaNearbyJet = Jet_btagDeepB[nearJet] < 0.4168;
                eleIso_mvaNearbyJet = Jet_btagDeepFlavB[nearJet] < deepFlavJetTTHCut;
            
            }
            //cout<<"\t"<<eleIso_pt <<"\t"<< eleIso_eta <<"\t"<< eleIso_dxy <<"\t"<< eleIso_dz <<"\t"<< 
            //   eleIso_sip3d <<"\t"<< eleIso_Ie <<"\t"<< eleIso_sieie <<"\t"<< eleIso_hoe <<"\t"<< 
            //   eleIso_eInvMinusPInv <<"\t"<< eleIso_convVeto <<"\t"<< eleIso_losthits <<"\t"<< 
            //   eleIso_mvaFall17V2Iso_WPL <<"\t"<< eleIso_mvaTTH <<"\t"<< eleIso_mvaNearbyJet<<endl;

            if(eleIso_pt && eleIso_eta && eleIso_dxy && eleIso_dz && 
               eleIso_sip3d && eleIso_Ie && eleIso_sieie && eleIso_hoe && 
               eleIso_eInvMinusPInv && eleIso_convVeto && eleIso_losthits && 
               eleIso_mvaFall17V2Iso_WPL && eleIso_mvaTTH && eleIso_mvaNearbyJet){
                   nGoodElectrons++;
                   vec_GoodElectrons.push_back(i);
               }
            if(eleIso_pt_veto && eleIso_eta && eleIso_dxy && eleIso_dz && 
               eleIso_sip3d && eleIso_Ie && eleIso_losthits_veto && 
               eleIso_mvaFall17V2Iso_WPL){
                   nVetoElectrons++;
                   vec_VetoElectrons.push_back(i);
               }
        }
        for(int i = 0; i < *nMuon; i++){
            bool muonIso_pt = Muon_pt[i] > 10;
            bool muonIso_pt_veto = Muon_pt[i] > 5;
            bool muonIso_eta = fabs(Muon_eta[i]) < 2.4;
            bool muonIso_dxy = fabs(Muon_dxy[i]) < 0.05;
            bool muonIso_dz = fabs(Muon_dz[i]) < 0.1;
            bool muonIso_sip3d = fabs(Muon_sip3d[i]) < 8;
            bool muonIso_Iu = Muon_miniPFRelIso_all[i] < (0.4);
            bool muonIso_muonMVA = Muon_mediumId[i];
            bool muonIso_muonMVA_veto = Muon_looseId[i];
            bool muonIso_mvaTTH = Muon_mvaTTH[i] > 0.85;

            bool muonIso_mvaNearbyJet = true;
            double minDR = 999.0;
            int nearJet = Muon_jetIdx[i];//-1;
            if(nearJet >= 0){
                //muonIso_mvaNearbyJet = Jet_btagDeepB[nearJet] < 0.4168;
                muonIso_mvaNearbyJet = Jet_btagDeepFlavB[nearJet] < deepFlavJetTTHCut;
            }
            if(muonIso_pt && muonIso_eta && muonIso_dxy && muonIso_dz && muonIso_sip3d && muonIso_Iu && muonIso_muonMVA && muonIso_mvaTTH && muonIso_mvaNearbyJet){
                nGoodMuons++;
                vec_GoodMuons.push_back(i);
            }
            if(muonIso_pt_veto && muonIso_eta && muonIso_dxy && muonIso_dz && muonIso_sip3d && muonIso_Iu && muonIso_muonMVA_veto){
                nVetoMuons++;
                vec_VetoMuons.push_back(i);
            }
        }
        nGoodLeptons = nGoodElectrons + nGoodMuons;
        nVetoLeptons = nVetoElectrons + nVetoMuons;

        int nPreSelectJets = 0;
        int nPreSelectFatJets =0;
        for(int iJet=0; iJet<*nJet; iJet++){
            if(Jet_pt[iJet] <= 20) continue;
            bool lepOverlap = false;
            for(int iM = 0; iM < nVetoMuons; iM++){
                int mu_idx = vec_VetoMuons.at(iM);
                double dR_mu = getDR(Muon_phi[mu_idx], Muon_eta[mu_idx], Jet_phi[iJet], Jet_eta[iJet]);
                if(dR_mu < 0.4) lepOverlap = true;
            }
            for(int iE = 0; iE < nVetoElectrons; iE++){
                int ele_idx = vec_VetoElectrons.at(iE);
                double dR_ele = getDR(Electron_phi[ele_idx], Electron_eta[ele_idx], Jet_phi[iJet], Jet_eta[iJet]);
                if(dR_ele < 0.4) lepOverlap = true;
            }
            if(!lepOverlap) nPreSelectJets++;
        }
        if(nPreSelectJets < 1) continue;
        RawEventNumber_PreSelect_Jets++;
        WeightedEventNumber_PreSelect_Jets+=eventScale;


        for(int iJet=0; iJet<*nFatJet; iJet++){
            if(FatJet_pt[iJet] <= 200) continue;
            if(FatJet_mass[iJet] <= 10) continue;
            if(FatJet_msoftdrop[iJet] <= 10) continue;
            bool lepOverlap = false;
            for(int iM = 0; iM < nVetoMuons; iM++){
                int mu_idx = vec_VetoMuons.at(iM);
                double dR_mu = getDR(Muon_phi[mu_idx], Muon_eta[mu_idx], FatJet_phi[iJet], FatJet_eta[iJet]);
                if(dR_mu < 0.8) lepOverlap = true;
            }
            for(int iE = 0; iE < nVetoElectrons; iE++){
                int ele_idx = vec_VetoElectrons.at(iE);
                double dR_ele = getDR(Electron_phi[ele_idx], Electron_eta[ele_idx], FatJet_phi[iJet], FatJet_eta[iJet]);
                if(dR_ele < 0.8) lepOverlap = true;
            }
            if(!lepOverlap) nPreSelectFatJets++;
        }
        if(nPreSelectFatJets < 1) continue;
        RawEventNumber_PreSelect_FatJets++;
        WeightedEventNumber_PreSelect_FatJets+=eventScale;

        bool base_leptonic_channel = false, base_semileptonic_channel = false;

        if(nVetoLeptons >= 2) base_leptonic_channel = true;
        if(nVetoLeptons >= 1) base_semileptonic_channel = true;

        if(!base_leptonic_channel && !base_semileptonic_channel) continue;

        if(base_leptonic_channel){
            RawEventNumber_PreSelect_LeptonicChannel++;
            WeightedEventNumber_PreSelect_LeptonicChannel+=eventScale;
        }
        if(base_semileptonic_channel){
            //cout<<"HEREE"<<endl;
            RawEventNumber_PreSelect_SemiLeptonicChannel++;
            WeightedEventNumber_PreSelect_SemiLeptonicChannel+=eventScale;
            test_sem+=eventScale;
        }
        //leptonic channels: 
        // Wm -> e + nu_e_bar, Wm -> mu + nu_mu_bar
        // Wp -> e_bar + nu_e, Wp -> mu_bar + nu_mu

        bool pass_singleEle = false, 
             pass_singleMuon = false, 
             pass_doubleEle = false, 
             pass_doubleMuon = false, 
             pass_muonEG = false,
             pass_jet = false,
             pass_hbb = false,
             pass_ucsd_l = false,
             pass_ucsd_sl = false,
             pass_ucsd_base = false;

        //my trigger
        if(*HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 || *HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2 || *HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) {
            
            pass_hbb = true;
            }

        if(*HLT_Ele115_CaloIdVT_GsfTrkIdT || *HLT_Ele27_WPTight_Gsf || *HLT_Ele28_WPTight_Gsf || 
           *HLT_Ele32_WPTight_Gsf || *HLT_Ele35_WPTight_Gsf || *HLT_Ele38_WPTight_Gsf ||
           *HLT_Ele40_WPTight_Gsf || *HLT_Ele32_WPTight_Gsf_L1DoubleEG || *HLT_Photon200) {
            
            pass_singleEle = true;
            }
        if(*HLT_IsoMu24 || *HLT_IsoMu27 || *HLT_IsoMu30 || *HLT_Mu50) {
            pass_singleMuon = true;
            }
        if(*HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || *HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL|| *HLT_DiEle27_WPTightCaloOnly_L1DoubleEG || 
           *HLT_DoubleEle33_CaloIdL_MW || *HLT_DoubleEle25_CaloIdL_MW || *HLT_DoubleEle27_CaloIdL_MW || *HLT_DoublePhoton70) {
            pass_doubleEle = true;
            }
        if(*HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || *HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || 
           *HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || *HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8) {
            pass_doubleMuon = true;
            }
        if(*HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || *HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL|| 
           *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || 
           *HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ|| *HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL) {
            pass_muonEG = true;
            }
        if(*HLT_AK8PFJet500 || *HLT_AK8PFJet360_TrimMass30 || *HLT_AK8PFJet380_TrimMass30 ||
           *HLT_AK8PFJet400_TrimMass30 || *HLT_AK8PFJet420_TrimMass30 || *HLT_AK8PFHT750_TrimMass50 ||
           *HLT_AK8PFHT800_TrimMass50 || *HLT_AK8PFHT850_TrimMass50 || *HLT_AK8PFHT900_TrimMass50 ||
           *HLT_PFHT1050){
               pass_jet = true;
           }

        

        if(*Flag_goodVertices
        && *Flag_HBHENoiseFilter
        && *Flag_HBHENoiseIsoFilter
        && *Flag_EcalDeadCellTriggerPrimitiveFilter
        && *Flag_BadPFMuonFilter)
        pass_ucsd_base = true;

        if(year == "17" || year == "18") {
            if(!*Flag_ecalBadCalibFilter) pass_ucsd_base = false;
        }

        if(!pass_ucsd_base) continue;

        bool continue_leptonic = false, continue_semileptonic = false;
        if(base_leptonic_channel){
            continue_leptonic = true;
            RawEventNumber_Trigger_LeptonicChannel++;
            WeightedEventNumber_Trigger_LeptonicChannel+=eventScale;
        }
        if(base_semileptonic_channel){
            continue_semileptonic = true;
            RawEventNumber_Trigger_SemiLeptonicChannel++;
            WeightedEventNumber_Trigger_SemiLeptonicChannel+=eventScale;
        }
        


                
        // --------------------------- EVENT SELECTION ---------------------------- //

        double invMass = -1.0;
        //W decay channels: leptonic, semileptonic, hadronic
        //various cuts:

        //vector< pair< pair< int, int >, float> > 
        std::vector< std::tuple< int, int, float, int, float, float, float > > lep_candidates_leptonic; // type, index, pt, charge, eta, phi, mass -> for easy sorting
        std::vector< std::tuple< int, int, float, int, float, float, float > > lep_candidates_semileptonic; // type, index, pt, charge, eta, phi, mass -> for easy sorting

        bool passed_leptonic = false, passed_semileptonic = false;
        bool LeptonSelection_LChannel_Pass = false;
        bool LeptonSelection_SLChannel_Pass=false;
        bool FatJetSelection_SLChannel_Pass=false;
        bool VBFSelection_LChannel_Pass = false;
        bool VBFSelection_SLChannel_Pass=false;

        int tmp_semileptonic_FatJet_idx = -1;
        std::tuple<int, int, float, int, float, float, float> lep_candidate_leptonic_lead_final;
        std::tuple<int, int, float, int, float, float, float> lep_candidate_leptonic_trail_final;
        std::tuple<int, int, float, int, float, float, float> lep_candidate_semileptonic_final;

        float lep_leptonic_mll_final;

        //select leptons
        //LEPTONIC
        //if(channel_leptonic){
        if(base_leptonic_channel){
            bool METCut_Leptonic_passSelection = (*MET_pt > 100);//300);

            nEvents_cuts_leptonicChannel[0]++; //Events entering leptonic channel consideration
            nEvents_weighted_cuts_leptonicChannel[0]+=eventScale; //Events entering leptonic channel consideration


            //find leading/trailing pt muons and electrons and taus
            //for(int iEle=0; iEle<*nElectron; iEle++){
            //    if(Electron_mvaFall17V2Iso_WP90[iEle] == false) continue;
            for(int ii=0; ii<vec_GoodElectrons.size(); ii++){
                int iEle = vec_GoodElectrons.at(ii);
                float sumPTCut = Electron_pt[iEle] * 0.1;
                if(Electron_pt[iEle] > 10 /*&& Electron_dr03TkSumPt[iEle] < sumPTCut && Electron_convVeto[iEle]*/) 
                    lep_candidates_leptonic.push_back(std::make_tuple(isEle, iEle, Electron_pt[iEle], Electron_charge[iEle], Electron_eta[iEle], Electron_phi[iEle], Electron_mass[iEle]));
            }
            //for(int iMuon=0; iMuon<*nMuon; iMuon++){
            //    if(Muon_tightId[iMuon] == false) continue;
            for(int ii=0; ii<vec_GoodMuons.size(); ii++){
                int iMuon = vec_GoodMuons.at(ii);
                float sumPTCut = Muon_pt[iMuon] * 0.1;
                if(Muon_pt[iMuon] > 10/* && Muon_tkRelIso[iMuon]<sumPTCut*/) 
                    lep_candidates_leptonic.push_back(std::make_tuple(isMuon, iMuon, Muon_pt[iMuon], Muon_charge[iMuon], Muon_eta[iMuon], Muon_phi[iMuon], Muon_mass[iMuon]));
            }
            sort(lep_candidates_leptonic.begin(), lep_candidates_leptonic.end(), TupleCompare<2>()); // sort candidates by pt

            std::vector< std::tuple< int, int, float, ROOT::Math::PtEtaPhiMVector> > lep_pairs;

            if(lep_candidates_leptonic.size() > 1){
                for(int iLead = 0; iLead < lep_candidates_leptonic.size(); iLead++){
                    if(get<2>(lep_candidates_leptonic.at(iLead)) < 20/*leptonPTCut_lead*/) continue;
                    for(int iTrail = iLead+1; iTrail < lep_candidates_leptonic.size(); iTrail++){
                        if((get<3>(lep_candidates_leptonic.at(iLead))*get<3>(lep_candidates_leptonic.at(iTrail))) >= 0) continue;

                        ROOT::Math::PtEtaPhiMVector leadLep_p4(get<2>(lep_candidates_leptonic.at(iLead)), get<4>(lep_candidates_leptonic.at(iLead)), get<5>(lep_candidates_leptonic.at(iLead)), get<6>(lep_candidates_leptonic.at(iLead)));
                        ROOT::Math::PtEtaPhiMVector trailLep_p4(get<2>(lep_candidates_leptonic.at(iTrail)), get<4>(lep_candidates_leptonic.at(iTrail)), get<5>(lep_candidates_leptonic.at(iTrail)), get<6>(lep_candidates_leptonic.at(iTrail)));
                        ROOT::Math::PtEtaPhiMVector lepSum = leadLep_p4+trailLep_p4;
                        double invMass_lep = sqrt(lepSum.Dot(lepSum));

                        if(invMass_lep < 93 && invMass_lep > 89) continue;

                        lep_pairs.push_back(std::make_tuple(iLead, iTrail, invMass_lep, lepSum));
            
                    }
                }
            }
            if(lep_pairs.size() > 0){
                sort(lep_pairs.begin(), lep_pairs.end(), TupleCompare<2>()); // sort by pt
                lep_candidate_leptonic_lead_final = lep_candidates_leptonic.at(get<0>(lep_pairs.at(0)));
                lep_candidate_leptonic_trail_final = lep_candidates_leptonic.at(get<1>(lep_pairs.at(0)));

                ROOT::Math::PtEtaPhiMVector leadLep_p4(get<2>(lep_candidates_leptonic.at(get<0>(lep_pairs.at(0)))), get<4>(lep_candidates_leptonic.at(get<0>(lep_pairs.at(0)))), get<5>(lep_candidates_leptonic.at(get<0>(lep_pairs.at(0)))), get<6>(lep_candidates_leptonic.at(get<0>(lep_pairs.at(0)))));
                ROOT::Math::PtEtaPhiMVector trailLep_p4(get<2>(lep_candidates_leptonic.at(get<1>(lep_pairs.at(0)))), get<4>(lep_candidates_leptonic.at(get<1>(lep_pairs.at(0)))), get<5>(lep_candidates_leptonic.at(get<1>(lep_pairs.at(0)))), get<6>(lep_candidates_leptonic.at(get<1>(lep_pairs.at(0)))));
                ROOT::Math::PtEtaPhiMVector lepSum = leadLep_p4+trailLep_p4;

                lep_leptonic_mll_final = get<2>(lep_pairs.at(0));
                //if(invMass_lep < 15) continue;

                //Candidate_Leptonic_MT = sqrt(2.0*(lepSum.Pt())*(*PuppiMET_pt) * (1 - cos((lepSum.Phi()) - (*PuppiMET_phi))));


                LeptonSelection_LChannel_Pass = true;
                RawEventNumber_Leptons_LeptonicChannel++;
                WeightedEventNumber_Leptons_LeptonicChannel+=eventScale;
                
            }
            else{
                continue_leptonic = false;
            }
        }

        //SEMI LEPTONIC
        if(base_semileptonic_channel){

            //for(int iEle=0; iEle<*nElectron; iEle++){
            //    if(Electron_mvaFall17V2Iso_WP90[iEle] == false) continue;
            for(int ii=0; ii<vec_GoodElectrons.size(); ii++){
                int iEle = vec_GoodElectrons.at(ii);
                float sumPTCut = Electron_pt[iEle] * 0.1;
                
                if(Electron_pt[iEle] > 10/*leptonPTCut_lead && Electron_dr03TkSumPt[iEle] < sumPTCut && Electron_convVeto[iEle]*/)
                    lep_candidates_semileptonic.push_back(std::make_tuple(isEle, iEle, Electron_pt[iEle], Electron_charge[iEle], Electron_eta[iEle], Electron_phi[iEle], Electron_mass[iEle]));
            }
            //for(int iMuon=0; iMuon<*nMuon; iMuon++){
            //    if(Muon_tightId[iMuon] == false) continue;
            for(int ii=0; ii<vec_GoodMuons.size(); ii++){
                int iMuon = vec_GoodMuons.at(ii);
                float sumPTCut = Muon_pt[iMuon] * 0.1;
                if(Muon_pt[iMuon] > 10 /*&& Muon_tkRelIso[iMuon]<sumPTCut*/) 
                    lep_candidates_semileptonic.push_back(std::make_tuple(isMuon, iMuon, Muon_pt[iMuon], Muon_charge[iMuon], Muon_eta[iMuon], Muon_phi[iMuon], Muon_mass[iMuon]));
            }
            //sort by PT
            sort(lep_candidates_semileptonic.begin(), lep_candidates_semileptonic.end(), TupleCompare<2>());
            if(lep_candidates_semileptonic.size() > 0){
                lep_candidate_semileptonic_final = lep_candidates_semileptonic.at(0);
                
                LeptonSelection_SLChannel_Pass=true;
                nEvents_cuts_semileptonicChannel[3]++; //Passed Lead Lepton PT Cut
                nEvents_weighted_cuts_semileptonicChannel[3]+=eventScale; //Passed Lead Lepton PT Cut

                RawEventNumber_Leptons_SemiLeptonicChannel++;
                WeightedEventNumber_Leptons_SemiLeptonicChannel+=eventScale;
            }
            else{
                continue_semileptonic = false;
            }
            
        }

        //Higgs FatJets
        //Just grabbing highest pt FatJet for now
        vector< std::tuple< int, float, float, int > > vec_higgsFatJet_candidates;
        //idx, pt, particleNetScore, jetID
        for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
            bool lepOverlap = false;
            for(int iM = 0; iM < nVetoMuons; iM++){
                int mu_idx = vec_VetoMuons.at(iM);
                double dR_mu = getDR(Muon_phi[mu_idx], Muon_eta[mu_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                if(dR_mu < 0.8) lepOverlap = true;
            }
            for(int iE = 0; iE < nVetoElectrons; iE++){
                int ele_idx = vec_VetoElectrons.at(iE);
                double dR_ele = getDR(Electron_phi[ele_idx], Electron_eta[ele_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                if(dR_ele < 0.8) lepOverlap = true;
            }
            if(!lepOverlap)
                vec_higgsFatJet_candidates.push_back(make_tuple(iFatJet,FatJet_pt[iFatJet], FatJet_particleNet_HbbvsQCD[iFatJet],FatJet_jetId[iFatJet] )); //sorting by deepTag Scores
        }
        vec_higgsFatJet_candidates.erase(std::remove_if( vec_higgsFatJet_candidates.begin(), vec_higgsFatJet_candidates.end(),
            [](const tuple< int, float, float, int >& x) -> bool{ 
                return get<1>(x) <= fatJetPTCut_Higgs; // put your condition here
            }), vec_higgsFatJet_candidates.end());

        vec_higgsFatJet_candidates.erase(std::remove_if( vec_higgsFatJet_candidates.begin(), vec_higgsFatJet_candidates.end(),
            [](const tuple< int, float, float, int >& x) -> bool{ 
                return get<3>(x) != 6; // put your condition here
            }), vec_higgsFatJet_candidates.end());

        vec_higgsFatJet_candidates.erase(std::remove_if( vec_higgsFatJet_candidates.begin(), vec_higgsFatJet_candidates.end(),
            [](const tuple< int, float, float, int >& x) -> bool{ 
                return get<2>(x) <= fatJetParticleNetCut_Higgs_cut; // put your condition here
            }), vec_higgsFatJet_candidates.end());

        if(vec_higgsFatJet_candidates.size() >= 1){
            sort(vec_higgsFatJet_candidates.begin(), vec_higgsFatJet_candidates.end(), TupleCompare<2>());
            CandidateHiggs_FatJet_idx = get<0>(vec_higgsFatJet_candidates.at(0));
        }
        if(CandidateHiggs_FatJet_idx < 0) continue;
        if(continue_leptonic){
            RawEventNumber_HiggsFatJet_LeptonicChannel++;
            WeightedEventNumber_HiggsFatJet_LeptonicChannel+=eventScale;
        }
        if(continue_semileptonic){
            RawEventNumber_HiggsFatJet_SemiLeptonicChannel++;
            WeightedEventNumber_HiggsFatJet_SemiLeptonicChannel+=eventScale;
        }


        int nBJet_loose_afterH = 0,
            nBJet_tight_afterH = 0;
        for(int iJet=0; iJet<*nJet; iJet++){
            double dR_bwithH = getDR(Jet_phi[iJet], Jet_eta[iJet], FatJet_phi[CandidateHiggs_FatJet_idx], FatJet_eta[CandidateHiggs_FatJet_idx]);
            if(dR_bwithH <= 0.8 || Jet_pt[iJet] < 30 || fabs(Jet_eta[iJet]) > 2.5) continue;

            if(Jet_btagDeepB[iJet] >= 0.4168) nBJet_tight_afterH++; //tight = 0.7100
        }
        if(nBJet_tight_afterH > 0) continue;
        if(continue_leptonic){
            RawEventNumber_BJetVeto_LeptonicChannel++;
            WeightedEventNumber_BJetVeto_LeptonicChannel+=eventScale;
        }
        if(continue_semileptonic){
            RawEventNumber_BJetVeto_SemiLeptonicChannel++;
            WeightedEventNumber_BJetVeto_SemiLeptonicChannel+=eventScale;
        }


        //SEMI LEPTONIC
        if(continue_semileptonic){
            vector< std::tuple< int, float > > fatJet_passedPT;
            for(int iFatJet=0; iFatJet<*nFatJet; iFatJet++){
                if(iFatJet == CandidateHiggs_FatJet_idx) continue;
                float dR_vbfJet_1 = -999, dR_vbfJet_2 = -999, dR_Higgs = -999;
                dR_Higgs = getDR(FatJet_phi[CandidateHiggs_FatJet_idx], FatJet_eta[CandidateHiggs_FatJet_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                        
                bool lepOverlap = false;
                for(int iM = 0; iM < nVetoMuons; iM++){
                    int mu_idx = vec_VetoMuons.at(iM);
                    double dR_mu = getDR(Muon_phi[mu_idx], Muon_eta[mu_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                    if(dR_mu < 0.8) lepOverlap = true;
                }
                for(int iE = 0; iE < nVetoElectrons; iE++){
                    int ele_idx = vec_VetoElectrons.at(iE);
                    double dR_ele = getDR(Electron_phi[ele_idx], Electron_eta[ele_idx], FatJet_phi[iFatJet], FatJet_eta[iFatJet]);
                    if(dR_ele < 0.8) lepOverlap = true;
                }
                if( FatJet_pt[iFatJet] > fatJetPTCut 
                    && dR_Higgs > 0.8 
                    && !lepOverlap
                    && FatJet_jetId[iFatJet] == 6 
                    //&& FatJet_deepTag_WvsQCD[iFatJet] > deepTag_WvsQCD[level_N])
                    && (FatJet_particleNet_WvsQCD[iFatJet] > fatJetParticleNetCut_WSL_cut || FatJet_particleNet_ZvsQCD[iFatJet] > fatJetParticleNetCut_WSL_cut )
                )
                    fatJet_passedPT.push_back(make_tuple(iFatJet, FatJet_pt[iFatJet]));
            }
            sort(fatJet_passedPT.begin(), fatJet_passedPT.end(), TupleCompare<1>());
            if(fatJet_passedPT.size() > 0){ 
                tmp_semileptonic_FatJet_idx = get<0>(fatJet_passedPT.at(0));
                FatJetSelection_SLChannel_Pass = true;
                RawEventNumber_WFatJet_SemiLeptonicChannel++;
                WeightedEventNumber_WFatJet_SemiLeptonicChannel+=eventScale;
            }
            else{
                continue_semileptonic = false;
            }
        }

        //Start New Jet
        vector<std::tuple<int,float>> jets_pos_l;
        vector<std::tuple<int,float>> jets_neg_l;
        vector<std::tuple<int,float>> jets_pos_sl;
        vector<std::tuple<int,float>> jets_neg_sl;

        for(int iJet_1=0; iJet_1<*nJet; iJet_1++){
            if(Jet_pt[iJet_1] < jetPT_lead_cut || Jet_jetId[iJet_1] != 6) continue;
            double dR_vbfJet_Higgs = getDR(Jet_phi[iJet_1], Jet_eta[iJet_1], FatJet_phi[CandidateHiggs_FatJet_idx], FatJet_eta[CandidateHiggs_FatJet_idx]);
            if(dR_vbfJet_Higgs < 0.8) continue;

            bool lepOverlap = false;
            for(int iM = 0; iM < nVetoMuons; iM++){
                int mu_idx = vec_VetoMuons.at(iM);
                double dR_mu = getDR(Muon_phi[mu_idx], Muon_eta[mu_idx], Jet_phi[iJet_1], Jet_eta[iJet_1]);
                if(dR_mu < 0.4) lepOverlap = true;
            }
            for(int iE = 0; iE < nVetoElectrons; iE++){
                int ele_idx = vec_VetoElectrons.at(iE);
                double dR_ele = getDR(Electron_phi[ele_idx], Electron_eta[ele_idx], Jet_phi[iJet_1], Jet_eta[iJet_1]);
                if(dR_ele < 0.4) lepOverlap = true;
            }
            if(lepOverlap) continue;

            ROOT::Math::PtEtaPhiMVector iJet_1_p4(Jet_pt[iJet_1], Jet_eta[iJet_1], Jet_phi[iJet_1], Jet_mass[iJet_1]);

            if(Jet_eta[iJet_1] >= 0) jets_pos_l.push_back(make_tuple(iJet_1, iJet_1_p4.P()));
            else jets_neg_l.push_back(make_tuple(iJet_1, iJet_1_p4.P()));

            if(continue_semileptonic){
                double dR_vbfJet_W = getDR(Jet_phi[iJet_1], Jet_eta[iJet_1], FatJet_phi[tmp_semileptonic_FatJet_idx], FatJet_eta[tmp_semileptonic_FatJet_idx]);
                if(dR_vbfJet_W < 0.8) continue;

                if(Jet_eta[iJet_1] >= 0) jets_pos_sl.push_back(make_tuple(iJet_1, iJet_1_p4.P()));
                else jets_neg_sl.push_back(make_tuple(iJet_1, iJet_1_p4.P()));

            }
        }
        int tmp_jet_lead_L = -1, tmp_jet_trail_L = -1,tmp_jet_lead_SL = -1, tmp_jet_trail_SL = -1;
        if(continue_leptonic){
            int nJetCands = jets_pos_l.size() + jets_neg_l.size();
            if(nJetCands >= 2){

                sort(jets_pos_l.begin(), jets_pos_l.end(), TupleCompare<1>());
                sort(jets_neg_l.begin(), jets_neg_l.end(), TupleCompare<1>());

                if(jets_pos_l.size() == 0){
                    tmp_jet_lead_L = get<0>(jets_neg_l.at(0));
                    tmp_jet_trail_L = get<0>(jets_neg_l.at(1));
                }
                else if(jets_neg_l.size() == 0){
                    tmp_jet_lead_L = get<0>(jets_pos_l.at(0));
                    tmp_jet_trail_L = get<0>(jets_pos_l.at(1));
                }
                else{
                    tmp_jet_lead_L = get<0>(jets_pos_l.at(0));
                    tmp_jet_trail_L = get<0>(jets_neg_l.at(0));
                }
                if(Jet_pt[tmp_jet_trail_L] > Jet_pt[tmp_jet_lead_L]) swap(tmp_jet_trail_L,tmp_jet_lead_L);

                ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[tmp_jet_lead_L], Jet_eta[tmp_jet_lead_L], Jet_phi[tmp_jet_lead_L], Jet_mass[tmp_jet_lead_L]);
                ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[tmp_jet_trail_L], Jet_eta[tmp_jet_trail_L], Jet_phi[tmp_jet_trail_L], Jet_mass[tmp_jet_trail_L]);
                ROOT::Math::PtEtaPhiMVector jetSum = leadJet_p4+trailJet_p4;
                invMass = sqrt(jetSum.Dot(jetSum));
                if(invMass >= jetInvMass_cut && fabs(Jet_eta[tmp_jet_lead_L] - Jet_eta[tmp_jet_trail_L]) >= jetEtaSep_cut){
                    RawEventNumber_VBFJets_LeptonicChannel++;
                    WeightedEventNumber_VBFJets_LeptonicChannel+=eventScale;
                    VBFSelection_LChannel_Pass = true;

                }
                else{ 
                    continue_leptonic = false;
                }
            }
            else{ 
                continue_leptonic = false;
            }
        }
        if(continue_semileptonic){
            int nJetCands = jets_pos_sl.size() + jets_neg_sl.size();
            if(nJetCands >= 2){ 

                sort(jets_pos_sl.begin(), jets_pos_sl.end(), TupleCompare<1>());
                sort(jets_neg_sl.begin(), jets_neg_sl.end(), TupleCompare<1>());


                if(jets_pos_sl.size() == 0){
                    tmp_jet_lead_SL = get<0>(jets_neg_sl.at(0));
                    tmp_jet_trail_SL = get<0>(jets_neg_sl.at(1));
                }
                else if(jets_neg_sl.size() == 0){
                    tmp_jet_lead_SL = get<0>(jets_pos_sl.at(0));
                    tmp_jet_trail_SL = get<0>(jets_pos_sl.at(1));
                }
                else{
                    tmp_jet_lead_SL = get<0>(jets_pos_sl.at(0));
                    tmp_jet_trail_SL = get<0>(jets_neg_sl.at(0));
                }

                if(Jet_pt[tmp_jet_trail_SL] > Jet_pt[tmp_jet_lead_SL]) swap(tmp_jet_trail_SL,tmp_jet_lead_SL);

                ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[tmp_jet_lead_SL], Jet_eta[tmp_jet_lead_SL], Jet_phi[tmp_jet_lead_SL], Jet_mass[tmp_jet_lead_SL]);
                ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[tmp_jet_trail_SL], Jet_eta[tmp_jet_trail_SL], Jet_phi[tmp_jet_trail_SL], Jet_mass[tmp_jet_trail_SL]);
                ROOT::Math::PtEtaPhiMVector jetSum = leadJet_p4+trailJet_p4;
                invMass = sqrt(jetSum.Dot(jetSum));
                if(invMass >= jetInvMass_cut && fabs(Jet_eta[tmp_jet_lead_SL] - Jet_eta[tmp_jet_trail_SL]) >= jetEtaSep_cut){
                    RawEventNumber_VBFJets_SemiLeptonicChannel++;
                    WeightedEventNumber_VBFJets_SemiLeptonicChannel+=eventScale;
                    VBFSelection_SLChannel_Pass = true;

                }
                else{ 
                    continue_semileptonic = false;
                }
            }
            else{ 
                continue_semileptonic = false;
            }
        }
        if(!VBFSelection_LChannel_Pass && !VBFSelection_SLChannel_Pass) continue;
        //cout<<"end vbf"<<endl;

        //End New Jet

        //st cuts
        double lt_l = -999.0;
        double st_l = -999.0;
        double lt_s = -999.0;
        double st_s = -999.0;
        double rpt_s = -999.0;
        double rpt_l = -999.0;
        double rpt_nomet_s = -999.0;
        double rpt_nomet_l = -999.0;
        if(continue_leptonic){
            ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[CandidateVBF_Lead_Jet_idx], Jet_eta[CandidateVBF_Lead_Jet_idx], Jet_phi[CandidateVBF_Lead_Jet_idx], Jet_mass[CandidateVBF_Lead_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[CandidateVBF_Trail_Jet_idx], Jet_eta[CandidateVBF_Trail_Jet_idx], Jet_phi[CandidateVBF_Trail_Jet_idx], Jet_mass[CandidateVBF_Trail_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector higgs_p4(FatJet_pt[CandidateHiggs_FatJet_idx], FatJet_eta[CandidateHiggs_FatJet_idx], FatJet_phi[CandidateHiggs_FatJet_idx], FatJet_mass[CandidateHiggs_FatJet_idx]);
            ROOT::Math::PtEtaPhiMVector lep1_p4(get<2>(lep_candidate_leptonic_lead_final), get<4>(lep_candidate_leptonic_lead_final), get<5>(lep_candidate_leptonic_lead_final), get<6>(lep_candidate_leptonic_lead_final));
            ROOT::Math::PtEtaPhiMVector lep2_p4(get<2>(lep_candidate_leptonic_trail_final), get<4>(lep_candidate_leptonic_trail_final), get<5>(lep_candidate_leptonic_trail_final), get<6>(lep_candidate_leptonic_trail_final));
            ROOT::Math::PtEtaPhiEVector met_p4(*MET_pt, 0.0, *MET_phi, *MET_sumEt);
                
            ROOT::Math::PtEtaPhiMVector p4_rpt = leadJet_p4 + trailJet_p4 + higgs_p4 + lep1_p4 + lep2_p4 + met_p4;
            ROOT::Math::PtEtaPhiMVector p4_rpt_nomet = leadJet_p4 + trailJet_p4 + higgs_p4 + lep1_p4 + lep2_p4;

            float pt_sum = Jet_pt[CandidateVBF_Lead_Jet_idx] + Jet_pt[CandidateVBF_Trail_Jet_idx] + FatJet_pt[CandidateHiggs_FatJet_idx] + 
                            get<2>(lep_candidate_leptonic_lead_final) + get<2>(lep_candidate_leptonic_trail_final) + *MET_pt;

            rpt_l = p4_rpt.Pt() / pt_sum;

            lt_l = get<2>(lep_candidate_leptonic_lead_final) + get<2>(lep_candidate_leptonic_trail_final) + *MET_pt;
            st_l = FatJet_pt[CandidateHiggs_FatJet_idx] + lt_l;
            if(lt_l > lt_l_cut && st_l > st_l_cut){
                passed_leptonic = true; 
                nEvents_cuts_leptonicChannel[5]++; //Total Passed Events
                nEvents_weighted_cuts_leptonicChannel[5]+=eventScale; //Total Passed Events
            }
        }
        //SEMI LEPTONIC
        if(continue_semileptonic){
            ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[CandidateVBF_Lead_Jet_idx], Jet_eta[CandidateVBF_Lead_Jet_idx], Jet_phi[CandidateVBF_Lead_Jet_idx], Jet_mass[CandidateVBF_Lead_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[CandidateVBF_Trail_Jet_idx], Jet_eta[CandidateVBF_Trail_Jet_idx], Jet_phi[CandidateVBF_Trail_Jet_idx], Jet_mass[CandidateVBF_Trail_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector higgs_p4(FatJet_pt[CandidateHiggs_FatJet_idx], FatJet_eta[CandidateHiggs_FatJet_idx], FatJet_phi[CandidateHiggs_FatJet_idx], FatJet_mass[CandidateHiggs_FatJet_idx]);
            ROOT::Math::PtEtaPhiMVector lep1_p4(get<2>(lep_candidate_semileptonic_final), get<4>(lep_candidate_semileptonic_final), get<5>(lep_candidate_semileptonic_final), get<6>(lep_candidate_semileptonic_final));
            ROOT::Math::PtEtaPhiMVector wak8_p4(FatJet_pt[CandidateW_SemiLeptonic_FatJet_idx], FatJet_eta[CandidateW_SemiLeptonic_FatJet_idx], FatJet_phi[CandidateW_SemiLeptonic_FatJet_idx], FatJet_mass[CandidateW_SemiLeptonic_FatJet_idx]);
            ROOT::Math::PtEtaPhiEVector met_p4(*MET_pt, 0.0, *MET_phi, *MET_sumEt);
                
            ROOT::Math::PtEtaPhiMVector p4_rpt = leadJet_p4 + trailJet_p4 + higgs_p4 + lep1_p4 + wak8_p4 + met_p4;
            ROOT::Math::PtEtaPhiMVector p4_rpt_nomet = leadJet_p4 + trailJet_p4 + higgs_p4 + lep1_p4 + wak8_p4;

            float pt_sum = Jet_pt[CandidateVBF_Lead_Jet_idx] + Jet_pt[CandidateVBF_Trail_Jet_idx] + FatJet_pt[CandidateHiggs_FatJet_idx] + 
                            get<2>(lep_candidate_leptonic_lead_final) + FatJet_pt[CandidateW_SemiLeptonic_FatJet_idx] + *MET_pt;

            rpt_s = p4_rpt.Pt() / pt_sum;
            lt_s = get<2>(lep_candidate_semileptonic_final) + FatJet_pt[tmp_semileptonic_FatJet_idx]  + *MET_pt;
            st_s = FatJet_pt[CandidateHiggs_FatJet_idx] + lt_s;
            if(lt_s > lt_sl_cut && st_s > st_sl_cut){
                passed_semileptonic = true;

                nEvents_cuts_semileptonicChannel[4]++; //Final Passed Events
                nEvents_weighted_cuts_semileptonicChannel[4]+=eventScale; //Final Passed Events
            }
        }









        int passed_temp = (passed_leptonic ? 1:0) + (passed_semileptonic ? 1:0);
        vector<int> to_fill;

        if(passed_temp == 0){
            RawEventNumber_Rejected++;
            WeightedEventNumber_Rejected+=eventScale;
            continue;
        }
        if(passed_temp > 0){
            if(passed_leptonic){
                CandidateVBF_Lead_Jet_idx = tmp_jet_lead_L;
                CandidateVBF_Trail_Jet_idx = tmp_jet_trail_L;
            }
            else if(passed_semileptonic){
                CandidateVBF_Lead_Jet_idx = tmp_jet_lead_SL;
                CandidateVBF_Trail_Jet_idx = tmp_jet_trail_SL;
            }
        }


        Events->GetEntry(myReader.GetCurrentEntry());

        if(passed_temp > 0){
            CandidateHiggs_FatJet_pt = FatJet_pt[CandidateHiggs_FatJet_idx];
            CandidateHiggs_FatJet_eta = FatJet_eta[CandidateHiggs_FatJet_idx];
            CandidateHiggs_FatJet_phi = FatJet_phi[CandidateHiggs_FatJet_idx];
            CandidateHiggs_FatJet_particlenetScore = FatJet_particleNet_HbbvsQCD[CandidateHiggs_FatJet_idx];
            CandidateHiggs_FatJet_mass = FatJet_mass[CandidateHiggs_FatJet_idx];
            CandidateHiggs_FatJet_msoftdrop = FatJet_msoftdrop[CandidateHiggs_FatJet_idx];

            CandidateVBF_Lead_Jet_pt = Jet_pt[CandidateVBF_Lead_Jet_idx];
            CandidateVBF_Lead_Jet_eta = Jet_eta[CandidateVBF_Lead_Jet_idx];
            CandidateVBF_Lead_Jet_phi = Jet_phi[CandidateVBF_Lead_Jet_idx];

            CandidateVBF_Trail_Jet_pt = Jet_pt[CandidateVBF_Trail_Jet_idx];
            CandidateVBF_Trail_Jet_eta = Jet_eta[CandidateVBF_Trail_Jet_idx];
            CandidateVBF_Trail_Jet_phi = Jet_phi[CandidateVBF_Trail_Jet_idx];


            ROOT::Math::PtEtaPhiMVector leadJet_p4(Jet_pt[CandidateVBF_Lead_Jet_idx], Jet_eta[CandidateVBF_Lead_Jet_idx], Jet_phi[CandidateVBF_Lead_Jet_idx], Jet_mass[CandidateVBF_Lead_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector trailJet_p4(Jet_pt[CandidateVBF_Trail_Jet_idx], Jet_eta[CandidateVBF_Trail_Jet_idx], Jet_phi[CandidateVBF_Trail_Jet_idx], Jet_mass[CandidateVBF_Trail_Jet_idx]);
            ROOT::Math::PtEtaPhiMVector jetSum = leadJet_p4+trailJet_p4;
            invMass = sqrt(jetSum.Dot(jetSum));

            CandidateVBF_Jet_invMass = invMass;
            CandidateVBF_Jet_etaSep = fabs(Jet_eta[CandidateVBF_Lead_Jet_idx] - Jet_eta[CandidateVBF_Trail_Jet_idx]);

            if(passed_leptonic){

                Candidate_Leptonic_Lead_Lepton_idx = get<1>(lep_candidate_leptonic_lead_final);
                Candidate_Leptonic_Lead_Lepton_type = get<0>(lep_candidate_leptonic_lead_final);
                Candidate_Leptonic_Lead_Lepton_pt = get<2>(lep_candidate_leptonic_lead_final);
                Candidate_Leptonic_Lead_Lepton_eta = get<4>(lep_candidate_leptonic_lead_final);
                Candidate_Leptonic_Lead_Lepton_phi = get<5>(lep_candidate_leptonic_lead_final);

                Candidate_Leptonic_Trail_Lepton_idx = get<1>(lep_candidate_leptonic_trail_final);
                Candidate_Leptonic_Trail_Lepton_type = get<0>(lep_candidate_leptonic_trail_final);
                Candidate_Leptonic_Trail_Lepton_pt = get<2>(lep_candidate_leptonic_trail_final);
                Candidate_Leptonic_Trail_Lepton_eta = get<4>(lep_candidate_leptonic_trail_final);
                Candidate_Leptonic_Trail_Lepton_phi = get<5>(lep_candidate_leptonic_trail_final);


                EventType_Leptonic = Candidate_Leptonic_Lead_Lepton_type + Candidate_Leptonic_Trail_Lepton_type;

                Candidate_Leptonic_Lepton_InvMass = lep_leptonic_mll_final;

                Candidate_Leptonic_ST = st_l;
                Candidate_Leptonic_LT = lt_l;
                Candidate_Leptonic_RpT = rpt_l;
            }
            else if(passed_semileptonic){

                CandidateW_SemiLeptonic_FatJet_idx = tmp_semileptonic_FatJet_idx;
                if(FatJet_particleNet_WvsQCD[CandidateW_SemiLeptonic_FatJet_idx] > FatJet_particleNet_ZvsQCD[CandidateW_SemiLeptonic_FatJet_idx]){
                    CandidateW_SemiLeptonic_FatJet_particlenetScore = FatJet_particleNet_WvsQCD[CandidateW_SemiLeptonic_FatJet_idx];
                    CandidateW_SemiLeptonic_FatJet_flavor = 0;
                }
                else{
                    CandidateW_SemiLeptonic_FatJet_particlenetScore = FatJet_particleNet_ZvsQCD[CandidateW_SemiLeptonic_FatJet_idx];
                    CandidateW_SemiLeptonic_FatJet_flavor = 1;
                }
                CandidateW_SemiLeptonic_FatJet_pt = FatJet_pt[CandidateW_SemiLeptonic_FatJet_idx];
                CandidateW_SemiLeptonic_FatJet_eta = FatJet_eta[CandidateW_SemiLeptonic_FatJet_idx];
                CandidateW_SemiLeptonic_FatJet_phi = FatJet_phi[CandidateW_SemiLeptonic_FatJet_idx];
                CandidateW_SemiLeptonic_FatJet_mass = FatJet_mass[CandidateW_SemiLeptonic_FatJet_idx];
                CandidateW_SemiLeptonic_FatJet_msoftdrop = FatJet_msoftdrop[CandidateW_SemiLeptonic_FatJet_idx];
                
                Candidate_SemiLeptonic_Lepton_idx = get<1>(lep_candidate_semileptonic_final);
                Candidate_SemiLeptonic_Lepton_type = get<0>(lep_candidate_semileptonic_final);
                Candidate_SemiLeptonic_Lepton_pt = get<2>(lep_candidate_semileptonic_final);
                Candidate_SemiLeptonic_Lepton_eta = get<4>(lep_candidate_semileptonic_final);
                Candidate_SemiLeptonic_Lepton_phi = get<5>(lep_candidate_semileptonic_final);

                EventType_SemiLeptonic = Candidate_SemiLeptonic_Lepton_type;

                Candidate_SemiLeptonic_ST = st_s;
                Candidate_SemiLeptonic_LT = lt_s;
                Candidate_SemiLeptonic_RpT = rpt_s;
            }
        }

        if(passed_temp > 0){
            Events_allPassed -> Fill();
            RawEventNumber_FinalPassed++;
            WeightedEventNumber_FinalPassed+=eventScale;
            if(passed_leptonic){
                Events_leptonic -> Fill();
                RawEventNumber_Passed_Leptonic++;
                WeightedEventNumber_Passed_Leptonic+=eventScale;
            }
            else if(passed_semileptonic){
                Events_semileptonic -> Fill();
                RawEventNumber_Passed_SemiLeptonic++;
                WeightedEventNumber_Passed_SemiLeptonic+=eventScale;
            }
        }

        if(passed_temp > 1){ 
            RawEventNumber_Passed_BothChannels++;
            WeightedEventNumber_Passed_BothChannels+=eventScale;
        }
        
    }//event loop
    CutFlow_tree->Fill();
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;


    outfile->cd();
    //Events_allPassed->Print(); 
    Events_allPassed->Write();
    //Events_leptonic->Print(); 
    Events_leptonic->Write();
    //Events_semileptonic->Print(); 
    Events_semileptonic->Write();
    //CutFlow_tree->Print(); 
    CutFlow_tree->Write();

    outfile->Close();
    cout<<setprecision(30)<<WeightedEventNumber_TotalEvents<<"\t"<<WeightedEventNumber_PreSelect_SemiLeptonicChannel<<endl;
    cout<<setprecision(30)<<test_total<<"\t"<<test_sem<<endl;
}

void nanoAODAnalyzer_SelectionTTreeMaker(std::string inputFile, std::string year, std::string name, float xsec, float sum_weights){
//main program
    ROOT::EnableImplicitMT();

    strcpy(Sample_name,name.c_str());
    strcpy(Year,year.c_str());
    SumWeights = sum_weights;
    XS = xsec;
    float lumi = 1.0;

    if(year == "18"){
        deepFlavJetTTHCut = 0.2783;
        lumi = 59.83;
    }
    else if(year == "16APV"){
        deepFlavJetTTHCut = 0.2598;
        lumi = 19.52;
    }
    else if(year == "16"){
        deepFlavJetTTHCut = 0.2489;
        lumi = 16.81;
    }
    else if(year == "17"){
        deepFlavJetTTHCut = 0.3040;
        lumi=41.48;
    }

    EventLoop_background(inputFile,year,xsec,sum_weights,lumi);

}
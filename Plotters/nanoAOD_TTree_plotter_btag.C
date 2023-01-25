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
#include "ROOT/RDataFrame.hxx"

#include<algorithm>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iostream>
#include<mutex>
#include<string>
#include<thread>
#include<vector>
#include <unordered_map>

#include<sys/stat.h>
#include<errno.h>

/* Variables */
    int nSamples=0;

    const int colors[] { 600,0,632, 416,616,432,840,900,880,860,820,800, 920, 893}; // kRed, kBlue, kGreen
    //const int colors[] { 600,632, 416,616,840,900,880,860,622,800,432}; // kRed, kBlue, kGreen
    /*kWhite  = 0,   
    kBlack  = 1,   
    kGray    = 920,
    kYellow = 400, */

    std::vector<TTree*> Background_trees;
    std::vector<TTree*> Signal_trees_vec;
    TTree* Signal_tree;

    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::string, double, TTree*>>> Signal_CatTrees;
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::string, double, TTree*>>> Background_CatTrees;

    //s sampleName, s groupName, f scale factor, s inFile
    std::tuple<std::string, std::string, float, std::string> signalInfo;
    vector<std::tuple<std::string, std::string, float, std::string>> backgroundInfo;
    vector<std::tuple<std::string, std::string, float, std::string>> signalInfo_vec;
    vector<string> groups;

    //const string sample_category[3] = {"bjetloose", "bjetmed", "bjettight"};
    //const string sample_category[2] = {"bjetloose", "bjetlooseinvert"};
    const string sample_category[2] = {"bjetmed", "bjetmedinvert"};
    //const string sample_category[2] = {"bjettight","bjetinvert"};
    const string event_category[4] = {"allPassed", "leptonic", "semileptonic","hadronic"};

    std::unordered_map<std::string, double> lumi_years = {{"16APV" , 19.52}, 
                                                            {"16" , 16.81}, 
                                                            {"17" , 41.48}, 
                                                            {"18" , 59.83}};

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

/*void DrawPlot_separated_background(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    //UNDORATIO TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.45,1.00,1.00);
    //UNDORATIO TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.45);
    
    //UNDORATIO cUp->SetBottomMargin(0.01); 
    //UNDORATIO cDown->SetTopMargin(0.01); 
    //UNDORATIO cDown->SetBottomMargin(0.2); 
        
    //UNDORATIO cUp->Draw();
    //UNDORATIO cDown->Draw();

    
    cout<<"CCCC"<<endl;
    TH1F* h_Signal = new TH1F(("h_Signal_"+outputName).c_str(), title.c_str(), bins, xmin, xmax);
    Signal_tree -> Project(("h_Signal_"+outputName).c_str(), varname.c_str(), selection.c_str());
    h_Signal->Scale(get<2>(signalInfo));
    cout<<"DDDD"<<endl;

    std::vector<TH1F*> h_Background;
    for(int i=0; i<backgroundInfo.size(); i++){
        h_Background.push_back(new TH1F(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Background_trees.at(i)->Project(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        h_Background.at(i)->Scale(get<2>(backgroundInfo.at(i)));
    }

    cout<<"BBBBB"<<endl;
    //Merge Background groups
    std::vector<TH1F*> h_Background_grouped;
    for(int i=0; i<groups.size(); i++){
        vector<int> included;
        for(int j=0; j<backgroundInfo.size(); j++){
            if(get<1>(backgroundInfo[j]) == groups[i]){
                included.push_back(j);
            }
        }
        TH1F* hists_merged = (TH1F*)h_Background[included[0]]->Clone( (varname+"_"+get<0>(backgroundInfo[included[0]])).c_str() );
        for(int i=1; i<included.size(); i++){
            hists_merged->Add(h_Background[included[i]]);
        }
        h_Background_grouped.push_back(hists_merged);
    }

    TH1F* h_Ratio = (TH1F*)h_Signal->Clone("h_Ratio");

    cout<<"AAAA"<<endl;
    std::vector<float> maxima;
    maxima.resize(h_Background_grouped.size() + 1);
    maxima[0] = h_Signal->GetMaximum();
    for(int i=1; i<h_Background_grouped.size()+1; i++){
        maxima[i] = h_Background_grouped[i-1]->GetMaximum(); 
    }
    std::sort(maxima.begin(),maxima.end()); 

    TLegend* legend = new TLegend(0.82, 0.5, 0.99, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04); 

    TPaveStats* st_hists_grouped[h_Background_grouped.size()+1];
    TH1F* h_ratios[h_Background_grouped.size()];

    h_Signal->SetMaximum(maxima.at(maxima.size()-1)*1.05);
    h_Signal->SetMinimum(0.00000000001);
    h_Signal->SetTitle(title.c_str());
    h_Signal-> SetLineColor(colors[0]+1);

    //UNDORATIO cUp->cd();
    //cUp->SetGrid();
    //UNDORATIO cUp->SetLogy();
    gPad->SetLogy();
    h_Signal->Draw("hist,sames");
    legend -> AddEntry(h_Signal, "SIGNAL", "F");
    for(int i=0; i<h_Background_grouped.size(); i++){
        if(h_Background_grouped[i]->Integral() <0.0000000000000000000001) continue;
        int color = colors[i+1] + 1;
        h_Background_grouped[i] -> SetLineColor(color);

        legend -> AddEntry(h_Background_grouped[i], groups[i].c_str(), "F");

        h_Background_grouped[i]->Draw("hist,sames");
        gPad -> Update();

        h_ratios[i] = (TH1F*)h_Ratio->Clone(("histo_ratio_"+to_string(i)).c_str());
        if(h_ratios[i]->GetSumw2N()<=0) h_ratios[i]->Sumw2();
        h_ratios[i] -> Divide(h_Background_grouped[i]);

        h_ratios[i] -> SetMarkerColor(color);
        h_ratios[i] -> SetLineColor(color);
        h_ratios[i] -> SetMarkerSize(0.5);
        h_ratios[i] ->SetTitle("");
        h_ratios[i] ->SetStats(0);
        

    }
    legend -> Draw("same");


    TPaveStats* st_ratio = new TPaveStats();
        
    //UNDORATIO cDown->cd();
    //UNDORATIO cDown->SetGrid();
    

    bool firstdone = false;
    for(int i=0; i<h_Background_grouped.size(); i++){
        if(h_Background_grouped[i]->Integral() <0.000000000000000001) continue;
        if(!firstdone){
    
            h_ratios[i]-> GetYaxis() -> SetTitle("Signal / Background");
            //h_ratios[1] -> SetMaximum(5000);
            //h_ratios[1] -> SetMinimum(-);
            h_ratios[i] -> GetXaxis() -> SetLabelSize(0.07);
            h_ratios[i] -> GetYaxis() -> SetLabelSize(0.07);
            h_ratios[i] -> GetXaxis() -> SetTitleSize(0.07);
            h_ratios[i] -> GetYaxis() -> SetTitleSize(0.07);
            h_ratios[i] -> GetYaxis() -> SetTitleOffset(0.7);
            firstdone = true;
        }
        //UNDORATIO h_ratios[i] -> Draw("e,sames");
    }
    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}

void DrawPlot_separated_background(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    //UNDORATIO TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.45,1.00,1.00);
    //UNDORATIO TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.45);
    
    //UNDORATIO cUp->SetBottomMargin(0.01); 
    //UNDORATIO cDown->SetTopMargin(0.01); 
    //UNDORATIO cDown->SetBottomMargin(0.2); 
        
    //UNDORATIO cUp->Draw();
    //UNDORATIO cDown->Draw();

    
    std::vector<TH1F*> h_Signal;
    for(int i=0; i<signalInfo_vec.size(); i++){
        h_Signal.push_back(new TH1F(("h_Signal2_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Signal_trees_vec.at(i)->Project(("h_Signal2_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        //h_Signal.at(i)->Scale(get<2>(signalInfo_vec.at(i)));
        h_Signal.at(i)->Scale(1/h_Signal.at(i)->Integral());
    }    

    std::vector<TH1F*> h_Background;
    for(int i=0; i<backgroundInfo.size(); i++){
        h_Background.push_back(new TH1F(("h_Background2_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Background_trees.at(i)->Project(("h_Background2_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        //h_Background.at(i)->Scale(get<2>(backgroundInfo.at(i)));
        h_Background.at(i)->Scale(1/h_Background.at(i)->Integral());
    }

    cout<<"BBBBB"<<endl;
    //Merge Background groups
    std::vector<TH1F*> h_Background_grouped;
    for(int i=0; i<groups.size(); i++){
        vector<int> included;
        for(int j=0; j<backgroundInfo.size(); j++){
            if(get<1>(backgroundInfo[j]) == groups[i]){
                included.push_back(j);
            }
        }
        TH1F* hists_merged = (TH1F*)h_Background[included[0]]->Clone( (varname+"_"+get<0>(backgroundInfo[included[0]])).c_str() );
        for(int i=1; i<included.size(); i++){
            hists_merged->Add(h_Background[included[i]]);
        }
        h_Background_grouped.push_back(hists_merged);
    }
    cout<<"AAAA"<<endl;

    std::vector<float> maxima;
    for(int i=0; i<signalInfo_vec.size(); i++){
        maxima.push_back( h_Signal[i]->GetMaximum()); 
    }
    for(int i=1; i<h_Background_grouped.size()+1; i++){
        maxima.push_back(h_Background_grouped[i-1]->GetMaximum()); 
    }
    std::sort(maxima.begin(),maxima.end()); 

    TLegend* legend = new TLegend(0.82, 0.5, 0.99, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04); 

    cout<<"DDDDD"<<endl;

    h_Signal[0]->SetMaximum(maxima.at(maxima.size()-1)*1.05);
    h_Signal[0]->SetMinimum(0.00000000001);
    h_Signal[0]->SetTitle(title.c_str());
    h_Signal[0]-> SetLineColor(colors[0]+1);

    cout<<"EEEE"<<endl;

    //UNDORATIO cUp->cd();
    //cUp->SetGrid();
    //UNDORATIO cUp->SetLogy();
    gPad->SetLogy();

    int color_counter = 0;
    for(int i=0; i<h_Signal.size(); i++){

        int color_sig = colors[color_counter++] + 1;
        h_Signal[i] -> SetLineColor(color_sig);
        legend -> AddEntry(h_Signal[i], (get<0>(signalInfo_vec.at(i))).c_str(), "F");
        h_Signal[i]->Draw("hist,sames");
        gPad -> Update();
    }

    for(int i=0; i<h_Background_grouped.size(); i++){
        if(h_Background_grouped[i]->Integral() <0.0000000000000000000001) continue;
        int color = colors[color_counter++] + 1;
        h_Background_grouped[i] -> SetLineColor(color);

        legend -> AddEntry(h_Background_grouped[i], groups[i].c_str(), "F");

        h_Background_grouped[i]->Draw("hist,sames");
        gPad -> Update();
        

    }
    legend -> Draw("same");

    cout<<"FFFF"<<endl;


    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}*/

void DrawPlot_merged_background(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    std::unordered_map<std::string, std::unordered_map<std::string, TH1F*>> h_Signal;
    std::unordered_map<std::string, std::vector<TH1F*>> h_Background;
    std::unordered_map<std::string, TH1F*> h_Background_Added;
    std::vector<float> maxima;
    for(auto cat : sample_category){

        //Create Signal TH1s
        for(auto sig_tree : Signal_CatTrees[cat]){
            string sig_nm = get<0>(sig_tree);
            string nm = "h_Signal_"+sig_nm+get<2>(sig_tree)+cat+outputName;
            h_Signal[cat][sig_nm] = new TH1F(nm.c_str(), title.c_str(), bins, xmin, xmax);
            get<4>(sig_tree) -> Project(nm.c_str(), varname.c_str(), selection.c_str());
            h_Signal[cat][sig_nm] -> Scale(get<3>(sig_tree));
            maxima.push_back( h_Signal[cat][sig_nm] -> GetMaximum()); 
        }

        //Create BKG TH1s
        for(auto bkg_tree : Background_CatTrees[cat]){
            string nm = "h_Background_"+get<0>(bkg_tree)+get<2>(bkg_tree)+cat+outputName;
            h_Background[cat].push_back(new TH1F(nm.c_str(), title.c_str(), bins, xmin, xmax));
            get<4>(bkg_tree)->Project(nm.c_str(), varname.c_str(), selection.c_str());
            h_Background[cat].back()->Scale(get<3>(bkg_tree));
        }

        h_Background_Added[cat] = (TH1F*)h_Background[cat].at(0)->Clone(("h_bkgAdded_"+cat).c_str());
        for(int i=1; i<h_Background[cat].size(); i++){
            h_Background_Added[cat]->Add(h_Background[cat].at(i));
        }
        maxima.push_back( h_Background_Added[cat]->GetMaximum()); 
    }

    std::sort(maxima.begin(),maxima.end()); 
    double max_val = maxima.at(maxima.size()-1)*1.05;

    TLegend* legend = new TLegend(0.82, 0.7, 0.9, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.03); 

    //Draw the TH1s
    gPad->SetLogy();

    bool base = true;
    int color_counter = 0;
    for(auto cat : sample_category){

        for(auto& h_sig_print : h_Signal[cat]){
            if(base){
                h_sig_print.second->SetMaximum(max_val );
                h_sig_print.second->SetMinimum(0.00000000001);
                h_sig_print.second->SetTitle(title.c_str());
                h_sig_print.second-> SetLineColor(colors[0]+1);
                base = false;
            }

            int color_sig = colors[color_counter++] + 1;
            h_sig_print.second -> SetLineColor(color_sig);
            legend -> AddEntry(h_sig_print.second, (h_sig_print.first + " " + cat).c_str(), "F");
            h_sig_print.second->Draw("hist,sames");
            gPad -> Update();
        }
        
        int color = colors[color_counter++] + 1;
        h_Background_Added[cat] -> SetLineColor(color);
        legend -> AddEntry(h_Background_Added[cat], ("Bkg " + cat).c_str(), "F");
        h_Background_Added[cat]->Draw("hist,sames");
    }

    legend -> Draw("same");

        
    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}

void DrawPlot_merged_background_merged_sig(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    std::unordered_map<std::string, std::vector<TH1F*>> h_Signal;
    std::unordered_map<std::string, std::vector<TH1F*>> h_Background;
    std::unordered_map<std::string, TH1F*> h_Background_Added;
    std::unordered_map<std::string, TH1F*> h_Signal_Added;
    std::vector<float> maxima;
    for(auto cat : sample_category){

        //Create Signal TH1s
        for(auto sig_tree : Signal_CatTrees[cat]){
            string sig_nm = get<0>(sig_tree);
            string nm = "h_Signal_"+sig_nm+get<2>(sig_tree)+cat+outputName+"sigmerged";
            h_Signal[cat].push_back(new TH1F(nm.c_str(), title.c_str(), bins, xmin, xmax));
            get<4>(sig_tree) -> Project(nm.c_str(), varname.c_str(), selection.c_str());
            h_Signal[cat].back() -> Scale(get<3>(sig_tree));
        }

        //Create BKG TH1s
        for(auto bkg_tree : Background_CatTrees[cat]){
            string nm = "h_Background_"+get<0>(bkg_tree)+get<2>(bkg_tree)+cat+outputName+"sigmerged";
            h_Background[cat].push_back(new TH1F(nm.c_str(), title.c_str(), bins, xmin, xmax));
            get<4>(bkg_tree)->Project(nm.c_str(), varname.c_str(), selection.c_str());
            h_Background[cat].back()->Scale(get<3>(bkg_tree));
        }

        h_Signal_Added[cat] = (TH1F*)h_Signal[cat].at(0)->Clone(("h_sigAdded_"+cat+"sigmerged").c_str());
        for(int i=1; i<h_Signal[cat].size(); i++){
            h_Signal_Added[cat]->Add(h_Signal[cat].at(i));
        }
        maxima.push_back( h_Signal_Added[cat]->GetMaximum()); 

        h_Background_Added[cat] = (TH1F*)h_Background[cat].at(0)->Clone(("h_bkgAdded_"+cat+"sigmerged").c_str());
        for(int i=1; i<h_Background[cat].size(); i++){
            h_Background_Added[cat]->Add(h_Background[cat].at(i));
        }
        maxima.push_back( h_Background_Added[cat]->GetMaximum()); 
    }

    std::sort(maxima.begin(),maxima.end()); 
    double max_val = maxima.at(maxima.size()-1)*1.05;

    TLegend* legend = new TLegend(0.82, 0.7, 0.9, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.03); 

    //Draw the TH1s
    gPad->SetLogy();

    bool base = true;
    int color_counter = 0;
    for(auto cat : sample_category){

        if(base){
            h_Signal_Added[cat]->SetMaximum(max_val );
            h_Signal_Added[cat]->SetMinimum(0.00000000001);
            h_Signal_Added[cat]->SetTitle((title+"_SigMerged").c_str());
            h_Signal_Added[cat]-> SetLineColor(colors[0]+1);
            base = false;
        }

        int color_sig = colors[color_counter++] + 1;
        h_Signal_Added[cat] -> SetLineColor(color_sig);
        legend -> AddEntry(h_Signal_Added[cat], ("Signal " + cat).c_str(), "F");
        h_Signal_Added[cat]->Draw("hist,sames");
        gPad -> Update();
        
        int color = colors[color_counter++] + 1;
        h_Background_Added[cat] -> SetLineColor(color);
        legend -> AddEntry(h_Background_Added[cat], ("Bkg " + cat).c_str(), "F");
        h_Background_Added[cat]->Draw("hist,sames");
    }

    legend -> Draw("same");

        
    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}

void SaveHistograms(){
    //CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9

    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT", "CandidateHiggs_FatJet_pt", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx && (Candidate_Leptonic_ST>1000 || Candidate_SemiLeptonic_ST>1000) && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT SemiLeptonic Channel", "CandidateHiggs_FatJet_pt_slep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(Candidate_SemiLeptonic_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT Leptonic Channel", "CandidateHiggs_FatJet_pt_lep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(Candidate_Leptonic_Lead_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background("FatJet_deepTagMD_HbbvsQCD", "Higgs Candidate FatJet DeepTag", "CandidateHiggs_FatJet_DeepTag", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx && (Candidate_Leptonic_ST>1000 || Candidate_SemiLeptonic_ST>1000) && CandidateVBF_Jet_etaSep>4)",8,0.9,1);
    DrawPlot_merged_background("Candidate_Leptonic_Lead_Lepton_pt", "Lepton Lead CandidatePT Leptonic", "Candidate_Leptonic_Lead_Lepton_pt","genWeight*(Candidate_Leptonic_Lead_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background("Candidate_Leptonic_Trail_Lepton_pt", "Lepton Trail CandidatePT Leptonic", "Candidate_Leptonic_Trail_Lepton_pt","genWeight*(Candidate_Leptonic_Trail_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background("Candidate_SemiLeptonic_Lepton_pt", "Lepton CandidatePT SemiLeptonic", "Candidate_SemiLeptonic_Lepton_pt","genWeight*(Candidate_SemiLeptonic_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background("CandidateW_SemiLeptonic_FatJet_pt", "FatJet CandidatePT SemiLeptonic", "CandidateW_SemiLeptonic_FatJet_pt","genWeight*(CandidateW_SemiLeptonic_FatJet_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,200,800);

    /*
    //Electron Isolation Plots
    string s_allElectrons_Leptonic = "(((Candidate_Leptonic_Lead_Lepton_type==0 && Iteration$==Candidate_Leptonic_Lead_Lepton_idx) || (Candidate_Leptonic_Lead_Lepton_idx != -1 && Candidate_Leptonic_Trail_Lepton_type==0 && Iteration$==Candidate_Leptonic_Trail_Lepton_idx)) && CandidateHiggs_FatJet_particlenetScore>0.9)";
    DrawPlot_merged_background("Electron_pt", "Electron Candidate Leptonic Channel PT", "CandidateElectron_pt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_phi", "Electron Candidate Leptonic Channel Phi", "CandidateElectron_phi_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_eta", "Electron Candidate Leptonic Channel Eta", "CandidateElectron_eta_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_dr03TkSumPt", "Electron Candidate Leptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate Leptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_leptonic", s_allElectrons_Leptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
    string s_allMuons_Leptonic = "(((Candidate_Leptonic_Lead_Lepton_type==1  && Iteration$==Candidate_Leptonic_Lead_Lepton_idx) ||(Candidate_Leptonic_Lead_Lepton_idx != -1 && Candidate_Leptonic_Trail_Lepton_type==1 && Iteration$==Candidate_Leptonic_Trail_Lepton_idx)) && CandidateHiggs_FatJet_particlenetScore>0.9)";
    DrawPlot_merged_background("Muon_pt", "Muon Candidate Leptonic Channel PT", "CandidateMuon_pt_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(),80,0,800);
    DrawPlot_merged_background("Muon_phi", "Muon Candidate Leptonic Channel Phi", "CandidateMuon_phi_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background("Muon_eta", "Muon Candidate Leptonic Channel Eta", "CandidateMuon_eta_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background("Muon_tkRelIso", "Muon Candidate Leptonic Channel tkRelIso", "CandidateMuon_kRelIso_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate Leptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_leptonic", s_allMuons_Leptonic,80,0,800,20, 0, 5);


    //Electron Isolation Plots
    string s_allElectrons_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_type==0 && Iteration$==Candidate_SemiLeptonic_Lepton_idx  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9)";
    DrawPlot_merged_background("Electron_pt", "Electron Candidate SemiLeptonic Channel PT", "CandidateElectron_pt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_phi", "Electron Candidate SemiLeptonic Channel Phi", "CandidateElectron_phi_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_eta", "Electron Candidate SemiLeptonic Channel Eta", "CandidateElectron_eta_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_dr03TkSumPt", "Electron Candidate SemiLeptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate SemiLeptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_semileptonic", s_allElectrons_SemiLeptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
    string s_allMuons_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_type==1 && Iteration$==Candidate_SemiLeptonic_Lepton_idx  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9)";
    DrawPlot_merged_background("Muon_pt", "Muon Candidate SemiLeptonic Channel PT", "CandidateMuon_pt_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(),80,0,800);
    DrawPlot_merged_background("Muon_phi", "Muon Candidate SemiLeptonic Channel Phi", "CandidateMuon_phi_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Muon_eta", "Muon Candidate SemiLeptonic Channel Eta", "CandidateMuon_eta_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Muon_tkRelIso", "Muon Candidate SemiLeptonic Channel tkRelIso", "CandidateMuon_kRelIso_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate SemiLeptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_semileptonic", s_allMuons_SemiLeptonic,80,0,800,20, 0, 5);
    */
    //MET
    string s_Leptonic = "(Candidate_Leptonic_Lead_Lepton_idx != -1  && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)";
    string s_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_idx != -1  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)";

    DrawPlot_merged_background("PuppiMET_pt", "PuppiMET pt SemiLeptonic Channel", "PuppiMET_pt_semileptonic", ("genWeight*"+s_SemiLeptonic).c_str(),40,0,1000);
    DrawPlot_merged_background("PuppiMET_pt", "PuppiMET pt Leptonic Channel", "PuppiMET_pt_leptonic", ("genWeight*"+s_Leptonic).c_str(),40,0,1000);

    DrawPlot_merged_background("Jet_pt", "VBF Candidate Lead Jet PT", "CandidateVBF_Lead_Jet_pt_semileptonic", "genWeight*(Iteration$==CandidateVBF_Lead_Jet_idx)*"+s_SemiLeptonic,40,0,400);
    DrawPlot_merged_background("Jet_pt", "VBF Candidate Trail Jet PT", "CandidateVBF_Trail_Jet_pt_semileptonic", "genWeight*(Iteration$==CandidateVBF_Trail_Jet_idx)*"+s_SemiLeptonic,40,0,400);
    DrawPlot_merged_background("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep_semileptonic", "genWeight*"+s_SemiLeptonic,20,5,10);
    DrawPlot_merged_background("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass_semileptonic", "genWeight*"+s_SemiLeptonic,100,0,8000);


    DrawPlot_merged_background("Jet_pt", "VBF Candidate Lead Jet PT", "CandidateVBF_Lead_Jet_pt_leptonic", "genWeight*(Iteration$==CandidateVBF_Lead_Jet_idx)*"+s_Leptonic,40,0,400);
    DrawPlot_merged_background("Jet_pt", "VBF Candidate Trail Jet PT", "CandidateVBF_Trail_Jet_pt_leptonic", "genWeight*(Iteration$==CandidateVBF_Trail_Jet_idx)*"+s_Leptonic,40,0,400);
    DrawPlot_merged_background("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep_leptonic", "genWeight*"+s_Leptonic,100,0,10);
    DrawPlot_merged_background("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass_leptonic", "genWeight*"+s_Leptonic,100,0,8000);

    DrawPlot_merged_background("Candidate_Leptonic_LT", "LT", "Candidate_Leptonic_LT", "genWeight*"+s_Leptonic,100,0,8000);
    DrawPlot_merged_background("Candidate_SemiLeptonic_LT", "LT", "Candidate_SemiLeptonic_LT", "genWeight*"+s_SemiLeptonic,100,0,8000);
    DrawPlot_merged_background("Candidate_Leptonic_ST", "ST", "Candidate_Leptonic_ST", "genWeight*"+s_Leptonic,100,0,8000);
    DrawPlot_merged_background("Candidate_SemiLeptonic_ST", "ST", "Candidate_SemiLeptonic_ST", "genWeight*"+s_SemiLeptonic,100,0,8000);
    
    //DrawPlot_separated_background("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep", "genWeight",20,5,10);
    //DrawPlot_separated_background("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass", "genWeight",100,0,8000);
    

    //CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9

    DrawPlot_merged_background_merged_sig("FatJet_pt", "Higgs Candidate FatJet PT", "CandidateHiggs_FatJet_pt", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx && (Candidate_Leptonic_ST>1000 || Candidate_SemiLeptonic_ST>1000) && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background_merged_sig("FatJet_pt", "Higgs Candidate FatJet PT SemiLeptonic Channel", "CandidateHiggs_FatJet_pt_slep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(Candidate_SemiLeptonic_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background_merged_sig("FatJet_pt", "Higgs Candidate FatJet PT Leptonic Channel", "CandidateHiggs_FatJet_pt_lep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(Candidate_Leptonic_Lead_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",100,0,10000);
    DrawPlot_merged_background_merged_sig("FatJet_deepTagMD_HbbvsQCD", "Higgs Candidate FatJet DeepTag", "CandidateHiggs_FatJet_DeepTag", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx && (Candidate_Leptonic_ST>1000 || Candidate_SemiLeptonic_ST>1000) && CandidateVBF_Jet_etaSep>4)",8,0.9,1);
    DrawPlot_merged_background_merged_sig("Candidate_Leptonic_Lead_Lepton_pt", "Lepton Lead CandidatePT Leptonic", "Candidate_Leptonic_Lead_Lepton_pt","genWeight*(Candidate_Leptonic_Lead_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background_merged_sig("Candidate_Leptonic_Trail_Lepton_pt", "Lepton Trail CandidatePT Leptonic", "Candidate_Leptonic_Trail_Lepton_pt","genWeight*(Candidate_Leptonic_Trail_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && Candidate_Leptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background_merged_sig("Candidate_SemiLeptonic_Lepton_pt", "Lepton CandidatePT SemiLeptonic", "Candidate_SemiLeptonic_Lepton_pt","genWeight*(Candidate_SemiLeptonic_Lepton_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,0,400);
    DrawPlot_merged_background_merged_sig("CandidateW_SemiLeptonic_FatJet_pt", "FatJet CandidatePT SemiLeptonic", "CandidateW_SemiLeptonic_FatJet_pt","genWeight*(CandidateW_SemiLeptonic_FatJet_idx!=-1 && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9 && Candidate_SemiLeptonic_ST>1000 && CandidateVBF_Jet_etaSep>4)",40,200,800);

    /*
    //Electron Isolation Plots
   // string s_allElectrons_Leptonic = "(((Candidate_Leptonic_Lead_Lepton_type==0 && Iteration$==Candidate_Leptonic_Lead_Lepton_idx) || (Candidate_Leptonic_Lead_Lepton_idx != -1 && Candidate_Leptonic_Trail_Lepton_type==0 && Iteration$==Candidate_Leptonic_Trail_Lepton_idx)) && CandidateHiggs_FatJet_particlenetScore>0.9)";
    DrawPlot_merged_background_merged_sig("Electron_pt", "Electron Candidate Leptonic Channel PT", "CandidateElectron_pt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_phi", "Electron Candidate Leptonic Channel Phi", "CandidateElectron_phi_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_eta", "Electron Candidate Leptonic Channel Eta", "CandidateElectron_eta_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_dr03TkSumPt", "Electron Candidate Leptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate Leptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_leptonic", s_allElectrons_Leptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
   // string s_allMuons_Leptonic = "(((Candidate_Leptonic_Lead_Lepton_type==1  && Iteration$==Candidate_Leptonic_Lead_Lepton_idx) ||(Candidate_Leptonic_Lead_Lepton_idx != -1 && Candidate_Leptonic_Trail_Lepton_type==1 && Iteration$==Candidate_Leptonic_Trail_Lepton_idx)) && CandidateHiggs_FatJet_particlenetScore>0.9)";
    DrawPlot_merged_background_merged_sig("Muon_pt", "Muon Candidate Leptonic Channel PT", "CandidateMuon_pt_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(),80,0,800);
    DrawPlot_merged_background_merged_sig("Muon_phi", "Muon Candidate Leptonic Channel Phi", "CandidateMuon_phi_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background_merged_sig("Muon_eta", "Muon Candidate Leptonic Channel Eta", "CandidateMuon_eta_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background_merged_sig("Muon_tkRelIso", "Muon Candidate Leptonic Channel tkRelIso", "CandidateMuon_kRelIso_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate Leptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_leptonic", s_allMuons_Leptonic,80,0,800,20, 0, 5);


    //Electron Isolation Plots
   // string s_allElectrons_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_type==0 && Iteration$==Candidate_SemiLeptonic_Lepton_idx  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9)";
    DrawPlot_merged_background_merged_sig("Electron_pt", "Electron Candidate SemiLeptonic Channel PT", "CandidateElectron_pt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_phi", "Electron Candidate SemiLeptonic Channel Phi", "CandidateElectron_phi_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_eta", "Electron Candidate SemiLeptonic Channel Eta", "CandidateElectron_eta_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background_merged_sig("Electron_dr03TkSumPt", "Electron Candidate SemiLeptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate SemiLeptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_semileptonic", s_allElectrons_SemiLeptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
   // string s_allMuons_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_type==1 && Iteration$==Candidate_SemiLeptonic_Lepton_idx  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9)";
    DrawPlot_merged_background_merged_sig("Muon_pt", "Muon Candidate SemiLeptonic Channel PT", "CandidateMuon_pt_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(),80,0,800);
    DrawPlot_merged_background_merged_sig("Muon_phi", "Muon Candidate SemiLeptonic Channel Phi", "CandidateMuon_phi_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background_merged_sig("Muon_eta", "Muon Candidate SemiLeptonic Channel Eta", "CandidateMuon_eta_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background_merged_sig("Muon_tkRelIso", "Muon Candidate SemiLeptonic Channel tkRelIso", "CandidateMuon_kRelIso_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate SemiLeptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_semileptonic", s_allMuons_SemiLeptonic,80,0,800,20, 0, 5);
    */
    //MET
   // string s_Leptonic = "(Candidate_Leptonic_Lead_Lepton_idx != -1  && CandidateHiggs_FatJet_particlenetScore>0.9)";
   // string s_SemiLeptonic = "(Candidate_SemiLeptonic_Lepton_idx != -1  && CandidateHiggs_FatJet_particlenetScore>0.9 && CandidateW_SemiLeptonic_FatJet_particlenetScore >0.9)";

    DrawPlot_merged_background_merged_sig("PuppiMET_pt", "PuppiMET pt SemiLeptonic Channel", "PuppiMET_pt_semileptonic", ("genWeight*"+s_SemiLeptonic).c_str(),40,0,1000);
    DrawPlot_merged_background_merged_sig("PuppiMET_pt", "PuppiMET pt Leptonic Channel", "PuppiMET_pt_leptonic", ("genWeight*"+s_Leptonic).c_str(),40,0,1000);

    DrawPlot_merged_background_merged_sig("Jet_pt", "VBF Candidate Lead Jet PT", "CandidateVBF_Lead_Jet_pt_semileptonic", "genWeight*(Iteration$==CandidateVBF_Lead_Jet_idx)*"+s_SemiLeptonic,40,0,400);
    DrawPlot_merged_background_merged_sig("Jet_pt", "VBF Candidate Trail Jet PT", "CandidateVBF_Trail_Jet_pt_semileptonic", "genWeight*(Iteration$==CandidateVBF_Trail_Jet_idx)*"+s_SemiLeptonic,40,0,400);
    DrawPlot_merged_background_merged_sig("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep_semileptonic", "genWeight*"+s_SemiLeptonic,20,5,10);
    DrawPlot_merged_background_merged_sig("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass_semileptonic", "genWeight*"+s_SemiLeptonic,100,0,8000);


    DrawPlot_merged_background_merged_sig("Jet_pt", "VBF Candidate Lead Jet PT", "CandidateVBF_Lead_Jet_pt_leptonic", "genWeight*(Iteration$==CandidateVBF_Lead_Jet_idx)*"+s_Leptonic,40,0,400);
    DrawPlot_merged_background_merged_sig("Jet_pt", "VBF Candidate Trail Jet PT", "CandidateVBF_Trail_Jet_pt_leptonic", "genWeight*(Iteration$==CandidateVBF_Trail_Jet_idx)*"+s_Leptonic,40,0,400);
    DrawPlot_merged_background_merged_sig("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep_leptonic", "genWeight*"+s_Leptonic,100,0,10);
    DrawPlot_merged_background_merged_sig("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass_leptonic", "genWeight*"+s_Leptonic,100,0,8000);

    DrawPlot_merged_background_merged_sig("Candidate_Leptonic_LT", "LT", "Candidate_Leptonic_LT", "genWeight*"+s_Leptonic,100,0,8000);
    DrawPlot_merged_background_merged_sig("Candidate_SemiLeptonic_LT", "LT", "Candidate_SemiLeptonic_LT", "genWeight*"+s_SemiLeptonic,100,0,8000);
    DrawPlot_merged_background_merged_sig("Candidate_Leptonic_ST", "ST", "Candidate_Leptonic_ST", "genWeight*"+s_Leptonic,100,0,8000);
    DrawPlot_merged_background_merged_sig("Candidate_SemiLeptonic_ST", "ST", "Candidate_SemiLeptonic_ST", "genWeight*"+s_SemiLeptonic,100,0,8000);
}

void nanoAOD_TTree_plotter_btag(){

    std::string eos_path = "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/AllYears/";
    std::unordered_map<std::string, double> bkg_xs = {{"DYJets", 6077.22 },
                                                    {"ttH", 0.6 },
                                                    {"QCD_Pt_300to470", 6826.0 },
                                                    {"QCD_Pt_470to600", 552.6 },
                                                    {"QCD_Pt_600to800", 156.6 },
                                                    {"QCD_Pt_800to1000", 26.25 },
                                                    {"QCD_Pt_1000to1400", 7.5 },
                                                    {"QCD_Pt_1400to1800", 0.6479 },
                                                    {"QCD_Pt_1800to2400", 0.08715 },
                                                    {"QCD_Pt_2400to3200", 0.005242 },
                                                    {"QCD_Pt_3200toInf", 0.0001349 },
                                                    {"WJets", 52850.0 },
                                                    {"WW", 75.95 },
                                                    {"WZ", 27.59 },
                                                    {"ZZ", 12.17 },
                                                    {"TTbb_TTToHadronic", 5.5 },
                                                    {"TTbb_TTTo2L2Nu", 4.0 },
                                                    {"TTZToLLNuNu", 0.2432 },
                                                    {"TTZToQQ", 0.5104 },
                                                    {"TTWJetsToLNu", 0.2161 },
                                                    {"TTWJetsToQQ", 0.4377 },
                                                    {"TTJets", 831.76 }};
    std::unordered_map<std::string, std::string> bkg_groups = {{"DYJets","DYJets"},
                                                            {"ttH","ttH"},
                                                            {"QCD_Pt_300to470","QCD"},
                                                            {"QCD_Pt_470to600","QCD"},
                                                            {"QCD_Pt_600to800","QCD"},
                                                            {"QCD_Pt_800to1000","QCD"},
                                                            {"QCD_Pt_1000to1400","QCD"},
                                                            {"QCD_Pt_1400to1800","QCD"},
                                                            {"QCD_Pt_1800to2400","QCD"},
                                                            {"QCD_Pt_2400to3200","QCD"},
                                                            {"QCD_Pt_3200toInf","QCD"},
                                                            {"WJets","WJets"},
                                                            {"WW","WW"},
                                                            {"WZ","WZ"},
                                                            {"ZZ","ZZ"},
                                                            {"TTbb_TTToHadronic","TTbb"},
                                                            {"TTbb_TTTo2L2Nu","TTbb"},
                                                            {"TTZToLLNuNu","TTZ"},
                                                            {"TTZToQQ","TTZ"},
                                                            {"TTWJetsToLNu","TTWJets"},
                                                            {"TTWJetsToQQ","TTWJets"},
                                                            {"TTJets","TTJets"}};

    std::unordered_map<std::string, double> sig_xs = {{"WWH_ucsd_C2V_Reweight" , 0.00243354}, 
                                                     {"WZH_ucsd_C2V_Reweight" , 0.00155167}};

    vector<TFile*> tree_infile;
    for(auto cat : sample_category){
        for(auto yr : lumi_years) {
            for(auto bkg : bkg_xs){
                string filename = eos_path + "BDT_Settings_loose11182022_" + cat + "/" + yr.first + "/" + bkg.first + "/"  + bkg.first + "_Full.root";
                cout<<filename<<endl;
                tree_infile.push_back(TFile::Open(filename.c_str()));

                auto df = ROOT::RDataFrame("CutFlow", filename.c_str());
                double sum_weights = df.Sum("WeightedEventNumber_TotalEvents").GetValue();
                float scale = 1000 * yr.second * bkg.second / sum_weights;

                Background_CatTrees[cat].push_back(std::make_tuple( bkg.first,
                                                                    bkg_groups[bkg.first],
                                                                    yr.first,
                                                                    scale,
                                                                    (TTree*)tree_infile.back()->Get("Events_allPassed")
                                                                    ));
            }
        }
    }

    vector<TFile*> tree_infile_sig;
    for(auto cat : sample_category){
        for(auto yr : lumi_years) {
            for(auto sig : sig_xs){
                string filename = eos_path + "BDT_Settings_loose11182022_" + cat + "/" + yr.first + "/" + sig.first + "/"  + sig.first + "_0.root";
                tree_infile_sig.push_back(TFile::Open(filename.c_str()));

                auto df = ROOT::RDataFrame("CutFlow", filename.c_str());
                double sum_weights = df.Sum("WeightedEventNumber_TotalEvents").GetValue();
                float scale = 1000 * yr.second * sig.second / sum_weights;

                Signal_CatTrees[cat].push_back(std::make_tuple( sig.first,
                                                                    "Sig",
                                                                    yr.first,
                                                                    scale,
                                                                    (TTree*)tree_infile_sig.back()->Get("Events_allPassed")
                                                                    ));
            }
        }
    }


    SaveHistograms();
}
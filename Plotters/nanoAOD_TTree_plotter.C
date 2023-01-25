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

    //s sampleName, s groupName, f scale factor, s inFile
    std::tuple<std::string, std::string, float, std::string> signalInfo;
    vector<std::tuple<std::string, std::string, float, std::string>> backgroundInfo;
    vector<std::tuple<std::string, std::string, float, std::string>> signalInfo_vec;
    vector<string> groups;

    const string event_category[4] = {"allPassed", "leptonic", "semileptonic","hadronic"};


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
}*/

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
}

void DrawPlot_merged_background(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);
    float scale_factor = 137.65*1000;


    std::vector<TH1F*> h_Signal;
    for(int i=0; i<signalInfo_vec.size(); i++){
        h_Signal.push_back(new TH1F(("h_Signal_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Signal_trees_vec.at(i)->Project(("h_Signal_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        h_Signal.at(i)->Scale(get<2>(signalInfo_vec.at(i))*scale_factor);
    }    

    std::vector<TH1F*> h_Background;
    for(int i=0; i<backgroundInfo.size(); i++){
        h_Background.push_back(new TH1F(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Background_trees.at(i)->Project(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        h_Background.at(i)->Scale(get<2>(backgroundInfo.at(i))*scale_factor);
    }

    TH1F* h_Background_Added = (TH1F*)h_Background.at(0)->Clone("h_Temp");
    for(int i=1; i<backgroundInfo.size(); i++){
        h_Background_Added->Add(h_Background.at(i));
    }



    std::vector<float> maxima;
    for(int i=0; i<signalInfo_vec.size(); i++){
        maxima.push_back( h_Signal[i]->GetMaximum()); 
    }
    maxima.push_back( h_Background_Added->GetMaximum()); 
    std::sort(maxima.begin(),maxima.end()); 

    TLegend* legend = new TLegend(0.82, 0.7, 0.9, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.03); 

    h_Signal[0]->SetMaximum(maxima.at(maxima.size()-1)*1.05);
    h_Signal[0]->SetMinimum(0.00000000001);
    h_Signal[0]->SetTitle(title.c_str());
    h_Signal[0]-> SetLineColor(colors[0]+1);

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

    cout<<"6"<<endl;
    int color = colors[color_counter++] + 1;
    h_Background_Added -> SetLineColor(color);
    cout<<"Background: "<<h_Background_Added->Integral()<<"Signal: "<<h_Signal[0]->Integral()<<endl;

    legend -> AddEntry(h_Background_Added, "BACKGROUND", "F");

    h_Background_Added->Draw("hist,sames");
    legend -> Draw("same");

        
    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}

void DrawPlot_2d(string varname, string title, string outputName, string selection = "", int xbins = 20, float xmin = 0.0, float xmax = 2000.0,int ybins = 20, float ymin = 0.0, float ymax = 2000.0){


    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    TH2F* h_Signal = new TH2F(("h_Signal_"+outputName).c_str(), title.c_str(), xbins, xmin, xmax,ybins,ymin,ymax);
    Signal_tree -> Project(("h_Signal_"+outputName).c_str(), varname.c_str(), selection.c_str());

    std::vector<TH2F*> h_Background;
    for(int i=0; i<backgroundInfo.size(); i++){
        h_Background.push_back(new TH2F(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), title.c_str(), xbins, xmin, xmax,ybins,ymin,ymax));
        Background_trees.at(i)->Project(("h_Background_"+get<0>(backgroundInfo.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
    }

    //Merge Background groups
    std::vector<TH2F*> h_Background_grouped;
    for(int i=0; i<groups.size(); i++){
        vector<int> included;
        for(int j=0; j<backgroundInfo.size(); j++){
            if(get<1>(backgroundInfo[j]) == groups[i]){
                included.push_back(j);
            }
        }
        TH2F* hists_merged = (TH2F*)h_Background[included[0]]->Clone( (varname+"_"+get<0>(backgroundInfo[included[0]])).c_str() );
        for(int i=1; i<included.size(); i++){
            hists_merged->Add(h_Background[included[i]]);
        }
        h_Background_grouped.push_back(hists_merged);
    }

    TLegend* legend = new TLegend(0.82, 0.5, 0.99, 0.9);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04); 


    h_Signal->SetTitle(title.c_str());
    h_Signal-> SetMarkerColor(colors[0]+1);
    h_Signal-> SetMarkerStyle(8);
    h_Signal-> SetMarkerSize(5);

    h_Signal->Draw("hist,sames");
    legend -> AddEntry(h_Signal, "SIGNAL", "P");
    for(int i=0; i<h_Background_grouped.size(); i++){
        int color = colors[i+1] + 1;
        h_Background_grouped[i] -> SetMarkerColor(color);
        h_Background_grouped[i] -> SetMarkerStyle(8);
        h_Background_grouped[i] -> SetMarkerSize(5);

        legend -> AddEntry(h_Background_grouped[i], groups[i].c_str(), "P");

        h_Background_grouped[i]->Draw("hist,sames");
        gPad -> Update();
        

    }
    legend -> Draw("same");


    TF1* f_const = new TF1("f_1", "0.1*x",xmin, xmax);
    f_const -> SetLineColor(kRed);
    f_const -> SetLineWidth(2);
    f_const->Draw("sames");

    TF1* f_2 = new TF1("f_2", "0.05*x",xmin, xmax);
    f_2 -> SetLineColor(kGreen);
    f_2 -> SetLineWidth(2);
    f_2->Draw("sames");

    TF1* f_3 = new TF1("f_2", "0.01*x",xmin, xmax);
    f_3 -> SetLineColor(kBlue);
    f_3 -> SetLineWidth(2);
    f_3->Draw("sames");

    c1->SaveAs((outdir+outputName+".png").c_str());
    c1->SaveAs((outdir+outputName+".pdf").c_str());

    c1->Clear();
}

void ReadInTrees(){

    for(int i=0; i<backgroundInfo.size(); i++){
        groups.push_back(get<1>(backgroundInfo[i]));
    }
    std::sort(groups.begin(),groups.end());
    groups.erase(unique(groups.begin(),groups.end()),groups.end());
    int nGroups = groups.size();


    //TFile *signal_infile = TFile::Open(get<3>(signalInfo).c_str());
    TFile *tree_infile[backgroundInfo.size()];
    TFile *tree_signal_infile[signalInfo_vec.size()];
    //Signal_tree = (TTree*)signal_infile->Get("Events_allPassed");

    for(int i=0; i<backgroundInfo.size(); i++){
        tree_infile[i] = TFile::Open(get<3>(backgroundInfo[i]).c_str());
        cout<<cur_time()<<"Reading in File: "<< (get<3>(backgroundInfo[i]))<<endl;

        Background_trees.push_back((TTree*)tree_infile[i]->Get("Events_allPassed"));
    }

    for(int i=0; i<signalInfo_vec.size(); i++){
        tree_signal_infile[i] = TFile::Open(get<3>(signalInfo_vec[i]).c_str());
        cout<<cur_time()<<"Reading in File: "<< (get<3>(signalInfo_vec[i]))<<endl;

        Signal_trees_vec.push_back((TTree*)tree_signal_infile[i]->Get("Events_allPassed"));
    }

    cout<<"Finished Reading in Trees"<<endl;
}

void SaveHistograms(){

    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT", "CandidateHiggs_FatJet_pt", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)",100,0,10000);
    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT SemiLeptonic Channel", "CandidateHiggs_FatJet_pt_slep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(CandidateW_SemiLeptonic_Lepton_idx!=-1)",100,0,10000);
    DrawPlot_merged_background("FatJet_pt", "Higgs Candidate FatJet PT Leptonic Channel", "CandidateHiggs_FatJet_pt_lep", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)*(CandidateW_Leptonic_Lead_Lepton_idx!=-1)",100,0,10000);
    DrawPlot_merged_background("FatJet_deepTagMD_HbbvsQCD", "Higgs Candidate FatJet DeepTag", "CandidateHiggs_FatJet_DeepTag", "genWeight*(Iteration$==CandidateHiggs_FatJet_idx)",8,0.9,1);
    DrawPlot_merged_background("CandidateW_Leptonic_Lead_Lepton_type", "Muons CandidatePT", "CandidateMuons_FatJet_type","",3,-1,2);
    DrawPlot_merged_background("CandidateW_Leptonic_Trail_Lepton_type", "Muons CandidatePT", "CandidateMuons_FatJet_type_trail","",3,-1,2);
    DrawPlot_merged_background("CandidateW_Leptonic_Lead_Lepton_pt", "Lepton Lead CandidatePT Leptonic", "CandidateW_Leptonic_Lead_Lepton_pt","genWeight*(CandidateW_Leptonic_Lead_Lepton_idx!=-1)",40,0,400);
    DrawPlot_merged_background("CandidateW_Leptonic_Trail_Lepton_pt", "Lepton Trail CandidatePT Leptonic", "CandidateW_Leptonic_Trail_Lepton_pt","genWeight*(CandidateW_Leptonic_Trail_Lepton_idx!=-1)",40,0,400);
    DrawPlot_merged_background("CandidateW_SemiLeptonic_Lepton_pt", "Lepton CandidatePT SemiLeptonic", "CandidateW_SemiLeptonic_Lepton_pt","genWeight*(CandidateW_SemiLeptonic_Lepton_idx!=-1)",40,0,400);
    DrawPlot_merged_background("CandidateW_SemiLeptonic_FatJet_pt", "FatJet CandidatePT SemiLeptonic", "CandidateW_SemiLeptonic_FatJet_pt","genWeight*(CandidateW_SemiLeptonic_FatJet_idx!=-1)",40,200,800);

    //Electron Isolation Plots
    string s_allElectrons_Leptonic = "((CandidateW_Leptonic_Lead_Lepton_type==0 && Iteration$==CandidateW_Leptonic_Lead_Lepton_idx) || (CandidateW_Leptonic_Lead_Lepton_idx != -1 && CandidateW_Leptonic_Trail_Lepton_type==0 && Iteration$==CandidateW_Leptonic_Trail_Lepton_idx))";
    DrawPlot_merged_background("Electron_pt", "Electron Candidate Leptonic Channel PT", "CandidateElectron_pt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_phi", "Electron Candidate Leptonic Channel Phi", "CandidateElectron_phi_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_eta", "Electron Candidate Leptonic Channel Eta", "CandidateElectron_eta_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str());
    DrawPlot_merged_background("Electron_dr03TkSumPt", "Electron Candidate Leptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_leptonic", ("genWeight*"+s_allElectrons_Leptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate Leptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_leptonic", s_allElectrons_Leptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
    string s_allMuons_Leptonic = "((CandidateW_Leptonic_Lead_Lepton_type==1  && Iteration$==CandidateW_Leptonic_Lead_Lepton_idx) ||(CandidateW_Leptonic_Lead_Lepton_idx != -1 && CandidateW_Leptonic_Trail_Lepton_type==1 && Iteration$==CandidateW_Leptonic_Trail_Lepton_idx))";
    DrawPlot_merged_background("Muon_pt", "Muon Candidate Leptonic Channel PT", "CandidateMuon_pt_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(),80,0,800);
    DrawPlot_merged_background("Muon_phi", "Muon Candidate Leptonic Channel Phi", "CandidateMuon_phi_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background("Muon_eta", "Muon Candidate Leptonic Channel Eta", "CandidateMuon_eta_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str());
    DrawPlot_merged_background("Muon_tkRelIso", "Muon Candidate Leptonic Channel tkRelIso", "CandidateMuon_kRelIso_leptonic", ("genWeight*"+s_allMuons_Leptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate Leptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_leptonic", s_allMuons_Leptonic,80,0,800,20, 0, 5);


    //Electron Isolation Plots
    string s_allElectrons_SemiLeptonic = "(CandidateW_SemiLeptonic_Lepton_type==0 && Iteration$==CandidateW_SemiLeptonic_Lepton_idx)";
    DrawPlot_merged_background("Electron_pt", "Electron Candidate SemiLeptonic Channel PT", "CandidateElectron_pt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_phi", "Electron Candidate SemiLeptonic Channel Phi", "CandidateElectron_phi_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_eta", "Electron Candidate SemiLeptonic Channel Eta", "CandidateElectron_eta_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Electron_dr03TkSumPt", "Electron Candidate SemiLeptonic Channel dr03TkSumPt", "CandidateElectron_dr03TkSumPt_semileptonic", ("genWeight*"+s_allElectrons_SemiLeptonic).c_str(), 100, 0, 10);
    //DrawPlot_2d("Electron_dr03TkSumPt:Electron_pt", "Electron Candidate SemiLeptonic Channel PT vs dr03TkSumP", "CandidateElectron_ptvsdr03TkSumP_semileptonic", s_allElectrons_SemiLeptonic,100,0,1500,100, 0, 10);

    //Muon Isolation Plots
    string s_allMuons_SemiLeptonic = "(CandidateW_SemiLeptonic_Lepton_type==1 && Iteration$==CandidateW_SemiLeptonic_Lepton_idx)";
    DrawPlot_merged_background("Muon_pt", "Muon Candidate SemiLeptonic Channel PT", "CandidateMuon_pt_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(),80,0,800);
    DrawPlot_merged_background("Muon_phi", "Muon Candidate SemiLeptonic Channel Phi", "CandidateMuon_phi_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Muon_eta", "Muon Candidate SemiLeptonic Channel Eta", "CandidateMuon_eta_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str());
    DrawPlot_merged_background("Muon_tkRelIso", "Muon Candidate SemiLeptonic Channel tkRelIso", "CandidateMuon_kRelIso_semileptonic", ("genWeight*"+s_allMuons_SemiLeptonic).c_str(), 10, 0, 5);
    //DrawPlot_2d("Muon_tkRelIso:Muon_pt", "Muon Candidate SemiLeptonic Channel PT vs tkRelIso", "CandidateMuon_ptvstkRelIso_semileptonic", s_allMuons_SemiLeptonic,80,0,800,20, 0, 5);

    //MET
    string s_Leptonic = "(CandidateW_Leptonic_Lead_Lepton_idx != -1)";
    string s_SemiLeptonic = "(CandidateW_SemiLeptonic_Lepton_idx != -1)";

    DrawPlot_merged_background("PuppiMET_pt", "PuppiMET pt SemiLeptonic Channel", "PuppiMET_pt_semileptonic", ("genWeight*"+s_SemiLeptonic).c_str(),40,0,1000);
    DrawPlot_merged_background("PuppiMET_pt", "PuppiMET pt Leptonic Channel", "PuppiMET_pt_leptonic", ("genWeight*"+s_Leptonic).c_str(),40,0,1000);

    DrawPlot_merged_background("Jet_pt", "VBF Candidate Lead Jet PT", "CandidateVBF_Lead_Jet_pt", "genWeight*(Iteration$==CandidateVBF_Lead_Jet_idx)",40,0,400);
    DrawPlot_merged_background("Jet_pt", "VBF Candidate Trail Jet PT", "CandidateVBF_Trail_Jet_pt", "genWeight*(Iteration$==CandidateVBF_Trail_Jet_idx)",40,0,400);
    DrawPlot_merged_background("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep", "genWeight",20,5,10);
    DrawPlot_merged_background("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass", "genWeight",100,0,8000);
    
    //DrawPlot_separated_background("CandidateVBF_Jet_etaSep", "VBF Candidate dEta", "CandidateVBF_Jet_etaSep", "genWeight",20,5,10);
    //DrawPlot_separated_background("CandidateVBF_Jet_invMass", "VBF Candidate InvMass", "CandidateVBF_Jet_invMass", "genWeight",100,0,8000);
    

}

void nanoAOD_TTree_plotter(){

    signalInfo = std::make_tuple("C2V_5","C2V_5",(0.0303768 / 457.691895), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_5/C2V_ucsd_5_0.root");
    //signalInfo_vec.push_back(std::make_tuple( "C2V_0p01","C2V_0p01",(0.000000120285913353 / 0.00609678355976939), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_0p01/C2V_0p01_0.root"));  
    //signalInfo_vec.push_back(std::make_tuple( "C2V_0p5","C2V_0p5",(0.00030913696902 / 15.2653474807739), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_0p5/C2V_0p5_0.root"));  
    //signalInfo_vec.push_back(std::make_tuple( "C2V_1","C2V_1",(0.00122269524783 / 60.9941749572754), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_1/C2V_1_0.root"));   
    //signalInfo_vec.push_back(std::make_tuple( "C2V_ucsd_3","C2V_ucsd_3",(0.005496 / 3433.06640625), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_3/C2V_ucsd_3_0.root"));  
    signalInfo_vec.push_back(std::make_tuple( "C2V_ucsd_4","C2V_ucsd_4",(0.01148 / 7175.4091796875), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_4/C2V_ucsd_4_0.root")); 
    //signalInfo_vec.push_back(std::make_tuple( "C2V_ucsd_4p5","C2V_ucsd_4p5",(0.01533 / 9567.375), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_4p5/C2V_ucsd_4p5_0.root")); 
    //signalInfo_vec.push_back(std::make_tuple( "C2V_ucsd_m1","C2V_ucsd_m1",(0.005189 / 3239.87109375), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_m1/C2V_ucsd_m1_0.root")); 
    //signalInfo_vec.push_back(std::make_tuple( "C2V_ucsd_m2","C2V_ucsd_m2",(0.01114 / 6955.7578125), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_m2/C2V_ucsd_m2_0.root")); 
    /*
    backgroundInfo.push_back(std::make_tuple("SM","SM",(0.00045489 / 115.632980), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/SM/SM_0.root"));
    backgroundInfo.push_back(std::make_tuple("DYJets","DYJets",(6077.22 / 103073120.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/DYJets/DYJets_Full.root"));
    backgroundInfo.push_back(std::make_tuple("ttH","ttH",(0.6*0.584 / 1000000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/ttH/ttH_Full.root"));
    backgroundInfo.push_back(std::make_tuple("WJets","WJets",(52850.0 / 6720919175168.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WJets/WJets_Full.root"));
    backgroundInfo.push_back(std::make_tuple("WW","WW",(75.95 / 7959266.5), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WW/WW_Full.root"));
    backgroundInfo.push_back(std::make_tuple("WZ","WZ",(27.59 / 4000000.000000), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WZ/WZ_Full.root"));
    backgroundInfo.push_back(std::make_tuple("ZZ","ZZ",(12.17 / 2000000.000000), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/ZZ/ZZ_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_300to470", "QCD", (6826.0 / 57616400.0),"/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_300to470/QCD_Pt_300to470_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_470to600", "QCD",(552.6 / 27343560.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_470to600/QCD_Pt_470to600_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_600to800", "QCD",(156.6 / 67447600), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_600to800/QCD_Pt_600to800_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_800to1000", "QCD",(26.25 / 36340300.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_800to1000/QCD_Pt_800to1000_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_1000to1400", "QCD",(7.5 / 19397100.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1000to1400/QCD_Pt_1000to1400_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_1400to1800", "QCD",(0.6479 / 5892300.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1400to1800/QCD_Pt_1400to1800_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_1800to2400", "QCD",(0.08715 / 2990400.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1800to2400/QCD_Pt_1800to2400_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_2400to3200", "QCD",(0.005242 / 1992800.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_2400to3200/QCD_Pt_2400to3200_Full.root"));
    backgroundInfo.push_back(std::make_tuple("QCD_Pt_3200toInf", "QCD",(0.0001349 / 722400.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/DYJets/DYJets_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTbb_TTToHadronic","TTbb",(5.5 / 160536624.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTbb_TTToHadronic/TTbb_TTToHadronic_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTbb_TTTo2L2Nu","TTbb",(4.0 / 22349624.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTbb_TTTo2L2Nu/TTbb_TTTo2L2Nu_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTZToLLNuNu", "TTZ",(0.2432 / 3402735.5), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTZToLLNuNu/TTZToLLNuNu_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTZToQQ", "TTZ",(0.5104 / 376016.19), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTZToQQ/TTZToQQ_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTWJetsToLNu", "TTWJets",(0.2161 / 3548132.75), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTWJetsToLNu/TTWJetsToLNu_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTWJetsToQQ", "TTWJets",(0.4377 / 655903.81), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTWJetsToQQ/TTWJetsToQQ_Full.root"));
    backgroundInfo.push_back(std::make_tuple("TTWets", "TTJets",(831.76 / 298162651136), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTJets/TTJets_Full.root"));
    */
    backgroundInfo.push_back(std::make_tuple( "DYJets", "DYJets",(6077.22 / 96233336.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/DYJets/DYJets_Full.root"));
    
    //
    backgroundInfo.push_back(std::make_tuple( "ttH", "ttH",(0.6*0.5269 / 4831113.5), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/ttH/ttH_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_300to470", "QCD",(6826.0 / 57868000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_300to470/QCD_Pt_300to470_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_470to600", "QCD",(552.6 / 52448116.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_470to600/QCD_Pt_470to600_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_600to800", "QCD",(156.6 / 66914000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_600to800/QCD_Pt_600to800_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_800to1000", "QCD",(26.25 / 36830000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_800to1000/QCD_Pt_800to1000_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_1000to1400", "QCD",(7.5 / 19664000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1000to1400/QCD_Pt_1000to1400_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_1400to1800", "QCD",(0.6479 / 10982000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1400to1800/QCD_Pt_1400to1800_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_1800to2400", "QCD",(0.08715 / 5491000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_1800to2400/QCD_Pt_1800to2400_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_2400to3200", "QCD",(0.005242 / 2931000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_2400to3200/QCD_Pt_2400to3200_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "QCD_Pt_3200toInf", "QCD",(0.0001349 / 1000000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/QCD_Pt_3200toInf/QCD_Pt_3200toInf_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "WJets", "WJets",(52850.0 / 1186012928.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WJets/WJets_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "WW", "WW",(75.95 / 15679123.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WW/WW_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "WZ", "WZ",(27.59 / 7940000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/WZ/WZ_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "ZZ", "ZZ",(12.17 / 3526000.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/ZZ/ZZ_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTZToLLNuNu", "TTZ",(0.2432 / 4793524.5), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTZToLLNuNu/TTZToLLNuNu_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTZToQQ", "TTZ",(0.5104 / 10140212.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTZToQQ/TTZToQQ_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTWJetsToLNu", "TTWJets",(0.2161 / 3520822.75), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTWJetsToLNu/TTWJetsToLNu_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTWJetsToQQ", "TTWJets",(0.4377 / 657418.8125), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTWJetsToQQ/TTWJetsToQQ_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTbb_TTToHadronic", "TTbb",(5.5 / 160721072.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTbb_TTToHadronic/TTbb_TTToHadronic_Full.root"));
    
    backgroundInfo.push_back(std::make_tuple( "TTbb_TTTo2L2Nu", "TTbb",(4.0 / 22022552.0), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTbb_TTTo2L2Nu/TTbb_TTTo2L2Nu_Full.root"));
    
    //backgroundInfo.push_back(std::make_tuple( "SM", "SM",(0.00045489 / 115.632980), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/SM/SM_0.root"));
    backgroundInfo.push_back(std::make_tuple( "SM", "SM",(0.0005609 / 3433.06640625), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/C2V_ucsd_1/C2V_ucsd_1_0.root"));
    backgroundInfo.push_back(std::make_tuple( "TTJets", "TTJets",(831.76 / 298162651136), "/eos/cms/store/user/lzygala/HVV/Selection_TTrees/Samples_VBFdEta5p5_HiggsPT350_ParticleNet0p9_LepLT800_SLepLT1000/TTJets/TTJets_Full.root"));
    
    nSamples=backgroundInfo.size();

    ReadInTrees();
    SaveHistograms();
}
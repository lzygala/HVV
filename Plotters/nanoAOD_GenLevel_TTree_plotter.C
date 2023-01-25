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
    //gPad->SetLogy();

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
    //gPad->SetLogy();

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

void DrawPlot_nobackground(string varname, string title, string outputName, string selection = "", int bins = 12, float xmin = 0.0, float xmax = 2000.0){

    gStyle->SetOptStat(0);
    string outdir = "../Output/histos/";
    setOutput(outdir);

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);
    float scale_factor = 137.65*1000;


    std::vector<TH1F*> h_Signal;
    for(int i=0; i<signalInfo_vec.size(); i++){
        h_Signal.push_back(new TH1F(("h_Signal_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), title.c_str(), bins, xmin, xmax));
        Signal_trees_vec.at(i)->Project(("h_Signal_"+get<0>(signalInfo_vec.at(i))+outputName).c_str(), varname.c_str(), selection.c_str());
        //h_Signal.at(i)->Scale(get<2>(signalInfo_vec.at(i))*scale_factor);
        if(h_Signal.at(i)->Integral() != 0) h_Signal.at(i)->Scale(1/h_Signal.at(i)->Integral());
    }    



    std::vector<float> maxima;
    for(int i=0; i<signalInfo_vec.size(); i++){
        maxima.push_back( h_Signal[i]->GetMaximum()); 
    }
    std::sort(maxima.begin(),maxima.end()); 

    TLegend* legend = new TLegend(0.75, 0.7, 0.9, 0.9);
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

    ///gPad->SetLogy();

    cout<<"CCC"<<endl;

    int color_counter = 0;
    for(int i=0; i<h_Signal.size(); i++){

        int color_sig = colors[color_counter++] + 1;
        h_Signal[i] -> SetLineColor(color_sig);
        legend -> AddEntry(h_Signal[i], (get<0>(signalInfo_vec.at(i))).c_str(), "F");
        h_Signal[i]->Draw("hist,sames");
        gPad -> Update();
    }

    int color = colors[color_counter++] + 1;


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

        Background_trees.push_back((TTree*)tree_infile[i]->Get("Events"));
    }

    for(int i=0; i<signalInfo_vec.size(); i++){
        tree_signal_infile[i] = TFile::Open(get<3>(signalInfo_vec[i]).c_str());
        cout<<cur_time()<<"Reading in File: "<< (get<3>(signalInfo_vec[i]))<<endl;

        Signal_trees_vec.push_back((TTree*)tree_signal_infile[i]->Get("Events"));
    }

    cout<<"Finished Reading in Trees"<<endl;
}

void SaveHistograms(){

    DrawPlot_nobackground("GenPart_pt", "Higgs GenParticle PT", "H_GenPart_pt", "genWeight*(GenMatching_Initial_FullMatch && Higgs_GenPart_idx != -1 && Iteration$==Higgs_GenPart_idx)",100,0,2000);
    DrawPlot_nobackground("GenPart_pt", "Z GenParticle PT", "Z_GenPart_pt", "genWeight*(GenMatching_Initial_FullMatch && Z_GenPart_idx != -1 && Iteration$==Z_GenPart_idx)",100,0,2000);
    DrawPlot_nobackground("GenPart_pt", "W GenParticle PT", "W_GenPart_pt", "genWeight*(GenMatching_Initial_FullMatch && ((Wp_GenPart_idx != -1 && Iteration$==Wp_GenPart_idx) || (Wm_GenPart_idx != -1 && Iteration$==Wm_GenPart_idx)))",100,0,2000);
    DrawPlot_nobackground("GenPart_pt", "Lead VBF GenParticle PT", "VBFL_GenPart_pt", "genWeight*(GenMatching_Initial_FullMatch && VBFLead_GenPart_idx != -1 && Iteration$==VBFLead_GenPart_idx)",100,0,500);
    DrawPlot_nobackground("GenPart_pt", "Trail VBF GenParticle PT", "VBFT_GenPart_pt", "genWeight*(GenMatching_Initial_FullMatch && VBFTrail_GenPart_idx != -1 && Iteration$==VBFTrail_GenPart_idx)",100,0,500);
    DrawPlot_nobackground("VBF_GenPart_Mjj", "VBF GenParticle Mjj", "VBF_GenPart_Mjj", "genWeight*(GenMatching_Initial_FullMatch)",100,0,6000);
    DrawPlot_nobackground("VBF_GenPart_dEta", "VBF GenParticle dEta", "VBF_GenPart_dEta", "genWeight*(GenMatching_Initial_FullMatch)",100,0,10);
    DrawPlot_nobackground("VBF_GenPart_dR", "VBF GenParticle dR", "VBF_GenPart_dR", "genWeight*(GenMatching_Initial_FullMatch)",100,0,10);


    DrawPlot_nobackground("FatJet_pt", "Higgs AK8 PT", "H_FatJet_pt", "genWeight*(POMatching_FullMatch && Higgs_Matched_AK8_idx != -1 && Iteration$==Higgs_Matched_AK8_idx)",100,0,2000);
    DrawPlot_nobackground("FatJet_particleNet_HbbvsQCD", "Higgs AK8 PNScore HbbvsQCD", "H_FatJet_particleNet_HbbvsQCD", "genWeight*(POMatching_FullMatch && Higgs_Matched_AK8_idx != -1 && Iteration$==Higgs_Matched_AK8_idx)",100,0,1);
    DrawPlot_nobackground("FatJet_pt", "Z AK8 PT", "Z_FatJet_pt", "genWeight*(POMatching_FullMatch && Z_Matched_AK8_idx != -1 && Iteration$==Z_Matched_AK8_idx)",100,0,2000);
    DrawPlot_nobackground("FatJet_particleNet_ZvsQCD", "Z AK8 PNScore ZvsQCD", "Z_FatJet_particleNet_ZvsQCD", "genWeight*(POMatching_FullMatch && Z_Matched_AK8_idx != -1 && Iteration$==Z_Matched_AK8_idx)",100,0,1);
    DrawPlot_nobackground("Jet_pt", "Lead VBF  Matched AK4 PT", "VBFL_Matched_AK4_pt", "genWeight*(POMatching_FullMatch && VBFLead_Matched_AK4_idx != -1 && Iteration$==VBFLead_Matched_AK4_idx)",100,0,500);
    DrawPlot_nobackground("Jet_pt", "Trail VBF  Matched AK4 PT", "VBFT_Matched_AK4_pt", "genWeight*(POMatching_FullMatch && VBFTrail_Matched_AK4_idx != -1 && Iteration$==VBFTrail_Matched_AK4_idx)",100,0,500);
    DrawPlot_nobackground("VBF_Matched_AK4_Mjj", "VBF  Matched AK4 Mjj", "VBF_Matched_AK4_Mjj", "genWeight*(POMatching_FullMatch)",100,0,6000);
    DrawPlot_nobackground("VBF_Matched_AK4_dEta", "VBF  Matched AK4 dEta", "VBF_Matched_AK4_dEta", "genWeight*(POMatching_FullMatch)",100,0,10);
    DrawPlot_nobackground("VBF_Matched_AK4_dR", "VBF Matched AK4 dR", "VBF_Matched_AK4_dR", "genWeight*(POMatching_FullMatch)",100,0,10);
    

}

void nanoAOD_GenLevel_TTree_plotter(){

    signalInfo_vec.push_back(std::make_tuple( "WWH_ucsd_C2V_0","WWH_ucsd_C2V_0",(0.00274 / 1054.77), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WWH_ucsd_C2V_0/WWH_ucsd_C2V_0_0.root"));
    signalInfo_vec.push_back(std::make_tuple( "WWH_ucsd_C2V_1","WWH_ucsd_C2V_1",(0.00115 / 336.548), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WWH_ucsd_C2V_1/WWH_ucsd_C2V_1_0.root")); 
    signalInfo_vec.push_back(std::make_tuple( "WWH_ucsd_C2V_4","WWH_ucsd_C2V_4",(0.01692 / 7175.41), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WWH_ucsd_C2V_4/WWH_ucsd_C2V_4_0.root"));  

    signalInfo_vec.push_back(std::make_tuple( "WZH_ucsd_C2V_0","WZH_ucsd_C2V_0",(0.00168 / 695.222), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WZH_ucsd_C2V_0/WZH_ucsd_C2V_0_0.root"));
    signalInfo_vec.push_back(std::make_tuple( "WZH_ucsd_C2V_1","WZH_ucsd_C2V_1",(0.006 / 236.081), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WZH_ucsd_C2V_1/WZH_ucsd_C2V_1_0.root")); 
    signalInfo_vec.push_back(std::make_tuple( "WZH_ucsd_C2V_4","WZH_ucsd_C2V_4",(0.01123 / 4754.75), "/eos/cms/store/user/lzygala/HVV/GenLevelInfo/18/WZH_ucsd_C2V_4/WZH_ucsd_C2V_4_0.root")); 
    
    nSamples=signalInfo_vec.size();

    ReadInTrees();
    SaveHistograms();
}
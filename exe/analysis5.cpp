#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>

#include "useful_functions.cpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

using namespace std;

void analysis5()
{   
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis1.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1110);
    gStyle->SetFitFormat("2.2e");
    

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    const char * path1 = "../data/magnetoresistenza.txt";
    int first_line = comment_lines(path1);
    ifstream file(path1);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH, sVH, B, sB;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path1) == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH.push_back(entry1);
            sVH.push_back(entry2);
            B.push_back(entry3);
            sB.push_back(entry4);
        }
    }
    else if(count_column(path1) == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH.push_back(entry1);
            sVH.push_back(entry2);
            B.push_back(entry3);
            sB.push_back(entry4);
        }
    }


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////

    const float V0 = 0.;
    const float sV0 = 0.1;
    const float I = 1.;
    const float sI = 0.1; 

    vector<float> R, sR;
    
    for(int i = 0; i < VH.size(); i++)   
    {
        entry1 = (VH.at(i) - V0)/I;
        entry2 = sqrt(pow(sVH.at(i)/I, 2) + pow(sV0/I, 2) + pow(sI*(VH.at(i)-V0)/(I*I), 2));
        
        R.push_back(entry1);
        sR.push_back(entry2);
    }

    string str1("\tR[Ohm]");
    string str2("\tsR[Ohm]");
    if(names.find(str1) == string::npos)
    {
        names += "\tR[Ohm]";
        append_column(path1, "R[Ohm]", R);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsR[Ohm]";   
        append_column(path1, "sR[Ohm]", sR);
    } 

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH.size(); i++)   
        cout << VH.at(i) << "\t" << sVH.at(i) << "\t" << B.at(i) << "\t" << sB.at(i) << "\t" << 
        R.at(i) << "\t\t" << sR.at(i) << endl;
    cout << endl;
    
    /////////////////////////////// FIT FORMA QUADRATICA //////////////////////////////////////////////////////
    
    TCanvas * canvas1 = new TCanvas("canvas1", "magnetoresistenza", 500, 5, 500, 600);
    canvas1->SetGrid();
    
    TF1 * tf1 = new TF1("tf1", "[0] + [1]*x*x", -15, 15);
    tf1->SetLineColor(38);
    tf1->SetParName(0, "R_{0}");

    TGraphErrors * graph1 = new TGraphErrors(B.size(), &B[0], &R[0], &sB[0], &sR[0]);
    graph1->SetTitle("#splitline{Magnetoresistenza}{R = R_{0} + p_{1} B^2};B [T];R [#Omega]");
    std_graph_settings(*graph1);
    
    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/magnetoresistenza1.jpg");
    canvas1->SaveAs("../graphs/magnetoresistenza1.pdf");
    
    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
    " (Probability: " << tf1->GetProb() << ")." << endl;

    /////////////////////////////// FIT FORMA FORNITA //////////////////////////////////////////////////////
    
    TCanvas * canvas2 = new TCanvas("canvas2", "magnetoresistenza", 500, 5, 500, 600);
    canvas1->SetGrid();
    
    TF1 * tf2 = new TF1("tf2", "[0]*(1+[1]*(x-[2]))", -15, 15);
    tf2->SetLineColor(38);
    tf2->SetParNames("R_{0}", "#mu^2", "B_{0}");

    TGraphErrors * graph2 = new TGraphErrors(B.size(), &B[0], &R[0], &sB[0], &sR[0]);
    graph2->SetTitle("#splitline{Magnetoresistenza}{R = R_{0} (1 + #mu^2 (B - B_{0})};B [T];R [#Omega]");
    std_graph_settings(*graph2);
    
    graph2->Fit(tf2, "ER");
    graph2->Draw("ap");
    canvas2->SaveAs("../graphs/magnetoresistenza2.jpg");
    canvas2->SaveAs("../graphs/magnetoresistenza2.pdf");
    
    cout << "Chi^2:" << tf2->GetChisquare() << ", number of DoF: " << tf2->GetNDF() << 
    " (Probability: " << tf2->GetProb() << ")." << endl;
    


}
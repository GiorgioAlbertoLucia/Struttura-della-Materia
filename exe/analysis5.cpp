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

    vector<float> VH, I, B, sVH, sB, sI;
    float entry1, entry2, entry3, entry4, entry5, entry6, entry7, entry8;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path1) == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH.push_back(entry1);
            I.push_back(entry2);
        }
    }
    else if(count_column(path1) == 8)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6 >> entry7 >> entry8)
        {
            VH.push_back(entry1);
            I.push_back(entry2);
        }
    }


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////

    const float V0 = 0.000000100441;    // mV
    const float sV0 = 0.00000213145;
    const float i = - 0.35;             // mA
    const float si = 0.03; 
    const float sB_sist = 4.205;        // mT

    // I -> B
    const float a_tilde = 0.20636;      // mT
    const float sa_tilde = 6.18027;
    const float b_tilde = 206.022;      // mT / A
    const float sb_tilde = 103.011;

    vector<float> R, sR;
    
    for(int j = 0; j < VH.size(); j++)   
    {
        entry1 = a_tilde + b_tilde * I.at(j);
        B.push_back(entry1);
        
        entry2 = (VH.at(j) - V0)/i;
        R.push_back(entry2);

        entry3 = abs(VH.at(j)) * 0.02;
        sVH.push_back(entry3);

        entry4 = 0.01;
        sI.push_back(entry4);

        entry5 = sqrt( (sa_tilde*sa_tilde + pow(b_tilde*sI.at(j), 2)) + sB_sist*sB_sist);
        //entry5 = sqrt( (sa_tilde*sa_tilde + pow(sb_tilde*I.at(j), 2) + pow(b_tilde*sI.at(j), 2)) + sB_sist*sB_sist);
        sB.push_back(entry5);

        entry6 = sqrt(pow(sVH.at(j)/i, 2) + pow(sV0/i, 2) + pow(si*(VH.at(j)-V0)/(i*i), 2));
        sR.push_back(entry6);
    }

    string str1("\tB[mT]"), str2("\tR[Ohm]"), str3("\tsVH[mV]"), str4("\tsI[A]"), str5("\tsB[mT]"), str6("\tsR[Ohm]");
    if(names.find(str1) == string::npos)
    {
        names += "\tB[mT]";
        append_column(path1, "B[mT]", B);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tR[Ohm]";   
        append_column(path1, "R[Ohm]", R);
    } 
    if(names.find(str3) == string::npos)
    {
        names += "\tsVH[mV]";
        append_column(path1, "sVH[mV]", sVH);
    }    
    if(names.find(str4) == string::npos)
    {
        names += "\tsI[A]";   
        append_column(path1, "sI[A]", sI);
    } 
    if(names.find(str5) == string::npos)
    {
        names += "\tsB[mT]";
        append_column(path1, "sB[mT]", sB);
    }    
    if(names.find(str6) == string::npos)
    {
        names += "\tsR[Ohm]";   
        append_column(path1, "sR[Ohm]", sR);
    } 
    

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int j = 0; j < VH.size(); j++)   
        cout << VH.at(j) << "\t" << I.at(j) << "\t" <<  B.at(j) << "\t"  << R.at(j) 
        << "\t" << sVH.at(j) << "\t" << sI.at(j) << "\t" << sB.at(j) << "\t\t" << sR.at(j) << endl;
    cout << endl;
    
    /////////////////////////////// FIT FORMA QUADRATICA //////////////////////////////////////////////////////
    
    TCanvas * canvas1 = new TCanvas("canvas1", "magnetoresistenza", 500, 5, 500, 600);
    canvas1->SetGrid();
    
    TF1 * tf1 = new TF1("tf1", "[0] + [1]*x*x", -400, 400);
    tf1->SetLineColor(38);
    tf1->SetParName(0, "R_{0}");

    TGraphErrors * graph1 = new TGraphErrors(B.size(), &B[0], &R[0], &sB[0], &sR[0]);
    graph1->SetTitle("#splitline{Magnetoresistenza}{R = R_{0} + p_{1} B^2};B [mT];R [#Omega]");
    std_graph_settings(*graph1);
    
    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/magnetoresistenza1.jpg");
    canvas1->SaveAs("../graphs/magnetoresistenza1.pdf");
    
    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
    " (Probability: " << tf1->GetProb() << ")." << endl;

    /////////////////////////////// FIT FORMA FORNITA //////////////////////////////////////////////////////
    
    TCanvas * canvas2 = new TCanvas("canvas2", "magnetoresistenza", 500, 5, 500, 600);
    canvas2->SetGrid();
    
    TF1 * tf2 = new TF1("tf2", "[0]*(1+[1]*(x-[2])*(x-[2]))", -400, 400);
    tf2->SetLineColor(38);
    tf2->SetParNames("R_{0}", "#mu^2", "B_{0}");

    TGraphErrors * graph2 = new TGraphErrors(B.size(), &B[0], &R[0], &sB[0], &sR[0]);
    graph2->SetTitle("#splitline{Magnetoresistenza}{R = R_{0} [1 + #mu^{2} (B - B_{0})^{2}]};B [mT];R [#Omega]");
    std_graph_settings(*graph2);
    
    graph2->Fit(tf2, "MR");
    graph2->Draw("ap");
    canvas2->SaveAs("../graphs/magnetoresistenza2.jpg");
    canvas2->SaveAs("../graphs/magnetoresistenza2.pdf");
    
    cout << "Chi^2:" << tf2->GetChisquare() << ", number of DoF: " << tf2->GetNDF() << 
    " (Probability: " << tf2->GetProb() << ")." << endl;
    


}
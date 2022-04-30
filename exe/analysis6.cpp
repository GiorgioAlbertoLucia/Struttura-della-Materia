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

void analysis_sample()
{  
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis1.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");
    
    /////////////////////////// CALIBRAZIONE ////////////////////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path1 = "../data/TempRiferimento.txt";
    int first_line = comment_lines(path1);
    ifstream file(path1);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n1');    

    vector<float> VH1, sVH1, T1, sT1;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    while (file >> entry1 >> entry2 >> entry3 >> entry4)
    {
        VH1.push_back(entry1);
        sVH1.push_back(entry2);
        T1.push_back(entry3);
        sT1.push_back(entry4);
    }
    

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH1.size(); i++)   
        cout << VH1.at(i) << "\t" << sVH1.at(i) << "\t" << T1.at(i) << "\t" << sT1.at(i) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    const int n1[] = {0, 10, 20};        // position where each dataset begins
    const int sets1 = 2;

    vector<float> chi, schi, omega, somega;

    TCanvas * canvas1 = new TCanvas("canvas1", "TempRiferimento", 500, 5, 500, 600);
    canvas1->SetGrid();
    canvas1->Divide(1, 2);
    
    for (int i = 0; i < sets1; i++)
    {
        vector<float> sub_VH1(VH1.begin()+n1[i], VH1.begin()+n1[i+1]);
        vector<float> sub_sVH1(sVH1.begin()+n1[i], sVH1.begin()+n1[i+1]);
        vector<float> sub_T1(T1.begin()+n1[i], T1.begin()+n1[i+1]);
        vector<float> sub_sT1(sT1.begin()+n1[i], sT1.begin()+n1[i+1]);

        TF1 * tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
        tf1->SetLineColor(38);

        TGraphErrors * graph1 = new TGraphErrors(sub_T1.size(), &sub_T1[0], &sub_VH1[0], &sub_sT1[0], &sub_sVH1[0]);
        graph1->SetTitle("#splitline{V_{H} vs T (B = 0)}{V_{H} = p_{0} + p_{1} T};T [°C];V_{H} [V]");
        std_graph_settings(*graph1);
    
        canvas1->cd(i+1);
        graph1->Fit(tf1, "ER");
        graph1->Draw("ap");
    
        cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
        " (Probability: " << tf1->GetProb() << ")." << endl;

        chi.push_back(tf1->GetParameter(0));
        schi.push_back(tf1->GetParError(0));
        omega.push_back(tf1->GetParameter(1));
        somega.push_back(tf1->GetParError(1));
    }
    canvas1->SaveAs("../graphs/tempriferimento.jpg");
    canvas1->SaveAs("../graphs/tempriferimento.pdf");

    
    














   ////////////////////////// MISURE CON B = const ////////////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path2 = "../data/TempBconst.txt";
    int first_line = comment_lines(path2);
    ifstream file(path2);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n1');    

    vector<float> VH2, sVH2, T2, sT2;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path2) == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            T2.push_back(entry3);
            sT2.push_back(entry4);
        }
    }
    else if(count_column(path2) == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            T2.push_back(entry3);
            sT2.push_back(entry4);
        }
    }
    


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    const float I = 8.;      // mA
    const float sI = 0.;

    const int n2[] = {0, 10, 20};

    vector<float> VH2_correct, sVH2_correct;

    for (int i = 0; i < VH.size(); i++)   
    {
        int j;
        if(i < n2[1])       j = 0;
        else if(i < n2[2])  j = 1;

        entry1 = VH2.at(i) - (chi.at(j) + omega.at(j)*I);
        entry2 = sqrt(sVH2.at(i)*sVH2.at(i) + schi.at(j)*schi.at(j) + I*I*somega.at(j)*somega.at(j) + omega.at(j)*omega.at(j)*sI*sI);

        VH2_correct.push_back(entry1);
        sVH2_correct.push_back(entry2);
    }

    string str1("\tVH_correct[V]");
    string str2("\tsVH_correct[V]");
    if(names.find(str1) == string::npos)
    {
        names += "\tVH_correct[V]";
        append_column(path1, "VH_correct[V]", VH2_correct);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH_correct[V]";   
        append_column(path1, "sVH_correct[V]", sVH2_correct);
    } 

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH.size(); i++)   
        cout << VH2.at(i) << "\t" << sVH2.at(i) << "\t" << T2.at(i) << "\t" << sT2.at(i) << "\t" << 
        VH2_correct.at(i) << "\t\t" << sVH2_correct.at(i) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas2 = new TCanvas("canvas2", "TempRiferimento", 500, 5, 500, 600);
    canvas2->SetGrid();
    canvas2->Divide(1, 2);
    
    for (int i = 0; i < sets1; i++)
    {
        vector<float> sub_VH2_correct(VH2_correct.begin()+n2[i], VH2_correct.begin()+n2[i+1]);
        vector<float> sub_sVH2_correct(sVH2_correct.begin()+n2[i], sVH2_correct.begin()+n2[i+1]);
        vector<float> sub_T2(T2.begin()+n2[i], T2.begin()+n2[i+1]);
        vector<float> sub_sT2(sT2.begin()+n1[i], sT2.begin()+n2[i+1]);

        TF1 * tf2 = new TF1("tf2", "[0]+[1]*x", -15, 15);
        tf2->SetLineColor(38);

        TGraphErrors * graph2 = new TGraphErrors(sub_T2.size(), &sub_T2[0], &sub_VH2_correct[0], &sub_sT2[0], &sub_sVH2_correct[0]);
        graph2->SetTitle("#splitline{V_{H} vs T (B = 0)}{V_{H} = p_{0} + p_{1} T};T [°C];V_{H} [V]");
        std_graph_settings(*graph2);
    
        canvas2->cd(i+1);
        graph2->Fit(tf2, "ER");
        graph2->Draw("ap");
    
        cout << "Chi^2:" << tf2->GetChisquare() << ", number of DoF: " << tf2->GetNDF() << 
        " (Probability: " << tf2->GetProb() << ")." << endl;

        chi.push_back(tf1->GetParameter(0));
        schi.push_back(tf1->GetParError(0));
        omega.push_back(tf1->GetParameter(1));
        somega.push_back(tf1->GetParError(1));
    }
    canvas2->SaveAs("../graphs/tempriferimento.jpg");
    canvas2->SaveAs("../graphs/tempriferimento.pdf");


}
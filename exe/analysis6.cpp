#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <Math/Interpolator.h>

#include "useful_functions.cpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

using namespace std;

void analysis6()
{  
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis1.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");
    
    /////////////////////////// INTERPOLAZIONE ////////////////////////////////////////////////

    const char * path0 = "../data/TempInterpolazione.txt";
    int first_line = comment_lines(path0);
    ifstream file(path0);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');

    vector<double> T0, R0; 
    string names;
    float entry0, entry1, entry2, entry3, entry4, entry5, entry6, entry7, entry8, entry9, entry10;
    getline(file,names);
    while(file >> entry0 >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6 >> entry7 >> entry8 >> entry9 >> entry10)
    {
        T0.push_back(entry0);
        R0.push_back(entry1);

        T0.push_back(entry0+1);
        R0.push_back(entry2);

        T0.push_back(entry0+2);
        R0.push_back(entry3);

        T0.push_back(entry0+3);
        R0.push_back(entry4);

        T0.push_back(entry0+4);
        R0.push_back(entry5);

        T0.push_back(entry0+5);
        R0.push_back(entry6);

        T0.push_back(entry0+6);
        R0.push_back(entry7);

        T0.push_back(entry0+7);
        R0.push_back(entry8);

        T0.push_back(entry0+8);
        R0.push_back(entry9);

        T0.push_back(entry0+9);
        R0.push_back(entry10);
    }

    cout << "check " << R0.size() << " " << T0.size() << endl;

    ROOT::Math::Interpolator interpolator(R0.size(), ROOT::Math::Interpolation::kCSPLINE);
    interpolator.SetData(R0.size(), &R0[0], &T0[0]);



    /////////////////////////// CALIBRAZIONE ////////////////////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path1 = "../data/TempRiferimento.txt";
    first_line = comment_lines(path1);
    file.open(path1);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH1, R1;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path1)==2)
    {
        while (file >> entry1 >> entry2)
        {
            VH1.push_back(entry1);
            R1.push_back(entry2);
        }
    }
    if(count_column(path1)==6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH1.push_back(entry1);          // mV
            R1.push_back(entry2);           // Ω
        }
    }

    cout << "check " << VH1.size() << endl;
    
    /////////////////////////////// ADD DATA ///////////////////////////////////////////////

    vector<float> T1, sVH1, sR1, sT1;
    
    for(int i = 0; i < VH1.size(); i++)   
    {
        entry1 = interpolator.Eval(R1.at(i)) + 273.15;      // K
        T1.push_back(entry1);
        
        entry2 = abs(VH1.at(i)) * 0.02;
        sVH1.push_back(entry2);

        entry3 = abs(R1.at(i)) * 0.008 + 0.1;
        sR1.push_back(entry3);

        entry4 = interpolator.Deriv(R1.at(i)) * sR1.at(i);
        sT1.push_back(entry4);
    }

    string str1("\tT[K]"), str2("\tsVH[mV]"), str3("\tsR[Ω]"), str4("\tsT[K]");
    if(names.find(str1) == string::npos)
    {
        names += "\tT[K]";
        append_column(path1, "T[K]", T1);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH[mV]";   
        append_column(path1, "sVH[mV]", sVH1);
    } 
    if(names.find(str3) == string::npos)
    {
        names += "\tsR[Ω]";
        append_column(path1, "sR[Ω]", sR1);
    }    
    if(names.find(str4) == string::npos)
    {
        names += "\tsT[K]";   
        append_column(path1, "sT[K]", sT1);
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH1.size(); i++)   
        cout << VH1.at(i) << "\t" << R1.at(i) << "\t" << T1.at(i) << "\t" 
             << sVH1.at(i) << "\t" << sR1.at(i) << "\t" << sT1.at(i) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    const int n1[] = {0, 22, 43};        // position where each dataset begins
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
        graph1->SetTitle("#splitline{V_{H} vs T (B = 0)}{V_{H} = p_{0} + p_{1} T};T [K];V_{H} [mV]");
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


    cout << "Dati Riferimento (V_H = chi + omega * T): " << endl;
    for(int i = 0; i < chi.size(); i++)
    {
        cout << "chi = (" << chi.at(i) << " ± " << schi.at(i) << ") mV" << endl;
        cout << "omega = (" << omega.at(i) << " ± " << somega.at(i) << ") mV/K" << endl;
    }
    
    const float chi_medio = (chi.at(0)+chi.at(1))/2;
    const float schi_medio = sqrt(schi.at(0)*schi.at(0) + schi.at(1)*schi.at(1));

    const float omega_medio = (omega.at(0)+omega.at(1))/2;
    const float somega_medio = sqrt(somega.at(0)*somega.at(0) + somega.at(1)*somega.at(1));














   ////////////////////////// MISURE CON B = const ////////////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path2 = "../data/TempBconst.txt";
    first_line = comment_lines(path2);
    file.open(path2);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH2, R2;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path2) == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH2.push_back(entry1);
            R2.push_back(entry2);
        }
    }
    else if(count_column(path2) == 8)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6 >> entry7 >> entry8)
        {
            VH2.push_back(entry1);
            R2.push_back(entry2);
        }
    }
    


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    const float I = 1.52;                           // A
    const float sI = 0.01;

    const int n2[] = {0, 21, 41};

    vector<float> T2, VH2_correct, sVH2, sR2, sT2, sVH2_correct;

    for (int i = 0; i < VH2.size(); i++)   
    {
        entry1 = interpolator.Eval(R2.at(i));       
        T2.push_back(entry1);                       // K

        entry2 = VH2.at(i) - (chi_medio + omega_medio*I);
        VH2_correct.push_back(entry2);              // mV

        entry3 = VH2.at(i) * 0.02;
        sVH2.push_back(entry3);

        entry4 = R2.at(i) * 0.008 + 0.1;
        sR2.push_back(entry4);

        entry5 = interpolator.Deriv(R2.at(i)) * sR2.at(i);
        sT2.push_back(entry5);

        entry6 = sqrt(pow(sVH2.at(i),2) + pow(schi_medio,2) + pow(I*somega_medio,2) + pow(omega_medio*sI,2));
        sVH2_correct.push_back(entry6);
    }

    str1 = "\tT[K]";
    str2 = "\tVH_corr[mV]";
    str3 = "\tsVH[mV]";
    str4 = "\tsR[Ω]";
    string str5("\tsT[K]"), str6("\tsVH_correct[mV]");
    if(names.find(str1) == string::npos)
    {
        names += str1;
        append_column(path1, "\tT[K]]", T2);
    }    
    if(names.find(str2) == string::npos)
    {
        names += str2;
        append_column(path1, "\tVH_correct[mV]", VH2_correct);
    }    
    if(names.find(str3) == string::npos)
    {
        names += str3;
        append_column(path1, "\tsVH[mV]", sVH2);
    }    
    if(names.find(str4) == string::npos)
    {
        names += str4;
        append_column(path1, "\tsR[Ω]", sR2);
    }    
    if(names.find(str5) == string::npos)
    {
        names += str5;
        append_column(path1, "\tsT[K]]", sT2);
    }    
    if(names.find(str6) == string::npos)
    {
        names += str6;   
        append_column(path1, "\tsVH_correct[mV]", sVH2_correct);
    } 

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH2.size(); i++)   
        cout << VH2.at(i) << "\t" << R2.at(i) << "\t" << T2.at(i) << "\t" << VH2_correct.at(i) << "\t" 
            << sVH2.at(i) << "\t" << sR2.at(i) << "\t" << sT2.at(i) << "\t" << sVH2_correct.at(i) << endl;
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

    }
    canvas2->SaveAs("../graphs/tempriferimento.jpg");
    canvas2->SaveAs("../graphs/tempriferimento.pdf");


}
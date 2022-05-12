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
#include <algorithm>

using namespace std;

void std_graph_settings_adv(TGraph& graph)
{
    gPad->SetTopMargin(0.15);
    graph.SetMarkerStyle(21);
    graph.SetMarkerColor(1);
    graph.SetMarkerSize(0.1);
    graph.SetLineColor(38);
    graph.SetLineWidth(4);
}

void analysis6()
{  
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis6.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1110);
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

    file.close();

    ROOT::Math::Interpolator interpolator(R0.size(), ROOT::Math::Interpolation::kCSPLINE);
    interpolator.SetData(R0.size(), &R0[0], &T0[0]);

    TCanvas * canvas0 = new TCanvas("canvas0", "TempInterpolazione", 500, 5, 500, 600);
    canvas0->SetGrid();

    TGraph * graph0 = new TGraph(R0.size(), &R0[0], &T0[0]);
    graph0->SetTitle("Interpolazione;R [#Omega];T [#circ C]");
    std_graph_settings(*graph0);
    graph0->SetLineColor(38);
    graph0->SetMarkerSize(0.2);
    graph0->SetMarkerColor(2);
    graph0->SetLineWidth(7);
    
    graph0->Draw("apl");
    canvas0->SaveAs("../graphs/interpolazione.pdf");

























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
    else if(count_column(path1)==6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH1.push_back(entry1);          // mV
            R1.push_back(entry2);           // Ω
        }
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();


    
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
        append_column(path1, "\tT[K]", T1);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH[mV]";   
        append_column(path1, "\tsVH[mV]", sVH1);
    } 
    if(names.find(str3) == string::npos)
    {
        names += "\tsR[Ω]";
        append_column(path1, "\tsR[Ω]", sR1);
    }    
    if(names.find(str4) == string::npos)
    {
        names += "\tsT[K]";   
        append_column(path1, "\tsT[K]", sT1);
    }

    
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < VH1.size(); i++)   
        cout << VH1.at(i) << "\t" << R1.at(i) << "\t" << T1.at(i) << "\t" 
             << sVH1.at(i) << "\t" << sR1.at(i) << "\t" << sT1.at(i) << endl;
    cout << endl;
    
    

    /////////////////////////////// INTERPOLAZIONE //////////////////////////////////////
    
    const int n1[] = {0, 22, 43};        // position where each dataset begins
    const int sets1 = 2;

    vector<float> chi, schi, omega, somega;

    TCanvas * canvas1 = new TCanvas("canvas1", "TempRiferimento", 500, 5, 500, 600);
    canvas1->Divide(1, 2);
    
    vector<double> sub1_VH1(VH1.begin(), VH1.begin()+22);
    vector<double> sub1_sVH1(sVH1.begin(), sVH1.begin()+22);
    vector<double> sub1_T1(T1.begin(), T1.begin()+22);
    vector<double> sub1_sT1(sT1.begin(), sT1.begin()+22);

    ROOT::Math::Interpolator interpolator1(sub1_VH1.size(), ROOT::Math::Interpolation::kCSPLINE);
    interpolator1.SetData(sub1_VH1.size(), &sub1_T1[0], &sub1_VH1[0]);


    TGraph * graph1 = new TGraph(sub1_T1.size(), &sub1_T1[0], &sub1_VH1[0]);
    graph1->SetTitle("Interpolazione V_{H};T [K];V_{H} [mV]");
    std_graph_settings_adv(*graph1);
    graph1->SetLineColor(38);

    canvas1->cd(1);
    gPad->SetGrid();
    graph1->Draw("apl");

//----------

    vector<double> sub2_VH1(VH1.begin()+22, VH1.end());
    vector<double> sub2_sVH1(sVH1.begin()+22, sVH1.end());
    vector<double> sub2_T1(T1.begin()+22, T1.end());
    vector<double> sub2_sT1(sT1.begin()+22, sT1.end());

    reverse(sub2_T1.begin(), sub2_T1.end());
    reverse(sub2_VH1.begin(), sub2_VH1.end());

    ROOT::Math::Interpolator interpolator2(sub2_T1.size(), ROOT::Math::Interpolation::kCSPLINE);
    interpolator2.SetData(sub2_T1.size(), &sub2_T1[0], &sub2_VH1[0]);

    TGraph * graph1_ = new TGraph(sub1_T1.size(), &sub1_T1[0], &sub1_VH1[0]);
    graph1_->SetTitle("Interpolazione V_{H};T [K];V_{H} [mV]");
    std_graph_settings_adv(*graph1_);
    graph1_->SetLineColor(38);

    canvas1->cd(2);
    gPad->SetGrid();
    graph1_->Draw("apl");

    canvas1->SaveAs("../graphs/interpolazione_VH1.png");
    canvas1->SaveAs("../graphs/interpolazione_VH1.pdf");
















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
    else if(count_column(path2) == 10)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6 >> entry7 >> entry8 >> entry9 >> entry10)
        {
            VH2.push_back(entry1);
            R2.push_back(entry2);
        }
    }
    


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    const float I = 1.01;                           // A
    const float sI = 0.01;
    const float B = 0.275176 + 205.948 * I;         // mT
    const float sB = sqrt((pow(6.18651, 2) + pow(0.281825*I, 2) + pow(205.948*sI, 2)) + 1.67*1.67);
    const float i_const = 8.00;                     // mA
    const float si_const = i_const * 0.006 + 0.03;
    const float t = 0.1;                            // cm

    const int n2[] = {0, 21, 41};

    cout << "check " << VH2.size();

    vector<float> T2, VH2_correct, RH2, sVH2, sR2, sT2, sVH2_correct, sRH2;

    for (int i = 0; i < VH2.size(); i++)   
    {
        entry1 = interpolator.Eval(R2.at(i)) + 273.15;       
        T2.push_back(entry1);                       // K

        if(i < 22)  entry2 = VH2.at(i) - interpolator1.Eval(T2.at(i));
        else        entry2 = VH2.at(i) - interpolator2.Eval(T2.at(i));
        VH2_correct.push_back(entry2);              // mV

        entry3 = VH2.at(i) * t * 1000/(i_const * B);    // V cm (A T)
        RH2.push_back(entry3);

        entry4 = VH2.at(i) * 0.02;
        sVH2.push_back(entry4);

        entry5 = R2.at(i) * 0.008 + 0.1;
        sR2.push_back(entry5);

        entry6 = interpolator.Deriv(R2.at(i)) * sR2.at(i);
        sT2.push_back(entry6);

        if(i < 22)  entry7 = sqrt(pow(sVH2.at(i),2) + pow(sT2.at(i)*interpolator1.Deriv(T2.at(i)),2));
        else        entry7 = sqrt(pow(sVH2.at(i),2) + pow(sT2.at(i)*interpolator2.Deriv(T2.at(i)),2));
        sVH2_correct.push_back(entry7);

        entry8 = sqrt(pow(sVH2.at(i)*t*1000/(i_const*B), 2) + pow(sB*VH2.at(i)*t*1000/(i_const*B*B), 2) + 
                        pow(VH2.at(i)*si_const*t*1000/(i_const*i_const*B), 2));
        sRH2.push_back(entry8);
    }

    str1 = "\tT[K]";
    str2 = "\tVH_correct[mV]";
    str3 = "\tRH[Vcm/AT]";
    str4 = "\tsVH[mV]";
    string str5("\tsR[Ω]"), str6("\tsT[K]"), str7("\tsVH_correct[mV]"), str8("\tsRH[Vcm/AT]");
    if(names.find(str1) == string::npos)
    {
        names += str1;
        append_column(path2, "\tT[K]", T2);
    }    
    if(names.find(str2) == string::npos)
    {
        names += str2;
        append_column(path2, "\tVH_correct[mV]", VH2_correct);
    }    
    if(names.find(str3) == string::npos)
    {
        names += str3;   
        append_column(path2, "\tRH[Vcm/AT]", RH2);
    }
    if(names.find(str4) == string::npos)
    {
        names += str4;
        append_column(path2, "\tsVH[mV]", sVH2);
    }    
    if(names.find(str5) == string::npos)
    {
        names += str5;
        append_column(path2, "\tsR[Ω]", sR2);
    }    
    if(names.find(str6) == string::npos)
    {
        names += str6;
        append_column(path2, "\tsT[K]", sT2);
    }    
    if(names.find(str7) == string::npos)
    {
        names += str7;   
        append_column(path2, "\tsVH_correct[mV]", sVH2_correct);
    } 
    if(names.find(str8) == string::npos)
    {
        names += str8;   
        append_column(path2, "\tsRH[Vcm/AT]", sRH2);
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
    
    TCanvas * canvas2 = new TCanvas("canvas2", "TempBconst", 500, 5, 500, 600);
    canvas2->Divide(1, 2);
    
    for (int i = 0; i < sets1; i++)
    {
        
        vector<float> sub_VH2_correct(VH2_correct.begin()+n2[i], VH2_correct.begin()+n2[i+1]);
        vector<float> sub_sVH2_correct(sVH2_correct.begin()+n2[i], sVH2_correct.begin()+n2[i+1]);
        vector<float> sub_T2(T2.begin()+n2[i], T2.begin()+n2[i+1]);
        vector<float> sub_sT2(sT2.begin()+n1[i], sT2.begin()+n2[i+1]);

        TGraphErrors * graph2 = new TGraphErrors(sub_T2.size(), &sub_T2[0], &sub_VH2_correct[0], &sub_sT2[0], &sub_sVH2_correct[0]);
        graph2->SetTitle("V_{H} vs T (B = const);T [K];V_{H} [mV]");
        std_graph_settings_adv(*graph2);
        graph2->SetLineColor(38);
    
        canvas2->cd(i+1);
        gPad->SetGrid();
        graph2->Draw("apl");
    }
    canvas2->SaveAs("../graphs/tempbconst.jpg");
    canvas2->SaveAs("../graphs/tempbconst.pdf");

    // RH coefficient
    
    TCanvas * canvas3 = new TCanvas("canvas3", "RHBconst", 500, 5, 500, 600);
    canvas3->Divide(1, 2);
    
    for (int i = 0; i < sets1; i++)
    {
        vector<float> sub_RH2_correct(RH2.begin()+n2[i], RH2.begin()+n2[i+1]);
        vector<float> sub_sRH2_correct(sRH2.begin()+n2[i], sRH2.begin()+n2[i+1]);
        vector<float> sub_T2(T2.begin()+n2[i], T2.begin()+n2[i+1]);
        vector<float> sub_sT2(sT2.begin()+n1[i], sT2.begin()+n2[i+1]);

        TGraphErrors * graph3 = new TGraphErrors(sub_T2.size(), &sub_T2[0], &sub_RH2_correct[0], &sub_sT2[0], &sub_sRH2_correct[0]);
        graph3->SetTitle("R_{H} vs T;T [K];R_{H} [#frac{V cm}{A T}]");
        std_graph_settings_adv(*graph3);
        graph3->SetLineColor(38);
        graph3->SetMarkerStyle(1);
        
    
        canvas3->cd(i+1);
        gPad->SetGrid();
        graph3->Draw("apl");
    }
    canvas3->SaveAs("../graphs/rh_va_temp.jpg");
    canvas3->SaveAs("../graphs/rh_vs_temp.pdf");


}
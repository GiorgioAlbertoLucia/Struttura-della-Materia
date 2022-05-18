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
#include <numeric>

using namespace std;

void std_graph_settings_adv(TGraph& graph)
{
    gPad->SetTopMargin(0.15);
    graph.SetMarkerColor(2);
    graph.SetMarkerStyle(1);
    graph.SetMarkerSize(0.5);
    graph.SetLineColor(38);
    graph.SetLineWidth(4);
}

void analysis3()
{

    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis3.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1110);
    gStyle->SetFitFormat("2.2e");

    ///////////////////// PART 1: REFERENCE MEASUREMENTS //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////

    int first_line = comment_lines("../data/VHriferimento.txt");
    ifstream file("../data/VHriferimento.txt");
    for (int i = 0; i < first_line; i++)
        file.ignore(10000, '\n');

    vector<float> VH0, i0;
    float entry1, entry2, entry3, entry4, entry5, entry6, entry7, entry8;
    string names;
    getline(file, names); // store the names of the variables

    if (count_column("../data/VHriferimento.txt") == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH0.push_back(entry1);
            i0.push_back(entry2);
        }
    }
    else if (count_column("../data/VHriferimento.txt") == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH0.push_back(entry1);
            i0.push_back(entry2);
        }
    }

    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////

    vector<float> sVH0, si0;
    for (int i = 0; i < VH0.size(); i++)
    {
        entry1 = abs(VH0.at(i)) * 0.02;
        entry2 = abs(i0.at(i)) * 0.006 + 0.02;

        sVH0.push_back(entry1);
        si0.push_back(entry2);
    }

    string str1("\tsVH0[mV]"), str2("\tsi0[mA]");
    if (names.find(str1) == string::npos)
    {
        names += "\tsVH0[mV]";
        append_column("../data/VHriferimento.txt", "sVH0[mV]", sVH0);
    }
    if (names.find(str2) == string::npos)
    {
        names += "\tsi0[mA]";
        append_column("../data/VHriferimento.txt", "si0[mA]", si0);
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////

    file.close();

    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////

    cout << endl
         << "Dati della misura di riferimento:" << endl
         << names << endl;
    for (int i = 0; i < VH0.size(); i++)
        cout << VH0.at(i) << "\t" << i0.at(i) << "\t" << sVH0.at(i) << "\t" << si0.at(i) << endl;
    cout << endl;

    /////////////////////////////// FIT //////////////////////////////////////////////////////

    TCanvas *canvas1 = new TCanvas("canvas1", "riferimento", 500, 5, 500, 500);
    canvas1->SetGrid();

    TF1 *tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
    tf1->SetParNames("#chi", "#omega");
    tf1->SetLineColor(38);

    TGraphErrors *graph1 = new TGraphErrors(i0.size(), &i0[0], &VH0[0], &si0[0], &sVH0[0]);
    graph1->SetTitle("#splitline{Misure di riferimento}{V = #chi + #omega i};i [mA];V_{H} [mV]");
    std_graph_settings(*graph1);

    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/tensione_hall0.pdf");

    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << " (Probability: " << tf1->GetProb() << ")." << endl;

    const float chi = tf1->GetParameter(0);
    const float schi = tf1->GetParError(0);

    const float omega = tf1->GetParameter(1);
    const float somega = tf1->GetParError(1);

    cout << endl
         << endl
         << "Parametri della retta di calibrazione del campo di Hall:" << endl
         << "chi (intercetta): (" << chi << " ± " << schi << ") mV" << endl
         << "omega (coeff angolare): (" << omega << " ± " << somega << ") Ω" << endl;   // mV / mA















































    ///////////////////// PART 2: B = const //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////

    first_line = comment_lines("../data/VHBconst_corretto.txt");
    file.open("../data/VHBconst_corretto.txt");
    for (int i = 0; i < first_line; i++)
        file.ignore(10000, '\n');

    vector<float> VH1, i1;
    getline(file, names); // store the names of the variables

    if (count_column("../data/VHBconst_corretto.txt") == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH1.push_back(entry1);
            i1.push_back(entry2);
        }
    }
    else if (count_column("../data/VHBconst_corretto.txt") == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH1.push_back(entry1);
            i1.push_back(entry2);
        }
    }

    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////

    vector<float> VH1_correct, sVH1, si1, sVH1_correct;
    for (int i = 0; i < VH1.size(); i++)
    {
        entry1 = VH1.at(i) - (chi + omega * i1.at(i));      // mV
        VH1_correct.push_back(entry1);

        entry2 = VH1.at(i) * 0.02;
        sVH1.push_back(entry2);

        entry3 = i1.at(i) * 0.006 + 0.03;                   // mA
        si1.push_back(entry3);

        entry4 = sqrt(sVH1.at(i) * sVH1.at(i) + schi * schi + i1.at(i) * i1.at(i) * somega * somega + omega * omega * si1.at(i) * si1.at(i));
        sVH1_correct.push_back(entry4);
    }

    str1 = "\tVH1_correct[mV]";
    str2 = "\tsVH1[mV]";
    string str3("\tsi1[mA]"), str4("\tsVH1_correct[mV]");
    if (names.find(str1) == string::npos)
    {
        names += "\tVH1_correct[mV]";
        append_column("../data/VHBconst_corretto.txt", "\tVH1_correct[mV]", VH1_correct);
    }
    if (names.find(str2) == string::npos)
    {
        names += "\tsVH1[mV]";
        append_column("../data/VHBconst_corretto.txt", "\tsVH1[mV]", sVH1);
    }
    if (names.find(str3) == string::npos)
    {
        names += "\tsi1[mA]";
        append_column("../data/VHBconst_corretto.txt", "\tsi1[mA]", si1);
    }
    if (names.find(str4) == string::npos)
    {
        names += "\tsVH1_correct[mV]";
        append_column("../data/VHBconst_corretto.txt", "\tsVH1_correct[mV]", sVH1_correct);
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////

    file.close();

    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////

    cout << endl
         << "Dati della misura con B = const:" << endl
         << names << endl;
    for (int i = 0; i < VH1.size(); i++)
        cout << VH1.at(i) << "\t" << i1.at(i) << "\t" << VH1_correct.at(i) << "\t" << sVH1.at(i) << "\t"
             << si1.at(i) << "\t" << sVH1_correct.at(i) << endl;
    cout << endl;

    //////////////////////////////// FIT FOR ALL SUB-DATASETS //////////////////////////////////////

    TCanvas *canvas2 = new TCanvas("canvas2", "costante di Hall", 500, 5, 500, 900);
    canvas2->Divide(2, 5);

    const int n1[] = {0, 9, 18, 27, 36, 45, 54, 63, 72, 81, 90}; // position where each subset begins
    const int sets1 = 10;                                        // number of sub-datasets

    vector<float> a1, sa1, b1, sb1;

    for (int i = 0; i < sets1; i++)
    {
        vector<float> sub_VH1correct(VH1_correct.begin() + n1[i], VH1_correct.begin() + n1[i + 1]);
        vector<float> sub_sVH1correct(sVH1_correct.begin() + n1[i], sVH1_correct.begin() + n1[i + 1]);
        vector<float> sub_i1(i1.begin() + n1[i], i1.begin() + n1[i + 1]);
        vector<float> sub_si1(si1.begin() + n1[i], si1.begin() + n1[i + 1]);

        TF1 *tf = new TF1("tf", "[0]+[1]*x", -15, 15);
        tf->SetLineColor(38);
        tf->SetParNames("V_{0}", "R_{H} B / t");
        tf->SetLineWidth(2);

        string title;
        if(i==0) title = "B = 62 mT;i [mA];V_{H} [mV]";
        if(i==1) title = "B = 123 mT;i [mA];V_{H} [mV]";
        if(i==2) title = "B = 185 mT;i [mA];V_{H} [mV]";
        if(i==3) title = "B = 247 mT;i [mA];V_{H} [mV]";
        if(i==4) title = "B = 309 mT;i [mA];V_{H} [mV]";
        if(i==5) title = "B = -62 mT;i [mA];V_{H} [mV]";
        if(i==6) title = "B = -123 mT;i [mA];V_{H} [mV]";
        if(i==7) title = "B = -185 mT;i [mA];V_{H} [mV]";
        if(i==8) title = "B = -247 mT;i [mA];V_{H} [mV]";
        if(i==9) title = "B = -309 mT;i [mA];V_{H} [mV]";

        TGraphErrors *graph = new TGraphErrors(sub_i1.size(), &sub_i1[0], &sub_VH1correct[0], &sub_si1[0], &sub_sVH1correct[0]);
        graph->SetTitle(title.c_str());

        canvas2->cd(i + 1);
        std_graph_settings_adv(*graph);
        gPad->SetGrid();
        graph->Fit(tf, "ER");
        graph->Draw("ap");

        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << " (Probability: " << tf->GetProb() << ")." << endl
             << endl;

        a1.push_back(tf->GetParameter(0));      // mV
        sa1.push_back(tf->GetParError(0));
        b1.push_back(tf->GetParameter(1));      // mV / mA
        sb1.push_back(tf->GetParError(1));
    }
    canvas2->SaveAs("../graphs/tensione_hall1.pdf");

    //////////////// CHECK

    cout << endl
         << "CHECK " << endl;
    for (int i = 0; i < a1.size(); i++)
    {
        cout << a1.at(i) << "\t" << sa1.at(i) << "\t" << b1.at(i) << "\t" << sb1.at(i) << endl;
    }

    /////////////////////////////////// HALL CONSTANT //////////////////////////////////////////////

    const float t = 0.1;                                                            // thickness of sample -- cm
    const float I[] = {0.3, 0.6, 0.9, 1.2, 1.5, -0.3, -0.6, -0.9, -1.2, -1.5};      // A
    float sI[sets1], B[sets1], sB[sets1];
    const float a_tilde = 0.275176;                                                 // mT
    const float sa_tilde = 6.18651;
    const float b_tilde = 205.948;                                                  // mT/A
    const float sb_tilde = 0.281825;
    const float sB_geom = 1.67;                                                     // mT

    vector<float> R_H, sR_H, inv_sRH_squared, RH_over_sRHsquared;

    for (int i = 0; i < sets1; i++)
    {
        sI[i] = 0.002 * I[i] + 0.003;
        B[i] = a_tilde + b_tilde * I[i];                                            // mT
        sB[i] = sqrt((sa_tilde * sa_tilde + pow(sb_tilde * I[i], 2) + pow(b_tilde * sI[i], 2)) + sB_geom * sB_geom);

        R_H.push_back(b1.at(i) * t / B[i]);                                         // V cm / (A mT)
        sR_H.push_back(sqrt(pow(sb1.at(i) * t / B[i], 2) + pow(b1.at(i) * t * sB[i] / (B[i] * B[i]), 2)));

        inv_sRH_squared.push_back(pow(sR_H.at(i), -2)); // will be used for the weighted average
        RH_over_sRHsquared.push_back(R_H.at(i) * inv_sRH_squared.at(i));
    }

    cout << endl << "IMPORTANTE: VALORI DEL CAMPO MAGNETICO" << endl;
    for(int i=0; i<sets1; i++) cout << B[i] << " ± " << sB[i] << endl;
    cout << endl;

    // reference for weighted average: http://web2.ba.infn.it/~palano//statistica/web/lab/chap2/node5_2.html

    const float R_H1 = accumulate(RH_over_sRHsquared.begin(), RH_over_sRHsquared.begin() + 5, 0.) /
                       accumulate(inv_sRH_squared.begin(), inv_sRH_squared.begin() + 5, 0.);
    const float sR_H1 = sqrt(1 / accumulate(inv_sRH_squared.begin(), inv_sRH_squared.begin() + 5, 0.));

    const float R_H2 = accumulate(RH_over_sRHsquared.begin() + 5, RH_over_sRHsquared.end(), 0.) /
                       accumulate(inv_sRH_squared.begin() + 5, inv_sRH_squared.end(), 0.);
    const float sR_H2 = sqrt(1 / accumulate(inv_sRH_squared.begin() + 5, inv_sRH_squared.end(), 0.));

    cout << "Z test: Hall constant R_H [mV cm / (mT mA)], B = const " << endl;
    z_test(R_H1, R_H2, sqrt(sR_H1 * sR_H1 + sR_H2 * sR_H2));

    /////////////////////////// HALL CONSTANT - FIT //////////////////////////////////

    TCanvas *canvas2_ = new TCanvas("canvas2_", "costante di Hall", 500, 5, 500, 500);

    const int iter1[] = {0, 5, 10};

    vector<float> R_H1_fit, sR_H1_fit;

    for (int i = 0; i < 2; i++)
    {
        vector<float> sub_b1, sub_sb1, sub_B, sub_sB;

        if (i == 0)
        {
            for (int j = 0; j < 5; j++)
            {
                sub_b1.push_back(b1.at(j));
                sub_sb1.push_back(sb1.at(j));
                sub_B.push_back(B[j]);
                sub_sB.push_back(sB[j]);
            }
        }
        if (i == 1)
        {
            for (int j = 5; j < 10; j++)
            {
                sub_b1.push_back(b1.at(j));
                sub_sb1.push_back(sb1.at(j));
                sub_B.push_back(B[j]);
                sub_sB.push_back(sB[j]);
            }
        }

        TF1 *tf = new TF1("tf", "[0]+[1]*x", -400, 400);
        tf->SetLineColor(38);
        tf->SetParNames("b_{0}", "R_{H} / t");

        TGraphErrors *graph = new TGraphErrors(sub_B.size(), &sub_B[0], &sub_b1[0], &sub_sB[0], &sub_sb1[0]);
        if(i==0)    graph->SetTitle("#splitline{Costante di Hall}{polarita 1};B [mT];b #left[#frac{mV}{mA}#right]");
        if(i==1)    graph->SetTitle("#splitline{Costante di Hall}{polarita 2};B [mT];b #left[#frac{mV}{mA}#right]");

        std_graph_settings(*graph);
        graph->GetYaxis()->SetTitleOffset(2.3);
        gPad->SetGrid();
        gPad->SetTopMargin(0.13);
        gPad->SetLeftMargin(0.20);
        graph->Fit(tf, "ER");
        graph->Draw("ap");

        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << " (Probability: " << tf->GetProb() << ")." << endl
             << endl;

        R_H1_fit.push_back(tf->GetParameter(1) * t);
        sR_H1_fit.push_back(tf->GetParError(1) * t);

        if(i==0)    canvas2_->SaveAs("../graphs/coeff_Hall1_1.pdf");
        if(i==1)    canvas2_->SaveAs("../graphs/coeff_Hall1_2.pdf");
    }

    cout << endl << "Z Test: Hall constant R_H [mV cm / (mT mA)] - fit results: " << endl;
    z_test(R_H1_fit.at(0), R_H1_fit.at(1), sqrt(pow(sR_H1_fit.at(0), 2) + pow(sR_H1_fit.at(1), 2)));

    const float RH1_mean = (R_H1_fit.at(0) + R_H1_fit.at(1))/2;
    const float sRH1_mean = sqrt(pow(sR_H1_fit.at(0), 2) + pow(sR_H1_fit.at(1), 2));

    cout << "Valore medio coefficiente di Hall (B=const): RH [mV cm / (mT mA)] = (" << RH1_mean << " ± " << sRH1_mean << ")" 
         << endl;











































    ///////////////////// PART 3: i = const //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////

    first_line = comment_lines("../data/VHiconst.txt");
    file.open("../data/VHiconst.txt");
    for (int i = 0; i < first_line; i++)    file.ignore(10000, '\n');

    vector<float> VH2, I2;
    getline(file, names); // store the names of the variables

    if (count_column("../data/VHiconst.txt") == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH2.push_back(entry1);
            I2.push_back(entry2);
        }
    }
    else if (count_column("../data/VHiconst.txt") == 8)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6 >> entry7 >> entry8)
        {
            VH2.push_back(entry1);
            I2.push_back(entry2);
        }
    }

    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////

    const int sets2 = 6; // number of sub-datasets

    const float i2[] = {1.5, 4.5, 7.5, -1.5, -4.5, -7.5}; // mA
    float si2[sets2];
    for (int i = 0; i < sets2; i++)     si2[i] = i2[i] * 0.006 + 0.02;

    const int n2[] = {0, 10, 20, 30, 40, 50, 60}; // position where each subset begins
    vector<float> B2, VH2_correct, sVH2, sI2, sB2, sVH2_correct;

    for (int i = 0; i < VH2.size(); i++)
    {
        int j;
        if (i < n2[1])          j = 0;
        else if (i < n2[2])     j = 1;
        else if (i < n2[3])     j = 2;
        else if (i < n2[4])     j = 3;
        else if (i < n2[5])     j = 4;
        else                    j = 5;

        entry1 = a_tilde + b_tilde * I2.at(i);          // mT
        B2.push_back(entry1);

        entry2 = VH2.at(i) - (chi + omega * i2[j]);     // mV
        VH2_correct.push_back(entry2);

        entry3 = VH2.at(i) * 0.02;
        sVH2.push_back(entry3);

        entry4 = I2.at(i) * 0.002 + 0.003;              // A
        sI2.push_back(entry4);

        entry5 = sqrt((sa_tilde * sa_tilde + pow(sb_tilde * I2.at(i), 2) + pow(b_tilde * sI2.at(i), 2)) + sB_geom * sB_geom);
        sB2.push_back(entry5);

        entry6 = sqrt(sVH2.at(i) * sVH2.at(i) + schi * schi + i2[j] * i2[j] * somega * somega + omega * omega * si2[j] * si2[j]);
        sVH2_correct.push_back(entry6);
    }

    str1 = "\tB2[mT]";
    str2 = "VH2_correct[mV]";
    str3 = "\tsVH2[mV]";
    str4 = "\tsI2[A]";
    string str5("\tsB2[mT]"), str6("\tsVH2_correct[mV]");
    if (names.find(str1) == string::npos)
    {
        names += "\tB2[mT]";
        append_column("../data/VHiconst.txt", "\tB2[mT]", B2);
    }
    if (names.find(str2) == string::npos)
    {
        names += "\tVH2_correct[mV]";
        append_column("../data/VHiconst.txt", "\tVH2_correct[mV]", VH2_correct);
    }
    if (names.find(str3) == string::npos)
    {
        names += "\tsVH2[mV]";
        append_column("../data/VHiconst.txt", "\tsVH2[mV]", sVH2);
    }
    if (names.find(str4) == string::npos)
    {
        names += "\tsI2[A]";
        append_column("../data/VHiconst.txt", "\tsI2[A]", sI2);
    }
    if (names.find(str5) == string::npos)
    {
        names += "\tsB2[mT]";
        append_column("../data/VHiconst.txt", "\tsB2[mT]", sB2);
    }
    if (names.find(str6) == string::npos)
    {
        names += "\tsVH2_correct[mV]";
        append_column("../data/VHiconst.txt", "\tsVH2_correct[mV]", sVH2_correct);
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////

    file.close();

    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////

    cout << endl
         << "Dati della misura con i = const:" << endl
         << names << endl;
    for (int i = 0; i < VH2.size(); i++)
        cout << VH2.at(i) << "\t" << I2.at(i) << "\t" << B2.at(i) << "\t" << VH2_correct.at(i) << "\t"
             << sVH2.at(i) << "\t" << sI2.at(i) << "\t" << sB2.at(i) << "\t" << sVH2_correct.at(i) << endl;
    cout << endl;

    //////////////////////////////// FIT ALL SUB-DATASETS //////////////////////////////////////

    TCanvas *canvas3 = new TCanvas("canvas3", "costante di Hall", 500, 5, 500, 900);
    canvas3->SetGrid();
    canvas3->Divide(2, 3);

    vector<float> a2, sa2, b2, sb2;

    for (int i = 0; i < sets2; i++)
    {
        vector<float> sub_VH2correct(VH2_correct.begin() + n2[i], VH2_correct.begin() + n2[i + 1]);
        vector<float> sub_sVH2correct(sVH2_correct.begin() + n2[i], sVH2_correct.begin() + n2[i + 1]);
        vector<float> sub_B2(B2.begin() + n2[i], B2.begin() + n2[i + 1]);
        vector<float> sub_sB2(sB2.begin() + n2[i], sB2.begin() + n2[i + 1]);

        TF1 *tf = new TF1("tf", "[0]+[1]*x", -400, 400);
        tf->SetLineColor(38);
        tf->SetParNames("V_{0}", "R_{H} B / t");
        tf->SetLineWidth(2);

        string title;
        if(i==0) title = "i = 1.5 mA;B [mT];V_{H} [mV]";
        if(i==1) title = "i = 4.5 mA;B [mT];V_{H} [mV]";
        if(i==2) title = "i = 7.5 mA;B [mT];V_{H} [mV]";
        if(i==3) title = "i = -1.5 mA;B [mT];V_{H} [mV]";
        if(i==4) title = "i = -4.5 mA;B [mT];V_{H} [mV]";
        if(i==5) title = "i = -7.5 mA;B [mT];V_{H} [mV]";


        TGraphErrors *graph = new TGraphErrors(sub_B2.size(), &sub_B2[0], &sub_VH2correct[0], &sub_sB2[0], &sub_sVH2correct[0]);
        graph->SetTitle(title.c_str());

        canvas3->cd(i + 1);
        std_graph_settings(*graph);
        gPad->SetGrid();
        graph->Fit(tf, "ER");
        graph->Draw("ap");

        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << " (Probability: " << tf->GetProb() << ")." << endl
             << endl;

        a2.push_back(tf->GetParameter(0));      // mV
        sa2.push_back(tf->GetParError(0));
        b2.push_back(tf->GetParameter(1));      // mV / mT
        sb2.push_back(tf->GetParError(1));
    }

    canvas3->SaveAs("../graphs/tensione_hall2.pdf");

    /////////////////////////////////// HALL CONSTANT //////////////////////////////////////////////

    vector<float> R_H_, sR_H_, inv_sRH_squared_, RH_over_sRHsquared_;

    for (int i = 0; i < sets2; i++)
    {
        R_H_.push_back(b2.at(i) * t / i2[i]);       // mV cm / (mA mT)
        sR_H_.push_back(sqrt(pow(sb2.at(i) * t / i2[i], 2) + pow(b2.at(i) * t * si2[i] / (i2[i] * i2[i]), 2)));

        inv_sRH_squared_.push_back(pow(sR_H_.at(i), -2)); // will be used for the weighted average
        RH_over_sRHsquared_.push_back(R_H_.at(i) * inv_sRH_squared_.at(i));
    }

    const float R_H1_ = accumulate(RH_over_sRHsquared_.begin(), RH_over_sRHsquared_.begin() + 3, 0.) /
                        accumulate(inv_sRH_squared_.begin(), inv_sRH_squared_.begin() + 3, 0.);
    const float sR_H1_ = sqrt(1 / accumulate(inv_sRH_squared_.begin(), inv_sRH_squared_.begin() + 3, 0.));

    const float R_H2_ = accumulate(RH_over_sRHsquared_.begin() + 3, RH_over_sRHsquared_.end(), 0.) /
                        accumulate(inv_sRH_squared_.begin() + 3, inv_sRH_squared_.end(), 0.);
    const float sR_H2_ = sqrt(1 / accumulate(inv_sRH_squared_.begin() + 3, inv_sRH_squared_.end(), 0.));

    cout << "Z test: Hall constant R_H [mV cm / (mT mA)], i = const " << endl;
    z_test(R_H1_, R_H2_, sqrt(sR_H1_ * sR_H1_ + sR_H2_ * sR_H2_));

    /////////////////////////// HALL CONSTANT - FIT //////////////////////////////////

    TCanvas *canvas3_ = new TCanvas("canvas3_", "costante di Hall", 500, 4, 500, 500);

    const int iter2[] = {0, sets2 / 2, sets2};
    vector<float> R_H2_fit, sR_H2_fit;

    for (int i = 0; i < 2; i++)
    {
        vector<float> sub_b2(b2.begin() + iter2[i], b2.begin() + iter2[i + 1]);
        vector<float> sub_sb2(sb2.begin() + iter2[i], sb2.begin() + iter2[i + 1]);
        vector<float> sub_i2, sub_si2;

        for (int j = iter2[i]; j < iter2[i + 1]; j++)
        {
            sub_i2.push_back(i2[j]);
            sub_si2.push_back(si2[j]);
        }

        TF1 *tf = new TF1("tf", "[0]+[1]*x", -400, 400);
        tf->SetLineColor(38);
        tf->SetParNames("b_{0}", "R_{H} / t");

        TGraphErrors *graph = new TGraphErrors(sub_i2.size(), &sub_i2[0], &sub_b2[0], &sub_si2[0], &sub_sb2[0]);
        if(i==0)    graph->SetTitle("#splitline{Costante di Hall}{polarita 1};i [mA];b #left[#frac{mV}{mT}#right]");
        if(i==1)    graph->SetTitle("#splitline{Costante di Hall}{polarita 2};i [mA];b #left[#frac{mV}{mT}#right]");

        std_graph_settings(*graph);
        graph->GetYaxis()->SetTitleOffset(2.3);
        gPad->SetGrid();
        gPad->SetTopMargin(0.13);
        gPad->SetLeftMargin(0.20);
        graph->Fit(tf, "ER");
        graph->Draw("ap");

        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << " (Probability: " << tf->GetProb() << ")." << endl
             << endl;

        R_H2_fit.push_back(tf->GetParameter(1) * t);
        sR_H2_fit.push_back(tf->GetParError(1) * t);

        if(i==0)    canvas3_->SaveAs("../graphs/coeff_Hall2_1.pdf");
        if(i==1)    canvas3_->SaveAs("../graphs/coeff_Hall2_2.pdf");
    }

    cout << endl << "Z Test: Hall constant R_H [mV cm / (mT mA)] - fit results: " << endl;
    z_test(R_H2_fit.at(0), R_H2_fit.at(1), sqrt(pow(sR_H2_fit.at(0), 2) + pow(sR_H2_fit.at(1), 2)));

    const float RH2_mean = (R_H2_fit.at(0) + R_H2_fit.at(1)) / 2;                 // 10^7 cm^3 / C
    const float sRH2_mean = sqrt(pow(sR_H2_fit.at(0), 2) + pow(sR_H2_fit.at(1), 2));
    cout << endl << "Valore medio del coefficiente di Hall (i=const): R_H  [mV cm / (mT mA)] = ("
         << RH2_mean << " ± " << sRH2_mean  << ") "<< endl;

    cout << "Z Test: costanti di Hall a confronto (metodo B=const e i=const): " << endl;
    z_test(RH1_mean, RH2_mean, sqrt(pow(sRH1_mean, 2) + pow(sRH2_mean, 2)));

    const float rh = (RH1_mean + RH2_mean) / 2;                             // 10^7 cm^3 / C
    const float srh = sqrt(pow(sRH1_mean, 2) + pow(sRH2_mean, 2));
    const float e = 1.60217663e-19;                                         // C
    const float p = 1/(rh * e);                                             // 10^-7 cm^-3
    const float sp = srh / (rh * rh * e);                                   // 10^-7 cm^-3

    cout << endl << "Costante di Hall definitiva: RH = (" << rh << " ± " << srh << ") mV cm / (mT mA)" << endl;
    cout << endl << "Numero dei portatori: p = (" << p << " ± " << sp << ") 10^-7 cm^-3 " << endl; 
}
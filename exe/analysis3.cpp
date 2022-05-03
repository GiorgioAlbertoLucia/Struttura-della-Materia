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
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH0, i0;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables
    
    if(count_column("../data/VHriferimento.txt") == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH0.push_back(entry1);
            i0.push_back(entry2);
        }
    }
    else if(count_column("../data/VHBconst.txt") == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH0.push_back(entry1);
            i0.push_back(entry2);
        }
    }


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    vector<float> sVH0, si0; 
    for(int i = 0; i < VH0.size(); i++)   
    {
        entry1 = VH0.at(i) * 0.02;
        entry2 = i0.at(i) * 0.001 + 0.03;
        sVH0.push_back(entry1);
        si0.push_back(entry2);
    }

    string str1("\tsVH0[mV]"), str2("\tsVH0[mV]");
    if(names.find(str1) == string::npos)
    {
        names += "\tsVH0[mV]";
        append_column("../data/VHBconst.txt", "sVH0[mV]", sVH0);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH1_correct[V]";   
        append_column("../data/VHBconst.txt", "sVH1_correct[V]", sVH1_correct);
    }    



    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati della misura di riferimento:" << endl << names << endl;
    for (int i = 0; i < VH0.size(); i++)   
        cout << VH0.at(i) << "\t" << sVH0.at(i) << "\t" << i0.at(i) << "\t" << si0.at(i) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas1 = new TCanvas("canvas1", "riferimento", 500, 5, 500, 600);
    canvas1->SetGrid();

    TF1 * tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
    tf1->SetLineColor(38);

    TGraphErrors * graph1 = new TGraphErrors(i0.size(), &i0[0], &VH0[0], &si0[0], &sVH0[0]);
    graph1->SetTitle("#splitline{Misure di riferimento}{V = p_{0} + p_{1} i};i [A];V_{H} [V]");
    std_graph_settings(*graph1);
    
    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/tensione_hall0.jpg");
    canvas1->SaveAs("../graphs/tensione_hall0.pdf");
    
    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
    " (Probability: " << tf1->GetProb() << ")." << endl;

    const float chi = tf1->GetParameter(0);
    const float schi = tf1->GetParError(0);

    const float omega = tf1->GetParameter(1);
    const float somega = tf1->GetParError(1);



























    ///////////////////// PART 2: B = const //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    first_line = comment_lines("../data/VHBconst.txt");
    file.open("../data/VHBconst.txt");
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH1, i1;
    getline(file, names);                                            // store the names of the variables

    if(count_column("../data/VHBconst.txt") == 2)
    {
        while (file >> entry1 >> entry2)
        {
            VH1.push_back(entry1);
            i1.push_back(entry2);
        }
    }
    else if(count_column("../data/VHBconst.txt") == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH1.push_back(entry1);
            i1.push_back(entry2);
        }
    }


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    vector<float> VH1_correct, sVH1, si1, sVH1_correct; 
    for(int i = 0; i < VH1.size(); i++)   
    {
        entry1 = VH1.at(i) - (chi + omega*i1.at(i));
        entry2 = VH1.at(i) * 0.02;
        entry3 = i1.at(i) * 0.001 + 0.03;
        entry4 = sqrt(sVH1.at(i)*sVH1.at(i) + schi*schi + i1.at(i)*i1.at(i)*somega*somega + omega*omega*si1.at(i)*si1.at(i));


        VH1_correct.push_back(entry1);
        sVH1.push_back(entry2);
        si1.push_back(entry3);
        sVH1_correct.push_back(entry4);
    }

    string str1("\tVH1_correct[V]"), str2("\tsVH1_correct[V]");
    if(names.find(str1) == string::npos)
    {
        names += "\tVH1_correct[V]";
        append_column("../data/VHBconst.txt", "VH1_correct[V]", VH1_correct);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH1_correct[V]";   
        append_column("../data/VHBconst.txt", "sVH1_correct[V]", sVH1_correct);
    }    

    


    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati della misura con B = const:" << endl << names << endl;
    for (int i = 0; i < VH1.size(); i++)   
        cout << VH1.at(i) << "\t" << sVH1.at(i) << "\t" << i1.at(i) << "\t" << si1.at(i) << "\t" 
        << VH1_correct.at(i) << "\t" << sVH1_correct.at(i) << endl;
    cout << endl;

    //////////////////////////////// FIT FOR ALL SUB-DATASETS //////////////////////////////////////

    TCanvas * canvas2 = new TCanvas("canvas2", "costante di Hall", 500, 5, 500, 600);
    canvas2->SetGrid();
    canvas2->Divide(2,5);

    const float B[] = {1., 2., 3., 4., 5., 1., 2., 3., 4., 5.};
    const float sB[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    const int n1[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};    // position where each subset begins
    const int sets1 = 10;                                           // number of sub-datasets

    vector<float> a1, sa1, b1, sb1;

    for (int i = 0; i < sets1; i++)
    {
        vector<float> sub_VH1correct(VH1_correct.begin()+n1[i], VH1_correct.begin()+n1[i+1]);
        vector<float> sub_sVH1correct(sVH1_correct.begin()+n1[i], sVH1_correct.begin()+n1[i+1]);
        vector<float> sub_i1(i1.begin()+n1[i], i1.begin()+n1[i+1]);
        vector<float> sub_si1(si1.begin()+n1[i], si1.begin()+n1[i+1]);

        TF1 * tf = new TF1("tf", "[0]+[1]*x", -15, 15);
        tf->SetLineColor(38);
        tf->SetParNames("V_{0}", "#frac{R_{H} B}{t}");

        TGraphErrors * graph = new TGraphErrors(sub_i1.size(), &sub_i1[0], &sub_VH1correct[0], &sub_si1[0], &sub_sVH1correct[0]);
        graph->SetTitle("#splitline{Costante di Hall}{V = V_{0} + #frac{R_{H} i B}{t}};i [A];V_{H} [V]");
        std_graph_settings(*graph);
    
        canvas2->cd(i+1);
        graph->Fit(tf, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << 
        " (Probability: " << tf->GetProb() << ")." << endl << endl;

        a1.push_back(tf->GetParameter(0));
        sa1.push_back(tf->GetParError(0));
        b1.push_back(tf->GetParameter(1));
        sb1.push_back(tf->GetParError(1));
    }

    canvas2->SaveAs("../graphs/tensione_hall1.jpg");
    canvas2->SaveAs("../graphs/tensione_hall1.pdf");

    /////////////////////////////////// HALL CONSTANT ////////////////////////////////////////////// 
    
    const float t = 1.;                              // thickness of sample
    const float st = 0.1;
    
    vector<float> R_H, sR_H, inv_sRH_squared, RH_over_sRHsquared;

    for (int i = 0; i < sets1; i++)
    {
        R_H.push_back(b1.at(i) * t / B[i]);
        sR_H.push_back(sqrt(pow(sb1.at(i)*t/B[i], 2) + pow(b1.at(i)*st/B[i], 2) + pow(b1.at(i)*t*sB[i]/(B[i]*B[i]), 2)));
        
        inv_sRH_squared.push_back(pow(sR_H.at(i), -2));  // will be used for the weighted average
        RH_over_sRHsquared.push_back(R_H.at(i)*inv_sRH_squared.at(i));
    }

    // reference for weighted average: http://web2.ba.infn.it/~palano//statistica/web/lab/chap2/node5_2.html 

    const float R_H1 = accumulate(RH_over_sRHsquared.begin(), RH_over_sRHsquared.begin()+5, 0.)/
                        accumulate(inv_sRH_squared.begin(), inv_sRH_squared.begin()+5, 0.);
    const float sR_H1 = sqrt(1/accumulate(inv_sRH_squared.begin(), inv_sRH_squared.begin()+5, 0.));

    const float R_H2 = accumulate(RH_over_sRHsquared.begin()+5, RH_over_sRHsquared.end(), 0.)/
                        accumulate(inv_sRH_squared.begin()+5, inv_sRH_squared.end(), 0.);
    const float sR_H2 = sqrt(1/accumulate(inv_sRH_squared.begin()+5, inv_sRH_squared.end(), 0.));

    cout << "Z test: Hall constant R_H, B = const " << endl;
    z_test(R_H1, R_H2, sqrt(sR_H1*sR_H1 + sR_H2*sR_H2));

    
    /////////////////////////// HALL CONSTANT - FIT //////////////////////////////////

    TCanvas * canvas2_ = new TCanvas("canvas2_", "costante di Hall", 500, 5, 500, 600);
    canvas2_->SetGrid();
    canvas2_->Divide(1,2);

    const int iter1[] = {0, sets1/2, sets1};

    vector<float> R_H1_fit, sR_H1_fit;

    for (int i = 0; i < 2; i++)
    {
        vector<float> sub_b1(b1.begin()+iter1[i], b1.begin()+iter1[i+1]);
        vector<float> sub_sb1(sb1.begin()+iter1[i], sb1.begin()+iter1[i+1]);
        vector<float> sub_B, sub_sB;
        
        for(int j = iter1[i]; j < iter1[i+1]; j++)
        {
            sub_B.push_back(B[j]);
            sub_sB.push_back(sB[j]);
        }

        TF1 * tf = new TF1("tf", "[0]+[1]*x", -15, 15);
        tf->SetLineColor(38);
        tf->SetParNames("b_{0}", "#frac{R_{H}}{t}");

        TGraphErrors * graph = new TGraphErrors(sub_B.size(), &sub_B[0], &sub_b1[0], &sub_sB[0], &sub_sb1[0]);
        graph->SetTitle("#splitline{Costante di Hall}{b = b_{0} + #frac{R_{H} B}{t}};B [T];b []");
        std_graph_settings(*graph);
    
        canvas2_->cd(i+1);
        graph->Fit(tf, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << 
        " (Probability: " << tf->GetProb() << ")." << endl << endl;

        R_H1_fit.push_back(tf->GetParameter(1) * t);
        sR_H1_fit.push_back(sqrt(pow(tf->GetParError(1)*t, 2) + pow(tf->GetParameter(1)*st, 2)));
    }

    canvas2_->SaveAs("../graphs/coeff_Hall1.jpg");
    canvas2_->SaveAs("../graphs/coeff_Hall1.pdf");

    cout << "Z Test: Hall constant - fit results: " << endl;
    z_test(R_H1_fit.at(0), R_H1_fit.at(1), sqrt(pow(sR_H1_fit.at(0), 2) + pow(sR_H1_fit.at(1), 2)));

























    ///////////////////// PART 3: i = const //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    first_line = comment_lines("../data/VHiconst.txt");
    file.open("../data/VHiconst.txt");
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH2, sVH2, B2, sB2;
    getline(file, names);                                            // store the names of the variables

    if(count_column("../data/VHiconst.txt") == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            B2.push_back(entry3);
            sB2.push_back(entry4);
        }
    }
    else if(count_column("../data/VHiconst.txt") == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            B2.push_back(entry3);
            sB2.push_back(entry4);
        }
    }


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    const float i2[] = {1., 2., 3., 1., 2., 3.};
    const float si2[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    const int n2[] = {0, 5, 10, 15, 20, 25, 30};            // position where each subset begins

    vector<float> VH2_correct, sVH2_correct; 
    for(int i = 0; i < VH2.size(); i++)   
    {
        int j;
        if (i < n2[1])          j = 0;
        else if (i < n2[2])     j = 1;
        else if (i < n2[3])     j = 2;
        else if (i < n2[4])     j = 3;
        else if (i < n2[5])     j = 4;
        else                    j = 5;

        entry1 = VH2.at(i) - (chi + omega*i2[j]);
        entry2 = sqrt(sVH2.at(i)*sVH2.at(i) + schi*schi + i2[j]*i2[j]*somega*somega + omega*omega*si2[j]*si2[j]);

        VH2_correct.push_back(entry1);
        sVH2_correct.push_back(entry2);
    }

    str1 = "\tVH2_correct[V]";
    str2 = "\tsVH2_correct[V]";
    if(names.find(str1) == string::npos)
    {
        names += "\tVH2_correct[V]";
        append_column("../data/VHiconst.txt", "VH2_correct[V]", VH2_correct);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH2_correct[V]";   
        append_column("../data/VHiconst.txt", "sVH2_correct[V]", sVH2_correct);
    }    

    


    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati della misura con B = const:" << endl << names << endl;
    for (int i = 0; i < VH2.size(); i++)   
        cout << VH2.at(i) << "\t" << sVH2.at(i) << "\t" << B2.at(i) << "\t" << sB2.at(i) << "\t" 
        << VH2_correct.at(i) << "\t" << sVH2_correct.at(i) << endl;
    cout << endl;

    //////////////////////////////// FIT ALL SUB-DATASETS //////////////////////////////////////

    TCanvas * canvas3 = new TCanvas("canvas3", "costante di Hall", 500, 5, 500, 600);
    canvas3->SetGrid();
    canvas3->Divide(2,3);

    const int sets2 = 6;                                            // number of sub-datasets

    vector<float> a2, sa2, b2, sb2;

    for (int i = 0; i < sets2; i++)
    {
        vector<float> sub_VH2correct(VH2_correct.begin()+n2[i], VH2_correct.begin()+n2[i+1]);
        vector<float> sub_sVH2correct(sVH2_correct.begin()+n2[i], sVH2_correct.begin()+n2[i+1]);
        vector<float> sub_B2(B2.begin()+n2[i], B2.begin()+n2[i+1]);
        vector<float> sub_sB2(sB2.begin()+n2[i], sB2.begin()+n2[i+1]);

        TF1 * tf = new TF1("tf", "[0]+[1]*x", -15, 15);
        tf->SetLineColor(38);
        tf->SetParNames("V_{0}", "#frac{R_{H} B}{t}");

        TGraphErrors * graph = new TGraphErrors(sub_B2.size(), &sub_B2[0], &sub_VH2correct[0], &sub_sB2[0], &sub_sVH2correct[0]);
        graph->SetTitle("#splitline{Costante di Hall}{V = V_{0} + #frac{R_{H} i B}{t}};B [T];V_{H} [V]");
        std_graph_settings(*graph);
    
        canvas3->cd(i+1);
        graph->Fit(tf, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << 
        " (Probability: " << tf->GetProb() << ")." << endl << endl;

        a2.push_back(tf->GetParameter(0));
        sa2.push_back(tf->GetParError(0));
        b2.push_back(tf->GetParameter(1));
        sb2.push_back(tf->GetParError(1));
    }

    canvas3->SaveAs("../graphs/tensione_hall2.jpg");
    canvas3->SaveAs("../graphs/tensione_hall2.pdf");

    /////////////////////////////////// HALL CONSTANT ////////////////////////////////////////////// 
   
    vector<float> R_H_, sR_H_, inv_sRH_squared_, RH_over_sRHsquared_;

    for (int i = 0; i < sets2; i++)
    {
        R_H_.push_back(b2.at(i) * t / i2[i]);
        sR_H_.push_back(sqrt(pow(sb2.at(i)*t/i2[i], 2) + pow(b2.at(i)*st/i2[i], 2) + pow(b2.at(i)*t*si2[i]/(i2[i]*i2[i]), 2)));
        
        inv_sRH_squared_.push_back(pow(sR_H_.at(i), -2));  // will be used for the weighted average
        RH_over_sRHsquared_.push_back(R_H_.at(i)*inv_sRH_squared_.at(i));
    }

    const float R_H1_ = accumulate(RH_over_sRHsquared_.begin(), RH_over_sRHsquared_.begin()+3, 0.)/
                        accumulate(inv_sRH_squared_.begin(), inv_sRH_squared_.begin()+3, 0.);
    const float sR_H1_ = sqrt(1/accumulate(inv_sRH_squared_.begin(), inv_sRH_squared_.begin()+3, 0.));

    const float R_H2_ = accumulate(RH_over_sRHsquared_.begin()+3, RH_over_sRHsquared_.end(), 0.)/
                        accumulate(inv_sRH_squared_.begin()+3, inv_sRH_squared_.end(), 0.);
    const float sR_H2_ = sqrt(1/accumulate(inv_sRH_squared_.begin()+3, inv_sRH_squared_.end(), 0.));

    cout << "Z test: Hall constant R_H, i = const " << endl;
    z_test(R_H1_, R_H2_, sqrt(sR_H1_*sR_H1_ + sR_H2_*sR_H2_));
    
    /////////////////////////// HALL CONSTANT - FIT //////////////////////////////////

    TCanvas * canvas3_ = new TCanvas("canvas3_", "costante di Hall", 500, 5, 500, 600);
    canvas3_->SetGrid();
    canvas3_->Divide(1,2);

    const int iter2[] = {0, sets2/2, sets2};
    vector<float> R_H2_fit, sR_H2_fit;

    for (int i = 0; i < 2; i++)
    {
        vector<float> sub_b2(b2.begin()+iter2[i], b2.begin()+iter2[i+1]);
        vector<float> sub_sb2(sb2.begin()+iter2[i], sb2.begin()+iter2[i+1]);
        vector<float> sub_i2, sub_si2;
        
        for(int j = iter2[i]; j < iter2[i+1]; j++)
        {
            sub_i2.push_back(B[j]);
            sub_si2.push_back(sB[j]);
        }

        TF1 * tf = new TF1("tf", "[0]+[1]*x", -15, 15);
        tf->SetLineColor(38);
        tf->SetParNames("b_{0}", "#frac{R_{H}}{t}");

        TGraphErrors * graph = new TGraphErrors(sub_i2.size(), &sub_i2[0], &sub_b2[0], &sub_si2[0], &sub_sb2[0]);
        graph->SetTitle("#splitline{Costante di Hall}{b = b_{0} + #frac{R_{H} i}{t}};i [A];b []");
        std_graph_settings(*graph);
    
        canvas3_->cd(i+1);
        graph->Fit(tf, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf->GetChisquare() << ", number of DoF: " << tf->GetNDF() << 
        " (Probability: " << tf->GetProb() << ")." << endl << endl;

        R_H2_fit.push_back(tf->GetParameter(1) * t);
        sR_H2_fit.push_back(sqrt(pow(tf->GetParError(1)*t, 2) + pow(tf->GetParameter(1)*st, 2)));
    }

    canvas3_->SaveAs("../graphs/coeff_Hall2.jpg");
    canvas3_->SaveAs("../graphs/coeff_Hall2.pdf");

    cout << "Z Test: Hall constant - fit results: " << endl;
    z_test(R_H2_fit.at(0), R_H2_fit.at(1), sqrt(pow(sR_H2_fit.at(0), 2) + pow(sR_H2_fit.at(1), 2)));

    

}
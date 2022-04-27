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
    
    ///////////////////// PART 1: REFERENCE MEASUREMENTS //////////////////////////////////////

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    int first_line = comment_lines("../data/VHriferimento.txt");
    ifstream file("../data/mappatura.txt");
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> VH0, sVH0, i0, si0;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables
    
    while (file >> entry1 >> entry2 >> entry3 >> entry4)
    {
        VH0.push_back(entry1);
        sVH0.push_back(entry2);
        i0.push_back(entry3);
        si0.push_back(entry4);
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

    vector<float> VH1, sVH1, i1, si1;
    getline(file, names);                                            // store the names of the variables

    if(count_column("../data/VHBconst.txt") == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH1.push_back(entry1);
            sVH1.push_back(entry2);
            i1.push_back(entry3);
            si1.push_back(entry4);
        }
    }
    else if(count_column("../data/VHBconst.txt") == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH1.push_back(entry1);
            sVH1.push_back(entry2);
            i1.push_back(entry3);
            si1.push_back(entry4);
        }
    }

    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    vector<float> VH1_correct, sVH1_correct; 
    for(int i = 0; i < VH1.size(); i++)   
    {
        entry1 = VH1.at(i) - (chi + omega*i1.at(i));
        entry2 = sqrt(sVH1.at(i)*sVH1.at(i) + schi*schi + i1.at(i)*i1.at(i)*somega*somega + omega*omega*si1.at(i)*si1.at(i));

        VH1_correct.push_back(entry1);
        sVH1_correct.push_back(entry2);
    }

    string str1(" VH1_correct[V] "), str2(" sVH1_correct[V] ");
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

    //////////////////////////////// DIVIDE TWO DATASETS //////////////////////////////////////

    const int n1 = 5;
    vector<float> sub1_VH1correct(VH1_correct.begin(), VH1_correct.begin()+n1-1);
    vector<float> sub1_sVH1correct(sVH1_correct.begin(), sVH1_correct.begin()+n1-1);
    vector<float> sub1_i1(i1.begin(), i1.begin()+n1-1);
    vector<float> sub1_si1(si1.begin(), si1.begin()+n1-1);

    vector<float> sub2_VH1correct(VH1_correct.begin()+n1, VH1_correct.end());
    vector<float> sub2_sVH1correct(sVH1_correct.begin()+n1, sVH1_correct.end());
    vector<float> sub2_i1(i1.begin()+n1, i1.end());
    vector<float> sub2_si1(si1.begin()+n1, si1.end());

    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas2 = new TCanvas("canvas2", "costante di Hall", 500, 5, 500, 600);
    canvas1->SetGrid();

    TF1 * tf2 = new TF1("tf2", "[0]+[1]*x", -15, 15);
    tf2->SetLineColor(38);
    tf2->SetParNames("V_{0}", "#frac{R_{H} B}{t}");

    TGraphErrors * graph2 = new TGraphErrors(sub1_i1.size(), &sub1_i1[0], &sub1_VH1correct[0], &sub1_si1[0], &sub1_sVH1correct[0]);
    graph2->SetTitle("#splitline{Determinazione della}{costante di Hall}{V = V_{0} + #frac{R_{H} i B}{t}};i [A];V_{H} [V]");
    std_graph_settings(*graph2);
    
    graph2->Fit(tf2, "ER");
    graph2->Draw("ap");
    canvas2->SaveAs("../graphs/tensione_hall1.jpg");
    canvas2->SaveAs("../graphs/tensione_hall1.pdf");
    
    cout << "Chi^2:" << tf2->GetChisquare() << ", number of DoF: " << tf2->GetNDF() << 
    " (Probability: " << tf2->GetProb() << ")." << endl;

    const float a1 = tf2->GetParameter(0);
    const float sa1 = tf2->GetParError(0);

    const float b1 = tf2->GetParameter(1);
    const float sb1 = tf2->GetParError(1);

    //-------------------------------------//

    TCanvas * canvas3 = new TCanvas("canvas3", "costante di Hall", 500, 5, 500, 600);
    canvas1->SetGrid();

    TF1 * tf3 = new TF1("tf3", "[0]+[1]*x", -15, 15);
    tf3->SetLineColor(38);
    tf3->SetParNames("V_{0}", "#frac{R_{H} B}{t}");

    TGraphErrors * graph3 = new TGraphErrors(sub2_i1.size(), &sub2_i1[0], &sub2_VH1correct[0], &sub2_si1[0], &sub2_sVH1correct[0]);
    graph3->SetTitle("#splitline{Determinazione della}{costante di Hall}{V = V_{0} + #frac{R_{H} i B}{t}};i [A];V_{H} [V]");
    std_graph_settings(*graph3);
    
    graph3->Fit(tf3, "ER");
    graph3->Draw("ap");
    canvas3->SaveAs("../graphs/tensione_hall2.jpg");
    canvas3->SaveAs("../graphs/tensione_hall2.pdf");
    
    cout << "Chi^2:" << tf3->GetChisquare() << ", number of DoF: " << tf3->GetNDF() << 
    " (Probability: " << tf3->GetProb() << ")." << endl;

    const float a2 = tf3->GetParameter(0);
    const float sa2 = tf3->GetParError(0);

    const float b2 = tf3->GetParameter(1);
    const float sb2 = tf3->GetParError(1);

    /////////////////////////////////// TEST Z //////////////////////////////////////////////
    z_test(a1, a2, sqrt(sa1*sa1 + sa2*sa2));
    z_test(b1, b2, sqrt(sb1*sb1 + sb2*sb2));

    const float b = (b1 + b2)/2;                    // average p1
    const float sb = sqrt(sb1*sb1 + sb2*sb2)/2;

    const float B = 1;
    const float sB = 0;

    const float t = 1;                              // thickness of sample
    const float st = 0.1;

    const float R_H = b * t / B;
    
    

}
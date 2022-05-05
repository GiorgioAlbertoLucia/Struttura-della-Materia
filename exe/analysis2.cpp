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

void analysis2()
{   
    ///////////////////// CURVA DI ISTERESI /////////////////////////////////////////////////////


    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    freopen("../output/analysis2.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");


    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path = "../data/isteresi2.txt";
    int first_line = comment_lines(path);
    ifstream file(path);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n'); 

    const int n_sets[] = {21, 21, 21, 21};       // how many points for each dataset
    int n = 0;

    vector<float> I1, I2, I3, I4, sI1, sI2, sI3, sI4, B1, B2, B3, B4, sB1, sB2, sB3, sB4;
    float entry1, entry2, entry3, entry4;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path) == 3)
    {
        while (file >> entry1 >> entry2 >> entry3)
        {
            if(n<n_sets[0])
            {
                B1.push_back(entry1);
                sB1.push_back(entry2);
                I1.push_back(entry3);
            }
            if(n>=n_sets[0] && n<(n_sets[0]+n_sets[1]))
            {
                B2.push_back(entry1);
                sB2.push_back(entry2);
                I2.push_back(entry3);
            }
            if(n>=(n_sets[0]+n_sets[1]) && n<(n_sets[0]+n_sets[1]+n_sets[2]))
            {
                B3.push_back(entry1);
                sB3.push_back(entry2);
                I3.push_back(entry3);
            }
            if(n>=(n_sets[0]+n_sets[1]+n_sets[2]) && n<(n_sets[0]+n_sets[1]+n_sets[2]+n_sets[3]))
            {
                B4.push_back(entry1);
                sB4.push_back(entry2);
                I4.push_back(entry3);
            }
            n++;
        }
    }
    else if(count_column(path) == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            if(n<n_sets[0])
            {
                B1.push_back(entry1);
                sB1.push_back(entry2);
                I1.push_back(entry3);
                sI1.push_back(entry4);
            }
            if(n>=n_sets[0] && n<(n_sets[0]+n_sets[1]))
            {
                B2.push_back(entry1);
                sB2.push_back(entry2);
                I2.push_back(entry3);
                sI2.push_back(entry4);
            }
            if(n>=(n_sets[0]+n_sets[1]) && n<(n_sets[0]+n_sets[1]+n_sets[2]))
            {
                B3.push_back(entry1);
                sB3.push_back(entry2);
                I3.push_back(entry3);
                sI3.push_back(entry4);
            }
            if(n>=(n_sets[0]+n_sets[1]+n_sets[2]) && n<(n_sets[0]+n_sets[1]+n_sets[2]+n_sets[3]))
            {
                B4.push_back(entry1);
                sB4.push_back(entry2);
                I4.push_back(entry3);
                sI4.push_back(entry4);
            }
            n++;
        }
    }

    
    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    for(int j = 0; j < 4; j++)
    {
        if(j == 0)  for(int i = 0; i < I1.size(); i++)   
                        {
                            entry1 = abs(I1.at(i)) * 0.002 + 0.003;
                            sI1.push_back(entry1);
                        }
        if(j == 1)  for(int i = 0; i < I2.size(); i++)   
                        {
                            entry1 = abs(I2.at(i)) * 0.002 + 0.003;
                            sI2.push_back(entry1);
                        }
        if(j == 2)  for(int i = 0; i < I3.size(); i++)   
                        {
                            entry1 = abs(I3.at(i)) * 0.002 + 0.003;
                            sI3.push_back(entry1);
                        }
        if(j == 3)  for(int i = 0; i < I4.size(); i++)   
                        {
                            entry1 = abs(I4.at(i)) * 0.002 + 0.003;
                            sI4.push_back(entry1);
                        }
    }
    
    vector<float> sI;
    for(int j = 0; j < 4; j++)
    {
        if(j == 0)  for(int i = 0; i < sI1.size(); i++)  sI.push_back(sI1.at(i));
        if(j == 1)  for(int i = 0; i < sI2.size(); i++)  sI.push_back(sI2.at(i)); 
        if(j == 2)  for(int i = 0; i < sI3.size(); i++)  sI.push_back(sI3.at(i)); 
        if(j == 3)  for(int i = 0; i < sI4.size(); i++)  sI.push_back(sI4.at(i));
    }
    string str1("\tsI[mA]");
    if(names.find(str1) == string::npos)
    {
        names += "\tsI[mA]";
        append_column(path, "sI[mA]", sI);
    }   


    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati prima discesa:" << endl << names << endl;
    for (int i = 0; i < I1.size(); i++)   
        cout << B1.at(i) << "\t" << sB1.at(i) << "\t" << I1.at(i) << "\t" << sI1.at(i) << endl;
    cout << endl << "Dati prima salita:" << endl << names << endl;
    for (int i = 0; i < I2.size(); i++)   
        cout << B2.at(i) << "\t" << sB2.at(i) << "\t" << I2.at(i) << "\t" << sI2.at(i) << endl;
    cout << endl << "Dati seconda discesa:" << endl << names << endl;
    for (int i = 0; i < I3.size(); i++)   
        cout << B3.at(i) << "\t" << sB3.at(i) << "\t" << I3.at(i) << "\t" << sI3.at(i) << endl;
    cout << endl << "Dati seconda salita:" << endl << names << endl;
    for (int i = 0; i < I4.size(); i++)   
        cout << B4.at(i) << "\t" << sB4.at(i) << "\t" << I4.at(i) << "\t" << sI4.at(i) << endl;
    cout << endl;
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas1 = new TCanvas("canvas1", "isteresi", 500, 5, 500, 600);
    canvas1->SetGrid();

    TF1 * tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
    tf1->SetLineColor(38);

    TGraphErrors * graph1 = new TGraphErrors(I1.size(), &I1[0], &B1[0], &sI1[0], &sB1[0]);
    graph1->SetTitle("#splitline{Ciclo di isteresi}{prima discesa};I [A];B [T]");
    std_graph_settings(*graph1);
    
    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/isteresi_1.jpg");
    canvas1->SaveAs("../graphs/isteresi_1.pdf");
    
    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
    " (Probability: " << tf1->GetProb() << ")." << endl;

    TCanvas * canvas2 = new TCanvas("canvas2", "isteresi", 500, 5, 500, 600);
    canvas2->SetGrid();

    TF1 * tf2 = new TF1("tf2", "[0]+[1]*x", -15, 15);
    tf2->SetLineColor(38);

    TGraphErrors * graph2 = new TGraphErrors(I2.size(), &I2[0], &B2[0], &sI2[0], &sB2[0]);
    graph2->SetTitle("#splitline{Ciclo di isteresi}{prima salita};I [A];B [T]");
    std_graph_settings(*graph2);
    
    graph2->Fit(tf2, "ER");
    graph2->Draw("ap");
    canvas2->SaveAs("../graphs/isteresi_2.jpg");
    canvas2->SaveAs("../graphs/isteresi_2.pdf");
    
    cout << "Chi^2:" << tf2->GetChisquare() << ", number of DoF: " << tf2->GetNDF() << 
    " (Probability: " << tf2->GetProb() << ")." << endl;

    TCanvas * canvas3 = new TCanvas("canvas3", "isteresi", 500, 5, 500, 600);
    canvas3->SetGrid();

    TF1 * tf3 = new TF1("tf3", "[0]+[1]*x", -15, 15);
    tf3->SetLineColor(38);

    TGraphErrors * graph3 = new TGraphErrors(I3.size(), &I3[0], &B3[0], &sI3[0], &sB3[0]);
    graph3->SetTitle("#splitline{Ciclo di isteresi}{seconda discesa};I [A];B [T]");
    std_graph_settings(*graph3);
    
    graph3->Fit(tf3, "ER");
    graph3->Draw("ap");
    canvas3->SaveAs("../graphs/isteresi_3.jpg");
    canvas3->SaveAs("../graphs/isteresi_3.pdf");
    
    cout << "Chi^2:" << tf3->GetChisquare() << ", number of DoF: " << tf3->GetNDF() << 
    " (Probability: " << tf3->GetProb() << ")." << endl;

    TCanvas * canvas4 = new TCanvas("canvas4", "isteresi", 500, 5, 500, 600);
    canvas4->SetGrid();

    TF1 * tf4 = new TF1("tf4", "[0]+[1]*x", -15, 15);
    tf4->SetLineColor(38);

    TGraphErrors * graph4 = new TGraphErrors(I4.size(), &I4[0], &B4[0], &sI4[0], &sB4[0]);
    graph4->SetTitle("#splitline{Ciclo di isteresi}{seconda discesa};I [A];B [T]");
    std_graph_settings(*graph4);
    
    graph4->Fit(tf4, "ER");
    graph4->Draw("ap");
    canvas4->SaveAs("../graphs/isteresi_4.jpg");
    canvas4->SaveAs("../graphs/isteresi_4.pdf");
    
    cout << "Chi^2:" << tf4->GetChisquare() << ", number of DoF: " << tf4->GetNDF() << 
    " (Probability: " << tf4->GetProb() << ")." << endl;


    ////////////////////////////////// ∆B /////////////////////////////////////////////

    float a1 = tf1->GetParameter(0);
    float a2 = tf2->GetParameter(0);
    float a3 = tf3->GetParameter(0);
    float a4 = tf4->GetParameter(0);

    float sa1 = tf1->GetParameter(0);
    float sa2 = tf2->GetParameter(0);
    float sa3 = tf3->GetParameter(0);
    float sa4 = tf4->GetParameter(0);

    float b1 = tf1->GetParameter(1);
    float b2 = tf2->GetParameter(1);
    float b3 = tf3->GetParameter(1);
    float b4 = tf4->GetParameter(1);

    float sb1 = tf1->GetParameter(1);
    float sb2 = tf2->GetParameter(1);
    float sb3 = tf3->GetParameter(1);
    float sb4 = tf4->GetParameter(1);

    // check compatibilità
    z_test(a1, a3, sqrt(sa1*sa1 + sa3*sa3));
    z_test(a2, a4, sqrt(sa2*sa2 + sa4*sa4));
    z_test(b1, b3, sqrt(sb1*sb1 + sb3*sb3));
    z_test(b2, b4, sqrt(sb2*sb2 + sb4*sb4));

    float c1 = (a1 + a3)/2;
    float c2 = (a2 + a4)/2;
    float sc1 = sqrt(sa1*sa1/4 + sa3*sa3/4);
    float sc2 = sqrt(sa2*sa2/4 + sa4*sa4/4);

    float d1 = (b1 + b3)/2;
    float d2 = (b2 + b4)/2;
    float sd1 = sqrt(sb1*sb1/4 + sb3*sb3/4);
    float sd2 = sqrt(sb2*sb2/4 + sb4*sb4/4);

    z_test(c1, c2, sqrt(sc1*sc1 + sc2*sc2));
    z_test(d1, d2, sqrt(sd1*sd1 + sd2*sd2));

    float e = (c1 + c2)/2;                  // parametro per calibrazione di B. Neglia appunti è a tilde
    float f = (d1 + d2)/2;                  // parametro per calibrazione di B. Neglia appunti è b tilde
    float se = abs((c1 - c2)/2);
    float sf = sqrt(sd1*sd1/4 + sd2*sd2/4);

    cout << "Calibrazione B:" << endl << "intercetta: q = " << e << " ± " << se << endl << "coeff angolare: m = " <<
            f << " ± " << sf << endl; 
    
    // per calcolare ∆B_cal nei prossimi file di analisi dati, tener conto del fatto che 
    // B = e + f * I -> ∆B = ...

}
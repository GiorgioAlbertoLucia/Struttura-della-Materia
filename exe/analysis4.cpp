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

void analysis4()
{   
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis4.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");


    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const float R_H = 0.;
    const float sR_H = 0.;

    const float t = 0.;     // z axis
    const float st = 0.;
    const float d = 0.;     // y axis
    const float sd = 0.;
    const float l = 0.;     // x axis
    const float sl = 0.;

    const char * path = "../data/mobilitàportatori.txt";
    int first_line = comment_lines(path);
    ifstream file(path);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> V, sV, i, si;
    float entry1, entry2, entry3, entry4;
    string names;
    getline(file, names);                                            // store the names of the variables

    while (file >> entry1 >> entry2 >> entry3 >> entry4)
    {
        V.push_back(entry1);
        sV.push_back(entry2);
        i.push_back(entry3);
        si.push_back(entry4);
    }
    

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int j = 0; j < V.size(); j++)   
        cout << V.at(j) << "\t" << sV.at(j) << "\t" << i.at(j) << "\t" << si.at(j) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas = new TCanvas("canvas", "mobilità", 500, 5, 500, 600);
    canvas->SetGrid();
    canvas->Divide(1,2);
    
    const int n[] = {0, 2, 4};    // position where each subset begins
    const int sets = 2;

    vector<float> V0, sV0, R, sR;

    for (int j = 0; j < sets; j++)
    {
        vector<float> sub_V(V.begin()+n[j], V.begin()+n[j+1]);
        vector<float> sub_sV(sV.begin()+n[j], sV.begin()+n[j+1]);
        vector<float> sub_i(i.begin()+n[j], i.begin()+n[j+1]);
        vector<float> sub_si(si.begin()+n[j], si.begin()+n[j+1]);

        TF1 * tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
        tf1->SetLineColor(38);

        TGraphErrors * graph = new TGraphErrors(sub_i.size(), &sub_i[0], &sub_V[0], &sub_si[0], &sub_sV[0]);
        graph->SetTitle("#splitline{Mobilità dei portatori}{V = V_{0} + R i};i [mA];V [V]");
        std_graph_settings(*graph);
    
        canvas->cd(j+1);
        graph->Fit(tf1, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
        " (Probability: " << tf1->GetProb() << ")." << endl;

        V0.push_back(tf1->GetParameter(0));
        sV0.push_back(tf1->GetParError(0));
        R.push_back(tf1->GetParameter(1));
        sR.push_back(tf1->GetParError(1));
    }
    
    canvas->SaveAs("../graphs/mobilità.jpg");
    canvas->SaveAs("../graphs/mobilità.pdf");

    ////////////////////////////////////// MOBILITÀ /////////////////////////////////////
    cout << "Z test, intercetta V0 compatibile con 0: " << endl;
    z_test(V0.at(0), 0, sV0.at(0));
    z_test(V0.at(1), 0, sV0.at(1));
    
    const float sigma1 = (t * d)/(R.at(0) * l);
    const float ssigma1 = sqrt(pow(st*d/(R.at(0)*l), 2) + pow(t*sd/(R.at(0)*l), 2) + pow(t*d*sl/(R.at(0)*l*l), 2) + 
                            pow(t*d*sR.at(0)/(R.at(0)*R.at(0)*l), 2));

    const float sigma2 = (t * d)/(R.at(1) * l);
    const float ssigma2 = sqrt(pow(st*d/(R.at(1)*l), 2) + pow(t*sd/(R.at(1)*l), 2) + pow(t*d*sl/(R.at(1)*l*l), 2) + 
                            pow(t*d*sR.at(1)/(R.at(1)*R.at(1)*l), 2));

    const float mu1 = R_H * sigma1;
    const float smu1 = sqrt(pow(sR_H*sigma1, 2) + pow(R_H*ssigma1, 2));
    const float mu2 = R_H * sigma2;
    const float smu2 = sqrt(pow(sR_H*sigma2, 2) + pow(R_H*ssigma2, 2));

    z_test(mu1, mu2, sqrt(smu1*smu1 + smu2*smu2));

    ///////////////////////////// DISALLINEAMENTO ///////////////////////////
    const float omega = 1;
    const float somega = 0.1;

    const float sigma = (sigma1 + sigma2)/2;
    const float ssigma = sqrt(ssigma1*ssigma1 + ssigma2*ssigma2)/2;

    const float dx = omega * sigma * d * t;
    const float sdx = sqrt(pow(somega*sigma*d*t, 2) + pow(omega*ssigma*d*t, 2) + pow(omega*sigma*sd*t, 2) + 
                        pow(omega*sigma*d*st, 2));

    cout << "Il disallineamento tra i contatti è di (" << dx << " +- " << sdx << ") cm." << endl;

}
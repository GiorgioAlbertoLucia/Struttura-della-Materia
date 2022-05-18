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

void std_graph_settings_adv(TGraph& graph)
{
    gPad->SetTopMargin(0.15);
    graph.SetMarkerStyle(21);
    graph.SetMarkerSize(0.3);
    graph.SetLineColor(38);
    graph.SetLineWidth(4);
}

void analysis4()
{   
    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis4.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1110);
    gStyle->SetFitFormat("2.2e");


    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const float R_H = 0.687329;                  // V T / (A cm)
    const float sR_H = 3.55239e-02;

    const float t = 0.1;                        // z axis cm
    const float d = 1.;                         // y axis cm
    const float l = 2.;                         // x axis cm

    const char * path = "../data/mobilitàportatori.txt";
    int first_line = comment_lines(path);
    ifstream file(path);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> V, sV, i, si;
    float entry1, entry2, entry3, entry4;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path) == 2)
    {
        while (file >> entry1 >> entry2)
        {
            V.push_back(entry1);
            i.push_back(entry2);
        }
    }
    else if(count_column(path) == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            V.push_back(entry1);
            i.push_back(entry2);
        }
    }
    
    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    for(int j = 0; j < V.size(); j++)   
    {
        entry1 = abs(V.at(j)) * 0.02;
        sV.push_back(entry1);

        entry2 = abs(i.at(j)) * 0.006 + 0.02;
        si.push_back(entry2);
    }

    string str1("\tsV[V]"), str2("\tsi[mA]");
    if(names.find(str1) == string::npos)
    {
        names += str1;
        append_column(path, "\tsV[V]", sV);
    }
    if(names.find(str2) == string::npos)
    {
        names += str2;   
        append_column(path, "\tsi[mA]", si);
    } 

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();
    


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int j = 0; j < V.size(); j++)   
        cout << V.at(j) << "\t" << i.at(j) << "\t" << sV.at(j) << si.at(j) << endl;
    cout << endl;
    
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    
    TCanvas * canvas = new TCanvas("canvas", "mobilità", 500, 5, 500, 500);
    
    const int n[] = {0, 9, 17};    // position where each subset begins
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
        tf1->SetParameter(1, 0.051);
        tf1->SetParNames("V_{0} [V]", "R [M #Omega]");

        TGraphErrors * graph = new TGraphErrors(sub_i.size(), &sub_i[0], &sub_V[0], &sub_si[0], &sub_sV[0]);
        graph->SetTitle("Mobilita dei portatori;i [mA];V [V]");
        std_graph_settings_adv(*graph);
    
        graph->GetYaxis()->SetTitleOffset(2.3);
        gPad->SetGrid();
        gPad->SetTopMargin(0.13);
        gPad->SetLeftMargin(0.20);
        graph->Fit(tf1, "ER");
        graph->Draw("ap");
    
        cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
        " (Probability: " << tf1->GetProb() << ")." << endl;

        if(j==0)    canvas->SaveAs("../graphs/mobilità_1.pdf");
        if(j==1)    canvas->SaveAs("../graphs/mobilità_2.pdf");

        V0.push_back(tf1->GetParameter(0));         // V
        sV0.push_back(tf1->GetParError(0));
        R.push_back(tf1->GetParameter(1)*1000);     // ohm
        sR.push_back(tf1->GetParError(1)*1000);
    }
    

    ////////////////////////////////////// MOBILITÀ /////////////////////////////////////
    cout << "Z test, intercetta V0 compatibile con 0: " << endl;
    z_test(V0.at(0), 0, sV0.at(0));
    z_test(V0.at(1), 0, sV0.at(1));

    cout << "V0 medio: (" << (V0.at(0)+V0.at(1))/2 << " ± " << sqrt(sV0.at(0)*sV0.at(0) + sV0.at(1)*sV0.at(1)) << ") V" << endl;

    cout << "Resistenza: " << endl;
    for (int i = 0; i < R.size(); i++)  cout << " (" << R.at(i) << " ± " << sR.at(i) << ") Ω" << endl;
    
    const float sigma1 = (l)/(R.at(0) * t * d);                 // cm^-1 Ω^-1
    const float ssigma1 = l*sR.at(0)/(R.at(0)*R.at(0)*t*d);

    const float sigma2 = (l)/(R.at(1) * t * d);                 // cm^-1 Ω^-1
    const float ssigma2 = l*sR.at(1)/(R.at(1)*R.at(1)*t*d);

    const float mu1 = R_H * sigma1;                             // m^2 / (V s)
    const float smu1 = sqrt(pow(sR_H*sigma1, 2) + pow(R_H*ssigma1, 2));
    const float mu2 = R_H * sigma2;                             // m^2 / (V s)
    const float smu2 = sqrt(pow(sR_H*sigma2, 2) + pow(R_H*ssigma2, 2));

    z_test(mu1, mu2, sqrt(smu1*smu1 + smu2*smu2)); cout << " m^2 / (V s) " << endl;
    cout << endl << "µ exp = (" << (mu1+mu2)/2 << " ± " << sqrt(smu1*smu1 + smu2*smu2)/2 << ") m^2 / (V s) " << endl; 

    ///////////////////////////// DISALLINEAMENTO ///////////////////////////
    const float omega = 0.244617;                               // Ω (dalla calibrazione in analysis3)
    const float somega = 0.0014403;

    const float sigma = (sigma1 + sigma2)/2;                    // cm^-1 Ω^-1 
    const float ssigma = sqrt(ssigma1*ssigma1 + ssigma2*ssigma2)/2;

    const float dx = omega * sigma * d * t;
    const float sdx = sqrt(pow(somega*sigma*d*t, 2) + pow(omega*ssigma*d*t, 2));

    cout << "La conducibilità del materiale è sigma = (" << sigma << " ± " << ssigma << ") cm^-1 Ω^-1" << endl;
    cout << "Il disallineamento tra i contatti è di (" << dx << " ± " << sdx << ") cm." << endl;

}
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

    hello_world();

    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    /*
    //freopen("../output/analysis1.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");
    */

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    /*
    int first_line = comment_lines("../data/mappatura.txt") + 1;
    ifstream file("../data/mappatura.txt");
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> x, sx, y, sy, B, sB;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
    {
        Vin.push_back(entry1);
        Vout.push_back(entry2);
        CHN.push_back(entry3);
        err_Vin.push_back(entry4);
        err_Vout.push_back(entry5);
        err_CHN.push_back(entry6);
    }
    */


    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    /*
    for(int i = 0; i < chn_centroid.size(); i++)   
    {
        float Ei = a + b * chn_centroid.at(i) + c *chn_centroid.at(i) *chn_centroid.at(i);
        float err_Ei = sqrt(err_a*err_a + chn_centroid.at(i)*err_b*err_b + chn_centroid.at(i)*chn_centroid.at(i)*err_c*err_c 
                        + (b + c*chn_centroid.at(i)*err_centroid.at(i)*err_centroid.at(i)));
        E.push_back(Ei);
        err_E.push_back(err_Ei);
    }

    // if there are no columns named E and err_E, respective values are appended to the file
    string str_E("E"), str_err_E("err_energy");
    if(line.find(str_E) == string::npos)        append_column(name, "E", E);
    if(line.find(str_err_E) == string::npos)    append_column(name, "err_energy", err_E);
    */

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    /*
    file.close();
    */


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    /*
    cout << endl << "Dati:" << endl << line << endl;
    for (int i = 0; i < Vin.size(); i++)   
        cout << Vin.at(i) << "\t" << Vout.at(i) << "\t" << CHN.at(i) << "\t" << err_Vin.at(i) << "\t" << 
        err_Vout.at(i) << "\t\t" << err_CHN.at(i) << endl;
    cout << endl;
    
    TCanvas * canvas1 = new TCanvas("Vin", "Vout", 500, 5, 500, 600);
    canvas1->SetGrid();
    */
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    /*
    TF1 * tf1 = new TF1("tf1", "[0]+[1]*x", -15, 15);
    tf1->SetLineColor(38);

    TGraphErrors * graph1 = new TGraphErrors(Vin.size(), &Vin[0], &Vout[0], &err_Vin[0], &err_Vout[0]);
    graph1->SetTitle("#splitline{Catena Elettronica}{y = p_{0} + p_{1} x};Vin [V];Vout [V]");
    std_graph_settings(*graph1);
    
    graph1->Fit(tf1, "ER");
    graph1->Draw("ap");
    canvas1->SaveAs("../graphs/linearità1.jpg");
    canvas1->SaveAs("../graphs/linearità1.pdf");
    
    cout << "Chi^2:" << tf1->GetChisquare() << ", number of DoF: " << tf1->GetNDF() << 
    " (Probability: " << tf1->GetProb() << ")." << endl;
    */


}
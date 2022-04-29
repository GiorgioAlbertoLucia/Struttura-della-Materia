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
    const char * path1 = "../data/mobilitàportatori.txt";
    int first_line = comment_lines(path1);
    ifstream file(path1);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> x, sx, y, sy, B, sB;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(count_column(path1) == 4)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            B2.push_back(entry3);
            sB2.push_back(entry4);
        }
    }
    else if(count_column(path1) == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            VH2.push_back(entry1);
            sVH2.push_back(entry2);
            B2.push_back(entry3);
            sB2.push_back(entry4);
        }
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

    str1 = "\tVH2_correct[V]";
    str2 = "\tsVH2_correct[V]";
    if(names.find(str1) == string::npos)
    {
        names += "\tVH2_correct[V]";
        append_column(path1, "VH2_correct[V]", VH2_correct);
    }    
    if(names.find(str2) == string::npos)
    {
        names += "\tsVH2_correct[V]";   
        append_column(path1, "sVH2_correct[V]", sVH2_correct);
    } 

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    /*
    file.close();
    */


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    /*
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < Vin.size(); i++)   
        cout << Vin.at(i) << "\t" << Vout.at(i) << "\t" << CHN.at(i) << "\t" << err_Vin.at(i) << "\t" << 
        err_Vout.at(i) << "\t\t" << err_CHN.at(i) << endl;
    cout << endl;
    */
    

    /////////////////////////////// FIT //////////////////////////////////////////////////////
    /*
    TCanvas * canvas1 = new TCanvas("Vin", "Vout", 500, 5, 500, 600);
    canvas1->SetGrid();
    
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
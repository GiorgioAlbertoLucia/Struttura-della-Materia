#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TStyle.h>
#include <TROOT.h>

#include "useful_functions.cpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

using namespace std;

void analysis1_2()
{   
    ///////////////////// MAPPATURA DEL MAGNETE /////////////////////////////////////////////////


    ///////////////////// SET OUTPUT IN A FILE AND OTHER GENERAL SETTINGS ///////////////////////

    //freopen("../output/analysis1.txt", "w", stdout);
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1100);
    gStyle->SetFitFormat("2.2e");
    

    ///////////////////// READ DATA FROM A FILE ////////////////////////////////////////////////
    
    const char * path = "../data/mappatura2.txt";
    int first_line = comment_lines(path);
    ifstream file(path);
    for(int i=0; i<first_line; i++) file.ignore(10000, '\n');    

    vector<float> x, sx, y, sy, B, sB;
    float entry1, entry2, entry3, entry4, entry5, entry6;
    string names;
    getline(file, names);                                            // store the names of the variables

    if(comment_lines(path) == 5)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5)
        {
            x.push_back(entry1);
            sx.push_back(entry2);
            y.push_back(entry3);
            sy.push_back(entry4);
            B.push_back(entry5);
        }
    }
    if(comment_lines(path) == 6)
    {
        while (file >> entry1 >> entry2 >> entry3 >> entry4 >> entry5 >> entry6)
        {
            x.push_back(entry1);
            sx.push_back(entry2);
            y.push_back(entry3);
            sy.push_back(entry4);
            B.push_back(entry5);
        }
    }

    ///////////////////////////// ADD DATA ///////////////////////////////////////////////////////
    
    for(int i = 0; i < B.size(); i++)   
    {
        entry1 = B.at(i) * 0.05 + 0.1;
        sB.push_back(entry1);
    }

    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < x.size(); i++)   
        cout << x.at(i) << "\t" << sx.at(i) << "\t" << y.at(i) << "\t" << sy.at(i) << "\t" << 
        B.at(i) << "\t\t" << sB.at(i) << endl;
    cout << endl;
    

    string str1("\tsB[mT]");
    if(names.find(str1) == string::npos)
    {
        names += str1;
        append_column(path, "\tsB[mT]", sB);
    }

    vector<float> sub_x, sub_sx, sub_y, sub_sy, Bx, sBx, By, sBy;

    for(int i = 0; i < x.size(); i++)   
    {
        if(i < 9)
        {
            sub_y.push_back(y.at(i));
            sub_sy.push_back(sy.at(i));
            By.push_back(B.at(i));
            sBy.push_back(sB.at(i));
        }
        else if(i < 18)
        {
            sub_x.push_back(x.at(i));
            sub_sx.push_back(sx.at(i));
            Bx.push_back(B.at(i));
            sBx.push_back(sB.at(i));
        }
    }

    ////////////////////////// CLOSE FILE ///////////////////////////////////////////////////
    
    file.close();


    ////////////////////////// PRINT OUT DATA /////////////////////////////////////////////////
    
    cout << endl << "Dati:" << endl << names << endl;
    for (int i = 0; i < sub_x.size(); i++)   
        cout << sub_x.at(i) << "\t" << sx.at(i) << "\t" << sub_y.at(i) << "\t" << sy.at(i) << "\t" << 
        B.at(i) << "\t\t" << sB.at(i) << endl;
    cout << endl;
    
    /////////////////////////////// GRAPHS //////////////////////////////////////////////////////
    TCanvas * canvas1 = new TCanvas("canvas1", "B", 500, 5, 500, 600);
    canvas1->SetGrid();

    TGraphErrors * graph1 = new TGraphErrors(sub_x.size(), &sub_x[0], &Bx[0], &sub_sx[0], &sBx[0]);
    graph1->SetTitle("#splitline{Mappatura del}{campo magnetico};x [mm];B [mT]");
    std_graph_settings(*graph1);
    
    graph1->Draw("apl");
    canvas1->SaveAs("../graphs/mappatura_x2.jpg");
    canvas1->SaveAs("../graphs/mappatura_x2.pdf");

    TCanvas * canvas2 = new TCanvas("canvas2", "B", 500, 5, 500, 600);
    canvas2->SetGrid();

    TGraphErrors * graph2 = new TGraphErrors(sub_y.size(), &sub_y[0], &By[0], &sub_sy[0], &sBy[0]);
    graph2->SetTitle("#splitline{Mappatura del}{campo magnetico};x [mm];B [mT]");
    std_graph_settings(*graph2);
    
    graph2->Draw("apl");
    canvas2->SaveAs("../graphs/mappatura_y2.jpg");
    canvas2->SaveAs("../graphs/mappatura_y2.pdf");

    TCanvas * canvas3 = new TCanvas("canvas3", "B", 500, 5, 500, 600);
    canvas3->SetGrid();

    TGraph2D * graph3 = new TGraph2D(x.size(), &x[0], &y[0], &B[0]);
    graph3->SetTitle("#splitline{Mappatura del}{campo magnetico};x [mm];y [mm];B [mT]");
    
    graph3->Draw("surf1");
    canvas3->SaveAs("../graphs/mappatura_xy2.jpg");
    canvas3->SaveAs("../graphs/mappatura_xy2.pdf");

    /*
    Int_t n = x.size();
    TGraph2DErrors * graph4 = new TGraph2DErrors(n, &x[0], &y[0], &B[0], &sx[0], &sy[0], &sB[0]);
    graph4->SetTitle("#splitline{Mappatura del}{campo magnetico};x [mm];y [mm];B [mT]");
    
    graph4->Draw("err");
    canvas->SaveAs("../graphs/mappatura_xyerr.jpg");
    canvas->SaveAs("../graphs/mappatura_xyerr.pdf");
    */

}
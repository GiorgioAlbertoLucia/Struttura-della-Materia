#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

////////////////////////////////// FILE MANAGEMENT ///////////////////////////////////////////////////////

/**
 * @brief Counts how many lines begin with '#' in a goven file (those lines will be used as comment lines).
 * To move to the desired line, just use file.seekg().
 * @param file_name 
 * @return int 
 */
int comment_lines(const char * file_name)
{
    int comment = 0;
    string line;
    ifstream file(file_name);
    if(file.is_open())  while(file >> line) if(line.at(0)=='#') comment++;
    file.close();
    return comment;
}

/**
 * @brief Counts the number of entries in a file (comment lines are not included).
 * 
 * @param file_name 
 * @return int 
 */
int entries(const char * file_name)
{
    ifstream file(file_name);

    if(file.is_open())
    {
        int first_line = comment_lines(file_name);
        for(int i=0; i < first_line; i++)   file.ignore(10000, '\n'); 

        int entries = 0;
        string line;
        while(getline(file, line))  entries++;

        file.close();
        return entries;
    }
    else
    {
        std::cerr << "Error: unable to open file" << std::endl;
        file.close();
        return 0;
    }    
}

/**
 * @brief This function add a column of floats (with their description as top line) to an existing .txt file
 * @param file_name path to the file
 * @param col_name description of the column data
 * @param column data to add in the column
 */
void append_column(const char * file_name, const char * col_name, vector<float> column)
{
    fstream file(file_name, ios::in);
    string line;
    vector<string> file_lines;
    while(getline(file, line))  file_lines.push_back(line);     // fill a vector with the file content
    file.close();

    int comment = comment_lines(file_name);

    file.open(file_name, ios::out);
    for (int i = 0; i < file_lines.size(); i++)
    {
        if (i < comment)        file << file_lines.at(i) << endl;
        else if (i == comment)  file << file_lines.at(i) << "\t\t" << col_name << endl;
        else                    file << file_lines.at(i) << "\t\t\t" << column.at(i-comment-1) << endl;                          
    }
    file.close();
}

/**
 * @brief Counts how many columns (text separated by a space, tab or comma) there are in a given string.
 * If different lines of the file have a different amount of columns, an error is displayed.
 * @return int number of columns.
 */
int count_column(const char * file_name)
{
    int columns = 0;
    int save_columns = 0;       // used to check if each line has the same number of columns

    std::ifstream file;
    file.open(file_name);

    if(file.is_open())
    {   
        int first_line = comment_lines(file_name);
        for(int i=0; i < first_line; i++)   file.ignore(10000, '\n'); 

        for (int j = 0; j < entries(file_name); j++)
        {
                std::string row, item;
                getline(file, row);
                std::istringstream iss(row);
                while(iss >> item)   if(item.length())  columns++;

                if (j > 0 && columns != save_columns)
                {
                    std::cerr << "Error: not all the lines have the same amount of columns. First " << j << " lines have "
                    << save_columns << " columns (this value will be returned)." << std::endl;
                    return save_columns;
                }
            
                save_columns = columns;
                columns = 0;
        }
    }
    return save_columns;
}

//////////////////////////////// GENERAL SETTINGS //////////////////////////////////////////////////////

/**
 * @brief Standard settings used for TGraph objects.
 * @param graph 
 */
void std_graph_settings(TGraph& graph)
{
    graph.GetYaxis()->SetTitleOffset(1.4);
    gPad->SetTopMargin(0.15);
    gPad->SetLeftMargin(0.15);
    graph.SetMarkerStyle(21);
    graph.SetMarkerSize(0.3);
    graph.SetLineColor(38);
    graph.SetLineWidth(4);
}

//////////////////////////////////// Z TEST /////////////////////////////////////////////////////////////

/**
 * @brief Gaussian test. The function prints the value of the variable Z and the result of the test
 * @param exp_value 
 * @param ref_value 
 * @param err 
 */
void z_test(const float exp_value, const float ref_value, const float err)
{
    float z = (exp_value - ref_value) / err;
    if (abs(z) < 1.96)  cout << exp_value << " +/- " << err << "\t" << ref_value << "\t"  << z << "\t" << "True" << endl;
    else                cout << exp_value << " +/- " << err << "\t" << ref_value << "\t"  << z << "\t" << "False" << endl;
}

/**
 * @brief This function takes a file in ../doc (no extension) where a TTree is stored. Executes a Z Test over all the values
 * stored in a specific branch in respect to the reference values stored in another branch (uncertainties are also stored in
 * some other branch).
 * 
 * @param file_name name of the .root file in ../doc (no extension)
 * @param exp_branch name of the branch where the experimental values are stored
 * @param ref_branch name of the branch where the reference values are stored
 * @param err_branch name of the branch where the experimental uncertainties are stored
 * @param name_branch name of the branch where the description of the data in the line is stored
 */
void z_test_branch(const char * file_name)
{
    ifstream file(file_name);

    string name;
    float entry1, entry2, entry3, ref, exp, err;

    string line;
    getline(file, line);    // skip first line

    cout << endl << "test z dal file " << file_name << endl;
    while (file >> name >> entry1 >> entry2 >> entry3 >> ref >> exp >> err)
    {
        cout << name << "\t";
        z_test(exp, ref, err);
    }
    file.close();
}

//////////////////////////////// HELLO WORLD ////////////////////////////////////////////////////////////
void hello_world()
{
    cout << "Hello world!" << endl;
}

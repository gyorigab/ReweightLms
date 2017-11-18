#ifndef LMS_SETTINGS_H
#define LMS_SETTINGS_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include <settings.h>


using namespace std;

template<typename Float=double, typename Index=std::size_t>
class LmsSettings: public Settings<>
{
public:

    LmsSettings()
    {
        file_A = new char[2048];
        file_l = new char[2048];
        file_p = new char[2048];

        sample_proc    = 95;
        remove_rows    = 5;
        maxiter        = 1000;
        extreme_weight = 10e-5;
        method      = "ExtremeWeight";
    }

    ~LmsSettings()
    {
        delete[] file_A;
        delete[] file_l;
        delete[] file_p;
    }

    int read(const char* file)
    {
        vector<ParVal> parval;

        if (read_file(file,parval)) return -1;

        for(Index i=Index(); i<parval.size(); i++)
        {
            if(parval[i].parameter == "DESIGN_MATRIX_FILE" )
            {
               // file_A = new char[parval[i].value.length()+1];
                strcpy(file_A,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "MEASUREMENT_VECTOR_FILE" )
            {
              //  file_l = new char[parval[i].value.length()+1];
                strcpy(file_l,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "WEIGHT_MATRIX" )
            {
              //  file_l = new char[parval[i].value.length()+1];
                strcpy(file_p,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "MAX_ITER" )
            {
                maxiter = atoi(parval[i].value.c_str());
            }
            else if (parval[i].parameter == "RANDOM_SAMPLE" )
            {
                sample_proc = atoi(parval[i].value.c_str());
            }
            else if (parval[i].parameter == "REMOVED_LINES" )
            {
                remove_rows = atoi(parval[i].value.c_str());
            }
            else if (parval[i].parameter == "METHOD" )
            {
                method = parval[i].value;
            }
            else if (parval[i].parameter == "EXTREME_WEIGHT" )
            {
                extreme_weight = atof(parval[i].value.c_str());
            }
            else
            {
                cerr << "Wrong settings parameter \n" << endl;
            }
        }
        return 0;
    }

    void print()
    {
        cout << "========LMS Settings======== \n";

        cout << "DESIGN_MATRIX_FILE = "      << file_A         << endl
             << "MEASUREMENT_VECTOR_FILE = " << file_l         << endl
             << "WEIGHT_MATRIX = "           << file_p         << endl
             << "MAX_ITER = "                << maxiter        << endl
             << "RANDOM_SAMPLE = "           << sample_proc    << endl
             << "REMOVED_LINES = "           << remove_rows    << endl
             << "METHOD = "                  << method         << endl
             <<  "EXTREME_WEIGHT = "         << extreme_weight << endl;

        cout << "=============================\n\n";

    }

    char* get_file_A(){ return file_A;}
    char* get_file_l(){ return file_l;}
    char* get_file_p(){ return file_p;}

    Index get_sample_proc()   { return sample_proc;}
    Index get_remove_rows()   { return remove_rows;}
    Index get_maxiter()       { return maxiter;}
    Float get_extreme_weight(){ return extreme_weight;}

    string get_method(){ return method;}

private:

    char* file_A;
    char* file_l;
    char* file_p;

    Index sample_proc;
    Index remove_rows;
    Index maxiter;

    string method;
    Float extreme_weight;
};

#endif // LMS_SETTINGS_H

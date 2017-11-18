#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

template<typename Float=double, typename Index=std::size_t>
class Settings
{
public:

    struct ParVal
    {
        ParVal():parameter(""),value(""){}

        string parameter;
        string value;

        void clear(){parameter=value="";}
    };

    Settings(){}
    ~Settings(){}

    int read_file(const char* file, vector<ParVal>& parval)
    {
        ifstream gin(file);
        string line;
        ParVal  pv;

        if(!gin)
        {
            cerr << "Can't read settings file \n\n";
            return -1;
        }

        while(getline(gin,line))
        {
            for(Index i=Index(); i < line.length(); i++)
            {
                if(line[i]=='=')
                {
                    for(Index j=i+1; j < line.length();j++,i++)
                    {
                        pv.value += line[j];
                    }
                }
                else
                {
                    pv.parameter += line[i];
                }
            }
            parval.push_back(pv);
            pv.clear();
        }
        return 0;
    }

private:

};

#endif // SETTINGS_H

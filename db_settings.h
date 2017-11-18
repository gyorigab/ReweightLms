#ifndef DB_SETTINGS_H
#define DB_SETTINGS_H

#include <vector>
#include <fstream>
#include <iostream>

#include <settings.h>

using namespace std;

template<typename Float=double, typename Index=std::size_t>
class DbSettings: public Settings<>
{
public:

    DbSettings()
    {
        host  = new char[2048];
        port  = new char[2048];
        db    = new char[2048];
        login = new char[2048];
        pass  = new char[2048];
    }

    ~DbSettings()
    {
        delete[] host;
        delete[] port;
        delete[] db;
        delete[] login;
        delete[] pass;
    }

    int read(const char* file)
    {
        vector<ParVal> parval;

        if ( read_file(file,parval)) return -1;

        for(Index i=Index(); i<parval.size(); i++)
        {
            if(parval[i].parameter == "HOST" )
            {
                strcpy(host,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "PORT" )
            {
                strcpy(port,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "DATABASE" )
            {
                strcpy(db,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "USER" )
            {
                strcpy(login,parval[i].value.c_str());
            }
            else if (parval[i].parameter == "PASSWORD" )
            {
                strcpy(pass,parval[i].value.c_str());
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
        cout << "========DB Settings======== \n";

        cout << "HOST = "     << host  << endl
             << "PORT = "     << port  << endl
             << "DATABASE = " << db    << endl
             << "LOGIN = "    << login << endl
             << "PASS = "     << pass  << endl;

        cout << "===========================\n";

    }

    char* get_host(){ return host; }
    char* get_port(){ return port; }
    char* get_db(){ return db; }
    char* get_login(){return login;}
    char* get_pass(){return pass;}

private:

     char* host;
     char* port;
     char* db;
     char* login;
     char* pass;

     vector<ParVal> parval;
};

#endif // DB_SETTINGS_H

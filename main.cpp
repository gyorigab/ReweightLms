#include <smatrix.h>
#include <iostream>
#include <algorithm>
#include <matvec.h>
#include <operators.h>
#include <vec.h>
#include <ctime>
#include <conjungate_gradients_nr.h>
#include <median.h>
#include <data.h>
#include <sample_model.h>
#include <general_cmp.h>
#include <settings.h>
#include <logger.h>
#include <db_results.h>
#include <stdlib.h>
#include <lms_settings.h>
#include <db_settings.h>
#include <lms_exact.h>

using namespace std;
using namespace GNU_gama;

int help();
int usage();
int version();

using namespace std;

int Logger::verb = 0;
ostream& Logger::out = std::cout;
/*
int main(int argc, char *argv[])
{
    srandom(time(0));

    Data<> d;
    LmsSettings<> set;

    if ( set.read("/home/gabriel/Doktstu/ds/repozitar/lms/settings/test1.set") ) return -1 ;
    set.print();

    //cout << "====== Nactene data =====" << endl;

    if(d.load_l(set.get_file_l())) return -1;
    if(d.load_P(set.get_file_p())) return -1;
    if(d.load_A(set.get_file_A())) return -1;

    cout << d.get_l() << endl;
    cout << d.get_P() << endl;
    cout << d.get_A() << endl;

    GNU_gama::SparseMatrix<>* ukA = d.get_A();

    double* lb = d.l_beg();
    double* le = d.l_end();

    double* Lb = d.P_beg();
    double* Le = d.P_end();

    LmsExact<> lmsexact(ukA, lb, le, Lb, Le);

    int n = 20;
    int k = 8;
    lmsexact.get_combinations(n,k);

    cout << "Min median: " << lmsexact.get_min_median() << endl;
}
*/


int main(int argc, char *argv[])
{
     srandom(time(0));

    time_t start;
    time_t end;
    double seconds;

    time(&start);

    DbResults db;

    if (argc < 2){return usage();}

    const std::string namearg(argv[1]);

    if      (namearg == "--help")    { return help();    }
    else if (namearg == "--version") { return version(); }

    const char* c;
    const char* argv_dbset;
    bool argv_dbsetb=false;

    for(int i=2; i<argc;i++)
    {
        c = argv[i];
        if (*c != '-') return help();
        if (*c && *c == '-') c++;
        if (*c && *c == '-') c++;

        const std::string name(c);

        c = argv[++i];

        if(name == "help") return help();
        else if(name == "dbset")  { argv_dbset=c; argv_dbsetb=true; }
        else return help();
    }

    DbSettings<> dbset;

    if(argv_dbsetb)
    {
        if( dbset.read(argv_dbset) ) return -1;

        dbset.print();
        db.connect(dbset);
    }
    else
    {
        db.connect("127.0.0.1","5432","lms","lms","lms");
    }

    if (namearg == "--dbinit")
    {
        return db.resetSchemas();
    }

    Data<> d;
    LmsSettings<> set;

    if ( set.read(argv[1]) ) return -1 ;
    set.print();

    db.insertSettings(set);

    //cout << "====== Nactene data =====" << endl;

    if(d.load_P(set.get_file_p())) return -1;
    if(d.load_l(set.get_file_l())) return -1;
    if(d.load_A(set.get_file_A())) return -1;

    //cout << d.get_P() << endl;
    //cout << d.get_l() << endl;
    //cout << d.get_A() << endl;

    GNU_gama::SparseMatrix<>* ukA = d.get_A();

    double* lb = d.l_beg();
    double* le = d.l_end();
    double* Pb = d.P_beg();
    double* Pe = d.P_end();

    GeneralCmp<> gc(ukA, lb, le, Pb, Pe);

    gc.homogenize();

    //cout << "====== Data po homogenizacii ======" << endl;
    //cout << d.get_A() << endl;
    //cout << d.get_l() << endl;

    int rows = d.get_rows();
    int cols = d.get_cols();

    GNU_gama::Vec<> x(cols);
    GNU_gama::Vec<> r(rows);

    x.set_zero();

    double* xb = x.begin();
    double* xe = x.end();

    double* rb = r.begin();
    double* re = r.end();

    //cout << "====== Vypocet medianu ======" << endl;

    gc.set_setting(&set);
    gc.set_database(&db);

    gc.solveLms(rb,re,xb,xe);

    db.closeConnection();

    time(&end);

    seconds = difftime(end,start);

    cout << " Duration of computation (sec) " << seconds << endl;


    return 0;
}


int help()
{
    cout << "Geodetic networks adjustment with Least Median Squares method: \n";
    cout << "\n" << endl;
    cout << "Usage: lms LMS_SETTINGS_FILE [options] \n";
    cout << "Usage: lms --help \n";
    cout << "Usage: lms --version \n";
    cout << "Usage: lms --dbinit [--dbset DB_SETTINGS_FILE] \n\n";
    cout << "Options: \n\n"
         << "--dbset DB_SETTINGS_FILE - db connection strings settings file \n"
         << "                           default: localhost, 5432, lms, lms, lms \n"
         << "\n";

    return 0;
}

int usage()
{
    cout << "Usage: lms SETTINGS_FILE [options] \n";
    cout << "Try 'lms --help' for more information. \n";
    return -1;
}

int version()
{
    cout << "Version: 1.0 Author: \"Gabriel Gyori 2014\" \n" << endl;

    return 0;
}

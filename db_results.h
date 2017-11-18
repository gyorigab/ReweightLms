#ifndef DB_RESULTS_H
#define DB_RESULTS_H

#include <lms_settings.h>
#include <sstream>
#include <stdlib.h>
#include <db_settings.h>

extern "C" {
#include <libpq-fe.h>
}

class DbResults
{
public:

    DbResults(){}
    ~DbResults(){}

    int connect(const char* host, const char* port, const char* db , const char* login, const char* pwd)
    {
        conn = PQsetdbLogin(host,port,NULL,NULL,db,login,pwd);

        if(PQstatus(conn) != CONNECTION_OK)
        {
            cout << "Konexi do databaze: " << db
                 << "\nUzivatel: " << login
                 << "\nHost: " << host
                 << "\nPort: " << port
                 << "\nsa nepodarilo naviazat " << endl;
            return -1;
        }
        else
        {
            cout << "===================================== \n";
            cout << "    Connection to database OK " << endl;
            cout << "===================================== \n";
            cout << "Database:           " << PQdb(conn) << endl;
            cout << "User:               " << PQuser(conn) << endl;
            cout << "Host:               " << PQhost(conn)  << endl;
            cout << "Port:               " << PQport(conn) << endl;
            cout << "Db server version:  " << PQserverVersion(conn) << endl;
            cout << "===================================== \n\n";

            return 0;
        }
    }

    int connect(DbSettings<>& set)
    {
        const char* host=set.get_host();
        const char* port=set.get_port();
        const char* db=set.get_db();
        const char* login=set.get_login();
        const char* pass=set.get_pass();

        connect(host,port,db,login,pass);

        return 0;
    }

    int createSchemas()
    {
        PGresult* res;

        const char* settings_tab="CREATE TABLE settings("
                " id SERIAL PRIMARY KEY, "
                " design_matrix_file VARCHAR(200),"
                " measurement_vector_file VARCHAR(200),"
                " weight_matrix_file VARCHAR(200),"
                " max_iter INTEGER,"
                " random_sample INTEGER,"
                " removed_lines INTEGER,"
                " method VARCHAR(50),"
                " extreme_weight FLOAT"
                " );";

        const char* median_cmp_tab="CREATE TABLE median_cmp("
                " settings_id INTEGER REFERENCES settings(id),"
                " cycle INTEGER NOT NULL,"
                " step INTEGER NOT NULL,"
                " median_before DOUBLE PRECISION NOT NULL,"
                " median_after DOUBLE PRECISION NOT NULL, "
                " max_residual DOUBLE PRECISION NOT NULL, "
                " min_residual DOUBLE PRECISION NOT NULL, "
                " abs_residual DOUBLE PRECISION NOT NULL "
                " );";

        const char* median_tab="CREATE TABLE median("
                " settings_id INTEGER REFERENCES settings(id),"
                " min_median   DOUBLE PRECISION NOT NULL, "
                " avg_median   DOUBLE PRECISION NOT NULL, "
                " cmp_time TIMESTAMP DEFAULT(now())"
                " );";

        const char* reweighted_rows_tab="CREATE TABLE reweighted_rows( "
                " settings_id INTEGER REFERENCES settings(id), "
                " removed_row INTEGER"
                " );";

        const char* estimated_unkonwns_tab= "CREATE TABLE estimated_unknowns( "
                " settings_id INTEGER REFERENCES settings(id), "
                " row INTEGER NOT NULL,"
                " value DOUBLE PRECISION NOT NULL"
                " );";

        const char* estimated_residuals_tab= "CREATE TABLE estimated_residuals( "
                " settings_id INTEGER REFERENCES settings(id), "
                " row INTEGER NOT NULL,"
                " value DOUBLE PRECISION NOT NULL"
                " );";

        const char* reweighted_rows_idx="CREATE INDEX rewighed_rrows ON reweighted_rows_tab COLUMN removed_row;";

        PQexec(conn,settings_tab);
        PQexec(conn,median_tab);
        PQexec(conn,median_cmp_tab);
        PQexec(conn,reweighted_rows_tab);
        PQexec(conn,reweighted_rows_idx);
        PQexec(conn,estimated_unkonwns_tab);
        PQexec(conn,estimated_residuals_tab);

        return 0;
    }

    int dropSchemas()
    {
        const char* settings_tab= "DROP TABLE settings CASCADE;";
        const char* median_cmp_tab= "DROP TABLE median_cmp CASCADE;";
        const char* median_tab= "DROP TABLE median CASCADE;";
        const char* reweighted_rows_tab= "DROP TABLE reweighted_rows CASCADE;";
        const char* estimated_unkonwns_tab="DROP TABLE estimated_unknowns CASCADE";
        const char* estimated_residuals_tab="DROP TABLE estimated_residuals CASCADE";

        PQexec(conn,settings_tab);
        PQexec(conn,median_tab);
        PQexec(conn,median_cmp_tab);
        PQexec(conn,reweighted_rows_tab);
        PQexec(conn,estimated_unkonwns_tab);
        PQexec(conn,estimated_residuals_tab);

        return 0;
    }


    int resetSchemas()
    {
        dropSchemas();
        createSchemas();
        return 0;
    }

    int insertSettings(LmsSettings<>& set)
    {
        vector<string> setvalues;
        vector<string>::iterator it;

        setvalues.push_back("\'"+string(set.get_file_A())+"\'");
        setvalues.push_back("\'"+string(set.get_file_l())+"\'");
        setvalues.push_back("\'"+string(set.get_file_p())+"\'");
        setvalues.push_back(toString(set.get_maxiter()));
        setvalues.push_back(toString(set.get_sample_proc()));
        setvalues.push_back(toString(set.get_remove_rows()));
        setvalues.push_back("\'"+string(set.get_method())+"\'");
        setvalues.push_back(toString(set.get_extreme_weight()));

        string insert_s="INSERT INTO settings("
                "design_matrix_file,"
                "measurement_vector_file,"
                "weight_matrix_file,"
                "max_iter,"
                "random_sample,"
                "removed_lines,"
                "method,"
                "extreme_weight"
                ") "
                "VALUES(";

        for(it=setvalues.begin(); it!=setvalues.end()-1;it++)
        {
            insert_s += *it+", ";
        }

        insert_s += *it+");";

        //cout << insert_s << endl;

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        const char * cmp_id_query="SELECT currval('settings_id_seq'::regclass);";

        PGresult* res;

        res = PQexec(conn,cmp_id_query);   // Vratim poslednu hodnotu sekvencie v aktualnej session

        //cout << "Pocet riadkov dotazu: " << PQntuples(res) << endl;
        //cout << "Pocet stlpcov dotazu: " << PQnfields(res) << endl;

        cmp_id = atoi(PQgetvalue(res,0,0));

        //cout << "cmpid " << cmp_id << endl;

        return 0;
    }

    int insertMedian(double min_median, double avg_median)
    {
        string insert_s="INSERT INTO median("
                "settings_id,"
                "min_median,"
                "avg_median"
                ") "
                "VALUES(";

        insert_s += toString(cmp_id) + ", ";
        insert_s += toString(min_median) + ", ";
        insert_s += toString(avg_median) + ");";

        //cout << insert_s << endl;

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        return 0;
    }

    int insertMedianCmp(int cycle, int step, double median_before, double median_after,
                        double max_residual, double min_residual, double abs_residual)
    {

        string insert_s="INSERT INTO median_cmp("
                "settings_id,"
                "cycle,"
                "step,"
                "median_before,"
                "median_after,"
                "max_residual,"
                "min_residual,"
                "abs_residual"
                ") "
                "VALUES(";

        insert_s += toString(cmp_id) + ", ";
        insert_s += toString(cycle) + ", ";
        insert_s += toString(step) + ", ";
        insert_s += toString(median_before) + ", ";
        insert_s += toString(median_after) + ", ";
        insert_s += toString(max_residual) + ", ";
        insert_s += toString(min_residual) + ", ";
        insert_s += toString(abs_residual) + ");";

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        return 0;
    }

    int insertReweightedRows(int reweighted_rows)
    {
        string insert_s="INSERT INTO reweighted_rows("
                "settings_id,"
                "removed_row"
                ") "
                "VALUES(";

        insert_s += toString(cmp_id) + ", ";
        insert_s += toString(reweighted_rows) + ");";

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        return 0;
    }

    int insertEstimatedUnknowns(int row ,double value)
    {
        string insert_s="INSERT INTO estimated_unknowns("
                "settings_id,"
                "row,"
                "value"
                ") "
                "VALUES(";

        insert_s += toString(cmp_id) + ", ";
        insert_s += toString(row) + ", ";
        insert_s += toString(value) + ");";

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        return 0;
    }

    int insertEstimatedResiduals(int row ,double value)
    {
        string insert_s="INSERT INTO estimated_residuals("
                "settings_id,"
                "row,"
                "value"
                ") "
                "VALUES(";

        insert_s += toString(cmp_id) + ", ";
        insert_s += toString(row) + ", ";
        insert_s += toString(value) + ");";

        const char* insert=insert_s.c_str();

        PQexec(conn,insert);

        return 0;
    }

    int closeConnection()
    {
        PQfinish(conn);
        return 0;
    }

    //Konverzni funkce
    template<typename T>
    string toString(T number)
    {
       std::stringstream ss;
       ss << number;
       return ss.str();
    }

private:
    PGconn *conn;          // pointer na jedno backend spojenie
    int cmp_id;            // computation id = current value id in table settings

};

#endif // DB_RESULTS_H

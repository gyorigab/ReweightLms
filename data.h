#ifndef DATA_H
#define DATA_H

#include <smatrix.h>
#include <fstream>
#include <matvec/matvec.h>

template<typename Float = double, typename Index=std::size_t>
class Data
{
public:
    Data()
    {
        rows   = 0;
        cols   = 0;
        nozero = 0;
    }
    ~Data()
    {
        // DOPLNIT DESTRUCTOR
    }

    int load_sparse_A(char* file)
    {

        if (check_sparse(file,rows,cols,nozero)) return -1;

        A = new SparseMatrix<Float, Index>(nozero,rows,cols);

        ifstream gin(file);

        read_sparse(gin,A);

        return 0;
    }

    int load_A(char* file)
    {

        if(check(file,rows,cols,nozero)) return -1 ;  // funkcia vrati dimenziu matice a pocet nenulovych prvkov
                                                         // matica je nacitana zo suboru
        A = new SparseMatrix<Float, Index>(nozero,rows,cols);

        ifstream gin(file);

        gin >> A;            // naplnim riedku maticu

        return 0;
    }

    int load_P(char* file)
    {
        ifstream gin(file);

        if(!gin)
        {
            cerr << "Can't read covariance matrix file \n\n";
            return -1;
        }

        gin >> P;

        return 0;
    }

    int load_l(char* file)
    {
        ifstream gin(file);

        if(!gin)
        {
            cerr << "Can't read measuremment vector file \n\n";
            return -1;
        }

        gin >> l;

        return 0;
    }

    GNU_gama::SparseMatrix<>* get_A(){ return A; }
    GNU_gama::Vec<> get_P()          { return P; }
    GNU_gama::Vec<> get_l()          { return l; }

    Float* P_beg(){ return P.begin();}
    Float* P_end(){ return P.end();  }

    Float* l_beg(){ return l.begin();}
    Float* l_end(){ return l.end();  }

    int get_rows()    {return rows;}
    int get_cols()    {return cols;}
    int get_nonzero() {return nozero;}

private:

    int rows;
    int cols;
    int nozero;

    GNU_gama::SparseMatrix<Float,Index>* A;
    GNU_gama::Vec<Float> P;
    GNU_gama::Vec<Float> l;
};

#endif // DATA_H

#include <stdlib.h>
#include "operators.h"

int check(char* subor,int &m, int &n, int &nonzero){

     ifstream dat(subor);            //otvori subor
     if(!dat)
     {
         cerr << "Can't read input matrix \n";
         return -1;
     }

     double d;
     m = n = nonzero = 0;

     dat >> m >> n;                          // nacita dimenziu do premennych (neuklada do struktury)

     while(dat >> d)
     {
         if(d != 0.0)
         {
             nonzero++;                      // spocita pocet nenulovych prvkov
         }
     }
   dat.close();

   return 0;
}

int check_sparse(const std::string subor, int &m, int &n, int &nonzero)
{

    ifstream dat(subor.c_str());            //otvori subor
    if(!dat)
    {
        cerr << "Can't read input matrix \n";
        return -1;
    }
    m = n = nonzero = 0;

    dat >> m >> n >> nonzero;

    dat.close();
    return 0;
}

std::istream& operator>>(std::istream &data, SparseMatrix<> *A)
{
    double val;
    int m,n;
    data >> m >> n;

    for(int i=1; i<=m; i++)
    {
        A->new_row();
        for(int j=1; j<=n; j++)
        {
         data >> val;

            if(val != 0.0)
            {
                A->add_element(val, j);
            }
        }
    }
    return data;
}

std::ostream& operator<<(std::ostream &data, SparseMatrix<> *A)
{
      data << std::endl;
      for (unsigned long k=1; k<=A->rows(); k++)
        {
          data << k << " : ";

          double* n = A->begin(k);
          double* e = A->end  (k);

          for(std::size_t* i=A->ibegin(k) ; n!=e; n++, i++)
            {
              data << *n << " [" << *i << "]  ";
            }
          data << std::endl;
        }
      return data;
}

std::istream& read_sparse(istream &data, SparseMatrix<> *A)
{
    int m,n,nonzero;
    int i,j;
    double val;

    int j0=0;

    data >> m >> n >> nonzero;

    while(data >> j >> i >> val)
    {
        if(j!=j0)
        {
            A->new_row();
            j0=j;
        }
        A->add_element(val,i);
    }
    return data;
}

#ifndef LMS_EXACT_H
#define LMS_EXACT_H

#include <smatrix.h>
#include <vector>
#include <conjungate_gradients_nr.h>

template<typename Float=double, typename Index=std::size_t,
          typename const_iterator=GNU_gama::Vec<>::const_iterator,
          typename       iterator=GNU_gama::Vec<>::iterator>
class LmsExact
{

public:

    enum{NoError, BadDimension, EmptySet, BadAlloc, RepeatedRollBack};

    LmsExact()
    {
        A    = 0;
        Acpy = 0;
        min  = cols;
        rows = cols = 0;
    }

    ~LmsExact()
    {
        if(Acpy)
        {
            delete Acpy;
        }
    }

    LmsExact( const GNU_gama::SparseMatrix<Float,Index>* Ap,
              const_iterator lbegin, const_iterator lend,
              const_iterator Lbegin, const_iterator Lend)
    {
        A = Ap;

        l_begin = lbegin;
        l_end   = lend;

        L_begin = Lbegin;
        L_end   = Lend;

        rows  = A->rows();
        cols  = A->columns();
        min = cols;

        Acpy  =  new SparseMatrix<Float, Index>(cols*cols,cols,cols);

        lcpy.clear();
        Lcpy.clear();
    }

    void pretty_print(const vector<int>& v)
    {
      static int count = 0;
      cout << "combination no " << (++count) << ": [ ";
      for (int i = 0; i < v.size(); ++i) { cout << v[i] << " "; }
      cout << "] " << endl;
    }

    int get_combination_vectors(int offset, int k)
    {
        if (k == 0)
        {
            pretty_print(combination);

            // Create sample of the vector l

            lcpy.clear();

            for (Index i=Index(); i < cols; i++)
            {
                lcpy.push_back(*(l_begin+combination[i]-1));
            }

            Lcpy.clear();

            for (Index i=Index(); i < cols; i++)
            {
                Lcpy.push_back(*(L_begin+combination[i]-1));
            }

            cout << "l sample" << "\n\n";
            for(int i=0; i<lcpy.size() ; i++)
                cout << lcpy[i] << "\n";

            cout << "L sample" << "\n\n";
            for(int i=0; i<Lcpy.size() ; i++)
                cout << Lcpy[i] << "\n";

            // Create sample of the matrix A

             Acpy->reset(cols*cols,cols,cols);

            for (int j=0; j< cols ; j++)
            {
                Float* n = A->begin(combination[j]);
                Float* e = A->end  (combination[j]);

                Acpy->new_row();

                for(Index* i = A->ibegin(combination[j]) ; n!=e; n++, i++)
                {
                    Acpy->add_element(*n, *i);
                }
            }

            cout << Acpy << "\n\n";

            vector<Float> x(cols);
            typename vector<Float>::iterator x_beg = x.begin();
            typename vector<Float>::iterator x_end = x.end();

            vector<Float> X(cols);
            typename vector<Float>::iterator X_beg = X.begin();
            typename vector<Float>::iterator X_end = X.end();

            ConjungateGradientsNr<Float, Index> cjg(Acpy);

            typename vector<Float>::const_iterator lcpy_beg = lcpy.begin();
            typename vector<Float>::const_iterator lcpy_end = lcpy.end();

            typename vector<Float>::const_iterator Lcpy_beg = Lcpy.begin();
            typename vector<Float>::const_iterator Lcpy_end = Lcpy.end();

            // solve linear equations by conjugate gradients
            cjg.solve(lcpy_beg, lcpy_end, x_beg, x_end);
            cjg.solve(Lcpy_beg, Lcpy_end, X_beg, X_end);


            cout << "Vector x:" << "\n\n";
            for(int i=0; i<x.size() ; i++)
                cout << x[i] << endl;

            cout << "Vector X:" << "\n\n";
            for(int i=0; i<x.size() ; i++)
                cout << X[i] << endl;

            vector<Float> r;

            // compute real residuals
            for(int i=0; i<x.size() ; i++) r.push_back((X[i]-x[i])*(X[i]-x[i]));

            /*
            cout << "Vector r:" << "\n\n";
            for(int i=0; i<r.size() ; i++)
                 cout << r[i] << endl;
            */

            typename vector<Float>::iterator r_beg = r.begin();
            typename vector<Float>::iterator r_end = r.end();

            Index mid = r.size()/2;

            Float median = 0.0;

            sort(r_beg, r_end);


            cout << "Vector r sorted:" << "\n\n";
            for(int i=0; i<r.size() ; i++)
                cout << r[i] << endl;


            if(r.size()%2 == 0) { median = (r[mid-1] + r[mid])/2; }
            else                { median = r[mid]; }


            cout << "Median: " << median << "\n\n";

            medians.push_back(median);

            return 0;
        }
        for (int i = offset; i <= elements.size() - k; ++i)
        {
            combination.push_back(elements[i]);
            get_combination_vectors(i+1, k-1);
            combination.pop_back();
        }
    }

    int get_combinations(int n, int k)
    {

        for (int i = 0; i < n; ++i) { elements.push_back(i+1); }
        get_combination_vectors(0, k);

        return 0;
    }

    double get_min_median()
    {

        double median = medians[0];

        for(int i = 1; i<medians.size() ; i++)
        {
            if(medians[i] < median) { median = medians[i]; }
        }

        return median;
    }

private:

    const GNU_gama::SparseMatrix<Float,Index>* A;
    GNU_gama::SparseMatrix<Float,Index>* Acpy;


    Index rows;                  // pocet riadkov matice planu
    Index cols;                  // pocet stlpcov matice planu
    Index min;                   // minimalny pocet merani pro vzorek (pocet merani bez nadbytocnych merani)

    vector<Float> lcpy;           // kopie prave strany
    vector<Float> Lcpy;           // kopie prave strany

    const_iterator l_begin;
    const_iterator l_end;

    const_iterator L_begin;
    const_iterator L_end;

    vector<int> elements;
    vector<int> combination;

    vector<Float> medians;

    int lasterror;

};


#endif // LMS_EXACT_H

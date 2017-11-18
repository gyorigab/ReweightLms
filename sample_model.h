#ifndef SAMPLE_MODEL_H
#define SAMPLE_MODEL_H

#include <smatrix.h>
#include <matvec/matvec.h>
#include <set>
#include <vector>
#include <logger.h>

using namespace std;

template<typename Float=double, typename Index=std::size_t,
          typename const_iterator=GNU_gama::Vec<>::const_iterator,
          typename       iterator=GNU_gama::Vec<>::iterator>
class SampleModel
{
public:

    enum{NoError, BadDimension, EmptySet, BadAlloc, RepeatedRollBack};

    enum STRATEGY{ExtremeWeight, MultipleWeight, IncreaseWeight};

    SampleModel()
    {
        A    = 0;
        Acpy = 0;
        min  = cols;
        srows = rows = cols = 0;
        matrix_norm = 0.0;
        extreme_weight = 10e-5;
        sample_size_proc = 95;
        sample_size = 0;

        sample.assign(rows,true);

        strat=ExtremeWeight;

        lasterror = NoError;
    }

    SampleModel(const GNU_gama::SparseMatrix<Float,Index>* Ap,
                const_iterator lbegin, const_iterator lend)
    {
        A     = Ap;
        Acpy  = A->replicate();

        l_begin = lbegin;
        l_end   = lend;

        rows  = A->rows();
        cols  = A->columns();
        srows = min = cols;

        lasterror = NoError;
        matrix_norm = matrixNorm();

        if (matrix_norm)  extreme_weight = 1/matrix_norm;
        else              extreme_weight = 10e-5;

        sample_size_proc = 95;
        sample_size =  Index(rows*(Float(sample_size_proc)/100));

        set_l(lbegin,lend);

        sample.assign(rows,true);
        strat=ExtremeWeight;
    }

    ~SampleModel()
    {
    }

    int init(const GNU_gama::SparseMatrix<Float,Index>* Ap,
             const_iterator lbegin, const_iterator lend)
    {
        A     = Ap;
        Acpy  = A->replicate();

        rows  = A->rows();
        cols  = A->columns();

        srows = min  = cols;

        matrix_norm = matrixNorm();

        if (matrix_norm)  extreme_weight = 1/matrix_norm;
        else              extreme_weight = 10e-5;

        sample_size_proc = 95;
        sample_size =  Index(rows*(Float(sample_size_proc)/100));

        set_l(lbegin,lend);

        sample.assign(rows,true);
        strat=ExtremeWeight;

        return lasterror;
    }

    int clear()
    {
        Acpy = A->replicate();

        lcpy.clear();
        drows.clear();

        set_l(l_begin,l_end);

        sample.clear();
        sample.assign(rows,true);

        return lasterror;
    }

    // Funkce vytvori kopiu vektoru pravej strany

    int set_l(const_iterator bbegin, const_iterator bend)
    {
        if((bend-bbegin) != rows) { return BadDimension; }

        for (Index i=Index(); i<rows;i++)
        {
            lcpy.push_back(*(bbegin+i));
        }

        return lasterror;
    }

    Float matrixNorm()
    {
        Logger::comment(Logger::DEBUG3,"Class: SampleModel, Function: matrixNorm");

        Float norm = Float();

        for(Index i=Index(); i<rows; i++)
        {
            Float* e = Acpy->end(i+1);

            for(Float* n = Acpy->begin(i+1); n!=e; n++)
            {
                norm += *n * *n;
            }
        }

        return sqrt(norm);
    }

    // Funkce vygeneruje mnozinu n percent merani

    int randomSample()
    {
        Logger::comment(Logger::DEBUG3,"Class: SampleModel, Function: randomSample");

        clear();

        Index item;
        Index size = rows - sample_size;

        //cout << "SIZE " << sample_size << endl;

        // Nahodne vybereme riadky ktore sa do vzorku nepouziju

        while(drows.size() < size)
        {
            item = Index(Float(rows)*rand()/(RAND_MAX+1.0));

            if(sample[item])
            {
                sample[item]=false;
                drows.push_back(item);
                //cout << item << endl;
            }
        }

        sample.clear();                // vektor som pouzil aby som zabranil viacnasobnemu vyberu
        sample.assign(rows,true);      // rovnakeho riadku teraz ho uz nepotrebujem

        reWeightRows(drows);

        drows.clear();

        return lasterror;
    }

    // Funkce prenasobi riedku maticu a vektor pravej strany extremnymi vahami

    int reWeightRows(vector<Index> &drows)
    {
        if(! lcpy.size()) { return BadDimension; }

        Float div = extreme_weight/drows.size();
        Float inc = div;

        if (strat == ExtremeWeight)
        {
            Logger::comment(Logger::DEBUG2,"Vybrana strategia extremnych vah");

            // nasobim len radky ktore chci zrusit
            for(IT k=drows.begin(); k!=drows.end(); k++)
            {
                // nenasobim riadky ktore uz boli extremnou vahou nasobene

                if(sample[*k])
                {
                    sample[*k]=false;

                    Float* e = Acpy->end(*k+1);

                    for(Float* n = Acpy->begin(*k+1); n!=e; n++)
                    {
                        *n *= extreme_weight;
                    }

                    lcpy[*k] *= extreme_weight;
                }
            }
        }
        else if (strat == IncreaseWeight)
        {
            Logger::comment(Logger::DEBUG2,"Vybrana strategia linearneho rozdelenia vah v zavislosti na velkosti rezidua");

            for(IT k=drows.begin(); k!=drows.end(); k++)
            {
                Float* e = Acpy->end(*k+1);

                for(Float* n = Acpy->begin(*k+1); n!=e; n++)
                {
                    *n *= inc;
                }

                lcpy[*k] *= inc;
                inc += div;
            }
        }
        else if (strat == MultipleWeight)
        {
            Logger::comment(Logger::DEBUG2,"Vybrana strategia nasobynch vah");

            for(IT k=drows.begin(); k!=drows.end(); k++)
            {
                Float* e = Acpy->end(*k+1);

                for(Float* n = Acpy->begin(*k+1); n!=e; n++)
                {
                    *n *= extreme_weight;
                }

                lcpy[*k] *= extreme_weight;
            }
        }
        else
        {
            Logger::comment(Logger::ERROR,"Unknown strategy");
            return -1;
        }

        return lasterror;
    }

    void  set_weight(Float w){ extreme_weight = w;}
    void  set_strategy(STRATEGY s){ strat = s;}

    void  set_proc(Float pro)
    {
        sample_size_proc = pro;
        sample_size =  Index(rows*(Float(sample_size_proc)/100));
    }

    Float get_weight()       { return extreme_weight;  }
    Float get_proc()         { return sample_size_proc;}
    Index get_sample_size()  { return sample_size;     }
    Index get_rows()         { return rows;            }
    Index get_cols()         { return cols;            }
    Index get_min()          { return min;             }
    Index get_srows()        { return srows;           }

    vector<bool>& get_sample_vec(){ return sample; }
    Index get_sample_vec_size(){ return sample.size(); }

    GNU_gama::SparseMatrix<>* get_A(){ return Acpy; }
    typename vector<Float>::const_iterator get_lbeg() { return lcpy.begin(); }
    typename vector<Float>::const_iterator get_lend() { return lcpy.end(); }

    void print_l()
    {
        for(Index i=Index(); i<lcpy.size();i++)
        {
            cout << lcpy[i] << endl;
        }
        cout <<"\n";
    }

private:

    const GNU_gama::SparseMatrix<Float,Index>* A;      // matica planu
          GNU_gama::SparseMatrix<Float,Index>* Acpy;   // kopie matice planu

    typedef typename vector<Index>::const_iterator CIT;    // konstanty iterator pres vektor
    typedef typename vector<Index>::iterator IT;           // iterator pres vektor


    Index  rows;                  // pocet riadkov matice planu
    Index  cols;                  // pocet stlpcov matice planu
    Index  min;                   // minimalny pocet merani pro vzorek (pocet merani bez nadbytocnych merani)

    vector<Index>  drows;         // cisla riadkov na prevahovanie
    vector<bool> sample;          // mnozina (ne)pouzitych riadkov

    vector<Float> lcpy;           // kopie prave strany

    const_iterator l_begin;
    const_iterator l_end;

    Float sample_size_proc;       // percentualny pocet pouzitych riadkov
    Index sample_size;            // pocet pouzitch riadkov

    Index srows;

    Float extreme_weight;
    Float matrix_norm;

    STRATEGY strat;

    int lasterror;
};

#endif // SAMPLE_MODEL_H

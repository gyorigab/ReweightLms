#ifndef GENREAL_CMP_H
#define GENREAL_CMP_H

#include <smatrix.h>
#include <matvec.h>
#include <cmath>
#include <vector>
#include <sample_model.h>
#include <conjungate_gradients_nr.h>
#include <median.h>
#include <logger.h>
#include <settings.h>
#include <db_results.h>

using namespace std;

template<typename Float=double, typename Index=std::size_t,
         typename const_iterator=GNU_gama::Vec<>::const_iterator,
         typename       iterator=GNU_gama::Vec<>::iterator>
class GeneralCmp
{
public:

    enum { noError, badMatrix, badDimension };

    GeneralCmp(GNU_gama::SparseMatrix<Float,Index>* Ap,
               iterator       l_bg, iterator       l_ed,
               const_iterator P_bg, const_iterator P_ed)
    {
        A = Ap;

        cols = A->columns();
        rows = A->rows();

        P_begin = P_bg;
        P_end   = P_ed;
        l_begin = l_bg;
        l_end   = l_ed;

        maxiter = 10;
        minmedian = 10e5;

        lasterror=noError;
    }

    ~GeneralCmp()
    {
    }

    int homogenize()
    {
        if(rows != (P_end-P_begin) ||
           rows != (l_end-l_begin)   )
        {
            return badDimension;
        }

        if(A->rows() < A->columns())
        {
            return badMatrix;
        }

        for(Index k=1;k<=rows;k++)
        {
            Float* e = A->end(k);

            for(Float* n = A->begin(k); n!=e; n++)
            {
                *n *= sqrt(*(P_begin + k - 1));
            }

            *(l_begin + k - 1) *= sqrt(*(P_begin + k - 1));
        }

        return lasterror;
    }

    int cmpSqrResiduals(const_iterator xbegin, const_iterator xend,
                        iterator       rbegin, iterator       rend)
    {
        if(cols != (xend-xbegin) ||
           rows != (rend-rbegin)   )
        {
            return badDimension;
        }

        for (Index j=Index(); j<rows; j++)
          {
            Float* n = A->begin(j+1);
            Float* e = A->end  (j+1);

            Float t = *(l_begin+j);

            for(Index* i = A->ibegin(j+1); n!=e; n++, i++)
              {
                t -= *n * *(xbegin + *i - 1);
              }

            *(rbegin + j) = t*t;
          }
        return lasterror;
    }    

    int solveLms(iterator rbeg, iterator rend,
                 iterator xbeg, iterator xend)
    {
        srandom(time(0));

        ConjungateGradientsNr<Float,Index> cjg(A);

        SampleModel<Float,Index,const_iterator,iterator> smpl(A,l_begin,l_end);

        GNU_gama::SparseMatrix<Float,Index>* Asmpl;

        typename vector<Float>::const_iterator lbeg;
        typename vector<Float>::const_iterator lend;

        vector<Index> drows;
        Median<Float, Index> med;

        Index minrows  = smpl.get_min();
        Index usedrows ;

        Float med_bfr =  0.0;
        Float med_aft =  0.0;
        Float max_residual = 0.0;
        Float min_residual = 0.0;
        Float abs_residual = 0.0;

        Index iter  = 0;
        Index citer = 0;

        vector<bool> worse_drows;               // pri ktorej bol spocitany najmensi median
        vector<Float> best_adjusted_residuals;  // najlepsie vyrovnane kvadraty rezidui
        vector<Float> best_adjusted_unkonwns;   // najlepsie vyrovnane nezname

        Float summedian=0;

        smpl.set_proc(set->get_sample_proc());
        med.set_drows_size(set->get_remove_rows());

        if(set->get_method() == "ExtremeWeight")       {smpl.set_strategy(SampleModel<>::ExtremeWeight);}
        else if (set->get_method() == "MultipleWeight"){smpl.set_strategy(SampleModel<>::MultipleWeight);}
        else if (set->get_method() == "IncreaseWeight"){smpl.set_strategy(SampleModel<>::IncreaseWeight);}

        while(maxiter)
        {
            usedrows = smpl.get_sample_size();

            Logger::comment(Logger::DEBUG2,"Generujem nahodny vzor");

            smpl.randomSample();                        // generujem nahondy vzorek

            Asmpl = smpl.get_A();                       // ukazatel na nahodny vzorek matice

            lbeg = smpl.get_lbeg();                     // iterator na zaciatok nahodne vzorku redukovanych merani
            lend = smpl.get_lend();

            cjg.reset_A( Asmpl );                       // nastavim vzorek matice gradientom

            Logger::comment(Logger::DEBUG2,"Pocitam vektor neznamych pomocou CJG");

            cjg.solve(lbeg,lend,xbeg,xend);             // spocitam vektor neznamych pomocou gradientov

            Logger::comment(Logger::DEBUG2,"Pocitam kavadraty rezidui");

            cmpSqrResiduals(xbeg,xend,rbeg,rend);       // spocitam kvadraty rezidui

            drows.clear();

            Logger::comment(Logger::DEBUG2,"Pocitam median");

            med.cmpMedian(rbeg,rend,drows);             // spocitam median a vratim vektor riadkov na odstranenie

            med_aft = med.get_median();

            summedian += med_aft;                       // pro vypocet prumerneho medianu
            iter++;
            citer++;

            do
            {
                med_bfr = med_aft;                      // median n+1 je mensi ako median n tak ho nastavim

                Logger::comment(Logger::DEBUG2,"Prevahujem maticu planu a vektor redukovanych merani");

                smpl.reWeightRows(drows);               // odstranim pozadovane riadky

                Asmpl = smpl.get_A();                   // zase vratim ukazatel na vzorek s odstranenymi riadkami
                lbeg  = smpl.get_lbeg();                // obdobne iterator na vektor redukovanych merani
                lend  = smpl.get_lend();

                // cout << Asmpl << endl;
                // for (int i=0; i < 63; i++) {cout <<"L2: " << *(lbeg+i) << endl;}

                cjg.reset_A( Asmpl );                   // nastavim novy vzorek matice gradientom

                Logger::comment(Logger::DEBUG2,"Pocitam vektor neznamych pomocou CJG");

                cjg.solve(lbeg,lend,xbeg,xend);         // spocitam vektor neznamych pomocou gradientov

                Logger::comment(Logger::DEBUG2,"Pocitam kavadraty rezidui");

                cmpSqrResiduals(xbeg,xend,rbeg,rend);   // spocitam kvadraty rezidui

                //for (int i=0; i < 22; i++) {cout <<"X: " << *(xbeg+i) << endl;}
                //for (int i=0; i < 63; i++) {cout <<"R: " << *(rbeg+i) << endl;}

                drows.clear();

                med.cmpMedian(rbeg,rend,drows);                  // spocitam median a vratim vektor riadkov na odstranenie

                //for (int i=0; i < drows.size(); i++){ cout <<"riadky na mazanie: " << drows[i] << endl; }

                Logger::comment(Logger::DEBUG2,"Pocitam median");

                med_aft=med.get_median();
                max_residual=med.get_max_residual();
                min_residual=med.get_min_residual();
                abs_residual=med.cmpAbsResidual(rbeg,rend);

                summedian += med_aft;                            // pro vypocet prumerneho medianu
                iter++;

                usedrows = smpl.get_sample_vec_size();

                //cout << "Median pred a po: " << med_bfr << " " << med_aft << endl;
                db->insertMedianCmp(citer,iter,med_bfr,med_aft,max_residual,min_residual,abs_residual);
            }
            while( med_bfr > med_aft && usedrows > minrows);

            //cout << "===========" << endl;

            // vysledky s najnizsiim medianom

            if(med_bfr < minmedian)
            {
                best_adjusted_residuals.clear();
                best_estimation(best_adjusted_residuals,rbeg,rend);

                best_adjusted_unkonwns.clear();
                best_estimation(best_adjusted_unkonwns,xbeg,xend);

                worse_drows.clear();
                worse_drows = smpl.get_sample_vec();

                minmedian = med_bfr;
            }

            maxiter-- ;
        }

        //cout << "Mnozina upravenych riadkov " << endl;

        for(Index i=0; i<worse_drows.size();  i++)
        {
            if(!worse_drows[i])
            {
                //cout << i << " ";
                db->insertReweightedRows(i);
            }
        }

        //cout << "\n Najlepsi odhad vyrovnanych rezidui: \n";

        Index id=1;

        for(typename vector<Float>::iterator i=best_adjusted_residuals.begin();
            i!=best_adjusted_residuals.end();i++)
        {
            //cout << *i << endl;
            db->insertEstimatedResiduals(id,*i);
            id++;
        }

        id=1;
        //cout << "\n Najlepsi odhad vyrovnanych neznamych: \n";

        for(typename vector<Float>::iterator i=best_adjusted_unkonwns.begin();
            i!=best_adjusted_unkonwns.end();i++)
        {
            //cout << *i << endl;
            db->insertEstimatedUnknowns(id,*i);
            id++;
        }

        //cout << endl;

        // najmensi median je median_bfr
        //cout << "Najmensi median: " << minmedian << endl;
        //cout << "Prumerny median: " << summedian/iter << endl;

        db->insertMedian(minmedian,summedian/iter);

        return lasterror;
    }

    Index best_estimation(vector<Float>& vec,const_iterator beg, const_iterator end)
    {
        for(;beg!=end;beg++)
        {
            vec.push_back(*beg);
        }
        return Index();
    }

    GNU_gama::SparseMatrix<>* get_A(){ return A; }

    void print(iterator l_bg, iterator l_ed)
    {
        for(Index i=Index(); i<rows; i++)
        {
            cout << *(l_bg + i) << endl;
        }
    }

    void set_setting(LmsSettings<Float,Index>* lset)
    {
        set=lset;
        maxiter = set->get_maxiter();
    }

    void set_database(DbResults* ldb){db=ldb;}

    void set_maxiter( Index fmaxiter)  { maxiter   = fmaxiter;}

    void set_log(Logger* flog){delete log; log=flog;}

    Index get_rows() { return rows; }
    Index get_cols() { return cols; }

private:

    GNU_gama::SparseMatrix<Float,Index>* A;

    LmsSettings<Float,Index> *set;
    DbResults* db;

    const_iterator P_begin;
    const_iterator P_end;

    iterator l_begin;
    iterator l_end;

    Index rows;
    Index cols;

    Index maxiter;
    Float minmedian;

    int lasterror;
};


#endif // GENREAL_CMP_H

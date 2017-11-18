#ifndef MEDIAN_H
#define MEDIAN_H

#include <matvec/matvec.h>
#include <algorithm>
#include <vector>
#include <set>
#include <sample_model.h>
#include <cmath>

using namespace std;

template<typename Float=double,typename Index=std::size_t>
class Median
{
public:

    enum { noError, badMatrix, badDimension };

    struct Pair
    {
        Float  val;
        Index index;
    };

    struct byVal
    {
        bool operator()(const Pair &left,const Pair &right)
        {
            return left.val < right.val;
        }
    };

    Median()
    {
        size            = 0;
        drows_size_proc = 10;
        drows_size      = Index(size*(Float(drows_size_proc)/100));
        median          = 0.0;
        lasterror       = noError;
    }

    // Funkcia spocita median a naplni mnozinu indexami odlahlych merani spocitanych z rezidui

    template<typename iterator>
    int cmpMedian(iterator rbeg, iterator rend, vector<Index> &drows)
    {
        std::vector<Pair> v;

        size = rend-rbeg;
        Index mid = size/2;

        Pair p;

        for(Index i=Index(); i<size; i++)
        {
            p.val   = *rbeg++;
            p.index = i;

            v.push_back(p);
        }

        // tiredim tak aby som vedel ako sa zmenili pozicie indexov

        sort(v.begin(),v.end(),byVal());

        //len pre vypis zmeny stvorcov reziudi
        //for(int i=0; i<size;i++)
        //{
        //    cout << v[i].val << endl;
        //}

        max_residual=v[size-1].val;
        min_residual=v[0].val;

        if(size%2 == 0) { median = (v[mid-1].val + v[mid].val)/2; }
        else            { median = v[mid].val; }

        for(Index i=size-1; i>(size-drows_size-1); i--)
        {
            drows.push_back(v[i].index);
        }

        return lasterror;
    }

    // skontrolovat... sucet stvorcov reziudi
    // v databaze stlpec abs_residual

    template<typename iterator>
    Float cmpAbsResidual(iterator rbeg, iterator rend)
    {
        size = rend-rbeg;
        Float p=0;

        for(Index i=Index(); i<size; i++)
        {
            p += *rbeg++;
        }

        return p;
    }

    Float get_median(){ return median; }
    Float get_max_residual(){ return max_residual;}
    Float get_min_residual(){ return min_residual;}

    void set_drows_size(Index drs){ drows_size = drs;}
    void set_drows_size_proc(Index drsp)
    {
        drows_size_proc = drsp;
        drows_size = Index(size*(Float(drows_size_proc)/100));
    }

private:

    Index size;             // velkost sustavy
    Float median;           // median
    Float max_residual;     // maximalny stvorec rezidua
    Float min_residual;     // minimalny stvorec rezidua
    Float abs_residual;     // absolutna hodnota rezidui

    Index drows_size_proc;
    Index drows_size;       // pocet odstranenych riadkov

    int lasterror;                 
};

#endif // MEDIAN_H

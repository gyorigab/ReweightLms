#ifndef CONJUGATE_GRADIENTS_NR_H
#define CONJUGATE_GRADIENTS_NR_H

#include <smatrix.h>
#include <logger.h>
#include <thread>
#include <atomic>
#include <mutex>

static std::mutex barrier;

template<typename Float=double, typename Index=std::size_t>
class ConjungateGradientsNr
{
public:

  enum { noError, badMatrix, badDimension };

  ConjungateGradientsNr(const GNU_gama::SparseMatrix<Float,Index>* Ap)
  {
    A = Ap;

    eps       = Float(1e-12);
    maxiter   = A->rows();

    rows      = A->rows();
    cols      = A->columns();
    itercount = Index();

    lasterror = noError;
    if (A->rows() < A->columns() || A->rows() <= Float())
    {
       lasterror = badMatrix;
    }

    threads_count = 2;
  }

  int reset_A(const GNU_gama::SparseMatrix<Float,Index>* Ap )
  {
      if (rows != Ap->rows() || cols != Ap->columns())
      {
        return badMatrix;
      }
      else { A = Ap; return lasterror; }
  }

  ~ConjungateGradientsNr()
  {
  }

  void setMax(int max) { maxiter = max;  }
  int  getMax() const  { return maxiter; }
  void setEps(int tol) { eps = tol;  }
  int  getEps() const  { return eps; }
  int  getIterations() const { return itercount; }
  void setThreadsCount(int threads_cnt) { threads_count = threads_cnt; }
  int getThreadsCount() { return threads_count; }

  // conjugate gradients for solution x = inv(A)*b

  template<typename const_iterator, typename iterator>
  int solve(const_iterator bbegin, const_iterator bend,
            iterator       xbegin, iterator       xend)
    throw()
  {
      Logger::comment(Logger::DEBUG3,"Class: ConjugateGradients, Function: solve");

    if (rows != Index(bend-bbegin) || cols != Index(xend-xbegin) ||
        rows <= Index())
      {
        lasterror = badDimension;
      }
    if (lasterror != noError) return lasterror;

    std::vector<std::thread> threads;
    std::vector<Index> bounds_cols = bound(threads_count, cols);
    std::vector<Index> bounds_rows = bound(threads_count, rows);


    Float* p  = new Float[cols];
    Float* r  = new Float[rows];
    Float* ap = new Float[rows];
    Float* s  = new Float[cols];

    Float alfa, r0, beta, temp;

    for(Index i=Index();i<cols;i++) s[i]=0;

    // initialization for residuals r = b - A*x
    // vector x is not initialized (not set to zero vector!)

    r0 = Float();
    for (Index j=Index(); j<rows; j++)
      {
        Float* n = A->begin(j+1);
        Float* e = A->end  (j+1);

        Float t = *(bbegin + j);
        for(Index* i = A->ibegin(j+1) ; n!=e; n++, i++)
          {
            t -= *n * *(xbegin + *i - 1);
          }

        r[j] = t;
        r0  += t*t;
      }

    // initialization s = A'r, p=s

    for(Index i=1; i<=rows; i++ )
    {
        Float* n = A->begin(i);
        Float* e = A->end  (i);

        for(Index* j = A->ibegin(i) ; n!=e; n++, j++)
        {
            s[*j-1] += *n * r[i-1];
        }
    }

    for(Index i=Index();i<cols;i++) p[i]=s[i];

    // main loop

    itercount = Index();
    for(Index i=Index(); i<maxiter; i++, itercount++)
      {
        //single-thread dot pruduct function (testing purposes)
        //Float temp1 = dot_product(s,s,cols);

        std::atomic<Float> dp(0.0);
        //Float dp = 0.0;

        for (Index i = 0; i < threads_count; i++)
        {
            threads.push_back(std::thread(dot_product_multithread_lockfree, s, s, std::ref(dp),bounds_cols[i], bounds_cols[i+1]));
            //threads.push_back(std::thread(dot_product_multithread, s, s,std::ref(dp),bounds_cols[i], bounds_cols[i+1]));
        }

        for(auto &t : threads) { t.join(); }

        threads.clear();

        temp = dp.load();

        if(temp <= r0*eps) break;

        // axp = A*p
        Float* axp = ap;
        for (Index k=1; k<=rows; k++, axp++)
          {
            Float* n = A->begin(k);
            Float* e = A->end  (k);

            *axp = Float();
            for(Index* i=A->ibegin(k) ; n!=e; n++, i++)
              {
                *axp += *n * p[*i - 1];
              }
          }

        dp.store(0.0);
        //dp = 0.0;

        for (Index i = 0; i < threads_count; i++)
        {
            threads.push_back(std::thread(dot_product_multithread_lockfree, ap, ap, std::ref(dp),bounds_rows[i], bounds_rows[i+1]));
            //threads.push_back(std::thread(dot_product_multithread, ap, ap,std::ref(dp),bounds_rows[i], bounds_rows[i+1]));
        }

        for(auto &t : threads) { t.join(); }

        threads.clear();

        // single-thread function dot pruduct (testing purposes)
        //Float alfa1 = temp/dot_product(ap,ap,rows);

        alfa = temp/dp.load();

        for(Index i=Index(); i<cols; i++)
        {
            *(xbegin+i) += alfa * p[i];
            s[i] = Float();
        }

        for(Index i=Index(); i<rows; i++)
        {
            r[i] -= alfa * ap[i];
        }

        for(Index i=1; i<=rows; i++ )
        {
            Float* n = A->begin(i);
            Float* e = A->end  (i);

            for(Index* j = A->ibegin(i) ; n!=e; n++, j++)
            {
                s[*j-1] += *n * r[i-1];
            }
        }

        dp.store(0.0);
        //dp = 0.0;

        for (Index i = 0; i < threads_count; i++)
        {
            threads.push_back(std::thread(dot_product_multithread_lockfree, s, s, std::ref(dp),bounds_cols[i], bounds_cols[i+1]));
            //threads.push_back(std::thread(dot_product_multithread, s, s,std::ref(dp),bounds_cols[i], bounds_cols[i+1]));
        }

        for(auto &t : threads) { t.join(); }

        threads.clear();

        // single-thread function dot pruduct (testing purposes)
        //Float beta1 = dot_product(s,s,cols)/temp;

        beta = dp.load()/temp;

        for(Index i=Index(); i<cols; i++)
          {
            p[i] = s[i] + beta * p[i];
          }
      }

    delete[] p;
    delete[] r;
    delete[] ap;
    delete[] s;

    return lasterror;
  }

private:

  const GNU_gama::SparseMatrix<Float,Index>* A;

  Index  rows;
  Index  cols;
  Index  itercount;
  Float  eps;
  Index  maxiter;
  int    lasterror;

  int threads_count;

  Float dot_product(const Float* a, const Float* b, Index size)
  {
    Float dp = Float();
    for (Index i=Index(); i<size; i++) dp += a[i]*b[i];
    return dp;
  }

  // multi-thread dot product lock free
  static void dot_product_multithread_lockfree(const Float *a, const Float *b, std::atomic<Float>& dp, Index left_bound, Index right_bound)
  {
      Float partial_dp = Float();

      for( Index i=left_bound; i < right_bound; i++)
      {
          partial_dp += a[i]*b[i];
      }

      auto dot_product = dp.load();
      while ( !dp.compare_exchange_weak(dot_product, dot_product + partial_dp) );
  }

  // multi-thread dot product
  static void dot_product_multithread(const Float *a, const Float *b, Float& dp, Index left_bound, Index right_bound)
  {
      Float partial_dp = Float();

      for( Index i=left_bound; i < right_bound; i++)
      {
          partial_dp += a[i]*b[i];
      }

      std::lock_guard<std::mutex> block_threads_until_finish_this_job(barrier);
      dp += partial_dp;
  }

  std::vector<Index> bound(Index threads, Index size) {

      std::vector<Index> bounds;

      Index delta = size / threads;
      Index reminder = size % threads;
      Index n1 = 0, n2 = 0;

      bounds.push_back(n1);

      for (Index i = Index(); i < threads; ++i)
      {
          n2 = n1 + delta;

          if (i == threads - 1) { n2 += reminder; }

          bounds.push_back(n2);
          n1 = n2;
      }

      return bounds;
  }


};

#endif // CONJUGATE_GRADIENTS_NR_H

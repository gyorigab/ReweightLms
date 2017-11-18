#ifndef CONJUNGATE_GRADIENTS_H
#define CONJUNGATE_GRADIENTS_H

#include <smatrix.h>

template<typename Float=double, typename Index=std::size_t>
class ConjungateGradients
{
public:

  enum { noError, badMatrix, badDimension };

  ConjungateGradients(const GNU_gama::SparseMatrix<Float,Index>* Ap)
  {
    A = Ap;

    eps       = Float(1e-12);
    maxiter   = A->rows();

    size      = A->rows();
    itercount = Index();

    lasterror = noError;
    if (A->rows() != A->columns() || A->rows() <= Float())
      {
        lasterror = badMatrix;
      }
  }

  ~ConjungateGradients()
  {
  }

  void setMax(int max) { maxiter = max;  }
  int  getMax() const  { return maxiter; }
  void setEps(int tol) { eps = tol;  }
  int  getEps() const  { return eps; }
  int  getIterations() const { return itercount; }

  // conjugate gradients for solution x = inv(A)*b

  template<typename const_iterator, typename iterator>
  int solve(const_iterator bbegin, const_iterator bend,
            iterator       xbegin, iterator       xend)
    throw()
  {
    if (size != Index(bend-bbegin) || size != Index(xend-xbegin) ||
        size <= Index())
      {
        lasterror = badDimension;
      }
    if (lasterror != noError) return lasterror;

    Float* p  = new Float[size];
    Float* r  = new Float[size];
    Float* ap = new Float[size];

    Float alfa, r0, beta, temp;

    // initialization for residuals r = b - A*x
    // vector x is not initialized (not set to zero vector!)

    r0 = Float();
    for (Index j=Index(), k=1; k<=size; k++, j++)
      {
        Float* n = A->begin(k);
        Float* e = A->end  (k);

        Float t = *(bbegin + j);
        for(Index* i=A->ibegin(k) ; n!=e; n++, i++)
          {
            t -= *n * *(xbegin + *i - 1);
          }

        r[j] = t;
        p[j] = t;
        r0  += t*t;
      }

    // main loop

    itercount = Index();
    for(Index i=Index(); i<maxiter; i++, itercount++)
      {
        temp = dot_product(r,r);

        if(temp <= r0*eps) break;

        // axp = A*p
        Float* axp = ap;
        for (Index k=1; k<=size; k++, axp++)
          {
            Float* n = A->begin(k);
            Float* e = A->end  (k);

            *axp = Float();
            for(Index* i=A->ibegin(k) ; n!=e; n++, i++)
              {
                *axp += *n * p[*i - 1];
              }
          }

        alfa = temp/dot_product(p,ap);

        for(Index i=Index(); i<size; i++)
          {
            *(xbegin+i) += alfa * p[i];
            r[i] -= alfa * ap[i];
          }

        beta = dot_product(r,r)/temp;

        for(Index i=Index(); i<size; i++)
          {
            p[i] = r[i] + beta * p[i];
          }
      }

    delete[] p;
    delete[] r;
    delete[] ap;

    return lasterror;
  }


private:

  const GNU_gama::SparseMatrix<Float,Index>* A;

  Index  size;
  Index  itercount;
  Float  eps;
  Index  maxiter;
  int    lasterror;

  Float dot_product(const Float* a, const Float* b)
  {
    Float dp = Float();
    for (Index i=Index(); i<size; i++) dp += a[i]*b[i];
    return dp;
  }

};

#endif // CONJUNGATE_GRADIENTS_H

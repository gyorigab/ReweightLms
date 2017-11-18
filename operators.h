#ifndef OPERATORS_H
#define OPERATORS_H

#include <smatrix.h>
#include <string>
#include <fstream>
#include <iostream>
#include <ios>

using namespace GNU_gama;
using namespace std;


int check(char *subor,int &m,int &n, int &nonzero);
int check_sparse(const std::string subor,int &m,int &n, int &nonzero);

std::istream& read_sparse(istream &data, SparseMatrix<> *A);

std::istream& operator>>(std::istream &is, SparseMatrix<>  *A);
std::ostream& operator<<(std::ostream &os, SparseMatrix<>  *A);

template<typename Float>
std::istream& read_vec(std::istream &is, Float* v);

template<typename Float>
std::ostream& write_vec(std::ostream &os, Float* v);


#endif // OPERATORS_H

#ifndef _UTILS_HPP
#define _UTILS_HPP

#include <fstream>
#include <iostream>
#include <complex>

#include "fermions_wd.hpp"

using namespace std;

typedef complex<double>       t_complex;

void saveEigenvaluesToBinary(const EigenResult& result, const int LS, const string& fname);
void saveEigenvectorsToBinary(const EigenResult& result, const int LS, const string& fname);
void printEigenvaluesAndSave(const EigenResult& result, const int LS, const string& eigval_fname, const string& eigvec_fname);

#endif // _UTILS_HPP
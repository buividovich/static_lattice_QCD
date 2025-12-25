#ifndef _LINALG_HPP_
#define _LINALG_HPP_

#include <cmath>
#include <functional>

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <complex>
#include <numbers>
#include <algorithm>
#include <lapacke.h>
#include <cblas.h>
#include <omp.h>
//#pragma omp declare reduction(ComplexAdd : t_complex : omp_out += omp_in) initializer(omp_priv = t_complex(0.0, 0.0))

using namespace std;
typedef complex<double>       t_complex;

constexpr double pi = numbers::pi;

void init_rng(int rng_engine_seed=0);

void rand_complex(t_complex* z, int n);
t_complex rand_complex();
void rand_complex_z4(t_complex* z, int n);

inline t_complex InnerProduct(const t_complex* phi, const t_complex* psi, const int n)
{
    //TODO: use BLAS functions here
    double sum_re = 0.0, sum_im = 0.0;
    #pragma omp parallel for reduction(+:sum_re,sum_im)
    for(int i=0; i<n; i++)
    {
        sum_re += phi[i].real()*psi[i].real() + phi[i].imag()*psi[i].imag();
        sum_im += phi[i].real()*psi[i].imag() - phi[i].imag()*psi[i].real();
    }
    return t_complex(sum_re, sum_im);
}

double   norm(t_complex* psi, uint n);
double   norm_diff(t_complex* psi1, t_complex* psi2, uint n);
t_complex trace(const t_complex* M, const int n);

void split_array(const t_complex* arr, int n, vector<t_complex>& pos, vector<t_complex>& neg, vector<t_complex>& zero);
bool match_arrays(
    vector<t_complex> zero1, vector<t_complex> zero2,
    vector<t_complex> pos1, vector<t_complex> pos2,
    vector<t_complex> neg1, vector<t_complex> neg2);


int diagonalize(t_complex* A, t_complex* evals, t_complex* revecs, t_complex* levecs, int n);
void conjugate_transpose(t_complex* A, int n);

//This implementation of the matrix exponential assumes A is Hermitian
int matrix_exponential(t_complex* res, const t_complex* A, const t_complex t = t_complex(0.0, 1.0), const int n=3);

t_complex element_product(const t_complex* A, const int n);
double element_log_sum(const t_complex* A, const int n);
void matrix_multiplication(t_complex* res, const t_complex* A, const t_complex* B, const int n);
void matrix_vector_mult(t_complex* res, const t_complex* A, const t_complex* v, const int n);

void    aA_plus_bB(t_complex* out, const t_complex a, const t_complex* A, const t_complex b, const t_complex* B, const uint n);
void   A_pluseq_bB(t_complex*   A, const t_complex b, const t_complex *B,                                        const uint n);
void	   rescale(t_complex*   A, const t_complex a,                                                            const uint n);

void Eigenstates(t_complex* eigenvalues, const t_complex* M, const int n);
void Eigenstates(t_complex* eigenvalues, t_complex* eigenvectors, const t_complex* M, const int n);
void EigenstatesHermitian(double* eigenvalues, t_complex* eigenvectors, const t_complex* M, const int n);

#endif
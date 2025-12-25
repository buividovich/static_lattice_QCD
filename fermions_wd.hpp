#ifndef _FERMIONS_WD_HPP
#define _FERMIONS_WD_HPP

#include <iostream>
#include <omp.h>

#include "arpackpp/arrgcomp.h"
#include "gauge_field.hpp"
#include "color_algebra.hpp"
#include "spinor_algebra.hpp"
#include "linalg.hpp"

using namespace std;

inline int get_4x4_spin(int spin_a, int spin_b) {
    return spin_a*NSPINOR + spin_b;
}

inline int get_color_spin(int color, int spinor) {
    return color*NSPINOR + spinor;
}

inline int get_index(int x, int color, int spinor) {
    return x*NCOLOR*NSPINOR + color*NSPINOR + spinor;
}

inline int get_tau_index(int tau, int x, int color, int spinor, int vol3D) {
    return tau*vol3D*NCOLOR*NSPINOR + x*NCOLOR*NSPINOR + color*NSPINOR + spinor;
}

//TODO: can we avoid the use of std::vector here? Plain pointers should be better for performance
struct EigenResult {
    vector<t_complex> eigenvalues;
    vector<t_complex> eigenvectors;
    int nconv;
    int NS3D;
    bool converged;
};


class FermionsWD : public GaugeField {
private:
public:
    //Summary of all parameters
    double xi = 1.0;               // bare gluonic antisotropy
    double gamma_f = 1.0;          // bare fermionic anisotropy
    double nu_ratio = 1.0;         // ratio of gluonic and fermionic anisotrophies
    double u_t = 1.0;              // temporal tadpole
    double u_s = 1.0;              // spatial tadpole
    double cT = 1.0;               // anisotropic temporal Clover parameter
    double cR = 1.0;               // anisotropic spatial Clover parameter
    double a_tau = 1.0;            // temporal lattice spacing
    double a_s = 1.0;              // spatial lattice spacing   
    double spacing_ratio = 1.0;    // a_tau / a_s
    double factors[4] = {1.0, 1.0, 1.0, 1.0};   // {1/u_t, 1/(gamma_f * u_t), 1/(gamma_f * u_t), 1/(gamma_f * u_t)}
    double mass = 0.0;             // fermion mass in lattice units
    double mass_term = 4.0;        // mass + 1/u_t + 3/(gamma_f * u_t)
    double mass_term_3D = 3.0;     // mass + 3/(gamma_f * u_t)
    int    NS = 0;                 // Total number of spinor components = vol * 3 * 4
    int    NS3D = 0;               // Total number of spinor components in 3D = vol3D * 3 * 4
    po::options_description wd_fermion_params{"Parameters of the Wilson-Dirac Fermions"};
    void parse_command_line_args(int argc, char* argv[]);
    void print_wd_fermion_params() const;
    void init_calculable_params();

    //Constructors
    FermionsWD(string fname, int argc, char* argv[]);
    FermionsWD(int argc, char* argv[]);

    //Methods that work with the fermion fields
    void CloverTerm(t_complex* res, const t_complex* psi, const int x);
    void CloverTerm(t_complex* res, const t_complex* psi, const int x, const int tau);
    void HamiltonianWD(t_complex* res, const t_complex* psi, const int tau = 0);
    void HamiltonianWDClover(t_complex* res, const t_complex* psi, const int tau = 0);
    void OperatorWD(t_complex* res, const t_complex* psi);
    void OperatorWDHamiltonian(t_complex* res, const t_complex* psi);
    void SpatialTerm(t_complex* res, const t_complex* psi, const int tau = 0);
    void Alpha(t_complex* res, const t_complex* psi, const int tau = 0);
    void Beta(t_complex* res, const t_complex* psi, const int tau = 0);
    void AlphaInvBeta(t_complex* res, const t_complex* psi);  // Placeholder
    void IdentityOperator3D(t_complex* res, const t_complex* psi);
    void ProjectionOperator(t_complex* res, const t_complex* psi);

    // Eigensolvers
    int  ExtremalEigenstates(const int nev, const int ncv, const double tol, double* evals, t_complex* evecs, int mode=0); //Find either the lowest or the highest eigenvalues of WD_clover
    EigenResult StdLowestEigenstates(const int nev, const int ncv, const double tol);
    void GenLowestEigenstates(int nev, int ncv, double tol);
    
    // Helper Functions
    void OperatorToMatrix(t_complex* res, void (FermionsWD::*O)(t_complex*, const t_complex*), const int n);
    void OperatorToMatrix(t_complex* res, void (FermionsWD::*O)(t_complex*, const t_complex*, const int), const int n, const int tau);
    void OperatorChaintoMatrix(t_complex* res, void (FermionsWD::*Ops[])(t_complex*, const t_complex*), int n, int num_ops);
    void MatrixInversion(t_complex* M, int n);
    bool CheckHermiticity(void (FermionsWD::*O)(t_complex*, const t_complex*, const int));
    bool CheckGamma5Hermiticity(void (FermionsWD::*O)(t_complex*, const t_complex*));

    // Momentum functions
    void HamiltonianWDMomentum(t_complex* res, int* n);
    void MomentumEigenstates(t_complex* eigenvalues, t_complex* mom_hamiltonian);
    void MomentumEigenstates(double* eigenvalues, t_complex* eigenvectors, t_complex* mom_hamiltonian);
    void EuclideanCorrelatorMomentum(t_complex& res, t_complex* op, const int tau);
    void RealTimeCorrelatorMomentum(t_complex& res, t_complex* op, const double t, const double* eigenvalues, const t_complex* eigenvectors);
};

#endif // _FERMIONS_WD_HPP
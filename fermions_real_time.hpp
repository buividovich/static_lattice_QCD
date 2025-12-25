#ifndef _FERMIONS_REAL_TIME_HPP
#define _FERMIONS_REAL_TIME_HPP

#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fermions_wd.hpp"
#include "minmax.hpp"

double FermiFactorFunction(double E, void* beta); // Fermi factor function 1/(exp(beta*E) + 1)
double EuclideanEvolutionFunction(double E, void* beta_tau); // Euclidean evolution function exp(-tau*E)/(1 + exp(-beta*E)), where beta_tau = {beta, tau}

class FermionsRealTime : public FermionsWD 
{
    public:
        //Parameters specific to real-time evolution
        int    step_taylor_expansion  = 12; // Number of terms in the Taylor expansion of the time-evolution operator
        int    static_projection_mode = 1; // Mode of static projection: 0 - no projection, 1 - project to diagonal form
        int    verbosity = 0;              // Verbosity level for real-time evolution output
        //Stout smearing parameters
        int    stout_nsteps = 2;        // Number of Stout smearing steps
        double stout_rho_s  = 0.14;     // Spatial Stout smearing parameter
        double stout_rho_t  = 0.0;      // Temporal Stout smearing parameter
        // Derived parameters
        double Emin = 0.0;               // Minimum eigenvalue of the Hamiltonian (for MinMax polynomial approximation)
        double Emax = 0.0;               // Maximum eigenvalue of the Hamiltonian
        // To avoid re-calculating the eigenspectrum each time, let's save the eigenvalues and eigenvectors in files in a separate folder
        bool write_eigenmodes = false;  // Whether to write the computed eigenvalues and eigenvectors to files
        bool read_eigenmodes  = false;  // Whether to read the eigenvalues and eigenvectors from files
        std::string eigen_folder = "./eigenmodes/"; // Folder to store eigenvalues and eigenvectors

        MinMaxPolynomial*  FermiFactor               = nullptr; // MinMax polynomial approximation of the Fermi factor function (real part, default)
        MinMaxPolynomial*  FermiFactorIm             = nullptr; // MinMax polynomial approximation of the Fermi factor function (imaginary part, in case of nontrivial global Polyakov loop)
        MinMaxPolynomial** EuclideanEvolutionFactors = nullptr; // MinMax polynomial approximation of the Euclidean evolution function exp(-tau*E)/(1 + exp(-beta*E)), where beta_tau = {beta, tau}

        //Temporary arrays for real-time evolution
        t_complex* tmp = nullptr;
        t_complex* tmp1 = nullptr;
        t_complex* tmp2 = nullptr;

        po::options_description rt_fermion_params{"Parameters of the Real-Time Evolution"};
        void parse_command_line_args(int argc, char* argv[]);
        void print_rt_fermion_params() const;
        void allocate_temporary_arrays();
        void init_hamiltonian();
        void init_euclidean_evolution(int argc, char* argv[]);
        void init_minkowski_evolution(int argc, char* argv[]);
        //Constructors
        FermionsRealTime(string fname, int argc, char* argv[]) : FermionsWD(fname, argc, argv)
        {
            parse_command_line_args(argc, argv);
            print_wd_fermion_params();  // Temporarily moved: should be called in main function
            print_rt_fermion_params();
            allocate_temporary_arrays();
        }
        FermionsRealTime(int argc, char* argv[]) : FermionsWD(argc, argv)
        {
            parse_command_line_args(argc, argv);
            print_wd_fermion_params();  // Temporarily moved: should be called in main function
            print_rt_fermion_params();
            allocate_temporary_arrays();
        }
        ~FermionsRealTime();
        //State preparation
        void ApplyMinMaxPolynomial(t_complex* res, const t_complex* psi, MinMaxPolynomial* P);
        void test_eigensystem(double* evals, t_complex* evecs, int nev);
        //Operators of physical observables
        void g0gi(t_complex* res, const t_complex* psi, int i); // gamma_0 * gamma_i
        //Methods for real-time evolution
        void TimeEvolutionStep(t_complex* res, const t_complex* psi, const double dt);
        int StepTaylorPrecision(int &step_taylor_expansion, const double dt, const double Emax, const double tmax, const double tol);
};

#endif // _FERMIONS_REAL_TIME_HPP
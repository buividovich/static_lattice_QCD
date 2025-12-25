#include <iostream>
#include <fstream>
#include <stdexcept>
#include <omp.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <vector>
#include <random>
using namespace std; 
namespace po = boost::program_options;

#include "fermions_wd.hpp"
#include "ansi_io.hpp"
#include "linalg.hpp"
#include "utils.hpp"
#include "fermions_real_time.hpp"

string ConfigDir = "./test_configs/";
string OutputDir = "./data/";
int    num_configs = 5;                 // Number of config files to use
int    num_stochastic_estimators = 10;  // Number of stochastic estimators to use
double dt = 0.9;                        // Real-time step size
double tmax = 1000.0;                   // Maximum evolution time (for real-time evolution
bool   euclidean_correlator = false;    // Whether to compute the Euclidean correlator
bool   real_time_correlator = false;    // Whether to compute the Real-Time correlator
bool   free_fermions        = false;    // Whether to run in free fermion mode (no gauge fields)
bool   momentum_space       = false;    // Whether to run in momentum space mode (no gauge fields)
bool   get_max_evals        = false;    // Whether to only compute min/max eigenvalues for config files

void EuclideanCorrelator(FermionsRealTime* WD)
{
    cout << ansi::green << "Starting calculation of the Euclidean correlator ..." << ansi::reset << endl;
    string OutputFileName = OutputDir + "/" + WD->ConfigFileName + ".GE";
    ofstream fout(OutputFileName, ios::binary | ios::app);
    if (!fout) { cerr << ansi::red << "Error: Cannot open output file " << OutputFileName << ansi::reset << endl; exit(1); }

    cout << ansi::green << "Output will be written to file " << ansi::cyan << OutputFileName << ansi::reset << endl;

    int LT2 = WD->LT/2 + 1; //Number of independent time separations due to time-reversal symmetry
    double* GE = new double[LT2]; // Euclidean correlator for each stochastic estimator

    t_complex* z4     = new t_complex[WD->NS3D]; //Stochastic source with Z(4) noise
    t_complex* rleft  = new t_complex[WD->NS3D]; //Resulting propagator
    t_complex* rright = new t_complex[WD->NS3D]; //Resulting propagator
    t_complex* tmp    = new t_complex[WD->NS3D]; //Temporary array

    cout << ansi::green << "Using " << ansi::yellow << num_stochastic_estimators << ansi::green << " stochastic estimators ..." << ansi::reset << endl << flush;

    for(int iest=0; iest<num_stochastic_estimators; iest++) //Number of stochastic estimators
    {
        cout << ansi::green << "Stochastic estimator " << ansi::yellow << (iest + 1) << "/" << num_stochastic_estimators << " " << ansi::magenta << flush;
        rand_complex_z4(z4, WD->NS3D);
        
        std::fill(GE, GE + LT2, 0.0);
        
        for(int it=0; it<LT2; it++)
        {
            // |rright> = exp(-(beta - tau)*H)/(1 + exp(-beta*H)) * |z4>
            WD->ApplyMinMaxPolynomial(rright, z4, WD->EuclideanEvolutionFactors[WD->LT - it]);
            for(int i=1; i<=3; i++)
            {
                WD->g0gi(tmp, z4, i); // |tmp> = gamma_0 * gamma_i |z4>
                WD->ApplyMinMaxPolynomial(rleft, tmp, WD->EuclideanEvolutionFactors[it]); // |rleft>  = exp(-tau*H)/(1 + exp(-beta*H)) * gamma_0 * gamma_i * |z4>
                
                WD->g0gi(tmp, rright, i); // gamma_0 * gamma_i
               
                GE[it] += InnerProduct(rleft, tmp, WD->NS3D).real()/3.0; // Average over i = 1, 2, 3
                cout << "." << std::flush;
            }
        }
        cout << ansi::green << " done." << ansi::reset << endl;

        // Write the results to the output file
        fout.write(reinterpret_cast<char*>(GE), LT2 * sizeof(double));
        fout.flush();
    }
    // Clean up
    delete[] z4;
    delete[] rleft;
    delete[] rright;
    delete[] tmp;
    delete[] GE;

    fout.close();
}

void RealTimeCorrelator(FermionsRealTime* WD)
{
    //To be implemented
    t_complex* z4     = new t_complex[WD->NS3D]; //Stochastic source with Z(4) noise
    t_complex* rright = new t_complex[WD->NS3D]; //Right thermal state
    t_complex* tmp    = new t_complex[WD->NS3D]; //Temporary

    t_complex* rleft[3]  = {new t_complex[WD->NS3D], new t_complex[WD->NS3D], new t_complex[WD->NS3D]}; //Left thermal state - one for each spatial direction

    int nsteps = (int)std::ceil(tmax/dt);
    std::cout << ansi::green << "Starting calculation of the Real-Time correlator, nsteps =" << ansi::yellow << nsteps << ansi::reset << endl;

    t_complex* GR = new t_complex[nsteps]; // Real-time correlator for each stochastic estimator

    string OutputFileName = OutputDir + "/" + WD->ConfigFileName + ".GR";
    ofstream fout(OutputFileName, ios::binary | ios::app);
    if (!fout) { cerr << ansi::red << "Error: Cannot open output file " << OutputFileName << ansi::reset << endl; exit(1); }
    std::cout << ansi::green << "Output will be written to file " << ansi::cyan << OutputFileName << ansi::reset << endl;

    for(int iest=0; iest<num_stochastic_estimators; iest++) //Number of stochastic estimators
    {
        std::cout << ansi::green << "Stochastic estimator " << ansi::yellow << (iest + 1) << "/" << num_stochastic_estimators << " " << ansi::magenta << flush;
        rand_complex_z4(z4, WD->NS3D);
        std::fill(GR, GR + nsteps, t_complex(0.0, 0.0));

        WD->ApplyMinMaxPolynomial(tmp, z4, WD->FermiFactor); // Right thermal state |tmp> = FermiFactor * |z4>
        std::copy(z4, z4 + WD->NS3D, rright); 
        A_pluseq_bB(rright, t_complex(-1.0, 0.0), tmp, WD->NS3D); // |rright> = |z4> - FermiFactor * |z4> = 1/(1 + exp(-beta*H)) * |z4>
        for(int i=0; i<3; i++)
        {
            WD->g0gi(tmp, z4, i+1); // gamma_0 * gamma_i * |z4>
            WD->ApplyMinMaxPolynomial(rleft[i], tmp, WD->FermiFactor); // Left thermal state |rleft> = FermiFactor * gamma_0 * gamma_i * |z4>
        }

        for(int istep=0; istep<nsteps; istep++)
        {
            double t = istep * dt;

            for(int i=0; i<3; i++)
            {
                //First, contract to get the correlator
                WD->g0gi(tmp, rright, i+1); // gamma_0 * gamma_i * |rright>
                GR[istep] += InnerProduct(rleft[i], tmp, WD->NS3D)/3.0; // Average over i = 1, 2, 3 <rleft_i| gamma_0 * gamma_i |rright>
                WD->TimeEvolutionStep(tmp, rleft[i], dt); // Evolve |rleft> to the next time step
                std::copy(tmp, tmp + WD->NS3D, rleft[i]);
            }
            
            WD->TimeEvolutionStep(tmp, rright, dt); // Evolve |rright> to the next time step
            std::copy(tmp, tmp + WD->NS3D, rright);

            cout << "." << std::flush;
        }
        cout << ansi::green << " done." << ansi::reset << endl;
        // Write the results to the output file
        fout.write(reinterpret_cast<char*>(GR), nsteps * sizeof(t_complex));
        fout.flush();
    }
    fout.close();

    delete[] z4;
    for (int i=0; i<3; i++) delete[] rleft[i];
    delete[] rright;
    delete[] tmp;
    delete[] GR;
}

void TestRealTimeEvolution(FermionsRealTime* WD);
void test_fermi_factor_application(FermionsRealTime* WD, int sm_nev);

int main(int argc, char** argv)
{
    po::options_description main_options("Main options");
    main_options.add_options()
        ("help,h",    "Print help message")
        ("ConfigDir",                 po::value<string>(&ConfigDir                )->default_value(ConfigDir                 ), "Directory to read configurations from")
        ("OutputDir",                 po::value<string>(&OutputDir                )->default_value(OutputDir                 ), "Directory to output data")
        ("num_configs",               po::value<int>(   &num_configs              )->default_value(num_configs               ), "Number of config files to use")
        ("num_stochastic_estimators", po::value<int>(   &num_stochastic_estimators)->default_value(num_stochastic_estimators ), "Number of stochastic estimators for trace estimation")
        ("dt",                        po::value<double>(&dt                       )->default_value(dt                        ), "Real-time step size for real-time evolution")
        ("tmax",                      po::value<double>(&tmax                     )->default_value(tmax                      ), "Maximum evolution time for real-time evolution")
        ("euclidean_correlator",      po::bool_switch(&euclidean_correlator       )->default_value(euclidean_correlator      ), "Whether to compute the Euclidean correlator")
        ("real_time_correlator",      po::bool_switch(&real_time_correlator       )->default_value(real_time_correlator      ), "Whether to compute the Real-Time correlator")
        ("free_fermions",             po::bool_switch(&free_fermions              )->default_value(free_fermions             ), "Whether to run in free fermion mode (no gauge fields)")
        ("momentum_space",            po::bool_switch(&momentum_space             )->default_value(momentum_space            ), "Whether to run in momentum space mode (no gauge fields)")
        ("get_max_evals",             po::bool_switch(&get_max_evals              )->default_value(get_max_evals             ), "Whether to only compute min/max eigenvalues for config files");

	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(main_options).allow_unregistered().run(), vm);
	po::notify(vm);

    // Check if help is needed and exit
    if (vm.count("help")) { cout << main_options << endl; return 0; };

    if (momentum_space)
    {
        FermionsWD* WD = new FermionsWD(argc, argv);
        WD->print_wd_fermion_params();

        if (euclidean_correlator) {
            cout << ansi::green << "Starting calculation of the Euclidean correlator ..." << ansi::reset << endl;
            string OutputFileName = OutputDir +
                "/momentum_exact_m" + std::to_string(WD->mass) +
                "_LT" + std::to_string(WD->LT) +
                "_LS" + std::to_string(WD->LS) +
                "_dt" + std::to_string(dt) +
                "_tmax" + std::to_string(tmax) +
                ".GE";
            ofstream fout(OutputFileName, ios::binary | ios::trunc);
            if (!fout) { cerr << ansi::red << "Error: Cannot open output file " << OutputFileName << ansi::reset << endl; exit(1); }

            cout << ansi::green << "Output will be written to file " << ansi::cyan << OutputFileName << ansi::reset << endl;

            int LT2 = (WD->LT/2) + 1;

            cout << ansi::green << "Using exact diagonalization in momentum space ..." << ansi::reset << endl << flush;

            t_complex GE[LT2];
            fill(GE, GE + LT2, t_complex(0.0, 0.0));

            for (int i = 1; i <= 3; i++) {
                t_complex g0gi[NSPINOR*NSPINOR];
                matrix_multiplication(g0gi, gamma_matrix(0), gamma_matrix(i), NSPINOR);

                for (int tau = 0; tau < LT2; tau++) {
                    t_complex correlator = t_complex(0.0, 0.0);
                    WD->EuclideanCorrelatorMomentum(correlator, g0gi, tau);
                    cout << "." << std::flush;
                    GE[tau] += correlator;
                }
            }

            cout << ansi::green << " done." << ansi::reset << endl;
            // Write the results to the output file
            fout.write(reinterpret_cast<char*>(GE), LT2 * sizeof(t_complex));
            fout.flush();
            fout.close();
        }

        if (real_time_correlator) {
            cout << ansi::green << "Starting calculation of the real-time correlator ..." << ansi::reset << endl;
            string OutputFileName = OutputDir +
                "/momentum_exact_m" + std::to_string(WD->mass) +
                "_LT" + std::to_string(WD->LT) +
                "_LS" + std::to_string(WD->LS) +
                "_dt" + std::to_string(dt) +
                "_tmax" + std::to_string(tmax) +
                ".GR";
            ofstream fout(OutputFileName, ios::binary | ios::trunc);
            if (!fout) { cerr << ansi::red << "Error: Cannot open output file " << OutputFileName << ansi::reset << endl; exit(1); }

            cout << ansi::green << "Output will be written to file " << ansi::cyan << OutputFileName << ansi::reset << endl;
            cout << ansi::green << "Using exact diagonalization in momentum space ..." << ansi::reset << endl << flush;

            int nsteps = std::ceil(tmax/dt);

            t_complex GR[nsteps];
            fill(GR, GR + nsteps, t_complex(0.0, 0.0));

            double* eigenvalues = new double[WD->vol3D*NSPINOR];
            t_complex* eigenvectors = new t_complex[WD->vol3D*NSPINOR*NSPINOR];
            t_complex mom_hamiltonian[NSPINOR*NSPINOR];

            WD->MomentumEigenstates(eigenvalues, eigenvectors, mom_hamiltonian);

            for (int i = 1; i <= 3; i++) {
                t_complex g0gi[NSPINOR*NSPINOR];
                matrix_multiplication(g0gi, gamma_matrix(0), gamma_matrix(i), NSPINOR);

                for (int step = 0; step < nsteps; step++) {
                    double t = step*dt;
                    t_complex correlator = t_complex(0.0, 0.0);
                    WD->RealTimeCorrelatorMomentum(correlator, g0gi, t, eigenvalues, eigenvectors);
                    cout << "." << std::flush;
                    GR[step] += correlator;
                }
            }

            cout << ansi::green << " done." << ansi::reset << endl;
            // Write the results to the output file
            fout.write(reinterpret_cast<char*>(GR), nsteps * sizeof(t_complex));
            fout.flush();
            fout.close();

            delete[] eigenvalues;
            delete[] eigenvectors;
        }

        delete WD;
        return 0;
    }
    if(free_fermions)
    {
        cout << ansi::cyan << "Running in free fermion mode ..." << ansi::reset << endl;
        FermionsRealTime* WD = new FermionsRealTime(argc, argv);
        WD->ConfigFileName = string("free_fermions_m") + std::to_string(WD->mass) + string("_LT") + std::to_string(WD->LT) + string("_LS") + std::to_string(WD->LS);
        if(euclidean_correlator)
        {
            WD->init_euclidean_evolution(argc, argv);
            EuclideanCorrelator(WD);
        }
        if(real_time_correlator)
        {
            WD->init_minkowski_evolution(argc, argv);
            RealTimeCorrelator(WD);
        }
        delete WD;
        return 0;
    };

    // requires #include <random>
    std::vector<std::filesystem::directory_entry> entries;
    for (const auto &e : std::filesystem::directory_iterator(ConfigDir))
        if (std::filesystem::is_regular_file(e.path()))
            entries.push_back(e);
    if (entries.empty()) { cerr << ansi::red << "No config files in " << ConfigDir << ansi::reset << endl; return 1; }

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, entries.size() - 1);

    for (int i = 0; i < num_configs; i++)
    {
        const auto &entry = entries[dist(rng)];
        FermionsRealTime* WD = new FermionsRealTime(entry.path(), argc, argv);
        auto filename = entry.path().filename().string();
        WD->print_wd_fermion_params(); std::cout << std::flush;

        cout << ansi::cyan << "Running in config mode..." << ansi::reset << endl;


        if(get_max_evals)
        {
            WD->init_minkowski_evolution(argc, argv);
        }
        if(euclidean_correlator)
        {
            WD->ConfigFileName = filename +
                            string("_m") + std::to_string(WD->mass) +
                            string("_LT") + std::to_string(WD->LT) +
                            string("_LS") + std::to_string(WD->LS);
            WD->init_euclidean_evolution(argc, argv);
            EuclideanCorrelator(WD);
        }
        if(real_time_correlator)
        {
            WD->ConfigFileName = filename +
                            string("_m") + std::to_string(WD->mass) +
                            string("_LT") + std::to_string(WD->LT) +
                            string("_LS") + std::to_string(WD->LS) +
                            string("_dt") + std::to_string(dt) +
                            string("_tmax") + std::to_string(tmax);
            WD->init_minkowski_evolution(argc, argv);
            RealTimeCorrelator(WD);
        }
        delete WD;
    }

    return 0;
}

//This piece of code was used to generate free gauge configurations and save them to files
// GaugeField* GF = new GaugeField(argc, argv);
// string output_file_name = "./free_configs/free_" + std::to_string(GF->LT) + "x" + std::to_string(GF->LS) + "n0";
// GF->WriteToFile(output_file_name);
// delete GF;

//cout << ansi::green << "Constructor finished ..." << ansi::reset << endl;
//WD->init_euclidean_evolution(argc, argv);
//std::cout << ansi::green << "Euclidean evolution initialized ..." << ansi::reset << std::endl;
//EuclideanCorrelator(WD);

// Testing functions used for development and debugging

void test_fermi_factor_application(FermionsRealTime* WD, int sm_nev)
{
    //Testing the Fermi factor application
    t_complex* psi = new t_complex[WD->NS3D];
    t_complex* res = new t_complex[WD->NS3D];
    rand_complex_z4(psi, WD->NS3D);
    double norm_psi = norm(psi, WD->NS3D);
    cout << ansi::green << "Norm of the input vector = " << ansi::yellow << norm_psi << " " << norm_psi*norm_psi << ansi::reset << endl;
    WD->ApplyMinMaxPolynomial(res, psi, WD->FermiFactor);
    double norm_res = norm(res, WD->NS3D);
    cout << ansi::green << "Norm of the output vector = " << ansi::yellow << norm_res << " " << ansi::reset << endl;
    //Now testing with eigenvectors
    int nev = 20;
    t_complex* evecs = new t_complex[nev*WD->NS3D];
    double*    evals = new double[nev];
    nev = WD->ExtremalEigenstates(nev, 7*nev, 1e-9, evals, evecs, 0);
    WD->test_eigensystem(evals, evecs, nev);

    for(int ie=0; ie<nev; ie++)
    {
        WD->ApplyMinMaxPolynomial(res, evecs + ie*WD->NS3D, WD->FermiFactor);
        double norm_evec = norm(evecs + ie*WD->NS3D, WD->NS3D);
        double ff = FermiFactorFunction(evals[ie], &(WD->beta));
        A_pluseq_bB(res, t_complex(-ff, 0.0), evecs + ie*WD->NS3D, WD->NS3D);
        std::cout << ansi::green << "Eigenvalue[" << ie << "] = " << ansi::yellow << evals[ie] << ansi::green << ", Fermi factor = " << ansi::yellow << ff << ansi::green << ", error after applying Fermi factor = " << ansi::yellow << norm(res, WD->NS3D) / norm_evec << ansi::reset << std::endl;    
    };

    delete[] psi;
    delete[] res;
    delete[] evals;
    delete[] evecs; 
}

void TestRealTimeEvolution(FermionsRealTime* WD)
{
    //Let's first test the real time evolution of a random state - whether the norm is preserved
    t_complex* z4 = new t_complex[WD->NS3D]; //Stochastic source with Z(4) noise
    t_complex* res = new t_complex[WD->NS3D]; //Resulting state after time evolution
    t_complex* tmp = new t_complex[WD->NS3D]; //Temporary array
    rand_complex_z4(z4, WD->NS3D); double norm_z4 = norm(z4, WD->NS3D);
    std::copy(z4, z4 + WD->NS3D, res);

    int nsteps = (int)std::ceil(tmax/dt);

    //Forward evolution
    for(int n=0; n<nsteps; n++)
    {
        WD->TimeEvolutionStep(tmp, res, dt);
        std::copy(tmp, tmp + WD->NS3D, res);

        double norm_res = norm(res, WD->NS3D);
        cout << ansi::green << "\t step " << ansi::yellow << n << ansi::green << ansi::green << ",\t norm(res)/norm(z4) = " << ansi::yellow << norm_res/norm_z4 << ansi::reset << endl;
    }

    //Forward evolution
    for(int n=0; n<nsteps; n++)
    {
        WD->TimeEvolutionStep(tmp, res, -dt);
        std::copy(tmp, tmp + WD->NS3D, res);

        double norm_res = norm(res, WD->NS3D);
        cout << ansi::green << "\t step " << ansi::yellow << n << ansi::green << ansi::green << ",\t norm(res)/norm(z4) = " << ansi::yellow << norm_res/norm_z4 << ansi::reset << endl;
    }

    cout << ansi::green << "Difference between initial and final state after forward and backward evolution = " << ansi::yellow << norm_diff(res, z4, WD->NS3D ) / norm_z4 << ansi::reset << endl;

    delete[] z4;
    delete[] res;
    delete[] tmp;
}

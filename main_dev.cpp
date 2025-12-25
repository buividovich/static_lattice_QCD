#include <iostream>
#include <fstream>
#include <stdexcept>
#include <omp.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <vector>
using namespace std; 
namespace po = boost::program_options;

#include "fermions_wd.hpp"
#include "ansi_io.hpp"
#include "linalg.hpp"
#include "utils.hpp"
    
int main(int argc, char** argv)
{
    cout << endl;
    string ConfigDir = "";
    string OutputDir = "";

    int nev = 8;
    int ncv = 32;
    double tol = 1e-8;  // Smallest eigenvalue gap of order 1e-5

    po::options_description main_options("Main options");
    main_options.add_options()
        ("help,h",    "Print help message")
        ("momentum",  "Momentum space eigenvalues")
        ("benchmark", "Benchmarking on Barkla for a range of nev/ncv values")
        ("ConfigDir", po::value<string>(&ConfigDir)->default_value("./test_configs/"),        "Directory to read configurations from")
        ("OutputDir", po::value<string>(&OutputDir)->default_value("./"             ),        "Directory to output data")
        ("nev",       po::value<int>(   &nev      )->default_value(                8),        "Number of eigenvalues to compute")
        ("ncv",       po::value<int>(   &ncv      )->default_value(               32),        "Number of Arnoldi vectors")
        ("tol",       po::value<double>(&tol      )->default_value(             1e-8),        "Tolerance of ARPACK eigensolver accuracy");
    
	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(main_options).allow_unregistered().run(), vm);
	po::notify(vm);

    // Check if help is needed and exit
    if (vm.count("help"))
    {
        cout << main_options << endl;
        return 0;
    }
    if (vm.count("momentum"))
    {
        FermionsWD* WD = new FermionsWD(argc, argv);

        t_complex mom_hamiltonian[NSPINOR*NSPINOR] = {};
        int nev = WD->NS3D;

        EigenResult result;
        result.eigenvalues.resize(nev);
        result.nconv = nev;
        string eigval_fname = "eigenvalues_momentum_" + to_string(WD->LS) + "x" + to_string(WD->LS) + "_" + to_string(nev) + ".bin";
        WD->MomentumEigenstates(result.eigenvalues.data(), mom_hamiltonian);

        saveEigenvaluesToBinary(result, WD->LS, eigval_fname);

        delete WD;      
        return 0;
    }

    if (vm.count("benchmark"))
    {
        init_rng();
        FermionsWD* WD = new FermionsWD(argc, argv);

        cout << ansi::green << "Running benchmark on random links (" << WD->LS << "^3):" << ansi::reset << endl;

        int nev_bench_arr[5] = {250, 500, 1000, 2000, 5000};
        int ncv_scale_arr[3] = {2, 3, 4};

        for (int nev_bench : nev_bench_arr) {
            for (int ncv_scale : ncv_scale_arr) {
                int ncv_bench = nev_bench*ncv_scale;
                double start;
                double end;

                string fname = "benchmark_nev_" + to_string(nev_bench) + "_ncv_" + to_string(ncv_bench) + ".bin";
                string fpath = OutputDir + fname;

                WD->RandomizeLinks();

                cout << "LS: " << WD->LS << endl;
                cout << "nev: " << nev_bench << ", ncv: " << ncv_bench << endl;

                start = omp_get_wtime();
                EigenResult result = WD->StdLowestEigenstates(nev_bench, ncv_bench, tol);
                end = omp_get_wtime();
                cout << "ARPACK solver time (seconds):" << (end - start) << endl;

                start = omp_get_wtime();
                saveEigenvaluesToBinary(result, WD->LS, fpath);
                saveEigenvectorsToBinary(result, WD->LS, fpath);
                end = omp_get_wtime();
                cout << "Write to file time (seconds):" << (end - start) << endl;

                size_t bytes = (result.eigenvalues.size() + result.eigenvectors.size()) * sizeof(t_complex);
                double gigabytes = static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
                cout << "Memory usage (GB): " << gigabytes << endl;
                cout << endl;
            }
        }

        delete WD;
        return 0;
    }

    if (vm.count("ConfigDir"))
    {
        cout << ansi::green << "ConfigDir = " << ansi::cyan << ConfigDir << ansi::reset << endl;

        init_rng();

        for (const auto & entry : filesystem::directory_iterator(ConfigDir))
        {
            FermionsWD* WD = new FermionsWD(entry.path(), argc, argv);

            string eigval_fname = "eigenvalues_" + entry.path().filename().string() + "_" + to_string(nev) + ".bin";
            string eigvec_fname = "eigenvectors_" + entry.path().filename().string() + "_" + to_string(nev) + ".bin";
            int ncv = min(nev*2 + 1, WD->NS3D - 1);

            EigenResult result = WD->StdLowestEigenstates(nev, ncv, tol);
            printEigenvaluesAndSave(result, WD->LS, eigval_fname, eigvec_fname);

            delete WD;
        }
    }
    else
    {
        cout << ansi::green << "Gauge links set to unity " << ansi::reset << endl;

        init_rng();

        FermionsWD* WD = new FermionsWD(argc, argv);

        string eigval_fname = "eigenvalues_unity_" + to_string(WD->LS) + "x" + to_string(WD->LS) + "_" + to_string(nev) + ".bin";
        string eigvec_fname = "eigenvectors_unity_"  + to_string(WD->LS) + "x" + to_string(WD->LS) + "_" + to_string(nev) + ".bin";
        int ncv = min(nev*2 + 1, WD->NS3D - 1);

        EigenResult result = WD->StdLowestEigenstates(nev, ncv, tol);
        printEigenvaluesAndSave(result, WD->LS, eigval_fname, eigvec_fname);

        delete WD;
    }

    return 0;
}
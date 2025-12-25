#include <iostream>
#include <fstream>
#include <stdexcept>
#include <omp.h>
#include <boost/program_options.hpp>
#include <filesystem>
using namespace std; 
namespace po = boost::program_options;

#include "ansi_io.hpp"
#include "gauge_field.hpp"
#include "linalg.hpp"
#include "utils.hpp"

int main(int argc, char** argv)
{
    string InputDir = "./test_configs/";
    string ConfigPrefix = "config_";
    string OutputDir = "./static_configs/";
    int StaticProjectionMode = 0;
    int FirstConfigIndex = 0;
    int LastConfigIndex = -1; // Process all files by default

    po::options_description main_options("Main options");
    main_options.add_options()
        ("help,h",    "Print help message")
        ("InputDir",             po::value<string>( &InputDir                   )->default_value(InputDir            ), "Directory to read input files from")
        ("ConfigPrefix",         po::value<string>( &ConfigPrefix               )->default_value(ConfigPrefix        ), "Prefix for input files")
        ("FirstConfigIndex",     po::value<int>(    &FirstConfigIndex           )->default_value(FirstConfigIndex    ), "Index of the first configuration to process")
        ("LastConfigIndex",      po::value<int>(    &LastConfigIndex            )->default_value(LastConfigIndex     ), "Index of the last configuration to process")
        ("OutputDir",            po::value<string>( &OutputDir                  )->default_value(OutputDir           ), "Directory to write output files to")
        ("StaticProjectionMode", po::value<int>(    &StaticProjectionMode       )->default_value(StaticProjectionMode), "Static projection mode (0 - set temporal links to center sector, copy spatial links from t=0; 1 - average spatial links over time and project to SU(3))");

	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(main_options).allow_unregistered().run(), vm);
	po::notify(vm);
    // Check if help is needed and exit
    if (vm.count("help")) { cout << main_options << endl; return 0; }

    //Echo the parameters before running
    cout << ansi::yellow << "This code performs static projection (mode " << StaticProjectionMode << ") of gauge configurations." << ansi::reset << endl;
    cout << ansi::green << "\tInputDir  = " << ansi::cyan << InputDir  << ansi::reset << endl;
    cout << ansi::green << "\tOutputDir = " << ansi::cyan << OutputDir << ansi::reset << endl;

    for(int iconf=FirstConfigIndex; iconf<=LastConfigIndex || LastConfigIndex<0; iconf++)
    {
        string in_fname = InputDir + "/" + ConfigPrefix + std::to_string(iconf);
        if(!filesystem::exists(in_fname)) { if(iconf==FirstConfigIndex) { cerr << ansi::red << "Error: No configuration file " << in_fname << " found!" << ansi::reset << endl; } break;}
        cout << ansi::green << "Processing configuration file " << ansi::cyan << in_fname << ansi::reset << endl;
        
        GaugeField* GF = new GaugeField(in_fname);
        cout << "Temporal and spatial plaquettes before the static projection: " << ansi::yellow << GF->MeanTemporalPlaquette() << ", " << GF->MeanSpatialPlaquette() << ansi::reset << endl;
        GF->StaticProjection(StaticProjectionMode, true);
        cout << "Temporal and spatial plaquettes after  the static projection: " << ansi::yellow << GF->MeanTemporalPlaquette() << ", " << GF->MeanSpatialPlaquette() << ansi::reset << endl;

        string out_fname = OutputDir + "/" + ConfigPrefix + std::to_string(iconf);
        cout << ansi::green << "Writing static configuration to " << ansi::cyan << out_fname << ansi::reset << endl;
        GF->WriteToFile(out_fname);
        cout << endl;
        //GF->StaticTests();
        delete GF;
    };

    return 0;
}

// // // Matrix exponential test
// std::ofstream bin_out("./data/matrix_exp_tests.bin", std::ios::binary);
// if (!bin_out) {
//     throw std::runtime_error("Failed to open 'matrix_exp_tests.bin' for binary writing");
// }
// for(int itrial=0; itrial<10; itrial++)
// {
//         rand_complex(A, 9);
//         for (int i = 0; i < 3; ++i) {
//             for (int j = i + 1; j < 3; ++j) {
//                 t_complex tmp = t_complex(0.5) * (A[i*3 + j] + std::conj(A[j*3 + i]));
//                 A[i*3 + j] = tmp;
//                 A[j*3 + i] = std::conj(tmp);
//             }
//             int d = i*3 + i;
//             A[d] = t_complex(std::real(A[d]), 0.0);
//         }
//         matrix_exponential(expA, A);

//         bin_out.write(reinterpret_cast<char*>(A),    sizeof(A));
//         bin_out.write(reinterpret_cast<char*>(expA), sizeof(expA));

//         matrix_multiplication(expA2_1, expA, expA, 3);

//         matrix_exponential(expA2_2, A, t_complex(0.0, 2.0));

//         double err = norm_diff(expA2_2, expA2_1, 9);
//         cout << "Trial " << itrial << ": ||exp(A)^2 - exp(2A)|| = " << err << endl;
//     }
//     bin_out.close();

//     return 0;
    
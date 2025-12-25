#ifndef _MINMAX_HPP
#define _MINMAX_HPP

#include <lapacke.h>
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "ansi_io.hpp"

class MinMaxPolynomial {
public:
    //Algorithm parameters
    int    max_remez_iterations  = 30;   // Maximum number of Remez iterations
    int    max_polynomial_degree = 200;  // Maximum allowed polynomial degree
    int    min_polynomial_degree = 10;   // Minimum polynomial degree to start iterations from. This might need to be lowered in case of very smooth/exactly polynomial functions
    int    degree_step           = 4;    // Step by which the polynomial degree is increased in the outer loop
    double absolute_precision    = 1e-5; // Desired precision of the approximation
    int    verbosity             = 0;    // Verbosity level: 0 - silent, 1 - final result, 2 - all iterations
    //Polynomial data
    double* cs       = nullptr;          // Polynomial coefficients
    double* xs       = nullptr;          // Critical points used for approximation
    double xmin      = -1.0, xmax = 1.0; // Approximation interval
    int degree       = 5;                // Polynomial degree
    double max_error = 1.0E+10;          // Maximum approximation error achieved
    //Rescaling coefficients to transform from [xmin, xmax] to [-1, 1]
    double a, b;                       // x = a + b*t maps [-1, 1] to [xmin, xmax]
    //Command line options interface
    po::options_description minmax_params{"Parameters of the MinMax Polynomial approximation"};
    void parse_command_line_args(int argc, char* argv[]);
    void print_minmax_params() const;
    // Basic Remez algorithm with a fixed degree
    void ConstructMinMaxPolynomial(double (*f)(double, void*), void* data, int degree);
    // Running Remez algorithm with increasing degree until absolute_precision is met
    void ConstructMinMaxPolynomial(double (*f)(double, void*), void* data);

    // Constructor: builds the minimax approximation so that max_error < absolute_precision. Uses default settings for the Remez algorithm
    MinMaxPolynomial(double (*f)(double, void*), double xmin, double xmax, void* data) : xmin(xmin), xmax(xmax)
    {
        ConstructMinMaxPolynomial(f, data);
    }

    //Constructor: reads in the parameters for the Remez algorithm from the command line arguments
    MinMaxPolynomial(int argc, char* argv[], double (*f)(double, void*), double xmin, double xmax, void* data) : xmin(xmin), xmax(xmax)
    {
        parse_command_line_args(argc, argv);
        ConstructMinMaxPolynomial(f, data);
        print_minmax_params();
    }
 
    // Destructor: cleanup allocated memory
    ~MinMaxPolynomial() {if (cs) delete[] cs; if (xs) delete[] xs;};
    
    // Evaluate polynomial at point x
    double operator()(double x) const;

    // Testing
    double TestPolynomial(double (*f)(double, void*), void* data, int num_test_points = 1000) const;
};

#endif // _MINMAX_HPP
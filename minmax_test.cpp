#include "minmax.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <fstream>

double sigmoid_function(double x, void* data) { return 1.0/(std::exp(16*x) + 1.0); }

int main(int argc, char* argv[]) 
{
    std::cout << std::fixed << std::setprecision(12);

    double xmin = -1.6;
    double xmax = 2.0;
    int    num_test_points = 10000;
    
    // Create the polynomial approximation
    MinMaxPolynomial poly(argc, argv, sigmoid_function, xmin, xmax, nullptr);

    std::cout << "Testing the polynomial approximation with " << num_test_points << " points in the range [" << xmin << ", " << xmax << "]" << std::endl;

    double max_error = poly.TestPolynomial(sigmoid_function, nullptr, num_test_points);
    std::cout << "Max error = " << max_error << std::endl;

    return 0;
}
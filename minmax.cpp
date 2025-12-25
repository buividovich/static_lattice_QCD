#include "minmax.hpp"
#include <numbers>

static constexpr double pi = std::numbers::pi;

void MinMaxPolynomial::parse_command_line_args(int argc, char* argv[])
{
    minmax_params.add_options()
        ("help,h", "Print help messages")
        ("minmax.max_remez_iterations",  po::value<int>(   &max_remez_iterations  )->default_value(max_remez_iterations ), "Maximum number of Remez iterations")
        ("minmax.max_polynomial_degree", po::value<int>(   &max_polynomial_degree )->default_value(max_polynomial_degree), "Maximum allowed polynomial degree")
        ("minmax.min_polynomial_degree", po::value<int>(   &min_polynomial_degree )->default_value(min_polynomial_degree), "Minimum polynomial degree to start iterations from. This might need to be lowered in case of very smooth/exactly polynomial functions")
        ("minmax.degree_step",           po::value<int>(   &degree_step           )->default_value(degree_step          ), "Step by which the polynomial degree is increased in the outer loop")
        ("minmax.absolute_precision",    po::value<double>(&absolute_precision    )->default_value(absolute_precision   ), "Desired precision of the approximation")
        ("minmax.verbosity",             po::value<int>(   &verbosity             )->default_value(verbosity            ), "Verbosity level: 0 - silent, 1 - final result, 2 - all iterations");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(minmax_params).allow_unregistered().run(), vm);
        po::notify(vm);
    }
    catch(const std::exception& e) { std::cerr << ansi::red << "Error parsing command line arguments for the MinMaxPolynomial class: " << e.what() << '\n' << ansi::reset; }   

    // Check if help is needed
    if(vm.count("help")) { std::cout << minmax_params << std::endl; }
}

void MinMaxPolynomial::print_minmax_params() const
{
    std::cout << ansi::cyan << "MinMax Polynomial parameters:" << ansi::reset << std::endl;
    std::cout << ansi::green << "\t max_remez_iterations         = " << ansi::yellow << max_remez_iterations     << ansi::reset << std::endl;
    std::cout << ansi::green << "\t max_polynomial_degree        = " << ansi::yellow << max_polynomial_degree    << ansi::reset << std::endl;
    std::cout << ansi::green << "\t min_polynomial_degree        = " << ansi::yellow << min_polynomial_degree    << ansi::reset << std::endl;
    std::cout << ansi::green << "\t degree_step                  = " << ansi::yellow << degree_step              << ansi::reset << std::endl;
    std::cout << ansi::green << "\t absolute_precision           = " << ansi::yellow << absolute_precision       << ansi::reset << std::endl;
    std::cout << ansi::green << "\t verbosity                    = " << ansi::yellow << verbosity                << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Argument range xmin ... xmax = " << ansi::yellow << xmin << " ... " << xmax  << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Actual polynomial degree     = " << ansi::yellow << degree                   << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Achieved max_error           = " << ansi::yellow << max_error                << ansi::reset << std::endl;
}

void MinMaxPolynomial::ConstructMinMaxPolynomial(double (*f)(double, void*), void* data, int degree)
{
    this->degree = std::max(degree, 2);
    a = 0.5*(xmax + xmin);
    b = 0.5*(xmax - xmin);

    if(xs!=nullptr) delete[] xs;
    if(cs!=nullptr) delete[] cs;
    xs = new double[degree + 2];
    cs = new double[degree + 2];

    // Initialize critical points with Chebyshev extrema
    for(int k=0; k<degree+2; k++) xs[k] = a + b*std::cos(pi*(degree+1-k)/(degree+1));

    //Allocate memory for the iterative refinement procedure ...
    double* A = new double[(degree+2)*(degree+2)];  lapack_int* ipiv = new lapack_int[degree+2];

    //Iterations of the Remez algorithm
    for(int iter=0; iter<max_remez_iterations; iter++)
    {
        if(verbosity>1)
        {
            std::cout << ansi::green << "Remez iteration " << ansi::yellow << iter << ansi::green << " xs: " << ansi::yellow;
            for (int k = 0; k < degree + 2; k++) std::cout << xs[k] << " ";
            std::cout << ansi::reset << std::endl;
        }

        //Prepare the matrix and the rhs vector and solve the linear equation for the new nodes
        for(int inode=0; inode<degree+2; inode++) 
        {
            double t = (xs[inode] - a)/b;
            // Fill matrix row: coefficients for c_0, c_1, ..., c_n, u
            //A[inode*(degree+2) + n] = T_n(t); We will perform the recursion directly on A's
            double* T = A + inode*(degree+2);
            T[0] = 1.0; T[1] = t; //Initializing the Chebyshev recursion
            for(int n=2; n<=degree; n++) T[n] = 2*t*T[n-1] - T[n-2];

            A[inode*(degree+2) + (degree+1)] = (inode%2==0? 1.0 : -1.0);

            cs[inode] = f(xs[inode], data);
        };

        // Solve using LU decomposition: A * x = b
        lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, degree+2, 1, A, degree+2, ipiv, cs, 1);
        if (info != 0) {std::cerr << ansi::red << "LAPACKE_dgesv failed with info = " << info << ansi::reset << std::endl;};

        double dx = xs[1] - xs[0];
        for(int k=1; k<degree+1; k++) dx = std::min(dx, xs[k+1]-xs[k]); //We assume the nodes are always sorted
        dx /= 100.0; //This will be the small extremum-finding step for the next part of the Remez algorithm
        int nx = int(std::ceil((xmax - xmin)/dx));
        dx = (xmax - xmin)/nx;
        max_error = 0.0;

        double hx0 = f(xmin, data) - (*this)(xmin); //Initial value of the error at the left boundary
        double hx1 = f(xmin + dx, data) - (*this)(xmin + dx); //Value of the error at the next point
        int nxs = 1; //xs[0] is already set to xmin, xs[degree+1] is already set to xmax
        for(int k=2; k<=nx; k++) //Loop over the interior points
        {
            double x = xmin + k*dx;
            double hx2 = f(x, data) - (*this)(x); //Value of the error at the next point
            if((hx1 - hx0)*(hx1 - hx2) > 0.0) //We have an extremum between x-dx and x+dx
            {
                if(nxs < degree + 1)
                {
                    xs[nxs] = xmin + (k-1)*dx; 
                    max_error = std::max(max_error, std::abs(hx1));
                }
                else if (verbosity>5) { std::cout << ansi::red << "Warning: more than degree+1 extrema found in Remez iteration " << iter << ", x = " << xmin + (k-1)*dx << ", dx = " << dx << ansi::reset << std::endl; std::exit(0); }
                nxs++;
            }
            hx0 = hx1; hx1 = hx2;
        }

        if(nxs < degree+1)
        { 
            std::cout << ansi::red << "Warning: only " << nxs << " extrema found in Remez iteration, terminating at iteration " << iter << ansi::reset << std::endl;
            std::exit(0);
        }
        if(verbosity>1) 
        {
            std::cout << ansi::green << "Found "          << ansi::yellow << nxs << " extrema,";
            std::cout << ansi::green << " max_error = "   << ansi::yellow << max_error << ",";
            std::cout << ansi::green << " u  = "          << ansi::yellow << cs[degree+1] << ",";
            std::cout << ansi::green << " max_error/u = " << ansi::yellow << max_error/cs[degree+1] << ansi::reset << std::endl;
        };

        if(std::abs(max_error/cs[degree+1] - 1.0) < 1e-5) break; //Converged
    }; //End of the Remez iteration loop
    
    delete[] A;  delete[] ipiv;
}

void MinMaxPolynomial::ConstructMinMaxPolynomial(double (*f)(double, void*), void* data)
{
    std::cout << "Initializing the minmax polynomial ..." << std::endl;

    // Start with degree 10 and increase until absolute_precision is met
    int adegree = min_polynomial_degree;
    while(max_error > absolute_precision && adegree < max_polynomial_degree)
    {
        if(verbosity > 0) std::cout << ansi::cyan << " Minmax polynomial tuning: trying degree " << ansi::magenta << adegree << ansi::cyan << "... "; 
        ConstructMinMaxPolynomial(f, data, adegree);
        if(max_error==0.0) {std::cout << ansi::red << " max_error = 0 for a polynomial of degree " << adegree << ansi::reset << std::endl; std::exit(0); }
        if(verbosity > 0) std::cout << ansi::green << "done. max_error = " << ansi::yellow << max_error << ansi::reset << std::endl;
        adegree += degree_step;
    };

    std::cout << ansi::green << "Minmax polynomial of degree " << ansi::yellow << degree << ansi::green << " constructed with max_error = " << ansi::yellow << max_error;
    if(adegree >= max_polynomial_degree) std::cout << ansi::red << " (maximum degree reached!)";
    std::cout << ansi::reset << std::endl;
}

double MinMaxPolynomial::operator()(double x) const 
{
    double t = (x - a)/b;
    double bn2 = 0.0, bn1 = 0.0, bn = 0.0;
    // Clenshaw recurrence: b_k = c_k + 2*t*b_{k+1} - b_{k+2}
    // Starting from the highest degree coefficient
    for(int k=degree; k>=1; k--)
    {
        bn = cs[k] + 2*t*bn1 - bn2;
        bn2 = bn1;
        bn1 = bn;
    }
    return cs[0] + t*bn1 - bn2;
}

double MinMaxPolynomial::TestPolynomial(double (*f)(double, void*), void* data, int num_test_points) const
{
    double max_err = 0.0;
    for(int i=0; i<num_test_points; i++)
    {
        double x = xmin + (xmax - xmin)*i/(num_test_points-1);
        double err = std::abs(f(x, data) - (*this)(x));
        max_err = std::max(max_err, err);
    }
    return max_err;
}
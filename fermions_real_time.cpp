#include "fermions_real_time.hpp"

double FermiFactorFunction(double E, void* beta)
{
    double b = *(double*)beta;
    return 1.0/(exp(b*E) + 1.0);
}

double EuclideanEvolutionFunction(double E, void* beta_tau)
{
    double* params = (double*)beta_tau;
    double  b = params[0];
    double  t = params[1];
    return exp(-t*E)/(1.0 + exp(-b*E));
}

void FermionsRealTime::parse_command_line_args(int argc, char* argv[])
{
    rt_fermion_params.add_options()
        ("FermionsRealTime.step_taylor_expansion",  po::value<int>(   &step_taylor_expansion  )->default_value(step_taylor_expansion  ), "Number of terms in the Taylor expansion of the time-evolution operator")
        ("FermionsRealTime.static_projection_mode", po::value<int>(   &static_projection_mode )->default_value(static_projection_mode ), "Mode of static projection: 0 - work with time slice 0, 1 - average links over all time slices")
        ("FermionsRealTime.write_eigenmodes",       po::bool_switch(  &write_eigenmodes       )->default_value(write_eigenmodes       ), "Whether to write the computed eigenvalues and eigenvectors to files")
        ("FermionsRealTime.read_eigenmodes",        po::bool_switch(  &read_eigenmodes        )->default_value(read_eigenmodes        ), "Whether to read the eigenvalues and eigenvectors from files")
        ("FermionsRealTime.eigen_folder",           po::value<string>(&eigen_folder           )->default_value(eigen_folder           ), "Folder to store eigenvalues and eigenvectors")
        ("FermionsRealTime.verbosity",              po::value<int>(   &verbosity              )->default_value(verbosity              ), "Verbosity level for real-time evolution output")
        ("StoutSmearing.nsteps",                    po::value<int>   (&stout_nsteps           )->default_value(stout_nsteps           ),  "Number of Stout smearing steps to apply to the gauge field")
        ("StoutSmearing.rho_s",                     po::value<double> (&stout_rho_s           )->default_value(stout_rho_s            ),  "Spatial Stout smearing parameter rho_s")
        ("StoutSmearing.rho_t",                     po::value<double> (&stout_rho_t           )->default_value(stout_rho_t            ),  "Temporal Stout smearing parameter rho_t");
        ("help,h", "Print help messages");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(rt_fermion_params).allow_unregistered().run(), vm);
        po::notify(vm);
    }
    catch(const std::exception& e)  { std::cerr << ansi::red << "Error parsing command line arguments for the FermionsRealTime class: " << e.what() << '\n' << ansi::reset; };

    // Check if help is needed
    if(vm.count("help")) { std::cout << rt_fermion_params << std::endl; }
}

void FermionsRealTime::print_rt_fermion_params() const
{
    std::cout << ansi::cyan << "Real-Time Fermion evolution parameters:" << ansi::reset << std::endl;
    std::cout << ansi::green << "\t step_taylor_expansion     = " << ansi::yellow << step_taylor_expansion     << ansi::reset << std::endl;
    std::cout << ansi::green << "\t static_projection_mode    = " << ansi::yellow << static_projection_mode    << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Saving eigenmodes?          " << ansi::yellow << (write_eigenmodes ? "Yes" : "No") << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Reading eigenmodes?         " << ansi::yellow << (read_eigenmodes  ? "Yes" : "No") << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Folder to save eigenmodes:  " << ansi::yellow << eigen_folder              << ansi::reset << std::endl;
}

void FermionsRealTime::init_hamiltonian()
{
    StaticProjection(static_projection_mode, (verbosity>0)); // Polyakov loop projected to the
    if(std::abs(global_Z3_phase - 1.0)>1.0e-6)
    {
        std::cout << ansi::red << "Warning: Nontrivial Polyakov loop after static projection: " << ansi::magenta << global_Z3_phase << ansi::reset << std::endl;
        std::cout << ansi::red << "         The imaginary part of the Fermi factor needs to be approximated as well." << ansi::reset << std::endl;
        std::cout << ansi::red << "         This functionality is not implemented yet!"                               << ansi::reset << std::endl;
        exit(1);
    };
    StoutSmearing(stout_nsteps, stout_rho_s, stout_rho_t);
    
    // Compute the extremal eigenvalues of the Hamiltonian
    double max_evals[10];
    int nev_res = 0;
    std::string evals_filename = eigen_folder + ConfigFileName + ".max.evals";
    if(read_eigenmodes)
    {
        std::ifstream evals_file(evals_filename, std::ios::binary);
        if (evals_file)
        {
            evals_file.read(reinterpret_cast<char*>(max_evals), sizeof(max_evals));
            nev_res = static_cast<int>(evals_file.gcount() / sizeof(double));
            evals_file.close();
            std::cout << ansi::green << "Read " << ansi::yellow << nev_res << ansi::green << " extremal eigenvalues of the Hamiltonian from file " << ansi::cyan << evals_filename << ansi::reset << std::endl;
        }
    }

    if (nev_res < 10)
    {
        nev_res = ExtremalEigenstates(10, 150, 1e-9, max_evals, nullptr, 1);
        if(nev_res>0 && write_eigenmodes)
        {
            std::ofstream out_evals_file(evals_filename, std::ios::binary);
            out_evals_file.write(reinterpret_cast<const char*>(max_evals), sizeof(max_evals));
            out_evals_file.close();
            std::cout << ansi::green << "Saved " << ansi::yellow << nev_res << ansi::green << " extremal eigenvalues of the Hamiltonian to file " << ansi::cyan << evals_filename << ansi::reset << std::endl;
        }
        else { std::cerr << ansi::red << "Error: Unable to compute the extremal eigenvalues of the Hamiltonian for MinMax polynomial approximation!" << ansi::reset << std::endl; exit(1); }
    }

    Emin = *std::min_element(max_evals, max_evals + nev_res);
    Emax = *std::max_element(max_evals, max_evals + nev_res);
    
    std::cout << ansi::cyan << "Real-time evolution initialization:" << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Minimum eigenvalue of the Hamiltonian Emin = " << ansi::yellow << Emin << ansi::reset << std::endl;
    std::cout << ansi::green << "\t Maximum eigenvalue of the Hamiltonian Emax = " << ansi::yellow << Emax << ansi::reset << std::endl;
    Emin = (Emin < 0.0) ? Emin*1.01 : Emin*0.99; //Make the interval a bit wider
    Emax = (Emax > 0.0) ? Emax*1.01 : Emax*0.99;
    std::cout << ansi::green << "\t Approximation range (after broadening) = " << ansi::yellow << Emin << " ... " << Emax << ansi::reset << std::endl;
}

//TODO: potential conflicts if both Minkowski and Euclidean evolutions are initialized

void FermionsRealTime::init_minkowski_evolution(int argc, char* argv[])
{
    init_hamiltonian();
    FermiFactor = new MinMaxPolynomial(argc, argv, FermiFactorFunction, Emin, Emax, &(this->beta));
}

void FermionsRealTime::init_euclidean_evolution(int argc, char* argv[])
{
    init_hamiltonian();
    if(EuclideanEvolutionFactors) return; // Already initialized
    EuclideanEvolutionFactors = new MinMaxPolynomial*[LT+1];
    double beta_tau[2] = {beta, 0.0};
    for(int it=0; it<=LT; it++)
    {
        beta_tau[1] = (double)it;
        EuclideanEvolutionFactors[it] = new MinMaxPolynomial(argc, argv, EuclideanEvolutionFunction, Emin, Emax, (void*)beta_tau);
        double max_err = EuclideanEvolutionFactors[it]->TestPolynomial(EuclideanEvolutionFunction, (void*)beta_tau);
        cout << ansi::green << "Euclidean evolution polynomial for it = " << ansi::yellow << it << ansi::green << " constructed with max_err = " << ansi::yellow << max_err << ansi::reset << endl;
    }
}

void FermionsRealTime::ApplyMinMaxPolynomial(t_complex* res, const t_complex* psi, MinMaxPolynomial* P)
{
    // Apply the MinMax polynomial approximation using Clenshaw iterations
    // Use tmp, tmp1, tmp2 as recurrence arrays
    // Use a and b from the P instance
    const double fr = 2.0/(P->b);
    const double fp = -2.0*(P->a)/(P->b);

    std::fill(tmp, tmp + NS3D, t_complex(0.0, 0.0));
    std::fill(tmp1, tmp1 + NS3D, t_complex(0.0, 0.0));
    std::fill(tmp2, tmp2 + NS3D, t_complex(0.0, 0.0));
    
    // Clenshaw recurrence
    for (int k = P->degree; k >= 1; k--) 
    {
        // tmp = (H - b) / a * tmp1 - tmp2 + coeffs[k] * psi
        HamiltonianWDClover(tmp, tmp1); 
        #pragma omp parallel for
        for (int i = 0; i < NS3D; ++i) tmp[i] = fr*tmp[i] + fp*tmp1[i]; // tmp = (2*(H-b)/a)*tmp1

        A_pluseq_bB(tmp, t_complex(-1.0, 0.0), tmp2, NS3D);
        A_pluseq_bB(tmp, t_complex(P->cs[k], 0.0), psi, NS3D);
        std::copy(tmp1, tmp1 + NS3D, tmp2);
        std::copy(tmp, tmp + NS3D, tmp1);
    }

    // Final step: res = coeffs[0] * psi + (H - b) / a * tmp1 - tmp2
    HamiltonianWDClover(res, tmp1);
    #pragma omp parallel for
    for (int i = 0; i < NS3D; ++i) res[i] = 0.5*(fr*res[i] + fp*tmp1[i]); // res = ((H-b)/a)*tmp1

    A_pluseq_bB(res, t_complex(-1.0, 0.0), tmp2, NS3D);
    A_pluseq_bB(res, t_complex(P->cs[0], 0.0), psi, NS3D);
}

void FermionsRealTime::g0gi(t_complex* res, const t_complex* psi, int i)
{
    // Apply the gamma_0 * gamma_i matrix to the spinor field psi
    // i = 1, 2, 3
    if(i<1 || i>3) { std::cerr << ansi::red << "Error: Invalid index for g0gi operator: " << i << ansi::reset << std::endl; exit(1); }

    #pragma omp parallel for
    for(int xc=0; xc<vol3D*NCOLOR; xc++)
    {
        t_complex tmp_spinor[NSPINOR];
        gamma_times_psi(      tmp_spinor, psi + xc*NSPINOR, i);
        for(int i=0; i<NSPINOR; i++) tmp_spinor[i] = t_complex(0.0, 1.0)*tmp_spinor[i]; // Multiply by gamma_0
        gamma_times_psi(res + xc*NSPINOR, tmp_spinor,       0);
    }
}

void FermionsRealTime::allocate_temporary_arrays()
{
    if(tmp)  delete[] tmp;
    if(tmp1) delete[] tmp1;
    if(tmp2) delete[] tmp2;

    tmp  = new t_complex[NS3D];
    tmp1 = new t_complex[NS3D];
    tmp2 = new t_complex[NS3D];
}

FermionsRealTime::~FermionsRealTime()
{
    if (FermiFactor) delete FermiFactor; FermiFactor = nullptr;
    if (FermiFactorIm) delete FermiFactorIm; FermiFactorIm = nullptr;

    if(EuclideanEvolutionFactors)
    {
        for(int it=0; it<=LT; it++)
            if(EuclideanEvolutionFactors[it]) delete EuclideanEvolutionFactors[it];
        delete[] EuclideanEvolutionFactors;
        EuclideanEvolutionFactors = nullptr;
    }

    if (tmp)  delete[] tmp; tmp = nullptr;
    if (tmp1) delete[] tmp1; tmp1 = nullptr;
    if (tmp2) delete[] tmp2; tmp2 = nullptr;
}

void FermionsRealTime::test_eigensystem(double* evals, t_complex* evecs, int nev)
{
	double max_evec_err = 0.0;
	std::cout << std::endl << "Checking the eigensystem... " << std::flush;
	//Eigenvector error
	for(uint ie=0; ie<nev; ie++)
	{
		HamiltonianWDClover(tmp, evecs + ie*NS3D);
		A_pluseq_bB(tmp, -evals[ie], evecs + ie*NS3D, NS3D);
		max_evec_err = std::max(max_evec_err, norm(tmp, NS3D));
	};
	std::cout << "\t Max. eigenstate error: " << max_evec_err << std::endl;
    double max_ortho_err = 0.0;
    //Orthogonality error
    for(uint ie1=0; ie1<nev; ie1++)
        for(uint ie2=ie1+1; ie2<nev; ie2++)
        {
            t_complex ip = InnerProduct(evecs + ie1*NS3D, evecs + ie2*NS3D, NS3D) - (ie1==ie2 ? t_complex(1.0, 0.0) : t_complex(0.0, 0.0));
            max_ortho_err = std::max(max_ortho_err, std::abs(ip));
        }
    
    std::cout << "\t Max. orthogonality error: " << max_ortho_err << std::endl;
}

void FermionsRealTime::TimeEvolutionStep(t_complex* res, const t_complex* psi, const double dt)
{
    // Apply one step of the time evolution operator using Taylor expansion
    // res = exp(i*H*dt) * psi
    t_complex dtc = t_complex(0.0, dt);
    t_complex coeff = dtc;

    std::copy(psi, psi + NS3D, res); // res = psi, zeros order term of the Taylor expansion
    HamiltonianWDClover(tmp, psi); // tmp = H * psi, first order term of the Taylor expansion
    rescale(tmp, coeff, NS3D); // tmp = i*dt * H * psi
    A_pluseq_bB(res, 1, tmp, NS3D); // res = (1 + i*dt*H) * psi

    for(int n=2; n<=step_taylor_expansion; n++)
    {
        std::copy(tmp, tmp + NS3D, tmp1); // tmp1 = tmp
        HamiltonianWDClover(tmp, tmp1); // tmp = H * tmp1 = H*tmp_prev
        coeff = dtc/(double)n;
        rescale(tmp, coeff, NS3D); // tmp = (i*dt)^n/n! * H^n * psi
        A_pluseq_bB(res, 1, tmp, NS3D); // res += tmp
    }
}

int FermionsRealTime::StepTaylorPrecision(int &step_taylor_expansion, const double dt, const double Emax, const double tmax, const double tol) {

    double error = 1.0;
    t_complex I = t_complex(0.0, 1.0);
    t_complex exact = std::exp(I * Emax * tmax);
    t_complex taylor_sum = t_complex(1.0, 0.0);
    t_complex term = t_complex(1.0, 0.0);

    step_taylor_expansion = 2;
    while (error > tol) {
        step_taylor_expansion++;
        term = t_complex(1.0, 0.0);
        taylor_sum = t_complex(1.0, 0.0);  // 0-th step
        for (int n = 1; n <= step_taylor_expansion; n++) {
            term *= (I * Emax * dt) / static_cast<t_complex>(n);
            taylor_sum += term;
        }
        taylor_sum = std::pow(taylor_sum, tmax/dt);
        error = std::abs(taylor_sum - exact);
        
        if (step_taylor_expansion > 30) {
            std::cout << "Warning: Taylor expansion does not converge." << std::endl;
            return 1;
        }
    }
    std::cout << "step_taylor_expansion: " << step_taylor_expansion << " | error: " << error << std::endl;
    return 0;
}
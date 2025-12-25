#include "fermions_wd.hpp"

using namespace std;

void FermionsWD::parse_command_line_args(int argc, char* argv[])
{
    wd_fermion_params.add_options()
        ("nu_ratio",             po::value<double>(&nu_ratio  )->default_value(nu_ratio), "Ratio of gluonic and fermionic anisotropies")
        ("xi",                   po::value<double>(&xi        )->default_value(xi), "Bare gluonic anisotropy")
        ("u_t",                  po::value<double>(&u_t       )->default_value(u_t), "Temporal tadpole improvement factor")
        ("u_s",                  po::value<double>(&u_s       )->default_value(u_s), "Spatial tadpole improvement factor")
        ("cT",                   po::value<double>(&cT        )->default_value(cT), "Anisotropic temporal Clover parameter")
        ("cR",                   po::value<double>(&cR        )->default_value(cR), "Anisotropic spatial Clover parameter")
        ("mass",                 po::value<double>(&mass      )->default_value(mass), "Fermion mass in lattice units");
    
    po::variables_map vm;
    try{
        //Parsing command line arguments
        po::store(po::command_line_parser(argc, argv).options(wd_fermion_params).allow_unregistered().run(), vm);
        po::notify(vm);
    }
    catch(std::exception& e) {
        std::cerr << "Error parsing command line arguments: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
   

    // Check if help is needed
    if(vm.count("help")) { std::cout << wd_fermion_params << std::endl; }

    //Parameters that are calculated in terms of input parameters
    gamma_f = xi / nu_ratio;                        // Bare fermionic antisotropy
    a_tau = u_t;                                    // lattice temporal spacing
    a_s = u_t * gamma_f;                            // lattice spatial spacing
    spacing_ratio = a_tau / a_s;
    factors[0] = 1.0/u_t;
    factors[1] = 1.0/(gamma_f*u_t);
    factors[2] = 1.0/(gamma_f*u_t);
    factors[3] = 1.0/(gamma_f*u_t);
    mass_term = mass + factors[0] + factors[1] + factors[2] + factors[3];
    mass_term_3D = mass + factors[1] + factors[2] + factors[3];
    NS = vol * NCOLOR * NSPINOR;                    // Total number of spinor components
    NS3D = vol3D * NCOLOR * NSPINOR;                // Total number of spinor components in 3D
}

void FermionsWD::print_wd_fermion_params() const
{
    std::cout << ansi::cyan << "Wilson-Dirac Fermion parameters:" << ansi::reset << std::endl;
    std::cout << ansi::green << "\t nu_ratio        = " << ansi::yellow << nu_ratio        << ansi::reset << std::endl;
    std::cout << ansi::green << "\t xi              = " << ansi::yellow << xi              << ansi::reset << std::endl;
    std::cout << ansi::green << "\t gamma_f         = " << ansi::yellow << gamma_f         << ansi::reset << std::endl;
    std::cout << ansi::green << "\t u_t             = " << ansi::yellow << u_t             << ansi::reset << std::endl;
    std::cout << ansi::green << "\t u_s             = " << ansi::yellow << u_s             << ansi::reset << std::endl;
    std::cout << ansi::green << "\t cT              = " << ansi::yellow << cT              << ansi::reset << std::endl;
    std::cout << ansi::green << "\t cR              = " << ansi::yellow << cR              << ansi::reset << std::endl;
    std::cout << ansi::green << "\t mass            = " << ansi::yellow << mass            << ansi::reset << std::endl;
    std::cout << ansi::green << "\t a_tau           = " << ansi::yellow << a_tau           << ansi::reset << std::endl;
    std::cout << ansi::green << "\t a_s             = " << ansi::yellow << a_s             << ansi::reset << std::endl;
    std::cout << ansi::green << "\t spacing_ratio   = " << ansi::yellow << spacing_ratio   << ansi::reset << std::endl;
    std::cout << ansi::green << "\t mass_term_3D    = " << ansi::yellow << mass_term_3D    << ansi::reset
              <<" (used in Hamiltonian)"<<ansi::reset<<std::endl;
    std::cout << ansi::green << "\t 4D linear space dimension = "  << ansi::yellow  << "vol * 3 * 4 = "  <<  NS   <<" (used in 4D computations)"<<ansi::reset<<std::endl;
    std::cout << ansi::green << "\t 3D linear space dimension = "  << ansi::yellow  << "vol3D * 3 * 4 = " << NS3D <<" (used in 3D computations)"<<ansi::reset<<std::endl;
}

void FermionsWD::init_calculable_params()
{
    gamma_f = xi / nu_ratio;                        // Bare fermionic antisotropy
    a_tau = u_t;                                    // lattice temporal spacing
    a_s = u_t * gamma_f;                            // lattice spatial spacing
    spacing_ratio = a_tau / a_s;
    factors[0] = 1.0/u_t;
    factors[1] = 1.0/(gamma_f*u_t);
    factors[2] = 1.0/(gamma_f*u_t);
    factors[3] = 1.0/(gamma_f*u_t);
    mass_term = mass + factors[0] + factors[1] + factors[2] + factors[3];
    mass_term_3D = mass + factors[1] + factors[2] + factors[3];
    NS = vol * NCOLOR * NSPINOR;                    // Total number of spinor components
    NS3D = vol3D * NCOLOR * NSPINOR;                // Total number of spinor components in 3D

}

FermionsWD::FermionsWD(int argc, char* argv[]) : GaugeField(argc, argv)
{ 
    parse_command_line_args(argc, argv);
    init_calculable_params();
}

FermionsWD::FermionsWD(string fname, int argc, char* argv[]) : GaugeField(fname)
{
    parse_command_line_args(argc, argv);
    init_calculable_params();
}

void FermionsWD::CloverTerm(t_complex* res, const t_complex* psi, const int x)
{
    t_complex tmp[NCOLOR*NSPINOR];
    t_complex tmp2[NCOLOR*NSPINOR];

    // Clover temporal term
    double t_term = (0.5 * cT) / (u_t*u_t * u_s*u_s);

    for (int nu = 1; nu < 4; nu++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &FieldStrength[get_fieldstrength_index(x, 0, nu, 0)], &psi[get_index(x, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            i_sigma_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], 0, nu);
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] -= t_term * tmp2[get_color_spin(color, spinor)];
            }
        }
    }

    // Clover spatial term
    double s_term = (0.5 * cR) / (gamma_f * nu_ratio * u_t * u_s*u_s*u_s);

    for (int nu = 2; nu < 4; nu++) {
        for (int mu = 1; mu < nu; mu++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &FieldStrength[get_fieldstrength_index(x, mu, nu, 0)], &psi[get_index(x, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                i_sigma_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu, nu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= s_term * tmp2[get_color_spin(color, spinor)];
                }
            }
        }
    }
}

void FermionsWD::CloverTerm(t_complex* res, const t_complex* psi, const int x, const int tau)
{
    t_complex tmp[NCOLOR*NSPINOR];
    t_complex tmp2[NCOLOR*NSPINOR];

    // Clover temporal term
    double t_term = (0.5 * cT) / (u_t*u_t * u_s*u_s);
    int x_t = x + tau*vol3D;

    for (int nu = 1; nu < 4; nu++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &FieldStrength[get_fieldstrength_index(x_t, 0, nu, 0)], &psi[get_index(x, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            i_sigma_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], 0, nu);
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] -= t_term * tmp2[get_color_spin(color, spinor)];
            }
        }
    }

    // Clover spatial term
    double s_term = (0.5 * cR) / (gamma_f * nu_ratio * u_t * u_s*u_s*u_s);

    for (int nu = 2; nu < 4; nu++) {
        for (int mu = 1; mu < nu; mu++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &FieldStrength[get_fieldstrength_index(x_t, mu, nu, 0)], &psi[get_index(x, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                i_sigma_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu, nu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= s_term * tmp2[get_color_spin(color, spinor)];
                }
            }
        }
    }
}

void FermionsWD::HamiltonianWD(t_complex* res, const t_complex* psi, const int tau)
{
    #pragma omp parallel for
    for (int x = 0; x < vol3D; x++) {
    
        // Mass term
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = mass_term_3D * psi[get_index(x, color, spinor)];
            }
        }
    
        // Hopping terms
        for (int mu = 1; mu < 4; mu++) {
            // Forwards spatial
            int nx_fwd = shift_fwd(x, mu);
            t_complex tmp[NSPINOR*NCOLOR];
            t_complex tmp2[NSPINOR*NCOLOR];

            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Link(tau, x, mu), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_minus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
            }
            A_pluseq_bB(&res[get_index(x, 0, 0)], -factors[mu], tmp2, NCOLOR*NSPINOR);

        
            // Backwards spatial
            int nx_bwd = shift_bwd(x, mu);
            t_complex Uconj[NCOLOR*NCOLOR];

            c3x3_conj(Uconj, Link(tau, nx_bwd, mu));
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Uconj, &psi[get_index(nx_bwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_plus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
            }
            A_pluseq_bB(&res[get_index(x, 0, 0)], -factors[mu], tmp2, NCOLOR*NSPINOR);
        }

        // Gamma_0
        for (int color = 0; color < NCOLOR; color++) {
            t_complex tmp[NSPINOR*NCOLOR];
            gamma_times_psi(&tmp[get_color_spin(color, 0)], &res[get_index(x, color, 0)], 0);
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = tmp[get_color_spin(color, spinor)];
            }
        }
    }
}

void FermionsWD::HamiltonianWDClover(t_complex* res, const t_complex* psi, const int tau)
{
    #pragma omp parallel for
    for (int x = 0; x < vol3D; x++) {
    
        // Mass term
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = mass_term_3D * psi[get_index(x, color, spinor)];
            }
        }
        CloverTerm(res, psi, x, tau);
    
        // Hopping terms
        for (int mu = 1; mu < 4; mu++) {
            // Forwards spatial
            int nx_fwd = shift_fwd(x, mu);
            t_complex tmp[NSPINOR*NCOLOR];
            t_complex tmp2[NSPINOR*NCOLOR];

            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Link(tau, x, mu), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_minus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
            }
            A_pluseq_bB(&res[get_index(x, 0, 0)], -factors[mu], tmp2, NCOLOR*NSPINOR);

        
            // Backwards spatial
            int nx_bwd = shift_bwd(x, mu);
            t_complex Uconj[NCOLOR*NCOLOR];

            c3x3_conj(Uconj, Link(tau, nx_bwd, mu));
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Uconj, &psi[get_index(nx_bwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_plus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
            }
            A_pluseq_bB(&res[get_index(x, 0, 0)], -factors[mu], tmp2, NCOLOR*NSPINOR);
        }

        // Gamma_0
        for (int color = 0; color < NCOLOR; color++) {
            t_complex tmp[NSPINOR*NCOLOR];
            gamma_times_psi(&tmp[get_color_spin(color, 0)], &res[get_index(x, color, 0)], 0);
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = tmp[get_color_spin(color, spinor)];
            }
        }
    }
}

void FermionsWD::OperatorWD(t_complex* res, const t_complex* psi)
{
    #pragma omp parallel for
    for (int x = 0; x < vol; x++) {
        // Mass term
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = mass_term * psi[get_index(x, color, spinor)];
            }
        }

        CloverTerm(res, psi, x);

        // Hopping terms
        for (int mu = 0; mu < 4; mu++) {
            // Forwards hopping
            int nx_fwd = shift_fwd(x, mu);
            t_complex tmp[NSPINOR*NCOLOR];
            t_complex tmp2[NSPINOR*NCOLOR];

            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Link(x, mu), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_minus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= factors[mu] * tmp2[get_color_spin(color, spinor)];
                }
            }

            // Backwards hopping
            int nx_bwd = shift_bwd(x, mu);
            t_complex Uconj[NCOLOR*NCOLOR];

            c3x3_conj(Uconj, Link(nx_bwd, mu));
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Uconj, &psi[get_index(nx_bwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_plus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= factors[mu] * tmp2[get_color_spin(color, spinor)];
                }
            }
        }
    }
}

void FermionsWD::OperatorWDHamiltonian(t_complex* res, const t_complex* psi)
{
    // Hamiltonian
    for (int tau = 0; tau < LT; tau++) {
        int tau_coord = tau*vol3D;
        t_complex tmp1[NS3D];
        t_complex tmp2[NS3D];

        HamiltonianWDClover(tmp1, &psi[get_index(tau_coord, 0, 0)], tau);

        for (int x = 0; x < vol3D; x++) {
            for (int color = 0; color < NCOLOR; color++) {
                gamma_times_psi(&tmp2[get_index(x, color, 0)], &tmp1[get_index(x, color, 0)], 0);
            }
            for (int color = 0; color < NCOLOR; color++) {
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x + tau_coord, color, spinor)] = tmp2[get_index(x, color, spinor)];
                }
            }
        }
    }

    // Additional mass term
    for (int x = 0; x < vol; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] += factors[0] * psi[get_index(x, color, spinor)];
            }
        }
    }

    // Forwards temporal
    for (int x = 0; x < vol; x++) {
        int nx_fwd = shift_fwd(x, 0);
        t_complex tmp[NSPINOR*NCOLOR];
        t_complex tmp2[NSPINOR*NCOLOR];

        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Link(x, 0), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_minus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], 0);
        }
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] -= factors[0] * tmp2[get_color_spin(color, spinor)];
            }
        }
    }

    // Backwards temporal
    for (int x = 0; x < vol; x++) {
        int nx_bwd = shift_bwd(x, 0);
        t_complex tmp[NSPINOR*NCOLOR];
        t_complex tmp2[NSPINOR*NCOLOR];
        t_complex Uconj[NCOLOR*NCOLOR];

        c3x3_conj(Uconj, Link(nx_bwd, 0));
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Uconj, &psi[get_index(nx_bwd, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_plus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], 0);
        }
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] -= factors[0] * tmp2[get_color_spin(color, spinor)];
            }
        }
    }
}

void FermionsWD::SpatialTerm(t_complex* res, const t_complex* psi, const int tau)
{
    #pragma omp parallel for
    for (int x = 0; x < vol3D; x++) {
        // Mass term
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = mass_term * psi[get_index(x, color, spinor)];
            }
        }

        CloverTerm(res, psi, x, tau);

        for (int mu = 1; mu < 4; mu++) {
            // Forwards spatial
            int nx_fwd = shift_fwd(x, mu);
            t_complex tmp[NSPINOR*NCOLOR];
            t_complex tmp2[NSPINOR*NCOLOR];

            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Link(tau, x, mu), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_minus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= factors[mu] * tmp2[get_color_spin(color, spinor)];
                }
            }

            // Backwards spatial
            int nx_bwd = shift_bwd(x, mu);
            t_complex Uconj[NCOLOR*NCOLOR];

            c3x3_conj(Uconj, Link(tau, nx_bwd, mu));
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                c3x3_times_c3(&tmp[get_color_spin(0, spinor)], Uconj, &psi[get_index(nx_bwd, 0, spinor)], NSPINOR);
            }
            for (int color = 0; color < NCOLOR; color++) {
                projection_mu_plus_times_psi(&tmp2[get_color_spin(color, 0)], &tmp[get_color_spin(color, 0)], mu);
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    res[get_index(x, color, spinor)] -= factors[mu] * tmp2[get_color_spin(color, spinor)];
                }
            }
        }
    }
}

void FermionsWD::Alpha(t_complex* res, const t_complex* psi, const int tau)
{
    t_complex* tmp_minus = new t_complex[NS3D];
    t_complex* tmp_plus = new t_complex[NS3D];
    t_complex* tmp_spatial = new t_complex[NS3D];

    for (int x = 0; x < vol3D; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_minus_times_psi(&tmp_minus[get_index(x, color, 0)], &psi[get_index(x, color, 0)], 0);
            projection_mu_plus_times_psi(&tmp_plus[get_index(x, color, 0)], &psi[get_index(x, color, 0)], 0);
        }
    }

    SpatialTerm(tmp_spatial, tmp_minus, tau);
    
    for (int x = 0; x < vol3D; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = tmp_spatial[get_index(x, color, spinor)] - factors[0] * tmp_plus[get_index(x, color, spinor)];
            }
        }
    }

    delete[] tmp_minus;
    delete[] tmp_plus;
    delete[] tmp_spatial;
}

void FermionsWD::Beta(t_complex* res, const t_complex* psi, const int tau)
{
    t_complex* Upsi = new t_complex[NS3D];
    t_complex* tmp_minus = new t_complex[NS3D];
    t_complex* tmp_plus = new t_complex[NS3D];
    t_complex* tmp_spatial = new t_complex[NS3D];

    for (int x = 0; x < vol3D; x++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            c3x3_times_c3(&Upsi[get_index(x, 0, spinor)], Link(tau, x, 0), &psi[get_index(x, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_minus_times_psi(&tmp_minus[get_index(x, color, 0)], &Upsi[get_index(x, color, 0)], 0);
            projection_mu_plus_times_psi(&tmp_plus[get_index(x, color, 0)], &Upsi[get_index(x, color, 0)], 0);
        }
    }
    
    SpatialTerm(tmp_spatial, tmp_plus, tau);

    for (int x = 0; x < vol3D; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = tmp_spatial[get_index(x, color, spinor)] - factors[0] * tmp_minus[get_index(x, color, spinor)];
            }
        }
    }

    delete[] Upsi;
    delete[] tmp_minus;
    delete[] tmp_plus;
    delete[] tmp_spatial;
}

void FermionsWD::AlphaInvBeta(t_complex* res, const t_complex* psi)
{
    //Placeholder
    for (int x = 0; x < vol3D; x++)
        for (int color = 0; color < NCOLOR; color++)
            for (int spinor = 0; spinor < NSPINOR; spinor++)
                res[get_index(x, color, spinor)] = psi[get_index(x, color, spinor)];
}

void FermionsWD::IdentityOperator3D(t_complex* res, const t_complex* psi)
{
    for (int x = 0; x < vol3D; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = psi[get_index(x, color, spinor)];
            }
        }
    }
}

void FermionsWD::ProjectionOperator(t_complex* res, const t_complex* psi)
{
    #pragma omp parallel for
    for (int x = 0; x < vol; x++)
    {
        t_complex tmp_minus[NCOLOR*NSPINOR];
        t_complex tmp_plus[NCOLOR*NSPINOR];
        t_complex Utmp[NCOLOR*NSPINOR];

        int nx_fwd = shift_fwd(x, 0);

        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_minus_times_psi(&tmp_minus[get_color_spin(color, 0)], &psi[get_index(x, color, 0)], 0);
        }

        for (int spinor = 0; spinor < NSPINOR; spinor++)
        {
            c3x3_times_c3(&Utmp[get_color_spin(0, spinor)], Link(x, 0), &psi[get_index(nx_fwd, 0, spinor)], NSPINOR);
        }
        for (int color = 0; color < NCOLOR; color++) {
            projection_mu_plus_times_psi(&tmp_plus[get_color_spin(color, 0)], &Utmp[get_color_spin(color, 0)], 0);
            for (int spinor = 0; spinor < NSPINOR; spinor++) {
                res[get_index(x, color, spinor)] = tmp_minus[get_color_spin(color, spinor)] + tmp_plus[get_color_spin(color, spinor)];
            }
        }
    }
}

int FermionsWD::ExtremalEigenstates(const int nev, const int ncv, const double tol, double* evals, t_complex* evecs, int mode)
{
    char which[3] = "SM"; // Both ends of the spectrum
    if(mode==1) sprintf(which, "LM"); // Largest magnitude

    ARrcCompStdEig<double> prob(NS3D, nev, which, ncv, tol, numeric_limits<int>::max());

    cout << ansi::green << "\tFinding " << nev << " eigenvalues with " << ncv << " Arnoldi vectors and tolerance " << tol << ansi::reset << endl;
    while (!prob.ArnoldiBasisFound()) 
    {
        prob.TakeStep();
        if ((prob.GetIdo() == 1) || (prob.GetIdo() == -1)) 
           HamiltonianWDClover(prob.PutVector(), prob.GetVector());
        cout << "." << flush;
    }
    cout << endl;

    prob.FindEigenvectors();

    int res_nev = prob.ConvergedEigenvalues();
    if(prob.EigenvaluesFound() && res_nev>=0 && res_nev<=nev)
    {
        for(int i=0; i<res_nev; i++) evals[i] = prob.Eigenvalue(i).real();
        if(evecs!=nullptr) std::copy(prob.RawEigenvectors(), prob.RawEigenvectors() + res_nev*NS3D, evecs);
    }

    return res_nev;  
}

EigenResult FermionsWD::StdLowestEigenstates(const int nev, const int ncv, const double tol)
{
    ARrcCompStdEig<double> prob(NS3D, nev, "SM", ncv, tol, numeric_limits<int>::max());

    while (!prob.ArnoldiBasisFound()) {
        prob.TakeStep();
        if ((prob.GetIdo() == 1) || (prob.GetIdo() == -1)) {
            HamiltonianWDClover(prob.PutVector(), prob.GetVector());
        }
    }

    prob.FindEigenvectors();

    EigenResult result;
    result.nconv = prob.ConvergedEigenvalues();
    result.NS3D = prob.GetN();
    result.converged = (prob.EigenvaluesFound() && result.nconv == nev);

    if (result.converged) {
        result.eigenvalues.resize(result.nconv);
        result.eigenvectors.resize(result.nconv * result.NS3D);
        for (int i = 0; i < result.nconv; i++) {
            result.eigenvalues[i] = prob.Eigenvalue(i);
            for (int j = 0; j < result.NS3D; j++) {
                result.eigenvectors[i*result.NS3D + j] = prob.Eigenvector(i, j);
            }
        }
    }
    
    return result;
}

void FermionsWD::GenLowestEigenstates(int nev, int ncv, double tol)
{

    ARrcCompGenEig<double> prob(NS3D, nev, "SM", ncv, tol, numeric_limits<int>::max());
	
    while (!prob.ArnoldiBasisFound())
    {
        prob.TakeStep();

        if ((prob.GetIdo() == 1) || (prob.GetIdo() == -1))
        {
            AlphaInvBeta(prob.PutVector(), prob.GetVector());
        }
        if (prob.GetIdo() == 2)
        {
            Alpha(prob.PutVector(), prob.GetVector());
        }
    }

    prob.FindEigenvectors();
    int nconv = prob.ConvergedEigenvalues();

    if (prob.EigenvaluesFound() && nconv==nev)
    {
        cout << ansi::green << "\tEigenvalues: " << ansi::reset << endl;
        for (int i = 0; i < nconv; ++i){
            cout << ansi::yellow << "\t\t\t" << prob.Eigenvalue(i) << ansi::reset << endl;
        }
    }
    else
    {
        cout << ansi::red << "\tDid not converge." << ansi::reset << endl;
        cout << ansi::green << "\tprob.ConvergedEigenvalues = " << ansi::reset;
        cout << ansi::yellow << nconv << ansi::reset << endl;
    }
}

void FermionsWD::OperatorToMatrix(t_complex* res, void (FermionsWD::*O)(t_complex*, const t_complex*), const int n)
{
    t_complex* basis = new t_complex[n];
    t_complex* O_basis = new t_complex[n];

    for (int i = 0; i < n; i++) {
        fill(basis, basis + n, 0.0);
        basis[i] = 1.0;

        (this->*O)(O_basis, basis);
        for (int j = 0; j < n; j++) {
            res[j*n + i] = O_basis[j];
        }
    }

    delete[] basis;
    delete[] O_basis;
}

void FermionsWD::OperatorToMatrix(t_complex* res, void (FermionsWD::*O)(t_complex*, const t_complex*, const int), const int n, const int tau)
{
    t_complex* basis = new t_complex[n];
    t_complex* O_basis = new t_complex[n];

    for (int i = 0; i < n; i++) {
        fill(basis, basis + n, 0.0);
        basis[i] = 1.0;

        (this->*O)(O_basis, basis, tau);
        for (int j = 0; j < n; j++) {
            res[j*n + i] = O_basis[j];
        }
    }

    delete[] basis;
    delete[] O_basis;
}

void FermionsWD::OperatorChaintoMatrix(t_complex* res, void (FermionsWD::*Ops[])(t_complex*, const t_complex*), int n, int num_ops)
{
    t_complex* temp1 = new t_complex[n];
    t_complex* temp2 = new t_complex[n];

    for (int i = 0; i < n; i++) {
        fill(temp1, temp1 + n, 0.0);
        temp1[i] = 1.0;

        t_complex* in = temp1;
        t_complex* out = temp2;

        for (int op_idx = 0; op_idx < num_ops; op_idx++) {
            (this->*Ops[op_idx])(out, in);

            // swap pointers for next operator
            t_complex* tmp = in;
            in = out;
            out = tmp;
        }

        for (int j = 0; j < n; j++) {
            res[j*n + i] = in[j];
        }
    }

    delete[] temp1;
    delete[] temp2;
}

void FermionsWD::MatrixInversion(t_complex* M, int n)
{
    lapack_complex_double* A = reinterpret_cast<lapack_complex_double*>(M);

    int info;
    int* ipiv = new int[n];

    // LU factorization
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
    if (info != 0) {
        std::cerr << "Error in zgetrf: info = " << info << std::endl;
        delete[] ipiv;
    }

    // Matrix inversion
    info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, A, n, ipiv);
    if (info != 0) {
        std::cerr << "Error in zgetri: info = " << info << std::endl;
        delete[] ipiv;
    }

    delete[] ipiv;
}

bool FermionsWD::CheckHermiticity(void (FermionsWD::*O)(t_complex*, const t_complex*, const int))
{
    bool check = false;
    int N = NS3D;
    double tol = 1e-12;

    t_complex* phi = new t_complex[N];
    t_complex* chi = new t_complex[N];
    t_complex* tmp = new t_complex[N];

    init_rng();
    rand_complex(phi, N);
    rand_complex(chi, N);

    (this->*O)(tmp, phi, 0);
    t_complex chi_Hphi = InnerProduct(chi, tmp, N);

    (this->*O)(tmp, chi, 0);
    t_complex Hchi_phi = InnerProduct(tmp, phi, N);

    double relative_diff = abs(chi_Hphi - Hchi_phi) / abs(chi_Hphi);

    if (relative_diff < tol) {
        check = true;
    }

    delete[] phi;
    delete[] chi;
    delete[] tmp;

    return check;
}

bool FermionsWD::CheckGamma5Hermiticity(void (FermionsWD::*O)(t_complex*, const t_complex*))
{
    bool check = false;
    int N = NS;
    double tol = 1e-8;

    t_complex* phi = new t_complex[N];
    t_complex* chi = new t_complex[N];
    t_complex* tmp = new t_complex[N];
    t_complex* comp = new t_complex[N];

    init_rng();
    rand_complex(phi, N);
    rand_complex(chi, N);

    (this->*O)(tmp, phi);
    for (int x = 0; x < vol; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            gamma_times_psi(&comp[get_index(x, color, 0)], &tmp[get_index(x, color, 0)], 5);
        }
    }
    t_complex chi_Hphi = InnerProduct(chi, comp, N);

    (this->*O)(tmp, chi);
    for (int x = 0; x < vol; x++) {
        for (int color = 0; color < NCOLOR; color++) {
            gamma_times_psi(&comp[get_index(x, color, 0)], &tmp[get_index(x, color, 0)], 5);
        }
    }
    t_complex Hchi_phi = InnerProduct(comp, phi, N);

    double relative_diff = abs(chi_Hphi - Hchi_phi) / abs(chi_Hphi);

    if (relative_diff < tol) {
        check = true;
    }

    delete[] phi;
    delete[] chi;
    delete[] tmp;

    return check;
}
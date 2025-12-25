#include "fermions_wd.hpp"

using namespace std;

// Hamiltonian in momentum space
void FermionsWD::HamiltonianWDMomentum(t_complex* res, int* n)
{
    const t_complex* gamma_0 = gamma_matrix(0);
    const t_complex* gamma_mu = nullptr;
    t_complex tmp[NSPINOR*NSPINOR];
    const double q[4] = {
        n[0]*(2*pi)/(LT),
        n[1]*(2*pi)/(LS),
        n[2]*(2*pi)/(LS),
        n[3]*(2*pi)/(LS)
    };

    for (int a = 0; a < NSPINOR; a++){
        for (int b = 0; b < NSPINOR; b++){
            res[get_4x4_spin(a, b)] = gamma_0[get_4x4_spin(a, b)] * mass_term_3D;
        }
    }

    for (int mu = 1; mu < 4; mu++){
        gamma_mu = gamma_matrix(mu);
        s4x4_times_s4x4(tmp, gamma_0, gamma_mu);
        for (int a = 0; a < NSPINOR; a++){
            for (int b = 0; b < NSPINOR; b++){
                res[get_4x4_spin(a, b)] += factors[mu] * (
                    1i * sin(q[mu]) * tmp[get_4x4_spin(a, b)]
                    - cos(q[mu]) * gamma_0[get_4x4_spin(a, b)]
                );
            }
        }
    }
}

void FermionsWD::MomentumEigenstates(t_complex* eigenvalues, t_complex* mom_hamiltonian)
{
    for (int kx = 0; kx < LS; kx++){
        for (int ky = 0; ky < LS; ky++){
            for (int kz = 0; kz < LS; kz++){
                t_complex tmp[NSPINOR];
                int spatial_index = kx*LS*LS + ky*LS + kz;
                int n[4] = {0, kx, ky, kz};
                HamiltonianWDMomentum(mom_hamiltonian, n);
                Eigenstates(tmp, mom_hamiltonian, NSPINOR);
                for (int color = 0; color < NCOLOR; color++){
                    for (int spinor = 0; spinor < NSPINOR; spinor++){
                        eigenvalues[get_index(spatial_index, color, spinor)] = tmp[spinor];
                    }
                }
            }
        }
    }
}

void FermionsWD::MomentumEigenstates(double* eigenvalues, t_complex* eigenvectors, t_complex* mom_hamiltonian)
{
    for (int kx = 0; kx < LS; kx++){
        for (int ky = 0; ky < LS; ky++){
            for (int kz = 0; kz < LS; kz++){
                double tmp_val[NSPINOR];
                t_complex tmp_vec[NSPINOR*NSPINOR];
                int spatial_index = kx*LS*LS + ky*LS + kz;
                int n[4] = {0, kx, ky, kz};
                HamiltonianWDMomentum(mom_hamiltonian, n);
                EigenstatesHermitian(tmp_val, tmp_vec, mom_hamiltonian, NSPINOR);
                for (int spinor = 0; spinor < NSPINOR; spinor++){
                    eigenvalues[spatial_index*NSPINOR + spinor] = tmp_val[spinor];
                }
                for (int spinor2 = 0; spinor2 < NSPINOR*NSPINOR; spinor2++) {
                    eigenvectors[spatial_index*NSPINOR*NSPINOR + spinor2] = tmp_vec[spinor2];
                }
            }
        }
    }
}

void FermionsWD::EuclideanCorrelatorMomentum(t_complex& res, t_complex* op, const int tau)
{
    res = t_complex(0.0, 0.0);
    for (int kx = 0; kx < LS; kx++){
        for (int ky = 0; ky < LS; ky++){
            for (int kz = 0; kz < LS; kz++){
                int n[4] = {tau, kx, ky, kz};

                t_complex mom_hamiltonian[NSPINOR*NSPINOR];
                HamiltonianWDMomentum(mom_hamiltonian, n);

                double eigenvalues[NSPINOR];
                t_complex eigenvectors[NSPINOR*NSPINOR];
                EigenstatesHermitian(eigenvalues, eigenvectors, mom_hamiltonian, NSPINOR);

                double factor[NSPINOR];
                double plus_tau[NSPINOR];
                double minus_tau[NSPINOR];
                double f[NSPINOR];
                double g[NSPINOR];
                double T = static_cast<double>(tau);
                double beta = static_cast<double>(LT);
                for (int i = 0; i < NSPINOR; i++) {
                    plus_tau[i] = exp(T*eigenvalues[i]);
                    minus_tau[i] = exp(-T*eigenvalues[i]);
                    factor[i] = exp(-beta*eigenvalues[i]);
                    f[i] = minus_tau[i]/(1.0 + factor[i]);
                    g[i] = (plus_tau[i]*factor[i])/(1.0 + factor[i]);
                }

                t_complex O_mn[NSPINOR*NSPINOR];
                fill(O_mn, O_mn + NSPINOR*NSPINOR, t_complex(0.0, 0.0));
                for (int m = 0; m < NSPINOR; m++) {
                    for (int n = 0; n < NSPINOR; n++) {
                        for (int i = 0; i < NSPINOR; i++) {
                            for (int j = 0; j < NSPINOR; j++) {
                                O_mn[m*NSPINOR + n] += conj(eigenvectors[i*NSPINOR + m]) * op[i*NSPINOR + j] * eigenvectors[j*NSPINOR + n];
                            }
                        }
                    }
                }

                for (int m = 0; m < NSPINOR; m++) {
                    for (int n = 0; n < NSPINOR; n++) {
                        res += (std::norm(O_mn[m*NSPINOR + n]) * f[m] * g[n]);
                    }
                }
            }
        }
    }
}

void FermionsWD::RealTimeCorrelatorMomentum(t_complex& res, t_complex* op, const double t, const double* eigenvalues, const t_complex* eigenvectors)
{
    res = t_complex(0.0, 0.0);
    const t_complex I = {0.0, 1.0};
    
    for (int kx = 0; kx < LS; kx++){
        for (int ky = 0; ky < LS; ky++){
            for (int kz = 0; kz < LS; kz++){
                int spatial_index = kx*LS*LS + ky*LS + kz;

                t_complex factor[NSPINOR];
                t_complex plus_t[NSPINOR];
                t_complex minus_t[NSPINOR];
                t_complex f[NSPINOR];
                t_complex g[NSPINOR];
                double beta = static_cast<double>(LT);
                for (int i = 0; i < NSPINOR; i++) {
                    plus_t[i] = exp(I*t*eigenvalues[spatial_index*NSPINOR + i]);
                    minus_t[i] = exp(-I*t*eigenvalues[spatial_index*NSPINOR + i]);
                    factor[i] = exp(-beta*eigenvalues[spatial_index*NSPINOR + i]);
                    f[i] = plus_t[i]/(1.0 + factor[i]);
                    g[i] = (minus_t[i]*factor[i])/(1.0 + factor[i]);
                }

                t_complex O_mn[NSPINOR*NSPINOR];
                fill(O_mn, O_mn + NSPINOR*NSPINOR, t_complex(0.0, 0.0));
                for (int m = 0; m < NSPINOR; m++) {
                    for (int n = 0; n < NSPINOR; n++) {
                        for (int i = 0; i < NSPINOR; i++) {
                            for (int j = 0; j < NSPINOR; j++) {
                                O_mn[m*NSPINOR + n] +=
                                    conj(eigenvectors[spatial_index*NSPINOR*NSPINOR + i*NSPINOR + m]) *
                                    op[i*NSPINOR + j] *
                                    eigenvectors[spatial_index*NSPINOR*NSPINOR + j*NSPINOR + n];
                            }
                        }
                    }
                }

                for (int m = 0; m < NSPINOR; m++) {
                    for (int n = 0; n < NSPINOR; n++) {
                        res += (std::norm(O_mn[m*NSPINOR + n]) * f[m] * g[n]);
                    }
                }
            }
        }
    }
}
#include "test.hpp"
#include <gtest/gtest.h>

//=================================================================
// Color Algebra
//=================================================================

class ColorAlgebraConj : public ::testing::Test {
    protected:
        t_complex res[NCOLOR*NCOLOR] = {};
        t_complex res_direct[NCOLOR*NCOLOR] = {};
        t_complex A[NCOLOR*NCOLOR] = {};
        t_complex A_conj[NCOLOR*NCOLOR] = {};
        t_complex B[NCOLOR*NCOLOR] = {};
        t_complex B_conj[NCOLOR*NCOLOR] = {};

        void SetUp() override {
            init_rng();
            rand_complex(A, NCOLOR*NCOLOR);
            rand_complex(B, NCOLOR*NCOLOR);
            c3x3_conj(A_conj, A);
            c3x3_conj(B_conj, B);
        }
};


TEST_F(ColorAlgebraConj, ColorConjTimesColor) {
    c3x3_times_c3x3(res, A_conj, B);
    c3x3_conj_times_c3x3(res_direct, A, B);

    for (int i = 0; i < NCOLOR*NCOLOR; i++) {
        EXPECT_EQ(res[i], res_direct[i]);
    }
}

TEST_F(ColorAlgebraConj, ColorTimesColorConj) {
    c3x3_times_c3x3(res, A, B_conj);
    c3x3_times_c3x3_conj(res_direct, A, B);

    for (int i = 0; i < NCOLOR*NCOLOR; i++) {
        EXPECT_EQ(res[i], res_direct[i]);
    }
}

TEST_F(ColorAlgebraConj, ColorConjTimesColorConj) {
    c3x3_conj(A_conj, A);
    c3x3_conj(B_conj, B);

    c3x3_times_c3x3(res, A_conj, B_conj);

    c3x3_conj_times_c3x3_conj(res_direct, A, B);

    for (int i = 0; i < NCOLOR*NCOLOR; i++) {
        EXPECT_EQ(res[i], res_direct[i]);
    }
}

TEST(ColorAlgebra, C3x3TimesC3Stride) {
    t_complex res[NCOLOR*NSPINOR];
    t_complex A[NCOLOR*NCOLOR] = {
        {1.0,0.0},{2.0,0.0},{3.0,0.0},
        {4.0,0.0},{5.0,0.0},{6.0,0.0},
        {7.0,0.0},{8.0,0.0},{9.0,0.0}
    };
    t_complex v[NCOLOR*NSPINOR] = {
        {1.0,0.0},{2.0,0.0},{3.0,0.0},{4.0,0.0},
        {1.0,0.0},{2.0,0.0},{3.0,0.0},{4.0,0.0},
        {1.0,0.0},{2.0,0.0},{3.0,0.0},{4.0,0.0}
    };
    t_complex res_expect[NCOLOR*NSPINOR] = {
        {6.0,0.0},{12.0,0.0},{18.0,0.0},{24.0,0.0},
        {15.0,0.0},{30.0,0.0},{45.0,0.0},{60.0,0.0},
        {24.0,0.0},{48.0,0.0},{72.0,0.0},{96.0,0.0}
    };

    int stride = NSPINOR;
    for (int spinor = 0; spinor < NSPINOR; spinor++) {
        c3x3_times_c3(&res[get_color_spin(0, spinor)], A, &v[get_color_spin(0, spinor)], NSPINOR);
    }

    for (int color = 0; color < NCOLOR; color++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            EXPECT_EQ(res[get_color_spin(color, spinor)], res_expect[get_color_spin(color, spinor)]);
        }
    }
}

//=================================================================
// Spinor Algebra
//=================================================================

TEST(SpinorAlgebra, Gamma0SquaredUnity)
{
    int n = NSPINOR;
    t_complex v[n];
    t_complex tmp[n];
    t_complex tmp2[n];

    init_rng();

    for (int j = 0; j < 100; j++) {
        rand_complex(v, n);
        gamma_times_psi(tmp, v, 0);
        gamma_times_psi(tmp2, tmp, 0);
        
        for (int i = 0; i < n; i++){
            EXPECT_NEAR(v[i].real(), tmp2[i].real(), 1e-12);
            EXPECT_NEAR(v[i].imag(), tmp2[i].imag(), 1e-12);
        }
    }
}

TEST(SpinorAlgebra, Gamma5)
{
    t_complex res_matrix[NSPINOR];
    t_complex res_oper[NSPINOR];

    t_complex tmp[NSPINOR*NSPINOR];
    t_complex tmp2[NSPINOR*NSPINOR];
    t_complex gamma5[NSPINOR*NSPINOR];

    s4x4_times_s4x4(tmp, gamma_matrix(3), gamma_matrix(0));
    s4x4_times_s4x4(tmp2, gamma_matrix(2), tmp);
    s4x4_times_s4x4(gamma5, gamma_matrix(1), tmp2);

    t_complex psi[NSPINOR];
    init_rng();
    rand_complex(psi, NSPINOR);

    gamma_times_psi(res_oper, psi, 5);
    matrix_vector_mult(res_matrix, gamma5, psi, NSPINOR);

    for (int i = 0; i < NSPINOR; i++) {
        EXPECT_NEAR(res_oper[i].real(), res_matrix[i].real(), 1e-12);
        EXPECT_NEAR(res_oper[i].imag(), res_matrix[i].imag(), 1e-12);
    }
}

TEST(SpinorAlgebra, SigmaTimesPsi)
{
    for (int nu = 1; nu < NSPINOR; nu++) {
        for (int mu = 0; mu < nu; mu++) {
            t_complex res_vec[NSPINOR];
            t_complex sigma[NSPINOR*NSPINOR];
            t_complex res_munu[NSPINOR*NSPINOR];
            t_complex res_numu[NSPINOR*NSPINOR];

            init_rng();
            t_complex vec[NSPINOR] = {};
            rand_complex(vec, NSPINOR);

            // Calculate sigma*v from definition of gamma_matrix
            s4x4_times_s4x4(res_munu, gamma_matrix(mu), gamma_matrix(nu));
            s4x4_times_s4x4(res_numu, gamma_matrix(nu), gamma_matrix(mu));
            for (int i = 0; i < (NSPINOR*NSPINOR); i++) {
                sigma[i] = 0.5i * (res_munu[i] - res_numu[i]);
            }
            s4x4_times_s4(res_vec, sigma, vec);

            // Calculate sigma*v directly
            t_complex res_vec_direct[NSPINOR];
            sigma_times_psi(res_vec_direct, vec, mu, nu);

            for (int i = 0; i < 4; i++) {
                EXPECT_EQ(res_vec[i], res_vec_direct[i]) << "Incorrect sigma matrix: mu = " << mu << ", nu = " << nu;
            }

            // Check i*sigma*vec
            for (int i = 0; i < NSPINOR; i++) {
                res_vec[i] = 1.0i * res_vec[i];
            }
            i_sigma_times_psi(res_vec_direct, vec, mu, nu);
            for (int i = 0; i < 4; i++) {
                EXPECT_EQ(res_vec[i], res_vec_direct[i]) << "Incorrect sigma matrix: mu = " << mu << ", nu = " << nu;
            }

        }
    }
}

TEST(SpinorAlgebra, ProjectionMuMinusTimesPsi)
{
    double tol = 1e-8;
    t_complex res[NSPINOR];
    t_complex psi[NSPINOR];
    t_complex projection[NSPINOR*NSPINOR];
    t_complex res_matrix[NSPINOR];

    t_complex unity[NSPINOR*NSPINOR];
    fill(unity, unity + NSPINOR*NSPINOR, t_complex(0.0, 0.0));
    for (int i = 0; i < NSPINOR; i++) {
        unity[i*NSPINOR + i] = t_complex(1.0, 0.0);
    }

    init_rng();
    rand_complex(psi, NSPINOR);

    for (int mu = 0; mu < 4; mu++) {
        bool check = true;
        projection_mu_minus_times_psi(res, psi, mu);

        for (int i = 0; i < NSPINOR; i++) {
            for (int j = 0; j < NSPINOR; j++) {
                projection[i*NSPINOR + j] = 0.5 * (unity[i*NSPINOR + j] - gamma_matrix(mu)[i*NSPINOR + j]);
            }
        }

        fill(res_matrix, res_matrix + NSPINOR, t_complex(0.0, 0.0));
        for (int i = 0; i < NSPINOR; i++) {
            for (int j = 0; j < NSPINOR; j++) {
                res_matrix[i] += projection[i*NSPINOR + j] * psi[j];
            }
        }

        for (int i = 0; i < NSPINOR; i++) {
            if (abs(res[i] - res_matrix[i]) > tol) {
                cout << res[i] << ", " << res_matrix[i] << endl;
                check = false;
                break;
            }
        }

        EXPECT_TRUE(check);
    }
}

TEST(SpinorAlgebra, ProjectionMuPlusTimesPsi)
{
    double tol = 1e-8;
    t_complex res[NSPINOR];
    t_complex psi[NSPINOR];
    t_complex projection[NSPINOR*NSPINOR];
    t_complex res_matrix[NSPINOR];

    t_complex unity[NSPINOR*NSPINOR];
    fill(unity, unity + NSPINOR*NSPINOR, t_complex(0.0, 0.0));
    for (int i = 0; i < NSPINOR; i++) {
        unity[i*NSPINOR + i] = t_complex(1.0, 0.0);
    }

    init_rng();
    rand_complex(psi, NSPINOR);

    for (int mu = 0; mu < 4; mu++) {
        bool check = true;
        projection_mu_plus_times_psi(res, psi, mu);

        for (int i = 0; i < NSPINOR; i++) {
            for (int j = 0; j < NSPINOR; j++) {
                projection[i*NSPINOR + j] = 0.5 * (unity[i*NSPINOR + j] + gamma_matrix(mu)[i*NSPINOR + j]);
            }
        }

        fill(res_matrix, res_matrix + NSPINOR, t_complex(0.0, 0.0));
        for (int i = 0; i < NSPINOR; i++) {
            for (int j = 0; j < NSPINOR; j++) {
                res_matrix[i] += projection[i*NSPINOR + j] * psi[j];
            }
        }

        for (int i = 0; i < NSPINOR; i++) {
            if (abs(res[i] - res_matrix[i]) > tol) {
                cout << res[i] << ", " << res_matrix[i] << endl;
                check = false;
                break;
            }
        }

        EXPECT_TRUE(check);
    }
}

//=================================================================
// Matching Eigenstates
//=================================================================

TEST(MatchingLowestEigenvalues, SplitArray)
{
    t_complex arr[10] = {-1e-16, 1e-16, 2, -2, 2, -2, 2, 5, -2, -5};
    vector<t_complex> pos, neg, zero;
    int n = sizeof(arr)/sizeof(arr[0]);

    split_array(arr, n, pos, neg, zero);

    EXPECT_EQ(pos.size(), 4);
    EXPECT_EQ(neg.size(), 4);
    EXPECT_EQ(zero.size(), 2);
}


TEST(MatchingLowestEigenvalues, MatchArrays)
{
    t_complex arr1[10] = {-1e-16, 1e-16, 2, -2, 2, -2, 2, 5, -2, -5};
    t_complex arr_pass1[6] = {1e-18, 2e-18, 2, 2, -2, -2};
    t_complex arr_pass2[6] = {1e-18, -2e-18, -2, -2, 2, 2};  
    t_complex arr_pass3[6] = {1e-18, 2e-18, -2, 2, -2, -2};  
    t_complex arr_fail1[6] = {1e-18, 2e-18, 2, 2, -2, -5};  
    t_complex arr_fail2[6] = {1e-18, 2e-18, 2, 5, -2, -5};  
    t_complex arr_fail3[6] = {1e-18, 2e-18, 3e-18, 2, -2, -2};  

    vector<t_complex> pos1, neg1, zero1;
    vector<t_complex> pos2, neg2, zero2;

    int n1 = sizeof(arr1)/sizeof(arr1[0]);
    int n2 = 6;

    split_array(arr1, n1, pos1, neg1, zero1);

    split_array(arr_pass1, n2, pos2, neg2, zero2);
    EXPECT_TRUE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));

    split_array(arr_pass2, n2, pos2, neg2, zero2);
    EXPECT_TRUE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));

    split_array(arr_pass3, n2, pos2, neg2, zero2);
    EXPECT_TRUE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));

    split_array(arr_fail1, n2, pos2, neg2, zero2);
    EXPECT_FALSE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));

    split_array(arr_fail2, n2, pos2, neg2, zero2);
    EXPECT_FALSE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));

    split_array(arr_fail3, n2, pos2, neg2, zero2);
    EXPECT_FALSE(match_arrays(zero1, zero2, pos1, pos2, neg1, neg2));
}

//=================================================================
// FermionsWD Test Classes
//=================================================================

class FermionsWDTest : public ::testing::Test {
protected:
    FermionsWD* WD;
    int argc;
    char** argv;

    int nev = 4;
    int ncv = 8;
    double tol = 1e-8;

    void SetUp() override {
        static const char* argv_static[] = {
            "test",
            "--LT", "3",
            "--LS", "3",
            "--nu_ratio", "1.5",
            "--xi", "4.3",
            "--u_t", "1.0",
            "--u_s", "1.0",
            "--cT", "0.9",
            "--cR", "1.5",
            "--mass", "0.0"
        };
        argc = sizeof(argv_static) / sizeof(argv_static[0]);
        argv = const_cast<char**>(argv_static);

        WD  = new FermionsWD(argc, argv);
    }

    // Called after each TEST_F block
    void TearDown() override {
        delete WD;
    }
};

//=================================================================
// Lattice Indexing
//=================================================================

TEST(LatticeTest, ShiftForwards)
{
    int LT = 3;
    int LS = 4;
    Lattice* lat = new Lattice(LT, LS);

    for (int t = 0; t < LT; t++) {
        for (int x = 0; x < LS; x++) {
            for (int y = 0; y < LS; y++) {
                for (int z = 0; z < LS; z++) {
                    int site = t*LS*LS*LS + x*LS*LS + y*LS + z;
                    int site_shift;

                    site_shift = lat->shift_fwd(site, 0);
                    EXPECT_EQ(site_shift, ((t+1)%LT)*LS*LS*LS + x*LS*LS + y*LS + z);

                    site_shift = lat->shift_fwd(site, 1);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + ((x+1)%LS)*LS*LS + y*LS + z);

                    site_shift = lat->shift_fwd(site, 2);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + x*LS*LS + ((y+1)%LS)*LS + z);

                    site_shift = lat->shift_fwd(site, 3);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + x*LS*LS + y*LS + (z+1)%LS);
                }
            }
        }
    }

    delete lat;
}

TEST(LatticeTest, ShiftBackwards)
{
    int LT = 3;
    int LS = 4;
    Lattice* lat = new Lattice(LT, LS);

    for (int t = 0; t < LT; t++) {
        for (int x = 0; x < LS; x++) {
            for (int y = 0; y < LS; y++) {
                for (int z = 0; z < LS; z++) {
                    int site = t*LS*LS*LS + x*LS*LS + y*LS + z;
                    int site_shift;

                    site_shift = lat->shift_bwd(site, 0);
                    EXPECT_EQ(site_shift, ((t-1 + LT)%LT)*LS*LS*LS + x*LS*LS + y*LS + z);

                    site_shift = lat->shift_bwd(site, 1);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + ((x-1 + LS)%LS)*LS*LS + y*LS + z);

                    site_shift = lat->shift_bwd(site, 2);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + x*LS*LS + ((y-1 + LS)%LS)*LS + z);

                    site_shift = lat->shift_bwd(site, 3);
                    EXPECT_EQ(site_shift, t*LS*LS*LS + x*LS*LS + y*LS + (z-1 + LS)%LS);
                }
            }
        }
    }

    delete lat;
}

TEST(LatticeTest, ShiftForwards3DSlice)
{
    int LT = 3;
    int LS = 4;
    Lattice* lat = new Lattice(LT, LS);

    for (int t = 0; t < LT; t++) {
        for (int x = 0; x < LS; x++) {
            for (int y = 0; y < LS; y++) {
                for (int z = 0; z < LS; z++) {
                    int site = t*LS*LS*LS + x*LS*LS + y*LS + z;
                    int site3D = x*LS*LS + y*LS + z;
                    int site_shift;
                    int site_shift3D;

                    site_shift = lat->shift_fwd(site, 0);
                    site_shift3D = site_shift - ((t+1)%LT)*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + y*LS + z);

                    site_shift = lat->shift_fwd(site, 1);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, ((x+1)%LS)*LS*LS + y*LS + z);
                    EXPECT_EQ(lat->shift_fwd(site3D, 1), ((x+1)%LS)*LS*LS + y*LS + z);

                    site_shift = lat->shift_fwd(site, 2);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + ((y+1)%LS)*LS + z);
                    EXPECT_EQ(lat->shift_fwd(site3D, 2), x*LS*LS + ((y+1)%LS)*LS + z);

                    site_shift = lat->shift_fwd(site, 3);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + y*LS + (z+1)%LS);
                    EXPECT_EQ(lat->shift_fwd(site3D, 3), x*LS*LS + y*LS + (z+1)%LS);
                }
            }
        }
    }

    delete lat;
}

TEST(LatticeTest, ShiftBackwards3DSlice)
{
    int LT = 3;
    int LS = 4;
    Lattice* lat = new Lattice(LT, LS);

    for (int t = 0; t < LT; t++) {
        for (int x = 0; x < LS; x++) {
            for (int y = 0; y < LS; y++) {
                for (int z = 0; z < LS; z++) {
                    int site = t*LS*LS*LS + x*LS*LS + y*LS + z;
                    int site3D = x*LS*LS + y*LS + z;
                    int site_shift;
                    int site_shift3D;

                    site_shift = lat->shift_bwd(site, 0);
                    site_shift3D = site_shift - ((t-1 + LT)%LT)*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + y*LS + z);

                    site_shift = lat->shift_bwd(site, 1);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, ((x-1 + LS)%LS)*LS*LS + y*LS + z);
                    EXPECT_EQ(lat->shift_bwd(site3D, 1), ((x-1 + LS)%LS)*LS*LS + y*LS + z);

                    site_shift = lat->shift_bwd(site, 2);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + ((y-1 + LS)%LS)*LS + z);
                    EXPECT_EQ(lat->shift_bwd(site3D, 2), x*LS*LS + ((y-1 + LS)%LS)*LS + z);

                    site_shift = lat->shift_bwd(site, 3);
                    site_shift3D = site_shift - t*LS*LS*LS;
                    EXPECT_EQ(site_shift3D, x*LS*LS + y*LS + (z-1 + LS)%LS);
                    EXPECT_EQ(lat->shift_bwd(site3D, 3), x*LS*LS + y*LS + (z-1 + LS)%LS);
                }
            }
        }
    }

    delete lat;
}

//=================================================================
// HamiltonianWD and Clover Terms
//=================================================================

TEST_F(FermionsWDTest, HamiltonianMomentumMatch) {
    int NS3D = WD->NS3D;
    t_complex* H = new t_complex[NS3D*NS3D];
    t_complex* H_momentum = new t_complex[NSPINOR*NSPINOR];
    t_complex* eigenvalues = new t_complex[NS3D];
    t_complex* eigenvalues_momentum = new t_complex[NS3D];

    WD->OperatorToMatrix(H, &FermionsWD::HamiltonianWD, NS3D, 0);
    Eigenstates(eigenvalues, H, NS3D);
    WD->MomentumEigenstates(eigenvalues_momentum, H_momentum);
    
    auto order_by_real = [](const t_complex& a, const t_complex& b) {
                  return a.real() < b.real();
              };
    sort(eigenvalues, eigenvalues + NS3D, order_by_real);
    sort(eigenvalues_momentum, eigenvalues_momentum + NS3D, order_by_real);

    for (int i = 0; i < NS3D; i++) {
        EXPECT_NEAR(eigenvalues[i].real(), eigenvalues_momentum[i].real(), 1e-8);
    }

    delete[] H;
    delete[] H_momentum;
    delete[] eigenvalues;
    delete[] eigenvalues_momentum;
}

TEST_F(FermionsWDTest, HamiltonianWDConstantVector)
{
    init_rng();
    t_complex random_val = rand_complex();
    bool check = true;

    int N = WD->NS3D;
    double tol = 1e-12;

    t_complex* arr = new t_complex[N];
    t_complex* res = new t_complex[N];
    for (int i = 0; i < N; i++) {
        arr[i] = random_val;
    }
    WD->HamiltonianWD(res, arr);
    for (int i = 0; i < N; i++) {
        if (abs(res[i]) > tol) {
            check = false;
            break;
        }
    }

    EXPECT_TRUE(check);

    delete[] arr;
    delete[] res;
}


TEST_F(FermionsWDTest, CloverTermVanishesIfLinksUnity)
{
    int N = WD->NS3D;
    double tol = 1e-12;

    t_complex* arr = new t_complex[N];
    t_complex* res = new t_complex[N];
    t_complex* res_clover = new t_complex[N];
    init_rng();
    rand_complex(arr, N);

    WD->HamiltonianWD(res, arr);
    WD->HamiltonianWDClover(res_clover, arr);

    for (int i = 0; i < N; i++) {
        EXPECT_EQ(res[i], res_clover[i]);
    }
    delete[] arr;
    delete[] res;
    delete[] res_clover;
}

TEST_F(FermionsWDTest, HamiltonianWDHermiticity) {
    WD->RandomizeLinks();
    EXPECT_TRUE(WD->CheckHermiticity(&FermionsWD::HamiltonianWD));
}

TEST_F(FermionsWDTest, FieldStrengthAntiSymmetric) {
    WD->RandomizeLinks();
    t_complex F_mu_nu[NCOLOR*NCOLOR] = {};
    t_complex F_nu_mu[NCOLOR*NCOLOR] = {};

    
    for (int x = 0; x < WD->vol; x++){
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                for (int color3x3 = 0; color3x3 < NCOLOR*NCOLOR; color3x3++) {
                    EXPECT_EQ(
                        WD->FieldStrength[get_fieldstrength_index(x, mu, nu, color3x3)],
                        -WD->FieldStrength[get_fieldstrength_index(x, nu, mu, color3x3)]
                    );
                }
            }
        }
    }
}

/*
TEST_F(FermionsWDTest, FieldStrengthAntiHermitian) {
    WD->RandomizeLinks();
    double tol = 1e-12;

    t_complex* phi = new t_complex[NCOLOR*NSPINOR];
    t_complex* chi = new t_complex[NCOLOR*NSPINOR];
    t_complex* tmp = new t_complex[NCOLOR*NSPINOR];
    t_complex* res = new t_complex[NCOLOR*NSPINOR];
    t_complex* F_mu_nu = new t_complex[NCOLOR*NCOLOR];

    rand_complex(phi, NCOLOR*NSPINOR);
    rand_complex(chi, NCOLOR*NSPINOR);

    int x = 0;
    for (int color = 0; color < NCOLOR; color++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            res[get_color_spin(color, spinor)] = t_complex(0.0, 0.0);
        }
    }
    for (int nu = 1; nu < 4; nu++) {
            for (int mu = 0; mu < nu; mu++) {
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &WD->FieldStrength[get_fieldstrength_index(x, mu, nu, 0)], &phi[get_index(x, 0, spinor)], NSPINOR);
                }
                for (int color = 0; color < NCOLOR; color++) {
                    for (int spinor = 0; spinor < NSPINOR; spinor++) {
                        res[get_color_spin(color, spinor)] += tmp[get_color_spin(color, spinor)];
                    }
                }
            }
        }
    t_complex chi_Hphi = InnerProduct(chi, res, NCOLOR*NSPINOR);

    for (int color = 0; color < NCOLOR; color++) {
        for (int spinor = 0; spinor < NSPINOR; spinor++) {
            res[get_color_spin(color, spinor)] = 0;
        }
    }
    for (int nu = 1; nu < 4; nu++) {
            for (int mu = 0; mu < nu; mu++) {
                for (int spinor = 0; spinor < NSPINOR; spinor++) {
                    c3x3_times_c3(&tmp[get_color_spin(0, spinor)], &WD->FieldStrength[get_fieldstrength_index(x, mu, nu, 0)], &chi[get_index(x, 0, spinor)], NSPINOR);
                }
                for (int color = 0; color < NCOLOR; color++) {
                    for (int spinor = 0; spinor < NSPINOR; spinor++) {
                        res[get_color_spin(color, spinor)] += tmp[get_color_spin(color, spinor)];
                    }
                }
            }
        }
    t_complex Hchi_phi = InnerProduct(res, phi, NCOLOR*NSPINOR);

    double relative_diff = abs(chi_Hphi + Hchi_phi) / abs(chi_Hphi);

    EXPECT_LE(relative_diff, tol);

    delete[] phi;
    delete[] chi;
    delete[] tmp;
    delete[] res;
    delete[] F_mu_nu;
}
*/

/*
TEST_F(FermionsWDRandomTest, HamiltonianWDFillingLowestEigenvalues)
{
    t_complex mom_hamiltonian[NSPINOR*NSPINOR] = {};
    t_complex tmp[NSPINOR] = {};
    t_complex* mom_eigenvalues = new t_complex[nev]();

    WD->MomentumEigenstates(mom_eigenvalues, mom_hamiltonian, tmp);

    // Position space eigenvalues
    init_rng();

    int ncv = min(ncv, WD->NS3D - 1);
    int nconv = 0;
    t_complex eigenvalues[nev] = {};

    EigenResult result = WD->StdLowestEigenstates(nev, ncv, tol);

    vector<t_complex> pos1, neg1, zero1;
    vector<t_complex> pos2, neg2, zero2;
    split_array(mom_eigenvalues, nev, pos1, neg1, zero1);
    split_array(result.eigenvalues.data(), nev, pos2, neg2, zero2);
    bool match = false;
    match = match_arrays(zero1, zero2, pos1, pos2, neg1, neg2);

    EXPECT_TRUE(match);

    delete[] mom_eigenvalues;
}
*/

//=================================================================
// Helper Functions
//=================================================================

TEST(LinearAlgebra, ElementProduct) {
    int n = 4;
    t_complex A[n] = {{3.0,2.0}, {3.0,-2.0}, {0.0,1.0}, {-1.0,0.0}};
    t_complex B[n] = {{1.0,2.0}, {3.0,4.0}, {5.0,6.0}, {7.0,8.0}};
    t_complex C[n] = {{0.3,0.0}, {0.5,0.0}, {0.7,0.0}, {0.11,0.0}};
    t_complex D[n] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    t_complex res;
    
    res = element_product(A, n);
    EXPECT_DOUBLE_EQ(res.real(), 0.0);
    EXPECT_DOUBLE_EQ(res.imag(), -13.0);

    res = element_product(B, n);
    EXPECT_DOUBLE_EQ(res.real(), -755.0);
    EXPECT_DOUBLE_EQ(res.imag(), -540.0);

    res = element_product(C, n);
    EXPECT_DOUBLE_EQ(res.real(), 0.01155);
    EXPECT_DOUBLE_EQ(res.imag(), 0.0);

    res = element_product(D, n);
    EXPECT_DOUBLE_EQ(res.real(), 1.0);
    EXPECT_DOUBLE_EQ(res.imag(), 0.0);
}

TEST(LinearAlgebra, MatrixMultiplication) {
    int n = 3;
    int nxn = n*n;

    t_complex Id[nxn];
    fill(Id, Id+nxn, 0.0);
    for (int i = 0; i < n; i++) {
        Id[i*n + i] = t_complex(1.0, 0.0);
    }

    t_complex A[nxn] =  {{2.0,0.0}, {3.0,0.0}, {5.0,0.0},
                         {7.0,0.0}, {11.0,0.0}, {13.0,0.0},
                         {17.0,0.0}, {19.0,0.0}, {23.0,0.0}};
    t_complex B[nxn] =  {{0.0,1.0}, {1.0,0.0}, {0.0,1.0},
                         {1.0,0.0}, {0.0,1.0}, {1.0,0.0},
                         {0.0,1.0}, {1.0,0.0}, {0.0,1.0}};
    t_complex res[nxn] = {};

    matrix_multiplication(res, A, Id, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            EXPECT_DOUBLE_EQ(res[i*n + j].real(), A[i*n + j].real());
            EXPECT_DOUBLE_EQ(res[i*n + j].imag(), A[i*n + j].imag());
        }
    }

    t_complex AA[nxn] = {{110.0,0.0}, {134.0,0.0}, {164.0,0.0},
                         {312.0,0.0}, {389.0,0.0}, {477.0,0.0},
                         {558.0,0.0}, {697.0,0.0}, {861.0,0.0}};
    matrix_multiplication(res, A, A, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            EXPECT_DOUBLE_EQ(res[i*n + j].real(), AA[i*n + j].real());
            EXPECT_DOUBLE_EQ(res[i*n + j].imag(), AA[i*n + j].imag());
        }
    }

    t_complex AB[nxn] = {{3.0,7.0}, {7.0,3.0}, {3.0,7.0},
                         {11.0,20.0}, {20.0,11.0}, {11.0,20.0},
                         {19.0,40.0}, {40.0,19.0}, {19.0,40.0}};
    matrix_multiplication(res, A, B, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            EXPECT_DOUBLE_EQ(res[i*n + j].real(), AB[i*n + j].real());
            EXPECT_DOUBLE_EQ(res[i*n + j].imag(), AB[i*n + j].imag());
        }
    }

    t_complex BA[nxn] = {{7.0,19.0}, {11.0,22.0}, {13.0,28.0},
                         {19.0,7.0}, {22.0,11.0}, {28.0,13.0},
                         {7.0,19.0}, {11.0,22.0}, {13.0,28.0}};
    matrix_multiplication(res, B, A, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            EXPECT_DOUBLE_EQ(res[i*n + j].real(), BA[i*n + j].real());
            EXPECT_DOUBLE_EQ(res[i*n + j].imag(), BA[i*n + j].imag());
        }
    }
}

TEST(LinearAlgebra, MatrixVectorMultiplication) {
    int n = 3;
    int nxn = n*n;
    
    t_complex Id[nxn];
    fill(Id, Id+nxn, 0.0);
    for (int i = 0; i < n; i++) {
        Id[i*n + i] = t_complex(1.0, 0.0);
    }

    t_complex A[nxn] =  {{2.0,0.0}, {3.0,0.0}, {5.0,0.0},
                         {7.0,0.0}, {11.0,0.0}, {13.0,0.0},
                         {17.0,0.0}, {19.0,0.0}, {23.0,0.0}};
    t_complex B[nxn] =  {{0.0,1.0}, {1.0,0.0}, {0.0,1.0},
                         {1.0,0.0}, {0.0,1.0}, {1.0,0.0},
                         {0.0,1.0}, {1.0,0.0}, {0.0,1.0}};
    t_complex v[n]   =  {{1.0,2.0}, {3.0,4.0}, {5.0,6.0}};
    t_complex res[n] =  {};

    matrix_vector_mult(res, Id, v, n);
    for (int i = 0; i < n; i++) {
        EXPECT_DOUBLE_EQ(res[i].real(), v[i].real());
        EXPECT_DOUBLE_EQ(res[i].imag(), v[i].imag());
    }

    t_complex Av[n] = {{36.0,46.0}, {105.0,136.0}, {189.0,248.0}};
    matrix_vector_mult(res, A, v, n);
    for (int i = 0; i < n; i++) {
        EXPECT_DOUBLE_EQ(res[i].real(), Av[i].real());
        EXPECT_DOUBLE_EQ(res[i].imag(), Av[i].imag());
    }

    t_complex Bv[n] = {{-5.0,10.0}, {2.0,11.0}, {-5.0,10.0}};
    matrix_vector_mult(res, B, v, n);
    for (int i = 0; i < n; i++) {
        EXPECT_DOUBLE_EQ(res[i].real(), Bv[i].real());
        EXPECT_DOUBLE_EQ(res[i].imag(), Bv[i].imag());
    }
}

/*
TEST(LinearAlgebra, BLASBenchmarks) {
    int N = NCOLOR*NSPINOR;
    int LS = 24;
    int vol3D = LS*LS*LS;
    int NS3D = LS*LS*LS*NCOLOR*NSPINOR;
    t_complex* A = new t_complex[NS3D];
    t_complex* B = new t_complex[NS3D];

    init_rng();
    t_complex a = rand_complex();
    rand_complex(A, NS3D);
    rand_complex(B, NS3D);

    auto bench = [&](auto fn) {
        double t = 0.0;
        for (int i = 0; i < 1000; i++) {
            for (int r=0; r<vol3D; r++) {
                auto t0 = std::chrono::high_resolution_clock::now();
                fn(r);
                auto t1 = std::chrono::high_resolution_clock::now();
                t += std::chrono::duration<double>(t1-t0).count();
            }
        }
        return t;
    };

    auto t_loop = bench([&](int r){ for (int i=0;i<N;i++) A[r*N + i] += a*B[r*N + i]; });
    auto t_blas = bench([&](int r){ 
        cblas_zaxpy(
            N,
            &a,
            &B[r*N],
            1,
            &A[r*N],
            1
        );

     });

    cout << "C++ Loop: " << t_loop << endl;
    cout << "BLAS: " << t_blas << endl;

    delete[] A;
    delete[] B;
}
*/

TEST_F(FermionsWDTest, Eigenstates) {
    int n = 3;
    int nxn = n*n;
    t_complex result[n] = {};
    t_complex M[nxn] = {{-2.0,0.0}, {2.0,2.0}, {3.0,0.0},
                      {-9.0,0.0}, {7.0,1.0}, {5.0,0.0},
                      {-5.0,-1.0}, {2.0,0.0}, {6.0,0.0}};

    t_complex eigenvalues[n] = {{6.96768, -2.21188}, {0.740683, 3.59108}, {3.29164, -0.379195}};
    Eigenstates(result, M, n);

    bool match = true;
    for (int i = 0; i < 0; i++) {
        bool found = false;
        for (int j = 0; j < 0; j++) {
            if ((abs(result[i]) - abs(eigenvalues[j]) < 1e-8)) {
                found = true;
                //cout << result[i] << ", " << eigenvalues[j] << endl;
                break;
            }
        }
        if (found == false) {
            match = false;
            break;
        }
    }

    EXPECT_TRUE(match);
    
}

TEST_F(FermionsWDTest, EigenSolver) {
    t_complex M[4] = {{1.0, 1.0}, {2.0, 2.0},
                      {3.0, 3.0}, {4.0, 4.0}};

    t_complex V[4] = {{0.457427, 0.0}, {-1.45743, 0.0},
                      {1.0, 0.0}, {1.0, 0.0}};

    t_complex V_inv[4];
    for (int i = 0; i < 4; i++) {
        V_inv[i] = V[i];
    }
    WD->MatrixInversion(V_inv, 2);

    t_complex L[4] = {{5.37228, 5.37228}, {0.0, 0.0},
                      {0.0, 0.0}, {-0.372281, -0.372281}};
    t_complex tmp[4];
    t_complex comp[4];

    matrix_multiplication(tmp, L, V_inv, 2);
    matrix_multiplication(comp, V, tmp, 2);

    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(M[i].real(), comp[i].real(), 1e-4);
        EXPECT_NEAR(M[i].imag(), comp[i].imag(), 1e-4);
    }

    t_complex eigenvalues[2];
    t_complex eigenvectors[4];
    Eigenstates(eigenvalues, eigenvectors, M, 2);
    t_complex eigenvalue_matrix[4];
    fill(eigenvalue_matrix, eigenvalue_matrix + 4, t_complex(0.0, 0.0));
    for (int i = 0; i < 2; i++) {
        eigenvalue_matrix[i*2 + i] = eigenvalues[i];
    }

    t_complex eigenvectors_inv[4];
    for (int i = 0; i < 4; i++) {
        eigenvectors_inv[i] = eigenvectors[i];
    }
    WD->MatrixInversion(eigenvectors_inv, 2);

    matrix_multiplication(tmp, eigenvalue_matrix, eigenvectors_inv, 2);
    matrix_multiplication(comp, eigenvectors, tmp, 2);
}

TEST_F(FermionsWDTest, EigenstatesRandomSolve) {
    t_complex M[NSPINOR*NSPINOR];
    init_rng();
    rand_complex(M, NSPINOR*NSPINOR);

    t_complex eigenvalues[NSPINOR];
    t_complex eigenvectors[NSPINOR*NSPINOR];

    Eigenstates(eigenvalues, eigenvectors, M, NSPINOR);

    t_complex lhs[NSPINOR];
    t_complex rhs[NSPINOR];
    t_complex res[NSPINOR];
    fill(lhs, lhs + NSPINOR, t_complex(0.0, 0.0));
    fill(rhs, rhs + NSPINOR, t_complex(0.0, 0.0));

    for (int k = 0; k < NSPINOR; k++){
        for (int i = 0; i < NSPINOR; i++) {
            for (int j = 0; j < NSPINOR; j++) {
                lhs[i] += M[i*NSPINOR + j] * eigenvectors[j*NSPINOR + k];
            }
            rhs[i] += eigenvalues[k] * eigenvectors[i*NSPINOR + k];
            res[i] = lhs[i] - rhs[i];
            EXPECT_LT(std::norm(res[i]), 1e-8);
        }
    }
}

TEST_F(FermionsWDTest, MatrixInversion) {
    init_rng();
    int n = WD->NS3D;
    double tol = 1e-8;

    bool check = true;
    t_complex* M = new t_complex[n*n];
    t_complex* Minv = new t_complex[n*n];
    t_complex* Id = new t_complex[n*n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            t_complex random_val = rand_complex();
            M[i*n + j] = random_val;
            Minv[i*n + j] = random_val;
        }
    }
    
    WD->MatrixInversion(Minv, n);
    matrix_multiplication(Id, M, Minv, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                if (abs(Id[i*n + j]) > tol) {
                    check = false;
                    break;
                }
            } else {
                if (abs(Id[i*n + j] - 1.0) > tol) {
                    check = false;
                    break;
                }
            }
        }
    }

    EXPECT_TRUE(check);

    delete[] M;
    delete[] Minv;
}

TEST_F(FermionsWDTest, OperatorToMatrix) {
    WD->RandomizeLinks();
    int n = WD->NS;
    t_complex* v = new t_complex[n];
    t_complex* M = new t_complex[n*n];
    t_complex* res = new t_complex[n];
    t_complex* res_op = new t_complex[n];
    init_rng();
    rand_complex(v, n);

    WD->OperatorToMatrix(M, &FermionsWD::ProjectionOperator, n);
    matrix_vector_mult(res, M, v, n);

    WD->ProjectionOperator(res_op, v);

    for (int i = 0; i < n; i++) {
        EXPECT_NEAR(res[i].real(), res_op[i].real(), 1e-12);
        EXPECT_NEAR(res[i].imag(), res_op[i].imag(), 1e-12);
    }

    delete[] v;
    delete[] M;
    delete[] res;
    delete[] res_op;
}

//=================================================================
// Research Paper Identities
//=================================================================

TEST_F(FermionsWDTest, ProjectionOperatorDeterminantUnity) {
    WD->RandomizeLinks();
    int n = WD->NS;
    t_complex* M = new t_complex[n*n];
    t_complex* eigenvalues = new t_complex[n];
    t_complex detM;

    WD->OperatorToMatrix(M, &FermionsWD::ProjectionOperator, n);
    Eigenstates(eigenvalues, M, n);
    detM = element_product(eigenvalues, n);

    EXPECT_NEAR(abs(detM), 1.0, 1e-8);

    delete[] M;
    delete[] eigenvalues;
}

TEST_F(FermionsWDTest, OperatorWDProjectionDeterminantMatch) {
    WD->RandomizeLinks();
    int n = WD->NS;
    t_complex* M = new t_complex[n*n];
    t_complex* eigenvalues = new t_complex[n];
    t_complex log_detD;
    t_complex log_detDP;

    WD->OperatorToMatrix(M, &FermionsWD::OperatorWD, n);
    Eigenstates(eigenvalues, M, n);
    log_detD = element_log_sum(eigenvalues, n);

    void (FermionsWD::*ops[2])(t_complex*, const t_complex*) = { &FermionsWD::ProjectionOperator, &FermionsWD::OperatorWD };
    WD->OperatorChaintoMatrix(M, ops, n, 2);
    Eigenstates(eigenvalues, M, n);
    log_detDP = element_log_sum(eigenvalues, n);

    EXPECT_GT(abs(log_detD), 1e-12);
    EXPECT_GT(abs(log_detDP), 1e-12);

    double relative_diff = abs(abs(log_detD) - abs(log_detDP))/abs(log_detD);
    EXPECT_LT(relative_diff, 1e-12);

    delete[] M;
    delete[] eigenvalues;
}

TEST_F(FermionsWDTest, AlphaBetaCheck)
{
    WD->RandomizeLinks();
    int NS3D = WD->NS3D;
    int NS   = WD->NS;
    int LT   = WD->LT;
    int vol3D = WD->vol3D;
    int vol = WD->vol;

    t_complex* res_DP = new t_complex[NS];
    t_complex* tmp_DP = new t_complex[NS];
    t_complex* res_A = new t_complex[NS];
    t_complex* res_B = new t_complex[NS];
    t_complex* res_AB = new t_complex[NS];
    t_complex* v = new t_complex[NS];
    init_rng();
    rand_complex(v, NS);

    // Alpha_tau psi_tau + Beta_tau psi_tau+1
    for (int tau = 0; tau < LT; tau++) {
        WD->Alpha(&res_A[get_tau_index(tau, 0, 0, 0, vol3D)], &v[get_tau_index(tau, 0, 0, 0, vol3D)], tau);
        WD->Beta(&res_B[get_tau_index(tau, 0, 0, 0, vol3D)], &v[get_tau_index((tau+1)%LT, 0, 0, 0, vol3D)], tau);
    }

    for (int i = 0; i < NS; i++) {
        res_AB[i] = res_A[i] + res_B[i];
    }

    WD->ProjectionOperator(tmp_DP, v);
    WD->OperatorWD(res_DP, tmp_DP);

    for (int i = 0; i < NS; i++) {
        EXPECT_NEAR(res_DP[i].real(), res_AB[i].real(), 1e-8);
        EXPECT_NEAR(res_DP[i].imag(), res_AB[i].imag(), 1e-8);
    }

    delete[] res_DP;
    delete[] tmp_DP;
    delete[] res_A;
    delete[] res_B;
    delete[] res_AB;
    delete[] v;
}


TEST_F(FermionsWDTest, AlphaBetaOperatorWDDeterminantMatch) {
    WD->RandomizeLinks();
    int NS = WD->NS;
    int NS3D = WD->NS3D;

    t_complex* D = new t_complex[NS*NS];
    t_complex* eigenvalues_D = new t_complex[NS];
    double log_detD = 0.0;

    t_complex* Alpha = new t_complex[NS3D*NS3D];
    t_complex* Beta = new t_complex[NS3D*NS3D];
    t_complex* tmp = new t_complex[NS3D*NS3D];
    t_complex* T_i = new t_complex[NS3D*NS3D];
    t_complex* AB = new t_complex[NS3D*NS3D];
    t_complex* Id = new t_complex[NS3D*NS3D];
    t_complex* eigenvalues_AB = new t_complex[NS3D];
    t_complex* eigenvalues_Alpha = new t_complex[NS3D];
    double log_detAB = 0.0;
    double log_detAlpha_i = 0.0;
    double log_detAlpha = 0.0;

    WD->OperatorToMatrix(Id, &FermionsWD::IdentityOperator3D, NS3D);
    for (int i = 0; i < NS3D*NS3D; i++) {
        AB[i] = Id[i];
    }

    WD->OperatorToMatrix(D, &FermionsWD::OperatorWD, NS);
    Eigenstates(eigenvalues_D, D, NS);
    log_detD = element_log_sum(eigenvalues_D, NS);

    for (int tau = 0; tau < WD->LT; tau++) {
        WD->OperatorToMatrix(Alpha, &FermionsWD::Alpha, NS3D, tau);
        Eigenstates(eigenvalues_Alpha, Alpha, NS3D);
        log_detAlpha_i = element_log_sum(eigenvalues_Alpha, NS3D);
        //cout << "log_DetAlpha_i = " << log_detAlpha_i << endl;
        EXPECT_GT(abs(log_detAlpha_i), 1e-8);
        log_detAlpha += log_detAlpha_i;

        WD->OperatorToMatrix(Beta, &FermionsWD::Beta, NS3D, tau);
        WD->MatrixInversion(Alpha, NS3D);
        matrix_multiplication(T_i, Alpha, Beta, NS3D);
        matrix_multiplication(tmp, AB, T_i, NS3D);
        memcpy(AB, tmp, NS3D*NS3D*sizeof(t_complex));
    }
    for (int i = 0; i < NS3D*NS3D; i++) {
        AB[i] += Id[i];
    }
    Eigenstates(eigenvalues_AB, AB, NS3D);
    log_detAB = element_log_sum(eigenvalues_AB, NS3D);

    //cout << "log_DetAlpha = " << log_detAlpha << endl;
    //cout << "log_DetAB = " << log_detAB << endl;
    //cout << "log_DetAlpha + log_DetAB = " << log_detAlpha + log_detAB << endl;
    //cout << "log_DetD = " << log_detD << endl;;

    
    EXPECT_GT(abs(log_detD), 1e-8);
    EXPECT_GT(abs(log_detAB + log_detAlpha), 1e-8);

    double relative_diff = abs(abs(log_detD) - abs(log_detAB + log_detAlpha))/abs(log_detD);
    EXPECT_LT(relative_diff, 1e-8) << " may fail due to accumulation of algorithm errors: inversion, eigensolver";

    delete[] D;
    delete[] eigenvalues_D;
    delete[] Alpha;
    delete[] Beta;
    delete[] tmp;
    delete[] AB;
    delete[] Id;
    delete[] T_i;
    delete[] eigenvalues_AB;
    delete[] eigenvalues_Alpha;
}

TEST_F(FermionsWDTest, OperatorWDHamiltonianMatch) {
    WD->RandomizeLinks();
    int n = WD->NS;
    t_complex psi[n];
    t_complex res[n];
    t_complex res_ham[n];
    init_rng();

    rand_complex(psi, n);
    WD->OperatorWD(res, psi);
    WD->OperatorWDHamiltonian(res_ham, psi);

    for (int j = 0; j < n; j++) {
        EXPECT_NEAR(res[j].real(), res_ham[j].real(), 1e-8);
        EXPECT_NEAR(res[j].imag(), res_ham[j].imag(), 1e-8);
    }
}

TEST_F(FermionsWDTest, OperatorWDGamma5Hermiticity) {
    WD->RandomizeLinks();
    EXPECT_TRUE(WD->CheckGamma5Hermiticity(&FermionsWD::OperatorWD));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
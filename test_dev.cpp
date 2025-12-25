#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>

#include "linalg.hpp"
#include "fermions_wd.hpp"
#include "utils.hpp"
#include <gtest/gtest.h>


using namespace std;

typedef complex<double>       t_complex;

// =======================================================================================
// CLASSES
// =======================================================================================

class FermionsWDTest : public ::testing::Test {
protected:
    FermionsWD* WD;
    int argc;
    char** argv;

    int nev = 4;
    int ncv = 8;
    double tol = 1e-8;

    // --LT 16 --LS 24 --nu_ratio 1.265 --xi 4.3 --u_t 1.0 --u_s 0.92674 --cR 1.5893 --cT 0.90278 --mass -0.084
    void SetUp() override {
        static const char* argv_static[] = {
            "test",
            "--LT", "16",
            "--LS", "24",
            "--nu_ratio", "1.265",
            "--xi", "4.3",
            "--u_t", "1.0",
            "--u_s", "0.92674",
            "--cT", "0.90278",
            "--cR", "1.5893",
            "--mass", "-0.084"  // "-0.084"
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


// =======================================================================================
// TESTS
// =======================================================================================

/*
TEST_F(FermionsWDTest, FieldStrengthPrecompute) {
    WD->RandomizeLinks();

    int LT = WD->LT;
    int vol3D = WD->vol3D;
    int vol = WD->vol;
    int NS3D = WD->NS3D;

    t_complex* res = new t_complex[NS3D];
    t_complex* res_pre = new t_complex[NS3D];
    t_complex* psi = new t_complex[NS3D];
    init_rng();
    rand_complex(psi, NS3D);

    cout << "F_tau size (Bytes): " << sizeof(t_complex)*vol*4*4*NCOLOR*NCOLOR << endl;
    cout << "F_tau size (MB): " << sizeof(t_complex)*vol*4*4*NCOLOR*NCOLOR/1e6 << endl;

    for (int tau = 0; tau < LT; tau++) {
        auto start = std::chrono::high_resolution_clock::now();
        WD->FieldStrengthCompute(WD->FieldStrength);
        auto checkpoint1 = std::chrono::high_resolution_clock::now();

        WD->HamiltonianWDCloverPrecompute(res_pre, psi, tau);
        auto checkpoint2 = std::chrono::high_resolution_clock::now();

        WD->HamiltonianWD(res, psi, tau);
        auto checkpoint3 = std::chrono::high_resolution_clock::now();

        WD->HamiltonianWDClover(res, psi, tau);
        auto checkpoint4 = std::chrono::high_resolution_clock::now();

        auto F_duration = duration_cast<std::chrono::microseconds>(checkpoint1 - start);
        auto H_pre_duration = duration_cast<std::chrono::microseconds>(checkpoint2 - checkpoint1);
        auto H_bare_duration = duration_cast<std::chrono::microseconds>(checkpoint3 - checkpoint2);
        auto H_duration = duration_cast<std::chrono::microseconds>(checkpoint4 - checkpoint3);

        cout << "F_tau: " << F_duration.count();
        cout << ", H_pre: " << H_pre_duration.count();
        cout << ", H: " << H_duration.count();
        cout << ", H_bare: " << H_bare_duration.count() << endl;

        for (int i = 0; i < NS3D; i++) {
            EXPECT_NEAR(res_pre[i].real(), res[i].real(), 1e-12);
            EXPECT_NEAR(res_pre[i].imag(), res[i].imag(), 1e-12);
        }
    }

    delete[] res;
    delete[] res_pre;
    delete[] psi;
}
*/


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
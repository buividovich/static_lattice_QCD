#include "utils.hpp"

void saveEigenvaluesToBinary(const EigenResult& result, const int LS, const string& fname)
{
    ofstream fout(fname, ios::binary);

    // Binary file header
    fout.write(reinterpret_cast<const char*>(&result.nconv), sizeof(int));
    fout.write(reinterpret_cast<const char*>(&LS), sizeof(int));

    for (int i = 0; i < result.nconv; i++) {
        double re = result.eigenvalues[i].real();
        double im = result.eigenvalues[i].imag();
        fout.write(reinterpret_cast<const char*>(&re), sizeof(double));
        fout.write(reinterpret_cast<const char*>(&im), sizeof(double));
    }

    fout.close();
}

void saveEigenvectorsToBinary(const EigenResult& result, const int LS, const string& fname)
{
    ofstream fout(fname, ios::binary);

    for (int i = 0; i < result.nconv; i++) {
        for (int j = 0; j < result.NS3D; j++) {
            double re = result.eigenvectors[i*result.NS3D + j].real();
            double im = result.eigenvectors[i*result.NS3D + j].imag();
            fout.write(reinterpret_cast<const char*>(&re), sizeof(double));
            fout.write(reinterpret_cast<const char*>(&im), sizeof(double));
        }   
    }
    fout.close();
}

void printEigenvaluesAndSave(const EigenResult& result, const int LS, const string& eigval_fname, const string& eigvec_fname)
{
    if (result.converged) {
        cout << ansi::green << "\tConverged with " << result.nconv << " eigenvalues:" << ansi::reset << endl;
        for (auto &eigenvalue : result.eigenvalues) {
            cout << ansi::yellow << "\t\t\t" << eigenvalue << ansi::reset << endl;
        }
        saveEigenvaluesToBinary(result, LS, eigval_fname);
        saveEigenvectorsToBinary(result, LS, eigvec_fname);
    } else {
        cout << ansi::red << "\tDid not converge." << ansi::reset << endl;
        cout << ansi::green << "Found " << ansi::reset;
        cout << ansi::yellow << result.nconv << " eigenvalues." << ansi::reset << endl;
    }
}
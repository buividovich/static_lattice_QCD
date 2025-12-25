#ifndef _GAUGE_FIELD_HPP
#define _GAUGE_FIELD_HPP

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "lattice.hpp"
#include "color_spinor.hpp"
#include "ansi_io.hpp"
#include "linalg.hpp"
#include "su3maximization.hpp"

inline int get_fieldstrength_index(int x, int mu, int nu, int c3x3) {
    return x*NCOLOR*NCOLOR*4*4 + mu*NCOLOR*NCOLOR*4 + nu*NCOLOR*NCOLOR + c3x3;
}

class GaugeField : public Lattice {
private:
    void init();
public:
    GaugeField(const Lattice& lat) : Lattice(lat) {init();}
    GaugeField(string fname);
    GaugeField(int argc, char* argv[]) : Lattice(argc, argv) {init();}
    ~GaugeField() {delete[] Links; delete[] FieldStrength;}
    string ConfigFileName = "NoConfiguration"; // Default configuration name, meaning no input file was used

    //TODO: remember to re-calculate the field strength when the gauge links are modified
    //TODO: we need to be mindful of the order of applying static projection, stout smearing, and field strength calculation
    //1) Apply static projection
    //2) Apply stout smearing
    //3) Compute field strengths
    //4) Use the gauge field and field strengths for fermion operations

    t_complex* Links;
    int        LinksSize;
    t_complex* Link(int x, uint mu) {return Links + NCOLOR2*(x*4 + mu);};
    t_complex* Link(int t, int xs, uint mu) {return Links + NCOLOR2*((t*vol3D + xs)*4 + mu);};
    t_complex  global_Z3_phase = t_complex(1.0, 0.0); // Global Z3 phase factor, obtained from volume-averaged Polyakov loop
    t_complex* FieldStrength = nullptr; // Field strength tensor array

    void       Plaquette(t_complex* res, int x, uint mu, uint nu);
    void       StaplesSum(t_complex* res, int x, int mu, int nu, double rw=1.0); //Sum of staples at the link (x, mu) in direction nu
    void       StaplesSum(t_complex* res, int x, int mu, double* rw); //Sum of staples at the link (x, mu) with different weights for each direction
    double     ReTrPlaquette(int x, uint mu, uint nu);
    double     MeanSpatialPlaquette();
    double     MeanTemporalPlaquette();
    double     MeanPlaquette();
    void       PolyakovLine(t_complex* res, int xs, int t0, int t1); //Parallel transporter in time direction from t0 to t1 at spatial position xs
    void       PolyakovLoop(t_complex* res, int xs);
    t_complex  MeanPolyakovLoop();
    void       WriteToFile(string fname);
    void       StaticGauge();
    void       LinkBasedStaticProjection(bool polyakov_lines, bool diagnostic_output=true);
    void       PlaqBasedStaticProjection(bool polyakov_lines, bool diagnostic_output=true, int ntrials=100);
    void       StaticProjection(int ProjectionMode=0, bool diagnostic_output=true);
    void       StoutSmearing(int nsteps, double rho_s, double rho_t);
    void       StaticTests();
    void       RandomizeLinks();
    void       SetLinks(const t_complex* A);
    void       CloverLeaf(t_complex* Q_mu_nu, const int x, const int mu, const int nu);
    void       FieldStrengthCompute(t_complex* FieldStrength);

};

/*std::ostream& operator<<(std::ostream& os, const GaugeField& U) 
{
    for (int x=0; x<U.lat->vol; x++) 
        for (uint mu=0; mu<4; mu++) //TODO: Link is a pointer, print out the content - create a function for this
            os << "\t x=" << x << ", mu=" << mu << ": " << ansi::cyan << U.Link(x, mu) << ansi::reset << std::endl;
    return os;
}*/

#endif // _GAUGE_FIELD_HPP
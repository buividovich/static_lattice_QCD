#include "gauge_field.hpp"

void GaugeField::init()
{
    LinksSize = vol*4*NCOLOR2;
    Links = new t_complex[LinksSize];
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
            c3x3_set_unity(Link(x, mu));
    int FieldStrengthSize = vol*16*NCOLOR2;
    FieldStrength = new t_complex[FieldStrengthSize];
    fill(FieldStrength, FieldStrength + FieldStrengthSize, t_complex(0.0,0.0));
};

void GaugeField::Plaquette(t_complex* res, int x, uint mu, uint nu)
{
    t_complex U1[NCOLOR2], U2[NCOLOR2];
    c3x3_times_c3x3(U1, Link(x, mu), Link(shift_fwd(x, mu), nu));
    c3x3_times_c3x3(U2, Link(x, nu), Link(shift_fwd(x, nu), mu));
    c3x3_times_c3x3_conj(res, U1, U2);
}

void GaugeField::StaplesSum(t_complex* res, int x, int mu, int nu, double rw)
{
    t_complex tmp[NCOLOR2], U[NCOLOR2];

    c3x3_times_c3x3(tmp, Link(x, nu), Link(shift_fwd(x, nu), mu));
    c3x3_times_c3x3_conj(U, tmp, Link(shift_fwd(x, mu), nu));
    A_pluseq_bB(res, t_complex(rw, 0.0), U, NCOLOR2);

    int xdn = shift_bwd(x, nu);

    c3x3_conj_times_c3x3(tmp, Link(xdn, nu), Link(xdn, mu));
    c3x3_times_c3x3(U, tmp, Link(shift_fwd(xdn, mu), nu));
    
    A_pluseq_bB(res, t_complex(rw, 0.0), U, NCOLOR2);
}

void GaugeField::StaplesSum(t_complex* res, int x, int mu, double* rw)
{
    std::fill(res, res+NCOLOR2, t_complex(0.0, 0.0));
    uint nu1 = (rw[0]==0.0? 1 : 0);
    for(uint nu=0; nu<4; nu++)
        if(nu!=mu) StaplesSum(res, x, mu, nu, rw[nu]);
}

double GaugeField::ReTrPlaquette(int x, uint mu, uint nu)
{
    t_complex res[NCOLOR2];
    Plaquette(res, x, mu, nu);
    return real(tr_c3x3(res));
}

double GaugeField::MeanSpatialPlaquette()
{
    double res = 0.0;
    for(int x=0; x<vol; x++)
        for(uint i=1; i<=3; i++)
            for(uint j=i+1; j<=3; j++)
                res += ReTrPlaquette(x, i, j);
    return res/(3.0*vol);
}

double GaugeField::MeanTemporalPlaquette()
{
    double res = 0.0;
    for(int x=0; x<vol; x++)
        for(uint i=1; i<=3; i++)
                res += ReTrPlaquette(x, i, 0);
    return res/(3.0*vol);
}

double GaugeField::MeanPlaquette()
{
    return 0.5*(MeanSpatialPlaquette() + MeanTemporalPlaquette());
}

void GaugeField::PolyakovLine(t_complex* res, int xs, int t0, int t1)
{
    c3x3_set_unity(res);
    for(int t=t0; t<t1; t++)
    {
        int x = (t%LT)*vol3D + xs;
        t_complex nP[NCOLOR2];
        c3x3_times_c3x3(nP, res, Link(x, 0));
        std::copy(nP, nP+NCOLOR2, res);
    };
}

void GaugeField::PolyakovLoop(t_complex* res, int xs)
{
    t_complex nP[NCOLOR2];
    c3x3_set_unity(res);
    for(int t=0; t<LT; t++)
    {
        int x = t*vol3D + xs;
        c3x3_times_c3x3(nP, res, Link(x, 0));
        std::copy(nP, nP+NCOLOR2, res);
    };
}

t_complex GaugeField::MeanPolyakovLoop()
{
    t_complex res = 0.0;
    for(int xs=0; xs<vol3D; xs++)
    {
        t_complex P[NCOLOR2];
        PolyakovLoop(P, xs);
        res += tr_c3x3(P)/(double)(NCOLOR);
    };
    return res/(double)(vol3D);
}

GaugeField::GaugeField(string fname) : Lattice(fname)
{
    // Getting just the configuration file name
    size_t pos = fname.find_last_of("/\\");
    ConfigFileName = (pos == std::string::npos) ? fname : fname.substr(pos + 1);

    std::ifstream fin(fname, ios::binary);
    int lat_size[4];
    double plaq0;
    fin.read((char*)lat_size, sizeof(lat_size));
    fin.read((char*)&plaq0, sizeof(plaq0));
    cout << "lat_size: " << lat_size[0] << " " << lat_size[1] << " " << lat_size[2] << " " << lat_size[3] << endl;

    init();

    for(int x=0; x<vol; x++)
    {
        int xc[4];
        get_coords(x, xc);
        if((xc[0] + xc[1] + xc[2] + xc[3])%2==1) //This is a definition of an odd lattice site
            for(uint mu=0; mu<4; mu++)
            {
                fin.read((char*)Link(               x, mu), sizeof(t_complex)*9);
                fin.read((char*)Link(shift_bwd(x, mu), mu), sizeof(t_complex)*9);
            };
    };
    fin.close();
        
    double max_unitarity_err  = 0.0;
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
            max_unitarity_err = max(max_unitarity_err, c3x3_unitarity_norm(Link(x, mu)));

    double plaq1 = MeanPlaquette();

    cout << setprecision(16) << "plaq1 = " << plaq1 << ", plaq0 = " << plaq0 << endl;   

    if(abs(plaq1 - plaq0)>1.0E-8 || max_unitarity_err>1.0E-6)
    {
        cerr << ansi::red << "Error reading gauge field from file " << ansi::cyan << fname << ":" << ansi::reset << endl;
        cerr << "\t" << ansi::red << "plaq0 = " << ansi::yellow << plaq0 << ansi::red << ", plaq1 = " << ansi::yellow << plaq1 << ansi::reset << endl;
        cerr << "\t" << ansi::red << "Unitarity error: " << ansi::yellow << max_unitarity_err << ansi::reset << endl;
        cerr << endl;
    }
    else
    {
        cout << ansi::green  << "Gauge field read from file " << ansi::cyan << fname << ansi::green << " successfully, ConfigFileName = " << ansi::cyan << ConfigFileName << ansi::reset << endl;
        cout << ansi::green  << "\tMean spatial and temporal plaquettes: ";
        cout << ansi::yellow << MeanSpatialPlaquette() << ansi::reset << ", ";
        cout << ansi::yellow << MeanTemporalPlaquette() << ansi::reset << endl;
        cout << ansi::green  << "\tMean Polyakov loop:                   ";
        cout << ansi::yellow << MeanPolyakovLoop() << ansi::reset << endl;
        cout << endl;
    };
    FieldStrengthCompute(FieldStrength);
}

void GaugeField::WriteToFile(string fname)
{
    std::ofstream fout(fname, ios::binary);
    if (!fout.is_open()) {
        cerr << ansi::red << "Error: Cannot open file " << ansi::cyan << fname << ansi::red << " for writing." << ansi::reset << endl;
        return;
    }
    
    // Write header: lattice size and mean plaquette
    int lat_size[4] = {LT, LS, LS, LS};
    double plaq0 = MeanPlaquette();
    fout.write((char*)lat_size, sizeof(lat_size));
    fout.write((char*)&plaq0, sizeof(plaq0));
    
    // Write gauge links in the same format as read
    for(int x=0; x<vol; x++)
    {
        int xc[4];
        get_coords(x, xc);
        if((xc[0] + xc[1] + xc[2] + xc[3])%2==1) // Only odd lattice sites
            for(uint mu=0; mu<4; mu++)
            {
                fout.write((char*)Link(               x, mu), sizeof(t_complex)*9);
                fout.write((char*)Link(shift_bwd(x, mu), mu), sizeof(t_complex)*9);
            };
    };
    
    fout.close();
    
    cout << ansi::green << "Gauge field written to file " << ansi::cyan << fname << ansi::green << " successfully." << ansi::reset << endl;
}

void GaugeField::StaticGauge()
{
    double max_gauge_err = 0.0;
    for(int xs=0; xs<vol3D; xs++)
    {
        t_complex P[NCOLOR2], tmp[NCOLOR2], evals[NCOLOR], levecs[NCOLOR2], revecs[NCOLOR2], U0A[NCOLOR2], evals_d[NCOLOR2];
        PolyakovLoop(P, xs);
        //Calculating the LT's root of P - this will be the value of the time-like links in the static gauge
        copy(P, P+NCOLOR2, tmp);
        diagonalize(tmp, evals, revecs, levecs, NCOLOR);
        fill(evals_d, evals_d+NCOLOR2, t_complex(0.0, 0.0));
        for(int i=0; i<NCOLOR; i++) evals_d[i*NCOLOR + i] = pow(evals[i], 1.0/((double)LT));
        c3x3_times_c3x3(tmp, evals_d, levecs); //tmp = D*L
        c3x3_times_c3x3(U0A, revecs, tmp); //Now U0A is the LT's root of P

        //Now gauge transform all time-like links at (t, xs) so that they become U0A
        //Space-like links are also transformed accordingly
        t_complex omega[NCOLOR2];
        c3x3_set_unity(omega); //We assume that the gauge transformation at t=0 is unity
        for(int t=0; t<LT; t++)
        {
            //Apply the current gauge transformation to all links at (t, xs)
            for(int i=1; i<=3; i++)
            {
                //Gauge-transforming the forward link attached to site xs
                std::copy(Link(t, xs, i), Link(t, xs, i)+NCOLOR2, tmp); //tmp = Link(t, xs, i)
                c3x3_times_c3x3(Link(t, xs, i), omega, tmp); //Link(t, xs, i) -> omega * Link(t, xs, i)
                int xsb = shift_bwd(t*vol3D + xs, i);
                //Gauge-transforming the backward link attached to site xs
                std::copy(Link(xsb, i), Link(xsb, i)+NCOLOR2, tmp); //tmp = Link(xsb, i)
                c3x3_times_c3x3_conj(Link(xsb, i), tmp, omega); //Link(xsb, i) -> Link(xsb, i) * omega^+
            };

            //Compute the gauge transformation matrix for the next time slice
            c3x3_times_c3x3(tmp, omega, Link(t, xs, 0)); //Link(t, xs, 0) is the time-like link at (t, xs)
            c3x3_conj_times_c3x3(omega, U0A, tmp);  //omega = U0A^+ * tmp
            std::copy(U0A, U0A+NCOLOR2, Link(t, xs, 0)); //Set the time-like link to U0A
        };
        c3x3_set_unity(tmp);
        double err = norm_diff(omega, tmp, NCOLOR2);
        max_gauge_err = max(max_gauge_err, err);
    };
    cout << ansi::green << "Max error of omega after static gauge fixing = " << ansi::magenta << max_gauge_err << ansi::reset << endl;
}

void GaugeField::LinkBasedStaticProjection(bool polyakov_lines, bool diagnostic_output)
{
    double max_det_err = 0.0;
    double max_unitarity_err = 0.0;
    double max_projection_err = 0.0;
    double min_overlap = 1.0E10, max_overlap = -1.0E10, mean_overlap = 0.0;
    int max_iter = 0; double mean_iter = 0.0;
    
    // For each spatial point and spatial direction, average over time
    for(int xs=0; xs<vol3D; xs++)
        for(int i=1; i<=3; i++) // Spatial directions only
        {
            t_complex avg_link[NCOLOR2], U[NCOLOR2], P1up[NCOLOR2], P2up[NCOLOR2], P1dn[NCOLOR2], P2dn[NCOLOR2], tmp[NCOLOR2], dlink[NCOLOR2];
            c3x3_set_unity(P1up); c3x3_set_unity(P2up); c3x3_set_unity(P1dn); c3x3_set_unity(P2dn);
            std::copy(Link(0, xs, i), Link(0, xs, i)+NCOLOR2, avg_link);

            // Sum all time slices
            int xup1 = xs, xup2 = shift_fwd(xs, i), xdn1 = xs, xdn2 = shift_fwd(xs, i);
            do
            {
                xdn1 = shift_bwd(xdn1, 0); xdn2 = shift_bwd(xdn2, 0);
                if(polyakov_lines)
                {
                    c3x3_times_c3x3(     tmp,           P1up, Link(xup1, 0)); std::copy(tmp, tmp+NCOLOR2, P1up); //P1up *= U0(xup1)
                    c3x3_times_c3x3(     tmp,           P2up, Link(xup2, 0)); std::copy(tmp, tmp+NCOLOR2, P2up); //P2up *= U0(xup2)
                    c3x3_times_c3x3_conj(tmp,           P1dn, Link(xdn1, 0)); std::copy(tmp, tmp+NCOLOR2, P1dn); //P1dn *= U0^+(xdn1)
                    c3x3_times_c3x3_conj(tmp,           P2dn, Link(xdn2, 0)); std::copy(tmp, tmp+NCOLOR2, P2dn); //P2dn *= U0^+(xdn2)
                }
                xup1 = shift_fwd(xup1, 0); xup2 = shift_fwd(xup2, 0);

                t_complex w = (xup1==xdn1) ? t_complex(0.5, 0.0) : t_complex(1.0, 0.0);

                c3x3_times_c3x3(        tmp, P1up, Link(xup1, i)); //tmp = P1up * Link(xup1, i)
                c3x3_times_c3x3_conj( dlink,  tmp,          P2up); //dlink = tmp * P2up^+
                A_pluseq_bB(avg_link, w, dlink, NCOLOR2);

                c3x3_times_c3x3(        tmp, P1dn, Link(xdn1, i)); //tmp = P1dn * Link(xdn1, i)
                c3x3_times_c3x3_conj( dlink,  tmp,          P2dn); //dlink = tmp * P2dn^+
                A_pluseq_bB(avg_link, w, dlink, NCOLOR2);
            }
            while(xup1!=xdn1);

            rescale(avg_link, t_complex(1.0/((double)LT), 0.0), NCOLOR2);

            // Project back to SU(3)
            int iter = MaximizeSU3Overlap(avg_link, U);
            if(iter<0) { cerr << ansi::red << "Error: MaximizeSU3Overlap failed at xs=" << xs << ", i=" << i << ", error code " << iter << ansi::reset << endl; continue;};

            if(diagnostic_output)
            {
                max_iter = max(max_iter, iter);
                mean_iter += (double)iter;

                //Check the projection quality
                t_complex tmp1[NCOLOR2], tmp2[NCOLOR2];
                c3x3_conj_times_c3x3(tmp1, U, avg_link);
                c3x3_conj_times_c3x3(tmp2, avg_link, U);
                t_complex diag = (tr_c3x3(tmp1) - tr_c3x3(tmp2))/3.0;
                double projection_err  = 0.0;
                for(int i=0; i<NCOLOR; i++)
                    for(int j=0; j<NCOLOR; j++)
                    {
                        t_complex d = tmp1[i*NCOLOR + j] - tmp2[i*NCOLOR + j];
                        if(i==j) d -= diag;
                        projection_err += norm(d);
                    }
                max_projection_err = max(max_projection_err, sqrt(projection_err));

                //And finally the actual value of the overlap
                double overlap = real(tr_c3x3(tmp1))/3.0;
                min_overlap = min(min_overlap, overlap);
                max_overlap = max(max_overlap, overlap);
                mean_overlap += overlap;

                t_complex det = c3x3_det(U);
                max_det_err = max(max_det_err, abs(1.0 - det));
                max_unitarity_err = max(max_unitarity_err, c3x3_unitarity_norm(U));
            };

            // Replace all time slices with the projected averaged link
            for(int t=0; t<LT; t++) copy(U, U+NCOLOR2, Link(t, xs, i));
        };

    if(diagnostic_output)
    {
        mean_overlap /= (double)(vol3D*3);
        mean_iter /= (double)(vol3D*3);

        std::cout << ansi::green << "Max. determinant error =  " << ansi::magenta << max_det_err << ansi::reset << endl;
        std::cout << ansi::green << "Max. unitarity error   =  " << ansi::magenta << max_unitarity_err << ansi::reset << endl;
        std::cout << ansi::green << "Max. projection error  =  " << ansi::magenta << max_projection_err << ansi::reset << endl;
        std::cout << ansi::green << "Range of overlap values:  " << ansi::magenta << min_overlap << " to " << max_overlap << " (mean = " << mean_overlap << ")" << ansi::reset << endl;
        std::cout << ansi::green << "Max. iterations:          " << ansi::magenta << max_iter << ansi::reset << endl;
        std::cout << ansi::green << "Mean iterations:          " << ansi::magenta << mean_iter << ansi::reset << endl;
    };
}

t_complex* S(t_complex* aS, int xs, int i, int j)
{
    return aS + NCOLOR2*(9*xs + (3*(i-1) + (j-1)));
}

void GaugeField::PlaqBasedStaticProjection(bool polyakov_lines, bool diagnostic_output, int ntrials)
{
    t_complex* aS = new t_complex[vol3D*9*NCOLOR2]; //S[NCOLOR2*9*xs + NCOLOR2*(3*i + j)] is for the spatial plaquette averages

    for(int xs=0; xs<vol3D; xs++)
        for(int i=1; i<=3; i++)
            for(int j=i+1; j<=3; j++)
            {
                t_complex* avg_plaq = S(aS, xs, i, j);
                Plaquette(avg_plaq, xs, i, j);

                t_complex Pup[NCOLOR2], Pdn[NCOLOR2];
                c3x3_set_unity(Pup); c3x3_set_unity(Pdn);

                int xup = xs, xdn = xs;
                do
                {
                    t_complex tmp[NCOLOR2], plaq[NCOLOR2];

                    xdn = shift_bwd(xdn, 0);
                    if(polyakov_lines)
                    {
                        c3x3_times_c3x3(     tmp,           Pup, Link(xup, 0)); std::copy(tmp, tmp+NCOLOR2, Pup); //P1up *= U0(xup1)
                        c3x3_times_c3x3_conj(tmp,           Pdn, Link(xdn, 0)); std::copy(tmp, tmp+NCOLOR2, Pdn); //P1dn *= U0^+(xdn1)
                    }
                    xup = shift_fwd(xup, 0);

                    t_complex w = (xup==xdn) ? t_complex(0.5, 0.0) : t_complex(1.0, 0.0);

                    Plaquette(plaq, xup, i, j);
                    c3x3_times_c3x3(        tmp,  Pup,  plaq); //tmp = Pup * plaq
                    c3x3_times_c3x3_conj(  plaq,  tmp,   Pup); //plaq = (Pup*plaq) * Pup^+
                    A_pluseq_bB(avg_plaq, w, plaq, NCOLOR2);

                    Plaquette(plaq, xdn, i, j);
                    c3x3_times_c3x3(        tmp,  Pdn, plaq); //tmp = Pdn * plaq
                    c3x3_times_c3x3_conj(  plaq,  tmp,  Pdn); //plaq = (Pdn * plaq) * Pdn^+
                    A_pluseq_bB(avg_plaq, w, plaq, NCOLOR2);

                } while (xup!=xdn);

                rescale(avg_plaq, t_complex(1.0/((double)LT), 0.0), NCOLOR2);

                // Conjugate plaquette - we need both orientations
                c3x3_conj(S(aS, xs, j, i), avg_plaq);
            };

    if(diagnostic_output)
    {
        // To double-check the correctness of the plaquette averages, one can compute the mean plaquette from S and compare with MeanSpatialPlaquette()
        double mean_plaq_check = 0.0;
        for(int xs=0; xs<vol3D; xs++)
            for(int i=1; i<=3; i++)
                for(int j=i+1; j<=3; j++)
                    mean_plaq_check += real(tr_c3x3(S(aS, xs, i, j)));
        
        mean_plaq_check /= (3.0*vol3D);

        cout << ansi::green << "Mean spatial plaquette from S:              " << ansi::magenta << mean_plaq_check << ansi::reset << endl;
        cout << ansi::green << "Mean spatial plaquette from link variables: " << ansi::magenta << MeanSpatialPlaquette() << ansi::reset << endl;
    };

    LinkBasedStaticProjection(polyakov_lines, false);
    //RandomizeLinks();

    // Now the optimization loop - link by link, optimize the overlap of spatial plaquettes with S
    // We use the t=0 slice to represent static links, and copy them to other time slices later
    ranlux48 rng;
    std::uniform_int_distribution<int> dist_xs(0, vol3D - 1);
    std::uniform_int_distribution<int> dist_k(1, 3);

    std::ofstream debug_file;
    if(diagnostic_output)  debug_file.open("plaq_based_projection_debug.dat");

    for(int itrial=0; itrial<ntrials*vol3D*3; itrial++)
    {
        int ys = dist_xs(rng);
        int k  = dist_k(rng);

        t_complex AB[NCOLOR2]; std::fill(AB, AB+NCOLOR2, t_complex(0.0, 0.0));

        for(int i=1; i<=3; i++)
        {
            if(i==k) continue;
            t_complex tmp[NCOLOR2], staple[NCOLOR2];

            // Staple in the positive i direction
            c3x3_times_c3x3_conj( staple, Link(shift_fwd(ys, k), i),  Link(shift_fwd(ys, i), k)); //tmp = Link(ys + k, i) * Link(ys+i, k)^+
            c3x3_times_c3x3_conj(    tmp,                    staple,                Link(ys, i)); //staple = tmp * Link(ys, i)^+
            c3x3_times_c3x3(      staple,                       tmp,            S(aS, ys, i, k)); //staple = Staple * Plaquette
            A_pluseq_bB(AB, t_complex(1.0, 0.0), staple, NCOLOR2);

            // Staple in the negative i direction
            c3x3_conj_times_c3x3_conj( staple, Link(shift_fwd(shift_bwd(ys, i), k), i),  Link(shift_bwd(ys, i), k)     ); //tmp = Link(ys - i +k, i)^+ * Link(ys - i, k)^+
            c3x3_times_c3x3(              tmp,                                  staple,  S(aS, shift_bwd(ys, i), k, i) ); //staple = Link(ys - i, k)^+ * Link(ys, k)
            c3x3_times_c3x3(           staple,                                     tmp,  Link(shift_bwd(ys, i), i)     ); //staple = tmp * Link(ys - i, k)
            A_pluseq_bB(AB, t_complex(1.0, 0.0), staple, NCOLOR2);
        };

        t_complex ABc[NCOLOR2];
        c3x3_conj(ABc, AB);

        // Project avg_staple to SU(3) and set as new link
        int err = MaximizeSU3Overlap(ABc, Link(0, ys, k));
        if(err<0) { cerr << ansi::red << "Error: MaximizeSU3Overlap failed during plaquette-based projection at xs=" << ys << ", k=" << k << ", error code " << err << ansi::reset << endl; continue;};

        if(diagnostic_output && itrial%(vol3D*3)==0)
        {
            //Explicitly compute the overlap value to monitor progress
            double overlap = 0.0, mean_plaq = 0.0;
            t_complex tmp[NCOLOR2], plaq[NCOLOR2];
            for(int xs=0; xs<vol3D; xs++)
                for(int i=1; i<=3; i++)
                    for(int j=i+1; j<=3; j++)
                    {
                        Plaquette(plaq, xs, i, j);
                        mean_plaq += real(tr_c3x3(plaq));
                        c3x3_times_c3x3_conj(tmp, plaq, S(aS, xs, i, j));
                        overlap += real(tr_c3x3(tmp))/(double)(NCOLOR);
                    };

            overlap /= (vol3D * 3.0);
            mean_plaq /= (3.0*vol3D);
            std::cout << ansi::green << "Overlap after trial        " << itrial << ": " << ansi::magenta << overlap << ansi::reset << std::endl;
            std::cout << ansi::green << "Mean plaquette after trial " << itrial << ": " << ansi::magenta << mean_plaq << ansi::reset << std::endl;
            std::cout << std::endl;

            debug_file << itrial/(vol3D*3) << "\t" << overlap << "\t" << mean_plaq << std::endl;
        }
    };

    if(diagnostic_output) debug_file.close();

    for(int xs=0; xs<vol3D; xs++)
        for(int i=1; i<=3; i++) // Spatial directions only
            for(int t=1; t<LT; t++) 
                std::copy(Link(0, xs, i), Link(0, xs, i) + NCOLOR2, Link(t, xs, i));

    delete[] aS;
}

void GaugeField::StaticProjection(int ProjectionMode, bool diagnostic_output)
{
    // First apply static gauge fixing
    if(ProjectionMode==0 || ProjectionMode==1 || ProjectionMode==3)
        StaticGauge();

    //Determine the Z3 center sector and replace the time-like links at the last time slice with the corresponding center element
    t_complex PolyakovLoopAvg = MeanPolyakovLoop();
    double angle = std::arg(PolyakovLoopAvg);
    angle = (abs(angle)<pi/3 ? 0 : 2*pi/3*angle/abs(angle)); //angle/abs(angle) is just the sign of angle
    global_Z3_phase = std::exp(t_complex(0.0, angle));

    std::cout << ansi::green << "Polyakov loop average:    " << ansi::magenta << PolyakovLoopAvg << ansi::reset << endl;
    std::cout << ansi::green << "Center sector projection: " << ansi::magenta << global_Z3_phase << ansi::reset << endl;

    if(ProjectionMode==0)
    {
        //All links are set to be equal to the links at t=0
        for(int it=1; it<LT; it++)
            for(int xs=0; xs<vol3D; xs++)
                for(int i=1; i<=3; i++) // Spatial directions only
                    std::copy(Link(0, xs, i), Link(0, xs, i)+NCOLOR2, Link(it, xs, i));
    };

    if(ProjectionMode==1) LinkBasedStaticProjection(false, diagnostic_output);
    if(ProjectionMode==2) LinkBasedStaticProjection(true,  diagnostic_output);
    if(ProjectionMode==3) PlaqBasedStaticProjection(false, diagnostic_output);
    if(ProjectionMode==4) PlaqBasedStaticProjection(true,  diagnostic_output);

    //Finally, set the time-like links at the last time slice to the global Z3 phase
    for(int xs=0; xs<vol3D; xs++)
        for(int it=0; it<LT; it++)
        {
            t_complex* link = Link(it, xs, 0);
            std::fill(link, link+NCOLOR2, t_complex(0.0, 0.0));
            for(int i=0; i<NCOLOR; i++)
                link[i*NCOLOR + i] = (it==LT-1? global_Z3_phase : t_complex(1.0, 0.0));
        };
    
    // Recalculate Field Strength
    FieldStrengthCompute(FieldStrength);
}

void GaugeField::StaticTests()
{
    StaticGauge();
    int LT2 = LT/2+1;
    double* GES = new double[LT2]; std::fill(GES, GES+LT2, 0.0);
    double* GEC = new double[LT2]; std::fill(GEC, GEC+LT2, 0.0);
    double* GPS = new double[LT2]; std::fill(GPS, GPS+LT2, 0.0);
    double* GPC = new double[LT2]; std::fill(GPC, GPC+LT2, 0.0);
    for(int xs=0; xs<vol3D; xs++)
        for(int t1=0; t1<LT; t1++)
            for(int dts=0; dts<LT2; dts++)
            {
                t_complex tmp[NCOLOR2], P[NCOLOR2];
                int t2 = (t1 + dts)%LT;
                PolyakovLine(P, xs, t1, t1 + dts);
                for(int i=1; i<=3; i++)
                {
                    //Link-link correlator without Polyakov lines
                    c3x3_conj_times_c3x3(tmp, Link(t1, xs, i), Link(t2, xs, i));
                    GES[dts] += real(tr_c3x3(tmp))/3.0;
                    //Link-link correlator with Polyakov lines
                    t_complex P_fwd[NCOLOR2], tmp1[NCOLOR2], res[NCOLOR2];
                    int xs_fwd = shift_fwd(xs, i);
                    PolyakovLine(P_fwd, xs_fwd, t1, t1 + dts);
                    c3x3_times_c3x3(tmp, Link(t1, xs, i), P_fwd);
                    c3x3_times_c3x3(tmp1, P, Link(t2, xs, i));
                    c3x3_conj_times_c3x3(res, tmp, tmp1);
                    GEC[dts] += real(tr_c3x3(res))/3.0;
                    for(int j=i+1; j<=3; j++)
                    {
                        //Plaquette-plaquette correlator without Polyakov lines
                        t_complex Plaq1[NCOLOR2], Plaq2[NCOLOR2];
                        Plaquette(Plaq1, t1*vol3D + xs, i, j);
                        Plaquette(Plaq2, t2*vol3D + xs, i, j);
                        c3x3_conj_times_c3x3(tmp, Plaq1, Plaq2);
                        GPS[dts] += real(tr_c3x3(tmp))/3.0;
                        //Plaquette-plaquette correlator with Polyakov lines
                        c3x3_times_c3x3(tmp, Plaq1, P);
                        c3x3_times_c3x3_conj(tmp1, tmp, Plaq2);
                        c3x3_times_c3x3_conj(res, tmp1, P);
                        GPC[dts] += real(tr_c3x3(res))/3.0;
                    }
                }
            }

    for(int dts=0; dts<LT2; dts++)
    {
        GES[dts] /= (double)(vol*3);
        GEC[dts] /= (double)(vol*3);
        GPS[dts] /= (double)(vol*3);
        GPC[dts] /= (double)(vol*3);
    };
    std::ofstream fout("./data/" + ConfigFileName + ".ST");
    if(!fout.is_open()){ cerr << ansi::red << "Error: Cannot open file ./data/static_tests.dat for writing." << ansi::reset << endl; return; };
    fout << std::setprecision(16);
    for(int i = 0; i < LT2; ++i)
        fout << GES[i] << " " << GEC[i] << " " << GPS[i] << " " << GPC[i] << "\n";
    fout.close();
    cout << ansi::green << "Wrote static tests to " << ansi::cyan << "./data/static_tests.dat" << ansi::reset << endl;
    delete[] GES;
    delete[] GEC;
    delete[] GPS;
    delete[] GPC;
}

void GaugeField::RandomizeLinks()
{
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
        {
            t_complex random_matrix[NCOLOR2];
            for(int i=0; i<NCOLOR2; i++)
                random_matrix[i] = t_complex(drand48(), drand48());
            int iter = MaximizeSU3Overlap(random_matrix, Link(x, mu));
            if(iter<0) { cerr << ansi::red << "Error: MaximizeSU3Overlap failed at x=" << x << ", mu=" << mu << ", error code " << iter << ansi::reset << endl; continue; }; 
        };
}

void GaugeField::SetLinks(const t_complex* A)
{
    for (int x = 0; x < vol; x++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int i = 0; i < NCOLOR2; i++) {
                Link(x, mu)[i] = A[i];
            }
        }
    }
}

void GaugeField::StoutSmearing(int nsteps, double rho_s, double rho_t)
{
    t_complex* new_links = new t_complex[vol * 4 * NCOLOR2];
    int mu1 = (rho_t==0.0? 1 : 0);

    for (int step = 0; step < nsteps; step++)
    { 
        for (int x = 0; x < vol; x++) 
            for (int mu = 0; mu<4; mu++) 
            {
                double rw[4] = {rho_t, rho_s, rho_s, rho_s};
                if(mu==0) {rw[0] = 0.0; rw[1] = rho_t; rw[2] = rho_t; rw[3] = rho_t; };

                t_complex staple[NCOLOR2], omega[NCOLOR2], tmp[NCOLOR2];
                StaplesSum(staple, x, mu, rw);
                c3x3_times_c3x3_conj(omega, staple, Link(x, mu));

                for(int i=0; i<NCOLOR; i++)
                    for(int j=i+1; j<NCOLOR; j++)
                    {
                        t_complex val = t_complex(0.0, 0.5)*(omega[i*NCOLOR + j] - std::conj(omega[j*NCOLOR + i]));
                        omega[i*NCOLOR + j] = val;
                        omega[j*NCOLOR + i] = std::conj(val);
                    }
                for(int i=0; i<NCOLOR; i++)
                    omega[i*NCOLOR + i] = t_complex(0.0, 0.5)*(omega[i*NCOLOR + i] - std::conj(omega[i*NCOLOR + i]));
                t_complex tr_omega = tr_c3x3(omega)/(double)(NCOLOR);
                for(int i=0; i<NCOLOR; i++)
                    omega[i*NCOLOR + i] -= tr_omega;

                matrix_exponential(tmp, omega, t_complex(0.0, -1.0), NCOLOR); //tmp = exp(omega)
                c3x3_times_c3x3(new_links + NCOLOR2*(x*4 + mu), tmp, Link(x, mu));
            }

        std::copy(new_links, new_links + vol * 4 * NCOLOR2, Links);
    }
    delete[] new_links;

    // Recalculate Field Strength
    FieldStrengthCompute(FieldStrength);
}

void GaugeField::CloverLeaf(t_complex* Q_mu_nu, const int x, const int mu, const int nu)
{
    t_complex tmp1[NCOLOR*NCOLOR];
    t_complex tmp2[NCOLOR*NCOLOR];

    for(int i = 0; i < NCOLOR*NCOLOR; i++) {
        Q_mu_nu[i] = 0;
    }
    
    // First plaquette: U(x,mu) U(x+mu,nu) U^dag(x+nu,mu) U^dag(x,nu)
    int x_mu_fwd = shift_fwd(x, mu);
    int x_nu_fwd = shift_fwd(x, nu);
    c3x3_times_c3x3(tmp1, Link(x, mu), Link(x_mu_fwd, nu));
    c3x3_times_c3x3_conj(tmp2, tmp1, Link(x_nu_fwd, mu));
    c3x3_times_c3x3_conj(tmp1, tmp2, Link(x, nu));
    for(int i = 0; i < NCOLOR*NCOLOR; i++) {
        Q_mu_nu[i] += tmp1[i];
    }
    
    // Second plaquette: U(x,nu) U^dag(x-mu+nu,mu) U^dag(x-mu,nu) U(x-mu,mu)
    int x_mu_bwd = shift_bwd(x, mu);
    int x_mu_bwd_nu_fwd = shift_fwd(x_mu_bwd, nu);
    c3x3_times_c3x3_conj(tmp1, Link(x, nu), Link(x_mu_bwd_nu_fwd, mu));
    c3x3_times_c3x3_conj(tmp2, tmp1, Link(x_mu_bwd, nu));
    c3x3_times_c3x3(tmp1, tmp2, Link(x_mu_bwd, mu));
    for(int i = 0; i < NCOLOR*NCOLOR; i++) {
        Q_mu_nu[i] += tmp1[i];
    }

    // Third plaquette: U^dag(x-mu,mu) U^dag(x-mu-nu,nu) U(x-mu-nu,mu) U(x-nu,nu)
    int x_nu_bwd = shift_bwd(x, nu);
    int x_mu_bwd_nu_bwd = shift_bwd(x_mu_bwd, nu);
    c3x3_conj_times_c3x3_conj(tmp1, Link(x_mu_bwd, mu), Link(x_mu_bwd_nu_bwd, nu));
    c3x3_times_c3x3(tmp2, tmp1, Link(x_mu_bwd_nu_bwd, mu));
    c3x3_times_c3x3(tmp1, tmp2, Link(x_nu_bwd, nu));
    for(int i = 0; i < NCOLOR*NCOLOR; i++) {
        Q_mu_nu[i] += tmp1[i];
    }
    
    // Fourth plaquette: U^dag(x-nu,nu) U(x-nu,mu) U(x+mu-nu,nu) U^dag(x,mu)
    int x_mu_fwd_nu_bwd = shift_fwd(x_nu_bwd, mu);
    c3x3_conj_times_c3x3(tmp1, Link(x_nu_bwd, nu), Link(x_nu_bwd, mu));
    c3x3_times_c3x3(tmp2, tmp1, Link(x_mu_fwd_nu_bwd, nu));
    c3x3_times_c3x3_conj(tmp1, tmp2, Link(x, mu));
    for(int i = 0; i < NCOLOR*NCOLOR; i++) {
        Q_mu_nu[i] += tmp1[i];
    }
}

void GaugeField::FieldStrengthCompute(t_complex* FieldStrength)
{
    #pragma omp parallel for
    for (int x = 0; x < vol; x++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                t_complex Q_mu_nu[NCOLOR*NCOLOR];
                t_complex Q_nu_mu[NCOLOR*NCOLOR];
                CloverLeaf(Q_mu_nu, x, mu, nu);
                CloverLeaf(Q_nu_mu, x, nu, mu);

                // F_mu_nu = (1/8)(Q_mu_nu - Q_nu_mu)
                for (int i = 0; i < NCOLOR*NCOLOR; i++) {
                    FieldStrength[get_fieldstrength_index(x, mu, nu, i)] = 0.125 * (Q_mu_nu[i] - Q_nu_mu[i]);
                }
            }
        }
    }
}


//These checks were a part of the StaticGauge function
//Checking if the diagonalization was successful
        // t_complex check[9], evals_d[9];
        // std::fill(evals_d, evals_d+9, t_complex(0.0, 0.0));
        // for(int i=0; i<3; i++)
        //     evals_d[3*i + i] = evals[i];
        
        // c3x3_times_c3x3(tmp, evals_d, levecs); //tmp = D*L
        // c3x3_times_c3x3(check, revecs, tmp);   //check = R*D*L
        // double err = norm_diff(P, check, 9);
        // max_diag_err = max(max_diag_err, err);

        //Checking the root calculation
        // t_complex check[NCOLOR2]; c3x3_set_unity(check);
        // for(int i=0; i<LT; i++)
        // {
        //     c3x3_times_c3x3(tmp, check, U0A); //tmp = check*U0A
        //     copy(tmp, tmp+NCOLOR2, check);
        // };

        // double err = norm_diff(P, check, NCOLOR2);
        // max_root_err = max(max_root_err, err);
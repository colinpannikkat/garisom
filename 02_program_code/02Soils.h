#ifndef SOILS_H
#define SOILS_H

#include <iostream>
#include <string>
#include <cmath> // math utility functions
#include <vector>
#include "04Component.h"

using namespace std;
class soils
{
public:
    // Van Genuchten functions
    void get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[2001][101], long &rowLR, long &colLR);
};

class RhizosphereComponent : public Component {

    private:
        
        double van_gen_alpha;
        double van_gen_n;
        double thetasat;

    public:

        void setVanGenAlpha(double alpha);
        double getVanGenAlpha();

        void setVanGenN(double n);
        double getVanGenN();

        void setThetasat(double thetasat);
        double getThetaSat();

        double vg(double pressure);

        void trapzd(const double &p1, const double &p2, double &s, const int &t, int &it);
        void qtrap(double &p1, double &p2, double &s);
        void calc_flow_rate(const double &p_inc, const double &k_min);

        /*

            Calculates hydraulic conductance, via a van Genuchten function.

            k < 0.0005 * k_max is physiological zero

            Inputs:
                - Absolute value of soil water potential: P_soil
                - Rhizosphere vulnerability curve: v
                - Maximum flow rate per pressure drop: K_max

            Outputs:
                - Hydraulic conductance of component: k
        */
        double& hydraulic_conductance(double &P_soil, double &v, double &K_max);

        /*
        
            Calculates vulnerability curve for Rhizosphere. 

            k < 0.0005 * k_max is physiological zero.

            Alpha and n are texture specific parameters.

            Inputs:
                - n
                - alpha

            Outputs:
                - Vulnerability (conductance)
        
        */
        double& rhizo_v(double &n, double &alpha);

        /*
        
            Calculates k_max.

            "We solved for the rhizosphere kmax from an inputed average percent 
            rhizosphere resistance. The percent of continuum resistance was 
            calculated from in-series vulnerability curves of rhizosphere, root, 
            stem, and leaf at same P" - Sperry et al. 2016

            This resistance is averaged over 0.1 MPa increments from P = 0 to
            P_crit.

            Inputs:
                - Average percent rhizosphere resistance: avg_res
        */
       double &calc_k_max(double &avg_res);

};

#endif
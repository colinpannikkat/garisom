#ifndef SOILS_H
#define SOILS_H

#include <iostream>
#include <string>
#include <cmath> // math utility functions
#include <vector>
#include "04Component.h"

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
        double rvg(double pressure);
        double swc(double pressure);

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

/* Used in GARISOM v2.0.5 */
class soils
{
public:
    // Van Genuchten functions
    /*Get Van Genuchten alpha // override if provided*/
    // This function obtains the Van Genuchten parameters for a given texture
    void get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[2001][101], long &rowLR, long &colLR) {
        std::string a, n, soilkmax, thetasat;
        if (texture == "sand"){
            a = "1479.5945"; 
            n = "2.68";
            soilkmax = "30305.88"; 
            thetasat = "0.43";
        } else if (texture == "loamy sand"){
            a = "1265.3084";
            n = "2.28";
            soilkmax = "14897.84";
            thetasat = "0.41";
        } else if (texture == "sandy loam") {
            a = "765.3075";
            n = "1.89";
            soilkmax = "4510.168";
            thetasat = "0.41";
        } else if (texture == "loam") {
            a = "367.3476";
            n = "1.56";
            soilkmax = "1061.216";
            thetasat = "0.43";
        } else if (texture == "silt") {
            a = "163.2656";
            n = "1.37";
            soilkmax = "255.1";
            thetasat = "0.46";
        } else if (texture == "silt loam") {
            a = "204.082";
            n = "1.41";
            soilkmax = "459.18";
            thetasat = "0.45";
        } else if (texture == "sandy clay loam") {
            a = "602.0419";
            n = "1.48";
            soilkmax = "1336.724";
            thetasat = "0.39";
        } else if (texture == "clay loam") {
            a = "193.8779";
            n = "1.31";
            soilkmax = "265.304";
            thetasat = "0.41";
        } else if (texture == "silty clay loam") {
            a = "102.041";
            n = "1.23";
            soilkmax = "71.428";
            thetasat = "0.43";
        } else if (texture == "sandy clay") {
            a = "275.5107";
            n = "1.23";
            soilkmax = "122.448";
            thetasat = "0.38";
        } else if (texture == "silty clay") {
            a = "51.0205";
            n = "1.09";
            soilkmax = "20.408";
            thetasat = "0.36";
        } else if (texture == "clay") {
            a = "81.6328";
            n = "1.09";
            soilkmax = "204.08";
            thetasat = "0.38";
        } else {
            std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
            std::cout << "SOURCE: Incorrect soil texture category" << std::endl;
            std::cout << "ACTION: Model stops " << std::endl;
            std::cout << std::endl;
            abort();
        }
        for(long k=1;k<=layers;k++){
            soillayersTable[rowLR + k][colLR + 3] = a; 
            soillayersTable[rowLR + k][colLR + 4] = n;
            soillayersTable[rowLR + k][colLR + 5] = soilkmax;
            soillayersTable[rowLR + k][colLR + 6] = thetasat;
        }
    }
};

#endif
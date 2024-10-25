#ifndef COMPONENT_H
#define COMPONENT_H

#include <stdio.h>
#include <cmath>    // math utility functions
#include <string>   // the C++ String Class, easier to deal with than char arrays for this application
#include <ctime>    // timers for performance testing
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <assert.h>
#include "01Utilities.h"

/*
    This implementation follows Sperry et al. 2016

    Implemented by Colin Pannikkat
    Last updated: 10/05/24
*/

class Component {
    
    protected:

        double  b_wb,           // b for weibull
                c_wb,           // c for weibull
                res_percent,    // percent resistances at ksat
                res_wb,         // weibull resistance, 1/cond_wb
                cond_wb,        // weibull conductance, 1/res_wb           
                k_max,          // kmax, conductance per basal area
                                // ksat used interchangably
                p_crit;         // critical pressure, of VC
        std::vector<double> wb_fatigue;
        double e_p[CURVE_MAX] = {0};    // E(P) curve for component
        double k[CURVE_MAX] = {0};      // K (conductivity) curve for component
        double k_v[CURVE_MAX] = {0};    // Virgin K (conductivity) curve for component
        double k_comp[CURVE_MAX] = {0}; // New(conductivity) curve for component based on other components

    public:

        /* Getters */
        double getBwb();
        double getCwb();
        double getResPercent();
        double getResWb();
        double getCondWb();
        double getKmax();
        double& getFatigue(int index);

        /* Setters */
        void setBwb(double value);
        void setCwb(double value);
        void setResPercent(double value);
        void setResWb(double value);
        void setCondWb(double value);
        void setKmax(double value);
        void setFatigue(int index, double value);

        /* 

            Calculates hydraulic conductance of component, via a weibull function.

            Inputs:
                - Maximum flow rate per pressure drop: K_max
                - Negative sap pressure (-MPa): P
                - Curve shift: b_wb
                - Shape: c_wb
                    Exponential corresponds to c <= 1, sigmoidal has c > 1

            Outputs:
                - Hydraulic conductance of component: k

        */
        double wb(const double &pressure);

        /*
        
            Calculates steady-state flow rate, E_i(P), for a component

            E_i = \int_{P_up}^{P_down}{k(P)_i}dP

            Inputs:
                - p_inc: Increment for pressure in Reimann integration
                - k_min: Minimum conductivity amount

            Outputs:
                - Steady-state flow rate: e_p
        */ 
        void calc_flow_rate(const double &p_inc, const double &k_min);

        void trapzd(const double &p1, const double &p2, double &s, const int &t, int &it);
        void virtual qtrap(double &p1, double &p2, double &s);

        void printCurveToFile(const double &p_inc, const std::string &filename) const {
            std::ofstream outFile(filename);
            if (!outFile) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
            }

            outFile << "p_inc,E(P)" << std::endl;
            for (int i = 0; p_inc * i <= this->p_crit; ++i) {
            outFile << p_inc * i << "," << e_p[i] << std::endl;
            }

            outFile.close();
        }
};

#endif